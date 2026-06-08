import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── Parameters ───────────────────────────────────────────────
params = {
    'k_in':           0.5,
    'k_cl':           0.02,
    'k_ros':          0.5,
    'k_sod':          5.0,
    'k_gsh':          2.0,
    'k_syn':          0.02,
    'k_deg':          0.001,
    'k_ox':           0.001,
    'k_ox_sustained': 0.0,
    'k_nrf2':         0.1,
    'k_nrf2deg':      0.01,
    'k_gshsyn':       0.1,
    'k_gshdeg':       0.02,
    'k_gshsyn_basal': 5.0,
    'k_gsh_turn':     0.001,
}

Km_gsh           = 500.0
k_damage_acute   = 0.1
k_repair_acute   = 0.001
k_damage_chronic = 0.0005
k_chronic_gsh    = 0.0002   # GSH-deficit driven chronic damage rate

# ── Pre-dox steady state ─────────────────────────────────────
Keap1_ss = params['k_syn'] / params['k_deg']
Nrf2_ss  = params['k_nrf2'] / (params['k_nrf2deg'] * Keap1_ss)
GSH_ss   = (params['k_gshsyn_basal'] + params['k_gshsyn'] * Nrf2_ss) / params['k_gsh_turn']

y0 = [0.0, 0.0, 0.0, Keap1_ss, Nrf2_ss, GSH_ss, 0.0, 0.0]

# ── Dosing schedule ───────────────────────────────────────────
cycle_days        = 21
n_cycles          = 6
dose_duration_min = 60
cycle_starts_min  = [c * cycle_days * 24 * 60 for c in range(n_cycles)]
total_time_min    = n_cycles * cycle_days * 24 * 60

def k_in_pulsed(t):
    for cs in cycle_starts_min:
        if cs <= t < cs + dose_duration_min:
            return params['k_in'] * 20
    return 0.0

# ── ODE system ───────────────────────────────────────────────
def odes_pulsed(t, y, p):
    Dox, O2, H2O2, Keap1, Nrf2, GSH, D_acute, D_chronic = y
    GSH       = max(GSH,       0.0)
    Keap1     = max(Keap1,     0.0)
    Nrf2      = max(Nrf2,      0.0)
    D_acute   = max(D_acute,   0.0)
    D_chronic = max(D_chronic, 0.0)
    gsh_sat   = GSH / (Km_gsh + GSH)
    kin       = k_in_pulsed(t)

    k_ox_total = p['k_ox'] + p['k_ox_sustained']
    GSH_deficit = max(GSH_ss - GSH, 0.0)

    dDox      = kin - p['k_cl'] * Dox
    dO2       = p['k_ros'] * Dox - p['k_sod'] * O2
    dH2O2     = p['k_sod'] * O2 - p['k_gsh'] * gsh_sat * H2O2
    dKeap1    = p['k_syn'] - p['k_deg'] * Keap1 \
                - k_ox_total * H2O2 * Keap1 \
                - p['k_ox_sustained'] * Keap1
    dNrf2     = p['k_nrf2'] - p['k_nrf2deg'] * Keap1 * Nrf2
    dGSH      = (p['k_gshsyn_basal'] + p['k_gshsyn'] * Nrf2
                 - p['k_gshdeg'] * gsh_sat * H2O2
                 - p['k_gsh_turn'] * GSH)
    dD_acute  = k_damage_acute * H2O2 - k_repair_acute * D_acute
    dD_chronic = k_damage_chronic * H2O2 + k_chronic_gsh * GSH_deficit

    return [dDox, dO2, dH2O2, dKeap1, dNrf2, dGSH, dD_acute, dD_chronic]

# ── Time evaluation ───────────────────────────────────────────
t_eval_segments = []
for cs in cycle_starts_min:
    t_eval_segments.append(np.linspace(cs, cs + dose_duration_min, 200))
    t_eval_segments.append(np.linspace(cs + dose_duration_min,
                                        cs + cycle_days * 24 * 60, 400))
t_eval = np.unique(np.concatenate(t_eval_segments))
t_eval = t_eval[t_eval <= total_time_min]

# ── Solve ─────────────────────────────────────────────────────
print("Solving baseline...")
sol = solve_ivp(odes_pulsed, (0, total_time_min), y0, t_eval=t_eval,
                args=(params,), method='RK45', max_step=120, rtol=1e-5, atol=1e-8)

params_int = params.copy()
params_int['k_ox_sustained'] = 0.002

print("Solving intervention (sustained sulforaphane)...")
sol_int = solve_ivp(odes_pulsed, (0, total_time_min), y0, t_eval=t_eval,
                    args=(params_int,), method='RK45', max_step=120, rtol=1e-5, atol=1e-8)

print("Done. Plotting...")

t_days          = sol.t     / (24 * 60)
t_days_int      = sol_int.t / (24 * 60)
cycle_day_marks = [cs / (24 * 60) for cs in cycle_starts_min]
D_total         = sol.y[6]     + sol.y[7]
D_total_int     = sol_int.y[6] + sol_int.y[7]

# ── Plot 1: Full variable panel ───────────────────────────────
labels_full = ['Dox', 'O2·⁻', 'H2O2', 'Keap1', 'Nrf2', 'GSH', 'D_acute', 'D_chronic']
fig, axes = plt.subplots(3, 3, figsize=(16, 12))
axes = axes.flatten()
for i in range(8):
    ax = axes[i]
    ax.plot(t_days,     sol.y[i],     linewidth=1.2, color='steelblue', label='Dox only')
    ax.plot(t_days_int, sol_int.y[i], linewidth=1.2, color='firebrick',
            linestyle='--', label='Dox + Sulforaphane')
    for cd in cycle_day_marks:
        ax.axvline(cd, color='gray', linewidth=0.6, linestyle=':', alpha=0.5)
    ax.set_title(labels_full[i])
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('nM' if i < 6 else 'AU')
    ax.legend(fontsize=6)
    ax.grid(True, alpha=0.3)
for j, cd in enumerate(cycle_day_marks):
    axes[0].annotate(f'C{j+1}', xy=(cd, axes[0].get_ylim()[1]*0.85),
                     fontsize=7, color='gray', ha='center')
axes[8].axis('off')
plt.suptitle('6-Cycle AC Regimen: All Variables\n(Dotted lines = cycle start)', fontsize=13)
plt.tight_layout()
plt.savefig('cumulative_dosing_full.png', dpi=150)
plt.show()
print("1/2 done")

# ── Plot 2: Clinical summary ──────────────────────────────────
fig2, axes2 = plt.subplots(1, 3, figsize=(15, 5))

axes2[0].plot(t_days,     D_total,     linewidth=2, color='steelblue', label='Dox only')
axes2[0].plot(t_days_int, D_total_int, linewidth=2, color='firebrick',
              linestyle='--', label='Dox + Sulforaphane')
for cd in cycle_day_marks:
    axes2[0].axvline(cd, color='gray', linewidth=0.8, linestyle=':', alpha=0.6)
axes2[0].set_title('Total Oxidative Damage\n(Acute + Irreversible)')
axes2[0].set_xlabel('Time (days)')
axes2[0].set_ylabel('Damage (AU)')
axes2[0].legend(fontsize=9)
axes2[0].grid(True, alpha=0.3)

axes2[1].plot(t_days,     sol.y[7],     linewidth=2.5, color='steelblue', label='Dox only')
axes2[1].plot(t_days_int, sol_int.y[7], linewidth=2.5, color='firebrick',
              linestyle='--', label='Dox + Sulforaphane')
for cd in cycle_day_marks:
    axes2[1].axvline(cd, color='gray', linewidth=0.8, linestyle=':', alpha=0.6)
axes2[1].set_title('Irreversible Damage Accumulation\n(Staircase = cumulative per cycle)')
axes2[1].set_xlabel('Time (days)')
axes2[1].set_ylabel('Irreversible Damage (AU)')
axes2[1].legend(fontsize=9)
axes2[1].grid(True, alpha=0.3)

cycle_chronic_dox, cycle_chronic_int = [], []
for cs in cycle_starts_min:
    ce = cs + cycle_days * 24 * 60
    mask_d = (sol.t     >= ce - 120) & (sol.t     <= ce)
    mask_i = (sol_int.t >= ce - 120) & (sol_int.t <= ce)
    cycle_chronic_dox.append(np.mean(sol.y[7][mask_d])     if mask_d.any() else 0)
    cycle_chronic_int.append(np.mean(sol_int.y[7][mask_i]) if mask_i.any() else 0)

x = np.arange(1, n_cycles + 1)
width = 0.35
axes2[2].bar(x - width/2, cycle_chronic_dox, width, color='steelblue', label='Dox only', alpha=0.85)
axes2[2].bar(x + width/2, cycle_chronic_int, width, color='firebrick', label='Dox + Sulforaphane', alpha=0.85)
axes2[2].set_title('Irreversible Damage at End of Each Cycle')
axes2[2].set_xlabel('Cycle number')
axes2[2].set_ylabel('Irreversible Damage (AU)')
axes2[2].set_xticks(x)
axes2[2].legend(fontsize=9)
axes2[2].grid(True, alpha=0.3, axis='y')

prot_c6 = (1 - cycle_chronic_int[-1] / cycle_chronic_dox[-1]) * 100 if cycle_chronic_dox[-1] > 0 else 0
print(f"\nSulforaphane protection at end of cycle 6: {prot_c6:.1f}% reduction in irreversible damage")
print(f"Dox only       — D_chronic cycle 6: {cycle_chronic_dox[-1]:.2f} AU")
print(f"+ Sulforaphane — D_chronic cycle 6: {cycle_chronic_int[-1]:.2f} AU")

plt.suptitle('Clinical Regimen Summary: 6-Cycle AC Chemotherapy', fontsize=13)
plt.tight_layout()
plt.savefig('cumulative_dosing_summary.png', dpi=150)
plt.show()
print("2/2 done.")