import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

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
k_chronic_gsh    = 0.0002

Keap1_ss = params['k_syn'] / params['k_deg']
Nrf2_ss  = params['k_nrf2'] / (params['k_nrf2deg'] * Keap1_ss)
GSH_ss   = (params['k_gshsyn_basal'] + params['k_gshsyn'] * Nrf2_ss) / params['k_gsh_turn']
y0 = [0.0, 0.0, 0.0, Keap1_ss, Nrf2_ss, GSH_ss]

t_span = (0, 600)
t_eval = np.linspace(0, 600, 2000)

def odes(t, y, p):
    Dox, O2, H2O2, Keap1, Nrf2, GSH = y
    GSH   = max(GSH,   0.0)
    Keap1 = max(Keap1, 0.0)
    Nrf2  = max(Nrf2,  0.0)
    gsh_sat     = GSH / (Km_gsh + GSH)
    k_ox_total  = p['k_ox'] + p['k_ox_sustained']
    GSH_deficit = max(GSH_ss - GSH, 0.0)

    dDox   = p['k_in'] - p['k_cl'] * Dox
    dO2    = p['k_ros'] * Dox - p['k_sod'] * O2
    dH2O2  = p['k_sod'] * O2 - p['k_gsh'] * gsh_sat * H2O2
    dKeap1 = p['k_syn'] - p['k_deg'] * Keap1 \
             - k_ox_total * H2O2 * Keap1 \
             - p['k_ox_sustained'] * Keap1
    dNrf2  = p['k_nrf2'] - p['k_nrf2deg'] * Keap1 * Nrf2
    dGSH   = (p['k_gshsyn_basal'] + p['k_gshsyn'] * Nrf2
              - p['k_gshdeg'] * gsh_sat * H2O2
              - p['k_gsh_turn'] * GSH)
    return [dDox, dO2, dH2O2, dKeap1, dNrf2, dGSH]

# ── Three conditions ──────────────────────────────────────────
# 1. Dox only
p_dox = params.copy()
sol_dox = solve_ivp(odes, t_span, y0, t_eval=t_eval, args=(p_dox,), method='RK45')

# 2. Dox + Sulforaphane (Nrf2 activator — sustained Keap1 degradation)
p_sulf = params.copy()
p_sulf['k_ox_sustained'] = 0.002
sol_sulf = solve_ivp(odes, t_span, y0, t_eval=t_eval, args=(p_sulf,), method='RK45')

# 3. Dox + Dexrazoxane (iron chelator — reduces ROS generation at source)
# Dexrazoxane works by chelating iron that catalyzes dox redox cycling
# Modeled as reduction in k_ros (less superoxide generated per dox molecule)
# and small reduction in k_in (some evidence of reduced dox uptake)
p_dex = params.copy()
p_dex['k_ros'] = params['k_ros'] * 0.4   # 60% reduction in ROS generation
sol_dex = solve_ivp(odes, t_span, y0, t_eval=t_eval, args=(p_dex,), method='RK45')

# 4. Combination: Dexrazoxane + Sulforaphane
p_combo = params.copy()
p_combo['k_ros']          = params['k_ros'] * 0.4
p_combo['k_ox_sustained'] = 0.002
sol_combo = solve_ivp(odes, t_span, y0, t_eval=t_eval, args=(p_combo,), method='RK45')

t_hours = t_eval / 60

# ── Plot 1: Four-condition comparison ────────────────────────
conditions = [
    (sol_dox,   'Dox only',                    'steelblue',  '-'),
    (sol_sulf,  'Dox + Sulforaphane (Nrf2)',   'firebrick',  '--'),
    (sol_dex,   'Dox + Dexrazoxane (source)',  'forestgreen','-.'),
    (sol_combo, 'Dox + Both',                  'purple',     ':'),
]

panels = [
    (2, 'H₂O₂ (nM)',    'ROS Level'),
    (4, 'Nrf2 (nM)',    'Nrf2 Activation'),
    (5, 'GSH (nM)',     'Antioxidant Reserve'),
]

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for ax, (idx, ylabel, title) in zip(axes, panels):
    for sol, label, color, ls in conditions:
        ax.plot(t_hours, sol.y[idx], linewidth=2, color=color,
                linestyle=ls, label=label)
    ax.set_title(title)
    ax.set_xlabel('Time (hours)')
    ax.set_ylabel(ylabel)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

plt.suptitle('Mechanism Comparison: Sulforaphane vs Dexrazoxane vs Combination\n'
             'Upstream (source reduction) vs Downstream (antioxidant upregulation)',
             fontsize=12)
plt.tight_layout()
plt.savefig('dexrazoxane_comparison.png', dpi=150)
print("Saved: dexrazoxane_comparison.png")

# ── Compute damage for all four conditions ───────────────────
def compute_damage(sol):
    H2O2 = sol.y[2]
    GSH  = sol.y[5]
    GSH_deficit = np.maximum(GSH_ss - GSH, 0)
    dt = np.diff(t_eval)
    d_acute_rate   = k_damage_acute  * H2O2
    d_chronic_rate = k_damage_chronic * H2O2 + k_chronic_gsh * GSH_deficit
    D_chronic = np.cumsum(d_chronic_rate[:-1] * dt)
    return D_chronic

dmg_dox   = compute_damage(sol_dox)
dmg_sulf  = compute_damage(sol_sulf)
dmg_dex   = compute_damage(sol_dex)
dmg_combo = compute_damage(sol_combo)
t_dmg     = t_hours[:-1]

fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

# Damage trajectories
ax1.plot(t_dmg, dmg_dox,   linewidth=2, color='steelblue',  linestyle='-',  label='Dox only')
ax1.plot(t_dmg, dmg_sulf,  linewidth=2, color='firebrick',  linestyle='--', label='+ Sulforaphane')
ax1.plot(t_dmg, dmg_dex,   linewidth=2, color='forestgreen',linestyle='-.', label='+ Dexrazoxane')
ax1.plot(t_dmg, dmg_combo, linewidth=2, color='purple',     linestyle=':',  label='+ Both')
ax1.set_title('Cumulative Irreversible Damage\nby Intervention Strategy')
ax1.set_xlabel('Time (hours)')
ax1.set_ylabel('Irreversible Damage (AU)')
ax1.legend(fontsize=9)
ax1.grid(True, alpha=0.3)

# Protection bar chart at t=600 min
final_dmg = [dmg_dox[-1], dmg_sulf[-1], dmg_dex[-1], dmg_combo[-1]]
prot = [(1 - d/dmg_dox[-1])*100 for d in final_dmg]
labels_bar = ['Dox only', 'Sulforaphane\n(Nrf2)', 'Dexrazoxane\n(source)', 'Combination']
colors_bar = ['steelblue', 'firebrick', 'forestgreen', 'purple']

bars = ax2.bar(labels_bar, prot, color=colors_bar, alpha=0.85, edgecolor='white', linewidth=1.2)
ax2.set_title('Protection at t=10h\n(% reduction in irreversible damage vs Dox only)')
ax2.set_ylabel('Damage reduction (%)')
ax2.set_ylim(0, 100)
ax2.grid(True, alpha=0.3, axis='y')
for bar, p in zip(bars, prot):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
             f'{p:.1f}%', ha='center', va='bottom', fontsize=10, fontweight='bold')

plt.suptitle('Dexrazoxane vs Sulforaphane: Complementary Mechanisms', fontsize=13)
plt.tight_layout()
plt.savefig('dexrazoxane_damage.png', dpi=150)
print("Saved: dexrazoxane_damage.png")

print(f"\nProtection at t=10h:")
print(f"  Sulforaphane:  {prot[1]:.1f}%")
print(f"  Dexrazoxane:   {prot[2]:.1f}%")
print(f"  Combination:   {prot[3]:.1f}%")

# ════════════════════════════════════════════════════════════
# PHASE PORTRAIT — H2O2 vs GSH
# ════════════════════════════════════════════════════════════

fig3, axes3 = plt.subplots(1, 2, figsize=(14, 6))

# Panel A: Phase portrait H2O2 vs GSH for all four conditions
ax = axes3[0]
for sol, label, color, ls in conditions:
    ax.plot(sol.y[2], sol.y[5], linewidth=2, color=color,
            linestyle=ls, label=label)
    # Mark start and end
    ax.plot(sol.y[2][0],  sol.y[5][0],  'o', color=color, markersize=8)
    ax.plot(sol.y[2][-1], sol.y[5][-1], 's', color=color, markersize=8)

ax.set_xlabel('H₂O₂ (nM)', fontsize=12)
ax.set_ylabel('GSH (nM)', fontsize=12)
ax.set_title('Phase Portrait: H₂O₂ vs GSH\n(Circle = start, Square = end at t=10h)', fontsize=11)
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
# Annotate steady state region
ax.axhline(GSH_ss, color='gray', linewidth=1, linestyle=':', alpha=0.5)
ax.text(0.1, GSH_ss + 2, 'Baseline GSH', fontsize=8, color='gray')

# Panel B: Phase portrait Keap1 vs Nrf2 — the stress sensor dynamics
ax2_p = axes3[1]
for sol, label, color, ls in conditions:
    ax2_p.plot(sol.y[3], sol.y[4], linewidth=2, color=color,
               linestyle=ls, label=label)
    ax2_p.plot(sol.y[3][0],  sol.y[4][0],  'o', color=color, markersize=8)
    ax2_p.plot(sol.y[3][-1], sol.y[4][-1], 's', color=color, markersize=8)

ax2_p.set_xlabel('Keap1 (nM)', fontsize=12)
ax2_p.set_ylabel('Nrf2 (nM)', fontsize=12)
ax2_p.set_title('Phase Portrait: Keap1 vs Nrf2\n(Stress sensor — antioxidant response axis)', fontsize=11)
ax2_p.legend(fontsize=8)
ax2_p.grid(True, alpha=0.3)
ax2_p.axvline(Keap1_ss, color='gray', linewidth=1, linestyle=':', alpha=0.5)
ax2_p.text(Keap1_ss + 0.2, 0.15, 'Baseline Keap1', fontsize=8, color='gray', rotation=90)

plt.suptitle('System Dynamics: Phase Portraits\nTrajectory of the ROS-Nrf2 System Under Dox Exposure',
             fontsize=13)
plt.tight_layout()
plt.savefig('phase_portraits.png', dpi=150)
print("Saved: phase_portraits.png")
print("\nAll figures complete.")
