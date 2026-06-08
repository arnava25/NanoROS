import numpy as np
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ── Base parameters ───────────────────────────────────────────
base_params = {
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

def compute_ss(p):
    Keap1_ss = p['k_syn'] / p['k_deg']
    Nrf2_ss  = p['k_nrf2'] / (p['k_nrf2deg'] * Keap1_ss)
    GSH_ss   = (p['k_gshsyn_basal'] + p['k_gshsyn'] * Nrf2_ss) / p['k_gsh_turn']
    return Keap1_ss, Nrf2_ss, GSH_ss

# ── FAST approach: simulate each cycle separately ─────────────
# Instead of one giant 126-day ODE, solve each 21-day cycle
# independently with just 2 phases: dose window + recovery.
# This is ~50x faster because we avoid integrating near-zero dynamics.

cycle_days        = 21
n_cycles          = 6
dose_min          = 60   # infusion duration

def simulate_patient(p, intervention=False):
    """
    Simulate 6 cycles. Returns D_chronic at end of each cycle and GSH nadir.
    """
    p = p.copy()
    if intervention:
        p['k_ox_sustained'] = 0.002

    Keap1_ss, Nrf2_ss, GSH_ss = compute_ss(p)
    # Start at steady state
    state = np.array([0.0, 0.0, 0.0, Keap1_ss, Nrf2_ss, GSH_ss, 0.0, 0.0])

    D_chronic_per_cycle = []
    gsh_min = GSH_ss

    def odes(t, y):
        Dox, O2, H2O2, Keap1, Nrf2, GSH, D_acute, D_chronic = y
        GSH   = max(GSH,   0.0)
        Keap1 = max(Keap1, 0.0)
        Nrf2  = max(Nrf2,  0.0)
        gsh_sat = GSH / (Km_gsh + GSH)
        k_ox_t  = p['k_ox'] + p['k_ox_sustained']
        GSH_def = max(GSH_ss - GSH, 0.0)

        dDox  = kin_rate - p['k_cl'] * Dox
        dO2   = p['k_ros'] * Dox - p['k_sod'] * O2
        dH2O2 = p['k_sod'] * O2 - p['k_gsh'] * gsh_sat * H2O2
        dK    = p['k_syn'] - p['k_deg'] * Keap1 \
                - k_ox_t * H2O2 * Keap1 \
                - p['k_ox_sustained'] * Keap1
        dN    = p['k_nrf2'] - p['k_nrf2deg'] * Keap1 * Nrf2
        dG    = (p['k_gshsyn_basal'] + p['k_gshsyn'] * Nrf2
                 - p['k_gshdeg'] * gsh_sat * H2O2
                 - p['k_gsh_turn'] * GSH)
        dDa   = k_damage_acute * H2O2 - k_repair_acute * D_acute
        dDc   = k_damage_chronic * H2O2 + k_chronic_gsh * GSH_def
        return [dDox, dO2, dH2O2, dK, dN, dG, dDa, dDc]

    for cycle in range(n_cycles):
        cycle_dur = cycle_days * 24 * 60

        # Phase 1: dose window (60 min, coarse)
        kin_rate = p['k_in'] * 20
        s1 = solve_ivp(odes, (0, dose_min), state,
                       method='RK45', max_step=10, rtol=1e-3, atol=1e-6,
                       dense_output=False)
        state = s1.y[:, -1]
        gsh_min = min(gsh_min, np.min(s1.y[5]))

        # Phase 2: recovery (rest of cycle, very coarse)
        kin_rate = 0.0
        recovery_dur = cycle_dur - dose_min
        s2 = solve_ivp(odes, (0, recovery_dur), state,
                       method='RK45', max_step=720, rtol=1e-3, atol=1e-6,
                       dense_output=False)
        state = s2.y[:, -1]
        gsh_min = min(gsh_min, np.min(s2.y[5]))
        D_chronic_per_cycle.append(state[7])

    return np.array(D_chronic_per_cycle), gsh_min

# ── Patient population ────────────────────────────────────────
np.random.seed(42)
n_patients = 200

k_ros_values       = np.random.lognormal(mean=np.log(0.5),  sigma=0.4, size=n_patients)
k_ros_values       = np.clip(k_ros_values, 0.1, 2.0)
k_gsh_basal_values = np.random.lognormal(mean=np.log(5.0),  sigma=0.3, size=n_patients)
k_gsh_basal_values = np.clip(k_gsh_basal_values, 1.0, 15.0)

print(f"Simulating {n_patients} patients (fast cycle-by-cycle method)...")
print(f"k_ros range:      {k_ros_values.min():.3f} — {k_ros_values.max():.3f}")
print(f"k_gsh_basal range:{k_gsh_basal_values.min():.2f} — {k_gsh_basal_values.max():.2f}")

t_start = time.time()
final_damage_dox = np.zeros(n_patients)
final_damage_int = np.zeros(n_patients)
gsh_nadir_dox    = np.zeros(n_patients)
gsh_nadir_int    = np.zeros(n_patients)
all_cycles_dox   = np.zeros((n_patients, n_cycles))
all_cycles_int   = np.zeros((n_patients, n_cycles))

for i in range(n_patients):
    p = base_params.copy()
    p['k_ros']          = k_ros_values[i]
    p['k_gshsyn_basal'] = k_gsh_basal_values[i]

    cyc_dox, gsh_d = simulate_patient(p, intervention=False)
    cyc_int, gsh_i = simulate_patient(p, intervention=True)

    all_cycles_dox[i]   = cyc_dox
    all_cycles_int[i]   = cyc_int
    final_damage_dox[i] = cyc_dox[-1]
    final_damage_int[i] = cyc_int[-1]
    gsh_nadir_dox[i]    = gsh_d
    gsh_nadir_int[i]    = gsh_i

    if (i + 1) % 10 == 0:
        elapsed   = time.time() - t_start
        rate      = (i + 1) / elapsed
        remaining = (n_patients - i - 1) / rate
        pct       = (i + 1) / n_patients * 100
        print(f"  [{pct:5.1f}%] {i+1}/{n_patients} | "
              f"{elapsed:.0f}s elapsed | ~{remaining:.0f}s remaining")

elapsed_total = time.time() - t_start
print(f"\nDone. Total time: {elapsed_total:.0f}s")

prot_pct = (1 - final_damage_int / np.maximum(final_damage_dox, 1e-9)) * 100

# ── Plot 1: Population summary ────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# A: Damage distribution
axes[0,0].hist(final_damage_dox, bins=30, color='steelblue', alpha=0.7,
               label=f'Dox only (median={np.median(final_damage_dox):.1f})')
axes[0,0].hist(final_damage_int, bins=30, color='firebrick', alpha=0.7,
               label=f'+ Sulforaphane (median={np.median(final_damage_int):.1f})')
axes[0,0].axvline(np.median(final_damage_dox), color='steelblue', linewidth=2, linestyle='--')
axes[0,0].axvline(np.median(final_damage_int), color='firebrick', linewidth=2, linestyle='--')
axes[0,0].set_title('Irreversible Damage Distribution at Cycle 6')
axes[0,0].set_xlabel('Irreversible Damage (AU)')
axes[0,0].set_ylabel('Number of patients')
axes[0,0].legend(fontsize=9)
axes[0,0].grid(True, alpha=0.3)

# B: Risk stratification
sc = axes[0,1].scatter(k_ros_values, final_damage_dox,
                        c=k_gsh_basal_values, cmap='RdYlGn', s=35, alpha=0.8)
axes[0,1].scatter(k_ros_values, final_damage_int,
                   c='firebrick', s=15, alpha=0.4, marker='^', label='+ Sulforaphane')
plt.colorbar(sc, ax=axes[0,1], label='Baseline GSH synthesis rate')
axes[0,1].set_title('Risk Stratification by Mitochondrial ROS Rate\n(Color = antioxidant capacity)')
axes[0,1].set_xlabel('k_ros (mitochondrial ROS generation)')
axes[0,1].set_ylabel('Irreversible Damage at Cycle 6 (AU)')
axes[0,1].legend(fontsize=9)
axes[0,1].grid(True, alpha=0.3)

# C: Protection distribution
axes[1,0].hist(prot_pct, bins=30, color='purple', alpha=0.7)
axes[1,0].axvline(np.median(prot_pct), color='black', linewidth=2,
                   linestyle='--', label=f'Median = {np.median(prot_pct):.1f}%')
axes[1,0].axvline(np.percentile(prot_pct, 25), color='gray', linewidth=1.5, linestyle=':')
axes[1,0].axvline(np.percentile(prot_pct, 75), color='gray', linewidth=1.5,
                   linestyle=':', label=f'IQR: {np.percentile(prot_pct,25):.1f}–{np.percentile(prot_pct,75):.1f}%')
axes[1,0].set_title('Distribution of Sulforaphane Protection\nacross Patient Population')
axes[1,0].set_xlabel('Damage reduction (%)')
axes[1,0].set_ylabel('Number of patients')
axes[1,0].legend(fontsize=9)
axes[1,0].grid(True, alpha=0.3)

# D: Who benefits most
axes[1,1].scatter(k_ros_values, prot_pct, c=k_gsh_basal_values,
                   cmap='RdYlGn', s=35, alpha=0.8)
axes[1,1].axhline(np.median(prot_pct), color='black', linewidth=1.5,
                   linestyle='--', label=f'Median: {np.median(prot_pct):.1f}%')
axes[1,1].set_title('Who Benefits Most?\n(Color = antioxidant capacity)')
axes[1,1].set_xlabel('k_ros (mitochondrial ROS generation)')
axes[1,1].set_ylabel('Protection (%)')
axes[1,1].legend(fontsize=9)
axes[1,1].grid(True, alpha=0.3)

plt.suptitle(f'Patient Variability Simulation (n={n_patients})', fontsize=13)
plt.tight_layout()
plt.savefig('patient_variability.png', dpi=150)
print("Saved: patient_variability.png")

# ── Plot 2: Damage accumulation per cycle across population ───
fig2, axes2 = plt.subplots(1, 2, figsize=(14, 6))
cycle_nums = np.arange(1, n_cycles + 1)

# Percentile bands
for data, ax, title in [
    (all_cycles_dox, axes2[0], 'Dox Only'),
    (all_cycles_int, axes2[1], 'Dox + Sulforaphane'),
]:
    p10 = np.percentile(data, 10, axis=0)
    p25 = np.percentile(data, 25, axis=0)
    p50 = np.percentile(data, 50, axis=0)
    p75 = np.percentile(data, 75, axis=0)
    p90 = np.percentile(data, 90, axis=0)

    ax.fill_between(cycle_nums, p10, p90, alpha=0.15, color='steelblue' if 'Only' in title else 'firebrick', label='10–90th %ile')
    ax.fill_between(cycle_nums, p25, p75, alpha=0.3,  color='steelblue' if 'Only' in title else 'firebrick', label='25–75th %ile')
    ax.plot(cycle_nums, p50, 'o-', linewidth=2.5,
            color='steelblue' if 'Only' in title else 'firebrick', label='Median')
    ax.set_title(f'{title}: Damage by Cycle\n(Population percentile bands)')
    ax.set_xlabel('Cycle number')
    ax.set_ylabel('Irreversible Damage (AU)')
    ax.set_xticks(cycle_nums)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

plt.suptitle(f'Population-Level Damage Accumulation (n={n_patients})', fontsize=13)
plt.tight_layout()
plt.savefig('patient_variability_trajectories.png', dpi=150)
print("Saved: patient_variability_trajectories.png")

# ── Summary stats ─────────────────────────────────────────────
print(f"\n── Population Summary ──────────────────────────────")
print(f"Median damage (Dox only):       {np.median(final_damage_dox):.1f} AU")
print(f"Median damage (+Sulforaphane):  {np.median(final_damage_int):.1f} AU")
print(f"Median protection:              {np.median(prot_pct):.1f}%")
print(f"IQR of protection:              {np.percentile(prot_pct,25):.1f}% — {np.percentile(prot_pct,75):.1f}%")
high = k_ros_values > np.median(k_ros_values)
print(f"High-risk (k_ros > median) n:   {np.sum(high)}")
print(f"Protection in high-risk:        {np.median(prot_pct[high]):.1f}%")
print(f"Protection in low-risk:         {np.median(prot_pct[~high]):.1f}%")