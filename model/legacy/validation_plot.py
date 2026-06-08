import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Parameters (same as final model) ────────────────────────
params = {
    'k_in':           0.5,
    'k_cl':           0.02,
    'k_ros':          0.5,
    'k_sod':          5.0,
    'k_gsh':          2.0,
    'k_syn':          0.1,
    'k_deg':          0.001,
    'k_ox':           0.001,
    'k_nrf2':         0.1,
    'k_nrf2deg':      0.01,
    'k_gshsyn':       0.1,
    'k_gshdeg':       0.02,
    'k_gshsyn_basal': 5.0,
    'k_gsh_turn':     0.001,
}

Km_gsh = 500.0

# ── Pre-dox steady state ─────────────────────────────────────
Keap1_ss = params['k_syn'] / params['k_deg']
Nrf2_ss  = params['k_nrf2'] / (params['k_nrf2deg'] * Keap1_ss)
GSH_ss   = (params['k_gshsyn_basal'] + params['k_gshsyn'] * Nrf2_ss) / params['k_gsh_turn']
y0 = [0.0, 0.0, 0.0, Keap1_ss, Nrf2_ss, GSH_ss]

# ── ODE system ───────────────────────────────────────────────
def odes(t, y, p):
    Dox, O2, H2O2, Keap1, Nrf2, GSH = y
    GSH   = max(GSH,   0.0)
    Keap1 = max(Keap1, 0.0)
    Nrf2  = max(Nrf2,  0.0)
    gsh_sat = GSH / (Km_gsh + GSH)

    dDox   = p['k_in'] - p['k_cl'] * Dox
    dO2    = p['k_ros'] * Dox - p['k_sod'] * O2
    dH2O2  = p['k_sod'] * O2 - p['k_gsh'] * gsh_sat * H2O2
    dKeap1 = p['k_syn'] - p['k_deg'] * Keap1 - p['k_ox'] * H2O2 * Keap1
    dNrf2  = p['k_nrf2'] - p['k_nrf2deg'] * Keap1 * Nrf2
    dGSH   = (p['k_gshsyn_basal'] + p['k_gshsyn'] * Nrf2
              - p['k_gshdeg'] * gsh_sat * H2O2
              - p['k_gsh_turn'] * GSH)
    return [dDox, dO2, dH2O2, dKeap1, dNrf2, dGSH]

# ── Solve over 24 hours ──────────────────────────────────────
t_span = (0, 1440)   # 0 to 1440 min = 24h
t_eval = np.linspace(0, 1440, 2000)

params_intervention = params.copy()
params_intervention['k_ox'] = 0.005

sol     = solve_ivp(odes, t_span, y0, t_eval=t_eval, args=(params,),              method='RK45')
sol_int = solve_ivp(odes, t_span, y0, t_eval=t_eval, args=(params_intervention,), method='RK45')

# ── Normalize H2O2 to % of baseline (t=0 value) ─────────────
# At t=0, H2O2=0 so we normalize to first non-zero value ~t=5min
# Instead use first timepoint where dox is just entering — normalize to mean of first 10 min
idx_10min = np.searchsorted(sol.t, 10)
H2O2_base = np.mean(sol.y[2][:idx_10min]) if np.mean(sol.y[2][:idx_10min]) > 0 else sol.y[2][50]

# Use the value at t=60min as baseline (1h, pre-significant-rise) for normalization
# This matches Ludke's "control = 100%" framing
idx_1h  = np.searchsorted(sol.t, 60)
H2O2_at_1h      = sol.y[2][idx_1h]
H2O2_int_at_1h  = sol_int.y[2][idx_1h]

# Normalize both trajectories so that t=1h = ~120% of a notional control
# Ludke's control is untreated cells; we set our t=1h dox value = 120 (their 1h reading)
scale_factor = 120.0 / H2O2_at_1h if H2O2_at_1h > 0 else 1.0

H2O2_norm     = sol.y[2]     * scale_factor
H2O2_int_norm = sol_int.y[2] * scale_factor

# Convert time to hours for x-axis
t_hours = sol.t / 60.0

# ── Ludke 2017 data points (read from Fig 2A) ────────────────
ludke_times     = [1, 3, 6, 12, 24]   # hours
ludke_dox       = [120, 145, 155, 180, 200]   # % control ± ~10
ludke_dox_sem   = [8,   10,  10,  12,  11]
ludke_vitc_dox  = [110, 130, 140, 150, 160]
ludke_vitc_sem  = [7,   9,   9,   10,  10]

# ── Validation plot ──────────────────────────────────────────
fig, ax = plt.subplots(figsize=(9, 6))

# Model trajectories
ax.plot(t_hours, H2O2_norm,     linewidth=2.5, color='steelblue',
        label='Model: Dox only')
ax.plot(t_hours, H2O2_int_norm, linewidth=2.5, color='firebrick',
        linestyle='--', label='Model: Dox + Nrf2 activator')

# Ludke data points with error bars
ax.errorbar(ludke_times, ludke_dox,
            yerr=ludke_dox_sem,
            fmt='o', color='steelblue', markersize=9,
            capsize=4, linewidth=2,
            label='Ludke 2017: Dox (10 µM)')
ax.errorbar(ludke_times, ludke_vitc_dox,
            yerr=ludke_vitc_sem,
            fmt='s', color='firebrick', markersize=9,
            capsize=4, linewidth=2, linestyle='none',
            label='Ludke 2017: Dox + Vit C')

# Reference line at 100% (untreated control)
ax.axhline(100, color='gray', linewidth=1.2, linestyle=':', label='Untreated control (100%)')

ax.set_xlabel('Time (hours)', fontsize=13)
ax.set_ylabel('ROS / H₂O₂ (% of control)', fontsize=13)
ax.set_title('Model Validation: H₂O₂ Trajectory vs Ludke et al. 2017\n'
             'Adult rat cardiomyocytes, Dox exposure', fontsize=12)
ax.set_xlim(0, 25)
ax.set_ylim(80, 230)
ax.legend(fontsize=9, loc='upper left')
ax.grid(True, alpha=0.3)

# Annotate the comparison
ax.annotate('Model trajectories normalized\nto match Ludke t=1h baseline',
            xy=(1, 120), xytext=(8, 195),
            fontsize=8, color='gray',
            arrowprops=dict(arrowstyle='->', color='gray', lw=1))

plt.tight_layout()
plt.savefig('validation_plot.png', dpi=150)
plt.show()
print("Validation plot saved.")
