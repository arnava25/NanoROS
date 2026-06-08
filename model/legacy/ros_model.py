import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Parameters ───────────────────────────────────────────────
params = {
    'k_in':           0.5,    # nM/min  — dox influx
    'k_cl':           0.02,   # min⁻¹  — dox clearance
    'k_ros':          0.5,    # min⁻¹  — superoxide generation per dox
    'k_sod':          5.0,    # min⁻¹  — SOD conversion rate
    'k_gsh':          2.0,    # min⁻¹  — GSH-mediated H2O2 clearance (MM)
    'k_syn':          0.1,    # nM/min — Keap1 synthesis
    'k_deg':          0.001,  # min⁻¹  — Keap1 basal degradation
    'k_ox':           0.001,  # nM⁻¹min⁻¹ — Keap1 oxidation by H2O2
    'k_nrf2':         0.1,    # nM/min — Nrf2 synthesis
    'k_nrf2deg':      0.01,   # nM⁻¹min⁻¹ — Nrf2 degradation via Keap1
    'k_gshsyn':       0.1,    # nM/min per nM Nrf2
    'k_gshdeg':       0.02,   # min⁻¹  — GSH consumption by H2O2 (MM)
    'k_gshsyn_basal': 5.0,    # nM/min — constitutive GSH synthesis
    'k_gsh_turn':     0.001,  # min⁻¹  — basal GSH turnover (GSSG export/recycling)
}

Km_gsh = 500.0  # nM — MM half-saturation for GSH

# ── Pre-dox steady state (analytical) ───────────────────────
# With Dox=0: H2O2=0, O2=0
# Keap1_ss = k_syn / k_deg
# Nrf2_ss  = k_nrf2 / (k_nrf2deg * Keap1_ss)
# GSH_ss   = k_gshsyn_basal + k_gshsyn * Nrf2_ss) / k_gsh_turn
p = params
Keap1_ss = p['k_syn'] / p['k_deg']
Nrf2_ss  = p['k_nrf2'] / (p['k_nrf2deg'] * Keap1_ss)
GSH_ss   = (p['k_gshsyn_basal'] + p['k_gshsyn'] * Nrf2_ss) / p['k_gsh_turn']

y0 = [0.0, 0.0, 0.0, Keap1_ss, Nrf2_ss, GSH_ss]
print(f"Pre-dox steady state:")
print(f"  Keap1 = {Keap1_ss:.1f} nM")
print(f"  Nrf2  = {Nrf2_ss:.4f} nM")
print(f"  GSH   = {GSH_ss:.0f} nM")

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

# ── ODE with damage variable ─────────────────────────────────
k_damage = 0.1
k_repair = 0.001

def odes_with_damage(t, y, p):
    Dox, O2, H2O2, Keap1, Nrf2, GSH, D = y
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
    dD     = k_damage * H2O2 - k_repair * D

    return [dDox, dO2, dH2O2, dKeap1, dNrf2, dGSH, dD]

# ── Time setup ───────────────────────────────────────────────
t_span = (0, 600)
t_eval = np.linspace(0, 600, 1000)
labels = ['Dox', 'O2·⁻', 'H2O2', 'Keap1', 'Nrf2', 'GSH']

# ── Plot 1: Baseline ─────────────────────────────────────────
sol = solve_ivp(odes, t_span, y0, t_eval=t_eval, args=(params,), method='RK45')

fig, axes = plt.subplots(2, 3, figsize=(14, 8))
axes = axes.flatten()
for i, ax in enumerate(axes):
    ax.plot(sol.t, sol.y[i], linewidth=2)
    ax.set_title(labels[i])
    ax.set_xlabel('Time (min)')
    ax.set_ylabel('Concentration (nM)')
    ax.grid(True, alpha=0.3)
plt.suptitle('ROS-Cardio ODE — Baseline Dox Exposure', fontsize=14)
plt.tight_layout()
plt.savefig('ros_model_output.png', dpi=150)
plt.show()
print("1/5 done")

# ── Plot 2: Intervention ─────────────────────────────────────
params_intervention = params.copy()
params_intervention['k_ox'] = 0.005  # 5x baseline — sulforaphane effect

sol_int = solve_ivp(odes, t_span, y0, t_eval=t_eval, args=(params_intervention,), method='RK45')

fig2, axes2 = plt.subplots(2, 3, figsize=(14, 8))
axes2 = axes2.flatten()
for i, ax in enumerate(axes2):
    ax.plot(sol.t,     sol.y[i],     linewidth=2, label='Dox only',         color='steelblue')
    ax.plot(sol_int.t, sol_int.y[i], linewidth=2, label='Dox + Sulforaphane',
            color='firebrick', linestyle='--')
    ax.set_title(labels[i])
    ax.set_xlabel('Time (min)')
    ax.set_ylabel('Concentration (nM)')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
plt.suptitle('ROS-Cardio ODE — Dox vs Dox + Sulforaphane', fontsize=13)
plt.tight_layout()
plt.savefig('ros_model_intervention.png', dpi=150)
plt.show()
print("2/5 done")

# ── Plot 3: Dose-Response ─────────────────────────────────────
kin_values = [0.1, 0.3, 0.5, 1.0, 2.0, 4.0]
colors = plt.cm.RdYlBu_r(np.linspace(0.1, 0.9, len(kin_values)))
fig3, axes3 = plt.subplots(1, 3, figsize=(15, 5))
gsh_nadirs, h2o2_peaks = [], []

for j, kin in enumerate(kin_values):
    p = params.copy()
    p['k_in'] = kin
    sol_dr = solve_ivp(odes, t_span, y0, t_eval=t_eval, args=(p,), method='RK45')
    gsh_nadirs.append(np.min(sol_dr.y[5]))
    h2o2_peaks.append(np.max(sol_dr.y[2]))
    label = f'k_in={kin}'
    axes3[0].plot(sol_dr.t, sol_dr.y[2], color=colors[j], label=label, linewidth=2)
    axes3[1].plot(sol_dr.t, sol_dr.y[5], color=colors[j], label=label, linewidth=2)
    axes3[2].plot(sol_dr.t, sol_dr.y[4], color=colors[j], label=label, linewidth=2)

for ax, title in zip(axes3, ['H2O2 (nM)', 'GSH (nM)', 'Nrf2 (nM)']):
    ax.set_title(title)
    ax.set_xlabel('Time (min)')
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)
plt.suptitle('Dose-Response: Varying Dox Influx Rate', fontsize=13)
plt.tight_layout()
plt.savefig('dose_response.png', dpi=150)
plt.show()
print("3/5 done")

# ── Plot 4: Sensitivity Analysis ──────────────────────────────
param_names = list(params.keys())
gsh_base  = sol.y[5][-1]          # use final GSH value now that it reaches SS
h2o2_base = max(np.max(sol.y[2]), 1e-9)
sa_gsh, sa_h2o2 = [], []

for pname in param_names:
    effects_gsh, effects_h2o2 = [], []
    for fold in [0.5, 2.0]:
        p = params.copy()
        p[pname] *= fold
        sol_sa = solve_ivp(odes, t_span, y0, t_eval=t_eval, args=(p,), method='RK45')
        effects_gsh.append(sol_sa.y[5][-1])
        effects_h2o2.append(np.max(sol_sa.y[2]))
    sa_gsh.append((effects_gsh[1] - effects_gsh[0]) / max(abs(gsh_base), 1.0))
    sa_h2o2.append((effects_h2o2[1] - effects_h2o2[0]) / h2o2_base)

fig4, (ax4a, ax4b) = plt.subplots(1, 2, figsize=(14, 6))
x = np.arange(len(param_names))

ax4a.barh(x, sa_gsh, color='steelblue')
ax4a.set_yticks(x); ax4a.set_yticklabels(param_names)
ax4a.set_title('Sensitivity: GSH at t=600 min')
ax4a.set_xlabel('Normalized sensitivity')
ax4a.axvline(0, color='black', linewidth=0.8)
ax4a.grid(True, alpha=0.3)

ax4b.barh(x, sa_h2o2, color='firebrick')
ax4b.set_yticks(x); ax4b.set_yticklabels(param_names)
ax4b.set_title('Sensitivity: Peak H2O2')
ax4b.set_xlabel('Normalized sensitivity')
ax4b.set_xlim(-5, 5)
ax4b.axvline(0, color='black', linewidth=0.8)
ax4b.grid(True, alpha=0.3)

plt.suptitle('Sensitivity Analysis — ±2x Parameter Perturbation', fontsize=13)
plt.tight_layout()
plt.savefig('sensitivity_analysis.png', dpi=150)
plt.show()
print("4/5 done")

# ── Plot 5: Damage vs Dose ────────────────────────────────────
y0_d = y0 + [0.0]
damage_dox, damage_int = [], []

for kin in kin_values:
    p_base = params.copy();              p_base['k_in'] = kin
    p_int  = params_intervention.copy(); p_int['k_in']  = kin
    s1 = solve_ivp(odes_with_damage, t_span, y0_d, t_eval=t_eval, args=(p_base,), method='RK45')
    s2 = solve_ivp(odes_with_damage, t_span, y0_d, t_eval=t_eval, args=(p_int,),  method='RK45')
    damage_dox.append(s1.y[6][-1])
    damage_int.append(s2.y[6][-1])

fig5, ax5 = plt.subplots(figsize=(8, 5))
ax5.plot(kin_values, damage_dox, 'o-',  linewidth=2, color='steelblue', label='Dox only')
ax5.plot(kin_values, damage_int, 's--', linewidth=2, color='firebrick', label='Dox + Sulforaphane')
ax5.set_title('Cumulative Oxidative Damage at t=600 min vs Dox Dose', fontsize=12)
ax5.set_xlabel('Dox influx rate k_in (nM/min)')
ax5.set_ylabel('Damage at t=600 min (AU)')
ax5.legend()
ax5.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('damage_variable.png', dpi=150)
plt.show()
print("5/5 done. All figures saved.")