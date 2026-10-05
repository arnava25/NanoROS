"""
ros_core_paper.py
=================
Single source of truth for the ROS-Cardio ODE model exactly as used to
produce the PLOS ONE manuscript (PONE-D-26-28903) results, refactored so the
sensitivity / uncertainty / goodness-of-fit analyses all call the same model.

The model is the 8-variable system (6 species + D_acute + D_chronic) from
model/legacy/. Intervention definitions are reproduced verbatim:
  - sulforaphane : k_ox_sustained = 0.002
  - dexrazoxane  : k_ros scaled by 0.4 (60% reduction)
  - combination  : both

NOTE ON k_syn: the headline protection numbers (32.7 / 52.6 / 71.0 %) were
generated with k_syn = 0.02 (Keap1_ss = 20 nM). Fig 1 / Results text use
k_syn = 0.1 (Keap1_ss = 100 nM). This module defaults to 0.02 to match the
results under review; set BASE['k_syn'] = 0.1 to test the alternative.
"""
import numpy as np
from scipy.integrate import solve_ivp

# base parameter vector (matches damage/population/comparator scripts).
# Damage constants and Km_gsh are included here so the sensitivity/uncertainty
# analyses can perturb them; the legacy scripts held them as module constants
# with these same values.
BASE = {
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
    'Km_gsh':           500.0,
    'k_damage_acute':   0.1,
    'k_repair_acute':   0.001,
    'k_damage_chronic': 0.0005,
    'k_chronic_gsh':    0.0002,
}


def steady_state(p):
    Keap1_ss = p['k_syn'] / p['k_deg']
    Nrf2_ss  = p['k_nrf2'] / (p['k_nrf2deg'] * Keap1_ss)
    GSH_ss   = (p['k_gshsyn_basal'] + p['k_gshsyn'] * Nrf2_ss) / p['k_gsh_turn']
    return Keap1_ss, Nrf2_ss, GSH_ss


def _rhs(t, y, p, kin_rate, GSH_ss):
    Dox, O2, H2O2, Keap1, Nrf2, GSH, D_acute, D_chronic = y
    GSH   = max(GSH,   0.0)
    Keap1 = max(Keap1, 0.0)
    Nrf2  = max(Nrf2,  0.0)
    gsh_sat = GSH / (p['Km_gsh'] + GSH)
    k_ox_t  = p['k_ox'] + p['k_ox_sustained']
    GSH_def = max(GSH_ss - GSH, 0.0)

    dDox  = kin_rate - p['k_cl'] * Dox
    dO2   = p['k_ros'] * Dox - p['k_sod'] * O2
    dH2O2 = p['k_sod'] * O2 - p['k_gsh'] * gsh_sat * H2O2
    dK    = (p['k_syn'] - p['k_deg'] * Keap1
             - k_ox_t * H2O2 * Keap1
             - p['k_ox_sustained'] * Keap1)
    dN    = p['k_nrf2'] - p['k_nrf2deg'] * Keap1 * Nrf2
    dG    = (p['k_gshsyn_basal'] + p['k_gshsyn'] * Nrf2
             - p['k_gshdeg'] * gsh_sat * H2O2
             - p['k_gsh_turn'] * GSH)
    dDa   = p['k_damage_acute'] * H2O2 - p['k_repair_acute'] * D_acute
    dDc   = p['k_damage_chronic'] * H2O2 + p['k_chronic_gsh'] * GSH_def
    return [dDox, dO2, dH2O2, dK, dN, dG, dDa, dDc]


def acute_run(p, t_end=600.0, n_eval=600):
    """Continuous-exposure acute run (constant influx). Returns (t, Y[8, n])."""
    Keap1_ss, Nrf2_ss, GSH_ss = steady_state(p)
    y0 = [0.0, 0.0, 0.0, Keap1_ss, Nrf2_ss, GSH_ss, 0.0, 0.0]
    t_eval = np.linspace(0, t_end, n_eval)
    sol = solve_ivp(_rhs, (0, t_end), y0, t_eval=t_eval,
                    args=(p, p['k_in'], GSH_ss), method='LSODA',
                    rtol=1e-6, atol=1e-9)
    return sol.t, sol.y


def peak_h2o2(p):
    _, Y = acute_run(p)
    return float(np.max(Y[2]))


def acute_chronic_damage(p, t_end=600.0):
    """D_chronic at t_end under continuous exposure (used by comparator fig)."""
    t, Y = acute_run(p, t_end=t_end, n_eval=2000)
    Keap1_ss, Nrf2_ss, GSH_ss = steady_state(p)
    H2O2 = Y[2]; GSH = Y[5]
    gsh_def = np.maximum(GSH_ss - GSH, 0.0)
    rate = p['k_damage_chronic'] * H2O2 + p['k_chronic_gsh'] * gsh_def
    dt = np.diff(t)
    return float(np.sum(rate[:-1] * dt))


def simulate_patient_cycle6(p, cycle_days=21, n_cycles=6, dose_min=60):
    """Fast cycle-by-cycle 6-cycle sim. Returns D_chronic at end of cycle 6."""
    Keap1_ss, Nrf2_ss, GSH_ss = steady_state(p)
    state = np.array([0.0, 0.0, 0.0, Keap1_ss, Nrf2_ss, GSH_ss, 0.0, 0.0])
    cycle_dur = cycle_days * 24 * 60
    for _ in range(n_cycles):
        # dose window
        s1 = solve_ivp(_rhs, (0, dose_min), state, method='LSODA',
                       rtol=1e-4, atol=1e-7,
                       args=(p, p['k_in'] * 20, GSH_ss))
        state = s1.y[:, -1]
        # recovery
        s2 = solve_ivp(_rhs, (0, cycle_dur - dose_min), state, method='LSODA',
                       rtol=1e-4, atol=1e-7,
                       args=(p, 0.0, GSH_ss))
        state = s2.y[:, -1]
    return float(state[7])


# intervention constructors -------------------------------------------------
def with_sulforaphane(p):
    q = p.copy(); q['k_ox_sustained'] = 0.002; return q

def with_dexrazoxane(p):
    q = p.copy(); q['k_ros'] = p['k_ros'] * 0.4; return q

def with_combination(p):
    q = p.copy(); q['k_ros'] = p['k_ros'] * 0.4; q['k_ox_sustained'] = 0.002; return q
