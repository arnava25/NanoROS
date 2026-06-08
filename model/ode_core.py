"""
NanoROS — ODE Core
===================
Thin wrapper around ros_model_extended that provides the run() interface
used by analysis and figure scripts.

All parameters come from parameters.py. All ODE math is in ros_model_extended.py.
"""

import numpy as np
from scipy.integrate import solve_ivp
from parameters import PARAMS


def steady_states(p=None):
    """Compute healthy pre-dox steady state values."""
    p = p or PARAMS
    Keap1_ss = p['k_syn'] / p['k_deg']
    Nrf2_ss  = p['k_nrf2'] / (p['k_nrf2deg'] * Keap1_ss)
    GSH_ss   = (p['k_gshsyn_b'] + p['k_gshsyn'] * Nrf2_ss) / p['k_gsh_turn']
    return Keap1_ss, Nrf2_ss, GSH_ss


def odes_base(t, y, p, k_sc=0.0, t_release=99999.0, payload='O2'):
    """
    7-state ODE RHS.
    payload='O2'   → MnTE-2-PyP intercepts O2.- upstream (adds to k_sod)
    payload='H2O2' → catalase-type scavenges H2O2 directly
    """
    Dox, O2, H2O2, Keap1, Nrf2, GSH, D = y
    GSH   = max(GSH,   0.0)
    Keap1 = max(Keap1, 0.0)
    Nrf2  = max(Nrf2,  0.0)
    H2O2  = max(H2O2,  0.0)
    O2    = max(O2,    0.0)

    gsh_sat = GSH / (p['Km_gsh'] + GSH)
    _, _, GSH_ss = steady_states(p)
    GSH_def = max(GSH_ss - GSH, 0.0)
    active  = (t >= t_release)

    k_sod_eff   = p['k_sod'] + (k_sc if active and payload == 'O2'   else 0.0)
    np_scav_h2  = k_sc * H2O2 if active and payload == 'H2O2' else 0.0

    dDox   = p['k_in']  - p['k_cl'] * Dox
    dO2    = p['k_ros'] * Dox - k_sod_eff * O2
    dH2O2  = k_sod_eff * O2 - p['k_gsh'] * gsh_sat * H2O2 - np_scav_h2
    dKeap1 = p['k_syn'] - p['k_deg'] * Keap1 - p['k_ox'] * H2O2 * Keap1
    dNrf2  = p['k_nrf2'] - p['k_nrf2deg'] * Keap1 * Nrf2
    dGSH   = (p['k_gshsyn_b'] + p['k_gshsyn'] * Nrf2
              - p['k_gshdeg'] * gsh_sat * H2O2
              - p['k_gsh_turn'] * GSH)
    dD     = p['k_damage'] * H2O2 + p['k_dam_gsh'] * GSH_def

    return [dDox, dO2, dH2O2, dKeap1, dNrf2, dGSH, dD]


def run(t_end=1440, n_pts=5761, k_sc=0.0, t_release=99999.0,
        payload='O2', params=None):
    """
    Run ODE simulation.

    Returns
    -------
    t : time array (min)
    y : solution array [7 x n_pts]
        y[0]=Dox, y[1]=O2-, y[2]=H2O2, y[3]=Keap1,
        y[4]=Nrf2, y[5]=GSH, y[6]=D
    """
    p = params or PARAMS
    Keap1_ss, Nrf2_ss, GSH_ss = steady_states(p)
    y0     = [0.0, 0.0, 0.0, Keap1_ss, Nrf2_ss, GSH_ss, 0.0]
    t_eval = np.linspace(0, t_end, n_pts)

    sol = solve_ivp(
        fun=lambda t, y: odes_base(t, y, p, k_sc, t_release, payload),
        t_span=(0, t_end),
        y0=y0,
        t_eval=t_eval,
        method='LSODA',
        rtol=1e-8, atol=1e-10,
        max_step=1.0
    )
    return sol.t, sol.y


def key_outputs(params=None):
    """Print key design constraint outputs."""
    p = params or PARAMS
    t, y = run(params=p)
    _, Nrf2_ss, GSH_ss = steady_states(p)

    print("=== NANOROS ODE — KEY OUTPUTS ===")
    print(f"[Dox] SS:              {p['k_in']/p['k_cl']:.1f} nM           (target: 25 nM)")
    print(f"Peak [H2O2]:           {y[2].max():.2f} nM at t={t[y[2].argmax()]:.0f} min  (target: 6.89 nM)")
    print(f"Peak [O2.-]:           {y[1].max():.2f} nM at t={t[y[1].argmax()]:.0f} min  (target: 2.51 nM)")

    idx_D = np.where(y[6] >= p['D_crit'])[0]
    if len(idx_D):
        print(f"Damage deadline D=50:  t={t[idx_D[0]]:.0f} min          (target: 123 min)")
    else:
        print("Damage deadline D=50:  not reached")

    nrf2_norm = y[4] / Nrf2_ss
    idx_nrf2 = np.where((nrf2_norm >= 2.0) & (t > 10))[0]
    if len(idx_nrf2):
        print(f"Nrf2 2x activation:    t={t[idx_nrf2[0]]:.0f} min         (target: 160 min)")

    print(f"GSH min depletion:     {(1 - y[5].min()/GSH_ss)*100:.1f}%              (target: 0.8%)")


if __name__ == '__main__':
    key_outputs()
