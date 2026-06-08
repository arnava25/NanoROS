"""
NanoROS — AND Gate Linker Kinetics
====================================
Arm A (dioxaborolane) and Arm B (O-alkylhydroxylamine) firing times.
References:
  Arm A: Bhattacharya & Chang, JACS 2014
  Arm B: Dirksen & Dawson, Bioconjug Chem 2008
"""

import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'models'))
from ode_core import run
from parameters import AND_GATE


def boronate_fire_time(k2=600.0, f_crit=0.001, params=None):
    t, y = run(params=params)
    H2O2_M = y[2] * 1e-9
    dt_sec = np.diff(t) * 60
    cum    = np.cumsum(H2O2_M[:-1] * dt_sec)
    B_ox   = 1 - np.exp(-k2 * cum)
    idx    = np.where(B_ox >= f_crit)[0]
    return float(t[idx[0]]) if len(idx) > 0 else None


def oxime_fire_time(HNE_uM, k2=1000.0, f_crit=0.15):
    if HNE_uM <= 0:
        return None
    return float(-np.log(1 - f_crit) / (k2 * HNE_uM * 1e-6) / 60.0)


def specificity_table():
    scenarios = [
        ("Dox cardiotox cycle 6",          7.5,  0.090, 600,  "FIRES — target"),
        ("Dox cardiotox cycle 3",          5.0,  0.050, 600,  "FIRES — target"),
        ("Post-MI reperfusion (sustained)",10.0, 0.060, 240,  "FIRES — also therapeutic"),
        ("Sepsis cardiomyopathy",           5.0, 0.030, 600,  "borderline"),
        ("HFpEF baseline",                  4.5, 0.020, 600,  "safe"),
        ("Trastuzumab cardiotox",           3.0, 0.008, 600,  "safe"),
        ("Prolonged exercise (2 hr)",       5.0, 0.013, 120,  "safe — 15% threshold"),
        ("Acute exercise (30 min)",         7.0, 0.005,  30,  "safe"),
        ("Neutrophil oxidative burst",    100.0, 0.002,  10,  "safe — no CL peroxidation"),
        ("Ischaemia-reperfusion",          15.0, 0.053, 240,  "safe (borderline)"),
        ("Healthy cardiomyocyte",           0.5, 0.001, 600,  "safe"),
    ]

    k2_B     = AND_GATE['k2_benzyloxyamine']
    k2_A     = AND_GATE['k2_boronate_pH80']
    f_bor    = AND_GATE['f_crit_boronate']
    f_oxime  = AND_GATE['f_crit_oxime']

    print(f"\n{'Scenario':<35} {'[4-HNE]µM':<11} {'t_B min':<10} {'AND gate':<10} {'Expected'}")
    print("-" * 88)

    for name, H2O2, HNE, dur, expected in scenarios:
        t_B     = oxime_fire_time(HNE, k2=k2_B, f_crit=f_oxime)
        fires_B = (t_B is not None and t_B <= dur)
        cum     = H2O2 * 1e-9 * dur * 60
        fires_A = (1 - np.exp(-k2_A * cum)) >= f_bor
        fires   = fires_A and fires_B
        t_B_str = f"{t_B:.0f}" if t_B else "∞"
        print(f"{name:<35} {HNE:<11.3f} {t_B_str:<10} {'FIRES' if fires else 'safe':<10} {expected}")


if __name__ == '__main__':
    print("=== ARM A (DIOXABOROLANE) ===\n")
    for f_pct in [0.1, 1.0, 5.0]:
        t = boronate_fire_time(k2=600.0, f_crit=f_pct/100)
        print(f"  f_crit={f_pct:.1f}%: t={t:.0f} min" if t else f"  f_crit={f_pct:.1f}%: not reached")

    print("\n=== ARM B (BENZYLOXYAMINE) ===\n")
    deadline = AND_GATE['damage_deadline']
    for HNE in [0.025, 0.034, 0.075, 0.10, 0.20, 0.46]:
        t = oxime_fire_time(HNE, k2=1000.0, f_crit=0.15)
        ok = "✓" if (t and t <= deadline) else "✗"
        print(f"  [4-HNE]={HNE:.3f} µM: t={t:.0f} min {ok}" if t else f"  {HNE}: N/A")

    print("\n=== AND GATE SPECIFICITY TABLE ===")
    specificity_table()
