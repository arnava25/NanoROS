"""
NanoROS — Figure 3: Comparator ODE Figure
==========================================
NanoROS vs dexrazoxane vs MitoQ vs no treatment on same damage axis.

Mechanisms modelled:
  No treatment:   baseline ODE
  Dexrazoxane:    40% reduction in k_ros (iron chelation reduces Complex I/III ROS)
                  Clinical: ~40% cardiotoxicity reduction (Swain 1997)
  MitoQ:          Net ~0 scavenging at 40% residual ΔΨm — pro-oxidant ceiling
                  at intramitochondrial concentrations from therapeutic plasma doses
                  (Doughan & Dikalov, Antioxid Redox Signal 2007)
  NanoROS:        MnTE-2-PyP payload (k_sc_O2 = 23,515 min-1)
                  Releases at t=22 min (boronate Arm A, f_crit=0.1%, pH8.0)
"""

import numpy as np
import copy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'models'))
from ode_core import run, steady_states
from parameters import PARAMS, AND_GATE, PAYLOAD


def make_figure(output_path='fig3_comparator.png'):
    # Baseline — no treatment
    t_base, y_base = run()
    D_no_NP = y_base[6, -1]

    # Dexrazoxane — 40% k_ros reduction
    p_dex = copy.deepcopy(PARAMS)
    p_dex['k_ros'] = PARAMS['k_ros'] * 0.60
    t_dex, y_dex = run(params=p_dex)

    # MitoQ — net ~0.5 min-1 scavenging (pro-oxidant offsets antioxidant)
    # At 40% residual ΔΨm: 9-fold accumulation × 100nM plasma = 0.9 µM
    # → approaching 1 µM pro-oxidant ceiling → net benefit minimal
    t_mq, y_mq = run(k_sc=0.5, t_release=0, payload='H2O2')

    # NanoROS — MnTE-2-PyP, boronate fires at t=22 min
    t_fire = AND_GATE['arm_A_fire_t']
    # Combined payload (MnTE-2-PyP + Mn cyclen): model H2O2 scavenging arm
    # SOD arm (MnTE-2-PyP) reduces O2- but does not change H2O2_ss (cancels algebraically)
    # Catalase arm (Mn cyclen) directly reduces H2O2 -- this drives the protection
    t_nr, y_nr = run(k_sc=PAYLOAD['k_sc_O2_at_EE'], t_release=t_fire, payload='H2O2')

    # ── Compute summary stats ──────────────────────────────────────────────────
    datasets = [
        ('No Treatment',  t_base, y_base, '#C00000', '--'),
        ('Dexrazoxane',   t_dex,  y_dex,  '#ED7D31', '-.'),
        ('MitoQ\n(partial ΔΨm)', t_mq, y_mq, '#7030A0', ':'),
        ('NanoROS',       t_nr,   y_nr,   '#2E75B6', '-'),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))
    fig.patch.set_facecolor('#F8F9FA')

    # ── Panel A: H2O2 trajectory ──
    ax = axes[0]
    for label, t, y, c, ls in datasets:
        ax.plot(t, y[2], color=c, ls=ls, lw=2.2, label=label)
    ax.axvline(t_fire, color='#2E75B6', ls=':', lw=1.5, alpha=0.6,
               label=f'NP release t={t_fire:.0f} min')
    ax.axvline(123, color='gray', ls='--', lw=1.2, alpha=0.5,
               label='Damage deadline (123 min)')
    ax.set_xlabel('Time (min)', fontsize=10)
    ax.set_ylabel('[H₂O₂] (nM)', fontsize=10)
    ax.set_title('A  H₂O₂ Trajectory', fontsize=11, fontweight='bold')
    ax.legend(fontsize=8, loc='upper right')
    ax.set_facecolor('white')
    ax.grid(True, alpha=0.3)

    # ── Panel B: Cumulative damage ──
    ax2 = axes[1]
    for label, t, y, c, ls in datasets:
        ax2.plot(t, y[6], color=c, ls=ls, lw=2.2, label=label)
    ax2.axhline(PARAMS['D_crit'], color='black', ls='--', lw=1.8,
                label=f"Clinical threshold (D={PARAMS['D_crit']} AU)")
    ax2.set_xlabel('Time (min)', fontsize=10)
    ax2.set_ylabel('Cumulative Damage D (AU)', fontsize=10)
    ax2.set_title('B  Cumulative Cardiac Damage', fontsize=11, fontweight='bold')
    ax2.legend(fontsize=8)
    ax2.set_facecolor('white')
    ax2.grid(True, alpha=0.3)

    # ── Panel C: Protection bar chart ──
    ax3 = axes[2]
    bar_labels = ['Dexrazoxane', 'MitoQ\n(partial ΔΨm)', 'NanoROS']
    d_finals = [y_dex[6,-1], y_mq[6,-1], y_nr[6,-1]]
    prots = [(1 - d/D_no_NP)*100 for d in d_finals]
    bar_colors = ['#ED7D31', '#7030A0', '#2E75B6']

    bars = ax3.bar(bar_labels, prots, color=bar_colors,
                   edgecolor='white', linewidth=2, width=0.5)

    # Clinical benchmark line for dexrazoxane
    ax3.axhline(40, color='#ED7D31', ls='--', lw=1.8, alpha=0.7,
                label='Dex clinical benchmark (~40%)')

    for bar, prot in zip(bars, prots):
        ax3.text(bar.get_x() + bar.get_width()/2,
                 bar.get_height() + 1.0,
                 f'{prot:.0f}%', ha='center', va='bottom',
                 fontsize=13, fontweight='bold')

    ax3.set_ylim(0, 100)
    ax3.set_ylabel('Per-Cycle Damage Reduction (%)', fontsize=10)
    ax3.set_title('C  Protection Comparison\n(vs no-treatment control)',
                  fontsize=11, fontweight='bold')
    ax3.legend(fontsize=9)
    ax3.set_facecolor('white')
    ax3.grid(True, alpha=0.3, axis='y')

    fig.suptitle(
        'NanoROS vs Existing Interventions — ODE Comparator\n'
        'Single doxorubicin cycle, 25 nM steady state | '
        'Validated ROS-Cardio model (Ludke et al. 2017)',
        fontsize=12, fontweight='bold', y=1.02
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='#F8F9FA')
    plt.close()
    print(f"Saved: {output_path}")

    print(f"\n=== PROTECTION SUMMARY ===")
    for label, prot in zip(bar_labels, prots):
        print(f"  {label.replace(chr(10),' '):<25}: {prot:.1f}%")
    print(f"\n  Dex clinical benchmark: ~40% (Swain 1997)")


if __name__ == '__main__':
    make_figure('fig3_comparator.png')
