"""Figure 1: ODE state variable trajectories."""
import numpy as np, matplotlib, sys, os
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'models'))
from ode_core import run, steady_states
from parameters import PARAMS

def make_figure(output_path='fig1_validation.png'):
    t, y = run(t_end=720)
    _, Nrf2_ss, GSH_ss = steady_states(PARAMS)
    Dox, O2, H2O2, Keap1, Nrf2, GSH, D = y
    C = {'H2O2':'#2E75B6','O2':'#ED7D31','Nrf2':'#70AD47',
         'GSH':'#FFC000','D':'#C00000','Dox':'#7030A0'}
    fig = plt.figure(figsize=(14,9)); fig.patch.set_facecolor('#F8F9FA')
    gs  = gridspec.GridSpec(2,3,figure=fig,hspace=0.45,wspace=0.35)
    for i,(ax,data,col,title,ylabel) in enumerate([
        (fig.add_subplot(gs[0,0]), H2O2, C['H2O2'], 'A  H₂O₂', '[H₂O₂] (nM)'),
        (fig.add_subplot(gs[0,1]), O2,   C['O2'],   'B  O₂·⁻',  '[O₂·⁻] (nM)'),
        (fig.add_subplot(gs[1,1]), D,    C['D'],    'E  Damage', 'Damage D (AU)'),
        (fig.add_subplot(gs[1,2]), Dox,  C['Dox'],  'F  Dox',    '[Dox] (nM)'),
    ]):
        ax.plot(t, data, color=col, lw=2.5); ax.set_facecolor('white')
        ax.set_xlabel('Time (min)'); ax.set_ylabel(ylabel)
        ax.set_title(title, fontweight='bold', fontsize=10); ax.grid(True, alpha=0.3)
    ax3 = fig.add_subplot(gs[0,2]); ax3.set_facecolor('white')
    ax3.plot(t, Nrf2/Nrf2_ss, color=C['Nrf2'], lw=2.5)
    ax3.axhline(2.0, color='gray', ls=':', lw=1.2)
    ax3.axvline(123, color=C['D'], ls=':', lw=1.2, label='Deadline t=123')
    idx = np.where((Nrf2/Nrf2_ss >= 2.0)&(t>10))[0]
    if len(idx): ax3.axvline(t[idx[0]], color=C['Nrf2'], ls='--', lw=1.2,
                              label=f'Nrf2 2x at t={t[idx[0]]:.0f}')
    ax3.set_title('C  Nrf2', fontweight='bold', fontsize=10)
    ax3.set_xlabel('Time (min)'); ax3.set_ylabel('Nrf2 (fold)'); ax3.legend(fontsize=8)
    ax3.grid(True, alpha=0.3)
    ax4 = fig.add_subplot(gs[1,0]); ax4.set_facecolor('white')
    ax4.plot(t, GSH/GSH_ss*100, color=C['GSH'], lw=2.5); ax4.axhline(100, color='gray',ls=':',lw=1)
    ax4.set_title(f'D  GSH (depl {(1-GSH.min()/GSH_ss)*100:.1f}%)', fontweight='bold', fontsize=10)
    ax4.set_xlabel('Time (min)'); ax4.set_ylabel('GSH (%)')
    ax4.grid(True, alpha=0.3)
    fig.suptitle('NanoROS — ODE Validation (Ludke 2017)', fontsize=13, fontweight='bold', y=1.01)
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='#F8F9FA')
    plt.close(); print(f"Saved: {output_path}")

if __name__ == '__main__':
    make_figure('fig1_validation.png')
