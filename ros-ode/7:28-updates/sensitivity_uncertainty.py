"""
sensitivity_uncertainty.py
==========================
Revision analyses for PONE-D-26-28903, addressing reviewer requests for:
  - global sensitivity analysis (Sobol variance decomposition + PRCC)   [R1.6, R3.12, R2.4]
  - uncertainty quantification / confidence bands                        [R2.4, R3.4]
  - goodness-of-fit metrics for the validation figure                    [R1.5, R3.5]

All analyses call ros_core_paper.py, which reproduces the published headline
numbers exactly (27.4 / 52.6 / 71.0 %, cycle-6 32.7 %, peak H2O2 6.868 nM).

Parameter ranges are multiplicative and tied to the Table 1 confidence tier:
  Strong   -> +/-25%   (base/1.25, base*1.25)
  Moderate -> +/-50%   (base/1.50, base*1.50)
  Weak     -> factor 2 (base/2,    base*2)
This makes the least-supported parameters (k_ros, k_ox, damage constants) the
widest, which is the honest choice and the one reviewers asked to see justified.

Outputs: figures in figures/revision/ and a text log in docs/revision_outputs.txt
"""
import os, time, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol as sobol_analyze
from scipy.stats import qmc

import ros_core_paper as M

RNG = np.random.default_rng(20260727)
FIGDIR = os.path.join(os.path.dirname(__file__), "..", "figures", "revision")
DOCDIR = os.path.join(os.path.dirname(__file__), "..", "docs")
os.makedirs(FIGDIR, exist_ok=True)
LOG = []
def log(s=""):
    print(s); LOG.append(str(s))

plt.rcParams.update({
    "font.family": "serif", "font.size": 10, "axes.grid": True,
    "grid.alpha": 0.25, "savefig.dpi": 300, "figure.dpi": 120,
})

# ---------------------------------------------------------------------------
# confidence tier -> multiplicative half-width
TIER = {  # from Table 1
    'k_in':'M','k_cl':'M','k_ros':'W','k_sod':'S','k_gsh':'M','Km_gsh':'M',
    'k_syn':'M','k_deg':'M','k_ox':'W','k_nrf2':'M','k_nrf2deg':'M',
    'k_gshsyn':'M','k_gshdeg':'M','k_gshsyn_basal':'M','k_gsh_turn':'M',
    'k_damage_acute':'W','k_repair_acute':'W','k_damage_chronic':'W','k_chronic_gsh':'W',
}
FACTOR = {'S':1.25,'M':1.5,'W':2.0}

def bounds_for(names):
    lo, hi = [], []
    for n in names:
        f = FACTOR[TIER[n]]; b = M.BASE[n]
        lo.append(b/f); hi.append(b*f)
    return np.array(lo), np.array(hi)

def make_p(names, x):
    p = dict(M.BASE)
    for n, v in zip(names, x):
        p[n] = float(v)
    return p

# parameters that influence peak H2O2 (exclude intervention switch + pure-damage constants)
H2O2_PARAMS = ['k_in','k_cl','k_ros','k_sod','k_gsh','Km_gsh','k_syn','k_deg',
               'k_ox','k_nrf2','k_nrf2deg','k_gshsyn','k_gshdeg','k_gshsyn_basal','k_gsh_turn']
# parameters that influence cycle-6 irreversible damage (add the two D_chronic constants)
DMG_PARAMS  = H2O2_PARAMS + ['k_damage_chronic','k_chronic_gsh']

# ===========================================================================
# 1. SOBOL on peak H2O2
# ===========================================================================
def run_sobol(names, output_fn, N, label, fname):
    lo, hi = bounds_for(names)
    problem = {'num_vars': len(names), 'names': names,
               'bounds': list(zip(lo, hi))}
    t0 = time.time()
    X = sobol_sample.sample(problem, N, calc_second_order=False, seed=1)
    Y = np.empty(X.shape[0])
    for i, row in enumerate(X):
        Y[i] = output_fn(make_p(names, row))
    Si = sobol_analyze.analyze(problem, Y, calc_second_order=False, seed=1)
    dt = time.time() - t0
    log(f"\n[Sobol: {label}]  N={N}  evals={X.shape[0]}  time={dt:.0f}s")
    order = np.argsort(Si['ST'])[::-1]
    log(f"  {'param':16s} {'S1':>8s} {'ST':>8s} {'ST_conf':>9s}")
    for k in order:
        log(f"  {names[k]:16s} {Si['S1'][k]:8.3f} {Si['ST'][k]:8.3f} {Si['ST_conf'][k]:9.3f}")
    interaction = float(np.sum(Si['ST']) - np.sum(Si['S1']))
    log(f"  sum(ST)-sum(S1) [interaction load] = {interaction:.3f}")

    # figure
    fig, ax = plt.subplots(figsize=(7, 5))
    yv = np.arange(len(names))[::-1]
    ax.barh(yv+0.18, Si['S1'][order], height=0.36, color="#3b6ea5",
            xerr=Si['S1_conf'][order], capsize=2, label="First-order $S_1$")
    ax.barh(yv-0.18, Si['ST'][order], height=0.36, color="#a53b3b",
            xerr=Si['ST_conf'][order], capsize=2, label="Total-order $S_T$")
    ax.set_yticks(yv); ax.set_yticklabels([names[k] for k in order])
    ax.set_xlabel("Sobol sensitivity index")
    ax.set_title(f"Global sensitivity (Sobol): {label}")
    ax.legend(loc="lower right", frameon=False)
    ax.axvline(0, color="k", lw=0.6)
    plt.tight_layout()
    path = os.path.join(FIGDIR, fname)
    plt.savefig(path); plt.close()
    log(f"  saved {fname}")
    return Si, order

# ===========================================================================
# 2. PRCC on cycle-6 damage (LHS)  -- cross-check, reviewer-named method
# ===========================================================================
def run_prcc(names, output_fn, N, label, fname):
    lo, hi = bounds_for(names)
    sampler = qmc.LatinHypercube(d=len(names), seed=7)
    U = sampler.random(N)
    X = qmc.scale(U, lo, hi)
    t0 = time.time()
    Y = np.array([output_fn(make_p(names, row)) for row in X])
    # PRCC: partial correlation of rank(X_j) vs rank(Y) controlling other ranks
    R = stats.rankdata(np.column_stack([X, Y]), axis=0).astype(float)
    C = np.corrcoef(R, rowvar=False)
    Cinv = np.linalg.pinv(C)
    prcc = np.empty(len(names)); pval = np.empty(len(names))
    yidx = len(names)
    for j in range(len(names)):
        pij = -Cinv[j, yidx] / np.sqrt(Cinv[j, j] * Cinv[yidx, yidx])
        prcc[j] = pij
        # t-test on PRCC
        dfree = N - len(names) - 2
        tstat = pij * np.sqrt(dfree / max(1 - pij**2, 1e-12))
        pval[j] = 2 * stats.t.sf(abs(tstat), dfree)
    dt = time.time() - t0
    log(f"\n[PRCC: {label}]  N={N}  time={dt:.0f}s  (df={dfree})")
    order = np.argsort(np.abs(prcc))[::-1]
    log(f"  {'param':18s} {'PRCC':>8s} {'p':>10s}")
    for k in order:
        log(f"  {names[k]:18s} {prcc[k]:8.3f} {pval[k]:10.2e}")

    fig, ax = plt.subplots(figsize=(7, 5))
    yv = np.arange(len(names))[::-1]
    colors = ["#a53b3b" if prcc[k] > 0 else "#3b6ea5" for k in order]
    ax.barh(yv, prcc[order], color=colors)
    ax.set_yticks(yv); ax.set_yticklabels([names[k] for k in order])
    ax.set_xlabel("Partial rank correlation coefficient (PRCC)")
    ax.set_title(f"Global sensitivity (PRCC): {label}")
    ax.axvline(0, color="k", lw=0.6); ax.set_xlim(-1, 1)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, fname)); plt.close()
    log(f"  saved {fname}")
    return prcc, X, Y

# ===========================================================================
# 3. UNCERTAINTY: full-parameter propagation to protection %
#    Key result: contrast full-parameter spread vs the published 2-parameter IQR
# ===========================================================================
def run_uq_protection(N, fname):
    names = DMG_PARAMS
    lo, hi = bounds_for(names)
    sampler = qmc.LatinHypercube(d=len(names), seed=11)
    X = qmc.scale(sampler.random(N), lo, hi)
    t0 = time.time()
    prot_s = np.empty(N); prot_d = np.empty(N); prot_c = np.empty(N)
    base_d = np.empty(N)
    for i, row in enumerate(X):
        p = make_p(names, row)
        d0 = M.simulate_patient_cycle6(p)
        base_d[i] = d0
        prot_s[i] = 100*(1 - M.simulate_patient_cycle6(M.with_sulforaphane(p))/d0)
        prot_d[i] = 100*(1 - M.simulate_patient_cycle6(M.with_dexrazoxane(p))/d0)
        prot_c[i] = 100*(1 - M.simulate_patient_cycle6(M.with_combination(p))/d0)
    dt = time.time() - t0
    def summ(a):
        return np.median(a), np.percentile(a,5), np.percentile(a,95), np.percentile(a,25), np.percentile(a,75)
    log(f"\n[UQ: cycle-6 protection, full {len(names)}-parameter uncertainty]  N={N}  time={dt:.0f}s")
    for nm, a in [("sulforaphane",prot_s),("dexrazoxane",prot_d),("combination",prot_c)]:
        med, p5, p95, p25, p75 = summ(a)
        log(f"  {nm:13s}: median {med:5.1f}%  90% CI [{p5:5.1f}, {p95:5.1f}]  IQR [{p25:.1f}, {p75:.1f}]")
    log(f"  baseline cycle-6 damage across sample: median {np.median(base_d):.1f} AU, "
        f"range [{base_d.min():.1f}, {base_d.max():.1f}]")

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    for a, c, nm in [(prot_s,"#a53b3b","Sulforaphane"),
                     (prot_d,"#2f7d4f","Dexrazoxane"),
                     (prot_c,"#5b3b8c","Combination")]:
        ax.hist(a, bins=40, alpha=0.55, color=c, label=nm)
        ax.axvline(np.median(a), color=c, lw=1.5, ls="--")
    # published 2-parameter IQR band for sulforaphane
    ax.axvspan(31.6, 33.7, color="#a53b3b", alpha=0.12)
    ax.text(32.65, ax.get_ylim()[1]*0.92, "published\n2-param IQR\n31.6–33.7%",
            ha="center", va="top", fontsize=7, color="#a53b3b")
    ax.set_xlabel("Cycle-6 damage reduction (%)")
    ax.set_ylabel("Count")
    ax.set_title("Protection under full-parameter uncertainty\n(vs published 2-parameter IQR)")
    ax.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, fname)); plt.close()
    log(f"  saved {fname}")
    return prot_s, prot_d, prot_c

# ===========================================================================
# 4. UNCERTAINTY bands on the acute H2O2 trajectory (for validation fig)
# ===========================================================================
def run_uq_trajectory(N, fname):
    names = H2O2_PARAMS
    lo, hi = bounds_for(names)
    X = qmc.scale(qmc.LatinHypercube(d=len(names), seed=13).random(N), lo, hi)
    t_end = 1440.0  # 24 h to overlay Ludke
    t_grid = np.linspace(0, t_end, 300)
    traces = np.empty((N, t_grid.size))
    for i, row in enumerate(X):
        p = make_p(names, row)
        Ke, Nr, Gs = M.steady_state(p)
        y0 = [0,0,0,Ke,Nr,Gs,0,0]
        from scipy.integrate import solve_ivp
        sol = solve_ivp(M._rhs,(0,t_end),y0,t_eval=t_grid,args=(p,p['k_in'],Gs),
                        method="LSODA",rtol=1e-6,atol=1e-9)
        h = sol.y[2]
        # normalize to model value at 60 min (Ludke normalization convention)
        h1 = np.interp(60.0, t_grid, h)
        traces[i] = 100.0 * h / max(h1, 1e-9)
    med = np.median(traces, axis=0)
    p5, p95 = np.percentile(traces, [5,95], axis=0)

    # Ludke 2017 Dox arm (digitized from Fig 2), % of t=1h control
    lud_t = np.array([1,3,6,12,24])*60.0
    lud_y = np.array([120,145,155,180,200.0])
    lud_e = np.array([8,10,10,12,12.0])  # approx SEM from figure

    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.fill_between(t_grid/60, p5, p95, color="#3b6ea5", alpha=0.25,
                    label="Model 90% band (full-param uncertainty)")
    ax.plot(t_grid/60, med, color="#1f3f66", lw=2, label="Model median")
    ax.errorbar(lud_t/60, lud_y, yerr=lud_e, fmt="o", color="#a53b3b",
                capsize=3, label="Ludke 2017 (Dox), mean$\\pm$SEM, n=5")
    ax.set_xlabel("Time (h)"); ax.set_ylabel("H$_2$O$_2$ (% of t=1 h)")
    ax.set_title("Model validation with uncertainty band")
    ax.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, fname)); plt.close()
    log(f"\n[UQ trajectory] saved {fname}")

# ===========================================================================
# 5. GOODNESS OF FIT vs Ludke (deterministic base model)
# ===========================================================================
def run_gof():
    p = dict(M.BASE)
    from scipy.integrate import solve_ivp
    Ke,Nr,Gs = M.steady_state(p); y0=[0,0,0,Ke,Nr,Gs,0,0]
    t_end=1440.0; t_grid=np.linspace(0,t_end,2000)
    sol=solve_ivp(M._rhs,(0,t_end),y0,t_eval=t_grid,args=(p,p['k_in'],Gs),
                  method="LSODA",rtol=1e-6,atol=1e-9)
    h=sol.y[2]; h1=np.interp(60.0,t_grid,h)
    model_pct = lambda th: 100.0*np.interp(th,t_grid,h)/h1
    lud_t = np.array([1,3,6,12,24])*60.0
    lud_y = np.array([120,145,155,180,200.0])
    pred  = model_pct(lud_t)

    def gof(idx, tag):
        o=lud_y[idx]; m=pred[idx]
        rmse=np.sqrt(np.mean((m-o)**2)); mae=np.mean(np.abs(m-o))
        ss_res=np.sum((o-m)**2); ss_tot=np.sum((o-np.mean(o))**2)
        r2 = 1-ss_res/ss_tot if ss_tot>0 else float('nan')
        log(f"  {tag:22s} n={len(idx)}  RMSE={rmse:5.1f}%  MAE={mae:5.1f}%  R2={r2:6.3f}")
        return rmse,mae,r2
    log("\n[Goodness of fit vs Ludke 2017 Dox arm]  (model normalized to t=1h)")
    log(f"  timepoints (h)      : {list((lud_t/60).astype(int))}")
    log(f"  observed (%)        : {list(lud_y)}")
    log(f"  model (%)           : {[round(x,1) for x in pred]}")
    gof(np.array([0,1,2]), "early phase t=1-6h")
    gof(np.arange(5),      "full window t=1-24h")

    # residual figure
    fig,(a1,a2)=plt.subplots(1,2,figsize=(10,4))
    a1.errorbar(lud_t/60,lud_y,yerr=[8,10,10,12,12],fmt="o",color="#a53b3b",capsize=3,label="Ludke (obs)")
    a1.plot(t_grid/60,model_pct(t_grid),color="#1f3f66",lw=2,label="Model")
    a1.axvspan(0,6,color="#2f7d4f",alpha=0.08); a1.text(3,95,"validation window",fontsize=7,ha="center",color="#2f7d4f")
    a1.set_xlabel("Time (h)"); a1.set_ylabel("H$_2$O$_2$ (% of t=1h)"); a1.set_title("Fit"); a1.legend(frameon=False,fontsize=8)
    res=pred-lud_y
    a2.axhline(0,color="k",lw=0.6)
    a2.stem(lud_t/60,res,linefmt="#3b6ea5",markerfmt="o",basefmt=" ")
    a2.set_xlabel("Time (h)"); a2.set_ylabel("Residual (model - obs, %)"); a2.set_title("Residuals")
    plt.tight_layout(); plt.savefig(os.path.join(FIGDIR,"figure_gof_residuals.png")); plt.close()
    log("  saved figure_gof_residuals.png")

# ===========================================================================
def append_log():
    path = os.path.join(DOCDIR, "revision_outputs.txt")
    with open(path, "a") as f:
        f.write("\n".join(LOG) + "\n")

if __name__ == "__main__":
    import sys
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"
    if stage in ("gof", "all"):
        run_gof()
    if stage in ("sobol_h2o2", "all"):
        run_sobol(H2O2_PARAMS, M.peak_h2o2, N=256, label="peak H2O2",
                  fname="figure_sobol_h2o2.png")
    if stage in ("sobol_dmg", "all"):
        run_sobol(DMG_PARAMS, M.simulate_patient_cycle6, N=128,
                  label="cycle-6 irreversible damage", fname="figure_sobol_damage.png")
    if stage in ("prcc", "all"):
        run_prcc(DMG_PARAMS, M.simulate_patient_cycle6, N=1200,
                 label="cycle-6 irreversible damage", fname="figure_prcc_damage.png")
    if stage in ("uq_prot", "all"):
        run_uq_protection(N=500, fname="figure_uq_protection.png")
    if stage in ("uq_traj", "all"):
        run_uq_trajectory(N=400, fname="figure_uq_trajectory.png")
    append_log()
