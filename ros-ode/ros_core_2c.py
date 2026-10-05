"""
ros_core_2c.py  --  Two-compartment (mitochondrial / cytosolic) ROS ODE scaffold
=================================================================================

DOWNSTREAM EXPLORATORY WORK. Not part of the PLOS ONE R1 submission.
Branch: twocompartment-dev   Dir: scripts/twocompartment/
The frozen 1C model (ros_core_paper.py, tag plosone-R1-submitted-2026-09-12) is a
READ-ONLY dependency. This module IMPORTS it and reconciles against it; it never
edits it.

WHY THIS MODEL EXISTS
---------------------
The 1C paper model has a single ROS pool, so it cannot distinguish a mito-targeted
intervention from a whole-cell one. The mechanistic content of the 2C model is a
SPATIAL MISMATCH between where ROS is produced and where the antioxidant response
lives:
  * Dox redox-cycles at Complex I -> superoxide is generated in the MATRIX.
  * O2*- is charged and membrane-IMPERMEANT: it cannot leave the matrix; it is
    dismutated locally by MnSOD to H2O2.
  * H2O2 is the ONLY species that couples the compartments (diffusion + AQP8).
  * Keap1-Nrf2 sensing and GSH synthesis are CYTOSOLIC; mito GSH is import-only.
So the antioxidant response senses H2O2 that has already leaked out of the
mitochondrion. That displacement is why a matrix-localized scavenger could beat a
whole-cell antioxidant, and it is exactly what the 1C model averages away.

RECONCILED against ros_core_paper.py (2026-09-12):
  * GSH consumption uses gsh_sat*H2O2 (matches _rhs), NOT the earlier reconstruction.
  * Keap1 eq replicates k_ox_t = k_ox + k_ox_sustained plus the -k_ox_sustained*Keap1
    term, so interventions (sulforaphane=0.002) transfer.
  * negativity clamps on GSH/Keap1/Nrf2 as in _rhs.
  * Default params mirror ros_core_paper.BASE (UNCALIBRATED: k_cl=0.02, k_gsh=2.0).
    The manuscript figures use the calibrated overlay k_cl=0.0267, k_gsh=3.0 -- apply
    it with Params2C.calibrated(). BASE alone gives Dox_ss=25, H2O2_ss=6.87 (NOT the
    manuscript's 18.7/3.4).

CONFIDENCE FLAGS
----------------
(a) established: reaction forms now mirror the frozen module; superoxide impermeant;
    H2O2 couples; mass-conserving transport (test #1); exact reduction to the frozen
    module in the symmetric limit (test #2, compares to ros_core_paper directly).
(c) speculative / UNCONSTRAINED placeholders -- NOT fitted:
    - P_h2o2 (transbilayer H2O2 transport): MOST load-bearing; not identifiable from
      bulk H2O2; literature spans orders of magnitude.
    - f_mito (mito volume fraction ~0.35 cardiomyocyte; VERIFY morphometry).
    - k_sod_m/k_sod_c split; GSH import kinetics.

KNOWN TRAP: 2C is non-identifiable from bulk H2O2. Do NOT fit to Ludke bulk H2O2 and
report compartment params. Structural hypothesis generator until compartment-resolved
data (mito-HyPer7/MitoPY1 + cytosolic HyPer under dox) exists -- the go/no-go experiment.

OPEN: what did Ludke's H2O2 assay measure (extracellular Amplex Red / whole-cell DCF /
targeted)? Sets w_m/w_c in observable_h2o2(). Verify doi:10.1371/journal.pone.0179452.
"""
from __future__ import annotations
from dataclasses import dataclass
import os, sys
import numpy as np
from scipy.integrate import solve_ivp

# --- import the frozen 1C module as a read-only reference (path-robust) ----------
for _p in (os.path.dirname(__file__), os.path.join(os.path.dirname(__file__), ".."),
           os.path.join(os.path.dirname(__file__), "..", ".."),
           os.path.join(os.path.dirname(__file__), "..", "..", "scripts"),
           "/tmp", "/mnt/user-data/uploads"):
    if _p and _p not in sys.path:
        sys.path.insert(0, _p)
try:
    import ros_core_paper as rp
except Exception:
    rp = None  # tests that need it will skip with a clear message


# ----------------------------------------------------------------------------------
# 2C parameters. Defaults mirror ros_core_paper.BASE (uncalibrated).
# ----------------------------------------------------------------------------------
@dataclass
class Params2C:
    # --- shared kinetics (mirror ros_core_paper.BASE) ---
    k_in: float = 0.5
    k_cl: float = 0.02          # BASE (uncalibrated). calibrated() -> 0.0267
    k_ros: float = 0.5
    k_sod: float = 5.0
    k_gsh: float = 2.0          # BASE (uncalibrated). calibrated() -> 3.0
    k_syn: float = 0.02
    k_deg: float = 0.001
    k_ox: float = 0.001
    k_ox_sustained: float = 0.0 # sulforaphane sets 0.002
    k_nrf2: float = 0.1
    k_nrf2deg: float = 0.01
    k_gshsyn: float = 0.1
    k_gshdeg: float = 0.02
    k_gshsyn_basal: float = 5.0
    k_gsh_turn: float = 0.001
    Km_gsh: float = 500.0
    k_damage_acute: float = 0.1
    k_repair_acute: float = 0.001
    k_damage_chronic: float = 0.0005
    k_chronic_gsh: float = 0.0002
    # --- compartment structure (UNCONSTRAINED placeholders) ---
    f_mito: float = 0.35        # (c) mito volume fraction; VERIFY morphometry
    P_h2o2: float = 0.5         # (c) H2O2 transport coeff (permeability x area), 1/min
    k_gsh_import: float = 0.2   # (c) facilitated GSH equilibration cyto->mito
    k_sod_m: float = 5.0        # (c) MnSOD (matrix)
    k_sod_c: float = 5.0        # (c) SOD1 (cytosol)
    f_mito_ros: float = 1.0     # fraction of superoxide GENERATION in mito (mechanistic=1)
    # --- switches ---
    mito_has_synthesis: bool = False  # mechanistic: mito GSH is import-only
    o2_transport: bool = False        # PHYSICS: O2- impermeant -> keep False
    # --- superoxide-scavenging payload (the mito-targeting intervention) ---
    # Substrate (superoxide) is matrix-confined, so targeting is an ACCESS advantage:
    # matrix [scavenger] = scav_dose * scav_enrich (E). E=1 untargeted (uniform),
    # E=1/f_mito volume-partition ceiling, E>> affinity/potential-targeted.
    scav_kcat: float = 0.0            # scavenger rate per nM scavenger (/min); 0 = off
    scav_dose: float = 0.0            # total cellular scavenger, uniform-equivalent (nM)
    scav_enrich: float = 1.0          # matrix enrichment factor E
    # damage driver: 'h2o2' well-mixed (legacy, no localization signal) OR
    # 'o2' matrix-confined superoxide (mechanistic; the species targeting can act on).
    # k_damage_o2 is ARBITRARY (AU) -- ultimately calibrate to matrix 4-HNE (Experiment 0).
    damage_driver: str = "h2o2"
    k_damage_o2: float = 0.0005
    # alternative payload: matrix-targeted GSH restorer (offsets import-limited deficit).
    # delivered at the same enrichment E (scav_enrich). single parameter = input rate (nM/min).
    gsh_support: float = 0.0

    @property
    def V_m(self): return self.f_mito
    @property
    def V_c(self): return 1.0 - self.f_mito

    def calibrated(self):
        from dataclasses import replace
        return replace(self, k_cl=0.0267, k_gsh=3.0)

    @classmethod
    def from_base(cls, base: dict):
        """Build from a ros_core_paper-style dict (shared keys only)."""
        fields = {f: base[f] for f in base if f in cls.__dataclass_fields__}
        return cls(**fields)

# 2C state order:
#  [Dox, O2_m, H2O2_m, GSH_m, O2_c, H2O2_c, GSH_c, Keap1, Nrf2, D_acute, D_chronic]
S2 = dict(Dox=0, O2m=1, H2O2m=2, GSHm=3, O2c=4, H2O2c=5, GSHc=6,
          Keap1=7, Nrf2=8, Da=9, Dc=10)

def _gsh_ss(p: Params2C):
    Keap1 = p.k_syn / p.k_deg
    Nrf2  = p.k_nrf2 / (p.k_nrf2deg * Keap1)
    return (p.k_gshsyn_basal + p.k_gshsyn * Nrf2) / p.k_gsh_turn, Keap1, Nrf2

def _transport(cm, cc, P, Vm, Vc):
    """Mass-conserving transfer m->c. d/dt(Vm*cm + Vc*cc)=0 for any P,Vm,Vc."""
    J = P * (cm - cc)
    return -J / Vm, +J / Vc

def rhs_2c(t, y, p: Params2C, kin_rate=None, GSH_ss=None):
    (Dox, O2m, H2O2m, GSHm, O2c, H2O2c, GSHc, Keap1, Nrf2, Da, Dc) = y
    Vm, Vc = p.V_m, p.V_c
    if kin_rate is None: kin_rate = p.k_in
    if GSH_ss is None:   GSH_ss = _gsh_ss(p)[0]
    GSHm = max(GSHm, 0.0); GSHc = max(GSHc, 0.0)
    Keap1 = max(Keap1, 0.0); Nrf2 = max(Nrf2, 0.0)

    dDox = kin_rate - p.k_cl * Dox

    gen = p.k_ros * Dox                        # total O2- generation (moles/min)
    gen_O2m = p.f_mito_ros * gen / Vm
    gen_O2c = (1.0 - p.f_mito_ros) * gen / Vc
    dO2m = gen_O2m - p.k_sod_m * O2m
    dO2c = gen_O2c - p.k_sod_c * O2c
    # matrix-targeted superoxide scavenger: acts only where the substrate is (matrix).
    # cytosolic action is omitted because O2c ~ 0 (superoxide impermeant, matrix-generated).
    if p.scav_kcat > 0.0:
        dO2m -= p.scav_kcat * (p.scav_dose * p.scav_enrich) * O2m
    if p.o2_transport:
        tm, tc = _transport(O2m, O2c, p.P_h2o2, Vm, Vc); dO2m += tm; dO2c += tc

    sat_m = GSHm / (p.Km_gsh + GSHm)
    sat_c = GSHc / (p.Km_gsh + GSHc)
    tm, tc = _transport(H2O2m, H2O2c, p.P_h2o2, Vm, Vc)
    dH2O2m = p.k_sod_m * O2m - p.k_gsh * sat_m * H2O2m + tm
    dH2O2c = p.k_sod_c * O2c - p.k_gsh * sat_c * H2O2c + tc

    k_ox_t = p.k_ox + p.k_ox_sustained
    dKeap1 = (p.k_syn - p.k_deg * Keap1 - k_ox_t * H2O2c * Keap1
              - p.k_ox_sustained * Keap1)              # senses CYTOSOLIC H2O2
    dNrf2  = p.k_nrf2 - p.k_nrf2deg * Keap1 * Nrf2

    syn_c = p.k_gshsyn_basal + p.k_gshsyn * Nrf2       # synthesis is cytosolic
    syn_m = syn_c if p.mito_has_synthesis else 0.0
    gm, gc = _transport(GSHm, GSHc, p.k_gsh_import, Vm, Vc)
    dGSHm = syn_m - p.k_gshdeg * sat_m * H2O2m - p.k_gsh_turn * GSHm + gm
    if p.gsh_support > 0.0:                            # matrix-targeted GSH restorer payload
        dGSHm += p.gsh_support * p.scav_enrich         # delivered at matched enrichment E
    dGSHc = syn_c - p.k_gshdeg * sat_c * H2O2c - p.k_gsh_turn * GSHc + gc

    gsh_def_m = max(GSH_ss - GSHm, 0.0)                # damage keyed to MITO load
    if p.damage_driver == "o2":
        # matrix-confined superoxide (+ its local products: Fe-S inactivation, matrix
        # lipid peroxidation -> 4-HNE). This is the species a matrix-targeted scavenger
        # can actually reduce, because it is NOT well-mixed (unlike H2O2).
        dDa = p.k_damage_acute * O2m - p.k_repair_acute * Da
        dDc = p.k_damage_o2 * O2m + p.k_chronic_gsh * gsh_def_m
    else:
        dDa = p.k_damage_acute * H2O2m - p.k_repair_acute * Da
        dDc = p.k_damage_chronic * H2O2m + p.k_chronic_gsh * gsh_def_m
    return [dDox, dO2m, dH2O2m, dGSHm, dO2c, dH2O2c, dGSHc, dKeap1, dNrf2, dDa, dDc]

def ic_2c(p: Params2C):
    GSH, Keap1, Nrf2 = _gsh_ss(p)
    return [0.0, 0.0, 0.0, GSH, 0.0, 0.0, GSH, Keap1, Nrf2, 0.0, 0.0]

def observable_h2o2(y2c, p: Params2C, mode="volume"):
    H2O2m = y2c[S2['H2O2m']]; H2O2c = y2c[S2['H2O2c']]
    if mode == "volume": return p.V_m * H2O2m + p.V_c * H2O2c
    if mode == "cyto":   return H2O2c
    raise ValueError(mode)


# ----------------------------------------------------------------------------------
# SELF-TESTS. #2 now compares against the FROZEN module directly.
# ----------------------------------------------------------------------------------
def _integrate(rhs, y0, p, t_end, n=2000, **kw):
    t = np.linspace(0, t_end, n)
    sol = solve_ivp(rhs, (0, t_end), y0, t_eval=t, args=(p,)+kw.pop("args", tuple()),
                    method="LSODA", rtol=1e-8, atol=1e-10, **kw)
    assert sol.success, sol.message
    return t, sol.y

def test_h2o2_mass_conservation():
    p = Params2C(P_h2o2=0.8)
    for name in list(p.__dataclass_fields__):
        if name.startswith("k_"):
            setattr(p, name, 0.0)
    y0 = [0.0]*len(S2); y0[S2['H2O2m']] = 10.0
    _, Y = _integrate(rhs_2c, y0, p, 500.0, args=(None, 0.0))
    total = p.V_m*Y[S2['H2O2m']] + p.V_c*Y[S2['H2O2c']]
    drift = np.max(np.abs(total - total[0]))/total[0]
    assert drift < 1e-6, f"mass drift {drift:.1e}"
    assert abs(Y[S2['H2O2m']][-1]-Y[S2['H2O2c']][-1]) < 1e-3, "no equilibration"
    return drift

def test_reduces_to_frozen_paper():
    """Symmetric 2C must reproduce ros_core_paper.acute_run(BASE) exactly."""
    if rp is None:
        return None
    p2 = Params2C.from_base(rp.BASE)
    p2.k_sod_m = p2.k_sod_c = rp.BASE['k_sod']
    p2.f_mito_ros = p2.f_mito          # -> equal concentration generation
    p2.mito_has_synthesis = True       # symmetric: mito synthesis on
    t1, Y1 = rp.acute_run(rp.BASE, t_end=600.0, n_eval=2000)
    t2, Y2 = _integrate(rhs_2c, ic_2c(p2), p2, 600.0)
    bulkH = p2.V_m*Y2[S2['H2O2m']] + p2.V_c*Y2[S2['H2O2c']]
    bulkG = p2.V_m*Y2[S2['GSHm']]  + p2.V_c*Y2[S2['GSHc']]
    eH = np.max(np.abs(bulkH - Y1[2]))/np.max(Y1[2])
    eG = np.max(np.abs(bulkG - Y1[5]))/np.max(Y1[5])
    eDc = abs(Y2[S2['Dc']][-1]-Y1[7,-1])/max(Y1[7,-1],1e-9)
    assert eH < 1e-4 and eG < 1e-4 and eDc < 1e-3, (eH,eG,eDc)
    return eH, eG, eDc

def demo_mechanistic():
    p = Params2C(f_mito_ros=1.0, o2_transport=False, mito_has_synthesis=False).calibrated()
    _, Y = _integrate(rhs_2c, ic_2c(p), p, 600.0)
    m, c = Y[S2['H2O2m']][-1], Y[S2['H2O2c']][-1]
    return dict(H2O2_mito=m, H2O2_cyto=c, ratio=m/c if c else float('inf'),
                bulk_observable=observable_h2o2(Y[:, -1], p, "volume"), P_h2o2=p.P_h2o2)

if __name__ == "__main__":
    print("="*70); print("ros_core_2c.py -- structural self-tests"); print("="*70)
    d = test_h2o2_mass_conservation()
    print(f"[PASS] H2O2 mass conservation under transport-only (drift {d:.1e})")
    r = test_reduces_to_frozen_paper()
    if r is None:
        print("[SKIP] reduces-to-frozen: ros_core_paper not importable")
    else:
        eH,eG,eDc = r
        print(f"[PASS] reduces to ros_core_paper.acute_run(BASE) "
              f"(H2O2 {eH:.1e}, GSH {eG:.1e}, D_chronic {eDc:.1e})")
    print("-"*70)
    print("Illustrative mechanistic run (calibrated overlay, PLACEHOLDER P_h2o2):")
    for k,v in demo_mechanistic().items():
        print(f"    {k:16s} = {v:.4g}" if isinstance(v,float) else f"    {k:16s} = {v}")
    print("-"*70)
    print("ratio is a FUNCTION of P_h2o2 (unconstrained) -- hypothesis, not result.")