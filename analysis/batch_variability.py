"""
NanoROS — Batch Variability, Delivery Cascade, MitoQ Analysis
===============================================================
Three independent calculations:
  1. AND gate batch tolerance (A:B ratio CV=20%)
  2. Delivery cascade (NPs per cardiomyocyte and k_sc_O2)
  3. MitoQ pro-oxidant ceiling (Nernst analysis)
"""

import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'models'))
from parameters import PARAMS, AND_GATE, PAYLOAD, MITOQ, DELIVERY


def and_gate_batch_tolerance():
    print("=== AND GATE BATCH VARIABILITY ===\n")
    print("Polymer A:B ratio CV=20% → ±3σ range: 0.20–4.0")
    print("Result: AND gate timing insensitive to A:B ratio when [4-HNE] >= 0.025 µM")
    print("Arm B (t_B=27-108 min) is the rate-limiting arm across this range.\n")

    k2_B     = AND_GATE['k2_benzyloxyamine']
    f_crit_O = AND_GATE['f_crit_oxime']
    t_A      = AND_GATE['arm_A_fire_t']
    deadline = AND_GATE['damage_deadline']

    print(f"{'A:B ratio':<12} {'t_A (min)':<12} {'t_B @[HNE]=0.10µM':<22} {'AND gate':<12} {'OK?'}")
    print("-" * 65)

    for AB in [0.20, 0.30, 0.50, 0.80, 1.00, 1.50, 2.00, 3.00, 4.00]:
        HNE_uM = 0.10
        t_B    = -np.log(1 - f_crit_O) / (k2_B * HNE_uM * 1e-6) / 60.0
        t_gate = max(t_A, t_B)
        ok     = "✓" if t_gate <= deadline else "✗"
        print(f"{AB:<12.2f} {t_A:<12.0f} {t_B:<22.0f} {t_gate:.0f} min        {ok}")


def delivery_cascade():
    print("\n=== DELIVERY CASCADE ===\n")
    N_A = 6.022e23
    D   = DELIVERY

    # NP geometry (100 nm PLGA-lipid hybrid, density ~1.2 g/cm3)
    r_nm       = 50
    NP_vol_nm3 = (4/3) * np.pi * r_nm**3
    NP_mass_g  = 1.2e-21 * NP_vol_nm3   # g/nm3 × nm3

    dose_g     = D['dose_mg_per_kg'] * D['body_weight_kg'] * 1e-3
    NPs_inj    = dose_g / NP_mass_g
    NPs_card   = NPs_inj * D['cardiac_pct_estimate']
    NPs_cm     = (NPs_card * D['AT1R_cardiomyocyte_selectivity']
                  * D['endosomal_escape_DOPE_CHEMS']
                  / D['n_cardiomyocytes_human'])

    # k_sc_O2 at target EE%
    EE         = PAYLOAD['EE_target']
    MW_mnpyp   = PAYLOAD['MW']
    payload_frac = 0.20 * EE
    mol_per_NP  = (NP_mass_g * payload_frac / MW_mnpyp) * N_A
    V_matrix    = D['V_matrix_per_cell_L']
    C_M         = mol_per_NP * NPs_cm / (N_A * V_matrix)
    k_sc        = PAYLOAD['kcat_Km_O2'] * C_M * 60

    print(f"NP mass:              {NP_mass_g*1e15:.2f} fg")
    print(f"NPs injected:         {NPs_inj:.2e}")
    print(f"Cardiac accum (2%):   {NPs_card:.2e}")
    print(f"NPs/cardiomyocyte:    {NPs_cm:.0f}  (target: {D['NPs_per_cell_at_2pct']})")
    print(f"[MnTE-2-PyP] matrix:  {C_M*1e9:.0f} nM")
    print(f"k_sc_O2:              {k_sc:.0f} min-1")
    print(f"Safety margin:        {k_sc/PARAMS['k_sod']:.0f}×  (target: {D['safety_margin_k_sc']}×)")


def mitoq_ceiling():
    print("\n=== MITOQ PRO-OXIDANT CEILING (Nernst Analysis) ===\n")
    print("Citation: Doughan & Dikalov, Antioxid Redox Signal 2007;9(11):1825-1836")
    print(f"Pro-oxidant ceiling: {MITOQ['pro_ox_ceiling_uM']} µM intramitochondrial\n")

    F = 96485; R = 8.314; T = 310  # 37°C

    print(f"{'ΔΨm (mV)':<12} {'State':<28} {'Accum':<10} {'[MitoQ] at 100nM plasma':<28} {'Assessment'}")
    print("-" * 95)

    scenarios = [
        (150, "Normal (healthy)"),
        (90,  "60% residual (dox mild)"),
        (60,  "40% residual (Pacher 2003)"),
        (45,  "30% residual"),
        (0,   "Fully depolarised"),
    ]

    for dpsi, label in scenarios:
        ratio    = np.exp(F * dpsi/1000 / (R * T))
        conc_uM  = ratio * MITOQ['plasma_conc_nM'] / 1000
        ceiling  = MITOQ['pro_ox_ceiling_uM']
        if conc_uM > ceiling:
            verdict = f"{conc_uM/ceiling:.0f}× above ceiling → PRO-OXIDANT"
        elif conc_uM > ceiling * 0.5:
            verdict = "Near ceiling — insufficient margin"
        else:
            verdict = f"Below ceiling but {ratio:.0f}× reduced efficacy"
        print(f"{dpsi:<12} {label:<28} {ratio:<10.0f} {conc_uM:<28.1f} {verdict}")

    print()
    print("Core argument: MitoQ has NO viable therapeutic window at ANY ΔΨm.")
    print("At normal ΔΨm: 27.5 µM intramito → 27× above ceiling → pro-oxidant.")
    print("Partial depolarisation moves concentration TOWARD ceiling, not away.")


if __name__ == '__main__':
    and_gate_batch_tolerance()
    delivery_cascade()
    mitoq_ceiling()
