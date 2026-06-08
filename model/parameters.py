"""
NanoROS — Model Parameters
===========================
All parameters with units and literature citations.
Validated ODE: Ludke et al. 2017 (J Mol Cell Cardiol)
Adult rat cardiomyocyte H2O2 dynamics under doxorubicin exposure.

Calibration notes:
  k_damage  = 0.090903  — derived analytically so D=50 AU at t=123 min
  k_dam_gsh = 0.005     — small GSH-deficit contribution
  k_syn     = 0.1       — matches original ros_model.py (Nrf2 2x at t=160 min)
  All other params match ros_model.py exactly.
"""

# ─── CORE ODE PARAMETERS ─────────────────────────────────────────────────────
PARAMS = {
    # Doxorubicin PK
    'k_in':   0.50,     # nM/min       Dox influx → [Dox]_ss = 25 nM
    'k_cl':   0.02,     # min-1        Dox clearance

    # ROS cascade
    'k_ros':  0.50,     # nM-1 min-1   O2.- production per nM Dox (Complex I/III)
    'k_sod':  5.00,     # min-1        SOD: O2.- → H2O2

    # GSH-mediated H2O2 scavenging (Michaelis-Menten)
    'k_gsh':  2.00,     # min-1        Rate constant (GSH-saturated)
    'Km_gsh': 500.0,    # nM           GSH half-saturation

    # Keap1/Nrf2 feedback
    'k_syn':      0.100,  # nM/min     Keap1 synthesis → [Keap1]_ss = 100 nM
    'k_deg':      0.001,  # min-1      Keap1 basal degradation
    'k_ox':       0.001,  # nM-1min-1  H2O2-mediated Keap1 oxidation
    'k_nrf2':     0.100,  # nM/min     Nrf2 nuclear translocation
    'k_nrf2deg':  0.010,  # nM-1min-1  Nrf2 degradation (Keap1-dependent)
    'k_gshsyn':   0.100,  # nM/min/nM  Nrf2-driven GSH synthesis per nM Nrf2
    'k_gshsyn_b': 5.000,  # nM/min     Basal GSH synthesis
    'k_gshdeg':   0.020,  # min-1      GSH oxidation per H2O2 saturation unit
    'k_gsh_turn': 0.001,  # min-1      Basal GSH turnover (GSSG export/recycling)

    # Damage accumulation — CALIBRATED to deadline t=123 min
    # k_damage derived: integral(H2O2, 0->123 min) * k_damage = 50 AU
    # Integral = 526.3 nM·min → k_damage = (50 - 0.005*GSHD_int) / 526.3
    'k_damage':  0.090903,  # AU nM-1 min-1  H2O2-driven irreversible damage
    'k_dam_gsh': 0.005,     # AU nM-1 min-1  GSH-deficit-driven damage
    'D_crit':    50.0,      # AU             Irreversibility threshold (t=123 min)
}

# ─── CARDIOLIPIN / 4-HNE MODULE ──────────────────────────────────────────────
# Extended model state variables: CL_OOH (nmol/mg), HNE (µM)
# Calibrated to Petrosillo 2003 (Biochem J): [CL-OOH]_ss = 2-8 nmol/mg
#
# NOTE: k_CL_ox corrected 28-fold downward from original estimate (0.001).
# Original produced [CL-OOH]_ss = 69 nmol/mg — 8-35x above Petrosillo range.
# Corrected value: 3.6e-5, giving [CL-OOH]_ss ≈ 5 nmol/mg.
CL_PARAMS = {
    'k_CL_ox':    3.6e-5,  # nmol mg-1 min-1 per nM H2O2  (Petrosillo 2003)
    'CL_total':   80.0,    # nmol/mg    total cardiac cardiolipin
    'k_CL_deg':   0.023,   # min-1      CL-OOH spontaneous decomposition (Antosiewicz 2006)
    'k_HNE_prod': 1.84,    # µM min-1 per nmol/mg CL-OOH  (Esterbauer 1991)
    'k_HNE_clr':  3.0,     # min-1      GSH Michael addition clearance
    # k_HNE_clr = k2_GSH(1e4 M-1s-1) × [GSH](5000 nM) × 60 = 3.0 min-1
    # GPx4 acts upstream on CL-OOH (not on 4-HNE directly)
    # GPx4 depletion by dox (Tanaka 2021, JCI Insight) increases CL_OOH_ss
    # → higher 4-HNE production per cycle (self-amplifying across cycles)
    '[4-HNE]_estimated_range': (0.034, 0.46),  # µM (90-50% GPx4 depletion)
    '[4-HNE]_model_peak':      0.12,           # µM (from extended model at ~400 min)
}

# ─── AND GATE CHEMISTRY ───────────────────────────────────────────────────────
AND_GATE = {
    # Arm A — Dioxaborolane (H2O2 temporal integrator)
    # Acid-stable, GSH-inert, fires at 0.1% shell boronate oxidation
    # Source: Bhattacharya & Chang, JACS 2014
    'k2_boronate_pH74':  156.0,  # M-1s-1  at pH 7.4
    'k2_boronate_pH80':  600.0,  # M-1s-1  at pH 8.0 (matrix) — 3.85x pH correction
    'f_crit_boronate':   0.001,  # 0.1% → amphiphilic shell disassembly
    'arm_A_fire_t':      22.0,   # min    at pH 8.0, f_crit=0.1%, ODE H2O2 trajectory

    # Arm B Branch A — Methoxyamine ([4-HNE] ≥ 0.075 µM)
    # Source: Dirksen & Dawson, Bioconjug Chem 2008
    'k2_methoxyamine':   300.0,  # M-1s-1
    'arm_B_fire_t_A':    120.0,  # min    at [4-HNE] = 0.075 µM
    'protection_A':      0.83,   # 83%    per-cycle (ODE model, ceiling)

    # Arm B Branch B — Benzyloxyamine (0.025 ≤ [4-HNE] < 0.075 µM)
    # Source: Dirksen & Dawson, Bioconjug Chem 2008
    'k2_benzyloxyamine': 1000.0, # M-1s-1
    'arm_B_fire_t_B':    (49, 108),    # min range at [4-HNE] = 0.075-0.025 µM
    'protection_B':      (0.83, 0.87),

    # Arm B Branch C — acrolein sensing OR Arm A only ([4-HNE] < 0.025 µM)
    # Updated from MDA (item 3): MDA bifunctionality crosslinking risk >> 1%
    # C1: if [acrolein] >= 0.1 µM → O-benzyloxyamine for acrolein (monofunctional)
    # C2: if [acrolein] < 0.1 µM → Arm A only, restrict to cycle 4+
    'arm_B_fire_t_C1':   27.0,   # min    benzyloxyamine at [acrolein] = 0.1 µM
    'arm_B_fire_t_C2':   22.0,   # min    Arm A only (boronate, H2O2)

    # Oxime formation threshold (updated from 5% to prevent prolonged-exercise FP)
    'f_crit_oxime':      0.15,   # 15% oxime formation for shell disassembly

    # Damage deadline (from ODE model)
    'damage_deadline':   123.0,  # min
}

# ─── PAYLOAD ─────────────────────────────────────────────────────────────────
# MnTE-2-PyP selected payload.
# Eliminated: MitoQ (no viable therapeutic window at any ΔΨm),
#             Mn(II) cyclen (matrix-corrected kcat/Km gives k_sc=1.5 min-1,
#             3x below requirement — phosphate competition per Kimura 2007),
#             catalase (240 kDa, cannot cross OMM via VDAC),
#             EUK-207 (catalase activity 400x lower than native)
PAYLOAD = {
    'name':              'MnTE-2-PyP',
    'MW':                900,        # Da — VDAC-permeable (< 5 kDa cutoff)
    'charge':            +5,         # at pH 7.4
    'kcat_Km_O2':        7e8,        # M-1s-1  Batinic-Haberle 2009/2012
    'EE_target':         0.05,       # 5% encapsulation efficiency (conservative)
    'k_sc_O2_at_EE':     23515.0,    # min-1   at 5% EE (4,700x above k_sod floor)
    'pro_ox_ceiling':    None,       # No semiquinone redox cycling side reaction
    'delta_psi_dep':     False,      # ΔΨm-independent (passive VDAC diffusion)
    'intervention':      'O2.-',     # Upstream interception before SOD → H2O2
    'NOAEL_mice_mgkg':   10.0,       # mg/kg/day × 18 days IV (Gad 2013, Int J Toxicol)
    'NOAEL_monkey_mgkg': 5.0,        # mg/kg/day × 14 days IV (Gad 2013)
}

# ─── MITOQ FAILURE ANALYSIS ──────────────────────────────────────────────────
# Primary citation: Doughan & Dikalov, Antioxid Redox Signal 2007;9(11):1825-1836
# Pro-oxidant ceiling ~1 µM intramitochondrial
# At normal ΔΨm (150 mV): 275-fold accum × 100 nM plasma = 27.5 µM → 27x above ceiling
# At 40% residual ΔΨm (60 mV): 9-fold accum → 0.9 µM → at ceiling, reduced efficacy
# MitoQ cannot be dosed to be simultaneously safe and effective across heterogeneous ΔΨm
MITOQ = {
    'pro_ox_ceiling_uM':  1.0,
    'dpsi_normal_mV':   150.0,
    'dpsi_dox_mV':       60.0,   # 40% residual (Pacher 2003; Bhatt 2013)
    'accum_normal':     275.0,   # fold (Nernst at 150 mV, 37°C)
    'accum_dox':          9.0,   # fold (Nernst at 60 mV)
    'plasma_conc_nM':   100.0,   # nM clinical therapeutic plasma concentration
}

# ─── DELIVERY CASCADE ────────────────────────────────────────────────────────
# At 20 mg/kg IV, CRPPR + ICAM-1 + CD47 + AT1R dual targeting
# Human cardiac accumulation: ~2% (estimated)
# CRPPR mouse: 4.2% (Dvir 2011) × 14x scale-down = 0.3-1.0% human
# × ICAM-1 4x improvement (Bhattacharya 2012) = 1.2-4.0%
# × CD47 ~1.5x = conservative estimate 2%
DELIVERY = {
    'dose_mg_per_kg':                  20.0,
    'body_weight_kg':                  70.0,
    'NP_diameter_nm':                  100.0,
    'cardiac_pct_estimate':            0.02,   # 2% (CRPPR + ICAM-1 + CD47)
    'AT1R_cardiomyocyte_selectivity':  0.60,   # 60% cardiomyocyte vs other cardiac cells
    'endosomal_escape_DOPE_CHEMS':     0.10,   # 10% (pH-triggered hexagonal transition)
    'n_cardiomyocytes_human':          7.5e9,
    'NPs_per_cell_at_2pct':            357,    # from delivery cascade table (Section 7)
    'V_matrix_per_cell_L':             12e-12, # 12 pL (5000 mito × 0.5 fL/mito)
    'safety_margin_k_sc':              4700,   # × above k_sod floor
}
