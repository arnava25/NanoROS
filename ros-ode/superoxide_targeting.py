"""
superoxide_targeting.py -- the mito-targeting question, on the correct species
==============================================================================
Downstream exploratory. Imports ros_core_2c. Saves a figure; no window (Agg).

CONTEXT (why this replaces the P_h2o2 sweep):
The P_h2o2 anchor showed H2O2 equilibrates across the inner membrane in ~ms
(P ~ 1e4 /min, deep well-mixed), so there is NO matrix H2O2 gradient to exploit.
But superoxide is charged and membrane-impermeant (P < 1e-7 cm/s), so it IS
matrix-confined. A superoxide-scavenging payload therefore only acts in the matrix
(its substrate exists nowhere else), and targeting is an ACCESS advantage measured
by the matrix enrichment factor E of the scavenger.

THE QUESTION:
  (A) how much matrix-damage reduction vs scavenger enrichment E?
  (B) at what membrane potential does a FIXED-E (potential-independent) particle
      overtake MitoQ, whose enrichment collapses as mitochondria depolarize?

CONFIDENCE:
(a) established: superoxide is matrix-confined; MitoQ enrichment is Nernstian and
    collapses toward 1x as dPsi -> 0; volume-partition ceiling for a potential-
    independent particle is 1/f_mito ~ 2.9x.
(c) arbitrary SCALING (flagged): the scavenger potency (scav_kcat*scav_dose) sets
    WHERE on the E axis you get a given reduction. Here it is normalized to endogenous
    MnSOD capacity (relative capacity 1.0 at E=1). The SHAPES and the CROSSOVER dPsi
    are the robust content; the absolute reduction magnitude is not.
    Also arbitrary: k_damage_o2 (matrix-superoxide damage, AU) -- calibrate to 4-HNE (Exp0).
"""
import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt

import ros_core_2c as M
from ros_core_2c import Params2C, rhs_2c, ic_2c, S2, _integrate

T_END   = 800.0
F_MITO  = 0.35
E_VOL   = 1.0 / F_MITO                 # volume-partition ceiling (potential-independent)
E_AFFIN = 10.0                         # example affinity/receptor-targeted enrichment
SCAV_REL_CAPACITY = 1.0                # scavenger capacity at E=1, relative to MnSOD (ARBITRARY)

def base_params(scav_enrich=1.0, scav_on=True):
    p = Params2C(f_mito_ros=1.0, o2_transport=False, mito_has_synthesis=False,
                 damage_driver="o2").calibrated()
    p.P_h2o2 = 1e4                     # anchored: H2O2 well-mixed
    if scav_on:
        p.scav_dose = 1.0
        p.scav_kcat = SCAV_REL_CAPACITY * p.k_sod_m   # capacity = rel*MnSOD per unit E
        p.scav_enrich = scav_enrich
    return p

def endpoint(p):
    _, Y = _integrate(rhs_2c, ic_2c(p), p, T_END)
    return Y[S2['Dc']][-1], Y[S2['O2m']][-1]

def mitoq_enrichment(dpsi_mV):
    """Nernstian matrix accumulation of a lipophilic cation vs membrane potential.
    dPsi is negative (matrix-negative). ~10-fold per 61.5 mV."""
    return 10.0 ** (-dpsi_mV / 61.5)

def main():
    Dc_base, O2m_base = endpoint(base_params(scav_on=False))

    # --- reduction(E) curve ---
    E_grid = np.logspace(np.log10(0.3), np.log10(1000), 90)
    red = np.empty_like(E_grid); o2m = np.empty_like(E_grid)
    for i, E in enumerate(E_grid):
        Dc, O2m = endpoint(base_params(scav_enrich=E))
        red[i] = 100 * (1 - Dc / Dc_base); o2m[i] = O2m
    red_of_E = lambda E: np.interp(E, E_grid, red)

    # --- MitoQ vs potential-independent particle across dPsi ---
    dpsi = np.linspace(0, -180, 120)
    E_mq = mitoq_enrichment(dpsi)
    red_mq  = red_of_E(E_mq)
    red_vol = red_of_E(E_VOL) * np.ones_like(dpsi)      # flat: potential-independent
    red_aff = red_of_E(E_AFFIN) * np.ones_like(dpsi)
    # crossover dPsi where MitoQ falls below the volume-targeted particle
    cross = None
    for i in range(1, len(dpsi)):
        if (red_mq[i-1]-red_vol[i-1])*(red_mq[i]-red_vol[i]) < 0:
            cross = dpsi[i]; break

    # ---- figure ----
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 4.8))

    ax1.plot(E_grid, red, color="#1a9850", lw=2)
    for E, lab, col in ((1.0, "untargeted (E=1)", "#666"),
                        (E_VOL, f"volume ceiling (E={E_VOL:.1f})", "#762a83"),
                        (E_AFFIN, f"affinity (E={E_AFFIN:.0f})", "#2166ac")):
        ax1.axvline(E, color=col, ls="--", lw=1)
        ax1.annotate(lab, (E, 4), rotation=90, fontsize=7, color=col,
                     va="bottom", ha="right")
    ax1.set_xscale("log"); ax1.set_xlabel("scavenger matrix enrichment  E")
    ax1.set_ylabel("cycle-end matrix damage reduction (%)")
    ax1.set_title("A. Damage reduction vs targeting enrichment")
    ax1.grid(alpha=.3, which="both"); ax1.set_ylim(0, 100)

    ax2.plot(-dpsi, red_mq, color="#b2182b", lw=2, label="MitoQ (Nernstian E)")
    ax2.plot(-dpsi, red_vol, color="#762a83", lw=2, ls="--",
             label=f"dPsi-independent, volume (E={E_VOL:.1f})")
    ax2.plot(-dpsi, red_aff, color="#2166ac", lw=2, ls=":",
             label=f"dPsi-independent, affinity (E={E_AFFIN:.0f})")
    if cross is not None:
        ax2.axvline(-cross, color="0.4", lw=0.8)
        ax2.annotate(f"MitoQ falls below\nvolume particle\nat dPsi={cross:.0f} mV",
                     (-cross, 20), fontsize=7.5, xytext=(6, 0),
                     textcoords="offset points", color="0.3")
    ax2.set_xlabel("mitochondrial depolarization  |dPsi|  (mV)   [healthy ~160 -> diseased ~0]")
    ax2.set_ylabel("matrix damage reduction (%)")
    ax2.set_title("B. Potential-independent particle vs MitoQ")
    ax2.invert_xaxis()                 # healthy (high |dPsi|) on left, diseased on right
    ax2.legend(fontsize=8, loc="lower left"); ax2.grid(alpha=.3); ax2.set_ylim(0, 100)

    fig.suptitle("Mito-targeting on the trapped superoxide pool "
                 "(potency scaling is ARBITRARY; shapes + crossover are the result)",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig("figure_superoxide_targeting.png", dpi=150); plt.close(fig)

    # ---- console ----
    print("=" * 70)
    print("Superoxide-targeting sweep (damage_driver='o2', calibrated, H2O2 well-mixed)")
    print("=" * 70)
    print(f"baseline (no scavenger): D_chronic={Dc_base:.3f}, matrix O2={O2m_base:.3f} nM")
    print(f"reduction at E=1   (untargeted)         : {red_of_E(1.0):5.1f}%")
    print(f"reduction at E={E_VOL:.1f} (volume ceiling)       : {red_of_E(E_VOL):5.1f}%")
    print(f"reduction at E={E_AFFIN:.0f}  (affinity example)     : {red_of_E(E_AFFIN):5.1f}%")
    print(f"reduction at E=400 (MitoQ, healthy dPsi) : {red_of_E(400):5.1f}%")
    print("-" * 70)
    print("MitoQ vs potential-independent particle:")
    for psi in (-160, -120, -80, -40, -20):
        Emq = mitoq_enrichment(psi)
        print(f"  dPsi={psi:4d} mV: MitoQ E={Emq:7.1f} -> {red_of_E(Emq):5.1f}% ; "
              f"volume particle (E={E_VOL:.1f}) -> {red_of_E(E_VOL):5.1f}%")
    if cross: print(f"\nCROSSOVER: below |dPsi|={-cross:.0f} mV the volume-targeted, "
                    f"potential-INDEPENDENT particle beats MitoQ.")
    print("saved: figure_superoxide_targeting.png")

if __name__ == "__main__":
    main()