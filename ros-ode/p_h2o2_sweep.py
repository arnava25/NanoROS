"""
p_h2o2_sweep.py -- sweep the (unconstrained) transbilayer H2O2 transport coefficient
====================================================================================
Downstream exploratory. Imports ros_core_2c. Saves a figure; no window (Agg backend).
Run:  python3 p_h2o2_sweep.py   (from scripts/twocompartment/, next to ros_core_2c.py)

WHAT THIS DOES
--------------
The mechanistic demo prints the mito/cyto H2O2 split at ONE placeholder P_h2o2. That
single point is meaningless alone. This sweeps P_h2o2 across its plausible range and
shows the whole hypothesis space:
  * fast transport  -> matrix and cytosol equilibrate -> ratio -> 1, bulk -> the 1C value
                       ("1C model is fine; mito targeting buys little")
  * slow transport  -> matrix traps H2O2 -> ratio climbs, matrix load >> 1C value
                       ("1C model hides a loaded matrix; mito targeting matters")
The scientific question is WHERE the plausible range lands relative to that crossover.

CONFIDENCE
----------
(a) established: the SHAPE of the curves and the two limiting regimes (model structure).
(c) speculative -- VERIFY before any conclusion enters a design doc:
    P_h2o2 is the model's EFFECTIVE transport coefficient (permeability x area / volume,
    lumped, ~1/min), NOT a raw membrane permeability (cm/s). Converting to a literature
    H2O2 permeability needs the inner-membrane area / matrix-volume ratio (morphometry).
    So P_RANGE and LIT_BAND below are EDITABLE PLACEHOLDERS, not primary-lit-anchored.
    The x-axis positioning against "the literature" is only as good as that anchor.
"""
import matplotlib
matplotlib.use("Agg")           # no window pop during script run
import numpy as np
import matplotlib.pyplot as plt

import ros_core_2c as M
from ros_core_2c import Params2C, rhs_2c, ic_2c, observable_h2o2, S2, _integrate

# --------------------------------------------------------------------------------
# EDITABLE PLACEHOLDERS -- (c)-speculative. Replace after primary-lit verification.
# --------------------------------------------------------------------------------
P_RANGE  = np.logspace(-3, 2, 80)     # effective transport coeff (1/min), log-spaced
LIT_BAND = (0.1, 10.0)                # PLACEHOLDER "plausible" band; NOT verified
T_END    = 800.0                      # min, integrate to steady state

def steady_2c(P):
    """Mechanistic config (O2- impermeant, matrix-localized generation, mito import-only,
    calibrated overlay). Returns (H2O2_mito, H2O2_cyto, bulk_observable) at steady state."""
    p = Params2C(f_mito_ros=1.0, o2_transport=False, mito_has_synthesis=False).calibrated()
    p.P_h2o2 = P
    _, Y = _integrate(rhs_2c, ic_2c(p), p, T_END)
    m, c = Y[S2['H2O2m']][-1], Y[S2['H2O2c']][-1]
    return m, c, observable_h2o2(Y[:, -1], p, "volume")

def h2o2_1c_reference():
    """1C steady-state H2O2 under the calibrated overlay -- the flat reference line."""
    if M.rp is None:
        return None
    q = dict(M.rp.BASE); q['k_cl'] = 0.0267; q['k_gsh'] = 3.0
    _, Y = M.rp.acute_run(q, t_end=T_END, n_eval=2000)
    return float(Y[2, -1])

def main():
    mito = np.empty_like(P_RANGE); cyto = np.empty_like(P_RANGE); bulk = np.empty_like(P_RANGE)
    for i, P in enumerate(P_RANGE):
        mito[i], cyto[i], bulk[i] = steady_2c(P)
    ratio = mito / cyto
    ref = h2o2_1c_reference()
    # matrix relative to the bulk-inferred (1C) value -- bounded, no denominator artifact
    enrich = mito / ref if ref else ratio

    # crossover: P where matrix enrichment drops through 2x the 1C value
    cross = None
    for i in range(1, len(P_RANGE)):
        if (enrich[i-1] - 2.0) * (enrich[i] - 2.0) < 0:
            lo, hi = np.log10(P_RANGE[i-1]), np.log10(P_RANGE[i])
            f = (2.0 - enrich[i-1]) / (enrich[i] - enrich[i-1])
            cross = 10 ** (lo + f * (hi - lo)); break

    # ---- figure ----
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))

    ax1.plot(P_RANGE, mito, color="#b2182b", lw=2, label="H2O2 matrix")
    ax1.plot(P_RANGE, cyto, color="#2166ac", lw=2, ls="--", label="H2O2 cytosol")
    ax1.plot(P_RANGE, bulk, color="#4d4d4d", lw=1.3, ls=":", label="H2O2 bulk (vol-weighted)")
    if ref is not None:
        ax1.axhline(ref, color="black", lw=1, ls=(0, (1, 1)),
                    label=f"1C calibrated ({ref:.2f} nM)")
    ax1.axvspan(*LIT_BAND, color="0.85", alpha=0.6, label="plausible band (PLACEHOLDER)")
    ax1.set_xscale("log"); ax1.set_xlabel("P_h2o2  (effective transport coeff, 1/min)")
    ax1.set_ylabel("H2O2 (nM)"); ax1.set_title("A. Compartment H2O2 vs transport")
    ax1.legend(fontsize=8, loc="upper right"); ax1.grid(alpha=0.3)

    ax2.plot(P_RANGE, enrich, color="#762a83", lw=2)
    ax2.axhline(1.0, color="black", lw=1, ls=(0, (1, 1)), label="matrix = 1C value (well-mixed)")
    ax2.axhline(2.0, color="0.5", lw=0.8, ls="--", label="2x the 1C value")
    ax2.axvspan(*LIT_BAND, color="0.85", alpha=0.6, label="plausible band (PLACEHOLDER)")
    if cross is not None:
        ax2.axvline(cross, color="#762a83", lw=0.8, ls=":")
        ax2.annotate(f"2x at\nP={cross:.2g}", (cross, 2.0), fontsize=8,
                     xytext=(8, 8), textcoords="offset points", color="#762a83")
    ax2.set_xscale("log")
    ax2.set_xlabel("P_h2o2  (effective transport coeff, 1/min)")
    ax2.set_ylabel("matrix H2O2 / 1C value  (enrichment)")
    ax2.set_title("B. Matrix enrichment vs transport")
    ax2.legend(fontsize=8, loc="upper right"); ax2.grid(alpha=0.3)

    fig.suptitle("P_h2o2 sweep: matrix loading vs transport rate "
                 "(x-axis anchor is a PLACEHOLDER -- verify vs primary lit)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = "figure_p_h2o2_sweep.png"
    fig.savefig(out, dpi=150); plt.close(fig)

    # ---- console readout ----
    print("=" * 68)
    print("P_h2o2 sweep  (mechanistic config, calibrated overlay)")
    print("=" * 68)
    print(f"1C calibrated H2O2 reference : {ref:.3f} nM" if ref else "1C reference: N/A (rp not found)")
    print(f"fast transport  (P={P_RANGE[-1]:.3g}) : matrix {mito[-1]:.3f} nM "
          f"= {enrich[-1]:.2f}x 1C, bulk {bulk[-1]:.3f}")
    print(f"slow transport  (P={P_RANGE[0]:.3g}) : matrix {mito[0]:.3f} nM "
          f"= {enrich[0]:.2f}x 1C, cyto {cyto[0]:.3f}")
    if cross: print(f"enrichment=2x crossover       : P_h2o2 = {cross:.3g} 1/min")
    print("-" * 68)
    print("Values across the PLACEHOLDER plausible band:")
    for edge, label in ((LIT_BAND[0], "band low "), (LIT_BAND[1], "band high")):
        m, c, b = steady_2c(edge)
        en = m / ref if ref else m / c
        print(f"  P={edge:6.3g}: matrix {m:6.3f}  cyto {c:6.3f}  "
              f"enrich {en:4.2f}x  bulk {b:5.3f}")
    print("-" * 68)
    print(f"saved: {out}")
    print("Note: the band is a PLACEHOLDER. Which regime you're in is not decided until")
    print("P_h2o2 is anchored to a verified inner-membrane H2O2 permeability + morphometry.")

if __name__ == "__main__":
    main()