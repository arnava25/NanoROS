# NanoROS

A computational design specification for a dual-trigger nanoparticle therapeutic targeting mitochondrial ROS in doxorubicin-induced cardiotoxicity.

**Arnav Amit | Independent Computational Researcher**
Specification v1, June 2026 · Status updated September 2026
[github.com/arnava25/NanoROS](https://github.com/arnava25/NanoROS)

---

## What This Is

NanoROS is a computationally parameterized nanoparticle design for precision mitochondrial ROS scavenging in doxorubicin (dox) cardiotoxicity. It is not a built particle. It is a design specification precise enough to hand to a wet lab, with every parameter traceable to a mechanistic model of the disease it aims to treat.

The core claim: existing mitochondria-targeted antioxidants (MitoQ) fail not because of poor mitochondrial uptake in diseased cells, but because they have no viable therapeutic window at any membrane potential. At therapeutic plasma concentrations, MitoQ accumulates to intramitochondrial levels ~27-fold above its pro-oxidant ceiling even in healthy mitochondria (Doughan & Dikalov, *Antioxid Redox Signal* 2007;9:1825–36). NanoROS replaces constitutive accumulation with a dual-trigger AND gate that activates only in the biochemical environment of the dox-damaged cardiomyocyte mitochondrion.

Full design rationale, novelty claims, quantitative parameters, experimental roadmap, and failure mode landscape: `docs/NanoROS_Research_Specification.pdf`.

> **Read this first.** The specification PDF is frozen at June 2026. Several of its locked design choices have since been superseded by the two-compartment modeling work and by primary-literature verification. The [Status](#status-september-2026) section below is authoritative where the two disagree.

---

## Status (September 2026)

**The project is experimentally gated, not knowledge-limited.** The bottleneck is Experiment 0, not further computation.

### Manuscript

The underlying kinetic model is written up as a standalone manuscript, **PONE-D-26-28903**, under review at PLOS ONE. Major revision submitted 12 September 2026; decision pending.

The revision added formal calibration rather than language-softening: weighted least-squares calibration of `k_cl` and `k_gsh` to the Ludke 2017 time course, with profile-likelihood identifiability profiling, global Sobol and PRCC sensitivity analysis, and full-parameter uncertainty propagation.

Code state at submission is frozen: **Zenodo v1.1-plos-revision, DOI [10.5281/zenodo.22718327](https://doi.org/10.5281/zenodo.22718327)**, git tag `v1.1-plos-revision`.

### Retired numbers

The following appeared in the June 2026 README and specification and should not be cited:

| Retired claim | Status |
|---|---|
| Dexrazoxane 52.6% protection | Acute endpoint at uncalibrated BASE. Cycle-6 calibrated value ≈ 55.1%; assumption-contingent across a 40–80% ROS-reduction sweep (35–77%). |
| Combination therapy 71.0% | Same provenance. Not carried forward. |
| Sulforaphane IQR 31.6–33.7% | Artifact of two-parameter sampling. Full-parameter UQ: central 32.5%, IQR 22.0–43.6%, 90% CI 11.9–60.8%. The central estimate is robust; the tight interval was not. |
| MnTE-2-PyP "~87% damage reduction" | A one-compartment ceiling estimate assuming complete upstream O2⁻ elimination. See *Payload* below. |
| "Zero false positives across eight modelled scenarios" | Extension to twelve clinical scenarios identified HFpEF as a possible false positive, mitigated by exclusion criteria rather than by gate design. |

### Calibration

Early-phase goodness of fit improved from R² = −0.30 to R² = 0.94 (RMSE 14.0% → 3.0%) after calibration. `k_ros` and `k_in` are structurally non-identifiable from normalized H2O2 data. `k_chronic_gsh` is the second most influential driver of cycle-6 damage (S_T = 0.19) despite being a weakly constrained parameter.

Note on provenance: `parameters_calibrated.py` holds the calibrated overlay (`k_cl` = 0.0267, `k_gsh` = 3.0) and corresponds to Table 1. The uncertainty-propagation script runs on the uncalibrated BASE (`k_cl` = 0.02, `k_gsh` = 2.0), so UQ distributions are centered on the wrong nominal. Re-centering is queued for the next revision.

### Two-compartment model

`ros_core_2c.py` extends the frozen single-compartment model into separate mitochondrial and cytosolic pools. It passes two structural self-tests: mass conservation under transport-only (drift 7.6e-16) and exact reduction to the frozen `ros_core_paper.acute_run(BASE)` (H2O2 error 1.4e-6). It is a representational gain, not a fit improvement.

Three results from it change the design:

**1. H2O2-gradient mitochondrial targeting is dead.** Anchoring transbilayer transport from inner-membrane permeability (~1e-3 cm/s) and cardiac cristae area-to-volume (~5e5 cm⁻¹) gives an effective P_h2o2 ≈ 1e4 min⁻¹ — roughly 1e4 above the well-mixed crossover, equilibrating in ~2 ms. H2O2 carries no matrix/cytosol concentration advantage. The targeting rationale cannot rest on it.

**2. Superoxide is the only genuinely matrix-confined species** (membrane permeability < 1e-7 cm/s). The targeting argument is therefore about payload *access*, not gradient exploitation.

**3. The payload decision is open, not locked.** In the current parameterization ~97% of chronic matrix damage is GSH-import-limited rather than superoxide-driven, so a superoxide scavenger returns ~3% protection even at 400-fold mitochondrial enrichment. Defining φ as the fraction of baseline irreversible matrix damage attributable to the superoxide pathway, the crossover between an SOD mimetic and a matrix-targeted GSH restorer sits at φ ≈ 0.5, insensitive to enrichment level. The model default sits at φ ≈ 0.02. The two payload classes are complementary, not substitutable.

A 24-hour diagnostic ruled out GSH exhaustion as the cause of late divergence from Ludke 2017 (GSH recovers to ~101% by 24 h via Nrf2). The structural gap is a missing feed-forward damage-to-ROS term, not a spatial-averaging error.

### Open design questions

| Question | Status |
|---|---|
| Payload identity | **Open.** MnTE-2-PyP is superseded as a locked choice. Resolution requires φ, which requires Experiment 0. Tadokoro 2020 and the 2024 *Redox Biology* mito-GPx4 work tilt toward a GSH/GPx4-axis payload. |
| SS-31 affinity for *oxidized* cardiolipin | **Unmeasured.** Birk 2013 (K_D 1.87 µM) and Mitchell 2020 (2.0–2.9 µM) measured native CL only. Novelty claim 3 (cardiolipin coherence) rests on inference by analogy to the cytochrome c / oxidized-CL affinity drop. This is the single largest unverified claim in the design. |
| Arm B linker chemistry | **Open.** The June spec's "acid-stable" oxime claim is contradicted by Kalia & Raines 2008 (oxime hydrolysis in aqueous/acidic media). Free 4-HNE is also cleared in seconds by 10–14 mM matrix GSH plus GSTA4/ALDH2 (Siems 1997, hepatocyte), making the kinetic race tight. Leads: built-in aniline catalysis (anthranilic-acid type) or electronically optimized low-pKa hydrazides. |
| Cardiomyocyte vs. immune selectivity | **Unresolved.** AT1R is functionally expressed on antigen-specific CD8 T cells and macrophages, so it does not discriminate as the spec assumes. CRPPR's ~300-fold cardiac homing is mouse in vivo; the human signal is HCAEC in vitro only. |
| Prior art | pB-DOX (boronate-caged ROS-responsive doxorubicin prodrug liposome) is the nearest competitor. NanoROS distinguishes on the AND gate and on cardioprotective rather than prodrug intent. |

### Settled since June

- **Trigger timing is comfortable.** The boronate integrator (Arm A) crosses at ~22 min and H2O2 reaches 90% of plateau at ~87 min, both inside the 123 min damage deadline. A 300-fold sweep of `k_ox_cl` puts 4-HNE half-rise at 15–38 min, so Arm B is not late. The **strict** AND gate is retained; the graded-gate contingency is unnecessary.
- **GSH is a payload, not a trigger.** Depletion is too small to carry signal, and the matrix deficit is largely static baseline rather than dox-driven.
- **MitoQ is out** on the pro-oxidant ceiling at any ΔΨm — not on depolarization-dependent uptake, which was the original and incorrect framing.

---

## The Engine: ROS-Cardio ODE Model

The nanoparticle design is not empirically tuned. Parameters trace back to a kinetic model of dox-induced ROS accumulation in cardiomyocytes.

State variables:

- Doxorubicin [Dox]
- Superoxide [O2⁻]
- Hydrogen peroxide [H2O2]
- Keap1 oxidation state
- Nrf2 nuclear concentration
- Glutathione [GSH]
- Cumulative irreversible damage D (acute + chronic)
- Extended module: CL-OOH and free 4-HNE

The Nrf2–GSH feedback axis captures the endogenous antioxidant response. The extended module adds cardiolipin peroxidation kinetics and 4-HNE production and clearance, required to parameterize AND gate Arm B.

**Validation.** Model H2O2 trajectories are validated against ROS time-course measurements in adult rat cardiomyocytes exposed to doxorubicin (Ludke et al., *PLoS ONE* 2017). The model recapitulates early-phase dynamics well after calibration (R² = 0.94); late-phase divergence remains, attributed to a missing feed-forward damage-to-ROS term.

**Scope limits worth stating plainly:** one dataset, one directly validated variable (H2O2), most state variables unvalidated, and damage reported in arbitrary units. Every design constraint below inherits those limits.

### Constraints that shape the design

| Constraint | Value | Design use |
|---|---|---|
| [Dox]_ss | 24.9 nM | Target cell context |
| Peak H2O2 | 6.88 nM at t = 587 min | Trigger range upper bound |
| Peak O2⁻ | 2.51 nM at t = 547 min | SOD mimic sizing |
| Damage deadline (D = 50 AU) | t = 123 min | NP must act before this |
| Nrf2 2× activation | t = 160 min | Endogenous defence too slow |
| GSH depletion | 0.8% † | Not a viable trigger signal |
| [4-HNE]_ss (estimated) | 0.034–0.46 µM | Arm B trigger range |
| AND gate fires | t = 22 min | Well within damage deadline |
| k_sc_O2 achieved | 23,515 min⁻¹ | 4,700× above k_sod floor |

† Discrepancy flagged: the later GSH-as-trigger ruling cites whole-cell maximum depletion of 0.16%, not 0.8%. The two figures come from different runs and have not been reconciled. The design conclusion (GSH is not a usable trigger) holds either way.

The damage deadline of 123 min is derived from the single-compartment model. If the feed-forward damage term is added, it moves. Current trigger-timing margins are wide enough that this is unlikely to flip a design decision, but the number is load-bearing and inherited.

---

## The Design: Four-Layer Architecture

**Layer 1 — Systemic delivery and cardiac accumulation.**
CRPPR peptide cardiac homing plus anti-ICAM-1 targeting. ICAM-1 is upregulated 3–5-fold on cardiac endothelium under dox stress (Hayashi 2004), so accumulation improves as disease severity increases. CD47 mimetic extends circulation half-life; the magnitude of that benefit in humans is uncertain. pH-labile PEG corona (hydrazone linker) sheds at endosomal pH 5.5.

**Layer 2 — Cardiomyocyte selectivity and endosomal escape.**
DOPE/CHEMS fusogenic lipid coat (7:3 molar) triggers membrane fusion at endosomal pH 5.5, 5–15% escape efficiency. SS-31 (elamipretide sequence) binds cardiolipin on the mitochondrial membrane — membrane-potential-independent, with the mechanism clinically validated by FDA accelerated approval of elamipretide (Forzinity, September 2025) for Barth syndrome. *The AT1R selectivity element is unresolved; see Open design questions.*

**Layer 3 — Dual-trigger AND gate release.**
Two mechanistically orthogonal polymer arms must co-activate for shell disassembly:

- **Arm A** — dioxaborolane boronate ester, H2O2 temporal integrator, k₂ = 600 M⁻¹s⁻¹ at matrix pH 8.0, GSH-inert, fires at t = 22 min. *Verified early and timely.*
- **Arm B** — 4-HNE-sensing linker as a lipid peroxidation marker, k₂ = 300–1,000 M⁻¹s⁻¹ depending on branch, GSH-inert. *Chemistry open; the June spec's O-alkylhydroxylamine oxime is the leading candidate but its claimed acid stability does not survive the primary literature.*

Arm B is a three-branch contingent design gated by Experiment 0 (LC-MS/MS measurement of free [4-HNE] in dox-treated cardiomyocyte mitochondrial fractions). All three branches fire within the 123-minute damage deadline.

**Layer 4 — Payload.**
*Open pending Experiment 0.* The June specification locked MnTE-2-PyP (Mn(III) meso-tetrakis(N-ethylpyridinium-2-yl)porphyrin; k_cat/K_m for O2⁻ = 7e8 M⁻¹s⁻¹; no pro-oxidant ceiling, no Fenton risk, no membrane-potential dependence, VDAC-permeant, Phase II safety data from the Duke group). The two-compartment result above means that choice is now contingent on φ rather than settled. The competing class is a matrix-targeted GSH/GPx4-axis restorer. Experiment 0 selects between them.

---

## Novelty Claims

1. **ODE-model-constrained trigger specification** — every nanoparticle parameter traces to a validated disease model output rather than chemical intuition. *Strongest claim; survives all revisions.*
2. **Orthogonal dual-signal AND gate** — co-requirement of H2O2 and 4-HNE achieves specificity no single-trigger design can match. *Holds, with HFpEF identified as a possible false positive handled by exclusion criteria.*
3. **Cardiolipin coherence** — a single lipid governs mitochondrial docking, Arm B trigger generation, and NP detachment via its oxidation state. *Conditional on the unmeasured SS-31 / oxidized-CL affinity drop. This is the claim most likely to fail.*
4. **ICAM-1 / 4-HNE self-amplifying therapeutic index** — both accumulation and release become more efficient as disease severity increases. *Holds.*
5. **4-HNE as a nanoparticle release trigger** — first such use. *Holds as a concept; the specific linker chemistry is open.*

---

## Experimental Roadmap

Five pre-specified gates determine which of three design branches is implemented. No experimental result renders the program undefined.

| Exp. | Name | Go criterion | If fail |
|---|---|---|---|
| **0** * | Free [4-HNE] LC-MS/MS in dox-treated CM mito fraction | [4-HNE] detectable ≥ 10 nM | Branch C: switch to MDA sensing |
| **1A** * | Payload encapsulation feasibility | EE > 5% AND boronate retention > 85% | Alternative: separate NP populations |
| 1B | FRET co-assembly homogeneity | FRET efficiency > 5% | Compatibiliser C or sequential addition |
| 2 | pH-responsive PEG shedding kinetics | > 80% PEG shed at pH 5.5 in 2 h | Switch to acetal linker |
| 3 | DOPE/CHEMS endosomal escape | > 5% calcein delivery vs macrophage | Add KALA peptide |

\* Pre-synthesis gates. All polymer synthesis is conditional on their outcomes.

**Experiment 0 is doubly decisive.** It selects the Arm B sensing chemistry *and* resolves φ, which selects the payload class. It is the single result that would most change the design.

**Known feasibility risk on Experiment 0:** free 4-HNE may be measurable only as adduct load rather than as a free-pool concentration, given clearance in seconds by matrix GSH. If so, the free-threshold premise needs rethinking before the branch structure is meaningful. This question should be answered before any sample is run.

Full roadmap (Experiments 0–8) in `docs/NanoROS_Research_Specification.pdf`.

---

## Repository Structure

```
NanoROS/
├── README.md
├── requirements.txt
├── run_all.py                              <- runs full pipeline
├── docs/
│   ├── NanoROS_Research_Specification.pdf  <- full design document (June 2026, see Status)
│   └── key_outputs.txt
├── model/
│   ├── parameters.py                       <- all parameters with citations
│   ├── parameters_calibrated.py            <- calibrated overlay; corresponds to Table 1
│   ├── ode_core.py                         <- ODE system + run() interface
│   ├── ros_core_paper.py                   <- frozen single-compartment module (manuscript)
│   ├── ros_core_2c.py                      <- two-compartment mito/cyto scaffold
│   ├── ros_model_extended.py               <- extended model: CL-OOH and 4-HNE
│   └── legacy/                             <- original exploratory scripts
│       ├── ros_model.py
│       ├── cumulative_dosing.py
│       ├── patient_variability.py
│       ├── dex_phase.py
│       └── validation_plot.py
├── analysis/
│   ├── linker_kinetics.py                  <- AND gate firing times + specificity table
│   ├── batch_variability.py                <- delivery cascade + batch tolerance
│   ├── sensitivity_uncertainty.py          <- Sobol / PRCC / full-parameter UQ
│   └── calibration.py                      <- WLS calibration + profile likelihood
└── figures/
    ├── fig1_validation.png
    ├── fig2_sensitivity_contour.png
    └── fig3_comparator.png
```

### Legacy scripts

`model/legacy/` contains the original ROS-Cardio ODE scripts written before the modular architecture. They are preserved as the documented reasoning chain behind the design, **not** as a current source of numbers — several headline values they produce are listed in the retired table above.

- **`ros_model.py`** — first complete implementation: baseline trajectory, intervention comparison, dose-response, sensitivity, two-component damage. Source of the damage deadline (123 min), Nrf2 lag (160 min), and GSH depletion constraints.
- **`cumulative_dosing.py`** — six-cycle AC regimen. Produced the original sulforaphane / dexrazoxane / combination comparison; the dexrazoxane and combination figures from this script are retired.
- **`patient_variability.py`** — 200-patient lognormal population simulation. Its narrow IQR is a two-parameter sampling artifact; superseded by full-parameter UQ.
- **`dex_phase.py`** — dexrazoxane (iron chelation) vs. sulforaphane (Nrf2/GSH) with phase portraits. Informed the original MnTE-2-PyP payload choice, which the two-compartment work has since reopened.
- **`validation_plot.py`** — fits model H2O2 trajectory to Ludke 2017.

---

## Quick Start

```bash
pip install -r requirements.txt
python3 run_all.py
```

---

## References

Ludke A, Akolkar G, Ayyappan P, Sharma AK, Singal PK. Time course of changes in oxidative stress and stress-induced proteins in cardiomyocytes exposed to doxorubicin and prevention by vitamin C. *PLoS ONE.* 2017;12(7):e0179452.

Doughan AK, Dikalov SI. Mitochondrial redox cycling of mitoquinone leads to superoxide production and cellular apoptosis. *Antioxid Redox Signal.* 2007;9(11):1825–36.

Hiemstra S, et al. Dynamic modeling of Nrf2 pathway activation in liver cells after toxicant exposure. *Sci Rep.* 2022. doi:10.1038/s41598-022-10857-x

Full reference list in `docs/NanoROS_Research_Specification.pdf`.