# NanoROS

**A computational design specification for a dual-trigger nanoparticle therapeutic targeting mitochondrial ROS in doxorubicin-induced cardiotoxicity.**

Arnav Amit | Independent Computational Researcher | June 2026
github.com/arnava25/NanoROS

---

## What This Is

NanoROS is a fully specified, computationally parameterized nanoparticle design for precision mitochondrial ROS scavenging in doxorubicin (dox) cardiotoxicity. It is not a built particle — it is a design specification precise enough to hand to a wet lab, with every parameter traceable to a validated mechanistic model of the disease it aims to treat.

The core claim: existing mitochondria-targeted antioxidants (MitoQ) fail not because of poor mitochondrial uptake in diseased cells, but because they have no viable therapeutic window at any membrane potential. At therapeutic plasma concentrations, MitoQ accumulates to intramitochondrial levels 27-fold above its pro-oxidant ceiling even in healthy mitochondria. NanoROS addresses this by replacing constitutive accumulation with a dual-trigger AND gate that only activates in the specific biochemical environment of the dox-damaged cardiomyocyte mitochondrion.

The full design rationale, five falsifiable novelty claims, quantitative parameters, experimental roadmap, and honest failure mode landscape are in `docs/NanoROS_Research_Specification.pdf`.

---

## The Engine: ROS-Cardio ODE Model

The nanoparticle design is not empirically tuned. Every parameter traces back to a validated kinetic model of dox-induced ROS accumulation in cardiomyocytes — the ROS-Cardio ODE model.

The model tracks eight state variables in a dox-exposed cardiomyocyte:

- Doxorubicin concentration [Dox]
- Superoxide [O2-]
- Hydrogen peroxide [H2O2]
- Keap1 oxidation state
- Nrf2 nuclear concentration
- Glutathione [GSH]
- Cumulative irreversible damage D
- Extended module: CL-OOH and free 4-HNE

The Nrf2-GSH feedback axis captures the endogenous antioxidant response. The extended module adds cardiolipin peroxidation kinetics and 4-HNE production and clearance, which are required to parameterize AND gate Arm B.

**Validation:** Model H2O2 trajectories are validated against ROS time-course measurements in adult rat cardiomyocytes exposed to doxorubicin (Ludke et al., PLoS ONE 2017). The model accurately recapitulates early-phase oxidative stress dynamics over a 6-hour exposure window.

**Key outputs that constrain the nanoparticle design:**

| Constraint | Value | Design Use |
|---|---|---|
| [Dox]_ss | 24.9 nM | Target cell context |
| Peak H2O2 | 6.88 nM at t = 587 min | Trigger range upper bound |
| Peak O2- | 2.51 nM at t = 547 min | SOD mimic sizing |
| Damage deadline (D = 50 AU) | t = 123 min | NP must act before this |
| Nrf2 2x activation | t = 160 min | Endogenous defence too slow |
| GSH depletion | 0.8% | NOT a viable trigger signal |
| [4-HNE]_ss (estimated) | 0.034-0.46 uM | Arm B trigger range |
| AND gate fires | t = 22 min | Well within damage deadline |
| k_sc_O2 achieved | 23,515 min-1 | 4,700x above k_sod floor |

---

## The Design: Four-Layer Architecture

**Layer 1 — Systemic delivery and cardiac accumulation**
CRPPR peptide cardiac homing + anti-ICAM-1 targeting. ICAM-1 is upregulated 3-5-fold on cardiac endothelium by dox stress — the NP accumulates more efficiently as disease severity increases. CD47 mimetic extends circulation half-life. pH-labile PEG corona (hydrazone linker) sheds at endosomal pH 5.5.

**Layer 2 — Cardiomyocyte selectivity and endosomal escape**
DOPE/CHEMS fusogenic lipid coat (7:3 molar ratio) triggers membrane fusion at endosomal pH 5.5, achieving 5-15% endosomal escape efficiency. AT1R targeting peptide provides cardiomyocyte vs. immune cell discrimination (Phase 2; omitted from Minimum Viable Design). SS-31 (elamipretide sequence) binds cardiolipin on the outer mitochondrial membrane — membrane-potential-independent, validated by FDA accelerated approval of elamipretide (Forzinity, September 2025) for Barth syndrome.

**Layer 3 — Dual-trigger AND gate release**
Two mechanistically orthogonal polymer arms must co-activate for shell disassembly:

- Arm A: Dioxaborolane boronate ester, H2O2 temporal integrator, k2 = 600 M-1s-1 at matrix pH 8.0, GSH-inert, fires at t = 22 min
- Arm B: O-alkylhydroxylamine oxime, 4-HNE lipid peroxidation marker, k2 = 300-1,000 M-1s-1 depending on branch, GSH-inert, acid-stable

Arm B is a three-branch contingent design gated by Experiment 0 (LC-MS/MS measurement of free [4-HNE] in dox-treated cardiomyocyte mitochondrial fractions). The branch determines linker chemistry; all three branches fire within the 123-minute damage deadline.

**Layer 4 — Payload: MnTE-2-PyP**
Manganese(III) meso-tetrakis(N-ethylpyridinium-2-yl)porphyrin. kcat/Km for O2- = 7x10^8 M-1s-1. No pro-oxidant ceiling. No Fenton risk. No membrane potential dependence. Crosses the outer mitochondrial membrane via VDAC passively. Phase II clinical safety data from Duke group. Delivers ceiling-level protection (~87% damage reduction) when released before the damage deadline.

---

## Repository Structure

```
NanoROS/
├── README.md
├── requirements.txt
├── run_all.py                        <- runs full pipeline
├── docs/
│   ├── NanoROS_Research_Specification.pdf   <- full design document
│   └── key_outputs.txt
├── model/
│   ├── parameters.py                 <- all parameters with citations
│   ├── ode_core.py                   <- ODE system + run() interface
│   ├── ros_model_extended.py         <- extended model: CL-OOH and 4-HNE
│   └── legacy/                       <- original exploratory scripts
│       ├── ros_model.py              <- baseline model, interventions, sensitivity
│       ├── cumulative_dosing.py      <- 6-cycle AC regimen simulation
│       ├── patient_variability.py    <- 200-patient population simulation
│       ├── dex_phase.py              <- dexrazoxane vs sulforaphane comparison
│       └── validation_plot.py        <- validation against Ludke et al. 2017
├── analysis/
│   ├── linker_kinetics.py            <- AND gate firing times + specificity table
│   └── batch_variability.py          <- delivery cascade + batch tolerance
└── figures/
    ├── fig1_validation.png
    ├── fig2_sensitivity_contour.png
    └── fig3_comparator.png
```

### Legacy Scripts

The `model/legacy/` directory contains the original ROS-Cardio ODE scripts written before the modular NanoROS architecture was built. These are the engine the design grew from:

- `ros_model.py` — the first complete implementation: baseline trajectory, MnTE-2-PyP and sulforaphane intervention comparison, dose-response curves, sensitivity analysis, and the two-component damage variable. This is where the key design constraints (damage deadline t = 123 min, Nrf2 lag t = 160 min, GSH depletion 0.8%) were first derived.
- `cumulative_dosing.py` — extends the model across a clinical 6-cycle AC chemotherapy regimen. Showed that sulforaphane (Nrf2 activator) reduces irreversible damage by 32.7% across cycles, dexrazoxane by 52.6%, and combination therapy by 71.0%. These results established that upstream O2- interception (MnTE-2-PyP mechanism) outperforms downstream iron chelation (dexrazoxane mechanism) for the mitochondrial ROS component of cardiotoxicity.
- `patient_variability.py` — 200-patient population simulation with lognormal parameter variation. Confirmed consistent protection across patient phenotypes (IQR 31.6-33.7%), establishing that the design is robust to biological variability in the target parameter space.
- `dex_phase.py` — mechanistic comparison of dexrazoxane (source reduction via iron chelation) vs. sulforaphane (Nrf2-mediated GSH upregulation) with phase portrait analysis. Informed the decision to use MnTE-2-PyP (direct O2- scavenging) rather than Nrf2 activation as the primary payload mechanism.
- `validation_plot.py` — fits model H2O2 trajectory to Ludke et al. 2017 experimental data. The calibrated parameters from this script propagate into all downstream design decisions.

The modular architecture in `model/` refactored these scripts into a clean ODE core with separated parameters, extended the model with CL-OOH and 4-HNE modules, and added the AND gate analysis. The legacy scripts are preserved as the documented reasoning chain behind the design.

---

## Quick Start

```bash
pip install -r requirements.txt
python3 run_all.py
```

---

## Five Novelty Claims

1. **ODE-model-constrained trigger specification** — every nanoparticle parameter traces to a validated disease model output, not chemical intuition
2. **Orthogonal dual-signal AND gate** — co-requirement of H2O2 and 4-HNE achieves specificity no single-trigger design can match; zero false positives across eight modelled biological scenarios
3. **Cardiolipin coherence** — single lipid governs mitochondrial docking, Arm B trigger generation, and NP detachment through its oxidation state
4. **ICAM-1 / 4-HNE self-amplifying therapeutic index** — both accumulation and release become more efficient as disease severity increases
5. **4-HNE oxime trigger as novel NP release chemistry** — first use of 4-HNE as a nanoparticle trigger; resolves membrane-contact, GSH competition, and acid instability problems that defeat direct CL-OOH sensing

---

## Experimental Roadmap

Five pre-specified gates determine which of three design branches is implemented. No experimental result renders the program undefined.

| Exp. | Name | Go Criterion | If Fail |
|---|---|---|---|
| 0 * | Free [4-HNE] LC-MS/MS in dox-treated CM mito fraction | [4-HNE] detectable >= 10 nM | Branch C: switch to MDA sensing |
| 1A * | MnTE-2-PyP encapsulation feasibility | EE > 5% AND boronate retention > 85% | Alternative: separate NP populations |
| 1B | FRET co-assembly homogeneity | FRET efficiency > 5% | Compatibiliser C or sequential addition |
| 2 | pH-responsive PEG shedding kinetics | > 80% PEG shed at pH 5.5 in 2 hr | Switch to acetal linker |
| 3 | DOPE/CHEMS endosomal escape | > 5% calcein delivery vs macrophage | Add KALA peptide |

*Experiments 0 and 1A are pre-synthesis gates. All polymer synthesis is conditional on their outcomes.*

Full experimental roadmap (Experiments 0-8) in `docs/NanoROS_Research_Specification.pdf`.

---

## Reference

Ludke A, Akolkar G, Ayyappan P, Sharma AK, Singal PK. Time course of changes in oxidative stress and stress-induced proteins in cardiomyocytes exposed to doxorubicin and prevention by vitamin C. PLoS ONE. 2017;12(7):e0179452.