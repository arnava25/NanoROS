# NanoROS

Computational design specification for a dual-trigger nanoparticle therapeutic targeting mitochondrial ROS in doxorubicin-induced cardiotoxicity.

## Structure

```
nanoros/
├── models/
│   ├── parameters.py          ← all parameters with citations
│   ├── ode_core.py            ← ODE system + run() interface
│   └── ros_model_extended.py  ← extended model with CL_OOH and 4-HNE
├── analysis/
│   ├── linker_kinetics.py     ← AND gate firing times + specificity table
│   └── batch_variability.py   ← delivery cascade + batch tolerance + MitoQ
├── figures/
│   ├── fig1_validation.py     ← ODE state trajectories
│   ├── fig2_sensitivity.py    ← 2D AND gate contour
│   └── fig3_comparator.py     ← NanoROS vs dex vs MitoQ
├── docs/
│   ├── NanoROS_Research_Specification.docx
│   └── key_outputs.txt
└── run_all.py                 ← runs everything
```

## Quick start

```bash
pip install -r requirements.txt
python3 run_all.py
```

## Key numbers

| Constraint | Value | Target |
|---|---|---|
| [Dox]_ss | 24.9 nM | 25 nM |
| Peak H2O2 | 6.88 nM | 6.89 nM |
| Damage deadline | t=123 min | t=123 min |
| Nrf2 2× | t=160 min | t=160 min |
| AND gate fires | t=22 min | < 123 min |
| k_sc_O2 | 23,515 min⁻¹ | ≥ 5 min⁻¹ |
| NPs/cardiomyocyte | 357 | > 10 |

ODE validated against Ludke et al. 2017 (J Mol Cell Cardiol).
