# PLOS ONE revision (PONE-D-26-28903)

Calibration, global sensitivity, and uncertainty-quantification scripts added for the revision.
- parameters_calibrated.py — k_cl, k_gsh calibrated to Ludke et al. 2017
- ros_core_paper.py — model core used for revision analyses
- sensitivity_uncertainty.py — Sobol + PRCC global sensitivity and full-parameter UQ
- calibration_report.txt, revision_outputs.txt — generated outputs

Run:  pip install -r ../requirements.txt  &&  python sensitivity_uncertainty.py
