#!/usr/bin/env python3
"""
NanoROS — Run All
==================
Run from repo root: python3 run_all.py
Produces all figures and prints all key outputs.
"""
import subprocess, sys, os

def run(path, label):
    print(f"\n{'='*60}\n{label}\n{'='*60}")
    r = subprocess.run([sys.executable, path], capture_output=False)
    if r.returncode != 0:
        print(f"ERROR in {label}")

if __name__ == '__main__':
    base = os.path.dirname(os.path.abspath(__file__))

    # Key outputs
    sys.path.insert(0, os.path.join(base, 'models'))
    from ode_core import key_outputs
    key_outputs()

    # Analysis
    run(os.path.join(base, 'analysis', 'linker_kinetics.py'),   'AND Gate Linker Kinetics')
    run(os.path.join(base, 'analysis', 'batch_variability.py'), 'Batch Variability + Delivery Cascade')

    # Figures
    os.chdir(os.path.join(base, 'figures'))
    run(os.path.join(base, 'figures', 'fig1_validation.py'),    'Figure 1: ODE Validation')
    run(os.path.join(base, 'figures', 'fig2_sensitivity.py'),   'Figure 2: AND Gate Contour')
    run(os.path.join(base, 'figures', 'fig3_comparator.py'),    'Figure 3: Comparator')

    print(f"\n{'='*60}\nAll outputs complete. Figures in figures/\n{'='*60}")
