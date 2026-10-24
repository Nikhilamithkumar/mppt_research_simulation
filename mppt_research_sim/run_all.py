"""Run validation and all MPPT simulation experiments in sequence.

This script is the top-level entry point for the project. It runs the quick
sanity checks first and then executes the full baseline experiment set.
"""
from __future__ import annotations
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
# Keep the project entry point simple by running the model checks before the
# larger MPPT experiment suite. The optional dynamic converter experiment is
# intentionally launched separately because it is a second, more detailed model.
SCRIPTS = ["validate_project.py", "experiments.py"]

for script in SCRIPTS:
    print(f"\n=== Running {script} ===", flush=True)
    subprocess.run([sys.executable, str(ROOT / script)], cwd=ROOT, check=True)

print("\nBaseline model checks and experiments completed. To run the optional averaged converter test separately, use: python dynamic_converter_experiment.py")
