# MPPT comparative simulation

This repository contains a simulation-only comparison of three maximum power
point tracking (MPPT) methods for a photovoltaic (PV) string:

- Perturb and Observe (P&O)
- Incremental Conductance
- Particle Swarm Optimization (PSO)

The project currently has two plant models. The **ideal-voltage benchmark**
applies an algorithm's requested PV voltage immediately, making it useful for
isolating and comparing MPPT logic. The separate **averaged boost-converter
experiment** adds input/output capacitors, an inductor, a resistive load, and
converter control dynamics. Neither model is hardware validation.

> **Status:** This is a research prototype and reproducible starting point, not
> yet a complete statistical comparison. The current ideal benchmark repeats
> PSO across 30 random seeds for one shading case, but does not provide equally
> broad repeated-run statistics for every method and scenario. See
> [Next steps for a defensible comparison](#next-steps-for-a-defensible-comparison)
> before treating the current results as a final ranking.

## What the simulation models

### PV array

`pv_model.py` implements a single-diode equivalent-circuit model for a
Kyocera KC200GT module. Its parameters reproduce the published standard test
condition (STC: 1000 W/m², 25 °C) nameplate points closely. The simulated array
uses three 54-cell modules in series. Each module is divided into three
18-cell substrings with an idealized 0.7 V bypass-diode clamp. This allows the
module-level and substring-level shading cases to create multiple peaks on the
array's power-voltage curve.

The off-STC temperature correction is semi-empirical. This is not a full
CEC/De Soto model or a manufacturer-certified characterization.

### MPPT algorithms

`mppt_algorithms.py` contains the controller implementations. Each requests a
PV voltage setpoint:

- **P&O** perturbs the voltage and reverses direction when measured power drops.
- **Incremental Conductance** estimates the power-voltage slope and moves the
  voltage according to which side of the MPP the operating point is on.
- **PSO** evaluates candidate voltages, tracks each particle's best result and
  the swarm's best result, and then holds the best candidate. Its search is
  stochastic and depends on its seed and settings.

### Simulation modes

1. **Ideal voltage actuation** (`experiments.py`): applies the requested voltage
   at each sample with no converter transient. This is a controlled algorithm
   benchmark, not a model of switching hardware.
2. **Averaged boost converter** (`boost_converter.py` and
   `dynamic_converter_experiment.py`): simulates converter state dynamics with
   an averaged continuous-conduction model. It does not resolve switching
   ripple, detailed semiconductor losses, or all parasitics.

## Quick start

Run these commands from the directory containing this README. Python 3.10 or
newer is recommended.

### Windows PowerShell

```powershell
py -3 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python validate_project.py
python run_all.py
python dynamic_converter_experiment.py
```

If PowerShell does not allow environment activation, use the environment's
interpreter directly:

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe validate_project.py
.\.venv\Scripts\python.exe run_all.py
.\.venv\Scripts\python.exe dynamic_converter_experiment.py
```

### macOS or Linux

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python validate_project.py
python run_all.py
python dynamic_converter_experiment.py
```

`run_all.py` first runs `validate_project.py`, then runs the baseline and
robustness experiments in `experiments.py`. The dynamic converter experiment
is separate and must be run explicitly. The scripts write generated files to
`outputs/`; they may overwrite files with the same names from a previous run.

## Experiments currently included

Running `python run_all.py` performs the following:

1. Checks the PV module's modeled STC points against its datasheet values.
2. Plots I-V and P-V curves for uniform irradiance, reduced irradiance,
   whole-module shading, and one shaded substring.
3. Compares all three methods over a 400-sample irradiance step from 1000 to
   600 W/m² and under fixed whole-module partial shading.
4. Compares the methods when one of the nine 18-cell substrings is shaded.
5. Runs 30 PSO seeds for the whole-module shading case.
6. Tests P&O step-size and Incremental Conductance dead-band sensitivity.
7. Sweeps cell temperature, tests current-sensor noise, and compares a gradual
   irradiance ramp.

Running `python dynamic_converter_experiment.py` separately compares the same
three methods in the averaged converter model for an irradiance step and
whole-module partial shading.

The baseline experiment uses 400 controller samples in its main comparisons.
These are **samples, not seconds**. The ideal-actuation baseline does not
represent a physical elapsed time. The dynamic experiment uses a 0.05 s MPPT
sample period, so its sample count can be converted to seconds using that
configured period.

## Generated outputs

The files below are generated under `outputs/` when the corresponding script
runs. Exact output names are defined in the Python scripts.

| Output | Produced by | What it contains |
|---|---|---|
| `module_stc_validation.json` | `experiments.py` | Modeled and datasheet module STC values |
| `pv_characteristics.png`, `pv_mpp_points.csv` | `experiments.py` | PV curves and detected local MPPs |
| `comparison_results.csv` | `experiments.py` | Summary metrics for baseline scenarios and algorithms |
| `irradiance_step_*.csv`, `partial_shading_*.csv`, `partial_substring_shading_*.csv` | `experiments.py` | Per-sample traces for each algorithm and scenario |
| `*_comparison.png` | `experiments.py` | Power and voltage plots for the main comparisons |
| `pso_30_seed_statistics.json` | `experiments.py` | Repeated-run PSO statistics for whole-module shading |
| `algorithm_sensitivity.csv`, `temperature_sweep.csv`, `current_noise_sensitivity.csv`, `irradiance_ramp_results.csv` | `experiments.py` | Sensitivity and robustness results |
| `run_summary.json` | `experiments.py` | Machine-readable summary of the baseline run |
| `averaged_boost_converter_comparison.png`, `averaged_boost_converter_results.json` | `dynamic_converter_experiment.py` | Dynamic converter plots and summary metrics |

The tracked [RESULTS_SUMMARY.md](RESULTS_SUMMARY.md) is a human-readable
snapshot of a previous/current example run. It is **not automatically updated**
when the scripts generate new outputs. If model parameters or experiment
settings change, regenerate the outputs and update the summary based on those
new files.

### Reading the metrics

- **Tracking efficiency** compares accumulated PV power to accumulated
  theoretical MPP power for the same profile. It is an aggregate across the
  complete run, not just the final operating point.
- **Global peak detected** checks whether the mean power over the last 30
  samples is at least 98% of the mean theoretical MPP over those samples. It is
  a threshold-based indicator, not proof of global optimization in every
  operating condition.
- **Settling time** in the baseline CSV is reported in samples after a change.
  A value of `-1` means the defined 2% tracking band was not maintained for a
  10-sample window before the run ended.
- **Ripple** and final-window metrics depend on the particular slices and
  thresholds used in `experiments.py`; document these definitions when
  presenting the results.

## Next steps for a defensible comparison

Use this checklist before presenting a final comparative study:

1. **Agree on the study question and protocol.** Decide which operating
   conditions, shading patterns, array configuration, and simulation mode
   answer the research question. Keep ideal-actuation results separate from
   dynamic-converter results.
2. **Validate the plant and assumptions.** Check PV parameters, bypass-diode
   assumptions, temperature behavior, and converter/load parameters against
   sources or independent measurements. Record any assumptions that cannot be
   validated. The references and model provenance are in [SOURCES.md](SOURCES.md).
3. **Set a fair tuning procedure.** Choose each algorithm's parameters before
   the final comparison. Tune using a documented set of development scenarios,
   then evaluate on separate held-out scenarios. Do not tune one method against
   the final test cases while leaving the others at defaults.
4. **Use repeated runs consistently.** Run every stochastic method over a
   predeclared set of seeds for every applicable scenario. Keep profiles,
   starting conditions, sample periods, and stopping criteria consistent.
   P&O and Incremental Conductance may be deterministic for fixed inputs, but
   repeated profiles can still test sensitivity to initial conditions/noise.
5. **Predefine metrics and aggregation.** In addition to tracking efficiency,
   consider final tracking error, settling time in seconds, steady-state
   oscillation, global-peak success rate, and computational cost. State how
   failures, incomplete settling, and repeated seeds are treated. Report
   distributions (for example mean and standard deviation or median and
   interquartile range), not only a best run.
6. **Make the implementation record reproducible.** Save the exact Python and
   dependency versions, algorithm parameters, random seeds, scenario
   definitions, simulation time steps, and source revision alongside results.
   Keep raw per-run traces as well as aggregate tables.
7. **Check numerical and measurement effects.** Repeat relevant cases with
   smaller simulation steps or finer PV curves, and test realistic sensor noise
   and converter response assumptions. Explain whether differences between
   methods are meaningful relative to numerical resolution and run-to-run
   variation.
8. **Regenerate and review outputs.** Run validation, run both experiment
   modes, inspect the CSVs and plots, then update `RESULTS_SUMMARY.md` and the
   study report from the generated data. Do not describe simulation results as
   hardware performance.

These are study recommendations, not additional behavior currently
implemented by the scripts. In particular, consistent multi-seed runs,
automated confidence intervals, metadata capture, and held-out tuning/test
scenarios would need to be added before claiming those properties.

## Project files

| File | Purpose |
|---|---|
| `pv_model.py` | PV module, series-string and substring-shading models |
| `mppt_algorithms.py` | P&O, Incremental Conductance, and PSO controllers |
| `experiments.py` | Ideal-actuation benchmark and robustness experiments |
| `boost_converter.py` | Averaged CCM boost-converter plant and voltage controller |
| `dynamic_converter_experiment.py` | MPPT experiments coupled to the converter model |
| `validate_project.py` | STC, string-curve, shading, and PSO sanity checks |
| `run_all.py` | Runs validation followed by the baseline experiment suite |
| `requirements.txt` | Python package dependencies |
| `SOURCES.md` | Model provenance and references |
| `RESULTS_SUMMARY.md` | Human-readable summary of a particular run |
| `outputs/` | Generated plots, traces, tables, and JSON summaries |

## Limitations

No physical hardware is modeled or tested. Ideal actuation assumes the PV
terminal voltage follows the requested setpoint immediately. The averaged
converter does not resolve switching ripple, detailed semiconductor losses, or
all parasitics. The PV array uses a fixed idealized bypass-diode drop and does
not model reverse breakdown, detailed cell-to-cell mismatch, spectral/angular
effects, or thermal dynamics. Results are therefore specific to these models,
scenarios, and settings; they should be treated as simulation results, not
universal algorithm rankings or hardware performance claims.
