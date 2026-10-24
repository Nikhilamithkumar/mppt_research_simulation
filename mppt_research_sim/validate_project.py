"""Lightweight reproducibility/physics sanity checks. Run: python validate_project.py

This script performs the quick checks that confirm the PV model is behaving as
expected before the larger experiment suite is run. It verifies datasheet-level
STC agreement and confirms that the shading model produces the expected number of
power peaks.
"""
import math
from pv_model import PVModule, make_string_iv_curve, make_substring_string_iv_curve, validate_module_stc
from experiments import run_profile, step_irradiance, step_temp, shaded_irradiance, metrics, MAX_V, N_MODULES, MODULE

# Confirm the STC module model is close to the published KC200GT performance
# numbers before using it as a benchmark for MPPT studies.
v = validate_module_stc()
assert abs(v["model_isc_a"] - v["datasheet_isc_a"]) / v["datasheet_isc_a"] < 0.01
assert abs(v["model_voc_v"] - v["datasheet_voc_v"]) / v["datasheet_voc_v"] < 0.01
assert abs(v["model_pmp_w"] - v["datasheet_pmp_w"]) / v["datasheet_pmp_w"] < 0.01

# The three-module series string should behave like approximately 600 W under
# uniform 1000 W/m^2 illumination.
uniform = make_string_iv_curve((1000,1000,1000), 25, points=6000)
assert abs(uniform.mpp()["power_w"] - 3 * 200.0) / (3*200.0) < 0.01
# Partial shading should create multiple maxima rather than a single smooth P-V
# curve. This is a key requirement for the MPPT comparisons in the project.
shaded = make_string_iv_curve((1000,1000,300), 25, points=6000)
peaks = shaded.mpps()
assert len(peaks) >= 2, f"Expected multiple peaks, got {peaks}"
assert peaks[0]["power_w"] > peaks[1]["power_w"]
within = make_substring_string_iv_curve([1000]*8 + [300], 25, points=6000)
assert len(within.mpps()) >= 2
# Ensure the PSO algorithm can navigate the multi-peak landscape and identify a
# global maximum rather than getting stuck on a local peak.
pso = run_profile("partial_shading", shaded_irradiance, step_temp, "PSO", n_samples=250,
                  start_v=0.90*N_MODULES*MODULE.v_oc_v, seed=4)
assert metrics(pso, step_index=250)["global_peak_detected"]
print("PASS: STC validation, series-string MPP, multi-peak shading curves, and PSO global-MPP test")
print(f"KC200GT STC model Pmp: {v['model_pmp_w']:.3f} W (nameplate {v['datasheet_pmp_w']:.1f} W)")
print(f"Three-module uniform MPP: {uniform.mpp()['power_w']:.2f} W")
print("Whole-module-shading peaks:", [(round(p['voltage_v'],2),round(p['power_w'],2)) for p in peaks])
