"""Run baseline MPPT comparisons and additional simulation-only experiments."""
from __future__ import annotations
import csv, json, math, sys, time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pv_model import PVModule, IVCurve, make_string_iv_curve, make_substring_string_iv_curve, validate_module_stc
from mppt_algorithms import make_algorithm

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)
MODULE = PVModule()
N_MODULES = 3
MAX_V = 1.05 * N_MODULES * (MODULE.v_oc_v + MODULE.beta_voc_v_per_k * (0.0 - 25.0))
GLOBAL_CURVE_CACHE: dict[tuple, IVCurve] = {}
GLOBAL_MPP_CACHE: dict[tuple, float] = {}


def curve_for(gs, temp, irradiance_mode="module"):
    # Build or reuse the static I-V curve for a given irradiance profile. The
    # same curve can be reused many times across the same experiment because MPPT
    # controllers repeatedly request the same plant behavior at each sample step.
    key = _curve_key(gs, temp, irradiance_mode)
    if key not in GLOBAL_CURVE_CACHE:
        if irradiance_mode == "module":
            curve = make_string_iv_curve(gs, temp, MODULE, bypass_v=0.7, points=8000)
        elif irradiance_mode == "substring":
            curve = make_substring_string_iv_curve(gs, temp, MODULE, bypass_v=0.7, points=8000)
        else:
            raise ValueError("irradiance_mode must be 'module' or 'substring'.")
        GLOBAL_CURVE_CACHE[key] = curve
        GLOBAL_MPP_CACHE[key] = curve.mpp()["power_w"]
    return GLOBAL_CURVE_CACHE[key]


def _curve_key(gs, temp, irradiance_mode="module"):
    return (irradiance_mode, tuple(round(float(g), 4) for g in gs), round(float(temp), 3))


def run_profile(name: str, irradiance_at, temperature_at, algorithm_name: str,
                n_samples: int = 400, start_v: float | None = None,
                seed: int = 0, po_step: float = 0.7,
                inc_deadband: float = 0.005, pso_particles: int = 6,
                pso_iterations: int = 12, current_noise_std_a: float = 0.0,
                change_indices: tuple[int, ...] = (),
                irradiance_mode: str = "module") -> dict:
    """Run an MPPT algorithm on a quasi-static PV string I-V map.

    At each controller sample the commanded PV voltage is applied immediately
    by an ideal voltage-controlled converter model. This intentionally isolates
    algorithm logic; dynamic converter validation is a separate model.
    """
    # Each profile emulates a sample-by-sample controller loop: the irradiance and
    # temperature set the operating curve, the measured current is computed with a
    # small amount of optional sensor noise, and the MPPT algorithm chooses the next
    # target voltage for the next iteration.
    cache: dict[tuple, IVCurve] = {}
    cache_mpp: dict[tuple, float] = {}
    refs_p = np.zeros(n_samples)
    true_p = np.zeros(n_samples)
    measured_p = np.zeros(n_samples)
    v_hist = np.zeros(n_samples)
    i_hist = np.zeros(n_samples)
    true_i_hist = np.zeros(n_samples)
    irr_hist = np.zeros((n_samples, len(tuple(float(x) for x in irradiance_at(0)))))
    temp_hist = np.zeros(n_samples)
    # Calculate an initial curve to pick a deterministic default start voltage.
    gs0 = tuple(float(x) for x in irradiance_at(0))
    t0 = float(temperature_at(0))
    c0 = curve_for(gs0, t0, irradiance_mode)
    key0 = _curve_key(gs0, t0, irradiance_mode)
    cache[key0] = c0
    cache_mpp[key0] = GLOBAL_MPP_CACHE[key0]
    if start_v is None:
        start_v = 0.58 * c0.voc_v
    alg = make_algorithm(algorithm_name, 0.0, MAX_V, start_v,
                         po_step=po_step, inc_deadband=inc_deadband,
                         pso_particles=pso_particles, pso_iterations=pso_iterations,
                         seed=seed)
    v = alg.initial_voltage()
    rng = np.random.default_rng(seed + 45000)
    for k in range(n_samples):
        # Evaluate one MPPT iteration: determine the irradiance/temperature state,
        # compute the exact and measured PV current at the current voltage, and then
        # ask the selected optimization algorithm for the next operating point.
        gs = tuple(float(x) for x in irradiance_at(k))
        temp = float(temperature_at(k))
        key = _curve_key(gs, temp, irradiance_mode)
        if key not in cache:
            cache[key] = curve_for(gs, temp, irradiance_mode)
            cache_mpp[key] = GLOBAL_MPP_CACHE[key]
        curve = cache[key]
        true_i = max(0.0, float(curve.current(v)))
        noise = float(rng.normal(0.0, current_noise_std_a)) if current_noise_std_a > 0 else 0.0
        meas_i = max(0.0, true_i + noise)
        p_actual = v * true_i
        p_measure = v * meas_i
        ref = cache_mpp[key]
        v_hist[k], i_hist[k], true_i_hist[k] = v, meas_i, true_i
        true_p[k], measured_p[k], refs_p[k] = p_actual, p_measure, ref
        irr_hist[k, :] = gs
        temp_hist[k] = temp
        v = float(alg.next_voltage(v, meas_i, p_measure))
    return {"name": name, "algorithm": algorithm_name, "voltage_v": v_hist,
            "current_a": i_hist, "true_current_a": true_i_hist, "power_w": true_p,
            "measured_power_w": measured_p, "mpp_power_w": refs_p,
            "irradiance_w_m2": irr_hist, "temp_c": temp_hist,
            "efficiency_pct": 100.0 * np.sum(true_p) / max(np.sum(refs_p), 1e-12),
            "seed": seed, "noise_std_a": current_noise_std_a,
            "change_indices": tuple(change_indices)}


def step_irradiance(k):
    g = 1000.0 if k < 200 else 600.0
    return (g, g, g)

def step_temp(k): return 25.0

def shaded_irradiance(k): return (1000.0, 1000.0, 300.0)

def ramp_irradiance(k, n=400):
    # A sampled staircase approximation of a linear ramp. Quantizing to 25 W/m²
    # keeps the repeated physical curve calculations efficient and reproducible.
    g = 1000.0 - 700.0 * (k / max(1, n - 1))
    g = float(np.clip(25.0 * round(g / 25.0), 300.0, 1000.0))
    return (g, g, g)

def fixed_temp(t):
    return lambda _k: float(t)

def metrics(run: dict, step_index: int = 200) -> dict[str, float | int | str]:
    # Summarize the MPPT run by comparing the measured PV power to the true MPP
    # for the current irradiance/temperature profile. The metrics are used both
    # for automation and for the JSON/CSV output tables in the outputs folder.
    p = run["power_w"]
    pref = run["mpp_power_w"]
    n = len(p)
    ripple_slice = p[160:200] if (step_index == 200 and n >= 200) else p[-40:]
    ripple = float(np.max(ripple_slice) - np.min(ripple_slice)) if len(ripple_slice) else float("nan")
    settle = -1
    if step_index < n:
        # First sample starting a 10-sample window for which every point stays
        # within 2% of the instantaneous MPP. -1 means not found.
        tol = np.maximum(0.02 * pref[step_index:], 1.0)
        for j in range(step_index, max(step_index, n - 9)):
            if j + 10 <= n and np.all(np.abs(p[j:j+10] - pref[j:j+10]) <= tol[j-step_index:j-step_index+10]):
                settle = j - step_index
                break
    final_n = min(30, n)
    return {"algorithm": run["algorithm"], "scenario": run["name"],
            "tracking_efficiency_pct": float(run["efficiency_pct"]),
            "steady_window_ripple_w": ripple,
            "settling_after_step_samples": int(settle),
            "final_power_mean_w": float(np.mean(p[-final_n:])),
            "final_mpp_mean_w": float(np.mean(pref[-final_n:])),
            "global_peak_detected": bool(np.mean(p[-final_n:]) >= 0.98 * np.mean(pref[-final_n:]))}


def _save_run_csv(run: dict, filename: str):
    path = OUT / filename
    n_irr = run["irradiance_w_m2"].shape[1]
    labels = "module" if n_irr == 3 else "substring"
    headers = ["sample", "voltage_v", "current_measured_a", "current_true_a", "power_true_w",
               "power_measured_w", "mpp_reference_w"]
    headers += [f"irradiance_{labels}{j+1}_w_m2" for j in range(n_irr)]
    headers += ["temperature_c"]
    with path.open("w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f); wr.writerow(headers)
        for k in range(len(run["power_w"])):
            wr.writerow([k, run["voltage_v"][k], run["current_a"][k], run["true_current_a"][k],
                         run["power_w"][k], run["measured_power_w"][k], run["mpp_power_w"][k],
                         *run["irradiance_w_m2"][k], run["temp_c"][k]])


def _save_results(rows, filename="comparison_results.csv"):
    if not rows: return
    with (OUT / filename).open("w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wr.writeheader(); wr.writerows(rows)


def plot_pv_curves():
    # Plot the static I-V and P-V curves for the most important irradiance states:
    # uniform illumination, reduced irradiance, whole-module shading, and a single
    # shaded substring. These visuals make the multi-peak behavior much easier to
    # reason about than inspecting raw numbers alone.
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
    profiles = [
        ("Uniform 1000 W/m²", curve_for((1000,1000,1000), 25.0)),
        ("Uniform 600 W/m²", curve_for((600,600,600), 25.0)),
        ("Whole module shaded: 1000/1000/300", curve_for((1000,1000,300), 25.0)),
        ("One substring shaded: 8×1000 + 300", make_substring_string_iv_curve([1000]*8+[300], 25.0, module=MODULE, points=8000)),
    ]
    peak_rows = []
    for label, c in profiles:
        axs[0].plot(c.voltage_v, c.current_a, label=label)
        axs[1].plot(c.voltage_v, c.power_w, label=label)
        for p in c.mpps():
            peak_rows.append({"profile":label, **p})
            axs[1].scatter([p["voltage_v"]], [p["power_w"]], s=20)
    axs[0].set(xlabel="PV string voltage (V)", ylabel="PV string current (A)", title="I–V characteristics")
    axs[1].set(xlabel="PV string voltage (V)", ylabel="PV string power (W)", title="P–V characteristics")
    for ax in axs:
        ax.grid(True, alpha=0.3); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(OUT / "pv_characteristics.png", dpi=180); plt.close(fig)
    _save_results(peak_rows, "pv_mpp_points.csv")
    return peak_rows


def plot_runs(runs, scenario, filename, title):
    fig, axs = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    colors = {"P&O": "tab:blue", "Incremental Conductance": "tab:orange", "PSO": "tab:green"}
    for run in runs:
        label = run["algorithm"]
        axs[0].plot(run["power_w"], label=label, color=colors.get(label))
        if label == "P&O":
            axs[0].plot(run["mpp_power_w"], "k--", lw=1.2, label="Theoretical MPP")
        axs[1].plot(run["voltage_v"], label=label, color=colors.get(label))
    axs[0].set_ylabel("PV power (W)"); axs[0].set_title(title)
    axs[1].set_ylabel("PV voltage (V)"); axs[1].set_xlabel("Controller sample")
    for ax in axs: ax.grid(True, alpha=0.3); ax.legend(fontsize=8, ncol=2)
    fig.tight_layout(); fig.savefig(OUT / filename, dpi=180); plt.close(fig)


def run_baseline_and_extensions():
    # This is the main experiment driver. It validates the PV model, runs the
    # baseline MPPT comparisons, then performs the sensitivity and robustness
    # studies used to characterize algorithm performance under different climate and
    # shading conditions.
    start_time = time.time()
    # Datasheet sanity check at STC.
    validation = validate_module_stc(MODULE)
    (OUT / "module_stc_validation.json").write_text(json.dumps(validation, indent=2), encoding="utf-8")
    plot_pv_curves()
    algs = ["P&O", "Incremental Conductance", "PSO"]
    step_runs, shade_runs = [], []
    rows = []
    for alg in algs:
        kwargs = {"seed": 3, "po_step": 0.7, "inc_deadband": 0.005,
                  "pso_particles": 6, "pso_iterations": 12}
        if alg == "PSO": kwargs["seed"] = 4
        run = run_profile("irradiance_step", step_irradiance, step_temp, alg,
                          n_samples=400, start_v=0.55 * MAX_V, change_indices=(200,), **kwargs)
        step_runs.append(run); _save_run_csv(run, f"irradiance_step_{alg.replace(' ','_').replace('&','and')}.csv")
        rows.append(metrics(run, step_index=200))
        run2 = run_profile("partial_shading", shaded_irradiance, step_temp, alg,
                           n_samples=400, start_v=0.90 * N_MODULES * MODULE.v_oc_v, **kwargs)
        shade_runs.append(run2); _save_run_csv(run2, f"partial_shading_{alg.replace(' ','_').replace('&','and')}.csv")
        rows.append(metrics(run2, step_index=400))
    plot_runs(step_runs, "step", "irradiance_step_comparison.png",
              "MPPT under a 1000 → 600 W/m² irradiance step (ideal voltage-actuation benchmark)")
    plot_runs(shade_runs, "shade", "partial_shading_comparison.png",
              "MPPT under whole-module partial shading (ideal voltage-actuation benchmark)")

    # Additional realistic pattern: one of the nine 18-cell substrings is shaded.
    subshade_profile = lambda k: tuple([1000.0] * 8 + [300.0])
    subshade_curve = curve_for(subshade_profile(0), 25.0, irradiance_mode="substring")
    subshade_runs = []
    for alg in algs:
        kwargs = {"seed": 5 if alg == "PSO" else 3, "po_step": 0.7,
                  "inc_deadband": 0.005, "pso_particles": 6, "pso_iterations": 12}
        run3 = run_profile("one_substring_shaded", subshade_profile, step_temp, alg,
                           n_samples=400, start_v=0.95*subshade_curve.voc_v,
                           irradiance_mode="substring", **kwargs)
        subshade_runs.append(run3)
        _save_run_csv(run3, f"partial_substring_shading_{alg.replace(' ','_').replace('&','and')}.csv")
        rows.append(metrics(run3, step_index=400))
    plot_runs(subshade_runs, "substring", "partial_substring_shading_comparison.png",
              "MPPT with one 18-cell substring shaded (ideal voltage-actuation benchmark)")
    _save_results(rows)

    # PSO stochastic repeatability test on whole-module partial shading.
    shade_success = []
    shade_eta = []
    for seed in range(30):
        run = run_profile("partial_shading", shaded_irradiance, step_temp, "PSO",
                          n_samples=400, start_v=0.90 * N_MODULES * MODULE.v_oc_v,
                          seed=seed, pso_particles=6, pso_iterations=12)
        m = metrics(run, step_index=400)
        shade_success.append(int(m["global_peak_detected"]))
        shade_eta.append(float(m["tracking_efficiency_pct"]))
    pso_stats = {"seeds": 30, "global_peak_successes": int(sum(shade_success)),
                 "global_peak_success_rate_pct": float(np.mean(shade_success) * 100),
                 "tracking_efficiency_mean_pct": float(np.mean(shade_eta)),
                 "tracking_efficiency_sd_pct": float(np.std(shade_eta, ddof=1))}
    (OUT / "pso_30_seed_statistics.json").write_text(json.dumps(pso_stats, indent=2), encoding="utf-8")

    # Sensitivity analysis: P&O voltage step and IncCond deadband.
    sens_rows = []
    for step in [0.2, 0.5, 0.7, 1.0, 2.0]:
        run = run_profile("uniform_1000", lambda k:(1000,1000,1000), step_temp,
                          "P&O", n_samples=300, start_v=0.50*MAX_V, po_step=step, seed=10)
        last = run["power_w"][-50:]
        sens_rows.append({"algorithm":"P&O", "parameter":"step_v", "value":step,
                          "tracking_efficiency_pct":run["efficiency_pct"],
                          "last50_ripple_w":float(last.max()-last.min())})
    for db in [0.0, 0.001, 0.005, 0.02, 0.05]:
        run = run_profile("uniform_1000", lambda k:(1000,1000,1000), step_temp,
                          "Incremental Conductance", n_samples=300, start_v=0.50*MAX_V,
                          po_step=0.7, inc_deadband=db, seed=10)
        last = run["power_w"][-50:]
        sens_rows.append({"algorithm":"Incremental Conductance", "parameter":"deadband_a_per_v", "value":db,
                          "tracking_efficiency_pct":run["efficiency_pct"],
                          "last50_ripple_w":float(last.max()-last.min())})
    _save_results(sens_rows, "algorithm_sensitivity.csv")

    # Temperature sweep under uniform STC irradiance, all modules at the same cell temperature.
    temp_rows = []
    for temp in [0, 25, 50, 75]:
        curve = curve_for((1000,1000,1000), temp)
        mpp = curve.mpp()
        temp_rows.append({"cell_temperature_c":temp, **mpp})
    _save_results(temp_rows, "temperature_sweep.csv")
    fig, ax = plt.subplots(figsize=(7,4.4))
    ax.plot([r["cell_temperature_c"] for r in temp_rows], [r["power_w"] for r in temp_rows], marker="o")
    ax.set(xlabel="Cell temperature (°C)", ylabel="Maximum power (W)", title="Temperature effect at 1000 W/m²")
    ax.grid(True, alpha=0.3); fig.tight_layout(); fig.savefig(OUT/"temperature_effect.png", dpi=180); plt.close(fig)

    # Current sensor-noise test on a uniform-irradiance profile.
    noise_rows = []
    for noise in [0.0, 0.02, 0.05, 0.10]:
        for alg in ["P&O", "Incremental Conductance"]:
            run = run_profile("uniform_noise", lambda k:(1000,1000,1000), step_temp,
                              alg, n_samples=300, start_v=0.50*MAX_V,
                              current_noise_std_a=noise, seed=11, po_step=0.7, inc_deadband=0.005)
            noise_rows.append({"algorithm":alg,"current_noise_std_a":noise,
                               "tracking_efficiency_pct":run["efficiency_pct"],
                               "final50_ripple_w":float(run["power_w"][-50:].max()-run["power_w"][-50:].min())})
    _save_results(noise_rows, "current_noise_sensitivity.csv")

    # Slowly changing irradiance; all modules see the same ramp.
    ramp_runs=[]; ramp_rows=[]
    for alg in algs:
        run = run_profile("irradiance_ramp", lambda k:ramp_irradiance(k, 400), step_temp,
                          alg, n_samples=400, start_v=0.55*MAX_V, seed=8,
                          po_step=0.7, inc_deadband=0.005, pso_particles=6, pso_iterations=12)
        ramp_runs.append(run)
        ramp_rows.append({"algorithm":alg,"tracking_efficiency_pct":run["efficiency_pct"],
                          "mean_absolute_power_error_w":float(np.mean(np.abs(run["mpp_power_w"]-run["power_w"])))})
    plot_runs(ramp_runs, "ramp", "irradiance_ramp_comparison.png",
              "MPPT under a gradual irradiance ramp (ideal voltage-actuation benchmark)")
    _save_results(ramp_rows, "irradiance_ramp_results.csv")

    summary = {
        "module": MODULE.name,
        "series_modules": N_MODULES,
        "model_scope": "single-diode electrical model + three substring bypass clamps per module + ideal voltage-controlled interface",
        "baseline_rows": rows,
        "pso_30_seed_partial_shading": pso_stats,
        "module_stc_validation": validation,
        "run_seconds": time.time() - start_time,
        "note": "The baseline algorithms use an ideal voltage-actuation interface; this is not a switching-level or averaged-converter simulation."
    }
    (OUT / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({"module_validation":validation,"pso_30_seed":pso_stats,"run_seconds":summary["run_seconds"]}, indent=2))

if __name__ == "__main__":
    run_baseline_and_extensions()
