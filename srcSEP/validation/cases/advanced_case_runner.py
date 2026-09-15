"""Application-level orchestration and scoring for validation cases CV06-CV12.

The linked executable owns every numerical sample.  This module only converts
reviewed JSON into native name/value arguments, launches ``amps --test CVxx``,
runs the case-local independent reference in a separate process, and evaluates
saved evidence.  Keeping those responsibilities separate prevents a Python
helper from being mistaken for the model under validation.
"""

from __future__ import annotations

import math
from pathlib import Path
import subprocess
import sys
import time
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from linked_case_common import (atomic_json, finish_result, load_input, metric,
                                read_csv, run_linked_model, sha256, write_csv)


def _csv(values: Iterable[Any]) -> str:
    """Serialize reviewed numeric arrays for the strict native line protocol."""
    return ",".join(str(value) for value in values)


def _reference(case_id: str, input_path: Path, output_path: Path,
               timeout: Optional[float]) -> None:
    """Execute the independent formula outside the linked model process."""
    script = Path(__file__).resolve().parent / case_id / "reference_solution.py"
    completed = subprocess.run(
        [sys.executable, str(script), "--input", str(input_path), "--output",
         str(output_path)], text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, check=False, timeout=timeout)
    if completed.returncode:
        raise RuntimeError(f"{case_id} reference failed: {completed.stdout.strip()}")


def _plot(case_id: str, output_dir: Path, title: str, x: Sequence[float],
          model: Sequence[float], exact: Sequence[float], xlabel: str,
          ylabel: str) -> List[Path]:
    """Render the same saved comparison series to PNG and vector EPS."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(2, 1, figsize=(8.2, 7.0), sharex=True,
                                gridspec_kw={"height_ratios": [3, 1]})
    axes[0].plot(x, exact, "k-", linewidth=2, label="independent reference")
    axes[0].plot(x, model, "o", markersize=4, label="linked srcSEP/AMPS")
    axes[0].set_title(f"{case_id}: {title}")
    axes[0].set_ylabel(ylabel)
    axes[0].grid(alpha=0.25)
    axes[0].legend()
    axes[1].axhline(0, color="black", linewidth=1)
    axes[1].plot(x, [a-b for a, b in zip(model, exact)], "o-", markersize=3)
    axes[1].set_xlabel(xlabel)
    axes[1].set_ylabel("residual")
    axes[1].grid(alpha=0.25)
    figure.tight_layout()
    paths: List[Path] = []
    for suffix in ("png", "eps"):
        path = output_dir / f"{case_id}_comparison.{suffix}"
        figure.savefig(path, dpi=180 if suffix == "png" else None)
        paths.append(path)
    plt.close(figure)
    return paths


def _order(errors: Dict[float, float]) -> float:
    """Calculate the endpoint convergence order on positive finite errors."""
    items = sorted((step, error) for step, error in errors.items()
                   if step > 0 and error > 0 and math.isfinite(error))
    if len(items) < 2:
        return 99.0 if items and items[0][1] < 1e-14 else 0.0
    fine, coarse = items[0], items[-1]
    return math.log(coarse[1] / fine[1]) / math.log(coarse[0] / fine[0])


def _mean(values: Sequence[float]) -> float:
    return sum(values) / len(values)


def _linear_slope(x: Sequence[float], y: Sequence[float]) -> float:
    xbar, ybar = _mean(x), _mean(y)
    denominator = sum((value-xbar)**2 for value in x)
    return sum((a-xbar)*(b-ybar) for a, b in zip(x, y)) / denominator


def _solution(case_id: str, output_dir: Path, x: Sequence[float],
              model: Sequence[float], exact: Sequence[float], units: str,
              title: str, xlabel: str, ylabel: str) -> List[Path]:
    """Persist the conventional analytical series consumed by run_tests.py."""
    path = output_dir / f"{case_id}_solution.csv"
    write_csv(path, ("coordinate", "numerical", "analytical", "residual", "units"),
              ({"coordinate": xx, "numerical": yy, "analytical": zz,
                "residual": yy-zz, "units": units}
               for xx, yy, zz in zip(x, model, exact)))
    return [path] + _plot(case_id, output_dir, title, x, model, exact,
                          xlabel, ylabel)


def _normal_cdf(value: float) -> float:
    return 0.5 * (1.0 + math.erf(value / math.sqrt(2.0)))


def _ig_cdf(time_s: float, boundary: float, drift: float, kappa: float) -> float:
    """CDF of first passage for ds=u dt+sqrt(2 kappa)dW to s=L."""
    if time_s <= 0:
        return 0.0
    scale = math.sqrt(2.0*kappa*time_s)
    return (_normal_cdf((drift*time_s-boundary)/scale) +
            math.exp(drift*boundary/kappa) *
            _normal_cdf(-(drift*time_s+boundary)/scale))


def _arguments(case_id: str, config: Dict[str, Any]) -> List[str]:
    """Map nested, human-readable inputs onto the native strict protocol."""
    p, n = config["physics"], config["numerics"]
    seed = ["--campaign-seed", str(n["campaign_seed"])]
    if case_id == "CV06":
        return ["--d0-per-s", str(p["d0_per_s"]), "--epsilon", str(p["initial_perturbation"]), "--modes", _csv(p["modes"]), "--final-time-s", str(n["final_time_s"]), "--time-steps-s", _csv(n["time_steps_s"]), "--sample-count", str(n["sample_count"]), "--particle-count", str(n["particle_count"]), "--seed-count", str(n["seed_count"])] + seed
    if case_id == "CV07":
        return ["--speed-m-per-s", str(p["speed_m_per_s"]), "--switching-rates-per-s", _csv(p["switching_rates_per_s"]), "--normalized-times", _csv(n["normalized_times"]), "--particle-count", str(n["particle_count"]), "--seed-count", str(n["seed_count"]), "--profile-bins", str(n["profile_bins"])] + seed
    if case_id == "CV08":
        return ["--boundary-m", str(p["boundary_m"]), "--kappa-m2-per-s", str(p["kappa_m2_per_s"]), "--drifts-m-per-s", _csv(p["drifts_m_per_s"]), "--time-steps-s", _csv(n["time_steps_s"]), "--maximum-time-s", str(n["maximum_time_s"]), "--particle-count", str(n["particle_count"]), "--seed-count", str(n["seed_count"])] + seed
    if case_id == "CV09":
        return ["--compression-ratios", _csv(p["compression_ratios"]), "--upstream-speed-m-per-s", str(p["upstream_speed_m_per_s"]), "--particle-speed-m-per-s", str(p["particle_speed_m_per_s"]), "--kappa-upstream-m2-per-s", str(p["kappa_upstream_m2_per_s"]), "--kappa-downstream-m2-per-s", str(p["kappa_downstream_m2_per_s"]), "--injection-momentum-kg-m-per-s", str(p["injection_momentum_kg_m_per_s"]), "--particle-count", str(n["particle_count"])] + seed
    if case_id == "CV10":
        return ["--resolutions", _csv(n["resolutions"]), "--spectral-bins", str(n["spectral_bins"]), "--duration-s", str(n["duration_s"]), "--cfl", str(n["cfl"])] + seed
    if case_id == "CV11":
        return ["--spectral-bins", str(n["spectral_bins"]), "--active-bin", str(n["active_bin"]), "--initial-energy-j", str(p["initial_energy_j"]), "--growth-rate-per-s", str(p["growth_rate_per_s"]), "--damping-rate-per-s", str(p["damping_rate_per_s"]), "--sinusoidal-rate-per-s", str(p["sinusoidal_rate_per_s"]), "--angular-frequency-per-s", str(p["angular_frequency_per_s"]), "--final-time-s", str(n["final_time_s"]), "--time-steps-s", _csv(n["time_steps_s"])] + seed
    if case_id == "CV12":
        return ["--particle-speed-m-per-s", str(p["particle_speed_m_per_s"]), "--alfven-speed-m-per-s", str(p["alfven_speed_m_per_s"]), "--initial-wave-energy-j", str(p["initial_wave_energy_j"]), "--target-exchange-j", str(p["target_exchange_j"]), "--time-steps-s", _csv(n["time_steps_s"]), "--particle-counts", _csv(n["particle_counts"]), "--steps", str(n["steps"])] + seed
    raise ValueError(f"unsupported advanced case: {case_id}")


def _cv06(config: Dict[str, Any], rows: List[Dict[str, str]], _: List[Dict[str, str]], out: Path):
    p, n, a = config["physics"], config["numerics"], config["acceptance"]
    d0, eps = float(p["d0_per_s"]), float(p["initial_perturbation"])
    finest = min(map(float, n["time_steps_s"]))
    rate_errors, leakage_z, iso_z, dt_errors = [], 0.0, 0.0, {}
    x, model, exact = [], [], []
    for mode in map(int, p["modes"]):
        for dt in map(float, n["time_steps_s"]):
            selected = [r for r in rows if int(r["initial_mode"]) == mode and int(r["measured_mode"]) == mode and float(r["dt_s"]) == dt]
            by_time = sorted({float(r["time_s"]) for r in selected})
            means = [_mean([float(r["coefficient"]) for r in selected if float(r["time_s"]) == t]) for t in by_time]
            slope = _linear_slope(by_time, [math.log(max(abs(v), 1e-30)) for v in means])
            expected_rate = mode*(mode+1)*d0
            rate_errors.append(abs((-slope)/expected_rate-1.0))
            endpoint = eps/(2*mode+1)*math.exp(-expected_rate*by_time[-1])
            dt_errors[dt] = max(dt_errors.get(dt, 0.0), abs(means[-1]-endpoint))
            if dt == finest:
                x.extend([mode+t/by_time[-1] for t in by_time])
                model.extend(means)
                exact.extend([eps/(2*mode+1)*math.exp(-expected_rate*t) for t in by_time])
        all_finest = [r for r in rows if int(r["initial_mode"]) == mode and float(r["dt_s"]) == finest]
        total_samples = int(n["particle_count"])*int(n["seed_count"])
        sigma = 1.0/math.sqrt(total_samples)
        leakage_z = max(leakage_z, *(abs(float(r["coefficient"]))/sigma for r in all_finest if int(r["measured_mode"]) not in (0, mode)))
        iso_z = max(iso_z, *(abs(float(r["coefficient"])-1.0)/sigma for r in all_finest if int(r["measured_mode"]) == 0))
    reflections = max(int(r["boundary_reflections"]) for r in rows)
    metrics = [metric("maximum_decay_rate_relative_error", max(rate_errors), a["decay_rate_relative_error_max"], "<=", "dimensionless"), metric("maximum_mode_leakage_standard_errors", leakage_z, a["leakage_standard_errors_max"], "<=", "standard-errors"), metric("isotropic_mode_standard_errors", iso_z, a["isotropic_mode_standard_errors_max"], "<=", "standard-errors"), metric("temporal_refinement_order", _order(dt_errors), a["refinement_order_min"], ">=", "dimensionless"), metric("production_mu_boundary_reflections", reflections, a["mu_boundary_reflections_min"], ">=", "count")]
    return metrics, _solution("CV06", out, x, model, exact, "legendre_coefficient", "Legendre modal decay", x and "mode + normalized time" or "time", "modal coefficient")


def _cv07(config: Dict[str, Any], rows: List[Dict[str, str]], _: List[Dict[str, str]], out: Path):
    p, n, a = config["physics"], config["numerics"], config["acceptance"]
    speed, particles = float(p["speed_m_per_s"]), int(n["particle_count"])
    moments = [r for r in rows if r["row_type"] == "moment"]
    support = max(float(r["outside_support"]) for r in moments)
    front_z = msd_z = event_rel = late_rel = 0.0
    x, model, exact = [], [], []
    for rate in map(float, p["switching_rates_per_s"]):
        for q in map(float, n["normalized_times"]):
            group = [r for r in moments if float(r["rate_per_s"]) == rate and float(r["normalized_time"]) == q]
            front = math.exp(-q)
            measured_front = _mean([float(r["front_mass"]) for r in group])
            front_se = math.sqrt(max(front*(1-front), 1e-30)/(particles*len(group)))
            front_z = max(front_z, abs(measured_front-front)/front_se)
            t = q/rate
            expected_msd = speed*speed*(t/rate-(1-math.exp(-2*rate*t))/(2*rate*rate))
            measured_msd = _mean([float(r["msd_m2"]) for r in group])
            seed_sd = math.sqrt(sum((float(r["msd_m2"])-measured_msd)**2 for r in group)/max(len(group)-1, 1))
            msd_z = max(msd_z, abs(measured_msd-expected_msd)/max(seed_sd/math.sqrt(len(group)), expected_msd/math.sqrt(particles*len(group))))
            event_rel = max(event_rel, abs(_mean([float(r["event_mean"]) for r in group])/q-1.0))
            if q >= 10:
                k_eff = measured_msd/(2*t)
                late_rel = max(late_rel, abs(k_eff/(speed*speed/(2*rate))-1.0))
            if rate == float(p["switching_rates_per_s"][0]):
                x.append(q); model.append(measured_msd); exact.append(expected_msd)
    metrics = [metric("probability_outside_causal_support", support, a["support_probability_max"], "<=", "probability"), metric("front_mass_standard_errors", front_z, a["front_mass_standard_errors_max"], "<=", "standard-errors"), metric("msd_standard_errors", msd_z, a["msd_standard_errors_max"], "<=", "standard-errors"), metric("late_diffusion_coefficient_relative_error", late_rel, a["late_diffusion_relative_error_max"], "<=", "dimensionless"), metric("poisson_event_mean_relative_error", event_rel, a["event_mean_relative_error_max"], "<=", "dimensionless")]
    return metrics, _solution("CV07", out, x, model, exact, "m2", "telegraph mean-square displacement", "normalized time nu*t", "MSD [m2]")


def _cv08(config: Dict[str, Any], rows: List[Dict[str, str]], _: List[Dict[str, str]], out: Path):
    p, n, a = config["physics"], config["numerics"], config["acceptance"]
    boundary, kappa, maximum = float(p["boundary_m"]), float(p["kappa_m2_per_s"]), float(n["maximum_time_s"])
    finest = min(map(float, n["time_steps_s"])); failures = 0; moment_z = 0.0
    overshoots: Dict[float, float] = {}; correct_ks = wrong_ks = 0.0
    x, model, exact = [], [], []
    for drift in map(float, p["drifts_m_per_s"]):
        for dt in map(float, n["time_steps_s"]):
            overshoots[dt] = max(overshoots.get(dt, 0.0), _mean([float(r["overshoot_m"]) for r in rows if float(r["drift_m_per_s"]) == drift and float(r["dt_s"]) == dt and r["arrived"] == "1"]))
            for seed in sorted({r["seed"] for r in rows}):
                group = [r for r in rows if float(r["drift_m_per_s"]) == drift and float(r["dt_s"]) == dt and r["seed"] == seed]
                arrivals = sorted(float(r["arrival_time_s"]) for r in group if r["arrived"] == "1")
                total = len(group); ks = 0.0
                for index, value in enumerate(arrivals, 1):
                    expected = _ig_cdf(value, boundary, drift, kappa)
                    ks = max(ks, abs(index/total-expected), abs((index-1)/total-expected))
                ks = max(ks, abs(len(arrivals)/total-_ig_cdf(maximum, boundary, drift, kappa)))
                if ks > math.sqrt(-0.5*math.log(a["ks_alpha"]/2.0)/total): failures += 1
                if dt == finest:
                    correct_ks = max(correct_ks, ks)
                    wrong = max(abs((index/total)-_ig_cdf(value, boundary, drift, 1.5*kappa)) for index, value in enumerate(arrivals, 1))
                    wrong_ks = max(wrong_ks, wrong)
                    if drift == max(map(float, p["drifts_m_per_s"])):
                        mean = _mean(arrivals); expected_mean = boundary/drift
                        expected_sd = math.sqrt(2*kappa*boundary/drift**3)
                        moment_z = max(moment_z, abs(mean-expected_mean)/(expected_sd/math.sqrt(len(arrivals))))
            if dt == finest:
                group = [r for r in rows if float(r["drift_m_per_s"]) == drift and float(r["dt_s"]) == dt]
                for value in [maximum*i/50 for i in range(1, 51)]:
                    x.append(value + 100*drift); model.append(sum(r["arrived"] == "1" and float(r["arrival_time_s"]) <= value for r in group)/len(group)); exact.append(_ig_cdf(value, boundary, drift, kappa))
    metrics = [metric("seed_level_ks_failures", failures, int(a["seeds_allowed_to_fail"])*len(p["drifts_m_per_s"])*len(n["time_steps_s"]), "<=", "count"), metric("arrival_mean_standard_errors", moment_z, a["moment_standard_errors_max"], "<=", "standard-errors"), metric("overshoot_refinement_order", _order(overshoots), a["overshoot_refinement_order_min"], ">=", "dimensionless"), metric("negative_control_ks_ratio", wrong_ks/max(correct_ks, 1e-30), a["negative_control_ratio_min"], ">=", "dimensionless")]
    return metrics, _solution("CV08", out, x, model, exact, "cdf", "first-passage CDF", "time with drift offset", "arrival CDF")


def _cv09(config: Dict[str, Any], rows: List[Dict[str, str]], _: List[Dict[str, str]], out: Path):
    p, n, a = config["physics"], config["numerics"], config["acceptance"]
    slope_error = time_error = accounting = 0.0; minimum_length = math.inf
    x, model, exact = [], [], []
    for ratio in map(float, p["compression_ratios"]):
        group = sorted([r for r in rows if float(r["compression_ratio"]) == ratio], key=lambda r: float(r["momentum_kg_m_per_s"]))
        lo, hi = n["fit_quantiles"]
        first, last = int(lo*len(group)), int(hi*len(group))
        subset = group[first:last]
        momenta = [float(r["momentum_kg_m_per_s"]) for r in subset]
        # Ranks are already known from the sorted slice.  Constructing them
        # arithmetically is both unambiguous for duplicate momenta and O(N),
        # unlike repeated list.index searches that would make this 50k-sample
        # scientific fit accidentally quadratic.
        survival = [(len(group)-index)/len(group)
                    for index in range(first, last)]
        alpha = -_linear_slope([math.log(v) for v in momenta], [math.log(v) for v in survival])
        qfit, qexact = alpha+3.0, 3.0*ratio/(ratio-1.0)
        slope_error = max(slope_error, abs(qfit-qexact))
        u1, u2 = float(p["upstream_speed_m_per_s"]), float(p["upstream_speed_m_per_s"])/ratio
        tacc = 3/(u1-u2)*(float(p["kappa_upstream_m2_per_s"])/u1+float(p["kappa_downstream_m2_per_s"])/u2)
        gain = 4*(u1-u2)/(3*float(p["particle_speed_m_per_s"]))
        observed = _mean([float(r["acceleration_time_s"])/max(math.log(float(r["momentum_kg_m_per_s"])/float(p["injection_momentum_kg_m_per_s"])), math.log1p(gain)) for r in group if int(r["cycles"]) > 0])
        time_error = max(time_error, abs(observed/tacc-1))
        accounting = max(accounting, abs(sum(float(r["escaped_weight"]) for r in group)/int(n["particle_count"])-1))
        minimum_length = min(minimum_length, *(float(r["upstream_diffusion_length_m"]) for r in group), *(float(r["downstream_diffusion_length_m"]) for r in group))
        x.append(ratio); model.append(qfit); exact.append(qexact)
    metrics = [metric("phase_space_slope_absolute_error", slope_error, a["momentum_slope_absolute_error_max"], "<=", "index"), metric("acceleration_time_relative_error", time_error, a["acceleration_time_relative_error_max"], "<=", "dimensionless"), metric("number_accounting_error", accounting, a["number_accounting_error_max"], "<=", "dimensionless"), metric("minimum_diffusion_length", minimum_length, a["resolved_diffusion_length_min_m"], ">=", "m")]
    return metrics, _solution("CV09", out, x, model, exact, "phase_space_index", "planar DSA spectrum", "compression ratio", "q")


def _sine_average(index: int, cells: int, shift: float) -> float:
    dx=1/cells; left=index*dx-shift; right=(index+1)*dx-shift
    return 1+0.25*(math.cos(2*math.pi*left)-math.cos(2*math.pi*right))/(2*math.pi*dx)


def _cv10(config: Dict[str, Any], rows: List[Dict[str, str]], _: List[Dict[str, str]], out: Path):
    n, a = config["numerics"], config["acceptance"]; finest=max(map(int,n["resolutions"])); errors={}; invariant=negative=0.0
    x=[]; model=[]; exact=[]
    for resolution in map(int,n["resolutions"]):
        error2=ref2=0.0
        for scenario in ("fixed","expanding-area","moving-grid"):
            for branch in ("plus","minus"):
                group=[r for r in rows if r["scenario"]==scenario and r["branch"]==branch and int(r["resolution"])==resolution and int(r["spectral_bin"])==0]
                cells=len(group); shift=float(n["duration_s"])*(1 if branch=="plus" else -1); fraction=1/(0.5*int(n["spectral_bins"])*(int(n["spectral_bins"])+1))
                for r in group:
                    ref=_sine_average(int(r["cell"]),cells,shift)/cells*fraction; value=float(r["wave_action_j"])
                    error2+=(value-ref)**2; ref2+=ref**2
                    if resolution==finest and scenario=="fixed": x.append(float(r["left_m"])); model.append(value); exact.append(ref)
                invariant=max(invariant,abs(float(group[0]["total_branch_energy_j"])/fraction-1))
                negative=max(negative,float(group[0]["negative_cells"]))
        errors[1/resolution]=math.sqrt(error2/max(ref2,1e-30))
    metrics=[metric("finest_relative_L2_error",errors[1/finest],a["relative_l2_error_max"],"<=","dimensionless"),metric("spatial_refinement_order",_order(errors),a["refinement_order_min"],">=","dimensionless"),metric("wave_action_relative_invariant_error",invariant,a["relative_invariant_error_max"],"<=","dimensionless"),metric("negative_cell_count",negative,a["negative_cell_count_max"],"<=","count")]
    order=sorted(range(len(x)),key=lambda i:x[i]); x=[x[i] for i in order]; model=[model[i] for i in order]; exact=[exact[i] for i in order]
    return metrics,_solution("CV10",out,x,model,exact,"J","spectral wave advection","cell left coordinate [m]","wave action [J]")


def _cv11(config: Dict[str, Any], rows: List[Dict[str, str]], _: List[Dict[str, str]], out: Path):
    p,n,a=config["physics"],config["numerics"],config["acceptance"]; active=int(n["active_bin"]); final=float(n["final_time_s"]); errors={}; constant=inactive=negative=cancel=0.0
    x=[];model=[];exact=[]
    def ref(s,t):
        if s=="growth": exponent=2*float(p["growth_rate_per_s"])*t
        elif s=="damping": exponent=-2*float(p["damping_rate_per_s"])*t
        elif s=="cancellation": exponent=0
        else: exponent=2*(float(p["sinusoidal_rate_per_s"])/float(p["angular_frequency_per_s"])*(1-math.cos(float(p["angular_frequency_per_s"])*t))-float(p["damping_rate_per_s"])*t)
        return float(p["initial_energy_j"])*math.exp(exponent)
    for dt in map(float,n["time_steps_s"]):
        endpoint=[r for r in rows if float(r["dt_s"])==dt and float(r["time_s"])==final and int(r["spectral_bin"])==active]
        errors[dt]=max(abs(float(r["wave_energy_j"])/ref(r["scenario"],final)-1) for r in endpoint if r["scenario"]=="sign-change")
        constant=max(constant,*(abs(float(r["wave_energy_j"])/ref(r["scenario"],final)-1) for r in endpoint if r["scenario"] in ("growth","damping")))
    inactive=max(abs(float(r["wave_energy_j"])) for r in rows if int(r["spectral_bin"])!=active); negative=sum(int(r["positive"])==0 for r in rows)
    cancel=max(abs(float(r["wave_energy_j"])/float(p["initial_energy_j"])-1) for r in rows if r["scenario"]=="cancellation" and int(r["spectral_bin"])==active)
    finest=min(map(float,n["time_steps_s"])); plot=[r for r in rows if r["scenario"]=="sign-change" and float(r["dt_s"])==finest and int(r["spectral_bin"])==active]
    for r in plot:x.append(float(r["time_s"]));model.append(float(r["wave_energy_j"]));exact.append(ref("sign-change",x[-1]))
    metrics=[metric("constant_rate_relative_error",constant,a["constant_rate_relative_error_max"],"<=","dimensionless"),metric("sign_change_relative_error",errors[finest],a["sign_change_relative_error_max"],"<=","dimensionless"),metric("temporal_refinement_order",_order(errors),a["temporal_refinement_order_min"],">=","dimensionless"),metric("inactive_bin_energy",inactive,a["inactive_bin_energy_max_j"],"<=","J"),metric("negative_cell_count",negative,a["negative_cell_count_max"],"<=","count"),metric("cancellation_relative_change",cancel,a["cancellation_relative_change_max"],"<=","dimensionless")]
    return metrics,_solution("CV11",out,x,model,exact,"J","time-dependent resonant growth/damping","time [s]","active-bin energy [J]")


def _cv12(config: Dict[str, Any], rows: List[Dict[str, str]], _: List[Dict[str, str]], out: Path):
    a=config["acceptance"]; coupled=ledger=uncoupled=0.0; negative=math.inf
    x=[];model=[];exact=[]
    for scenario in ("particle-only","wave-only","coupled","negative-control"):
        groups={(r["particle_count"],r["dt_s"]) for r in rows if r["scenario"]==scenario}
        for key in groups:
            group=sorted([r for r in rows if r["scenario"]==scenario and (r["particle_count"],r["dt_s"])==key],key=lambda r:int(r["step"])); initial=float(group[0]["total_energy_j"]); drift=abs(float(group[-1]["total_energy_j"])/initial-1)
            if scenario=="coupled":
                coupled=max(coupled,drift); ledger=max(ledger,*(abs(float(r["step_residual_j"]))/initial for r in group));
                if key==max(groups):
                    x.extend(float(r["time_s"]) for r in group);model.extend(float(r["total_energy_j"]) for r in group);exact.extend(initial for r in group)
            elif scenario in ("particle-only","wave-only"): uncoupled=max(uncoupled,drift)
            else: negative=min(negative,drift)
    metrics=[metric("coupled_total_relative_drift",coupled,a["coupled_total_relative_drift_max"],"<=","dimensionless"),metric("ledger_step_residual_relative",ledger,a["ledger_step_residual_relative_max"],"<=","dimensionless"),metric("uncoupled_relative_drift",uncoupled,a["uncoupled_relative_drift_max"],"<=","dimensionless"),metric("negative_control_relative_drift",negative,a["negative_control_relative_drift_min"],">=","dimensionless")]
    return metrics,_solution("CV12",out,x,model,exact,"J","closed particle-wave energy","time [s]","total energy [J]")


SCORERS = {"CV06":_cv06,"CV07":_cv07,"CV08":_cv08,"CV09":_cv09,"CV10":_cv10,"CV11":_cv11,"CV12":_cv12}


def run_advanced_case(case_id: str, *, source_root: Path, input_path: Path,
                      output_dir: Path, executable: Path,
                      timeout: Optional[float]) -> Dict[str, Any]:
    """Run one CV06-CV12 linked model/reference/acceptance pipeline."""
    started=time.monotonic(); config=load_input(input_path,case_id); output_dir.mkdir(parents=True,exist_ok=True)
    resolved=output_dir/"resolved_input.json"; atomic_json(resolved,config)
    native=run_linked_model(case_id=case_id,arguments=_arguments(case_id,config),source_root=source_root,output_dir=output_dir,executable=executable,timeout=timeout)
    reference=output_dir/f"{case_id}_reference.csv"; _reference(case_id,resolved,reference,timeout)
    metrics,derived=SCORERS[case_id](config,read_csv(native["model"]),read_csv(reference),output_dir)
    provenance=output_dir/"provenance.json"; atomic_json(provenance,{"schema":"srcsep-validation-provenance-v1","case_id":case_id,"execution":"linked-srcsep-amps","sha256":{"executable":sha256(executable),"input":sha256(resolved),"model":sha256(native["model"]),"reference":sha256(reference)},"reference_implementation":str(Path(__file__).resolve().parent/case_id/"reference_solution.py")})
    artifacts=[resolved,native["manifest"],native["model"],native["report"],native["junit"],native["log"],reference,provenance]+derived
    return finish_result(case_id=case_id,started=started,seed=int(config["numerics"]["campaign_seed"]),input_path=input_path,executable=executable,metrics=metrics,artifacts=artifacts,message=f"{case_id} linked-application validation passed")
