#!/usr/bin/env python3
"""Production-linked M0/M1 execution and independent acceptance checks.

The executable writes raw SI-valued receipts by calling srcMoon production
entry points.  This module owns orchestration and independent reference
calculations only: it never supplies forces or states to the production mover.
That separation prevents a copied test formula from making an implementation
agree with itself.  Every run preserves the command stream, receipt, result,
and a manifest that identifies the executable, source, inputs, kernels, random
seed, and rank/thread layout.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
from typing import Callable

REPO = Path(__file__).resolve().parents[2]
TEST_ROOT = Path(__file__).resolve().parent
SPICE_ROOT = Path("/home/vtenishe/SPICE/Kernels")
SPICE_FILES = (
    "NAIF/naif0012.tls", "NAIF/de421.bsp", "NAIF/pck00011.tpc",
    "NAIF/moon_pa_de421_1900-2050.bpc", "NAIF/moon_080317.tf",
    "NAIF/moon_assoc_me.tf", "OTHER/Moon.LSO.tf",
    "OTHER/GSE.tf",
)
G = 6.67300e-11
M_MOON = 7.3477e22
M_EARTH = 5.9736e24
M_SUN = 1.9891e30
R_MOON = 1737.10e3
PASS, FAIL, ERROR, SKIPPED = "PASS", "FAIL", "ERROR", "SKIPPED"


class MissingPrerequisite(RuntimeError):
    """A declared campaign prerequisite/configuration is unavailable."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_receipt(path: Path) -> dict[str, object]:
    """Decode the deliberately simple, locale-independent C++ receipt."""
    values: dict[str, object] = {}
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if "=" not in line:
            raise ValueError(f"receipt line {number} has no '=': {line!r}")
        key, text = line.split("=", 1)
        try:
            values[key] = ([float(item) for item in text.split(",")]
                           if "," in text else float(text))
        except ValueError:
            values[key] = text
    return values


def norm(vector: list[float]) -> float:
    return math.sqrt(sum(value * value for value in vector))


def sub(left: list[float], right: list[float]) -> list[float]:
    return [left[i] - right[i] for i in range(3)]


def cross(left: list[float], right: list[float]) -> list[float]:
    return [left[1] * right[2] - left[2] * right[1],
            left[2] * right[0] - left[0] * right[2],
            left[0] * right[1] - left[1] * right[0]]


def force_without_radiation(x: list[float], v: list[float], receipt: dict[str, object]) -> list[float]:
    """Independent SI evaluator for the frozen-state I06/I10 ODE.

    The expression is specified in the roadmap and U03/U04 contracts.  It is
    evaluated from the raw SPICE vectors emitted by the executable, not by
    calling Moon::OrbitalDynamics.  The fixture remains in lunar shadow, so no
    empirical Na radiation-pressure table enters this reference trajectory.
    """
    radius = norm(x)
    acceleration = [-G * M_MOON * item / radius**3 for item in x]
    sun = [float(receipt["moon_heliocentric_distance_m"]), 0.0, 0.0]
    earth = list(receipt["earth_position_lso_m"])  # type: ignore[arg-type]
    for body, mass in ((sun, M_SUN), (earth, M_EARTH)):
        delta = sub(body, x)
        body_radius, delta_radius = norm(body), norm(delta)
        for component in range(3):
            acceleration[component] += G * mass * (
                delta[component] / delta_radius**3 -
                body[component] / body_radius**3)
    omega = list(receipt["lso_angular_velocity_rad_s"])  # type: ignore[arg-type]
    centrifugal = cross(omega, cross(omega, x))
    coriolis = cross(omega, v)
    return [acceleration[i] - centrifugal[i] - 2.0 * coriolis[i]
            for i in range(3)]


def rk4(initial_x: list[float], initial_v: list[float], duration: float,
        receipt: dict[str, object], steps: int = 32768) -> tuple[list[float], list[float]]:
    """High-resolution independent RK4 reference for the frozen-state ODE."""
    state = initial_x + initial_v
    dt = duration / steps

    def derivative(candidate: list[float]) -> list[float]:
        return candidate[3:] + force_without_radiation(candidate[:3], candidate[3:], receipt)

    for _ in range(steps):
        k1 = derivative(state)
        k2 = derivative([state[i] + 0.5 * dt * k1[i] for i in range(6)])
        k3 = derivative([state[i] + 0.5 * dt * k2[i] for i in range(6)])
        k4 = derivative([state[i] + dt * k3[i] for i in range(6)])
        state = [state[i] + dt * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) / 6.0
                 for i in range(6)]
    return state[:3], state[3:]


def rk4_uniform_lorentz(receipt: dict[str, object], steps: int = 32768) -> tuple[list[float], list[float]]:
    """Independent Na+ reference including uniform q(E+v×B)/m and Moon gravity."""
    electric = list(receipt["electric_field_v_m"])  # type: ignore[arg-type]
    magnetic = list(receipt["magnetic_field_t"])  # type: ignore[arg-type]
    charge_to_mass = float(receipt["ion_charge_c"]) / float(receipt["ion_mass_kg"])
    state = list(receipt["initial_x_m"]) + list(receipt["initial_v_m_s"])  # type: ignore[arg-type]
    dt = float(receipt["total_time_s"]) / steps

    def derivative(candidate: list[float]) -> list[float]:
        x, v = candidate[:3], candidate[3:]
        radius = norm(x)
        vx_b = cross(v, magnetic)
        acceleration = [charge_to_mass * (electric[i] + vx_b[i]) -
                        G * M_MOON * x[i] / radius**3 for i in range(3)]
        return v + acceleration

    for _ in range(steps):
        k1 = derivative(state)
        k2 = derivative([state[i] + 0.5*dt*k1[i] for i in range(6)])
        k3 = derivative([state[i] + 0.5*dt*k2[i] for i in range(6)])
        k4 = derivative([state[i] + dt*k3[i] for i in range(6)])
        state = [state[i] + dt*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6.0
                 for i in range(6)]
    return state[:3], state[3:]


def scaled_state_error(x: list[float], v: list[float], reference: tuple[list[float], list[float]]) -> float:
    x_ref, v_ref = reference
    return math.sqrt((norm(sub(x, x_ref)) / max(R_MOON, norm(x_ref)))**2 +
                     (norm(sub(v, v_ref)) / max(1.0, norm(v_ref)))**2)


def relative_change(left: float, right: float) -> float:
    """Symmetric scale used by the frozen I11 successive-grid contract."""
    scale = max(abs(left), abs(right))
    return 0.0 if scale == 0.0 else abs(left - right) / scale


def moon_setting(path: Path, key: str) -> str:
    """Return one explicit scalar Moon setting or fail closed.

    The production C++ parser remains authoritative for the complete receipt.
    This small preflight prevents expensive linked tests from being launched
    with a known-ineligible profile or full-volume output request; it does not
    duplicate or infer any resolution values.
    """
    pattern = rf"(?mi)^\s*{re.escape(key)}\s*=\s*([^!\s]+)\s*$"
    matches = re.findall(pattern,
                         path.read_text(encoding="utf-8"))
    if len(matches) != 1:
        raise ValueError(f"expected exactly one {key} in {path}")
    return matches[0].lower()


def mesh_profile(path: Path) -> str:
    """Return the explicit QUICK/FULL/CUSTOM provenance label."""
    return moon_setting(path, "mesh_profile")


def derive_i11_input(template: Path, destination: Path,
                     near_surface_km: int) -> None:
    """Create one archived I11 receipt without editing the reviewed template.

    The 10 km template is the frozen physical case. Near-surface, spherical-
    boundary, angular, and polar targets are multiplied together because they
    jointly determine boundary-adjacent resolution. Far-field and outer-domain
    targets remain fixed: scaling the complete 400-lunar-radius tail would mix
    a near-surface convergence question with a separate far-domain sensitivity
    question and produces an intractable three-dimensional block count.
    Physical extents, altitude comparison points, epoch, sources, and
    observation geometry remain identical. Every expected assignment must
    occur exactly once, otherwise derivation fails closed.
    """
    factor = near_surface_km / 10.0
    text = template.read_text(encoding="utf-8")
    scalar_keys = ("surface_mesh_resolution_m",
                   "mesh_near_surface_cell_size_m",
                   "mesh_polar_resolution_m")
    for key in scalar_keys:
        pattern = rf"(?m)^({re.escape(key)}\s*=\s*)([0-9]+(?:\.[0-9]*)?)\s*$"
        matches = re.findall(pattern, text)
        if len(matches) != 1:
            raise ValueError(f"expected exactly one {key} in {template}")
        scaled = float(matches[0][1]) * factor
        text = re.sub(pattern, lambda match: match.group(1) + f"{scaled:.17g}",
                      text)

    angular_pattern = r"(?m)^(mesh_angular_resolution_m_by_sza\s*=\s*)(.+)$"
    angular_matches = re.findall(angular_pattern, text)
    if len(angular_matches) != 1:
        raise ValueError("expected one mesh_angular_resolution_m_by_sza table")
    pairs: list[str] = []
    for field in angular_matches[0][1].split(","):
        angle, value = field.strip().split(":", 1)
        pairs.append(f"{angle}:{float(value) * factor:.17g}")
    text = re.sub(angular_pattern,
                  lambda match: match.group(1) + ",".join(pairs), text)
    destination.write_text(text, encoding="utf-8")


def read_qualified_metrics(path: Path) -> dict[str, dict[str, float | str]]:
    """Load the production-sampled vector and reject duplicate/nonfinite rows."""
    metrics: dict[str, dict[str, float | str]] = {}
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            name = row["metric"]
            if name in metrics:
                raise ValueError(f"duplicate I11 metric {name!r} in {path}")
            numeric = {key: float(row[key]) for key in
                       ("mean_si", "relative_ci95_half_width", "lag1")}
            if not all(math.isfinite(value) for value in numeric.values()):
                raise ValueError(f"nonfinite I11 metric {name!r} in {path}")
            metrics[name] = {"classification": row["classification"], **numeric}
    if not metrics:
        raise ValueError(f"empty I11 metric file {path}")
    return metrics


def descendant_pids(root_pid: int) -> set[int]:
    """Return the live Linux process subtree used by mpiexec and AMPS."""
    descendants = {root_pid}
    changed = True
    while changed:
        changed = False
        for entry in Path("/proc").iterdir():
            if not entry.name.isdigit() or int(entry.name) in descendants:
                continue
            try:
                fields = (entry / "stat").read_text(encoding="utf-8").split()
                parent = int(fields[3])
            except (OSError, ValueError, IndexError):
                continue
            if parent in descendants:
                descendants.add(int(entry.name)); changed = True
    return descendants


def proc_rss_kib(pid: int) -> int:
    try:
        for line in Path(f"/proc/{pid}/status").read_text(
                encoding="utf-8").splitlines():
            if line.startswith("VmRSS:"):
                return int(line.split()[1])
    except (OSError, ValueError, IndexError):
        pass
    return 0


def host_available_kib() -> int:
    try:
        for line in Path("/proc/meminfo").read_text(encoding="utf-8").splitlines():
            if line.startswith("MemAvailable:"):
                return int(line.split()[1])
    except (OSError, ValueError, IndexError):
        pass
    return -1


def run_monitored(command: list[str], cwd: Path, environment: dict[str, str],
                  stdout_path: Path, stderr_path: Path,
                  resource_path: Path) -> int:
    """Run without a resource/time cap while recording the complete process tree."""
    with stdout_path.open("w", encoding="utf-8") as stdout, \
            stderr_path.open("w", encoding="utf-8") as stderr, \
            resource_path.open("w", encoding="utf-8") as resources:
        resources.write("wall_seconds,process_tree_rss_kib,host_mem_available_kib\n")
        started = time.monotonic()
        process = subprocess.Popen(command, cwd=cwd, env=environment,
                                   stdout=stdout, stderr=stderr)
        try:
            while process.poll() is None:
                rss = sum(proc_rss_kib(pid) for pid in descendant_pids(process.pid))
                resources.write(f"{time.monotonic()-started:.6f},{rss},"
                                f"{host_available_kib()}\n")
                resources.flush()
                time.sleep(1.0)
        except KeyboardInterrupt:
            # A resource operator may interrupt a run before the host OOM
            # killer acts. Terminate the complete mpiexec process group through
            # its launcher and return a normal nonzero code so the caller can
            # write result.json=ERROR instead of losing the structured status.
            process.terminate()
            try:
                process.wait(timeout=30.0)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
            return 130
        return process.returncode


def check(test_id: str, receipt: dict[str, object], work: Path) -> tuple[str, dict[str, object], str]:
    """Apply only predeclared physical/numerical criteria for one frozen ID."""
    metrics: dict[str, object] = {}
    if test_id == "I01":
        return PASS, metrics, "all required provenance fields and hashes were captured"
    if test_id == "I02":
        conditions = {
            "neutral_na_only": receipt.get("species_count") == 1.0 and receipt.get("na_species_index") == 0.0 and float(receipt.get("na_plus_species_index", 0.0)) < 0.0,
            "surface_callbacks": receipt.get("sphere_present") == 1.0 and receipt.get("sphere_injection_callback") == 1.0 and receipt.get("sphere_interaction_callback") == 1.0,
            "source_positive": float(receipt.get("total_source_rate_s", -1.0)) >= float(receipt.get("impact_source_rate_s", 0.0)) > 0.0,
            "photo_gate": receipt.get("photo_lit_allowed") == 1.0 and receipt.get("photo_shadow_allowed") == 0.0,
        }
        metrics["conditions"] = conditions
        return (PASS if all(conditions.values()) else FAIL), metrics, "production dispatch audit"
    if test_id == "I03":
        mass_error = abs(float(receipt["species_0_mass_kg"]) - 3.82e-26) / 3.82e-26
        conditions = {"mass_relative_error": mass_error,
                      "neutral_charge": receipt["species_0_charge_c"] == 0.0,
                      "positive_timestep": float(receipt["species_0_global_dt_s"]) > 0.0,
                      "positive_weight": float(receipt["species_0_global_weight"]) > 0.0,
                      "reference_particles": receipt["max_reference_injected_particles"] == 20.0}
        metrics.update(conditions)
        ok = mass_error <= 1.0e-14 and all(bool(value) for key, value in conditions.items() if key != "mass_relative_error")
        return (PASS if ok else FAIL), metrics, "species/weight/time-step initialization"
    if test_id == "I04":
        diff = work / "test_Moon.diff"
        size = diff.stat().st_size if diff.exists() else -1
        metrics["test_Moon_diff_bytes"] = size
        return (PASS if size == 0 else FAIL), metrics, "canonical sampler regression"
    if test_id == "I05":
        x0, v0 = list(receipt["initial_x_m"]), list(receipt["initial_v_m_s"])  # type: ignore[arg-type]
        duration = float(receipt["total_time_s"])
        mu = G * M_MOON
        angular_rate = math.sqrt(mu / norm(x0)**3)
        angle = angular_rate * duration
        reference = ([x0[0]*math.cos(angle), x0[0]*math.sin(angle), 0.0],
                     [-v0[1]*math.sin(angle), v0[1]*math.cos(angle), 0.0])
        error = scaled_state_error(list(receipt["final_x_m"]), list(receipt["final_v_m_s"]), reference)  # type: ignore[arg-type]
        metrics["scaled_state_error"] = error
        return (PASS if error <= 1.0e-3 else FAIL), metrics, "closed-form circular lunar orbit"
    if test_id in ("I06", "I10"):
        initial_x = list(receipt["initial_x_m"]); initial_v = list(receipt["initial_v_m_s"])  # type: ignore[arg-type]
        reference = rk4(initial_x, initial_v, float(receipt["total_time_s"]), receipt)
        if test_id == "I06":
            error = scaled_state_error(list(receipt["final_x_m"]), list(receipt["final_v_m_s"]), reference)  # type: ignore[arg-type]
            metrics["scaled_state_error"] = error
            return (PASS if error <= 2.0e-4 else FAIL), metrics, "independent frozen-SPICE RK4 trajectory"
        errors = []
        for level in range(4):
            errors.append(scaled_state_error(list(receipt[f"level_{level}_final_x_m"]), list(receipt[f"level_{level}_final_v_m_s"]), reference))  # type: ignore[arg-type]
        ratios = [errors[i] / errors[i+1] for i in range(3)]
        metrics.update({"scaled_state_errors": errors, "successive_error_ratios": ratios,
                        "observed_orders": [math.log(value, 2.0) for value in ratios]})
        ok = all(errors[i+1] < errors[i] for i in range(3)) and all(value >= 3.0 for value in ratios) and errors[-1] <= 0.01
        return (PASS if ok else FAIL), metrics, "second-order convergence and <=1% finest-state error"
    if test_id == "I07":
        shadow_expected = force_without_radiation(list(receipt["shadow_initial_x_m"]), list(receipt["shadow_initial_v_m_s"]), receipt)  # type: ignore[arg-type]
        shadow_error = norm(sub(list(receipt["shadow_acceleration_m_s2"]), shadow_expected))  # type: ignore[arg-type]
        lit_expected = force_without_radiation(list(receipt["lit_initial_x_m"]), list(receipt["lit_initial_v_m_s"]), receipt)  # type: ignore[arg-type]
        radiation = sub(list(receipt["lit_acceleration_m_s2"]), lit_expected)  # type: ignore[arg-type]
        metrics.update({"shadow_acceleration_abs_error_m_s2": shadow_error,
                        "lit_radiation_residual_m_s2": norm(radiation),
                        "lit_radiation_x_component_m_s2": radiation[0]})
        ok = shadow_error <= 1.0e-12 and radiation[0] < 0.0 and abs(radiation[1]) <= 1.0e-10 and abs(radiation[2]) <= 1.0e-10
        return (PASS if ok else FAIL), metrics, "shadow gate and sun-line radiation residual"
    if test_id == "I08":
        # The pure uniform-field gyro solution is the analytic limiting oracle;
        # srcMoon correctly retains lunar gravity, so the scored endpoint uses
        # a high-resolution independent integration of the exact declared
        # q(E+v cross B)/m - GM r/r^3 ODE rather than deleting gravity in a
        # test-only production mode.
        reference = rk4_uniform_lorentz(receipt)
        errors = []
        for level in range(4):
            errors.append(scaled_state_error(list(receipt[f"level_{level}_final_x_m"]), list(receipt[f"level_{level}_final_v_m_s"]), reference))  # type: ignore[arg-type]
        ratios = [errors[i]/errors[i+1] for i in range(3)]
        metrics.update({"scaled_state_errors": errors, "successive_ratios": ratios,
                        "gyrofrequency_rad_s": receipt["gyrofrequency_rad_s"]})
        ok = all(errors[i+1] < errors[i] for i in range(3)) and all(value >= 3.0 for value in ratios) and errors[-1] <= 0.01
        return (PASS if ok else FAIL), metrics, "uniform-field Lorentz plus lunar-gravity reference convergence"
    if test_id == "I09":
        samples = int(receipt["samples"]); accounted = int(receipt["stuck"] + receipt["reemitted"] + receipt["other_status"])
        metrics.update({"samples": samples, "accounted": accounted,
                        "inward_reemission": int(receipt["inward_reemission"]),
                        "physical_escape_energy_j_kg": receipt["physical_escape_specific_energy_j_kg"],
                        "numerical_loss_energy_j_kg": receipt["numerical_domain_loss_specific_energy_j_kg"],
                        "physical_escape_mover_status": receipt["physical_escape_mover_status"],
                        "numerical_loss_mover_status": receipt["numerical_domain_loss_mover_status"],
                        "compiled_left_domain_status": receipt["particle_left_domain_status"],
                        "compiled_motion_finished_status": receipt["particle_motion_finished_status"]})
        # Both trajectories must reach the same production domain-loss return
        # path; the independently meaningful energy signs distinguish their
        # physical budget labels.  Compare against the constants emitted by
        # the same compiled AMPS executable: copying a numeric macro value into
        # this runner would silently assume that all configurations agree.
        left_domain = receipt["particle_left_domain_status"]
        motion_finished = receipt["particle_motion_finished_status"]
        ok = (accounted == samples and receipt["other_status"] == 0.0 and
              receipt["inward_reemission"] == 0.0 and
              float(receipt["physical_escape_specific_energy_j_kg"]) > 0.0 and
              float(receipt["numerical_domain_loss_specific_energy_j_kg"]) <= 0.0 and
              receipt["physical_escape_mover_status"] == left_domain and
              receipt["numerical_domain_loss_mover_status"] == left_domain and
              left_domain != motion_finished)
        return (PASS if ok else FAIL), metrics, "surface accounting plus energy-based escape/domain classification"
    if test_id == "I12":
        # Reconstruct the declared piecewise-linear selector solely from raw
        # input values echoed by the executable.  This Python evaluator does
        # not call the C++ production kernel, while the realized-cell half of
        # the gate is measured from the linked AMPS tree.
        knot_count = int(receipt["mesh_knot_count"])
        knots = [(float(receipt[f"mesh_knot_{index}_sza_deg"]),
                  float(receipt[f"mesh_knot_{index}_resolution_m"]))
                 for index in range(knot_count)]
        terminator_width = float(receipt["mesh_terminator_half_width_deg"])
        polar_start = float(receipt["mesh_polar_cap_start_latitude_deg"])
        polar_resolution = float(receipt["mesh_polar_resolution_m"])
        profile = str(receipt["mesh_profile"]).lower()

        def expected_resolution(sza: float, latitude: float) -> float:
            right = next((index for index, item in enumerate(knots)
                          if item[0] >= sza), len(knots) - 1)
            if knots[right][0] == sza or right == 0:
                value = knots[right][1]
            else:
                left = right - 1
                fraction = ((sza - knots[left][0]) /
                            (knots[right][0] - knots[left][0]))
                value = knots[left][1] + fraction * (
                    knots[right][1] - knots[left][1])
            terminator_value = next(value_m for angle, value_m in knots
                                    if angle == 90.0)
            if abs(sza - 90.0) <= terminator_width:
                value = min(value, terminator_value)
            if abs(latitude) >= polar_start:
                value = min(value, polar_resolution)
            return value

        samples: list[dict[str, float | bool]] = []
        for index in range(6):
            sza = float(receipt[f"sza_{index}_angle_deg"])
            latitude = float(receipt[f"sza_{index}_absolute_body_latitude_deg"])
            requested = float(receipt[f"sza_{index}_requested_resolution_m"])
            realized = float(receipt[f"sza_{index}_realized_cell_size_m"])
            expected = expected_resolution(sza, latitude)
            formula_error = abs(requested - expected) / expected
            realization_error = abs(realized - requested) / requested
            # AMPS refines a Cartesian parent into power-of-two children. A
            # requested target can therefore lie between two available sizes
            # even when the callback is exact. The roadmap explicitly permits
            # documented AMR-constraint exceptions. The first is the nearest
            # dyadic level. The second is deliberately narrower: at an input
            # knot adjoining a finer target, AMPS's enforced one-level neighbor
            # balance may refine the containing block once beyond that nearest
            # level. Its dyadic parent must still satisfy the original 25%
            # gate. This cannot excuse a coarse cell, a non-transition sample,
            # or unbounded over-refinement.
            nearest_dyadic = abs(math.log2(realized / requested)) <= 0.5
            discrete_exception = (realization_error > 0.25 and nearest_dyadic)
            knot_index = next(
                (knot_index for knot_index, (angle, _) in enumerate(knots)
                 if math.isclose(angle, sza, rel_tol=0.0, abs_tol=1.0e-12)),
                None,
            )
            finer_neighbor = False
            if knot_index is not None:
                if knot_index > 0:
                    finer_neighbor = knots[knot_index - 1][1] < requested
                if knot_index + 1 < len(knots):
                    finer_neighbor = (
                        finer_neighbor or
                        knots[knot_index + 1][1] < requested
                    )
            parent_realization_error = abs(2.0 * realized - requested) / requested
            balance_exception = (
                realization_error > 0.25 and
                realized < requested and
                finer_neighbor and
                parent_realization_error <= 0.25
            )
            samples.append({"sza_deg": sza,
                            "absolute_body_latitude_deg": latitude,
                            "requested_resolution_m": requested,
                            "expected_resolution_m": expected,
                            "realized_cell_size_m": realized,
                            "formula_relative_error": formula_error,
                            "realization_relative_error": realization_error,
                            "formula_pass": formula_error <= 1.0e-14,
                            "amr_discrete_level_exception": discrete_exception,
                            "amr_balance_transition_exception": balance_exception,
                            "dyadic_parent_realization_relative_error":
                                parent_realization_error,
                            "realization_pass": (realization_error <= 0.25 or
                                                 discrete_exception or
                                                 balance_exception)})
        pole_sza = float(receipt["north_pole_sza_deg"])
        pole_requested = float(receipt["north_pole_requested_resolution_m"])
        pole_realized = float(receipt["north_pole_realized_cell_size_m"])
        pole_expected = expected_resolution(pole_sza, 90.0)
        pole_realization_error = abs(pole_realized-pole_requested)/pole_requested
        pole_discrete_exception = (pole_realization_error > 0.25 and
                                   abs(math.log2(pole_realized /
                                                 pole_requested)) <= 0.5)
        pole = {"sza_deg": pole_sza,
                "requested_resolution_m": pole_requested,
                "expected_resolution_m": pole_expected,
                "realized_cell_size_m": pole_realized,
                "formula_relative_error": abs(pole_requested-pole_expected)/pole_expected,
                "realization_relative_error": pole_realization_error,
                "amr_discrete_level_exception": pole_discrete_exception}
        metrics.update({"mesh_profile": profile, "samples": samples,
                        "north_pole": pole})
        ok = (all(bool(item["formula_pass"]) and
                  bool(item["realization_pass"]) for item in samples) and
              pole["formula_relative_error"] <= 1.0e-14 and
              (pole["realization_relative_error"] <= 0.25 or
               pole_discrete_exception))
        if not ok:
            return FAIL, metrics, (
                "input-defined angular formula or linked realized-cell "
                "tolerance/declared dyadic AMR constraint failed")
        if profile != "full":
            return SKIPPED, metrics, (
                f"{profile.upper()} linked wiring checks passed, but only "
                "mesh_profile=full qualifies the M1 I12 gate")
        return PASS, metrics, (
            "FULL input-defined angular formula and linked realized-cell "
            "tolerance/declared dyadic AMR constraint passed")
    if test_id == "I13":
        # Geometry has a first-principles oracle and is always evaluated.  The
        # science baseline is deliberately not inferred from an initialized,
        # empty particle population: it may be frozen only after I11 sampling
        # passes and a human reviews the candidate CSV.
        area = float(receipt["sphere_area_m2"])
        radius = float(receipt["sphere_radius_m"])
        analytical_area = 4.0 * math.pi * radius * radius
        area_error = abs(area - analytical_area) / analytical_area
        metrics.update({
            "sphere_area_relative_error": area_error,
            "sphere_surface_element_count": int(receipt["sphere_surface_element_count"]),
            "amr_leaf_node_count": int(receipt["amr_leaf_node_count"]),
            "amr_allocated_leaf_block_count": int(receipt["amr_allocated_leaf_block_count"]),
            "species_count": int(receipt["species_count"]),
            "simulation_environment_mode": str(
                receipt["simulation_environment_mode"]),
            "mesh_near_surface_cell_size_m": float(
                receipt["mesh_near_surface_cell_size_m"]),
            "mesh_near_surface_extent_m": float(
                receipt["mesh_near_surface_extent_m"]),
            "mesh_far_field_cell_size_m": float(
                receipt["mesh_far_field_cell_size_m"]),
            "mesh_far_field_extent_m": float(
                receipt["mesh_far_field_extent_m"]),
            "mesh_outer_cell_size_m": float(
                receipt["mesh_outer_cell_size_m"]),
            "na_total_source_rate_s": float(receipt["na_total_source_rate_s"]),
            "na_impact_source_rate_s": float(receipt["na_impact_source_rate_s"]),
        })
        if area_error > 1.0e-10:
            return FAIL, metrics, "production sphere area exceeds the I13 analytical tolerance"
        baseline = TEST_ROOT / "integration-tests/I13/reference/I13_archived_baseline_metrics.csv"
        qualified = work.parent / "I11/5km/qualified_metrics.csv"
        i11_result = work.parent / "I11/result.json"
        i11_finest_input = work.parent / "I11/5km/application.in"
        if not qualified.is_file() or not i11_result.is_file():
            metrics.update({"baseline_path": str(baseline),
                            "qualified_i11_path": str(qualified)})
            return SKIPPED, metrics, (
                "sphere invariant passed; a passing finest-grid I11 artifact is unavailable")
        i11_status = json.loads(i11_result.read_text(encoding="utf-8")).get("status")
        if i11_status != PASS:
            return SKIPPED, metrics, (
                f"sphere invariant passed; I11 status is {i11_status}, not PASS")
        if (not i11_finest_input.is_file() or
                sha256(i11_finest_input) != sha256(work / "application.in")):
            return ERROR, metrics, (
                "I13 must use the exact passing I11 5 km derived receipt; "
                "the preserved input is missing or has a different hash")
        candidate = read_qualified_metrics(qualified)
        candidate_path = work / "I13_candidate_baseline_metrics.csv"
        shutil.copy2(qualified, candidate_path)
        # The roadmap freezes more than sampled science values: the realized
        # mesh counts, compiled species table, source totals, and exact input
        # configuration must also remain identifiable. Keep those deterministic
        # controls in a separate JSON object because species symbols and hashes
        # cannot be represented honestly as SI-valued rows in the science CSV.
        species_count = int(receipt["species_count"])
        configuration_candidate = {
            "schema": "srcMoon-I13-configuration-v1",
            "application_input_sha256": sha256(work / "application.in"),
            "surface_geometry": "sphere",
            "simulation_environment_mode": str(
                receipt["simulation_environment_mode"]),
            "sphere_radius_m": float(receipt["sphere_radius_m"]),
            "sphere_surface_element_count": int(
                receipt["sphere_surface_element_count"]),
            "amr_leaf_node_count": int(receipt["amr_leaf_node_count"]),
            "amr_allocated_leaf_block_count": int(
                receipt["amr_allocated_leaf_block_count"]),
            "mesh": {
                "near_surface_cell_size_m": float(
                    receipt["mesh_near_surface_cell_size_m"]),
                "near_surface_extent_m": float(
                    receipt["mesh_near_surface_extent_m"]),
                "far_field_cell_size_m": float(
                    receipt["mesh_far_field_cell_size_m"]),
                "far_field_extent_m": float(
                    receipt["mesh_far_field_extent_m"]),
                "outer_cell_size_m": float(receipt["mesh_outer_cell_size_m"]),
            },
            "species": [{
                "index": index,
                "symbol": str(receipt[f"species_{index}_symbol"]),
                "mass_kg": float(receipt[f"species_{index}_mass_kg"]),
                "charge_c": float(receipt[f"species_{index}_charge_c"]),
                "global_timestep_s": float(
                    receipt[f"species_{index}_global_timestep_s"]),
                "global_particle_weight": float(
                    receipt[f"species_{index}_global_particle_weight"]),
            } for index in range(species_count)],
            "sources": {
                "na_total_source_rate_s": float(receipt["na_total_source_rate_s"]),
                "na_impact_source_rate_s": float(
                    receipt["na_impact_source_rate_s"]),
            },
        }
        configuration_candidate_path = work / "I13_candidate_configuration.json"
        configuration_candidate_path.write_text(
            json.dumps(configuration_candidate, indent=2, sort_keys=True) + "\n",
            encoding="utf-8")
        configuration_baseline = (TEST_ROOT /
            "integration-tests/I13/reference/I13_archived_configuration.json")
        metrics.update({"candidate_path": str(candidate_path),
                        "candidate_sha256": sha256(candidate_path),
                        "configuration_candidate_path": str(
                            configuration_candidate_path),
                        "configuration_candidate_sha256": sha256(
                            configuration_candidate_path),
                        "baseline_path": str(baseline),
                        "configuration_baseline_path": str(
                            configuration_baseline)})
        if not baseline.is_file() or not configuration_baseline.is_file():
            return SKIPPED, metrics, (
                "sphere invariant and I11 qualification passed; science and configuration candidates await independent review and freeze")

        frozen = read_qualified_metrics(baseline)
        metrics["baseline_sha256"] = sha256(baseline)
        frozen_configuration = json.loads(
            configuration_baseline.read_text(encoding="utf-8"))
        metrics["configuration_baseline_sha256"] = sha256(
            configuration_baseline)
        if configuration_candidate != frozen_configuration:
            return FAIL, metrics, (
                "realized mesh/species/source/input configuration differs from "
                "the reviewed I13 configuration baseline")
        if set(candidate) != set(frozen):
            return FAIL, metrics, "frozen and current I13 scoring vectors differ"
        failures: list[dict[str, object]] = []
        comparisons: list[dict[str, object]] = []
        for name in sorted(candidate):
            current, reference = candidate[name], frozen[name]
            if current["classification"] != reference["classification"]:
                return FAIL, metrics, f"classification changed for I13 metric {name}"
            numerical = 0.10 if current["classification"] == "pointwise" else 0.05
            allowance = (numerical +
                         float(current["relative_ci95_half_width"]) +
                         float(reference["relative_ci95_half_width"]))
            difference = relative_change(float(current["mean_si"]),
                                         float(reference["mean_si"]))
            row = {"metric": name, "relative_change": difference,
                   "combined_allowance": allowance,
                   "pass": difference <= allowance}
            comparisons.append(row)
            if not row["pass"]:
                failures.append(row)
        metrics.update({"comparison_count": len(comparisons),
                        "failed_comparisons": failures})
        return (PASS if not failures else FAIL), metrics, (
            "sphere invariant and frozen qualified-science baseline comparison")
    return ERROR, metrics, f"no acceptance implementation for {test_id}"


def git_value(*arguments: str) -> str:
    completed = subprocess.run(["git", *arguments], cwd=REPO, text=True,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    return completed.stdout.strip() if completed.returncode == 0 else "UNAVAILABLE"


def manifest(test_id: str, executable: Path, command: list[str], seed: object | None,
             application_input: Path | None = None,
             mpi_ranks: int = 1) -> dict[str, object]:
    """Capture deterministic controls; manifests contain no timestamps by design.

    ``mpi_ranks`` is passed explicitly by every configurable linked runner so
    the retained layout is the layout actually supplied to ``mpiexec``.  The
    default preserves the frozen one-rank M0 behavior for existing Python
    callers that do not request a decomposition study.
    """
    files = [REPO / "input/moon.input", REPO / "input/species.input",
             REPO / "srcMoon/development_and_validation.md",
             Path("/data/vtenishe/moon_validation_data/development_and_validation.md")]
    if test_id == "I08":
        files.append(REPO / "input/moon-i08.input")
    if application_input is not None:
        files.append(application_input)
    if test_id == "I13":
        frozen_baseline = (TEST_ROOT /
            "integration-tests/I13/reference/I13_archived_baseline_metrics.csv")
        if frozen_baseline.is_file():
            files.append(frozen_baseline)
        frozen_configuration = (TEST_ROOT /
            "integration-tests/I13/reference/I13_archived_configuration.json")
        if frozen_configuration.is_file():
            files.append(frozen_configuration)
    # Only the seven kernels declared by the frozen DE421 M0 configuration are
    # dependencies.  Hashing unrelated files elsewhere in the installation
    # would make the manifest host-dependent without improving provenance.
    kernels = ([] if test_id in ("I05", "I08") else
               [SPICE_ROOT / relative for relative in SPICE_FILES
                if (SPICE_ROOT / relative).is_file()])
    compiler = subprocess.run(["mpicxx", "--version"], text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              check=False).stdout.splitlines()
    return {
        "schema": "srcMoon-linked-run-manifest-v2", "test_id": test_id,
        "source_revision": git_value("rev-parse", "HEAD"),
        "source_dirty_paths": git_value("status", "--short").splitlines(),
        "executable": {"path": str(executable.resolve()), "sha256": sha256(executable)},
        "command": command, "mpi_ranks": mpi_ranks, "openmp_threads": 1,
        "compiler": compiler[0] if compiler else "UNAVAILABLE",
        "build_configuration_hashes": {
            str(path): sha256(path) for path in (REPO / ".amps.conf", REPO / "Makefile.conf")
            if path.is_file()
        },
        "random_seed": seed, "species_fixture": "Na,NA_PLUS" if test_id == "I08" else "Na",
        "configuration": ("no-SPICE" if test_id in ("I05", "I08") else
                          ("canonical-test_Moon" if test_id == "I04" else "SPICE-DE421")),
        "inputs": {str(path): sha256(path) for path in files if path.is_file()},
        "external_data": {str(path.relative_to(SPICE_ROOT)): sha256(path) for path in kernels},
    }


def run_i11_campaign(executable: Path, output_root: Path,
                     application_input: Path | None,
                     mpi_ranks: int = 1,
                     campaign_mode: str = "qualification") -> int:
    """Run and score the frozen 40/20/10/5 km production campaign."""
    output = output_root / "I11"
    output.mkdir(parents=True, exist_ok=True)
    result: dict[str, object] = {"id": "I11", "status": ERROR,
                                "metrics": {}, "artifacts": {}}
    try:
        if application_input is None:
            raise ValueError("I11 requires --input with the reviewed M1 receipt")
        if campaign_mode == "regression":
            # I12 owns the single bounded linked production run in the default
            # regression campaign.  I11's adaptive statistical loop has no
            # defensible fixed runtime, so starting it here would make a normal
            # developer regression unexpectedly unbounded. Preserve a complete
            # non-execution record and never label it PASS.
            shutil.copy2(application_input.resolve(), output / "application.in")
            (output / "stdout.log").write_text(
                "I11 FULL convergence was not launched in REGRESSION mode.\n",
                encoding="utf-8")
            (output / "stderr.log").write_text("", encoding="utf-8")
            skipped_manifest = {
                "schema": "srcMoon-I11-regression-skip-v1",
                "source_revision": git_value("rev-parse", "HEAD"),
                "campaign_mode": campaign_mode,
                "execution": "not-run",
                "reason": "FULL I11 is reserved for qualification mode",
                "input": {"path": str(application_input.resolve()),
                          "sha256": sha256(application_input.resolve())},
                "mpi_ranks_requested": mpi_ranks,
            }
            (output / "run_manifest.json").write_text(
                json.dumps(skipped_manifest, indent=2, sort_keys=True) + "\n",
                encoding="utf-8")
            result.update(
                status=SKIPPED,
                reason=("REGRESSION mode does not launch the unbounded FULL "
                        "I11 convergence campaign; use --campaign-mode "
                        "qualification for M1 evidence"),
                metrics={"linked_model_executed": False,
                         "mpi_ranks_requested": mpi_ranks},
                artifacts={"input_snapshot": str(output / "application.in"),
                           "run_manifest": str(output / "run_manifest.json"),
                           "stdout": str(output / "stdout.log"),
                           "stderr": str(output / "stderr.log")})
            raise StopIteration
        if mesh_profile(application_input.resolve()) != "full":
            raise MissingPrerequisite(
                "I11 requires mesh_profile=full; QUICK is a wiring profile "
                "and cannot be used as numerical-convergence evidence")
        if moon_setting(application_input.resolve(),
                        "mesh_output_tecplot") != "false":
            raise MissingPrerequisite(
                "I11 requires mesh_output_tecplot=false; full-volume text "
                "visualization is not a science artifact and can exhaust disk")
        levels = (40, 20, 10, 5)
        shutil.copy2(application_input.resolve(), output / "application.in")
        (output / "stdout.log").write_text(
            "I11 aggregate runner: per-resolution model stdout is retained "
            "under 40km/, 20km/, 10km/, and 5km/.\n", encoding="utf-8")
        (output / "stderr.log").write_text(
            "I11 aggregate runner: per-resolution model stderr is retained "
            "beside each derived receipt.\n", encoding="utf-8")
        initial_manifest = {
            "schema": "srcMoon-I11-campaign-manifest-v1",
            "source_revision": git_value("rev-parse", "HEAD"),
            "template": {"path": str(application_input.resolve()),
                         "sha256": sha256(application_input.resolve())},
            "random_seed": 731993, "mpi_ranks": mpi_ranks,
            "openmp_threads": 1,
            "levels": {}}
        (output / "run_manifest.json").write_text(
            json.dumps(initial_manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8")
        level_metrics: dict[int, dict[str, dict[str, float | str]]] = {}
        level_receipts: dict[int, dict[str, object]] = {}
        level_summaries: dict[str, object] = {}
        environment = dict(os.environ); environment["OMP_NUM_THREADS"] = "1"

        for resolution in levels:
            level = output / f"{resolution}km"
            level.mkdir(parents=True, exist_ok=True)
            derived_input = level / "application.in"
            derive_i11_input(application_input.resolve(), derived_input,
                             resolution)
            receipt = level / "receipt.txt"
            stdout, stderr = level / "stdout.log", level / "stderr.log"
            command = ["mpiexec", "-n", str(mpi_ranks),
                       str(executable.resolve()),
                       "--moon-integration-test", "I11",
                       "--integration-output", str(receipt.resolve()),
                       "-input", str(derived_input.resolve())]
            level_manifest = manifest("I11", executable, command, 731993,
                                      derived_input, mpi_ranks)
            (level / "run_manifest.json").write_text(
                json.dumps(level_manifest, indent=2, sort_keys=True) + "\n",
                encoding="utf-8")
            runtime = level / "runtime"
            runtime.mkdir(exist_ok=True)
            returncode = run_monitored(
                command, runtime, environment, stdout, stderr,
                level / "process_resource_metrics.csv")
            if returncode != 0 or not receipt.is_file():
                raise RuntimeError(
                    f"I11 {resolution} km executable exited {returncode}; "
                    f"see {stdout} and {stderr}")
            parsed = parse_receipt(receipt)
            shutil.copy2(receipt, output / "receipt.txt")
            if parsed.get("test_id") != "I11":
                raise ValueError(f"I11 {resolution} km receipt has wrong ID")
            qualified = level / "qualified_metrics.csv"
            level_metrics[resolution] = read_qualified_metrics(qualified)
            level_receipts[resolution] = parsed
            initial_manifest["levels"][f"{resolution}km"] = {  # type: ignore[index]
                "input_sha256": sha256(derived_input),
                "manifest_sha256": sha256(level / "run_manifest.json"),
                "qualified_metrics_sha256": sha256(qualified),
            }
            (output / "run_manifest.json").write_text(
                json.dumps(initial_manifest, indent=2, sort_keys=True) + "\n",
                encoding="utf-8")

            vector_names = ("sun_position_lso_m", "earth_position_lso_m",
                            "lso_angular_velocity_rad_s")
            vector_drift = {name: norm(sub(list(parsed[f"step_2_{name}"]),
                                           list(parsed[f"step_1_{name}"])))
                            for name in vector_names}
            et_drift = abs(float(parsed["step_2_ephemeris_time_tdb_s"]) -
                           float(parsed["step_1_ephemeris_time_tdb_s"]))
            wiring_ok = (parsed.get("simulation_environment_mode") == "frozen" and
                         et_drift == 0.0 and
                         all(value == 0.0 for value in vector_drift.values()) and
                         float(parsed["particle_global_timestep_s"]) > 0.0 and
                         int(parsed["random_seed"]) == 731993 and
                         int(parsed["retained_independent_batch_count"]) >= 9)

            unresolved = []
            for name, metric in level_metrics[resolution].items():
                gate = 0.0333 if metric["classification"] == "pointwise" else 0.0167
                if float(metric["relative_ci95_half_width"]) > gate:
                    unresolved.append(name)
            level_status = PASS if wiring_ok and not unresolved else (
                SKIPPED if wiring_ok else FAIL)
            level_result = {
                "id": "I11", "resolution_km": resolution,
                "status": level_status,
                "reason": ("equilibrated independent production blocks satisfy CI gates"
                           if level_status == PASS else
                           "required scoring bins remain statistically unresolved"
                           if level_status == SKIPPED else
                           "frozen-environment or deterministic-control invariant failed"),
                "metrics": {"metric_count": len(level_metrics[resolution]),
                            "unresolved_bins": unresolved,
                            "ephemeris_time_drift_s": et_drift,
                            "external_vector_absolute_drifts": vector_drift,
                            "base_blocks": int(parsed["base_blocks_executed"]),
                            "batch_size_base_blocks": int(
                                parsed["retained_batch_size_base_blocks"]),
                            "independent_batch_count": int(
                                parsed["retained_independent_batch_count"])},
            }
            (level / "result.json").write_text(
                json.dumps(level_result, indent=2, sort_keys=True) + "\n",
                encoding="utf-8")
            level_summaries[f"{resolution}km"] = level_result
            if level_status != PASS:
                result.update(status=level_status,
                              reason=level_result["reason"],
                              metrics={"levels": level_summaries})
                raise StopIteration

        names = set(level_metrics[levels[0]])
        if any(set(level_metrics[level]) != names for level in levels[1:]):
            raise ValueError("I11 scoring-vector names differ across resolutions")
        convergence_rows: list[dict[str, object]] = []
        convergence_ok = True
        unresolved_pairs: list[str] = []
        for coarse, fine in zip(levels, levels[1:]):
            for name in sorted(names):
                left, right = level_metrics[coarse][name], level_metrics[fine][name]
                if left["classification"] != right["classification"]:
                    raise ValueError(f"I11 classification changed for {name}")
                classification = str(left["classification"])
                ci_gate = 0.0333 if classification == "pointwise" else 0.0167
                resolved = (float(left["relative_ci95_half_width"]) <= ci_gate and
                            float(right["relative_ci95_half_width"]) <= ci_gate)
                change = relative_change(float(left["mean_si"]),
                                         float(right["mean_si"]))
                tolerance = 0.10 if classification == "pointwise" else 0.05
                passed = resolved and change < tolerance
                if not resolved:
                    unresolved_pairs.append(f"{coarse}->{fine}:{name}")
                elif not passed:
                    convergence_ok = False
                convergence_rows.append({"coarse_km": coarse, "fine_km": fine,
                    "metric": name, "classification": classification,
                    "coarse_mean_si": left["mean_si"],
                    "fine_mean_si": right["mean_si"],
                    "relative_change": change, "tolerance": tolerance,
                    "both_ci_resolved": resolved, "pass": passed})
        with (output / "convergence_metrics.csv").open(
                "w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(
                convergence_rows[0].keys()))
            writer.writeheader(); writer.writerows(convergence_rows)

        status = SKIPPED if unresolved_pairs else (PASS if convergence_ok else FAIL)
        reason = ("required cross-resolution bins remain statistically unresolved"
                  if unresolved_pairs else
                  "all integrated and pointwise successive-resolution gates passed"
                  if convergence_ok else
                  "one or more evaluated successive-resolution gates failed")
        shutil.copy2(output / "5km" / "receipt.txt", output / "receipt.txt")
        # The campaign manifest binds the four child manifests and derived
        # receipts without hiding their individual commands or resource logs.
        campaign_manifest = {**initial_manifest,
            "levels": {f"{level}km": {
                "input_sha256": sha256(output / f"{level}km/application.in"),
                "manifest_sha256": sha256(output / f"{level}km/run_manifest.json"),
                "qualified_metrics_sha256": sha256(
                    output / f"{level}km/qualified_metrics.csv")}
                for level in levels}}
        (output / "run_manifest.json").write_text(
            json.dumps(campaign_manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8")
        result.update(status=status, reason=reason,
                      metrics={"levels": level_summaries,
                               "comparison_count": len(convergence_rows),
                               "unresolved_pairs": unresolved_pairs,
                               "failed_comparisons": sum(
                                   not bool(row["pass"]) and
                                   bool(row["both_ci_resolved"])
                                   for row in convergence_rows)},
                      artifacts={"convergence_metrics": str(
                          output / "convergence_metrics.csv")})
    except StopIteration:
        pass
    except MissingPrerequisite as exception:
        result.update(status=SKIPPED, reason=str(exception))
    except (OSError, ValueError, RuntimeError, KeyError, ZeroDivisionError) as exception:
        result.update(status=ERROR, reason=str(exception))
    (output / "result.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"{result['status']} I11: {result.get('reason', 'campaign stopped')}")
    return {PASS: 0, FAIL: 1, ERROR: 2, SKIPPED: 77}[str(result["status"])]


def run(test_id: str, executable: Path, output_root: Path,
        application_input: Path | None = None,
        mpi_ranks: int = 1,
        campaign_mode: str = "qualification") -> int:
    # Reject invalid layouts before constructing a command or manifest.  This
    # guard also protects direct Python callers, which bypass argparse below.
    if mpi_ranks < 1:
        raise ValueError("mpi_ranks must be a positive integer")
    if test_id == "I11":
        return run_i11_campaign(executable, output_root, application_input,
                                mpi_ranks, campaign_mode)
    output = output_root / test_id
    output.mkdir(parents=True, exist_ok=True)
    effective_input = application_input
    if test_id == "I13" and application_input is not None:
        # I13 freezes the numerically qualified finest I11 case, not the 10 km
        # human-authored template. Derive it with the same audited transformer
        # used by I11, preserve it as the I13 input snapshot, and later require
        # its hash to equal I11/5km/application.in before baseline scoring.
        effective_input = output / "application.in"
    receipt_path = output / "receipt.txt"
    stdout_path, stderr_path = output / "stdout.log", output / "stderr.log"
    command = ["mpiexec", "-n", str(mpi_ranks), str(executable.resolve()),
               "--moon-integration-test", test_id,
               "--integration-output", str(receipt_path.resolve())]
    if effective_input is not None:
        command.extend(["-input", str(effective_input.resolve())])
    # I09 freezes its callback fixture directly.  I11-I13 instead receive the
    # same value through the parsed production application receipt; retaining
    # it here makes that deterministic control visible in every linked manifest.
    seed = 731993 if test_id in ("I09", "I11", "I12", "I13") else None
    result: dict[str, object] = {"id": test_id, "status": ERROR, "metrics": {},
                                "artifacts": {"stdout": str(stdout_path), "stderr": str(stderr_path), "receipt": str(receipt_path)}}
    try:
        if test_id == "I13" and campaign_mode == "regression":
            if application_input is None:
                raise ValueError("I13 regression skip requires --input provenance")
            shutil.copy2(application_input.resolve(), output / "application.in")
            stdout_path.write_text(
                "I13 baseline freeze was not launched in REGRESSION mode.\n",
                encoding="utf-8")
            stderr_path.write_text("", encoding="utf-8")
            skip_manifest = {
                "schema": "srcMoon-I13-regression-skip-v1",
                "source_revision": git_value("rev-parse", "HEAD"),
                "campaign_mode": campaign_mode,
                "execution": "not-run",
                "reason": "I13 requires passing FULL I11 qualification evidence",
                "input": {"path": str(application_input.resolve()),
                          "sha256": sha256(application_input.resolve())},
                "mpi_ranks_requested": mpi_ranks,
            }
            (output / "run_manifest.json").write_text(
                json.dumps(skip_manifest, indent=2, sort_keys=True) + "\n",
                encoding="utf-8")
            result.update(
                status=SKIPPED,
                reason=("REGRESSION mode cannot freeze or score I13 without "
                        "passing FULL I11 evidence"),
                metrics={"linked_model_executed": False,
                         "mpi_ranks_requested": mpi_ranks},
                artifacts={"input_snapshot": str(output / "application.in"),
                           "run_manifest": str(output / "run_manifest.json"),
                           "stdout": str(stdout_path),
                           "stderr": str(stderr_path)})
            raise StopIteration
        if test_id in ("I12", "I13") and application_input is not None:
            if moon_setting(application_input.resolve(),
                            "mesh_output_tecplot") != "false":
                raise MissingPrerequisite(
                    f"{test_id} requires mesh_output_tecplot=false; optional "
                    "full-volume visualization is outside its test contract")
        if test_id == "I13" and application_input is not None:
            if mesh_profile(application_input.resolve()) != "full":
                raise MissingPrerequisite(
                    "I13 requires mesh_profile=full because its 5 km case "
                    "must be derived from the exact I11 FULL template")
            derive_i11_input(application_input.resolve(), effective_input, 5)  # type: ignore[arg-type]
        run_manifest = manifest(test_id, executable, command, seed,
                                effective_input, mpi_ranks)
        (output / "run_manifest.json").write_text(json.dumps(run_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        if application_input is not None and test_id != "I13":
            # Preserve the exact human-readable receipt beside its hash. This
            # is a run artifact, never a replacement for the reviewed source
            # template and never an automatically updated reference solution.
            shutil.copy2(application_input, output / "application.in")
        if test_id == "I01":
            # R0.1 requires both reproducibility and sensitivity.  Recreate the
            # same manifest from the same controls, then change only the thread
            # control.  Timestamps are intentionally absent, so equality is an
            # exact deterministic comparison rather than a filtered one.
            repeated = manifest(test_id, executable, command, seed,
                                application_input, mpi_ranks)
            changed = json.loads(json.dumps(repeated))
            changed["openmp_threads"] = 2
            (output / "run_manifest_repeat.json").write_text(json.dumps(repeated, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            (output / "run_manifest_changed_control.json").write_text(json.dumps(changed, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            metrics = {"identical_control_equal": repeated == run_manifest,
                       "changed_control_detected": changed != repeated,
                       "hashed_external_files": len(run_manifest["external_data"])}
            status = PASS if (metrics["identical_control_equal"] and
                              metrics["changed_control_detected"] and
                              metrics["hashed_external_files"] == len(SPICE_FILES)) else FAIL
            reason = "manifest repeatability and changed-control sensitivity"
        elif test_id == "I04":
            status, metrics, reason = check(test_id, {}, output)
        else:
            environment = dict(os.environ); environment["OMP_NUM_THREADS"] = "1"
            # AMPS caches an approximately 650 MB cell-measure file in its
            # working directory.  Tests are serial and share identical mesh
            # geometry within a configuration, so use one campaign-owned cache
            # instead of regenerating/copying it for every I ID.  Receipts and
            # logs still remain in their per-test artifact directories.
            runtime_label = "ion-nospice" if test_id == "I08" else (
                "neutral-nospice" if test_id == "I05" else "neutral-spice")
            runtime_directory = output_root / "_runtime_cache" / runtime_label
            runtime_directory.mkdir(parents=True, exist_ok=True)
            # Retain the same one-second process-tree RSS and host-availability
            # history for I12/I13 that I11 records. This is observation only:
            # it imposes no timeout or memory threshold. Actual C++ allocation
            # failure remains responsible for the detailed AMPS fatal exit.
            returncode = run_monitored(
                command, runtime_directory, environment, stdout_path,
                stderr_path, output / "process_resource_metrics.csv")
            if returncode != 0 or not receipt_path.is_file():
                raise RuntimeError(f"linked executable exited {returncode}; see {stdout_path} and {stderr_path}")
            receipt = parse_receipt(receipt_path)
            if receipt.get("test_id") != test_id:
                raise ValueError(f"receipt ID {receipt.get('test_id')!r} does not match {test_id}")
            status, metrics, reason = check(test_id, receipt, output)
        result.update(status=status, metrics=metrics, reason=reason)
    except StopIteration:
        pass
    except MissingPrerequisite as exception:
        result.update(status=SKIPPED, reason=str(exception))
    except (OSError, ValueError, RuntimeError, KeyError, ZeroDivisionError) as exception:
        result.update(status=ERROR, reason=str(exception))
    (output / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"{result['status']} {test_id}: {result['reason']}")
    return {PASS: 0, FAIL: 1, ERROR: 2, SKIPPED: 77}[str(result["status"])]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    # Individual `IXX/test.py` wrappers intentionally need no redundant ID on
    # their command line.  Direct calls to this shared module may still provide
    # one explicitly, and the frozen parent-directory name supplies the default.
    inferred_id = Path(sys.argv[0]).resolve().parent.name
    parser.add_argument("test_id", nargs="?", default=inferred_id,
                        choices=[f"I{number:02d}" for number in range(1, 14)])
    parser.add_argument("--amps", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--input", type=Path,
                        help="validated srcMoon .in receipt for I11-I13")
    parser.add_argument(
        "--mpi-ranks", type=int, default=1,
        help=("positive MPI process count; the exact value is used by mpiexec "
              "and recorded in run_manifest.json (default: 1)"))
    parser.add_argument(
        "--campaign-mode", choices=("regression", "qualification"),
        default="qualification",
        help=("REGRESSION permits only the bounded linked wiring path; "
              "QUALIFICATION runs the formal M1 gate (default: qualification)"))
    arguments = parser.parse_args()
    if arguments.mpi_ranks < 1:
        parser.error("--mpi-ranks must be a positive integer")
    return run(arguments.test_id, arguments.amps, arguments.output_dir.resolve(),
               arguments.input, arguments.mpi_ranks, arguments.campaign_mode)


if __name__ == "__main__":
    raise SystemExit(main())
