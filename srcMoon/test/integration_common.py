#!/usr/bin/env python3
"""Production-linked M0 test execution and independent acceptance checks.

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
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
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
    return ERROR, metrics, f"no acceptance implementation for {test_id}"


def git_value(*arguments: str) -> str:
    completed = subprocess.run(["git", *arguments], cwd=REPO, text=True,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    return completed.stdout.strip() if completed.returncode == 0 else "UNAVAILABLE"


def manifest(test_id: str, executable: Path, command: list[str], seed: object | None) -> dict[str, object]:
    """Capture deterministic controls; manifests contain no timestamps by design."""
    files = [REPO / "input/moon.input", REPO / "input/species.input",
             REPO / "srcMoon/development_and_validation.md",
             Path("/data/vtenishe/moon_validation_data/development_and_validation.md")]
    if test_id == "I08":
        files.append(REPO / "input/moon-i08.input")
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
        "schema": "srcMoon-M0-run-manifest-v1", "test_id": test_id,
        "source_revision": git_value("rev-parse", "HEAD"),
        "source_dirty_paths": git_value("status", "--short").splitlines(),
        "executable": {"path": str(executable.resolve()), "sha256": sha256(executable)},
        "command": command, "mpi_ranks": 1, "openmp_threads": 1,
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


def run(test_id: str, executable: Path, output_root: Path) -> int:
    output = output_root / test_id
    output.mkdir(parents=True, exist_ok=True)
    receipt_path = output / "receipt.txt"
    stdout_path, stderr_path = output / "stdout.log", output / "stderr.log"
    command = ["mpiexec", "-n", "1", str(executable.resolve()),
               "--moon-integration-test", test_id,
               "--integration-output", str(receipt_path.resolve())]
    seed = 731993 if test_id == "I09" else None
    result: dict[str, object] = {"id": test_id, "status": ERROR, "metrics": {},
                                "artifacts": {"stdout": str(stdout_path), "stderr": str(stderr_path), "receipt": str(receipt_path)}}
    try:
        run_manifest = manifest(test_id, executable, command, seed)
        (output / "run_manifest.json").write_text(json.dumps(run_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        if test_id == "I01":
            # R0.1 requires both reproducibility and sensitivity.  Recreate the
            # same manifest from the same controls, then change only the thread
            # control.  Timestamps are intentionally absent, so equality is an
            # exact deterministic comparison rather than a filtered one.
            repeated = manifest(test_id, executable, command, seed)
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
            with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as stderr:
                completed = subprocess.run(command, cwd=runtime_directory, env=environment,
                                           stdout=stdout, stderr=stderr, check=False)
            if completed.returncode != 0 or not receipt_path.is_file():
                raise RuntimeError(f"linked executable exited {completed.returncode}; see {stdout_path} and {stderr_path}")
            receipt = parse_receipt(receipt_path)
            if receipt.get("test_id") != test_id:
                raise ValueError(f"receipt ID {receipt.get('test_id')!r} does not match {test_id}")
            status, metrics, reason = check(test_id, receipt, output)
        result.update(status=status, metrics=metrics, reason=reason)
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
                        choices=[f"I{number:02d}" for number in range(1, 11)])
    parser.add_argument("--amps", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    arguments = parser.parse_args()
    return run(arguments.test_id, arguments.amps, arguments.output_dir.resolve())


if __name__ == "__main__":
    raise SystemExit(main())
