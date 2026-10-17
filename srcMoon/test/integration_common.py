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
import struct
import subprocess
import sys
import time
from typing import Callable

REPO = Path(__file__).resolve().parents[2]
TEST_ROOT = Path(__file__).resolve().parent
SPICE_ROOT = Path("/home/vtenishe/SPICE/Kernels")
DEFAULT_VALIDATION_DATA_ROOT = Path("/data/vtenishe/moon_validation_data")
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
I11_SAMPLE_LENGTH_STEPS = 500
I11_MAXIMUM_PARTICLE_STEPS = 5000
I11_MAXIMUM_COMPLETE_SAMPLES = (
    I11_MAXIMUM_PARTICLE_STEPS // I11_SAMPLE_LENGTH_STEPS
)


class MissingPrerequisite(RuntimeError):
    """A declared campaign prerequisite/configuration is unavailable."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validation_data_root() -> Path:
    """Return the explicit external-data root used by this campaign.

    The canonical path comes from the srcMoon validation contract.  The
    environment override is retained for a reviewed replica, but the resolved
    path is always written to the run evidence; silently searching another
    directory would make the scientific input irreproducible.
    """
    return Path(os.environ.get(
        "MOON_VALIDATION_DATA_ROOT",
        str(DEFAULT_VALIDATION_DATA_ROOT))).resolve()


def dataset_package_preflight(
    dataset_id: str,
) -> tuple[bool, str, dict[str, object]]:
    """Verify package-level readiness before reading processed observations.

    Raw-file presence is deliberately insufficient.  AGENTS.md requires the
    package README, provenance record, and QA report to exist and to agree on
    an explicit READY state.  Missing or non-READY metadata returns a normal
    scientific prerequisite miss (SKIPPED); malformed JSON or contradictory
    status records raise ``ValueError`` and are therefore classified ERROR as
    data/infrastructure corruption by the caller.
    """
    package_name = {
        "D01": "lola", "D02": "lola", "D03": "diviner",
    }.get(dataset_id)
    if package_name is None:
        raise ValueError(f"no package mapping is declared for {dataset_id}")

    package = validation_data_root() / package_name
    required = {
        "readme": package / "README.md",
        "provenance": package / "provenance.json",
        "qa_report": package / "qa_report.json",
    }
    evidence: dict[str, object] = {
        "dataset_id": dataset_id,
        "package": str(package),
        "files": {},
    }
    missing = [str(path) for path in required.values() if not path.is_file()]
    evidence_files = evidence["files"]
    if not isinstance(evidence_files, dict):
        raise AssertionError("dataset evidence file table was not initialized")
    for label, path in required.items():
        if path.is_file():
            evidence_files[label] = {
                "path": str(path), "sha256": sha256(path),
            }
    if missing:
        return False, (
            f"{dataset_id} package metadata is incomplete; missing " +
            ", ".join(missing)), evidence

    # Read the README as a required human contract even though readiness is
    # taken only from the two machine-readable records below.
    required["readme"].read_text(encoding="utf-8")
    provenance = json.loads(required["provenance"].read_text(encoding="utf-8"))
    qa_report = json.loads(required["qa_report"].read_text(encoding="utf-8"))
    provenance_status = (
        provenance.get("datasets", {}).get(dataset_id, {}).get("status"))
    qa_status = qa_report.get("dataset_status", {}).get(dataset_id)
    evidence["provenance_status"] = provenance_status
    evidence["qa_status"] = qa_status
    if provenance_status is None or qa_status is None:
        return False, (
            f"{dataset_id} package does not declare its status in both "
            "provenance.json and qa_report.json"), evidence
    if provenance_status != qa_status:
        raise ValueError(
            f"{dataset_id} package status conflict: provenance="
            f"{provenance_status!r}, QA={qa_status!r}")
    if provenance_status != "READY":
        return False, (
            f"{dataset_id} package status is {provenance_status}, not READY"), \
            evidence
    return True, f"{dataset_id} package is READY", evidence


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


def ldem4_elevation_m(
    image: bytes, latitude_deg: float, longitude_east_deg: float
) -> float:
    """Independently decode/interpolate the frozen LDEM_4 native raster.

    The dimensions, edge coordinates, little-endian signed sample type, and
    0.5 m/DN scale are declared by LDEM_4.LBL and frozen by U08.  This Python
    oracle intentionally does not call ``Moon::Surface::LolaDem``.
    """
    lines, samples = 720, 1440
    longitude = longitude_east_deg % 360.0
    latitude = max(-90.0, min(90.0, latitude_deg))
    column = longitude / (360.0 / samples) - 0.5
    row = (90.0 - latitude) / (180.0 / lines) - 0.5
    column_floor = math.floor(column)
    column_weight = column - column_floor
    row_clamped = max(0.0, min(lines - 1.0, row))
    row0 = math.floor(row_clamped)
    row1 = min(row0 + 1, lines - 1)
    row_weight = row_clamped - row0

    def value(row_index: int, column_index: int) -> float:
        wrapped = column_index % samples
        offset = 2 * (row_index * samples + wrapped)
        return 0.5 * struct.unpack_from("<h", image, offset)[0]

    north = ((1.0 - column_weight) * value(row0, column_floor) +
             column_weight * value(row0, column_floor + 1))
    south = ((1.0 - column_weight) * value(row1, column_floor) +
             column_weight * value(row1, column_floor + 1))
    return (1.0 - row_weight) * north + row_weight * south


def render_comparison(work: Path, title: str) -> dict[str, str]:
    """Render an already-written row-level comparison without changing it."""
    comparison = work / "comparison.csv"
    png = work / "comparison.png"
    eps = work / "comparison.eps"
    command = [
        sys.executable,
        str(TEST_ROOT / "common/write_comparison_artifacts.py"),
        "--input", str(comparison), "--png", str(png), "--eps", str(eps),
        "--title", title,
    ]
    completed = subprocess.run(
        command, cwd=REPO, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            "comparison plot renderer failed: " + completed.stderr.strip())
    return {"comparison_csv": str(comparison), "comparison_png": str(png),
            "comparison_eps": str(eps)}


def render_i11_stage_diagnostic(
    campaign: Path, completed_resolutions_km: list[int], formal_status: str,
) -> dict[str, object]:
    """Render raw diagnostics after a completed I11 resolution stage.

    The dedicated renderer reads production ``block_metrics.csv`` files and
    explicitly marks every artifact as unqualified.  It is deliberately
    separate from ``qualified_metrics.csv`` and ``convergence_metrics.csv``:
    producing a meeting/demo figure cannot affect the formal I11 decision.
    """
    command = [
        sys.executable, "-B",
        str(TEST_ROOT / "integration-tests/I11/write_stage_diagnostics.py"),
        "--campaign-dir", str(campaign),
        "--completed-resolutions-km",
        *[str(value) for value in completed_resolutions_km],
        "--formal-level-status", formal_status,
    ]
    completed = subprocess.run(
        command, cwd=REPO, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            "I11 stage diagnostic renderer failed: " +
            (completed.stderr.strip() or completed.stdout.strip()))
    try:
        metadata = json.loads(completed.stdout)
    except json.JSONDecodeError as exception:
        raise RuntimeError(
            "I11 stage diagnostic renderer returned malformed metadata") from exception
    if metadata.get("accepted_convergence_input") is not False:
        raise RuntimeError("I11 diagnostic artifact was not marked unqualified")
    return metadata


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


def moon_setting(path: Path, key: str, preserve_case: bool = False) -> str:
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
    # Enum comparisons are case insensitive, but filesystem paths are not.
    # In particular the staged D03-RDR-I17 directory must retain its case.
    return matches[0] if preserve_case else matches[0].lower()


def mesh_profile(path: Path) -> str:
    """Return the explicit QUICK/FULL/CUSTOM provenance label."""
    return moon_setting(path, "mesh_profile")


def replace_moon_setting(text: str, key: str, value: str, source: Path) -> str:
    """Replace exactly one flat Moon assignment in a derived run receipt."""
    pattern = rf"(?m)^({re.escape(key)}\s*=\s*)([^!\n]+?)(\s*)$"
    matches = list(re.finditer(pattern, text))
    if len(matches) != 1:
        raise ValueError(f"expected exactly one {key} in {source}")
    return re.sub(pattern, lambda match: match.group(1) + value +
                  match.group(3), text)


def prepare_na_ar_input(template: Path, destination: Path) -> None:
    """Derive the runtime half of the reviewed two-species M2 receipt.

    The compile-time half is produced by ``build_m2_qualification.py``.  This
    transformer enables the explicit radiogenic-Ar switch and coarsens only
    the *volume* AMR ownership mesh to the bounded regression targets.  I18
    evaluates surface-impact and source distributions, not transport
    convergence: retaining the FULL 5--10 km transport mesh would allocate
    hundreds of GiB without changing any sampled probability law.  The LOLA
    surface remains at the approved 100 km resolution, and source rate,
    geography, time dependence, thermal state, surface parameters, seed, and
    every acceptance threshold remain unchanged.  ``mesh_profile=custom``
    prevents this bounded fixture from being mistaken for M1/M2 spatial-
    convergence evidence.
    """
    text = template.read_text(encoding="utf-8")
    text = replace_moon_setting(
        text, "argon_source_enabled", "true", template)
    bounded_volume_settings = {
        "mesh_profile": "custom",
        "mesh_near_surface_cell_size_m": "120000",
        "mesh_far_field_cell_size_m": "480000",
        "mesh_outer_cell_size_m": "208452000",
        "mesh_angular_resolution_m_by_sza":
            "0:120000,30:120000,60:60000,90:60000,120:120000,180:120000",
        "mesh_polar_resolution_m": "60000",
    }
    for key, value in bounded_volume_settings.items():
        text = replace_moon_setting(text, key, value, template)
    # Surface serialization is part of normal LOLA initialization.  Point the
    # two products into this test's immutable artifact directory so a focused
    # I18 run cannot overwrite another campaign's scientifically meaningful
    # mesh files.  These paths affect output ownership only, not geometry.
    text = replace_moon_setting(
        text, "surface_cea_file",
        str((destination.parent / "lola_surface.cea").resolve()), template)
    text = replace_moon_setting(
        text, "surface_tecplot_file",
        str((destination.parent / "lola_surface.dat").resolve()), template)
    text = (
        "! Generated I18 Na+Ar receipt: the 100 km LOLA surface and all "
        "physics are unchanged; only the non-scoring volume ownership mesh "
        "uses bounded regression targets.\n" + text
    )
    destination.write_text(text, encoding="utf-8")


def parse_temperature_probability_table(path: Path, key: str) -> list[tuple[float, float]]:
    """Parse one explicitly configured ``temperature:probability`` table."""
    raw = moon_setting(path, key)
    rows: list[tuple[float, float]] = []
    for field in raw.split(","):
        parts = field.split(":")
        if len(parts) != 2:
            raise ValueError(f"invalid {key} entry {field!r}")
        temperature, probability = map(float, parts)
        if (not math.isfinite(temperature) or temperature <= 0.0 or
                not math.isfinite(probability) or not 0.0 <= probability <= 1.0):
            raise ValueError(f"invalid {key} entry {field!r}")
        rows.append((temperature, probability))
    if not rows or any(rows[index][0] <= rows[index - 1][0]
                       for index in range(1, len(rows))):
        raise ValueError(f"{key} temperatures must be strictly increasing")
    return rows


def interpolate_bounded_table(rows: list[tuple[float, float]], value: float) -> float:
    """Independent linear interpolation with constant endpoint extension."""
    if value <= rows[0][0]:
        return rows[0][1]
    if value >= rows[-1][0]:
        return rows[-1][1]
    for index in range(1, len(rows)):
        if value <= rows[index][0]:
            x0, y0 = rows[index - 1]
            x1, y1 = rows[index]
            fraction = (value - x0) / (x1 - x0)
            return (1.0 - fraction) * y0 + fraction * y1
    raise AssertionError("bounded interpolation did not find an interval")


def one_sample_ks(values: list[float], cdf: Callable[[float], float]) -> float:
    """Return the exact one-sample empirical-CDF Kolmogorov distance."""
    ordered = sorted(values)
    if not ordered:
        raise ValueError("a KS comparison requires at least one sample")
    count = len(ordered)
    distance = 0.0
    for index, value in enumerate(ordered, 1):
        expected = min(1.0, max(0.0, cdf(value)))
        distance = max(distance, abs(index / count - expected),
                       abs((index - 1) / count - expected))
    return distance


def prepare_i16_input(template: Path, destination: Path,
                      output_root: Path) -> Path:
    """Process D03 and derive the auditable DIVINER_MAP linked receipt."""
    data_directory = output_root / "_m2_validation_data" / "D03-GCP"
    table = data_directory / "temperature_table.csv"
    benchmark = data_directory / "benchmark_points.csv"
    qa_path = data_directory / "d03_gcp_qa.json"
    if not (table.is_file() and benchmark.is_file() and qa_path.is_file()):
        raw = validation_data_root() / "diviner/raw/gcp"
        if len(list(raw.glob("global_cumul_avg_cyl_*_002.tab"))) != 18:
            raise MissingPrerequisite(
                "I16 requires all 18 authoritative D03 GCP tables")
        command = [sys.executable, "-B", str(TEST_ROOT /
            "data_processing/prepare_d03_gcp.py"), "--output-dir",
            str(data_directory)]
        completed = subprocess.run(
            command, cwd=REPO, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False)
        (data_directory / "processor.stdout.log").write_text(
            completed.stdout, encoding="utf-8")
        (data_directory / "processor.stderr.log").write_text(
            completed.stderr, encoding="utf-8")
        if completed.returncode != 0:
            raise RuntimeError(
                "D03 GCP processor failed; see " +
                str(data_directory / "processor.stderr.log"))
    qa = json.loads(qa_path.read_text(encoding="utf-8"))
    if qa.get("status") != "READY":
        raise MissingPrerequisite("processed D03 GCP package is not READY")
    if (qa.get("processed", {}).get("temperature_table", {}).get("sha256") !=
            sha256(table) or
            qa.get("processed", {}).get("benchmark_points", {}).get("sha256") !=
            sha256(benchmark)):
        raise ValueError("processed D03 GCP hash does not match its QA record")
    text = template.read_text(encoding="utf-8")
    text = replace_moon_setting(
        text, "surface_temperature_mode", "diviner-map", template)
    text = replace_moon_setting(
        text, "diviner_temperature_table_file", str(table.resolve()), template)
    destination.write_text(text, encoding="utf-8")
    return qa_path


def prepare_i17_input(template: Path, destination: Path,
                      output_root: Path) -> Path:
    """Stage a qualified frozen RDR product and derive an auditable receipt.

    Processed validation products belong below the canonical dataset's
    ``processed/`` directory, while linked-run artifacts belong below the
    campaign output.  The runner therefore never performs an implicit 11-GB
    raw-data transformation.  It requires the separately generated canonical
    product, verifies its QA hashes, and copies the exact evidence into this
    run before AMPS sees it.  The processor has no AMPS-output argument and one
    pre-holdout contract is applied to both dates, so December residuals cannot
    alter membership, native-flag cuts, binning, channel rules, or uncertainty.
    """
    canonical_directory = validation_data_root() / "diviner/processed/i17"
    canonical_history = canonical_directory / "eclipse_thermal_history.csv"
    canonical_qa = canonical_directory / "d03_rdr_qa.json"
    canonical_contract = canonical_directory / "selection_contract.json"
    contract = (TEST_ROOT /
        "integration-tests/I17/reference/selection_contract.json").resolve()
    if not (canonical_history.is_file() and canonical_qa.is_file() and
            canonical_contract.is_file()):
        raise MissingPrerequisite(
            "I17 requires the committed processor to be run separately into "
            f"{canonical_directory}; linked tests never create or edit the "
            "canonical validation-data package implicitly")
    qa = json.loads(canonical_qa.read_text(encoding="utf-8"))
    if qa.get("status") != "READY":
        raise MissingPrerequisite("processed D03 RDR eclipse product is not READY")
    processor = (TEST_ROOT /
        "data_processing/prepare_d03_rdr_eclipse.py").resolve()
    if (qa.get("processed", {}).get("sha256") != sha256(canonical_history) or
            qa.get("contract", {}).get("sha256") !=
            sha256(canonical_contract) or
            sha256(canonical_contract) != sha256(contract) or
            qa.get("processor", {}).get("sha256") != sha256(processor)):
        raise ValueError(
            "D03 RDR history/contract/processor hash mismatch")

    # Preserve the exact observational product and QA receipt inside the
    # linked campaign.  AMPS reads this immutable staged copy, so later edits
    # to a shared processed directory cannot change an already recorded run.
    data_directory = output_root / "_m2_validation_data" / "D03-RDR-I17"
    data_directory.mkdir(parents=True, exist_ok=True)
    history = data_directory / canonical_history.name
    training_history = data_directory / "training_eclipse_thermal_history.csv"
    qa_path = data_directory / canonical_qa.name
    contract_snapshot = data_directory / contract.name

    def stage(source: Path, target: Path) -> None:
        # A campaign directory may be inspected/restarted, but scientifically
        # meaningful evidence is never overwritten with different bytes.
        if target.is_file():
            if sha256(target) != sha256(source):
                raise ValueError(
                    f"refusing to overwrite changed I17 artifact {target}")
            return
        shutil.copy2(source, target)

    stage(canonical_history, history)
    stage(canonical_qa, qa_path)
    stage(canonical_contract, contract_snapshot)
    for name in ("processor.stdout.log", "processor.stderr.log"):
        source = canonical_directory / name
        if source.is_file():
            stage(source, data_directory / name)

    # The first-comparison campaign is allowed to inspect June training rows
    # only.  Materialize a byte-stable training view before launching AMPS so
    # neither C++ nor the comparison scorer can accidentally read a December
    # temperature while numerical thresholds are still unfrozen.  Filtering
    # uses only the already-frozen ``selection`` label; it never examines a
    # model residual and never edits the canonical processed product.
    temporary = training_history.with_suffix(".csv.tmp")
    with history.open(newline="", encoding="utf-8") as source, \
            temporary.open("w", newline="", encoding="utf-8") as target:
        reader = csv.DictReader(source)
        if reader.fieldnames is None or "selection" not in reader.fieldnames:
            raise ValueError("I17 processed history has no selection field")
        writer = csv.DictWriter(target, fieldnames=reader.fieldnames,
                                lineterminator="\n")
        writer.writeheader()
        count = 0
        for row in reader:
            if row["selection"] == "training":
                writer.writerow(row)
                count += 1
            elif row["selection"] != "holdout":
                raise ValueError(
                    f"I17 history has unknown selection {row['selection']!r}")
    if count == 0:
        temporary.unlink(missing_ok=True)
        raise ValueError("I17 processed history has no June training rows")
    if training_history.is_file():
        if sha256(training_history) != sha256(temporary):
            temporary.unlink(missing_ok=True)
            raise ValueError(
                f"refusing to overwrite changed I17 artifact {training_history}")
        temporary.unlink()
    else:
        temporary.replace(training_history)
    text = template.read_text(encoding="utf-8")
    # This receipt moves to a campaign directory. Resolve source input paths
    # against the original template and give generated outputs unique owners;
    # otherwise relative paths silently change meaning during staging.
    for key in ("spice_path", "lola_image_file", "lola_label_file"):
        original = Path(moon_setting(template, key, preserve_case=True))
        absolute = (original if original.is_absolute() else template.parent / original)
        text = replace_moon_setting(text, key, str(absolute.resolve()), template)
    for key, name in (("surface_cea_file", "surface.cea"),
                      ("surface_tecplot_file", "surface.dat"),
                      ("surface_inventory_checkpoint_file", "inventory.chk")):
        text = replace_moon_setting(text, key,
            str((destination.parent / name).resolve()), template)
    text = replace_moon_setting(
        text, "diviner_eclipse_history_file",
        str(training_history.resolve()), template)
    text = replace_moon_setting(
        text, "diviner_eclipse_selection_contract_file",
        str(contract_snapshot.resolve()), template)
    text = (
        "! Generated I17 receipt: RDR evidence and its frozen contract are "
        "observation inputs only and never initialize or tune the thermal state.\n" +
        text)
    destination.write_text(text, encoding="utf-8")
    return qa_path


def qualification_levels(template: Path) -> tuple[int, ...]:
    """Select an approved campaign, never promote an arbitrary CUSTOM mesh.

    Lengths below are metres. The 2026-10-14 memory campaign was approved
    before its science results were examined. Its finest input has measured
    allocation evidence, but this whitelist is NOT convergence evidence.
    FULL retains the historical campaign for explicitly provisioned hosts.
    """
    profile = mesh_profile(template)
    if profile == "full":
        return (40, 20, 10, 5)
    expected = {
        "surface_geometry": "sphere",
        "surface_mesh_resolution_m": "100000",
        "mesh_near_surface_cell_size_m": "60000",
        "mesh_near_surface_extent_m": "200000",
        "mesh_far_field_cell_size_m": "480000",
        "mesh_far_field_extent_m": "10000000",
        "mesh_outer_cell_size_m": "104226000",
        "mesh_angular_resolution_m_by_sza":
            "0:60000,30:60000,60:30000,90:30000,120:60000,180:60000",
        "mesh_terminator_half_width_deg": "15",
        "mesh_polar_cap_start_latitude_deg": "75",
        "mesh_polar_resolution_m": "30000",
        "simulation_environment_mode": "frozen",
        "mesh_output_tecplot": "false",
    }
    if profile == "custom" and all(
            moon_setting(template, key) == value for key, value in expected.items()):
        return (480, 240, 120, 60)
    raise MissingPrerequisite(
        "qualification requires the approved m1_memory200.in mesh settings "
        "or the legacy FULL profile; arbitrary CUSTOM/QUICK meshes are not qualified")


def derive_i11_input(template: Path, destination: Path,
                     near_surface_km: int) -> None:
    """Create one archived I11 receipt without editing the reviewed template.

    The selected template is the frozen physical case. Near-surface, spherical-
    boundary, angular, and polar targets are multiplied together because they
    jointly determine boundary-adjacent resolution. Far-field and outer-domain
    targets remain fixed: scaling the complete 400-lunar-radius tail would mix
    a near-surface convergence question with a separate far-domain sensitivity
    question and produces an intractable three-dimensional block count.
    Physical extents, altitude comparison points, epoch, sources, and
    observation geometry remain identical. Every expected assignment must
    occur exactly once, otherwise derivation fails closed.
    """
    # Normalize against the actual template, not a hidden 10 km assumption:
    # the memory campaign's template is 60 km and its finest factor is one.
    base_m = float(moon_setting(template, "mesh_near_surface_cell_size_m"))
    if not math.isfinite(base_m) or base_m <= 0 or near_surface_km <= 0:
        raise ValueError("I11 resolution targets must be finite and positive")
    # Mirror the production input constraint BEFORE MPI launch. In particular
    # the approved 480 km coarse level cannot use a fixed 240 km far target.
    # Do not silently clamp individual levels: that would change two spatial
    # controls at once and obscure the intended near-surface convergence test.
    if near_surface_km * 1000.0 > float(moon_setting(
            template, "mesh_far_field_cell_size_m")):
        raise ValueError("I11 near-surface target exceeds fixed far-field target")
    factor = near_surface_km * 1000.0 / base_m
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
            try:
                qualification_levels(work / "application.in")
            except MissingPrerequisite:
                return SKIPPED, metrics, (
                    f"{profile.upper()} linked wiring checks passed; mesh is "
                    "not an approved qualification profile")
        return PASS, metrics, (
            "approved input-defined angular formula and linked realized-cell "
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
        finest = qualification_levels(work / "application.in")[-1]
        qualified = work.parent / f"I11/{finest}km/qualified_metrics.csv"
        i11_result = work.parent / "I11/result.json"
        i11_finest_input = work.parent / f"I11/{finest}km/application.in"
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
                "I13 must use the exact passing I11 finest derived receipt; "
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
    if test_id == "I14":
        package = Path(os.environ.get(
            "MOON_VALIDATION_DATA_ROOT",
            "/data/vtenishe/moon_validation_data")) / "lola"
        image_path = package / "raw/LDEM_4.IMG"
        label_path = package / "raw/LDEM_4.LBL"
        qa_path = package / "qa_report.json"
        provenance_path = package / "provenance.json"
        required = [image_path, label_path, qa_path, provenance_path]
        missing = [str(path) for path in required if not path.is_file()]
        if missing:
            raise MissingPrerequisite(
                "I14 requires the qualified D01 package; missing: " +
                ", ".join(missing))
        qa = json.loads(qa_path.read_text(encoding="utf-8"))
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
        if (qa.get("dataset_status", {}).get("D01") != "READY" or
                provenance.get("datasets", {}).get("D01", {}).get("status") !=
                "READY"):
            raise MissingPrerequisite(
                "I14 requires D01 READY in both QA and provenance metadata")
        u08_contract = json.loads((TEST_ROOT /
            "U08_lola_geometry/reference/acceptance.json").read_text(
                encoding="utf-8"))
        for filename, expected in u08_contract["data_hashes"].items():
            path = package / "raw" / filename
            if sha256(path) != expected:
                raise ValueError(f"D01 SHA-256 mismatch for {path}")

        image = image_path.read_bytes()
        if len(image) != 720 * 1440 * 2:
            raise ValueError("LDEM_4 byte count violates its detached label")
        geometry_path = work / "receipt.txt.I14_geometry.csv"
        with geometry_path.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        expected_header = [
            "vertex_index", "latitude_deg", "longitude_east_deg",
            "vertex_radius_m", "ray_radius_m", "ray_face",
        ]
        if not rows or list(rows[0]) != expected_header:
            raise ValueError("I14 linked geometry table is empty or malformed")
        comparison_rows: list[dict[str, object]] = []
        errors_m: list[float] = []
        ray_errors_m: list[float] = []
        for row in rows:
            latitude = float(row["latitude_deg"])
            longitude = float(row["longitude_east_deg"])
            model_radius = float(row["vertex_radius_m"])
            ray_radius = float(row["ray_radius_m"])
            if not all(math.isfinite(value) for value in
                       (latitude, longitude, model_radius, ray_radius)):
                raise ValueError("I14 geometry contains a non-finite value")
            reference_radius = 1_737_400.0 + ldem4_elevation_m(
                image, latitude, longitude)
            errors_m.append(abs(model_radius - reference_radius))
            ray_errors_m.append(abs(ray_radius - model_radius))
            comparison_rows.append({
                "comparison_id": "I14_D01_vertex_radius",
                "x_value": row["vertex_index"],
                "x_unit": "vertex-index",
                "reference_value": f"{reference_radius:.17g}",
                "model_value": f"{model_radius:.17g}",
                "y_unit": "m",
                # The roadmap freezes max(one vertical quantization unit,
                # 50 m). LDEM_4 uses 0.5 m/DN, so 50 m controls this gate.
                "uncertainty_value": "50",
                "selection": "verification",
            })
        comparison_path = work / "comparison.csv"
        with comparison_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(
                stream, fieldnames=list(comparison_rows[0]),
                lineterminator="\n")
            writer.writeheader()
            writer.writerows(comparison_rows)
        artifacts = render_comparison(
            work, "I14 production LOLA vertex radius versus D01")
        metrics.update({
            "comparison_count": len(rows),
            "maximum_control_point_error_m": max(errors_m),
            "rms_control_point_error_m": math.sqrt(
                sum(value * value for value in errors_m) / len(errors_m)),
            "maximum_vertex_ray_difference_m": max(ray_errors_m),
            "degenerate_face_count": int(receipt["degenerate_face_count"]),
            "inverted_face_count": int(receipt["inverted_face_count"]),
            "radial_query_count": int(receipt["radial_query_count"]),
            "radial_hit_count": int(receipt["radial_hit_count"]),
            "mover_fixture_rank_count": int(
                receipt["mover_fixture_rank_count"]),
            "mover_collision_rank_count": int(
                receipt["mover_collision_rank_count"]),
            "data_hashes": {"LDEM_4.IMG": sha256(image_path),
                            "LDEM_4.LBL": sha256(label_path)},
            "artifacts": artifacts,
        })
        ok = (
            max(errors_m) <= 50.0 and
            max(ray_errors_m) <= float(
                receipt.get("surface_collision_distance_tolerance_m", 1.0e-4)) and
            int(receipt["degenerate_face_count"]) == 0 and
            int(receipt["inverted_face_count"]) == 0 and
            int(receipt["radial_hit_count"]) ==
                int(receipt["radial_query_count"]) and
            int(receipt["mover_fixture_rank_count"]) > 0 and
            int(receipt["mover_collision_rank_count"]) ==
                int(receipt["mover_fixture_rank_count"])
        )
        return (PASS if ok else FAIL), metrics, (
            "production LOLA control points, topology, radial rays, and "
            "normal mover collision path")
    if test_id == "I16":
        table_path = Path(str(receipt["diviner_temperature_table_file"]))
        benchmark_path = table_path.parent / "benchmark_points.csv"
        qa_path = table_path.parent / "d03_gcp_qa.json"
        model_path = work / "receipt.txt.I16_model.csv"
        for path in (table_path, benchmark_path, qa_path, model_path):
            if not path.is_file():
                raise ValueError(f"I16 required artifact is missing: {path}")
        qa = json.loads(qa_path.read_text(encoding="utf-8"))
        if qa.get("status") != "READY":
            raise MissingPrerequisite("I16 processed D03 GCP QA is not READY")
        if (qa["processed"]["temperature_table"]["sha256"] !=
                sha256(table_path) or
                qa["processed"]["benchmark_points"]["sha256"] !=
                sha256(benchmark_path)):
            raise ValueError("I16 D03 processed evidence hash mismatch")
        reference_by_location: dict[
            tuple[float, float], dict[int, dict[str, str]]
        ] = {}
        with benchmark_path.open(newline="", encoding="utf-8") as stream:
            for row in csv.DictReader(stream):
                if (row["quality_flag"] != "GOOD" or
                        row["selection"] != "hold-out"):
                    continue
                key = (round(float(row["latitude_deg"]), 8),
                       round(float(row["longitude_east_deg"]), 8))
                time_index = int(round(
                    float(row["local_time_hour"]) * 4.0 - 0.5)) % 96
                if time_index in reference_by_location.setdefault(key, {}):
                    raise ValueError("duplicate I16 hold-out coordinate/time bin")
                reference_by_location[key][time_index] = row
        with model_path.open(newline="", encoding="utf-8") as stream:
            model_rows = list(csv.DictReader(stream))
        if len(model_rows) != int(receipt["benchmark_query_count"]):
            raise ValueError("I16 model row count disagrees with linked receipt")
        comparisons: list[dict[str, object]] = []
        residuals: list[float] = []
        for index, model in enumerate(model_rows):
            key = (round(float(model["latitude_deg"]), 8),
                   round(float(model["longitude_east_deg"]), 8))
            time_rows = reference_by_location.get(key, {})
            query_time = float(model["query_local_time_hour"])
            phase = (query_time * 4.0 - 0.5) % 96.0
            lower_index = math.floor(phase)
            upper_index = (lower_index + 1) % 96
            fraction = phase - lower_index
            lower = time_rows.get(lower_index)
            upper = time_rows.get(upper_index)
            if lower is None or (fraction > 1.0e-12 and upper is None):
                # A missing native hold-out bin is a real archive gap.  It is
                # not crossed by interpolation and is excluded independently
                # of the model value.
                continue
            if fraction <= 1.0e-12:
                observed = float(lower["temperature_k"])
            else:
                observed = ((1.0 - fraction) * float(lower["temperature_k"]) +
                            fraction * float(upper["temperature_k"]))  # type: ignore[index]
            predicted = float(model["production_temperature_k"])
            if not (math.isfinite(observed) and observed > 0.0 and
                    math.isfinite(predicted) and predicted > 0.0):
                raise ValueError("I16 comparison contains invalid temperature")
            residuals.append(predicted - observed)
            comparisons.append({
                "comparison_id": "I16_Diviner_GCP_holdout",
                "x_value": len(comparisons),
                "x_unit": "sample-index",
                "reference_value": f"{observed:.17g}",
                "model_value": f"{predicted:.17g}",
                "y_unit": "K",
                # The archive supplies no per-bin uncertainty.  Empty is
                # preserved through CSV and plotting rather than replaced by
                # an assumed zero or a fit-derived value.
                "uncertainty_value": "",
                "selection": "hold-out",
            })
        if not comparisons:
            raise MissingPrerequisite(
                "I16 found no hold-out local-time bin at the linked epoch")
        comparison_path = work / "comparison.csv"
        with comparison_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(
                stream, fieldnames=list(comparisons[0]), lineterminator="\n")
            writer.writeheader(); writer.writerows(comparisons)
        artifacts = render_comparison(
            work, "I16 Diviner GCP hold-out temperature")
        rmse = math.sqrt(sum(value * value for value in residuals) /
                         len(residuals))
        mean_bias = sum(residuals) / len(residuals)
        metrics.update({
            "comparison_count": len(comparisons),
            "rmse_k": rmse,
            "mean_bias_k": mean_bias,
            "maximum_absolute_residual_k": max(map(abs, residuals)),
            "temperature_table_sha256": sha256(table_path),
            "benchmark_points_sha256": sha256(benchmark_path),
            "uncertainty_metadata": "unavailable in native GCP",
            "artifacts": artifacts,
        })
        return (PASS if rmse <= 15.0 else FAIL), metrics, (
            "production surface callback versus frozen GOOD D03 GCP hold-out; "
            "RMSE gate is 15 K")
    if test_id == "I17":
        # The approved first comparison is June training only.  C++ refuses a
        # non-training row, and this scorer consumes only the model table
        # emitted after observation-free periodic spin-up.  It deliberately
        # creates a threshold *candidate*, never a PASS and never a frozen
        # hold-out criterion.
        model_path = work / "receipt.txt.I17_training_model.csv"
        if not model_path.is_file():
            raise ValueError(f"I17 linked model table is missing: {model_path}")
        with model_path.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        expected = int(receipt["training_comparison_row_count"])
        if not rows or len(rows) != expected:
            raise ValueError(
                "I17 model row count disagrees with the linked receipt")
        if int(receipt.get("holdout_rows_accessed", -1)) != 0:
            raise ValueError("I17 first comparison accessed hold-out rows")
        if receipt.get("capability_status") != "EXPERIMENTAL / NOT VALIDATED":
            raise ValueError("I17 receipt lacks its experimental status label")

        comparisons: list[dict[str, object]] = []
        diagnostics: list[dict[str, object]] = []
        residuals: list[float] = []
        normalized: list[float] = []
        first_epoch = float(rows[0]["ephemeris_time_tdb_s"])
        for row in rows:
            observed = float(
                row["observed_bolometric_brightness_temperature_K"])
            predicted = float(
                row["model_bolometric_brightness_temperature_K"])
            uncertainty = float(row["observed_uncertainty_K"])
            if not (math.isfinite(observed) and observed > 0.0 and
                    math.isfinite(predicted) and predicted > 0.0 and
                    math.isfinite(uncertainty) and uncertainty > 0.0):
                raise ValueError("I17 comparison contains invalid temperature")
            residual = predicted - observed
            residuals.append(residual)
            normalized.append(residual / uncertainty)
            comparisons.append({
                "comparison_id": "I17_Diviner_June_training_Tbol",
                "x_value": (float(row["ephemeris_time_tdb_s"]) -
                            first_epoch) / 60.0,
                "x_unit": "minutes-from-first-training-sample",
                "reference_value": f"{observed:.17g}",
                "model_value": f"{predicted:.17g}",
                "y_unit": "K",
                "uncertainty_value": f"{uncertainty:.17g}",
                "selection": "training-experimental",
            })
            for channel in range(6, 10):
                observed_radiance = float(
                    row[f"observed_c{channel}_radiance_W_m-2_sr-1"])
                model_radiance = float(
                    row[f"model_c{channel}_radiance_W_m-2_sr-1"])
                observed_brightness = float(
                    row[f"observed_c{channel}_brightness_temperature_K"])
                model_brightness = float(
                    row[f"model_c{channel}_brightness_temperature_K"])
                if not all(math.isfinite(value) for value in
                           (observed_radiance, model_radiance,
                            observed_brightness, model_brightness)):
                    raise ValueError("I17 C6-C9 diagnostic is non-finite")
                diagnostics.append({
                    "site_time_bin_id": row["site_time_bin_id"],
                    "utc": row["utc"],
                    "channel": channel,
                    "observed_radiance_W_m-2_sr-1": observed_radiance,
                    "model_radiance_W_m-2_sr-1": model_radiance,
                    "radiance_residual_W_m-2_sr-1": (
                        model_radiance - observed_radiance),
                    "observed_brightness_temperature_K": observed_brightness,
                    "model_brightness_temperature_K": model_brightness,
                    "brightness_temperature_residual_K": (
                        model_brightness - observed_brightness),
                    "status": "diagnostic-only",
                })
        comparison_path = work / "comparison.csv"
        with comparison_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(
                stream, fieldnames=list(comparisons[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(comparisons)
        diagnostic_path = work / "channel_diagnostics.csv"
        with diagnostic_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(
                stream, fieldnames=list(diagnostics[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(diagnostics)
        artifacts = render_comparison(
            work,
            "I17 EXPERIMENTAL June training: model vs Diviner T_bol")
        rmse = math.sqrt(sum(value * value for value in residuals) /
                         len(residuals))
        bias = sum(residuals) / len(residuals)
        weight_sum = sum(1.0 / float(row["observed_uncertainty_K"]) ** 2
                         for row in rows)
        weighted_bias = sum(
            residual / float(row["observed_uncertainty_K"]) ** 2
            for residual, row in zip(residuals, rows)) / weight_sum
        weighted_rmse = math.sqrt(sum(
            residual * residual /
            float(row["observed_uncertainty_K"]) ** 2
            for residual, row in zip(residuals, rows)) / weight_sum)
        candidate = {
            "schema": "srcMoon-I17-June-threshold-candidate-v1",
            "status": "REVIEW_REQUIRED_NOT_FROZEN",
            "selection": "training",
            "comparison_count": len(rows),
            "rmse_k": rmse,
            "mean_bias_k": bias,
            "weighted_rmse_k": weighted_rmse,
            "weighted_bias_k": weighted_bias,
            "normalized_residual_rms": math.sqrt(sum(
                value * value for value in normalized) / len(normalized)),
            "maximum_absolute_residual_k": max(map(abs, residuals)),
            "warning": (
                "These June diagnostics are not acceptance thresholds. "
                "Investigator review and a separately committed frozen "
                "artifact are required before December rows may be scored."),
            "assumptions": {
                "solar_disc": receipt["solar_disc_radiance_model"],
                "earth_atmosphere": receipt[
                    "earth_atmosphere_transmission_model"],
                "capability": receipt["capability_status"],
            },
        }
        candidate_path = work / "june_threshold_candidate.json"
        candidate_path.write_text(
            json.dumps(candidate, indent=2, sort_keys=True) + "\n",
            encoding="utf-8")
        metrics.update({
            **candidate,
            "spinup_completed_cycles": int(receipt["spinup_completed_cycles"]),
            "spinup_maximum_node_temperature_change_k": float(receipt[
                "spinup_maximum_node_temperature_change_k"]),
            "spinup_maximum_column_enthalpy_change_j_m2": float(receipt[
                "spinup_maximum_column_enthalpy_change_j_m2"]),
            "maximum_production_update_temperature_change_k": float(receipt[
                "maximum_production_update_temperature_change_k"]),
            "channel_diagnostics": str(diagnostic_path),
            "threshold_candidate": str(candidate_path),
            "artifacts": artifacts,
        })
        return SKIPPED, metrics, (
            "June training comparison and EXPERIMENTAL plots were produced; "
            "the candidate threshold is not frozen and December hold-out "
            "rows remain inaccessible, so I17 is NOT VALIDATED")
    if test_id == "I18":
        # These row-level samples are emitted by the ordinary production impact
        # and radiogenic-source kernels on every MPI rank.  Python supplies the
        # independent analytic distribution oracle and never feeds values back
        # into the model.
        sample_paths = sorted(work.glob("receipt.txt.I18.rank*.csv"))
        expected_ranks = int(receipt["mpi_ranks"])
        if len(sample_paths) != expected_ranks:
            raise ValueError(
                f"I18 retained {len(sample_paths)} rank tables; expected "
                f"{expected_ranks}")
        samples: list[dict[str, str]] = []
        for path in sample_paths:
            with path.open(newline="", encoding="utf-8") as stream:
                rows = list(csv.DictReader(stream))
            samples.extend(rows)
        sodium = [row for row in samples if row["species"] == "Na"]
        argon = [row for row in samples if row["species"] == "Ar"]
        expected_sodium = (expected_ranks *
                           int(receipt["sample_count_per_species_per_rank"]))
        linked_argon_count = int(receipt["argon_generated_count"])
        if (len(sodium) != expected_sodium or len(argon) != linked_argon_count
                or len(argon) <
                int(receipt["sample_count_per_species_per_rank"])):
            raise ValueError(
                "I18 Na/Ar sample counts disagree with the linked receipt or "
                "the fixed-attempt Ar ensemble is too small")

        application = work / "application.in"
        sticking_table = parse_temperature_probability_table(
            application, "sodium_sticking_probability_by_temperature_k")
        sodium_accommodation = float(moon_setting(
            application, "sodium_accommodation_coefficient"))
        if moon_setting(application, "sodium_accommodation_speed_model") != \
                "kinetic-energy":
            raise MissingPrerequisite(
                "I18 analytic oracle is frozen to the approved Na "
                "kinetic-energy accommodation law")
        if not 0.0 < sodium_accommodation <= 1.0:
            raise ValueError("invalid I18 Na accommodation coefficient")

        boltzmann = 1.380649e-23
        na_mass = float(receipt["na_mass_kg"])
        ar_mass = float(receipt["argon_mass_kg"])
        na_reemitted = [row for row in sodium
                        if int(float(row["status"])) != 0]
        na_expected_reemitted = 0.0
        na_binomial_variance = 0.0
        for row in sodium:
            probability = 1.0 - interpolate_bounded_table(
                sticking_table, float(row["temperature_k"]))
            na_expected_reemitted += probability
            na_binomial_variance += probability * (1.0 - probability)
        na_half_width_count = 1.96 * math.sqrt(na_binomial_variance)
        na_binomial_ok = (abs(len(na_reemitted) - na_expected_reemitted) <=
                          na_half_width_count)

        na_energy: list[float] = []
        na_cosine: list[float] = []
        for row in na_reemitted:
            temperature = float(row["temperature_k"])
            incident = float(row["incident_speed_m_s"])
            output = float(row["output_speed_m_s"])
            cosine = float(row["outward_cosine"])
            numerator = (output * output -
                         (1.0 - sodium_accommodation) * incident * incident)
            denominator = (sodium_accommodation * 2.0 * boltzmann *
                           temperature / na_mass)
            value = numerator / denominator
            if value < -1.0e-10 or not 0.0 <= cosine <= 1.0:
                raise ValueError("I18 Na production sample is non-physical")
            na_energy.append(max(0.0, value))
            na_cosine.append(cosine)

        ar_energy: list[float] = []
        ar_cosine: list[float] = []
        for row in argon:
            temperature = float(row["temperature_k"])
            speed = float(row["output_speed_m_s"])
            cosine = float(row["outward_cosine"])
            value = ar_mass * speed * speed / (2.0 * boltzmann * temperature)
            if value < 0.0 or not 0.0 <= cosine <= 1.0:
                raise ValueError("I18 Ar production sample is non-physical")
            ar_energy.append(value)
            ar_cosine.append(cosine)

        def maxwell_energy_cdf(value: float) -> float:
            if value <= 0.0:
                return 0.0
            root = math.sqrt(value)
            return math.erf(root) - 2.0 * root * math.exp(-value) / math.sqrt(math.pi)

        def flux_energy_cdf(value: float) -> float:
            return 0.0 if value <= 0.0 else 1.0 - (1.0 + value) * math.exp(-value)

        distributions: list[tuple[str, list[float], Callable[[float], float]]] = [
            ("Na_energy", na_energy, maxwell_energy_cdf),
            ("Na_outward_cosine", na_cosine,
             lambda value: min(1.0, max(0.0, value))),
            ("Ar_flux_energy", ar_energy, flux_energy_cdf),
            ("Ar_outward_cosine", ar_cosine,
             lambda value: min(1.0, max(0.0, value)) ** 2),
        ]
        ks_metrics: dict[str, dict[str, float | bool]] = {}
        comparison_rows: list[dict[str, object]] = []
        distributions_ok = True
        for name, values, cdf in distributions:
            distance = one_sample_ks(values, cdf)
            threshold = 1.36 / math.sqrt(len(values))
            passed = distance <= threshold
            distributions_ok = distributions_ok and passed
            ks_metrics[name] = {"sample_count": len(values),
                                "ks_distance": distance,
                                "ks_95_threshold": threshold,
                                "pass": passed}
            ordered = sorted(values)
            stride = max(1, len(ordered) // 512)
            for index in range(stride - 1, len(ordered), stride):
                value = ordered[index]
                comparison_rows.append({
                    "comparison_id": f"I18_{name}_CDF",
                    "x_value": f"{value:.17g}",
                    "x_unit": "dimensionless",
                    "reference_value": f"{cdf(value):.17g}",
                    "model_value": f"{(index + 1) / len(ordered):.17g}",
                    "y_unit": "cumulative-probability",
                    "uncertainty_value": "",
                    "selection": "verification",
                })
        actual_fraction = len(na_reemitted) / len(sodium)
        expected_fraction = na_expected_reemitted / len(sodium)
        comparison_rows.append({
            "comparison_id": "I18_Na_reemission_fraction",
            "x_value": f"{float(sodium[0]['temperature_k']):.17g}",
            "x_unit": "K",
            "reference_value": f"{expected_fraction:.17g}",
            "model_value": f"{actual_fraction:.17g}",
            "y_unit": "probability",
            "uncertainty_value": f"{na_half_width_count / len(sodium):.17g}",
            "selection": "verification",
        })

        # The linked C++ adapter temporarily installs a known atom population
        # in one ordinary logical-surface element, invokes the production
        # thermal-desorption callback for both species, and restores the
        # original inventory.  Evaluate the cited residence laws independently
        # from explicit input values; no production output is used to construct
        # either reference rate.
        residence_temperature = float(
            receipt["residence_fixture_temperature_k"])
        residence_atoms = float(receipt["residence_fixture_atoms"])
        electron_volt_j = 1.602176634e-19
        planck_j_s = 6.62607015e-34
        gas_constant_j_mol_k = 8.31446261815324
        na_activation_j = float(moon_setting(
            application, "sodium_activation_energy_ev")) * electron_volt_j
        na_attempt_hz = float(moon_setting(
            application, "sodium_attempt_frequency_hz"))
        ar_activation_j_mol = float(moon_setting(
            application, "argon_desorption_energy_j_mol"))
        na_residence_s = math.exp(
            na_activation_j / (boltzmann * residence_temperature)) / \
            na_attempt_hz
        ar_residence_s = (
            planck_j_s / (boltzmann * residence_temperature) *
            math.exp(ar_activation_j_mol /
                     (gas_constant_j_mol_k * residence_temperature)))
        expected_residence_rates = {
            "Na": residence_atoms / na_residence_s,
            "Ar": residence_atoms / ar_residence_s,
        }
        linked_residence_rates = {
            "Na": float(receipt["linked_na_desorption_rate_atoms_s"]),
            "Ar": float(receipt["linked_argon_desorption_rate_atoms_s"]),
        }
        residence_relative_errors: dict[str, float] = {}
        residence_ok = True
        for species_index, species in enumerate(("Na", "Ar")):
            expected_rate = expected_residence_rates[species]
            linked_rate = linked_residence_rates[species]
            relative_error = abs(linked_rate - expected_rate) / expected_rate
            residence_relative_errors[species] = relative_error
            residence_ok = residence_ok and relative_error <= 1.0e-12
            comparison_rows.append({
                "comparison_id": "I18_linked_residence_desorption_rate",
                "x_value": species_index,
                "x_unit": "species-index-0-Na-1-Ar",
                "reference_value": f"{expected_rate:.17g}",
                "model_value": f"{linked_rate:.17g}",
                "y_unit": "atoms-s^-1",
                "uncertainty_value": "",
                "selection": "verification",
            })
        comparison_path = work / "comparison.csv"
        with comparison_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(
                stream, fieldnames=list(comparison_rows[0]), lineterminator="\n")
            writer.writeheader(); writer.writerows(comparison_rows)
        artifacts = render_comparison(
            work, "I18 linked Na/Ar surface-distribution verification")
        source_rate = float(receipt["argon_source_rate_atoms_s"])
        source_rate_ok = math.isclose(source_rate, 3.8e20, rel_tol=1.0e-15)
        metrics.update({
            "na_total_samples": len(sodium),
            "na_reemitted_samples": len(na_reemitted),
            "na_expected_reemitted_samples": na_expected_reemitted,
            "na_reemission_95_half_width_samples": na_half_width_count,
            "na_reemission_fraction": actual_fraction,
            "na_expected_reemission_fraction": expected_fraction,
            "na_binomial_95_pass": na_binomial_ok,
            "argon_total_samples": len(argon),
            "argon_attempt_count_per_rank": int(
                receipt["argon_attempt_count_per_rank"]),
            "rank_count_with_owned_surface_facet": int(
                receipt["rank_count_with_owned_surface_facet"]),
            "argon_source_rate_atoms_s": source_rate,
            "argon_source_rate_pass": source_rate_ok,
            "residence_fixture_temperature_k": residence_temperature,
            "residence_fixture_atoms": residence_atoms,
            "expected_desorption_rate_atoms_s": expected_residence_rates,
            "linked_desorption_rate_atoms_s": linked_residence_rates,
            "residence_rate_relative_error": residence_relative_errors,
            "residence_rate_pass": residence_ok,
            "ks": ks_metrics,
            "random_seed_base": int(receipt["random_seed_base"]),
            "artifacts": artifacts,
        })
        passed = (na_binomial_ok and distributions_ok and source_rate_ok and
                  residence_ok)
        return (PASS if passed else FAIL), metrics, (
            "linked production Na/Ar residence rates, Na sticking/"
            "accommodation, and Ar Maxwell-flux distributions against frozen "
            "analytical and 95% stochastic gates")
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
    external_data = {
        str(path.relative_to(SPICE_ROOT)): sha256(path) for path in kernels
    }
    if test_id == "I14":
        d01 = Path(os.environ.get(
            "MOON_VALIDATION_DATA_ROOT",
            "/data/vtenishe/moon_validation_data")) / "lola/raw"
        for name in ("LDEM_4.IMG", "LDEM_4.LBL"):
            path = d01 / name
            if path.is_file():
                external_data[str(path)] = sha256(path)
    if test_id == "I16" and application_input is not None:
        table_matches = re.findall(
            r"(?mi)^\s*diviner_temperature_table_file\s*=\s*([^!\s]+)\s*$",
            application_input.read_text(encoding="utf-8"))
        if len(table_matches) != 1:
            raise ValueError(
                f"expected one Diviner table path in {application_input}")
        table_text = table_matches[0]
        table = Path(table_text)
        for path in (table, table.parent / "benchmark_points.csv",
                     table.parent / "d03_gcp_qa.json"):
            if path.is_file():
                external_data[str(path)] = sha256(path)
    if test_id == "I17" and application_input is not None:
        # Hash the exact processed evidence and the pre-holdout contract named
        # by the derived application receipt.  The processor logs/QA are also
        # retained when present; absence of the required pair fails closed.
        for key in ("diviner_eclipse_history_file",
                    "diviner_eclipse_selection_contract_file"):
            path = Path(moon_setting(application_input, key, preserve_case=True))
            if not path.is_file():
                raise MissingPrerequisite(
                    f"I17 configured {key} is unavailable: {path}")
            external_data[str(path)] = sha256(path)
        history = Path(moon_setting(
            application_input, "diviner_eclipse_history_file", preserve_case=True))
        for path in (history.parent / "d03_rdr_qa.json",
                     history.parent / "processor.stdout.log",
                     history.parent / "processor.stderr.log"):
            if path.is_file():
                external_data[str(path)] = sha256(path)
    # Bind the package-level qualification records as well as the numerical
    # products.  This prevents a copied processed table from appearing
    # qualified after its README/provenance/QA contract is removed or changed.
    dataset_dependencies = {
        "I14": ("D01",), "I15": ("D02",), "I16": ("D03",),
        "I17": ("D03",), "I20": ("D02", "D03"),
    }
    for dataset_id in dataset_dependencies.get(test_id, ()):
        ready, reason, evidence = dataset_package_preflight(dataset_id)
        if not ready:
            raise MissingPrerequisite(reason)
        evidence_files = evidence.get("files")
        if not isinstance(evidence_files, dict):
            raise ValueError(f"{dataset_id} evidence has no file table")
        for record in evidence_files.values():
            if not isinstance(record, dict):
                raise ValueError(f"{dataset_id} evidence contains a bad file record")
            external_data[str(record["path"])] = str(record["sha256"])
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
        "random_seed": seed,
        "species_fixture": ("Na,NA_PLUS" if test_id == "I08" else
                            ("Na,Ar" if test_id in ("I18", "I19", "I20")
                             else "Na")),
        "configuration": ("no-SPICE" if test_id in ("I05", "I08") else
                          ("canonical-test_Moon" if test_id == "I04" else "SPICE-DE421")),
        "inputs": {str(path): sha256(path) for path in files if path.is_file()},
        "external_data": external_data,
    }


def run_i11_campaign(executable: Path, output_root: Path,
                     application_input: Path | None,
                     mpi_ranks: int = 1,
                     campaign_mode: str = "qualification") -> int:
    """Run the approved memory-budget or legacy FULL production campaign."""
    output = output_root / "I11"
    output.mkdir(parents=True, exist_ok=True)
    result: dict[str, object] = {"id": "I11", "status": ERROR,
                                "metrics": {}, "artifacts": {}}
    try:
        if application_input is None:
            raise ValueError("I11 requires --input with the reviewed M1 receipt")
        if campaign_mode == "regression":
            # I12 owns the single bounded linked production run in the default
            # regression campaign. I11 now has a reviewed 5000-step ceiling,
            # but four production meshes remain too costly for an ordinary
            # developer regression. Preserve a complete non-execution record
            # and never label it PASS.
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
                "reason": "I11 qualification is reserved for qualification mode",
                "input": {"path": str(application_input.resolve()),
                          "sha256": sha256(application_input.resolve())},
                "mpi_ranks_requested": mpi_ranks,
            }
            (output / "run_manifest.json").write_text(
                json.dumps(skipped_manifest, indent=2, sort_keys=True) + "\n",
                encoding="utf-8")
            result.update(
                status=SKIPPED,
                reason=("REGRESSION mode does not launch the bounded "
                        "I11 convergence campaign; use --campaign-mode "
                        "qualification for M1 evidence"),
                metrics={"linked_model_executed": False,
                         "mpi_ranks_requested": mpi_ranks},
                artifacts={"input_snapshot": str(output / "application.in"),
                           "run_manifest": str(output / "run_manifest.json"),
                           "stdout": str(output / "stdout.log"),
                           "stderr": str(output / "stderr.log")})
            raise StopIteration
        levels = qualification_levels(application_input.resolve())
        if moon_setting(application_input.resolve(),
                        "mesh_output_tecplot") != "false":
            raise MissingPrerequisite(
                "I11 requires mesh_output_tecplot=false; full-volume text "
                "visualization is not a science artifact and can exhaust disk")
        shutil.copy2(application_input.resolve(), output / "application.in")
        (output / "stdout.log").write_text(
            "I11 aggregate runner: per-resolution model stdout is retained "
            "under " + ", ".join(f"{level}km/" for level in levels) +
            ".\n", encoding="utf-8")
        (output / "stderr.log").write_text(
            "I11 aggregate runner: per-resolution model stderr is retained "
            "beside each derived receipt.\n", encoding="utf-8")
        initial_manifest = {
            "schema": "srcMoon-I11-campaign-manifest-v2",
            "source_revision": git_value("rev-parse", "HEAD"),
            "template": {"path": str(application_input.resolve()),
                         "sha256": sha256(application_input.resolve())},
            "random_seed": 731993, "mpi_ranks": mpi_ranks,
            "openmp_threads": 1,
            "sample_length_particle_steps": I11_SAMPLE_LENGTH_STEPS,
            "maximum_particle_steps_per_resolution":
                I11_MAXIMUM_PARTICLE_STEPS,
            "maximum_complete_samples_per_resolution":
                I11_MAXIMUM_COMPLETE_SAMPLES,
            "particle_restart_writes": False,
            "resolution_sequence_km": list(levels),
            "levels": {}}
        (output / "run_manifest.json").write_text(
            json.dumps(initial_manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8")
        level_metrics: dict[int, dict[str, dict[str, float | str]]] = {}
        level_receipts: dict[int, dict[str, object]] = {}
        level_summaries: dict[str, object] = {}
        bounded_failures: list[str] = []
        completed_diagnostic_levels: list[int] = []
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
            level_receipts[resolution] = parsed
            level_manifest_entry = {  # type: ignore[var-annotated]
                "input_sha256": sha256(derived_input),
                "manifest_sha256": sha256(level / "run_manifest.json"),
            }
            initial_manifest["levels"][f"{resolution}km"] = level_manifest_entry  # type: ignore[index]
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
            bounded_controls_ok = (
                int(parsed["base_block_length_steps"]) ==
                    I11_SAMPLE_LENGTH_STEPS and
                int(parsed["maximum_particle_steps"]) ==
                    I11_MAXIMUM_PARTICLE_STEPS and
                int(parsed["production_particle_steps_executed"]) <=
                    I11_MAXIMUM_PARTICLE_STEPS and
                int(parsed["particle_restart_writes_enabled"]) == 0)
            wiring_ok = (bounded_controls_ok and
                         parsed.get("simulation_environment_mode") == "frozen" and
                         et_drift == 0.0 and
                         all(value == 0.0 for value in vector_drift.values()) and
                         float(parsed["particle_global_timestep_s"]) > 0.0 and
                         int(parsed["random_seed"]) == 731993)
            if not wiring_ok:
                raise ValueError(
                    f"I11 {resolution} km bounded-control, frozen-environment, "
                    "or deterministic wiring invariant failed")

            criteria_met = int(parsed["convergence_criteria_met"]) == 1
            if not criteria_met:
                steps = int(parsed["production_particle_steps_executed"])
                blocks = int(parsed["base_blocks_executed"])
                cap_reached = int(parsed["iteration_limit_reached"]) == 1
                if (not cap_reached or steps != I11_MAXIMUM_PARTICLE_STEPS or
                        blocks != I11_MAXIMUM_COMPLETE_SAMPLES or
                        int(parsed["convergence_metrics_available"]) != 0):
                    raise ValueError(
                        f"I11 {resolution} km stopped without either qualified "
                        "metrics or the declared ten-sample iteration-cap "
                        "receipt")
                completed_diagnostic_levels.append(resolution)
                diagnostic = render_i11_stage_diagnostic(
                    output, completed_diagnostic_levels, FAIL)
                level_result = {
                    "id": "I11", "resolution_km": resolution,
                    "status": FAIL,
                    "reason": (
                        "5000-particle-step limit reached before the unchanged "
                        "independence, equilibrium, and CI criteria were met"),
                    "metrics": {
                        "sample_length_particle_steps":
                            I11_SAMPLE_LENGTH_STEPS,
                        "maximum_particle_steps":
                            I11_MAXIMUM_PARTICLE_STEPS,
                        "maximum_complete_samples":
                            I11_MAXIMUM_COMPLETE_SAMPLES,
                        "particle_steps_executed": steps,
                        "base_blocks": blocks,
                        "minimum_batches_required": 9,
                        "particle_restart_writes_enabled": False,
                        "ephemeris_time_drift_s": et_drift,
                        "external_vector_absolute_drifts": vector_drift,
                    },
                    "artifacts": {
                        "block_metrics": str(level / "block_metrics.csv"),
                        "resource_metrics": str(level / "resource_metrics.csv"),
                        "process_resource_metrics": str(
                            level / "process_resource_metrics.csv"),
                        "stage_diagnostic": diagnostic["artifacts"],
                    },
                }
                (level / "result.json").write_text(
                    json.dumps(level_result, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")
                # Bind every diagnostic needed to audit a capped scientific
                # FAIL. Missing raw metrics are infrastructure corruption and
                # therefore make sha256() raise into the ERROR path instead of
                # leaving a superficially complete campaign manifest.
                level_manifest_entry.update({
                    "receipt_sha256": sha256(receipt),
                    "result_sha256": sha256(level / "result.json"),
                    "block_metrics_sha256": sha256(
                        level / "block_metrics.csv"),
                    "resource_metrics_sha256": sha256(
                        level / "resource_metrics.csv"),
                    "process_resource_metrics_sha256": sha256(
                        level / "process_resource_metrics.csv"),
                    "stage_diagnostic_csv_sha256": sha256(
                        level / "stage_diagnostic.csv"),
                    "stage_diagnostic_png_sha256": sha256(
                        level / "stage_diagnostic.png"),
                    "stage_diagnostic_eps_sha256": sha256(
                        level / "stage_diagnostic.eps"),
                    "stage_diagnostic_metadata_sha256": sha256(
                        level / "stage_diagnostic.json"),
                })
                (output / "run_manifest.json").write_text(
                    json.dumps(initial_manifest, indent=2, sort_keys=True) +
                    "\n", encoding="utf-8")
                level_summaries[f"{resolution}km"] = level_result
                bounded_failures.append(f"{resolution}km")
                # Run every requested resolution even though the overall
                # campaign is already FAIL. This honors the per-test 5000-step
                # cap and preserves a comparable bounded diagnostic set.
                continue

            qualified = level / "qualified_metrics.csv"
            level_metrics[resolution] = read_qualified_metrics(qualified)
            level_manifest_entry["qualified_metrics_sha256"] = sha256(qualified)
            (output / "run_manifest.json").write_text(
                json.dumps(initial_manifest, indent=2, sort_keys=True) + "\n",
                encoding="utf-8")
            if int(parsed["retained_independent_batch_count"]) < 9:
                raise ValueError(
                    f"I11 {resolution} km claims convergence with fewer than "
                    "nine retained independent batches")

            unresolved = []
            for name, metric in level_metrics[resolution].items():
                gate = 0.0333 if metric["classification"] == "pointwise" else 0.0167
                if float(metric["relative_ci95_half_width"]) > gate:
                    unresolved.append(name)
            level_status = PASS if wiring_ok and not unresolved else (
                SKIPPED if wiring_ok else FAIL)
            completed_diagnostic_levels.append(resolution)
            diagnostic = render_i11_stage_diagnostic(
                output, completed_diagnostic_levels, level_status)
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
                                parsed["retained_independent_batch_count"]),
                            "stage_diagnostic_notice": diagnostic["notice"]},
                "artifacts": {"stage_diagnostic": diagnostic["artifacts"]},
            }
            (level / "result.json").write_text(
                json.dumps(level_result, indent=2, sort_keys=True) + "\n",
                encoding="utf-8")
            level_manifest_entry.update({
                "receipt_sha256": sha256(receipt),
                "result_sha256": sha256(level / "result.json"),
                "block_metrics_sha256": sha256(level / "block_metrics.csv"),
                "resource_metrics_sha256": sha256(
                    level / "resource_metrics.csv"),
                "process_resource_metrics_sha256": sha256(
                    level / "process_resource_metrics.csv"),
                "stage_diagnostic_csv_sha256": sha256(
                    level / "stage_diagnostic.csv"),
                "stage_diagnostic_png_sha256": sha256(
                    level / "stage_diagnostic.png"),
                "stage_diagnostic_eps_sha256": sha256(
                    level / "stage_diagnostic.eps"),
                "stage_diagnostic_metadata_sha256": sha256(
                    level / "stage_diagnostic.json"),
            })
            (output / "run_manifest.json").write_text(
                json.dumps(initial_manifest, indent=2, sort_keys=True) + "\n",
                encoding="utf-8")
            level_summaries[f"{resolution}km"] = level_result
            if level_status != PASS:
                result.update(status=level_status,
                              reason=level_result["reason"],
                              metrics={"levels": level_summaries})
                raise StopIteration

        if bounded_failures:
            result.update(
                status=FAIL,
                reason=("iteration cap reached before statistical qualification "
                        "at: " + ", ".join(bounded_failures)),
                metrics={"levels": level_summaries,
                         "sample_length_particle_steps":
                             I11_SAMPLE_LENGTH_STEPS,
                         "maximum_particle_steps_per_resolution":
                             I11_MAXIMUM_PARTICLE_STEPS,
                         "particle_restart_writes": False},
                artifacts={"campaign_manifest": str(
                    output / "run_manifest.json")})
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
        shutil.copy2(output / f"{levels[-1]}km" / "receipt.txt", output / "receipt.txt")
        # The campaign manifest binds the four child manifests and derived
        # receipts without hiding their individual commands or resource logs.
        # ``initial_manifest`` has been extended after every completed stage
        # with both accepted-science and explicitly unqualified diagnostic
        # hashes.  Preserve those entries here; rebuilding a smaller levels
        # dictionary would silently detach the requested meeting artifacts
        # from the campaign provenance record.
        campaign_manifest = initial_manifest
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
        # I13 freezes the numerically qualified finest I11 case, not a coarser
        # human-authored template. Derive it with the same audited transformer
        # used by I11, preserve it as the I13 input snapshot, and later require
        # its hash to equal I11/<finest>km/application.in before baseline scoring.
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
    seed = 731993 if test_id in (
        "I09", "I11", "I12", "I13", "I17", "I18", "I19", "I20") else None
    result: dict[str, object] = {"id": test_id, "status": ERROR, "metrics": {},
                                "artifacts": {"stdout": str(stdout_path), "stderr": str(stderr_path), "receipt": str(receipt_path)}}
    try:
        if test_id in {"I18", "I19", "I20"}:
            if application_input is None:
                raise ValueError(
                    f"{test_id} requires its reviewed M2 input template")
            effective_input = output / "application.in"
            prepare_na_ar_input(application_input.resolve(), effective_input)
            # The command is assembled before per-ID receipt derivation so the
            # common I01--I17 path stays byte-for-byte stable.  Replace only
            # the final input argument with the preserved two-species receipt.
            command[-1] = str(effective_input.resolve())
        if 14 <= int(test_id[1:]) <= 20:
            # M2 is deliberately layered after the FULL M1 freeze.  A direct
            # I14--I20 invocation cannot infer that evidence from a local unit
            # pass or from an arbitrary old output directory; it must find the
            # I13 result in this same preserved campaign.
            predecessor = output_root / "I13" / "result.json"
            predecessor_status = None
            if predecessor.is_file():
                predecessor_status = json.loads(
                    predecessor.read_text(encoding="utf-8")).get("status")
            if predecessor_status != PASS:
                if application_input is None:
                    raise ValueError(f"{test_id} requires its M2 --input receipt")
                if effective_input.resolve() != (output / "application.in").resolve():
                    shutil.copy2(effective_input.resolve(), output / "application.in")
                stdout_path.write_text(
                    f"{test_id} linked model was not launched because I13 in "
                    "this campaign is not PASS.\n", encoding="utf-8")
                stderr_path.write_text("", encoding="utf-8")
                skip_manifest = {
                    "schema": "srcMoon-M2-predecessor-skip-v1",
                    "source_revision": git_value("rev-parse", "HEAD"),
                    "test_id": test_id,
                    "execution": "not-run",
                    "predecessor": {"id": "I13", "status": predecessor_status,
                                    "path": str(predecessor)},
                    "input": {"path": str((output / "application.in").resolve()),
                              "sha256": sha256(output / "application.in")},
                    "mpi_ranks_requested": mpi_ranks,
                }
                (output / "run_manifest.json").write_text(
                    json.dumps(skip_manifest, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")
                result.update(
                    status=SKIPPED,
                    reason=(f"{test_id} requires a PASS from I13 in the same "
                            "approved qualification campaign; later-layer M2 "
                            "execution cannot replace the unresolved M1 gate"),
                    metrics={"linked_model_executed": False,
                             "predecessor_status": predecessor_status},
                    artifacts={"input_snapshot": str(output / "application.in"),
                               "run_manifest": str(output / "run_manifest.json"),
                               "stdout": str(stdout_path),
                               "stderr": str(stderr_path)})
                raise StopIteration
        # A linked observational/environment comparison may inspect its data
        # only after the preceding numerical gate passes and only when the
        # external package itself is qualified.  Processed CSV files generated
        # in a prior work directory are not a substitute for package metadata.
        dataset_dependencies = {
            "I14": ("D01",), "I15": ("D02",), "I16": ("D03",),
            "I17": ("D03",), "I20": ("D02", "D03"),
        }
        data_evidence: dict[str, object] = {}
        data_blockers: list[str] = []
        for dataset_id in dataset_dependencies.get(test_id, ()):
            ready, reason, evidence = dataset_package_preflight(dataset_id)
            data_evidence[dataset_id] = evidence
            if not ready:
                data_blockers.append(reason)
        if data_blockers:
            if application_input is None:
                raise ValueError(f"{test_id} requires its M2 input receipt")
            if effective_input.resolve() != (output / "application.in").resolve():
                shutil.copy2(effective_input.resolve(), output / "application.in")
            reason = "; ".join(data_blockers)
            stdout_path.write_text(
                f"{test_id} linked model was not launched: {reason}.\n",
                encoding="utf-8")
            stderr_path.write_text("", encoding="utf-8")
            data_manifest = {
                "schema": "srcMoon-M2-data-skip-v1",
                "source_revision": git_value("rev-parse", "HEAD"),
                "test_id": test_id,
                "execution": "not-run",
                "reason": reason,
                "datasets": data_evidence,
                "input": {
                    "path": str((output / "application.in").resolve()),
                    "sha256": sha256(output / "application.in"),
                },
                "mpi_ranks_requested": mpi_ranks,
            }
            (output / "run_manifest.json").write_text(
                json.dumps(data_manifest, indent=2, sort_keys=True) + "\n",
                encoding="utf-8")
            result.update(
                status=SKIPPED, reason=reason,
                metrics={"linked_model_executed": False,
                         "dataset_status": "NOT READY"},
                artifacts={"input_snapshot": str(output / "application.in"),
                           "run_manifest": str(output / "run_manifest.json"),
                           "stdout": str(stdout_path),
                           "stderr": str(stderr_path)})
            raise StopIteration
        if test_id == "I16":
            if application_input is None:
                raise ValueError("I16 requires its reviewed M2 input template")
            effective_input = output / "application.in"
            prepare_i16_input(
                application_input.resolve(), effective_input, output_root)
        if test_id == "I17":
            if application_input is None:
                raise ValueError("I17 requires its reviewed M2 input template")
            effective_input = output / "application.in"
            prepare_i17_input(
                application_input.resolve(), effective_input, output_root)
        if test_id in ("I16", "I17"):
            # These tests derive an observation-bearing receipt after the
            # initial command was assembled. Bind the executable to that same
            # receipt before writing its manifest: launching the original
            # template would use `none` paths (I17), or the wrong thermal mode
            # (I16), despite provenance describing the staged observation data.
            # The final argument is the path following the existing -input.
            command[-1] = str(effective_input.resolve())
        unavailable_m2_capability = {
            "I15": (
                "the production terrain service currently evaluates a "
                "point-Sun ray; D02 requires the published extended-solar-disc "
                "18.6-year/one-hour operator before observational scoring"),
            "I19": (
                "production-wide weighted counters for external injection, "
                "desorption, adsorption, escape, and numerical domain loss "
                "are not yet all wired for synchronized Na/Ar closure"),
            "I20": (
                "I20 requires passing linked I15, I17, I18, and I19 evidence; "
                "the prepared D02/PRP tables do not substitute for those "
                "production capability gates"),
        }
        if test_id in unavailable_m2_capability:
            if application_input is None:
                raise ValueError(f"{test_id} requires its M2 input receipt")
            if effective_input.resolve() != (output / "application.in").resolve():
                shutil.copy2(effective_input.resolve(), output / "application.in")
            reason = unavailable_m2_capability[test_id]
            stdout_path.write_text(
                f"{test_id} linked model was not launched: {reason}.\n",
                encoding="utf-8")
            stderr_path.write_text("", encoding="utf-8")
            capability_manifest = {
                "schema": "srcMoon-M2-capability-skip-v1",
                "source_revision": git_value("rev-parse", "HEAD"),
                "test_id": test_id,
                "execution": "not-run",
                "reason": reason,
                "input": {"path": str((output / "application.in").resolve()),
                          "sha256": sha256(output / "application.in")},
                "mpi_ranks_requested": mpi_ranks,
            }
            (output / "run_manifest.json").write_text(
                json.dumps(capability_manifest, indent=2, sort_keys=True) + "\n",
                encoding="utf-8")
            result.update(
                status=SKIPPED, reason=reason,
                metrics={"linked_model_executed": False,
                         "capability_status": "EXPERIMENTAL / NOT QUALIFIED"},
                artifacts={"input_snapshot": str(output / "application.in"),
                           "run_manifest": str(output / "run_manifest.json"),
                           "stdout": str(stdout_path),
                           "stderr": str(stderr_path)})
            raise StopIteration
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
                "reason": "I13 requires passing approved I11 qualification evidence",
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
            finest = qualification_levels(application_input.resolve())[-1]
            # Do not build another expensive mesh when the science predecessor
            # is already known to be absent/failed. Never manufacture baseline
            # evidence by running geometry alone after a failed I11 campaign.
            predecessor = output_root / "I11/result.json"
            if not predecessor.is_file() or json.loads(
                    predecessor.read_text(encoding="utf-8")).get("status") != PASS:
                shutil.copy2(application_input, effective_input)
                stdout_path.write_text("I13 not launched: I11 is not PASS.\n")
                stderr_path.write_text("")
                (output / "run_manifest.json").write_text(json.dumps({
                    "execution": "not-run", "source_revision": git_value("rev-parse", "HEAD"),
                    "input_sha256": sha256(effective_input),
                    "resolution_km": finest, "mpi_ranks": mpi_ranks,
                    "reason": "I11 prerequisite not PASS"}, indent=2) + "\n")
                raise MissingPrerequisite("I13 requires same-campaign I11 PASS before launch")
            derive_i11_input(application_input.resolve(), effective_input, finest)  # type: ignore[arg-type]
        run_manifest = manifest(test_id, executable, command, seed,
                                effective_input, mpi_ranks)
        (output / "run_manifest.json").write_text(json.dumps(run_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        if (application_input is not None and test_id != "I13" and
                effective_input == application_input):
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
            if test_id == "I17":
                # Numerical qualification precedes all residual computation.
                # Each refined execution independently spins up the identical
                # full production state using the same physical parameters.
                # Changing only dt includes initialization sensitivity in the
                # error estimate rather than reusing a coarse initial state.
                from i17_convergence import compare_time_steps
                base_step = float(moon_setting(effective_input, "thermal_time_step_s"))
                if not math.isfinite(base_step) or base_step <= 0.0:
                    raise ValueError("I17 thermal time step must be finite and positive")
                directories = [output]
                receipts = [receipt]
                for divisor in (2, 4):
                    refined = output / f"dt_div_{divisor}"
                    refined.mkdir(exist_ok=False)
                    refined_input = refined / "application.in"
                    text = replace_moon_setting(
                        effective_input.read_text(encoding="utf-8"),
                        "thermal_time_step_s", format(base_step / divisor, ".17g"),
                        effective_input)
                    # Keep each run's generated surface and inventory outputs
                    # separate. All science inputs in the staged I17 receipt
                    # are absolute and remain identical across the three runs.
                    for key, name in (("surface_cea_file", "surface.cea"),
                                      ("surface_tecplot_file", "surface.dat"),
                                      ("surface_inventory_checkpoint_file", "inventory.chk")):
                        text = replace_moon_setting(text, key,
                            str((refined / name).resolve()), effective_input)
                    refined_input.write_text(text, encoding="utf-8")
                    refined_receipt = refined / "receipt.txt"
                    refined_command = command[:]
                    refined_command[refined_command.index("-input") + 1] = str(refined_input.resolve())
                    refined_command[refined_command.index("--integration-output") + 1] = str(refined_receipt.resolve())
                    identity = manifest("I17", executable, refined_command, seed,
                                        refined_input, mpi_ranks)
                    (refined / "run_manifest.json").write_text(
                        json.dumps(identity, indent=2, sort_keys=True) + "\n")
                    code = run_monitored(refined_command, runtime_directory,
                        environment, refined / "stdout.log", refined / "stderr.log",
                        refined / "process_resource_metrics.csv")
                    if code != 0 or not refined_receipt.is_file():
                        raise RuntimeError(f"I17 dt/{divisor} exited {code}; see {refined}")
                    finer = parse_receipt(refined_receipt)
                    if (finer.get("test_id") != "I17" or
                            finer.get("holdout_rows_accessed") != 0 or
                            float(finer["thermal_time_step_s"]) != base_step / divisor):
                        raise ValueError("I17 refined receipt disagrees with its configuration")
                    directories.append(refined)
                    receipts.append(finer)
                passed, convergence = compare_time_steps(
                    directories, output / "thermal_time_step_convergence.csv")
                # The existing <50 K/update rule applies to every execution,
                # including the coarse run, and is never inferred from dt alone.
                jumps = [float(item["maximum_production_update_temperature_change_k"])
                         for item in receipts]
                updates_ok = all(math.isfinite(value) and 0.0 <= value < 50.0
                                 for value in jumps)
                convergence["maximum_update_changes_K"] = jumps
                convergence["time_steps_s"] = [base_step / d for d in (1, 2, 4)]
                if passed and updates_ok:
                    status, metrics, reason = check(test_id, receipts[-1], directories[-1])
                    metrics["thermal_convergence"] = convergence
                else:
                    status, metrics, reason = FAIL, convergence, (
                        "I17 thermal numerical gate failed; no observational residuals or plots produced")
            else:
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
                        choices=[f"I{number:02d}" for number in range(1, 21)])
    parser.add_argument("--amps", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--input", type=Path,
                        help="validated srcMoon .in receipt for I11-I20")
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
