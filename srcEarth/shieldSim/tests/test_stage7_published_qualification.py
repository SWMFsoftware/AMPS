#!/usr/bin/env python3
"""Stage-7 published-component representative-matrix qualification.

This gate uses every bin of the downloaded AMS-02 A=1 proton component.  It is
intentionally labelled an incomplete GCR *component* qualification.  It must
not be interpreted as a complete GCR environment.  The GOES-16 product is
reloaded and checked here, but remains ingestion/event-integration only because
its sensor-to-boundary transform and angular response are not available.

For the two approved Al/HDPE layer orders, the test compares three independent
paths for silicon-slab ionizing deposition:

1. production sampling from exact published bin edges and isotropic p(mu)=2mu;
2. a reviewable energy/angle quadrature sampled directly by Geant4; and
3. the same quadrature folded through monoenergetic response matrices.

Five fixed seeds are pooled for each direct estimator.  Matrix response cells
use disjoint deterministic seeds.  The water rows execute the same downloaded
source through the validated TID material boundary but are not folded into the
silicon response.  All generated cases and logs remain below RUN_DIR.
"""

from __future__ import annotations

import copy
import concurrent.futures
import hashlib
import json
import math
import os
import pathlib
import shutil
import subprocess
import sys


if len(sys.argv) != 10:
    raise SystemExit(
        "usage: test_stage7_published_qualification.py SHIELDSIM SOURCE_TOOL "
        "RESPONSE_TOOL EXPOSURE_TOOL WEIGHT_TOOL SOURCE_ROOT EVIDENCE_DIR RUN_DIR PYTHON"
    )
EXE = pathlib.Path(sys.argv[1]).resolve()
SOURCE_TOOL = pathlib.Path(sys.argv[2]).resolve()
RESPONSE_TOOL = pathlib.Path(sys.argv[3]).resolve()
EXPOSURE_TOOL = pathlib.Path(sys.argv[4]).resolve()
WEIGHT_TOOL = pathlib.Path(sys.argv[5]).resolve()
ROOT = pathlib.Path(sys.argv[6]).resolve()
EVIDENCE_DIR = pathlib.Path(sys.argv[7]).resolve()
RUN_DIR = pathlib.Path(sys.argv[8]).resolve()
PYTHON = pathlib.Path(sys.argv[9]).resolve()
FIXTURES = ROOT / "tests/reference_spectra"
PLAN_PATH = ROOT / "tests/reference_stage7/stage7_qualification_plan.v1.json"
SEEDS = [710003, 710021, 710047, 710083, 710111]
# TID behind the mixed stack has a heavy nuclear/secondary tail. A 40,000-
# history diagnostic demonstrated that the earlier 4,000-history ensemble was
# not stable enough for the registered 10%/95% comparison, so the five-seed
# direct estimators each receive 8,000 histories (40,000 total per row).
DIRECT_HISTORIES = 8000
PILOT_CELL_HISTORIES = 256
# The production response budget is allocated after an independent 256-history
# pilot.  Each of five independent production replicates receives the same
# frozen Neyman allocation: a 512-history floor protects every cell, and the
# remaining histories follow |source weight|*pilot standard deviation.  A
# single-stream precision study showed that doubling one stream could still
# miss the mixed-stack nuclear tail and understate its sampled variance.  Five
# disjoint response seed families therefore provide an empirical
# between-replicate standard error; the conservative larger of that error and
# the propagated within-matrix error controls the comparison.  Pilot moments
# remain excluded, so data-dependent allocation does not feed a pilot
# fluctuation back into the estimator.
RESPONSE_PRODUCTION_HISTORIES = 73728
RESPONSE_MINIMUM_PER_CELL = 512
RESPONSE_SEED_BASES = [820001, 920001, 1020001, 1120001, 1220001]
# The second adaptive energy refinement adds the 11 interleaved published-bin
# representatives to the original 12-node grid.  Its 99 angle-energy strata
# receive an independent pilot/allocation.  Keeping the base allocation frozen
# preserves all prior unbiased samples; refinement adds information rather than
# changing their stopping counts after their outcomes have been observed.
REFINEMENT_PRODUCTION_HISTORIES = 56320
REFINEMENT_MINIMUM_PER_CELL = 512
# The 23-node refinement reduced the reversed-stack relative difference below
# 10%, but its independent-comparison statistic remained 2.078, just outside
# the registered two-sided 95% limit of 1.96.  The acceptance limit is not
# relaxed.  Instead, a final independent refinement adds the remaining 22
# published energy representatives.  Its history budget is twice the first
# refinement budget because it contains twice as many angle-energy strata;
# this keeps the same 512-history floor and approximately the same mean budget
# per new cell without changing any already-observed stopping count.
EXACT_REFINEMENT_PRODUCTION_HISTORIES = 112640
EXACT_REFINEMENT_MINIMUM_PER_CELL = 512

# Three-point Gauss-Legendre quadrature on [0,1], multiplied by p(mu)=2mu.
# Azimuth uses a separate three-point Gauss-Legendre rule on the irreducible
# square-symmetry interval [0,pi/4]. A two-endpoint development rule biased the
# HDPE->Al continuous-source comparison by 12%; the higher-order rule is fixed
# here and in the registered plan rather than relaxing the acceptance gate.
MU_NODES = [0.1127016653792583, 0.5, 0.8872983346207417]
MU_BASE_WEIGHTS = [5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0]
MU_WEIGHTS = [2.0 * mu * weight for mu, weight in zip(MU_NODES, MU_BASE_WEIGHTS)]
PHI_GAUSS_ABSCISSA = math.sqrt(3.0 / 5.0)
PHI_NODES = [math.pi / 8.0 * (1.0 - PHI_GAUSS_ABSCISSA),
             math.pi / 8.0,
             math.pi / 8.0 * (1.0 + PHI_GAUSS_ABSCISSA)]
PHI_WEIGHTS = [5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0]
# The initial 12-node grid resolves every fourth published bin.  A qualification
# run found a statistically resolved 7.09% TID difference for HDPE->Al, so the
# accepted fine grid interleaves every second published bin.  This is an actual
# response refinement; the comparison threshold is unchanged.
BASE_ENERGY_INDICES = [0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44]
COARSE_ENERGY_INDICES = BASE_ENERGY_INDICES
FINE_ENERGY_INDICES = list(range(0, 45, 2))
REFINEMENT_ENERGY_INDICES = [index for index in FINE_ENERGY_INDICES
                             if index not in BASE_ENERGY_INDICES]
EXACT_ENERGY_INDICES = list(range(45))
EXACT_REFINEMENT_ENERGY_INDICES = [index for index in EXACT_ENERGY_INDICES
                                   if index not in FINE_ENERGY_INDICES]
SELECTED_RESPONSE_SCORE_IDS = [
    "absorber/silicon_slab/ionizing_deposited_energy_MeV",
    "shield_exit/primary/proton/crossing_count",
]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def configured_worker_count() -> int:
    """Return the explicitly bounded number of concurrent shieldSim cases.

    A worker is one independent *serial* shieldSim process.  This is not
    event-level Geant4 threading: every case retains its registered random
    seed, history count, private working directory, and private output prefix.
    The default remains one worker so ordinary CTest behavior is conservative.
    Qualification runs can opt in with SHIELDSIM_STAGE7_WORKERS=N.
    """
    raw = os.environ.get("SHIELDSIM_STAGE7_WORKERS", "1")
    try:
        workers = int(raw)
    except ValueError as error:
        raise ValueError("SHIELDSIM_STAGE7_WORKERS must be an integer") from error
    try:
        available = len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        available = os.cpu_count() or 1
    require(workers >= 1, "SHIELDSIM_STAGE7_WORKERS must be at least 1")
    require(workers <= available,
            f"SHIELDSIM_STAGE7_WORKERS={workers} exceeds {available} available CPUs")
    return workers


def canonical_write(path: pathlib.Path, document: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(document, sort_keys=True, separators=(",", ":")) + "\n",
                    encoding="utf-8")


def campaign_cache_identity() -> dict:
    """Hash every executable, rule, and fixed input that can affect this run.

    Stage-7 qualification is intentionally long.  A killed CTest process must
    be restartable without accepting stale results from different code or
    inputs.  The marker includes this driver, every invoked executable/script,
    the frozen plan, base case/material library, and all published fixtures.
    """
    paths = [
        pathlib.Path(__file__).resolve(), EXE, SOURCE_TOOL, RESPONSE_TOOL,
        EXPOSURE_TOOL, WEIGHT_TOOL, PLAN_PATH,
        ROOT / "examples/stage4_tid_case.json",
        ROOT / "data/material_library.v1.json",
        FIXTURES / "gcr/ams02_20131118_hapi_info.json",
        FIXTURES / "gcr/ams02_20131118_hapi_data.csv",
        FIXTURES / "gcr/ams02_20131118_interpretation.v1.json",
        FIXTURES / "sep/goes16_sgps_gle72_timeseries.csv",
        FIXTURES / "sep/goes16_gle72_event_definition.v1.json",
    ]
    files = []
    for path in paths:
        require(path.is_file(), f"campaign identity input is missing: {path}")
        files.append({"path": str(path),
                      "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    return {"schema_version": "shieldSim.stage7_campaign_cache.v1", "files": files}


def command(argv: list[str], cwd: pathlib.Path, log_name: str,
            success: bool = True) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(
        argv, cwd=cwd, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, check=False,
    )
    (cwd / log_name).write_text(completed.stdout, encoding="utf-8")
    if success:
        require(completed.returncode == 0,
                f"command failed ({completed.returncode}): {' '.join(argv)}\n"
                f"log: {cwd / log_name}\n{completed.stdout[-6000:]}")
    else:
        require(completed.returncode != 0, f"command unexpectedly succeeded: {' '.join(argv)}")
    return completed


def file_sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return "sha256:" + digest.hexdigest()


def compact_response_case(
        directory: pathlib.Path, response_path: pathlib.Path,
        result_path: pathlib.Path, result: dict) -> dict:
    """Validate a full response cell, retain registered scores, and checkpoint.

    A Stage-7 Geant4 result and its unabridged covariance response cell carry
    the same large diagnostic score ledger.  Keeping both for thousands of
    nodes exhausted the test filesystem.  The production C++ response tool
    first validates the complete cell and all of its covariance pairs.  Only
    then is the cell reduced to the two preregistered primitive scores and
    their covariance block.  The compact cell remains a valid v1 response cell
    and is validated again whenever a matrix is assembled.

    The full result is replaced by a content-hashed completion marker because
    response-cell callers need only event closure and score presence.  Direct
    published/quadrature cases retain their full result documents.
    """
    histories = int(result["events"])
    result_score_ids = {item["score_id"] for item in result["raw_transport_scores"]}
    for score_id in SELECTED_RESPONSE_SCORE_IDS:
        require(score_id in result_score_ids,
                f"response result is missing selected score {score_id}")

    # Release the redundant large result before the validator writes its small
    # temporary matrix.  If validation fails, no completion marker is written
    # and the next run regenerates this case from its canonical input.
    result_path.unlink()
    validation_matrix = directory / "compact_validation_matrix.json"
    validation_csv = validation_matrix.with_suffix(".csv")
    command([str(RESPONSE_TOOL), "assemble-selected", str(validation_matrix),
             ",".join(SELECTED_RESPONSE_SCORE_IDS), str(response_path)],
            directory, "compact_validation.log")

    full_cell = json.loads(response_path.read_text(encoding="utf-8"))
    selected = set(SELECTED_RESPONSE_SCORE_IDS)
    full_cell["scores"] = [
        item for item in full_cell["scores"] if item["score_id"] in selected
    ]
    full_cell["score_covariance"] = [
        item for item in full_cell["score_covariance"]
        if item["score_id_i"] in selected and item["score_id_j"] in selected
    ]
    full_cell.pop("secondary_production_by_layer", None)
    require({item["score_id"] for item in full_cell["scores"]} == selected,
            "response-cell compaction lost a selected score")
    require(len(full_cell["score_covariance"]) == len(selected) * (len(selected) + 1) // 2,
            "response-cell compaction did not retain the selected covariance triangle")
    compact_path = directory / "response_cell.compact.tmp"
    canonical_write(compact_path, full_cell)
    os.replace(compact_path, response_path)

    validation_matrix.unlink(missing_ok=True)
    validation_csv.unlink(missing_ok=True)
    completion = {
        "schema_version": "shieldSim.stage7_response_case_completion.v1",
        "events": histories,
        "score_ids": SELECTED_RESPONSE_SCORE_IDS,
        "response_cell_sha256": file_sha256(response_path),
        "retention": "full_cell_validated_then_registered_scores_compacted",
    }
    canonical_write(directory / "response_case_completion.json", completion)
    return {
        "events": histories,
        "raw_transport_scores": [{"score_id": score_id}
                                 for score_id in SELECTED_RESPONSE_SCORE_IDS],
    }


def set_fixed_histories(case: dict, histories: int) -> None:
    case["accuracy"] = {
        "target": 0.10,
        "scope": "statistical",
        "coverage": {
            "standard_k": 1.0,
            "standard_probability": 0.68,
            "interval_probability": 0.95,
            "interval_method": "student_t_mc_plus_normal_systematics_quadrature",
        },
        "systematic_budget": {"status": "incomplete"},
        "quantity_overrides": {},
        "sampling": {
            "min_histories": histories,
            "max_histories": histories,
            "batch_size": histories,
        },
    }


def base_case(stack_id: str, absorber: str, seed: int, histories: int) -> dict:
    case = json.loads((ROOT / "examples/stage4_tid_case.json").read_text(encoding="utf-8"))
    case["case_id"] = f"stage7_{stack_id}_{absorber}"
    case["purpose"] = "verification_test"
    case["material_library"]["path"] = str(ROOT / "data/material_library.v1.json")
    # The response boundary is a plane immediately upstream of the shield, as
    # required for an incident environment fold.  The Stage-4 example places
    # it much farther away for a separate geometry test; retaining that gap
    # here would turn the energy-grid convergence test into a rare geometric-
    # acceptance test for grazing rays.
    case["geometry"]["source_plane"]["transform"]["translation"]["z"] = {
        "value": -22.0, "unit": "mm"
    }
    materials = ("Al", "HDPE") if stack_id == "al_then_hdpe" else ("HDPE", "Al")
    thickness = {"Al": 1.0, "HDPE": 2.0}
    case["geometry"]["shield_stack"]["layers"] = [
        {"layer_id": "layer_0", "material": materials[0],
         "thickness": {"type": "geometric", "value": thickness[materials[0]], "unit": "mm"}},
        {"layer_id": "layer_1", "material": materials[1],
         "thickness": {"type": "geometric", "value": thickness[materials[1]], "unit": "mm"}},
    ]
    item = case["geometry"]["absorbers"][0]
    item["absorber_id"] = absorber
    item["material"] = "Si" if absorber == "silicon_slab" else "Water"
    item["geometry"]["size_x"] = {"value": 50.0, "unit": "mm"}
    item["geometry"]["size_y"] = {"value": 50.0, "unit": "mm"}
    item["geometry"]["thickness"] = {
        # Stage 7 qualifies TID using the already validated Stage-4 1 mm slab.
        # The 0.1 mm silicon fixture remains appropriate for Stage-6 LET but
        # has a rare-deposition tail too large for a 2% response-grid monitor.
        "value": 1.0, "unit": "mm"
    }
    case["quantities"] = ["TID"]
    case["transport"]["random_seed"] = seed
    case["transport"]["production_cut"] = {"value": 0.1, "unit": "mm"}
    case["transport"]["max_step"] = {"value": 0.1, "unit": "mm"}
    set_fixed_histories(case, histories)
    return case


def published_source() -> dict:
    return {
        "mode": "isotropic",
        "spectrum": {
            "format": "ams02_hapi_component",
            "info_path": str(FIXTURES / "gcr/ams02_20131118_hapi_info.json"),
            "data_path": str(FIXTURES / "gcr/ams02_20131118_hapi_data.csv"),
            "interpretation_path": str(FIXTURES / "gcr/ams02_20131118_interpretation.v1.json"),
            "component_id": "proton",
        },
        "energy_limits": {
            "proton": {"min": {"value": 400.0, "unit": "MeV"},
                       "max": {"value": 100.0, "unit": "GeV"}},
            "alpha": {"min": {"value": 1.0, "unit": "MeV"},
                      "max": {"value": 100.0, "unit": "GeV"}},
        },
    }


def set_response_incident(case: dict, incident_id: str, energy: float,
                          mu: float, phi: float) -> None:
    case["source"] = {
        "mode": "directional_plane",
        "spectrum": {
            "format": "monoenergetic_response_cell",
            "incident_id": incident_id,
            "particle": {"name": "proton", "Z": 1, "A": 1},
            "energy": {"coordinate": "total_kinetic_energy", "value": energy, "unit": "MeV"},
            "direction": {"mu": mu, "phi": {"value": phi, "unit": "rad"}},
        },
    }


def run_case(case: dict, directory: pathlib.Path, name: str) -> tuple[dict, pathlib.Path]:
    directory.mkdir(parents=True, exist_ok=True)
    document = copy.deepcopy(case)
    document["case_id"] = name
    document["output"]["prefix"] = str(directory / "run")
    path = directory / "case.json"
    absorber = document["geometry"]["absorbers"][0]["absorber_id"]
    prefix = directory / f"run__{absorber}"
    result_path = pathlib.Path(str(prefix) + "_result.json")
    response_path = pathlib.Path(str(prefix) + "_response_cell.json")
    completion_path = directory / "response_case_completion.json"
    requested_histories = int(document["accuracy"]["sampling"]["max_histories"])
    response_required = (document["source"]["spectrum"]["format"] ==
                         "monoenergetic_response_cell")

    # Reuse a completed case only when its canonical input is exactly the
    # desired document and its event count closes to the fixed history count.
    # Compact response cases use a content-hashed completion marker; direct
    # cases retain and validate the ordinary full result document.
    if response_required and path.is_file() and completion_path.is_file() and \
            response_path.is_file():
        try:
            previous_case = json.loads(path.read_text(encoding="utf-8"))
            completion = json.loads(completion_path.read_text(encoding="utf-8"))
            response = json.loads(response_path.read_text(encoding="utf-8"))
            if (previous_case == document and
                    completion["schema_version"] ==
                    "shieldSim.stage7_response_case_completion.v1" and
                    int(completion["events"]) == requested_histories and
                    int(response["histories"]) == requested_histories and
                    completion["score_ids"] == SELECTED_RESPONSE_SCORE_IDS and
                    completion["response_cell_sha256"] == file_sha256(response_path)):
                return {
                    "events": requested_histories,
                    "raw_transport_scores": [{"score_id": score_id}
                                             for score_id in SELECTED_RESPONSE_SCORE_IDS],
                }, prefix
        except (OSError, ValueError, KeyError, json.JSONDecodeError):
            pass

    # Migrate complete pre-compaction response artifacts in place. The C++
    # response parser validates the full cell before the large result is
    # discarded. Any incomplete legacy artifact falls through and is rerun.
    if path.is_file() and result_path.is_file():
        try:
            previous_case = json.loads(path.read_text(encoding="utf-8"))
            previous_result = json.loads(result_path.read_text(encoding="utf-8"))
            if (previous_case == document and
                    int(previous_result["events"]) == requested_histories and
                    (not response_required or response_path.is_file())):
                if response_required:
                    return compact_response_case(
                        directory, response_path, result_path, previous_result), prefix
                return previous_result, prefix
        except (AssertionError, OSError, ValueError, KeyError, json.JSONDecodeError):
            pass

    canonical_write(path, document)
    command([str(EXE), f"--case={path}"], directory, "shieldSim.log")
    require(result_path.is_file(), f"missing result {result_path}")
    result = json.loads(result_path.read_text(encoding="utf-8"))
    require(int(result["events"]) == requested_histories,
            f"case event count differs from fixed history count: {result_path}")
    if response_required:
        require(response_path.is_file(), f"missing response cell {response_path}")
        return compact_response_case(directory, response_path, result_path, result), prefix
    return result, prefix


def run_case_batch(
        jobs: list[tuple[dict, pathlib.Path, str]], workers: int,
        label: str) -> list[tuple[dict, pathlib.Path]]:
    """Run independent serial cases concurrently and preserve input ordering.

    ThreadPoolExecutor is used only as a subprocess launcher.  Python does not
    execute transport work in these threads; each job starts an independent
    shieldSim process in a unique directory.  Returning results in job order
    keeps pooling, evidence, and failure messages deterministic regardless of
    process completion order.
    """
    if not jobs:
        return []
    directories = [str(directory.resolve()) for _, directory, _ in jobs]
    require(len(directories) == len(set(directories)),
            f"{label}: parallel jobs must have unique working directories")
    active_workers = min(workers, len(jobs))
    print(f"Stage-7 {label}: {len(jobs)} independent cases, "
          f"{active_workers} worker(s)", flush=True)
    if active_workers == 1:
        return [run_case(case, directory, name)
                for case, directory, name in jobs]

    ordered: list[tuple[dict, pathlib.Path] | None] = [None] * len(jobs)
    with concurrent.futures.ThreadPoolExecutor(
            max_workers=active_workers,
            thread_name_prefix="stage7-case") as executor:
        futures = {
            executor.submit(run_case, case, directory, name): index
            for index, (case, directory, name) in enumerate(jobs)
        }
        for future in concurrent.futures.as_completed(futures):
            index = futures[future]
            # Calling result propagates the original command/assertion error.
            # The executor then waits for already-running processes, preventing
            # orphaned Geant4 jobs from writing after CTest reports failure.
            ordered[index] = future.result()
    require(all(item is not None for item in ordered),
            f"{label}: parallel result ledger is incomplete")
    return [item for item in ordered if item is not None]


def scores(document: dict, field: str = "raw_transport_scores") -> dict[str, dict]:
    return {item["score_id"]: item for item in document[field]}


def pooled_moment(documents: list[dict], score_id: str) -> dict[str, float | str]:
    histories = sum(float(document["events"]) for document in documents)
    total = 0.0
    squares = 0.0
    unit = ""
    for document in documents:
        row = scores(document).get(score_id)
        if row is not None:
            total += float(row["sum"])
            squares += float(row["sum_squares"])
            unit = row["unit"]
    mean = total / histories
    variance = max(0.0, (squares - total * total / histories) / (histories - 1.0))
    return {"mean": mean, "standard_error": math.sqrt(variance / histories),
            "histories": histories, "unit": unit}


def replicated_direct_moment(documents: list[dict], score_id: str) -> dict:
    """Combine equal-history direct seeds with conservative replicate error.

    The direct TID estimator samples the same rare nuclear/secondary tail as
    each response cell.  Treating the five seed files as one undifferentiated
    event stream can report less uncertainty than the observed seed-to-seed
    dispersion.  Use the same registered conservative rule as the response
    ensemble: the larger of propagated within-seed and empirical between-seed
    standard error.  This changes neither the mean nor either acceptance limit.
    """
    require(len(documents) >= 2, "direct replicate ensemble requires at least two members")
    rows = [scores(document)[score_id] for document in documents]
    histories = [int(document["events"]) for document in documents]
    require(len(set(histories)) == 1,
            "direct replicate uncertainty requires equal fixed history counts")
    values = [float(row["mean_per_primary"]) for row in rows]
    errors = [float(row["standard_error"]) for row in rows]
    count = len(rows)
    mean = sum(values) / count
    within = math.sqrt(sum(error * error for error in errors)) / count
    between = math.sqrt(sum((value - mean) ** 2 for value in values) /
                        (count * (count - 1)))
    pooled = pooled_moment(documents, score_id)
    require(math.isclose(mean, float(pooled["mean"]), rel_tol=1.0e-14,
                         abs_tol=1.0e-15),
            "equal-history direct replicate mean differs from pooled mean")
    return {
        "mean": mean,
        "standard_error": max(within, between),
        "within_replicate_standard_error": within,
        "between_replicate_standard_error": between,
        "pooled_event_standard_error": pooled["standard_error"],
        "replicate_values": values,
        "replicate_count": count,
        "histories": sum(histories),
        "unit": rows[0]["unit"],
    }


def equivalence(first: dict, second: dict, label: str) -> dict:
    delta = float(second["mean"]) - float(first["mean"])
    se = math.hypot(float(first["standard_error"]), float(second["standard_error"]))
    scale = max(abs(float(first["mean"])), abs(float(second["mean"])), 1.0e-300)
    relative = abs(delta) / scale
    z = abs(delta) / se if se > 0.0 else (0.0 if delta == 0.0 else math.inf)
    require(relative <= 0.10, f"{label}: relative difference {relative:.3%} exceeds 10%")
    require(z <= 1.96, f"{label}: difference z={z:.3f} exceeds two-sided 95% criterion")
    return {"id": label, "relative_difference": relative, "difference_z": z,
            "relative_limit": 0.10, "z_limit": 1.96}


def replicated_fold_moment(folds: list[dict], score_id: str) -> dict:
    """Combine equal-budget, independent folded response replicates.

    Response TID has a rare nuclear/secondary tail for which a single stream's
    event-level variance can be optimistic.  The replicate mean is unbiased
    because every replicate uses the same fixed allocation and an independent
    seed family.  The reported standard error is deliberately conservative:
    it is the larger of propagated within-replicate matrix uncertainty and the
    empirical standard error of the replicate means.
    """
    rows = [{item["score_id"]: item for item in fold["results"]}[score_id]
            for fold in folds]
    values = [float(row["value"]) for row in rows]
    errors = [float(row["matrix_statistical_standard_uncertainty"])
              for row in rows]
    count = len(rows)
    require(count >= 2, "response replicate ensemble requires at least two members")
    mean = sum(values) / count
    within = math.sqrt(sum(error * error for error in errors)) / count
    between = math.sqrt(sum((value - mean) ** 2 for value in values) /
                        (count * (count - 1)))
    return {
        "mean": mean,
        "standard_error": max(within, between),
        "within_replicate_standard_error": within,
        "between_replicate_standard_error": between,
        "replicate_values": values,
        "replicate_count": count,
        "unit": rows[0]["unit"],
    }


def source_nodes(component: dict) -> tuple[list[float], list[float], list[dict]]:
    energies = [0.5 * (float(low) + float(high))
                for low, high in zip(component["bin_lower"], component["bin_upper"])]
    unnormalized = [float(value) * (float(high) - float(low))
                    for value, low, high in zip(component["value"], component["bin_lower"],
                                                component["bin_upper"])]
    total = sum(unnormalized)
    energy_weights = [value / total for value in unnormalized]
    nodes = []
    for energy_index, (energy, energy_weight) in enumerate(zip(energies, energy_weights)):
        for mu_index, (mu, mu_weight) in enumerate(zip(MU_NODES, MU_WEIGHTS)):
            for phi_index, (phi, phi_weight) in enumerate(zip(PHI_NODES, PHI_WEIGHTS)):
                nodes.append({
                    "particle": "proton", "Z": 1, "A": 1,
                    "total_kinetic_energy_MeV": energy, "mu": mu, "phi_rad": phi,
                    "probability": energy_weight * mu_weight * phi_weight,
                    "incident_id": f"source_e{energy_index}_m{mu_index}_p{phi_index}",
                })
    probability_sum = sum(item["probability"] for item in nodes)
    nodes[-1]["probability"] += 1.0 - probability_sum
    require(abs(sum(item["probability"] for item in nodes) - 1.0) < 1.0e-14,
            "source quadrature probabilities do not close")
    return energies, energy_weights, nodes


def direct_quadrature_source(nodes: list[dict]) -> dict:
    return {
        "mode": "directional_plane",
        "spectrum": {
            "format": "discrete_response_nodes",
            "nodes": [
                {
                    "incident_id": item["incident_id"],
                    "particle": {"name": "proton", "Z": 1, "A": 1},
                    "energy": {"coordinate": "total_kinetic_energy",
                               "value": item["total_kinetic_energy_MeV"], "unit": "MeV"},
                    "direction": {"mu": item["mu"],
                                  "phi": {"value": item["phi_rad"], "unit": "rad"}},
                    "weight": item["probability"],
                }
                for item in nodes
            ],
        },
    }


def output_by_name(result: dict, name: str) -> dict:
    return {item["name"]: item for item in result["outputs"]}[name]


def main() -> int:
    workers = configured_worker_count()
    identity = campaign_cache_identity()
    marker = RUN_DIR / "campaign_cache_identity.json"
    adopt_existing = os.environ.get("SHIELDSIM_STAGE7_ADOPT_EXISTING", "") == "1"
    preserve = False
    if marker.is_file():
        try:
            preserve = json.loads(marker.read_text(encoding="utf-8")) == identity
        except (OSError, ValueError, json.JSONDecodeError):
            preserve = False
    # Adoption exists solely to recover artifacts produced immediately before
    # hash-guarded restart support was introduced. Normal CTest runs never set
    # it; without an identical marker they remove the entire run directory.
    if adopt_existing and RUN_DIR.is_dir():
        preserve = True
    if RUN_DIR.exists() and not preserve:
        shutil.rmtree(RUN_DIR)
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    canonical_write(marker, identity)
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    plan = json.loads(PLAN_PATH.read_text(encoding="utf-8"))
    require(plan["acceptance"]["independent_seed_count"] == len(SEEDS), "seed contract changed")
    execution = plan["execution"]
    require(execution["direct_seeds"] == SEEDS and
            execution["direct_histories_per_seed"] == DIRECT_HISTORIES,
            "direct ensemble execution contract changed")
    require(execution["response_pilot_histories_per_cell"] == PILOT_CELL_HISTORIES and
            execution["response_production_replicates"] == len(RESPONSE_SEED_BASES) and
            execution["response_production_histories_per_replicate"] ==
            RESPONSE_PRODUCTION_HISTORIES and
            execution["response_production_total_histories_per_stack"] ==
            len(RESPONSE_SEED_BASES) * RESPONSE_PRODUCTION_HISTORIES and
            execution["response_production_minimum_histories_per_cell"] ==
            RESPONSE_MINIMUM_PER_CELL and
            execution["response_production_seed_bases"] == RESPONSE_SEED_BASES,
            "two-phase response allocation contract changed")
    require(execution["response_refinement_total_histories_per_replicate"] ==
            REFINEMENT_PRODUCTION_HISTORIES and
            execution["response_refinement_minimum_histories_per_cell"] ==
            REFINEMENT_MINIMUM_PER_CELL and
            execution["response_refinement_energy_indices"] ==
            REFINEMENT_ENERGY_INDICES,
            "response refinement allocation contract changed")
    require(execution["response_exact_refinement_total_histories_per_replicate"] ==
            EXACT_REFINEMENT_PRODUCTION_HISTORIES and
            execution["response_exact_refinement_minimum_histories_per_cell"] ==
            EXACT_REFINEMENT_MINIMUM_PER_CELL and
            execution["response_exact_refinement_energy_indices"] ==
            EXACT_REFINEMENT_ENERGY_INDICES,
            "exact-node response refinement allocation contract changed")
    require(plan["acceptance"]["response_replicate_relative_standard_error_max"] == 0.10,
            "response replicate precision gate changed")
    require(plan["acceptance"]["direct_replicate_relative_standard_error_max"] == 0.10 and
            plan["acceptance"]["direct_replicate_uncertainty_rule"] ==
            "max(propagated_within_replicate_standard_error, "
            "empirical_between_replicate_standard_error)",
            "direct replicate uncertainty contract changed")
    require(execution["base_published_energy_indices"] == BASE_ENERGY_INDICES and
            execution["coarse_published_energy_indices"] == COARSE_ENERGY_INDICES and
            execution["fine_published_energy_indices"] == FINE_ENERGY_INDICES and
            execution["exact_published_energy_indices"] == EXACT_ENERGY_INDICES and
            execution["mu_nodes"] == MU_NODES and
            execution["phi_nodes_rad"] == PHI_NODES,
            "response quadrature contract changed")

    # Canonical published products are regenerated offline from committed,
    # checksum-pinned inputs. This does not contact an upstream service.
    ams_path = RUN_DIR / "ams02.json"
    command([
        str(SOURCE_TOOL), "ams02",
        f"--info={FIXTURES / 'gcr/ams02_20131118_hapi_info.json'}",
        f"--data={FIXTURES / 'gcr/ams02_20131118_hapi_data.csv'}",
        f"--interpretation={FIXTURES / 'gcr/ams02_20131118_interpretation.v1.json'}",
        f"--output={ams_path}",
    ], RUN_DIR, "ams02_adapter.log")
    ams = json.loads(ams_path.read_text(encoding="utf-8"))
    component = {item["component_id"]: item for item in ams["components"]}["proton"]
    energies, _, quadrature_nodes = source_nodes(component)
    require(len(energies) == 45 and component["mass_number"] == 1,
            "AMS A=1 proton component contract")

    goes_path = RUN_DIR / "goes16.json"
    command([
        str(SOURCE_TOOL), "goes16",
        f"--csv={FIXTURES / 'sep/goes16_sgps_gle72_timeseries.csv'}",
        f"--definition={FIXTURES / 'sep/goes16_gle72_event_definition.v1.json'}",
        f"--output={goes_path}",
    ], RUN_DIR, "goes16_adapter.log")
    goes = json.loads(goes_path.read_text(encoding="utf-8"))
    directions = {channel["spacecraft_axis"] for window in goes["windows"]
                  for channel in window["channels"]}
    require(directions == {"+X", "-X"}, "GOES sensor directions were not retained")
    require(any(not channel["complete"] for window in goes["windows"]
                for channel in window["channels"]), "GOES gaps were hidden")

    comparisons: list[dict] = []
    acceptance: list[dict] = []
    matrix_rows: list[dict] = []
    event_allocations: list[dict] = []
    published_results: dict[tuple[str, str], list[dict]] = {}

    # Every shield x absorber row is executed with the complete 45-bin proton
    # fixture and a five-seed ensemble. Only already-validated TID is requested.
    published_keys: list[tuple[str, str, int]] = []
    published_jobs: list[tuple[dict, pathlib.Path, str]] = []
    for stack_id in ("al_then_hdpe", "hdpe_then_al"):
        for absorber in ("silicon_slab", "water_slab"):
            for seed_index, seed in enumerate(SEEDS):
                case = base_case(stack_id, absorber, seed, DIRECT_HISTORIES)
                case["source"] = published_source()
                published_keys.append((stack_id, absorber, seed_index))
                published_jobs.append((
                    case,
                    RUN_DIR / "published" / stack_id / absorber /
                    f"seed{seed_index + 1}",
                    f"stage7_published_{stack_id}_{absorber}_seed{seed_index + 1}",
                ))
    for (stack_id, absorber, _), (result, _) in zip(
            published_keys,
            run_case_batch(published_jobs, workers, "published-source rows")):
        require(result["source_component_id"] == "proton", "published component changed")
        require(result["source_id"] == ams["dataset_id"], "published dataset identity changed")
        tid = output_by_name(result, "TID")
        require(tid["value"] > 0.0 and tid["source_normalized_value"] > 0.0,
                f"{stack_id}/{absorber} produced no TID")
        published_results.setdefault((stack_id, absorber), []).append(result)

    for stack_id in ("al_then_hdpe", "hdpe_then_al"):
        for absorber in ("silicon_slab", "water_slab"):
            matrix_rows.append({
                "shield_stack_id": stack_id, "absorber_id": absorber,
                "histories": len(SEEDS) * DIRECT_HISTORIES,
                "source_id": ams["dataset_id"],
                "environment_completeness": "incomplete_not_full_gcr",
                "requested_quantities": ["TID"], "status": "executed",
            })

    # Response-matrix comparisons are made for the silicon rows. The exact
    # same 45 energy-bin representatives and angle weights drive direct and
    # folded paths; only response interpolation differs.
    score_id = "absorber/silicon_slab/ionizing_deposited_energy_MeV"
    grid_score_id = "shield_exit/primary/proton/crossing_count"
    for stack_index, stack_id in enumerate(("al_then_hdpe", "hdpe_then_al")):
        pilot_cell_paths: dict[tuple[int, int, int], pathlib.Path] = {}
        pilot_keys: list[tuple[int, int, int]] = []
        pilot_jobs: list[tuple[dict, pathlib.Path, str]] = []
        for energy_index in BASE_ENERGY_INDICES:
            for mu_index, mu in enumerate(MU_NODES):
                for phi_index, phi in enumerate(PHI_NODES):
                    incident_id = f"matrix_e{energy_index}_m{mu_index}_p{phi_index}"
                    seed = 720001 + stack_index * 10000 + energy_index * 100 + mu_index * 10 + phi_index
                    case = base_case(stack_id, "silicon_slab", seed, PILOT_CELL_HISTORIES)
                    set_response_incident(case, incident_id, energies[energy_index], mu, phi)
                    pilot_keys.append((energy_index, mu_index, phi_index))
                    pilot_jobs.append((
                        case, RUN_DIR / "response" / stack_id / incident_id,
                        f"stage7_{stack_id}_{incident_id}",
                    ))
        for key, (result, prefix) in zip(
                pilot_keys, run_case_batch(pilot_jobs, workers,
                                           f"{stack_id} base-grid pilots")):
            incident_id = f"matrix_e{key[0]}_m{key[1]}_p{key[2]}"
            require(score_id in scores(result), f"{incident_id} missing silicon score")
            cell_path = pathlib.Path(str(prefix) + "_response_cell.json")
            require(cell_path.is_file(), f"missing response cell {cell_path}")
            pilot_cell_paths[key] = cell_path

        source_document = {
            "schema_version": "shieldSim.response_source_nodes.v1",
            "normalization": "unit probability from AMS proton J_i*DeltaE_i and angular quadrature",
            "quadrature_policy": "published exact edges; arithmetic energy representatives; "
                                 "three-point p(mu)=2mu; three-point Gauss phi on the "
                                 "square-symmetry interval [0,pi/4]",
            "nodes": [{key: value for key, value in node.items() if key != "incident_id"}
                      for node in quadrature_nodes],
        }
        source_path = RUN_DIR / "response" / stack_id / "source_nodes.json"
        canonical_write(source_path, source_document)

        # Build weights against the pilot grid, then invoke the production
        # allocation tool with each cell's source weight and measured pilot
        # standard deviation. This is the real Stage-7 two-phase allocation
        # path, not merely a synthetic unit test of the integer allocator.
        directory = RUN_DIR / "response" / stack_id
        pilot_matrix_path = directory / "matrix_pilot.json"
        pilot_selected = [str(pilot_cell_paths[(energy_index, mu_index, phi_index)])
                          for energy_index in BASE_ENERGY_INDICES
                          for mu_index in range(len(MU_NODES))
                          for phi_index in range(len(PHI_NODES))]
        command([str(RESPONSE_TOOL), "assemble-selected", str(pilot_matrix_path),
                 f"{score_id},{grid_score_id}", *pilot_selected],
                directory, "assemble_pilot.log")
        pilot_weights_path = directory / "weights_pilot.json"
        command([str(PYTHON), str(WEIGHT_TOOL), str(pilot_matrix_path), str(source_path),
                 str(pilot_weights_path)], directory, "weights_pilot.log")
        pilot_weights = json.loads(pilot_weights_path.read_text(encoding="utf-8"))
        source_weight = {item["incident_id"]: float(item["weight"])
                         for item in pilot_weights["weights"]}
        strata = []
        for key in sorted(pilot_cell_paths):
            cell = json.loads(pilot_cell_paths[key].read_text(encoding="utf-8"))
            incident_id = cell["incident"]["incident_id"]
            pilot_score = next(item for item in cell["scores"]
                               if item["score_id"] == score_id)
            pilot_deviation = float(pilot_score["standard_error"]) * math.sqrt(
                float(pilot_score["histories"])
            )
            strata.append({"stratum_id": incident_id,
                            "source_weight": source_weight.get(incident_id, 0.0),
                            "pilot_standard_deviation": pilot_deviation})
        allocation_input = {
            "schema_version": "shieldSim.neyman_allocation.v1",
            "policy": "pilot_neyman_no_transport_biasing",
            "total_histories": RESPONSE_PRODUCTION_HISTORIES,
            "minimum_per_stratum": RESPONSE_MINIMUM_PER_CELL,
            "strata": strata,
        }
        allocation_input_path = directory / "allocation_input.json"
        allocation_output_path = directory / "allocation_output.json"
        canonical_write(allocation_input_path, allocation_input)
        command([str(EXPOSURE_TOOL), "allocate", str(allocation_input_path),
                 str(allocation_output_path)], directory, "allocation.log")
        allocation = json.loads(allocation_output_path.read_text(encoding="utf-8"))
        histories_by_id = {item["stratum_id"]: int(item["histories"])
                           for item in allocation["allocation"]}
        require(sum(histories_by_id.values()) == RESPONSE_PRODUCTION_HISTORIES,
                f"{stack_id}: allocated history total changed")
        require(min(histories_by_id.values()) >= RESPONSE_MINIMUM_PER_CELL,
                f"{stack_id}: allocated history floor was violated")
        require(allocation["physics_biasing"] is False,
                f"{stack_id}: allocation enabled physics biasing")
        event_allocations.append({
            "allocation_phase": "base_grid",
            "shield_stack_id": stack_id,
            "pilot_histories_per_cell": PILOT_CELL_HISTORIES,
            "production_replicates": len(RESPONSE_SEED_BASES),
            "production_histories_per_replicate": RESPONSE_PRODUCTION_HISTORIES,
            "production_total_histories": (len(RESPONSE_SEED_BASES) *
                                           RESPONSE_PRODUCTION_HISTORIES),
            "minimum_histories_per_cell": RESPONSE_MINIMUM_PER_CELL,
            "production_seed_bases": RESPONSE_SEED_BASES,
            "physics_biasing": allocation["physics_biasing"],
            "allocation_path": str(allocation_output_path),
        })

        # The first completed qualification exposed a statistically resolved
        # response-grid bias for the reversed stack. Add the interleaved energy
        # nodes with a separate pilot/allocation so the already observed base
        # histories remain frozen. Exact published-node probabilities are safe
        # allocation weights here: allocation changes variance only and never
        # changes a primary's Geant4 weight or interaction physics.
        refinement_pilot_paths: dict[tuple[int, int, int], pathlib.Path] = {}
        node_probabilities = {item["incident_id"]: float(item["probability"])
                              for item in quadrature_nodes}
        refinement_pilot_keys: list[tuple[int, int, int]] = []
        refinement_pilot_jobs: list[tuple[dict, pathlib.Path, str]] = []
        for energy_index in REFINEMENT_ENERGY_INDICES:
            for mu_index, mu in enumerate(MU_NODES):
                for phi_index, phi in enumerate(PHI_NODES):
                    incident_id = f"matrix_e{energy_index}_m{mu_index}_p{phi_index}"
                    seed = (720001 + stack_index * 10000 + energy_index * 100 +
                            mu_index * 10 + phi_index)
                    case = base_case(stack_id, "silicon_slab", seed,
                                     PILOT_CELL_HISTORIES)
                    set_response_incident(case, incident_id, energies[energy_index], mu, phi)
                    refinement_pilot_keys.append((energy_index, mu_index, phi_index))
                    refinement_pilot_jobs.append((
                        case, RUN_DIR / "response_refinement" / stack_id / "pilot" /
                        incident_id,
                        f"stage7_refinement_pilot_{stack_id}_{incident_id}",
                    ))
        for key, (result, prefix) in zip(
                refinement_pilot_keys,
                run_case_batch(refinement_pilot_jobs, workers,
                               f"{stack_id} 23-node-refinement pilots")):
            incident_id = f"matrix_e{key[0]}_m{key[1]}_p{key[2]}"
            require(score_id in scores(result),
                    f"refinement pilot {incident_id} missing silicon score")
            cell_path = pathlib.Path(str(prefix) + "_response_cell.json")
            require(cell_path.is_file(),
                    f"missing refinement pilot response cell {cell_path}")
            refinement_pilot_paths[key] = cell_path

        refinement_strata = []
        for (energy_index, mu_index, phi_index), cell_path in sorted(
                refinement_pilot_paths.items()):
            cell = json.loads(cell_path.read_text(encoding="utf-8"))
            incident_id = cell["incident"]["incident_id"]
            pilot_score = next(item for item in cell["scores"]
                               if item["score_id"] == score_id)
            pilot_deviation = float(pilot_score["standard_error"]) * math.sqrt(
                float(pilot_score["histories"]))
            published_node_id = f"source_e{energy_index}_m{mu_index}_p{phi_index}"
            refinement_strata.append({
                "stratum_id": incident_id,
                "source_weight": node_probabilities[published_node_id],
                "pilot_standard_deviation": pilot_deviation,
            })
        refinement_allocation_input = {
            "schema_version": "shieldSim.neyman_allocation.v1",
            "policy": "pilot_neyman_no_transport_biasing",
            "total_histories": REFINEMENT_PRODUCTION_HISTORIES,
            "minimum_per_stratum": REFINEMENT_MINIMUM_PER_CELL,
            "strata": refinement_strata,
        }
        refinement_directory = RUN_DIR / "response_refinement" / stack_id
        refinement_input_path = refinement_directory / "allocation_input.json"
        refinement_output_path = refinement_directory / "allocation_output.json"
        canonical_write(refinement_input_path, refinement_allocation_input)
        command([str(EXPOSURE_TOOL), "allocate", str(refinement_input_path),
                 str(refinement_output_path)], refinement_directory, "allocation.log")
        refinement_allocation = json.loads(
            refinement_output_path.read_text(encoding="utf-8"))
        refinement_histories_by_id = {
            item["stratum_id"]: int(item["histories"])
            for item in refinement_allocation["allocation"]
        }
        require(sum(refinement_histories_by_id.values()) ==
                REFINEMENT_PRODUCTION_HISTORIES,
                f"{stack_id}: refinement allocated history total changed")
        require(min(refinement_histories_by_id.values()) >=
                REFINEMENT_MINIMUM_PER_CELL,
                f"{stack_id}: refinement allocated history floor was violated")
        require(refinement_allocation["physics_biasing"] is False,
                f"{stack_id}: refinement allocation enabled physics biasing")
        event_allocations.append({
            "allocation_phase": "interleaved_energy_refinement",
            "shield_stack_id": stack_id,
            "pilot_histories_per_cell": PILOT_CELL_HISTORIES,
            "production_replicates": len(RESPONSE_SEED_BASES),
            "production_histories_per_replicate":
                REFINEMENT_PRODUCTION_HISTORIES,
            "production_total_histories": (len(RESPONSE_SEED_BASES) *
                                           REFINEMENT_PRODUCTION_HISTORIES),
            "minimum_histories_per_cell": REFINEMENT_MINIMUM_PER_CELL,
            "production_seed_bases": RESPONSE_SEED_BASES,
            "physics_biasing": refinement_allocation["physics_biasing"],
            "allocation_path": str(refinement_output_path),
        })

        # The 23-node campaign still missed the unchanged two-sided 95% gate
        # for the reversed stack (z=2.078).  Complete the registered corrective
        # action with the remaining odd published-bin representatives.  This
        # second supplement has its own pilot and frozen allocation, so neither
        # the base cells nor the first refinement are resampled or re-budgeted
        # after observing their outcomes.
        exact_pilot_paths: dict[tuple[int, int, int], pathlib.Path] = {}
        exact_pilot_keys: list[tuple[int, int, int]] = []
        exact_pilot_jobs: list[tuple[dict, pathlib.Path, str]] = []
        for energy_index in EXACT_REFINEMENT_ENERGY_INDICES:
            for mu_index, mu in enumerate(MU_NODES):
                for phi_index, phi in enumerate(PHI_NODES):
                    incident_id = f"matrix_e{energy_index}_m{mu_index}_p{phi_index}"
                    seed = (720001 + stack_index * 10000 + energy_index * 100 +
                            mu_index * 10 + phi_index)
                    case = base_case(stack_id, "silicon_slab", seed,
                                     PILOT_CELL_HISTORIES)
                    set_response_incident(case, incident_id, energies[energy_index], mu, phi)
                    exact_pilot_keys.append((energy_index, mu_index, phi_index))
                    exact_pilot_jobs.append((
                        case, RUN_DIR / "response_exact_refinement" / stack_id /
                        "pilot" / incident_id,
                        f"stage7_exact_pilot_{stack_id}_{incident_id}",
                    ))
        for key, (result, prefix) in zip(
                exact_pilot_keys,
                run_case_batch(exact_pilot_jobs, workers,
                               f"{stack_id} 45-node-refinement pilots")):
            incident_id = f"matrix_e{key[0]}_m{key[1]}_p{key[2]}"
            require(score_id in scores(result),
                    f"exact-refinement pilot {incident_id} missing silicon score")
            cell_path = pathlib.Path(str(prefix) + "_response_cell.json")
            require(cell_path.is_file(),
                    f"missing exact-refinement pilot response cell {cell_path}")
            exact_pilot_paths[key] = cell_path

        exact_strata = []
        for (energy_index, mu_index, phi_index), cell_path in sorted(
                exact_pilot_paths.items()):
            cell = json.loads(cell_path.read_text(encoding="utf-8"))
            incident_id = cell["incident"]["incident_id"]
            pilot_score = next(item for item in cell["scores"]
                               if item["score_id"] == score_id)
            pilot_deviation = float(pilot_score["standard_error"]) * math.sqrt(
                float(pilot_score["histories"]))
            published_node_id = f"source_e{energy_index}_m{mu_index}_p{phi_index}"
            exact_strata.append({
                "stratum_id": incident_id,
                "source_weight": node_probabilities[published_node_id],
                "pilot_standard_deviation": pilot_deviation,
            })
        exact_allocation_input = {
            "schema_version": "shieldSim.neyman_allocation.v1",
            "policy": "pilot_neyman_no_transport_biasing",
            "total_histories": EXACT_REFINEMENT_PRODUCTION_HISTORIES,
            "minimum_per_stratum": EXACT_REFINEMENT_MINIMUM_PER_CELL,
            "strata": exact_strata,
        }
        exact_directory = RUN_DIR / "response_exact_refinement" / stack_id
        exact_input_path = exact_directory / "allocation_input.json"
        exact_output_path = exact_directory / "allocation_output.json"
        canonical_write(exact_input_path, exact_allocation_input)
        command([str(EXPOSURE_TOOL), "allocate", str(exact_input_path),
                 str(exact_output_path)], exact_directory, "allocation.log")
        exact_allocation = json.loads(exact_output_path.read_text(encoding="utf-8"))
        exact_histories_by_id = {
            item["stratum_id"]: int(item["histories"])
            for item in exact_allocation["allocation"]
        }
        require(sum(exact_histories_by_id.values()) ==
                EXACT_REFINEMENT_PRODUCTION_HISTORIES,
                f"{stack_id}: exact refinement allocated history total changed")
        require(min(exact_histories_by_id.values()) >=
                EXACT_REFINEMENT_MINIMUM_PER_CELL,
                f"{stack_id}: exact refinement history floor was violated")
        require(exact_allocation["physics_biasing"] is False,
                f"{stack_id}: exact refinement enabled physics biasing")
        event_allocations.append({
            "allocation_phase": "all_published_energy_nodes_refinement",
            "shield_stack_id": stack_id,
            "pilot_histories_per_cell": PILOT_CELL_HISTORIES,
            "production_replicates": len(RESPONSE_SEED_BASES),
            "production_histories_per_replicate":
                EXACT_REFINEMENT_PRODUCTION_HISTORIES,
            "production_total_histories": (len(RESPONSE_SEED_BASES) *
                                           EXACT_REFINEMENT_PRODUCTION_HISTORIES),
            "minimum_histories_per_cell": EXACT_REFINEMENT_MINIMUM_PER_CELL,
            "production_seed_bases": RESPONSE_SEED_BASES,
            "physics_biasing": exact_allocation["physics_biasing"],
            "allocation_path": str(exact_output_path),
        })

        # Five equal-budget production matrices use independent seed families
        # and exclude pilot moments. Replication is essential for this TID
        # score: it exposes rare-tail variance that one event stream can miss.
        # All replicates retain the same pilot-derived allocation, so their
        # arithmetic mean estimates the same response without allocation bias.
        folds_by_label: dict[str, list[dict]] = {
            "coarse": [], "fine": [], "exact": []
        }
        for replicate_index, seed_base in enumerate(RESPONSE_SEED_BASES, start=1):
            cell_paths: dict[tuple[int, int, int], pathlib.Path] = {}
            base_keys: list[tuple[int, int, int]] = []
            base_jobs: list[tuple[dict, pathlib.Path, str]] = []
            for energy_index in BASE_ENERGY_INDICES:
                for mu_index, mu in enumerate(MU_NODES):
                    for phi_index, phi in enumerate(PHI_NODES):
                        incident_id = f"matrix_e{energy_index}_m{mu_index}_p{phi_index}"
                        histories = histories_by_id[incident_id]
                        seed = (seed_base + stack_index * 10000 + energy_index * 100 +
                                mu_index * 10 + phi_index)
                        case = base_case(stack_id, "silicon_slab", seed, histories)
                        set_response_incident(case, incident_id, energies[energy_index],
                                              mu, phi)
                        base_keys.append((energy_index, mu_index, phi_index))
                        base_jobs.append((
                            case,
                            RUN_DIR / "response_production" / stack_id /
                            f"replicate{replicate_index}" / incident_id,
                            f"stage7_production_{stack_id}_replicate{replicate_index}_"
                            f"{incident_id}",
                        ))
            for key, (result, prefix) in zip(
                    base_keys,
                    run_case_batch(base_jobs, workers,
                                   f"{stack_id} replicate {replicate_index} base cells")):
                incident_id = f"matrix_e{key[0]}_m{key[1]}_p{key[2]}"
                require(score_id in scores(result),
                        f"production {incident_id} missing silicon score")
                cell_path = pathlib.Path(str(prefix) + "_response_cell.json")
                require(cell_path.is_file(),
                        f"missing production response cell {cell_path}")
                cell_paths[key] = cell_path

            # Supplemental production cells use the same replicate seed base
            # but previously unused energy indices, so all random seeds remain
            # disjoint from the base grid and from every other replicate.
            refinement_keys: list[tuple[int, int, int]] = []
            refinement_jobs: list[tuple[dict, pathlib.Path, str]] = []
            for energy_index in REFINEMENT_ENERGY_INDICES:
                for mu_index, mu in enumerate(MU_NODES):
                    for phi_index, phi in enumerate(PHI_NODES):
                        incident_id = f"matrix_e{energy_index}_m{mu_index}_p{phi_index}"
                        histories = refinement_histories_by_id[incident_id]
                        seed = (seed_base + stack_index * 10000 + energy_index * 100 +
                                mu_index * 10 + phi_index)
                        case = base_case(stack_id, "silicon_slab", seed, histories)
                        set_response_incident(case, incident_id, energies[energy_index],
                                              mu, phi)
                        refinement_keys.append((energy_index, mu_index, phi_index))
                        refinement_jobs.append((
                            case, RUN_DIR / "response_refinement" / stack_id /
                            f"replicate{replicate_index}" / incident_id,
                            f"stage7_refinement_{stack_id}_replicate{replicate_index}_"
                            f"{incident_id}",
                        ))
            for key, (result, prefix) in zip(
                    refinement_keys,
                    run_case_batch(refinement_jobs, workers,
                                   f"{stack_id} replicate {replicate_index} "
                                   "23-node supplement")):
                incident_id = f"matrix_e{key[0]}_m{key[1]}_p{key[2]}"
                require(score_id in scores(result),
                        f"refinement {incident_id} missing silicon score")
                cell_path = pathlib.Path(str(prefix) + "_response_cell.json")
                require(cell_path.is_file(),
                        f"missing refinement response cell {cell_path}")
                cell_paths[key] = cell_path

            # Final supplement: all odd energy indices.  These jobs are the
            # expensive new work after the 23-node comparison; executing them
            # concurrently is safe because every directory and seed is unique.
            exact_keys: list[tuple[int, int, int]] = []
            exact_jobs: list[tuple[dict, pathlib.Path, str]] = []
            for energy_index in EXACT_REFINEMENT_ENERGY_INDICES:
                for mu_index, mu in enumerate(MU_NODES):
                    for phi_index, phi in enumerate(PHI_NODES):
                        incident_id = f"matrix_e{energy_index}_m{mu_index}_p{phi_index}"
                        histories = exact_histories_by_id[incident_id]
                        seed = (seed_base + stack_index * 10000 + energy_index * 100 +
                                mu_index * 10 + phi_index)
                        case = base_case(stack_id, "silicon_slab", seed, histories)
                        set_response_incident(case, incident_id, energies[energy_index],
                                              mu, phi)
                        exact_keys.append((energy_index, mu_index, phi_index))
                        exact_jobs.append((
                            case, RUN_DIR / "response_exact_refinement" / stack_id /
                            f"replicate{replicate_index}" / incident_id,
                            f"stage7_exact_{stack_id}_replicate{replicate_index}_"
                            f"{incident_id}",
                        ))
            for key, (result, prefix) in zip(
                    exact_keys,
                    run_case_batch(exact_jobs, workers,
                                   f"{stack_id} replicate {replicate_index} "
                                   "45-node supplement")):
                incident_id = f"matrix_e{key[0]}_m{key[1]}_p{key[2]}"
                require(score_id in scores(result),
                        f"exact refinement {incident_id} missing silicon score")
                cell_path = pathlib.Path(str(prefix) + "_response_cell.json")
                require(cell_path.is_file(),
                        f"missing exact-refinement response cell {cell_path}")
                cell_paths[key] = cell_path

            replicate_directory = directory / f"replicate{replicate_index}"
            replicate_directory.mkdir(parents=True, exist_ok=True)
            for label, indices in (("coarse", COARSE_ENERGY_INDICES),
                                   ("fine", FINE_ENERGY_INDICES),
                                   ("exact", EXACT_ENERGY_INDICES)):
                matrix_path = replicate_directory / f"matrix_{label}.json"
                selected = [str(cell_paths[(energy_index, mu_index, phi_index)])
                            for energy_index in indices
                            for mu_index in range(len(MU_NODES))
                            for phi_index in range(len(PHI_NODES))]
                # Keep only the registered TID observable and high-population
                # grid monitor. Full cell validation still occurs first.
                command([str(RESPONSE_TOOL), "assemble-selected", str(matrix_path),
                         f"{score_id},{grid_score_id}", *selected],
                        replicate_directory, f"assemble_{label}.log")
                weights_path = replicate_directory / f"weights_{label}.json"
                command([str(PYTHON), str(WEIGHT_TOOL), str(matrix_path), str(source_path),
                         str(weights_path)], replicate_directory, f"weights_{label}.log")
                fold_path = replicate_directory / f"fold_{label}.json"
                command([str(RESPONSE_TOOL), "fold", str(matrix_path), str(weights_path),
                         str(fold_path)], replicate_directory, f"fold_{label}.log")
                folds_by_label[label].append(
                    json.loads(fold_path.read_text(encoding="utf-8")))

        folded_moment = replicated_fold_moment(folds_by_label["exact"], score_id)
        exact_grid_moment = replicated_fold_moment(folds_by_label["exact"], grid_score_id)
        fine_grid_moment = replicated_fold_moment(folds_by_label["fine"], grid_score_id)
        coarse_grid_moment = replicated_fold_moment(folds_by_label["coarse"], grid_score_id)
        # Grid refinement needs a statistically resolved primitive. Silicon
        # ionizing deposition retains a rare nuclear tail and is therefore
        # judged by the direct-vs-matrix 10%/95% gate below, not misused as a
        # 2% deterministic interpolation oracle. Primary transmission is the
        # registered high-population grid monitor for this qualification.
        grid_change = abs(float(exact_grid_moment["mean"]) -
                          float(fine_grid_moment["mean"])) / max(
            abs(float(exact_grid_moment["mean"])), 1.0e-300
        )
        require(grid_change <= 0.02,
                f"{stack_id}: final response-grid change {grid_change:.3%} exceeds 2%")
        acceptance.append({"id": f"{stack_id}_response_grid_refinement",
                           "score_id": grid_score_id,
                           "relative_change": grid_change, "limit": 0.02})
        response_relative_standard_error = (float(folded_moment["standard_error"]) /
                                            max(abs(float(folded_moment["mean"])),
                                                1.0e-300))
        require(response_relative_standard_error <= 0.10,
                f"{stack_id}: replicated response relative standard error "
                f"{response_relative_standard_error:.3%} exceeds 10%")
        acceptance.append({
            "id": f"{stack_id}_response_replicate_precision",
            "replicate_count": folded_moment["replicate_count"],
            "replicate_values": folded_moment["replicate_values"],
            "within_replicate_standard_error":
                folded_moment["within_replicate_standard_error"],
            "between_replicate_standard_error":
                folded_moment["between_replicate_standard_error"],
            "used_standard_error": folded_moment["standard_error"],
            "relative_standard_error": response_relative_standard_error,
            "relative_limit": 0.10,
        })
        comparisons.extend([
            {"id": f"{stack_id}_coarse_vs_fine_response_grid",
             "expected": float(fine_grid_moment["mean"]),
             "actual": float(coarse_grid_moment["mean"]),
             "unit": fine_grid_moment["unit"]},
            {"id": f"{stack_id}_fine_vs_exact_response_grid",
             "expected": float(exact_grid_moment["mean"]),
             "actual": float(fine_grid_moment["mean"]),
             "unit": exact_grid_moment["unit"]},
        ])

        quadrature_jobs: list[tuple[dict, pathlib.Path, str]] = []
        for seed_index, seed in enumerate(SEEDS):
            case = base_case(stack_id, "silicon_slab", seed + 1000, DIRECT_HISTORIES)
            case["source"] = direct_quadrature_source(quadrature_nodes)
            quadrature_jobs.append((
                case, RUN_DIR / "quadrature_direct" / stack_id / f"seed{seed_index + 1}",
                f"stage7_quadrature_{stack_id}_seed{seed_index + 1}",
            ))
        quadrature_documents = [result for result, _ in run_case_batch(
            quadrature_jobs, workers, f"{stack_id} direct quadrature")]

        quadrature_moment = replicated_direct_moment(quadrature_documents, score_id)
        published_moment = replicated_direct_moment(
            published_results[(stack_id, "silicon_slab")], score_id)
        for source_path_id, moment in (("quadrature", quadrature_moment),
                                       ("published", published_moment)):
            relative_standard_error = (float(moment["standard_error"]) /
                                       max(abs(float(moment["mean"])), 1.0e-300))
            require(relative_standard_error <= 0.10,
                    f"{stack_id}: {source_path_id} direct replicate relative "
                    f"standard error {relative_standard_error:.3%} exceeds 10%")
            acceptance.append({
                "id": f"{stack_id}_{source_path_id}_direct_replicate_precision",
                "replicate_count": moment["replicate_count"],
                "replicate_values": moment["replicate_values"],
                "within_replicate_standard_error":
                    moment["within_replicate_standard_error"],
                "between_replicate_standard_error":
                    moment["between_replicate_standard_error"],
                "pooled_event_standard_error": moment["pooled_event_standard_error"],
                "used_standard_error": moment["standard_error"],
                "relative_standard_error": relative_standard_error,
                "relative_limit": 0.10,
            })
        acceptance.append(equivalence(folded_moment, quadrature_moment,
                                      f"{stack_id}_direct_vs_matrix"))
        acceptance.append(equivalence(published_moment, quadrature_moment,
                                      f"{stack_id}_published_vs_quadrature"))
        comparisons.extend([
            {"id": f"{stack_id}_direct_vs_matrix", "expected": folded_moment["mean"],
             "actual": quadrature_moment["mean"], "unit": folded_moment["unit"]},
            {"id": f"{stack_id}_published_vs_quadrature", "expected": published_moment["mean"],
             "actual": quadrature_moment["mean"], "unit": folded_moment["unit"]},
        ])

    # Convert each published matrix row into an explicit exposure over the
    # actual 27-day AMS averaging interval. Source covariance is unpublished,
    # so every result must remain not_assessable rather than gaining a diagonal
    # covariance assumption.
    exposure_rows = []
    for (stack_id, absorber), documents in published_results.items():
        rates = [float(output_by_name(document, "TID")["source_normalized_value"])
                 for document in documents]
        rate = sum(rates) / len(rates)
        input_document = {
            "schema_version": "shieldSim.time_resolved_exposure.v1",
            "case_id": f"stage7_{stack_id}_{absorber}_br2460_exposure",
            "time_system": "UTC",
            "source_qualification": {
                "environment": "ams02_a1_proton_component",
                "environment_completeness": "incomplete_not_full_gcr",
                "goes_transport_status": "disabled_missing_boundary_mapping",
                "transport_admission": "verification_test_only",
                "limitations": [
                    "AMS A=1 proton component is not a complete GCR environment.",
                    "Published source covariance is unavailable; uncertainty is incomplete.",
                ],
            },
            "integration_policy": {
                "gap_policy": "preserve_no_fill",
                "cadence_policy": "native_intervals_no_implicit_rebin",
                "configuration_policy": "split_before_fold",
                "extrapolation_policy": "refuse",
            },
            "anomaly_timestamps": [],
            "intervals": [{
                "interval_id": "BR2460_average_support",
                "start_utc": ams["observation_start_utc"],
                "end_utc": ams["observation_stop_utc"],
                "original_time_coordinate": "AMS Bartels-rotation average support",
                "source_id": ams["dataset_id"],
                "source_hash": documents[0]["source_hash"],
                "response_id": f"{stack_id}_{absorber}_direct_ensemble",
                "response_hash": documents[0].get("case_hash", "case-hash-not-exported"),
                "shield_stack_id": stack_id, "absorber_id": absorber,
                "attitude_id": "isotropic_simplification",
                "quality_state": "valid", "coverage_fraction": 1.0,
                "rates": [{"series_id": f"{absorber}/TID", "quantity": "TID",
                           "rate_unit": "Gy/s", "integrated_unit": "Gy", "value": rate}],
            }],
            "uncertainty": {
                "status": "incomplete",
                "limitation": "AMS total per-bin errors have no published covariance for this fold.",
                "independence_basis": "No unprovided source covariance is assumed.",
                "covariance_components": [], "sampled_ensembles": [],
            },
        }
        directory = RUN_DIR / "exposures" / stack_id / absorber
        directory.mkdir(parents=True, exist_ok=True)
        input_path = directory / "input.json"
        output_path = directory / "result.json"
        canonical_write(input_path, input_document)
        command([str(EXPOSURE_TOOL), "integrate", str(input_path), str(output_path)],
                directory, "exposure.log")
        result = json.loads(output_path.read_text(encoding="utf-8"))
        integrated = result["integrated_exposures"][0]
        require(integrated["accuracy_status"] == "not_assessable" and
                integrated["standard_uncertainty"] is None,
                "unpublished source covariance acquired an accuracy claim")
        expected = rate * 27.0 * 86400.0
        require(math.isclose(integrated["value"], expected, rel_tol=1.0e-14),
                "27-day exposure integration")
        exposure_rows.append({"shield_stack_id": stack_id, "absorber_id": absorber,
                              "rate_Gy_per_s": rate, "exposure_Gy": integrated["value"],
                              "uncertainty_status": integrated["uncertainty_status"]})

    evidence = {
        "schema_version": "shieldSim.validation_evidence.v1",
        "test_case": "stage7_published_qualification",
        "validation_status": "stage7_software_framework_and_component_qualification_implemented_"
                             "incomplete_physical_gcr_goes_qualification",
        "plan_sha256": hashlib.sha256(PLAN_PATH.read_bytes()).hexdigest(),
        "source_boundaries": {
            "gcr": "AMS A=1 proton component only; incomplete_not_full_gcr",
            "sep": "GOES directions retained; transport disabled_missing_boundary_mapping",
        },
        "independent_seed_count": len(SEEDS),
        "representative_matrix": matrix_rows,
        "comparisons": comparisons,
        "acceptance_details": acceptance,
        "event_allocations": event_allocations,
        "time_integrated_exposures": exposure_rows,
        "source_uncertainty_status": "incomplete_published_covariance_unavailable",
        "parallel_execution": {
            "mode": "independent_serial_shieldSim_processes",
            "workers": workers,
            "per_case_seed_and_history_contract": "unchanged",
            "shared_output_directories": False,
            "event_level_mpi_or_geant4_mt": False,
        },
        "run_directory": str(RUN_DIR),
    }
    canonical_write(EVIDENCE_DIR / "stage7_published_qualification.json", evidence)
    canonical_write(RUN_DIR / "acceptance_summary.json", evidence)
    print(f"PASS Stage-7 published component qualification; artifacts: {RUN_DIR}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AssertionError, OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        print(f"FAIL Stage-7 published qualification: {error}", file=sys.stderr)
        raise SystemExit(1)
