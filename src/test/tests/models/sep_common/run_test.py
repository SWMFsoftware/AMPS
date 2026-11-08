#!/usr/bin/env python3
"""Import every authoritative ``src/models/sep_common`` verification scope.

This adapter deliberately owns no physical expected value.  It builds and runs
the five model-local verification entry points in an isolated copy, validates
their native registries/reports, and only then translates their typed outcomes
to the global schema.  Isolation is important here: the perpendicular library
builds the parallel library as a dependency, while several reference scripts
write generated audit products beside their inputs.  Running in the checkout
would therefore make concurrent nightlies race and would mutate a read-mostly
source tree.

The pitch-angle Monte-Carlo acceptance campaign takes several hours.  Routine
T11 imports use the local runner's explicit ``--skip-acceptance`` mode, then
query the authoritative acceptance binary with ``--list`` and publish every
one of its 52 identities as SKIP.  This retains the complete qualification
boundary; it cannot be mistaken for an acceptance PASS or silently collapsed
to the runner's single aggregate ``AP`` record.
"""

from __future__ import annotations

import json
import re
import shlex
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple


TEST_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(TEST_ROOT))

from lib.runner_support import (  # noqa: E402
    RunnerError, case, command, error_result, load_context, main_guard,
    run_command, write_result)


PARALLEL_CPP_IDS = [
    "PD13-UNCONFIGURED", "PD13-SECTION-GRAMMAR", "PD13-SECTION-ERRORS",
    "PD01-UNITS", "PD01-INVALID-PARTICLE", "PD02-CONSTANT-LAMBDA",
    "PD02-CONSTANT-KAPPA", "PD02-POWER-REFERENCE", "PD02-POWER-SCALING",
    "PD02-POWER-SPECIES", "PD02-POWER-VARIANTS",
    "PD02-POWER-DISABLED-FIELD", "PD02-POWER-OPTIONAL-INPUT",
    "PD02-BROKEN-REFERENCE", "PD02-BROKEN-SINGLE",
    "PD02-BROKEN-DERIVATIVE", "PD02-BROKEN-ASYMPTOTES",
    "PD02-BROKEN-SPECIES", "PD02-BOHM", "PD02-SELECTIVE-FIELD",
    "PD02-BOHM-FIELD-DEFINITION", "PD01-DISPATCH", "PD05-PARSER",
    "PD01-SCHEMA", "PD01-PARSER-BRIDGE", "PD01-REGISTRY-PROVENANCE",
    "PD13-BOUND-REGISTRY", "PD13-BOUND-DISPATCH", "PD13-HOST-GATE",
    "PD02-KSTAR-CONVERTER",
]

PARALLEL_ADVANCED_IDS = [
    "PD04-08-10-API", "PD04-QLT-0.001", "PD04-QLT-0.01",
    "PD04-QLT-1", "PD04-QLT-10",
] + [
    "PD05-06-{}-{}".format(model, point)
    for point in "ABCDEFG"
    for model in ("NLGCE_F", "NLGCE_N", "NLGC_E")
] + [
    "PD06-NLPA-B:kx_1e-3", "PD06-NLPA-C:kx_1e-3",
    "PD06-NLPA-D:kx_1e-3", "PD06-NLPA-F:kx_1e-3",
    "PD07-BROADENED-gaussian:0.1", "PD07-BROADENED-lorentzian:0.1",
]

PITCH_BASE_IDS = ["build"] + [
    "PAD{:02d}".format(index) for index in range(1, 13)
] + ["PAR{:02d}".format(index) for index in range(1, 12)]

# This is a registry contract, not a list of expected physical values.  The
# adapter independently asks the just-built local binary for --list and
# rejects any insertion, deletion or reordering before it reports a result.
PITCH_ACCEPTANCE_IDS = [
    "AP-ISO-01.em.f0.01", "AP-ISO-01.em.f0.003", "AP-ISO-01.em.f0.001",
    "AP-ISO-01.milstein.f0.01", "AP-ISO-01.milstein.f0.003",
    "AP-ISO-01.milstein.f0.001", "AP-ISO-01.rotation.f0.01",
    "AP-ISO-01.rotation.f0.003", "AP-ISO-01.rotation.f0.001",
    "AP-ISO-02.milstein.f0.003.t50", "AP-ISO-02.milstein.f0.001.t10",
    "AP-ISO-02.rotation.f0.003.t10", "AP-ISO-02.milstein.f0.001.t50",
    "AP-TIME-01.em", "AP-TIME-01.milstein", "AP-TIME-01.rotation",
    "AP-ZERO-02.eps0.01.milstein.f0.01",
    "AP-ZERO-02.eps0.01.milstein.f0.003",
    "AP-ZERO-02.eps0.01.milstein.f0.001",
    "AP-ZERO-02.eps0.048.milstein.f0.01",
    "AP-ZERO-02.eps0.048.milstein.f0.003",
    "AP-ZERO-02.eps0.048.milstein.f0.001",
    "AP-ZERO-02.eps1.milstein.f0.01",
    "AP-ZERO-02.eps1.milstein.f0.003",
    "AP-ZERO-02.eps1.milstein.f0.001", "AP-ZERO-02.eps0.em",
    "AP-ZERO-02.eps0.milstein", "AP-FOC-01.milstein.xi0.f0.003",
    "AP-FOC-01.milstein.xi0.1.f0.003",
    "AP-FOC-01.milstein.xi-0.1.f0.003",
    "AP-FOC-01.milstein.xi1.5.f0.003",
    "AP-FOC-01.milstein.xi-1.5.f0.003",
    "AP-FOC-01.milstein.xi5.f0.003",
    "AP-FOC-01.milstein.xi-5.f0.003",
    "AP-FOC-01.milstein.xi0.f0.001",
    "AP-FOC-01.milstein.xi0.1.f0.001",
    "AP-FOC-01.milstein.xi-0.1.f0.001",
    "AP-FOC-01.milstein.xi1.5.f0.001",
    "AP-FOC-01.milstein.xi-1.5.f0.001",
    "AP-FOC-01.milstein.xi5.f0.001",
    "AP-FOC-01.milstein.xi-5.f0.001", "AP-ZERO-01", "AP-RNG-01",
    "AP-CLOCK-01", "AP-FOC-02", "AP-TENSOR-01", "AP-FRAME-01",
    "AP-SPEC-01", "AP-WAVE-01", "AP-ENERGY-01", "AP-ENERGY-02",
    "AP-GEOM-01..04",
]


def normalize(local_id: str) -> str:
    """Make a collision-free global suffix while retaining the local ID in messages."""

    return re.sub(r"-+", "-", re.sub(r"[^A-Za-z0-9]+", "-", local_id)).strip("-").upper()


def global_id(scope: str, local_id: str) -> str:
    return "SEP-COMMON-{}-{}".format(scope, normalize(local_id))


EXPECTED_CASES = (
    [global_id("ARCHIVE", "SEP_COMMON01"), global_id("MFP", "VERIFY")]
    + [global_id("PARALLEL", "GENERATOR")]
    + [global_id("PARALLEL", value) for value in PARALLEL_CPP_IDS]
    + [global_id("PARALLEL-ADV", value) for value in PARALLEL_ADVANCED_IDS]
    + [global_id("PERP", "VERIFY"), global_id("PERP", "HIGH-PRECISION")]
    + [global_id("PITCH", value) for value in PITCH_BASE_IDS]
    + [global_id("PITCH", value) for value in PITCH_ACCEPTANCE_IDS]
)


def copy_source(source: Path, target: Path) -> None:
    """Copy only source/reference inputs; generated local state is never evidence."""

    def ignored(_directory: str, names: Sequence[str]) -> set:
        return {name for name in names if name in {
            "build", "output", "__pycache__", "sep_common.a",
            "libparallel_diffusion.a", "libperpendicular_diffusion.a",
            "libpitch_angle_diffusion.a"}}

    shutil.copytree(str(source), str(target), ignore=ignored)

    # The delivered perpendicular bundle's SHA256SUMS names the model Markdown
    # file as a bundle member, while this checkout stores the byte-identical
    # file one directory above the bundle so documentation links remain
    # natural. Recreate that declared packaging layout only in the disposable
    # copy. The verifier still authenticates the bytes against its frozen
    # SHA-256; no expected value or model content is manufactured here.
    perpendicular = target / "perpendicular_diffusion"
    model_document = perpendicular / "PERPENDICULAR_DIFFUSION_COEFFICIENT_MODEL.md"
    bundle_document = (perpendicular / "perpendicular_diffusion_model_data" /
                       model_document.name)
    if model_document.is_file() and not bundle_document.exists():
        shutil.copy2(str(model_document), str(bundle_document))


def load_json(path: Path) -> Mapping:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise RunnerError("cannot read authoritative report {}: {}".format(path, error))
    if not isinstance(payload, dict):
        raise RunnerError("authoritative report {} is not an object".format(path))
    return payload


def exact_ids(observed: Sequence[str], expected: Sequence[str], label: str) -> None:
    if list(observed) != list(expected) or len(set(observed)) != len(observed):
        raise RunnerError("{} registry mismatch: expected {} observed {}".format(
            label, list(expected), list(observed)))


def make_command(context: Mapping, directory: Path, target: str) -> List[str]:
    jobs = int(context["execution"].get("build_jobs", 1))
    return command(context, "make") + ["-C", str(directory), "-j{}".format(jobs), target]


def save_transcript(artifacts: Path, name: str, text: str) -> str:
    path = artifacts / name
    path.write_text(text or "", encoding="utf-8")
    return str(path)


def status_from_passed(value: object) -> str:
    return "PASS" if value is True else "FAIL"


def run_all(context: Mapping) -> Tuple[Dict[str, Dict[str, str]], List[str]]:
    repo = Path(str(context["repo_root"]))
    work = Path(str(context["work_dir"]))
    artifacts = Path(str(context["artifact_dir"]))
    timeout = int(context["test"].get("timeout_seconds", 7200))
    source = repo / "src" / "models" / "sep_common"
    model = work / "sep_common"
    artifacts.mkdir(parents=True, exist_ok=True)
    copy_source(source, model)
    evidence: List[str] = []
    results: Dict[str, Dict[str, str]] = {}

    # Pass the configured compiler vector through Make as one command-valued
    # variable.  No shell is involved in the adapter; Make performs its normal
    # recipe expansion and thus preserves machine-profile compiler wrappers.
    environment = {"CXX": shlex.join(command(context, "cxx")),
                   "PYTHONDONTWRITEBYTECODE": "1"}

    archive = run_command(make_command(context, model, "verify"), model,
                          timeout=timeout, env=environment)
    evidence.append(save_transcript(artifacts, "archive.txt", archive.stdout or ""))
    archive_ok = (archive.returncode == 0 and
                  (archive.stdout or "").count("SEP_COMMON01 PASS") == 1)
    results[global_id("ARCHIVE", "SEP_COMMON01")] = case(
        global_id("ARCHIVE", "SEP_COMMON01"), "PASS" if archive_ok else "FAIL",
        "local SEP_COMMON01 exact archive/PIC/MPI boundary" if archive_ok else
        "authoritative make verify failed or omitted its unique marker")

    mfp = run_command(make_command(context, model / "mean_free_path", "verify"),
                      model / "mean_free_path", timeout=timeout, env=environment)
    evidence.append(save_transcript(artifacts, "mean-free-path.txt", mfp.stdout or ""))
    mfp_ok = mfp.returncode == 0 and "host-neutral build and tests PASS" in (mfp.stdout or "")
    results[global_id("MFP", "VERIFY")] = case(
        global_id("MFP", "VERIFY"), "PASS" if mfp_ok else "FAIL",
        "complete local algebraic/reference/host-neutral verification" if mfp_ok else
        "authoritative mean-free-path verify failed or omitted its terminal marker")

    parallel_dir = model / "parallel_diffusion"
    parallel = run_command(make_command(context, parallel_dir, "verify"), parallel_dir,
                           timeout=timeout, env=environment)
    evidence.append(save_transcript(artifacts, "parallel-diffusion.txt",
                                    parallel.stdout or ""))
    cpp_report = load_json(parallel_dir / "build" / "test-report.json")
    advanced_report = load_json(parallel_dir / "build" / "advanced-test-report.json")
    if cpp_report.get("schema") != "parallel-diffusion-test-report-v1":
        raise RunnerError("parallel C++ report has an unsupported schema")
    cpp_rows = cpp_report.get("tests")
    if not isinstance(cpp_rows, list):
        # Older report revisions use results while retaining the same schema.
        cpp_rows = cpp_report.get("results")
    if not isinstance(cpp_rows, list):
        raise RunnerError("parallel C++ report has no tests/results array")
    cpp_ids = [str(row.get("id")) for row in cpp_rows if isinstance(row, dict)]
    exact_ids(cpp_ids, PARALLEL_CPP_IDS, "parallel C++")
    if advanced_report.get("schema") != "parallel_diffusion_advanced_test_report_v1":
        raise RunnerError("parallel advanced report has an unsupported schema")
    advanced_rows = advanced_report.get("records")
    if not isinstance(advanced_rows, list):
        raise RunnerError("parallel advanced report has no records array")
    advanced_ids = [str(row.get("name")) for row in advanced_rows if isinstance(row, dict)]
    exact_ids(advanced_ids, PARALLEL_ADVANCED_IDS, "parallel advanced")
    generator_ok = "generate_nlgce_coefficients.py --check" in (parallel.stdout or "")
    results[global_id("PARALLEL", "GENERATOR")] = case(
        global_id("PARALLEL", "GENERATOR"),
        "PASS" if generator_ok and parallel.returncode == 0 else "FAIL",
        "coefficient source regeneration check" if generator_ok else
        "authoritative generator check missing or verify failed")
    for row in cpp_rows:
        local = str(row["id"])
        passed = row.get("passed")
        if passed is None:
            passed = row.get("status") == "PASS"
        results[global_id("PARALLEL", local)] = case(
            global_id("PARALLEL", local), status_from_passed(passed),
            "local_id={}: {}".format(local, row.get("message", row.get("detail", ""))))
    for row in advanced_rows:
        local = str(row["name"])
        results[global_id("PARALLEL-ADV", local)] = case(
            global_id("PARALLEL-ADV", local), status_from_passed(row.get("passed")),
            "local_id={}: {}".format(local, row.get("detail", "")))
    for report in (parallel_dir / "build" / "test-report.json",
                   parallel_dir / "build" / "advanced-test-report.json"):
        target = artifacts / report.name
        shutil.copy2(str(report), str(target)); evidence.append(str(target))

    perp_dir = model / "perpendicular_diffusion"
    perp = run_command(make_command(context, perp_dir, "verify"), perp_dir,
                       timeout=timeout, env=environment)
    evidence.append(save_transcript(artifacts, "perpendicular-diffusion.txt",
                                    perp.stdout or ""))
    optional_skip = "SKIP review_additions_verification.py" in (perp.stdout or "")
    results[global_id("PERP", "VERIFY")] = case(
        global_id("PERP", "VERIFY"), "PASS" if perp.returncode == 0 else "FAIL",
        "complete mandatory local verification" if perp.returncode == 0 else
        "authoritative perpendicular-diffusion verify failed")
    results[global_id("PERP", "HIGH-PRECISION")] = case(
        global_id("PERP", "HIGH-PRECISION"),
        "SKIP" if optional_skip else ("PASS" if perp.returncode == 0 else "FAIL"),
        "optional mpmath audit unavailable" if optional_skip else
        "local review_additions high-precision audit executed")

    pitch_dir = model / "pitch_angle_diffusion"
    pitch_report = artifacts / "pitch-angle.json"
    pitch = run_command(command(context, "python") + [
        str(pitch_dir / "run_all_tests.py"), "--skip-acceptance",
        "--threads", str(max(1, int(context["execution"].get("build_jobs", 1)))),
        "--json", str(pitch_report),
    ], pitch_dir, timeout=timeout, env=environment)
    evidence.append(save_transcript(artifacts, "pitch-angle.txt", pitch.stdout or ""))
    payload = load_json(pitch_report)
    pitch_rows = payload.get("records")
    if not isinstance(pitch_rows, list):
        raise RunnerError("pitch-angle report has no records array")
    exact_ids([str(row.get("id")) for row in pitch_rows], PITCH_BASE_IDS + ["AP"],
              "pitch-angle base")
    for row in pitch_rows[:-1]:
        local = str(row["id"])
        local_status = str(row.get("status", ""))
        if local_status not in {"PASS", "FAIL", "SKIP"}:
            raise RunnerError("invalid pitch-angle status for {}".format(local))
        results[global_id("PITCH", local)] = case(
            global_id("PITCH", local), local_status,
            "local_id={}: {}".format(local, row.get("message", "")))
    # The JSON file may be complete before a later control-plane failure.
    # Reconcile it with the authoritative runner's documented exit contract
    # (1 iff any record is FAIL, otherwise 0) so a stale-looking report cannot
    # conceal an abnormal post-publication exit.  The final AP record is the
    # intentional --skip-acceptance SKIP in this T11 profile.
    reported_failure = any(
        isinstance(row, dict) and row.get("status") == "FAIL"
        for row in pitch_rows)
    expected_pitch_exit = 1 if reported_failure else 0
    if pitch.returncode != expected_pitch_exit:
        raise RunnerError(
            "pitch-angle exit/report mismatch: exit {} expected {}".format(
                pitch.returncode, expected_pitch_exit))
    acceptance_binary = pitch_dir / "build" / "test_pitch_angle_montecarlo"
    listed = run_command([str(acceptance_binary), "--list"], pitch_dir,
                         timeout=timeout, env=environment)
    evidence.append(save_transcript(artifacts, "pitch-acceptance-registry.txt",
                                    listed.stdout or ""))
    if listed.returncode != 0:
        raise RunnerError("pitch-angle acceptance --list failed")
    local_acceptance = [line.strip() for line in (listed.stdout or "").splitlines()
                        if line.strip()]
    exact_ids(local_acceptance, PITCH_ACCEPTANCE_IDS, "pitch acceptance")
    for local in local_acceptance:
        results[global_id("PITCH", local)] = case(
            global_id("PITCH", local), "SKIP",
            "local_id={}: full Monte-Carlo acceptance explicitly not run in routine T11 import".format(local))

    exact_ids(list(results), EXPECTED_CASES, "global SEP-common")
    return results, evidence


def main() -> int:
    context = load_context("Import complete SEP-common local test registries")
    try:
        # A local registry check cannot by itself detect that the central
        # manifest deleted a case and thereby stopped asking the adapter to
        # report it. Freeze the ordered cross-boundary contract before
        # applying an individual-case selection. Keep this preflight inside
        # the receipt boundary: registry drift is an infrastructure ERROR for
        # every requested identity, not an unstructured process exit.
        configured = list(context["test"].get("expected_cases", ()))
        if configured != EXPECTED_CASES:
            raise RunnerError(
                "SEP-common manifest registry drift: configured {} ordered cases, "
                "adapter requires {}".format(len(configured), len(EXPECTED_CASES)))
        requested = list(context["test"].get(
            "selected_cases", configured))
        if (not requested or len(set(requested)) != len(requested) or
                any(identifier not in EXPECTED_CASES for identifier in requested)):
            raise RunnerError("invalid SEP-common case selection: {}".format(requested))
        all_results, evidence = run_all(context)
        selected = [all_results[identifier] for identifier in requested]
        severity = {"PASS": 0, "SKIP": 1, "FAIL": 2, "ERROR": 3}
        status = max(selected, key=lambda row: severity[row["status"]])["status"]
        return write_result(
            context, status,
            "{} selected of {} validated SEP-common registry cases; acceptance SKIPs are not qualification".format(
                len(selected), len(EXPECTED_CASES)),
            cases=selected, evidence=evidence,
            metrics={"complete_registry_cases": len(EXPECTED_CASES),
                     "pitch_acceptance_cases": len(PITCH_ACCEPTANCE_IDS)})
    except Exception as error:
        return error_result(context, error)


if __name__ == "__main__":
    main_guard(main)
