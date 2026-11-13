#!/usr/bin/env python3
"""Single top-level runner for the layered srcMoon verification campaign.

The scientific work remains in the existing per-ID runners.  This driver is
deliberately an orchestrator: it audits every implemented test contract,
delegates configuration-separated M0 work to ``run_integration_tests.py``,
delegates U/M1/M2 work to ``run_all_test.py``, and combines their structured
results without reinterpreting a metric or weakening a threshold.

I21--I46 are frozen roadmap IDs but do not yet have executable runners.  They
are listed explicitly and appear as SKIPPED in ``release`` scope.  A missing
production capability, independent oracle, or authoritative dataset must not
be disguised by a generic runner that returns PASS.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
from typing import Any


TEST_ROOT = Path(__file__).resolve().parent
REPO_ROOT = TEST_ROOT.parents[1]
SOURCE_ROADMAP = REPO_ROOT / "srcMoon/development_and_validation.md"
DATA_ROOT = Path("/data/vtenishe/moon_validation_data")
DATA_ROADMAP = DATA_ROOT / "development_and_validation.md"
DEFAULT_OUTPUT_PARENT = REPO_ROOT / "test_output/srcMoon/global-campaign"
VALID_STATUSES = ("PASS", "FAIL", "SKIPPED", "ERROR")

# AGENTS.md freezes these fields for every test specification.  Extra fields
# are permitted, but an implemented runner cannot silently omit one of the
# common scientific/provenance questions.
REQUIRED_ACCEPTANCE_FIELDS = {
    "id", "name", "purpose", "production_capability",
    "production_entry_points", "configuration", "inputs",
    "independent_reference", "units", "frame_time_convention",
    "deterministic_stochastic_controls", "procedure", "metrics",
    "acceptance_criteria", "expected_pass_result",
    "expected_failure_modes", "required_data", "data_hashes",
    "predecessor_gates", "artifacts", "status_semantics",
}

PLANNED_I_TESTS: dict[str, str] = {
    "I21": "Sodium photoionization survival and Na-to-Na+ conversion",
    "I22": "Photoionization for He, Ne, Ar, and future H2O/OH",
    "I23": "Electron-impact ionization rate coefficients",
    "I24": "Charge-exchange reactions with solar-wind ions",
    "I25": "Measured solar-wind and magnetospheric forcing",
    "I26": "Helium implantation, re-emission, and endogenous source",
    "I27": "Neon source and noncondensable accommodation",
    "I28": "Radiogenic 40Ar geography and transient outgassing",
    "I29": "Sodium PSD/TD/sputtering/impact source partition",
    "I30": "Time-dependent meteoroid forcing and impact vaporization",
    "I31": "H2O/OH surface migration and photochemistry",
    "I32": "Kaguya TVIS sodium limb/exosphere brightness validation",
    "I33": "Ground-based sodium tail and Earth-focusing validation",
    "I34": "LADEE NMS helium local-time and altitude validation",
    "I35": "LRO/LAMP helium response with ARTEMIS forcing",
    "I36": "LADEE NMS neon local-time validation",
    "I37": "LADEE NMS argon diurnal/selenographic validation",
    "I38": "Apollo 17 LACE and LADEE argon synthesis",
    "I39": "LADEE NMS meteoroid-driven water releases",
    "I40": "Kaguya PACE pickup-ion and ion-escape comparison",
    "I41": "Frozen cross-species cross-instrument validation",
    "I42": "Automated one-command Moon V&V release gate",
    "I43": "Sensitivity, uncertainty, and identifiability",
    "I44": "Monte Carlo and random-seed robustness",
    "I45": "MPI/OpenMP reproducibility and scaling",
    "I46": "Release-readiness science scorecard",
}


def sha256(path: Path) -> str:
    """Return the byte-level identity recorded by the aggregate receipt."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_revision() -> str:
    """Return the current source revision without changing the worktree."""
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    return completed.stdout.strip() if completed.returncode == 0 else "unknown"


def default_output() -> Path:
    """Return a unique campaign directory so evidence is never overwritten."""
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return DEFAULT_OUTPUT_PARENT / stamp


def implemented_layout() -> dict[str, Path]:
    """Map every implemented frozen ID to its directory."""
    layout: dict[str, Path] = {}
    for path in sorted(TEST_ROOT.glob("U[0-9][0-9]_*")):
        if path.is_dir():
            layout[path.name[:3]] = path
    for number in range(1, 21):
        test_id = f"I{number:02d}"
        layout[test_id] = TEST_ROOT / "integration-tests" / test_id
    return layout


def audit_architecture() -> list[str]:
    """Return structural contract errors without compiling or running physics."""
    errors: list[str] = []
    layout = implemented_layout()
    expected = ({f"U{number:02d}" for number in range(1, 27)} |
                {f"I{number:02d}" for number in range(1, 21)})
    if set(layout) != expected:
        errors.append(
            "implemented layout differs from frozen U01-U26/I01-I20 set: "
            f"missing={sorted(expected-set(layout))}, "
            f"extra={sorted(set(layout)-expected)}")

    for test_id in sorted(expected):
        directory = layout.get(test_id)
        if directory is None:
            continue
        required = {
            "runner": directory / "test.py",
            "readme": directory / "README.md",
            "acceptance": directory / "reference/acceptance.json",
        }
        for label, path in required.items():
            if not path.is_file():
                errors.append(f"{test_id}: missing {label}: {path}")
        acceptance_path = required["acceptance"]
        if not acceptance_path.is_file():
            continue
        try:
            contract = json.loads(acceptance_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exception:
            errors.append(f"{test_id}: unreadable acceptance.json: {exception}")
            continue
        if contract.get("id") != test_id:
            errors.append(
                f"{test_id}: acceptance ID is {contract.get('id')!r}")
        missing_fields = sorted(REQUIRED_ACCEPTANCE_FIELDS - contract.keys())
        if missing_fields:
            errors.append(
                f"{test_id}: acceptance contract lacks " +
                ", ".join(missing_fields))

        # I14-I20 are the currently implemented comparison-facing M2 layer.
        # If their criterion is evaluated, the same rows must be preserved as
        # CSV and rendered to both raster and vector formats.  A skipped test
        # produces none of these artifacts, as documented by each runner.
        if test_id.startswith("I") and 14 <= int(test_id[1:]) <= 20:
            artifacts = set(contract.get("artifacts", []))
            comparison = {"comparison.csv", "comparison.png", "comparison.eps"}
            if not comparison <= artifacts:
                errors.append(
                    f"{test_id}: comparison contract lacks " +
                    ", ".join(sorted(comparison-artifacts)))

    planned = set(PLANNED_I_TESTS)
    expected_planned = {f"I{number:02d}" for number in range(21, 47)}
    if planned != expected_planned:
        errors.append(
            "planned registry must contain exactly I21-I46: "
            f"missing={sorted(expected_planned-planned)}, "
            f"extra={sorted(planned-expected_planned)}")
    return errors


def run_logged(command: list[str], log: Path, stream: bool) -> int:
    """Execute one phase runner and preserve its combined outer transcript."""
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8") as output:
        output.write("command: " + " ".join(command) + "\n")
        output.flush()
        if not stream:
            completed = subprocess.run(
                command, cwd=REPO_ROOT, stdout=output,
                stderr=subprocess.STDOUT, check=False)
            return completed.returncode
        process = subprocess.Popen(
            command, cwd=REPO_ROOT, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, text=True, bufsize=1)
        assert process.stdout is not None
        for line in process.stdout:
            output.write(line)
            output.flush()
            print(line, end="", flush=True)
        return process.wait()


def read_json(path: Path) -> dict[str, Any]:
    """Read a required delegated-runner result with a precise error."""
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exception:
        raise RuntimeError(f"delegated runner did not create {path}") from exception
    except (OSError, json.JSONDecodeError) as exception:
        raise RuntimeError(f"cannot read delegated result {path}: {exception}") from exception


def normalized_results(payload: dict[str, Any], source: str) -> list[dict[str, Any]]:
    """Validate and annotate per-ID records from an authoritative phase runner."""
    results = payload.get("results")
    if not isinstance(results, list):
        raise RuntimeError(f"{source} summary has no results list")
    normalized: list[dict[str, Any]] = []
    for row in results:
        if not isinstance(row, dict):
            raise RuntimeError(f"{source} summary contains a non-object result")
        test_id = row.get("id")
        status = row.get("status")
        if (not isinstance(test_id, str) or
                status not in VALID_STATUSES):
            raise RuntimeError(
                f"{source} result lacks a valid id/status: {row!r}")
        record = dict(row)
        record["phase_summary"] = source
        normalized.append(record)
    return normalized


def skipped_record(test_id: str, name: str, reason: str) -> dict[str, Any]:
    """Create an explicit not-run roadmap record; never a placeholder PASS."""
    return {
        "id": test_id, "name": name, "status": "SKIPPED",
        "execution": "not-run", "runner": None, "reason": reason,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--phase", choices=("audit", "local", "m0", "frontier",
                             "implemented", "release", "m2-development"), default="audit",
        help=("audit structure; run U tests; run M0; run I11-I20; run all "
              "implemented layers; include planned release gaps; or run "
              "bounded supporting M2 development without waiving I13"))
    parser.add_argument("--list", action="store_true",
                        help="list all frozen U01-U26/I01-I46 IDs and exit")
    parser.add_argument("--output-dir", type=Path,
                        help="new/empty global campaign artifact directory")
    parser.add_argument("--campaign-mode", choices=("regression", "qualification"),
                        default="regression")
    parser.add_argument("--mpi-ranks", type=int, default=8)
    parser.add_argument("--m2-amps",type=Path,
                        help="explicit M2 executable with matching generated build for m2-development")
    parser.add_argument("--spice-path", type=Path,
                        default=Path("/home/vtenishe/SPICE/cspice"))
    parser.add_argument("--spice-kernels", type=Path,
                        default=Path("/home/vtenishe/SPICE/Kernels"))
    parser.add_argument("--skip-build", action="store_true",
                        help="reuse the current build for local/frontier phases only")
    parser.add_argument("--stream", action="store_true")
    parser.add_argument("--strict-skips", action="store_true")
    arguments = parser.parse_args()
    if arguments.mpi_ranks < 1:
        parser.error("--mpi-ranks must be positive")
    if arguments.phase == "m2-development" and arguments.skip_build:
        parser.error("m2-development reuse requires explicit --m2-amps and its matching generated tree")

    layout = implemented_layout()
    if arguments.list:
        for test_id in sorted(layout):
            contract = read_json(layout[test_id] / "reference/acceptance.json")
            print(f"{test_id}  implemented  {contract.get('name', '')}")
        for test_id, name in PLANNED_I_TESTS.items():
            print(f"{test_id}  planned      {name}")
        return 0

    output = (arguments.output_dir.resolve() if arguments.output_dir
              else default_output().resolve())
    if output.exists() and any(output.iterdir()):
        parser.error(
            f"output directory is not empty: {output}; prior evidence is immutable")
    output.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    architecture_errors = audit_architecture()
    results: list[dict[str, Any]] = []
    campaign_errors: list[str] = list(architecture_errors)
    phase_runs: list[dict[str, Any]] = []
    development_checks: list[dict[str, Any]] = []

    if not architecture_errors and arguments.phase != "audit":
        if arguments.phase == "m2-development":
            phase_output=output/"m2-development"
            command=[sys.executable,"-B",str(TEST_ROOT/"run_m2_development.py"),
                     "--output-dir",str(phase_output),"--mpi-ranks",str(arguments.mpi_ranks),
                     "--spice-path",str(arguments.spice_path),"--spice-kernels",str(arguments.spice_kernels)]
            if arguments.m2_amps:command.extend(["--amps",str(arguments.m2_amps.resolve())])
            log=output/"m2-development.log"
            code=run_logged(command,log,arguments.stream)
            phase_runs.append({"phase":"m2-development","command":command,"returncode":code,"log":str(log)})
            try:
                child=read_json(phase_output/"development_summary.json")
                results.extend(child["results"])
                development_checks.extend(child["supporting_checks"])
            except (OSError,ValueError,KeyError) as error:campaign_errors.append(str(error))
            if code not in (0,1,2):campaign_errors.append(f"M2 development runner exited {code}")
        if arguments.phase == "local":
            phase_output = output / "local"
            command = [sys.executable, str(TEST_ROOT / "run_all_test.py"),
                       "--output-dir", str(phase_output)]
            for number in range(1, 27):
                command.extend(["--only", f"U{number:02d}"])
            if arguments.skip_build:
                command.append("--skip-build")
            if arguments.stream:
                command.append("--stream")
            code = run_logged(command, output / "logs/local.log", arguments.stream)
            phase_runs.append({"phase": "local", "command": command,
                               "returncode": code})
            try:
                results.extend(normalized_results(
                    read_json(phase_output / "run_all_summary.json"), "local"))
            except RuntimeError as exception:
                campaign_errors.append(str(exception))

        if arguments.phase in {"m0", "implemented", "release"}:
            m0_output = output / "m0"
            command = [
                sys.executable, str(TEST_ROOT / "run_integration_tests.py"),
                "--output-dir", str(m0_output),
                "--spice-path", str(arguments.spice_path.resolve()),
                "--spice-kernels", str(arguments.spice_kernels.resolve()),
            ]
            code = run_logged(command, output / "logs/m0.log", arguments.stream)
            phase_runs.append({"phase": "m0", "command": command,
                               "returncode": code})
            try:
                baseline = read_json(m0_output / "baseline_summary.json")
                local = read_json(m0_output / "u-series/run_all_summary.json")
                results.extend(normalized_results(local, "m0/u-series"))
                results.extend(normalized_results(baseline, "m0/I01-I10"))
                m0_pass = baseline.get("m0_gate_status") == "PASS"
            except RuntimeError as exception:
                campaign_errors.append(str(exception))
                m0_pass = False
        else:
            m0_pass = arguments.phase == "frontier"

        if arguments.phase in {"frontier", "implemented", "release"}:
            if m0_pass:
                frontier_output = output / "frontier"
                command = [
                    sys.executable, str(TEST_ROOT / "run_all_test.py"),
                    "--output-dir", str(frontier_output),
                    "--campaign-mode", arguments.campaign_mode,
                    "--mpi-ranks", str(arguments.mpi_ranks),
                    "--spice-path", str(arguments.spice_path.resolve()),
                    "--spice-kernels", str(arguments.spice_kernels.resolve()),
                ]
                for number in range(11, 21):
                    command.extend(["--only", f"I{number:02d}"])
                # A completed M0 run restores a canonical Moon/SPICE build.
                # Qualification may additionally need the dedicated Na+Ar M2
                # executable, so let run_all_test build both in that mode.
                if (arguments.skip_build or
                        (arguments.phase in {"implemented", "release"} and
                         arguments.campaign_mode == "regression")):
                    command.append("--skip-build")
                if arguments.stream:
                    command.append("--stream")
                code = run_logged(
                    command, output / "logs/frontier.log", arguments.stream)
                phase_runs.append({"phase": "frontier", "command": command,
                                   "returncode": code})
                try:
                    results.extend(normalized_results(
                        read_json(frontier_output / "run_all_summary.json"),
                        "frontier/I11-I20"))
                except RuntimeError as exception:
                    campaign_errors.append(str(exception))
            else:
                for number in range(11, 21):
                    test_id = f"I{number:02d}"
                    contract = read_json(
                        layout[test_id] / "reference/acceptance.json")
                    results.append(skipped_record(
                        test_id, str(contract.get("name", test_id)),
                        "not executed because the preceding M0 aggregate gate did not PASS"))

    if arguments.phase == "release":
        for test_id, name in PLANNED_I_TESTS.items():
            results.append(skipped_record(
                test_id, name,
                "frozen roadmap ID has no executable production-linked runner yet; "
                "see TESTING_PROCEDURE.md and STATUS.md"))

    # Duplicate IDs would make a later layer overwrite evidence from an
    # earlier layer.  Treat that as runner ERROR rather than choosing one.
    ids = [str(row.get("id")) for row in results]
    duplicates = sorted({test_id for test_id in ids if ids.count(test_id) > 1})
    if duplicates:
        campaign_errors.append("duplicate aggregate result IDs: " + ", ".join(duplicates))

    counts = {status: sum(row.get("status") == status for row in results)
              for status in VALID_STATUSES}
    counts["executed_or_explicitly_skipped"] = len(results)
    campaign_status = "PASS"
    if campaign_errors or counts["ERROR"] or any(r.get("status")=="ERROR" for r in development_checks):
        campaign_status = "ERROR"
    elif counts["FAIL"] or any(r.get("status")=="FAIL" for r in development_checks):
        campaign_status = "FAIL"
    elif counts["SKIPPED"] or any(r.get("status")=="SKIPPED" for r in development_checks):
        campaign_status = "SKIPPED"

    summary = {
        "schema": "srcMoon-global-test-campaign-v1",
        "phase": arguments.phase,
        "campaign_mode": arguments.campaign_mode,
        "campaign_status": campaign_status,
        "strict_skips": arguments.strict_skips or arguments.phase == "release",
        "source_revision": git_revision(),
        "roadmaps": {
            "source": {"path": str(SOURCE_ROADMAP),
                       "sha256": sha256(SOURCE_ROADMAP)},
            "validation_data": {"path": str(DATA_ROADMAP),
                                "sha256": sha256(DATA_ROADMAP)},
        },
        "roadmaps_identical": sha256(SOURCE_ROADMAP) == sha256(DATA_ROADMAP),
        "architecture_audit": {
            "status": "PASS" if not architecture_errors else "ERROR",
            "errors": architecture_errors,
            "implemented_ids": sorted(layout),
            "planned_ids": list(PLANNED_I_TESTS),
        },
        "phase_runs": phase_runs,
        "development_checks": development_checks,
        "development_counts": {s:sum(r.get("status")==s for r in development_checks) for s in VALID_STATUSES},
        "counts": counts,
        "campaign_errors": campaign_errors,
        "duration_seconds": time.monotonic() - started,
        "results": sorted(results, key=lambda row: str(row.get("id"))),
    }
    summary_path = output / "test_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("srcMoon global campaign summary")
    print(f"  phase:    {arguments.phase}")
    print(f"  status:   {campaign_status}")
    print(f"  records:  {len(results)}")
    for status in VALID_STATUSES:
        print(f"  {status.lower():7s}: {counts[status]}")
    if development_checks:
        print("  supporting development: "+json.dumps(summary["development_counts"],sort_keys=True))
    # A supporting failure must be as discoverable as a formal one, without
    # adding supporting records to the permanent-ID result namespace.
    for row in results + development_checks:
        if row.get("status") in ("FAIL", "ERROR"):
            print(f"  {row['status']} {row.get('check', row.get('id', 'unknown'))}")
            print(f"    runner: {row.get('runner') or 'see phase_runs in test_summary.json'}")
            print(f"    log:    {row.get('log') or row.get('artifacts', {}).get('stdout') or output}")
    if campaign_errors:
        print("  campaign errors:")
        for error in campaign_errors:
            print(f"    - {error}")
    print(f"  result:   {summary_path}")

    if campaign_errors or counts["ERROR"] or any(r.get("status")=="ERROR" for r in development_checks):
        return 2
    if counts["FAIL"] or any(r.get("status")=="FAIL" for r in development_checks):
        return 1
    if (arguments.strict_skips or arguments.phase == "release") and (
            counts["SKIPPED"] or any(r.get("status")=="SKIPPED" for r in development_checks)):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
