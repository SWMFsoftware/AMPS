#!/usr/bin/env python3
"""Build configuration-separated executables and run the M0 I01-I10 campaign.

AMPS generates source under build/, so SPICE, no-SPICE, and Na/Na+ objects must
never be mixed.  This runner removes build/ before every configuration, keeps
each complete log and executable, and restores the normal SPICE Moon build at
the end.  It does not update a regression reference or weaken a failed gate.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

from integration_common import ERROR, FAIL, PASS, REPO, run


def command_log(command: list[str], log: Path) -> int:
    """Run one build command at repository root and append its exact stream."""
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("a", encoding="utf-8") as stream:
        stream.write("COMMAND: " + " ".join(command) + "\n")
        completed = subprocess.run(command, cwd=REPO, stdout=stream,
                                   stderr=subprocess.STDOUT, check=False)
        stream.write(f"EXIT: {completed.returncode}\n")
    return completed.returncode


def clean_configure_build(application: str, spice_path: str,
                          kernels: str, output: Path, label: str) -> Path:
    """Create one internally consistent generated tree and retain its binary."""
    build = REPO / "build"
    if build.exists():
        shutil.rmtree(build)
    log = output / "build_logs" / f"{label}.log"
    config = [str(REPO / "Config.pl"), f"-application={application}",
              f"-spice-path={spice_path}", f"-spice-kernels={kernels}"]
    if command_log(config, log) != 0 or command_log(["make", "-j"], log) != 0:
        raise RuntimeError(f"{label} configuration/build failed; see {log}")
    source = REPO / "amps"
    destination = output / "executables" / label / "amps"
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    return destination


def load_result(output: Path, test_id: str) -> dict[str, object]:
    path = output / test_id / "result.json"
    if not path.is_file():
        return {"id": test_id, "status": ERROR,
                "reason": f"runner did not create {path}"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--spice-path", default="/home/vtenishe/SPICE/cspice")
    parser.add_argument("--spice-kernels", default="/home/vtenishe/SPICE/Kernels")
    parser.add_argument("--no-restore-build", action="store_true",
                        help="leave the final test fixture configured instead of restoring the normal SPICE Moon build")
    arguments = parser.parse_args()
    output = arguments.output_dir.resolve()
    if output.exists() and any(output.iterdir()):
        parser.error(f"output directory must be new or empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, object]] = []
    prerequisite_artifacts: dict[str, object] = {}
    generated_i08 = REPO / "input/moon-i08.input"

    try:
        spice_executable = clean_configure_build(
            "moon", arguments.spice_path, arguments.spice_kernels,
            output, "moon-spice")
        # Execute every frozen U ID against the same clean SPICE build.  Later
        # capabilities remain honest SKIPPED; M0 uses only the documented local
        # prerequisite subset.  U06 is rerun below in its dedicated ion build.
        u_output = output / "u-series"
        u_command = [sys.executable, str(REPO / "srcMoon/test/run_all_test.py"),
                     "--skip-build", "--output-dir", str(u_output)]
        u_log = output / "build_logs" / "u-series.log"
        prerequisite_artifacts["u_series_returncode"] = command_log(u_command, u_log)
        prerequisite_artifacts["u_series_summary"] = str(u_output / "run_all_summary.json")
        for test_id in ("I01", "I02", "I03", "I06", "I07", "I09", "I10"):
            run(test_id, spice_executable, output)
            results.append(load_result(output, test_id))

        neutral_executable = clean_configure_build(
            "moon", "nospice", "nospice", output, "moon-nospice")
        run("I05", neutral_executable, output)
        results.append(load_result(output, "I05"))

        # The generator enforces that only the declared species and SPICE
        # controls differ from the canonical production input.
        generation_log = output / "build_logs" / "i08-input.log"
        if command_log([sys.executable, str(REPO / "srcMoon/test/generate_i08_input.py"),
                        "--source", str(REPO / "input/moon.input"),
                        "--output", str(generated_i08)], generation_log) != 0:
            raise RuntimeError(f"I08 input generation failed; see {generation_log}")
        ion_executable = clean_configure_build(
            "moon-i08", "nospice", "nospice", output, "moon-i08-nospice")
        u06_output = output / "u06-ion-fixture"
        u06_command = [sys.executable, str(REPO / "srcMoon/test/U06_lorentz/test.py"),
                       "--build-dir", str(REPO / "build"),
                       "--output-dir", str(u06_output)]
        u06_log = output / "build_logs" / "u06-ion-fixture.log"
        prerequisite_artifacts["u06_ion_returncode"] = command_log(u06_command, u06_log)
        prerequisite_artifacts["u06_ion_result"] = str(u06_output / "U06/result.json")
        run("I08", ion_executable, output)
        results.append(load_result(output, "I08"))

        # I04 deliberately uses the repository's canonical regression target.
        # Its reference is never regenerated by this campaign.
        build = REPO / "build"
        if build.exists():
            shutil.rmtree(build)
        regression_log = output / "build_logs" / "test_Moon.log"
        regression_code = command_log(["make", "-j", "test_Moon"], regression_log)
        i04_dir = output / "I04"; i04_dir.mkdir(parents=True, exist_ok=True)
        diff = REPO / "test_Moon.diff"
        if diff.is_file():
            shutil.copy2(diff, i04_dir / "test_Moon.diff")
        regression_executable = REPO / "run_test_Moon/amps"
        if regression_code == 0 and regression_executable.is_file():
            run("I04", regression_executable, output)
        else:
            result = {"id": "I04", "status": ERROR,
                      "reason": f"canonical regression failed; see {regression_log}"}
            (i04_dir / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        results.append(load_result(output, "I04"))
    except (OSError, RuntimeError, ValueError) as exception:
        results.append({"id": "campaign", "status": ERROR, "reason": str(exception)})
    finally:
        # This file is a deterministic generated test input, never a source
        # artifact.  Remove only the exact path owned by this runner.
        if generated_i08.is_file():
            generated_i08.unlink()
        if not arguments.no_restore_build:
            try:
                clean_configure_build("moon", arguments.spice_path,
                                      arguments.spice_kernels, output,
                                      "moon-spice-restored")
            except (OSError, RuntimeError) as exception:
                results.append({"id": "restore-build", "status": ERROR,
                                "reason": str(exception)})

    # Resolve the phase-appropriate local gate explicitly.  A later U test may
    # remain SKIPPED, but every kernel used directly by I01-I10 must PASS.
    prerequisite_ids = ("U01", "U02", "U03", "U04", "U05", "U06",
                        "U12", "U14", "U22")
    local_by_id: dict[str, dict[str, object]] = {}
    u_summary_path = output / "u-series/run_all_summary.json"
    if u_summary_path.is_file():
        for item in json.loads(u_summary_path.read_text(encoding="utf-8")).get("results", []):
            local_by_id[str(item.get("id"))] = item
    u06_path = output / "u06-ion-fixture/U06/result.json"
    if u06_path.is_file():
        local_by_id["U06"] = json.loads(u06_path.read_text(encoding="utf-8"))
    local_results = [local_by_id.get(test_id, {"id": test_id, "status": ERROR,
                                              "reason": "missing prerequisite result"})
                     for test_id in prerequisite_ids]
    local_gate_status = PASS if all(item.get("status") == PASS for item in local_results) else FAIL

    # Sort by frozen ID for the human and machine-readable summaries even
    # though configuration grouping determines execution order.
    results.sort(key=lambda item: str(item["id"]))
    counts = {status: sum(item.get("status") == status for item in results)
              for status in (PASS, FAIL, "SKIPPED", ERROR)}
    integration_gate_status = PASS if counts[PASS] == 10 and not any(
        counts[status] for status in (FAIL, "SKIPPED", ERROR)) else FAIL
    m0_gate_status = PASS if integration_gate_status == PASS and local_gate_status == PASS else FAIL
    summary = {"counts": counts, "results": results,
               "local_prerequisite_results": local_results,
               "local_prerequisite_status": local_gate_status,
               "integration_gate_status": integration_gate_status,
               "m0_gate_status": m0_gate_status,
               "local_prerequisite_artifacts": prerequisite_artifacts}
    (output / "baseline_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("srcMoon I01-I10 baseline summary")
    for status in (PASS, FAIL, "SKIPPED", ERROR):
        print(f"  {status.lower():7s}: {counts[status]}")
    for result in results:
        print(f"  {result.get('status', ERROR):7s} {result['id']}: {result.get('reason', '')}")
    print(f"  local prerequisites: {local_gate_status}")
    print(f"  M0 gate: {m0_gate_status}")
    return 2 if counts[ERROR] else (1 if m0_gate_status != PASS else 0)


if __name__ == "__main__":
    raise SystemExit(main())
