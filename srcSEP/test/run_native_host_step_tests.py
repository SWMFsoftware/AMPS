#!/usr/bin/env python3
"""Run actual host/PIC birth-clock diagnostics in separate native processes.

The required fixture supplies explicit controlled numerical/source parameters.
No missing physical profile, thermal closure, guiding-center boundary law or
production qualification is supplied by this software timing campaign. Every
native job owns its cwd and files; PIC private first-call/statistics state is
isolated by process termination, not an incomplete in-process restore claim.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
import sys
import time

sys.dont_write_bytecode = True
# Reuse the existing native process-group timeout/hash implementation. This
# file lives beside that source-owned runner in source and copied packages.
from run_fresh_process_checkpoint_tests import binary_hash, run_role

SOURCE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SOURCE / "validation"))
try:
    from sep_native_execution import create_native_work_directory
finally:
    sys.path.pop(0)

CASES = {"HOST01": "parker", "HOST02": "fte-dmumu", "HOST03": "fte-mfp"}


def load_fixture(path: Path) -> dict:
    """Read a complete flat contract; no nested path or physics inference."""
    def unique(pairs):
        output = {}
        for key, value in pairs:
            if key in output:
                raise ValueError("duplicate host fixture key: " + key)
            output[key] = value
        return output
    value = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=unique)
    fields = {"schema", "steps", "source_macros", "source_energy_mev", "source_mu", "constant_dmumu_per_s"}
    if not isinstance(value, dict) or set(value) != fields or value["schema"] != "srcsep-native-host-step-fixture-v1":
        raise ValueError("host fixture requires its complete version-one field set")
    for key, maximum in (("steps", 2), ("source_macros", 8)):
        if type(value[key]) is not int or not 1 <= value[key] <= maximum:
            raise ValueError(key + " is outside its explicit bounded integer range")
    for key in ("source_energy_mev", "source_mu", "constant_dmumu_per_s"):
        if type(value[key]) not in (int, float) or not math.isfinite(value[key]) or value[key] <= 0:
            raise ValueError(key + " must be finite and positive")
    if value["source_mu"] > 1:
        raise ValueError("source_mu exceeds one")
    return value


def fixture_protocol(values: dict) -> str:
    """Serialize exactly the reviewed flat controls, without path/default rules."""
    return "SEP_NATIVE_HOST_STEP 1\n" + "".join(
        key + " " + format(values[key], ".17g") + "\n" for key in
        ("steps", "source_macros", "source_energy_mev", "source_mu", "constant_dmumu_per_s"))


def inspect_traces(directory: Path, identifier: str, mover: str, ranks: int,
                   expected_steps: int, require_valid_phase: bool = True) -> list[dict]:
    """Check software timing/coverage; never generate a trajectory reference.

    Native C++ checks the existing analytic governing law and completed shock
    geometry. This reader checks retained clocks/durations/identity and that
    every actual step-start birth had exactly one pending-first mover call.
    """
    sources, first_moves, source_steps, checks = [], [], set(), []
    for rank in range(ranks):
        local = directory if rank == 0 else directory / ("rank-" + str(rank))
        paths = sorted(local.glob("host-step-trace-*.json"))
        if not paths:
            raise ValueError("native host immutable trace snapshots missing: " + str(local))
        previous, document, counter_policy = [], {}, None
        for sequence, path in enumerate(paths):
            document = json.loads(path.read_text(encoding="utf-8"))
            expected_name = "host-step-trace-" + format(sequence, "06d") + ".json"
            if (not re.fullmatch(r"host-step-trace-[0-9]{6}\.json", path.name) or path.name != expected_name or
                    document.get("schema") != "srcsep-native-host-step-trace-v2" or
                    document.get("case_id") != identifier or document.get("mover") != mover or
                    document.get("rank") != rank or type(document.get("counter_at_end")) is not bool or
                    type(document.get("snapshot_sequence")) is not int or document["snapshot_sequence"] != sequence or
                    type(document.get("complete")) is not bool or
                    document.get("complete") != (sequence == len(paths)-1) or
                    type(document.get("requested_host_steps")) is not int or document["requested_host_steps"] != expected_steps or
                    type(document.get("completed_host_steps")) is not int or not 0 <= document["completed_host_steps"] <= expected_steps or
                    not isinstance(document.get("events"), list)):
                raise ValueError("missing/mismatched native immutable trace contract: " + str(path))
            if sequence and document["counter_at_end"] != counter_policy:
                raise ValueError("compiled PIC counter policy changed between trace snapshots")
            counter_policy = document["counter_at_end"]
            # Snapshots are immutable cumulative prefixes, so a later file may
            # append actual events but cannot rewrite or discard early evidence.
            if document["events"][:len(previous)] != previous:
                raise ValueError("native host cumulative trace changed or lost earlier events")
            previous = document["events"]
            checks.append({"status": "PASS", "message": "immutable cumulative rank trace retained", "artifact": str(path)})
        if document["completed_host_steps"] != expected_steps:
            raise ValueError("complete final native trace does not include every requested host step")
        for event in document["events"]:
            if (not isinstance(event, dict) or type(event.get("particle")) is not int or event["particle"] <= 0 or
                    type(event.get("step")) is not int or not 0 <= event["step"] < expected_steps or
                    event.get("phase") not in ("source", "mover") or type(event.get("pending_first")) is not bool):
                raise ValueError("invalid native host event identity/phase")
            for key in ("clock_s", "parent_start_s", "parent_end_s", "first_duration_s", "requested_s", "source_radius_m"):
                if type(event.get(key)) not in (int, float) or not math.isfinite(event[key]):
                    raise ValueError("nonfinite/missing native host event " + key)
            if event["first_duration_s"] <= 0 or event["parent_end_s"] < event["parent_start_s"] or event["source_radius_m"] <= 0:
                raise ValueError("invalid declared host event interval/radius")
            # Keep explicit duration separate from subtraction of clock epochs:
            # the native birth ABI deliberately preserves that information.
            if require_valid_phase and (event["phase"] == "source" or event["pending_first"]):
                expected_clock = event["parent_start_s"] if document["counter_at_end"] else event["parent_end_s"]
                if event["clock_s"] != expected_clock:
                    raise ValueError("source/first mover clock differs from the actual compiled host policy")
            if event["phase"] == "source":
                if type(event.get("cartesian_matches")) is not bool or not event["pending_first"]:
                    raise ValueError("source phase/Cartesian flag is missing")
                if require_valid_phase and event.get("cartesian_matches") is not True:
                    raise ValueError("source PB::X/coordinate mismatch or absent pending first birth")
                sources.append(event["particle"])
                source_steps.add(event["step"])
            elif event["pending_first"]:
                if require_valid_phase and event["requested_s"] != event["first_duration_s"]:
                    raise ValueError("actual first mover duration differs from the birth record")
                first_moves.append(event["particle"])
    if not sources or len(sources) != len(set(sources)) or sorted(sources) != sorted(first_moves):
        raise ValueError("birth identity coverage is empty, duplicated, or lacks exactly one actual first movement")
    if source_steps != set(range(expected_steps)):
        raise ValueError("actual birth traces do not exercise every requested host step")
    checks.append({"status": "PASS", "message": "every source identity has exactly one first mover call", "births": len(sources)})
    return checks


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amps", type=Path, required=True)
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mpi-ranks", type=int, default=2)
    parser.add_argument("--mpi-exec", default="mpiexec")
    parser.add_argument("--timeout", type=float, default=120)
    args = parser.parse_args()
    amps, fixture, output = args.amps.resolve(), args.fixture.resolve(), args.output.resolve()
    try:
        values = load_fixture(fixture)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    if not amps.is_file() or args.mpi_ranks < 1 or not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("linked executable and positive MPI count/finite timeout are required")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        parser.error("choose a new or empty output directory; earlier evidence is retained")
    output.mkdir(parents=True, exist_ok=True)
    initial_hash, fixture_hash = binary_hash(amps), binary_hash(fixture)
    request = output / "fixture.args"
    request.write_text(fixture_protocol(values), encoding="utf-8")
    request_hash = binary_hash(request)
    preflights, results = [], []
    # Negative controls prove the PIC-private-state boundary before any host
    # job. An explicit protocol may not share a process or repeat a callback.
    for name, extra, expected_code, expected_message in (
            ("single-case", ["--test", "ARBN01"], 1,
             "end-to-end test input requires exactly one explicit --test ID"),
            ("repeat-one", ["--test-repeat", "2"], 2, "fresh process")):
        directory = output / ("preflight-" + name);directory.mkdir()
        work = create_native_work_directory(directory, "exclusive")
        command = [args.mpi_exec, "-n", str(args.mpi_ranks), str(amps), "--test", "HOST01",
                   "--test-input", str(request), "--test-output-dir", str(directory / "artifacts")] + extra
        log = directory / "native.log"
        try:
            code = run_role(command, work, log, args.timeout)
            text = log.read_text(encoding="utf-8")
            # The established parser rejects multiple explicit input IDs with
            # exit1 before campaign configuration. The additive native-only
            # repeat guard rejects with exit2. Pin both exact rejection paths
            # rather than accepting any early failure (e.g. a missing paired
            # output option) as proof of process isolation.
            passed = code == expected_code and expected_message in text and "The mesh is generated" not in text
            passed = (passed and binary_hash(amps) == initial_hash and binary_hash(fixture) == fixture_hash and
                      binary_hash(request) == request_hash)
            message = "explicit exclusive-process preflight rejected before mesh" if passed else "exclusive-process guard/exit/binary check failed"
        except (OSError, ValueError) as error:
            code, passed, message = 2, False, str(error)
        preflights.append({"name": name, "status": "PASS" if passed else "ERROR", "message": message,
                           "expected_exit_code": expected_code, "expected_diagnostic": expected_message,
                           "exit_code": code, "command": command, "log": str(log), "working_directory": str(work)})
    for identifier, mover in CASES.items():
        directory = output / identifier;directory.mkdir()
        work = create_native_work_directory(directory, identifier)
        report, log = directory / "srcsep-tests.json", directory / "native.log"
        command = [args.mpi_exec, "-n", str(args.mpi_ranks), str(amps), "--test", identifier,
                   "--test-input", str(request), "--test-output-dir", str(directory / "artifacts"),
                   "--test-json", str(report), "--test-junit", str(directory / "srcsep-tests.xml")]
        code, status, checks, message = None, "ERROR", [], "exclusive-process preflight failed"
        started = time.monotonic()
        if all(row["status"] == "PASS" for row in preflights):
            try:
                code = run_role(command, work, log, args.timeout)
                document = json.loads(report.read_text(encoding="utf-8"))
                rows = document.get("results", [])
                if len(rows) != 1 or rows[0].get("id") != identifier:
                    raise ValueError("native host report must contain exactly the requested case")
                native_status = rows[0].get("status")
                if (native_status, code) not in (("PASS", 0), ("FAIL", 1)):
                    raise ValueError("native PASS/FAIL must match expected zero/one exit; SKIP or exit contradictions are errors")
                metrics = {item.get("name"): item.get("value") for item in rows[0].get("metrics", [])}
                if (metrics.get("completed_host_steps") != values["steps"] or
                        not isinstance(metrics.get("global_first_birth_moves"), (int, float)) or metrics["global_first_birth_moves"] <= 0):
                    raise ValueError("native host metrics do not prove all requested actual steps and first births")
                failures = metrics.get("assertion_failures")
                if (type(failures) not in (int, float) or not math.isfinite(failures) or
                        (native_status == "PASS" and failures != 0) or (native_status == "FAIL" and failures <= 0)):
                    raise ValueError("native assertion metrics contradict the reported outcome")
                if (binary_hash(amps) != initial_hash or binary_hash(fixture) != fixture_hash or
                        binary_hash(request) != request_hash):
                    raise ValueError("linked executable, fixture or actual native protocol changed during the host campaign")
                # An actual assertion FAIL with exit1 is retained as FAIL after
                # complete immutable evidence/coverage validation. Its recorded
                # physical phase discrepancy must not be relabelled as a runner
                # error; missing traces or contradictory exits remain ERROR.
                checks = inspect_traces(directory / "artifacts", identifier, mover, args.mpi_ranks,
                                        values["steps"], require_valid_phase=native_status == "PASS")
                status = native_status
                message = ("actual host/PIC timing case and retained per-rank traces passed" if status == "PASS" else
                           "actual native assertions failed: " + str(rows[0].get("message", "")))
            except (OSError, ValueError, AttributeError, TypeError) as error:
                message = str(error)
                checks.append({"status": "FAIL", "message": message})
        results.append({"id": identifier, "status": status, "message": message,
                        "exit_code": code, "duration_s": time.monotonic() - started, "command": command,
                        "working_directory": str(work), "report": str(report), "log": str(log), "postchecks": checks})
        print(identifier, results[-1]["status"], message, flush=True)
    summary = {"schema": "srcsep-native-host-step-campaign-v1", "executable": str(amps),
               "executable_sha256": initial_hash, "fixture": str(fixture), "fixture_sha256": fixture_hash,
               "protocol_input": str(request), "protocol_input_sha256": request_hash,
               "mpi_ranks": args.mpi_ranks, "preflight_checks": preflights, "results": results,
               "scope": "actual-PIC-host-steps;controlled-fixed-analytic-source;production-physical-closure-unqualified"}
    temporary = output / "summary.json.tmp";temporary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    temporary.replace(output / "summary.json")
    return 0 if all(row["status"] == "PASS" for row in preflights + results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
