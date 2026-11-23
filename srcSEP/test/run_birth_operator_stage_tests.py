#!/usr/bin/env python3
"""Qualify native birth lifecycles with independent frozen-operator processes.

BIRTHN01–05 roles initialize one unchanged .01-second global timestep.
The direct control calls the existing native diagnostic operator with h=.0025;
it does not perform a shorter physical host step or sample a fractional interval.
Every role retains the same curved, nonzero-wind departure realization, and the
native checks compare complete particle/RNG images before and after the same
full material geometry transaction. BIRTHN06 uses the existing qualified reduced
checkpoint startup instead of pretending the diagnostic owner has a codec.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
from run_fresh_process_checkpoint_tests import binary_hash, create_native_work_directory, quoted, run_role
from run_native_coefficient_radius_tests import replace_key

SOURCE = Path(__file__).resolve().parents[1]
CASES = tuple(f"BIRTHN0{i}" for i in range(1, 7))


def analytic_fixture(template: str) -> tuple[str, dict[str, str]]:
    """Edit only explicit numerical/support-domain startup keys.

    The complete maintained equatorial Parker input supplies the curved field,
    actual400km/s radial wind, source/thermal metadata and all required schema6
    keys. Keep all4001 native vertices and its2e11m arc. The larger support cube
    contains that whole arc; an unrefined Eulerian support mesh does not change
    the segment-attached particle equations or any lifecycle test oracle.
    """
    changes = {"time_step_s": ".01", "maximum_level": "0", "outer_radius_m": "3.0e11"}
    text = template
    for key, value in changes.items():
        text = replace_key(text, key, value)
    # Refinement is inactive at level0; leave its complete declared profile
    # present so the shared scientific parser still validates the same record.
    return text, changes


def expected_case_seed(campaign: int, identifier: str, replicate: int) -> int:
    """Independent metadata oracle for the existing registry's seed schema1.

    This is NOT a random generator or a new stream law. It reproduces only
    DeriveTestCaseSeed's documented little-endian FNV serialization so a stale
    native report cannot claim that replicate1 used replicate0's case key.
    Future schema changes must update this explicit verification contract.
    """
    value = 14695981039346656037
    def append(payload: bytes) -> None:
        nonlocal value
        for byte in payload:
            value = ((value ^ byte) * 1099511628211) & ((1 << 64)-1)
    for word in (0x5445535443415345, 1, campaign):
        append(word.to_bytes(8, "little"))
    append(identifier.upper().encode("ascii"))
    append(replicate.to_bytes(8, "little"))
    return value


def validate_seed_metadata(native: dict, campaign: int | None, case: str, replicate: int) -> None:
    """Every seeded role must report the exact original campaign case key."""
    if campaign is None:
        return
    configuration = native.get("configuration", [])
    keys = [item for item in configuration if item.startswith("test_case_seed=")]
    schemas = [item for item in configuration if item.startswith("test_seed_schema=")]
    expected = "test_case_seed="+str(expected_case_seed(campaign, case, replicate))
    if keys != [expected] or schemas != ["test_seed_schema=1"]:
        raise ValueError("native seed metadata differs from shared registry campaign/ID/replicate schedule")


def selected_result(report: Path, case: str, code: int) -> dict:
    """A role needs an actual matching native PASS and zero process exit."""
    document = json.loads(report.read_text(encoding="utf-8"))
    rows = document.get("results", [])
    if len(rows) != 1 or rows[0].get("id") != case:
        raise ValueError("native role report lacks exactly its requested case")
    row = rows[0]
    status = row.get("status")
    if status not in ("PASS", "FAIL", "ERROR") or code != ({"PASS": 0, "FAIL": 1, "ERROR": 2}.get(status)):
        raise ValueError(f"contradictory native {status} exit={code}: {row.get('message')}")
    metrics = {metric.get("name"): metric.get("value") for metric in row.get("metrics", [])}
    if status == "PASS" and metrics.get("assertion_failures") != 0:
        raise ValueError("native assertion metric does not establish zero failures")
    if status == "FAIL" and not (isinstance(metrics.get("assertion_failures"), (int, float)) and metrics["assertion_failures"] > 0):
        raise ValueError("native FAIL does not retain a positive assertion count")
    return row


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amps", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mpi-ranks", type=int, default=2)
    parser.add_argument("--mpi-exec", default="mpiexec")
    parser.add_argument("--timeout", type=float, default=1200)
    parser.add_argument("--test-id", choices=CASES)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--repeat-count", "--test-repeat", type=int, default=1)
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    if args.mpi_ranks < 1 or args.repeat_count < 1 or args.timeout <= 0 or (args.seed is not None and not 0 <= args.seed < 2**64):
        parser.error("positive ranks/repeats/timeout and an optional uint64 seed are required")
    amps, output = args.amps.resolve(strict=True), args.output.resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        parser.error("choose a fresh or empty output directory")
    output.mkdir(parents=True, exist_ok=True)
    executable_hash = binary_hash(amps)
    runner_hash = binary_hash(Path(__file__).resolve())
    template = SOURCE / "examples/sep_parker_mesh.in"
    template_text = template.read_text(encoding="utf-8")
    fixture_text, changes = analytic_fixture(template_text)
    startup = output / "radial-operator-startup.in"
    startup.write_text(fixture_text, encoding="utf-8")
    startup_hash = binary_hash(startup)
    reduced_startup = SOURCE / "examples/reduced-shock/field_line.in"
    reduced_event = SOURCE / "examples/reduced-shock/positive_1au.event"
    # Pin the files ACTUALLY consumed by each physical startup. The analytic
    # owned deck is not evidence for BIRTHN06's independent reduced authority.
    protected_inputs = {str(path): binary_hash(path) for path in
                        (template, startup, reduced_startup, reduced_event)}
    results = []
    for case in ([args.test_id] if args.test_id else CASES):
        repetitions = []
        for repetition in range(args.repeat_count):
            directory = output / case / f"replicate-{repetition}"
            directory.mkdir(parents=True)
            prefix = directory / "birth"
            role_records = []
            # Each repetition is a NEW process group. Preserve the exact
            # registry DeriveTestCaseSeed(campaign,id,replicate) schedule:
            # C++ consumes the explicit replicate index before its first
            # TestSeed call. Never invent a seed+replicate RNG schedule.
            seed = args.seed if args.seed is not None else (0 if args.repeat_count > 1 else None)
            roles = ("direct", "full", "actual") if case <= "BIRTHN03" else (("routine",) if case == "BIRTHN06" else ("actual",))
            for role in roles:
                role_dir = directory / role
                role_dir.mkdir()
                work = create_native_work_directory(role_dir, role)
                request = role_dir / "fixture.args"
                if case == "BIRTHN06":
                    text = ("SEP_CHECKPOINT_FRESH_PROCESS 2\nrole \"routine\"\nprefix " + quoted(str(prefix)) +
                            "\nstartup_input " + quoted(str(reduced_startup)) +
                            "\nevent " + quoted(str(reduced_event)) + "\n")
                else:
                    text = ("SEP_BIRTH_OPERATOR_STAGE 1\nrole " + quoted(role) + "\nprefix " + quoted(str(prefix)) +
                            "\nstartup_input " + quoted(str(startup)) + "\n")
                text += "campaign_replicate_index " + str(repetition) + "\n"
                request.write_text(text, encoding="utf-8")
                report, log = role_dir / "srcsep-tests.json", role_dir / "native.log"
                command = [args.mpi_exec, "-n", str(args.mpi_ranks), "--bind-to", "none", str(amps),
                           "--test", case, "--test-input", str(request),
                           "--test-output-dir", str(role_dir / "artifacts"),
                           "--test-json", str(report), "--test-junit", str(role_dir / "srcsep-tests.xml")]
                if case != "BIRTHN06":
                    command.extend(("--shock-model", "analytical", "--turbulence-source", "prescribed",
                                    "--coefficient-source", "prescribed", "--no-particle-coupling",
                                    "--no-cascade", "--no-reflection"))
                if seed is not None:
                    command.extend(("--test-seed", str(seed)))
                code, status, message, native = None, "PREPARED", "", None
                started = time.monotonic()
                if not args.prepare_only:
                    try:
                        code = run_role(command, work, log, args.timeout)
                        native = selected_result(report, case, code)
                        if native["status"] == "PASS":
                            validate_seed_metadata(native, seed, case, repetition)
                        if binary_hash(Path(__file__).resolve()) != runner_hash:
                            raise ValueError("birth role runner source changed during campaign")
                        if binary_hash(amps) != executable_hash or any(binary_hash(Path(path)) != digest
                                                                       for path, digest in protected_inputs.items()):
                            raise ValueError("executable or actual scientific startup/template/event changed during campaign")
                        status = native["status"]
                        message = "" if status == "PASS" else str(native.get("message", ""))
                        if status == "PASS" and case <= "BIRTHN03":
                            for rank in range(args.mpi_ranks):
                                for phase in ("R0", "operator", "completed"):
                                    artifact = Path(str(prefix)+f".{role}.{phase}.rank{rank}.bin")
                                    if not artifact.is_file():
                                        raise ValueError("missing native immutable role artifact: "+str(artifact))
                    except (OSError, ValueError, TypeError) as error:
                        status, message = "ERROR", str(error)
                role_records.append({"role": role, "status": status, "exit_code": code, "message": message,
                                     "command": command, "working_directory": str(work), "report": str(report),
                                     "log": str(log), "duration_s": time.monotonic()-started, "native_result": native,
                                     "physical_input_hashes": {str(path): protected_inputs[str(path)] for path in
                                                               ((reduced_startup, reduced_event) if case == "BIRTHN06" else (startup,))}})
                if status in ("FAIL", "ERROR"):
                    break
            passed = len(role_records) == len(roles) and all(row["status"] == "PASS" for row in role_records)
            failed_status = "ERROR" if any(row["status"] == "ERROR" for row in role_records) else "FAIL"
            repetitions.append({"replicate": repetition, "seed": seed, "status": "PREPARED" if args.prepare_only else "PASS" if passed else failed_status, "roles": role_records})
        passed = all(row["status"] == "PASS" for row in repetitions)
        # Retain actual native metrics/config/seed, rather than replacing
        # them with a synthetic aggregate PASS. Independent reference roles
        # remain mandatory nested evidence for every requested repetition.
        actual_native = None
        if repetitions and repetitions[0]["roles"]:
            actual_native = repetitions[0]["roles"][-1].get("native_result")
        result = dict(actual_native or {})
        failed_status = "ERROR" if any(row["status"] == "ERROR" for row in repetitions) else "FAIL"
        result.update({"id": case, "status": "PREPARED" if args.prepare_only else "PASS" if passed else failed_status,
                  "message": "fresh-process operator-stage lifecycle gate; unchanged fixed physical dt; real completed material map" if passed else "see independent role reports",
                  "repetitions": repetitions})
        result.setdefault("configuration", []).append("requested_fresh_process_repetitions="+str(args.repeat_count))
        actual_seeds = [row["roles"][-1].get("native_result", {}).get("seed")
                        for row in repetitions if row["roles"][-1].get("native_result")]
        # Match Registry::RunCampaign: an aggregate with different invocation
        # seeds has no single seed. Each original native invocation retains its
        # own exact seed/metrics; neither the first seed nor its duration stands
        # in for the complete requested multi-process campaign.
        result["seed"] = actual_seeds[0] if actual_seeds and len(set(actual_seeds)) == 1 else None
        result["elapsed_seconds"] = sum(role["duration_s"] for row in repetitions for role in row["roles"])
        result["invocations"] = [{"replicate_index": row["replicate"], "result": row["roles"][-1].get("native_result") or
                                 {"id": case, "status": row["status"], "message": "independent role failed before actual"}}
                                for row in repetitions]
        results.append(result)
        print(result["status"], case, flush=True)
    summary = {"schema": "srcsep-birth-operator-stage-v1", "executable": str(amps), "executable_sha256": executable_hash,
               "runner": str(Path(__file__).resolve()), "runner_sha256": runner_hash,
               "template": str(template), "template_sha256": hashlib.sha256(template_text.encode()).hexdigest(),
               "owned_startup": str(startup), "owned_startup_sha256": startup_hash, "explicit_changes": changes,
               "protected_input_hashes": protected_inputs,
               "mpi_ranks": args.mpi_ranks, "repeat_count": args.repeat_count,
               "scope": "frozen-departure-native-operator-parity-and-completed-material-publication;no-new-source-physics",
               "results": results}
    temporary = output / "summary.json.tmp"
    temporary.write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    temporary.replace(output / "summary.json")
    # Standard native-runner surface consumed by the existing routine and
    # orchestrator readers: exactly the requested stable IDs, no missing/Skip
    # records, and original actual-role metrics/configuration retained.
    standard = {"schema": "srcsep-component-tests-v1", "results": results,
                "exit_code": 2 if any(row["status"] == "ERROR" for row in results) else 1 if any(row["status"] == "FAIL" for row in results) else 0,
                "not_run_invocation_count": 0, "not_run": [],
                "totals": {key: sum(row["status"] == status for row in results)
                           for key, status in (("passed", "PASS"), ("failed", "FAIL"),
                                               ("skipped", "SKIP"), ("errors", "ERROR"))}}
    (output / "srcsep-tests.json").write_text(json.dumps(standard, indent=2)+"\n", encoding="utf-8")
    suite = ET.Element("testsuite", name="native-birth-operator-stage", tests=str(len(results)),
                       failures=str(sum(row["status"] == "FAIL" for row in results)), errors=str(sum(row["status"] == "ERROR" for row in results)))
    for row in results:
        case_xml = ET.SubElement(suite, "testcase", name=row["id"], classname="birth-window-native")
        if row["status"] not in ("PASS", "PREPARED"):
            ET.SubElement(case_xml, "failure" if row["status"] == "FAIL" else "error", message=row["message"]).text=json.dumps(row["repetitions"])
    ET.ElementTree(suite).write(output / "srcsep-tests.xml", encoding="utf-8", xml_declaration=True)
    return 0 if args.prepare_only else standard["exit_code"]


if __name__ == "__main__":
    raise SystemExit(main())
