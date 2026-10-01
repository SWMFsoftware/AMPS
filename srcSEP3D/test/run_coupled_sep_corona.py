#!/usr/bin/env python3
"""Run all shared SEP/corona gates and every live native sep-corona descriptor.

From the AMPS root:
  python3 srcSEP3D/test/run_coupled_sep_corona.py --amps ./amps --ranks 4 --test-input srcSEP3D/examples/sep3d_analytic_parker_active_tube.in --test-steps 0

The shared runner includes Stage 11 C++ and Stage 12 offline Python tests.
The linked executable owns native suite membership. Counts/IDs are discovered
from both authorities on every invocation; this script maintains no copied
test catalogue. Results retain their scope: shared physics/preprocessing
verification is distinct from live AMPS/MPI initialization evidence.
Child output streams live. Phase/test counters and elapsed-time heartbeats
use --progress-interval (15 seconds by default); JSON remains authoritative.
The final failure summary identifies cases, diagnostics and fresh log paths.
"""
from __future__ import annotations

import argparse
import codecs
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import selectors
import shlex
import signal
import subprocess
import sys
import time
import uuid
import xml.etree.ElementTree as ET

AMPS_ROOT = Path(__file__).resolve().parents[2]
STATUSES = ("PASS", "FAIL", "SKIP", "ERROR")


class SuiteError(RuntimeError):
    """Invalid/incomplete evidence, separate from a scientific FAIL."""


def require(condition, message):
    if not condition:
        raise SuiteError(message)


def read_json(path):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, "duplicate JSON key: " + key)
            result[key] = value
        return result
    def nonfinite(value):
        raise SuiteError("nonfinite JSON value: " + value)
    def floating(value):
        number = float(value)
        require(math.isfinite(number), "out-of-range JSON number: " + value)
        return number
    try:
        return json.loads(Path(path).read_text(), object_pairs_hook=pairs,
                          parse_constant=nonfinite, parse_float=floating)
    except (OSError, ValueError) as error:
        raise SuiteError("missing/malformed fresh report: " + str(error)) from error


def sha256(path):
    checksum = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


class PhaseProgress:
    """Console observations are progress, never a substitute for checked JSON.

    Only IDs in the discovered registry count toward completion, and repeated
    result messages count once. A failed test still completes an execution.
    The final report readers, rather than these stdout messages, own PASS/FAIL.
    """
    def __init__(self, phase, started, registry):
        self.phase, self.started = phase, started
        self.ids = None if registry is None else {row["id"] for row in registry}
        self.seen = {}
        self.pending = ""

    def update(self, text, final=False):
        if self.ids is None:
            return
        self.pending += text
        lines = self.pending.split("\n")
        self.pending = "" if final else lines.pop()
        for line in lines:
            matched = re.match(r"^\[([A-Z][A-Z0-9]+)\]\s+(PASS|FAIL|SKIP|ERROR)(?:\s|$)", line)
            if matched and matched[1] in self.ids and matched[1] not in self.seen:
                self.seen[matched[1]] = matched[2]
                self.emit("result="+matched[1]+" "+matched[2])
        # Non-result output can have very long lines (compiler diagnostics or
        # progress dots). Keep a bounded tail; it cannot be a canonical result
        # line if its beginning exceeded this bound without matching a record.
        if len(self.pending) > 65536:
            self.pending = self.pending[-65536:]

    def emit(self, message):
        count = ""
        if self.ids is not None:
            count = " completed="+str(len(self.seen))+"/"+str(len(self.ids))
            if self.ids:
                count += " ("+format(100.0*len(self.seen)/len(self.ids), ".1f")+"%)"
        print("[progress] "+self.phase+count+" elapsed="+
              format(time.monotonic()-self.started, ".1f")+"s "+message, flush=True)


def stop_process(process):
    """Reap a timed-out launcher and its local children, including pipe writers.

    The child owns a new POSIX session. Killing only the launcher could leave a
    compiler/rank holding stdout open forever after the execution deadline.
    Remote MPI cleanup is the launcher's responsibility; its loss is an ERROR.
    """
    if process is None:
        return
    try:
        if os.name == "posix":
            os.killpg(process.pid, signal.SIGKILL)
        elif process.poll() is None:
            process.kill()
    except ProcessLookupError:
        pass
    process.wait(timeout=5)


def execute(command, cwd, log, timeout, operations, phase, progress_interval,
            registry=None, stream=True):
    """Stream subprocess output and retain diagnostics even on launch failure.

    Each phase has a distinct fresh log. No success is inferred from exit zero:
    report validation below must also recover every discovered test exactly once.
    A selectable pipe gives live output plus heartbeats even when the child is
    silent. Raw child bytes are written unchanged to the log. Python children
    inherit unbuffered stdout so each completed shared test is visible at once.
    """
    started = time.monotonic()
    log.parent.mkdir(parents=True, exist_ok=True)
    error = None
    process = None
    progress = PhaseProgress(phase, started, registry)
    progress.emit("START log="+str(log))
    print("RUN " + " ".join(shlex.quote(part) for part in command), flush=True)
    decoder = codecs.getincrementaldecoder("utf-8")(errors="replace")
    selector = selectors.DefaultSelector()
    last_output, next_progress = started, started+progress_interval
    try:
        environment = dict(os.environ, PYTHONUNBUFFERED="1")
        with log.open("wb") as output:
            process = subprocess.Popen(command, cwd=str(cwd), stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, env=environment, bufsize=0,
                start_new_session=os.name == "posix")
            selector.register(process.stdout, selectors.EVENT_READ)
            while process.poll() is None or selector.get_map():
                now = time.monotonic()
                if now-started >= timeout:
                    raise subprocess.TimeoutExpired(command, timeout)
                wait = max(0.0, min(0.2, next_progress-now, timeout-(now-started)))
                for key, _ in selector.select(wait):
                    chunk = os.read(key.fileobj.fileno(), 65536)
                    if chunk:
                        output.write(chunk); output.flush()
                        text = decoder.decode(chunk)
                        if stream:
                            print(text, end="", flush=True)
                        progress.update(text)
                        last_output = time.monotonic()
                    else:
                        selector.unregister(key.fileobj)
                        tail = decoder.decode(b"", final=True)
                        if stream and tail:
                            print(tail, end="", flush=True)
                        progress.update(tail, final=True)
                now = time.monotonic()
                if now >= next_progress:
                    progress.emit("RUNNING output_idle="+format(now-last_output, ".1f")+"s")
                    next_progress = now+progress_interval
            code = process.wait()
    except (OSError, subprocess.TimeoutExpired) as failure:
        code, error = 2, str(failure)
        stop_process(process)
    except BaseException:
        stop_process(process)
        raise
    finally:
        selector.close()
        if process is not None and process.stdout is not None:
            process.stdout.close()
    progress.emit("END exit="+str(code)+( " ERROR="+error if error else ""))
    text = log.read_text(errors="replace") if log.exists() else ""
    operation = {"command": command, "returncode": code,
                 "seconds": time.monotonic()-started, "log": str(log), "error": error,
                 "phase": phase, "console_results_seen": len(progress.seen)}
    operations.append(operation)
    return operation, text


def discover_shared(text):
    records = []
    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) != 3 or not fields[1].startswith("stage="):
            continue
        require(re.fullmatch(r"[A-Z][A-Z0-9]+", fields[0]), "invalid shared test ID")
        stage = int(fields[1][6:])
        require(stage >= 0 and fields[2], "invalid shared stage/backend")
        records.append({"id": fields[0], "stage": stage, "backend": fields[2]})
    require(records and len({r["id"] for r in records}) == len(records),
            "empty/duplicate shared test registry")
    return records


def discover_native(text):
    records, seen = [], set()
    for line in text.splitlines():
        fields = [field.strip() for field in line.split("|")]
        if len(fields) != 4 or not fields[3].startswith("suite="):
            continue
        require(re.fullmatch(r"[A-Z][A-Z0-9]+", fields[0]) and fields[0] not in seen,
                "invalid/duplicate native registry descriptor")
        seen.add(fields[0])
        if fields[3] == "suite=sep-corona":
            records.append({"id": fields[0], "name": fields[1]})
    require(records, "executable advertises no native sep-corona suite")
    return records


def exact_rows(rows, registry):
    require(isinstance(rows, list) and all(isinstance(r, dict) for r in rows),
            "report results must be an array of records")
    ids = [r.get("id") for r in rows]
    require(all(isinstance(i, str) for i in ids) and len(set(ids)) == len(ids) and
            set(ids) == {r["id"] for r in registry},
            "report omits, duplicates or adds a discovered test")
    return {r["id"]: r for r in rows}


def read_shared(report, registry, code):
    data = read_json(report)
    require(data.get("suite") == "sep_coronal_cme", "wrong shared report suite")
    known = exact_rows(data.get("tests"), registry)
    rows = []
    for test in registry:
        row = known[test["id"]]
        require(type(row.get("passed")) is bool and row.get("stage") == test["stage"] and
                isinstance(row.get("seconds"), (int, float)) and math.isfinite(row["seconds"]) and row["seconds"] >= 0,
                "invalid shared result/stage/duration")
        status=row.get("status","PASS" if row["passed"] else "FAIL")
        require(status in {"PASS","FAIL","SKIP"} and row["passed"]==(status=="PASS"),"shared status/boolean disagreement")
        rows.append(dict(row, scope="shared-model", status=status, executed=status!="SKIP"))
    passed = sum(r["status"] == "PASS" for r in rows)
    failed = sum(r["status"] == "FAIL" for r in rows);skipped=sum(r["status"]=="SKIP" for r in rows)
    require(data.get("total") == len(rows) and data.get("passed") == passed and
            data.get("failed") == failed and data.get("skipped",0)==skipped, "shared report totals disagree with rows")
    require(code == (1 if failed or (data.get("require_no_skips",False) and skipped) else 0), "shared report/process exit disagree")
    return rows


def read_native(report, registry, code, ranks):
    data = read_json(report)
    require(data.get("schema") == "srcsep-component-tests-v1", "unsupported native report schema")
    require(type(data.get("mpi_ranks")) is int and data["mpi_ranks"] == ranks,
            "native report MPI rank count differs from requested count")
    require(isinstance(data.get("configuration_fingerprint"), str) and data["configuration_fingerprint"],
            "native report has no configured-state identity")
    providers = data.get("providers", {})
    require(type(providers.get("input_schema_version")) is int and providers["input_schema_version"] > 0 and
            all(isinstance(providers.get(k), str) and providers[k] not in {"", "unreported"} for k in ("background", "shock")),
            "native provider identity missing; rebuild amps with the updated native report boundary")
    require(type(data.get("completed_steps")) is int and data["completed_steps"] >= 0,
            "native completed-step evidence missing")
    known = exact_rows(data.get("results"), registry)
    rows = []
    for test in registry:
        row = known[test["id"]]
        require(row.get("status") in STATUSES and isinstance(row.get("metrics"), list) and
                isinstance(row.get("artifacts"), list), "invalid native status/metrics/artifacts")
        rows.append(dict(row, scope="native-amps", executed=True))
    expected = 2 if any(r["status"] == "ERROR" for r in rows) else (1 if any(r["status"] == "FAIL" for r in rows) else 0)
    # A suite exit represents ALL rows. A PASS beside another case's FAIL is
    # still a valid PASS; comparing each row separately with the suite exit
    # would incorrectly turn every successful row into ERROR.
    require(code == expected, "native suite report/process exit disagree")
    return rows, {"providers": providers, "configuration_fingerprint": data["configuration_fingerprint"],
                  "mpi_ranks": ranks, "completed_steps": data["completed_steps"]}


def failure_summary(run, results, errors, operations, require_no_skips):
    """Make checked failures reviewable without searching the streamed output.

    Shared cases retain their full captured output in a dedicated diagnostic
    file. Native cases expose their captured message, metrics and artifacts;
    their MPI stdout is the suite log, not separable per-case output. Console
    reasons are bounded, but neither diagnostic files nor raw logs are clipped.
    Include infrastructure errors even when discovery yielded no case IDs.
    """
    bad = [row for row in results if row["status"] in {"FAIL", "ERROR"}]
    strict_skips = [row for row in results if row["status"] == "SKIP"] if require_no_skips else []
    counts = Counter(row["status"] for row in bad)
    lines = ["failed_test_summary: fail="+str(counts["FAIL"])+" error="+str(counts["ERROR"])+
             " runner_errors="+str(len(errors))]
    if strict_skips:
        lines.append("strict_skip_summary: skip="+str(len(strict_skips))+" (--require-no-skips)")
    # The actual executed phases own log paths. Discovery/build failures may
    # have no execution log; never imply that such a file was generated.
    logs = {operation["phase"]: operation["log"] for operation in operations
            if Path(operation["log"]).is_file()}
    execution_phase = {"shared-model": "shared-model-tests", "native-amps": "native-amps-tests"}
    for row in bad+strict_skips:
        captured = str(row.get("output") or "")
        message = str(row.get("message") or "")
        # Assertion/traceback conclusions usually appear at the end of shared
        # output. The full body stays available if this short reason is insufficient.
        reason = message or next((line.strip() for line in reversed(captured.splitlines()) if line.strip()),
                                 "No diagnostic message recorded; inspect the phase log and source report.")
        reason = " ".join(reason.split())
        if len(reason) > 500:
            reason = reason[:497]+"..."
        directory = run/"failures"/row["scope"]
        directory.mkdir(parents=True, exist_ok=True)
        diagnostic = directory/(row["id"]+".log")
        row["diagnostic_log"] = str(diagnostic)
        phase_log = logs.get(execution_phase[row["scope"]])
        contents = [row["status"]+" "+row["scope"]+"/"+row["id"],
                    "run_directory="+str(run), "phase_log="+(phase_log or "not generated")]
        if message:
            contents.extend(["", "Message:", message])
        if captured:
            contents.extend(["", "Captured test output:", captured])
        if "metrics" in row or "artifacts" in row:
            contents.extend(["", "Native metrics/artifact references:",
                             json.dumps({key: row[key] for key in ("metrics", "artifacts") if key in row}, indent=2)])
        diagnostic.write_text("\n".join(contents)+"\n", encoding="utf-8")
        lines.extend(["  "+row["status"]+" "+row["scope"]+"/"+row["id"]+": "+reason,
                      "    diagnostic_log="+str(diagnostic)])
        if phase_log:
            lines.append("    execution_log="+phase_log)
    for error in errors:
        lines.append("  runner_error "+error["scope"]+": "+" ".join(error["message"].split()))
    lines.append("test_logs_directory="+str(run))
    for phase, path in logs.items():
        lines.append("phase_log "+phase+"="+path)
    for scope, path in (("shared-model", run/"shared/results.json"), ("native-amps", run/"native/native.json")):
        if path.is_file():
            lines.append("source_report "+scope+"="+str(path))
    lines.append("failure_summary="+str(run/"failures.txt"))
    return "\n".join(lines)+"\n"


def write_reports(output, report, failure_text):
    output.mkdir(parents=True, exist_ok=True)
    temporary = output/("summary.json."+uuid.uuid4().hex+".tmp")
    temporary.write_text(json.dumps(report, indent=2)+"\n")
    os.replace(temporary, output/"summary.json")
    cases = report["results"]
    # An unknown/failed registry has no scientific rows to mark ERROR. Expose
    # its infrastructure error explicitly in JUnit so XML consumers cannot
    # mistake the surviving phase's passing cases for a complete suite.
    infrastructure = [error for error in report["runner_errors"] if not any(
        row["scope"] == error["scope"] and row["status"] == "ERROR" for row in cases)]
    suite = ET.Element("testsuite", name="sep-corona-all", tests=str(len(cases)+len(infrastructure)),
                       failures=str(report["counts"]["fail"]), errors=str(report["counts"]["error"]+len(infrastructure)),
                       skipped=str(report["counts"]["skip"]))
    for row in cases:
        case = ET.SubElement(suite, "testcase", name=row["id"], classname=row["scope"],
                             time=str(row.get("seconds", 0)))
        if row["status"] != "PASS":
            tag = {"FAIL": "failure", "ERROR": "error", "SKIP": "skipped"}[row["status"]]
            ET.SubElement(case, tag, message=str(row.get("message", ""))).text = row.get("output", "")
    for error in infrastructure:
        case = ET.SubElement(suite, "testcase", name="runner-discovery", classname="runner."+error["scope"])
        ET.SubElement(case, "error", message=error["message"])
    temporary = output/("junit.xml."+uuid.uuid4().hex+".tmp")
    ET.ElementTree(suite).write(temporary, encoding="utf-8", xml_declaration=True)
    os.replace(temporary, output/"junit.xml")
    # Keep a convenient latest text summary as well as the immutable run copy.
    temporary = output/("failures.txt."+uuid.uuid4().hex+".tmp")
    temporary.write_text(failure_text, encoding="utf-8")
    os.replace(temporary, output/"failures.txt")


def positive(value):
    result = int(value)
    if result <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--amps", type=Path, default=AMPS_ROOT/"amps", help="rebuilt linked srcSEP3D executable")
    parser.add_argument("--test-input", type=Path, default=AMPS_ROOT/"srcSEP3D/examples/sep3d_analytic_parker_active_tube.in")
    parser.add_argument("--ranks", type=positive, default=4)
    parser.add_argument("--test-steps", type=int, default=0, help="native horizon; zero checks initialization only")
    parser.add_argument("--launcher", default="mpiexec -n {ranks}", help="MPI argv template, e.g. 'srun -n {ranks}'")
    parser.add_argument("--jobs", type=positive, default=8)
    parser.add_argument("--timeout", type=float, default=1800, help="seconds allowed per phase")
    parser.add_argument("--progress-interval", type=float, default=15,
                        help="seconds between running/elapsed heartbeats (default: 15); child output and completed-test updates stream live")
    parser.add_argument("--output-dir", type=Path, default=AMPS_ROOT/"test_output/coupled-sep-corona")
    parser.add_argument("--model-dir", type=Path, default=AMPS_ROOT/"src/models/sep_coronal_cme")
    parser.add_argument("--model-only", action="store_true", help="explicitly omit live AMPS; report shared-only scope")
    parser.add_argument("--list", action="store_true", help="discover both registries without building/running cases")
    parser.add_argument("--require-no-skips", action="store_true", help="return failure if any selected case is SKIP")
    args = parser.parse_args(argv)
    parser.error("--test-steps must be nonnegative") if args.test_steps < 0 else None
    parser.error("--timeout must be positive and finite") if not math.isfinite(args.timeout) or args.timeout <= 0 else None
    parser.error("--progress-interval must be positive and finite") if not math.isfinite(args.progress_interval) or args.progress_interval <= 0 else None
    model, executable, deck, output = (p.resolve() for p in (args.model_dir, args.amps, args.test_input, args.output_dir))
    # Fresh run-owned paths prevent stale reports from satisfying a crashed
    # invocation, while a latest summary remains convenient for CI/operators.
    run = output/"runs"/(datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"-"+uuid.uuid4().hex[:12])
    run.mkdir(parents=True, exist_ok=False)
    operations, results, errors, registries, native_state, ownership = [], [], [], {}, {}, {}
    runner = model/"test/run_tests.py"
    for scope in ("shared-model", "native-amps"):
        if scope == "native-amps" and args.model_only:
            continue
        registry = []
        try:
            if scope == "shared-model":
                operation, text = execute([sys.executable, str(runner), "--list"], AMPS_ROOT,
                                          run/"shared-discovery.log", args.timeout, operations,
                                          "shared-discovery", args.progress_interval, stream=False)
                require(operation["returncode"] == 0, "shared registry discovery failed: " + text[-1500:])
                registry = discover_shared(text)
                registries[scope] = registry
                print("[progress] shared registry selected="+str(len(registry))+" tests", flush=True)
                if args.list:
                    continue
                build, text = execute(["make", "-C", str(model), "-j"+str(args.jobs),
                                       "build/sep_coronal_cme_tests", "check-adapters"], AMPS_ROOT,
                                      run/"shared-build.log", args.timeout, operations,
                                      "shared-build", args.progress_interval)
                require(build["returncode"] == 0, "shared build failed: " + text[-1500:])
                operation, text = execute([sys.executable, str(runner), "--all", "--output-dir", str(run/"shared")],
                                          AMPS_ROOT, run/"shared.log", args.timeout, operations,
                                          "shared-model-tests", args.progress_interval, registry)
                require(operation["error"] is None, str(operation["error"]))
                results.extend(read_shared(run/"shared/results.json", registry, operation["returncode"]))
                print("[progress] shared-model report verified="+str(len(registry))+" cases", flush=True)
            else:
                require(executable.is_file() and os.access(executable, os.X_OK), "linked --amps executable is missing/nonexecutable")
                ownership["executable_sha256"] = sha256(executable)
                operation, text = execute([str(executable), "--list-tests"], AMPS_ROOT,
                                          run/"native-discovery.log", args.timeout, operations,
                                          "native-discovery", args.progress_interval, stream=False)
                require(operation["returncode"] == 0, "native registry discovery failed: " + text[-1500:])
                registry = discover_native(text)
                registries[scope] = registry
                print("[progress] native registry selected="+str(len(registry))+" tests; aggregate selected="+
                      str(sum(len(rows) for rows in registries.values())), flush=True)
                if args.list:
                    continue
                require(deck.is_file(), "native --test-input is missing")
                ownership["input_sha256"] = sha256(deck)
                require("{ranks}" in args.launcher, "--launcher must contain {ranks}")
                launcher = shlex.split(args.launcher.format(ranks=args.ranks))
                require(launcher, "empty MPI launcher")
                command = launcher+[str(executable), "--test-suite", "sep-corona", "--test-input", str(deck),
                    "--test-steps", str(args.test_steps), "--expect-mpi-ranks", str(args.ranks),
                    "--test-json", str(run/"native/native.json"), "--artifact-directory", str(run/"native/artifacts")]
                operation, text = execute(command, AMPS_ROOT, run/"native.log", args.timeout, operations,
                                          "native-amps-tests", args.progress_interval, registry)
                require(operation["error"] is None, str(operation["error"]))
                require(sha256(executable) == ownership["executable_sha256"] and sha256(deck) == ownership["input_sha256"],
                        "executable/input bytes changed during native invocation")
                rows, native_state = read_native(run/"native/native.json", registry, operation["returncode"], args.ranks)
                results.extend(rows)
                print("[progress] native-amps report verified="+str(len(registry))+" cases", flush=True)
        except (SuiteError, OSError, ValueError, TypeError, KeyError) as error:
            message = str(error)
            errors.append({"scope": scope, "message": message})
            # These ERROR records describe invalid/missing execution evidence,
            # not failed physics. Known cases cannot quietly disappear from
            # the aggregate when their build/launcher/report fails.
            if not args.list:
                results.extend(dict(test, scope=scope, status="ERROR", executed=False, message=message) for test in registry)
            print(scope+" ERROR: "+message, file=sys.stderr, flush=True)
    if args.list:
        for scope, registry in registries.items():
            for row in registry:
                print(scope+"\t"+row["id"]+("\tstage="+str(row["stage"]) if "stage" in row else ""))
        print("DISCOVERY: "+", ".join(scope+"="+str(len(rows)) for scope, rows in registries.items()))
        return 2 if errors else 0
    counts = Counter(row["status"].lower() for row in results)
    totals = dict(total=len(results), **{status.lower(): counts[status.lower()] for status in STATUSES})
    code = 2 if errors or counts["error"] else (1 if counts["fail"] or (args.require_no_skips and counts["skip"]) else 0)
    failure_text = failure_summary(run, results, errors, operations, args.require_no_skips)
    report = {"schema": "sep-corona-aggregate-v1", "run_directory": str(run),
              "evidence_kind": "shared-software-and-generic-native-host",
              "production_release_qualified": False,
              "observational_campaign_qualified": False,
              "scope": "shared-model-only" if args.model_only else "shared-model-and-native-amps",
              "counts": totals, "registries": registries, "results": results, "runner_errors": errors,
              "operations": operations, "ownership": ownership, "native_state": native_state,
              "requested_native_steps": None if args.model_only else args.test_steps,
              "require_no_skips": args.require_no_skips, "exit_code": code,
              "failure_summary": str(run/"failures.txt"),
              "coronal_runtime_provider_qualification": "not-established-by-generic-native-initialization"}
    print("[progress] writing aggregate JSON/JUnit reports; total="+str(len(results)), flush=True)
    write_reports(run, report, failure_text)
    write_reports(output, report, failure_text)
    print("coupled_test_summary: "+" ".join(k+"="+str(v) for k, v in totals.items())+
          " runner_errors="+str(len(errors))+" scope="+report["scope"], flush=True)
    print(failure_text, end="", flush=True)
    print("aggregate_json="+str(output/"summary.json"), flush=True)
    print("aggregate_junit="+str(output/"junit.xml"), flush=True)
    if native_state:
        print("native_providers="+json.dumps(native_state["providers"], sort_keys=True), flush=True)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
