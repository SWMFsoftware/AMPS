#!/usr/bin/env python3
"""Bounded native MPI harness checks; none is a physical-model validation."""

import argparse
import json
import pathlib
import subprocess
import tempfile
import traceback
import xml.etree.ElementTree as ET


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amps", type=pathlib.Path, required=True)
    parser.add_argument("--mpiexec", default="mpiexec")
    parser.add_argument("--output", type=pathlib.Path)
    args = parser.parse_args()
    executable = args.amps.resolve(strict=True)
    root = (args.output or pathlib.Path(tempfile.mkdtemp(prefix="srcsep-native-cli."))).resolve()
    root.mkdir(parents=True, exist_ok=True)
    records = []
    planned = ["seed-A", "seed-reordered", "seed-B-keep", "rank1-fault", "invalid-fault",
               "fixed-contract", "unknown-selector", "output-preflight", "rank1-artifact-error"]
    active = None

    def finish():
        if active is not None and active["status"] == "RUNNING":
            active["status"] = "PASS"
            active["message"] = "expected process outcome and all recorded postchecks passed"

    def check(condition, message):
        if active is not None:
            active["postchecks"].append({"status": "PASS" if condition else "FAIL", "message": message})
        if not condition:
            raise AssertionError(message)

    def run(name, arguments, expected=0):
        nonlocal active
        finish()
        command = [args.mpiexec, "-n", "2", str(executable)] + arguments
        active = {"id": "NCLI%02d" % (planned.index(name) + 1), "name": name,
                  "status": "RUNNING", "message": "process/postchecks in progress",
                  "command": command, "exit_code": None, "expected_exit_code": expected,
                  "log": str(root / (name + ".log")), "postchecks": []}
        records.append(active)
        # Every MPI job owns its fixed PIC scratch names, even though these
        # controlled infrastructure cases normally return before mesh setup.
        # Unique leaves also make two callers reusing an output parent safe.
        work = pathlib.Path(tempfile.mkdtemp(prefix=name+"-", dir=str(root)))
        active["working_directory"] = str(work)
        try:
            result = subprocess.run(command, cwd=str(work), stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, text=True, timeout=90)
        except (OSError, subprocess.TimeoutExpired) as error:
            captured = getattr(error, "stdout", "") or ""
            if isinstance(captured, bytes):
                captured = captured.decode("utf-8", errors="replace")
            (root / (name + ".log")).write_text(captured + "\n" + str(error))
            raise
        (root / (name + ".log")).write_text(result.stdout)
        active["exit_code"] = result.returncode
        check(result.returncode == expected, "%s actual exit=%s expected=%s; see %s" % (
            name, result.returncode, expected, root / (name + ".log")))
        return result

    def summary(name):
        result = json.loads((root / name / "results.json").read_text())
        check(result["schema"] == "srcsep-component-tests-v1", 'postcheck: result["schema"] == "srcsep-component-tests-v1"')
        ET.parse(str(root / name / "results.xml"))
        return result

    def draws(name, count=2):
        result = []
        for replicate in range(count):
            directory = root / name / "CLIEX01" / ("replicate-%d" % replicate)
            pair = [json.loads((directory / suffix / "CLIEX01_draw.json").read_text())
                    for suffix in (".", "rank-1")]
            check([item["rank"] for item in pair] == [0, 1], 'postcheck: [item["rank"] for item in pair] == [0, 1]')
            check(all(item["replicate"] == replicate for item in pair), 'postcheck: all(item["replicate"] == replicate for item in pair)')
            for suffix in (".", "rank-1"):
                check((directory / suffix / "CLIEX01_immutable_reference.txt").read_text() == \
                    "CLIEX01 immutable reference sentinel\n", 'postcheck: (directory / suffix / "CLIEX01_immutable_reference.txt").read_text() == \\ "CLIEX01 immutable reference sentinel\\n"')
            # This infrastructure stream intentionally has one physical key on
            # both ranks; the rank is only the writer/fault provenance.
            for key in ("seed", "uniform", "normal", "integer_state"):
                check(pair[0][key] == pair[1][key], 'postcheck: pair[0][key] == pair[1][key]')
            result.append(tuple(pair[0][key] for key in ("seed", "uniform", "normal", "integer_state")))
        return result

    try:
        run("seed-A", ["--test", "CLIEX01", "--test-seed", "123", "--test-repeat", "2",
                       "--test-format", "json", "--test-output", str(root / "seed-A")])
        first = summary("seed-A")
        check(len(first["results"]) == 1 and first["results"][0]["id"] == "CLIEX01", 'postcheck: len(first["results"]) == 1 and first["results"][0]["id"] == "CLIEX01"')
        executions = first["results"][0]["invocations"]
        check(len(executions) == 2 and [item["replicate_index"] for item in executions] == [0, 1], 'postcheck: len(executions) == 2 and [item["replicate_index"] for item in executions] == [0, 1]')
        for invocation in executions:
            child = invocation["result"]
            check(child["status"] == "PASS", 'postcheck: child["status"] == "PASS"')
            for rank in (0, 1):
                check(any(value.startswith("mpi_rank_%d=" % rank) for value in child["configuration"]), 'postcheck: any(value.startswith("mpi_rank_%d=" % rank) for value in child["configuration"])')
        a = draws("seed-A")
        check(a[0][0] != a[1][0] and a[0][1:3] != a[1][1:3], 'postcheck: a[0][0] != a[1][0] and a[0][1:3] != a[1][1:3]')
        check(not list((root / "seed-A").rglob(".cli-probe-unlisted.tmp")), 'postcheck: not list((root / "seed-A").rglob(".cli-probe-unlisted.tmp"))')
        check(not list((root / "seed-A").rglob("scratch.*")), 'postcheck: not list((root / "seed-A").rglob("scratch.*"))')

        run("seed-reordered", ["--tests", "FTED02,CLIEX01", "--test-seed=123", "--test-repeat=2",
                               "--test-output", str(root / "seed-reordered")])
        check(draws("seed-reordered") == a, 'postcheck: draws("seed-reordered") == a')
        check(summary("seed-reordered")["totals"]["passed"] == 2, 'postcheck: summary("seed-reordered")["totals"]["passed"] == 2')

        run("seed-B-keep", ["--test", "CLIEX01", "--test-seed", "456", "--test-repeat", "2",
                            "--test-format", "junit", "--keep-test-artifacts=on",
                            "--test-output", str(root / "seed-B-keep")])
        b = draws("seed-B-keep")
        check(all(x[0] != y[0] and x[1:3] != y[1:3] for x, y in zip(a, b)), 'postcheck: all(x[0] != y[0] and x[1:3] != y[1:3] for x, y in zip(a, b))')
        check(len(list((root / "seed-B-keep").rglob(".cli-probe-unlisted.tmp"))) == 4, 'postcheck: len(list((root / "seed-B-keep").rglob(".cli-probe-unlisted.tmp"))) == 4')
        junit = ET.parse(str(root / "seed-B-keep" / "results.xml")).getroot()
        check(junit.attrib["tests"] == "2", 'postcheck: junit.attrib["tests"] == "2"')
        check([item.attrib["name"] for item in junit.findall("testcase")] == [
            "CLIEX01[replicate=0]", "CLIEX01[replicate=1]"], 'postcheck: [item.attrib["name"] for item in junit.findall("testcase")] == [ "CLIEX01[replicate=0]", "CLIEX01[replicate=1]"]')

        fixture = pathlib.Path(__file__).resolve().parent / "fixtures" / "cli_fault_rank1_rep0.txt"
        run("rank1-fault", ["--test", "CLIEX01", "--test-repeat", "3", "--test-seed", "123",
                            "--fail-fast", "--test-input", str(fixture),
                            "--test-output-dir", str(root / "rank1-fault"),
                            "--test-output", str(root / "rank1-fault")], expected=1)
        fault = summary("rank1-fault")
        check(fault["totals"]["failed"] == 1 and fault["not_run_invocation_count"] == 2, 'postcheck: fault["totals"]["failed"] == 1 and fault["not_run_invocation_count"] == 2')
        check(fault["not_run"] == ["CLIEX01#1..2"], 'postcheck: fault["not_run"] == ["CLIEX01#1..2"]')
        invocation = fault["results"][0]["invocations"][0]["result"]
        check(invocation["status"] == "FAIL", 'postcheck: invocation["status"] == "FAIL"')
        config = invocation["configuration"]
        check(any(item.startswith("mpi_rank_0=status=PASS;") for item in config), 'postcheck: any(item.startswith("mpi_rank_0=status=PASS;") for item in config)')
        check(any(item.startswith("mpi_rank_1=status=FAIL;") for item in config), 'postcheck: any(item.startswith("mpi_rank_1=status=FAIL;") for item in config)')
        check(not (root / "rank1-fault" / "CLIEX01" / "replicate-1").exists(), 'postcheck: not (root / "rank1-fault" / "CLIEX01" / "replicate-1").exists()')
        provenance = root / "rank1-fault" / "CLIEX01" / "replicate-0" / "rank-1" / "CLIEX01_fault_input.txt"
        check(provenance.read_bytes() == fixture.read_bytes(), 'postcheck: provenance.read_bytes() == fixture.read_bytes()')

        invalid = root / "invalid-fault.txt"
        invalid.write_text("srcsep-cli-fault-v1\nfault_rank=1\nfault_rank=1\nfault_replicate=0\n")
        run("invalid-fault", ["--test", "CLIEX01", "--test-input", str(invalid),
                              "--test-output-dir", str(root / "invalid-fault"),
                              "--test-output", str(root / "invalid-fault")], expected=2)
        check(summary("invalid-fault")["totals"]["errors"] == 1, 'postcheck: summary("invalid-fault")["totals"]["errors"] == 1')

        run("fixed-contract", ["--test", "CV01", "--test-seed", "1",
                               "--test-output", str(root / "fixed-contract")], expected=2)
        fixed = summary("fixed-contract")
        check(fixed["totals"]["errors"] == 1, 'postcheck: fixed["totals"]["errors"] == 1')
        check("unsupported" in fixed["results"][0]["message"], 'postcheck: "unsupported" in fixed["results"][0]["message"]')
        run("unknown-selector", ["--tests", "not-a-registered-selector"], expected=2)
        blocked = root / "report-target-is-directory"
        blocked.mkdir()
        run("output-preflight", ["--test", "CLIEX01", "--test-json", str(blocked),
                                 "--test-output", str(root / "output-preflight")], expected=2)
        check(not list((root / "output-preflight").rglob("CLIEX01_draw.json")), 'postcheck: not list((root / "output-preflight").rglob("CLIEX01_draw.json"))')
        rank_blocked = root / "rank1-artifact-error" / "CLIEX01" / "replicate-0" / "rank-1"
        rank_blocked.parent.mkdir(parents=True)
        rank_blocked.write_text("explicit rank-local directory obstruction\n")
        run("rank1-artifact-error", ["--test", "CLIEX01", "--test-seed", "123",
                                     "--test-output", str(root / "rank1-artifact-error")], expected=2)
        check(summary("rank1-artifact-error")["totals"]["errors"] == 1, 'postcheck: summary("rank1-artifact-error")["totals"]["errors"] == 1')
        print("PASS NATIVE-CLI: applied draws, repeats/order, JSON/JUnit, rank writers, artifact retention, non-root fail-fast and strict negative controls")
        return 0
    except Exception as error:
        if active is not None:
            active["status"] = "FAIL" if isinstance(error, AssertionError) else "ERROR"
            active["message"] = "%s: %s" % (type(error).__name__, error)
            active["postchecks"].append({"status": active["status"], "message": traceback.format_exc()})
        raise
    finally:
        finish()
        completed = {item["name"] for item in records}
        for index, name in enumerate(planned):
            if name not in completed:
                records.append({"id": "NCLI%02d" % (index + 1), "name": name, "status": "SKIP",
                                "message": "blocked by preceding native CLI infrastructure failure",
                                "command": [], "exit_code": None, "expected_exit_code": None,
                                "log": str(root / (name + ".log")), "postchecks": []})
        payload = {"schema": "srcsep-native-cli-campaign-v1",
                   "evidence_level": "infrastructure-probe;not-physics-validation",
                   "totals": {status: sum(item["status"] == status for item in records)
                              for status in ("PASS", "FAIL", "ERROR", "SKIP")},
                   "results": records, "runs": records}
        for leaf in ("summary.json", "native-cli-campaign-record.json"):
            temporary = root / (leaf + ".tmp")
            temporary.write_text(json.dumps(payload, indent=2) + "\n")
            temporary.replace(root / leaf)
        for item in records:
            print("[%s] %s %s: %s" % (item["id"], item["status"], item["name"], item["message"]))
        print("Native CLI campaign evidence: %s" % root)


if __name__ == "__main__":
    raise SystemExit(main())
