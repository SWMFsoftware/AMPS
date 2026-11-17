#!/usr/bin/env python3
"""Native production-CLI software contracts and explicit remaining gates."""

import argparse
import hashlib
import json
import pathlib
import re
import subprocess
import tempfile
import traceback


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amps", required=True, type=pathlib.Path)
    parser.add_argument("--mpiexec", default="mpiexec")
    parser.add_argument("--output", type=pathlib.Path)
    parser.add_argument("--source-root", type=pathlib.Path,
                        default=pathlib.Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    executable = args.amps.resolve(strict=True)
    root = (args.output or pathlib.Path(tempfile.mkdtemp(prefix="srcsep-production-cli."))).resolve()
    root.mkdir(parents=True, exist_ok=True)
    results = []
    active = None

    def require(condition, message):
        if active is not None:
            active["postchecks"].append({"status": "PASS" if condition else "FAIL", "message": message})
        if not condition:
            raise AssertionError(message)

    def record(identifier, message):
        nonlocal active
        if active is not None and active["status"] == "RUNNING":
            active["status"] = "PASS"
        active = {"id": identifier, "status": "RUNNING", "message": message,
                  "command": [], "exit_code": None, "log": "", "postchecks": [], "artifacts": []}
        results.append(active)
        return active

    def process(identifier, command, cwd, expected, message, check_exit=True):
        item = record(identifier, message)
        item["command"] = command
        item["expected_exit_code"] = expected
        item["log"] = str(root / (identifier + ".log"))
        try:
            child = subprocess.run(command, cwd=str(cwd), stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True, timeout=90)
        except (OSError, subprocess.TimeoutExpired) as error:
            pathlib.Path(item["log"]).write_text(str(error))
            raise
        pathlib.Path(item["log"]).write_text(child.stdout)
        item["exit_code"] = child.returncode
        if check_exit:
            require(child.returncode == expected, "actual exit=%d expected=%d" % (child.returncode, expected))
        return child

    try:
        native = root / "native"
        native_process=process("PCLIN01", [args.mpiexec, "-n", "2", str(executable), "--test-group", "production-cli",
                            "--test-repeat", "2", "--test-output", str(native)], root, 0,
                "actual linked production globals and scoped repetition on two MPI ranks", check_exit=False)
        payload = json.loads((native / "results.json").read_text())
        for case in payload["results"]:
            case["log"]=str(root / "PCLIN01.log")
            case["artifacts"] = [str(native / "results.json")] + case.get("artifacts", [])
            results.append(case)
        require(native_process.returncode == 0, "native production CLI group must exit successfully")
        require(payload["totals"] == {"passed": 12, "failed": 0, "skipped": 0, "errors": 0},
                "all twelve linked software contracts must pass")
        require([item["id"] for item in payload["results"]] == ["PCLI%02d" % i for i in range(1, 13)],
                "canonical production CLI IDs must be complete and unique")
        active["artifacts"].extend([str(native / "results.json"), str(native / "results.xml")])
        for item in payload["results"]:
            require(len(item["invocations"]) == 2, item["id"] + " must retain both repeated invocations")
            for invocation in item["invocations"]:
                child = invocation["result"]
                require(child["status"] == "PASS", item["id"] + " repeat failed")
                for rank in (0, 1):
                    require(any(value.startswith("mpi_rank_%d=status=PASS;" % rank) for value in child["configuration"]),
                            item["id"] + " must retain successful evidence for rank%d" % rank)

        hostile = root / "hostile-default-input"
        hostile.mkdir(exist_ok=True)
        (hostile / "amps.in").write_text("deliberately invalid input: help/errors must not initialize it\n")
        for identifier, flag in (("PCLIH01", "--help"), ("PCLIH02", "-h")):
            child = process(identifier, [str(executable), flag], hostile, 0,
                            "fresh executable help exits without reading default input or creating model output")
            require("--particle-mover" in child.stdout and "--test-seed" in child.stdout,
                    "fresh help must advertise production and campaign options")
            require(sorted(item.name for item in hostile.iterdir()) == ["amps.in"],
                    "help must not create a model/evidence file")
        child = process("PCLIE01", [str(executable), "--unregistered-production-cli-option"], hostile, 1,
                        "unknown production argument preserves nonzero compatibility exit before initialization")
        require("unknown SEP command-line option" in child.stdout, "unknown option must be diagnosed by the production parser")
        require(sorted(item.name for item in hostile.iterdir()) == ["amps.in"],
                "unknown option must not initialize default input")

        # This is a limited source callsite audit, explicitly separate from an
        # actual coupled SWMF execution. Remove comments/literals so prose is
        # not counted as an invocation. It does not substitute for a C++ AST.
        def code(path):
            text = path.read_text()
            return re.sub(r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'', "", text, flags=re.S)

        source = args.source_root.resolve(strict=True)
        main_source, coupled_source = code(source / "main.cpp"), code(source / "main_lib.cpp")
        record("PCLIS01", "source callsites preserve input-before-CLI application order")
        legacy = main_source.index("SEP::Parser::ReadFile(PIC::PostCompileInputFileName)")
        apply = main_source.index("SEP::Util::CLI::ApplyTurbulenceOptions(cli_options)")
        layer = main_source.index("SEP::Util::CLI::ApplyInputFileTurbulence(")
        reparse = main_source.index("SEP::Util::CLI::ParseCommandLine(argc, argv, layered")
        require(legacy < apply and layer < reparse < apply, "input parsing/layering must precede CLI application")
        active["artifacts"] = [str(source / "main.cpp")]
        record("PCLIS02", "limited coupled-library standalone-CLI callsite audit; no SWMF-run claim")
        for symbol in ("ParseCommandLine", "TestManager", "RunSelectedComponentTests", "ConfigureComponentTestCampaign"):
            require(re.search(r"\b" + symbol + r"\s*\(", coupled_source) is None,
                    "coupled main_lib must not call " + symbol)
        require(re.search(r"\bTestManager\s*\(", main_source) is not None,
                "standalone driver must retain explicit legacy TestManager ownership")
        active["artifacts"] = [str(source / "main_lib.cpp"), str(source / "main.cpp")]
        active["scope"] = "source-callsite-audit"
        audit = root / "ownership-source-audit.json"
        audit.write_text(json.dumps({"scope": "limited-source-callsite-audit;not-full-SWMF-execution",
                                     "sources": [{"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                                                 for path in (source / "main.cpp", source / "main_lib.cpp")],
                                     "forbidden_coupled_calls": ["ParseCommandLine", "TestManager", "RunSelectedComponentTests",
                                                                "ConfigureComponentTestCampaign"]}, indent=2) + "\n")
        active["artifacts"].append(str(audit))
        active["status"] = "PASS"
        active = None
        results.extend([
            {"id": "PCLI11-PRODUCTION", "status": "SKIP", "message":
             "192 native applied-state cases do not establish post-allocation production evolution; polarized closure remains unqualified",
             "log": str(root / "PCLIN01.log"), "postchecks": [], "artifacts": []},
            {"id": "PCLI12-SWMF", "status": "SKIP", "message":
             "source ownership and private parser repetition do not establish a complete coupled SWMF execution",
             "log": str(root / "PCLIN01.log"), "postchecks": [], "artifacts": []}])
        return 0
    except Exception as error:
        if active is not None:
            active["status"] = "FAIL" if isinstance(error, AssertionError) else "ERROR"
            active["message"] = "%s: %s" % (type(error).__name__, error)
            active["postchecks"].append({"status": active["status"], "message": traceback.format_exc()})
        raise
    finally:
        if active is not None and active["status"] == "RUNNING":
            active["status"] = "PASS"
        report = {"schema": "srcsep-production-cli-compatibility-v1", "qualification": "partial-software-contract-evidence",
                  "results": results,
                  "totals": {status: sum(item["status"] == status for item in results)
                             for status in ("PASS", "FAIL", "ERROR", "SKIP")}}
        temporary = root / "summary.json.tmp"
        temporary.write_text(json.dumps(report, indent=2) + "\n")
        temporary.replace(root / "summary.json")
        for item in results:
            print("[%s] %s %s" % (item["id"], item["status"], item["message"]))
        print("Production CLI software evidence: %s" % root)


if __name__ == "__main__":
    raise SystemExit(main())
