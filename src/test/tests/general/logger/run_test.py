#!/usr/bin/env python3
"""Build and run bounded/concurrent diagnostic logger tests."""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path


TEST_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(TEST_ROOT))

from lib.runner_support import (  # noqa: E402
    case, command, error_result, load_context, main_guard, run_command,
    write_result)


EXPECTED_CASES = [
    "GENERAL-LOGGER-FORMAT-LONG-001",
    "GENERAL-LOGGER-CONCURRENT-001",
    "GENERAL-LOGGER-ROTATE-FLUSH-001",
    "GENERAL-LOGGER-FATAL-001",
]


def prepare_header_mirror(repo: Path, work: Path) -> Path:
    source = repo / "src" / "general"
    target = work / "header-mirror" / "src" / "general"
    target.mkdir(parents=True, exist_ok=True)
    for item in source.iterdir():
        if item.is_file() and (item.suffix in (".h", ".dfn") or
                               item.name.endswith(".inc")):
            shutil.copyfile(str(item), str(target / item.name))
    (work / "header-mirror" / ".general.conf").write_text("", encoding="utf-8")
    return target


def main() -> int:
    context = load_context("Test AMPS synchronous diagnostic logging")
    repo = Path(context["repo_root"])
    work = Path(context["work_dir"])
    artifacts = Path(context["artifact_dir"])
    source = Path(__file__).with_name("test_general_logger.cpp")
    executable = work / "test_general_logger"
    timeout = int(context["test"].get("timeout_seconds", 180))
    requested = list(context["test"].get(
        "selected_cases", context["test"].get("expected_cases", EXPECTED_CASES)))
    try:
        include = prepare_header_mirror(repo, work)
        compiled = run_command(command(context, "mpi_cxx") + [
            "-std=c++17", "-O2", "-Wall", "-Wextra", "-Wpedantic",
            "-DOMPI_SKIP_MPICXX=1", "-pthread", "-I", str(include),
            "-I", str(TEST_ROOT / "native"),
            str(TEST_ROOT / "native" / "test_main.cpp"), str(source),
            "-o", str(executable),
        ], repo, timeout=timeout)
        if compiled.returncode != 0:
            return write_result(
                context, "FAIL", "diagnostic logger fixture did not compile",
                cases=[case(identifier, "FAIL", "compiler exit {}".format(
                    compiled.returncode)) for identifier in requested])

        results = []
        transcripts = []
        for identifier in requested:
            run = run_command([str(executable), "--case", identifier,
                               "--order", "normal"], work, timeout=timeout)
            transcripts.append("===== {} =====\n{}".format(
                identifier, run.stdout or ""))
            found = re.findall(
                r"^\[([A-Z][A-Z0-9-]+)\]\s+(PASS|FAIL)\s*(.*)$",
                run.stdout or "", flags=re.MULTILINE)
            exact = [entry for entry in found if entry[0] == identifier]
            foreign = [entry[0] for entry in found if entry[0] != identifier]
            if len(exact) == 1 and not foreign:
                _, status, message = exact[0]
                if run.returncode != (0 if status == "PASS" else 1):
                    status = "FAIL"
                    message = "case status/exit mismatch: exit {}".format(
                        run.returncode)
                results.append(case(identifier, status, message or "ok"))
            else:
                results.append(case(
                    identifier, "FAIL",
                    "expected one exact record; exact={} foreign={} exit={}".format(
                        len(exact), foreign, run.returncode)))

        artifacts.mkdir(parents=True, exist_ok=True)
        transcript = artifacts / "general-logger-output.txt"
        transcript.write_text("\n".join(transcripts), encoding="utf-8")
        evidence = [str(transcript)]
        for log in sorted(work.glob("amps-test-logger*.log*")):
            target = artifacts / log.name
            shutil.copyfile(str(log), str(target))
            evidence.append(str(target))
        status = "PASS" if all(item["status"] == "PASS" for item in results) else "FAIL"
        return write_result(
            context, status, "{} diagnostic logger cases executed".format(
                len(results)), cases=results, evidence=evidence)
    except Exception as error:
        return error_result(context, error)


if __name__ == "__main__":
    main_guard(main)
