"""Run per-ID infrastructure checks without changing scientific gate status.

These checks exercise temporary metadata/arithmetic fixtures only. Success
allows the ordinary shared runner to decide PASS/FAIL/SKIPPED; it can never
qualify an observational comparison or replace a predecessor gate. Failure is
test infrastructure ERROR, including when the science run would be skipped.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys


def run_protocol_checks(test_id: str, scripts: tuple[Path, ...]) -> bool:
    """Preserve each subordinate command/log and return whether all passed.

    The individual runner's existing --output-dir remains authoritative. The
    aggregate runner already gives every campaign a fresh directory, so its
    per-ID protocol_checks/ files cannot replace earlier science evidence.
    With no output option (for example --help), diagnostics remain on stdout
    and stderr; the normal argument parser still owns CLI validity.
    """
    output_root = None
    for index, argument in enumerate(sys.argv[1:], start=1):
        if argument.startswith("--output-dir="):
            output_root = Path(argument.split("=", 1)[1])
        elif argument == "--output-dir" and index + 1 < len(sys.argv):
            output_root = Path(sys.argv[index + 1])
    output = None if output_root is None else output_root / test_id
    directory = None if output is None else output / "protocol_checks"
    if directory is not None:
        directory.mkdir(parents=True, exist_ok=True)
    checks = []
    for index, script in enumerate(scripts, start=1):
        script = script.resolve()
        command = [sys.executable, "-B", str(script)]
        try:
            completed = subprocess.run(command, text=True, capture_output=True,
                                       check=False)
            code, stdout, stderr = completed.returncode, completed.stdout, completed.stderr
        except OSError as error:
            code, stdout, stderr = None, "", str(error) + "\n"
        # Preserve diagnostics in the aggregate runner log as well as the
        # individual files, so a failed check is discoverable from either CLI.
        print(f"{test_id} infrastructure-only contract: {script}", flush=True)
        print(stdout, end="", flush=True)
        print(stderr, end="", file=sys.stderr, flush=True)
        item = {"script": str(script), "command": command,
                "scope": "infrastructure-only-not-scientific-validation",
                "status": "PASS" if code == 0 else "ERROR", "returncode": code}
        if directory is not None:
            label = f"{index:02d}-{script.stem}"
            stdout_path, stderr_path = directory / (label + ".stdout.log"), directory / (label + ".stderr.log")
            stdout_path.write_text(stdout, encoding="utf-8")
            stderr_path.write_text(stderr, encoding="utf-8")
            item.update(stdout=str(stdout_path), stderr=str(stderr_path))
        if script.is_file():
            item["script_sha256"] = hashlib.sha256(script.read_bytes()).hexdigest()
        checks.append(item)
    passed = all(item["status"] == "PASS" for item in checks)
    summary = {"schema": "srcMoon-per-ID-infrastructure-checks-v1", "id": test_id,
               "scope": "infrastructure-only-not-scientific-validation",
               "status": "PASS" if passed else "ERROR", "checks": checks}
    if directory is not None:
        (directory / "result.json").write_text(json.dumps(summary, indent=2) + "\n")
        if not passed:
            # main() will not execute after a failed supporting check. Supply
            # its required per-ID ERROR artifact instead of allowing a missing
            # result to be confused with a skipped scientific prerequisite.
            (output / "result.json").write_text(json.dumps({
                "id": test_id, "status": "ERROR", "execution": "not-run",
                "metrics": {}, "reason": "infrastructure contract failed; see protocol_checks logs",
                "infrastructure_checks": summary,
            }, indent=2) + "\n")
    return passed
