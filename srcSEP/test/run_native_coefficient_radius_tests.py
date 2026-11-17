#!/usr/bin/env python3
"""Run the native radius gate through opted-in immutable schema-6 --test-input.

The complete shipped schema6 input is the reviewed base. Only explicit
geometry/numerical initialization settings and the declared wind are changed.
There is no in-process Initialization override, second scientific parser or
source/particle timestep. Every output belongs to fresh isolated scratch.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import time

sys.dont_write_bytecode = True
from run_fresh_process_checkpoint_tests import binary_hash, create_native_work_directory, run_role

SOURCE = Path(__file__).resolve().parents[1]


def replace_key(text: str, key: str, value: str) -> str:
    """Edit one known template key, requiring exactly one reviewed occurrence."""
    expression = re.compile(r"^([ \t]*" + re.escape(key) + r"[ \t]*=[ \t]*)[^!\r\n]*", re.MULTILINE)
    result, count = expression.subn(lambda match: match.group(1) + value, text)
    if count != 1:
        raise ValueError(f"reviewed input key {key!r} has {count} occurrences")
    return result


def fixture_input(template: str) -> tuple[str, dict[str, str]]:
    """Preserve the actual complete physics record; declare a bounded geometry.

    The unchanged2e11m arc fits strictly inside the declared3e11m outer sphere:
    every point's radial distance is <= initial radius+arc length. Existing
    observer/source/thermal fields remain explicit. A coarse support mesh and
    17vertices suffice for this read-only accessor, not transport convergence.
    """
    changes = {
        "origin_x_m": "2.0e10", "origin_y_m": "-1.0e10", "origin_z_m": "3.0e9",
        "initial_x_m": "3.3914e10", "initial_y_m": "-1.0e10", "initial_z_m": "3.0e9",
        "ambient.wind_speed": "725 km/s", "point_count": "17", "maximum_level": "0",
        "outer_radius_m": "3.0e11",
    }
    text = template
    for key, value in changes.items():
        text = replace_key(text, key, value)
    for kind in ("solar", "tube"):
        expression = re.compile(r"(^[ \t]*#begin refinement:[ \t]*" + kind +
                                r"[^\r\n]*\n)(.*?)(^[ \t]*#end refinement\b)", re.MULTILINE | re.DOTALL)
        matches = list(expression.finditer(text))
        if len(matches) != 1:
            raise ValueError(f"reviewed refinement:{kind} block is not unique")
        match = matches[0]
        body = replace_key(match.group(2), "enabled", "false")
        text = text[:match.start(2)] + body + text[match.end(2):]
        changes[f"mesh.refinement:{kind}.enabled"] = "false"
    return text, changes


def library_fixture_input(template: str) -> tuple[str, dict[str, str]]:
    """Declare the reviewed library's explicit speed/radius diagnostic law.

    The complete shipped library example remains the scientific base. Only
    its one identified parallel-diffusion section is replaced by explicit SI
    model parameters, all accepted and installed by the actual library parser.
    This is a deterministic fixture editor, not a second scientific parser.
    No nucleon, field decomposition or externally supplied factor is assumed.
    """
    text, changes = fixture_input(template)
    expression = re.compile(r"(^[ \t]*#begin physics:[ \t]*parallel_diffusion[^\r\n]*\n)"
                            r"(.*?)(^[ \t]*#end physics\b)", re.MULTILINE | re.DOTALL)
    values = ("  ! Explicit read-only native radius diagnostic, not a fitted transport law.\n"
              "  model = power_law_lambda\n"
              "  lambda0_m = 3.0e9\n"
              "  independent_variable = speed\n"
              "  speed0_m_per_s = 2.0e6\n"
              "  independent_exponent = 0.75\n"
              "  radius0_m = 1.495978707e11\n"
              "  radial_exponent = 1.5\n")
    text, count = expression.subn(lambda match: match.group(1)+values+match.group(3), text)
    if count != 1:
        raise ValueError(f"reviewed parallel_diffusion block has {count} occurrences")
    changes["physics.parallel_diffusion"] = values
    return text, changes


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amps", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mpi-ranks", type=int, default=2)
    parser.add_argument("--mpi-exec", default="mpiexec")
    parser.add_argument("--timeout", type=float, default=1200)
    parser.add_argument("--prepare-only", action="store_true", help="write reviewable input/manifest without executing; never reports PASS")
    parser.add_argument("--case", choices=("RADIUSN01", "RADIUSN02"), default="RADIUSN01",
                        help="QLT1 scalar radius or independently installed speed/radius library diagnostic")
    args = parser.parse_args()
    amps, output = args.amps.resolve(), args.output.resolve()
    if not amps.is_file() or args.mpi_ranks < 1 or args.timeout <= 0:
        parser.error("a linked executable, positive MPI count and timeout are required")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        parser.error("choose a new/empty output directory; previous evidence is preserved")
    output.mkdir(parents=True, exist_ok=True)
    template_name = "sep_parker_mesh.in" if args.case == "RADIUSN01" else "sep_parker_mesh_parallel_diffusion.in"
    template_path = SOURCE / "examples" / template_name
    template = template_path.read_text(encoding="utf-8")
    text, changes = fixture_input(template) if args.case == "RADIUSN01" else library_fixture_input(template)
    owned_input = output / "translated-radius.in"
    owned_input.write_text(text, encoding="utf-8")
    executable_hash = binary_hash(amps)
    summary = {
        "schema": "srcsep-native-coefficient-radius-v1", "id": args.case,
        "executable": str(amps), "executable_sha256": executable_hash,
        "template": str(template_path), "template_sha256": hashlib.sha256(template.encode()).hexdigest(),
        "owned_input": str(owned_input), "owned_input_sha256": hashlib.sha256(text.encode()).hexdigest(),
        "explicit_changes": changes, "mpi_ranks": args.mpi_ranks,
        "scope": "actual-immutable-initialization-and-read-only-coefficient;no-hoststep;no-production-qualification",
    }
    if args.prepare_only:
        summary["status"] = "PREPARED"
        print("PREPARED", owned_input)
    else:
        work = create_native_work_directory(output, "radius")
        report = output / "srcsep-tests.json"
        # Only the selected radius descriptor's explicit capability authorizes the
        # standalone driver to load this complete reviewed record. The same
        # production LoadFile/Install and file-then-CLI layering are used;
        # production --input remains incompatible with test mode. The template
        # declares relative output leaves, all confined to this owned native
        # cwd. No initialization-only output flag is smuggled into test mode.
        command = [args.mpi_exec, "-n", str(args.mpi_ranks), str(amps),
                   "--test-input", str(owned_input), "--shock-model", "analytical",
                   "--particle-mover", "parker",
                   "--test", args.case, "--test-output-dir", str(output / "artifacts"),
                   "--test-json", str(report), "--test-junit", str(output / "srcsep-tests.xml")]
        if args.case == "RADIUSN02":
            command.extend(("--spatial-diffusion-provider", "parallel-diffusion-library",
                            "--coefficient-source", "prescribed"))
        started = time.monotonic()
        try:
            code = run_role(command, work, output / "native.log", args.timeout)
            document = json.loads(report.read_text(encoding="utf-8")) if report.is_file() else {}
            rows = [row for row in document.get("results", []) if row.get("id") == args.case]
            passed = code == 0 and len(rows) == 1 and rows[0].get("status") == "PASS"
            passed = passed and binary_hash(amps) == executable_hash
            passed = passed and hashlib.sha256(owned_input.read_bytes()).hexdigest() == summary["owned_input_sha256"]
            message = "" if passed else f"native exit={code}; selected result={[row.get('status') for row in rows]}; binary/input identity pinned"
        except (OSError, ValueError) as error:
            code, passed, message = 2, False, str(error)
        summary.update(status="PASS" if passed else "FAIL", exit_code=code, message=message,
                       duration_s=time.monotonic()-started, report=str(report), log=str(output/"native.log"),
                       working_directory=str(work), command=command)
        print(args.case, summary["status"], message, flush=True)
    temporary = output / "summary.json.tmp"
    temporary.write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    os.replace(temporary, output/"summary.json")
    return 0 if summary["status"] in ("PASS", "PREPARED") else 1


if __name__ == "__main__":
    raise SystemExit(main())
