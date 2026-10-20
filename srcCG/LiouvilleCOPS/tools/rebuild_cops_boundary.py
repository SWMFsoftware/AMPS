#!/usr/bin/env python3
"""Regenerate a static COPS fitting base with selected boundary and shape.

Uses only the Python standard library; runs the already compiled liouville_cops.
Relative output paths resolve under <AMPS>/cops_output (or COPS_OUTPUT_ROOT).
Other relative CLI and input-data paths resolve from the application root,
NOT the input file's directory. ``--shape-model`` and
``--shape-model-format`` override the corresponding source-input values in the
generated, archived input.  Existing output paths are always refused.
See README.md, section 6.6, for physics, safeguards, and limitations.
"""

import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import struct
import subprocess
import sys

from cops_paths import output_path, product_input
from cops_io import file_identity, load_manifest, RUN_SCHEMA


ROOT = Path(__file__).resolve().parents[1]
PATH_KEYS = {"ShapeModel", "ObservationFile", "KernelDirectory",
             "TemperatureRGTable", "ActivityFile"}


def rooted(value, base=None):
    path = Path(value)
    return Path(os.path.abspath(path if path.is_absolute() else (base or ROOT) / path))


def safe_input_path(path):
    value = str(path)
    if any(c.isspace() for c in value) or any(c in value for c in "#\"'$~"):
        raise ValueError("Input paths must not contain whitespace, quotes, #, $, or ~: " + value)
    return value


def parse_input(text):
    """Match C++ comment and repeated-key handling; retain all repeated values."""
    values = {}
    for n, line in enumerate(text.splitlines(), 1):
        active = line.split("#", 1)[0].strip()
        if not active:
            continue
        if "=" not in active:
            raise ValueError("Input line %d is not Key = Value" % n)
        key, value = (s.strip() for s in active.split("=", 1))
        if not key or not value:
            raise ValueError("Empty key/value on input line %d" % n)
        values.setdefault(key, []).append(value)
    return values


def last(values, key, default=None):
    return values.get(key, [default])[-1]


def integer(values, key, default, minimum):
    raw = last(values, key, str(default))
    try:
        number = float(raw)
        result = int(number)
    except (ValueError, OverflowError):
        raise ValueError(key + " must be an integer") from None
    if result != number or result < minimum:
        raise ValueError("%s must be an integer >= %d" % (key, minimum))
    return result


def prepare_input(text, boundary, out, threads, shape_model=None,
                  shape_model_format=None):
    """Return a normalized solver input with explicit rebuild overrides.

    ``shape_model`` is expected to be an already resolved path.  Keeping path
    resolution in ``main`` makes the function deterministic in unit tests and
    ensures that the generated input—not mutable shell state—records the exact
    geometry used to build every table and transport kernel.
    """
    values = parse_input(text)
    required = {"Mode": "observations", "FlightTimeRotation": "none",
                "TemperatureModel": "analytical"}
    for key, expected in required.items():
        if last(values, key, expected) != expected:
            raise ValueError("Fitting-base rebuild requires %s = %s (not changed automatically)"
                             % (key, expected))
    if float(last(values, "SurfaceTemperatureCorrectionFactor", "1")) != 1.0:
        raise ValueError("Fitting-base rebuild requires SurfaceTemperatureCorrectionFactor = 1")
    for key in ("ShapeModel", "ObservationFile", "KernelDirectory"):
        # A command-line mesh is a complete replacement for a missing source
        # ShapeModel.  Observation and kernel roots have no CLI replacement.
        if key not in values and not (key == "ShapeModel" and shape_model is not None):
            raise ValueError("Missing required input key: " + key)
    if not values.get("Kernel"):
        raise ValueError("No Kernel entries in input")
    first = integer(values, "FirstPoint", 0, 0)
    if "LastPoint" in values and integer(values, "LastPoint", 0, 0) < first:
        raise ValueError("LastPoint must be >= FirstPoint")
    integer(values, "PointStep", 1, 1)
    if threads is None and "Threads" in values:
        integer(values, "Threads", 1, 1)

    # The wrapper itself creates a brand-new directory before launching the
    # solver so it can archive inputs and a live log.  This explicit switch
    # permits that owned directory while the C++ solver still refuses any
    # pre-existing scientific products inside it.
    overrides = {"SurfaceBC": boundary, "WriteKernels": "on",
                 "AllowExistingOutputDirectory": "on",
                 "OutputDirectory": safe_input_path(out)}
    if threads is not None:
        overrides["Threads"] = str(threads)
    if shape_model is not None:
        overrides["ShapeModel"] = safe_input_path(shape_model)
    if shape_model_format is not None:
        if shape_model_format not in ("nastran", "dsk"):
            raise ValueError("ShapeModelFormat must be nastran or dsk")
        overrides["ShapeModelFormat"] = shape_model_format
    kernel_root = rooted(last(values, "KernelDirectory"))
    lines = ["# Generated by tools/rebuild_cops_boundary.py; keep with this base."]
    for line in text.splitlines():
        active, mark, comment = line.partition("#")
        if not active.strip():
            lines.append(line)
            continue
        key, value = (s.strip() for s in active.split("=", 1))
        if key in overrides:
            lines.append("# Replaced by rebuild tool: " + line)
            continue
        if key in PATH_KEYS or key == "Kernel":
            base = kernel_root if key == "Kernel" else ROOT
            path = product_input(value, app_root=ROOT) if key in ("TemperatureRGTable", "ActivityFile") else rooted(value, base)
            value = safe_input_path(path)
            line = key + " = " + value + ("  #" + comment if mark else "")
        lines.append(line)
    lines.append("\n# Boundary rebuild overrides (temperature and selection are preserved)")
    lines.extend(key + " = " + value for key, value in overrides.items())
    generated = "\n".join(lines) + "\n"
    return generated, parse_input(generated)


def preflight(values, executable):
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise ValueError("Executable missing/not executable; compile first with make: " + str(executable))
    files = [last(values, "ShapeModel"), last(values, "ObservationFile")]
    files.extend(values["Kernel"])
    if last(values, "TemperatureRGMode", "constant") == "hk_table":
        files.append(last(values, "TemperatureRGTable"))
    if last(values, "ActivityModel", "shap5") == "file":
        files.append(last(values, "ActivityFile"))
    for name in files:
        if not name or not Path(name).is_file():
            raise ValueError("Required data file is missing: " + str(name))
    # SPICE coverage, data syntax, unknown keys, and model validity are checked by C++.


def verify_outputs(out, values, boundary):
    """Check completion and the version-1 binary layout, not scientific accuracy."""
    log = (out / "rebuild.log").read_text()
    count = re.search(r"^Observations .*: (\d+) points,", log, re.M)
    if count is None:
        raise ValueError("Solver log has no observation count; cannot verify completeness")
    nobs = int(count.group(1))
    first = integer(values, "FirstPoint", 0, 0)
    final = min(integer(values, "LastPoint", nobs - 1, 0), nobs - 1)
    expected = list(range(first, final + 1, integer(values, "PointStep", 1, 1)))
    if not expected:
        raise ValueError("Requested observation range is empty")
    table = (out / "COPS_model.dat").read_text()
    headers = "\n".join(line for line in table.splitlines() if line.startswith("#"))
    if not re.search(r"\bSurfaceBC=" + boundary + r"\b", headers):
        raise ValueError("Output table has the wrong/missing SurfaceBC")
    if "FlightTimeRotation=none" not in headers:
        raise ValueError("Output table is not a static run")
    actual = [int(line.split()[0]) for line in table.splitlines()
              if line.strip() and not line.startswith("#")]
    if actual != expected:
        raise ValueError("Output table observation indices are incomplete or unexpected")
    elements = out / "elements.bin"
    with elements.open("rb") as stream:
        header = stream.read(12)
    if len(header) != 12:
        raise ValueError("Truncated elements header")
    version, nfaces, nbasis = struct.unpack("<iii", header)
    if (version != 1 or nfaces <= 0 or nbasis != 25
            or elements.stat().st_size != 12 + nfaces * 272):
        raise ValueError("Unsupported/truncated elements.bin")
    for ip in expected:
        path = out / "kernels" / ("point=%04d.bin" % ip)
        with path.open("rb") as stream:
            header = stream.read(92)
        if len(header) != 92:
            raise ValueError("Truncated kernel header: " + str(path))
        version, saved_ip = struct.unpack_from("<ii", header)
        nk, = struct.unpack_from("<i", header, 88)
        if (version != 1 or saved_ip != ip or not 0 <= nk <= nfaces
                or path.stat().st_size != 92 + 40 * nk):
            raise ValueError("Unsupported/truncated/mismatched kernel: " + str(path))
    return len(expected)


def timestamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", default="fit_base.in", help="source input (default: fit_base.in)")
    parser.add_argument("--surface-bc", required=True, choices=("flat", "rough"), help="new angular emission law")
    parser.add_argument("--out", required=True, help="NEW directory under cops_output (or absolute path outside LiouvilleCOPS); existing paths refused")
    parser.add_argument("--threads", type=int, help="positive OpenMP thread count; otherwise preserve input/environment")
    parser.add_argument("--shape-model",
                        help="surface mesh overriding ShapeModel in the source input; relative paths use the application root")
    parser.add_argument("--shape-model-format", choices=("nastran", "dsk"),
                        help="mesh reader overriding ShapeModelFormat; preserve the source value/default when omitted")
    parser.add_argument("--executable", default="liouville_cops", help="compiled solver (default: application liouville_cops)")
    parser.add_argument("--dry-run", action="store_true", help="print input and command; create nothing, do not check external files")
    args = parser.parse_args(argv)
    out = None
    manifest = None
    try:
        if args.threads is not None and args.threads <= 0:
            raise ValueError("--threads must be positive")
        source, executable = (rooted(v) for v in (args.input, args.executable))
        if not source.exists():
            source = Path(product_input(args.input, app_root=ROOT))
        shape_model = rooted(args.shape_model) if args.shape_model else None
        out = output_path(args.out).absolute()
        safe_input_path(source)
        safe_input_path(out)
        if os.path.lexists(out):
            raise ValueError("Output already exists; choose a fresh directory: " + str(out))
        original_bytes = source.read_bytes()
        text = original_bytes.decode("utf-8")
        generated, values = prepare_input(text, args.surface_bc, out, args.threads,
                                          shape_model, args.shape_model_format)
        input_path = out / "boundary.in"
        command = [str(executable), str(input_path)]
        print("Source: %s\nNew boundary: %s\nShape model: %s (%s)\nOutput: %s\nCommand: %s" %
              (source, args.surface_bc, last(values, "ShapeModel"),
               last(values, "ShapeModelFormat", "nastran"), out,
               shlex.join(command)), flush=True)
        if args.dry_run:
            print("Dry run: external file availability and SPICE coverage are NOT checked.\n")
            print(generated)
            return 0
        preflight(values, executable)
        out.mkdir(parents=True, exist_ok=False)
        (out / "source.in").write_bytes(original_bytes)
        input_path.write_text(generated)
        manifest = {"status": "running", "started_utc": timestamp(),
                    "source_input": str(source),
                    "source_sha256": hashlib.sha256(original_bytes).hexdigest(),
                    "generated_input": str(input_path),
                    "surface_bc": args.surface_bc,
                    "shape_model": last(values, "ShapeModel"),
                    "shape_model_format": last(values, "ShapeModelFormat", "nastran"),
                    "shape_model_overridden": args.shape_model is not None,
                    "shape_model_format_overridden": args.shape_model_format is not None,
                    "command": command, "cwd": str(ROOT)}
        manifest_path = out / "rebuild.json"
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        with (out / "rebuild.log").open("w") as log:
            process = subprocess.Popen(command, cwd=ROOT, stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT, text=True, bufsize=1)
            try:
                for line in process.stdout:
                    log.write(line)
                    log.flush()
                    print(line, end="", flush=True)
                code = process.wait()
            except BaseException:
                process.terminate()
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
                raise
            finally:
                process.stdout.close()
        manifest["solver_returncode"] = code
        if code:
            raise ValueError("Forward solver failed with exit status %d; see rebuild.log" % code)
        manifest["verified_points"] = verify_outputs(out, values, args.surface_bc)
        run_manifest_path = out / "run_manifest.json"
        if run_manifest_path.is_file():
            run_manifest = load_manifest(run_manifest_path, RUN_SCHEMA)
            if run_manifest.get("states", {}).get("process") != "complete":
                raise ValueError("C++ run manifest does not mark the process complete")
            if run_manifest.get("states", {}).get("cache") != "verified":
                raise ValueError("C++ run manifest does not mark the cache verified")
            manifest["run_manifest"] = file_identity(run_manifest_path)
        else:
            # Fixture/legacy executables can still exercise this wrapper, but
            # the absence is explicit and cannot be confused with W01 output.
            manifest["run_manifest"] = None
        manifest.update(status="complete", finished_utc=timestamp())
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        print("Complete: regenerated and checked %d static point kernels." % manifest["verified_points"])
        print("Next: python3 tools/fit_cops.py %s --data-only --out %s" %
              (shlex.quote(str(out)), shlex.quote(str(out) + "_selection")))
        print("Choose NG exclusion/saturation switches explicitly before fitting; see quick_guide.md.")
        return 0
    except (OSError, ValueError, KeyboardInterrupt) as exc:
        message = str(exc) or "Interrupted"
        if manifest is not None:
            manifest.update(status="failed", finished_utc=timestamp(), error=message)
            (out / "rebuild.json").write_text(json.dumps(manifest, indent=2) + "\n")
        print("Error: " + message, file=sys.stderr)
        if manifest is not None:
            print("Partial output retained for diagnosis; do not fit it. Retry with a new --out.", file=sys.stderr)
        return 130 if isinstance(exc, KeyboardInterrupt) else 1


if __name__ == "__main__":
    sys.exit(main())
