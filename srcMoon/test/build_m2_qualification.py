#!/usr/bin/env python3
"""Build the dedicated Na+Ar srcMoon M2 qualification executable.

The canonical ``input/moon.input`` is the frozen Na-only M0/M1 application and
must not be edited merely to run an M2 species test.  This script derives a
temporary Config.pl application with reviewed compile-time changes: add
neutral Ar to ``SpeciesList``, register the dedicated ``RadiogenicArgon``
source, disable photolytic chemistry for the I19 accounting contract, and
enable individual particle weights for fractional finite-release packets.
It verifies that each applicable edit occurs once,
performs a clean SPICE build, copies the resulting executable into the
campaign directory, and removes only the three temporary input files it made.
It also removes the uniquely prefixed root copies and `.Assembled*`
intermediates that Config.pl derives from those inputs.

Runtime activation remains separate: I18--I20 derive a preserved ``.in``
receipt with ``argon_source_enabled=true``.  This two-key arrangement makes it
impossible for an accidental Na-only build or an accidental runtime switch to
silently masquerade as a two-species qualification run.

The optional purpose ``sampling`` applies only the first two species/source
changes, retaining canonical chemistry and individual-weight controls,
for I11's deterministic storage regression. Its executable is named separately
and the build receipt explicitly omits the I19 chemistry contract. The default
M2 purpose still requires the original strict I19 generated-configuration gate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time


REPO = Path(__file__).resolve().parents[2]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def replace_once(text: str, pattern: str, replacement: str,
                 description: str) -> str:
    updated, count = re.subn(pattern, replacement, text, count=0,
                             flags=re.MULTILINE)
    if count != 1:
        raise ValueError(
            f"expected exactly one {description}; found {count}")
    return updated


def derive_compile_input(source: Path, application_name: str,
                         purpose: str = "m2") -> str:
    """Return the reviewed two-species Config.pl input without writing it."""
    text = source.read_text(encoding="utf-8")
    text = replace_once(
        text, r"^(\s*SpeciesList\s*=\s*)Na(\s*!.*)$",
        r"\1Na,Ar\2", "active canonical SpeciesList=Na declaration")

    # Anchor the source edit to its unique process code.  A broad replacement
    # of ``Source:UserDefined=off`` could enable the historical MySource block,
    # which the M0 decision explicitly keeps disabled.
    block_pattern = (
        r"(^\s*Source:UserDefined\s*=\s*)off"
        r"(?=(?:[^\n]*\n){1,3}\s*SourceProcessCode\s*=\s*"
        r"RadiogenicArgon\b)")
    text = replace_once(
        text, block_pattern, r"\1on", "RadiogenicArgon source block")
    if purpose == "m2":
        # I19 explicitly excludes chemistry. The generated macro check below
        # protects this physical contract; make the corresponding source
        # transformation rather than failing after an expensive wrong build.
        text = replace_once(text, r"^(\s*PhotolyticReactions\s*=\s*)on\b",
                            r"\1off", "canonical photolytic reaction switch")
        text = replace_once(text, r"^(\s*ParticleWeightCorrectionMode\s*=\s*)off\b",
                            r"\1on", "finite-release individual particle weights")

    # Config.pl copies files matching the selected application basename into
    # the repository root.  Point the derived input at its equally temporary
    # parser/kernel companions so it never depends on a stale generated
    # ``moon.input.pl`` from an earlier canonical configuration.
    text = replace_once(
        text, r"^#block\s+\.\./moon\.input\.pl\s*.*$",
        f"#block ../{application_name}.input.pl ! generated M2 qualification parser",
        "moon.input.pl block reference")
    text = replace_once(
        text, r"^\s*#include\s+moon\.input\.spicekernels\s*$",
        f"  #include {application_name}.input.spicekernels",
        "moon.input.spicekernels include")
    return text


def run_logged(command: list[str], log: Path) -> tuple[int, float]:
    started = time.monotonic()
    with log.open("a", encoding="utf-8") as stream:
        stream.write("command: " + " ".join(command) + "\n")
        stream.flush()
        completed = subprocess.run(
            command, cwd=REPO, stdout=stream, stderr=subprocess.STDOUT,
            check=False)
    return completed.returncode, time.monotonic() - started


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--spice-path", type=Path, required=True)
    parser.add_argument("--spice-kernels", type=Path, required=True)
    parser.add_argument("--purpose", choices=("m2", "sampling"), default="m2",
                        help="sampling builds a two-species storage fixture without claiming the I19 chemistry contract")
    args = parser.parse_args()

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    executable = output / ("amps-moon-m2-na-ar" if args.purpose == "m2" else
                           "amps-moon-sampling-na-ar")
    log = output / "build-m2-na-ar.log"
    result_path = output / "build-m2-na-ar.json"
    # Sandbox/container launches can reuse the same small PID, so a PID alone
    # is not a safe application identity.  Bind the generated name to the
    # absolute campaign path; the runner already requires that path to be new
    # or empty, making the first 12 SHA-256 digits deterministic and collision
    # resistant without introducing an unrecorded random value.
    application_token = hashlib.sha256(
        str(output).encode("utf-8")).hexdigest()[:12]
    application_name = f"moon_m2_qualification_{application_token}"
    source_input = REPO / "input/moon.input"
    staged = {
        REPO / f"input/{application_name}.input": None,
        REPO / f"input/{application_name}.input.pl":
            REPO / "input/moon.input.pl",
        REPO / f"input/{application_name}.input.spicekernels":
            REPO / "input/moon.input.spicekernels",
    }
    # Config.pl copies the selected application's three inputs from input/
    # into the repository root.  They are generated build intermediates, not
    # source files.  Track and remove these exact uniquely named copies as
    # well as the staged inputs so M2 qualification leaves no Python/build
    # debris for a later commit.
    generated_copies = tuple(
        REPO / f"{application_name}{suffix}"
        for suffix in (".input", ".input.pl", ".input.spicekernels")
    )
    assembled_pattern = f"{application_name}.input.Assembled*"
    result: dict[str, object] = {
        "schema": "srcMoon-M2-Na-Ar-build-v1",
        "purpose": args.purpose,
        "status": "ERROR",
        "source_input": {"path": str(source_input),
                         "sha256": sha256(source_input)},
        "compile_time_changes": [
            "SpeciesList=Na -> SpeciesList=Na,Ar",
            "RadiogenicArgon Source:UserDefined=off -> on",
        ],
        "commands": [],
        "log": str(log),
    }
    # Cleanup ownership starts only after every temporary name is verified
    # absent. If a previous/user-created file causes the preflight to refuse,
    # finally must preserve it rather than remove the very collision reported.
    cleanup_owned = False
    if args.purpose == "m2":
        result["compile_time_changes"].extend([
            "PhotolyticReactions=on -> off (I19 chemistry contract)",
            "ParticleWeightCorrectionMode=off -> on (exact partial release packets)"])
    try:
        for path in (*staged, *generated_copies,
                     *tuple(REPO.glob(assembled_pattern))):
            if path.exists():
                raise FileExistsError(
                    f"refusing to overwrite temporary application file {path}")
        cleanup_owned = True
        derived = derive_compile_input(source_input, application_name,args.purpose)
        # The staged files are build inputs, not validation evidence.  Their
        # exact hashes are retained before cleanup in the machine-readable
        # build receipt.
        (REPO / f"input/{application_name}.input").write_text(
            derived, encoding="utf-8")
        for destination, source in staged.items():
            if source is not None:
                shutil.copy2(source, destination)
        result["staged_input_hashes"] = {
            path.name: sha256(path) for path in staged
        }

        build = REPO / "build"
        if build.exists():
            shutil.rmtree(build)
        commands = [
            # Select production mode explicitly: Config.pl otherwise retains
            # a previous nightly build's sorting/debug/iteration policies.
            [str(REPO / "Config.pl"), f"-application={application_name}", "-amps-test=off",
             f"-spice-path={args.spice_path.resolve()}",
             f"-spice-kernels={args.spice_kernels.resolve()}"],
            # The current core native owners use sixteen build jobs. Keep
            # debug symbols/optimization explicit while bounding compilers.
            ["make", "-j16", "OPT3=-O3 -g"],
        ]
        log.write_text("", encoding="utf-8")
        for command in commands:
            returncode, duration = run_logged(command, log)
            result["commands"].append({
                "command": command, "returncode": returncode,
                "duration_seconds": duration})
            if returncode != 0:
                raise RuntimeError(
                    f"M2 qualification build command failed ({returncode}): "
                    + " ".join(command))

        # I19's approved M2 accounting case excludes chemistry so its
        # interspecies conversion terms have a physical expectation of
        # exactly zero.  Verify the *generated* contract instead of assuming
        # that a human-readable Config.pl input spelling maps to a particular
        # macro.  A future configuration change that enables photochemistry
        # must fail this dedicated build rather than silently alter I19's
        # closure equation.
        generated_definitions = REPO / "build/pic/picGlobal.dfn"
        definitions_text = generated_definitions.read_text(encoding="utf-8")
        chemistry_pattern = (
            r"^#define\s+_PIC_PHOTOLYTIC_REACTIONS_MODE_\s+"
            r"_PIC_PHOTOLYTIC_REACTIONS_MODE_OFF_\s*$"
        )
        if (args.purpose == "m2" and
                re.search(chemistry_pattern, definitions_text,
                          flags=re.MULTILINE) is None):
            raise RuntimeError(
                "dedicated M2 build violates the approved I19 contract: "
                "photolytic chemistry is not compiled OFF")
        if args.purpose == "m2" and re.search(
                r"^#define\s+_INDIVIDUAL_PARTICLE_WEIGHT_MODE_\s+_INDIVIDUAL_PARTICLE_WEIGHT_ON_\s*$",
                definitions_text,flags=re.MULTILINE) is None:
            raise RuntimeError("finite M2 release requires generated individual particle weights ON")
        # A sampling-only fixture needs the real two-species storage layout,
        # but performs no chemistry or transport advances. It must not claim
        # I19's chemistry-off configuration when the normal build enables it.
        result["i19_chemistry_contract"] = None if args.purpose == "sampling" else {
            "photolytic_reactions": "OFF",
            "converted_in_atoms": 0.0,
            "converted_out_atoms": 0.0,
            "generated_definitions": str(generated_definitions),
            "generated_definitions_sha256": sha256(generated_definitions),
        }
        shutil.copy2(REPO / "amps", executable)
        result.update(
            status="PASS",
            executable={"path": str(executable),
                        "sha256": sha256(executable)},
        )
    except (OSError, ValueError, RuntimeError) as error:
        result["error"] = str(error)
    finally:
        cleanup_errors: list[str] = []
        # Config.pl determines the exact set of `.Assembled*` derivatives, so
        # enumerate only this uniquely named application's prefix after the
        # build. No broad repository cleanup or unrelated generated file is
        # touched.
        cleanup_paths = ((*staged, *generated_copies,
                          *tuple(REPO.glob(assembled_pattern)))
                         if cleanup_owned else ())
        for path in cleanup_paths:
            try:
                path.unlink(missing_ok=True)
            except OSError as error:
                cleanup_errors.append(f"{path}: {error}")
        if cleanup_errors:
            result["status"] = "ERROR"
            result["cleanup_errors"] = cleanup_errors
        result_path.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n",
            encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == "PASS" else 2

if __name__ == "__main__":
    raise SystemExit(main())
