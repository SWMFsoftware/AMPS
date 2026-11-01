#!/usr/bin/env python3
"""Validate native Enceladus CLI/configuration gates without entering MPI."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import run_tests as suite


def execute(command: list[str], cwd: Path) -> subprocess.CompletedProcess:
    return subprocess.run(command, cwd=cwd, text=True, capture_output=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path,
                        default=suite.ROOT.parent)
    parser.add_argument("--output", type=Path,
                        default=suite.ROOT / "test_output" /
                        "srcEnceladus" / "native_startup")
    args = parser.parse_args()
    executable = args.executable.resolve()
    repo_root = args.repo_root.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    outcomes: list[dict[str, str]] = []

    help_result = execute([str(executable), "--help"], repo_root)
    help_pass = (help_result.returncode == 0 and
                 "Default run deck: input/enceladus.in" in help_result.stdout)
    outcomes.append({"test": "N00_help", "status":
                     "PASS" if help_pass else "FAIL",
                     "detail": help_result.stdout + help_result.stderr})

    default_result = execute([str(executable), "--validate-input"], repo_root)
    default_text = default_result.stdout + default_result.stderr
    default_pass = (default_result.returncode == 0 and
                    "input_validation=PASS" in default_text and
                    "source_distribution=e3_h2o_co2_equal_startup" in
                    default_text and
                    "solar_activity=active" in default_text and
                    "sample_length_steps=1000" in default_text and
                    "total_iterations=10000" in default_text and
                    "flyby_density_output=on" in default_text and
                    "plasma_background_mode=none" in default_text)
    outcomes.append({"test": "N01_valid_default", "status":
                     "PASS" if default_pass else "FAIL",
                     "detail": default_text})

    with tempfile.TemporaryDirectory(prefix="enceladus_native_startup_") as name:
        temporary = Path(name)
        source = temporary / "sources.in"
        exported = execute([
            sys.executable,
            str(suite.ROOT / "standalone" / "export_forward_sources.py"),
            "--manifest", str(suite.TEST / "data" / "shape_manifest.json"),
            "--output", str(source),
        ], suite.ROOT)
        # Corrupt the projected receipt inward by one percent to prove the
        # native adapter independently enforces the declared surface mapping.
        inward_source = temporary / "inward_sources.in"
        inward_lines = []
        for line in source.read_text(encoding="utf-8").splitlines():
            if line.startswith(("x_m=", "y_m=", "z_m=")):
                key, value = line.split("=", 1)
                line = f"{key}={0.99 * float(value):.17g}"
            inward_lines.append(line)
        inward_source.write_text(
            "\n".join(inward_lines) + "\n", encoding="utf-8")
        runtime_input = temporary / "enceladus.in"
        suite.write_transport_only_input(runtime_input, inward_source.name)
        geometry_result = execute(
            [str(executable), "--input", str(runtime_input)], repo_root)
        geometry_text = geometry_result.stdout + geometry_result.stderr
        geometry_pass = (exported.returncode == 0 and
                         geometry_result.returncode != 0 and
                         "does not lie on the configured spherical boundary "
                         "252100 m" in
                         geometry_text and
                         "local_vertical_to_reference_sphere" in geometry_text)
        outcomes.append({"test": "N02_stripe_surface_gate", "status":
                         "PASS" if geometry_pass else "FAIL",
                         "detail": geometry_text})

        cache_input = temporary / "missing_cache_directory.in"
        suite.write_transport_only_input(cache_input, source.name)
        cache_input.write_text(
            cache_input.read_text(encoding="utf-8").replace(
                "plasma_cache_directory=.",
                "plasma_cache_directory=missing_plasma_cache"),
            encoding="utf-8")
        cache_result = execute(
            [str(executable), "--input", str(cache_input)], repo_root)
        cache_text = cache_result.stdout + cache_result.stderr
        cache_pass = (cache_result.returncode != 0 and
                      "plasma_cache_directory is not an accessible directory" in
                      cache_text and "missing_plasma_cache" in cache_text)
        outcomes.append({"test": "N03_plasma_cache_directory_gate", "status":
                         "PASS" if cache_pass else "FAIL",
                         "detail": cache_text})

    generated = repo_root / "build" / "pic" / "picGlobal.dfn"
    mover_detail = "generated picGlobal.dfn is missing"
    mover_pass = False
    if generated.is_file():
        definitions = [line.strip() for line in
                       generated.read_text(encoding="utf-8").splitlines()
                       if line.startswith(
                           "#define _PIC_PARTICLE_MOVER_LEGACY_SETTINGS_")]
        mover_detail = "\n".join(definitions) + "\n"
        mover_pass = bool(definitions) and definitions[-1].endswith(
            "_PIC_MODE_OFF_")
    outcomes.append({"test": "N04_pointer_mover_generated", "status":
                     "PASS" if mover_pass else "FAIL", "detail": mover_detail})

    chemistry_detail = "generated picGlobal.dfn is missing"
    chemistry_pass = False
    if generated.is_file():
        lines = [line.strip() for line in
                 generated.read_text(encoding="utf-8").splitlines()]
        modes = [line for line in lines if line.startswith(
            "#define _PIC_USER_PARTICLE_PROCESSING__MODE_")]
        callbacks = [line for line in lines if line.startswith(
            "#define _PIC_USER_PARTICLE_PROCESSING__FUNCTION_")]
        chemistry_detail = "\n".join(modes[-1:] + callbacks[-1:]) + "\n"
        chemistry_pass = (bool(modes) and modes[-1].endswith("_PIC_MODE_ON_") and
                          bool(callbacks) and
                          "ForwardModel::AMPS::ProcessParticle" in callbacks[-1])
    outcomes.append({"test": "N05_chemistry_processing_generated", "status":
                     "PASS" if chemistry_pass else "FAIL",
                     "detail": chemistry_detail})

    passed = all(item["status"] == "PASS" for item in outcomes)
    receipt = {
        "runner": str(Path(__file__).resolve()),
        "executable": str(executable),
        "repo_root": str(repo_root),
        "status": "PASS" if passed else "FAIL",
        "results": outcomes,
    }
    (output / "test_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    for item in outcomes:
        print(f"{item['status']} {item['test']}")
        if item["status"] != "PASS":
            print(item["detail"].rstrip())
    print(f"SUMMARY status={receipt['status']} receipt={output/'test_receipt.json'}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
