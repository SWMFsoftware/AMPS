#!/usr/bin/env python3
"""Compile/run portable native-boundary and corridor/solar regressions.

This runner uses actual native evaluators and shared geometry, with controlled
input records. It provides component evidence; it cannot qualify linked AMPS
allocation or MPI behavior. It works with the supplied partial source package.
"""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess

ROOT = Path(__file__).resolve().parents[1]

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cxx", default=os.environ.get("CXX", "g++"))
    parser.add_argument("--model-dir", type=Path,
                        default=ROOT.parent / "src/models/sep_coronal_cme")
    parser.add_argument("--output-dir", type=Path,
                        default=ROOT / "test_output/native-boundary")
    args = parser.parse_args()
    model = args.model_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    compiler = shlex.split(args.cxx)
    if not compiler:
        parser.error("--cxx must name a compiler")
    # This archive is independently buildable, unlike the omitted full AMPS
    # framework and older SEP shared coefficient/kernel build files.
    subprocess.run(["make", "-C", str(model), "-j8", "lib",
                    "CXX=" + args.cxx], check=True)
    binary = output / "native_boundary_regression"
    report = output / "native-boundary.json"
    if report.exists():
        report.unlink()
    # Section garbage collection removes unused mesh/storage methods that
    # belong to the omitted complete runtime. The functions exercised here
    # compile directly from their production sources without replacement APIs.
    command = compiler + ["-O2", "-std=c++17", "-Wall", "-Wextra",
        "-Wpedantic", "-Werror", "-ffunction-sections", "-fdata-sections",
        "-I" + str(model / "include"), "-I" + str(model.parent / "sep_common"),
        str(ROOT / "test/native_boundary_regression.cpp"),
        str(ROOT / "validation/coronal_cme_application_test.cpp"),
        str(ROOT / "mesh/mesh_model.cpp"),
        str(ROOT / "core/parker_geometry.cpp"),
        str(ROOT / "runtime/standalone_command_line.cpp"),
        str(model / "build/libsep_coronal_cme.a"),
        "-Wl,--gc-sections", "-o", str(binary)]
    subprocess.run(command, check=True)
    subprocess.run([str(binary), str(report)], check=True)
    data = json.loads(report.read_text(encoding="utf-8"))
    expected = {"mode": "parker-tube", "plan_installed": True,
                "pruning_applied": True, "allocation_verified": True,
                "planned_active_leaves": 100, "planned_inactive_leaves": 900,
                "planned_solar_interior_leaves": 8, "allocated_blocks": 100}
    if data.get("active_region") != expected:
        raise RuntimeError("native JSON lost active-region qualification evidence")
    if data.get("schema") != "srcsep-component-tests-v1" or \
            data["results"][0]["status"] != "PASS":
        raise RuntimeError("native JSON schema/status regression")
    print("native JSON active-region evidence: PASS")
    return 0

if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, subprocess.CalledProcessError) as error:
        print("native boundary regression ERROR:", error)
        raise SystemExit(2)
