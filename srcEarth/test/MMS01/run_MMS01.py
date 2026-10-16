#!/usr/bin/env python3
"""Run the independent QUICK or COMPLETE MMS1/HPCA product validation."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from compare_MMS01 import compare
from mms01_common import PROFILES, atomic_json, sha256_file, utc_now
from plot_MMS01 import plot
from prepare_MMS01 import prepare


HERE = Path(__file__).resolve().parent
NATIVE_SOURCE = HERE / "native_hpca_products.cpp"


def _new_campaign(parent: Path, profile: str) -> Path:
    """Allocate a unique evidence root; never overwrite a prior campaign."""
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    candidate = parent.resolve() / f"mms01_{profile.lower()}_{stamp}_pid{os.getpid()}"
    if candidate.exists():
        raise FileExistsError(f"campaign directory already exists: {candidate}")
    return candidate


def calculate(output_root: Path, amps_root: Path, cxx: str) -> dict[str, Any]:
    """Compile and execute the native adapter against the immutable prepared table."""
    prepared = json.loads((output_root / "prepared_manifest.json").read_text(encoding="utf-8"))
    native_dir = output_root / "native"
    native_dir.mkdir(exist_ok=True)
    executable = native_dir / "native_hpca_products"
    header = amps_root.resolve() / "srcEarth" / "util" / "FluxNumerics.h"
    if not header.is_file():
        raise FileNotFoundError(f"AMPS FluxNumerics header not found: {header}")
    compile_command = [
        cxx, "-std=c++11", "-O2", "-Wall", "-Wextra", "-Werror", "-pedantic",
        str(NATIVE_SOURCE), "-o", str(executable),
    ]
    compile_result = subprocess.run(compile_command, cwd=str(amps_root),
                                    text=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, check=False)
    (native_dir / "compile.log").write_text(compile_result.stdout, encoding="utf-8")
    if compile_result.returncode != 0:
        raise RuntimeError(f"native MMS01 compilation failed; see {native_dir / 'compile.log'}")
    directional = native_dir / "directional_products.csv"
    spectrum = native_dir / "spectrum_products.csv"
    run_command = [str(executable), str(output_root / "input" / "hpca_psd_native_input.dat"),
                   str(directional), str(spectrum)]
    run_result = subprocess.run(run_command, cwd=str(output_root), text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    (native_dir / "run.log").write_text(run_result.stdout, encoding="utf-8")
    if run_result.returncode != 0:
        raise RuntimeError(f"native MMS01 calculation failed; see {native_dir / 'run.log'}")
    artifacts = (directional, spectrum, native_dir / "compile.log", native_dir / "run.log")
    manifest = {
        "schema": "amps-mms01-native-run-v2", "created_utc": utc_now(),
        "profile": prepared["profile"], "return_code": run_result.returncode,
        "compiler": cxx, "compile_command": compile_command,
        "run_command": run_command,
        "production_header": {"path": str(header), "sha256": sha256_file(header)},
        "adapter_source": {"path": str(NATIVE_SOURCE), "sha256": sha256_file(NATIVE_SOURCE)},
        "prepared_input_sha256": sha256_file(output_root / "input" / "hpca_psd_native_input.dat"),
        "artifacts": {path.name: {"path": str(path.resolve()), "sha256": sha256_file(path)}
                      for path in artifacts},
    }
    atomic_json(output_root / "native_run_manifest.json", manifest)
    return manifest


def run(args: argparse.Namespace) -> tuple[Path, dict[str, Any]]:
    """Execute requested stages while preserving evidence needed by later stages."""
    stage = args.stage.upper()
    if stage == "ALL":
        output_root = (args.output_root.resolve() if args.output_root else
                       _new_campaign(args.output_parent, args.profile))
    else:
        if args.output_root is None:
            raise ValueError(f"--output-root is required for stage {stage}")
        output_root = args.output_root.resolve()
        prepared_path = output_root / "prepared_manifest.json"
        if not prepared_path.is_file():
            raise FileNotFoundError(f"missing preparation evidence: {prepared_path}")
        prepared_profile = json.loads(prepared_path.read_text(encoding="utf-8")).get("profile")
        if prepared_profile != args.profile:
            raise ValueError(
                f"requested profile {args.profile} does not match saved profile {prepared_profile}"
            )
    summary: dict[str, Any] = {"status": "NOT_COMPARED"}
    if stage in ("ALL", "PREPARE"):
        prepare(args.data_root, output_root, args.profile)
    if stage in ("ALL", "CALCULATE"):
        calculate(output_root, args.amps_root, args.cxx)
    if stage in ("ALL", "COMPARE"):
        summary = compare(output_root)
    if stage in ("ALL", "PLOT"):
        plot(output_root, tuple(args.formats))
        if stage == "PLOT":
            summary = json.loads((output_root / "comparison_summary.json").read_text())
    return output_root, summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path,
                        default=Path("/data/vtenishe/GEO/observational_validation_datasets_20261013"))
    parser.add_argument("--amps-root", type=Path, default=Path.cwd())
    parser.add_argument("--output-parent", type=Path, default=Path("test_output/MMS01"))
    parser.add_argument("--output-root", type=Path,
                        help="exact evidence root; required for a single non-ALL stage")
    parser.add_argument("--profile", choices=sorted(PROFILES), default="QUICK")
    parser.add_argument("--stage", choices=("ALL", "PREPARE", "CALCULATE", "COMPARE", "PLOT"),
                        default="ALL")
    parser.add_argument("--formats", nargs="+", default=("png", "eps"))
    parser.add_argument("--cxx", default=os.environ.get("CXX", "c++"))
    args = parser.parse_args()
    try:
        output_root, summary = run(args)
        workflow = {
            "schema": "amps-mms01-workflow-v1", "created_utc": utc_now(),
            "stage": args.stage, "profile": args.profile,
            "status": "PASS" if summary.get("status") in ("PASS", "NOT_COMPARED") else "FAIL",
            "comparison_status": summary.get("status"), "output_root": str(output_root),
        }
        # Preserve the full-campaign receipt when a reviewer later regenerates only
        # comparison tables or figures.  Stage-specific receipts document that
        # derived-evidence action without rewriting the original ALL provenance.
        receipt = (output_root / "workflow_summary.json" if args.stage == "ALL" else
                   output_root / f"workflow_{args.stage.lower()}_summary.json")
        atomic_json(receipt, workflow)
        print(f"[MMS01] stage={args.stage} profile={args.profile} status={workflow['status']}")
        print(f"[MMS01] evidence={output_root}")
        return 0 if workflow["status"] == "PASS" else 1
    except Exception as error:
        print(f"ERROR: {error}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
