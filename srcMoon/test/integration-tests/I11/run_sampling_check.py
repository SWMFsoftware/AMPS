#!/usr/bin/env python3
"""Execute I11's deterministic production sampling-window prerequisite.

See SAMPLING_CONTRACT.md for the independent constant-window identity. This
runner records a storage-wiring result, never a spatial-convergence result.
The supplied mesh is the approved 480 km child derived from the M1 template.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from integration_common import (REPO, derive_i11_input, manifest,
                                parse_receipt, run_monitored, sha256)


def run(executable: Path, template: Path, output: Path, ranks: int,
        build_receipt: Path | None = None) -> int:
    """Preserve one linked storage check and its exact executable/input hashes."""
    output = output.resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"refusing to overwrite sampling evidence: {output}")
    output.mkdir(parents=True, exist_ok=True)
    runtime = output / "runtime"
    runtime.mkdir()
    input_path = output / "application.in"
    derive_i11_input(template.resolve(), input_path, 480)
    surface_template = None
    if build_receipt is not None:
        # A compiled Ar species needs its documented positive source rate for
        # ordinary AMPS weight initialization, even though this storage check
        # performs no transport. Reuse the complete approved M2 receipt and
        # its literature parameters; do not invent a dummy weight/source.
        # Sphere geometry and simple-cosine thermal mode avoid loading LOLA or
        # spinning up thermal columns in a storage-only regression.
        surface_template = REPO / "srcMoon/examples/lola_surface.in"
        derive_i11_input(surface_template, input_path, 480)
        text = input_path.read_text()
        for key, value in {"surface_geometry": "sphere",
                           "surface_temperature_mode": "simple-cosine",
                           "argon_source_enabled": "true",
                           "simulation_environment_mode": "frozen",
                           "surface_inventory_checkpoint_file": "none"}.items():
            text, count = re.subn(rf"(?m)^({key}\s*=\s*).+$",
                                 lambda match: match[1] + value, text)
            if count != 1:
                raise ValueError(f"storage fixture expected exactly one {key}")
        input_path.write_text(text)
    receipt = output / "receipt.txt"
    command = ["mpiexec", "-n", str(ranks), str(executable.resolve()),
               "--moon-integration-test", "I11", "--i11-sampling-check",
               "--integration-output", str(receipt), "-input", str(input_path)]
    evidence = manifest("I11", executable, command, 731993, input_path, ranks)
    evidence["purpose"] = "sampling-window-storage-prerequisite"
    evidence["species_fixture"] = "compiled layout; count verified by runtime receipt"
    if surface_template is not None:
        evidence["inputs"][str(surface_template)] = sha256(surface_template)
    if build_receipt is not None:
        build = json.loads(build_receipt.read_text())
        if (build.get("status") != "PASS" or
                build.get("purpose") != "sampling" or
                build["executable"]["sha256"] != sha256(executable)):
            raise ValueError("sampling-only build receipt/executable identity mismatch")
        evidence["compiled_build_receipt"] = {
            "path": str(build_receipt.resolve()), "sha256": sha256(build_receipt),
            "build": build}
    (output / "run_manifest.json").write_text(
        json.dumps(evidence, indent=2, sort_keys=True) + "\n")
    result = {"id": "I11", "scope": "sampling-window-storage-prerequisite",
              "formal_i11_status": "NOT EVALUATED", "status": "ERROR"}
    try:
        environment = dict(os.environ, OMP_NUM_THREADS="1")
        code = run_monitored(command, runtime, environment,
                             output / "stdout.log", output / "stderr.log",
                             output / "process_resource_metrics.csv")
        if code != 0:
            raise RuntimeError(f"linked sampling fixture exited {code}")
        data = parse_receipt(receipt)
        errors = {name: float(data[name]) for name in (
            "sampling_check_completed_error", "sampling_check_collecting_error",
            "sampling_check_callback_error")}
        passed = (data.get("test_id") == "I11" and
                  data.get("mpi_ranks") == ranks and
                  data.get("sampling_check_windows") == 9 and
                  int(data["sampling_check_species"]) == (2 if build_receipt else 1) and
                  all(value == 0 for value in errors.values()))
        result.update(status="PASS" if passed else "FAIL", metrics=data,
                      receipt_sha256=sha256(receipt))
        evidence["compiled_species_count"] = int(data["sampling_check_species"])
        evidence["receipt_sha256"] = sha256(receipt)
        (output / "run_manifest.json").write_text(
            json.dumps(evidence, indent=2, sort_keys=True) + "\n")
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        result["reason"] = str(error)
    (output / "result.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"{result['status']} I11 sampling storage: {output}", flush=True)
    return {"PASS": 0, "FAIL": 1, "ERROR": 2}[result["status"]]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amps", type=Path, required=True)
    parser.add_argument("--input", type=Path,
                        default=REPO / "srcMoon/examples/m1_memory200.in")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--mpi-ranks", type=int, default=8)
    parser.add_argument("--build-receipt", type=Path,
                        help="hashed sampling-only Na+Ar build receipt, when using that fixture")
    args = parser.parse_args()
    if args.mpi_ranks < 1:
        parser.error("--mpi-ranks must be positive")
    raise SystemExit(run(args.amps, args.input, args.output_dir, args.mpi_ranks,
                         args.build_receipt))
