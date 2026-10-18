#!/usr/bin/env python3
"""Validate LiouvilleCOPS caches or register C++ reproduction of a fit.

Examples
--------
Validate a static fitting base without changing it::

  python3 tools/validate_cops_run.py output_fit_base

After running the generated ``cops_fit.in``, compare its aggregate pressures
with the Python fit prediction and update only the fit's verification state::

  python3 tools/validate_cops_run.py --fit output_fit_new \
      --verification output_fit_new_verify

This tool checks numerical/product identity.  It never changes
``scientific_review`` and does not assert calibration or physical validity.
"""

import argparse
import math
from pathlib import Path
import sys

from cops_paths import product_input
from cops_io import (FIT_SCHEMA, RUN_SCHEMA, ProductError, file_identity,
                     load_manifest, read_model_table, utc_now,
                     validate_fit_base, verify_identity, write_json_atomic)


def _dependencies(manifest):
    """Return ordered (role, hash) identities for cross-run comparison."""
    return [(entry.get("role"), entry.get("file", {}).get("sha256"))
            for entry in manifest.get("dependencies", [])]


def _verify_declared_outputs(fit_directory, manifest):
    """Detect edits to top-level fit outputs made after manifest creation."""
    for name, identity in manifest.get("outputs", {}).items():
        if name == "lambda_directories" or not isinstance(identity, dict):
            continue
        # Manifests store absolute paths, but require that the declared file is
        # still the corresponding member of this fit directory.
        expected = (fit_directory / name).resolve()
        if Path(identity.get("path", "")).resolve() != expected:
            raise ProductError("manifest output path for %s is not inside the fit" % name)
        verify_identity(identity, "fit output " + name)


def validate_base(run):
    audit = validate_fit_base(run)
    print("base: %s" % audit["run"])
    print("manifest: %s" % ("legacy/missing" if audit["legacy_manifest"] else "valid"))
    print("table points: %d" % len(audit["point_ids"]))
    print("elements: %d faces, %d basis functions" %
          (audit["elements"]["nfaces"], audit["elements"]["nbasis"]))
    print("kernels: %d files, %d records, %d bytes" %
          (audit["kernels"]["count"], audit["kernels"]["records"],
           audit["kernels"]["size_bytes"]))
    print("cache: verified")


def register_verification(fit, verification, tolerance):
    fit, verification = Path(fit).resolve(), Path(verification).resolve()
    manifest_path = fit / "fit_manifest.json"
    manifest = load_manifest(manifest_path, FIT_SCHEMA)
    if manifest.get("states", {}).get("process") != "complete":
        raise ProductError("fit process is not complete")
    _verify_declared_outputs(fit, manifest)

    verification_manifest_path = verification / "run_manifest.json"
    verification_manifest = load_manifest(verification_manifest_path, RUN_SCHEMA)
    if verification_manifest.get("states", {}).get("process") != "complete":
        raise ProductError("verification run process is not complete")
    if verification_manifest.get("states", {}).get("cache") not in ("verified", "not_applicable"):
        raise ProductError("verification run cache is not verified")

    # When the parent run was manifested, require the complete ordered set of
    # shape/observation/SPICE/HK dependencies to be byte-identical.  Activity
    # coefficients legitimately differ and are excluded by role.
    parent = manifest.get("parent", {})
    geometry_state = "not_verifiable_from_legacy_parent"
    if parent.get("kind") == "manifested_run":
        verify_identity(parent["manifest"], "parent run manifest")
        parent_manifest = load_manifest(parent["manifest"]["path"], RUN_SCHEMA)
        excluded_roles = {"activity_coefficients"}
        original = [item for item in _dependencies(parent_manifest)
                    if item[0] not in excluded_roles]
        reproduced = [item for item in _dependencies(verification_manifest)
                      if item[0] not in excluded_roles]
        if original != reproduced:
            raise ProductError("verification run does not use the parent's shape/observation/SPICE/HK bytes")
        geometry_state = "verified"

    required = ("iPoint", "UTC", "inside", "pNG", "pRG")
    _, _, predicted, pc = read_model_table(fit / "COPS_model.dat", required)
    _, _, checked, cc = read_model_table(verification / "COPS_model.dat", required)
    checked_by_id = {int(row[cc["iPoint"]]): row for row in checked}
    if len(checked_by_id) != len(checked):
        raise ProductError("duplicate points in verification table")
    maximum = {"NG": 0.0, "RG": 0.0}
    compared = 0
    for row in predicted:
        point = int(row[pc["iPoint"]])
        other = checked_by_id.pop(point, None)
        if other is None:
            raise ProductError("verification table is missing iPoint %d" % point)
        if row[pc["UTC"]] != other[cc["UTC"]]:
            raise ProductError("UTC mismatch at iPoint %d" % point)
        if int(row[pc["inside"]]):
            continue
        for gauge in ("NG", "RG"):
            predicted_pressure = float(row[pc["p" + gauge]])
            checked_pressure = float(other[cc["p" + gauge]])
            if not (math.isfinite(predicted_pressure) and math.isfinite(checked_pressure)
                    and predicted_pressure > 0.0 and checked_pressure > 0.0):
                raise ProductError("non-positive/non-finite %s pressure at iPoint %d" %
                                   (gauge, point))
            difference = abs(predicted_pressure / checked_pressure - 1.0)
            maximum[gauge] = max(maximum[gauge], difference)
        compared += 1
    if checked_by_id:
        raise ProductError("verification table has points absent from the fit prediction")
    if max(maximum.values()) > tolerance:
        raise ProductError("C++ pressure reproduction exceeds tolerance: NG %.3e, RG %.3e > %.3e" %
                           (maximum["NG"], maximum["RG"], tolerance))

    manifest["states"]["cpp_prediction_agreement"] = "verified"
    manifest["updated_utc"] = utc_now()
    manifest["cpp_verification"] = {
        "run_manifest": file_identity(verification_manifest_path),
        "model_table": file_identity(verification / "COPS_model.dat"),
        "compared_outside_points": compared,
        "relative_tolerance": tolerance,
        "max_relative_difference": maximum,
        "parent_dependency_identity": geometry_state,
    }
    write_json_atomic(manifest_path, manifest)
    print("fit: %s" % fit)
    print("verification: %s" % verification)
    print("C++ prediction agreement: verified (%d outside points; NG %.3e, RG %.3e)" %
          (compared, maximum["NG"], maximum["RG"]))
    print("scientific review: %s (unchanged)" % manifest["states"]["scientific_review"])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run", nargs="?", help="static fitting-base directory")
    parser.add_argument("--fit", help="fit directory containing fit_manifest.json")
    parser.add_argument("--verification", help="C++ run generated from the fit's cops_fit.in")
    parser.add_argument("--tolerance", type=float, default=1e-6,
                        help="maximum relative NG/RG pressure difference (default: %(default)g)")
    args = parser.parse_args(argv)
    try:
        if args.fit or args.verification:
            if not args.fit or not args.verification or args.run:
                parser.error("--fit and --verification must be used together, without RUN")
            if args.tolerance <= 0.0:
                parser.error("--tolerance must be positive")
            register_verification(product_input(args.fit), product_input(args.verification), args.tolerance)
        elif args.run:
            validate_base(product_input(args.run))
        else:
            parser.error("provide RUN, or --fit FIT --verification RUN")
        return 0
    except (OSError, ProductError, ValueError) as exc:
        print("Error: %s" % exc, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
