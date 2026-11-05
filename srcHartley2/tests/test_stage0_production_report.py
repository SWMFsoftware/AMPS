#!/usr/bin/env python3
"""Validate an instrumented Stage-0 comet sampling report.

The production directory is supplied explicitly because the focused unit
suite must remain runnable in a checkout that does not retain validation
archives.  This checker deliberately recomputes the C2 number-rate and
species-flux identities instead of trusting the booleans written by C++.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def close(left: float, right: float) -> bool:
    scale = max(abs(left), abs(right), 1.0)
    return abs(left - right) <= 5.0e-13 * scale


def validate(root: Path, require_observed_c2: bool = False) -> None:
    diagnostics = root / "PT" / "diagnostics"
    measure_path = diagnostics / "cell_measures.csv"
    if not measure_path.is_file():
        raise AssertionError(f"missing exact cell-measure audit: {measure_path}")

    with measure_path.open(newline="", encoding="utf-8") as stream:
        measures = list(csv.DictReader(stream))
    if not measures:
        raise AssertionError("cell-measure audit is empty")
    radii = [float(row["radius_m"]) for row in measures]
    volumes = [float(row["measure_m3"]) for row in measures]
    if [int(row["cell_index"]) for row in measures] != list(range(len(measures))):
        raise AssertionError("cell-measure indices are not exact traversal order")
    if any(not math.isfinite(value) or value <= 0.0 for value in volumes):
        raise AssertionError("cell-measure audit contains a non-positive volume")
    if any(right <= left for left, right in zip(radii, radii[1:])):
        raise AssertionError("cell-measure audit is not radially ordered")

    json_paths = sorted(diagnostics.glob("comet_stage0.out=*.json"))
    csv_paths = sorted(diagnostics.glob("comet_stage0.out=*.csv"))
    if not json_paths or len(json_paths) != len(csv_paths):
        raise AssertionError("sampling report JSON/CSV file counts differ")

    observed_c2_selection = False
    observed_c2_carrier = False
    observed_c2_destruction = False
    for json_path, csv_path in zip(json_paths, csv_paths):
        report = json.loads(json_path.read_text(encoding="utf-8"))
        if report.get("schema") != "amps.comet.stage0-diagnostics.v1":
            raise AssertionError(f"unexpected schema in {json_path}")
        if report["sampling_length"] <= 0:
            raise AssertionError("non-positive diagnostic sampling length")

        with csv_path.open(newline="", encoding="utf-8") as stream:
            csv_rows = list(csv.DictReader(stream))
        if any(None in row for row in csv_rows):
            raise AssertionError(f"misaligned columns in {csv_path}")
        kinds = {row["record_type"] for row in csv_rows}
        if not {"metadata", "reaction", "product", "species_flux"} <= kinds:
            raise AssertionError(f"incomplete record types in {csv_path}")

        reactions = report["reactions"]
        identifiers = [reaction["id"] for reaction in reactions]
        if len(identifiers) != len(set(identifiers)):
            raise AssertionError("reaction identifiers are not unique")
        if "C2_2C" not in identifiers:
            raise AssertionError("the C2 destruction path is not exposed")

        c2_products = []
        for reaction in reactions:
            if reaction["reaction_selections"] > reaction["eligible_parent_visits"]:
                raise AssertionError("reaction selections exceed eligible channel visits")
            for product in reaction["products"]:
                if product["species"] == "C2":
                    c2_products.append((reaction, product))
                if product["allocated_slots"] + product["reused_slots"] != product["model_carriers"]:
                    raise AssertionError("product slot/carrier accounting does not close")
                difference = abs(
                    product["represented_number_rate_s-1"]
                    - product["expected_number_rate_s-1"]
                )
                accepted = difference <= product["production_tolerance_s-1"]
                if accepted != product["production_agreement"]:
                    raise AssertionError("stored product acceptance disagrees with its tolerance")

        c2_paths = {reaction["id"] for reaction, _ in c2_products}
        required_paths = {"C2H2_C2_H2", "C2H_C2_H"}
        if not required_paths <= c2_paths:
            raise AssertionError(f"missing C2 creation paths: {required_paths - c2_paths}")

        flux_by_species = {entry["species"]: entry for entry in report["species_flux"]}
        c2_flux = flux_by_species["C2"]
        expected = sum(product["expected_number_rate_s-1"] for _, product in c2_products)
        represented = sum(product["represented_number_rate_s-1"] for _, product in c2_products)
        if not close(expected, c2_flux["chemical_expected_production_s-1"]):
            raise AssertionError("C2 expected reaction products do not reconcile with species flux")
        if not close(represented, c2_flux["chemical_represented_production_s-1"]):
            raise AssertionError("C2 represented reaction products do not reconcile with species flux")
        if expected > 0.0:
            observed_c2_selection = True
            expected_tolerance = max(
                5.0 * c2_flux["chemical_production_sigma_s-1"],
                0.01 * abs(expected),
            )
            if not close(expected_tolerance, c2_flux["production_tolerance_s-1"]):
                raise AssertionError("C2 tolerance is not max(5 sigma, 1% expected)")
            if not c2_flux["production_agreement"]:
                raise AssertionError("C2 represented production fails the approved acceptance")
        if any(product["model_carriers"] > 0 for _, product in c2_products):
            observed_c2_carrier = True

        for species, entry in flux_by_species.items():
            if not close(entry["inner_incident_s-1"], entry["inner_reemitted_s-1"]):
                raise AssertionError(f"inner thermalization flux does not close for {species}")
            expected_net = (
                entry["sampled_surface_source_s-1"]
                + entry["chemical_expected_production_s-1"]
                + entry["inner_reemitted_s-1"]
                - entry["chemical_destruction_s-1"]
                - entry["inner_incident_s-1"]
                - entry["outer_escape_s-1"]
            )
            represented_net = (
                entry["sampled_surface_source_s-1"]
                + entry["chemical_represented_production_s-1"]
                + entry["inner_reemitted_s-1"]
                - entry["chemical_destruction_s-1"]
                - entry["inner_incident_s-1"]
                - entry["outer_escape_s-1"]
            )
            if not close(expected_net, entry["net_expected_flux_s-1"]):
                raise AssertionError(f"expected flux identity fails for {species}")
            if not close(represented_net, entry["net_represented_flux_s-1"]):
                raise AssertionError(f"represented flux identity fails for {species}")

        c2_loss = next(reaction for reaction in reactions if reaction["id"] == "C2_2C")
        if c2_flux["chemical_destruction_s-1"] < 0.0 or c2_flux["outer_escape_s-1"] < 0.0:
            raise AssertionError("C2 loss/escape flux is negative")
        if c2_loss["parent_destructions"] < 0:
            raise AssertionError("C2 destruction event count is negative")
        if c2_loss["parent_destructions"] > 0:
            observed_c2_destruction = True

    if require_observed_c2:
        if not observed_c2_selection:
            raise AssertionError("no physical C2-producing reaction was selected")
        if not observed_c2_carrier:
            raise AssertionError("no C2 model carrier was allocated or reused")

    print(
        f"PASS Stage-0 production reports={len(json_paths)} "
        f"cell-measures={len(measures)} c2-selected={observed_c2_selection} "
        f"c2-carrier={observed_c2_carrier} "
        f"c2-destruction={observed_c2_destruction} root={root}"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("production_root", type=Path)
    parser.add_argument(
        "--require-observed-c2",
        action="store_true",
        help="require at least one selected C2-production channel and represented carrier",
    )
    arguments = parser.parse_args()
    validate(arguments.production_root.resolve(), arguments.require_observed_c2)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
