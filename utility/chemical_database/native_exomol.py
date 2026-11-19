"""Lossless lexical ExoMol projections, separate from physical event models.

The official definitions file controls OH field names/count and source version.
Published wavelength/cross-section columns control the selected H2O tables.
Original numeric tokens/order survive exactly. Negative measurement samples,
source infinities/qualifiers and missing state/formation conventions are never
silently corrected. These CSVs provide source access, not a reaction sampler.
"""
from __future__ import annotations

import argparse
import bz2
import csv
from decimal import Decimal, InvalidOperation
import hashlib
import json
from pathlib import Path


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def project(source, output, kind, definition=None):
    """Project exactly declared source tokens into a NEW CSV and QA receipt.

    Quantities are never converted: photo wavelengths stay nm, cross sections
    cm2/molecule, OH energies cm-1, and source lifetimes s. ``definition`` must
    be the actual pinned JSON for the states file; it cannot create molecular
    formation energies, products or a physical rest mass crosswalk.
    """
    source, output = Path(source).resolve(), Path(output).resolve()
    if output.exists():
        raise ValueError("refuse to overwrite an existing ExoMol projection")
    raw = source.read_bytes()
    if source.suffix == ".bz2":
        raw = bz2.decompress(raw)
    lines = raw.decode("utf-8").splitlines()
    if not lines or any(not line.strip() for line in lines):
        raise ValueError("empty source/blank row requires explicit source interpretation")
    rows = [line.split() for line in lines]
    units, limitations = {}, []
    if kind == "photo":
        columns = ["wavelength_nm", "cross_section_cm2_molecule"]
        units = {"wavelength_nm": "nm", "cross_section_cm2_molecule": "cm2/molecule"}
        if any(len(row) != 2 for row in rows):
            raise ValueError("published two-column photo schema does not match native source")
        try:
            numbers = [[Decimal(token) for token in row] for row in rows]
        except InvalidOperation as error:
            raise ValueError("non-numeric source photo token") from error
        if any(not x.is_finite() for row in numbers for x in row):
            raise ValueError("nonfinite source photo value")
        if any(row[0] <= 0 for row in numbers):
            raise ValueError("nonpositive published wavelength")
        if any(b[0] <= a[0] for a, b in zip(numbers, numbers[1:])):
            raise ValueError("nonincreasing native wavelengths; do not sort source data")
        negative = sum(row[1] < 0 for row in numbers)
        limitations += ["Published absorption/cross-section samples do not supply channel products or daughter kinetics.",
                        "Native signed samples are preserved without clipping, smoothing or uncertainty assignment."]
        qa = {"negative_sample_count": negative, "zero_sample_count": sum(row[1] == 0 for row in numbers),
              "first_wavelength_native": rows[0][0], "last_wavelength_native": rows[-1][0],
              "nonnegative_probability_table": negative == 0}
    elif kind == "states":
        if definition is None:
            raise ValueError("OH states require the original source definition JSON")
        definition = Path(definition).resolve()
        definitions = json.loads(definition.read_text())
        state = definitions["dataset"]["states"]
        columns = [item["name"] for item in state["states_file_fields"]]
        if len(set(columns)) != len(columns) or any(len(row) != len(columns) for row in rows):
            raise ValueError("native OH columns disagree with supplied definition; no guessed column offsets")
        if len(rows) != state["number_of_states"]:
            raise ValueError("native OH state count disagrees with supplied definition")
        if "ID" not in columns or "E" not in columns:
            raise ValueError("required native state ID/energy fields absent")
        ids = [row[columns.index("ID")] for row in rows]
        if any(not value.isdigit() or int(value) <= 0 for value in ids) or len(ids) != len(set(ids)):
            raise ValueError("native OH state IDs must be distinct positive integers")
        try:
            energy = [Decimal(row[columns.index("E")]) for row in rows]
        except InvalidOperation as error:
            raise ValueError("invalid OH source energy token") from error
        if any(not value.is_finite() for value in energy):
            raise ValueError("nonfinite OH source energy")
        units = {"ID": "1", "E": "cm-1", "gtot": "1", "J": "1", "unc": "cm-1",
                 "tau": "s", "gfactor": "1", "v": "1", "hunda:Lambda": "1",
                 "hunda:Sigma": "1", "hunda:Omega": "1", "Auxiliary:Ecal": "cm-1"}
        units = {key: value for key, value in units.items() if key in columns}
        qa = {"definition_sha256": _sha(definition), "source_dataset_version": definitions["dataset"]["version"],
              "declared_states": state["number_of_states"], "observed_states": len(rows),
              "minimum_energy_native_cm-1": str(min(energy)), "maximum_energy_native_cm-1": str(max(energy))}
        # The current source definition and actual E-column upper extent differ.
        # Preserve both: do not silently crop states, alter the header or assume
        # this metadata field supplies a validated physical energy cutoff.
        if "max_energy" in state:
            declared = Decimal(str(state["max_energy"]))
            qa["definition_max_energy_native_cm-1"] = str(declared)
            qa["maximum_energy_matches_definition"] = declared == max(energy)
            if declared != max(energy):
                limitations.append("Definition max_energy differs from the native E-column maximum; no physical cutoff is inferred.")
        limitations += ["E is the source spectroscopic excitation coordinate, not molecular formation energy.",
                        "Catalogue mass is not admitted as an isotope/state physical-rest-mass crosswalk.",
                        "Source lifetime/uncertainty tokens are preserved; no missing value or complete cascade is inferred."]
    elif kind == "partition":
        columns = ["temperature_K", "partition_function_native"]
        units = {"temperature_K": "K", "partition_function_native": "1"}
        if any(len(row) != 2 for row in rows):
            raise ValueError("partition-function source schema differs from two-column definition")
        try:
            pairs = [[Decimal(token) for token in row] for row in rows]
        except InvalidOperation as error:
            raise ValueError("non-numeric partition-function token") from error
        if any(not value.is_finite() or value <= 0 for row in pairs for value in row):
            raise ValueError("nonpositive/nonfinite native partition quantity")
        if any(b[0] <= a[0] for a, b in zip(pairs, pairs[1:])):
            raise ValueError("nonincreasing source partition temperatures; no sorting")
        qa = {}
        limitations += ["The source population/degeneracy convention is retained; no initial-state population is selected."]
    else:
        raise ValueError("unsupported ExoMol lexical projection")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)
    # Re-read the actual bytes through an independent CSV decoder before any
    # source-access table can be published. This proves lexical identity only.
    with output.open(newline="", encoding="utf-8") as stream:
        reread = list(csv.reader(stream))
    if reread != [columns] + rows:
        raise ValueError("ExoMol lexical projection changed original source tokens")
    return {"status": "PASS", "scope": "native source/schema/lexical QA only", "kind": kind,
            "source_sha256": _sha(source), "output_sha256": _sha(output), "row_count": len(rows),
            "columns": columns, "units": units, "limitations": limitations, "physical_activation": False,
            "native_token_identity": True, "source_checks": qa,
            "recipe_sha256": _sha(__file__)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--kind", choices=("photo", "states", "partition"), required=True)
    parser.add_argument("--definition", type=Path)
    args = parser.parse_args()
    result = project(args.source, args.output, args.kind, args.definition)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
