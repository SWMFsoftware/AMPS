#!/usr/bin/env python3
"""Export the exact r14 sodium model point for the Moon C++ particle provider."""
import argparse
import csv
import io
from pathlib import Path
import sys

MODEL = "hplus-na3s-ground-capture-4kev-u-isotropic-v1"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", required=True, type=Path, help="installed r14 database-lock.json")
    parser.add_argument("--output", required=True, type=Path, help="new CSV path")
    args = parser.parse_args()
    lock = args.lock.resolve()
    # Import the actual installed data-access API. It validates the selected
    # release/profile/export and checks identities on queries. The C++ consumer
    # pins these projected bytes independently before using the single knot.
    sys.path.insert(0, str(lock.parent / "python"))
    from amps_chemical_database import Database
    db = Database.from_lock(lock)
    model = db.model(MODEL)
    binding = model["source_bindings"]
    source, = (row for row in db.describe(binding["dataset_id"])["files"]
               if row["id"] == binding["source_id"])
    point = source
    for key in binding["metadata_path"]:
        point = point[key]
    # No binary floating-point conversion during export. The original frame,
    # native decimal precision, and the angular and recoil policy IDs survive.
    rows = [("model", "source", "energy", "energy_unit", "energy_frame",
             "cross_section", "area_unit", "angular_model", "recoil_model"),
            (model["id"], source["id"], point["projectile_energy"], point["energy_unit"],
             point["energy_frame"], point[binding["area_field"]], binding["area_unit"],
             binding["angular_model_id"], binding["recoil_model_id"])]
    contents = io.StringIO(newline="")
    csv.writer(contents, lineterminator="\n").writerows(rows)
    with args.output.open("x", encoding="ascii", newline="") as stream:
        stream.write(contents.getvalue())


if __name__ == "__main__":
    main()
