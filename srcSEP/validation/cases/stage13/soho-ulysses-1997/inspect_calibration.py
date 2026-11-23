#!/usr/bin/env python3
"""Verify SolarSoft radiance/WCS exports and write their provenance receipt."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.io import readsav

from reproduce import digest, read_inputs, require

NATIVE_SHAPE = (143, 30, 21)  # SciPy reverses IDL's array axis order.
RADIANCE_UNITS = b"ERGS/CM^2/SEC/STERAD"


def inspect(source, output):
    wcs = output / (source.stem + ".wcs.csv")
    sav = output / (source.stem + ".sav")
    saved = readsav(str(sav), python_dict=True)
    require(Path(saved["input_path"].decode()).name == source.name,
            f"{sav}: recorded source is a different raster")
    x, y = saved["x"], saved["y"]
    require(x.shape == y.shape == NATIVE_SHAPE and np.isfinite(x).all() and np.isfinite(y).all(),
            f"{sav}: invalid wavelength-dependent geometry")
    with wcs.open(newline="") as stream:
        entries = list(csv.DictReader(stream))
    # Check against the full native wavelength-dependent array. The ASCII
    # grid used for field tracing is the central Mg IX pixel, rounded to six
    # places in arcseconds by the exporter. That tolerance is serialization
    # only. It is not an instrument pointing error or an acceptance gate.
    native = set()
    for entry in entries:
        ix, iy = int(entry["x_index_zero"]), int(entry["y_index_zero"])
        require(0 <= ix < 30 and 0 <= iy < 143 and (ix, iy) not in native,
                f"{wcs}: invalid or repeated native pixel")
        native.add((ix, iy))
        u, v = float(entry["x_arcsec"]), float(entry["y_arcsec"])
        require(abs(x[iy, ix, 10] - u) <= 5.1e-7 and abs(y[iy, ix, 10] - v) <= 5.1e-7,
                f"{wcs}: does not reproduce the central wavelength coordinates")
    require(len(native) == 30 * 143, f"{wcs}: incomplete spatial grid")

    # Preserve both the uncalibrated detector measurements and converted
    # physical samples. A successful export alone is insufficient if the
    # pointer returned the raw array or a different radiometry convention.
    qlds = saved["qlds"]
    windows = qlds["DETDESC"][0]
    raw_windows = saved["original_windows"]
    require(len(windows) == len(raw_windows) > 0
            and all(label.strip() == RADIANCE_UNITS for label in windows["UNITS"])
            and all(label.strip() == b"ADC" for label in raw_windows["UNITS"]),
            f"{sav}: invalid radiance or raw detector units")
    calibrated = qlds["DETDATA"][0]
    raw = saved["original_adc"]
    require(calibrated.shape == raw.shape and not np.array_equal(calibrated, raw),
            f"{sav}: detector and calibrated arrays were not both processed")
    return {"source_sha256": digest(source), "wcs_sha256": digest(wcs),
            "radiance_sha256": digest(sav), "wcs_pixels": len(native),
            "science_windows": len(windows), "original_and_processed_shape": list(raw.shape),
            "units": "erg cm^-2 s^-1 sr^-1 per spectral pixel",
            "idl_compatibility_version": saved["runtime_version"].decode()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--calibration-dir", type=Path, required=True)
    parser.add_argument("--cds-data-dir", type=Path, required=True,
                        help="SolarSoft CDS data tree read by the archived run")
    parser.add_argument("--software", required=True, help="SolarSoft release and runtime identification")
    parser.add_argument("--archive-receipt", type=Path,
                        help="optional JSON with the assessed SolarSoft archive URLs and hashes")
    args = parser.parse_args()
    receipt = args.calibration_dir / "calibration_receipt.json"
    require(not receipt.exists(), "choose a calibration directory without an existing receipt")
    pinned = read_inputs(args.input_dir)
    archives = json.loads(args.archive_receipt.read_bytes()) if args.archive_receipt else None
    inventory = {}
    # The calibrated QLDS depends on pointing, wavelength and VDS/NIS
    # response. Record the exact installed inventory rather than guessing
    # which date lookup selected each underlying record.
    for relative in ("database", "vds/cal/int", "vds/cal/wav"):
        folder = args.cds_data_dir / relative
        files = sorted(item for item in folder.rglob("*") if item.is_file())
        require(files, f"missing SolarSoft calibration inventory {folder}")
        inventory.update({str(item.relative_to(args.cds_data_dir)): digest(item) for item in files})
    rasters = {name: inspect(args.input_dir / name, args.calibration_dir)
               for name in pinned["products"] if name.startswith("s7997")}
    result = {"software": args.software, "software_source_archives": archives,
              "database_sha256": inventory, "rasters": rasters}
    temporary = receipt.with_suffix(".json.partial")
    temporary.write_text(json.dumps(result, indent=2) + "\n")
    temporary.replace(receipt)
    print(f"Verified two radiance arrays and {sum(r['wcs_pixels'] for r in rasters.values()):,} WCS rows.")


if __name__ == "__main__":
    main()
