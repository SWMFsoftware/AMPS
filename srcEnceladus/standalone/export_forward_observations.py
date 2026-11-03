#!/usr/bin/env python3
"""Export traceable E7 CO2 and E21 native-count forward-model inputs.

This exporter deliberately keeps two observation contracts separate.  E7
uses the calibrated mass-44 density vector distributed in ``Flybys2.mat`` and
is matched count-for-count to official PDS Level-1A geometry.  E21 retains the
published open-source H2O/CO2 detector counts in their native observable; it
does not invent a scalar count-to-density conversion.  The latter is suitable
for a labelled context panel until the velocity/attitude-dependent SIMION
response is accumulated from the AMPS phase-space distribution.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.io import loadmat


DEFAULT_VALIDATION_ROOT = Path("/data/vtenishe/cassini_validation_data")
DEFAULT_OUTPUT = Path(__file__).resolve().parent / "observations"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty observation table: {path}")
    with path.open("w", newline="", encoding="utf-8") as stream:
        # AMPS' deliberately small runtime CSV reader treats a carriage return
        # as data in the final header token.  Emit the repository's canonical
        # LF records explicitly so the release is identical on every host and
        # the last uncertainty column retains its declared name.
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]),
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")


def pds_mass44_rows(raw_root: Path) -> tuple[list[dict[str, object]], list[dict[str, str]]]:
    selected: list[dict[str, object]] = []
    lineage: list[dict[str, str]] = []
    for path in sorted((raw_root / "E7").glob("*.CSV")):
        lineage.append({"path": str(path.resolve()), "sha256": sha256(path)})
        with path.open(newline="", encoding="utf-8") as stream:
            for row in csv.DictReader(stream):
                try:
                    time_s = float(row["time_ca"]) / 1000.0
                    mass = float(row["mass_per_charge"])
                except (TypeError, ValueError):
                    continue
                if (row["source"].strip().lower() != "csn" or
                        abs(mass - 44.0) > 1.0e-8 or
                        not -70.0 <= time_s <= 70.0):
                    continue
                selected.append({
                    "time_ca_s": time_s,
                    "counts": float(row["c1counts"]),
                    "x_m": 1000.0 * float(row["sc_pos_t_x"]),
                    "y_m": 1000.0 * float(row["sc_pos_t_y"]),
                    "z_m": 1000.0 * float(row["sc_pos_t_z"]),
                    "ram_angle_deg": float(row["ram_angle_t"]),
                    "source_file": path.name,
                })
    selected.sort(key=lambda item: float(item["time_ca_s"]))
    if not selected:
        raise RuntimeError("no E7 closed-source mass-44 PDS rows found")
    return selected, lineage


def closest_indices(reference: np.ndarray, query: np.ndarray) -> np.ndarray:
    index = np.searchsorted(reference, query)
    index = np.clip(index, 1, len(reference) - 1)
    previous = np.abs(reference[index - 1] - query) < np.abs(reference[index] - query)
    return index - previous.astype(int)


def export_e7(raw_inms: Path, flybys_mat: Path, output_dir: Path) -> dict[str, object]:
    flybys = loadmat(flybys_mat, squeeze_me=True,
                     struct_as_record=False)["Flybys"]
    e7 = next(item for item in np.atleast_1d(flybys)
              if str(item.drctry) == "E7")
    mass44 = next(item for item in np.atleast_1d(e7.datVls)
                  if int(item.mass) == 44)
    published_time = np.asarray(mass44.dat[:, 0], dtype=float)
    counts = np.asarray(mass44.dat[:, 1], dtype=float)
    density = np.asarray(mass44.Dens, dtype=float)
    if (published_time.shape != counts.shape or counts.shape != density.shape or
            len(counts) != 43 or np.any(counts < 0.0) or np.any(density < 0.0)):
        raise ValueError("invalid E7 mass-44 vectors in Flybys2.mat")

    pds, lineage = pds_mass44_rows(raw_inms)
    raw_time = np.asarray([float(row["time_ca_s"]) for row in pds])
    raw_counts = np.asarray([float(row["counts"]) for row in pds])
    candidates = np.linspace(-0.8, 0.8, 1601)
    scores = []
    for lag in candidates:
        index = closest_indices(raw_time, published_time + lag)
        scores.append((
            -float(np.mean(raw_counts[index] == counts)),
            float(np.mean(np.abs(raw_counts[index] - counts))),
            float(np.sqrt(np.mean((raw_time[index] - published_time - lag) ** 2))),
            float(lag),
        ))
    lag = min(scores)[3]
    index = closest_indices(raw_time, published_time + lag)
    if not np.array_equal(raw_counts[index], counts):
        raise RuntimeError("E7 Flybys2/PDS mass-44 count sequences do not match")
    residual = raw_time[index] - published_time - lag
    if np.max(np.abs(residual)) > 0.02 or len(set(index.tolist())) != len(index):
        raise RuntimeError("E7 Flybys2/PDS geometry match is not unique within 0.02 s")

    nonzero = counts > 0.0
    sensitivity = np.empty_like(density)
    sensitivity[nonzero] = density[nonzero] / counts[nonzero]
    sensitivity[~nonzero] = np.interp(
        published_time[~nonzero], published_time[nonzero], sensitivity[nonzero])
    sigma = sensitivity * np.sqrt(np.maximum(counts, 1.0))
    rows: list[dict[str, object]] = []
    for number, raw_index in enumerate(index):
        raw = pds[int(raw_index)]
        rows.append({
            "flyby": "E7",
            "time_ca_s": f"{published_time[number]:.6f}",
            "mass_per_charge_u": "44",
            "species_interpretation": "CO2",
            "density_m3": f"{density[number]:.12g}",
            "sigma_counting_m3": f"{sigma[number]:.12g}",
            "counts_per_integration": f"{counts[number]:.12g}",
            "spacecraft_x_m": f"{float(raw['x_m']):.6f}",
            "spacecraft_y_m": f"{float(raw['y_m']):.6f}",
            "spacecraft_z_m": f"{float(raw['z_m']):.6f}",
            "ram_angle_deg": f"{float(raw['ram_angle_deg']):.9g}",
            "density_source": "SwRI_Flybys2_mass44_Dens",
            "source_file": str(raw["source_file"]),
        })
    destination = output_dir / "e7_co2_density.csv"
    write_csv(destination, rows)
    metadata = {
        "event_id": "E7",
        "observable": "CO2 (mass-44) local number density",
        "rows": len(rows),
        "density_source": str(flybys_mat.resolve()),
        "density_source_sha256": sha256(flybys_mat),
        "density_field": "Flybys[E7].datVls[mass==44].Dens",
        "position_source": "exact-count-matched official PDS INMS Level-1A rows",
        "published_to_pds_time_lag_s": lag,
        "published_to_pds_max_abs_time_residual_s": float(np.max(np.abs(residual))),
        "point_sigma": "Poisson counting sigma propagated through local density/count factor",
        "common_calibration_fractional_sigma_not_in_point_sigma": 0.15,
        "pds_lineage": lineage,
        "reference": "Perry et al. (2015), Icarus 257, 139-162",
    }
    metadata_path = destination.with_suffix(".metadata.json")
    write_json(metadata_path, metadata)
    # Product names are relative to the manifest directory so a checked-in
    # release stays usable after the repository is moved to another host.
    return {"path": destination.name, "sha256": sha256(destination),
            "metadata": metadata_path.name, "rows": len(rows)}


def export_e21_species(swri: Path, species: str,
                       trajectory: np.ndarray, output_dir: Path) -> dict[str, object]:
    source = swri / f"E21{species}osnb.mat"
    data = np.asarray(loadmat(source, squeeze_me=True)["SignalMapDat"], dtype=float)
    if data.ndim != 2 or data.shape[1] != 2 or np.any(data[:, 1] < 0.0):
        raise ValueError(f"invalid E21 {species} count table")
    # The supplied trajectory spans -1360:+2238 s. H2O counts begin earlier;
    # retain only rows with explicitly supplied geometry rather than
    # extrapolating a spacecraft path.
    keep = ((data[:, 0] >= trajectory[0, 0]) &
            (data[:, 0] <= trajectory[-1, 0]))
    data = data[keep]
    positions = np.column_stack([
        np.interp(data[:, 0], trajectory[:, 0], trajectory[:, component])
        for component in (1, 2, 3)
    ])
    if np.any(np.linalg.norm(positions, axis=1) <= 252100.0):
        raise RuntimeError(f"E21 {species} interpolation produced a subsurface point")
    mass = 18 if species == "H2O" else 44
    rows = []
    for (time_s, counts), position in zip(data, positions):
        rows.append({
            "event_id": "E21",
            "time_ca_s": f"{time_s:.7f}",
            "spacecraft_x_m": f"{position[0]:.9f}",
            "spacecraft_y_m": f"{position[1]:.9f}",
            "spacecraft_z_m": f"{position[2]:.9f}",
            "species": species,
            "mass_per_charge_u": mass,
            "counts_per_integration": f"{counts:.12g}",
            "poisson_sigma_counts": f"{math.sqrt(max(counts, 1.0)):.12g}",
        })
    destination = output_dir / f"e21_{species.lower()}_osnb_counts.csv"
    write_csv(destination, rows)
    trajectory_path = swri / "E21_trajectoryB.mat"
    metadata = {
        "event_id": "E21",
        "species": species,
        "mass_per_charge_u": mass,
        "observable": "INMS open-source neutral-beaming counts per integration",
        "comparison_role": "native-count context; no density conversion and no goodness-of-fit",
        "rows": len(rows),
        "source": str(source.resolve()),
        "source_sha256": sha256(source),
        "trajectory_source": str(trajectory_path.resolve()),
        "trajectory_source_sha256": sha256(trajectory_path),
        "geometry": "linear interpolation of E21_trajectory columns 1:3 in IAU_ENCELADUS metres",
        "excluded_count_rows_without_trajectory": int(np.count_nonzero(~keep)),
        "required_for_quantitative_comparison": (
            "accumulate the species phase-space distribution through the supplied "
            "velocity-, voltage-, and attitude-dependent SIMION response"
        ),
    }
    metadata_path = destination.with_suffix(".metadata.json")
    write_json(metadata_path, metadata)
    return {"path": destination.name, "sha256": sha256(destination),
            "metadata": metadata_path.name, "rows": len(rows)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation-root", type=Path,
                        default=DEFAULT_VALIDATION_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    root = args.validation_root.resolve()
    raw_inms = root / "data" / "raw" / "inms"
    swri = (root / "data" / "raw" / "external" /
            "swri_inms_modeling" / "extracted")
    flybys = swri / "Flybys2.mat"
    trajectory_path = swri / "E21_trajectoryB.mat"
    trajectory = np.asarray(loadmat(trajectory_path, squeeze_me=True)["E21_trajectory"],
                            dtype=float)
    if (trajectory.ndim != 2 or trajectory.shape[1] < 4 or
            np.any(np.diff(trajectory[:, 0]) <= 0.0)):
        raise ValueError("invalid E21 trajectory table")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    products = [export_e7(raw_inms, flybys, args.output_dir)]
    products.extend(export_e21_species(swri, species, trajectory, args.output_dir)
                    for species in ("H2O", "CO2"))
    manifest = {
        "schema_version": 1,
        "contract": (
            "E7 is a calibrated local-density comparison; E21 products retain "
            "native counts and must not be scored against AMPS density"
        ),
        "products": products,
    }
    manifest_path = args.output_dir / "release_manifest.json"
    write_json(manifest_path, manifest)
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
