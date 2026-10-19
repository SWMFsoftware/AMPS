"""Small, self-consistent LiouvilleCOPS base runs for unit tests.

``make_run`` writes the four products a base run consists of, in the exact
binary/text formats of the C++ application:

* ``COPS_model.dat``: header with every key ``fit_cops.Model`` parses, and
  rows whose geometry equals the kernel headers;
* ``elements.bin``: version 1, 25 basis values per facet;
* ``kernels/point=NNNN.bin``: version 1;
* ``run_manifest.json``: digests computed like ``cops_io``.

The values are random but valid. Each point sees a random subset of facets in
a random order, which exercises the record-order preservation. One point sees
nothing, which exercises empty CSR rows. The physics is meaningless; these
runs only test formats, identities, layout and the backend contract. Each
corruption helper damages one invariant, so a test can check that exactly
that invariant is detected.
"""

import hashlib
import json
from pathlib import Path

import numpy as np

AU = 149597870700.0
HEADER_LINES = (
    "LiouvilleCOPS model of the COPS gauge pressures; synthetic test run",
    "ProductKind=forward TableSchema=liouville_cops_forward_v1",
    "ShapeModel=synthetic FlightTimeRotation=none FitInterval=all",
    "AnalyticalTmax=182.09922 AnalyticalT75=136.11149 AnalyticalTmin=133 "
    "HeliocentricModel=shap5 ActivityFile=none",
    "SurfaceBC=flat TemperatureModel=analytical ActivityModel=shap5 SZAMode=continuous "
    "CastShadows=on SpacecraftVelocity=on",
    "NGAcceptance=srccg_shield NGHeadOnWeight=on RGAcceptance=hemisphere "
    "RGPressureModel=equilibrium_sphere TemperatureNG=293 TemperatureRGMode=hk_table "
    "SpeciesSensitivity=on (H2O 0.893, CO2 0.704) ReferenceDataUnit=mbar",
)
COLUMNS = ("iPoint", "UTC", "et", "x", "y", "z", "altitude", "inside", "rh_AU", "VdotB_RG",
           "nVisible", "Q_H2O", "Q_CO2", "n_H2O", "n_CO2", "PhiRG_H2O", "PhiRG_CO2",
           "pNG_H2O", "pNG_CO2", "pNG", "pRG_H2O", "pRG_CO2", "pRG", "refNG", "refRG",
           "NGBelowDeckFraction", "T_RG", "useNG", "useRG")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fibonacci_sphere(n, radius=2000.0):
    """Nearly uniform facet centroids so harmonic Gram matrices are well posed."""
    k = np.arange(n) + 0.5
    polar = np.arccos(1.0 - 2.0 * k / n)
    azimuth = np.pi * (1.0 + 5.0 ** 0.5) * k
    unit = np.column_stack((np.sin(polar) * np.cos(azimuth),
                            np.sin(polar) * np.sin(azimuth), np.cos(polar)))
    return radius * unit, unit


def write_elements(path, table):
    with Path(path).open("wb") as stream:
        stream.write(np.array([1, len(table), 25], dtype="<i4").tobytes())
        stream.write(np.ascontiguousarray(table, dtype="<f8").tobytes())


def write_kernel(path, version, point, geometry, records, count=None):
    header = np.zeros(1, dtype=[("version", "<i4"), ("point", "<i4"),
                                ("geometry", "<f8", (10,)), ("count", "<i4")])
    header["version"], header["point"] = version, point
    header["geometry"] = geometry
    header["count"] = len(records) if count is None else count
    with Path(path).open("wb") as stream:
        stream.write(header.tobytes())
        stream.write(records.tobytes())


RECORD = np.dtype([("face", "<i4"), ("Kn", "<f8"), ("KPhi", "<f8"), ("KV", "<f8"),
                   ("cos", "<f8"), ("shadow", "<i4")])


def make_run(directory, nfaces=240, npoints=6, seed=7, empty_point=2):
    """Write a synthetic base run; return a description dictionary."""
    rng = np.random.default_rng(seed)
    run = Path(directory)
    (run / "kernels").mkdir(parents=True)
    centroid, normal = fibonacci_sphere(nfaces)
    table = np.empty((nfaces, 34))
    table[:, 0:3] = centroid
    table[:, 3:6] = normal
    table[:, 6] = rng.uniform(1.0e3, 3.0e3, nfaces)
    table[:, 7:32] = rng.normal(0.0, 1.0, (nfaces, 25))
    table[:, 7] = 1.0                             # Y_00 keeps the prior positive
    table[:, 32] = rng.uniform(0.5, 1.5, nfaces)
    table[:, 33] = rng.uniform(0.0, 0.2, nfaces)
    write_elements(run / "elements.bin", table)

    point_ids = [10 * i + 3 for i in range(npoints)]   # ids are not row indices
    rows, kernels = [], {}
    for index, point in enumerate(point_ids):
        et = 528465672.0 + 60.0 * index
        sc = np.array([9000.0 - 50 * index, -6400.0 + 40 * index, 15000.0 - 30 * index])
        sun = np.array([0.6, 0.0, 0.8])
        rh = 3.831 * AU
        altitude = float(np.linalg.norm(sc) - 2000.0)
        vdotb = 0.28 + 0.01 * index
        geometry = np.array([et, vdotb, sc[0], sc[1], sc[2], sun[0], sun[1], sun[2], rh, altitude])
        count = 0 if index == empty_point else int(rng.integers(nfaces // 3, nfaces))
        records = np.zeros(count, dtype=RECORD)
        records["face"] = rng.permutation(nfaces)[:count]
        records["Kn"] = rng.uniform(0.0, 1e-6, count)
        records["KPhi"] = np.where(rng.random(count) < 0.2, 0.0, rng.uniform(0.0, 1e-6, count))
        records["KV"] = np.where(records["KPhi"] == 0.0, 0.0, rng.uniform(0.0, 1e-7, count))
        records["cos"] = rng.uniform(-1.0, 1.0, count)
        records["shadow"] = (rng.random(count) < 0.1).astype("<i4")
        path = run / "kernels" / ("point=%04d.bin" % point)
        write_kernel(path, 1, point, geometry, records)
        kernels[point] = (geometry, records)
        rows.append([str(point), "2016-09-30T00:%02d:04" % index, "%.6f" % et,
                     "%.6e" % sc[0], "%.6e" % sc[1], "%.6e" % sc[2], "%.6e" % altitude, "0",
                     "%.8f" % (rh / AU), "%.6e" % vdotb, str(count)] +
                    ["1.0e+20"] * 6 + ["2.0e-07", "2.0e-08", "2.2e-07",
                                       "6.0e-07", "6.0e-08", "6.6e-07",
                                       "%.6e" % (2.2e-7 * (1 + 0.05 * index)),
                                       "%.6e" % (6.6e-7 * (1 - 0.03 * index)),
                                       "0.0", "250.082", "1", "1"])
    with (run / "COPS_model.dat").open("w") as stream:
        for line in HEADER_LINES:
            stream.write("# %s\n" % line)
        stream.write("#" + " ".join(COLUMNS) + "\n")
        for row in rows:
            stream.write(" ".join(row) + "\n")
    write_manifest(run, point_ids)
    return {"run": run, "nfaces": nfaces, "point_ids": point_ids, "kernels": kernels,
            "elements": table}


def write_manifest(run, point_ids):
    """Run manifest with cops_io digest definitions (overwrites)."""
    run = Path(run)
    set_digest, records, size = hashlib.sha256(), 0, 0
    for point in point_ids:
        name = "point=%04d.bin" % point
        path = run / "kernels" / name
        set_digest.update(name.encode("ascii"))
        set_digest.update(sha256(path).encode("ascii"))
        records += (path.stat().st_size - 92) // 40
        size += path.stat().st_size
    manifest = {"schema": "liouville_cops.run_manifest", "schema_version": 1,
                "states": {"process": "complete", "cache": "verified"},
                "outputs": {
                    "model_table": {"path": str(run / "COPS_model.dat"),
                                    "size_bytes": (run / "COPS_model.dat").stat().st_size,
                                    "sha256": sha256(run / "COPS_model.dat")},
                    "elements": {"path": str(run / "elements.bin"),
                                 "size_bytes": (run / "elements.bin").stat().st_size,
                                 "sha256": sha256(run / "elements.bin")},
                    "kernels": {"version": 1, "count": len(point_ids), "records": records,
                                "size_bytes": size, "set_sha256": set_digest.hexdigest()}}}
    (run / "run_manifest.json").write_text(json.dumps(manifest, indent=1))


def kernel_path(run, point):
    return Path(run) / "kernels" / ("point=%04d.bin" % point)
