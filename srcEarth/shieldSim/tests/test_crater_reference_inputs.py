#!/usr/bin/env python3
"""Offline integrity and semantic checks for recovered CRaTER inputs.

This test is intentionally narrower than Stage-8 physical validation.  It
proves that the reviewed source package has not changed and that the two most
important machine-readable products retain the dimensions, units, dates, and
species coverage needed by a future instrument harness.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import pathlib
import zipfile


PROJECT = pathlib.Path(__file__).resolve().parents[1]
REFERENCE = PROJECT / "tests" / "reference_crater"


def sha256(path: pathlib.Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def check_manifest() -> None:
    """Verify every byte-level input, including the shared LET table."""
    manifest = REFERENCE / "manifest.sha256"
    entries = 0
    for line_number, raw in enumerate(manifest.read_text(encoding="utf-8").splitlines(), 1):
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        try:
            expected, relative = raw.split("  ", 1)
        except ValueError as error:
            raise AssertionError(f"malformed manifest line {line_number}") from error
        assert len(expected) == 64 and all(c in "0123456789abcdef" for c in expected)
        path = (REFERENCE / relative).resolve()
        # Manifest paths may reference the shared Stage-2 observation fixture,
        # but they may never escape the project source tree.
        try:
            path.relative_to(PROJECT)
        except ValueError as error:
            raise AssertionError(f"manifest path escapes project: {relative}") from error
        assert path.is_file(), f"missing reference input: {relative}"
        actual = sha256(path)
        assert actual == expected, f"SHA-256 mismatch for {relative}: {actual}"
        entries += 1
    assert entries == 40, f"unexpected manifest entry count: {entries}"


def check_bon_spectrum(path: pathlib.Path, expected_dates: str,
                       expected_phi: float) -> None:
    """Check the GLACE/BON arrays without assigning unpublished isotopes."""
    value = json.loads(path.read_text(encoding="utf-8"))
    assert value["BON_CODE_VERSION"] == "BON GCR_2020_Cpp 2.2,kw,ts,19Mar2021"
    assert value["RUN_DATES"] == expected_dates
    assert value["PHI_UNITS"] == "MV"
    assert value["ENERGY_UNITS"] == "MeV/nucleon"
    assert value["FLUX_UNITS"] == "ions / (MeV/nuc m^2 s sr)"
    assert math.isclose(float(value["RUN_PHI"]), expected_phi, rel_tol=0.0, abs_tol=1e-12)

    charges = value["GCR_Z"]
    energies = value["GCR_ENERGIES"]
    fluxes = value["GCR_FLUXES"]
    assert charges == list(range(1, 29))
    assert len(energies) == 100
    assert all(math.isfinite(float(x)) and float(x) > 0.0 for x in energies)
    assert all(float(b) > float(a) for a, b in zip(energies, energies[1:]))

    # The distributed JSON is energy-major: one 28-element row per energy.
    # Checking the actual serialization prevents a future adapter from silently
    # following the transposed shorthand used in the GLACE documentation.
    assert len(fluxes) == len(energies)
    assert all(len(row) == len(charges) for row in fluxes)
    assert all(
        math.isfinite(float(item)) and float(item) >= 0.0
        for row in fluxes for item in row
    )


def check_observed_let() -> None:
    path = (
        PROJECT / "tests" / "reference_spectra" / "observations"
        / "crater_gcr_let_spectra_v1_3.csv"
    )
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.reader(stream))
    assert len(rows) > 100
    assert rows[0][0] == "LET (keV / micron)"
    assert len(rows[0]) == 35  # one LET coordinate plus 34 time intervals
    assert rows[0][1] == "2009259-2010065"
    assert rows[0][-1] == "2026066-2026247"
    assert all(len(row) == len(rows[0]) for row in rows[1:])
    numeric = [[float(item) for item in row] for row in rows[1:]]
    assert all(math.isfinite(item) and item >= 0.0 for row in numeric for item in row)
    assert all(b[0] > a[0] for a, b in zip(numeric, numeric[1:]))


def check_pds_pipeline() -> None:
    """Prove that the pinned PDS ZIP is source, not an opaque executable."""
    with zipfile.ZipFile(REFERENCE / "pds" / "PIPE2.ZIP") as archive:
        names = set(archive.namelist())
        required = {
            "PIPE2/README.TXT",
            "PIPE2/pipe2.cpp",
            "PIPE2/crater/L1.hpp",
            "PIPE2/crater/L2.hpp",
            "PIPE2/conf/crat_hk_table_v02.txt",
        }
        assert required <= names
        level1 = archive.read("PIPE2/crater/L1.hpp").decode("utf-8")
        level2 = archive.read("PIPE2/crater/L2.hpp").decode("utf-8")
        assert "The MIT License" in level1
        assert "detCorOffset[j] + a *" in level1
        assert "thisLET = p1->energy[i][j]/detThickness[j]" in level2
        assert "p1->amplitude[i][j] > 3890" in level2


def check_document_boundaries() -> None:
    for relative in (
        "publications/CRaTER_SSR2010.pdf",
        "publications/Looper2020_CRaTER_GCR_LET.pdf",
        "publications/swe21079-sup-0001-2020SW002543-ts1.pdf",
        "zenodo/000README_GCR_LET_v1_3.pdf",
        "glace/000README_GLACE_v1_0.pdf",
        "pds/VOLSIS.PDF",
    ):
        assert (REFERENCE / relative).read_bytes().startswith(b"%PDF-")

    withheld = (REFERENCE / "pds" / "CRAT_CAL_PLAN.PDF").read_bytes()
    assert withheld.startswith(b"%PDF-")
    # The checksum pins the small PDF notice itself.  The companion text is
    # machine-readable evidence that the actual plan was removed; extracting
    # the notice sentence is documented in README.md but would add a poppler
    # dependency to this otherwise standard-library-only test.
    assert len(withheld) == 8265
    calibration_info = (REFERENCE / "pds" / "CALINFO.TXT").read_text()
    assert "Updated to reflect removal of CRaTER" in calibration_info

    label = (REFERENCE / "pds" / "CRAT_L2_PRI_2009259_V01.LBL").read_text()
    assert "FILE_RECORDS                 = 5287264" in label
    assert "RECORD_BYTES                 = 240" in label
    assert "START_TIME                   = 2009-09-16T00:00:00.26" in label


def check_acquisition_status() -> None:
    status = json.loads((REFERENCE / "acquisition.v1.json").read_text(encoding="utf-8"))
    assert status["schema_version"] == "shieldSim.crater_stage8_acquisition.v1"
    assert "not_stage8_validation_complete" in status["validation_status"]
    downloaded = {item["id"] for item in status["downloaded_inputs"]}
    assert {
        "observed_gcr_let_spectra_v1_3",
        "time_matched_bon2020_h_to_ni",
        "published_background_rejection_logic",
        "pds_level2_reconstruction_software",
    } <= downloaded
    unresolved = {item["id"] for item in status["unresolved_inputs_or_decisions"]}
    assert {
        "complete_crater_geant4_engineering_mass_model",
        "full_calibration_plan",
        "observational_bin_covariance",
        "calibration_and_blind_intervals",
        "acceptance_metric_and_uncertainty_model",
        "independent_stage8_review",
    } <= unresolved


def main() -> int:
    check_manifest()
    check_bon_spectrum(
        REFERENCE / "glace" / "solmin_gcr.json",
        "2009-09-16 to 2010-03-06 (2009/259 to 2010/065)",
        411.143,
    )
    check_bon_spectrum(
        REFERENCE / "glace" / "solmax_gcr.json",
        "2014-03-07 to 2014-09-04 (2014/066 to 2014/247)",
        816.77,
    )
    check_observed_let()
    check_pds_pipeline()
    check_document_boundaries()
    check_acquisition_status()
    print("CRaTER reference-input checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
