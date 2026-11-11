#!/usr/bin/env python3
"""Stage 0 tests for the read-only DIM1 minor-species analyzer."""

from __future__ import annotations

import json
import os
import sys
import tempfile
import unittest
from pathlib import Path


SOURCE_ROOT = Path(__file__).resolve().parents[2]
VISUALIZATION = SOURCE_ROOT / "srcHartley2" / "visualization"
sys.path.insert(0, str(VISUALIZATION))

import analyze_minor_species as analyzer  # noqa: E402


VARIABLES = (
    'VARIABLES="X", "Number Density[1/m^3]", "Particle Number", '
    '"Local Time Step", "log10(Local Particle Weight)", "Cell Measure [m^3]"\n'
)


def write_table(path: Path, rows, include_measure: bool = True) -> None:
    variables = VARIABLES
    if not include_measure:
        variables = variables.replace(', "Cell Measure [m^3]"', "")
    with path.open("w", encoding="utf-8") as stream:
        stream.write(variables)
        for row in rows:
            selected = row if include_measure else row[:-1]
            stream.write(" ".join(str(value) for value in selected) + "\n")


class MinorSpeciesAnalyzerTest(unittest.TestCase):
    def test_exact_measures_and_comparison_contract(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            plots = root / "plots"
            plots.mkdir()
            previous = [
                (1000.0, 5.0, 0.0, 1.0, 10.0, 1.0),
                (2000.0, 25.0, 1.0, 2.0, 20.0, 2.0),
                (3000.0, 0.0, 3.0, 4.0, 40.0, 3.0),
            ]
            current = [
                (1000.0, 10.0, 0.0, 1.0, 10.0, 1.0),
                (2000.0, 20.0, 2.0, 2.0, 20.0, 2.0),
                (3000.0, 0.0, 4.0, 4.0, 40.0, 3.0),
            ]
            write_table(plots / "pic.C2.s=4.out=0.dat", previous)
            write_table(plots / "pic.C2.s=4.out=1.dat", current)
            comparison = root / "comparison.csv"
            comparison.write_text(
                "radius_km,AMPS_C2_density_cm-3,Rubin2011_C2_density_cm-3\n"
                "1,10,\n"
                "2,20,10\n",
                encoding="utf-8",
            )
            run_log = root / "run.log"
            run_log.write_text(
                "AMPS: radial convergence window=2 iteration=20 "
                "max-relative-L1=2.727273e-01 limiting-species=C2 "
                "consecutive=0/3 tolerance=5.000000e-02\n",
                encoding="utf-8",
            )
            prefix = root / "report"
            status = analyzer.main(
                [
                    str(plots), "--species", "C2", "--output-index", "1",
                    "--sampling-length", "10", "--radii-km", "1", "2",
                    "--run-log", str(run_log), "--comparison-csv", str(comparison),
                    "--output-prefix", str(prefix),
                ]
            )
            self.assertEqual(status, 0)
            report = json.loads(prefix.with_suffix(".json").read_text(encoding="utf-8"))
            species = report["species"][0]
            self.assertEqual(species["occupancy"]["zero_cell_count"], 1)
            self.assertAlmostEqual(species["occupancy"]["mean_all_cells"], 2.0)
            # Difference norm=|10-5|*1+|20-25|*2=15; reference norm=max(50,55).
            self.assertAlmostEqual(
                species["integrated_relative_l1_change"]["value"], 15.0 / 55.0
            )
            self.assertEqual(
                species["integrated_relative_l1_change"]["measure_source"],
                "amps_output_column",
            )
            self.assertIsNone(species["requested_radii"][0]["amps_to_reference_ratio"])
            radius = species["requested_radii"][1]
            self.assertAlmostEqual(radius["relative_window_change"], 0.2)
            self.assertEqual(radius["amps_to_reference_ratio"], 2.0)
            self.assertEqual(report["convergence"]["status"], "nonconverged")

    def test_legacy_output_does_not_invent_cell_measures(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            plots = Path(temporary)
            rows = [
                (1000.0, 1.0, 1.0, 1.0, 10.0, 1.0),
                (2000.0, 2.0, 1.0, 1.0, 10.0, 2.0),
            ]
            write_table(plots / "pic.C2.s=4.out=0.dat", rows, include_measure=False)
            write_table(plots / "pic.C2.s=4.out=1.dat", rows, include_measure=False)
            prefix = plots / "legacy"
            self.assertEqual(
                analyzer.main(
                    [str(plots), "--output-index", "1", "--output-prefix", str(prefix)]
                ),
                0,
            )
            report = json.loads(prefix.with_suffix(".json").read_text(encoding="utf-8"))
            integrated = report["species"][0]["integrated_relative_l1_change"]
            self.assertIsNone(integrated["value"])
            self.assertIn("exact AMPS cell measures", integrated["status"])


class ArchivedBaselineTest(unittest.TestCase):
    def test_preserved_c2_final_window(self) -> None:
        # The complete science archive remains outside Git because it contains
        # the executable and every species/output table.  The portable
        # orchestrator instead supplies a compact, checksummed manifest under
        # tests/data.  Direct invocations retain the complete-archive manifest
        # as their default so this test cannot silently redefine science
        # provenance when a user points it at the historical archive.
        baseline_override = os.environ.get("HARTLEY2_STAGE0_BASELINE_MANIFEST")
        baseline_path = (
            Path(baseline_override)
            if baseline_override
            else SOURCE_ROOT
            / "srcHartley2"
            / "halley2011"
            / "reference"
            / "halley1d_stage0_baseline.json"
        )
        baseline = json.loads(baseline_path.read_text(encoding="utf-8"))
        archive = os.environ.get("HARTLEY2_STAGE0_BASELINE_ROOT")
        if archive is None:
            candidate = SOURCE_ROOT / baseline["archive_relative_to_repository"]
            if candidate.is_dir():
                archive = str(candidate)
        if archive is None:
            self.skipTest(
                "set HARTLEY2_STAGE0_BASELINE_ROOT to verify the external validation archive"
            )
        archive_path = Path(archive)

        for relative, expected_digest in baseline["artifacts"].items():
            artifact = archive_path / relative
            self.assertTrue(artifact.is_file(), f"missing baseline artifact: {artifact}")
            self.assertEqual(analyzer.file_sha256(artifact), expected_digest)

        production = archive_path / "production"
        prefix = Path(tempfile.mkdtemp(prefix="hartley2-stage0-baseline-")) / "report"
        try:
            status = analyzer.main(
                [
                    str(production / "PT" / "plots"), "--species", "C2",
                    "--output-index", str(baseline["final_output_index"]),
                    "--sampling-length", str(baseline["sampling_length"]),
                    "--run-log", str(production / "run.log"),
                    "--comparison-csv",
                    str(production / "PT" / "plots" / "rubin2011_comparison" / "comparison_profiles.csv"),
                    "--output-prefix", str(prefix),
                ]
            )
            self.assertEqual(status, 0)
            report = json.loads(prefix.with_suffix(".json").read_text(encoding="utf-8"))
        finally:
            # TemporaryDirectory cannot be used because the analyzer writes
            # two sibling files; explicit unlink keeps failure artifacts out
            # of the repository while leaving no hidden test dependency.
            for suffix in (".json", ".csv"):
                path = prefix.with_suffix(suffix)
                if path.exists():
                    path.unlink()
            prefix.parent.rmdir()

        self.assertEqual(report["output_index"], baseline["final_output_index"])
        self.assertEqual(report["convergence"]["status"], baseline["scientific_status"])
        actual = report["species"][0]
        expected = baseline["c2"]
        self.assertEqual(actual["row_count"], expected["row_count"])
        self.assertEqual(actual["occupancy"]["zero_cell_count"], expected["zero_cell_count"])
        for actual_value, expected_name in (
            (actual["occupancy"]["mean_all_cells"], "mean_occupancy_all_cells"),
            (actual["occupancy"]["quantiles_positive_cells"]["p50"], "median_occupancy_positive_cells"),
            (actual["occupancy"]["maximum"], "maximum_occupancy"),
            (actual["occupancy"]["zero_cell_fraction"], "zero_cell_fraction"),
            (actual["occupancy"]["first_nonzero_radius_m"], "first_nonzero_radius_m"),
            (actual["occupancy"]["last_nonzero_radius_m"], "last_nonzero_radius_m"),
            (actual["carrier"]["weight_over_timestep_s-1_range"][0], "weight_over_timestep_min_s-1"),
            (actual["carrier"]["weight_over_timestep_s-1_range"][1], "weight_over_timestep_max_s-1"),
        ):
            self.assertAlmostEqual(
                actual_value / expected[expected_name], 1.0, places=12, msg=expected_name
            )
        radius_10000 = next(
            row for row in actual["requested_radii"] if row["radius_km"] == 10000.0
        )
        self.assertAlmostEqual(
            radius_10000["density_m-3"] / expected["density_at_10000_km_m-3"],
            1.0,
            places=12,
        )
        self.assertAlmostEqual(
            radius_10000["amps_to_reference_ratio"]
            / expected["amps_to_rubin_at_10000_km"],
            1.0,
            places=9,
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
