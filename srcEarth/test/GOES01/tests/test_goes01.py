#!/usr/bin/env python3
"""Small deterministic tests for calculations that do not require an AMPS run."""

from __future__ import annotations

import math
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from goes01_common import (  # noqa: E402
    atomic_json, power_law_bin_average, ratio_metrics, read_spectrum_zones,
    sha256_file, write_csv,
)
from plot_GOES01 import plot  # noqa: E402
from prepare_GOES01 import TEMPLATE, shift_boundary_to_midpoints  # noqa: E402


class SpectrumFoldingTests(unittest.TestCase):
    def test_constant_intensity_is_unchanged(self) -> None:
        self.assertAlmostEqual(power_law_bin_average((80.0, 90.0, 100.0),
                                                      (7.0, 7.0, 7.0), 83.0, 99.0), 7.0)

    def test_power_law_uses_analytic_bin_integral(self) -> None:
        energy = (80.0, 90.0, 100.0)
        intensity = tuple(value ** -2 for value in energy)
        expected = (1.0 / 83.0 - 1.0 / 99.0) / (99.0 - 83.0)
        self.assertAlmostEqual(power_law_bin_average(energy, intensity, 83.0, 99.0),
                               expected, places=15)

    def test_negative_native_intensity_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            power_law_bin_average((80.0, 100.0), (1.0, -1.0), 83.0, 99.0)

    def test_boundary_relabel_changes_only_time(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "source.dat"
            destination = Path(directory) / "midpoints.dat"
            source.write_text(
                "# source header\n"
                "2017-09-11T10:30:00Z 1.000000e+00 2.000000e+00\n"
                "2017-09-11T10:35:00Z 3.000000e+00 4.000000e+00\n",
                encoding="ascii")
            shift_boundary_to_midpoints(source, destination)
            data = [line for line in destination.read_text(encoding="ascii").splitlines()
                    if line and not line.startswith("#")]
            self.assertEqual(data, [
                "2017-09-11T10:32:30Z 1.000000e+00 2.000000e+00",
                "2017-09-11T10:37:30Z 3.000000e+00 4.000000e+00",
            ])


class MetricTests(unittest.TestCase):
    def test_physical_zero_stays_in_factor_two_denominator(self) -> None:
        metrics = ratio_metrics((
            {"model": 1.0, "reference": 1.0},
            {"model": 0.0, "reference": 1.0},
        ))
        self.assertEqual(metrics["factor_two_fraction"], 0.5)
        self.assertIsNone(metrics["median_abs_log10_ratio"])

    def test_factor_two_boundary_is_inclusive(self) -> None:
        metrics = ratio_metrics(({"model": 2.0, "reference": 1.0},))
        self.assertEqual(metrics["factor_two_fraction"], 1.0)
        self.assertAlmostEqual(metrics["median_abs_log10_ratio"], math.log10(2.0))


class NativeSpectrumReaderTests(unittest.TestCase):
    def test_snapshot_epoch_and_declared_rows_are_preserved(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "spectrum.dat"
            path.write_text(
                'TITLE="test"\n'
                'AUXDATA SNAPSHOT_EPOCH_UTC="2017-09-11T10:32:30Z"\n'
                'VARIABLES="E_MeV" "J_local_perMeV"\n'
                'ZONE T="point", I=2, F=POINT\n'
                '82.9 1.0\n99.8 2.0\n', encoding="utf-8")
            zones = read_spectrum_zones(path)
            self.assertEqual(zones[0].utc, "2017-09-11T10:32:30Z")
            self.assertEqual(zones[0].columns, ("E_MeV", "J_local_perMeV"))
            self.assertEqual(zones[0].rows.shape, (2, 2))

    def test_truncated_declared_zone_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "spectrum.dat"
            path.write_text(
                'AUXDATA SNAPSHOT_EPOCH_UTC="2017-09-11T10:32:30Z"\n'
                'VARIABLES="E_MeV" "J_local_perMeV"\n'
                'ZONE T="point", I=2, F=POINT\n82.9 1.0\n', encoding="utf-8")
            with self.assertRaises(ValueError):
                read_spectrum_zones(path)


class InputTemplateTests(unittest.TestCase):
    def test_every_section_is_supported_by_the_production_parser(self) -> None:
        """Prevent another typo/fictional section from reaching an expensive run."""
        supported = {
            "#RUN_INFO", "#CALCULATION_MODE", "#CUTOFF_RIGIDITY",
            "#PARTICLE_SPECIES", "#BACKGROUND_FIELD", "#ELECTRIC_FIELD",
            "#DOMAIN_BOUNDARY", "#OUTPUT_DOMAIN", "#NUMERICAL", "#MODE3D_MESH",
            "#DENSITY_3D", "#PARTICLE_TRAJECTORY", "#DENSITY_SPECTRUM",
            "#SPECTRUM", "#BOUNDARY_ANISOTROPY", "#OUTPUT_OPTIONS", "#TEMPORAL",
            "#ENERGY_CHANNELS", "#DETECTOR_RESPONSES",
        }
        sections = {line.strip() for line in TEMPLATE.read_text().splitlines()
                    if line.startswith("#") and line.strip() != "#END"}
        self.assertEqual(sections - supported, set())


class PlotTests(unittest.TestCase):
    def test_all_evidence_plots_retain_a_physical_zero(self) -> None:
        """Plotting must not crash or silently require every AMPS value to be positive."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            comparison = root / "comparison" / "p8a_intensity_comparison.csv"
            base = {
                "record_start_utc": "2017-09-11T10:30:00Z",
                "reference": 1.0, "reference_uncertainty": 0.1,
                "two_head_envelope_separation": 0.2,
            }
            write_csv(comparison, (
                {**base, "comparison_midpoint_utc": "2017-09-11T10:32:30Z",
                 "head": "WEST_MINUS_X", "model": 1.1, "model_observed_ratio": 1.1},
                {**base, "comparison_midpoint_utc": "2017-09-11T10:32:30Z",
                 "head": "EAST_PLUS_X", "model": 0.0, "model_observed_ratio": 0.0},
                {**base, "comparison_midpoint_utc": "2017-09-11T11:32:30Z",
                 "head": "WEST_MINUS_X", "model": 0.9, "model_observed_ratio": 0.9},
                {**base, "comparison_midpoint_utc": "2017-09-11T11:32:30Z",
                 "head": "EAST_PLUS_X", "model": 1.0, "model_observed_ratio": 1.0},
            ))
            ledger = root / "selection.csv"
            write_csv(ledger, ({
                "comparison_midpoint_utc": "2017-09-11T10:32:30Z",
                "two_head_envelope_separation": 0.2, "eligible": 1,
                "exclusion_reasons": "",
            },))
            (root / "GOES16").mkdir()
            atomic_json(root / "GOES16" / "prepared_run.json",
                        {"selection_ledger_file": str(ledger)})
            atomic_json(root / "comparison_summary.json", {
                "comparison_table": str(comparison),
                "comparison_table_sha256": sha256_file(comparison),
            })
            products = plot(root, ("png",))
            self.assertEqual(len(products), 4)
            self.assertTrue(all(path.stat().st_size > 0 for path in products))


if __name__ == "__main__":
    unittest.main()
