#!/usr/bin/env python3
"""Check I17 plumbing without producing model or observational evidence.

Temporary values below test CSV syntax only. No numerical physics is mocked
and called verified: the launcher test deliberately stops at manifest creation,
before MPI, using MissingPrerequisite. These checks never create a plot.
"""
import csv
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

TEST_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(TEST_ROOT))
sys.path.insert(0, str(TEST_ROOT / "common"))
import integration_common as integration
import write_comparison_artifacts as plotting
from i17_convergence import compare_time_steps


class InfrastructureContract(unittest.TestCase):
    def test_case_sensitive_path_is_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "application.in"
            path.write_text("diviner_eclipse_history_file = /tmp/D03-RDR-I17/June.csv\n")
            self.assertEqual(integration.moon_setting(
                path, "diviner_eclipse_history_file", preserve_case=True),
                "/tmp/D03-RDR-I17/June.csv")

    def test_numerical_budget_and_decreasing_difference(self):
        # These are arithmetic gate fixtures, never observations or thermal
        # solutions. They independently distinguish the two required rules:
        # a small but increasing difference must fail, as must a decreasing
        # difference that still exceeds the approved uncertainty budget.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            directories = [root / str(i) for i in range(3)]
            for directory in directories:
                directory.mkdir()
            def write_values(values):
                for directory, value in zip(directories, values):
                    row = {"site_time_bin_id": "arithmetic-fixture",
                           "utc": "fixture", "ephemeris_time_tdb_s": "0",
                           "latitude_deg": "0", "longitude_east_deg": "0",
                           "observed_uncertainty_K": "1",
                           "observed_bolometric_brightness_temperature_K": "100",
                           "model_bolometric_brightness_temperature_K": str(value)}
                    with (directory / "receipt.txt.I17_training_model.csv").open("w", newline="") as stream:
                        writer = csv.DictWriter(stream, fieldnames=list(row))
                        writer.writeheader()
                        writer.writerow(row)
            for values, expected in (([100, 101, 101.05], True),
                                     ([100, 100.01, 100.05], False),
                                     ([100, 101, 101.2], False),
                                     ([100, 100, 100], True)):
                write_values(values)
                passed, metrics = compare_time_steps(directories, root / "gate.csv")
                self.assertEqual(passed, expected)
                self.assertEqual(metrics["compared_bins"], 1)

    def test_derived_input_reaches_manifest_and_mpi(self):
        # Exercise the actual orchestration through input derivation. Fake gate
        # records live only inside TemporaryDirectory and cannot qualify a real
        # campaign. The manifest hook stops before launching an executable.
        for test_id in ("I16", "I17"):
            with self.subTest(test_id=test_id), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                template = root / "template.in"
                template.write_text("template sentinel\n")
                (root / "I13").mkdir()
                (root / "I13/result.json").write_text('{"status":"PASS"}')
                captured = {}

                def prepare(source, destination, output_root):
                    destination.write_text("derived sentinel\n")

                def stop_at_manifest(tid, executable, command, seed, receipt, ranks):
                    captured.update(command=command[:], receipt=receipt, ranks=ranks)
                    raise integration.MissingPrerequisite("infrastructure-only stop")

                with patch.object(integration, "dataset_package_preflight",
                                  return_value=(True, "fixture", {})), \
                     patch.object(integration, "prepare_" + test_id.lower() + "_input",
                                  side_effect=prepare), \
                     patch.object(integration, "manifest", side_effect=stop_at_manifest), \
                     patch.object(integration, "run_monitored") as launch:
                    code = integration.run(test_id, root / "unused-executable",
                                           root, template, 8, "qualification")
                self.assertEqual(code, 77)
                launch.assert_not_called()
                expected = (root / test_id / "application.in").resolve()
                self.assertEqual(captured["receipt"], expected)
                command = captured["command"]
                self.assertEqual(Path(command[command.index("-input") + 1]), expected)
                self.assertEqual(command[command.index("-n") + 1], "8")

    def test_experimental_selection_is_accepted_but_unknown_is_rejected(self):
        # Parser fixture only: these numbers carry no scientific significance
        # and never enter render(), a production run, or validation artifacts.
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "parser_fixture.csv"
            row = dict(zip(plotting.REQUIRED_COLUMNS,
                           ("parser-fixture", "0", "s", "1", "1", "K",
                            "0.1", "training-experimental")))
            def write():
                with path.open("w", newline="") as stream:
                    writer = csv.DictWriter(stream, plotting.REQUIRED_COLUMNS)
                    writer.writeheader()
                    writer.writerow(row)
            write()
            self.assertEqual(plotting.read_comparison(path)[0]["selection"],
                             "training-experimental")
            row["selection"] = "unknown"
            write()
            with self.assertRaises(ValueError):
                plotting.read_comparison(path)


if __name__ == "__main__":
    unittest.main()
