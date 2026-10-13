"""Synthetic forensic-contract checks; never count these as physical SA03 evidence."""
from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import subprocess
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from analyze_diagnostics import field_audits, repeatability_evidence, complete_direction_population, compare
from diagnose_SA03 import launch_command, override_directions_per_energy
from sa03_common import CHANNELS, sha256_file


class RepeatabilityTests(unittest.TestCase):
    """Ensure missing/changed fields cannot be hidden by successful trajectory closure."""

    def test_equal_products_without_field_evidence_are_incomplete(self):
        """An integrated difference of zero is not proof of a controlled experiment."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            epoch = "2017-09-06T10:00:00Z"
            for label in ("a", "b"):
                branch = root / label
                (branch / "input").mkdir(parents=True)
                deck = branch / "input/AMPS_PARAM.in"
                deck.write_text("synthetic immutable deck\n")
                (branch / "diagnostic_status.json").write_text(json.dumps({
                    "status": "COMPLETE", "return_code": 0, "complete_products": True,
                    "input_sha256": sha256_file(deck), "source_deck_sha256": "fixture",
                    "epochs": [epoch]}))
                (branch / "amps.log").write_text("no field audits in this legacy fixture\n")
                flux_columns = ["F_tot_m2s1"] + [f"F_{name}_m2s1" for name, _, _ in CHANNELS]
                for kind, columns in (("flux", flux_columns), ("density", ["N_cm^-3"])):
                    (branch / f"mode3d_points_{kind}_snapshot_0.dat").write_text(
                        'VARIABLES="TimeUTC" ' + " ".join(f'"{column}"' for column in columns) +
                        f'\nZONE I=1 F=POINT\n"{epoch}" ' + " ".join("1" for _ in columns) + "\n")
            result = compare(root / "a", root / "b", root / "comparison.json")
            self.assertEqual(result["integrated_difference_gate"], "PASS")
            self.assertEqual(result["status"], "INCOMPLETE")

    def test_production_unique_writer_plan(self):
        """Compile the production helper and reproduce the old alias/order failure."""
        with tempfile.TemporaryDirectory() as directory:
            executable = Path(directory) / "check-field-writers"
            source = Path(__file__).with_name("check_field_writer_plan.cpp")
            subprocess.run(["c++", "-std=c++14", "-O2", str(source), "-o", str(executable)], check=True)
            subprocess.run([str(executable)], check=True)

    def test_independent_concurrency_controls(self):
        """A serial field test must not silently retain the parallel initializer."""
        command = launch_command(Path("/example/amps"), 1, 1, "SERIAL", "SERIAL", "HC4")
        self.assertNotIn("-mode3d-parallel-field-init", command)
        self.assertEqual(command[command.index("-density-parallel") + 1], "SERIAL")
        command = launch_command(Path("/example/amps"), 2, 8, "OPENMP", "PARALLEL", "HC4")
        self.assertIn("-mode3d-parallel-field-init", command)
        # Binding is a recorded, explicit launcher choice, not a silent default
        # change that would confound comparisons with preserved older experiments.
        self.assertNotIn("--bind-to", command)
        unbound = launch_command(Path("/example/amps"), 1, 1, "SERIAL", "SERIAL", "HC4", "NONE")
        self.assertEqual(unbound[unbound.index("--bind-to") + 1], "none")

    def test_direction_override_changes_only_population_cap(self):
        """Angular refinement must preserve the declared energy and parent grids."""
        deck = """#DENSITY_SPECTRUM
  DS_DIRECTION_NMU 24
  DS_DIRECTION_NPHI 48
  DS_NINTERVALS 24
  DS_MAX_PARTICLES 2400
"""
        refined, evidence = override_directions_per_energy(deck, 384)
        self.assertEqual(evidence["energy_nodes"], 25)
        self.assertEqual(evidence["source_parent_direction_count"], 1152)
        self.assertEqual(evidence["requested_parent_direction_count"], 1152)
        self.assertEqual(evidence["source_directions_per_energy"], 96)
        self.assertEqual(evidence["requested_directions_per_energy"], 384)
        self.assertIn("DS_MAX_PARTICLES 9600", refined)
        self.assertIn("DS_DIRECTION_NMU 24", refined)
        self.assertIn("DS_DIRECTION_NPHI 48", refined)
        self.assertIn("DS_NINTERVALS 24", refined)

        expanded, evidence = override_directions_per_energy(deck, 1536, 32, 48)
        self.assertEqual(evidence["source_direction_nmu"], 24)
        self.assertEqual(evidence["requested_direction_nmu"], 32)
        self.assertEqual(evidence["requested_parent_direction_count"], 1536)
        self.assertIn("DS_DIRECTION_NMU 32", expanded)
        self.assertIn("DS_DIRECTION_NPHI 48", expanded)
        self.assertIn("DS_MAX_PARTICLES 38400", expanded)

    def test_direction_override_rejects_ambiguous_or_oversized_population(self):
        """Do not manufacture refinement evidence from a malformed source deck."""
        malformed = """DS_DIRECTION_NMU 2
DS_DIRECTION_NPHI 4
DS_NINTERVALS 2
DS_MAX_PARTICLES 10
"""
        with self.assertRaisesRegex(ValueError, "not an integral direction population"):
            override_directions_per_energy(malformed, 4)
        valid = malformed.replace("DS_MAX_PARTICLES 10", "DS_MAX_PARTICLES 12")
        with self.assertRaisesRegex(ValueError, "exceeds"):
            override_directions_per_energy(valid, 9)
        with self.assertRaisesRegex(ValueError, "requires two positive dimensions"):
            override_directions_per_energy(valid, 4, 8, 0)

    def test_direction_population_is_independent_of_log_agreement(self):
        """Two equally truncated logs must not pass; allow only energy roundoff."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            epoch = "2017-09-06T10:00:00Z"
            (root / "mode3d_points_spectrum_snapshot_0.dat").write_text(
                f'AUXDATA SNAPSHOT_EPOCH_UTC="{epoch}"\n'
                'VARIABLES="E_MeV" "N_sampled"\nZONE T="loc_0" I=1 F=POINT\n18.5 2\n')
            directions = {str(i): {"epoch": epoch, "particle_energy_MeV": "18.50000000000001",
                                   "direction": str(i)} for i in range(2)}
            self.assertTrue(complete_direction_population(root, directions))
            self.assertFalse(complete_direction_population(root, {"0": directions["0"]}))
            directions["1"]["particle_energy_MeV"] = "18.51"
            self.assertFalse(complete_direction_population(root, directions))

    def test_audit_requires_complete_exact_evidence(self):
        """Missing audits are ERROR; any changed array or sampled value is FAIL."""
        epoch = "2017-09-06T10:00:00Z"
        record = {epoch: {"content": "standalone-field-v1-0123", "meshRevision": "mesh-geometry-v1-4567",
                          "owner": {"sampled": 8, "differing": 0},
                          "interpolation": {"sampled": 12, "differing": 0}}}
        self.assertEqual(repeatability_evidence(record, record, [epoch])["status"], "PASS")
        for section, key, value in ((None, "content", "standalone-field-v1-9876"),
                                    ("owner", "differing", 1), ("interpolation", "differing", 1)):
            changed = copy.deepcopy(record)
            target = changed[epoch] if section is None else changed[epoch][section]
            target[key] = value
            self.assertEqual(repeatability_evidence(record, changed, [epoch])["status"], "FAIL")
        with self.assertRaisesRegex(ValueError, "incomplete"):
            repeatability_evidence(record, {}, [epoch])

    def test_real_log_format_and_missing_evidence(self):
        """Parse the production log grammar, including optional UTC trailing Z."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "amps.log"
            path.write_text(
                "[Mode3D] Prepared compact global cell-centered fields: snapshot=id, "
                "content=standalone-field-v1-abcd, meshRevision=mesh-geometry-v1-ef01, "
                "epoch=2017-09-06T10:00:00Z, frame=GSM\n"
                "[field-audit] owner epoch=2017-09-06T10:00:00 sampled=8 differing=0 maximum_relative_error=0\n"
                "[field-audit] interpolation epoch=2017-09-06T10:00:00Z sampled=12 differing=0\n")
            records = field_audits(root)
            self.assertEqual(repeatability_evidence(records, records, list(records))["status"], "PASS")
            path.write_text("ordinary run without forensic instrumentation\n")
            self.assertEqual(field_audits(root), {})


if __name__ == "__main__":
    unittest.main()
