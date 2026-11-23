#!/usr/bin/env python3
"""Negative controls for the material campaign reader, not native evidence."""
import contextlib
import io
import json
from pathlib import Path
import struct
import sys
import tempfile
import unittest

sys.dont_write_bytecode = True
import run_material_checkpoint_tests as gate


class MaterialEvidenceControls(unittest.TestCase):
    def test_signed_energy_and_birth_are_independent(self):
        # Manufactured exact open-system budget: four represented particles,
        # birth energy 10 J, 2 J background heating, retained event energy 12 J.
        values = dict(epoch_s=12300, born=4, born_kinetic_j=10,
                      retained=4, returned=0, domain_exit=0, geometry_retired=0,
                      retained_birth_j=10, returned_birth_j=0, domain_birth_j=0, retired_birth_j=0,
                      retained_kinetic_j=12, returned_event_j=0, domain_event_j=0, retired_event_j=0,
                      background_work_j=2, wave_work_j=0)
        def text():
            return "REDUCED_SOURCE_LEDGER " + " ".join(f"{key}={value}" for key, value in values.items())
        self.assertEqual(gate.ledger(text(), [1])[1], 0)
        values["retained_birth_j"] = 9
        with self.assertRaisesRegex(AssertionError, "budget"):
            gate.ledger(text(), [1])
        values["retained_birth_j"] = 10
        values["background_work_j"] = -2
        with self.assertRaisesRegex(AssertionError, "budget"):
            gate.ledger(text(), [1])

    def test_no_particle_history_and_duplicate_ids_reject(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory/"material-particles-rank0.bin").write_bytes(struct.pack("<Q", 4))
            with self.assertRaisesRegex(ValueError, "empty population"):
                gate.particles(directory, 1)
            record = struct.pack("<QQ", 4, 123) + b"abcd"
            for rank in range(2):
                (directory/f"material-particles-rank{rank}.bin").write_bytes(record)
            with self.assertRaisesRegex(ValueError, "duplicate"):
                gate.particles(directory, 2)

    def test_global_denominator_and_unmeasured_negative_controls(self):
        with tempfile.TemporaryDirectory() as tmp:
            reference, actual = Path(tmp)/"reference.dat", Path(tmp)/"actual.dat"
            reference.write_text("TITLE=verification\n1.000000e+00 4.000000e+00\n")
            actual.write_bytes(reference.read_bytes())
            metadata = "0 -1 12300 1 1 300 300 1\n"
            for path in (reference, actual):
                path.with_suffix(".dat.sampling").write_text(metadata)
            self.assertEqual(gate.compare_output(reference, actual, True), 0)
            # A new segment's true one-sample average divided by the old
            # two-sample global window would be half its physical density.
            actual.write_text("TITLE=verification\n1.000000e+00 2.000000e+00\n")
            with self.assertRaisesRegex(AssertionError, "text precision"):
                gate.compare_output(reference, actual, True)
            actual.write_bytes(reference.read_bytes())
            actual.with_suffix(".dat.sampling").write_text("0 -1 12300 1 0 0 300 0\n")
            with self.assertRaisesRegex(AssertionError, "exposure"):
                gate.compare_output(reference, actual, True)

    def test_missing_receipt_and_failed_process_are_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            job = dict(id="MCP07", role="reject", ranks=1, cwd=str(directory), argv=["test-only"])
            plan = dict(hashes={}, jobs=[job])
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(gate.collect(plan, directory), 1)
            self.assertEqual(json.loads((directory/"summary.json").read_text())["results"][0]["status"], "ERROR")
            # Some MPI launchers can shut down gracefully after a timeout
            # signal. A zero exit code must not erase the deadline violation.
            (directory/"process-result.json").write_text(json.dumps(
                dict(argv=job["argv"], exit_code=0, launch_error="declared deadline exceeded")))
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(gate.collect(plan, directory), 1)
            self.assertEqual(json.loads((directory/"summary.json").read_text())["results"][0]["status"], "ERROR")
            (directory/"native.rank0.log").write_text("[MCP07] PASS\n")
            (directory/"process-result.json").write_text(json.dumps(dict(argv=job["argv"], exit_code=2)))
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(gate.collect(plan, directory), 1)
            self.assertEqual(json.loads((directory/"summary.json").read_text())["results"][0]["status"], "ERROR")

    def test_changed_source_identity_is_not_accepted(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            source = directory/"physics.cpp"
            source.write_text("original physical model")
            plan = dict(hashes={str(source): gate.digest(source)}, jobs=[])
            source.write_text("changed physical model")
            with self.assertRaisesRegex(ValueError, "identity changed"):
                gate.collect(plan, directory)

    def test_late_failure_must_stage_particles_and_retry(self):
        # These are collector negative controls, not manufactured MPI evidence.
        # The real C++ case independently compares full live images and links.
        with tempfile.TemporaryDirectory() as tmp:
            directory=Path(tmp)
            job=dict(id='MCP12',role='late-reject',ranks=8,cwd=str(directory),argv=['unit-only'])
            plan=dict(hashes={},jobs=[job])
            (directory/'process-result.json').write_text(json.dumps(dict(argv=job['argv'],exit_code=0)))
            for prepared,failed,retry in [(0,7,4),(32,0,4),(32,7,0)]:
                (directory/'native.rank0.log').write_text(
                    f'MATERIAL_LATE_REJECTION prepared_particles={prepared} failed_rank={failed} retry_step={retry}\n[MCP12] PASS\n')
                with contextlib.redirect_stdout(io.StringIO()):self.assertEqual(gate.collect(plan,directory),1)
                self.assertEqual(json.loads((directory/'summary.json').read_text())['results'][0]['status'],'FAIL')


if __name__ == "__main__":
    unittest.main()
