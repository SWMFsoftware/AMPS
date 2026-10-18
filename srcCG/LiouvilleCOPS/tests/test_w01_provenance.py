"""Fast W01 provenance/cache tests; no mission data or CSPICE required."""

import hashlib
import json
from pathlib import Path
import struct
import sys
import tempfile
import unittest


APP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP / "tools"))
from cops_io import (MODEL_COLUMNS, ProductError, file_identity, load_manifest,
                     read_model_table, sha256_file, validate_fit_base,
                     write_json_atomic)


class ProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="cops_w01_")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_streaming_sha256_and_atomic_json(self):
        payload = self.root / "payload"
        payload.write_bytes(bytes(range(256)) * 4097)
        self.assertEqual(sha256_file(payload), hashlib.sha256(payload.read_bytes()).hexdigest())
        identity = file_identity(payload, "fixture")
        self.assertEqual(identity["size_bytes"], payload.stat().st_size)
        manifest = self.root / "manifest.json"
        value = {"schema": "fixture", "schema_version": 1, "identity": identity}
        write_json_atomic(manifest, value)
        self.assertEqual(load_manifest(manifest)["identity"], identity)
        self.assertFalse((self.root / "manifest.json.tmp").exists())

    def test_model_table_rejects_duplicate_points_and_width(self):
        table = self.root / "COPS_model.dat"
        table.write_text("#iPoint UTC pNG\n0 t0 1\n0 t1 2\n")
        with self.assertRaisesRegex(ProductError, "duplicate iPoint"):
            read_model_table(table)
        table.write_text("#iPoint UTC pNG\n0 t0\n")
        with self.assertRaisesRegex(ProductError, "fields"):
            read_model_table(table)

    def make_base(self):
        run = self.root / "run"
        (run / "kernels").mkdir(parents=True)
        values = {name: "0" for name in MODEL_COLUMNS}
        values.update({"iPoint": "0", "UTC": "2016-09-30T0:00:04", "et": "1",
                       "x": "3", "y": "4", "z": "5", "altitude": "10",
                       "rh_AU": "2", "VdotB_RG": "2", "pNG": "1e-8",
                       "pRG": "2e-8", "refNG": "1e-8", "refRG": "2e-8",
                       "T_RG": "293", "useNG": "1", "useRG": "1"})
        (run / "COPS_model.dat").write_text(
            "# fixture\n#" + " ".join(MODEL_COLUMNS) + "\n" +
            " ".join(values[name] for name in MODEL_COLUMNS) + "\n")
        element = [0.0] * 34
        element[6] = 1.0
        (run / "elements.bin").write_bytes(
            struct.pack("<iii", 1, 1, 25) + struct.pack("<34d", *element))
        header = struct.pack("<ii10di", 1, 0, 1.0, 2.0, 3.0, 4.0, 5.0,
                             0.0, 0.0, 1.0, 2.0 * 149597870700.0, 10.0, 1)
        record = struct.pack("<i4di", 0, 1.0, 2.0, 3.0, 0.5, 0)
        (run / "kernels/point=0000.bin").write_bytes(header + record)
        return run

    def test_complete_static_cache_validation(self):
        run = self.make_base()
        audit = validate_fit_base(run)
        self.assertTrue(audit["legacy_manifest"])
        self.assertEqual(audit["elements"]["nfaces"], 1)
        self.assertEqual(audit["kernels"]["count"], 1)
        self.assertEqual(audit["kernels"]["records"], 1)

    def test_cache_rejects_truncation_face_range_and_stale_kernel(self):
        run = self.make_base()
        kernel = run / "kernels/point=0000.bin"
        original = kernel.read_bytes()
        kernel.write_bytes(original[:-1])
        with self.assertRaisesRegex(ProductError, "size mismatch"):
            validate_fit_base(run)
        kernel.write_bytes(original[:92] + struct.pack("<i4di", 1, 1., 2., 3., .5, 0))
        with self.assertRaisesRegex(ProductError, "out-of-range face"):
            validate_fit_base(run)
        kernel.write_bytes(original)
        (run / "kernels/point=0001.bin").write_bytes(original)
        with self.assertRaisesRegex(ProductError, "unexpected"):
            validate_fit_base(run)


if __name__ == "__main__":
    unittest.main()
