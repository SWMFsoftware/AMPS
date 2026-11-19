"""Actual shared model-metadata API controls; no observational test fixtures.

The authority graphs contain administrative labels only. Mathematical kernels
are independently checked in test_reduced_models.py, not by copying these
descriptions into an expected physical result.
"""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "utility/chemical_database/python"))
from amps_chemical_database import Database, DatabaseError
from amps_chemical_database import model_definitions as models
from test_protocol import Fixture, digest


class ModelMetadata(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.fixture = Fixture(self.root / "administrative")

    def open(self):
        f = self.fixture.fields
        return Database.open(self.fixture.root, f["release_id"], f["release_manifest_sha256"],
            f["profile_id"], f["profile_sha256"], "export.in", digest(self.fixture.root / "export.in"))

    def replace(self, values):
        path = self.fixture.root / "catalogue.json"
        catalogue = json.loads(path.read_text())
        catalogue["models"] = values
        path.write_text(json.dumps(catalogue))
        self.fixture.publish()

    def test_historical_catalogue_has_no_implicit_models(self):
        database = self.open()
        self.assertEqual(database.models(), ())
        with self.assertRaises(DatabaseError) as context:
            database.model("polanyi-wigner-fixed-sites-v1")
        self.assertEqual(context.exception.code, "Unavailable")

    def test_explicit_definitions_are_immutable_and_describable(self):
        self.replace(models.definitions())
        database = self.open()
        value = database.model("opal-energy-marginal-v1")
        self.assertEqual(value, database.describe(value["id"]))
        self.assertFalse(value["physical_activation"])
        with self.assertRaises(TypeError):
            value["inputs"][0] = "changed"
        with self.assertRaises(TypeError):
            value["physical_activation"] = True

    def test_model_metadata_does_not_grant_runtime_capabilities(self):
        self.replace(models.definitions())
        database = self.open()
        for capability in ("scalar-rates", "complete-daughter-kinetics", "surface-neutralization"):
            with self.subTest(capability=capability), self.assertRaises(DatabaseError) as context:
                database.require(capability)
            self.assertEqual(context.exception.code, "Unavailable")

    def test_invalid_model_query_has_typed_input_failure(self):
        database = self.open()
        for identifier in (None, [], {}, "", "not an ID"):
            with self.subTest(identifier=identifier), self.assertRaises(DatabaseError) as context:
                database.model(identifier)
            self.assertEqual(context.exception.code, "InvalidInput")

    def test_duplicate_or_ambiguous_model_id_is_corruption(self):
        value = models.definitions()[0]
        for values in ((value, value), (dict(value, id="fixture-a"),), (dict(value, id="source-z"),),
                       (dict(value, id="table-z"),)):
            with self.subTest(values=values):
                self.replace(values)
                with self.assertRaises(DatabaseError) as context:
                    self.open()
                self.assertEqual(context.exception.code, "CorruptSource")

    def test_activation_or_malformed_identity_is_corruption(self):
        value = models.definitions()[0]
        for bad in (dict(value, physical_activation=True), dict(value, physical_activation=None),
                    dict(value, id=""), dict(value, id="bad model"), "not a model object"):
            with self.subTest(bad=bad):
                self.replace([bad])
                with self.assertRaises(DatabaseError) as context:
                    self.open()
                self.assertEqual(context.exception.code, "CorruptSource")

    def test_registry_returns_defensive_copies(self):
        before = models.definitions()
        changed = models.definitions()
        changed[0]["inputs"].append("ADMINISTRATIVE CHANGE")
        self.assertEqual(models.definitions(), before)

    def test_authored_csv_json_are_lossless_software_definitions(self):
        dataset, values = models.materialize(self.root / "publication", "administrative", ROOT)
        originals = {row["id"]: row for row in dataset["files"]}
        table = originals["amps-models:definitions-csv"]
        with Path(table["path"]).open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual([row["model_id"] for row in rows], [value["id"] for value in values])
        for row, value in zip(rows, values):
            for key in ("equations", "inputs", "outputs", "assumptions", "limitations", "citations"):
                self.assertEqual(json.loads(row[key + "_json"]), value[key])
        for original in originals.values():
            self.assertEqual(original["sha256"], digest(original["path"]))
            self.assertNotIn("http_status", original)
            self.assertNotIn("retrieved_utc", original)
        with self.assertRaises(FileExistsError):
            models.materialize(self.root / "publication", "administrative", ROOT)

    def test_missing_software_licence_does_not_publish_definitions(self):
        with self.assertRaisesRegex(ValueError, "licence"):
            models.materialize(self.root / "publication", "administrative", self.root)


if __name__ == "__main__":
    unittest.main(verbosity=2)
