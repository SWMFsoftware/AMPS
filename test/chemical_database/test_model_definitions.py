"""Actual shared model-metadata API controls; no observational test fixtures.

The authority graphs contain administrative labels only. Mathematical kernels
are independently checked in test_reduced_models.py, not by copying these
descriptions into an expected physical result.
"""
import csv
from decimal import Decimal
import json
import os
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "utility/chemical_database/python"))
from amps_chemical_database import Database, DatabaseError
from amps_chemical_database import model_definitions as models
import test_protocol
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

    def test_scattering_review_and_sensitivity_model_have_distinct_evidence(self):
        self.replace(models.definitions())
        database = self.open()
        review = database.model("proton-sodium-ground-capture-evidence-v1")
        isotropic = database.model("isotropic-com-two-body-scattering-v1")
        self.assertEqual(review["evidence_kind"], "published source-scope review")
        self.assertEqual(isotropic["evidence_kind"], "explicit experimental sensitivity model")
        self.assertEqual(review["equations"], ())
        self.assertIn("dP/dOmega = 1/(4*pi) [sr^-1]", isotropic["equations"])
        for entry in (review, isotropic):
            self.assertEqual(entry["implementation"], "descriptor-only")
            self.assertFalse(entry["physical_activation"])

    def test_reference_cases_require_their_reviewed_sources(self):
        # Protocol-only control. No physical data are manufactured here. A
        # familiar packet name or an unreviewed edition is insufficient to
        # publish a source-bound model with unresolved numerical inputs.
        empty = {"id": "poppe2022-lunar-pickup-rates", "source_packet": {}, "files": []}
        other_edition = {"id": "sodium-ground-capture-benchmark-20261124", "source_packet": {},
                         "files": [{"id": "shingal1986-na:ground-partial-benchmark", "sha256": "0" * 64}]}
        dataset, values = models.materialize(self.root / "publication", "control", ROOT, [empty, other_edition])
        self.assertFalse(any("source_bindings" in value for value in values))
        self.assertNotIn(models.BINDINGS_TABLE_ID, [entry.get("table", {}).get("id") for entry in dataset["files"]])


class SodiumReferencePublication(unittest.TestCase):
    requires_authoritative_source = True
    pickup_id = "poppe2022-na-pickup-median-v1"
    collision_id = "hplus-na3s-ground-capture-4kev-u-isotropic-v1"

    @classmethod
    def setUpClass(cls):
        if test_protocol.CONTEXT is None:
            raise unittest.SkipTest("shared database test context required")
        cls.db = test_protocol.CONTEXT.python_database()
        cls.available = {model["id"] for model in cls.db.models()}
        cls.require_both = os.environ.get("AMPS_CHECK_NA_REFERENCE_MODELS") == "1"

    def require_case(self, identifier):
        if identifier not in self.available and not self.require_both:
            self.skipTest("selected release omits " + identifier)

    def test_pickup_binding_reproduces_published_medians_and_birth_frame(self):
        self.require_case(self.pickup_id)
        data = self.db.model(self.pickup_id)["source_bindings"]
        row, = (row for row in self.db.table(data["table_id"]).dictionaries()
                if all(row[key] == value for key, value in data["row_match"].items()))
        # The raw ASCII source is independent of the projection and its binding.
        # Table 1 in the paper independently checks the roles and rounded units.
        deposited = self.db.source_path("poppe2022:ionization-rates").read_text().splitlines()[15].split()
        self.assertEqual(deposited[0], "Na")
        for process, index, published in (("photoionization", 2, "7.59E-6"),
                                         ("solar_wind_proton_charge_exchange", 5, "8.40E-7"),
                                         ("solar_wind_electron_impact", 8, "5.00E-7")):
            token = row[data["columns"][process]]
            self.assertEqual(token, deposited[index])
            self.assertEqual(Decimal(format(Decimal(token), ".3g")), Decimal(published))
        self.assertEqual((data["birth_velocity_m_s"], data["velocity_frame"]), ((0, 0, 0), "lunar"))
        self.assertEqual((data["unit"], data["heliocentric_distance_au"]), ("s^-1", 1))

    def test_ground_capture_binding_selects_theory_point_and_recoil(self):
        self.require_case(self.collision_id)
        data = self.db.model(self.collision_id)["source_bindings"]
        point, = (item for item in self.db.describe(data["dataset_id"])["files"] if item["id"] == data["source_id"])
        self.assertEqual(point["policy"], "metadata-only")
        for key in data["metadata_path"]:
            point = point[key]
        # Use the native area (Å²) in the independent unit conversion. The
        # table is indexed by projectile energy per atomic mass unit.
        self.assertEqual(Decimal(point[data["area_field"]]), Decimal("0.22") * Decimal("1e-16"))
        self.assertEqual((data["area_unit"], point["energy_unit"], data["projectile_energy_kev_u"]),
                         ("cm^2", "keV/u", point["projectile_energy"]))
        self.assertEqual((point["method"], data["energy_frame"], data["off_energy_policy"]),
                         ("A34 theory", "target-rest laboratory", "undefined"))
        self.assertIsNone(point["uncertainty"])
        self.assertEqual(data["product_states"], ("H(1s)", "Na+(2p6)"))
        self.assertIn("dP/dOmega = 1/(4*pi) [sr^-1]", self.db.model(data["angular_model_id"])["equations"])
        self.assertTrue(self.db.model(data["recoil_model_id"])["implementation"].endswith("::TwoBodyProducts"))

    def test_scenario_selectors_agree_in_shared_python_and_cpp_readers(self):
        from amps_chemical_database.api import plain
        if not self.require_both and self.available.isdisjoint((self.pickup_id, self.collision_id)):
            self.skipTest("selected release has no sodium scenario bindings")
        identity = dict(self.db.identity)
        fixture = SimpleNamespace(root=self.db.root, fields=identity, export_path=identity["export_path"])
        result = test_protocol.CONTEXT.native(fixture, "table", models.BINDINGS_TABLE_ID, "table-access")
        self.assertTrue(result["query_ok"], result)
        native = {model_id: json.loads(bindings) for model_id, bindings in result["rows"]}
        for model_id in (self.pickup_id, self.collision_id):
            if not self.require_both and model_id not in self.available:
                continue
            model = self.db.model(model_id)
            self.assertEqual(native[model_id], plain(model["source_bindings"]))
            self.assertFalse(model["physical_activation"])

    def test_each_source_packet_alone_selects_its_own_case(self):
        for packet_id, model_id in (("poppe2022-lunar-pickup-rates", self.pickup_id),
                                   ("sodium-ground-capture-benchmark-20261124", self.collision_id)):
            if not self.require_both and model_id not in self.available:
                continue
            with self.subTest(model=model_id):
                selections = models.sodium_scenarios([self.db.describe(packet_id)])
                self.assertEqual([entry["id"] for entry in selections], [model_id])
        if not self.require_both and self.available.isdisjoint((self.pickup_id, self.collision_id)):
            self.skipTest("selected release has no sodium scenario bindings")


class ScatteringPublication(unittest.TestCase):
    requires_authoritative_source = True

    @classmethod
    def setUpClass(cls):
        if os.environ.get("AMPS_CHECK_CEX_SCATTERING_RELEASE") != "1" or test_protocol.CONTEXT is None:
            raise unittest.SkipTest("explicit scattering-catalogue release test not selected")
        cls.database = test_protocol.CONTEXT.python_database()

    def test_native_integral_table_carries_review_without_value_edits(self):
        table = self.database.table("collisiondb_119990_embedded_arrays.csv")
        self.assertEqual(table.row_count, 12)
        self.assertEqual(table.rows[0], ("4000.0", "6.34e-15"))
        self.assertEqual(table.metadata["sha256"],
                         "2a887ad43497a7f9e2b3d7b6173ae1b2597cc893d897926762c181946c92db93")
        metadata = self.database.describe(table.metadata["table_id"])
        review = self.database.model(metadata["source_scope_review_id"])
        self.assertEqual(review["id"], "proton-sodium-ground-capture-evidence-v1")
        self.assertIn("4 keV/u angular data are unresolved", metadata["notes"])
        self.assertTrue(any("6.34e-15" in item and "Figure 5.6" in item for item in review["assumptions"]))
        with self.assertRaises(test_protocol.CONTEXT.python_api().DatabaseError) as refused:
            self.database.source_path("chemistry-archive-20261124:knoop-h-na")
        self.assertEqual(refused.exception.code, "PermissionRequired")

    def test_scattering_definitions_read_through_native_catalogue(self):
        table_id = "amps-reduced-model-definitions-v1"
        table = self.database.table(table_id)
        identity = dict(self.database.identity)
        fixture = SimpleNamespace(root=self.database.root, fields=identity, export_path=identity["export_path"])
        native = test_protocol.CONTEXT.native(fixture, "table", table_id, "table-access")
        self.assertTrue(native["query_ok"], native)
        self.assertEqual(native["rows"], [list(row) for row in table.rows])
        by_id = {row["model_id"]: row for row in table.dictionaries()}
        for key in ("proton-sodium-ground-capture-evidence-v1", "isotropic-com-two-body-scattering-v1"):
            with self.subTest(model=key):
                entry = self.database.model(key)
                self.assertFalse(entry["physical_activation"])
                self.assertEqual(by_id[key]["implementation"], "descriptor-only")
                for field in ("equations", "inputs", "outputs", "assumptions", "limitations", "citations"):
                    self.assertEqual(json.loads(by_id[key][field + "_json"]), list(entry[field]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
