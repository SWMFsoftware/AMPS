"""Administrative packet-integrity tests, never observations or physical QA.

The tiny labelled files below exercise parser/publication admission only. Their
URLs and HTTP fields are administrative fixtures, not claims of actual external
acquisition. Real downloaded source tests are recorded separately by the owner.
"""
from pathlib import Path
import copy
import hashlib
import importlib.util
import json
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("tested_source_packets", ROOT / "utility/chemical_database/source_packets.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class SourcePackets(unittest.TestCase):
    """Source permission, identity and provenance boundaries of the real importer."""

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.native = self.root / "native.txt"
        self.native.write_text("ADMINISTRATIVE FIXTURE ONLY\n")
        self.licence = self.root / "licence.txt"
        self.licence.write_text("ADMINISTRATIVE PERMISSION FIXTURE ONLY\n")
        self.manifest = self.root / "packet.json"
        self.policy = {"storage_permission": "allowed", "evaluation_permission": "allowed",
                       "redistribution": "administrative-test", "evidence": [], "conditions": "not scientific data"}
        self.policies = {"metadata-only": dict(self.policy, storage_permission="unknown", evaluation_permission="unknown"),
                         "exomol-cc-by-sa-4.0": self.policy,
                         "nasa-cea-apache-2.0": self.policy}
        self.row = {"id": "administrative:source", "logical_name": "native.txt", "path": "native.txt",
                    "sha256": digest(self.native), "source_uri": "https://www.exomol.com/ADMINISTRATIVE-TEST",
                    "retrieved_utc": "2026-11-19T00:00:00+00:00", "http_status": 200,
                    "role": "source-metadata", "policy": "exomol-cc-by-sa-4.0", "quantity_metadata": {},
                    "rights_evidence": [{"path": "licence.txt", "sha256": digest(self.licence),
                                         "source_uri": "https://www.exomol.com/data/licence/"}]}
        self.packet = {"schema": "amps-chemical-source-packet-v1", "id": "administrative-packet",
                       "title": "ADMINISTRATIVE FIXTURE", "provider": "ADMINISTRATIVE FIXTURE",
                       "version": "test-only", "description": "Parser exercise, no numerical observations", "sources": [self.row]}

    def load(self, value=None):
        self.manifest.write_text(json.dumps(self.packet if value is None else value))
        return module.load_packets([self.manifest], self.policies)

    def test_exact_native_identity_and_immutable_input(self):
        before = digest(self.native), digest(self.licence)
        record = self.load()[0]["files"][0]
        self.assertEqual(record["sha256"], before[0])
        self.assertEqual(record["packet_acquisition_records"][0]["retrieved_utc"], self.row["retrieved_utc"])
        self.assertFalse(record["physical_activation"])
        self.assertEqual(before, (digest(self.native), digest(self.licence)))

    def test_unknown_policy_does_not_grant_storage(self):
        self.row["policy"] = "metadata-only"
        self.assertEqual(self.load()[0]["files"][0]["storage_permission"], "unknown")

    def test_arbitrary_permission_policy_refused(self):
        self.row["policy"] = "all-downloads-are-approved"
        with self.assertRaisesRegex(ValueError, "unreviewed"):
            self.load()

    def test_source_hash_corruption_refused(self):
        self.native.write_text("CHANGED ADMINISTRATIVE BYTES\n")
        with self.assertRaisesRegex(ValueError, "corrupt"):
            self.load()

    def test_rights_evidence_corruption_refused(self):
        self.licence.write_text("CHANGED EVIDENCE\n")
        with self.assertRaisesRegex(ValueError, "corrupt"):
            self.load()

    def test_missing_rights_evidence_refused(self):
        self.row["rights_evidence"] = []
        with self.assertRaisesRegex(ValueError, "rights evidence"):
            self.load()

    def test_provider_domain_mismatch_refused(self):
        self.row["source_uri"] = "https://example.invalid/source"
        with self.assertRaisesRegex(ValueError, "provider"):
            self.load()

    def test_nasa_policy_cannot_clear_another_repository(self):
        self.row.update(policy="nasa-cea-apache-2.0", source_uri="https://github.com/another/project/data")
        with self.assertRaisesRegex(ValueError, "another repository"):
            self.load()

    def test_new_policies_remain_provider_and_evidence_bound(self):
        # These administrative strings exercise guards, not actual rights
        # evidence. Real source/article notices are independently packet-pinned.
        for policy, host in (("chianti-use-attribution", "sohoftp.nascom.nasa.gov"),
                             ("springer-cc-by-4.0", "static-content.springer.com")):
            with self.subTest(policy=policy):
                self.policies[policy] = self.policy
                self.row.update(policy=policy, source_uri="https://" + host + "/ADMINISTRATIVE-TEST")
                self.assertEqual(self.load()[0]["files"][0]["storage_permission"], "allowed")
                self.row["source_uri"] = "https://example.invalid/ADMINISTRATIVE-TEST"
                with self.assertRaisesRegex(ValueError, "provider"):
                    self.load()

    def test_ntrs_determination_cannot_clear_another_record(self):
        policy = "ntrs-19700032941-public-use"
        self.policies[policy] = self.policy
        self.row.update(policy=policy, source_uri="https://ntrs.nasa.gov/api/citations/19700032941/downloads/19700032941.pdf")
        self.assertEqual(self.load()[0]["files"][0]["storage_permission"], "allowed")
        self.row["source_uri"] = "https://ntrs.nasa.gov/api/citations/197000329410"
        with self.assertRaisesRegex(ValueError, "specific to record"):
            self.load()

    def test_insecure_uri_refused(self):
        self.row["source_uri"] = "http://www.exomol.com/source"
        with self.assertRaisesRegex(ValueError, "HTTPS"):
            self.load()

    def test_http_failure_is_not_a_source(self):
        self.row["http_status"] = 500
        with self.assertRaisesRegex(ValueError, "HTTP200"):
            self.load()

    def test_missing_explicit_timestamp_refused(self):
        self.row["retrieved_utc"] = "2026-11-19T00:00:00"
        with self.assertRaisesRegex(ValueError, "UTC"):
            self.load()

    def test_absolute_member_refused(self):
        self.row["path"] = str(self.native)
        with self.assertRaisesRegex(ValueError, "escapes"):
            self.load()

    def test_traversal_refused(self):
        self.row["path"] = "../native.txt"
        with self.assertRaisesRegex(ValueError, "escapes"):
            self.load()

    def test_external_symlink_refused(self):
        external = self.root.parent / (self.root.name + "-outside")
        external.write_text("ADMINISTRATIVE OUTSIDE OBJECT")
        self.addCleanup(external.unlink)
        (self.root / "escape").symlink_to(external)
        self.row.update(path="escape", sha256=digest(external))
        with self.assertRaisesRegex(ValueError, "escapes"):
            self.load()

    def test_duplicate_source_ids_refused(self):
        self.packet["sources"].append(copy.deepcopy(self.row))
        with self.assertRaisesRegex(ValueError, "duplicate"):
            self.load()

    def test_duplicate_packet_ids_refused(self):
        self.load()
        with self.assertRaisesRegex(ValueError, "duplicate"):
            module.load_packets([self.manifest, self.manifest], self.policies)

    def test_false_physical_activation_refused(self):
        for target in (self.packet, self.row):
            target["physical_activation"] = True
            with self.assertRaisesRegex(ValueError, "activate"):
                self.load()
            target.pop("physical_activation")

    def test_duplicate_json_keys_refused(self):
        self.manifest.write_text('{"schema":"a","schema":"b"}')
        with self.assertRaisesRegex(ValueError, "duplicate"):
            module.load_packets([self.manifest], self.policies)

    def test_nonfinite_json_refused(self):
        self.manifest.write_text('{"schema":NaN}')
        with self.assertRaisesRegex(ValueError, "nonfinite"):
            module.load_packets([self.manifest], self.policies)

    def test_missing_quantity_metadata_refused(self):
        self.row.pop("quantity_metadata")
        with self.assertRaisesRegex(ValueError, "quantity"):
            self.load()

    def derived(self):
        recipe = self.root / "recipe.py"
        recipe.write_text("# Administrative fixture; never scientific processing.\n")
        row = copy.deepcopy(self.row)
        row.update(id="administrative:projection", role="derived-table", derived_from=[self.row["id"]],
                   recipe={"path": "recipe.py", "sha256": digest(recipe)})
        self.packet["sources"].append(row)
        return row

    def test_derived_source_preserves_native_parent_policy(self):
        row = self.derived()
        result = self.load()[0]["files"][1]
        self.assertEqual(result["derived_from"], [self.row["id"]])
        self.assertEqual(result["storage_permission"], "allowed")

    def test_derived_source_cannot_clear_held_parent(self):
        self.derived()
        self.row["policy"] = "metadata-only"
        with self.assertRaisesRegex(ValueError, "cannot clear"):
            self.load()

    def test_derived_source_requires_actual_parent(self):
        self.derived()["derived_from"] = ["not-declared"]
        with self.assertRaisesRegex(ValueError, "identity is missing"):
            self.load()

    def test_derived_recipe_corruption_refused(self):
        self.derived()
        (self.root / "recipe.py").write_text("changed")
        with self.assertRaisesRegex(ValueError, "corrupt"):
            self.load()


if __name__ == "__main__":
    unittest.main(verbosity=2)
