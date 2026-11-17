#!/usr/bin/env python3
"""Shared M2 contract dispatch tests; temporary metadata, no science values."""
import json
from pathlib import Path
import tempfile
import unittest

import m2_contracts as contracts


class ContractSelection(unittest.TestCase):
    def fixture(self, root):
        child=root/"named child.json"
        child.write_text('{"schema":"srcMoon-I15-reference-contract-v1","fixture":"dispatch only"}')
        value=dict(schema="srcMoon-M2-formal-contract-set-v1",status="APPROVED",
            review=dict(approved=True,reviewer="PROTOCOL FIXTURE",approved_utc="2026-01-01T00:00:00Z"),
            contracts={"I15":{"path":child.name,"sha256":contracts.sha256(child)}})
        bundle=root/"bundle.json";bundle.write_text(json.dumps(value))
        return child,bundle,value

    def test_relative_child_path_and_both_document_hashes_are_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            child,bundle,_=self.fixture(Path(tmp))
            selected,evidence=contracts.resolve_contract(bundle,"I15")
            self.assertEqual(selected,child.resolve())
            self.assertEqual(evidence["input_document"]["sha256"],contracts.sha256(bundle))
            self.assertEqual(evidence["selected_document"]["sha256"],contracts.sha256(child))

    def test_missing_id_preserves_callers_unavailable_status(self):
        class MissingScientificInput(RuntimeError):pass
        with tempfile.TemporaryDirectory() as tmp:
            _,bundle,_=self.fixture(Path(tmp))
            with self.assertRaises(MissingScientificInput):
                contracts.resolve_contract(bundle,"I19",unavailable=MissingScientificInput)

    def test_unreviewed_set_never_selects_even_existing_child(self):
        with tempfile.TemporaryDirectory() as tmp:
            _,bundle,value=self.fixture(Path(tmp));value["review"]["approved"]=False
            bundle.write_text(json.dumps(value))
            with self.assertRaises(contracts.Unavailable):contracts.resolve_contract(bundle,"I15")

    def test_duplicate_keys_and_nonfinite_metadata_are_corruption(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/"bad.json"
            for text in ('{"schema":"one","schema":"two"}','{"limit":NaN}'):
                path.write_text(text)
                with self.assertRaises(ValueError):contracts.resolve_contract(path,"I15")

    def test_changed_child_bytes_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            child,bundle,_=self.fixture(Path(tmp));child.write_text('{"changed":true}')
            with self.assertRaisesRegex(ValueError,"SHA-256 mismatch"):contracts.resolve_contract(bundle,"I15")

    def test_path_only_or_bad_hash_descriptor_is_not_accepted(self):
        with tempfile.TemporaryDirectory() as tmp:
            _,bundle,value=self.fixture(Path(tmp))
            for bad in ({"path":"named child.json"},{"path":"named child.json","sha256":"not-a-hash"}):
                value["contracts"]["I15"]=bad;bundle.write_text(json.dumps(value))
                with self.assertRaises(ValueError):contracts.resolve_contract(bundle,"I15")

    def test_nested_bundle_is_rejected_without_recursing(self):
        with tempfile.TemporaryDirectory() as tmp:
            child,bundle,value=self.fixture(Path(tmp))
            child.write_text('{"schema":"srcMoon-M2-formal-contract-set-v1"}')
            value["contracts"]["I15"]["sha256"]=contracts.sha256(child);bundle.write_text(json.dumps(value))
            with self.assertRaisesRegex(ValueError,"nested"):contracts.resolve_contract(bundle,"I15")

    def test_wrong_document_container_types_are_parser_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for field in ("review","contracts"):
                child,bundle,value=self.fixture(root);value[field]=[];bundle.write_text(json.dumps(value))
                with self.assertRaises(ValueError):contracts.resolve_contract(bundle,"I15")
            child,bundle,value=self.fixture(root);child.write_text('[]')
            value["contracts"]["I15"]["sha256"]=contracts.sha256(child);bundle.write_text(json.dumps(value))
            with self.assertRaises(ValueError):contracts.resolve_contract(bundle,"I15")

    def test_direct_document_is_selected_but_not_certified(self):
        with tempfile.TemporaryDirectory() as tmp:
            child,_,_=self.fixture(Path(tmp))
            selected,evidence=contracts.resolve_contract(child,"I15")
            self.assertEqual(selected,child)
            self.assertNotIn("status",evidence)
            self.assertNotIn("approved",evidence)

    def test_explicit_timezone_is_required_for_review(self):
        with tempfile.TemporaryDirectory() as tmp:
            _,bundle,value=self.fixture(Path(tmp));value["review"]["approved_utc"]="2026-01-01T00:00:00"
            bundle.write_text(json.dumps(value))
            with self.assertRaisesRegex(ValueError,"time zone"):contracts.resolve_contract(bundle,"I15")


if __name__=="__main__":unittest.main(verbosity=2)
