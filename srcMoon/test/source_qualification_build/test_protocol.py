#!/usr/bin/env python3
"""Administrative path/mapping guards only; no native build or fake physics."""
import json
import shutil
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import build_m3_source_qualification as production


class Protocol(unittest.TestCase):
    def test_actual_canonical_compiled_ion_symbol(self):
        text=(production.ROOT/'input/species.input').read_text()
        self.assertIn('#component=NA_PLUS',text)
        self.assertEqual(production.EXPECTED_MAPPING,{'Na':'23Na:ground','NA_PLUS':'23Na+:ground'})

    def test_exact_explicit_mapping(self):
        with tempfile.TemporaryDirectory()as d:
            p=Path(d)/'mapping.json';p.write_text(json.dumps(production.EXPECTED_MAPPING))
            self.assertEqual(production.validate_mapping(p),production.EXPECTED_MAPPING)

    def test_three_species_capture_mapping(self):
        # This recipe owns two independently weighted transported daughters.
        # H is the third species. An H+ packet would count the incoming bath
        # twice and would conflate a charged projectile with the outgoing atom.
        with tempfile.TemporaryDirectory() as directory:
            mapping = Path(directory) / 'mapping.json'
            mapping.write_text(json.dumps(production.GROUND_CAPTURE_MAPPING))
            self.assertEqual(production.validate_mapping(mapping), production.GROUND_CAPTURE_MAPPING)
            for state in ('1H+:ground', '2H:ground'):
                mapping.write_text(json.dumps({**production.EXPECTED_MAPPING, 'H': state}))
                with self.subTest(state=state), self.assertRaises(ValueError):
                    production.validate_mapping(mapping)

    def test_wrong_case_unknown_state_or_reordered_mapping_refused(self):
        with tempfile.TemporaryDirectory()as d:
            p=Path(d)/'mapping.json'
            for value in ({'Na':'23Na:ground','Na_PLUS':'23Na+:ground'},
                {'Na':'23Na:ground','NA_PLUS':'23Na+:unresolved'},
                {'NA_PLUS':'23Na+:ground','Na':'23Na:ground'}):
                p.write_text(json.dumps(value))
                with self.subTest(value=value),self.assertRaises(ValueError):production.validate_mapping(p)

    def test_duplicate_or_nonfinite_json_refused(self):
        with tempfile.TemporaryDirectory()as d:
            p=Path(d)/'mapping.json'
            for text in ('{"Na":1,"Na":2}','{"value":NaN}'):
                p.write_text(text)
                with self.subTest(text=text),self.assertRaises(ValueError):production.strict(p)

    def test_active_checkout_refused_before_inventory_or_build(self):
        with patch.object(production,'checkout_identity',side_effect=AssertionError('must not inspect/build active root')):
            with self.assertRaises(ValueError):production.validate_checkout(production.ROOT,production.ROOT)

    def test_nested_checkout_refused_before_inventory_or_build(self):
        with patch.object(production,'checkout_identity',side_effect=AssertionError('must not inspect/build nested root')):
            with self.assertRaises(ValueError):production.validate_checkout(production.ROOT,production.ROOT/'test_output/private')

    def test_missing_private_checkout_refused(self):
        with tempfile.TemporaryDirectory()as d:
            with self.assertRaises(ValueError):production.validate_checkout(production.ROOT,Path(d)/'missing')

    def fixture_checkout(self,directory):
        # Administrative files only; no fake executable or PASS receipt exists.
        target=Path(directory)/'private';target.mkdir();(target/'Config.pl').write_text('ADMINISTRATIVE FIXTURE ONLY\n')
        for relative in production.WORKFLOW:
            p=target/relative;p.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(production.ROOT/relative,p)
        return target

    def test_stale_source_map_refused_before_build(self):
        with tempfile.TemporaryDirectory()as d:
            target=self.fixture_checkout(d)
            with patch.object(production,'checkout_identity',return_value={'changed.cpp':'different'}):
                with self.assertRaises(ValueError):production.validate_checkout(production.ROOT,target,{'source.cpp':'pinned'})

    def test_unrebased_site_path_refused_before_build(self):
        with tempfile.TemporaryDirectory()as d:
            target=self.fixture_checkout(d);(target/'Makefile').write_text('PATH_TO_ACTIVE='+str(production.ROOT)+'\n')
            with patch.object(production,'checkout_identity',return_value={}):
                with self.assertRaises(ValueError):production.validate_checkout(production.ROOT,target,{})

    def test_missing_generated_site_makefile_refused_before_build(self):
        with tempfile.TemporaryDirectory()as d:
            target=self.fixture_checkout(d);(target/'Makefile.def').write_text('include Makefile.test\n')
            with patch.object(production,'checkout_identity',return_value={}):
                with self.assertRaisesRegex(FileNotFoundError,'Config.pl -install'):
                    production.validate_checkout(production.ROOT,target,{})


if __name__=='__main__':unittest.main(verbosity=2)
