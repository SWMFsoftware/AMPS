#!/usr/bin/env python3
"""Independent lexical/native controls; no synthetic observational evidence."""
from pathlib import Path
import bz2
import csv
import importlib.util
import json
import os
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('new_molecular', ROOT/'utility/chemical_database/new_molecular_exomol.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)


class Administrative(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='molecular-parser-fixture-'); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.paths = {k:self.root/n for k,n in [('definition','def.json'),('definition_text','native.def'),
                    ('states','states.bz2'),('transitions','trans.bz2')]}
        self.value = {'isotopologue':{'iso_formula':'ADMINISTRATIVE NOT PHYSICAL'},
            'dataset':{'name':'PARSER-FIXTURE','version':1,
                'states':{'number_of_states':2,'max_energy':1,
                    'states_file_fields':[{'name':n} for n in ('ID','E','gtot','J','unresolved')]},
                'transitions':{'number_of_transitions':1}}}
        self.paths['definition'].write_text(json.dumps(self.value))
        self.paths['definition_text'].write_text('1 # Version number with format YYYYMMDD\n1 # Higher energy with complete set of transitions (in cm-1)\n')
        self.states = ['1 0.000000 1 0 SOURCE_A','2 2.000000 3 1 SOURCE_B']
        self.transitions = ['2 1 1.2300E-02 2.000000']
        self.write_rows();self.output=self.root/'view'
    def write_rows(self):
        self.paths['states'].write_bytes(bz2.compress(('\n'.join(self.states)+'\n').encode()))
        self.paths['transitions'].write_bytes(bz2.compress(('\n'.join(self.transitions)+'\n').encode()))
    def project(self, hashes=None):
        return m.project(**self.paths,output=self.output,
            expected_hashes=hashes or {k:m.sha(p) for k,p in self.paths.items()})
    def refused(self):
        with self.assertRaises((ValueError,EOFError,OSError)):self.project()
        self.assertFalse(self.output.exists())
    def test_literal_tokens_and_order_preserved(self):
        result=self.project()
        with (self.output/'states.csv').open(newline='')as stream:rows=list(csv.reader(stream))
        self.assertEqual(rows[1:],[r.split() for r in self.states])
        self.assertEqual(rows[0],['ID','E','gtot','J','native_token_5'])
        with (self.output/'transitions.csv').open(newline='')as stream:rows=list(csv.reader(stream))
        self.assertEqual(rows[1:],[r.split() for r in self.transitions])
        self.assertEqual(result['scientific_status'],'PARTIAL_SOURCE_INPUT_NOT_QUALIFIED')
    def test_completeness_energy_is_not_maximum_energy(self):
        r=self.project();self.assertEqual(r['actual_maximum_state_energy_cm-1'],'2.000000')
        self.assertEqual(r['native_transition_completeness_energy']['native_value'],'1')
        self.assertEqual(r['status'],'PASS')
    def test_count_discrepancy_preserved_without_filling(self):
        self.value['dataset']['states']['number_of_states']=3
        self.paths['definition'].write_text(json.dumps(self.value));r=self.project()
        self.assertFalse(r['declared_state_count_matches']);self.assertEqual(r['actual_state_count'],2)
    def test_native_unknown_completeness_is_not_invented(self):
        self.paths['definition_text'].write_text('1 # Version number with format YYYYMMDD\nNA # Higher energy with complete set of transitions (in cm-1)\n')
        r=self.project();self.assertFalse(r['native_transition_completeness_energy']['numeric_value_available'])
        self.assertEqual(r['native_transition_completeness_energy']['native_value'],'NA')
    def test_optional_state_labels_are_not_guessed(self):
        self.value['dataset']['states']['states_file_fields'][4]['name']='v'
        self.paths['definition'].write_text(json.dumps(self.value));r=self.project()
        self.assertEqual(r['state_columns'][4],'native_token_5')
    def test_changed_original_sha_refused(self):
        hashes={k:m.sha(p)for k,p in self.paths.items()};self.states[0]+=' X';self.write_rows()
        with self.assertRaises(ValueError):self.project(hashes)
        self.assertFalse(self.output.exists())
    def test_missing_sha_refused(self):
        with self.assertRaises(ValueError):self.project({'states':m.sha(self.paths['states'])})
        self.assertFalse(self.output.exists())
    def test_mixed_definition_versions_refused(self):
        self.value['dataset']['version']=2;self.paths['definition'].write_text(json.dumps(self.value));self.refused()
    def test_duplicate_state_id_refused(self):
        self.states[1]='1 2.000000 3 1 SOURCE_B';self.write_rows();self.refused()
    def test_missing_transition_state_refused(self):
        self.transitions=['3 1 1e-2 2'];self.write_rows();self.refused()
    def test_negative_A_not_clipped(self):
        self.transitions=['2 1 -1e-2 2'];self.write_rows();self.refused()
    def test_nonfinite_energy_refused(self):
        self.states[1]='2 NaN 3 1 SOURCE_B';self.write_rows();self.refused()
    def test_ragged_row_not_repaired(self):
        self.states[1]='2 2.000000 3 1';self.write_rows();self.refused()
    def test_truncated_bzip2_refused(self):
        self.paths['states'].write_bytes(self.paths['states'].read_bytes()[:12]);self.refused()
    def test_existing_output_preserved(self):
        self.project();before=(self.output/'states.csv').read_bytes()
        with self.assertRaises(ValueError):self.project()
        self.assertEqual(before,(self.output/'states.csv').read_bytes())


class ActualSources(unittest.TestCase):
    requires_authoritative_source = True
    @classmethod
    def setUpClass(cls):
        selected=os.environ.get('AMPS_NEW_MOLECULAR_SOURCE_ROOT')
        if not selected: raise unittest.SkipTest('explicit AMPS_NEW_MOLECULAR_SOURCE_ROOT not selected')
        cls.sources=Path(selected).resolve()
        if not cls.sources.is_dir(): raise unittest.SkipTest('actual ExoMol originals unavailable')
        cls.receipts={r['path']:r for r in json.loads((cls.sources.parent/'acquisition_manifest.json').read_text())['records']}
    def test_actual_original_tokens_and_transition_references(self):
        expected={'16O-1H_p__MoLLIST':(815,12044),'16O-1H__MoLLIST-OH':(1878,54276),
                  '1H2__ARLR':(2051,32915),'1H2__RACPPK':(302,4712)}
        for stem,counts in expected.items():
            with self.subTest(source=stem),tempfile.TemporaryDirectory()as d:
                paths={k:self.sources/(stem+'.'+ext)for k,ext in
                    [('definition','def.json'),('definition_text','def'),('states','states.bz2'),('transitions','trans.bz2')]}
                hashes={}
                for key,path in paths.items():
                    receipt=self.receipts['exomol/'+path.name]
                    self.assertEqual(receipt['http_status'],200);self.assertTrue(receipt['complete'])
                    self.assertEqual(m.sha(path),receipt['sha256']);hashes[key]=receipt['sha256']
                result=m.project(**paths,output=Path(d)/'view',expected_hashes=hashes)
                self.assertEqual((result['actual_state_count'],result['actual_transition_count']),counts)
                for name,key in [('states.csv','states'),('transitions.csv','transitions')]:
                    with (Path(d)/'view'/name).open(newline='')as stream:csv_rows=list(csv.reader(stream))[1:]
                    native=[line.split()for line in bz2.decompress(paths[key].read_bytes()).decode('ascii').splitlines()]
                    self.assertEqual(csv_rows,native)
                self.assertEqual(result['missing_transition_state_references'],0)
    def test_h2_ground_state_is_explicit_source_zero(self):
        rows=[line.split()for line in bz2.decompress((self.sources/'1H2__ARLR.states.bz2').read_bytes()).decode().splitlines()]
        self.assertEqual(rows[0],['1','0.000000','1','0','X','0'])


if __name__=='__main__':unittest.main(verbosity=2)
