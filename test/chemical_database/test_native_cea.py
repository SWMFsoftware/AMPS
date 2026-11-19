#!/usr/bin/env python3
"""NASA source-index verification, with labelled parser fixtures and real QA.

Fixtures below exercise text layout/integrity only; none are observations or
physical state profiles. Optional --source (or AMPS_CEA_SOURCE) admits the
actual hash-pinned acquired NASA file for additional lexical/source checks.
The preserved eleven condensed-phase domain issues are expected source facts,
not repaired data or relaxed kinetic/validation criteria.
"""
import argparse
import csv
from collections import Counter
from dataclasses import FrozenInstanceError
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'utility/chemical_database'))
import native_cea as production

SOURCE=os.environ.get('AMPS_CEA_SOURCE')
EXPECTED='fa7746572952d74e249e818a82a35c113829742fb421a308e167185528884363'
NAMES=('e-','H','H+','H2','H2+','H2O','H2O+','OH','OH+','O','O+','He','He+','Ne','Ne+','Na','Na+','Ar','Ar+')


def fixture_record(name='H', *, n=1, lo=200., hi=1000., delta=6197.428,
                   phase=0, ignored=' '*16):
    """Construct an explicit administrative FORMAT-layout fixture, not data."""
    header=f'{name:<15}   DECLARED PARSER-ONLY FIXTURE\n'
    formula=f'{"H":<2}{1.:6.2f}'+''.join(f'{"":2}{0.:6.2f}'for _ in range(4))
    metadata=f'{n:2d} {"g 1/00":6s} {formula}{phase:2d}{1.00794:13.7f}{0.:15.3f}\n'
    exponents=''.join(f'{v:5.1f}'for v in(-2,-1,0,1,2,3,4,0))
    row=f'{lo:11.3f}{hi:11.3f}7{exponents}  {delta:15.3f}\n'
    number=lambda value:f'{value:16.9E}'.replace('E','D')
    left=''.join(number(v)for v in(0,0,2.5,0,0))+'\n'
    right=''.join(number(v)for v in(0,0))+ignored+''.join(number(v)for v in(-745.375,-1.1))+'\n'
    return header+metadata+(row+left+right)*n if n else header+metadata+f'{298.15:11.3f}\n'


def fixture(*records, reactants=''):
    schedule=''.join(f'{t:10.2f}'for t in(200,1000,6000,20000))+'fixture-date\n'
    return ('! ADMINISTRATIVE PARSER FIXTURE, NOT OBSERVATIONS\nthermo\n'+schedule+
        ''.join(records or(fixture_record(),))+'END PRODUCTS\n'+reactants+'END REACTANTS\n').encode('ascii')


class LexicalProtocol(unittest.TestCase):
    def test_finite_native_tokens_and_nine_coefficients(self):
        c=production.parse_thermo_bytes(fixture());r=c.records[0];i=r['intervals'][0]
        self.assertEqual(r['molecular_weight_native_u'],'1.0079400')
        self.assertEqual(i['coefficients_a1_to_a7_b1_b2_native'][2],'2.500000000D+00')
        self.assertEqual(len(i['coefficients_a1_to_a7_b1_b2_native']),9)

    def test_ignored_slot_preserved_without_new_polynomial_term(self):
        c=production.parse_thermo_bytes(fixture(fixture_record(ignored='IGNORED FIXTURE!')))
        i=c.records[0]['intervals'][0]
        self.assertEqual(i['native_compiler_ignored_columns_33_to_48'],'IGNORED FIXTURE!')
        self.assertEqual(len(i['coefficients_a1_to_a7_b1_b2_native']),9)

    def test_unavailable_zero_marker_is_not_a_physical_zero(self):
        i=production.parse_thermo_bytes(fixture(fixture_record(delta=0))).records[0]['intervals'][0]
        self.assertEqual(i['h298_minus_h0_native_J_mol'],'0.000')
        self.assertEqual(i['h298_minus_h0_status'],'UNAVAILABLE_NATIVE_ZERO_MARKER')

    def test_native_reversed_domain_is_preserved_and_fails_domain_quality(self):
        d=production.parse_thermo_bytes(fixture(fixture_record(lo=400,hi=300))).as_dict()
        self.assertEqual(d['native_domain_quality_status'],'FAIL')
        self.assertEqual(d['records'][0]['intervals'][0]['temperature_min_K'],'400.000')
        self.assertFalse(d['records'][0]['intervals'][0]['temperature_limits_increasing'])

    def test_reactant_only_record_preserved(self):
        c=production.parse_thermo_bytes(fixture(reactants=fixture_record('R',n=0)))
        self.assertEqual(c.records[1]['interval_count'],0)
        self.assertEqual(c.records[1]['reactant_reference_temperature_K'],'298.150')
        self.assertEqual(c.records[1]['section'],'reactants')

    def test_immutable_snapshot_and_defensive_nested_copy(self):
        c=production.parse_thermo_bytes(fixture());r=c.records[0]
        r['intervals'][0]['temperature_min_K']='CORRUPTED RETURNED COPY'
        self.assertEqual(c.records[0]['intervals'][0]['temperature_min_K'],'200.000')
        with self.assertRaises(FrozenInstanceError):c._snapshot='{}'

    def test_duplicate_native_names_are_not_merged(self):
        c=production.parse_thermo_bytes(fixture(fixture_record(),fixture_record()))
        self.assertEqual(len(c.records),2)
        with self.assertRaises(LookupError):production.selected_gas_record(c,'H')

    def test_missing_native_name_refused(self):
        with self.assertRaises(LookupError):production.selected_gas_record(production.parse_thermo_bytes(fixture()),'UNAVAILABLE')

    def test_selection_order_preserved(self):
        c=production.parse_thermo_bytes(fixture(fixture_record('A'),fixture_record('B')))
        text=production.selected_species_csv(c,['B','A'])
        self.assertEqual([s.split(',')[0]for s in text.splitlines()[1:]],['B','A'])

    def test_duplicate_projection_request_refused(self):
        with self.assertRaises(ValueError):production.selected_species_csv(production.parse_thermo_bytes(fixture()),['H','H'])

    def test_coefficient_projection_retains_native_fields_and_domain(self):
        c=production.parse_thermo_bytes(fixture(fixture_record(lo=400,hi=300)))
        rows=list(csv.DictReader(io.StringIO(production.coefficients_csv(c,['H']))))
        self.assertEqual(len(rows),1)
        self.assertEqual(rows[0]['a3_native'],'2.500000000D+00')
        self.assertEqual(rows[0]['native_domain_quality'],'FAIL')
        self.assertEqual(rows[0]['temperature_min_K'],'400.000')

    def test_coefficient_projection_preserves_declared_species_order(self):
        c=production.parse_thermo_bytes(fixture(fixture_record('A'),fixture_record('B')))
        rows=list(csv.DictReader(io.StringIO(production.coefficients_csv(c,['B','A']))))
        self.assertEqual([row['name']for row in rows],['B','A'])

    def test_exact_explicit_resource_limit(self):
        data=fixture();self.assertEqual(production.parse_thermo_bytes(data,maximum_bytes=len(data)).source_sha256,hashlib.sha256(data).hexdigest())
        with self.assertRaises(ValueError):production.parse_thermo_bytes(data,maximum_bytes=len(data)-1)

    def test_invalid_resource_limits_refused(self):
        for value in(True,0,-1,'10'):
            with self.subTest(value=value),self.assertRaises(ValueError):production.parse_thermo_bytes(fixture(),maximum_bytes=value)

    def test_nonbytes_refused(self):
        with self.assertRaises(ValueError):production.parse_thermo_bytes(fixture().decode('ascii'))

    def test_nonascii_refused(self):
        with self.assertRaises(ValueError):production.parse_thermo_bytes(fixture()+b'\xff')

    def test_no_thermo_sentinel_refused(self):
        with self.assertRaises(ValueError):production.parse_thermo_bytes(fixture().replace(b'thermo\n',b'other\n'))

    def test_duplicate_thermo_sentinel_refused(self):
        with self.assertRaises(ValueError):production.parse_thermo_bytes(fixture()+b'thermo\n')

    def test_unknown_prefix_refused(self):
        with self.assertRaises(ValueError):production.parse_thermo_bytes(b'UNEXPLAINED\n'+fixture())

    def test_missing_products_end_refused(self):
        with self.assertRaises(ValueError):production.parse_thermo_bytes(fixture().replace(b'END PRODUCTS\n',b''))

    def test_missing_reactants_end_refused(self):
        with self.assertRaises(ValueError):production.parse_thermo_bytes(fixture().replace(b'END REACTANTS\n',b''))

    def test_unknown_trailing_record_refused(self):
        with self.assertRaises(ValueError):production.parse_thermo_bytes(fixture()+b'UNEXPLAINED\n')

    def test_truncated_interval_refused(self):
        rows=fixture().splitlines(keepends=True);rows.pop(7)
        with self.assertRaises(ValueError):production.parse_thermo_bytes(b''.join(rows))

    def test_nonnumeric_mass_refused(self):
        rows=fixture().splitlines(keepends=True);row=rows[4];rows[4]=row[:52]+b'XXXXXXXXXXXXX'+row[65:]
        with self.assertRaises(ValueError):production.parse_thermo_bytes(b''.join(rows))

    def test_nonfinite_thermal_increment_refused(self):
        rows=fixture().splitlines(keepends=True);row=rows[5];rows[5]=row[:65]+b'            NaN\n'
        with self.assertRaises(ValueError):production.parse_thermo_bytes(b''.join(rows))

    def test_wrong_cp_count_refused(self):
        rows=fixture().splitlines(keepends=True);row=rows[5];rows[5]=row[:22]+b'8'+row[23:]
        with self.assertRaises(ValueError):production.parse_thermo_bytes(b''.join(rows))


class ActualSource(unittest.TestCase):
    # The global CORE protocol promises no site-specific originals. Full
    # source verification still executes these methods or reports SKIPPED.
    requires_authoritative_source = True
    @classmethod
    def setUpClass(cls):
        if not SOURCE:raise unittest.SkipTest('actual NASA source path not explicitly supplied; no replacement data manufactured')
        cls.native=Path(SOURCE).read_bytes()
        if hashlib.sha256(cls.native).hexdigest()!=EXPECTED:raise RuntimeError('actual pinned NASA source bytes changed')
        cls.catalogue=production.parse_thermo_bytes(cls.native)
        cls.data=cls.catalogue.as_dict()

    def test_actual_2111_records_and_complete_native_boundaries(self):
        self.assertEqual(len(self.data['records']),2111)
        self.assertEqual(self.data['native_line_count'],15802)
        self.assertEqual(self.data['native_global_date'],'9/8/2021')

    def test_all_native_blocks_byte_roundtrip(self):
        rows=self.native.splitlines(keepends=True)
        for record in self.data['records']:
            self.assertEqual(hashlib.sha256(b''.join(rows[record['native_start_line']-1:record['native_end_line']])).hexdigest(),record['native_block_sha256'])

    def test_real_nineteen_selected_gas_domains_are_ordered(self):
        for name in NAMES:
            record=production.selected_gas_record(self.catalogue,name)
            self.assertTrue(all(i['temperature_limits_increasing']for i in record['intervals']))
            self.assertTrue(all(i['h298_minus_h0_status']=='REPORTED'for i in record['intervals']))

    def test_native_mass_convention_preserved_not_silently_replaced(self):
        self.assertEqual(production.selected_gas_record(self.catalogue,'Na')['molecular_weight_native_u'],'22.9897700')
        self.assertEqual(production.selected_gas_record(self.catalogue,'Na+')['molecular_weight_native_u'],'22.9892214')

    def test_source_condensed_domain_issues_remain_explicit_fail(self):
        self.assertEqual(self.data['native_domain_quality_status'],'FAIL')
        self.assertEqual(len(self.data['native_quality_issues']),11)
        self.assertEqual(self.data['native_quality_issues'][0]['name'],'Br2(cr)')

    def test_actual_coefficient_projection_count_and_exact_tokens(self):
        rows=list(csv.DictReader(io.StringIO(production.coefficients_csv(self.catalogue,NAMES))))
        expected=[]
        for name in NAMES:
            record=production.selected_gas_record(self.catalogue,name)
            for interval in record['intervals']:expected.append(interval)
        self.assertEqual(len(rows),len(expected))
        for row,interval in zip(rows,expected):
            self.assertEqual([row[name+'_native']for name in('a1','a2','a3','a4','a5','a6','a7','b1','b2')],interval['coefficients_a1_to_a7_b1_b2_native'])


class Recorded(unittest.TextTestResult):
    def __init__(self,*args,**kwargs):super().__init__(*args,**kwargs);self.rows=[]
    def addSuccess(self,test):super().addSuccess(test);self.rows.append({'test':str(test),'status':'PASS'})
    def addFailure(self,test,error):super().addFailure(test,error);self.rows.append({'test':str(test),'status':'FAIL'})
    def addError(self,test,error):super().addError(test,error);self.rows.append({'test':str(test),'status':'ERROR'})
    def addSkip(self,test,reason):super().addSkip(test,reason);self.rows.append({'test':str(test),'status':'SKIPPED','reason':reason})


def main():
    global SOURCE
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args();out=args.output_dir.resolve()
    out.mkdir(parents=True,exist_ok=False)
    SOURCE=str(args.source.resolve())if args.source else SOURCE
    started=time.monotonic();suite=unittest.TestSuite(unittest.defaultTestLoader.loadTestsFromTestCase(owner)for owner in(LexicalProtocol,ActualSource))
    with(out/'unittest.log').open('w')as stream:result=unittest.TextTestRunner(stream=stream,verbosity=2,resultclass=Recorded).run(suite)
    counts={status:sum(row['status']==status for row in result.rows)for status in('PASS','FAIL','ERROR','SKIPPED')}
    status='ERROR'if counts['ERROR']else'FAIL'if counts['FAIL']else'PASS'if counts['PASS']else'SKIPPED'
    value={'schema':'amps-native-cea-source-tests-v1','status':status,'counts':counts,'checks':result.rows,
        'scope':'native lexical/parser integrity verification only; NOT individual physical state qualification',
        'native_full_source_domain_quality':'FAIL (eleven preserved condensed-phase limits)'if SOURCE else'NOT EVALUATED',
        'source_revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in(Path(__file__).resolve(),Path(production.__file__).resolve())},
        'actual_input':{'path':SOURCE,'sha256':hashlib.sha256(Path(SOURCE).read_bytes()).hexdigest()}if SOURCE else None,
        'elapsed_s':time.monotonic()-started,'runner':str(Path(__file__).resolve()),'log':str(out/'unittest.log')}
    (out/'result.json').write_text(json.dumps(value,indent=2)+'\n');print(json.dumps(value['counts']),status)
    return {'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[status]


if __name__=='__main__':raise SystemExit(main())
