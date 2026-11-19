#!/usr/bin/env python3
"""Shared CHIANTI lexical/API verification, never observational validation.

Administrative fixed-column fixtures below test parser/ownership/failure
semantics only. Optional real --levels/--radiative inputs are independently
hash-pinned O I originals; absence is reported SKIPPED, not replaced by a
fixture called data. The production parser is the sole conversion path.
"""
import argparse
import csv
from dataclasses import FrozenInstanceError
from fractions import Fraction
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'utility/chemical_database'))
import native_chianti as production

LEVELS = os.environ.get('AMPS_CHIANTI_LEVELS')
RADIATIVE = os.environ.get('AMPS_CHIANTI_RADIATIVE')
LEVEL_SHA = '230bbdcd95c269fc0eea04a1fd37f534c1a42979d2e29c1b346e49dd71c688da'
RATE_SHA = 'f24efc2bd0fa5120baca6dce06ffb7831ae2a59e17cc186c233a43b5145866a8'


def level(index, observed='0.000', theory='0.000', *, config='PARSER FIXTURE',
          label='', spin=1, term='S', j='0.0'):
    """Declared fixed-format administrative row, not an atomic measurement."""
    return f'{index:7d}{config:<30}{label:<5}{spin:5d}{term:>5}{j:>5}{observed:>15}{theory:>15}\n'


def rate(lower=1, upper=2, wave='100.000', gf='0.000e+00', a='1.000e+00', label='PARSER FIXTURE'):
    return f'{lower:5d}{upper:5d}{wave:>15}{gf:>15}{a:>15}   {label}\n'


def fixture(*, levels=None, rates=None, footer='% ADMINISTRATIVE FIXTURE; NOT OBSERVATIONS\n'):
    return ((levels or level(1)+level(2,'100.000','110.000'))+' -1\n'+footer).encode('ascii'), \
           ((rates or rate())+' -1\n'+footer+' -1\n').encode('ascii')


class LexicalProtocol(unittest.TestCase):
    def test_native_tokens_and_units_preserved(self):
        c = production.parse_ion_bytes(*fixture())
        self.assertEqual(c.levels[1]['observed_energy_cm_inverse_native'], '100.000')
        self.assertEqual(c.transitions[0]['decay_rate_s_inverse_native'], '1.000e+00')
        self.assertEqual(c.transitions[0]['gf_native'], '0.000e+00')

    def test_source_row_hashes_include_original_newlines(self):
        left, right = fixture(); c = production.parse_ion_bytes(left, right)
        self.assertEqual(c.levels[0]['native_row_sha256'], hashlib.sha256(left.splitlines(keepends=True)[0]).hexdigest())
        self.assertEqual(c.transitions[0]['native_row_sha256'], hashlib.sha256(right.splitlines(keepends=True)[0]).hexdigest())

    def test_immutable_owned_defensive_snapshot(self):
        c = production.parse_ion_bytes(*fixture()); rows = c.levels
        rows[0]['configuration_native'] = 'MUTATED RETURN VALUE'
        self.assertEqual(c.levels[0]['configuration_native'], 'PARSER FIXTURE')
        d=c.as_dict(); d['transitions'][0]['decay_rate_s_inverse_native']='99'
        self.assertEqual(c.transitions[0]['decay_rate_s_inverse_native'], '1.000e+00')
        with self.assertRaises(FrozenInstanceError): c._snapshot='{}'

    def test_observed_missingness_never_filled_from_theory(self):
        c=production.parse_ion_bytes(*fixture(levels=level(1)+level(2,'-1.000','100.000')))
        r=c.levels[1]
        self.assertEqual(r['observed_energy_cm_inverse_native'],'-1.000')
        self.assertEqual(r['observed_energy_status'],'UNAVAILABLE_NATIVE_MINUS_ONE')
        self.assertEqual(r['theoretical_energy_cm_inverse_native'],'100.000')

    def test_two_missing_energies_retained_as_source_quality_fail(self):
        d=production.parse_ion_bytes(*fixture(levels=level(1)+level(2,'-1.000','-1.000'))).as_dict()
        self.assertEqual(d['native_domain_quality_status'],'FAIL')
        self.assertEqual(d['native_quality_issues'][0]['code'],'BOTH_LEVEL_ENERGIES_UNAVAILABLE')

    def test_negative_wavelength_provenance_preserved(self):
        r=production.parse_ion_bytes(*fixture(rates=rate(wave='-100.000'))).transitions[0]
        self.assertEqual(r['signed_wavelength_angstrom_native'],'-100.000')
        self.assertEqual(r['wavelength_status'],'THEORETICAL_OR_UNOBSERVED_ENDPOINT')

    def test_zero_wavelength_is_not_invented_single_photon(self):
        r=production.parse_ion_bytes(*fixture(rates=rate(wave='0.000'))).transitions[0]
        self.assertEqual(r['rate_kind'],'UNRESOLVED_TWO_PHOTON_OR_AUTOIONIZATION')

    def test_zero_rate_not_deleted_and_no_invented_uncertainty(self):
        c=production.parse_ion_bytes(*fixture(rates=rate(a='0.000e+00')))
        self.assertEqual(len(c.transitions),1)
        self.assertEqual(c.transitions[0]['decay_rate_s_inverse_native'],'0.000e+00')
        self.assertTrue(all(x['uncertainty_status']=='NOT_SUPPLIED_IN_FILE' for x in c.levels+c.transitions))

    def test_fortran_D_token_preserved(self):
        r=production.parse_ion_bytes(*fixture(rates=rate(a='1.234D-06'))).transitions[0]
        self.assertEqual(r['decay_rate_s_inverse_native'],'1.234D-06')

    def test_freeform_comments_retained_as_documented(self):
        footer='arbitrary original credit\n123 numeric-looking footer is a comment\n'
        d=production.parse_ion_bytes(*fixture(footer=footer)).as_dict()
        self.assertIn('123 numeric-looking footer is a comment',d['level_comments'])
        self.assertEqual(len(d['transitions']),1)

    def test_index_order_is_not_imposed_as_energy_order(self):
        d=production.parse_ion_bytes(*fixture(levels=level(2)+level(1,'100.000','100.000'),rates=rate(2,1))).as_dict()
        self.assertEqual([x['level_index'] for x in d['levels']],[2,1])
        self.assertEqual(d['native_domain_quality_status'],'PASS')

    def test_native_energy_order_defect_flagged_without_repair(self):
        c=production.parse_ion_bytes(*fixture(levels=level(1,'100.000','100.000')+level(2)))
        self.assertEqual(c.as_dict()['native_domain_quality_status'],'FAIL')
        self.assertEqual(c.transitions[0]['signed_wavelength_angstrom_native'],'100.000')

    def test_explicit_resource_limit_inclusive(self):
        a,b=fixture(); self.assertEqual(len(production.parse_ion_bytes(a,b,maximum_bytes=len(a)+len(b)).levels),2)
        with self.assertRaises(ValueError):production.parse_ion_bytes(a,b,maximum_bytes=len(a)+len(b)-1)

    def test_bad_resource_limits_refused(self):
        for x in (True,0,-1,'100'):
            with self.subTest(x=x),self.assertRaises(ValueError):production.parse_ion_bytes(*fixture(),maximum_bytes=x)

    def test_mutable_or_text_inputs_refused(self):
        a,b=fixture()
        for x in (a.decode(),bytearray(a),None):
            with self.subTest(x=type(x).__name__),self.assertRaises(ValueError):production.parse_ion_bytes(x,b)

    def test_nonascii_refused(self):
        a,b=fixture()
        with self.assertRaises(ValueError):production.parse_ion_bytes(a+b'\xff',b)

    def test_missing_sentinel_refused(self):
        a,b=fixture()
        for x,y in ((a.replace(b' -1\n',b''),b),(a,b.replace(b' -1\n',b''))):
            with self.subTest(),self.assertRaises(ValueError):production.parse_ion_bytes(x,y)

    def test_empty_native_data_refused(self):
        with self.assertRaises(ValueError):production.parse_ion_bytes(b'-1\n',fixture()[1])

    def test_truncated_elvlc_refused(self):
        a,b=fixture()
        with self.assertRaises(ValueError):production.parse_ion_bytes(a.replace(b'0.000\n',b'0.00\n',1),b)

    def test_truncated_wgfa_refused(self):
        a,b=fixture()
        with self.assertRaises(ValueError):production.parse_ion_bytes(a,b[:30]+b'\n-1\n')

    def test_duplicate_level_refused(self):
        with self.assertRaises(ValueError):production.parse_ion_bytes(*fixture(levels=level(1)+level(1)))

    def test_duplicate_transition_refused(self):
        with self.assertRaises(ValueError):production.parse_ion_bytes(*fixture(rates=rate()+rate()))

    def test_unknown_endpoint_refused(self):
        with self.assertRaises(ValueError):production.parse_ion_bytes(*fixture(rates=rate(1,9)))

    def test_self_transition_refused(self):
        with self.assertRaises(ValueError):production.parse_ion_bytes(*fixture(rates=rate(1,1)))

    def test_nonfinite_rate_refused(self):
        for token in ('NaN','inf','-Infinity'):
            with self.subTest(token=token),self.assertRaises(ValueError):production.parse_ion_bytes(*fixture(rates=rate(a=token)))

    def test_negative_rate_or_gf_refused(self):
        for kw in ({'a':'-1.0'},{'gf':'-0.1'}):
            with self.subTest(kw=kw),self.assertRaises(ValueError):production.parse_ion_bytes(*fixture(rates=rate(**kw)))

    def test_invalid_quantum_numbers_refused(self):
        for kw in ({'j':'0.3'},{'j':'-0.5'},{'spin':0}):
            with self.subTest(kw=kw),self.assertRaises(ValueError):production.parse_ion_bytes(*fixture(levels=level(1,**kw)+level(2,'100','100')))

    def test_unrecognized_negative_energy_refused(self):
        with self.assertRaises(ValueError):production.parse_ion_bytes(*fixture(levels=level(1)+level(2,'-2.000','100')))

    def test_explicit_selection_order_preserved(self):
        c=production.parse_ion_bytes(*fixture())
        rows=list(csv.DictReader(io.StringIO(production.levels_csv(c,[2,1]))))
        self.assertEqual([x['level_index'] for x in rows],['2','1'])

    def test_invalid_selections_refused(self):
        c=production.parse_ion_bytes(*fixture())
        for x in ([],[1,1],[True],[1.0],['1']):
            with self.subTest(x=x),self.assertRaises(ValueError):production.levels_csv(c,x)

    def test_unavailable_selection_refused(self):
        with self.assertRaises(LookupError):production.levels_csv(production.parse_ion_bytes(*fixture()),[9])

    def test_partial_downward_cascade_refused_not_renormalized(self):
        c=production.parse_ion_bytes(*fixture())
        with self.assertRaises(ValueError):production.transitions_csv(c,[2])
        self.assertEqual(len(list(csv.DictReader(io.StringIO(production.transitions_csv(c,[1,2]))))),1)


class ActualSource(unittest.TestCase):
    requires_authoritative_source=True
    @classmethod
    def setUpClass(cls):
        if not LEVELS and not RADIATIVE:raise unittest.SkipTest('actual O I source inputs not explicitly supplied; no replacement observations')
        if not LEVELS or not RADIATIVE:raise RuntimeError('both actual source files must be supplied together')
        cls.left,cls.right=Path(LEVELS).read_bytes(),Path(RADIATIVE).read_bytes()
        if hashlib.sha256(cls.left).hexdigest()!=LEVEL_SHA or hashlib.sha256(cls.right).hexdigest()!=RATE_SHA:
            raise RuntimeError('actual CHIANTI O I source SHA256 mismatch')
        cls.c=production.parse_ion_bytes(cls.left,cls.right)

    def test_full_native_inventory(self):
        self.assertEqual((len(self.c.levels),len(self.c.transitions)),(7,16))
        self.assertEqual(self.c.as_dict()['native_domain_quality_status'],'PASS')

    def test_independent_whitespace_native_energy_tokens(self):
        for row,text in zip(self.c.levels,self.left.decode().splitlines()[:7]):
            tokens=text.split()
            self.assertEqual(row['observed_energy_cm_inverse_native'],tokens[-2])
            self.assertEqual(row['theoretical_energy_cm_inverse_native'],tokens[-1])

    def test_independent_whitespace_native_radiative_tokens(self):
        for row,text in zip(self.c.transitions,self.right.decode().splitlines()[:16]):
            tokens=text.split()
            self.assertEqual((row['lower_index'],row['upper_index']),(int(tokens[0]),int(tokens[1])))
            self.assertEqual([row[n] for n in ('signed_wavelength_angstrom_native','gf_native','decay_rate_s_inverse_native')],tokens[2:5])

    def test_source_row_byte_identities(self):
        for rows,native in ((self.c.levels,self.left),(self.c.transitions,self.right)):
            original=native.splitlines(keepends=True)
            for row in rows:self.assertEqual(row['native_row_sha256'],hashlib.sha256(original[row['native_line']-1]).hexdigest())

    def test_five_level_projection_and_nine_original_links(self):
        levels=list(csv.DictReader(io.StringIO(production.levels_csv(self.c,[1,2,3,4,5]))))
        links=list(csv.DictReader(io.StringIO(production.transitions_csv(self.c,[1,2,3,4,5]))))
        self.assertEqual((len(levels),len(links)),(5,9))
        pairs={(int(r['lower_index']),int(r['upper_index'])) for r in links}
        self.assertEqual(pairs,{(1,2),(1,3),(2,3),(1,4),(2,4),(3,4),(1,5),(2,5),(4,5)})

    def test_weak_exits_not_lost(self):
        rates={(r['lower_index'],r['upper_index']):r['decay_rate_s_inverse_native'] for r in self.c.transitions}
        self.assertEqual(rates[(3,4)],'6.181e-07')
        self.assertEqual(rates[(1,5)],'3.494e-04')

    def test_exact_source_lifetime_and_branch_arithmetic(self):
        source={tuple(map(int,line.split()[:2])):Fraction(line.split()[4]) for line in self.right.decode().splitlines()[:9]}
        parsed={(r['lower_index'],r['upper_index']):Fraction(r['decay_rate_s_inverse_native']) for r in self.c.transitions if r['upper_index']<=5}
        self.assertEqual(source,parsed)
        self.assertEqual(sum(v for (lo,up),v in parsed.items() if up==4),Fraction('0.0086046181'))
        self.assertEqual(sum(v for (lo,up),v in parsed.items() if up==5),Fraction('1.3387694'))
        self.assertEqual(parsed[(4,5)]/parsed[(2,5)],Fraction(1_260_000,78_420))

    def test_native_theory_and_observed_credits_retained(self):
        d=self.c.as_dict()
        self.assertTrue(any('Froese Fischer' in line for line in d['radiative_comments']))
        self.assertTrue(any('NIST' in line for line in d['level_comments']))
        self.assertNotEqual(self.c.levels[3]['observed_energy_cm_inverse_native'],self.c.levels[3]['theoretical_energy_cm_inverse_native'])

    def test_no_error_bars_inferred(self):
        self.assertTrue(all(r['uncertainty_status']=='NOT_SUPPLIED_IN_FILE' for r in self.c.levels+self.c.transitions))


class Recorded(unittest.TextTestResult):
    def __init__(self,*args,**kw):super().__init__(*args,**kw);self.rows=[]
    def addSuccess(self,test):super().addSuccess(test);self.rows.append({'test':str(test),'status':'PASS'})
    def addFailure(self,test,error):super().addFailure(test,error);self.rows.append({'test':str(test),'status':'FAIL'})
    def addError(self,test,error):super().addError(test,error);self.rows.append({'test':str(test),'status':'ERROR'})
    def addSkip(self,test,reason):super().addSkip(test,reason);self.rows.append({'test':str(test),'status':'SKIPPED','reason':reason})


def main():
    global LEVELS,RADIATIVE
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--levels',type=Path);p.add_argument('--radiative',type=Path)
    p.add_argument('--output-dir',type=Path,required=True)
    args=p.parse_args();out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
    LEVELS=str(args.levels.resolve()) if args.levels else LEVELS
    RADIATIVE=str(args.radiative.resolve()) if args.radiative else RADIATIVE
    started=time.monotonic()
    suite=unittest.TestSuite(unittest.defaultTestLoader.loadTestsFromTestCase(c) for c in (LexicalProtocol,ActualSource))
    with(out/'unittest.log').open('w') as stream:r=unittest.TextTestRunner(stream=stream,verbosity=2,resultclass=Recorded).run(suite)
    counts={s:sum(x['status']==s for x in r.rows) for s in ('PASS','FAIL','ERROR','SKIPPED')}
    status='ERROR' if counts['ERROR'] else 'FAIL' if counts['FAIL'] else 'PASS' if counts['PASS'] else 'SKIPPED'
    value={'schema':'amps-native-chianti-source-tests-v1','status':status,'counts':counts,'checks':r.rows,
           'scope':'native lexical/parser/API checks; NOT physical cascade or reaction qualification',
           'source_revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
           'source_sha256':{str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in (Path(__file__).resolve(),Path(production.__file__).resolve())},
           'actual_inputs':[{'path':x,'sha256':hashlib.sha256(Path(x).read_bytes()).hexdigest() if Path(x).is_file() else None} for x in (LEVELS,RADIATIVE) if x],
           'elapsed_s':time.monotonic()-started,'runner':str(Path(__file__).resolve()),'log':str(out/'unittest.log')}
    (out/'result.json').write_text(json.dumps(value,indent=2)+'\n');print(json.dumps(counts),status)
    return {'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[status]


if __name__=='__main__':raise SystemExit(main())
