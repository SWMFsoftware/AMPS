#!/usr/bin/env python3
"""Independent P02 arithmetic/source/provenance controls, never fake PIC runs."""
import copy
import csv
from decimal import Decimal,localcontext
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import verify_p02 as p

REPO=HERE.parents[2]
PROFILE=REPO/'test_output/srcMoon/m3-completion-20261119T091742Z/state_profiles/controls-final/profile/profile.json'
BUDGET=REPO/'srcMoon/test/M3_P02_NUMERICAL_BUDGET.json'

class Mathematics(unittest.TestCase):
    """Labelled analytical scalar/linear-algebra fixtures; no observations."""
    def test_distinct_mass_denominators(self):
        self.assertEqual(p.mass_errors(F(2),F(3)),(F(1,3),F(1,2)))
    def test_equal_masses_exact_zero(self):
        self.assertEqual(p.mass_errors(F(7,3),F(7,3)),(F(),F()))
    def test_mass_nonpositive_refused(self):
        for a,b in [(F(0),F(1)),(F(1),F(-1))]:
            with self.subTest(pair=(a,b)):
                with self.assertRaises(ValueError):p.mass_errors(a,b)
    def test_kinetic_exact_zero_boundary(self):
        self.assertTrue(p.kinetic_pass(F(),F(1),F()))
    def test_kinetic_subnormal_speed_not_false_zero(self):
        v2=F(1,2**22000)
        self.assertFalse(p.kinetic_pass(v2,F(299792458),F()))
        self.assertTrue(p.kinetic_pass(v2,F(299792458),F(1,10**100)))
    def test_kinetic_exact_nonzero_boundary(self):
        # beta=3/5,gamma=5/4: exact Kclass=.18mc² and Krel=.25mc².
        self.assertTrue(p.kinetic_pass(F(9,25),F(1),F(7,25)))
        self.assertFalse(p.kinetic_pass(F(9,25),F(1),F(7,25)-F(1,10**30)))
        self.assertTrue(p.kinetic_pass(F(9,25),F(1),F(7,25)+F(1,10**30)))
    def test_kinetic_independent_direct_decimal_energy_oracle(self):
        with localcontext()as ctx:
            ctx.prec=110
            for beta in ['0.000001','0.001','0.01','0.1','0.6','0.99']:
                with self.subTest(beta=beta):
                    b=Decimal(beta);gamma=1/(1-b*b).sqrt()
                    direct=((gamma-1)-b*b/2)/(gamma-1)
                    actual=Decimal(p.kinetic_diagnostic(F(b*b),F(1)))
                    self.assertLess(abs(actual-direct),Decimal('1e-75'))
                    margin=F(1,10**40)
                    central=F(direct)
                    self.assertTrue(p.kinetic_pass(F(b*b),F(1),central+margin))
                    self.assertFalse(p.kinetic_pass(F(b*b),F(1),central-margin))
    def test_invalid_kinetic_controls(self):
        for a,b,c in [(F(-1),F(1),F()),(F(),F(),F()),(F(),F(1),F(1))]:
            with self.assertRaises(ValueError):p.kinetic_pass(a,b,c)
    def test_luminal_and_superluminal_refusal(self):
        self.assertFalse(p.kinetic_pass(F(1),F(1),F(9,10)))
        self.assertFalse(p.kinetic_pass(F(2),F(1),F(9,10)))
    def test_exact_affine_inverse(self):
        r=list(map(F,[0,-1,0,1,0,0,0,0,1]))
        self.assertEqual(p.invert_affine(r,[F(2),F(3),F(4)],[F(-3),F(8),F(11)]),[F(5),F(5),F(7)])
    def test_inverse_not_rounded_transpose(self):
        epsilon=F(1,10**15);r=[1+epsilon,F(),F(),F(),F(1),F(),F(),F(),F(1)]
        self.assertEqual(p.invert_affine(r,[F(),F(),F()],[2*(1+epsilon),F(3),F(4)]),[F(2),F(3),F(4)])
    def test_native_caps_match_explicit_budget(self):
        caps={p.CAPS[0]:F(1,10**8),p.CAPS[1]:F(1,10**8),p.CAPS[2]:F(1,10**4),p.CAPS[3]:F(3000000)}
        controls={'maximum_classical_speed_ms':{'value':'3e6'}}
        for state in ('parent','ion','electron'):
            for field,value in [('inertial_mass','1e-8'),('acceleration','1e-8'),('relativistic_kinetic','1e-4')]:
                controls[state+'_maximum_relative_'+field+'_error']={'value':value}
        self.assertEqual(len(p.native_budget_checks(controls,caps)),10)
        self.assertTrue(all(a==b for _,_,a,b in p.native_budget_checks(controls,caps)))
        controls['electron_maximum_relative_inertial_mass_error']['value']='1e-7'
        self.assertEqual(sum(a!=b for _,_,a,b in p.native_budget_checks(controls,caps)),1)
    def test_missing_actual_native_cap_unavailable(self):
        with self.assertRaises(p.ledger.Unavailable):p.native_budget_checks({},dict(zip(p.CAPS,[F(1,10)]*4)))
    def test_actual_rounded_rotation_uses_unit_norm_scale(self):
        # Actual commissioning SPICE coefficients, source-derived geometry
        # fixture only; these are neither observations nor selected residuals.
        r=list(map(F,['0.560211314016633821211','-0.760052659393080420891','-0.329367938022605444282',
            '0.828347105332689293','0.515022425867632138718','0.220438140847663277544',
            '0.00208727926048927385841','-0.396322918594075035426','0.918108810252010965058']))
        limit=F(1,10**15);rows=p.rotation_checks(r,F(),limit)
        self.assertTrue(all(abs(residual)<=bound for _,_,residual,bound in rows))
        off=[row for row in rows if row[0]in {'proper-rotation-01','proper-rotation-02','proper-rotation-10','proper-rotation-12','proper-rotation-20','proper-rotation-21'}]
        self.assertEqual(len(off),6)
        self.assertTrue(all(bound==limit for _,_,_,bound in off))
        self.assertTrue(any(abs(residual)>limit*abs(actual)for _,actual,residual,_ in off))
    def test_rotation_offdiagonal_exact_boundary_and_violation(self):
        limit=F(1,10**15)
        for epsilon,expected in [(limit,True),(limit+F(1,10**30),False)]:
            r=[F(1),epsilon,F(),F(),F(1),F(),F(),F(),F(1)]
            rows={name:(actual,residual,bound)for name,actual,residual,bound in p.rotation_checks(r,F(),limit)}
            actual,residual,bound=rows['proper-rotation-01']
            self.assertEqual(abs(residual)<=bound,expected);self.assertEqual(bound,limit)
            self.assertEqual(rows['proper-rotation-00'][0],1) # columnGram, not rowGram
    def test_improper_rotation_still_fails_unit_determinant(self):
        rows=p.rotation_checks(list(map(F,[1,0,0,0,1,0,0,0,-1])),F(),F(1,10**15))
        self.assertTrue(all(abs(residual)<=bound for name,_,residual,bound in rows if name!='proper-rotation-determinant'))
        self.assertTrue(any(abs(residual)>bound for name,_,residual,bound in rows if name=='proper-rotation-determinant'))
    def test_singular_affine_refusal(self):
        with self.assertRaises(ValueError):p.invert_affine([F()]*9,[F()]*3,[F()]*3)

class NativeGrammar(unittest.TestCase):
    def test_native_include_and_child_scope(self):
        with tempfile.TemporaryDirectory()as d:
            r=Path(d);(r/'part.in').write_text('inertial_frame = J2000\n')
            (r/'native.in').write_text('#begin sodium_qualification: TEST\n#include "part.in"\ninertial_origin = SUN\n#begin source: ONE\npath = "original.txt"\n#end source\n#end sodium_qualification\n')
            values,docs=p.native_controls(r/'native.in')
            self.assertEqual(values['inertial_frame']['value'],'J2000');self.assertNotIn('path',values)
            self.assertEqual(set(docs),{r/'native.in',r/'part.in'})
    def test_duplicate_scalar_refused(self):
        with tempfile.TemporaryDirectory()as d:
            f=Path(d)/'bad.in';f.write_text('#begin sodium_qualification: X\nx = 1\nx = 2\n#end sodium_qualification\n')
            with self.assertRaises(ValueError):p.native_controls(f)
    def test_child_cannot_close_owner_in_include(self):
        with tempfile.TemporaryDirectory()as d:
            r=Path(d);(r/'child').write_text('#end sodium_qualification\n')
            (r/'main').write_text('#begin sodium_qualification: X\n#include child\n#end sodium_qualification\n')
            with self.assertRaises(ValueError):p.native_controls(r/'main')
    def test_include_cycle_refused(self):
        with tempfile.TemporaryDirectory()as d:
            r=Path(d);(r/'main').write_text('#include main\n')
            with self.assertRaises(ValueError):p.native_controls(r/'main')

class Sources(unittest.TestCase):
    requires_authoritative_source=True
    @classmethod
    def setUpClass(cls):
        if not PROFILE.is_file()or not BUDGET.is_file():raise unittest.SkipTest('actual source/budget unavailable')
    def test_real_source_physical_reference_masses(self):
        states,c=p.source_masses(PROFILE,{'profile_sha256':p.ledger.digest(PROFILE),'parent':'23Na:ground','ion':'23Na+:ground','electron':'electron:free'}, {})
        self.assertEqual(c,F(299792458));self.assertEqual(states['parent']['physical_mass'],states['parent']['mass'])
        self.assertGreater(states['ion']['physical_mass'],states['ion']['mass'])
        self.assertEqual(states['electron']['physical_mass'],states['electron']['mass'])
    def test_real_source_mass_caps_and_declared_speed(self):
        ids={};caps=p.load_budget({'path':str(BUDGET),'sha256':p.ledger.digest(BUDGET)},HERE,ids)
        states,c=p.source_masses(PROFILE,{'profile_sha256':p.ledger.digest(PROFILE),'parent':'23Na:ground','ion':'23Na+:ground','electron':'electron:free'},ids)
        for state in states.values():
            ie,ae=p.mass_errors(state['mass'],state['physical_mass']);self.assertLessEqual(ie,caps[p.CAPS[0]]);self.assertLessEqual(ae,caps[p.CAPS[1]])
        self.assertTrue(p.kinetic_pass(caps[p.CAPS[3]]**2,c,caps[p.CAPS[2]]))
    def test_source_hash_corruption_error(self):
        with self.assertRaises(ValueError):p.source_masses(PROFILE,{'profile_sha256':'0'*64}, {})
    def test_missing_ion_remains_unavailable(self):
        with self.assertRaises(p.ledger.Unavailable):p.source_masses(PROFILE,{'profile_sha256':p.ledger.digest(PROFILE),'parent':'23Na:ground','ion':'23Na+:unresolved','electron':'electron:free'}, {})
    def test_budget_hash_change_refused(self):
        with self.assertRaises(ValueError):p.load_budget({'path':str(BUDGET),'sha256':'0'*64},HERE,{})

class Admission(unittest.TestCase):
    def test_missing_binding_skipped(self):
        with tempfile.TemporaryDirectory()as d:self.assertEqual(p.verify(Path(d)/'missing')['status'],'SKIPPED')
    def test_malformed_binding_error(self):
        with tempfile.TemporaryDirectory()as d:
            path=Path(d)/'binding';path.write_text('{}\n');self.assertEqual(p.verify(path)['status'],'ERROR')
    def test_no_owner_source_only_skipped(self):
        if not PROFILE.exists():self.skipTest('real profile unavailable')
        with tempfile.TemporaryDirectory()as d:
            path=Path(d)/'binding';path.write_text(json.dumps({'schema':p.BINDING_SCHEMA,'profile':{'path':str(PROFILE),'sha256':p.ledger.digest(PROFILE)},'budget':{'path':str(BUDGET),'sha256':p.ledger.digest(BUDGET)},'owners':[]}))
            self.assertEqual(p.verify(path)['status'],'SKIPPED')
    def test_wrong_owner_file_bytes_error(self):
        if not PROFILE.exists():self.skipTest('real profile unavailable')
        with tempfile.TemporaryDirectory()as d:
            path=Path(d)/'binding';owner=Path(d)/'owner';owner.write_text('not a completed owner')
            path.write_text(json.dumps({'schema':p.BINDING_SCHEMA,'profile':{'path':str(PROFILE),'sha256':p.ledger.digest(PROFILE)},'budget':{'path':str(BUDGET),'sha256':p.ledger.digest(BUDGET)},'owners':[{'owner_manifest':{'path':str(owner),'sha256':'0'*64},'owner_result':{'path':str(owner),'sha256':'0'*64},'mpi_ranks':8}]}))
            self.assertEqual(p.verify(path)['status'],'ERROR')

if __name__=='__main__':unittest.main()
