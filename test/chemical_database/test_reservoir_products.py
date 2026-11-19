"""Production shared-kernel moment tests, NOT observations or Na/He validation.

The C++ probe uses declared nominal masses/energies to isolate conservation and
source/domain guards. Python's independent scalar/vector expressions supply
analytic oracles. No rate, daughter PDF, experimental error bar or source
parameter is fitted here. All native commands/logs/source identities are kept.
"""
from pathlib import Path
import argparse
from datetime import datetime,timezone
from decimal import Decimal,localcontext
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[2]
OUTPUT=None


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class ReservoirProducts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        global OUTPUT
        if shutil.which('g++')is None:raise unittest.SkipTest('C++17 compiler unavailable')
        if OUTPUT is None:
            base=ROOT/'test_output/chemical_database';base.mkdir(parents=True,exist_ok=True)
            OUTPUT=Path(tempfile.mkdtemp(prefix='reservoir-products-',dir=base))
        else:OUTPUT.mkdir(parents=True,exist_ok=False)
        cls.output=OUTPUT;cls.exe=OUTPUT/'reservoir-probe'
        sources=[ROOT/'test/chemical_database/reservoir_products_probe.cpp',ROOT/'src/models/chemical_database/Runtime.cpp',ROOT/'src/models/chemical_database/DataIO.cpp']
        dependencies=sources+[ROOT/'src/models/chemical_database'/name for name in ('ReservoirProducts.h','ReducedModels.h','Runtime.h','DataIO.h')]+[ROOT/'src/general/input/amps_input.h']
        cls.before={str(path):digest(path)for path in dependencies}
        command=['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror','-pedantic']+[str(p)for p in sources]+['-o',str(cls.exe)]
        r=subprocess.run(command,capture_output=True,text=True,timeout=120)
        (OUTPUT/'compile.stdout.log').write_text(r.stdout);(OUTPUT/'compile.stderr.log').write_text(r.stderr)
        (OUTPUT/'build_manifest.json').write_text(json.dumps({'command':command,'returncode':r.returncode,'source_hashes':cls.before,'compiler':subprocess.check_output(['g++','--version'],text=True).splitlines()[0],'executable_sha256':digest(cls.exe)if cls.exe.exists()else None,'scope':'Analytical shared production-kernel fixture only; NOT observations'},indent=2)+'\n')
        if r.returncode:raise RuntimeError('shared reservoir probe compilation failed; see '+str(OUTPUT))
    @classmethod
    def tearDownClass(cls):
        changed=[name for name,expected in cls.before.items()if digest(name)!=expected]
        if changed:raise RuntimeError('tested production source changed during execution '+str(changed))
    def probe(self,mode):
        """Retain actual executable output; a process failure is an ERROR."""
        r=subprocess.run([str(self.exe),mode],capture_output=True,text=True,timeout=10)
        stem=self.output/(self._testMethodName+'-'+mode)
        stem.with_suffix('.stdout.log').write_text(r.stdout);stem.with_suffix('.stderr.log').write_text(r.stderr)
        if r.returncode:raise RuntimeError('native probe failed '+str(r.returncode))
        rows={}
        for line in r.stdout.splitlines():
            key,_,value=line.partition(' ')
            if key in rows:raise ValueError('duplicate native diagnostic '+key)
            rows[key]=value
        return rows
    def number(self,r,key):return float(r[key])
    def vector(self,r,key):return [float(x)for x in r[key].split()]
    def close(self,a,b):
        # This is an explicit numerical oracle budget, not physical uncertainty:
        # at most32 elementary binary64 operations in these independent oracles,
        # with factor4 allowance for sqrt/decimal-to-binary comparison paths.
        python_oracle_roundoff_factor=128
        self.assertAlmostEqual(a,b,delta=python_oracle_roundoff_factor*sys.float_info.epsilon*max(1,abs(a),abs(b)))
    def accepted(self,mode):
        r=self.probe(mode);self.assertEqual(r['code'],'Evaluated',r);self.assertEqual(r['sentinel'],'0');return r
    def refused(self,mode,code):
        r=self.probe(mode);self.assertEqual(r['code'],code,r);self.assertEqual(r['sentinel'],'1','failed composition changed output');return r
    def test_constant_cross_section_collision_mean_is_two_kT(self):
        r=self.probe('mean');self.assertEqual(r['code'],'Evaluated');self.assertEqual(r['sentinel'],'1');self.close(self.number(r,'mean'),4)
    def test_thresholded_constant_sigma_mean_matches_independent_integral(self):
        # Integral E^2 exp(-E/theta) / integral E exp(-E/theta), E>=I.
        # This is an analytical test source, not an assigned lunar parameter.
        r=self.probe('threshold-mean');I=3;theta=2
        self.assertEqual(r['code'],'Evaluated');self.close(self.number(r,'mean'),(I*I+2*I*theta+2*theta*theta)/(I+theta))
    def test_EI_aggregate_counts_and_explicit_charge(self):
        r=self.accepted('ei');self.assertEqual(r['incoming-number'],'1');self.assertEqual(r['outgoing-number'],'2');self.close(self.number(r,'net-charge')/-1.602176634e-19,1)
    def test_EI_incoming_and_outgoing_moments_are_separate(self):
        r=self.accepted('ei');self.assertEqual(self.vector(r,'incoming-P'),[1,2,0]);self.assertEqual(self.vector(r,'outgoing-P'),[2,4,0]);self.close(self.number(r,'incoming-K'),8.3);self.close(self.number(r,'outgoing-K'),7.8)
    def test_EI_signed_bath_and_state_Q_conserve(self):
        r=self.accepted('ei');self.close(self.number(r,'Q'),-3);self.close(self.number(r,'net-E'),-.5);[self.close(value,0)for value in self.vector(r,'p-residual')];self.close(self.number(r,'energy-residual'),0)
    def test_EI_recoil_bound_has_independent_moment_oracle(self):
        r=self.accepted('ei');impulse=math.sqrt(2*1*5.8)+math.sqrt(4*1*2.8)
        self.close(self.number(r,'speed-error'),impulse/10);self.close(self.number(r,'energy-error'),math.sqrt(5)*impulse+impulse**2/20)
    def test_actual_velocity_time_and_weight_are_preserved(self):
        r=self.accepted('ei');self.assertEqual(self.vector(r,'v'),[1,2,0]);self.close(self.number(r,'epoch'),100.375);self.close(self.number(r,'elapsed'),.375);self.close(self.number(r,'remaining'),.625);self.close(self.number(r,'weight'),7.5);self.close(self.number(r,'weighted-E'),-3.75)
    def test_rotation_covariance_of_EI_moments(self):
        a=self.accepted('ei');b=self.accepted('rotated');self.assertEqual(self.vector(b,'net-P'),[0,1,2]);self.close(self.number(a,'speed-error'),self.number(b,'speed-error'));self.close(self.number(a,'energy-error'),self.number(b,'energy-error'))
    def test_EI_boost_transforms_both_bath_sides(self):
        r=self.accepted('boosted');self.assertEqual(self.vector(r,'v'),[4,2,0]);self.assertEqual(self.vector(r,'incoming-P'),[4,2,0]);self.assertEqual(self.vector(r,'outgoing-P'),[8,4,0]);self.close(self.number(r,'net-E'),7);self.close(self.number(r,'weighted-E'),52.5)
    def test_partial_CEX_has_independent_closed_form_velocity(self):
        r=self.accepted('cex')
        with localcontext()as c:
            c.prec=70;vx=Decimal(21)/13-(Decimal(12120).sqrt()/13)/10
        self.close(self.vector(r,'v')[0],float(vx));self.close(self.vector(r,'v')[1],2);self.close(self.number(r,'Q'),2);self.close(self.number(r,'incoming-relative-K'),176/13)
    def test_CEX_consumed_plasma_and_fast_neutral_are_distinct(self):
        r=self.accepted('cex');self.assertEqual(r['incoming-id'],'B+');self.assertEqual(r['outgoing-id'],'B');self.assertEqual(r['incoming-number'],'1');self.assertEqual(r['outgoing-number'],'1');self.close(self.number(r,'net-charge')/-1.602176634e-19,1)
    def test_CEX_zero_declared_angle_has_zero_direction_error(self):
        r=self.accepted('cex');self.close(self.number(r,'speed-error'),0);self.close(self.number(r,'energy-error'),0)
    def test_CEX_explicit_all_angle_bound_is_not_an_angle_measurement(self):
        r=self.accepted('cex-angle');self.close(self.number(r,'speed-error'),2*math.sqrt(12120)/130)
    def test_CEX_rotation_and_boost(self):
        a=self.accepted('cex');b=self.accepted('cex-rotated');c=self.accepted('cex-boosted');v=self.vector(a,'v');self.close(self.vector(b,'v')[1],v[0]);self.close(self.vector(b,'v')[2],v[1]);self.close(self.vector(c,'v')[0],v[0]+3)
    def test_CEX_all_momentum_and_energy_close(self):
        r=self.accepted('cex');self.close(self.number(r,'energy-residual'),0);[self.close(value,0)for value in self.vector(r,'p-residual')];self.close(self.number(r,'weighted-E'),7.5*self.number(r,'net-E'))
    def test_nonfinite_supplied_derivative_is_invalid_not_unavailable(self):self.refused('infinite-derivative','InvalidInput')
    def test_missing_profile_refuses_transactionally(self):self.refused('missing','Unavailable')
    def test_corrupt_hash_refuses(self):self.refused('bad-hash','InvalidInput')
    def test_mismatched_energy_reference_refuses(self):self.refused('bad-reference','InvalidInput')
    def test_missing_state_Q_refuses(self):self.refused('missing-Q','Unavailable')
    def test_unsupported_Maxwellian_drift_refuses(self):self.refused('drift','Unavailable')
    def test_actual_nonzero_Maxwellian_drift_refuses(self):self.refused('drift-vector','OutsideDomain')
    def test_missing_collision_rotation_assumption_refuses(self):self.refused('anisotropic','Unavailable')
    def test_mean_below_threshold_refuses(self):self.refused('threshold','OutsideDomain')
    def test_source_threshold_state_crosswalk_refuses(self):self.refused('source-mismatch','OutsideDomain')
    def test_temperature_outside_derivative_domain_refuses(self):self.refused('temperature','OutsideDomain')
    def test_velocity_model_budget_cannot_be_relaxed(self):self.refused('budget','OutsideDomain')
    def test_energy_model_budget_refuses(self):self.refused('energy-budget','OutsideDomain')
    def test_physical_mass_budget_refuses(self):self.refused('mass-budget','OutsideDomain')
    def test_integer_nuclei_mismatch_refuses(self):self.refused('nuclei','InvalidInput')
    def test_integer_charge_mismatch_refuses(self):self.refused('charge','InvalidInput')
    def test_reference_material_mass_mismatch_refuses(self):self.refused('mass','InvalidInput')
    def test_actual_event_clock_mismatch_refuses(self):self.refused('clock','InvalidInput')
    def test_zero_represented_number_refuses(self):self.refused('N0','InvalidInput')
    def test_represented_moment_overflow_refuses(self):self.refused('Noverflow','NumericalFailure')
    def test_represented_moment_underflow_refuses(self):self.refused('Nunderflow','NumericalFailure')
    def test_CLI_unknown_argument_is_not_silently_ignored(self):
        target=self.output/'unknown-argument-output'
        r=subprocess.run([sys.executable,'-B',str(Path(__file__).resolve()),'--output-dir',str(target),'--unknown-source-selection'],capture_output=True,text=True,timeout=10)
        self.assertEqual(r.returncode,2);self.assertFalse(target.exists())
    def test_CLI_existing_output_is_refused_without_overwriting(self):
        existing=self.output/'existing-evidence';existing.mkdir();sentinel=existing/'result.json';sentinel.write_text('EXISTING USER EVIDENCE\n')
        r=subprocess.run([sys.executable,'-B',str(Path(__file__).resolve()),'--output-dir',str(existing)],capture_output=True,text=True,timeout=10)
        self.assertEqual(r.returncode,2);self.assertEqual(sentinel.read_text(),'EXISTING USER EVIDENCE\n');self.assertFalse((existing/'build_manifest.json').exists())
    def test_zero_CEX_relative_direction_refuses(self):self.refused('cex-zero','OutsideDomain')
    def test_negative_CEX_available_energy_refuses(self):self.refused('cex-impossible','OutsideDomain')
    def test_actual_CEX_source_energy_domain_refuses(self):self.refused('cex-domain','OutsideDomain')
    def test_missing_CEX_angle_evidence_refuses(self):self.refused('cex-missing-angle','Unavailable')
    def test_missing_fast_neutral_state_refuses(self):self.refused('cex-missing-neutral','Unavailable')
    def test_CEX_nuclei_mismatch_refuses(self):self.refused('cex-bad-nuclei','OutsideDomain')
    def test_CEX_requires_actual_nucleus_carriers(self):self.refused('cex-no-nuclei','OutsideDomain')
    def test_elastic_nontransfer_cannot_be_labelled_CEX(self):self.refused('cex-not-transfer','OutsideDomain')
    def test_integer_nuclei_accumulation_overflow_refuses(self):self.refused('cex-nuclei-overflow','NumericalFailure')
    def test_CEX_angle_outside_physical_range_refuses(self):self.refused('cex-bad-angle','InvalidInput')


def main():
    global OUTPUT
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args();OUTPUT=args.output_dir.resolve()
    if OUTPUT.exists():parser.error('output directory already exists; choose a new evidence directory')
    started=datetime.now(timezone.utc)
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(ReservoirProducts)
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    status='ERROR'if result.errors else'FAIL'if result.failures else'SKIPPED'if result.skipped else'PASS'
    record={'schema':'amps-reservoir-products-verification-v1','status':status,'counts':{'executed':result.testsRun,'TOTAL':max(result.testsRun,len(result.failures)+len(result.errors)+len(result.skipped)),'PASS':max(0,result.testsRun-len(result.failures)-len(result.errors)-len(result.skipped)),'FAIL':len(result.failures),'ERROR':len(result.errors),'SKIPPED':len(result.skipped)},'scope':'Production shared mathematical/moment software only; analytical nominal masses, NOT measured Na/He or physical source qualification','started_utc':started.isoformat(),'finished_utc':datetime.now(timezone.utc).isoformat(),'test_source_sha256':digest(__file__),'failures':[str(a)+'\n'+b for a,b in result.failures],'errors':[str(a)+'\n'+b for a,b in result.errors]}
    OUTPUT.mkdir(parents=True,exist_ok=True);(OUTPUT/'result.json').write_text(json.dumps(record,indent=2)+'\n')
    return {'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[status]

if __name__=='__main__':raise SystemExit(main())
