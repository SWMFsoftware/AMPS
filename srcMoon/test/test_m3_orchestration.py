#!/usr/bin/env python3
"""M3 orchestration protocol controls only: no builds, MPI or observations.

All output trees are temporary. Mocked launchers describe CONTROL FIXTURES,
never actual scientific approvals, predecessor qualification or physical data.
The tests exercise production orchestration/status/restoration logic itself;
they do not replace the native coefficient/callback verification campaign.
"""
from __future__ import annotations
import contextlib
import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import run_m3_development as development
import run_m3_rate_callback as callback

STATES=('PASS','FAIL','SKIPPED','ERROR')


def summary(status='SKIPPED',local='PASS'):
    records=[{'id':'U16','status':local,'scope':'CONTROL FIXTURE ONLY'}]
    records += [{'id':'I'+str(n),'status':'SKIPPED','scope':'formal-chemistry-not-qualified'} for n in range(21,25)]
    support=[{'check':'CONTROL FIXTURE ONLY','status':'PASS'}]
    return {'schema':'srcMoon-initial-M3-development-v1','status':status,
        'results':records,'supporting_checks':support,'error':None,
        'formal_M3_status':'SKIPPED / NOT QUALIFIED','observational_scoring':False,
        'counts':{s:sum(x['status']==s for x in records) for s in STATES},
        'supporting_counts':{s:sum(x['status']==s for x in support) for s in STATES}}


class WorkerProtocol(unittest.TestCase):
    def execute(self,*,baseline='PASS',local='PASS',callback_status='PASS',
                callback_error=None,restore_error=None,corrupt=None,
                baseline_error_count=0):
        temporary=tempfile.TemporaryDirectory();self.addCleanup(temporary.cleanup)
        output=Path(temporary.name)/'campaign'
        def launch(command,**kwargs):
            destination=Path(command[command.index('--output-dir')+1]);destination.mkdir(parents=True)
            if corrupt!='missing':
                baseline_record={'m0_gate_status':baseline}
                if baseline_error_count:
                    # The existing M0 owner can retain a coarse FAIL gate
                    # while its structured records identify build/restore
                    # ERRORs. The dependent owner must honor that stronger
                    # infrastructure evidence, not call it a physics miss.
                    baseline_record.update(counts={'ERROR':baseline_error_count},
                        results=[{'id':'CONTROL_BUILD_ONLY','status':'ERROR',
                                  'scope':'CONTROL FIXTURE; no model execution'}])
                text='BROKEN JSON CONTROL FIXTURE' if corrupt=='json' else json.dumps(baseline_record)
                (destination/'baseline_summary.json').write_text(text)
            (destination/'u-series').mkdir()
            rows=[] if local is None else [{'id':'U16','status':local,'scope':'CONTROL FIXTURE ONLY'}]
            (destination/'u-series/run_all_summary.json').write_text(json.dumps({'results':rows}))
            return SimpleNamespace(returncode=0 if baseline=='PASS' else 1)
        native=patch.object(development,'run_callback',side_effect=callback_error,
                            return_value={'status':callback_status,'scope':'CONTROL FIXTURE ONLY'})
        restore=patch.object(development,'clean_configure_build',side_effect=restore_error)
        with patch.object(development.subprocess,'run',side_effect=launch),native as call,restore as restored:
            with contextlib.redirect_stdout(io.StringIO()):
                code=development.worker(output,8,100,'CONTROL_TOOLKIT','CONTROL_KERNELS',10)
        return code,json.loads((output/'development_summary.json').read_text()),call,restored

    def formal_skips(self,value):
        formal={row['id']:row['status'] for row in value['results'] if row['id'].startswith('I')}
        self.assertEqual(formal,{'I'+str(n):'SKIPPED' for n in range(21,25)})
        self.assertEqual(value['formal_M3_status'],'SKIPPED / NOT QUALIFIED')
        self.assertFalse(value['observational_scoring'])

    def test_baseline_fail_stops_dependent_callback(self):
        code,value,native,restore=self.execute(baseline='FAIL')
        self.assertEqual((code,value['status']),(1,'FAIL'));native.assert_not_called();self.formal_skips(value)

    def test_baseline_gate_fail_with_build_error_is_infrastructure_error(self):
        code,value,native,restore=self.execute(baseline='FAIL',baseline_error_count=2)
        self.assertEqual((code,value['status']),(2,'ERROR'))
        native.assert_not_called();self.formal_skips(value)
        self.assertEqual(value['supporting_checks'][0]['status'],'ERROR')

    def test_local_scientific_fail_is_not_hidden_as_skip(self):
        code,value,native,restore=self.execute(local='FAIL')
        self.assertEqual((code,value['status']),(1,'FAIL'));native.assert_not_called();self.formal_skips(value)

    def test_local_skip_stops_callback_and_keeps_formal_skips(self):
        code,value,native,restore=self.execute(local='SKIPPED')
        self.assertEqual((code,value['status']),(0,'SKIPPED'));native.assert_not_called();self.formal_skips(value)

    def test_callback_scientific_fail_preserved_and_restored(self):
        code,value,native,restore=self.execute(callback_status='FAIL')
        self.assertEqual((code,value['status']),(1,'FAIL'));restore.assert_called_once();self.formal_skips(value)
        self.assertEqual(restore.call_args[0][0],'moon')

    def test_callback_infrastructure_error_restores(self):
        code,value,native,restore=self.execute(callback_error=RuntimeError('CONTROL child runtime error'))
        self.assertEqual((code,value['status']),(2,'ERROR'));restore.assert_called_once();self.formal_skips(value)

    def test_restore_error_overrides_scientific_fail(self):
        code,value,native,restore=self.execute(callback_status='FAIL',restore_error=RuntimeError('CONTROL restore failure'))
        self.assertEqual((code,value['status']),(2,'ERROR'));self.assertIn('restore failed',value['error'])

    def test_missing_baseline_output_is_error(self):
        code,value,native,restore=self.execute(corrupt='missing')
        self.assertEqual((code,value['status']),(2,'ERROR'));native.assert_not_called()

    def test_corrupt_baseline_output_is_error(self):
        code,value,native,restore=self.execute(corrupt='json')
        self.assertEqual((code,value['status']),(2,'ERROR'));native.assert_not_called()

    def test_missing_local_record_is_error(self):
        code,value,native,restore=self.execute(local=None)
        self.assertEqual((code,value['status']),(2,'ERROR'));native.assert_not_called()

    def test_support_success_never_closes_formal_chemistry(self):
        code,value,native,restore=self.execute()
        self.assertEqual((code,value['status']),(0,'SKIPPED'));restore.assert_called_once();self.formal_skips(value)


class SupervisorProtocol(unittest.TestCase):
    def execute(self,child=None,code=0,raw=None,reason=None):
        temporary=tempfile.TemporaryDirectory();self.addCleanup(temporary.cleanup)
        output=Path(temporary.name)/'campaign'
        def launch(command,cwd,log,*,output):
            if raw is not None:(output/'development_summary.json').write_text(raw)
            elif child is not None:(output/'development_summary.json').write_text(json.dumps(child))
            return {'returncode':code,'reason':reason,'scope':'CONTROL FIXTURE ONLY'}
        with patch.object(sys,'argv',['run_m3_development.py','--output-dir',str(output),
                '--wall-seconds','100','--cleanup-reserve-seconds','10']),\
                patch.object(development.CampaignBudget,'run',side_effect=launch):
            returned=development.main()
        return returned,json.loads((output/'development_summary.json').read_text()),output

    def test_scientific_fail_exit_one_is_preserved(self):
        code,value,output=self.execute(summary('FAIL','FAIL'),1)
        self.assertEqual((code,value['status']),(1,'FAIL'))

    def test_formal_skip_exit_zero_remains_skip(self):
        code,value,output=self.execute(summary(),0)
        self.assertEqual((code,value['status']),(0,'SKIPPED'))

    def test_missing_summary_even_after_zero_exit_is_error(self):
        code,value,output=self.execute()
        self.assertEqual((code,value['status']),(2,'ERROR'))

    def test_corrupt_summary_preserved_before_error(self):
        text='CORRUPT CONTROL FIXTURE'
        code,value,output=self.execute(raw=text)
        self.assertEqual((code,value['status']),(2,'ERROR'))
        self.assertEqual((output/'invalid-worker-summary.json').read_text(),text)

    def test_exit_status_disagreement_is_error(self):
        code,value,output=self.execute(summary('FAIL','FAIL'),0)
        self.assertEqual((code,value['status']),(2,'ERROR'))

    def test_counts_disagreement_is_error(self):
        child=summary();child['counts']['PASS']=999
        code,value,output=self.execute(child)
        self.assertEqual((code,value['status']),(2,'ERROR'))

    def test_record_failure_cannot_hide_under_coherent_skip_counts(self):
        code,value,output=self.execute(summary('SKIPPED','FAIL'),0)
        self.assertEqual((code,value['status']),(2,'ERROR'))

    def test_budget_termination_cannot_be_success(self):
        code,value,output=self.execute(summary(),124,reason='CONTROL budget expired')
        self.assertEqual((code,value['status']),(2,'ERROR'))


class CallbackProtocol(unittest.TestCase):
    def execute(self,case,*,native_code=0,configure_error=False):
        temporary=tempfile.TemporaryDirectory();self.addCleanup(temporary.cleanup)
        root=Path(temporary.name);repo=root/'isolated-control-repo';(repo/'input').mkdir(parents=True)
        canonical='SpeciesList = Na ! CONTROL ONLY\n#block ../moon.input.pl\n#include moon.input.spicekernels\n'
        for name,text in [('moon.input',canonical),('moon.input.pl','CONTROL ONLY'),('moon.input.spicekernels','CONTROL ONLY'),('keep.in','UNRELATED')]:
            (repo/'input'/name).write_text(text)
        (repo/'keep-root.in').write_text('UNRELATED')
        coefficients=root/'control-only-coefficients.dat'
        coefficients.write_text('CONTROL FIXTURE; not authoritative atomic data')
        # Actual native block/assignment grammar, with intentionally minimal
        # control fields. No complete physical input/parser approval is claimed
        # and the native executable remains mocked in this infrastructure test.
        application=root/'input.in';application.write_text(
            '#begin application: moon\nelectron_impact_rate_model = maxwellian-voronov\n'
            'electron_impact_coefficients_file = '+str(coefficients)+'\n#end application\n')
        output=root/'result'
        derived_species=[]
        old=str(Path(callback.__file__).resolve().parent/'common/production_kernel_probe.cpp')
        def build_command(build,executable):return ['CONTROL-CXX',old,'-o',str(executable)]
        def launch(command,cwd,log,*,output):
            label=Path(log).stem
            Path(log).write_text('CONTROL FIXTURE ONLY; no process launched\n')
            if label=='configure':
                # Inspect the actual derived compile-input bytes before the
                # mocked compiler sees them. Display labels such as Na+ are
                # not the native species registry's NA_PLUS identifier.
                key=next(value for value in command if value.startswith('-application='))
                staged=repo/'input'/(key.split('=',1)[1]+'.input')
                for line in staged.read_text().splitlines():
                    if line.startswith('SpeciesList'):
                        derived_species.extend(line.split('=',1)[1].split('!',1)[0].strip().split(','))
                if configure_error:return {'returncode':2,'reason':'CONTROL configure failed'}
                for relative in ('pic/pic.h','pic/picGlobal.dfn','pic/picSpeciesMacro.dfn','general/rnd.h'):
                    path=repo/'build'/relative;path.parent.mkdir(parents=True,exist_ok=True);path.write_text('CONTROL ABI BYTES')
            if label=='compile-probe':Path(command[-1]).write_text('CONTROL ONLY; not executable model')
            if label=='stdout':
                Path(log).write_text('M3 native rate callback rank=0 status='+('FAIL' if native_code==1 else 'PASS')+'\n')
                for rank in range(2):
                    if case=='missing-rank' and rank==1:continue
                    rows=[]
                    for index,name in enumerate(('He','Na','Ne','Ar')):
                        row=dict(rank=rank,atom=name,neutral_index=2*index,expected_ion_index=2*index+1,
                                 actual_ion_index=2*index+1,callback_rate_s='0',scalar_rate_s='0',
                                 # Arbitrary positive SI metadata fixtures,
                                 # not measured or asserted atomic masses.
                                 neutral_mass_kg=(index+1)*1e-26,ion_mass_kg=(index+1)*1e-26,
                                 neutral_charge_c=0,ion_charge_c=1.602176634e-19)
                        if case=='missing-field':row.pop('actual_ion_index')
                        if case=='bad-rank':row['rank']='INVALID CONTROL BYTE'
                        if case=='scientific-fail':row['actual_ion_index']=99
                        if case=='zero-ion-mass':row['ion_mass_kg']=0
                        rows.append(row)
                    with (output/('callback.rank'+str(rank)+'.csv')).open('w',newline='') as stream:
                        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
                return {'returncode':native_code,'reason':None}
            return {'returncode':0,'reason':None}
        def digest(path):
            path=Path(path)
            return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() and str(path).startswith(str(root)) else 'CONTROL-NOT-SCIENTIFIC-HASH'
        with patch.object(callback,'REPO',repo),patch.object(callback,'_probe_compile_command',side_effect=build_command),\
                patch.object(callback.CampaignBudget,'run',side_effect=launch),\
                patch.object(callback,'manifest',return_value={'external_data':{},'scope':'CONTROL FIXTURE ONLY'}),\
                patch.object(callback,'production_source_identity',return_value={'scope':'CONTROL FIXTURE ONLY'}),\
                patch.object(callback,'sha256',side_effect=digest),\
                patch('integration_common.moon_setting',side_effect=lambda path,key:'0' if key=='electron_density_m3' else '10'):
            value=callback.run(output,application,mpi_ranks=2,wall_seconds=100)
        self.assertEqual((repo/'input/moon.input').read_text(),canonical)
        self.assertEqual((repo/'input/keep.in').read_text(),'UNRELATED')
        self.assertEqual((repo/'keep-root.in').read_text(),'UNRELATED')
        self.assertEqual(list((repo/'input').glob('moon_m3_rates_*')),[])
        self.assertEqual(json.loads((output/'result.json').read_text())['status'],value['status'])
        # Errors before MPI launch still need the stream paths advertised by
        # the global runner, rather than a pointer to a nonexistent test log.
        self.assertTrue((output/'stdout.log').is_file())
        self.assertTrue((output/'stderr.log').is_file())
        value['_control_derived_species_list']=derived_species
        return value

    def test_complete_control_rows_do_not_advertise_formal_qualification(self):
        value=self.execute('complete')
        self.assertEqual(value['status'],'PASS');self.assertEqual(value['formal_I23_status'],'SKIPPED')

    def test_completed_native_scientific_fail_stays_fail(self):
        value=self.execute('scientific-fail',native_code=1)
        self.assertEqual(value['status'],'FAIL')

    def test_derived_species_recipe_uses_native_ion_identifiers(self):
        value=self.execute('complete')
        self.assertEqual(value['_control_derived_species_list'],
            ['Na','NA_PLUS','He','HE_PLUS','Ne','NE_PLUS','Ar','AR_PLUS'])

    def test_zero_ion_mass_is_independently_scored_as_scientific_fail(self):
        # A nominally successful mocked launcher must not hide invalid mass
        # metadata. This specifically exercises Python's independent scorer,
        # rather than relying on a native exit1/FAIL fallback diagnostic.
        self.assertEqual(self.execute('zero-ion-mass')['status'],'FAIL')

    def test_missing_rank_output_is_error(self):
        self.assertEqual(self.execute('missing-rank')['status'],'ERROR')

    def test_missing_required_csv_field_is_error(self):
        self.assertEqual(self.execute('missing-field')['status'],'ERROR')

    def test_corrupt_numeric_output_is_error(self):
        self.assertEqual(self.execute('bad-rank')['status'],'ERROR')

    def test_build_failure_is_error_and_only_owned_temporary_inputs_removed(self):
        self.assertEqual(self.execute('complete',configure_error=True)['status'],'ERROR')


class RateManifestProtocol(unittest.TestCase):
    def fixture(self,included=False,absolute=False):
        temporary=tempfile.TemporaryDirectory();self.addCleanup(temporary.cleanup)
        root=Path(temporary.name)
        directory=root/'included-owner' if included else root
        (directory/'data').mkdir(parents=True)
        coefficients=directory/'data/cfit.dat'
        coefficients.write_bytes(b'CONTROL FIXTURE ONLY; not a scientific coefficient table\n')
        literal=str(coefficients) if absolute else 'data/cfit.dat'
        application=root/'application.in'
        text='#begin application: moon\nelectron_impact_rate_model = maxwellian-voronov\n'
        if included:
            owner=directory/'rates.in'
            owner.write_text('electron_impact_coefficients_file = '+literal+'\n')
            text+='#include "included-owner/rates.in"\n'
            # This decoy catches resolution against the root input, instead
            # of the actual included assignment's native source location.
            (root/'data').mkdir();(root/'data/cfit.dat').write_text('WRONG CONTROL DECOY')
        else:
            owner=application
            text+='electron_impact_coefficients_file = '+literal+'\n'
        application.write_text(text+'#end application\n')
        generic={'uses_spice_kernels':True,'configuration':'SPICE-DE421',
            'species_fixture':'Na','external_data':{'NAIF/control-only.bsp':'CONTROL-INSTALLED-HASH'},
            'inputs':{}}
        with patch.object(callback,'manifest',return_value=generic),\
                patch.object(callback,'production_source_identity',return_value={'scope':'CONTROL FIXTURE ONLY'}):
            result=callback.rate_manifest(root/'not-a-native-executable',['CONTROL COMMAND'],application,8)
        return result,coefficients,application,owner

    def test_absolute_canonical_style_path_hashes_selected_bytes(self):
        result,coefficients,application,owner=self.fixture(absolute=True)
        self.assertEqual(result['external_data'],{str(coefficients):hashlib.sha256(coefficients.read_bytes()).hexdigest()})
        self.assertEqual(result['coefficient_input']['resolved_file'],str(coefficients))

    def test_relative_coefficient_path_uses_input_owner(self):
        result,coefficients,application,owner=self.fixture()
        self.assertEqual(result['coefficient_input']['resolved_file'],str(coefficients))
        self.assertEqual(result['coefficient_input']['assignment_file'],str(owner))

    def test_included_coefficient_path_uses_included_owner_not_decoy(self):
        result,coefficients,application,owner=self.fixture(included=True)
        self.assertEqual(result['external_data'],{str(coefficients):hashlib.sha256(coefficients.read_bytes()).hexdigest()})
        self.assertEqual(result['coefficient_input']['assignment_file'],str(owner))
        for path in (application,owner):
            self.assertEqual(result['inputs'][str(path)],hashlib.sha256(path.read_bytes()).hexdigest())

    def test_unused_spice_inventory_cannot_be_reported_as_consumed_data(self):
        result,coefficients,application,owner=self.fixture()
        self.assertFalse(result['uses_spice_kernels'])
        self.assertNotEqual(result['configuration'],'SPICE-DE421')
        self.assertEqual(result['configuration_preflight_inventory']['files'],
            {'NAIF/control-only.bsp':'CONTROL-INSTALLED-HASH'})
        self.assertIn('not consumed',result['configuration_preflight_inventory']['role'])
        self.assertNotIn('NAIF/control-only.bsp',result['external_data'])
        self.assertEqual(result['species_fixture'],'Na,NA_PLUS,He,HE_PLUS,Ne,NE_PLUS,Ar,AR_PLUS')
        self.assertEqual(result['formal_I23_status'],'SKIPPED')
        self.assertFalse(result['observational_scoring'])


if __name__=='__main__':unittest.main()
