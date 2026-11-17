#!/usr/bin/env python3
"""Protocol/oracle checks; artificial bookkeeping is never observations.

Registered automatically by the global M2 development runner. Spawned sleepers
exercise real process-tree termination, including MPI-like independent sessions;
no science capability receives PASS from these infrastructure tests.
"""
from __future__ import annotations
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch
from types import SimpleNamespace

import campaign_budget
from campaign_budget import CampaignBudget
from m2_protocol import cases, grid_pairs, write_case, evaluate_cohorts

ROOT=Path(__file__).resolve().parents[2]

class Protocol(unittest.TestCase):
    def test_complete_design_and_seed_pairing(self):
        development=cases();cluster=cases('cluster')
        self.assertEqual(len(development),36);self.assertEqual(len(cluster),288)
        self.assertEqual(len({r['key'] for r in cluster}),288)
        self.assertEqual(len(grid_pairs()[1]),18)
        for r in development:self.assertIn(r['seed'],(2000003,2001012))
        self.assertTrue({r['seed'] for r in development}.isdisjoint({r['seed'] for r in cluster}))

    def test_derived_input_is_explicit_and_pins_controls(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'case.in';case=cases()[0]
            reference={'coarse_timestep_s':.5,'reference_na_weight_atoms':1e20,'reference_ar_weight_atoms':2e19}
            write_case(ROOT/'srcMoon/examples/lola_surface.in',path,case,reference=reference)
            text=path.read_text()
            for setting in ('m2_verification_case = source-fed','m2_verification_horizon_s = short',
                            'm2_verification_coarse_timestep_s = 0.5','random_seed = 2000003'):
                self.assertIn(setting,text)
            self.assertIn('surface_inventory_checkpoint_file = none',text)
            self.assertEqual(text.count('#end application'),1)

    def test_frozen_sites_override_template_and_preserve_shared_support(self):
        """The declared native mapping must control stock and normalization.

        These are input-contract fixtures only. Shared logical support is
        deliberately retained; it is never counted as independent observation.
        """
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'case.in'
            write_case(ROOT/'srcMoon/examples/m2_development.in',path,cases()[0],
                       elements=(12,12,34))
            self.assertIn('m2_verification_surface_elements = 12,12,34',path.read_text())

    def test_baseline_adapter_receipt_retains_physics_without_numerical_dispatch(self):
        from m2_protocol import without_verification_controls
        text=(ROOT/'srcMoon/examples/m2_development.in').read_text()
        physical=without_verification_controls(text)
        self.assertNotIn('m2_verification_case =',physical)
        self.assertIn('surface_geometry = lola',physical)
        self.assertIn('surface_inventory_initial_state = zero',physical)
        self.assertEqual(without_verification_controls(physical),physical)

    def test_explicit_first_case_references_are_retained(self):
        """Input controls must survive the first-case launcher transformation."""
        from integration_common import replace_moon_setting,moon_setting
        with tempfile.TemporaryDirectory() as tmp:
            template=Path(tmp)/'template.in';output=Path(tmp)/'case.in'
            text=(ROOT/'srcMoon/examples/m2_development.in').read_text()
            controls={'m2_verification_coarse_timestep_s':'0.5',
                      'm2_verification_na_weight_atoms':'1e21',
                      'm2_verification_ar_weight_atoms':'1e20'}
            for key,value in controls.items():text=replace_moon_setting(text,key,value,template)
            template.write_text(text)
            write_case(template,output,cases()[0])
            for key,value in controls.items():self.assertEqual(float(moon_setting(output,key)),float(value))

    def test_cluster_preserves_reviewed_physics_while_development_stays_bounded(self):
        """Contrasting input controls are parser fixtures, never lunar data.

        A formal launch must not obtain a passing identity by first changing
        the reviewed problem to the inexpensive I18 ownership fixture.
        """
        from integration_common import replace_moon_setting, moon_setting
        from m2_initial_state import physics_identity_sha256
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);template=root/'template.in'
            text=(ROOT/'srcMoon/examples/m2_development.in').read_text()
            declared={'mesh_near_surface_cell_size_m':'36000',
                      'mesh_polar_resolution_m':'18000',
                      'argon_source_rate_atoms_s':'4.1e20',
                      'thermal_time_step_s':'77',
                      'simulation_start_utc':'2009-01-25T01:02:03'}
            for key,value in declared.items():
                text=replace_moon_setting(text,key,value,template)
            template.write_text(text);formal=root/'formal/case.in'
            write_case(template,formal,cases('cluster')[0],tier='cluster')
            for key,value in declared.items():
                self.assertEqual(moon_setting(formal,key,preserve_case=True),value)
            self.assertEqual(physics_identity_sha256(formal),physics_identity_sha256(template))
            self.assertIn('m2_verification_record_initial_thermal_state = true',formal.read_text())
            development=root/'development.in'
            write_case(template,development,cases()[0])
            self.assertEqual(moon_setting(development,'mesh_near_surface_cell_size_m'),'120000')

    def test_cluster_missing_physical_prerequisite_skips_without_forcing_it(self):
        from integration_common import MissingPrerequisite,replace_moon_setting
        for key,value in [('argon_source_enabled','false'),
                          ('simulation_environment_mode','frozen'),
                          ('surface_inventory_initial_state','checkpoint')]:
            with self.subTest(key=key),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);template=root/'template.in';output=root/'case.in'
                original=(ROOT/'srcMoon/examples/m2_development.in').read_text()
                template.write_text(replace_moon_setting(original,key,value,template))
                before=template.read_bytes()
                with self.assertRaises(MissingPrerequisite):
                    write_case(template,output,cases('cluster')[0],tier='cluster')
                self.assertFalse(output.exists());self.assertEqual(template.read_bytes(),before)

    def test_formal_path_relocation_uses_included_assignment_owner(self):
        """Native include/continuation rules, not shell path expansion."""
        from m2_initial_state import stage_cluster_template,assignments,physics_identity
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);parts=root/'Parts';parts.mkdir()
            fragment=parts/'Physics Fragment.inc'
            fragment.write_text('lola_image_file = "../DEM/Case Data.IMG"\n'
                'sodium_surface_interaction_reference = "literal\\path"\n'
                'mesh_angular_resolution_m_by_sza = 0:100, \\\n    180:200 ! continuation boundary\n')
            template=root/'source.in'
            template.write_text('#begin application: moon\nargon_source_enabled = true\n'
                'simulation_environment_mode = evolving\nsurface_inventory_initial_state = zero\n'
                '#include "Parts/Physics Fragment.inc"\n#end application\n')
            destination=root/'elsewhere/case.in';stage_cluster_template(template,destination)
            after=assignments(destination)
            self.assertEqual(after['lola_image_file']['value'],str(root/'DEM/Case Data.IMG'))
            self.assertEqual(after['sodium_surface_interaction_reference']['value'],'literal\\path')
            self.assertEqual(after['mesh_angular_resolution_m_by_sza']['value'],'0:100, 180:200')
            self.assertEqual(physics_identity(template),physics_identity(destination))

    def test_physics_identity_ignores_only_named_controls_and_keeps_unknown_fields(self):
        from m2_initial_state import physics_identity_sha256
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);path=root/'case.in'
            def write(extra):
                path.write_text('#begin application: moon\nthermal_time_step_s = 77\n'+extra+
                                '\n#end application\n')
                return physics_identity_sha256(path)
            baseline=write('random_seed = 1\nm2_verification_dt_divisor = 1')
            self.assertEqual(baseline,write('random_seed = 99\nm2_verification_dt_divisor = 4\n'
                'surface_cea_file = anywhere\nmesh_output_tecplot = false\n'
                'thermal_initialization_mode = declared-uniform\n'
                'm2_polar_queries_file = geometry-only.csv\n'
                'm2_illumination_queries_file = other-queries.csv\n'
                'm2_formal_contract_file = separately-reviewed.json'))
            self.assertNotEqual(baseline,write('thermal_initialization_mode = periodic-spinup'))
            self.assertNotEqual(baseline,write('m2_verification_future_physical_switch = enabled'))
            self.assertNotEqual(write('future_surface_parameter = 1'),
                                write('future_surface_parameter = 2'))

    def test_cluster_invalid_native_boolean_is_error_instead_of_a_switch(self):
        from integration_common import replace_moon_setting
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);template=root/'template.in'
            text=(ROOT/'srcMoon/examples/m2_development.in').read_text()
            template.write_text(replace_moon_setting(text,'argon_source_enabled','on',template))
            # Moon's native ParseBoolean accepts true/false only. Broader
            # shell-style Boolean aliases cannot silently change its grammar.
            with self.assertRaises(ValueError):
                write_case(template,root/'case.in',cases('cluster')[0],tier='cluster')

    def test_development_keeps_explicit_initial_state_export_request(self):
        for value in ('true','false'):
            with self.subTest(value=value),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);template=root/'template.in'
                text=(ROOT/'srcMoon/examples/m2_development.in').read_text()
                text=text.replace('#end application',
                    'm2_verification_record_initial_thermal_state = '+value+'\n#end application')
                template.write_text(text);output=root/'case.in'
                write_case(template,output,cases()[0])
                self.assertIn('m2_verification_record_initial_thermal_state = '+value,
                              output.read_text())

    def test_cluster_unknown_verification_field_cannot_disappear_during_staging(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);template=root/'template.in'
            text=(ROOT/'srcMoon/examples/m2_development.in').read_text()
            text=text.replace('#end application',
                'm2_verification_future_physical_switch = enabled\n#end application')
            template.write_text(text)
            with self.assertRaisesRegex(ValueError,'changed frozen physics'):
                write_case(template,root/'case.in',cases('cluster')[0],tier='cluster')
            self.assertEqual(template.read_text(),text)

    def test_canonical_input_rejects_duplicate_cycle_and_dangling_continuation(self):
        from m2_initial_state import physics_identity
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);path=root/'case.in'
            for text in ['#begin application: moon\nx = 1\nx = 2\n#end application\n',
                         '#include case.in\n',
                         '#begin application: moon\nx = 1\\\n']:
                with self.subTest(text=text):
                    path.write_text(text)
                    with self.assertRaises(ValueError):physics_identity(path)

    def test_cluster_never_overwrites_existing_input_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'case.in';path.write_text('retained earlier receipt\n')
            before=path.read_bytes()
            with self.assertRaises(ValueError):
                write_case(ROOT/'srcMoon/examples/m2_development.in',path,
                           cases('cluster')[0],tier='cluster')
            self.assertEqual(path.read_bytes(),before)

    def test_cluster_manifest_binds_actual_physics_before_executor_launch(self):
        """Intercept the launch boundary; no native job or model data exists."""
        from datetime import datetime
        from m2_initial_state import physics_identity_sha256
        from m2_protocol import execute_design
        class BeforeNativeLaunch(Exception):
            pass
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            template=ROOT/'srcMoon/examples/m2_development.in'
            def inspect(jobs,*args,**kwargs):
                self.assertEqual(len(jobs),1)
                case=next(iter(jobs.values()));directory=case['directory']
                identity=json.loads((directory/'run_manifest.json').read_text())
                self.assertIsNotNone(datetime.fromisoformat(
                    identity['started_utc'].replace('Z','+00:00')).tzinfo)
                self.assertEqual(identity['physics_identity_sha256'],
                                 physics_identity_sha256(directory/'application.in'))
                self.assertEqual(identity['physics_identity_sha256'],
                                 physics_identity_sha256(template))
                self.assertIn('m2_verification_record_initial_thermal_state = true',
                              (directory/'application.in').read_text())
                raise BeforeNativeLaunch()
            with patch('integration_common.manifest',return_value={'scope':'protocol-only'}), \
                 patch('m1_execution.execute_levels',side_effect=inspect):
                with self.assertRaises(BeforeNativeLaunch):
                    execute_design(root/'not-an-executable',template,root/'cases',
                        tier='cluster',mpi_ranks=8,workers=1,remaining_seconds=lambda:1)

    def test_normalization_staging_distinguishes_formal_and_development_physics(self):
        """Stop at the launch boundary; never execute the placeholder binary."""
        from prepare_m2_normalization import prepare
        from integration_common import replace_moon_setting,moon_setting
        from m2_initial_state import physics_identity_sha256
        class BeforeNativeLaunch(Exception):
            pass
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);template=root/'template.in'
            text=(ROOT/'srcMoon/examples/m2_development.in').read_text()
            template.write_text(replace_moon_setting(text,'mesh_near_surface_cell_size_m','36000',template))
            for tier,expected in [('development','120000'),('cluster','36000')]:
                def inspect(command,*args,**kwargs):
                    application=Path(command[-1])
                    self.assertEqual(moon_setting(application,'mesh_near_surface_cell_size_m'),expected)
                    self.assertEqual(moon_setting(application,'m2_verification_horizon_s'),'2')
                    self.assertEqual(moon_setting(application,'m2_verification_profile'),'cluster')
                    self.assertEqual(physics_identity_sha256(application)==physics_identity_sha256(template),
                                     tier=='cluster')
                    raise BeforeNativeLaunch()
                with self.subTest(tier=tier), \
                     patch('prepare_m2_normalization.verify_retained_executable'), \
                     patch('prepare_m2_normalization.manifest',return_value={'scope':'protocol-only'}), \
                     patch.object(CampaignBudget,'run',side_effect=inspect):
                    with self.assertRaises(BeforeNativeLaunch):
                        prepare(root/'not-an-executable',template,root/tier,tier=tier,horizon_s=2)

    def test_normalization_ready_metadata_preserves_actual_weight_provenance(self):
        """A temporary bookkeeping receipt tests metadata, never validation.

        The fake executable and explicitly non-native manifest cannot satisfy
        the formal evidence validator. No frozen reference or observation is
        created, and no physical model or cluster job runs in this test.
        """
        from prepare_m2_normalization import prepare
        from m2_initial_state import physics_identity_sha256
        from integration_common import sha256,replace_moon_setting
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);executable=root/'protocol-placeholder';executable.write_text('not a native executable\n')
            template=root/'template.in';output=root/'normalization'
            text=(ROOT/'srcMoon/examples/m2_development.in').read_text()
            controls=[]
            for key,name in [('lola_image_file','NOT-A-DEM.control'),
                             ('lola_label_file','NOT-A-DEM-LABEL.control')]:
                file=root/name;file.write_text('PROTOCOL FIXTURE; NOT OBSERVATIONAL OR NATIVE DATA\n')
                text=replace_moon_setting(text,key,str(file),template);controls.append(file)
            template.write_text(replace_moon_setting(text,'thermal_layer_count','2',template))
            kernel=root/'NOT-A-SPICE-KERNEL.control';kernel.write_text('PROTOCOL FIXTURE ONLY\n');controls.append(kernel)
            model_identity={'scope':'protocol-only-not-native',
                            'external_data':{str(file):sha256(file) for file in controls}}
            def bookkeeping(command,*args,**kwargs):
                receipt=Path(command[command.index('--integration-output')+1])
                state=receipt.with_name('PROTOCOL-ONLY-initial-state.csv')
                state.write_text('logical_surface_element,depth_index,depth_m,kinetic_temperature_K\n'
                    '0,-1,0,250\n0,0,0.1,250\n0,1,0.2,250\n')
                receipt.write_text('supplier=PROTOCOL-FIXTURE-NOT-NATIVE\nduration_s=2\n'
                    'expected_external_na_atoms=100\nexpected_external_ar_atoms=7.6e20\n'
                    'initial_na_cohort_atoms_per_site=64\ninitial_ar_cohort_atoms_per_site=128\n'
                    'thermal_initialization_mode=declared-uniform\nthermal_initialization_et_s=0\n'
                    'initial_thermal_state_csv='+str(state)+'\ninitial_thermal_state_rows=3\n'
                    'initial_thermal_state_element_count=1\ninitial_thermal_state_layer_count=2\n'
                    'thermal_initialization_spice_kernel_count=1\n'
                    'thermal_initialization_spice_kernel_0_file='+str(kernel)+'\n')
                return {'returncode':0,'reason':None}
            with patch('prepare_m2_normalization.verify_retained_executable') as verifier, \
                 patch('prepare_m2_normalization.manifest',return_value=model_identity), \
                 patch.object(CampaignBudget,'run',side_effect=bookkeeping):
                result=prepare(executable,template,output,tier='cluster',horizon_s=2)
            verifier.assert_called_once_with(executable,require_current_generated_tree=False)
            self.assertEqual(result['status'],'READY')
            self.assertEqual(result['production_input_tier'],'cluster')
            self.assertEqual(result['physics_identity_sha256'],physics_identity_sha256(template))
            for name in ('actual_input','native_manifest','native_receipt'):
                self.assertEqual(sha256(Path(result[name]['path'])),result[name]['sha256'])
            provenance=result['native_weight_provenance']
            self.assertEqual(provenance['cohort_packets'],32)
            self.assertEqual(provenance['reference_weight_atoms'],{'Na':2,'Ar':4})
            candidate=result['candidate_initial_state']
            self.assertIs(candidate['qualified'],False)
            self.assertNotIn('review',candidate)
            self.assertEqual(candidate['state_shape']['rows'],3)
            self.assertEqual(candidate['physical_external_data']['ordered_actual_spice_kernels'],
                             [{'path':str(kernel),'sha256':sha256(kernel)}])

    def test_requested_normalization_state_must_exist_and_match_native_shape(self):
        """Missing, duplicated and inconsistent control CSVs remain ERROR."""
        from prepare_m2_normalization import prepare
        from integration_common import replace_moon_setting
        for fault in ('missing','duplicate','count'):
            with self.subTest(fault=fault),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);template=root/'template.in';executable=root/'not-native'
                executable.write_text('PROTOCOL ONLY\n')
                text=(ROOT/'srcMoon/examples/m2_development.in').read_text()
                template.write_text(replace_moon_setting(text,'thermal_layer_count','2',template))
                def bookkeeping(command,*args,**kwargs):
                    application=Path(command[-1])
                    self.assertIn('m2_verification_record_initial_thermal_state = true',application.read_text())
                    receipt=Path(command[command.index('--integration-output')+1]);state=receipt.parent/'control-state.csv'
                    if fault!='missing':
                        state.write_text('logical_surface_element,depth_index,depth_m,kinetic_temperature_K\n'
                            '0,-1,0,250\n0,0,0.1,250\n'+
                            ('0,0,0.2,250\n' if fault=='duplicate' else '0,1,0.2,250\n'))
                    receipt.write_text('duration_s=2\nexpected_external_na_atoms=100\n'
                        'expected_external_ar_atoms=7.6e20\ninitial_na_cohort_atoms_per_site=64\n'
                        'initial_ar_cohort_atoms_per_site=128\nthermal_initialization_mode=declared-uniform\n'
                        'thermal_initialization_et_s=0\ninitial_thermal_state_csv='+str(state)+'\n'
                        'initial_thermal_state_rows='+('4' if fault=='count' else '3')+'\n'
                        'initial_thermal_state_element_count=1\ninitial_thermal_state_layer_count=2\n')
                    return {'returncode':0,'reason':None}
                with patch('prepare_m2_normalization.verify_retained_executable'), \
                     patch('prepare_m2_normalization.manifest',return_value={'scope':'protocol-only-not-native'}), \
                     patch.object(CampaignBudget,'run',side_effect=bookkeeping):
                    result=prepare(executable,template,root/'normalization',horizon_s=2,record_initial_state=True)
                self.assertEqual(result['status'],'ERROR')
                self.assertNotIn('candidate_initial_state',result)

    def test_m2_rng_selection_does_not_change_canonical_or_sampling_input(self):
        from build_m2_qualification import derive_compile_input
        source=ROOT/'input/moon.input';original=source.read_text()
        self.assertIn('RND=MERSENNETWISTER',derive_compile_input(source,'fixture','m2'))
        self.assertNotIn('RND=MERSENNETWISTER',derive_compile_input(source,'fixture','sampling'))
        self.assertEqual(source.read_text(),original)

    def test_zero_and_finite_hazard_keep_censored_residents(self):
        # Exact independent first-order bookkeeping fixtures, not model output
        # and never evidence about real lunar temperature or abundance.
        import math
        for rate in (0.0,.5):
            with self.subTest(rate=rate),tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);H=2.;survival=math.exp(-rate*H)
                mean=H if rate==0 else -math.expm1(-rate*H)/rate
                row={'elapsed_s':H,'species':'Ar','cohort':0,'initial_atoms':100,
                     'first_survival':survival,'restricted_first_residence_s':mean,
                     'surface_fraction':survival,'gas_fraction':1-survival,
                     'escape_fraction':0,'numerical_loss_fraction':0,'tag_residual_atoms':0}
                hazard={'begin_elapsed_s':0,'dt_s':H,'species':'Ar','cohort':0,'total_release_hazard_s-1':rate}
                for name,value in (('cohorts.csv',row),('hazards.csv',hazard)):
                    with (root/name).open('w',newline='') as stream:
                        writer=csv.DictWriter(stream,fieldnames=list(value));writer.writeheader();writer.writerow(value)
                self.assertEqual(evaluate_cohorts(root/'cohorts.csv',root/'hazards.csv')['status'],'PASS')
                row['restricted_first_residence_s']=0
                with (root/'cohorts.csv').open('w',newline='') as stream:
                    writer=csv.DictWriter(stream,fieldnames=list(row));writer.writeheader();writer.writerow(row)
                self.assertEqual(evaluate_cohorts(root/'cohorts.csv',root/'hazards.csv')['status'],'FAIL')

    def test_equivalence_requires_precision_and_zero_support(self):
        from m2_equivalence import equivalence_interval
        self.assertTrue(equivalence_interval([.001+(-1)**i*.0001 for i in range(16)],.01)['pass'])
        self.assertFalse(equivalence_interval([(-1)**i*.05 for i in range(16)],.01)['pass'])
        self.assertFalse(equivalence_interval([0.0]*16,.01)['pass'])
        self.assertTrue(equivalence_interval([0.0]*16,.01,deterministic=True)['pass'])
        self.assertFalse(equivalence_interval([.02]*16,.01,deterministic=True)['pass'])
        with self.assertRaises(ValueError):equivalence_interval([0]*15,.01)

    def test_lost_cohort_artifacts_cannot_shrink_the_statistical_family(self):
        from m2_equivalence import load_primary
        # Zero-valued transfer rows are protocol fixtures, not observations or
        # a model qualification run. The fault under test is a lost CSV file.
        with tempfile.TemporaryDirectory() as tmp:
            directory=Path(tmp);(directory/'receipt.txt').write_text('steps=8\n')
            fields=['completed_steps','species','gas_atoms','surface_atoms',
                    'emission_atoms','surface_release_atoms','adsorption_atoms','escape_atoms']
            with (directory/'m2_budget.csv').open('w',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
                for species in ('Na','Ar'):
                    for step in range(9):
                        row=dict.fromkeys(fields,0)
                        row.update(completed_steps=step,species=species)
                        writer.writerow(row)
            with self.assertRaises(FileNotFoundError):load_primary(directory,require_cohorts=True)

    def test_cluster_does_not_promote_support_or_accept_false_receipt_identity(self):
        """Protocol fixtures exercise false-green handling, never lunar data."""
        sys.path.insert(0,str(Path(__file__).parent/'data_processing'))
        from run_m2_cluster import preflight
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);directory=root/'I18';directory.mkdir()
            (directory/'result.json').write_text(json.dumps({'status':'PASS','scope':'supporting-only'}))
            identity={'executable':{'sha256':'protocol-fixture-not-an-executable'},
                      'mpi_ranks':8,'formal_test_status':'SKIPPED'}
            (directory/'run_manifest.json').write_text(json.dumps(identity))
            (directory/'receipt.txt').write_text('test_id=I18\nmpi_ranks=8\n')
            args=SimpleNamespace(execute=False,qualification_root=root,amps=None,
                sites_file=None,forecast_file=None,normalizations_file=None,
                statistical_support_file=None,mpi_ranks=8,core_hours=300,
                input=ROOT/'srcMoon/examples/m2_development.in')
            with patch('audit_m2_packages.inspect',return_value={'records':[]}):
                gate=preflight(args)
                self.assertTrue(any('I18 evidence is explicitly supporting' in g for g in gate['gaps']))
                (directory/'result.json').write_text(json.dumps({'status':'PASS','scope':'formal-qualification'}))
                identity.pop('formal_test_status')
                (directory/'run_manifest.json').write_text(json.dumps(identity))
                (directory/'receipt.txt').write_text('test_id=I18\nmpi_ranks=1\n')
                gate=preflight(args)
                self.assertEqual(gate['status'],'ERROR')
                self.assertTrue(any('I18 native receipt identity/layout' in e for e in gate['errors']))

    def test_campaign_clock_is_not_reset_by_children(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);budget=CampaignBudget(1.2,.3)
            for i in range(2):
                result=budget.run([sys.executable,'-c','import time;time.sleep(.55)'],ROOT,
                    root/f'{i}.log',output=root/f'monitor-{i}')
            self.assertEqual(result['returncode'],124)
            self.assertLess(result['wall_seconds'],1.7)

    def test_allocation_includes_preflight_and_cpu_is_labelled_honestly(self):
        """Charge earlier allocation time; sampled CPU is not batch accounting."""
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);budget=CampaignBudget(2.,.2)
            budget.started-=.2  # explicit already-charged coordinator interval
            command=[sys.executable,'-c',
                'import time;end=time.monotonic()+.35\nwhile time.monotonic()<end: pass']
            result=budget.run(command,ROOT,root/'cpu.log',output=root/'monitor',allocated_cores=8)
            self.assertEqual(result['returncode'],0)
            self.assertGreaterEqual(result['allocated_core_hours'],8*.5/3600)
            self.assertGreater(result['sampled_cpu_core_hours_lower_bound'],0)
            self.assertIsNone(result['actual_scheduler_cpu_core_hours'])

    def test_detached_orphan_terminated_unrelated_child_preserved(self):
        sentinel=subprocess.Popen([sys.executable,'-c','import time;time.sleep(30)'])
        try:
            with tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);pidfile=root/'pid.txt'
                code=('import subprocess,sys,pathlib;'
                      'p=subprocess.Popen([sys.executable,"-c","import time;time.sleep(60)"],start_new_session=True);'
                      'pathlib.Path(sys.argv[1]).write_text(str(p.pid))')
                budget=CampaignBudget(1.1,.3)
                result=budget.run([sys.executable,'-c',code,str(pidfile)],ROOT,root/'child.log',output=root/'monitor')
                self.assertEqual(result['returncode'],124)
                pid=int(pidfile.read_text());stat=Path('/proc')/str(pid)/'stat'
                self.assertTrue(not stat.exists() or stat.read_text().rsplit(')',1)[1].split()[0]=='Z')
                self.assertIsNone(sentinel.poll())
        finally:
            sentinel.terminate();sentinel.wait()

    def test_reap_refresh_finds_fork_after_proc_enumeration(self):
        """Real fork handshake pins the non-atomic /proc/reap race exactly.

        The child waits until a real pre-fork table has been captured, then
        forks a new-session sleeper and exits. The first observation retains
        that old PID list but reads the launcher's real final zombie state,
        exactly as non-atomic /proc enumeration can. No process identity is
        invented. A successful reap must be followed by a NEW adoption scan;
        the old code returned zero and left the actual sleeper running.
        """
        real_table=campaign_budget.process_table
        sentinel=subprocess.Popen([sys.executable,'-c','import time;time.sleep(30)'])
        owned={}
        try:
            with tempfile.TemporaryDirectory() as tmp:
                root=Path(tmp);ready=root/'ready';release=root/'release';pidfile=root/'orphan'
                code=('import os,pathlib,subprocess,sys,time\n'
                      'ready,release,pidfile=map(pathlib.Path,sys.argv[1:])\n'
                      'temp=ready.with_suffix(".tmp");temp.write_text(str(os.getpid()));temp.replace(ready)\n'
                      'while not release.exists(): time.sleep(.001)\n'
                      'p=subprocess.Popen([sys.executable,"-c","import time;time.sleep(60)"],start_new_session=True)\n'
                      'temp=pidfile.with_suffix(".tmp");temp.write_text(str(p.pid));temp.replace(pidfile)\n')
                pinned=[False]
                def first_incomplete_snapshot():
                    snapshot=real_table()
                    if pinned[0] or not ready.is_file():return snapshot
                    launcher=int(ready.read_text())
                    if launcher not in snapshot:return snapshot
                    owned[launcher]=snapshot[launcher][2]
                    pinned[0]=True
                    release.write_text('fork now')
                    # Handshake is based on actual kernel state, not a sleep
                    # that assumes the fork or reparenting finished in time.
                    stop=time.monotonic()+.5
                    while time.monotonic()<stop:
                        current=real_table()
                        if pidfile.is_file():
                            orphan=int(pidfile.read_text())
                            if orphan in current:owned[orphan]=current[orphan][2]
                        if launcher in current and current[launcher][3]=='Z':
                            snapshot[launcher]=current[launcher]
                            return snapshot
                        time.sleep(.001)
                    self.fail('controlled launcher did not reach its real zombie handshake')
                budget=CampaignBudget(1.1,.3)
                try:
                    with patch.object(campaign_budget,'process_table',side_effect=first_incomplete_snapshot):
                        result=budget.run([sys.executable,'-c',code,str(ready),str(release),str(pidfile)],
                            ROOT,root/'child.log',output=root/'monitor')
                    self.assertTrue(pinned[0])
                    self.assertEqual(result['returncode'],124)
                    self.assertLess(result['wall_seconds'],2.)
                    pid=int(pidfile.read_text());current=real_table()
                    self.assertTrue(pid not in current or current[pid][3]=='Z')
                    self.assertIsNone(sentinel.poll())
                finally:
                    # Preserve the unrelated pre-existing child. Cleanup is
                    # restricted to the exact birth identities recorded from
                    # this fixture, even when an assertion deliberately fails.
                    for pid,birth in owned.items():
                        current=real_table()
                        if pid in current and current[pid][2]==birth and current[pid][3]!='Z':
                            try:os.kill(pid,15)
                            except ProcessLookupError:pass
                        try:os.waitpid(pid,os.WNOHANG)
                        except ChildProcessError:pass
        finally:
            sentinel.terminate();sentinel.wait()

if __name__=='__main__':unittest.main(buffer=True)
