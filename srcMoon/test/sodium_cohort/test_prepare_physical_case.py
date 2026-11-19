#!/usr/bin/env python3
"""Administrative parser/provenance controls; fixtures are NOT physical data.

The deliberately fictional one-eV table and energy records below test only
strict source identity, controlled serialization and refusal behavior. They are
never admitted to AMPS, compared to observations, or presented as a Na profile.
"""
import contextlib
import csv
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import shutil

sys.path.insert(0, str(Path(__file__).resolve().parent))
import prepare_physical_case as P


class Preparation(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='moon-administrative-preparation-')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.profile_path = self.root/'profile.json'
        self.crosswalk = self.root/'crosswalk.csv'
        self.case_path = self.root/'controls.json'
        self.application = self.root/'application.in'
        self.output = self.root/'prepared'
        self.reference = self.root/'reference.txt'
        self.reference.write_text('FICTIONAL ADMINISTRATIVE REFERENCE, NOT PHYSICAL EVIDENCE\n')
        self.ion_source = self.root/'ion.txt'
        self.ion_source.write_text('FICTIONAL FORMATION SOURCE, NOT PHYSICAL EVIDENCE\n')
        self.reference_id = 'REF@sha256:'+P.sha(self.reference)+'#declared-zero'
        self.ion_id = 'ION@sha256:'+P.sha(self.ion_source)+'#formation'
        self.thresholds = self.root/'thresholds.csv'
        self.thresholds.write_text('reaction_id,threshold_eV\nPHIDRATES:Na:00:na_e,1.0\n')
        self.rates = self.root/'rates.csv'
        self.rates.write_text('reaction_id,rate_1au_s-1\nPHIDRATES:Na:00:na_e,1.0\n')
        self.evidence = self.root/'administrative-policy.txt'
        self.evidence.write_text('ARITHMETIC FIXTURE ONLY; NOT MODEL AUTHORIZATION\n')
        self.row = dict(zip(P.CROSSWALK_COLUMNS, ('PHIDRATES:Na:00:na_e',
            '23Na:ground','23Na+:ground','electron:free',self.reference_id,
            'source-effective-photon','1.0',self.ion_id)))
        self.write_crosswalk()
        self.meta = {**self.row,'schema':'srcMoon-M3-Na-source-energy-crosswalk-v1',
                     'source_files':[self.pin(self.thresholds),self.pin(self.rates)]}
        self.meta_path = self.root/'crosswalk.json'
        self.write_json(self.meta_path,self.meta)
        fields = ['state_id','nuclei','charge','formation_J','formation_source_id',
                  'formation_source_record','excitation_J','excitation_source_id',
                  'excitation_source_record','status']
        self.states = self.root/'states.csv'
        with self.states.open('w',newline='') as stream:
            writer=csv.writer(stream);writer.writerow(fields)
            writer.writerow(['23Na:ground','23Na:1','0','0','REF','declared-zero','0','REF','declared-zero','available'])
            writer.writerow(['23Na+:ground','23Na:1','1','1','ION','formation','0','REF','declared-zero','available'])
            writer.writerow(['electron:free','','-1','0','REF','declared-zero','0','REF','declared-zero','available'])
        self.isotopes=self.root/'isotopes.csv';self.isotopes.write_text('fixture\nNOT_SCIENCE\n')
        self.constants=self.root/'constants.csv';self.constants.write_text('fixture\nNOT_SCIENCE\n')
        self.profile={'schema':'srcMoon-M3-state-profile-v1','physical_activation':False,
            'isotopes':self.pin(self.isotopes),'constants':self.pin(self.constants),
            'states':self.pin(self.states),'sources':[
                {'id':'REF',**self.pin(self.reference),'role':'reference-definition'},
                {'id':'ION',**self.pin(self.ion_source),'role':'energy'}],
            'atomic_extension':{'crosswalk':self.pin(self.meta_path),
                                'crosswalk_csv':self.pin(self.crosswalk)}}
        self.write_json(self.profile_path,self.profile)
        self.fragment=self.root/'fragment.in';self.fragment.write_text('! EMPTY ADMINISTRATIVE FRAGMENT\n')
        self.application.write_text('#begin application: moon\n#include "fragment.in"\n'
            'particle_tracing_mover = force-cut-surface-second-order\n'
            'tracing_clock_mode = common-time\nsimulation_environment_mode = frozen\n'
            'm2_environment = disabled\n#end application\n')
        controls={k:'0' for k in P.NUMERIC_KEYS}
        controls.update({k:'1' for k in P.UINT_KEYS})
        controls.update(case_kind='homogeneous-fixed-position',inertial_frame='J2000',
            inertial_origin='SUN',source_reference_distance_identity='ADMINISTRATIVE-ONE-METRE',
            earth_occultation='off',model_identity='ADMINISTRATIVE-NOT-PHYSICS',
            solar_endpoint='quiet',transport_policy_identity='ADMINISTRATIVE-NOT-AUTHORIZATION',
            reference_weight_mode='production',
            transport_policy_evidence_file=str(self.evidence),
            transport_policy_evidence_sha256=P.sha(self.evidence),parent_compiled_species='0',
            ion_compiled_species='1',physical_step_s='1',physical_horizon_s='1',
            represented_number_per_parent='1',source_reference_distance_m='1',
            maximum_classical_speed_ms='1',exposure_root_time_tolerance_s='0.001',
            exposure_absolute_depth_tolerance='0.001')
        self.case={'schema':P.SCHEMA,'cohort_id':'administrative-fixture','controls':controls}
        self.write_json(self.case_path,self.case)

    @staticmethod
    def pin(path):return {'path':str(path),'sha256':P.sha(path)}
    @staticmethod
    def write_json(path,value):path.write_text(json.dumps(value,indent=2)+'\n')
    def write_case(self):self.write_json(self.case_path,self.case)
    def write_crosswalk(self):
        with self.crosswalk.open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=P.CROSSWALK_COLUMNS)
            writer.writeheader();writer.writerow(self.row)
    def refresh_metadata(self):
        self.write_json(self.meta_path,self.meta)
        self.profile['atomic_extension']={'crosswalk':self.pin(self.meta_path),
                                          'crosswalk_csv':self.pin(self.crosswalk)}
        self.write_json(self.profile_path,self.profile)
    def prepare(self):return P.prepare(self.profile_path,self.crosswalk,self.application,self.case_path,self.output)
    def refused(self,exception=ValueError):
        with self.assertRaises(exception):self.prepare()
        self.assertFalse(self.output.exists())

    def test_complete_case_preserves_explicit_controls_and_projection(self):
        result=self.prepare();self.assertEqual(result['status'],'PASS')
        text=(self.output/'native.in').read_text()
        for key,value in self.case['controls'].items():
            quoted='"'+value+'"' if key.endswith('_file') else value
            self.assertIn('  '+key+' = '+quoted+'\n',text)
        self.assertIn('photon_energy_policy = source-effective-photon',text)
        self.assertNotIn('source_threshold_eV =',text)  # CSV stays the authority.
        self.assertNotIn('qualified =',text)
        projection=self.root/'independent-projection.in'
        P.project(self.profile_path,projection)
        self.assertEqual(projection.read_bytes(),(self.output/'profile.native.in').read_bytes())
        manifest=json.loads((self.output/'run_manifest.json').read_text())
        self.assertEqual([r['path'] for r in manifest['application_documents']],
                         [str(self.application),str(self.fragment)])
        self.assertFalse(manifest['physical_activation'])
        for record in manifest['outputs']:
            self.assertEqual(record['sha256'],P.sha(record['path']))

    def test_actual_native_parser_preserves_raw_text_values_and_quoted_paths(self):
        """Owner.Text reads RAW values; only Owner.Path removes outer quotes.

        Compile the production generic parser, without PIC/AMPS configuration.
        This catches accidental JSON quoting and inline-# truncation using the
        exact native Assignment.value returned to the existing owner.
        """
        compiler=shutil.which('g++')
        if not compiler:self.skipTest('g++ unavailable for maintained native grammar check')
        self.prepare()
        header=Path(P.__file__).resolve().parents[3]/'src/general/input/amps_input.h'
        source=self.root/'raw_parser.cpp'
        source.write_text('#include "'+str(header)+'"\n#include <iostream>\n'
            'int main(int argc,char**argv){AMPS::Input::Document d;auto s=AMPS::Input::ParseFile(argv[1],&d);'
            'if(!s.ok()){std::cerr<<s.message;return 2;}'
            'if(d.root.blocks.size()!=1)return 3;const auto&b=d.root.blocks[0];'
            'if(b.kind!="sodium_qualification"||b.name!="administrative-fixture")return 4;'
            'for(int i=2;i+1<argc;i+=2){auto*a=b.Find(argv[i]);'
            'if(!a||a->value!=argv[i+1]){std::cerr<<argv[i]<<" raw value mismatch";return 5;}}return 0;}\n')
        executable=self.root/'raw_parser'
        compile_run=subprocess.run([compiler,'-std=c++17','-Wall','-Wextra','-Werror','-pedantic',
            str(source),'-o',str(executable)],text=True,capture_output=True,timeout=30)
        self.assertEqual(compile_run.returncode,0,compile_run.stderr)
        command=[str(executable),str(self.output/'native.in')]
        for key in ('parent_state_id','ion_state_id','electron_state_id',
                    'energy_reference_identity','photon_energy_policy'):
            command.extend((key,self.row[key]))
        command.extend(('solar_endpoint','quiet','application_input','"'+str(self.application)+'"'))
        parsed=subprocess.run(command,text=True,capture_output=True,timeout=10)
        self.assertEqual(parsed.returncode,0,parsed.stderr)

    def test_missing_control_is_skipped(self):
        del self.case['controls']['energy_absolute_bound'];self.write_case()
        self.refused(P.MissingConfiguration)
    def test_template_nulls_cannot_publish(self):
        self.write_json(self.case_path,P.template());self.refused(P.MissingConfiguration)
    def test_unknown_control_is_error(self):
        self.case['controls']['qualified']=True;self.write_case();self.refused()
    def test_float_loses_lexical_precision_and_is_refused(self):
        self.case['controls']['energy_absolute_bound']=0.0;self.write_case();self.refused()
    def test_negative_bound(self):
        self.case['controls']['rate_absolute_bound']='-1';self.write_case();self.refused()
    def test_model_error_domain_is_not_widened(self):
        self.case['controls']['ion_maximum_relative_inertial_mass_error']='1'
        self.write_case();self.refused()
    def test_zero_root_tolerance_is_not_runnable(self):
        self.case['controls']['exposure_root_time_tolerance_s']='0'
        self.write_case();self.refused()
    def test_no_positive_exposure_bound(self):
        self.case['controls']['exposure_absolute_depth_tolerance']='0'
        self.write_case();self.refused()
    def test_text_whitespace_is_not_silently_changed(self):
        self.case['controls']['inertial_origin']=' SUN'
        self.write_case();self.refused()
    def test_noncanonical_uint(self):
        self.case['controls']['campaign_seed']='01';self.write_case();self.refused()
    def test_horizon_never_truncates(self):
        self.case['controls'].update(physical_steps='3',physical_step_s='0.1',physical_horizon_s='0.3')
        self.write_case();self.refused()
    def test_exact_dyadic_partition(self):
        self.case['controls'].update(physical_steps='4',physical_step_s='0.25',physical_horizon_s='1')
        self.write_case();self.assertEqual(self.prepare()['status'],'PASS')
    def test_nonzero_homogeneous_velocity(self):
        self.case['controls']['source_velocity_lso_ms_2']='1';self.write_case();self.refused()
    def test_explicit_unequal_weights_are_preserved(self):
        self.case['controls'].update(reference_weight_mode='explicit',parent_reference_weight='1',ion_reference_weight='8')
        self.write_case();self.prepare()
        text=(self.output/'native.in').read_text()
        self.assertIn('parent_reference_weight = 1',text);self.assertIn('ion_reference_weight = 8',text)
    def test_explicit_weight_mode_requires_both_weights(self):
        self.case['controls'].update(reference_weight_mode='explicit',parent_reference_weight='1')
        self.write_case();self.refused(P.MissingConfiguration)
    def test_production_weight_mode_forbids_explicit_override(self):
        self.case['controls']['parent_reference_weight']='1';self.write_case();self.refused()
    def test_nonpositive_reference_weight_is_refused(self):
        self.case['controls'].update(reference_weight_mode='explicit',parent_reference_weight='1',ion_reference_weight='0')
        self.write_case();self.refused()
    def test_malformed_profile_hash(self):
        self.profile['states']['sha256']='0'*64;self.write_json(self.profile_path,self.profile);self.refused()
    def test_source_mutation(self):
        self.ion_source.write_text('MODIFIED');self.refused()
    def test_crosswalk_mutation(self):
        self.crosswalk.write_text('CORRUPTED');self.refused()
    def test_crosswalk_wrong_state(self):
        self.row['ion_state_id']='unknown';self.meta.update(self.row)
        self.write_crosswalk();self.refresh_metadata();self.refused(P.MissingConfiguration)
    def test_wrong_formation_source_identity(self):
        self.row['physical_threshold_source_identity']=self.reference_id;self.meta.update(self.row)
        self.write_crosswalk();self.refresh_metadata();self.refused()
    def test_wrong_native_threshold(self):
        self.row['source_threshold_eV']='1.0001';self.meta.update(self.row)
        self.write_crosswalk();self.refresh_metadata();self.refused()
    def test_missing_atomic_crosswalk(self):
        del self.profile['atomic_extension'];self.write_json(self.profile_path,self.profile)
        self.refused(P.MissingConfiguration)
    def test_policy_payload_hash_mismatch(self):
        self.evidence.write_text('CHANGED');self.refused()
    def test_newline_cannot_inject_native_assignment(self):
        self.case['controls']['model_identity']='fixture\nqualified = true';self.write_case();self.refused()
    def test_comment_cannot_cut_off_native_scalar(self):
        self.case['controls']['inertial_frame']='J2000!';self.write_case();self.refused()
    def test_output_nonoverwrite(self):
        self.prepare();before={p.name:p.read_bytes() for p in self.output.iterdir()}
        with self.assertRaises(ValueError):self.prepare()
        self.assertEqual(before,{p.name:p.read_bytes() for p in self.output.iterdir()})
    def test_broken_output_symlink_preserved(self):
        self.output.symlink_to(self.root/'missing')
        with self.assertRaises(ValueError):self.prepare()
        self.assertTrue(self.output.is_symlink())
    def test_unknown_application_temporal_coupling(self):
        self.application.write_text(self.application.read_text().replace('frozen','evolving'))
        self.refused()
    def test_duplicate_json_key(self):
        self.case_path.write_text('{"schema":"x","schema":"y"}')
        self.refused()
    def test_cli_missing_configuration_is_77(self):
        self.write_json(self.case_path,P.template())
        result=subprocess.run([sys.executable,'-B',str(Path(P.__file__)),
            '--profile',str(self.profile_path),'--crosswalk',str(self.crosswalk),
            '--application-input',str(self.application),'--case-json',str(self.case_path),
            '--output-dir',str(self.output)],text=True,capture_output=True)
        self.assertEqual(result.returncode,77,result.stderr)
        self.assertEqual(json.loads(result.stdout)['status'],'SKIPPED')
        self.assertFalse(self.output.exists())
    def test_cli_malformed_case_is_2(self):
        self.case['controls']['qualified']=True;self.write_case()
        result=subprocess.run([sys.executable,'-B',str(Path(P.__file__)),
            '--profile',str(self.profile_path),'--crosswalk',str(self.crosswalk),
            '--application-input',str(self.application),'--case-json',str(self.case_path),
            '--output-dir',str(self.output)],text=True,capture_output=True)
        self.assertEqual(result.returncode,2,result.stderr)
        self.assertEqual(json.loads(result.stdout)['status'],'ERROR')
        self.assertFalse(self.output.exists())
    def test_cli_template_is_explicitly_incomplete(self):
        target=self.root/'template.json'
        result=subprocess.run([sys.executable,'-B',str(Path(P.__file__)),
            '--write-template',str(target)],text=True,capture_output=True)
        self.assertEqual(result.returncode,77,result.stderr)
        self.assertEqual(json.loads(result.stdout)['status'],'SKIPPED')
        template=json.loads(target.read_text())
        self.assertEqual(set(template['controls']),P.REQUIRED)
        self.assertTrue(all(value is None for value in template['controls'].values()))


if __name__=='__main__':unittest.main(verbosity=2)
