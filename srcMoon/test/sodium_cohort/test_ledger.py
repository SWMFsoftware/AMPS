#!/usr/bin/env python3
"""Independent-verifier software controls using labelled toy bookkeeping.

The toy profile/events are deterministic arithmetic fixtures, not sodium
physics, observations, or a qualification/reference release.
"""
import copy,csv,json,sys,tempfile,unittest,subprocess
from pathlib import Path
from fractions import Fraction as F
from verify_ledger import verify,digest,binary,STATE_HEADER,load_photon_source


def fixture(root):
    root=Path(root);raw=root/'analytical-source.txt';raw.write_text('ANALYTICAL SOFTWARE FIXTURE ONLY; NOT PHYSICAL DATA\n')
    sources=[{'id':role,'path':str(raw),'sha256':digest(raw),'role':role} for role in ('isotope-mass','constants','energy')]
    def table(name,head,rows):
        path=root/(name+'.csv')
        with path.open('w',newline='') as stream:w=csv.writer(stream);w.writerow(head);w.writerows(rows)
        return {'path':str(path),'sha256':digest(path)}
    isotopes=table('isotopes',['isotope','neutral_mass_kg','source_id','source_record'],[['TOY','2','isotope-mass','fixture:1']])
    constants=table('constants',['name','value','unit','source_id','source_record'],[[name,value,unit,'constants','fixture:1'] for name,value,unit in [('electron mass','1','kg'),('elementary charge','1.602176634e-19','C'),('speed of light in vacuum','299792458','m s^-1')]])
    c=F(299792458);formation=c-F(1,2)
    states=table('states',STATE_HEADER,[
        ['toy-neutral','toy neutral','ground','TOY:1','0','0','energy','fixture:zero','0','energy','fixture:zero','2','isotope-mass','fixture:1','available',''],
        ['toy-ion','toy ion','ground','TOY:1','1',str(float(formation)),'energy','fixture:formation','0','energy','fixture:zero','1','energy','fixture:mass','available',''],
        ['toy-electron','toy free electron','free-electron','','-1','0','energy','fixture:zero','0','energy','fixture:zero','1','constants','fixture:1','available','']])
    profile=root/'profile.json';profile.write_text(json.dumps({'schema':'srcMoon-M3-state-profile-v1','scope':'analytical-software-fixture','sources':sources,'isotopes':isotopes,'constants':constants,'states':states}))
    mapping={'parent':'toy-neutral','ion':'toy-ion','electron':'toy-electron','profile_sha256':digest(profile),'ion_compiled_species':1,'ion_catalog_species':1,'event_frame':'TOY-EVENT'}
    bounds={'schema':'srcMoon-M3-event-ledger-bounds-v1','scope':'analytical-software-fixture','event_frame':'TOY-EVENT',**{k:{'absolute':'0','relative':'0'} for k in ('mass_kg','momentum_kg_m_s','energy_j','number','rotation')}}
    # A 90 degree rotation with velocity offset (2,0,0) independently exercises
    # the signed reservoir coordinate transformation. Incoming photon p=(1,0,0),
    # E=c; ion inertial v=(1,0,0), outgoing electron v=0, and source formation
    # c-1/2 close both kinetic+internal energy and additive-reference mass.
    values={'parent_number':3,'ion_number':3,'parent_reference_weight':2,'parent_individual_correction':F(3,2),'ion_reference_weight':4,'ion_individual_correction':F(3,4),'parent_internal_j':0,'ion_internal_j':formation,'photon_inertial_j':c,'bath_event_j':2-c,'bath_mass_kg':1,'electron_state_j':0,'source_threshold_j':formation,'source_mean_excess_j':F(1,2),'source_q_j':-formation,'ledger_parent_mass_kg':6,'ledger_ion_mass_kg':3,'ledger_bath_mass_kg':3,'ledger_tracked_de_j':3*(c-2),'ledger_bath_de_j':3*(2-c),'tracked_charge_change_e':1,'bath_charge_change_e':-1,'exposure_error_bound':0,'clock_storage_error_bound':0,'actual_ion_reference_weight':4,'actual_ion_individual_correction':F(3,4),'actual_ion_internal_j':formation,'actual_ion_compiled_index':1,'actual_ion_catalog_species':1,'actual_ion_child_count':1,'event_epoch_et_s':8,'actual_ion_epoch_et_s':8,'planned_interval_start_s':0,'planned_interval_duration_s':10,'planned_interval_elapsed_s':8,'actual_ion_interval_start_s':0,'actual_ion_interval_duration_s':10,'actual_ion_interval_elapsed_s':8}
    for name,vector in {'parent_event_v':[2,0,0],'ion_event_v':[2,1,0],'electron_event_v':[2,0,0],'parent_event_x':[0,0,0],'actual_ion_x':[0,0,0],'actual_ion_v':[2,1,0],'parent_inertial_v':[0,0,0],'electron_inertial_v':[0,0,0],'photon_inertial_p':[1,0,0],'bath_event_p':[2,-1,0],'velocity_offset_event':[2,0,0],'ledger_tracked_dp':[-6,3,0],'ledger_bath_dp':[6,-3,0]}.items():
        values.update({name+str(i):F(value) for i,value in enumerate(vector)})
    values.update({'actual_parent_compiled_mass_kg':F(2),'actual_ion_compiled_mass_kg':F(1),'actual_parent_charge_c':F(0),'actual_ion_charge_c':F.from_float(1.602176634e-19)})
    values.update({'inertial_to_event_r'+str(i):F(value) for i,value in enumerate([0,-1,0,1,0,0,0,0,1])})
    return profile,mapping,bounds,{k:F(v) for k,v in values.items()}


def write_events(root,values,name='events.hex.csv'):
    path=Path(root)/name
    with path.open('w',newline='') as stream:
        w=csv.writer(stream);w.writerow(['original_parent_id','key','value'])
        for key,value in values.items():
            # All toy values here are dyadic and exactly binary64 representable.
            if F.from_float(float(value))!=value:raise ValueError('fixture value not exact binary64')
            w.writerow(['73',key,float(value).hex()])
    return path

class LedgerTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.profile,self.mapping,self.bounds,self.values=fixture(self.root)
    def score(self,values=None,mapping=None,bounds=None):
        return verify([write_events(self.root,self.values if values is None else values)],self.profile,self.mapping if mapping is None else mapping,self.bounds if bounds is None else bounds)
    def test_independent_toy_frame_signed_bath_all_criteria_close(self):
        result=self.score();self.assertEqual(result['status'],'PASS',result);self.assertGreater(result['counts']['PASS'],60);self.assertEqual(result['events_evaluated'],1)
    def test_each_actual_child_property_has_an_exact_readback_gate(self):
        for key in ['actual_ion_v0','actual_ion_x1','actual_ion_reference_weight','actual_ion_individual_correction','actual_ion_internal_j','actual_ion_compiled_index','actual_ion_catalog_species','actual_ion_child_count','actual_ion_epoch_et_s','actual_ion_interval_elapsed_s']:
            with self.subTest(key=key):
                v=dict(self.values);v[key]+=1;self.assertEqual(self.score(v)['status'],'FAIL')
    def test_native_boolean_cannot_hide_wrong_bath_sign(self):
        v=dict(self.values);v['bath_event_j']=-v['bath_event_j'];v['kinematicConservationVerified']=1
        self.assertEqual(self.score(v)['status'],'FAIL')
    def test_profile_formation_Q_is_recomputed_not_taken_from_native_receipt(self):
        v=dict(self.values);v['source_q_j']+=1;self.assertEqual(self.score(v)['status'],'FAIL')
    def test_missing_actual_initial_readback_is_skipped(self):
        v=dict(self.values);del v['actual_ion_v0'];self.assertEqual(self.score(v)['status'],'SKIPPED')
    def test_evaluated_failure_is_not_hidden_by_a_later_missing_readback(self):
        v=dict(self.values);v['source_q_j']+=1;del v['actual_ion_v0'];r=self.score(v);self.assertEqual(r['status'],'FAIL');self.assertEqual(r['counts']['SKIPPED'],1)
    def test_changed_profile_source_is_infrastructure_error(self):
        (self.root/'analytical-source.txt').write_text('changed');self.assertEqual(self.score()['status'],'ERROR')
    def test_profile_json_identity_cannot_be_replaced(self):
        m=dict(self.mapping);m['profile_sha256']='0'*64;self.assertEqual(self.score(mapping=m)['status'],'ERROR')
    def test_missing_frozen_bound_is_error_not_inferred_tolerance(self):
        b=copy.deepcopy(self.bounds);del b['momentum_kg_m_s'];self.assertEqual(self.score(bounds=b)['status'],'ERROR')
    def test_negative_or_nonfinite_bound_is_error(self):
        for value in ['-1','NaN','Infinity']:
            b=copy.deepcopy(self.bounds);b['mass_kg']['absolute']=value;self.assertEqual(self.score(bounds=b)['status'],'ERROR')
    def test_duplicate_and_nonhex_native_rows_are_errors(self):
        path=write_events(self.root,self.values)
        with path.open('a') as s:s.write('73,parent_number,0x1.8p+1\n')
        self.assertEqual(verify([path],self.profile,self.mapping,self.bounds)['status'],'ERROR')
        path.write_text('original_parent_id,key,value\n73,parent_number,nan\n')
        self.assertEqual(verify([path],self.profile,self.mapping,self.bounds)['status'],'ERROR')
    def test_duplicate_rank_ownership_and_empty_event_identity_are_distinct(self):
        a=write_events(self.root,self.values,'a.csv');b=write_events(self.root,self.values,'b.csv')
        self.assertEqual(verify([a,b],self.profile,self.mapping,self.bounds)['status'],'ERROR')
        a.write_text('original_parent_id,key,value\n');empty=verify([a],self.profile,self.mapping,self.bounds);self.assertEqual(empty['status'],'PASS');self.assertEqual(empty['events_evaluated'],0);self.assertFalse(empty['event_validation_evaluated']);self.assertEqual(empty['counts']['PASS'],1)
        self.assertEqual(verify([],self.profile,self.mapping,self.bounds)['status'],'SKIPPED')
    def test_represented_weight_and_integer_charge_independently_checked(self):
        for key in ['ion_number','tracked_charge_change_e','bath_charge_change_e','ledger_tracked_dp1','ledger_bath_de_j']:
            v=dict(self.values);v[key]+=1;self.assertEqual(self.score(v)['status'],'FAIL')
    def test_frame_reflection_is_refused_by_proper_rotation_gate(self):
        v=dict(self.values);v['inertial_to_event_r8']=-1;self.assertEqual(self.score(v)['status'],'FAIL')
    def test_complete_binary_exponent_span_is_not_narrowed(self):
        self.assertEqual(binary('0x1.000000000000001p-16000'),F((1<<60)+1,1<<(16060)))
    def test_current_real_profile_missing_Na_ion_is_honestly_skipped(self):
        root=Path(__file__).resolve().parents[3]
        p=root/'test_output/srcMoon/m3-physical-completion-20261118T223145Z/state-profiles/attempt3/profile/profile.json'
        if not p.is_file():self.skipTest('retained real source package absent')
        m={'profile_sha256':digest(p),'parent':'23Na:ground','ion':'23Na+:unresolved','electron':'electron:free','event_frame':'LSO','ion_compiled_species':1,'ion_catalog_species':1}
        b=copy.deepcopy(self.bounds);b['event_frame']='LSO'
        result=verify([],p,m,b);self.assertEqual(result['status'],'SKIPPED');self.assertIn('unavailable',result['reason'])
    def test_actual_compiled_mass_and_SI_charge_are_not_taken_from_plan(self):
        for key in ['actual_parent_compiled_mass_kg','actual_ion_compiled_mass_kg','actual_parent_charge_c','actual_ion_charge_c']:
            v=dict(self.values);v[key]=v[key]*2 if key=='actual_ion_charge_c' else v[key]+1;self.assertEqual(self.score(v)['status'],'FAIL')
    def test_real_frozen_D04_endpoint_is_independently_selected_without_state_activation(self):
        base=Path('/data/vtenishe/moon_validation_data/photochemistry/processed/phidrates')
        if not base.is_dir():self.skipTest('canonical D04 package unavailable')
        hashes={'photo_rates_1au.csv':'3aa54915a08258371dbe975949e6afb2ab74618164fadca8bceebe876c926847','thresholds.csv':'a3335b3b37c3d275a4f72e92ae8e5b1f6ecfcbc4e2d84735cb20694dc93b48f1'}
        before={str(base/name):(digest(base/name),(base/name).stat().st_mtime_ns) for name in hashes}
        for activity,rate,mean in [('0','0.00000726','0.972'),('1','0.00000791','2.98')]:
            source={'reaction_id':'PHIDRATES:Na:00:na_e','solar_activity':activity,'rates':{'path':str(base/'photo_rates_1au.csv'),'sha256':hashes['photo_rates_1au.csv']},'thresholds':{'path':str(base/'thresholds.csv'),'sha256':hashes['thresholds.csv']}}
            actual=load_photon_source({'photon_source':source},base,{}, {'elementary charge':(F('1.602176634e-19'),'C')},'source-bound-experimental')
            self.assertEqual(actual['rate'],F(rate));self.assertEqual(actual['mean'],F(mean)*F('1.602176634e-19'));self.assertEqual(actual['threshold'],F('5.13907097104346')*F('1.602176634e-19'))
            source['rates']['sha256']='0'*64
            with self.assertRaises(ValueError):load_photon_source({'photon_source':source},base,{}, {'elementary charge':(F('1.602176634e-19'),'C')},'source-bound-experimental')
        self.assertEqual(before,{str(base/name):(digest(base/name),(base/name).stat().st_mtime_ns) for name in hashes})
    def test_cli_malformed_metadata_is_error_and_existing_artifacts_are_immutable(self):
        events=write_events(self.root,self.values);mapping=self.root/'mapping.json';bounds=self.root/'bounds.json';mapping.write_text('{broken');bounds.write_text(json.dumps(self.bounds))
        output=self.root/'cli-error';command=[sys.executable,'-B',str(Path(__file__).with_name('verify_ledger.py')),'--events',str(events),'--profile',str(self.profile),'--mapping',str(mapping),'--bounds',str(bounds),'--output-dir',str(output)]
        child=subprocess.run(command,capture_output=True,text=True,timeout=10);self.assertEqual(child.returncode,2);self.assertEqual(json.loads((output/'result.json').read_text())['status'],'ERROR')
        before={str(p):(digest(p),p.stat().st_mtime_ns) for p in output.rglob('*') if p.is_file()}
        child=subprocess.run(command,capture_output=True,text=True,timeout=10);self.assertEqual(child.returncode,2)
        self.assertEqual(before,{str(p):(digest(p),p.stat().st_mtime_ns) for p in output.rglob('*') if p.is_file()})
    def test_result_counts_exactly_match_all_retained_rows(self):
        v=dict(self.values);v['source_q_j']+=1;del v['actual_ion_v0'];r=self.score(v)
        self.assertEqual(r['counts'],{s:sum(row['status']==s for row in r['checks']) for s in r['counts']})
if __name__=='__main__':unittest.main()
