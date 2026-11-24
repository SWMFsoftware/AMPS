#!/usr/bin/env python3
"""Analytical and administrative controls; NEVER native/physical cohort evidence."""
from __future__ import annotations
import argparse,copy,csv,io,json,math,struct,subprocess,sys,tempfile,time,traceback,unittest
from decimal import Decimal,localcontext
from fractions import Fraction
from pathlib import Path
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).with_name('sodium_cohort')
sys.path.insert(0,str(HERE))
import verify_numerical as numerical
import audit
import m3_qualification as q
import m3_qualification_executor as executor
from test_m3_sodium_ledger import Recorded

COHORT='ANALYTICAL-FIRST-HIT-NOT-OBSERVATIONS'
CHANNEL='PHIDRATES:Na:00:na_e'


def reference(visible=1):
    return {'solar_endpoint':'quiet','total_acceleration_macro':'TotalParticleAcceleration_default(a,s,p,x,v,n);',
        'point_sun_visible':visible,'sun_position_lso_m':['1','0','0'],
        'source_position_lso_m':['0','0','0'],'source_reference_distance_m':'1',
        'source_reference_distance_identity':'ANALYTICAL-UNIT-DISTANCE-NOT-SPICE',
        'geometry_basis_epoch_et_s':'0'}


def plan(visible=1):
    rate=float('7.26e-6') if visible else 0.
    return {'schema':'srcMoon-M3-cohort-plan-v1','approved_addendum_sha256':q.APPROVED_SHA256,
        'test_id':'I21','cohort_id':COHORT,'original_parents':16,'initial_represented_number':.125,
        'horizon_s':262144.,'reference_total_rate_s':rate,
        'source_channel_rates_s':[{'id':CHANNEL,'rate_s':rate}],'branch_claims':1,'seed':731993,
        'primary_estimator':'lifetime' if visible else 'zero-hazard',
        'test_kind':'homogeneous-fixed-position' if visible else 'zero-hazard-invariant',
        'source_profile':{'path':'MISSING_REAL_PHYSICAL_PROFILE.json','sha256':'0'*64}}


def first_hit(parent,configuration):
    # Independent Decimal reference, retaining the actual source scalar's
    # binary64 representation. Exact HMAC uniform is checked against literals
    # separately; this fixture contains no observed or model-generated truth.
    if configuration['reference_total_rate_s']==0:return None
    with localcontext() as c:
        c.prec=150
        u=Decimal.from_float(numerical.source_uniform(configuration['seed'],configuration['cohort_id'],parent))
        return -u.ln()/Decimal.from_float(configuration['reference_total_rate_s'])


def write_level(root,configuration,steps,ranks,error,visible=1):
    folder=Path(root)/(str(ranks)+'mpi-'+str(steps));folder.mkdir()
    ref=reference(visible)
    receipt={'schema':'srcMoon-native-sodium-input-receipt-v1','fixture_scope':'ANALYTICAL-NOT-NATIVE-NOT-OBSERVATIONS',
        'homogeneous_reference':ref,'geometry_basis_et_s':'0','case_kind':'homogeneous-fixed-position',
        'original_parents':configuration['original_parents'],'seed':configuration['seed'],
        'cohort_id':configuration['cohort_id'],'physical_horizon_s':str(configuration['horizon_s']),
        'physical_step_s':str(configuration['horizon_s']/steps),'physical_steps':steps,
        'initial_represented_number':str(configuration['initial_represented_number'])}
    (folder/'input_receipt.json').write_text(json.dumps(receipt))
    events={}
    with (folder/'parent_records.csv').open('w',newline='') as stream:
        w=csv.writer(stream);w.writerow(q.PARENT_COLUMNS)
        for parent in range(configuration['original_parents']):
            event=first_hit(parent,configuration);reacted=event is not None and event<=Decimal.from_float(configuration['horizon_s'])
            actual=event+Decimal(str(error))*Decimal.from_float(configuration['horizon_s']) if reacted else None
            if reacted:
                # Project the analytic fixture onto a binary64 clock and
                # reproduce the native parent CSV output format. Tests below
                # also supply a finer native bit that the parent CSV loses.
                events[str(parent)]={'relative_event_time_s':Fraction.from_float(float(actual)),
                    'accepted_rate_s':Fraction.from_float(configuration['reference_total_rate_s'])}
            w.writerow([parent,'REACTED' if reacted else 'SURVIVED',
                format(float(actual),'.17g') if reacted else '',CHANNEL if reacted else '',configuration['initial_represented_number']])
    return folder,ranks,events


def request_fixture(root):
    root=Path(root);producer=root/'ADMINISTRATIVE-producer.json';producer.write_text('{"scope":"ADMINISTRATIVE ONLY; NOT A NATIVE BUILD"}')
    included=root/'included.in';included.write_text('! administrative include\n')
    main=root/'request.in';main.write_text('#include "included.in"\n')
    p=plan();p['native_adapter']={'schema':'srcMoon-M3-normal-sodium-execution-v1','producer_id':'normal-source-bound-sodium-cohort-v1',
        'producer_manifest':{'path':str(producer),'sha256':q.sha(producer)},'state_mapping':{'parent':'not-a-real-parent','ion':'not-a-real-ion','electron':'not-a-real-electron'},
        'runs':[{'mpi_ranks':ranks,'physical_steps':steps,'input':{'path':str(main),'sha256':q.sha(main)},
                 'input_documents':[{'path':str(x),'sha256':q.sha(x)} for x in [main,included]]} for ranks in (1,8) for steps in (8,16,32)]}
    return p,main,included


def owner_identity_fixture(root,receipt,documents,steps=8,ranks=1):
    """Administrative source-document identity ONLY; never physics evidence."""
    root=Path(root);native=root/'native';native.mkdir(exist_ok=True)
    (native/'input_receipt.json').write_text(json.dumps(receipt))
    (native/'parent_records.csv').write_text('ADMINISTRATIVE BYTES ONLY; NOT PARENT SCIENCE\n')
    for rank in range(ranks):(native/('event_metrics.rank%d.hex.csv'%rank)).write_text('original_parent_id,key,value\n')
    input_descriptor=documents[0]
    manifest_path=root/'ADMINISTRATIVE-owner-manifest.json';result_path=root/'ADMINISTRATIVE-owner-result.json'
    manifest={'scope':'ADMINISTRATIVE IDENTITY CONTROL; NOT A NATIVE LAUNCH','execution':'native-run','mode':'run','mpi_ranks':ranks,'openmp_threads':1,'input':input_descriptor,
        'artifacts':{str(p.relative_to(root)):q.sha(p) for p in native.iterdir() if p.is_file()}}
    manifest_path.write_text(json.dumps(manifest));result_path.write_text(json.dumps({'status':'PASS','scope':'ADMINISTRATIVE IDENTITY ONLY; NO PHYSICAL RUN'}))
    return {'output':str(native),'mpi_ranks':ranks,'physical_steps':steps,'input_documents':documents,'input':input_descriptor,
        'owner_manifest':{'path':str(manifest_path),'sha256':q.sha(manifest_path)},'owner_result':{'path':str(result_path),'sha256':q.sha(result_path)}}


def administrative_launch_matrix(root,configuration,first):
    """Keep the frozen six launch identities intact in corrupted-graph tests.

    The owning audit checks the complete launch matrix before reading a
    single input graph. Supply the other five administrative identities too,
    so an expected graph error exercises the intended check.
    """
    launches=[first]
    for index,case in enumerate(configuration['native_adapter']['runs'][1:],start=1):
        folder=root/('ADMIN-launch-%d'%index)
        folder.mkdir()
        docs=case['input_documents']
        receipt={'input_documents':docs,'physical_steps':case['physical_steps']}
        launches.append(owner_identity_fixture(folder,receipt,docs,case['physical_steps'],case['mpi_ranks']))
    return launches

class NumericalControls(unittest.TestCase):
    def test_finite_lunar_shadow_geometry(self):
        # Analytic diametral crossing, off-axis grazing, and a line which
        # intersects behind either endpoint distinguish a finite segment from
        # an infinite line or an unbounded ray in both directions.
        segment=numerical.lunar_segment_geometry
        self.assertEqual(segment([-2,0,0],[6,0,0],1)['closest_segment_fraction'],'1/4')
        self.assertEqual(segment([-2,0,0],[6,0,0],1)['intersection'],'interior')
        self.assertEqual(segment([-2,1,0],[6,1,0],1)['intersection'],'tangent')
        for start,end,closest in (([2,0,0],[4,0,0],'0'),([-4,0,0],[-2,0,0],'1')):
            result=segment(start,end,1)
            self.assertEqual(result['intersection'],'exterior')
            self.assertEqual(result['closest_segment_fraction'],closest)
            self.assertEqual(result['minimum_distance_m'],2)
        with self.assertRaises(ValueError):segment([2,0,0],[2,0,0],1)

    def test_fixed_shadow_refuses_mutated_radius_and_detects_wrong_visibility(self):
        # Administrative controls bind the analytical unit sphere to a
        # synthetic producer manifest. These files never enter a native plan.
        header_key='src/general/constants.PlanetaryData.h'
        header=self.root/header_key
        header.parent.mkdir(parents=True)
        header.write_text('#define _MOON__RADIUS_ 1.0E0\n')
        producer=self.root/'producer.json'
        producer.write_text(json.dumps({'source_root':str(self.root),
            'production_source_hashes':{header_key:q.sha(header)}}))
        configuration={'native_adapter':{'producer_manifest':{'path':str(producer),'sha256':q.sha(producer)}}}
        app=self.root/'sphere.in'
        app.write_text('#begin application: moon\nsurface_geometry = sphere\nsimulation_environment_mode = frozen\n#end application\n')
        native=self.root/'native.in'
        native.write_text('earth_occultation = off\n')
        launch={'input':{'path':str(native)},'input_documents':[{'path':str(x)} for x in (native,app)]}
        ref=reference(0)
        ref.update(source_position_lso_m=['-2','0','0'],sun_position_lso_m=['2','0','0'])
        receipt={'case_kind':'homogeneous-fixed-position','schema':'srcMoon-native-sodium-input-receipt-v1',
                 'homogeneous_reference':ref}
        self.assertEqual(numerical.fixed_moon_shadow(configuration,launch,receipt)['status'],'PASS')
        ref['point_sun_visible']=1
        self.assertEqual(numerical.fixed_moon_shadow(configuration,launch,receipt)['status'],'FAIL')
        header.write_text('#define _MOON__RADIUS_ 1.5\n')
        with self.assertRaisesRegex(ValueError,'lunar-radius source differs'):
            numerical.fixed_moon_shadow(configuration,launch,receipt)

    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
    def levels(self,errors=(.0001,.00005,.00001),visible=1):
        p=plan(visible)
        return p,[write_level(self.root,p,steps,ranks,error,visible) for ranks in (1,8) for steps,error in zip((8,16,32),errors)]
    def test_source_hmac_golden_open_uniform_vectors(self):
        expected={0:'0x1.808c32cde1ed6p-2',1:'0x1.3d010e1e2ae46p-2',7:'0x1.3cd835abd3ecap-2'}
        for parent,value in expected.items():self.assertEqual(numerical.source_uniform(731993,COHORT,parent).hex(),value)
        self.assertNotEqual(numerical.source_uniform(731993,COHORT,0),numerical.source_uniform(731993,COHORT,1))
    def test_native_nested_geometry_schema_is_used_without_inferred_visibility(self):
        r={'schema':'srcMoon-native-sodium-input-receipt-v1','homogeneous_reference':reference(),'geometry_basis_et_s':'0'}
        self.assertEqual(numerical.reference_receipt(r),reference())
        r['point_sun_visible']=0
        with self.assertRaises(ValueError):numerical.reference_receipt(r)
        del r['point_sun_visible'];del r['homogeneous_reference']['point_sun_visible']
        with self.assertRaises((ValueError,KeyError)):numerical.reference_receipt(r)
    def test_real_D04_binary_rate_and_inverse_square_are_independent(self):
        r=reference();ids={};actual,channel=numerical.source_rate(r,ids)
        self.assertEqual(actual,Decimal.from_float(7.26e-6));self.assertEqual(channel,CHANNEL);self.assertGreaterEqual(len(ids),14)
        r['sun_position_lso_m']=['2','0','0'];scaled,_=numerical.source_rate(r,{})
        with localcontext() as context:
            context.prec=80;expected=Decimal.from_float(7.26e-6)/4
        self.assertEqual(scaled,expected)
    def test_declared_zero_visibility_is_exact_zero_hazard(self):
        self.assertEqual(numerical.source_rate(reference(0),{})[0],0)
        p,runs=self.levels((0,0,0),0);actual=numerical.verify(runs,p)
        self.assertEqual(actual['status'],'PASS');self.assertTrue(all(row['reference_event_s'] is None for row in actual['comparisons']))
    def test_Decimal_first_hits_all_original_parents_and_layouts(self):
        p,runs=self.levels();actual=numerical.verify(runs,p)
        self.assertEqual(actual['status'],'PASS');self.assertEqual(len(actual['comparisons']),p['original_parents']*6)
        self.assertTrue(any(r['expected_status']=='SURVIVED' for r in actual['comparisons']))
        self.assertTrue(any(r['expected_status']=='REACTED' for r in actual['comparisons']))
        self.assertEqual(actual['scores']['1']['approved_relative_bias_limit'],.002)
    def test_strict_decrease_and_original_numeric_allocations_are_preserved(self):
        self.assertEqual(q.NUMERICAL_LIMIT,{'I21':.002,'I22':.002,'I23':.003,'I24':.003})
        p,runs=self.levels((.0001,.0001,.0001));self.assertEqual(numerical.verify(runs,p)['status'],'FAIL')
    def test_source_clock_precision_binding_and_partial_diagnostics(self):
        p,runs=self.levels()
        p['optical_clock_precision']='compensated-extended-log'
        with self.assertRaisesRegex(ValueError,'source-clock precision'):numerical.verify(runs,p)
        for path,_,_ in runs:
            key=path/'input_receipt.json';receipt=json.loads(key.read_text())
            receipt['optical_clock_precision']='compensated-extended-log';key.write_text(json.dumps(receipt))
        scored=numerical.verify(runs,p)
        self.assertEqual(scored['status'],'PASS')
        for row in scored['homogeneous_timing_breakdown']:
            self.assertEqual(row['optical_clock_precision'],p['optical_clock_precision'])
            self.assertNotIn('source_clock_rounding_s',row['component_extrema'])
            self.assertIn('native_clock_and_transport_remainder_s',row['component_extrema'])

    def test_score_uses_native_bits_despite_identical_binary64_parent_projection(self):
        p,runs=self.levels()
        native=runs[0][2]
        key=next(iter(native))
        base=native[key]['relative_event_time_s']
        # Successor of this normal binary64 value on a 64-significand-bit
        # grid, as on the x87 native host. Its projection to binary64 is unchanged.
        shifted=base+Fraction(2)**(math.frexp(float(base))[1]-64)
        self.assertEqual(float(base),float(shifted))
        before=numerical.verify(runs,p)['comparisons']
        native[key]['relative_event_time_s']=shifted
        after=numerical.verify(runs,p)['comparisons']
        parent=int(key)
        self.assertEqual(before[parent]['original_parent_id'],parent)
        self.assertNotEqual(before[parent]['normalized_error'],after[parent]['normalized_error'])
        # Conversely, printing the exact decimal expansion of the double is
        # the very same CSV projection and must leave physical timing intact.
        path=runs[0][0]/'parent_records.csv'
        with path.open(newline='') as stream:table=list(csv.reader(stream))
        table[parent+1][2]=str(Decimal.from_float(float(base)))
        with path.open('w',newline='') as stream:csv.writer(stream).writerows(table)
        self.assertEqual(after,numerical.verify(runs,p)['comparisons'])
        del native[key]['relative_event_time_s']
        with self.assertRaises((KeyError,ValueError)):numerical.verify(runs,p)
    def test_clock_rate_transport_and_reporting_components_have_independent_signs(self):
        # All reference endpoints below are terminating binary fractions.
        # At unit optical depth the +2^-54 increment falls inside a binary64
        # half-ULP. A 2/s accepted rate halves the rounded reference time.
        with localcontext() as c:
            c.prec=80
            depth=Decimal(1)+Decimal(2)**-54
            rate=Fraction(2)
            seconds=Fraction(3,4)+Fraction(1,2**54)
            displayed='0.75000000000000001'
            values=numerical.homogeneous_timing_breakdown(depth,Decimal(1),
                {'accepted_rate_s':rate,'relative_event_time_s':seconds},displayed)
            self.assertEqual(values['source_clock_rounding_s'],-Decimal(2)**-54)
            self.assertEqual(values['accepted_rate_effect_s'],Decimal('-.5'))
            self.assertEqual(values['native_clock_and_transport_remainder_s'],Decimal('.25')+Decimal(2)**-54)
            self.assertEqual(values['parent_binary64_projection_s'],-Decimal(2)**-54)
            self.assertEqual(values['parent_csv_text_rendering_s'],Decimal('1e-17'))
            self.assertEqual(values['native_minus_reference_s'],Decimal('-.25'))
            self.assertEqual(sum(list(values.values())[:5]),values['csv_text_minus_reference_s'])
    def test_decreasing_but_over_budget_is_still_failed(self):
        p,runs=self.levels((.02,.01,.0021));self.assertEqual(numerical.verify(runs,p)['status'],'FAIL')
    def test_wrong_first_hit_membership_is_evaluated_failure(self):
        p,runs=self.levels();path=runs[0][0]/'parent_records.csv'
        with path.open(newline='') as s:rows=list(csv.reader(s))
        row=next(r for r in rows[1:] if r[1]=='REACTED');row[1:4]=['SURVIVED','','']
        with path.open('w',newline='') as s:csv.writer(s).writerows(rows)
        with self.assertRaises(ArithmeticError):numerical.verify(runs,p)
    def test_corrupt_parent_CSV_is_infrastructure_not_membership_failure(self):
        p,runs=self.levels();(runs[0][0]/'parent_records.csv').write_text('wrong,header\n')
        with self.assertRaises(ValueError):numerical.verify(runs,p)
    def test_variable_hazard_cannot_enter_homogeneous_scorer(self):
        p,runs=self.levels();path=runs[0][0]/'input_receipt.json';v=json.loads(path.read_text());v['case_kind']='trajectory-integrated-hazard';path.write_text(json.dumps(v))
        with self.assertRaises(ValueError):numerical.verify(runs,p)
    def test_changed_seed_geometry_or_original_parent_count_is_error(self):
        p,runs=self.levels();path=runs[0][0]/'input_receipt.json';original=json.loads(path.read_text())
        for key,value in [('seed',104729),('original_parents',15),('physical_horizon_s','131072')]:
            changed=copy.deepcopy(original);changed[key]=value;path.write_text(json.dumps(changed))
            with self.assertRaises(ValueError):numerical.verify(runs,p)
        path.write_text(json.dumps(original));changed=copy.deepcopy(original);changed['homogeneous_reference']['geometry_basis_epoch_et_s']='1';changed['geometry_basis_et_s']='1';path.write_text(json.dumps(changed))
        with self.assertRaises(ValueError):numerical.verify(runs,p)
    def test_nan_boolean_and_bad_decimal_geometry_are_errors_not_scientific_failures(self):
        for value in ['NaN','Infinity','not-a-number',True]:
            r=reference();r['sun_position_lso_m'][0]=value
            with self.assertRaises(ValueError):numerical.source_rate(r,{})
    def test_incomplete_or_mismatched_refinement_matrix_is_error(self):
        p,runs=self.levels()
        with self.assertRaises(ValueError):numerical.verify(runs[:-1],p)
        with self.assertRaises(ValueError):numerical.verify(runs+[runs[0]],p)
    def test_nonzero_force_macro_refuses_constant_position_oracle(self):
        for macro in ('Moon::TotalParticleAcceleration','PIC::Mover::TotalParticleAcceleration_default','TotalParticleAcceleration_default','TotalParticleAcceleration_default(a,s,p,x,v,n);otherForce(a,s,p,x,v,n);','TotalParticleAcceleration_default(a,s,p,x,v,n);suffix'):
            r=reference();r['total_acceleration_macro']=macro
            with self.assertRaises(ValueError):numerical.source_rate(r,{})
        for macro in ('TotalParticleAcceleration_default(a,s,p,x,v,n);','PIC::Mover::TotalParticleAcceleration_default(a,s,p,x,v,n)'):
            r=reference();r['total_acceleration_macro']=macro;self.assertGreater(numerical.source_rate(r,{})[0],0)
    def test_existing_inverse_square_gate_checks_both_actual_rate_fields(self):
        correct={'73':{'source_rate_event_s':Fraction(1000),'accepted_rate_s':Fraction(1000)}}
        result=numerical.inverse_square_checks(correct,Decimal(1000))
        self.assertEqual(result['status'],'PASS');self.assertEqual(len(result['checks']),2)
        self.assertEqual(numerical.I21_INVERSE_SQUARE_RELATIVE_LIMIT,Decimal('.001'))
        for key in ('source_rate_event_s','accepted_rate_s'):
            changed=copy.deepcopy(correct);changed['73'][key]=Fraction(1002);changed['73']['loose_caller_bound']=Fraction(1000000)
            self.assertEqual(numerical.inverse_square_checks(changed,Decimal(1000))['status'],'FAIL')
    def test_inverse_square_exact_boundary_and_next_binary_value(self):
        for value in (Fraction(1001),Fraction(999)):
            events={'73':{'source_rate_event_s':value,'accepted_rate_s':value}}
            self.assertEqual(numerical.inverse_square_checks(events,Decimal(1000))['status'],'PASS')
        value=Fraction(1001)+Fraction(1,2**40)
        events={'73':{'source_rate_event_s':value,'accepted_rate_s':value}}
        self.assertEqual(numerical.inverse_square_checks(events,Decimal(1000))['status'],'FAIL')
    def test_inverse_square_zero_and_empty_are_explicit_identities(self):
        result=numerical.inverse_square_checks({},Decimal(0));self.assertEqual(result['status'],'PASS');self.assertEqual(result['events_evaluated'],0);self.assertFalse(result['checks']);self.assertIn('empty',result['scope'])
        zeros={'73':{'source_rate_event_s':Fraction(0),'accepted_rate_s':Fraction(0)}}
        self.assertEqual(numerical.inverse_square_checks(zeros,Decimal(0))['status'],'PASS')
        zeros['73']['accepted_rate_s']=Fraction(1,2**1074)
        self.assertEqual(numerical.inverse_square_checks(zeros,Decimal(0))['status'],'FAIL')
    def test_inverse_square_missing_actual_metric_is_infrastructure_error(self):
        with self.assertRaises(ValueError):numerical.inverse_square_checks({'73':{'accepted_rate_s':Fraction(1000)}},Decimal(1000))
    def test_pinned_request_checks_real_bytes_and_root_graph_membership(self):
        p,main,included=request_fixture(self.root);request=executor.native_request(p,self.root)
        self.assertEqual(len(request['runs']),6);self.assertTrue(all(included in r['document_paths'] for r in request['runs']))
        included.write_text('changed')
        with self.assertRaises(ValueError):executor.native_request(p,self.root)
    def test_missing_duplicate_or_rootless_input_graph_is_rejected(self):
        p,_,_=request_fixture(self.root)
        for edit in ('missing','duplicate','rootless'):
            v=copy.deepcopy(p)
            if edit=='missing':del v['native_adapter']['runs'][0]['input_documents']
            elif edit=='duplicate':v['native_adapter']['runs'][0]['input_documents'].append(v['native_adapter']['runs'][0]['input_documents'][0])
            else:v['native_adapter']['runs'][0]['input_documents']=v['native_adapter']['runs'][0]['input_documents'][1:]
            with self.assertRaises(ValueError):executor.native_request(v,self.root)
    def test_adapter_request_cannot_choose_arbitrary_command_or_mpi_layout(self):
        p,_,_=request_fixture(self.root)
        for key,value in [('mpi_ranks',2),('physical_steps',True),('shell_command','echo fake')]:
            v=copy.deepcopy(p);v['native_adapter']['runs'][0][key]=value
            with self.assertRaises(ValueError):executor.native_request(v,self.root)
    def test_exact_native_input_graph_mismatch_is_retained_as_ERROR(self):
        p,main,included=request_fixture(self.root);profile=self.root/'ADMIN-profile.json';profile.write_text('{"scope":"ADMINISTRATIVE ONLY"}');p['source_profile']={'path':str(profile),'sha256':q.sha(profile)}
        planpath=self.root/'plan.json';planpath.write_text(json.dumps(p));native=self.root/'native';native.mkdir()
        docs=p['native_adapter']['runs'][0]['input_documents'];launch=owner_identity_fixture(self.root,{'input_documents':docs[:1],'physical_steps':8},docs)
        registry=self.root/'runs.json';registry.write_text(json.dumps(administrative_launch_matrix(self.root,p,launch)))
        result=audit.run(planpath,registry,self.root/'audit');self.assertEqual(result['status'],'ERROR');self.assertIn('include graph',result['checks'][-1]['reason'])
    def test_duplicate_actual_input_graph_is_not_silently_collapsed(self):
        p,_,_=request_fixture(self.root);profile=self.root/'ADMIN-profile.json';profile.write_text('{}');p['source_profile']={'path':str(profile),'sha256':q.sha(profile)}
        planpath=self.root/'plan.json';planpath.write_text(json.dumps(p));native=self.root/'native';native.mkdir();docs=p['native_adapter']['runs'][0]['input_documents']
        launch=owner_identity_fixture(self.root,{'input_documents':docs+[docs[0]],'physical_steps':8},docs)
        registry=self.root/'runs.json';registry.write_text(json.dumps(administrative_launch_matrix(self.root,p,launch)))
        result=audit.run(planpath,registry,self.root/'audit');self.assertEqual(result['status'],'ERROR');self.assertIn('duplicate',result['checks'][-1]['reason'].lower())
    def test_portable_adjacent_binary64_matches_independent_literal_boundaries(self):
        # These literal IEEE encodings are independent of the implementation's
        # struct stepping, and include binade asymmetry and signed zero.
        cases=[('0x0.0p+0',True,'0x0.0000000000001p-1022'),('-0x0.0p+0',True,'0x0.0000000000001p-1022'),
            ('0x0.0p+0',False,'-0x0.0000000000001p-1022'),('-0x0.0p+0',False,'-0x0.0000000000001p-1022'),
            ('0x0.0000000000001p-1022',False,'0x0.0p+0'),('-0x0.0000000000001p-1022',True,'-0x0.0p+0'),
            ('0x1.0p+0',True,'0x1.0000000000001p+0'),('0x1.0p+0',False,'0x1.fffffffffffffp-1'),
            ('-0x1.0p+0',True,'-0x1.fffffffffffffp-1'),('-0x1.0p+0',False,'-0x1.0000000000001p+0'),
            ('0x1.0p-1',False,'0x1.fffffffffffffp-2'),('0x1.fffffffffffffp+1023',True,'inf'),
            ('-0x1.fffffffffffffp+1023',False,'-inf'),('0x1.fffffffffffffp+1023',False,'0x1.ffffffffffffep+1023')]
        for value,upward,expected in cases:
            self.assertEqual(numerical.adjacent_binary64(float.fromhex(value),upward).hex(),float.fromhex(expected).hex())
    def test_portable_adjacent_invalid_or_nonfinite_epoch_envelope_refuses(self):
        for value,direction in [(math.nan,True),(math.inf,False),(1,True),(1.,1)]:
            with self.assertRaises(ValueError):numerical.adjacent_binary64(value,direction)
        maximum=float.fromhex('0x1.fffffffffffffp+1023')
        with self.assertRaises(ValueError):numerical.event_time_links({'73':{'relative_event_time_s':Fraction(0),'event_epoch_et_s':Fraction.from_float(maximum)}},[{'original_parent_id':'73','status':'REACTED','relative_event_time_s':'0'}],str(maximum))
    def test_time_link_preserves_long_relative_and_binary64_CSV_projection(self):
        relative=Fraction(1)+Fraction(1,2**70);origin=1000000000.
        events={'73':{'relative_event_time_s':relative,'event_epoch_et_s':Fraction.from_float(float(Fraction.from_float(origin)+relative))}}
        parents=[{'original_parent_id':'73','status':'REACTED','relative_event_time_s':format(float(relative),'.17g')}]
        result=numerical.event_time_links(events,parents,str(origin));self.assertEqual(result['status'],'PASS');self.assertEqual(len(result['checks']),2)
    def test_time_link_wrong_CSV_relative_and_actual_epoch_fail(self):
        events={'73':{'relative_event_time_s':Fraction(1,4),'event_epoch_et_s':Fraction(1000000000)+Fraction(1,4)}}
        parents=[{'original_parent_id':'73','status':'REACTED','relative_event_time_s':'.5'}]
        self.assertEqual(numerical.event_time_links(events,parents,'1000000000')['status'],'FAIL')
        parents[0]['relative_event_time_s']='.25';events['73']['event_epoch_et_s']+=1
        self.assertEqual(numerical.event_time_links(events,parents,'1000000000')['status'],'FAIL')
    def test_time_link_existing_two_adjacent_ULPs_allowed_three_fail(self):
        origin=1000000000.;relative=Fraction(1,4);projected=origin+.25
        parents=[{'original_parent_id':'73','status':'REACTED','relative_event_time_s':'.25'}]
        for direction in (-math.inf,math.inf):
            epoch=projected
            for distance in range(1,4):
                # Independent integer representation of the known positive
                # binary64 origin fixture; no dependency on math.nextafter
                # (unavailable in supported system Python3.8).
                bits=struct.unpack('>Q',struct.pack('>d',projected))[0]
                epoch=struct.unpack('>d',struct.pack('>Q',bits+(distance if direction>0 else -distance)))[0]
                events={'73':{'relative_event_time_s':relative,'event_epoch_et_s':Fraction.from_float(epoch)}}
                self.assertEqual(numerical.event_time_links(events,parents,str(origin))['status'],'PASS' if distance<=2 else 'FAIL')
    def test_time_link_geometry_basis_is_not_an_assumed_transport_origin(self):
        events={'73':{'relative_event_time_s':Fraction(1,4),'event_epoch_et_s':Fraction(1000000000)+Fraction(1,4),'geometry_basis_et_s':Fraction(-4000000000)}}
        parents=[{'original_parent_id':'73','status':'REACTED','relative_event_time_s':'.25'}]
        self.assertEqual(numerical.event_time_links(events,parents,'1000000000')['status'],'PASS')
    def test_time_link_missing_origin_metrics_or_membership_is_error(self):
        full={'73':{'relative_event_time_s':Fraction(1,4),'event_epoch_et_s':Fraction(1000000000)+Fraction(1,4)}}
        parents=[{'original_parent_id':'73','status':'REACTED','relative_event_time_s':'.25'}]
        for origin in (None,'NaN',True):
            with self.assertRaises(ValueError):numerical.event_time_links(full,parents,origin)
        for key in ('relative_event_time_s','event_epoch_et_s'):
            events=copy.deepcopy(full);del events['73'][key]
            with self.assertRaises(ValueError):numerical.event_time_links(events,parents,'1000000000')
        with self.assertRaises(ValueError):numerical.event_time_links(full,[],'1000000000')
    def test_original_owner_artifact_binding_accepts_identity_only(self):
        p,_,_=request_fixture(self.root);docs=p['native_adapter']['runs'][0]['input_documents']
        launch=owner_identity_fixture(self.root,{'input_documents':docs,'physical_steps':8},docs)
        self.assertEqual(audit.authenticate_launch(launch),self.root/'native')
    def test_changed_consumed_CSV_cannot_be_rehashed_into_an_original_launch(self):
        p,_,_=request_fixture(self.root);docs=p['native_adapter']['runs'][0]['input_documents']
        launch=owner_identity_fixture(self.root,{'input_documents':docs,'physical_steps':8},docs)
        (self.root/'native/parent_records.csv').write_text('ALTERED VALID-LOOKING MODEL OUTPUT IS NOT ORIGINAL\n')
        with self.assertRaisesRegex(ValueError,'artifact changed'):audit.authenticate_launch(launch)
    def test_replaced_owner_manifest_or_result_is_error(self):
        p,_,_=request_fixture(self.root);docs=p['native_adapter']['runs'][0]['input_documents']
        launch=owner_identity_fixture(self.root,{'input_documents':docs,'physical_steps':8},docs)
        for key in ('owner_manifest','owner_result'):
            path=Path(launch[key]['path']);before=path.read_bytes();path.write_text('{}')
            with self.assertRaisesRegex(ValueError,'owner identity changed'):audit.authenticate_launch(launch)
            path.write_bytes(before)
    def test_wrong_owner_layout_or_frozen_root_input_is_error(self):
        p,_,_=request_fixture(self.root);docs=p['native_adapter']['runs'][0]['input_documents']
        launch=owner_identity_fixture(self.root,{'input_documents':docs,'physical_steps':8},docs)
        wrong=copy.deepcopy(launch);wrong['mpi_ranks']=8
        with self.assertRaisesRegex(ValueError,'this actual native'):audit.authenticate_launch(wrong)
        wrong=copy.deepcopy(launch);wrong['input']['sha256']='0'*64
        with self.assertRaisesRegex(ValueError,'frozen native root input'):audit.authenticate_launch(wrong)
    def test_consumed_event_file_must_be_in_original_owner_map(self):
        p,_,_=request_fixture(self.root);docs=p['native_adapter']['runs'][0]['input_documents']
        launch=owner_identity_fixture(self.root,{'input_documents':docs,'physical_steps':8},docs)
        path=Path(launch['owner_manifest']['path']);data=json.loads(path.read_text());del data['artifacts']['native/event_metrics.rank0.hex.csv'];path.write_text(json.dumps(data));launch['owner_manifest']['sha256']=q.sha(path)
        with self.assertRaisesRegex(ValueError,'not bound'):audit.authenticate_launch(launch)
    def test_dispatch_failure_provenance_never_becomes_not_run_or_physicsPASS(self):
        p=plan();planpath=self.root/'plan.json';planpath.write_text(json.dumps(p))
        ready={'status':'PASS','identities':[],'blocking_reasons':[],'scope':'ADMINISTRATIVE-DISPATCH-FIXTURE-NOT-SCIENCE'}
        attempted={'status':'ERROR','reason':'CONTROLLED DISPATCH FAILURE; NO PHYSICS RUN','execution':'native-attempted','criteria_evaluated':False,'commands':[{'command':['ADMINISTRATIVE-NOT-EXECUTED'],'returncode':2}]}
        with patch.object(executor,'readiness',return_value=ready),patch.object(executor,'execute_cohort',return_value=attempted):
            result=executor.run(planpath,self.root,self.root/'attempt')
        self.assertEqual(result['status'],'ERROR');manifest=q.strict_object(self.root/'attempt/run_manifest.json')
        self.assertEqual(manifest['execution'],'native-attempted');self.assertFalse(manifest['criteria_evaluated'])


def run(output):
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    manifest={'scope':'analytical timing and administrative protocol software ONLY; no actual native/physical cohort','command':[sys.executable,'-B',str(Path(__file__).resolve()),'--output-dir',str(output)],'cwd':str(Path.cwd())}
    sources=[Path(__file__),HERE/'verify_numerical.py',HERE/'audit.py',HERE/'NUMERICAL_TEST_CONTRACT.md',Path(executor.__file__),Path(q.__file__),Path(__file__).with_name('test_m3_sodium_ledger.py')]
    manifest['source_sha256']={str(p.resolve()):q.sha(p) for p in sources}
    stream=io.StringIO();suite=unittest.defaultTestLoader.loadTestsFromTestCase(NumericalControls)
    actual=unittest.TextTestRunner(stream=stream,verbosity=2,resultclass=Recorded).run(suite)
    (output/'unittest.log').write_text(stream.getvalue());rows=actual.rows
    try:
        fixture=output/'ANALYTICAL-SCORER-FIXTURE-NOT-NATIVE';fixture.mkdir()
        configuration=plan();runs=[write_level(fixture,configuration,steps,ranks,error) for ranks in (1,8) for steps,error in zip((8,16,32),(.0001,.00005,.00001))]
        diagnostic=numerical.verify(runs,configuration)
        diagnostic['scope']='GENERATED ANALYTICAL SCORER FIXTURE; NOT native/model/observational evidence'
        (fixture/'plan.json').write_text(json.dumps(configuration,indent=2)+'\n')
        (fixture/'scorer_result.json').write_text(json.dumps(diagnostic,indent=2)+'\n')
        with (fixture/'comparison.csv').open('w',newline='') as target:
            writer=csv.DictWriter(target,fieldnames=list(diagnostic['comparisons'][0]));writer.writeheader();writer.writerows(diagnostic['comparisons'])
        import matplotlib
        matplotlib.use('Agg');import matplotlib.pyplot as plt
        figure,axis=plt.subplots(figsize=(8,4))
        for ranks in ('1','8'):
            axis.plot([8,16,32],diagnostic['scores'][ranks]['relative_errors'],marker='o',label=ranks+' layout arithmetic fixture')
        axis.axhline(q.NUMERICAL_LIMIT['I21'],color='black',linestyle='--',label='unchanged approved I21 numerical allocation')
        axis.set_xlabel('Declared analytical refinement steps');axis.set_ylabel('Fixture normalized first-hit timing error');axis.legend(framealpha=1)
        figure.suptitle('SCORER SOFTWARE CONTROL — generated arithmetic fixture');figure.tight_layout()
        for suffix in ('png','eps'):figure.savefig(fixture/('comparison.'+suffix))
        plt.close(figure);rows.append({'name':'retained-analytical-firsthit-artifacts','status':diagnostic['status']})
    except q.MissingPrerequisite as error:rows.append({'name':'required-source-artifacts','status':'SKIPPED','reason':str(error)})
    except ImportError as error:rows.append({'name':'required-plot-dependency','status':'SKIPPED','reason':str(error)})
    except Exception as error:
        rows.append({'name':'analytical-artifact-infrastructure','status':'ERROR','reason':str(error)});(output/'artifact_exception.log').write_text(traceback.format_exc())
    counts={s:sum(r['status']==s for r in rows) for s in ('PASS','FAIL','ERROR','SKIPPED')}
    status=next(s for s in ('ERROR','FAIL','SKIPPED','PASS') if counts[s])
    result={'schema':'srcMoon-M3-sodium-numerical-controls-v1','status':status,'counts':counts,'results':rows,'scope':manifest['scope'],'formal_I21_I24_status':'SKIPPED / NOT QUALIFIED','elapsed_wall_seconds':time.monotonic()-started}
    manifest['source_changed_during_run']=[p for p,h in manifest['source_sha256'].items() if q.sha(p)!=h]
    if manifest['source_changed_during_run']:
        rows.append({'name':'source-identity-during-run','status':'ERROR'});result['status']='ERROR';result['counts']['ERROR']+=1;result['source_identity_note']='parent-owned implementation changed during review; retain this attempt and rerun on freeze'
    (output/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');(output/'result.json').write_text(json.dumps(result,indent=2)+'\n');return result
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output-dir',required=True,type=Path);args=parser.parse_args()
    if args.output_dir.exists():parser.error('output must not exist')
    r=run(args.output_dir);print(json.dumps(r,indent=2));raise SystemExit({'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[r['status']])
