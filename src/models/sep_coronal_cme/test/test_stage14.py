#!/usr/bin/env python3
"""Manufactured verification for the thirteen non-release Stage-14 IDs.

These tests call public kernels and immutable offline producers. EVT/XMD/SLM
verify their campaign protocols on synthetic data; they do not certify an
observed campaign. Full AMPS/MPI adapter, moving-frame, and event-data gates
remain separately required and cannot be inferred from this software result.
"""
from __future__ import annotations
import argparse
import copy
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from preprocessing.core import PreprocessingError, digest, freeze_record, verify_frozen
from preprocessing.research import (CAPABILITIES,research_configuration,WAVE_COMMON_BINDINGS,WAVE_COUPLED_BINDINGS,
    validate_wave_asset,advance_resonant_waves,resonant_diffusion,streaming_limit,family_bundle,transfer_campaign)
from preprocessing.research_products import (kinetic,renewal_checkpoint,renew,drift_exposure,read_family,publish_family)
from preprocessing.wind_envelope import gate_wind,certified_extrema,outward
from preprocessing.impulsive import impulsive_plan
from preprocessing.nonradial import SteadyPiolaProvider,norm,dot,cross,plus,scaled,matvec,inverse
from preprocessing.protocols import (freeze_transfer_protocol,bind_transfer_run,TRANSFER_AUTHORITIES,
    compare_backgrounds,MATCHED_FIELDS,BOUNDARY_FIELDS,VARIABLES)
IDS=['CPL3D10','TUR3D08','MFP3D08','SLM3D01','FTE3D10','FTE3D11','ELL3D11',
     'SRC3D20','SRC3D21','SRC3D22','WND3D21','EVT3D01','XMD3D01']


def refreeze(value):
    out=copy.deepcopy(value);out.pop('identity',None);return freeze_record(out)


def cpp(identifier):
    result=subprocess.run([str(ROOT/'build/sep_coronal_cme_tests'),'--test',identifier],
        stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
    if result.returncode: raise AssertionError(result.stdout)


class CPL3D10(unittest.TestCase):
    def test_characteristic_piola_flux_interfaces_and_inverse(self):
        axis=scaled([1.,2.,3.],1/math.sqrt(14))
        def velocity(x):
            # u_r=1. This smooth steady flow has an independently known map:
            # F=r exp([axis]x * .005*(r-1)^2) n. It deflects theta and phi.
            r=norm(x);return plus(scaled(x,1/r),cross(axis,x),.01*(r-1))
        def b0(x):return scaled(x,1/norm(x)**3)
        provider=SteadyPiolaProvider(1,5,.02,1e-4,velocity,b0,digest('steady-flow'),2e-6)
        reference=[2.,1.,1.];r=norm(reference);angle=.005*(r-1)**2
        exact=plus(plus(scaled(reference,math.cos(angle)),cross(axis,reference),math.sin(angle)),
                   scaled(axis,dot(axis,reference)),1-math.cos(angle))
        value=provider.evaluate(reference,3,-1,[0,0,1])
        self.assertLess(norm(plus(exact,value['position_m'],-1)),1e-8)
        self.assertGreater(value['jacobian'],0);self.assertEqual(value['sector'],-1)
        self.assertAlmostEqual(value['mass_density_kg_m3']*norm(value['velocity_m_per_s'])/norm(value['magnetic_field_t']),3,places=12)
        inverse_point=provider.inverse_position(exact,reference,1e-10)
        self.assertLess(norm(plus(inverse_point,reference,-1)),1e-9)
        # Cofactor surface-flux identity, evaluated with independently formed
        # tangent area vectors. This also tests inverse-transpose event normals.
        a,b=[1.,0.,0.],[0.,1.,0.]
        area=cross(matvec(value['deformation'],a),matvec(value['deformation'],b))
        self.assertAlmostEqual(dot(value['magnetic_field_t'],area),dot(b0(reference),cross(a,b)),places=10)
        self.assertLess(norm(cross(value['interface_normal'],area)),1e-10)
        h=2e-4;divergence=0.;grad=[]
        for j in range(3):
            step=[0.,0.,0.];step[j]=h
            xp,xm=plus(exact,step),plus(exact,step,-1)
            ap=provider.inverse_position(xp,reference,1e-10);am=provider.inverse_position(xm,reference,1e-10)
            bp=provider.evaluate(ap,3,1)['magnetic_field_t'];bm=provider.evaluate(am,3,1)['magnetic_field_t']
            divergence+=(bp[j]-bm[j])/(2*h);grad.append((norm(bp)-norm(bm))/(2*h))
        self.assertLess(abs(divergence),2e-6)
        # Line evaluation and 3-D inverse evaluation agree at a physical point;
        # focusing uses the transformed gradient, never a radial-only stencil.
        line=provider.evaluate(reference,3,1);volume=provider.evaluate(inverse_point,3,1)
        self.assertLess(norm(plus(line['magnetic_field_t'],volume['magnetic_field_t'],-1)),1e-9)
        focusing=-dot(scaled(line['magnetic_field_t'],1/norm(line['magnetic_field_t'])),grad)/norm(line['magnetic_field_t'])
        self.assertTrue(math.isfinite(focusing));self.assertGreater(abs(focusing),.1)
        radial=SteadyPiolaProvider(1,5,.03,1e-4,lambda x:scaled(x,1/norm(x)),b0,digest('radial'),1e-6)
        reduced=radial.evaluate(reference,3,1)
        self.assertLess(norm(plus(reduced['position_m'],reference,-1)),1e-12)
        self.assertLess(norm(plus(reduced['magnetic_field_t'],b0(reference),-1)),1e-9)
        with self.assertRaises(PreprocessingError): inverse([[1,0,0],[0,1,0],[0,0,-1]])
        turning=SteadyPiolaProvider(1,5,.03,1e-4,lambda x:scaled(x,-1),b0,digest('turning'),1e-6)
        with self.assertRaises(PreprocessingError):turning.forward(reference)
        with self.assertRaises(PreprocessingError):provider.forward([6,0,0])
        mismatch=SteadyPiolaProvider(1,5,.03,1e-4,velocity,b0,digest('badtrace'),1e-6,lambda x:[2,0,0])
        with self.assertRaises(PreprocessingError):mismatch.evaluate(reference,3,1)


def wave_controls():
    return dict(ds_m=1.,dlogk=.5,alfven_speed_m_per_s=1.,tube_area_m2=2.,
        front_coordinate_support_m=[2,6],reference_distance_m=1.,shock_propagation_approximation='frozen-shock-over-coupling-step',
        flow_speed_m_per_s=0.,logk_advection_per_s=0.,cascade_per_s=0.,damping_per_s=0.,focusing_per_s=0.,
        gyrofrequency_per_s=1.,cr_to_ion_density=.01,particle_mass_density=.001,growth_enabled=True,
        spatial_boundary='periodic',spectral_boundary='outflow',conservation_relative_tolerance=1e-10)


def wave_asset():
    bindings={key:digest(key) for key in WAVE_COMMON_BINDINGS|WAVE_COUPLED_BINDINGS}
    bindings.update(background_generation=3,shock_history_generation=7,shock_frame='inertial',
        wave_frame_signed_direction='inertial:+/-along-oriented-B',wave_number_sign_convention='signed-k-in-oriented-B',external_wave_field='not-applicable')
    bindings['table_content']=digest([1.,2.]);bindings['coverage_mask']=digest([True,True])
    return freeze_record(dict(schema='sccm-foreshock-wave-asset-v6',bindings=bindings,
        coupling_kind='self-consistent-particle-wave-iteration',coverage_mask=[True,True],spectrum=[1.,2.]))


class TUR3D08(unittest.TestCase):
    def test_growth_balances_sinks_positivity_and_replay(self):
        state=dict(wave_energy=[[[.0001,.0001] for _ in range(4)] for _ in range(4)],particle_energy=[1.]*4,particle_momentum=[.002]*4)
        controls=wave_controls();dt=.1;result=advance_resonant_waves(state,controls,dt)
        gamma=math.pi/4*.01 # v_stream/v_A=2, independently evaluated growth.
        self.assertAlmostEqual(result['wave_energy'][0][0][0],.0001*(1+2*gamma*dt),places=15)
        self.assertAlmostEqual(result['wave_energy'][0][0][1],.0001,places=15)
        self.assertGreater(result['ledger']['particle_wave_energy'],0)
        self.assertLess(abs(result['ledger']['energy_residual_j']),1e-12)
        self.assertLess(abs(result['ledger']['momentum_residual_kg_m_per_s']),1e-12)
        # Separately activate every explicit sink/background-work channel.
        controls.update(logk_advection_per_s=.05,cascade_per_s=.03,damping_per_s=.1,focusing_per_s=.02)
        result=advance_resonant_waves(state,controls,dt)
        self.assertGreater(result['ledger']['damping_heat'],0);self.assertGreater(result['ledger']['spectral_boundary_sink'],0)
        self.assertGreater(result['ledger']['background_work'],0)
        controls=wave_controls();controls['growth_enabled']=False
        self.assertEqual(advance_resonant_waves(state,controls,dt)['wave_energy'],state['wave_energy'])
        with self.assertRaises(PreprocessingError):advance_resonant_waves(state,controls,2)
        asset=wave_asset();active=copy.deepcopy(asset['bindings']);validate_wave_asset(asset,active)
        for key in WAVE_COMMON_BINDINGS|WAVE_COUPLED_BINDINGS:
            changed=copy.deepcopy(active);changed[key]='different'
            with self.assertRaises(PreprocessingError):validate_wave_asset(asset,changed)
        external=copy.deepcopy(asset);external['coupling_kind']='prescribed-external-wave-field'
        for key in WAVE_COUPLED_BINDINGS:external['bindings'][key]='not-applicable'
        external['bindings']['external_wave_field']=digest('external');external=refreeze(external)
        validate_wave_asset(external,external['bindings'])
        changed=copy.deepcopy(external);changed['bindings']['species_weights']=digest('invented')
        with self.assertRaises(PreprocessingError):validate_wave_asset(refreeze(changed),changed['bindings'])
        bad=copy.deepcopy(asset);bad['coverage_mask'][0]=False;bad['bindings']['coverage_mask']=digest(bad['coverage_mask'])
        with self.assertRaises(PreprocessingError):validate_wave_asset(refreeze(bad),bad['bindings'])

    def test_independent_advection_and_coupling_cadence_convergence(self):
        errors=[]
        for n in (16,32,64):
            ds=2*math.pi/n;dt=ds/4;steps=round(.8/dt);dt=.8/steps
            controls=wave_controls();controls.update(ds_m=ds,growth_enabled=False)
            state=dict(wave_energy=[[[1+.2*math.sin(j*ds),1.] for _ in range(4)] for j in range(n)],particle_energy=[1.]*n,particle_momentum=[0.]*n)
            for _ in range(steps):state=advance_resonant_waves(state,controls,dt)
            errors.append(math.sqrt(sum((state['wave_energy'][j][0][0]-(1+.2*math.sin(j*ds-.8)))**2 for j in range(n))/n))
        self.assertGreater(errors[0]/errors[1],1.6);self.assertGreater(errors[1]/errors[2],1.6)
        # Constant spectrum gives D=.5*nu*(1-mu²), hence lambda=v/nu.
        speed=2.;omega=3.;b=1e-3;power=1e-5;nu=math.pi/2*omega*(4*math.pi*1e-7)*power/(b*b)
        estimates=[]
        for n in (128,256):
            integral=math.fsum((1-(-1+(j+.5)*2/n)**2)**2/resonant_diffusion(-1+(j+.5)*2/n,speed,omega,b,[-5,5],[power,power],.1)*2/n for j in range(n))
            estimates.append(3*speed/8*integral)
        self.assertLess(abs(estimates[1]-speed/nu),abs(estimates[0]-speed/nu))
        self.assertAlmostEqual(estimates[1]/(speed/nu),1,places=4)
        with self.assertRaises(PreprocessingError):resonant_diffusion(.5,speed,omega,b,[-10,-9],[power,power],.1)


class MFP3D08(unittest.TestCase):
    def test_proxy_cpp_and_independent_capability_switches(self):
        cpp('MFP3D08')
        flags={name:False for name in CAPABILITIES};flags['foreshock-distance-proxy']=True
        authority={'foreshock-distance-proxy':dict(algorithm_version='compact-c2-v1',verification_gate='MFP3D08',domain='upstream-resolved',validation_claim='coefficient-sensitivity')}
        record=research_configuration(flags,authority)
        self.assertFalse(record['schema5_selector_available']);self.assertFalse(record['capability_flags']['self-generated-waves'])
        with self.assertRaises(PreprocessingError):research_configuration(flags,{})


class SLM3D01(unittest.TestCase):
    def test_response_radial_covariance_all_strata_and_typed_zeros(self):
        registration=freeze_record(dict(schema='sccm-streaming-limit-registration-v6',frozen_before_sep=True,
            source_asset_sha256=digest('limit'),calibration_asset_sha256=digest('calibration'),response_sha256=digest('response'),
            required_strata=['observer-A/proton/E1/mover1/cadence1/closure1','observer-B/proton/E1/mover1/cadence1/closure1'],
            response_weights=[.25,.75],units='m^-2 s^-1 sr^-1 J^-1',window=[0,10],radius_support_m=[1,3],
            radial_transform=dict(law='registered-power-law-sensitivity',exponent=-2.,target_radius_m=2.,limit_reference_radius_m=2.,authority_sha256=digest('radial-law'))))
        rows=[dict(stratum=s,units=registration['units'],window=[0,10],radius_m=1.,model=[8.,8.],limit=[1.,1.],
            variance_log_model=.01,variance_log_limit=.04,log_covariance=.005,state='valid') for s in registration['required_strata']]
        rows[1].update(model=[0.,0.],state='censored')
        result=streaming_limit(rows,registration)
        self.assertAlmostEqual(result['rows'][0]['q_SL'],math.log10(2),places=13)
        self.assertIsNone(result['rows'][1]['q_SL']);self.assertEqual(result['rows'][1]['state'],'censored')
        self.assertFalse(result['observational_agreement_gate']);self.assertFalse(result['rows'][0]['feedback_qualified'])
        self.assertGreater(result['rows'][0]['q_interval'][1],result['rows'][0]['q_SL'])
        with self.assertRaises(PreprocessingError):streaming_limit(rows[:1],registration)
        rows[0]['radius_m']=10
        with self.assertRaises(PreprocessingError):streaming_limit(rows,registration)


class FTE3D10(unittest.TestCase):
    def test_independent_full_orbit_reference(self):cpp('FTE3D10')


class FTE3D11(unittest.TestCase):
    def test_unified_ownership_disabled_recovery_and_exposure(self):
        cpp('FTE3D11')
        segments=[dict(velocity_m_per_s=[1.,0.,0.],dt_s=2.,frame='inertial',generation=1),
                  dict(velocity_m_per_s=[-1.,0.,0.],dt_s=2.,frame='inertial',generation=1)]
        widths=dict(tube_radius_m=4.,rms_1d_perpendicular_width_m=2.,rms_plane_perpendicular_width_m=math.sqrt(8),diffusion_registration_sha256=digest('disabled-diffusion'))
        result=drift_exposure(segments,widths,'inertial',1)
        self.assertEqual(result['net_magnitude_m'],0);self.assertEqual(result['path_exposure_m'],4)
        self.assertEqual(result['exposure_ratios']['tube_radius_m'],1)
        self.assertAlmostEqual(result['exposure_ratios']['rms_plane_perpendicular_width_m'],math.sqrt(2))
        segments[1]['generation']=2
        with self.assertRaises(PreprocessingError):drift_exposure(segments,widths,'inertial',1)


class ELL3D11(unittest.TestCase):
    def test_translation_rotation_shape_normal_speed(self):cpp('ELL3D11')


def impulsive_asset():
    return freeze_record(dict(schema='sccm-impulsive-source-v6',intent='attribution-sensitivity',
        observation_response_sha256=digest('response'),mapping_covariance_sha256=digest('mapping'),frame='inertial',
        random_stream='impulsive-only',normalization_authority='independent-impulsive',
        footprint=[dict(id='open-A',open_mapping=True,supported=True,probability=.4),dict(id='open-B',open_mapping=True,supported=True,probability=.6)],
        time_response=[dict(time_s=0.,probability=.25),dict(time_s=1.,probability=.75)],
        momentum_spectrum=[dict(momentum_si=1e-20,probability=1.)],pitch_distribution=[dict(mu=.5,probability=1.)],
        species_mixture=[dict(species_id='proton',mass_kg=1.6726219e-27,probability=1.)]))


class SRC3D20(unittest.TestCase):
    def test_disjoint_source_measure_and_front_encounters(self):
        cpp('SRC3D20');asset=impulsive_asset();plan=impulsive_plan(asset,100,7,'upstream')
        self.assertEqual(len(plan['births']),4);self.assertAlmostEqual(sum(b['number'] for b in plan['births']),100)
        self.assertAlmostEqual(plan['births'][0]['number'],10);self.assertGreater(plan['births'][0]['kinetic_energy_j'],0)
        self.assertEqual(plan['shock_first_passage_number'],0);self.assertFalse(plan['shock_calibration_mutated'])
        with self.assertRaises(PreprocessingError):impulsive_plan(asset,100,7,'downstream')
        bad=copy.deepcopy(asset);bad['footprint'][0]['open_mapping']=False
        with self.assertRaises(PreprocessingError):impulsive_plan(refreeze(bad),100,7,'upstream')
        bad=copy.deepcopy(asset);bad['random_stream']='shock-source'
        with self.assertRaises(PreprocessingError):impulsive_plan(refreeze(bad),100,7,'upstream')


def family_fixture():
    calibration=freeze_record(dict(schema='sccm-reference-calibration-v6',authorities={key:digest(key) for key in
        ('coefficient_model','mean_free_path_authority','mover','return_policy','calibration_horizon','geometry')}))
    common=dict(species_id='proton',patch_id=1,time_lower_s=0.,time_upper_s=10.,geometry_generation=7,
        source_is_separable=False,integrated_peclet=2.,joint_measure=3.,joint_measure_derivation_sha256=digest('normal-coarea-source'),
        geometry_certificate={key:True for key in ('positive_jacobian','unique_normal_root','reach','overlap','mask','solar_clearance')})
    members=[dict(common,momentum_lower_si=1.,momentum_upper_si=2.,offset_m=4.),dict(common,momentum_lower_si=2.,momentum_upper_si=3.,offset_m=2.)]
    return calibration,members


class SRC3D21(unittest.TestCase):
    def test_joint_family_bundle_major_support_and_restart(self):
        cpp('SRC3D21');calibration,members=family_fixture();bundle=family_bundle(members,calibration,2.)
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'family.json';publish_family(bundle,path)
            self.assertEqual(read_family(path),bundle)
            with self.assertRaises(PreprocessingError):read_family(path,supported_major=3)
            with self.assertRaises(PreprocessingError):publish_family(bundle,path)
        for key,value in [('momentum_lower_si',2.1),('geometry_generation',8),('source_is_separable',True)]:
            bad=copy.deepcopy(members);bad[1][key]=value
            with self.assertRaises(PreprocessingError):family_bundle(bad,calibration,2.)
        bad=copy.deepcopy(members);bad[1]['geometry_certificate']['solar_clearance']=False
        with self.assertRaises(PreprocessingError):family_bundle(bad,calibration,2.)


def renewal_fixture(number=10.):
    mass=1.6726219e-27;incoming=[1e-20,0.,0.];outgoing=[2e-20,0.,0.]
    cohort=dict(number=number,mass_kg=mass,momentum_si=incoming,cycle=0,time_s=1.,ancestry_id='firstpass-1',species_id='proton',frame_id='inertial')
    condition=dict(momentum_si=incoming,pitch_cosine=.5,surface_position_m=[1.,0.,0.],time_s=1.,species_id='proton',frame_id='inertial',
        wave_authority=digest('ambient'),front_generation=7,reference_family_identity=digest('family'))
    branches=[dict(outcome='Absorbed',probability=.3,residence_time_s=0.,momentum_si=incoming,shock_work_j=0.,shock_impulse_si=[0.,0.,0.]),
        dict(outcome='Downstream',probability=.2,residence_time_s=0.,momentum_si=incoming,shock_work_j=0.,shock_impulse_si=[0.,0.,0.]),
        dict(outcome='ReReleased',probability=.5,residence_time_s=3.,momentum_si=outgoing,shock_work_j=kinetic(mass,outgoing)-kinetic(mass,incoming),shock_impulse_si=[1e-20,0.,0.])]
    kernel=freeze_record(dict(schema='sccm-conditional-renewal-kernel-v6',validated_work_authority=True,work_validation_sha256=digest('conditional-work'),conditioning=condition,branches=branches,relative_tolerance=1e-10))
    first=freeze_record(dict(schema='immutable-first-passage-v1',number=10.,momentum_spectrum_sha256=digest('original-g')))
    return first,kernel,cohort


class SRC3D22(unittest.TestCase):
    def test_conditional_renewal_ancestry_split_restart_and_work(self):
        cpp('SRC3D22');first,kernel,cohort=renewal_fixture();checkpoint=renewal_checkpoint(first,kernel,cohort);result=renew(checkpoint)
        self.assertEqual(result['ledger']['outcome_number'],dict(Absorbed=3.,Downstream=2.,ReReleased=5.))
        self.assertEqual(result['continuations'][0]['cycle'],1);self.assertEqual(result['continuations'][0]['time_s'],4)
        self.assertEqual(result['continuations'][0]['ancestry_id'],'firstpass-1');self.assertEqual(result['first_passage_identity'],first['identity'])
        # JSON restart is content identical; deterministic cohort splitting closes
        # before/after a repartition, without rank-local random draws.
        self.assertEqual(renew(json.loads(json.dumps(checkpoint))),result)
        halves=[renew(renewal_checkpoint(first,kernel,dict(cohort,number=5.))) for _ in range(2)]
        self.assertAlmostEqual(sum(r['ledger']['shock_work_j'] for r in halves),result['ledger']['shock_work_j'],places=25)
        self.assertEqual(sum(r['ledger']['outcome_number']['ReReleased'] for r in halves),5)
        # A new cycle must condition on the outgoing state, never original g(p).
        next_cohort=result['continuations'][0];next_kernel=copy.deepcopy(kernel)
        next_kernel['conditioning'].update(momentum_si=next_cohort['momentum_si'],time_s=next_cohort['time_s'])
        next_kernel['branches']=[dict(kernel['branches'][0],probability=1.,momentum_si=next_cohort['momentum_si'])]
        second=renew(renewal_checkpoint(first,refreeze(next_kernel),next_cohort,result['transactions']))
        self.assertEqual(len(second['transactions']),2);self.assertEqual(second['ledger']['cycle'],1)
        self.assertEqual(second['first_passage_identity'],first['identity'])
        bad=copy.deepcopy(kernel);bad['branches'][2]['shock_work_j']*=2
        with self.assertRaises(PreprocessingError):renew(renewal_checkpoint(first,refreeze(bad),cohort))
        bad=copy.deepcopy(kernel);bad['branches'][2]['residence_time_s']=0
        with self.assertRaises(PreprocessingError):renew(renewal_checkpoint(first,refreeze(bad),cohort))
        with self.assertRaises(PreprocessingError):renew(renewal_checkpoint(first,kernel,dict(cohort,momentum_si=[2e-20,0.,0.])))


def wind_fixture(coefficients=(2.,1.),speed_bounds=(1.,4.),accel_bounds=(-1.,4.)):
    channels=[]
    for identifier,quantity,units,bounds in [('velocity','field-aligned-speed','m/s',speed_bounds),('acceleration','quasi-steady-field-aligned-advective-acceleration','m/s^2',accel_bounds)]:
        channels.append(dict(id=identifier,quantity=quantity,units=units,frame='inertial',support=[0.,1.],interpolant='polynomial-in-s-SI-v1',lower_bound=bounds[0],upper_bound=bounds[1],enclosure_tolerance_si=1e-5))
    asset=freeze_record(dict(schema='sccm-wind-envelope-v6',data_use_role='qualification',covariance_sha256=digest('covariance'),inference_provenance_sha256=digest('wind-inference'),channels=channels,
        gates=[dict(selector='open:tube-A',velocity_channel='velocity',acceleration_channel='acceleration')]))
    profile=dict(selector='open:tube-A',quasi_steady=True,frame='inertial',support=[0.,1.],speed_coefficients=list(coefficients),tube_id='A',segment_id='A:0',unsigned_flux_wb=3.)
    return asset,[profile]


class WND3D21(unittest.TestCase):
    def test_continuous_paired_channel_extrema_and_hidden_overshoot(self):
        asset,profiles=wind_fixture();self.assertTrue(gate_wind(asset,profiles)['passed'])
        # Endpoints both u=1 and u*u'=0, yet the interior violates speed and
        # advective-acceleration bounds. A node-only check would miss both.
        asset,profiles=wind_fixture([1,0,16,-32,16],(.5,1.5),(-.5,.5))
        result=gate_wind(asset,profiles);self.assertFalse(result['passed']);row=result['rows'][0]
        self.assertEqual(row['rejected_magnetic_flux_wb'],3)
        self.assertGreater(row['channels'][0]['upper_excursion_si'],0)
        self.assertGreater(row['channels'][1]['upper_excursion_si'],0)
        cert=certified_extrema([0,4,-4],[0,1],1e-5)
        self.assertLessEqual(cert['maximum_interval'][0],1);self.assertGreaterEqual(cert['maximum_interval'][1],1)
        self.assertEqual(outward(1.,math.inf),1.0000000000000002)
        profiles[0]['quasi_steady']=False
        with self.assertRaises(PreprocessingError):gate_wind(asset,profiles)
        bad=copy.deepcopy(asset);bad['channels'][1]['quantity']='field-aligned-speed'
        with self.assertRaises(PreprocessingError):gate_wind(refreeze(bad),wind_fixture()[1])


def transfer_fixture():
    protocol=dict(schema='sep-event-transfer-protocol-v1',calibration_event='event-D',transfer_event='2020-05-29',held_out_group='2020-05-27/2020-06-02',
        frozen_authorities={key:'frozen-v1' for key in TRANSFER_AUTHORITIES},transferable_constants={'source_efficiency':1e-4},event_specific_inputs=['magnetogram','front'],
        freeze_before_withheld_sep=True,front_reconstruction_independent=True,classification='primary-transfer',
        observer_paths=[dict(observer_id=x,icme_screening_asset_sha256=digest(x+'screening'),background_coverage_complete=True,sep_coverage_complete=True,
            prior_magnetic_cloud=False,cloud_background_validated=False) for x in ('PSP','STEREO-A')])
    protocol=freeze_transfer_protocol(protocol)
    event={key:dict(asset_id=key+'-E',content_sha256=digest(key+'-E'),inference_procedure='frozen-v1') for key in ('magnetogram','front')}
    release={key:protocol[key] for key in ('transferable_constants','frozen_authorities')}
    run=bind_transfer_run(protocol,event,release)
    metrics=dict(onset_time=1.,peak_time=2.,log_intensity=3.,fluence_ratio=4.,spectral_index=2.,anisotropy=.3,uncertainty_coverage=.9)
    obs=freeze_record(dict(data_use_role='withheld-validation',held_out_group=protocol['held_out_group'],pure_radial_same_flux_tube_claim=False,
        rows=[dict(observer_id=x,metrics=metrics) for x in ('PSP','STEREO-A')]))
    registration=freeze_record(dict(protocol_identity=protocol['identity'],frozen_before_observations=True,calibration_event_fingerprint=digest('D-run'),tolerances={k:.1 for k in metrics}))
    return protocol,run,obs,[dict(observer_id=r['observer_id'],metrics=copy.deepcopy(metrics)) for r in obs['rows']],registration,event,release


class EVT3D01(unittest.TestCase):
    def test_frozen_transfer_all_metrics_group_holdout_and_failed_values(self):
        p,r,o,pred,reg,event,release=transfer_fixture();result=transfer_campaign(p,r,o,pred,reg)
        self.assertTrue(result['passed']);self.assertEqual(result['independent_event_count'],1)
        pred[0]['metrics']['peak_time']=100
        failed=transfer_campaign(p,r,o,pred,reg);self.assertFalse(failed['passed'])
        self.assertEqual(len(failed['rows'][0]['metrics']),7);self.assertEqual(failed['rows'][0]['metrics']['peak_time']['value'],98)
        changed=copy.deepcopy(event);changed['front']['content_sha256']=digest('new-front')
        self.assertNotEqual(bind_transfer_run(p,changed,release)['release_calibration_fingerprint'],r['release_calibration_fingerprint'])
        bad=copy.deepcopy(o);bad['pure_radial_same_flux_tube_claim']=True
        with self.assertRaises(PreprocessingError):transfer_campaign(p,r,refreeze(bad),pred,reg)
        bad=copy.deepcopy(reg);bad['frozen_before_observations']=False
        with self.assertRaises(PreprocessingError):transfer_campaign(p,r,o,pred,refreeze(bad))
        stress=copy.deepcopy(p);stress['transfer_event']='2012-05-17';stress.pop('identity')
        with self.assertRaises(PreprocessingError):freeze_transfer_protocol(stress)
        stress['classification']='stress-test';freeze_transfer_protocol(stress)


def mhd_fixture():
    metadata={key:'same' for key in MATCHED_FIELDS|BOUNDARY_FIELDS}
    metadata.update(epoch_utc='2000-01-01T00:00:00Z',frame='inertial',cadence_s=1.,spatial_support=dict(minimum_m=[0,-1,-1],maximum_m=[2,1,1],time_s=[0,1]),masks=[True],interpolation='none-analytic-at-nodes',
        magnetogram_sha256=digest('mag'),preprocessing_sha256=digest('pre'),front_history_sha256=digest('front'),
        magnetic_normalization=dict(value=1.,units='T',definition='inner-B'),open_flux_normalization=dict(value=1.,units='Wb',definition='unsigned-open-flux'),
        variable_definitions={k:'synthetic-'+k for k in VARIABLES},units=dict(number_density_m3='m^-3',mass_density_kg_m3='kg/m^3',temperature_K='K',pressure_Pa='Pa',
        magnetic_field_T='T',velocity_m_per_s='m/s',alfven_speed_m_per_s='m/s',fast_speed_m_per_s='m/s',signed_open_flux_wb='Wb',unsigned_open_flux_wb='Wb',D2_first_fast_height_m='m',D2_first_supercritical_height_m='m'))
    variables={k:dict(values=[1.,2.,3.] if k in {'magnetic_field_T','velocity_m_per_s'} else [1.],
        covariance=[[.01 if i==j else 0. for j in range(3 if k in {'magnetic_field_T','velocity_m_per_s'} else 1)] for i in range(3 if k in {'magnetic_field_T','velocity_m_per_s'} else 1)]) for k in VARIABLES}
    return freeze_record(dict(schema='sep-offline-background-samples-v1',metadata=metadata,variables=variables,sample_coordinates=[[1,0,0,0]],source_uri='synthetic://background',source_sha256=digest('synthetic-background'),model_version='v1',run_id='manufactured',residual_tuning_used=False,
        topology=dict(definition='open-closed',values=['open']),connectivity=dict(definition='stable-footpoint',values=['A'])))


class XMD3D01(unittest.TestCase):
    def test_offline_matched_samples_uncertainty_and_boundary_mismatch(self):
        analytic=mhd_fixture();mhd=copy.deepcopy(analytic);mhd['variables']['temperature_K']['values']=[1.1];mhd=refreeze(mhd)
        result=compare_backgrounds(analytic,mhd)
        self.assertAlmostEqual(result['metrics']['temperature_K']['rms_difference'],.1)
        self.assertFalse(result['truth_validation']);self.assertFalse(result['runtime_imported_provider_qualified'])
        self.assertEqual(len(result['metrics']),len(VARIABLES))
        bad=copy.deepcopy(mhd);bad['metadata']['magnetogram_sha256']=digest('other-mag');bad=refreeze(bad)
        with self.assertRaises(PreprocessingError):compare_backgrounds(analytic,bad)
        self.assertEqual(compare_backgrounds(analytic,bad,True)['classification'],'combined-boundary-plus-model-discrepancy')
        for mutation in ('masks','frame','residual-tuning'):
            bad=copy.deepcopy(mhd)
            if mutation=='masks':bad['metadata']['masks']=[False]
            elif mutation=='frame':bad['metadata']['frame']='other'
            else:bad['residual_tuning_used']=True
            with self.assertRaises(PreprocessingError):compare_backgrounds(analytic,refreeze(bad),True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--test',choices=IDS);args=parser.parse_args()
    classes=[globals()[args.test]] if args.test else [globals()[name] for name in IDS]
    suite=unittest.TestSuite(unittest.defaultTestLoader.loadTestsFromTestCase(case) for case in classes)
    success=unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful()
    # The specification reserves these IDs for actual campaign evidence. The
    # synthetic contract tests above are useful verification, but cannot make
    # a canonical observed-event/MHD/streaming campaign PASS. The owning runner
    # retains that verification separately and reports the missing campaign.
    campaign_ids={'EVT3D01','XMD3D01','SLM3D01'}
    if success:
        for identifier in ([args.test] if args.test else IDS):
            if identifier in campaign_ids:
                print('research_gate_status='+json.dumps(dict(id=identifier,status='SKIP',verification_passed=True,
                    reason='Synthetic contract checks passed; independently registered actual campaign evidence was not supplied.')))
    raise SystemExit(0 if success else 1)
