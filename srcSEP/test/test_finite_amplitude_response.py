"""Bounded independent math/protocol checks; no broad amplitude scan.

The Lorentz/Kubo normalization was already verified by RWN. These tests add
relativistic finite-time references and preserved unresolved/provenance paths.
They do not assert a Markov diffusion law or choose production domain limits.
"""
import copy
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
PHYSICS=ROOT/"srcSEP/validation/physics"
sys.path.insert(0,str(PHYSICS))
sys.dont_write_bytecode=True
import finite_amplitude_response as field
import run_finite_amplitude_study as runner
import plot_finite_amplitude_study as plotter


def control():
    value=json.loads((PHYSICS/"controlled_finite_amplitude_pilot.json").read_text())
    # Explicit short software/math campaign, separately recorded from pilot.
    value["provenance"]="bounded test-only numerical verification coordinates"
    value["amplitudes_delta_b_over_b0"]=["0", "0.0001"]
    value["observation_gyroperiods"]=["0.5"]
    value["ensemble"]["resonant_k_ratios"]=["1"]
    return value


class LorentzReferences(unittest.TestCase):
    def probe(self,helicity=1):
        # Dimensionless controlled q=m=c=B0=1, p=(.6,0,.4), signedu=.1.
        # This explicitly relativistic gamma~1.23 is not a solar-wind input.
        return field.Probe(1,1,1,1,.1,helicity,(.6,0,.4))

    def test_uniform_field_exact_gyro_and_fourth_order(self):
        p=self.probe();k=np.array([field.resonance_scale_per_m(p)])
        errors=[]
        for steps in (32,64,128):
            orbit=field.integrate_rk4(p,k,0,4,steps,["1"],10000)
            final=orbit["states"][0][0,1:]
            errors.append(float(np.linalg.norm(final-np.array(p.momentum))))
        self.assertLess(errors[2],errors[1]);self.assertGreater(math.log(errors[1]/errors[2],2),3.8)
        self.assertLess(errors[2],2e-6)

    def test_analytic_force_integral_resonance_and_wrong_helicity(self):
        p=self.probe();k=np.array([field.resonance_scale_per_m(p)]);a=1e-5;tau=math.pi*p.gamma
        phases=2*math.pi*np.arange(8)/8
        actual=field.linear_pitch(p,k,a,tau,phases)[0]
        beta=p.momentum_norm/p.gamma;mu=p.mean_pitch
        # At exact resonance sinc=1; independently integrate sin(phi) in time.
        coefficient=-(1/p.gamma)*math.sqrt(1-mu*mu)*(1-p.beta_u*mu/beta)*a*tau
        # The nominal resonant k is binary64, so its actual detuning can be
        # nonzero at roundoff. Pointwise relative error is undefined at the
        # sine zeros; normalize absolute error by the physical response
        # amplitude and retain the same2e-14 numerical accuracy criterion.
        self.assertLess(float(np.max(np.abs(actual-coefficient*np.sin(phases))))/abs(coefficient),2e-14)
        wrong=self.probe(-1)
        off=field.linear_pitch(wrong,k,a,tau,phases)
        self.assertLess(float(np.max(np.abs(off))),1e-18)

    def test_weak_full_orbit_recovers_finite_time_linear_response(self):
        p=self.probe();k=np.array([field.resonance_scale_per_m(p)]);weights=np.array([1.0]);a=1e-5
        orbit=field.integrate_rk4(p,k,a,16,256,["0.5"],10000)
        values=field.moments(p,orbit,k,weights,a)[0]
        theory=float(field.analytic_variance(p,k,weights,a,orbit["times"][0],60))
        self.assertLess(abs(values["pitch_change_variance"]/theory-1),2e-7)
        self.assertLess(values["maximum_sampled_invariant_residual_over_initial_kinetic"],1e-9)

    def test_finite_amplitude_dop853_and_exact_energy_momentum_invariant(self):
        p=self.probe();k=np.array([field.resonance_scale_per_m(p)]);a=.1
        orbit=field.integrate_rk4(p,k,a,8,512,["0.5","1"],10000)
        cross=field.dop853_comparison(p,k,a,orbit,1e-12,1e-14,4)
        self.assertLess(cross["maximum_absolute_pitch_difference"],1e-8)
        values=field.moments(p,orbit,k,np.array([1.0]),a)
        self.assertLess(max(x["maximum_sampled_invariant_residual_over_initial_kinetic"] for x in values),1e-9)
        self.assertGreater(values[-1]["rms_pitch_linear_residual"],1e-3)

    def test_phase_and_mode_quadrature_are_separate_ensembles(self):
        p=self.probe();scale=field.resonance_scale_per_m(p);a=1e-4;tau=2*math.pi*p.gamma
        continuous=field.analytic_variance(p,np.array([scale]),np.array([1.0]),a,tau,60,(.5*scale,2*scale))
        errors=[]
        for count in (4,8,16):
            nodes,weights=np.polynomial.legendre.leggauss(count)
            ks=1.25*scale+.75*scale*nodes
            discrete=field.analytic_variance(p,ks,weights/2,a,tau,60)
            errors.append(abs(float(discrete/continuous)-1))
        self.assertLess(errors[-1],1e-12);self.assertLess(errors[1],errors[0])
        # Quadratic first-order phase average is exact for N>=3. Full nonlinear
        # phase moments are independently refined by the study driver.
        for count in (4,8,16):
            values=field.linear_pitch(p,np.array([scale]),a,tau,2*math.pi*np.arange(count)/count)[0]
            theory=float(field.analytic_variance(p,np.array([scale]),np.array([1.0]),a,tau,60))
            self.assertLess(abs(float(np.mean(values*values))/theory-1),1e-13)

    def test_gap_has_no_k_floor_or_implicit_bridge(self):
        p=self.probe();equal_speed=field.Probe(1,1,1,1,p.momentum[2]/p.gamma,1,p.momentum)
        with self.assertRaises(field.Unresolved):field.resonance_scale_per_m(equal_speed)
        with self.assertRaises(ValueError):field.integrate_rk4(p,np.array([0.0]),.1,8,64,["1"],1000)
        # Wrong-helicity finite-time force may be nonzero at a generic time;
        # its response is not interpreted as broadband gap diffusion.
        values=field.linear_pitch(self.probe(-1),np.array([field.resonance_scale_per_m(p)]),.01,1.0,np.array([.3]))
        self.assertGreater(abs(float(values[0,0])),0)

    def test_decimal_underflow_rejects_without_zero_substitution(self):
        self.assertEqual(field.explicit_number("0","amplitude"),0)
        with self.assertRaises(field.Unresolved):field.explicit_number("1e-400","amplitude")
        with self.assertRaises(field.Unresolved):field.explicit_number("-1e-400","signed input")

    def test_both_branches_charges_and_helicities(self):
        for charge in (-1,1):
            for speed in (-.1,.1):
                for helicity in (-1,1):
                    p=field.Probe(1,charge,1,1,speed,helicity,(.6,0,.4))
                    k=np.array([field.resonance_scale_per_m(p)]);a=1e-6
                    orbit=field.integrate_rk4(p,k,a,8,256,["0.25"],10000)
                    actual=field.moments(p,orbit,k,np.array([1.0]),a)[0]["pitch_change_variance"]
                    expected=float(field.analytic_variance(p,k,np.array([1.0]),a,orbit["times"][0],60))
                    self.assertLess(abs(actual/expected-1),5e-7)

    def test_zero_momentum_pitch_is_unresolved_without_normalization(self):
        p=self.probe();fake={"phases":np.arange(4),"times":[1.0],"states":[np.zeros((4,4))]}
        with self.assertRaisesRegex(field.Unresolved,"pitch is undefined"):
            field.moments(p,fake,np.array([1.0]),np.array([1.0]),.1)


class ProtocolChecks(unittest.TestCase):
    def prepare(self,directory,controls=None):
        path=Path(directory)/"controls.json";path.write_text(json.dumps(control() if controls is None else controls))
        return runner.prepare(PHYSICS/"controlled_ep_maxwellian_baseline.json",path,Path(directory)/"bundle")

    def test_prepare_is_explicit_and_archived_replay_is_portable(self):
        with tempfile.TemporaryDirectory() as directory:
            protocol=self.prepare(directory);record=json.loads(protocol.read_text())
            self.assertFalse((protocol.parent/"study.json").exists())
            self.assertFalse(record["physical_limits_selected"])
            command=[sys.executable,str(protocol.parent/record["sources"]["runner"]["archive"]),"--execute-protocol",str(protocol)]
            result=subprocess.run(command,cwd="/tmp",capture_output=True,text=True)
            self.assertIn(result.returncode,(0,1),result.stderr)
            self.assertEqual(result.stderr,"")
            self.assertTrue((protocol.parent/"study.json").exists(),result.stderr)
            study=json.loads((protocol.parent/"study.json").read_text())
            self.assertEqual(len(study["rows"]),2)
            zero=next(row for row in study["rows"] if row["amplitude_delta_b_over_b0"]=="0")
            self.assertEqual(zero["linear_phase_variance"],0)
            self.assertIsNone(zero["relative_variance_difference"])
            self.assertFalse(study["physical_production_qualified"])

    def test_tampered_source_and_coordinates_reject(self):
        with tempfile.TemporaryDirectory() as directory:
            path=self.prepare(directory);record=json.loads(path.read_text())
            core=path.parent/record["sources"]["core"]["archive"];core.write_text(core.read_text()+"\n# changed\n")
            with self.assertRaisesRegex(ValueError,"identity mismatch"):runner.read_protocol(path)
        with tempfile.TemporaryDirectory() as directory:
            path=self.prepare(directory);record=json.loads(path.read_text())
            record["resolved_inputs"]["actual_probe"]["field_t"]="0x1p+1";path.write_text(json.dumps(record))
            with self.assertRaisesRegex(ValueError,"coordinates differ"):runner.read_protocol(path)

    def test_unknown_physics_and_nonaligned_observation_grid_reject(self):
        with tempfile.TemporaryDirectory() as directory:
            value=control();value["alfven_speed_source"]="guessed"
            with self.assertRaises(ValueError):self.prepare(directory,value)
        with tempfile.TemporaryDirectory() as directory:
            value=control();value["amplitudes_delta_b_over_b0"]=["1e-400"]
            # The pinned core is loaded under its own module identity, so
            # its typed Unresolved has a distinct Python class object. Match
            # the shared RuntimeError contract and explicit underflow reason.
            with self.assertRaisesRegex(RuntimeError,"underflows binary64"):self.prepare(directory,value)
        with tempfile.TemporaryDirectory() as directory:
            value=control();value["observation_gyroperiods"]=["0.123456789"]
            with self.assertRaises(ValueError):self.prepare(directory,value)

    def test_resource_failure_preserves_complete_unresolved_coverage(self):
        with tempfile.TemporaryDirectory() as directory:
            value=control();value["observation_gyroperiods"]=["0.5","1"]
            value["numerical_controls"]["maximum_total_trajectory_updates"]=1
            path=self.prepare(directory,value);record=json.loads(path.read_text())
            result=subprocess.run([sys.executable,str(path.parent/record["sources"]["runner"]["archive"]),"--execute-protocol",str(path)],capture_output=True,text=True)
            self.assertEqual(result.returncode,1,result.stderr)
            study=json.loads((path.parent/"study.json").read_text())
            self.assertEqual(len(study["rows"]),4);self.assertEqual(study["unresolved_rows"],4)
            self.assertEqual(len({(r["amplitude_delta_b_over_b0"],r["observation_gyroperiods"]) for r in study["rows"]}),4)

    def test_nonzero_analytic_response_underflow_stays_unresolved(self):
        with tempfile.TemporaryDirectory() as directory:
            value=control();value["amplitudes_delta_b_over_b0"]=["1e-200"]
            path=self.prepare(directory,value);record=json.loads(path.read_text())
            result=subprocess.run([sys.executable,str(path.parent/record["sources"]["runner"]["archive"]),"--execute-protocol",str(path)],capture_output=True,text=True)
            self.assertEqual(result.returncode,1,result.stderr)
            row=json.loads((path.parent/"study.json").read_text())["rows"][0]
            self.assertEqual(row["status"],"UNRESOLVED")
            self.assertIn("underflows",row["unresolved_reason"])
            self.assertIsNone(row["linear_phase_variance"])
            self.assertIsNone(row["relative_variance_difference"])
            self.assertGreater(float(row["linear_phase_variance_high_precision"].split('e')[0]),0)

    def test_archived_plot_and_coordinate_binding(self):
        with tempfile.TemporaryDirectory() as directory:
            value=control();value["amplitudes_delta_b_over_b0"]=["0.0001"]
            path=self.prepare(directory,value);record=json.loads(path.read_text())
            calculated=subprocess.run([sys.executable,str(path.parent/record["sources"]["runner"]["archive"]),"--execute-protocol",str(path)],capture_output=True,text=True)
            self.assertIn(calculated.returncode,(0,1),calculated.stderr)
            prefix=Path(directory)/"figure"
            command=[sys.executable,str(path.parent/record["sources"]["plotter"]["archive"]),"--protocol",str(path),"--output-prefix",str(prefix)]
            plotted=subprocess.run(command,cwd="/tmp",capture_output=True,text=True)
            self.assertEqual(plotted.returncode,0,plotted.stderr)
            self.assertTrue(Path(str(prefix)+".png").exists());self.assertTrue(Path(str(prefix)+".pdf").exists())
            # Explicit newer rendering source may replot old scientific data.
            # It must archive/report its actual code, never claim the old
            # protocol plotter identity. This performs no new orbit solve.
            newer=Path(directory)/"newer_plotter.py"
            newer.write_text((path.parent/record["sources"]["plotter"]["archive"]).read_text()+"\n# rendering provenance negative control\n")
            upgraded=subprocess.run([sys.executable,str(newer),"--protocol",str(path),"--output-prefix",str(prefix)],cwd="/tmp",capture_output=True,text=True)
            self.assertEqual(upgraded.returncode,0,upgraded.stderr)
            metadata=json.loads(Path(str(prefix)+".plot.json").read_text())
            self.assertFalse(metadata["uses_original_protocol_plotter"])
            self.assertNotEqual(metadata["plotter_sha256"],metadata["original_protocol_plotter_sha256"])
            self.assertEqual(runner.digest((prefix.parent/metadata["actual_plotter_archive"]).read_bytes()),metadata["plotter_sha256"])
            study_path=path.parent/"study.json";study=json.loads(study_path.read_text())
            study["rows"][0]["actual_k_per_m_hex"]=[float(123).hex()];study_path.write_text(json.dumps(study))
            changed=subprocess.run(command,cwd="/tmp",capture_output=True,text=True)
            self.assertEqual(changed.returncode,2);self.assertIn("k/weights differ",changed.stderr)

    def test_plot_positive_underflow_rejects(self):
        with self.assertRaisesRegex(ValueError,"underflow"):plotter.number("1e-400","positive finite response")


if __name__=="__main__":unittest.main()
