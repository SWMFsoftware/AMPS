#!/usr/bin/env python3
"""Build and compare distinct srcSEP and srcSEP3D production-core drivers.

The executable hashes and source paths are evidence.  The gate rejects a
single binary used for both sides, preventing a shared-helper self-comparison.
Statistical tolerances are declared here before execution and cover Parker
mean/variance plus focused pitch mean/variance under matched SI inputs.
"""
import argparse, hashlib, json, math, pathlib, subprocess, sys

ROOT=pathlib.Path(__file__).resolve().parents[1]
def run(cmd): subprocess.run(cmd,check=True,cwd=ROOT)
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser(); p.add_argument('--sep1d-root',type=pathlib.Path,required=True,help='explicit srcSEP root; never inferred from a sibling'); p.add_argument('--output-dir',type=pathlib.Path,required=True); a=p.parse_args()
 out=a.output_dir.resolve(); out.mkdir(parents=True,exist_ok=True); cxx='g++'; sep1d=a.sep1d_root.resolve(); amps=sep1d.parent
 if not (sep1d/'util/sep_parker_core.cpp').is_file(): raise RuntimeError('--sep1d-root is not a srcSEP source tree')
 common=amps/'src/models/sep_common/sep_common.a'
 run(['make','-C',str(common.parent)])
 one=out/'sep1d-parity'; three=out/'sep3d-parity'
 run([cxx,'-std=c++17','-O2','-I'+str(sep1d/'util'),'-I'+str(amps/'src/models/sep_common'),
      str(ROOT/'test/parity/sep1d_parity_driver.cpp'),str(sep1d/'util/sep_parker_core.cpp'),
      str(sep1d/'util/sep_focused_transport_core.cpp'),str(common),'-o',str(one)])
 run([cxx,'-std=c++17','-O2','-I'+str(ROOT/'transport'),'-I'+str(ROOT/'core'),
      str(ROOT/'test/parity/sep3d_parity_driver.cpp'),str(ROOT/'transport/keyed_random.cpp'),
      str(ROOT/'transport/perpendicular_transport.cpp'),str(ROOT/'transport/parker_transport.cpp'),
      str(ROOT/'transport/focused_transport.cpp'),'-o',str(three)])
 if one.resolve()==three.resolve() or sha(one)==sha(three): raise RuntimeError('parity producers are not distinct')
 j1=out/'srcsep.json'; j3=out/'srcsep3d.json'; run([str(one),str(j1)]); run([str(three),str(j3)])
 x=json.loads(j1.read_text()); y=json.loads(j3.read_text()); limits={'parker_mean_m':0.06,'parker_variance_m2':0.06,'focused_mu_mean':0.02,'focused_mu_variance':0.06}
 # The symmetric pitch distribution has exact mean zero, so a relative error
 # is undefined/noise-amplifying there. Its declared acceptance is absolute;
 # the positive Parker moments and pitch variance use relative error.
 errors={k:(abs(x[k]-y[k]) if k=='focused_mu_mean' else abs(x[k]-y[k])/max(abs(x[k]),abs(y[k]),1e-12)) for k in limits}
 status='PASS' if all(errors[k]<=limits[k] for k in limits) else 'FAIL'
 report={'schema':'srcsep3d-cross-model-parity-v1','status':status,'executables':{'srcSEP':sha(one),'srcSEP3D':sha(three)},'records':{'srcSEP':sha(j1),'srcSEP3D':sha(j3)},'relative_errors':errors,'limits':limits}
 (out/'parity-report.json').write_text(json.dumps(report,indent=2)+'\n')
 print(f"V02 {status}: distinct srcSEP/srcSEP3D Parker and focused records")
 return 0 if status=='PASS' else 1
if __name__=='__main__': sys.exit(main())
