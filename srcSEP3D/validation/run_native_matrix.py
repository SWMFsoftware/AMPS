#!/usr/bin/env python3
"""Execute V03 against a real linked AMPS executable.

This runner never simulates MPI. It launches the production binary for every
declared rank/thread pair, records the exact argv/environment/executable hash,
and preserves each native JSON report. Absence of MPI or the executable is an
error because this command is the qualification gate, not a source-only test.
"""
import argparse, hashlib, json, os, pathlib, shlex, subprocess, sys, time
ROOT=pathlib.Path(__file__).resolve().parents[1]
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser(); p.add_argument('--amps',type=pathlib.Path,required=True); p.add_argument('--profile',choices=['small','medium','production'],default='small'); p.add_argument('--launcher',default='mpiexec -n {ranks}'); p.add_argument('--output-dir',type=pathlib.Path,required=True); p.add_argument('--timeout',type=float,default=3600); a=p.parse_args()
 exe=a.amps.resolve(); out=a.output_dir.resolve(); out.mkdir(parents=True,exist_ok=True)
 if not exe.is_file() or not os.access(exe,os.X_OK): raise RuntimeError('linked AMPS executable is absent or not executable')
 profiles=json.loads((ROOT/'validation/native_profiles.json').read_text())['profiles']; profile=next(x for x in profiles if x['name']==a.profile)
 records=[]
 for ranks in profile['ranks']:
  for threads in profile['threads']:
   for case in profile['required_cases']:
    directory=out/f'r{ranks}-t{threads}'/case; directory.mkdir(parents=True,exist_ok=True); report=directory/'native.json'
    cmd=[part.format(ranks=ranks) for part in shlex.split(a.launcher)]+[str(exe),'--test',case,'--test-json',str(report),'--artifact-directory',str(directory)]
    env=dict(os.environ); env['OMP_NUM_THREADS']=str(threads); start=time.monotonic()
    completed=subprocess.run(cmd,cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=a.timeout,env=env)
    status='ERROR'; message='native JSON report absent'
    if report.is_file():
     payload=json.loads(report.read_text()); row=next((r for r in payload.get('results',[]) if r.get('id','').upper()==case),{}); status=row.get('status','ERROR').upper(); message=row.get('message','')
    records.append({'case':case,'ranks':ranks,'threads':threads,'status':status,'message':message,'returncode':completed.returncode,'elapsed_seconds':time.monotonic()-start,'command':cmd,'stdout_tail':completed.stdout[-2000:]})
 summary={'schema':'srcsep3d-native-matrix-v1','profile':a.profile,'executable':str(exe),'executable_sha256':digest(exe),'records':records}
 (out/'native-matrix.json').write_text(json.dumps(summary,indent=2)+'\n')
 return 0 if all(r['status']=='PASS' and r['returncode']==0 for r in records) else 1
if __name__=='__main__':
 try: sys.exit(main())
 except Exception as e: print('ERROR:',e,file=sys.stderr); sys.exit(2)
