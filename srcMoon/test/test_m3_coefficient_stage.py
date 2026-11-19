#!/usr/bin/env python3
"""Exact native zero-force macro admission controls; no physical observations."""
import argparse,csv,json,os,shutil,subprocess,time,traceback
from pathlib import Path
CASES=(('TotalParticleAcceleration_default(a,s,p,x,v,n)',True),('TotalParticleAcceleration_default(a,s,p,x,v,n);',True),('PIC::Mover::TotalParticleAcceleration_default(a,s,p,x,v,n)',True),('PIC :: Mover :: TotalParticleAcceleration_default ( a, s, p, x, v, n ) ;',True),('TotalParticleAcceleration_default',False),('PIC::Mover::TotalParticleAcceleration_default',False),('TotalParticleAcceleration_default(a,s,p,x,v,n);other()',False),('PIC::Mover::TotalParticleAcceleration_default(a,s,p,x,v,n);;',False),('Moon::TotalParticleAcceleration(a,s,p,x,v,n)',False),('prefixTotalParticleAcceleration_default(a,s,p,x,v,n)',False),('TotalParticleAcceleration_default(a,s,p,x,v,n)+1',False),('',False))
def run(output):
 output=Path(output).resolve();output.mkdir(parents=True,exist_ok=False);start=time.monotonic();result={'status':'ERROR','scope':'native input macro/recipe admission only; no force/chemistry validation','formal_M3_status':'SKIPPED / NOT QUALIFIED'}
 try:
  compiler=shutil.which('g++')
  if not compiler:result.update(status='SKIPPED',reason='C++17 compiler unavailable');return result
  header=Path(__file__).resolve().parent/'sodium_cohort/zero_force_recipe.h';source=output/'probe.cpp';source.write_text('#include '+json.dumps(str(header))+'\n#include <iostream>\nint main(){\n'+''.join('std::cout<<Moon::SodiumCohortInput::ExactZeroAccelerationMacro('+json.dumps(value)+')<<"\\n";\n'for value,_ in CASES)+'}\n')
  command=[compiler,'-std=c++17',str(source),'-o',str(output/'probe')];compiled=subprocess.run(command,text=True,capture_output=True,timeout=30);(output/'compile.log').write_text(compiled.stdout+compiled.stderr)
  if compiled.returncode:raise RuntimeError('actual thin-owner zero-force guard failed compilation')
  native=subprocess.run([str(output/'probe')],text=True,capture_output=True,timeout=10);(output/'native.log').write_text(native.stdout+native.stderr)
  if native.returncode or len(native.stdout.splitlines())!=len(CASES):raise RuntimeError('actual guard returned incomplete or failed native output')
  checks=[]
  for (expression,expected),actual in zip(CASES,native.stdout.splitlines()):
   if actual not in('0','1'):raise ValueError('malformed native guard result')
   checks.append({'expression':expression,'expected':expected,'actual':actual=='1','status':'PASS'if(actual=='1')==expected else'FAIL'})
  with(output/'comparison.csv').open('w',newline='')as stream:
   writer=csv.DictWriter(stream,fieldnames=list(checks[0]));writer.writeheader();writer.writerows(checks)
  # Actual MPI consensus controls are ADMINISTRATIVE fixtures, not physical
  # profiles. Each rank always participates in every case; ERROR outranks an
  # unavailable prerequisite before any rank may enter physical startup.
  mpicxx=shutil.which('mpicxx');mpiexec=shutil.which('mpiexec')
  if not mpicxx or not mpiexec:result.update(status='SKIPPED',reason='MPI consensus compiler/launcher unavailable');return result
  admission=header.with_name('collective_admission.h');mpi_source=output/'admission.cpp'
  mpi_source.write_text('#include '+json.dumps(str(admission))+r"""
#include <iostream>
int main(int argc,char**argv){MPI_Init(&argc,&argv);int rank=0,size=0,failed=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
for(int k=0;k<7;++k){int local=k==0?0:k==1?1:k==2?2:k==3?(rank==size-1?1:0):k==4?(rank==size-1?2:0):k==5?(rank==size-1?2:rank==0?1:0):3;
int expected=k==0?0:k==1||k==3?1:2;
auto r=Moon::SodiumCohortInput::CollectiveAdmission(MPI_COMM_WORLD,local,"ADMINISTRATIVE FIXTURE ONLY","parser-admission");
bool ok=r.severity==expected && (expected==0?r.diagnostic.empty():!r.diagnostic.empty());int value=ok,all=0;MPI_Allreduce(&value,&all,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);if(!all)++failed;
if(rank==0)std::cout<<k<<','<<all<<'\n';}
MPI_Finalize();return failed?1:0;}
""")
  compile_mpi=[mpicxx,'-std=c++17',str(mpi_source),'-o',str(output/'admission')]
  c=subprocess.run(compile_mpi,text=True,capture_output=True,timeout=30);(output/'mpi_compile.log').write_text(c.stdout+c.stderr)
  if c.returncode:raise RuntimeError('actual MPI admission helper compile failed')
  for ranks in(1,8):
   r=subprocess.run([mpiexec,'-n',str(ranks),str(output/'admission')],text=True,capture_output=True,timeout=20,env=dict(os.environ,OMP_NUM_THREADS='1'))
   (output/('mpi_np'+str(ranks)+'.log')).write_text(r.stdout+r.stderr)
   rows=r.stdout.splitlines()
   if r.returncode not in(0,1)or len(rows)!=7:raise RuntimeError('actual MPI consensus returned missing/signal/incomplete result')
   for k,row in enumerate(rows):
    if row not in(str(k)+',0',str(k)+',1'):raise ValueError('actual MPI consensus output malformed')
    checks.append({'expression':'administrativeMPI-np'+str(ranks)+'-case'+str(k),'expected':True,'actual':row.endswith(',1'),'status':'PASS'if row.endswith(',1')else'FAIL'})
  result.update(status='FAIL'if any(c['status']=='FAIL'for c in checks)else'PASS',counts={s:sum(c['status']==s for c in checks)for s in('PASS','FAIL','ERROR','SKIPPED')},checks=len(checks))
 except Exception as e:(output/'error.log').write_text(traceback.format_exc());result.update(status='ERROR',reason=str(e))
 finally:
  result['elapsed_wall_seconds']=time.monotonic()-start;(output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
 return result
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',required=True);a=p.parse_args();r=run(a.output_dir);print(json.dumps(r,indent=2));raise SystemExit({'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[r['status']])
