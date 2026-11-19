#!/usr/bin/env python3
"""Compile/run the conditional normal physical Na owner, retaining honest77.

It uses a freshly source-matched ordinary producer and actual native .in
owner, never a guessed AMPS CLI. Source-only/preflight PASS is not physical
execution or formal qualification. The actual current Na+ gap returns77 before
SPICE/mesh/particle initialization. More than an hour requires cluster planning.
"""
import argparse,json,os,re,shutil,subprocess,sys,time,traceback
from pathlib import Path
from campaign_budget import CampaignBudget
from build_m2_qualification import production_source_identity,sha256
from integration_common import REPO
from run_m3_moon_force_tracers import native_sources
sys.path.insert(0,str(Path(__file__).resolve().parent/'common'))
from moon_testlib import _probe_compile_command

EXIT={'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}
def run(output,producer_manifest,input_file,mode='source-preflight',mpi_ranks=8,wall_seconds=3600):
 output=Path(output).resolve();output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
 result={'schema':'srcMoon-M3-conditional-Na-owner-v1','status':'ERROR','scope':'conditional source-bound normal owner; formal qualification separate','formal_M3_status':'SKIPPED / NOT QUALIFIED','observational_scoring':False,'commands':[]}
 manifest={'execution':'not-run','mode':mode,'mpi_ranks':mpi_ranks,'openmp_threads':1,'source_revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),'random_seed':'recorded from explicit native input if physical phase admitted','source_root':str(REPO)}
 manifest['dirty_tracked_paths']=subprocess.check_output(['git','status','--porcelain','--untracked-files=no'],cwd=REPO,text=True,env=dict(os.environ,GIT_OPTIONAL_LOCKS='0')).splitlines()
 budget=CampaignBudget(wall_seconds,min(120,wall_seconds*.1))
 try:
  if mode not in('source-preflight','preflight','run')or type(mpi_ranks)is not int or mpi_ranks not in(1,8):raise ValueError('explicit native phase and1/8MPI required')
  input_file=Path(input_file).resolve();producer_manifest=Path(producer_manifest).resolve()
  if not input_file.is_file():raise ValueError('native .in request is unavailable')
  shutil.copy2(input_file,output/'input.snapshot.in');manifest['input']={'path':str(input_file),'sha256':sha256(input_file)}
  if not producer_manifest.is_file():result.update(status='SKIPPED',reason='fresh source-matched native producer unavailable');return result
  producer=json.loads(producer_manifest.read_text());normal=producer.get('schema')=='srcMoon-actual-normal-restore-producer-v1'
  if normal:
   # This separate, explicit normal-build receipt is backed by actual successful
   # configuration/build commands and retained bytes. A bare status flag or
   # archive presence is never sufficient, and no Na+ capability is inferred.
   if producer.get('status')!='PASS':result.update(status='SKIPPED',reason='normal producer build did notPASS');return result
   for role in('build_log','executable'):
    record=producer[role]
    if sha256(Path(record['path']))!=record['sha256']:raise ValueError('normal producer '+role+' bytes changed')
   text=Path(producer['build_log']['path']).read_text()
   if re.findall(r'^EXIT: (-?\d+)$',text,re.M)!=['0','0']or len(re.findall(r'^COMMAND: ',text,re.M))!=2:raise ValueError('normal producer lacks actual successful configure/build command pair')
  else:
   producer_result=producer_manifest.with_name('result.json')
   if json.loads(producer_result.read_text()).get('status')!='PASS':result.update(status='SKIPPED',reason='native producer predecessor did notPASS');return result
  library=Path(producer['library']['path']);build=library.parent
  if Path(producer['source_root']).resolve()!=REPO or sha256(library)!=producer['library']['sha256']:raise ValueError('actual native producer root/library identity mismatch')
  snapshots=producer.get('retained_generated_snapshots')if normal else producer['generated_snapshots']
  if not isinstance(snapshots,dict)or not snapshots:raise ValueError('producer has no retained generated ABI snapshots')
  for path,digest in snapshots.items():
   relative=Path(path).relative_to(producer_manifest.parent/'generated')
   if sha256(Path(path))!=digest or sha256(build/relative)!=digest:raise ValueError('producer generated ABI differs from actual preserved/current headers')
  if normal:
   for path,digest in producer['generated_snapshots'].items():
    if sha256(Path(path))!=digest:raise ValueError('normal producer current generated header changed')
  if native_sources(producer['production_source_hashes'])!=native_sources(production_source_identity()):raise ValueError('production source differs from actually compiled producer')
  manifest.update(producer_manifest={'path':str(producer_manifest),'sha256':sha256(producer_manifest)},library=producer['library'],generated_snapshots=producer['generated_snapshots'],production_source_hashes=producer['production_source_hashes'])
  here=Path(__file__).resolve().parent;source=here/'sodium_cohort/native.cpp'
  manifest['test_sources']={str(p):sha256(p)for p in(source,Path(__file__).resolve(),here/'sodium_cohort/profile_projection.py',here/'sodium_cohort/prepare_source_preflight.py',here/'sodium_cohort/README.md',here/'sodium_cohort/zero_force_recipe.h',here/'sodium_cohort/collective_admission.h')}
  executable=output/'native-sodium-cohort';command=_probe_compile_command(build,executable);command.insert(1,'-DMPI_ON');command[command.index(str(here/'common/production_kernel_probe.cpp'))]=str(source)
  def execute(argv,label,cwd,cores):
   folder=output/(label+'-resources');folder.mkdir()
   value=budget.run(argv,cwd,output/(label+'.log'),output=folder,allocated_cores=cores,environment=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
   result['commands'].append(value)
   if value.get('reason'):raise RuntimeError('resource guard: '+str(value['reason']))
   return value
  compiled=execute(command,'compile',REPO,1)
  if compiled['returncode']!=0:raise RuntimeError('conditional actual native owner failed to compile/link')
  manifest['executable']={'path':str(executable),'sha256':sha256(executable)};manifest['execution']='native-'+mode
  work=output/'native-work';work.mkdir();(work/'PT').mkdir();native_output=output/'native-output'
  receipt=execute(['mpiexec','-n',str(mpi_ranks),str(executable),'--input',str(input_file),'--output-dir',str(native_output),'--mode',mode],'native',work,mpi_ranks)
  target=native_output/'result.json'
  if not target.is_file():raise RuntimeError('actual native completed result missing; retain native log')
  native=json.loads(target.read_text());status=native.get('status')
  if status not in EXIT or native.get('native_returncode')!=EXIT[status]or receipt['returncode']!=EXIT[status]:raise ValueError('native status/summary/exit contradiction')
  if mode=='run'and status=='PASS':
   required=[native_output/'parent_records.csv']+[native_output/(f'event_metrics.rank{r}.hex.csv')for r in range(mpi_ranks)]
   if any(not p.is_file()for p in required):raise ValueError('actual physical execution omitted parent/live-product/bath ledger artifacts')
   result['physical_execution_status']='COMPLETED_NOT_YET_QUALIFIED'
  result.update(status=status,reason=native.get('reason'),native_result={'path':str(target),'sha256':sha256(target)})
  manifest['actual_native_phase_result']=native
 except Exception as error:
  (output/'infrastructure_exception.log').write_text(traceback.format_exc());result.update(status='ERROR',reason=str(error),exception_type=type(error).__name__)
 finally:
  result['elapsed_wall_seconds']=time.monotonic()-started;result['counts']={s:int(result['status']==s)for s in EXIT}
  manifest['artifacts']={str(p.relative_to(output)):sha256(p)for p in output.rglob('*')if p.is_file()and p.name not in('result.json','run_manifest.json')}
  (output/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');(output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
 return result

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',required=True);p.add_argument('--producer-manifest',required=True);p.add_argument('--input',required=True);p.add_argument('--mode',choices=('source-preflight','preflight','run'),default='source-preflight');p.add_argument('--mpi-ranks',type=int,choices=(1,8),default=8);p.add_argument('--wall-seconds',type=float,default=3600);a=p.parse_args();r=run(a.output_dir,a.producer_manifest,a.input,a.mode,a.mpi_ranks,a.wall_seconds);print(json.dumps(r,indent=2));raise SystemExit(EXIT[r['status']])
