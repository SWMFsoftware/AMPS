#!/usr/bin/env python3
"""Actual normal-SPICE coefficient stage checks; no chemistry activation."""
import argparse,collections,csv,json,os,re,shutil,sys,time,traceback
from pathlib import Path
from campaign_budget import CampaignBudget
from build_m2_qualification import production_source_identity,sha256
from integration_common import REPO
from run_m3_moon_force_tracers import native_sources
from test_m3_runtime_rates import change
sys.path.insert(0,str(Path(__file__).resolve().parent/'common'))
from moon_testlib import _probe_compile_command
EXIT={'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}
NAMES=('same_epoch_Sun_Earth_positions_SI','same_epoch_Sun_Earth_velocities_SI','begin_heliocentric_state_SI','free_velocity_midpoint_rotation','independent_Rdot_Rtranspose_angular_velocity','actual_body_rotation_and_declared_derivative_policy','begin_middle_end_epoch_policy','stage_does_not_inject_or_sample')
def run(output,build_evidence,wall_seconds=600):
 output=Path(output).resolve();output.mkdir(parents=True,exist_ok=False);start=time.monotonic();counts=collections.Counter();metrics=[]
 result={'status':'ERROR','scope':'actual normal frozen/evolving environment stage wiring; no observations or physical chemistry qualification','formal_M3_status':'SKIPPED / NOT QUALIFIED','commands':[]}
 manifest={'mpi_ranks':[1,8],'frozen_orders':[1,2],'evolving_mover':'legacy','OMP_NUM_THREADS':1,'external_observations':[],'source_root':str(REPO)};budget=CampaignBudget(wall_seconds,min(30,wall_seconds*.1))
 try:
  if build_evidence is None:result.update(status='SKIPPED',reason='authentic normal producer unavailable');return result
  proof=json.loads(Path(build_evidence).read_text())if isinstance(build_evidence,(str,Path))else build_evidence
  if proof.get('status')!='PASS'or Path(proof['source_root']).resolve()!=REPO:raise ValueError('actual normal producer status/root mismatch')
  for key in('library','executable','build_log'):
   if sha256(Path(proof[key]['path']))!=proof[key]['sha256']:raise ValueError('normal producer '+key+' bytes changed')
  if re.findall(r'^EXIT: (-?\d+)$',Path(proof['build_log']['path']).read_text(),re.M)!=['0','0']:raise ValueError('producer was not actually configured/built successfully')
  for path,digest in proof['generated_snapshots'].items():
   if sha256(Path(path))!=digest:raise ValueError('normal generated ABI differs from actual producer')
  if native_sources(proof['production_source_hashes'])!=native_sources(production_source_identity()):raise ValueError('stage production sources differ from actual compiled producer')
  manifest['producer']=proof;here=Path(__file__).resolve().parent;source=here/'m3_native_coefficient_stage_probe.cpp';manifest['test_sources']={str(p):sha256(p)for p in(source,Path(__file__).resolve(),here/'coefficient_stage/CONTRACT.md')}
  library=Path(proof['library']['path']);exe=output/'native-stage';command=_probe_compile_command(library.parent,exe);command.insert(1,'-DMPI_ON');command[command.index(str(here/'common/production_kernel_probe.cpp'))]=str(source)
  def execute(argv,label,cwd,ranks):
   resource=output/(label+'-resources');resource.mkdir(parents=True)
   value=budget.run(argv,cwd,output/(label+'.log'),output=resource,allocated_cores=ranks,environment=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
   result['commands'].append(value)
   if value.get('reason')or value['returncode']not in(0,1):raise RuntimeError(label+' infrastructure failure: '+str(value.get('reason')or value['returncode']))
   return value
  if execute(command,'compile',REPO,1)['returncode']:raise RuntimeError('actual coefficient-stage probe failed to compile/link')
  manifest['executable']={'path':str(exe),'sha256':sha256(exe)};kernels={};layouts=[]
  for mode,order,ranks in[(m,o,r)for m,orders in [('frozen',('first','second')),('evolving',('legacy',))]for o in orders for r in(1,8)]:
   label=mode+'-'+order+'-np'+str(ranks);folder=output/label;folder.mkdir();work=folder/'runtime';(work/'PT/plots').mkdir(parents=True)
   template=REPO/('srcMoon/examples/m1_regression.in'if order=='legacy'else'srcMoon/examples/m3_common_time_'+order+'_order.in');text=template.read_text();text=change(text,'simulation_environment_mode',mode)
   for key,extension in(('surface_cea_file','cea'),('surface_tecplot_file','dat')):text=change(text,key,str(folder/('sphere.'+extension)))
   application=folder/'application.in';application.write_text(text);manifest.setdefault('inputs',[]).append({'path':str(application),'sha256':sha256(application)})
   execution=execute(['mpiexec','-n',str(ranks),str(exe),str(application),str(folder/'stage'),mode,order],label,work,ranks)
   with(folder/'stage.checks.csv').open(newline='')as stream:
    reader=csv.DictReader(stream)
    if reader.fieldnames!=['check','status']:raise ValueError('native stage named-check schema changed')
    rows=list(reader)
   expected={'step'+str(step)+'_'+name for step in range(3)for name in NAMES}|{'three_stage_external_epoch_evolution'}
   if mode=='frozen':expected.add('actual_normal_Na_force_uses_initialized_coefficients_on_allocated_owner')
   summary=re.findall(r'M3 coefficient stage mode=(frozen|evolving) checks=(\d+) status=(PASS|FAIL)',(output/(label+'.log')).read_text())
   if len(rows)!=len(expected)or{r['check']for r in rows}!=expected or any(r['status']not in('PASS','FAIL')for r in rows)or len(summary)!=1 or summary[0][0]!=mode or int(summary[0][1])!=len(expected):raise ValueError('stage completed membership/summary incomplete or corrupt')
   failed=any(r['status']=='FAIL'for r in rows)
   if execution['returncode']!=int(failed)or(summary[0][2]=='FAIL')!=failed:raise ValueError('native stage exit/summary/criteria contradiction')
   counts.update(r['status']for r in rows);layouts.append({'mode':mode,'order':order,'mpi_ranks':ranks,'checks':len(rows)})
   with(folder/'stage.metrics.csv').open(newline='')as stream:
    for row in csv.DictReader(stream):metrics.append(dict(row,layout=label))
   files=(folder/'stage.spice_files.txt').read_text().splitlines()
   if not files:raise ValueError('actual furnished kernel inventory empty')
   for name in files:kernels[name]=sha256(Path(name))
  with(output/'comparison.csv').open('w',newline='')as stream:
   writer=csv.DictWriter(stream,fieldnames=['step','quantity','index','actual','expected','layout']);writer.writeheader();writer.writerows(metrics)
  import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
  figure,axes=plt.subplots(1,2,figsize=(12,4));selected=[r for r in metrics if r['quantity']=='angular_velocity_rad_s'and r['index']=='2']
  axes[0].plot([float(r['expected'])for r in selected],label='independent CSPICE derivative',marker='x');axes[0].plot([float(r['actual'])for r in selected],label='normal production stage',linestyle='none',marker='o',fillstyle='none');axes[0].set_ylabel('LSO angular velocity z [rad/s]');axes[0].legend(framealpha=1)
  axes[1].plot([float(r['actual'])-float(r['expected'])for r in selected],'.');axes[1].set_ylabel('Production minus independent [rad/s]')
  for axis in axes:axis.set_xlabel('Declared mode/order/rank/step record');axis.grid(True)
  figure.suptitle('Normal environment stage wiring — no observations/physical Na qualification');figure.tight_layout()
  for extension in('png','eps'):figure.savefig(output/('comparison.'+extension))
  plt.close(figure);manifest.update(actual_furnished_kernels=kernels,layouts=layouts)
  result.update(status='FAIL'if counts['FAIL']else'PASS',counts={s:counts[s]for s in EXIT},checks=sum(counts.values()))
 except Exception as e:
  (output/'infrastructure_exception.log').write_text(traceback.format_exc());result.update(status='ERROR',reason=str(e),exception_type=type(e).__name__)
 finally:
  result['elapsed_wall_seconds']=time.monotonic()-start;manifest['artifact_sha256']={str(p.relative_to(output)):sha256(p)for p in output.rglob('*')if p.is_file()and p.name not in('result.json','run_manifest.json')};(output/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');(output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
 return result
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',required=True);p.add_argument('--normal-build-evidence',required=True);p.add_argument('--wall-seconds',type=float,default=600);a=p.parse_args();r=run(a.output_dir,a.normal_build_evidence,a.wall_seconds);print(json.dumps(r,indent=2));raise SystemExit(EXIT[r['status']])
