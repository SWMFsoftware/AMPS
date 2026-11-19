#!/usr/bin/env python3
"""Run real force-tracer analytical controls; never generate observations.

This is the normal nightly/core target's thin adapter. It does not configure
the active source tree. Every MPI invocation has its own new working directory
and preserves actual stdout/report/metrics plus executable/source identity.
Fatal MPI_ABORT may return1, so a complete native criterion report is required
before treating exit1 as scientific FAIL rather than runtime ERROR.
"""
from __future__ import annotations
import argparse,csv,datetime,hashlib,json,math,os,re,shlex,shutil,subprocess,sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'srcMoon/test'))
from campaign_budget import CampaignBudget

ROW=re.compile(r'^(PASS|FAIL) order=([12]) (.*?) value=(\S+) expected=(\S+)$')
SUMMARY=re.compile(r'^SUMMARY checks=(\d+) passed=(\d+) failed=(\d+)$')
FATAL='quadratic arc predicts entry but integrated event velocity is leaving'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
class CriterionFailure(RuntimeError):pass

def score_native_report(report,stdout,code,mode):
    """Interpret the actual native protocol; do not infer physics from exit1.

    This pure reporting boundary also permits explicitly labelled temporary
    protocol fixtures to exercise corruption/failure handling without calling
    their invented rows physical tests or observational evidence.
    """
    try:lines=Path(report).read_text().splitlines()
    except OSError:return 'ERROR','native report is unavailable',[]
    metrics=[]
    try:
        for line in lines:
            match=ROW.fullmatch(line)
            if match:metrics.append(dict(status=match[1],order=int(match[2]),check=match[3],value=float(match[4]),reference=float(match[5])))
    except ValueError:return 'ERROR','malformed native metric',[]
    summary=SUMMARY.fullmatch(lines[-1]) if lines else None
    expected={'fixed':124,'adaptive':134}.get(mode)
    complete=expected is not None and summary is not None and int(summary[1])==expected and len(metrics)==expected and int(summary[2])+int(summary[3])==expected and all(math.isfinite(row['value']) and math.isfinite(row['reference']) for row in metrics)
    if not complete or 'AMPS: exit:' in stdout or code not in (0,1):return 'ERROR','missing/truncated report or fatal runtime/build failure',metrics
    failures=sum(row['status']=='FAIL' for row in metrics)
    if failures!=int(summary[3]) or code!=(1 if failures else 0):return 'ERROR','native report and process status contradict',metrics
    return ('FAIL' if failures else 'PASS'),'declared native analytical criteria evaluated',metrics

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--amps',type=Path,required=True,help='already built native force-cut fixture executable')
    parser.add_argument('--output-root',type=Path,default=Path.cwd(),help='parent for a NEW timestamped evidence directory')
    parser.add_argument('--mpi-ranks',type=int,default=8,help='parallel rank count; serial1 is also always tested')
    parser.add_argument('--mpiexec',default=os.environ.get('MPIEXEC','mpiexec'),help='launcher, optionally with arguments')
    parser.add_argument('--wall-seconds',type=float,default=3600,help='one clock including every MPI case and reporting')
    args=parser.parse_args()
    if args.mpi_ranks<1:parser.error('--mpi-ranks must be positive')
    executable=args.amps.resolve();base=args.output_root.resolve();base.mkdir(parents=True,exist_ok=True)
    # The latest pointer is merely an index; meaningful results are always
    # unique and never overwritten. Pre-write a nonempty sentinel so an
    # interrupted launcher cannot leave a stale empty nightly diff as PASS.
    latest=base/'force-cut-latest.diff';latest.write_text('ERROR force-tracer verification incomplete\n')
    stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S.%fZ')
    output=base/('force-cut-verification-'+stamp);output.mkdir()
    (base/'force-cut-latest-result.txt').write_text(str(output/'result.json')+'\n')
    budget=CampaignBudget(wall_s=args.wall_seconds,reserve_s=min(120,args.wall_seconds/4),memory_bytes=200_000_000_000)
    try:
        manifest=dict(schema='force-cut-core-verification-v1',scope='analytical production-path verification; NOT lunar observational validation',executable=str(executable),executable_sha256=sha(executable),mpi_ranks=sorted(set((1,args.mpi_ranks))),omp_threads=1,external_observations_used=False,source_inventory_role='workspace bytes; retain normal configuration/build logs separately to bind these to the executable',sources={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'src/pic/pic_mover_force_cut_surface.cpp',ROOT/'src/pic/pic_mover_force_cut_surface_math.h',Path(__file__).resolve(),ROOT/'test/srcForceCutSurfaceMoverTest/main.cpp']})
    except OSError as error:
        # A missing executable/source is infrastructure ERROR with concrete
        # artifacts even before any native child can be launched.
        record=dict(status='ERROR',reason=str(error),runner=str(Path(__file__).resolve()))
        (output/'result.json').write_text(json.dumps(record,indent=2)+'\n');latest.write_text('ERROR '+str(error)+'\n');print(str(error),file=sys.stderr);return 2
    versions={}
    for tool,argv in [('cxx',shlex.split(os.environ.get('CXX','g++'))),('mpi',shlex.split(args.mpiexec))]:
        if not argv or shutil.which(argv[0]) is None:
            record=dict(status='SKIPPED',reason='required '+tool+' executable is unavailable',runner=str(Path(__file__).resolve()))
            (output/'result.json').write_text(json.dumps(record,indent=2)+'\n');latest.write_text('SKIPPED '+record['reason']+'\n');return 77
        try:
            version=subprocess.run([argv[0],'--version'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=5)
            versions[tool]=dict(executable=shutil.which(argv[0]),returncode=version.returncode,version=version.stdout.splitlines()[:3])
        except (OSError,subprocess.TimeoutExpired) as error:versions[tool]=dict(executable=shutil.which(argv[0]),version_unavailable_reason=str(error))
    manifest['tool_versions']=versions
    revision=subprocess.run(['git','-C',str(ROOT),'rev-parse','HEAD'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    manifest['source_revision']=revision.stdout.strip() if revision.returncode==0 else None
    (output/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    results=[]
    try:
        protocol=output/'reporting-protocol';protocol.mkdir()
        control=budget.run([sys.executable,'-B',str(ROOT/'test/srcForceCutSurfaceMoverTest/test_runner_protocol.py')],ROOT,protocol/'stdout.log',output=protocol,environment=os.environ.copy())
        if control['returncode']!=0:raise RuntimeError('reporting-only protocol checks failed: '+str(protocol/'stdout.log'))
        results.append(dict(test='reporting-protocol',status='PASS',reason='10 infrastructure controls; no physical evidence inferred',runner=str(Path(__file__).resolve()),log=str(protocol/'stdout.log'),evaluated_native_checks=0))
        case=output/'numerics';case.mkdir()
        numerical=case/'numerics_test'
        compiler=shlex.split(os.environ.get('CXX','g++'))
        compile_run=budget.run(compiler+['-std=c++17','-O2','-I',str(ROOT/'src/pic'),str(ROOT/'test/srcForceCutSurfaceMoverTest/numerics_test.cpp'),'-o',str(numerical)],ROOT,case/'build.log',output=case/'build-resources',environment=os.environ.copy())
        if compile_run['returncode']!=0:raise RuntimeError('production numerical primitive adapter failed to build: '+str(case/'build.log'))
        execution=budget.run([str(numerical)],case,case/'stdout.log',output=case/'run-resources',environment=os.environ.copy())
        text=(case/'stdout.log').read_text()
        if execution['returncode']==1 and 'FAIL production numerical primitive criterion:' in text:
            raise CriterionFailure('declared production numerical primitive criterion failed: '+str(case/'stdout.log'))
        if execution['returncode']!=0 or 'PASS production numerical primitives checks=27' not in text:
            raise RuntimeError('production numerical primitive criteria failed: '+str(case/'stdout.log'))
        results.append(dict(test='numerics',status='PASS',reason='27 analytical production-primitive checks',runner=str(Path(__file__).resolve()),log=str(case/'stdout.log'),evaluated_native_checks=27))
        rotation=[]
        for match in re.finditer(r'PASS order=([12]) refinement_ratio=(\S+) coarse_error=(\S+) fine_error=(\S+)',text):
            for steps,error in ((64,float(match[3])),(128,float(match[4]))):rotation.append(dict(order=int(match[1]),steps=steps,step_s=1/steps,velocity_error_m_s=error))
        if len(rotation)!=4 or any(not math.isfinite(row['velocity_error_m_s']) or row['velocity_error_m_s']<=0 for row in rotation):raise RuntimeError('missing/corrupt independent rotation error metrics')
        with (case/'rotation_comparison.csv').open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rotation[0]));writer.writeheader();writer.writerows(rotation)
        try:
            import matplotlib;matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except ImportError:
            visualization=dict(status='SKIPPED',reason='optional matplotlib is unavailable; native criteria remain separate')
        else:
            figure,axis=plt.subplots(figsize=(6,4))
            for order in (1,2):
                subset=[r for r in rotation if r['order']==order]
                axis.loglog([r['step_s'] for r in subset],[r['velocity_error_m_s'] for r in subset],'o-',label='production order '+str(order))
            axis.set(xlabel='Fixed step [s]',ylabel='Velocity error against exact rotation [m/s]',title='Analytical unit-frequency rotation; no observations')
            axis.grid(True,which='both');axis.legend();figure.tight_layout()
            for extension in ('png','eps'):figure.savefig(case/('rotation_comparison.'+extension),dpi=180)
            plt.close(figure);visualization=dict(status='PASS',scope='diagnostic visualization of independently evaluated numerical error',files=['rotation_comparison.csv','rotation_comparison.png','rotation_comparison.eps'])
        (case/'visualization.json').write_text(json.dumps(visualization,indent=2)+'\n')
        for mode in ('fixed','adaptive','inconsistent-impact'):
            ranks_list=(1,) if mode=='inconsistent-impact' else sorted(set((1,args.mpi_ranks)))
            for ranks in ranks_list:
                name=mode+'-np'+str(ranks);case=output/name;(case/'PT/plots').mkdir(parents=True)
                argv=shlex.split(args.mpiexec)+['-n',str(ranks),str(executable)]+([] if mode=='fixed' else [mode])
                run=budget.run(argv,case,case/'stdout.log',output=case,allocated_cores=ranks,environment=dict(os.environ,OMP_NUM_THREADS='1'))
                text=(case/'stdout.log').read_text(errors='replace');code=run['returncode'];metrics=[]
                if mode=='inconsistent-impact':
                    # Expected refusal is a local safety-contract PASS, never
                    # a completed physical trajectory or chemistry ensemble.
                    status='PASS' if code==1 and FATAL in text and 'Segmentation fault' not in text else 'ERROR'
                    reason='expected fixed-step inconsistent-impact refusal' if status=='PASS' else 'missing expected refusal or unrelated runtime failure'
                else:
                    report=case/('force-cut-surface.np'+str(ranks)+'.report')
                    status,reason,metrics=score_native_report(report,text,code,mode)
                    if metrics:
                        with (case/'metrics.csv').open('w',newline='') as stream:
                            writer=csv.DictWriter(stream,fieldnames=['status','order','check','value','reference']);writer.writeheader();writer.writerows(metrics)
                item=dict(test=name,status=status,reason=reason,native_returncode=code,runner=str(Path(__file__).resolve()),log=str(case/'stdout.log'),evaluated_native_checks=len(metrics),artifact_directory=str(case))
                results.append(item);(case/'result.json').write_text(json.dumps(item,indent=2)+'\n')
                print(status+' '+name+' log='+str(case/'stdout.log'),flush=True)
                if status=='ERROR':break
            if results[-1]['status']=='ERROR':break
    except CriterionFailure as error:
        results.append(dict(test='numerics',status='FAIL',reason=str(error),runner=str(Path(__file__).resolve()),log=str(output/'numerics/stdout.log')))
        results.append(dict(test='native-tracers',status='SKIPPED',reason='local numerical predecessor failed; linked execution cannot substitute for it',runner=str(Path(__file__).resolve()),log=str(output/'numerics/stdout.log')))
    except (OSError,ValueError,RuntimeError) as error:
        results.append(dict(test='runner',status='ERROR',reason=str(error),runner=str(Path(__file__).resolve()),log=str(output)))
    status='ERROR' if any(x['status']=='ERROR' for x in results) else 'FAIL' if any(x['status']=='FAIL' for x in results) else 'PASS'
    summary=dict(schema='force-cut-core-verification-v1',status=status,scope=manifest['scope'],results=results,counts={s:sum(x['status']==s for x in results) for s in ['PASS','FAIL','ERROR','SKIPPED']})
    (output/'result.json').write_text(json.dumps(summary,indent=2)+'\n')
    latest.write_text('' if status=='PASS' else '\n'.join(x['status']+' '+x['test']+' '+x['reason']+' log='+x['log'] for x in results if x['status']!='PASS')+'\n')
    return {'PASS':0,'FAIL':1,'ERROR':2}[status]
if __name__=='__main__':raise SystemExit(main())
