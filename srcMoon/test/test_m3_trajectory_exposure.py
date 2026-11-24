#!/usr/bin/env python3
"""Analytical production-integrator controls; no physical/observational promotion."""
import argparse
import csv
from datetime import datetime,timezone
from decimal import Decimal, localcontext
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'srcMoon/test'))
from m3_qualification import score_refinement

FILES=['srcMoon/species/TrajectoryExposure.h','srcMoon/species/TrajectoryExposure.cpp',
 'srcMoon/species/PiecewiseHazard.h','srcMoon/species/PiecewiseHazard.cpp',
 'srcMoon/species/TrajectoryIntervals.h','srcMoon/test/trajectory_exposure/probe.cpp',
 'srcMoon/test/trajectory_exposure/README.md','srcMoon/test/M3_TRAJECTORY_EXPOSURE_CONTRACT.md',
 'srcMoon/test/test_m3_trajectory_exposure.py',
 'srcMoon/test/integration-tests/I21/VARYING_EXPOSURE_CONTRACT.md',
 'srcMoon/test/m3_qualification.py','srcMoon/test/M3_QUALIFICATION_ADDENDUM_20261118.json']

def run(output):
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=False)
    result={'status':'ERROR','scope':'production numerical integration only; NOT physical M3 qualification',
        'counts':{s:0 for s in ('PASS','FAIL','ERROR','SKIPPED')},'results':[]}
    manifest={'created_utc':datetime.now(timezone.utc).isoformat(),'commands':[],
        'source_revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'source_hashes':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FILES},
        'external_data':[],'seed':None,'MPI':'none','OpenMP':'none','comparison_rounding_relative':256*sys.float_info.epsilon}
    comparisons=[];radial=[]
    try:
        compiler=shutil.which('g++')
        if not compiler:
            result.update(status='SKIPPED',reason='C++17 compiler unavailable');result['counts']['SKIPPED']=1
            return result
        executable=output/'probe'
        cmd=[compiler,'-std=c++17','-O2','-Wall','-Wextra','-Werror','-pedantic',
            str(ROOT/FILES[5]),str(ROOT/FILES[1]),str(ROOT/FILES[3]),'-o',str(executable)]
        child=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True,timeout=120)
        (output/'compile.stdout.log').write_text(child.stdout);(output/'compile.stderr.log').write_text(child.stderr)
        manifest['commands'].append({'argv':cmd,'returncode':child.returncode})
        if child.returncode:raise RuntimeError('production compilation failed; see compile.stderr.log')
        manifest['executable_sha256']=hashlib.sha256(executable.read_bytes()).hexdigest()
        def query(mode):
            cmd=[str(executable),mode];r=subprocess.run(cmd,capture_output=True,text=True,timeout=30)
            (output/(mode+'.stdout.log')).write_text(r.stdout);(output/(mode+'.stderr.log')).write_text(r.stderr)
            manifest['commands'].append({'argv':cmd,'returncode':r.returncode})
            if r.returncode:raise RuntimeError('native adapter failed: '+mode)
            return json.loads(r.stdout)
        def check(name,passed,**details):
            result['results'].append({'name':name,'status':'PASS' if passed else 'FAIL',**details})
        # H(t)=3*t/(4+t); solve H(t)=1/2 independently of the native
        # interpolation and exact binary clock. Decimal retains the output
        # error below binary64 projection, as in the original I21 auditor.
        with localcontext() as context:
            context.prec=80
            for steps in (8,16,32):
                values=query(f'radial-{steps}')
                depth=Decimal(values['depth'])
                depth_error=abs(depth-1)
                rounding=Decimal(256)*Decimal(2)**-52
                check(f'radial-{steps}-depth-envelope',
                      depth_error<=Decimal(values['bound'])+rounding,
                      absolute_depth_error=str(depth_error),declared_bound=values['bound'])
                check(f'radial-{steps}-event-and-survivor',
                      values['steps']==steps and Decimal(values['accepted_rate_s'])>0
                      and values['survivor_clock_updated'])
                radial.append({'steps':steps,'actual':values['event_s'],'reference':'0.8','scale':'2',
                               'actual_depth':values['depth'],'reference_depth':'1',
                               'normalized_time_error':str(abs(Decimal(values['event_s'])-Decimal('0.8'))/2),
                               'segments':values['segments']})
        score=score_refinement('I21',radial)
        result['radial_refinement']=score
        check('radial-independent-first-hit-refinement',score['status']=='PASS')
        expected={'constant':2.3,'linear':200*(1/9-1/10)+.3,
            'accelerated':200*(1/(2*10*9.5)+math.atanh(math.sqrt(.5/10))/(2*10*math.sqrt(5)))+.3,
            'constants':.3}
        for mode,radius in [('sphere',1),('thin-sphere',.001)]:
            tangent=100*radius/math.sqrt(2500-radius*radius)
            expected[mode]=(200/3)*(math.atan(.03)-math.atan(tangent/100))+.3
        # For a screen halfway to the Sun the shadow is |y|<=0.2. Integrate
        # 200/(100+y²) dy/2 outside it using the independent arctan primitive.
        expected['terrain']=20*(math.atan(.1)-math.atan(.02))+.3
        for mode,oracle in expected.items():
            actual=query(mode);check(mode+'-evaluated',actual['evaluated'])
            if not actual['evaluated']:continue
            numerical=actual['bound']+actual['shadow_bound'];rounding=256*sys.float_info.epsilon*max(1,abs(oracle))
            delta=abs(actual['depth']-oracle)
            check(mode+'-independent-integral',delta<=numerical+rounding,absolute_error=delta,bound=numerical,rounding=rounding)
            check(mode+'-declared-total-budget',numerical<=1e-10+1e-8*oracle)
            comparisons.append({'case':mode,'actual_depth':actual['depth'],'independent_depth':oracle,
                'absolute_error':delta,'interpolation_bound':actual['bound'],'shadow_bound':actual['shadow_bound'],'segments':actual['segments']})
        far=query('far-sun');check('far-sun-evaluated',far['evaluated'])
        if far['evaluated']:
            envelope=far['bound']+far['shadow_bound']+256*sys.float_info.epsilon
            check('far-diagonal-sun-independent-shadow-range',far['independent_lower']-envelope<=far['depth']<=far['independent_upper']+envelope)
            check('far-sun-silhouette-roots-retained',far['roots']>=2)
        close=query('close-roots')
        if close['evaluated']:
            tangent=100*4.5e-7/math.sqrt(2500-(4.5e-7)**2)
            exact=(200/.006)*(math.atan(.003/100)-math.atan(tangent/100))
            envelope=close['bound']+close['shadow_bound']
            check('close-roots-correct-integral-or-explicit-refusal',abs(close['depth']-exact)<=envelope+256*sys.float_info.epsilon and envelope<=2.5e-4)
        else:check('close-roots-correct-integral-or-explicit-refusal',close['unchanged'])
        projection=query('position-projection')
        check('position-projection-evaluated',projection['evaluated'])
        if projection['evaluated']:
            # Exact line-plane support gives t_edge=2*h/v. Reference-distance
            # normalization differs from unity by <3e-11 throughout this
            # fixture; this analytic bound is independent of callback samples.
            visible=1-2*1e-12/1e-11
            independent_lo=visible*(4e12/(4e12+100+1e-4))
            independent_hi=visible
            envelope=projection['bound']+projection['shadow_bound']
            check('position-projection-full-cell-error-covered',
                  independent_lo-envelope<=projection['depth']<=independent_hi+envelope,
                  independent_lower=independent_lo,independent_upper=independent_hi,
                  actual_depth=projection['depth'],bound=envelope)
            check('position-projection-declared-budget',envelope<=2)
            comparisons.append({'case':'position-projection','actual_depth':projection['depth'],
                'independent_depth':visible,'absolute_error':abs(projection['depth']-visible),
                'interpolation_bound':projection['bound'],'shadow_bound':projection['shadow_bound'],
                'segments':projection['segments']})
        restrictive=query('position-projection-budget')
        check('position-projection-restrictive-budget-refusal',not restrictive['evaluated'] and restrictive['unchanged'])
        for mode in ('negative','duplicate','nan','no-frame','missing-callback','root-budget','resource','callback-error'):
            actual=query(mode);check(mode+'-transactional-refusal',not actual['evaluated'] and actual['unchanged'])
        for row in result['results']:result['counts'][row['status']]+=1
        result['status']='FAIL' if result['counts']['FAIL'] else 'PASS'
        # Plot the same unmodified numerical comparisons retained in CSV.
        # The legend explicitly distinguishes an analytical oracle from data;
        # no observational value or uncertainty is generated by this runner.
        os.environ.setdefault('MPLCONFIGDIR',str(output/'matplotlib-cache'))
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(2,1,figsize=(9,7),sharex=True)
        x=list(range(len(comparisons)));names=[r['case'] for r in comparisons]
        axes[0].plot(x,[r['independent_depth'] for r in comparisons],'o',label='independent analytic integral')
        axes[0].plot(x,[r['actual_depth'] for r in comparisons],'x',label='production integral')
        axes[0].set_ylabel('dimensionless optical depth');axes[0].legend(framealpha=1)
        axes[1].plot(x,[r['absolute_error'] for r in comparisons],'o',label='absolute numerical error')
        axes[1].plot(x,[r['interpolation_bound']+r['shadow_bound'] for r in comparisons],'x',label='declared analytic error envelope')
        axes[1].set_ylabel('absolute depth error');axes[1].set_xticks(x);axes[1].set_xticklabels(names)
        axes[1].legend(framealpha=1);fig.suptitle('Analytical numerical fixtures — no physical validation');fig.tight_layout()
        for extension in ('png','eps'):fig.savefig(output/('comparison.'+extension))
        plt.close(fig)
        fig,axis=plt.subplots(figsize=(7,4))
        axis.loglog([row['steps'] for row in radial],[float(row['normalized_time_error']) for row in radial],
                    'o-',label='independent first-hit error')
        axis.axhline(.002,color='black',linestyle='--',label='approved I21 numerical ceiling')
        axis.set_xlabel('Equal steps in the 2 s horizon')
        axis.set_ylabel('Absolute time error / horizon')
        axis.legend(framealpha=1)
        fig.tight_layout()
        for extension in ('png','eps'):fig.savefig(output/('radial_refinement.'+extension))
        plt.close(fig)
    except ImportError as error:
        result.update(status='SKIPPED',reason='required plotting dependency unavailable: '+str(error));result['counts']['SKIPPED']+=1
    except (OSError,ValueError,RuntimeError,subprocess.TimeoutExpired) as error:
        result.update(status='ERROR',reason=str(error));result['counts']['ERROR']+=1
    finally:
        if comparisons:
            with (output/'comparison.csv').open('w',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=list(comparisons[0]));writer.writeheader();writer.writerows(comparisons)
        if radial:
            with (output/'radial_refinement.csv').open('w',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=list(radial[0]));writer.writeheader();writer.writerows(radial)
        (output/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
        (output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output-dir',required=True,type=Path)
    args=parser.parse_args();r=run(args.output_dir);print(json.dumps(r,indent=2));sys.exit({'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[r['status']])
