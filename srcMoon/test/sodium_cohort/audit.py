#!/usr/bin/env python3
"""Audit actual complete cohort launches; no observations or provider stubs.

The enclosing resource supervisor owns the original campaign deadline. This
child performs independent source-ledger, original-parent statistics and
first-hit convergence audits together, retaining CSV/PNG/EPS and all failures.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path
import sys
import traceback
import subprocess

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import m3_qualification as q
import verify_ledger as ledger
import verify_numerical as numerical

EXIT={'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}


def input_documents(rows):
    """Validate the complete receipt graph BEFORE mapping paths to identities."""
    if not isinstance(rows,list) or not rows:raise ValueError('complete nonempty native input graph required')
    result={}
    for row in rows:
        if not isinstance(row,dict) or set(row)!={'path','sha256'} or not isinstance(row['path'],str) or not row['path'] or not q.digest_syntax(row['sha256']):
            raise ValueError('native input graph requires explicit path/SHA256 records')
        path=str(Path(row['path']).resolve())
        if path in result:raise ValueError('duplicate native input document: '+path)
        result[path]=row['sha256']
    return result


def authenticate_launch(launch):
    """Bind consumed native bytes to their ORIGINAL completed owner receipt.

    Rehashing only at scoring time cannot detect an altered-but-plausible CSV.
    The executor pins these two owner documents immediately on return; their
    artifact map authenticates every consumed input/parent/exact-event record.
    This checks actual provenance, never creates a physical/readiness approval.
    """
    def descriptor(value):
        if not isinstance(value,dict) or set(value)!={'path','sha256'} or not isinstance(value['path'],str) or not q.digest_syntax(value['sha256']):
            raise ValueError('actual owner requires explicit path/SHA256 descriptor')
        path=Path(value['path']).resolve()
        if q.sha(path)!=value['sha256']:raise ValueError('actual completed owner identity changed: '+str(path))
        return path
    manifest_path=descriptor(launch['owner_manifest']);result_path=descriptor(launch['owner_result'])
    manifest=q.strict_object(manifest_path);result=q.strict_object(result_path)
    if result.get('status')!='PASS' or manifest.get('execution')!='native-run' or manifest.get('mode')!='run' or manifest.get('mpi_ranks')!=launch['mpi_ranks'] or manifest.get('openmp_threads')!=1:
        raise ValueError('owner did not complete this actual native physical-phase launch')
    requested=launch['input'];actual=manifest.get('input')
    if not isinstance(actual,dict) or Path(actual['path']).resolve()!=Path(requested['path']).resolve() or actual['sha256']!=requested['sha256']:
        raise ValueError('actual owner did not use the frozen native root input')
    native=Path(launch['output']).resolve()
    paths=([native/'input_receipt.json',native/'parent_records.csv']+
        [native/('event_metrics.rank%d.hex.csv'%rank) for rank in range(launch['mpi_ranks'])])
    artifacts=manifest.get('artifacts')
    if not isinstance(artifacts,dict):raise ValueError('actual owner artifact map unavailable')
    for path in paths:
        relative=str(path.relative_to(manifest_path.parent))
        expected=artifacts.get(relative)
        if not q.digest_syntax(expected) or q.sha(path)!=expected:raise ValueError('completed native artifact changed or was not bound: '+str(path))
    return native


def compare_parent_layouts(first, second):
    """Compare original-parent endpoints independently of MPI order.

    Exact keyed fates and channels are invariant at an identical timestep.
    The original first-hit oracle scores timing against the frozen numerical
    criterion. Report the largest serialized inter-layout timing difference
    in seconds in addition to that criterion.
    """
    by_id = lambda rows: {row['original_parent_id']: row for row in rows}
    a, b = by_id(first), by_id(second)
    missing = sorted(set(a) ^ set(b), key=int)
    mismatches = []
    max_difference = 0.0
    for parent in sorted(set(a) & set(b), key=int):
        lhs, rhs = a[parent], b[parent]
        if (lhs['status'], lhs['channel_id']) != (rhs['status'], rhs['channel_id']):
            mismatches.append(parent)
        elif lhs['status'] == 'REACTED':
            max_difference = max(max_difference, abs(float(lhs['relative_event_time_s']) -
                                                       float(rhs['relative_event_time_s'])))
    return {'status': 'FAIL' if missing or mismatches else 'PASS',
            'matched_original_parents': len(set(a) & set(b)),
            'unpaired_parent_ids': missing, 'changed_fate_or_channel_parent_ids': mismatches,
            'maximum_serialized_time_difference_s': max_difference}


def run(plan_path,runs_path,output,*,development=False):
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=False)
    records=[];manifest={'scope':'actual cohort audits, NOT whole I21-I24/M3 qualification',
        'development':development,
        'source_revision':None,'input_sha256':{str(Path(plan_path).resolve()):q.sha(plan_path),
            str(Path(runs_path).resolve()):q.sha(runs_path)}}
    try:
        plan=q.strict_object(plan_path);q.validate_plan(plan)
        launches=json.loads(Path(runs_path).read_text())
        if not isinstance(launches,list) or not launches:raise ValueError('actual native launch registry required')
        # The independent audit must score the entire frozen matrix. An
        # otherwise authentic run from a different level or input cannot
        # substitute for a failed trial that was part of this plan.
        frozen = plan['native_adapter']['runs']
        def launch_identity(record):
            return (record['mpi_ranks'],record['physical_steps'],
                    str(Path(record['input']['path']).resolve()),record['input']['sha256'])
        if list(map(launch_identity,launches))!=list(map(launch_identity,frozen)):
            raise ValueError('actual launch registry differs from the complete frozen case matrix')
        varying=plan.get('test_kind')=='trajectory-integrated-hazard'
        if varying and not development:
            raise ValueError('straight-line trajectory scoring is a numerical development scope')
        verified_runs=[]
        unmatched_layouts={}
        profile=Path(plan['source_profile']['path'])
        if not profile.is_absolute():profile=Path(plan_path).resolve().parent/profile
        if q.sha(profile)!=plan['source_profile']['sha256']:raise ValueError('frozen physical profile changed')
        for index,launch in enumerate(launches):
            native=authenticate_launch(launch);receipt=q.strict_object(native/'input_receipt.json')
            expected_docs=input_documents(launch['input_documents'])
            actual_docs=input_documents(receipt['input_documents'])
            if expected_docs!=actual_docs:raise ValueError('actual native include graph differs from frozen input documents')
            if receipt['physical_steps']!=launch['physical_steps']:raise ValueError('actual common step count differs from frozen refinement')
            for path,digest in expected_docs.items():
                if q.sha(path)!=digest:raise ValueError('input changed while native cohort was executing')
            parent_path=native/'parent_records.csv';parents=q.read_original_parents(parent_path,plan)
            level=launch['physical_steps']
            if level in unmatched_layouts:
                other_ranks,other_parents=unmatched_layouts.pop(level)
                if other_ranks==launch['mpi_ranks']:
                    raise ValueError('duplicate MPI layout at matched physical step count')
                rank_result=compare_parent_layouts(other_parents,parents)
                (output/f'parent_mpi_steps_{level}.json').write_text(json.dumps(rank_result,indent=2)+'\n')
                records.append({'name':f'original-parent-MPI-match-{level}', 'status':rank_result['status']})
            else:
                unmatched_layouts[level]=(launch['mpi_ranks'],parents)
            reacted={row['original_parent_id'] for row in parents if row['status']=='REACTED'}
            paths=[native/('event_metrics.rank%d.hex.csv'%rank) for rank in range(launch['mpi_ranks'])]
            actual_events=ledger.read_events(paths)
            if set(actual_events)!=reacted:raise ValueError('actual per-rank event ledgers omit/add original parent removals')
            verified_runs.append((native,launch['mpi_ranks'],actual_events))
            time_links=numerical.event_time_links(actual_events,parents,receipt['transport_origin_et_s'])
            (output/('event_time_links-%02d.json'%index)).write_text(json.dumps(time_links,indent=2)+'\n')
            records.append({'name':'actual-parent-native-time-links-%02d'%index,'status':time_links['status']})
            mapping=receipt['mapping']
            if mapping.get('profile_sha256')!=plan['source_profile']['sha256'] or receipt['profile']['sha256']!=mapping['profile_sha256']:
                raise ValueError('native state mapping is not bound to frozen profile')
            if {k:mapping.get(k) for k in ('parent','ion','electron')}!=plan['native_adapter']['state_mapping']:
                raise ValueError('native state selections differ from frozen mapping')
            reference=numerical.reference_receipt(receipt)
            if plan.get('test_kind')=='zero-hazard-invariant':
                # The visibility flag recorded by production is itself derived
                # from the photon rate. Require a separate finite-segment
                # geometric proof before crediting a fixed lunar dark case.
                shadow=numerical.fixed_moon_shadow(plan,launch,receipt)
                (output/('lunar_shadow-%02d.json'%index)).write_text(json.dumps(shadow,indent=2)+'\n')
                records.append({'name':'independent-fixed-lunar-shadow-%02d'%index,'status':shadow['status']})
            reference_rate,channel=numerical.source_rate(reference,manifest['input_sha256'])
            motion=numerical.StraightLineExposure(plan,receipt,reference_rate) if varying else None
            inverse_square=numerical.inverse_square_checks(actual_events,motion.point_rate if motion else reference_rate)
            if not actual_events and reference_rate > 0:
                # An exposed reference still needs actual accepted-event rate
                # evidence. An empty set only proves the zero-rate identity.
                inverse_square.update(status='SKIPPED',reason='no accepted events for positive source rate')
            (output/('inverse_square-%02d.json'%index)).write_text(json.dumps(inverse_square,indent=2)+'\n')
            records.append({'name':'fixed-I21-inverse-square-%02d'%index,'status':inverse_square['status']})
            base=numerical.D04/'processed/phidrates'
            mapping={**mapping,'photon_source':{'reaction_id':channel,
                'solar_activity':'0' if reference['solar_endpoint']=='quiet' else '1',
                **{key:{'path':str(base/name),'sha256':q.sha(base/name)} for key,name in
                    [('rates','photo_rates_1au.csv'),('thresholds','thresholds.csv')]}}}
            if 'rate_s' not in receipt['bounds']:raise ValueError('actual native rate_s numerical bound missing')
            bounds={**receipt['bounds'],'schema':'srcMoon-M3-event-ledger-bounds-v1',
                'scope':'source-bound-experimental','event_frame':mapping['event_frame']}
            reconstructed=ledger.verify(paths,profile,mapping,bounds)
            ledger.write_artifacts(reconstructed,output/('ledger-%02d'%index))
            records.append({'name':'actual-ledger-%02d'%index,'status':reconstructed['status']})
            if development:
                # Small, explicitly selected development cohorts diagnose the
                # clock and ledger. The observed endpoint fraction is retained
                # without applying a statistical gate to an underpowered run.
                survivors=len(parents)-len(reacted)
                score={'status':'NOT_EVALUATED','original_parents':len(parents),
                       'survivors':survivors,'first_removals':len(reacted),
                       'observed_survival_fraction':survivors/len(parents)}
                if motion:
                    score['reference_survival_probability_interval']=[
                        math.exp(-float(motion.end_depth_upper)),math.exp(-float(motion.end_depth_lower))]
                    score['reference_end_optical_depth_interval']=[
                        str(motion.end_depth_lower),str(motion.end_depth_upper)]
                else:
                    score['reference_survival_probability']=math.exp(-float(reference_rate)*float(plan['horizon_s']))
            else:
                score=q.score_cohort(parent_path,plan)
                records.append({'name':'independent-statistics-%02d'%index,'status':score['status']})
            (output/('statistics-%02d.json'%index)).write_text(json.dumps(score,indent=2)+'\n')
        if unmatched_layouts:
            raise ValueError('unpaired MPI layout in the frozen step matrix')
        convergence=numerical.verify(verified_runs,plan)
        (output/'numerical.json').write_text(json.dumps(convergence,indent=2)+'\n')
        comparisons=convergence['comparisons']
        if comparisons:
            with (output/'numerical_comparison.csv').open('w',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=list(comparisons[0]));writer.writeheader();writer.writerows(comparisons)
        records.append({'name':'independent-numerical-refinement','status':convergence['status']})
        import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
        figure,axis=plt.subplots(figsize=(7,4))
        for ranks,score in convergence['scores'].items():
            axis.plot(range(len(score['relative_errors'])),score['relative_errors'],marker='o',label=ranks+'MPI')
        axis.axhline(q.NUMERICAL_LIMIT['I21'],label='approved numerical bias budget',color='black',linestyle='--')
        axis.set_xlabel('Frozen successive time refinement');axis.set_ylabel('Maximum normalized first-hit time error');axis.legend(framealpha=1)
        figure.tight_layout()
        for extension in ('png','eps'):figure.savefig(output/('numerical_comparison.'+extension))
        plt.close(figure)
    except (q.MissingPrerequisite,ledger.Unavailable) as error:
        records.append({'name':'required-evidence','status':'SKIPPED','reason':str(error)})
    except ArithmeticError as error:
        records.append({'name':'independent-first-hit','status':'FAIL','reason':str(error)})
    except Exception as error:
        (output/'exception.log').write_text(traceback.format_exc());records.append({'name':'infrastructure','status':'ERROR','reason':str(error)})
    counts={s:sum(r['status']==s for r in records) for s in EXIT}
    status=next(s for s in ('ERROR','FAIL','SKIPPED','PASS') if counts[s])
    result={'schema':'srcMoon-M3-actual-cohort-audit-v1','status':status,'checks':records,'counts':counts,
        'reason':'actual cohort source/budget/statistical/numerical checks evaluated where available',
        'statistical_acceptance':'NOT_EVALUATED' if development else 'SEE_PER_LAUNCH_RESULTS',
        'formal_M3_status':'NOT QUALIFIED BY ONE COHORT; complete frozen protocol matrix remains',
        'observational_scoring':False}
    manifest['source_revision']=subprocess.check_output(['git','rev-parse','HEAD'],cwd=q.TEST.parents[1],text=True).strip()
    manifest['verifier_sha256']={str(path):q.sha(path) for path in (Path(__file__),Path(ledger.__file__),Path(numerical.__file__),Path(q.__file__))}
    (output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    (output/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',required=True,type=Path);parser.add_argument('--runs',required=True,type=Path)
    parser.add_argument('--output-dir',required=True,type=Path)
    parser.add_argument('--development',action='store_true',help='evaluate clock and ledger with diagnostic survival only')
    args=parser.parse_args()
    result=run(args.plan,args.runs,args.output_dir,development=args.development)
    print(json.dumps(result,indent=2));raise SystemExit(EXIT[result['status']])
