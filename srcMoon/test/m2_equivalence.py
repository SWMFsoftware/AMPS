"""Independent fixed-family I19/I20 resolution scoring, never observation fitting.

Paired Student-t intervals are approximate unless their marginal/stream
assumptions are justified. The caller must supply frozen normalization and
support evidence; zero empirical variance alone is not a precision claim.
"""
from __future__ import annotations
import csv
import json
import math
from pathlib import Path
import statistics

from m2_protocol import cases, grid_pairs


def equivalence_interval(values, margin, *, family_size=8640, deterministic=False):
    """Evaluate all 16 paired seeds with the approved simultaneous budget."""
    if len(values)!=16 or not all(math.isfinite(x) for x in values):
        raise ValueError('exactly sixteen finite paired differences are required')
    if margin<=0:raise ValueError('equivalence margin must be positive')
    from scipy.stats import t
    mean=statistics.mean(values);sd=statistics.stdev(values)
    critical=float(t.ppf(1-.05/(2*family_size),15));half=critical*sd/4
    supported=sd>0 or deterministic
    return {'mean':mean,'standard_deviation':sd,'half_width':half,
        'lower':mean-half,'upper':mean+half,'margin':margin,
        'supported':supported,'confidence':'simultaneous approximate 95%, conditional on declared assumptions',
        'pass':supported and abs(mean)+half<=margin and half<=margin/2}


def load_primary(directory,*,require_cohorts=False):
    """Select indexed common physical checkpoints; never interpolate a ledger."""
    from integration_common import parse_receipt
    receipt=parse_receipt(directory/'receipt.txt');steps=int(receipt['steps'])
    rows=list(csv.DictReader((directory/'m2_budget.csv').open()))
    picked={}
    for species in ('Na','Ar'):
        group={int(r['completed_steps']):r for r in rows if r['species']==species}
        for j in range(1,9):
            row=group[steps*j//8];values={
                'gas':float(row['gas_atoms']),'surface':float(row['surface_atoms']),
                'external':float(row['emission_atoms'])-float(row['surface_release_atoms']),
                'adsorption':float(row['adsorption_atoms']),
                'release':float(row['surface_release_atoms']),'escape':float(row['escape_atoms'])}
            picked.update({(species,j,key):value for key,value in values.items()})
    path=directory/'m2_cohorts.csv'
    # A lost cohort artifact must not silently shrink the declared comparison
    # family to the easier global metrics after native runs have completed.
    if require_cohorts and not path.is_file():
        raise FileNotFoundError('required native cohort diagnostics unavailable: '+str(path))
    if not require_cohorts and path.exists():
        raise ValueError('unexpected cohort diagnostics in source-fed case: '+str(path))
    if path.is_file():
        cohorts=list(csv.DictReader(path.open()));end=max(float(r['elapsed_s']) for r in cohorts)
        for r in cohorts:
            if float(r['elapsed_s'])!=end:continue
            for field in ('first_survival','surface_fraction','restricted_first_residence_s'):
                value=float(r[field]);name=field
                if field=='restricted_first_residence_s':value/=end;name='restricted_residence_over_H'
                picked[(r['species'],'site-'+r['cohort'],name)]=value
    return receipt,picked


def score_campaign(root,normalizations,support,*,application_input=None,executable=None,support_base=None,sites=None):
    """Score a complete predeclared matrix; missing support remains SKIPPED.

    The normalization manifest is frozen before execution and supplies A_s per
    physical case/species. It is never inferred from the residual being tested.
    Deterministic-zero declarations require independent analytical evidence,
    not a post-hoc claim that sixteen samples happened to be equal.
    """
    root=Path(root)
    if support.get('status')!='READY' or not support.get('marginal_intervals_justified') or not support.get('rng_stream_policy_reviewed'):
        return {'status':'SKIPPED','reason':'statistical interval/stream qualification unavailable','formal_i19_status':'SKIPPED','formal_i20_status':'SKIPPED'}
    selected=cases('cluster');data={};failures=[];records=[]
    # Authenticate ALL 288 applied initial conditions before opening a single
    # numerical difference or interval. Equal but unqualified uniform starts
    # cannot be promoted by a generic periodic proof or a small residual.
    # Re-hashing each native full-state file proves it has the complete grid
    # parsed from the reviewed pilot; no production cached-state load exists.
    from m2_initial_state_evidence import (Unavailable,load_qualification,verify_case_initial_state,bind_frozen_sites)
    if application_input is None or executable is None or support_base is None:
        return {'status':'SKIPPED','reason':'explicit physical input/executable/initial-state proof owner required',
                'formal_i19_status':'SKIPPED','formal_i20_status':'SKIPPED'}
    initial_states={}
    try:
        qualification=load_qualification(support,application_input,executable,base=support_base)
        qualification=bind_frozen_sites(qualification,sites)
        for case in selected:
            initial_states[case['key']]=verify_case_initial_state(root/case['key'],qualification)
    except Unavailable as error:
        return {'status':'SKIPPED','reason':str(error),'scope':'applied initial-state prerequisite; numerical comparisons not evaluated',
                'formal_i19_status':'SKIPPED','formal_i20_status':'SKIPPED'}
    except (OSError,ValueError,KeyError,TypeError) as error:
        return {'status':'ERROR','reason':str(error),'scope':'applied initial-state authentication; numerical comparisons not evaluated',
                'formal_i19_status':'SKIPPED','formal_i20_status':'SKIPPED'}
    grid,pairs=grid_pairs();deterministic=set(support.get('analytically_verified_deterministic_metrics',[]))
    for case in selected:
        directory=root/case['key'];result=json.loads((directory/'result.json').read_text())
        if result['status']!='PASS':
            failures.append({'case':case['key'],'reason':'selected native case did not pass'})
        receipt,metrics=load_primary(directory,require_cohorts=case['case']=='cohort')
        data[(case['case'],case['grid_index'],case['seed_index'])]=(receipt,metrics)
    for name in ('source-fed','cohort'):
        for left,right in pairs:
            keys=data[(name,right,0)][1].keys()
            for key in keys:
                species,checkpoint,metric=key;is_site=isinstance(checkpoint,str)
                scale=1.0 if is_site else float(normalizations[name][species])
                if not math.isfinite(scale) or scale<=0:raise ValueError('missing positive frozen atom normalization')
                differences=[]
                for seed in range(16):
                    a,ma=data[(name,left,seed)];b,mb=data[(name,right,seed)]
                    if not math.isclose(float(a['duration_s']),float(b['duration_s']),rel_tol=1e-14):
                        raise ValueError('refinement physical durations differ')
                    differences.append((ma[key]-mb[key])/scale)
                label=f'{name}/{species}/{checkpoint}/{metric}'
                interval=equivalence_interval(differences,.01,deterministic=label in deterministic)
                record={'case':name,'species':species,'checkpoint':checkpoint,'metric':metric,
                        'left_grid':left,'right_grid':right,**interval}
                records.append(record)
                if not interval['pass']:failures.append({'comparison':label,'left':left,'right':right,'reason':'equivalence/precision/statistical support'})
        # Numerical loss remains distinct from physical escape even when the
        # signed accounting closes. Check every seed/run, not its ensemble mean.
        for case in [c for c in selected if c['case']==name]:
            directory=root/case['key']
            for r in csv.DictReader((directory/'m2_budget.csv').open()):
                if float(r['numerical_removal_atoms'])/float(normalizations[name][r['species']])>.001:
                    failures.append({'case':case['key'],'reason':'numerical-domain loss exceeds 0.1% budget'})
    # Report resolved growth along the dt axis at fixed weight. No second-order
    # assumption or monotone single-seed weight means is imposed.
    by_key={(r['case'],r['species'],str(r['checkpoint']),r['metric'],r['left_grid'],r['right_grid']):r for r in records}
    for name,species,checkpoint,metric in {k[:4] for k in by_key}:
        for w in (1,2,4):
            ids=[grid.index((dt,w)) for dt in (1,2,4)]
            coarse=by_key.get((name,species,checkpoint,metric,ids[0],ids[1]))
            fine=by_key.get((name,species,checkpoint,metric,ids[1],ids[2]))
            if coarse and fine and max(0,abs(fine['mean'])-fine['half_width'])>abs(coarse['mean'])+coarse['half_width']:
                failures.append({'comparison':f'{name}/{species}/{checkpoint}/{metric}','reason':'resolved timestep sensitivity increases'})
    return {'status':'FAIL' if failures else 'PASS','scope':'fixed full numerical design; not observational agreement',
        'formal_i19_status':'SKIPPED','formal_i20_status':'SKIPPED',
        'formal_reason':'numerical design score alone does not supply formal predecessor/data/operator evidence',
        'comparisons':records,'failures':failures,'family_size_bound':8640,
        'applied_initial_state':{'qualification_sha256':qualification['qualification_sha256'],
            'state_sha256':qualification['state_sha256'],'initialization_mode':qualification['initialization_mode'],
            'physics_identity_sha256':qualification['physics_identity_sha256'],'native_case_count':len(initial_states),
            'cases':initial_states,'scope':'qualified applied model initial condition; not independent observations'}}
