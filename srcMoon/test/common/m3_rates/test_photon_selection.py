#!/usr/bin/env python3
"""Frozen U15 statistical control of the actual D04 production selector.

This is local table/sampler verification, never daughter-kinetics or linked
chemistry qualification. Constants and exact-binomial criteria are declared
in PHOTON_SELECTION_CONTRACT.md before any sampled output is inspected.
"""
import argparse
from collections import Counter
import csv
from decimal import Decimal, localcontext
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[3]
sys.path.insert(0,str(REPO/'srcMoon/test'))
from chemical_database_build import shared_sources, shared_dependencies
PARENTS=("He","Na","Ne","Ar","H2O","OH_exp")
COUNTS={"He":1,"Na":1,"Ne":1,"Ar":1,"H2O":7,"OH_exp":4}
DRAWS=100000
MASTER_SEED=731993
CLAIMS=30
FAMILY_ALPHA=0.01
FROZEN={
 "reaction_channels.csv":"2d9fa1c8eb72d3bbba0b91ab4c9a912e6ebce2abc80b251096df2f41a4354061",
 "photo_rates_1au.csv":"3aa54915a08258371dbe975949e6afb2ab74618164fadca8bceebe876c926847",
 "thresholds.csv":"a3335b3b37c3d275a4f72e92ae8e5b1f6ecfcbc4e2d84735cb20694dc93b48f1",
 "cross_sections.csv":"e900eb5514b616dcbbde7e3ca124b5c06d4d22a28a3cd429ed58a14c1933c860"}

class Missing(Exception):
    """A real dependency/source prerequisite is unavailable, not failed physics."""

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def unique_object(pairs):
    value={}
    for key,item in pairs:
        if key in value:raise ValueError("duplicate JSON key: "+key)
        value[key]=item
    return value

def read_json(path):
    return json.loads(Path(path).read_text(),object_pairs_hook=unique_object,
                      parse_constant=lambda token:(_ for _ in ()).throw(ValueError("nonfinite JSON: "+token)))

def checked_file(package,relative,expected,evidence):
    # Metadata cannot use traversal or symlinks to silently retarget a source.
    raw=Path(relative)
    if raw.is_absolute() or ".." in raw.parts:raise ValueError("unsafe metadata path: "+relative)
    path=package/raw
    if not path.is_file():raise Missing("authoritative source absent: "+str(path))
    if path.resolve()!=path.absolute():raise ValueError("metadata source path is symlinked: "+str(path))
    digest=sha(path)
    if digest!=expected:raise ValueError("authoritative source hash changed: "+str(path))
    evidence[relative]={"sha256":digest,"bytes":path.stat().st_size}
    return path

def data_preflight(root):
    """Verify D04 independently of D05/CEX files and preserve native metadata."""
    package=root/"photochemistry"
    required=("README.md","provenance.json","qa_report.json","raw_inventory.json",
              "processed/phidrates/processing_summary.json")
    for name in required:
        if not (package/name).is_file():raise Missing("D04 metadata unavailable: "+str(package/name))
    inventory=root/"INVENTORY.md"
    if not inventory.is_file():raise Missing("native validation-data inventory unavailable")
    records=[line.split('|') for line in inventory.read_text().splitlines() if line.startswith('| D04 |')]
    if len(records)!=1:raise ValueError("missing/duplicate native D04 inventory row")
    if records[0][4].strip()!="Complete/data-ready":raise Missing("D04 native inventory is not Complete/data-ready")
    provenance=read_json(package/"provenance.json");qa=read_json(package/"qa_report.json")
    if qa.get('status') in ('PARTIAL','INCOMPLETE','SKIPPED'):raise Missing("D04 package QA is "+qa['status'])
    if qa.get('status')!='PASS' or qa.get('errors') or qa.get('unresolved_issues'):
        raise ValueError("D04 QA is not an unambiguous source PASS")
    if any(c.get('status')!='PASS' for c in qa.get('checks',[])) or not qa.get('checks'):
        raise ValueError("D04 QA checks are incomplete/nonpassing")
    if provenance.get('schema_version')!=1:raise ValueError("unsupported D04 provenance schema")
    evidence={name:{'sha256':sha(package/name),'bytes':(package/name).stat().st_size} for name in required}
    evidence['../INVENTORY.md']={'sha256':sha(inventory),'native_D04_status':records[0][4].strip()}
    raw=read_json(package/"raw_inventory.json")
    indexed={row['path']:row for row in raw['files']}
    if len(indexed)!=len(raw['files']):raise ValueError("duplicate raw-inventory identity")
    checked_file(package,provenance['raw_inventory']['path'],provenance['raw_inventory']['sha256'],evidence)
    used=[rel for rel in provenance['raw_files_actually_used'] if rel.startswith('raw/phidrates/')]
    if len(used)!=30 or len(set(used))!=30:raise ValueError("D04 source raw membership changed")
    for relative in used:
        record=indexed[relative]
        if record.get('immutable') is not True:raise ValueError("D04 raw source not declared immutable")
        checked_file(package,relative,record['sha256'],evidence)
    for record in provenance['outputs']:
        if record['path'].startswith('processed/phidrates/'):
            checked_file(package,record['path'],record['sha256'],evidence)
    # Photon-specific processing scripts plus the shared pipeline/QA scripts
    # bind the documented original processing; D05 numerical files are unused.
    for record in provenance['processing_scripts']:
        if 'phidrates' in record['path'] or record['path'].endswith(('run_photochemistry_pipeline.py','validate_photochemistry.py')):
            checked_file(package,record['path'],record['sha256'],evidence)
    photo=package/'processed/phidrates'
    for name,digest in FROZEN.items():checked_file(package,'processed/phidrates/'+name,digest,evidence)
    with (photo/'photo_rates_1au.csv').open(newline='') as stream:rates=list(csv.DictReader(stream))
    with (photo/'reaction_channels.csv').open(newline='') as stream:catalogue=list(csv.DictReader(stream))
    reactive={r['reaction_id'] for r in catalogue if r['reaction_type']!='total_photoabsorption'}
    if len(rates)!=30 or len(reactive)!=15:raise ValueError("D04 reactive family membership changed")
    cases=[]
    for parent in PARENTS:
        for solar in (0,1):
            rows=[r for r in rates if r['species']==parent and Decimal(r['solar_activity'])==solar]
            if len(rows)!=COUNTS[parent] or len({r['reaction_id'] for r in rows})!=len(rows):
                raise ValueError("wrong source family: "+parent)
            with localcontext() as context:
                context.prec=80
                values=[Decimal(r['rate_1au_s-1']) for r in rows]
                if any(not v.is_finite() or v<=0 for v in values):raise ValueError("nonpositive reactive hazard")
                total=sum(values,Decimal(0));probabilities=[v/total for v in values]
            channels=[]
            for row,p in zip(rows,probabilities):
                if row['reaction_id'] not in reactive or Decimal(row['optical_depth'])!=0 or Decimal(row['distance_au'])!=1:
                    raise ValueError("source identity/solar normalization changed")
                if row['source_file'] not in used:raise ValueError("rate row is not bound to actual raw membership")
                channels.append({'channel':row['reaction_id'],'rate_1au_s':row['rate_1au_s-1'],
                                 'probability_decimal':str(p),'probability':float(p)})
            cases.append({'case_index':len(cases),'parent':parent,'solar':solar,'channels':channels})
    return photo,evidence,cases

def write_json(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def exact_bounds(distribution,n,p,alpha):
    """Inclusive exact-binomial quantiles with analytic degenerate support.

    Some installed SciPy versions return0 for a low quantile at p=1. The
    distribution itself is exactly X=N in that case, independent of a special
    function implementation. Check quantile minimality as well as rejection
    tails: checking tails alone would incorrectly accept the broad[0,N] band.
    """
    if type(n) is not int or n<0 or not (0<=p<=1) or not (0<alpha<1):
        raise ValueError('invalid binomial control parameters')
    if p==1:
        lower=upper=n
    elif p==0:
        lower=upper=0
    else:
        lower=int(distribution.ppf(alpha/2,n,p))
        upper=int(distribution.ppf(1-alpha/2,n,p))
    q=alpha/2
    if not (0<=lower<=upper<=n):raise ValueError('binomial quantile outside count support')
    for quantile,probability in ((lower,q),(upper,1-q)):
        if distribution.cdf(quantile,n,p)<probability or distribution.cdf(quantile-1,n,p)>=probability:
            raise ValueError('binomial quantile minimality check failed')
    left=float(distribution.cdf(lower-1,n,p));right=float(distribution.sf(upper,n,p))
    if left>q or right>q:raise ValueError('binomial quantile tail self-check failed')
    return lower,upper,left,right

def oracle_controls(distribution):
    """Independent small/degenerate laws and a deliberately broken oracle.

    For five fair Bernoulli draws, the two excluded extreme outcomes each
    have exact probability1/32; at alpha0.1 the interval is[1,4]. These are
    arithmetic fixtures of the statistical scorer, never invented D04 data.
    """
    if exact_bounds(distribution,5,0.5,0.1)!=(1,4,1/32,1/32):
        raise ValueError('independent five-draw binomial identity failed')
    if exact_bounds(distribution,DRAWS,1.0,FAMILY_ALPHA/CLAIMS)!=(DRAWS,DRAWS,0.0,0.0):
        raise ValueError('degenerate p=1 source support failed')
    if exact_bounds(distribution,17,0.0,0.1)!=(0,0,0.0,0.0):
        raise ValueError('degenerate p=0 oracle control failed')
    class BrokenQuantile:
        ppf=staticmethod(lambda q,n,p:0)
        cdf=staticmethod(distribution.cdf)
        sf=staticmethod(distribution.sf)
    try:
        exact_bounds(BrokenQuantile,5,0.5,0.1)
    except ValueError:
        return 4
    raise ValueError('incorrect quantile oracle was not refused')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',required=True,type=Path)
    parser.add_argument('--data-root',default=Path('/data/vtenishe/moon_validation_data'),type=Path)
    parser.add_argument('--cxx',default='/usr/bin/g++')
    args=parser.parse_args();out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=False)
    started=time.monotonic();report={'status':'ERROR','scope':'U15 local production D04 table/selector statistics only',
                                  'formal_I22_status':'SKIPPED','daughter_kinetics_status':'NOT_QUALIFIED'}
    try:
        try:
            import scipy
            from scipy.stats import binom
            import numpy as np
            os.environ['MPLCONFIGDIR']=str(out/'matplotlib-cache')
            os.environ['XDG_CACHE_HOME']=str(out/'cache')
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except ImportError as error:raise Missing('statistical/plot dependency unavailable: '+str(error))
        if not shutil.which(args.cxx):raise Missing('C++ compiler unavailable: '+args.cxx)
        photo,data,cases=data_preflight(args.data_root.resolve())
        write_json(out/'data_hashes.json',data)
        alpha=FAMILY_ALPHA/CLAIMS
        protocol_controls=oracle_controls(binom)
        for case in cases:
            case['seed_words']=[MASTER_SEED,case['case_index']]
            for claim in case['channels']:
                p=claim['probability'];lower,upper,left,right=exact_bounds(binom,DRAWS,p,alpha)
                claim.update(lower_count=lower,upper_count=upper,rejection_left_tail=left,rejection_right_tail=right)
        if sum(len(c['channels']) for c in cases)!=CLAIMS:raise ValueError('predeclared claim count changed')
        sources=[HERE/'photon_selection_probe.cpp',REPO/'srcMoon/species/Photochemistry.cpp']
        sources += shared_sources(REPO,'photochemistry')
        identities={str(path):sha(path) for path in sources+[REPO/'srcMoon/species/Photochemistry.h',Path(__file__),HERE/'PHOTON_SELECTION_CONTRACT.md']+shared_dependencies(REPO,'photochemistry')}
        frozen={'schema':'srcMoon-U15-photon-selector-statistical-contract-v1','frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
                'draws_per_case':DRAWS,'master_seed':MASTER_SEED,'familywise_alpha':FAMILY_ALPHA,'positive_claims':CLAIMS,
                'alpha_per_claim':alpha,'cases':cases,'source_sha256':identities,'data_hashes_sha256':sha(out/'data_hashes.json'),
                'statistical_oracle_controls_passed':protocol_controls,
                'scipy_version':scipy.__version__,'numpy_version':np.__version__,'matplotlib_version':matplotlib.__version__}
        # This file is complete before compilation or sampling starts. Scoring
        # re-reads its exact bytes; it does not recompute favorable intervals.
        write_json(out/'frozen_bounds.json',frozen);frozen_digest=sha(out/'frozen_bounds.json')
        executable=out/'photon_selection_probe'
        command=[args.cxx,'-std=c++17','-O2','-g','-Wall','-Wextra','-Werror','-pedantic','-I',str(REPO/'srcMoon/species')]
        command += [str(p) for p in sources]+['-o',str(executable)]
        with (out/'compile.log').open('w') as stream:
            result=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,timeout=45)
        write_json(out/'build.json',{'command':command,'returncode':result.returncode,'source_sha256':identities,
                    'compiler':subprocess.check_output([args.cxx,'--version'],text=True).splitlines()[0],
                    'executable_sha256':sha(executable) if executable.is_file() else None})
        if result.returncode:raise RuntimeError('native photon selector probe compilation failed')
        command=[str(executable),str(photo)]
        result=subprocess.run(command,text=True,capture_output=True,timeout=10)
        (out/'native.stdout').write_text(result.stdout);(out/'native.stderr').write_text(result.stderr)
        write_json(out/'native_call.json',{'command':command,'returncode':result.returncode})
        if result.returncode:raise RuntimeError('actual production photon selector failed')
        if frozen_digest!=sha(out/'frozen_bounds.json'):raise ValueError('pre-output bounds were changed')
        frozen=read_json(out/'frozen_bounds.json')
        actual=[json.loads(line,object_pairs_hook=unique_object) for line in result.stdout.splitlines()]
        if len(actual)!=12:raise ValueError('native family count differs from frozen12')
        comparison=[];failed=[]
        for native,case in zip(actual,frozen['cases']):
            for key in ('case_index','parent','solar','seed_words'):
                if native.get(key)!=case[key]:raise ValueError('native identity mismatch: '+key)
            if native.get('draws')!=DRAWS or not (0<native['uniform_min']<=native['uniform_max']<1):raise ValueError('invalid native draw/support metadata')
            counts=native['counts'];expected={c['channel'] for c in case['channels']}
            if set(counts)!=expected or any(type(n) is not int or n<0 for n in counts.values()) or sum(counts.values())!=DRAWS:
                raise ValueError('native selected channel identity/count is malformed')
            for claim in case['channels']:
                count=counts[claim['channel']];passed=claim['lower_count']<=count<=claim['upper_count']
                row={'parent':case['parent'],'solar':case['solar'],'channel':claim['channel'],'count':count,
                     'expected_count':DRAWS*claim['probability'],'source_probability':claim['probability'],
                     'lower_count':claim['lower_count'],'upper_count':claim['upper_count'],'status':'PASS' if passed else 'FAIL'}
                comparison.append(row)
                if not passed:failed.append(row)
        with (out/'branch_statistics.csv').open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(comparison[0]));writer.writeheader();writer.writerows(comparison)
        write_json(out/'branch_statistics.json',comparison)
        # Preserve the same frozen acceptance bands in both export formats;
        # plotted points are numerical selector outcomes, never observations.
        fig,axes=plt.subplots(6,2,figsize=(13,17),squeeze=False)
        for case,axis in zip(frozen['cases'],axes.flat):
            subset=[r for r in comparison if r['parent']==case['parent'] and r['solar']==case['solar']]
            x=np.arange(len(subset));centers=np.array([r['expected_count'] for r in subset])
            lower=np.array([r['lower_count'] for r in subset]);upper=np.array([r['upper_count'] for r in subset])
            axis.errorbar(x,centers,yerr=[centers-lower,upper-centers],fmt='o',color='black',capsize=3,label='source + frozen binomial interval')
            axis.scatter(x,[r['count'] for r in subset],marker='x',color='tab:blue',label='production sampler')
            axis.set_xticks(x);axis.set_xticklabels([r['channel'].split(':')[-1] for r in subset],rotation=30,ha='right',fontsize=7)
            axis.set_title(case['parent']+' / '+('quiet' if case['solar']==0 else 'active'))
            axis.set_ylabel('selected events / 100000');axis.grid(axis='y',linestyle=':')
        axes[0,0].legend(fontsize=7,framealpha=1.0)
        fig.suptitle('U15 LOCAL SELECTOR VERIFICATION — NOT OBSERVATIONAL VALIDATION\n30 predeclared claims, familywise alpha=0.01',fontsize=12)
        fig.tight_layout(rect=[0,0,1,0.96])
        for extension in ('png','eps'):fig.savefig(out/('branch_statistics.'+extension),dpi=160)
        plt.close(fig)
        report.update(status='FAIL' if failed else 'PASS',executed=CLAIMS,passed=CLAIMS-len(failed),failed=len(failed),errors=0,skipped=0,
                      families_executed=12,families_passed=12-len({(r['parent'],r['solar']) for r in failed}),
                      statistical_oracle_controls_passed=protocol_controls,
                      frozen_bounds_sha256=frozen_digest,draws_per_family=DRAWS,total_draws=12*DRAWS,failed_claims=failed,
                      plots=['branch_statistics.png','branch_statistics.eps'])
    except Missing as error:
        report.update(status='SKIPPED',executed=0,passed=0,failed=0,errors=0,skipped=1,reason=str(error))
    except (OSError,ValueError,KeyError,TypeError,RuntimeError,subprocess.SubprocessError) as error:
        report.update(status='ERROR',executed=0,passed=0,failed=0,errors=1,skipped=0,reason=str(error))
    report['elapsed_s']=time.monotonic()-started
    write_json(out/'result.json',report);print(json.dumps(report,indent=2))
    return {'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[report['status']]

if __name__=='__main__':raise SystemExit(main())
