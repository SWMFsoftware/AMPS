#!/usr/bin/env python3
"""Actual accepted-flight controller checks, reusing a fresh native producer.

The declared spatial hazard and zero-Q packet are analytical controls, never
physical Na products or observational validation. No configuration/build is
changed here; the caller owns its single deadline and ordinary restoration.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import re
import sys

from campaign_budget import CampaignBudget
from build_m2_qualification import production_source_identity,sha256
from integration_common import REPO
from run_m3_moon_force_tracers import native_sources
sys.path.insert(0,str(Path(__file__).resolve().parent/'common'))
from moon_testlib import _probe_compile_command

def run(output,producer_manifest,wall_seconds=600):
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=False)
    producer_manifest=Path(producer_manifest).resolve()
    result={'status':'ERROR','scope':'actual PIC variable-exposure mechanics only; NOT physical M3 qualification',
        'counts':{s:0 for s in ('PASS','FAIL','ERROR','SKIPPED')},'results':[],'commands':[]}
    manifest={'MPI_layouts':[1,8],'orders':[1,2],'OMP_NUM_THREADS':1,
        'hazard':'explicit lambda(x)=1+(x-x0)/730 s^-1; zero-force mathematical control',
        'kinematics':'nominal equal-heavy-mass zero-Q analytical packet, not Na physics',
        'external_observations':[],'random_seed':None}
    budget=CampaignBudget(wall_seconds,min(30,wall_seconds*.1))
    try:
        if not producer_manifest.is_file():
            result.update(status='SKIPPED',reason='fresh actual all-species native producer unavailable');result['counts']['SKIPPED']=1;return result
        producer=json.loads(producer_manifest.read_text());prior=json.loads(producer_manifest.with_name('result.json').read_text())
        if prior.get('status')!='PASS':
            result.update(status='SKIPPED',reason='actual native backend predecessor is not PASS');result['counts']['SKIPPED']=1;return result
        if Path(producer['source_root']).resolve()!=REPO:raise ValueError('native producer belongs to another checkout')
        library=Path(producer['library']['path']);build=library.parent
        if sha256(library)!=producer['library']['sha256']:raise ValueError('actual native library changed')
        for path,digest in producer['generated_snapshots'].items():
            relative=Path(path).relative_to(producer_manifest.parent/'generated')
            if sha256(Path(path))!=digest or sha256(build/relative)!=digest:raise ValueError('producer/current generated ABI mismatch')
        if native_sources(producer['production_source_hashes'])!=native_sources(production_source_identity()):
            raise ValueError('trajectory controller sources differ from actually built native producer')
        manifest.update(producer_manifest={'path':str(producer_manifest),'sha256':sha256(producer_manifest)},
            library=producer['library'],production_source_hashes=producer['production_source_hashes'],generated_snapshots=producer['generated_snapshots'])
        here=Path(__file__).resolve().parent
        sources=[Path(__file__).resolve(),here/'m3_native_trajectory_exposure_probe.cpp',here/'m3_native_reaction_controller_probe.cpp',here/'M3_TRAJECTORY_EXPOSURE_CONTRACT.md']
        manifest['test_sources']={str(p):sha256(p) for p in sources}
        executable=output/'probe';command=_probe_compile_command(build,executable);command.insert(1,'-DMPI_ON')
        old=next(x for x in command if x.endswith('/common/production_kernel_probe.cpp'));command[command.index(old)]=str(sources[1])
        def execute(argv,label,cwd,ranks):
            receipt=budget.run(argv,cwd,output/(label+'.log'),output=output,allocated_cores=ranks,
                environment=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
            result['commands'].append(receipt)
            if receipt.get('reason') or receipt.get('returncode') not in (0,1):raise RuntimeError(label+' infrastructure failure')
            return receipt
        if execute(command,'compile',REPO,1)['returncode']:raise RuntimeError('actual native trajectory probe compilation failed')
        manifest['executable']={'path':str(executable),'sha256':sha256(executable)}
        for ranks in (1,8):
            for order in (1,2):
                label=f'np{ranks}-order{order}';directory=output/label;directory.mkdir();(directory/'PT').mkdir()
                receipt=execute(['mpiexec','-n',str(ranks),str(executable),str(order),str(directory/'trajectory')],label,directory,ranks)
                with (directory/'trajectory.checks.csv').open(newline='') as stream:rows=list(csv.DictReader(stream))
                text=(output/(label+'.log')).read_text();match=re.findall(r'M3 trajectory exposure order=(\d+) checks=(\d+) status=(PASS|FAIL)',text)
                if len(rows)!=14 or len({r['check'] for r in rows})!=14 or len(match)!=1 or int(match[0][0])!=order or int(match[0][1])!=14:
                    raise ValueError('actual native trajectory check membership/summary incomplete')
                if any(row['status'] not in ('PASS','FAIL') for row in rows):raise ValueError('native trajectory result status invalid')
                failed=any(row['status']=='FAIL' for row in rows)
                if receipt['returncode']!=int(failed) or match[0][2]!=('FAIL' if failed else 'PASS'):raise ValueError('native status/exit contradiction')
                for row in rows:
                    result['counts'][row['status']]+=1
                    result['results'].append(dict(row,mpi_ranks=ranks,order=order,log=str(output/(label+'.log'))))
        result['status']='FAIL' if result['counts']['FAIL'] else 'PASS'
    except (OSError,ValueError,RuntimeError,KeyError) as error:
        result.update(status='ERROR',reason=str(error));result['counts']['ERROR']+=1
    finally:
        (output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
        (output/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',required=True,type=Path)
    p.add_argument('--producer-manifest',required=True,type=Path);p.add_argument('--wall-seconds',type=float,default=600)
    args=p.parse_args();value=run(args.output_dir,args.producer_manifest,args.wall_seconds);print(json.dumps(value,indent=2))
    sys.exit({'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[value['status']])
