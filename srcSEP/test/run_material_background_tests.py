#!/usr/bin/env python3
"""Strict/sanitized shared material-contract and real canonical-adapter gates.

These source tests exercise actual core queries and interface error paths.
They do not supply new background evolution laws or replace native MPI gates.
Every build/runtime/report failure emits the complete twenty-case ERROR catalog.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time
import xml.etree.ElementTree as ET

REPO=Path(__file__).resolve().parents[2]
IDS=['MBC%02d'%i for i in range(1,21)]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path)
    parser.add_argument('--sanitizers',action='store_true')
    args=parser.parse_args()
    reports=Path(os.environ['SRCSEP_REPORT_DIR']).resolve() if os.environ.get('SRCSEP_REPORT_DIR') else None
    if reports:reports.mkdir(parents=True,exist_ok=True)
    out=args.output_dir.resolve() if args.output_dir else Path(tempfile.mkdtemp(
        prefix='material-background-',dir=str(reports) if reports else None))
    # Never overwrite an earlier successful or failed scientific record.
    if args.output_dir and out.exists():
        print('ERROR: choose a fresh material-background evidence directory',flush=True)
        return 1
    out.mkdir(parents=True,exist_ok=True)
    root=out/'source';commands=[];hashes={};rows=[];started=time.monotonic()
    def execute(argv,cwd,name,timeout=600):
        # Record the attempted command before launch. A missing compiler or
        # timeout must retain its real argv and error, with no invented exit.
        log=out/name
        receipt=dict(argv=argv,cwd=str(cwd),exit_code=None,log=str(log))
        commands.append(receipt)
        launch=time.monotonic()
        try:
            with log.open('w') as stream:
                result=subprocess.run(argv,cwd=cwd,stdout=stream,stderr=subprocess.STDOUT,
                    timeout=timeout,env=dict(os.environ,ASAN_OPTIONS='detect_leaks=0') if args.sanitizers else None)
            receipt['exit_code']=result.returncode
            return result.returncode,log
        except (OSError,subprocess.TimeoutExpired) as error:
            receipt['error']=str(error)
            raise
        finally:
            receipt['seconds']=time.monotonic()-launch
    try:
        ignore=shutil.ignore_patterns('build','test_output','__pycache__','*.o','*.a','*.d','*.pyc')
        for name in ('sep_common','sep_coronal_cme','sep_corona_swcme','swcme'):
            shutil.copytree(REPO/'src/models'/name,root/'src/models'/name,ignore=ignore)
        shutil.copytree(REPO/'src/general/input',root/'src/general/input',ignore=ignore)
        shutil.copytree(REPO/'srcSEP/adapters',root/'srcSEP/adapters',ignore=ignore)
        (root/'srcSEP/util').mkdir()
        for name in ('sep_common_header_path.h','sep_reduced_geometry.h','sep_reduced_geometry.cpp','sep_geometric_spacing.h'):
            shutil.copy2(REPO/'srcSEP/util'/name,root/'srcSEP/util'/name)
        (root/'src/pic').mkdir()
        shutil.copy2(REPO/'src/pic/pic_field_line_geometry.h',root/'src/pic')
        (root/'srcSEP/test').mkdir()
        shutil.copy2(REPO/'srcSEP/test/test_material_background.cpp',root/'srcSEP/test')
        shutil.copytree(REPO/'srcSEP/examples/reduced-shock',root/'srcSEP/examples/reduced-shock',ignore=ignore)
        event=root/'srcSEP/examples/reduced-shock/positive_1au.event'
        original=event.read_text();needle='plasma.alpha_to_proton_number_ratio = 0\n'
        if original.count(needle)!=1:raise ValueError('explicit composition fixture changed; no guessed replacement')
        alpha=event.with_name('material-alpha-fixture.event')
        alpha.write_text(original.replace(needle,'plasma.alpha_to_proton_number_ratio = 0.1\n'))
        # Hash all actual copied inputs BEFORE any generated objects exist.
        # This distinguishes immutable originals from the labelled composition
        # verification fixture and binds the actual compiled source closure.
        hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest()
                for p in root.rglob('*') if p.is_file()}
        (out/'source-sha256.json').write_text(json.dumps(hashes,indent=2)+'\n')
        for name in ('sep_coronal_cme','sep_corona_swcme'):
            code,log=execute(['make','-j16','lib'],root/'src/models'/name,name+'-build.log')
            if code:raise RuntimeError('canonical library build failed: '+str(log))
        includes=['srcSEP/adapters','srcSEP/util','src/pic','src/models/sep_common','src/models/swcme',
                  'src/models/sep_coronal_cme/include','src/models/sep_corona_swcme/include',
                  'src/models/sep_corona_swcme/shock_front']
        # A neutral public header must compile by itself, without an earlier
        # PIC/model include accidentally supplying its required declarations.
        code,log=execute([os.environ.get('CXX','c++'),'-std=c++17','-Wall','-Wextra',
            '-Wpedantic','-Werror','-x','c++','-fsyntax-only','-include',
            str(root/'src/models/sep_common/sep_material_background.h'),'/dev/null'],
            root,'public-header-build.log')
        if code:raise RuntimeError('self-contained public material header failed: '+str(log))
        argv=[os.environ.get('CXX','c++'),'-std=c++17','-O2','-Wall','-Wextra','-Wpedantic','-Werror']
        if args.sanitizers:argv+=['-fsanitize=address,undefined','-fno-omit-frame-pointer']
        argv+=['-I'+str(root/p) for p in includes]
        argv+=[str(root/p) for p in ['srcSEP/test/test_material_background.cpp',
             'srcSEP/adapters/reduced_shock_background_adapter.cpp','srcSEP/adapters/swcme1d_adapter.cpp',
             'srcSEP/util/sep_reduced_geometry.cpp','src/models/sep_common/sep_background_snapshot.cpp']]
        argv+=['-Wl,--start-group',str(root/'src/models/sep_corona_swcme/build/libsep_corona_swcme.a'),
               str(root/'src/models/sep_coronal_cme/build/libsep_coronal_cme.a'),'-Wl,--end-group',
               '-o',str(out/'material-background-test')]
        code,log=execute(argv,root,'fixture-build.log')
        if code:raise RuntimeError('strict contract fixture compilation failed: '+str(log))
        code,log=execute([str(out/'material-background-test'),str(event),str(alpha)],root,'fixture-run.log',120)
        records=re.findall(r'^\[(MBC\d+)\] (PASS|FAIL)(.*)$',log.read_text(),re.M)
        if [x[0] for x in records]!=IDS:raise RuntimeError('missing/foreign/duplicate material contract cases')
        # Import actual per-case measurements from the C++ assertions. Do not
        # divide total build time among cases or invent zero residuals for
        # cases without numerical reference comparisons. The bound is the
        # unchanged value used by the numerical test, not a runner override.
        measured=re.findall(r'^\[(MBC\d+)\] METRICS assertions=(\d+) near_calls=(\d+) maximum_relative_error=(\S+) relative_tolerance=(\S+) seconds=(\S+)$',log.read_text(),re.M)
        if [x[0] for x in measured]!=IDS:raise RuntimeError('missing/foreign/duplicate material case measurements')
        rows=[]
        for (identity,status,message),metric in zip(records,measured):
            _,checks,near_calls,residual,bound,seconds=metric
            checks=int(checks);near_calls=int(near_calls);residual=float(residual);bound=float(bound);seconds=float(seconds)
            if checks<1 or not math.isfinite(seconds) or seconds<0 or not math.isfinite(bound) or bound<=0:
                raise RuntimeError('invalid actual material assertion count/timing/bound: '+identity)
            rows.append(dict(test_id=identity,status=status,message=message.strip() or
                'shared contract and canonical adapter assertion passed',seed=0,log=str(log),seconds=seconds,
                metrics=dict(assertions=checks,near_calls=near_calls,
                    maximum_relative_error=residual if near_calls and math.isfinite(residual) else None),
                tolerance=dict(metadata='exact',SI_reference_relative=bound if near_calls else None),
                configuration=dict(compiled_source_root=str(root),source_hashes=str(out/'source-sha256.json'),
                    evidence_scope='software contracts and canonical adapter verification; no new physics admission')))
        if code!=int(any(row['status']=='FAIL' for row in rows)):
            raise RuntimeError('contract test exit disagrees with per-case assertions')
        # Detect a changed copy AND a changed primary source during execution.
        # The sole generated physical fixture is not claimed as an original.
        for name,sha in hashes.items():
            if hashlib.sha256((root/name).read_bytes()).hexdigest()!=sha:
                raise RuntimeError('compiled material source changed during execution: '+name)
            if name.endswith('material-alpha-fixture.event'):continue
            original_path=REPO/name
            if original_path.is_file() and hashlib.sha256(original_path.read_bytes()).hexdigest()!=sha:
                raise RuntimeError('checkout material source changed during execution: '+name)
    except (OSError,ValueError,RuntimeError,subprocess.TimeoutExpired) as error:
        rows=[dict(test_id=identity,status='ERROR',message=str(error)) for identity in IDS]
    payload=dict(schema='srcsep-material-background-v1',scope='stage11 software contracts and canonical adapters; no new backend physics admission',
                 sanitizers=args.sanitizers,results=rows,commands=commands,source_sha256=hashes,seconds=time.monotonic()-started)
    (out/'summary.json').write_text(json.dumps(payload,indent=2)+'\n')
    suite=ET.Element('testsuite',name='material-background',tests=str(len(rows)))
    for row in rows:
        print('[%s] %s %s'%(row['test_id'],row['status'],row['message']),flush=True)
        case=ET.SubElement(suite,'testcase',name=row['test_id'])
        if 'seconds' in row:case.set('time',str(row['seconds']))
        properties=ET.SubElement(case,'properties')
        for name in ('seed','metrics','tolerance','configuration'):
            if name in row:ET.SubElement(properties,'property',name=name,value=json.dumps(row[name],sort_keys=True))
        if row['status']!='PASS':ET.SubElement(case,'error' if row['status']=='ERROR' else 'failure',message=row['message'])
    ET.ElementTree(suite).write(out/'junit.xml',encoding='unicode')
    if reports:(reports/'material-background-results.json').write_text(json.dumps({'results':[dict(row,id=row['test_id'],group='material-background') for row in rows]},indent=2)+'\n')
    return int(any(row['status']!='PASS' for row in rows))

if __name__=='__main__':
    raise SystemExit(main())
