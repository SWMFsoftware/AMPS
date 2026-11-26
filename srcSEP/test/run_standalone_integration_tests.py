#!/usr/bin/env python3
"""Run production core adapters in an isolated, clean repository-shaped tree.

This runner owns numerical/software verification only. It never creates an SEP
injection prescription or promotes a pure test into native/MPI qualification.
JSON is written on build failure as well as assertion failure so orchestration
cannot mistake missing output for successful physical evaluation.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
import time
import xml.etree.ElementTree as ET

REPO = Path(__file__).resolve().parents[2]
IDS = ['SINT01', 'SINT02', 'SINT03', 'SINT04', 'SINT05', 'SINT06', 'SINT07', 'SINT07R', 'SINT10', 'SINT11', 'SINT12', 'SINT13', 'SINT14', 'SINT15', 'SINT16']
IDS += ['SINT17', 'SINT18', 'SINT19', 'SINT20', 'SINT21', 'SINT22', 'SINT23', 'SINT24', 'SINT26', 'SINT27', 'SINT28', 'SINT29', 'SINT30', 'SINT31', 'SINT32']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--sanitizers', action='store_true', help='instrument adapter, geometry, transport and changed core source/jump/ambient paths with ASan/UBSan')
    args = parser.parse_args()
    reports=Path(os.environ['SRCSEP_REPORT_DIR']).resolve() if os.environ.get('SRCSEP_REPORT_DIR') else None
    if reports:reports.mkdir(parents=True,exist_ok=True)
    output=args.output_dir.resolve() if args.output_dir else Path(tempfile.mkdtemp(
        prefix='standalone-core-',dir=str(reports) if reports else None))
    output.mkdir(parents=True, exist_ok=True)
    root = output/'source'
    if root.exists():
        parser.error('output source tree already exists; use a fresh evidence directory')
    started = time.monotonic()
    commands = []
    rows = []
    ignore = shutil.ignore_patterns('build', '__pycache__', '*.o', '*.a', '*.d', '*.pyc', 'test_output', 'output')

    def copy(source, target):
        shutil.copytree(source, target, ignore=ignore)

    def execute(argv, cwd, name):
        log = output/name
        commands.append({'argv': argv, 'cwd': str(cwd), 'log': str(log)})
        with log.open('w') as stream:
            result = subprocess.run(argv, cwd=cwd, stdout=stream, stderr=subprocess.STDOUT, timeout=600,
                                    env=dict(os.environ, ASAN_OPTIONS='detect_leaks=0') if args.sanitizers else None)
        return result.returncode, log

    try:
        for name in ['sep_common', 'sep_coronal_cme', 'sep_corona_swcme', 'swcme']:
            copy(REPO/'src/models'/name, root/'src/models'/name)
        copy(REPO/'src/general/input', root/'src/general/input')
        copy(REPO/'srcSEP/adapters', root/'srcSEP/adapters')
        (root/'srcSEP/util').mkdir()
        shutil.copy2(REPO/'srcSEP/util/sep_common_header_path.h', root/'srcSEP/util')
        for name in ['sep_reduced_geometry.h', 'sep_reduced_geometry.cpp', 'sep_geometric_spacing.h',
                     'sep_reduced_particle_options.h', 'sep_reduced_source.h', 'sep_reduced_source.cpp',
                     'sep_absorbing_bridge.h', 'sep_absorbing_bridge.cpp',
                     'sep_parker_core.h', 'sep_parker_core.cpp',
                     'sep_focused_transport_core.h', 'sep_focused_transport_core.cpp',
                     'sep_focused_transport_mfp_core.h', 'sep_focused_transport_mfp_core.cpp',
                     'sep_transport_coefficients.h']:
            shutil.copy2(REPO/'srcSEP/util'/name, root/'srcSEP/util'/name)
        (root/'src/pic').mkdir()
        shutil.copy2(REPO/'src/pic/pic_field_line_geometry.h', root/'src/pic')
        (root/'srcSEP/test').mkdir()
        shutil.copy2(REPO/'srcSEP/test/test_standalone_core_integration.cpp', root/'srcSEP/test')
        copy(REPO/'srcSEP/examples/reduced-shock', root/'srcSEP/examples/reduced-shock')
        fixture = root/'srcSEP/examples/reduced-shock/positive_1au.event'
        alpha = fixture.with_name('alpha-test.event')
        original = fixture.read_text()
        needle = 'plasma.alpha_to_proton_number_ratio = 0\n'
        if original.count(needle) != 1:
            raise RuntimeError('composition fixture no longer has its declared pure-H input')
        alpha.write_text(original.replace(needle, 'plasma.alpha_to_proton_number_ratio = 0.1\n'))
        # Capture exactly the immutable bytes being compiled, before generation
        # of any model archives. Generated objects never serve as source input.
        hashes = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                  for p in root.rglob('*') if p.is_file()}
        for model in ['sep_coronal_cme', 'sep_corona_swcme']:
            code, log = execute(['make', '-j4', 'lib'], root/'src/models'/model, model+'-build.log')
            if code:
                raise RuntimeError(f'{model} library build failed: {log}')
        includes = ['srcSEP/adapters', 'srcSEP/util', 'src/pic', 'src/models/sep_corona_swcme/include',
                    'src/models/sep_corona_swcme/shock_front', 'src/models/sep_coronal_cme/include',
                    'src/models/sep_common']
        argv = [os.environ.get('CXX', 'c++'), '-std=c++17', '-O2', '-Wall', '-Wextra', '-Wpedantic', '-Werror']
        if args.sanitizers:
            argv += ['-fsanitize=address,undefined', '-fno-omit-frame-pointer']
        argv += ['-I'+str(root/p) for p in includes]
        argv += [str(root/'srcSEP/adapters/reduced_shock_background_adapter.cpp'),
                 str(root/'srcSEP/util/sep_reduced_geometry.cpp'),
                 str(root/'srcSEP/util/sep_reduced_source.cpp'),
                 str(root/'srcSEP/util/sep_absorbing_bridge.cpp'),
                 str(root/'src/models/sep_common/sep_background_snapshot.cpp'),
                 str(root/'src/models/sep_common/sep_transport_common.cpp'),
                 str(root/'srcSEP/util/sep_parker_core.cpp'),
                 str(root/'srcSEP/util/sep_focused_transport_core.cpp'),
                 str(root/'srcSEP/util/sep_focused_transport_mfp_core.cpp'),
                 str(root/'srcSEP/test/test_standalone_core_integration.cpp')]
        if args.sanitizers:
            # Put the changed core translation units before their ordinary
            # archives. The linker resolves their exported symbols from these
            # instrumented objects and does not extract the same archive
            # members. This covers the actual source CDF, RH search and ambient
            # shortcut without claiming instrumentation of every unchanged
            # core member or of the separately linked PIC driver.
            argv += [str(root/'src/models/sep_coronal_cme/src/particle_source.cpp'),
                     str(root/'src/models/sep_coronal_cme/src/mhd_jump_solver.cpp'),
                     str(root/'src/models/sep_corona_swcme/src/ambient_state.cpp')]
        argv += ['-Wl,--start-group', str(root/'src/models/sep_corona_swcme/build/libsep_corona_swcme.a'),
                 str(root/'src/models/sep_coronal_cme/build/libsep_coronal_cme.a'),
                 '-Wl,--end-group', '-o', str(output/'integration-test')]
        code, log = execute(argv, root, 'fixture-build.log')
        if code:
            raise RuntimeError(f'fixture build failed: {log}')
        code, log = execute([str(output/'integration-test'), str(fixture), str(alpha)], root, 'fixture-run.log')
        records = re.findall(r'^\[(SINT\d+R?)\] (PASS|FAIL) (.*)$', log.read_text(), re.M)
        if [r[0] for r in records] != IDS:
            raise RuntimeError('missing, duplicated or reordered integration result identities')
        metrics={identity:[] for identity in IDS}
        for identity,fields in re.findall(r'^METRIC (SINT\d+R?) (.*)$',log.read_text(),re.M):
            metrics[identity].append(dict(token.split('=',1) for token in fields.split() if '=' in token))
        rows = [{'scope': 'production-adapter', 'test_id': identity, 'status': status,
                 'message': message, 'log': str(log), 'metrics':metrics[identity],
                 'seed':{'SINT14':4289,'SINT16':5297,'SINT22':773,'SINT23':997}.get(identity),
                 'configuration':{'units':'SI/HCI','event':str(fixture),'alpha_event':str(alpha),
                     'source_sha256':str(output/'source-sha256.json'),'sanitizers':args.sanitizers},
                 'additional_seeds':[20261122] if identity=='SINT23' else [],
                 'acceptance':'Frozen analytical/statistical assertions in source/srcSEP/test/test_standalone_core_integration.cpp',
                 'artifacts':[str(log),str(output/'commands.json')]} for identity, status, message in records]
        if code != int(any(r['status'] == 'FAIL' for r in rows)):
            raise RuntimeError('test executable exit disagrees with its case results')
        for path, digest in hashes.items():
            if hashlib.sha256((root/path).read_bytes()).hexdigest() != digest:
                raise RuntimeError('compiled input changed during execution: '+path)
        (output/'source-sha256.json').write_text(json.dumps(hashes, indent=2)+'\n')
    except Exception as error:
        rows = [{'scope': 'production-adapter', 'test_id': identity, 'status': 'ERROR',
                 'message': str(error), 'log': str(output/'fixture-build.log')} for identity in IDS]
    for row in rows:
        print(f"[{row['test_id']}] {row['status']} {row['message']}", flush=True)
    (output/'commands.json').write_text(json.dumps(commands, indent=2)+'\n')
    (output/'summary.json').write_text(json.dumps({'schema': 'srcsep-standalone-integration-v1',
        'scope': 'production-adapter; not native/MPI or observational qualification',
        'seconds': time.monotonic()-started, 'results': rows}, indent=2)+'\n')
    suite=ET.Element('testsuite', name='standalone-core-integration', tests=str(len(rows)))
    for row in rows:
        node=ET.SubElement(suite,'testcase',name=row['test_id'])
        if row['status']!='PASS':
            ET.SubElement(node,'error' if row['status']=='ERROR' else 'failure',message=row['message'])
    ET.ElementTree(suite).write(output/'junit.xml',encoding='unicode')
    if reports:
        imported=[dict(row,id=row['test_id'],group='standalone-core-integration',artifacts=[str(output/'summary.json')]) for row in rows]
        (reports/'standalone-integration-results.json').write_text(json.dumps({'schema':'srcsep-component-tests-v1','results':imported},indent=2)+'\n')
    print('Evidence: '+str(output/'summary.json'), flush=True)
    return int(any(row['status'] in ['FAIL', 'ERROR'] for row in rows))


if __name__ == '__main__':
    raise SystemExit(main())
