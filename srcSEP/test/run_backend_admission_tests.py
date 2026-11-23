#!/usr/bin/env python3
"""Verify existing backend boundaries without qualifying unsupported transport.

Source mode compares actual canonical SWCME queries and neutral exchange state
with the production admission checks. Native mode separately requires all three
standalone solvers to reject SWCME before allocation on two real MPI ranks.
These are rejection/compatibility results, not successful backend propagation.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import time
import xml.etree.ElementTree as ET

REPO = Path(__file__).resolve().parents[2]
SOURCE_IDS = ['SBACK%02d' % i for i in range(1, 6)]
NATIVE_IDS = ['SBACK%02d' % i for i in range(6, 9)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--sanitizers', action='store_true')
    parser.add_argument('--native', action='store_true')
    parser.add_argument('--amps', type=Path)
    action = parser.add_mutually_exclusive_group()
    action.add_argument('--prepare-only', action='store_true', help='native MPI jobs only; no PASS is emitted')
    action.add_argument('--collect', action='store_true', help='collect externally executed native jobs with actual process receipts')
    args = parser.parse_args()
    reports = Path(os.environ['SRCSEP_REPORT_DIR']).resolve() if os.environ.get('SRCSEP_REPORT_DIR') else None
    if reports:
        reports.mkdir(parents=True, exist_ok=True)
    root = args.output_dir.resolve() if args.output_dir else Path(tempfile.mkdtemp(
        prefix='backend-admission-', dir=str(reports) if reports else None))
    # An explicitly reused directory can contain valid prior evidence. Refuse
    # preparation/execution without replacing its verdict or starting MPI.
    # The implicit mkdtemp directory above is newly owned; --collect explicitly
    # re-evaluates its existing plan and remains a separate operation.
    if args.output_dir and root.exists() and not args.collect:
        print('ERROR: choose a fresh backend evidence directory: '+str(root), file=sys.stderr)
        return 1
    root.mkdir(parents=True, exist_ok=True)
    identities = NATIVE_IDS if args.native else SOURCE_IDS
    commands, hashes, rows = [], {}, []

    def execute(argv, name, cwd, timeout):
        """Capture actual return codes and commands, including expected rejects."""
        log = root/name
        started = time.monotonic()
        if args.collect:
            receipt = json.loads((root/(name+'.process.json')).read_text())
            if receipt.get('argv') != argv or receipt.get('cwd') != str(cwd) or type(receipt.get('exit_code')) is not int:
                raise ValueError('missing or foreign native process receipt: '+name)
            code = receipt['exit_code']
        else:
            with log.open('w') as stream:
                result = subprocess.run(argv, cwd=cwd, stdout=stream, stderr=subprocess.STDOUT,
                    timeout=timeout, env=dict(os.environ, ASAN_OPTIONS='detect_leaks=0') if args.sanitizers else None)
            code = result.returncode
        commands.append({'argv': argv, 'cwd': str(cwd), 'exit_code': code,
                         'seconds': time.monotonic()-started, 'log': str(log)})
        return code, log.read_text()

    try:
        if (args.prepare_only or args.collect) and not args.native:
            raise ValueError('external preparation/collection requires --native')
        paths = [Path(__file__), REPO/'srcSEP/util/sep_backend_admission.h',
                 REPO/'srcSEP/adapters/field_line_bundle_adapter.h',
                 REPO/'srcSEP/adapters/field_line_bundle_adapter.cpp',
                 REPO/'srcSEP/adapters/swcme1d_adapter.h', REPO/'srcSEP/adapters/swcme1d_adapter.cpp',
                 REPO/'srcSEP/util/sep_common_header_path.h', REPO/'srcSEP/test/test_backend_admission.cpp']
        paths += sorted((REPO/'src/models/swcme').glob('*.hpp'))
        paths += sorted((REPO/'src/models/sep_common').glob('*.h'))
        sources = [REPO/p for p in ['srcSEP/test/test_backend_admission.cpp',
            'srcSEP/adapters/swcme1d_adapter.cpp', 'srcSEP/adapters/field_line_bundle_adapter.cpp',
            'src/models/sep_common/sep_background_snapshot.cpp',
            'src/models/sep_common/sep_field_line_exchange.cpp', 'src/models/sep_common/sep_field_line_bundle_io.cpp']]
        paths += sources
        if args.native:
            if not args.amps:
                raise ValueError('--native requires a rebuilt --amps executable')
            executable = args.amps.resolve()
            paths += [executable, REPO/'srcSEP/main.cpp', REPO/'srcSEP/transport_common.cpp',
                      REPO/'srcSEP/examples/reduced-shock/field_line.in']
        hashes = {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
        if args.native:
            jobs = []
            if args.collect:
                frozen = json.loads((root/'prepared.json').read_text())
                if frozen['hashes'] != hashes:
                    raise ValueError('native admission source/executable identity changed')
            for identity, mode in zip(identities, ('parker', 'fte-dmumu', 'fte-mfp')):
                directory = root/mode
                directory.mkdir(exist_ok=args.collect)
                argv = ['mpiexec', '-n', '2', '--bind-to', 'none', str(executable),
                    '--input', str(REPO/'srcSEP/examples/reduced-shock/field_line.in'),
                    '--shock-model', 'swcme1d', '--particle-mover', mode, '--total-iterations', '1',
                    '--no-test-manager', '--turbulence-source', 'prescribed',
                    '--turbulence-coupling-policy', 'disabled', '--coupling', 'off', '--cascade', 'off',
                    '--reflection', 'off', '--turbulence-advection', 'off', '--shock-injection', 'off']
                jobs.append({'test_id': identity, 'argv': argv, 'cwd': str(directory),
                             'log': str(root/(identity+'.log')), 'receipt': str(root/(identity+'.log.process.json'))})
                if args.prepare_only:
                    continue
                code, text = execute(argv, identity+'.log', directory, 60)
                good = code == 1 and 'SWCME1D moving SEP transport is unsupported:' in text
                good = good and 'Iteration:' not in text and not (directory/'PT').exists()
                # If the intended diagnostic was never reached (MPI startup,
                # parser/build/input failure), there was no evaluated backend
                # rejection. Keep that infrastructure outcome as ERROR.
                evaluated = 'SWCME1D moving SEP transport is unsupported:' in text
                rows.append({'test_id': identity, 'status': 'PASS' if good else 'FAIL' if evaluated else 'ERROR',
                    'message': 'actual two-rank '+mode+' rejects unsupported material backend before allocation',
                    'actual_exit_code': code, 'log': str(root/(identity+'.log')), 'seed': 0})
            if args.prepare_only:
                (root/'prepared.json').write_text(json.dumps({'hashes': hashes, 'jobs': jobs}, indent=2)+'\n')
                print('PREPARED, not executed: '+str(root/'prepared.json'), flush=True)
                return 0
            if args.collect and frozen['jobs'] != jobs:
                raise ValueError('native admission launch catalog changed')
        else:
            executable = root/'backend-test'
            argv = [os.environ.get('CXX', 'c++'), '-std=c++17', '-O2', '-Wall', '-Wextra', '-Wpedantic', '-Werror']
            if args.sanitizers:
                argv += ['-fsanitize=address,undefined', '-fno-omit-frame-pointer']
            argv += ['-I'+str(REPO/p) for p in ['srcSEP/util', 'srcSEP/adapters', 'src/models/sep_common', 'src/models/swcme']]
            argv += [str(path) for path in sources]+['-o', str(executable)]
            code, _ = execute(argv, 'build.log', REPO, 300)
            if code:
                raise RuntimeError('strict backend fixture compilation failed; see build.log')
            code, text = execute([str(executable), str(root)], 'native-kernel.log', root, 60)
            records = re.findall(r'^\[(SBACK\d+)\] (PASS|FAIL)(.*)$', text, re.M)
            if [row[0] for row in records] != identities:
                raise RuntimeError('backend fixture has missing/foreign/duplicate cases')
            rows = [{'test_id': identity, 'status': status, 'message': message.strip() or
                'canonical query/exchange and material admission boundary checked',
                'log': str(root/'native-kernel.log'), 'seed': 0} for identity, status, message in records]
            if code != int(any(row['status'] == 'FAIL' for row in rows)):
                raise RuntimeError('fixture exit disagrees with actual assertions')
        for name, before in hashes.items():
            if hashlib.sha256(Path(name).read_bytes()).hexdigest() != before:
                raise RuntimeError('backend source/input changed during execution: '+name)
    except (OSError, ValueError, KeyError, RuntimeError, subprocess.TimeoutExpired) as error:
        rows = [{'test_id': identity, 'status': 'ERROR', 'message': str(error)} for identity in identities]
    payload = {'schema': 'srcsep-backend-admission-v1', 'native': args.native,
               'scope': 'backend rejection and diagnostic compatibility; no unsupported transport qualification',
               'results': rows, 'commands': commands, 'source_sha256': hashes}
    (root/'summary.json').write_text(json.dumps(payload, indent=2)+'\n')
    suite = ET.Element('testsuite', name='backend-admission', tests=str(len(rows)))
    for row in rows:
        print('[%s] %s %s' % (row['test_id'], row['status'], row['message']), flush=True)
        node = ET.SubElement(suite, 'testcase', name=row['test_id'])
        if row['status'] != 'PASS':
            ET.SubElement(node, 'error' if row['status'] == 'ERROR' else 'failure', message=row['message'])
    ET.ElementTree(suite).write(root/'junit.xml', encoding='unicode')
    if reports:
        (reports/'backend-admission-results.json').write_text(json.dumps({'results': [dict(row,
            id=row['test_id'], group='backend-admission') for row in rows]}, indent=2)+'\n')
    return int(any(row['status'] != 'PASS' for row in rows))


if __name__ == '__main__':
    raise SystemExit(main())
