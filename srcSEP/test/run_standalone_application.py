#!/usr/bin/env python3
"""Run the maintained 400-step, three-mover reduced-core application gate.

The physical deck, 300-s timestep, 400-sample output window, 1800-s deadline
and independent verifier are the existing stage-6 contract. Preparation only
freezes commands and inputs; it never reports a simulation PASS. External MPI
launch is supported for hosts which prohibit MPI sockets in Python children.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time

import verify_standalone_application as audit

REPO = Path(__file__).resolve().parents[2]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, payload):
    """Publish complete evidence atomically; retain every failed run directory."""
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(payload, indent=2) + '\n')
    temporary.replace(path)


def jobs_for(root):
    """Freeze the previously qualified physical campaign, not a tunable smoke run.

    All shell arguments are individually quoted. The only shell expansion is
    Open MPI's rank in the output filename; no user-provided text is executed
    as shell code. No step-count, timestep or tolerance override is offered.
    """
    common = [str(root/'amps'), '--input', str(REPO/'srcSEP/examples/reduced-shock/field_line.in'),
              '--reduced-shock-event', str(REPO/'srcSEP/examples/reduced-shock/positive_1au.event'),
              '--reduced-shock-particles', 'on', '--shock-release-distance-m', '1e8',
              '--shock-return-policy', 'delete', '--shock-source-normalization', 'net-upstream-escape',
              '--total-iterations', '400', '--no-test-manager', '--turbulence-source', 'prescribed',
              '--turbulence-coupling-policy', 'disabled', '--coupling', 'off', '--cascade', 'off',
              '--reflection', 'off', '--turbulence-advection', 'off', '--shock-injection', 'off',
              '--coefficient-source', 'prescribed', '--pitch-angle-diffusion-provider', 'constant',
              '--constant-dmumu', '0.001', '--spatial-diffusion-provider', 'from-dmumu',
              '--mean-free-path-provider', 'from-spatial']
    jobs = []
    for mode in audit.MODES:
        child = 'exec timeout 1800s ' + shlex.join(common + ['--particle-mover', mode])
        child += ' > application.rank${OMPI_COMM_WORLD_RANK}.log 2>&1'
        argv = ['mpiexec', '-n', '1', '--bind-to', 'none', '/bin/sh', '-c', child]
        jobs.append({'mode': mode, 'cwd': str(root/mode), 'argv': argv, 'command': shlex.join(argv)})
    return jobs


def check_prepared(root):
    """Reject stale/partial/foreign plans before launch and again at collection.

    Exact command and catalog identity prevents a short or differently selected
    job from satisfying this gate. Hash checks bind the actual native closure,
    not simply a later executable with the same filename.
    """
    plan = json.loads((root/'launch-plan.json').read_text())
    if plan.get('schema') != 'srcsep-full-application-plan-v1' or plan.get('jobs') != jobs_for(root):
        raise ValueError('full application plan/catalog differs from the maintained physical campaign')
    hashes = json.loads((root/'build-sha256.json').read_text())
    if not hashes:
        raise ValueError('missing source/build identity closure')
    for name, before in hashes.items():
        if digest(Path(name)) != before:
            raise ValueError('prepared application source/build/input changed: ' + name)
    return plan


def prepare(root, build):
    """Use the existing native dependency audit before admitting a real driver.

    Its prepared fixture supplies compiler dependency/object/archive identity
    and checkout/generated-source equivalence. It does not run the application
    and its successful preparation is never substituted for the 400-step gate.
    The ordinary standalone main is additionally checked because fixture mains
    deliberately do not link that object.
    """
    root.mkdir(parents=True)
    fixture = root/'native-prerequisite'
    command = [sys.executable, str(REPO/'srcSEP/test/run_native_standalone_integration.py'),
               '--build-root', str(build), '--output-dir', str(fixture), '--ranks', '1', '--prepare-only']
    with (root/'prepare.log').open('w') as stream:
        completed = subprocess.run(command, cwd=REPO, stdout=stream, stderr=subprocess.STDOUT, timeout=600)
    if completed.returncode:
        raise RuntimeError('native source/build prerequisite failed; see prepare.log')
    hashes = json.loads((fixture/'prepared.json').read_text())['hashes']
    source, generated = REPO/'srcSEP/main.cpp', build/'build/main/main.cpp'
    normalize = lambda text: re.sub(r'\s+', '', re.sub(
        r'/\* AMPS_CHECKMACRO_BEGIN \*/.*?/\* AMPS_CHECKMACRO_END \*/', '', text, flags=re.S).replace('$PREFIX:', 'AMPS:'))
    if normalize(source.read_text()) != normalize(generated.read_text()):
        raise ValueError('native standalone main differs from checkout')
    executable, main_object = build/'amps', build/'build/main/main.o'
    for dependency in (generated, build/'build/main/makefile'):
        if main_object.stat().st_mtime_ns < dependency.stat().st_mtime_ns:
            raise ValueError('standalone main object is stale: ' + str(dependency))
    if any(executable.stat().st_mtime_ns < path.stat().st_mtime_ns
           for path in (main_object, build/'build/libAMPS.a')):
        raise ValueError('standalone executable predates its main/archive')
    # Freeze a private executable without changing the configured build or
    # historical test artifacts. Every source asset remains separately hashed.
    shutil.copy2(executable, root/'amps')
    paths = [source, generated, main_object, executable, root/'amps', Path(__file__),
             Path(audit.__file__), REPO/'srcSEP/test/test_standalone_application_verifier.py']
    paths += [path for path in (REPO/'srcSEP/examples/reduced-shock').rglob('*') if path.is_file()]
    hashes.update({str(path): digest(path) for path in paths})
    for job in jobs_for(root):
        Path(job['cwd']).mkdir()
    write_json(root/'build-sha256.json', hashes)
    write_json(root/'launch-plan.json', {'schema': 'srcsep-full-application-plan-v1',
        'scope': 'synthetic prescribed reduced SSE/Parker; fixed 300-s physical dt; no observational claim',
        'seed': 0, 'jobs': jobs_for(root)})


def execute(root):
    """Retain actual launcher exits and fail closed on timeout or startup error.

    Each mover runs independently in its owned directory. A failure does not
    discard the remaining mover evidence or overwrite a previous receipt.
    """
    plan = check_prepared(root)
    for job in plan['jobs']:
        directory = Path(job['cwd'])
        if (directory/'process-result.json').exists():
            raise ValueError('refusing to overwrite an existing process receipt: ' + job['mode'])
        started = time.monotonic()
        receipt = {'argv': job['argv'], 'cwd': job['cwd']}
        try:
            with (directory/'launcher.log').open('w') as stream:
                result = subprocess.run(job['argv'], cwd=directory, stdout=stream,
                                        stderr=subprocess.STDOUT, timeout=1830)
            receipt['exit_code'] = result.returncode
        except (OSError, subprocess.TimeoutExpired) as error:
            receipt.update(exit_code=None, error=str(error))
        receipt['wall_time_seconds'] = time.monotonic()-started
        write_json(directory/'process-result.json', receipt)


def collect(root):
    """Bind each completed receipt to its planned actual launcher invocation."""
    for job in check_prepared(root)['jobs']:
        receipt = json.loads((Path(job['cwd'])/'process-result.json').read_text())
        if receipt.get('argv') != job['argv'] or receipt.get('cwd') != job['cwd']:
            raise ValueError('process receipt belongs to another application invocation: ' + job['mode'])
    return subprocess.call([sys.executable, str(Path(audit.__file__)), '--campaign-dir', str(root)], cwd=REPO)


def report_error(root, message):
    """Missing build/MPI/output is ERROR for every required mover, never PASS."""
    import xml.etree.ElementTree as ET
    rows = [{'test_id': 'SAPP%02d' % i, 'mode': mode, 'status': 'ERROR', 'message': message}
            for i, mode in enumerate(audit.MODES, 1)]
    write_json(root/'verified-results.json', {'schema': 'srcsep-stage6-application-v1', 'results': rows})
    suite = ET.Element('testsuite', name='standalone-full-application', tests='3', errors='3')
    for row in rows:
        case = ET.SubElement(suite, 'testcase', name=row['test_id'])
        ET.SubElement(case, 'error', message=message)
    temporary = root/'junit.xml.tmp'
    ET.ElementTree(suite).write(temporary, encoding='unicode')
    temporary.replace(root/'junit.xml')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-root', type=Path, default=os.environ.get('AMPS_FIELD_LINE_NATIVE_BUILD'))
    parser.add_argument('--output-dir', type=Path, required=True)
    action = parser.add_mutually_exclusive_group()
    action.add_argument('--prepare-only', action='store_true')
    action.add_argument('--collect', action='store_true')
    args = parser.parse_args()
    root = args.output_dir.resolve()
    # Preparation must never overwrite a prior qualification or its failure
    # evidence. In particular, catching mkdir's FileExistsError below and then
    # publishing ERROR in that SAME directory would destroy its old verdict.
    # Collection is the explicit operation permitted to re-evaluate a plan.
    if root.exists() and not args.collect:
        print('ERROR: choose a fresh application evidence directory: '+str(root), file=sys.stderr)
        return 1
    try:
        if not args.collect:
            if args.build_root is None:
                raise ValueError('a matching configured --build-root is required')
            prepare(root, Path(args.build_root).resolve())
        check_prepared(root)
        if args.prepare_only:
            print('PREPARED, not executed: ' + str(root/'launch-plan.json'), flush=True)
            return 0
        if not args.collect:
            execute(root)
        # The original independent checker owns numerical acceptance. Its own
        # process emits every SAPP result and JUnit after a native failure.
        return collect(root)
    except (OSError, ValueError, KeyError, RuntimeError, subprocess.TimeoutExpired) as error:
        root.mkdir(parents=True, exist_ok=True)
        write_json(root/'preflight-error.json', {'status': 'ERROR', 'message': str(error)})
        report_error(root, str(error))
        print('ERROR: ' + str(error), file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
