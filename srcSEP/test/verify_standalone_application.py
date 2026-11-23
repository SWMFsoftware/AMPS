#!/usr/bin/env python3
"""Audit the frozen stage-6 full-driver campaign, including a complete output.

This reads actual process receipts and production output; it never manufactures
a reference solution, reruns a kernel to grade itself, or changes acceptance
limits. The verification case fixes 400 physical steps of 300 s. A short run,
missing receipt, crash, changed executable/input or absent output cannot pass.
Launch commands and the compiled source/build hashes remain campaign artifacts.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import shlex
import xml.etree.ElementTree as ET

MODES = ('parker', 'fte-dmumu', 'fte-mfp')
STEPS, DT, TOLERANCE = 400, 300.0, 1e-10


class GateFailure(RuntimeError):
    """An evaluated numerical/software acceptance condition was not met."""


def require(condition, message):
    if not condition:
        raise GateFailure(message)


def rows(text, prefix):
    """Parse numerical production records without inventing missing fields."""
    return [dict((key, float(value)) for key, value in
                 (token.split('=', 1) for token in line.split()[1:]))
            for line in text.splitlines() if line.startswith(prefix + ' ')]


def balance(terms, scale, description):
    # fsum is independent of the production reduction order. Keep the SAME
    # relative conservation bound, including exact zero when the scale is zero;
    # a numerical floor would hide a budget error in an empty population.
    residual = math.fsum(terms)
    require(math.isfinite(residual) and abs(residual) <= TOLERANCE * scale,
            description + ' does not close')
    return abs(residual) / scale if scale else 0.0


def verify_mode(root, mode, job, hashes):
    directory = root / mode
    receipt = json.loads((directory / 'process-result.json').read_text())
    if type(receipt.get('exit_code')) is not int:
        raise ValueError('missing completed process return code')
    if receipt['exit_code']:
        raise RuntimeError('application exited with code ' + str(receipt['exit_code']))
    # The shell is an external MPI-launch compatibility mechanism. Recover the
    # actual executable/options from its recorded command, then verify the
    # immutable executable and physical inputs rather than a later live build.
    outer = shlex.split(job['command'])
    require(outer[:3] == ['mpiexec', '-n', '1'], 'stage 6 requires actual one-rank MPI')
    child = shlex.split(outer[-1])
    require(child[:3] == ['exec', 'timeout', '1800s'], 'campaign deadline changed')
    argv = child[3:child.index('>')]
    executable = Path(argv[0])
    expected = {'--total-iterations': '400', '--particle-mover': mode,
                '--reduced-shock-particles': 'on', '--shock-release-distance-m': '1e8',
                '--shock-return-policy': 'delete', '--shock-source-normalization': 'net-upstream-escape',
                '--turbulence-source': 'prescribed', '--coefficient-source': 'prescribed',
                '--constant-dmumu': '0.001', '--turbulence-coupling-policy': 'disabled',
                '--coupling': 'off', '--cascade': 'off', '--reflection': 'off',
                '--turbulence-advection': 'off', '--shock-injection': 'off',
                '--pitch-angle-diffusion-provider': 'constant',
                '--spatial-diffusion-provider': 'from-dmumu', '--mean-free-path-provider': 'from-spatial'}
    for flag, value in expected.items():
        require(argv[argv.index(flag) + 1] == value, 'verification assignment changed: ' + flag)
    for path in [executable, Path(argv[argv.index('--input') + 1]),
                 Path(argv[argv.index('--reduced-shock-event') + 1])]:
        require(hashlib.sha256(path.read_bytes()).hexdigest() == hashes[str(path)],
                'executable or physical input changed: ' + str(path))
    # This is also a CODE-validation gate. An unchanged old executable/output
    # cannot certify an edited checkout or a subsequently rebuilt archive.
    # Verify every recorded source/header/library identity, not just inputs.
    # A new build therefore requires fresh matching campaign evidence.
    for name, digest in hashes.items():
        require(hashlib.sha256(Path(name).read_bytes()).hexdigest() == digest,
                'recorded source or build input changed: ' + name)

    text = (directory / 'application.rank0.log').read_text()
    require('AMPS: exit:' not in text, 'native host error appears in output')
    steps = re.findall(r'Iteration: (\d+)\s+\(current sample length:(\d+), (\d+) iterations', text)
    require(steps == [(str(i), '400', str(400-i)) for i in range(STEPS)],
            'physical step sequence or 400-step sampling cadence changed')
    ledgers = rows(text, 'REDUCED_SOURCE_LEDGER')
    geometry = rows(text, 'REDUCED_GEOMETRY')
    require(len(ledgers) == STEPS and len(geometry) == STEPS, 'incomplete physical campaign')
    maximum = 0.0
    for i, ledger in enumerate(ledgers):
        require(all(math.isfinite(v) for v in ledger.values()), 'nonfinite SI ledger')
        require(ledger['epoch_s'] == (i+1)*DT, 'physical endpoint changed')
        n, e = ledger['born'], ledger['born_kinetic_j']
        maximum = max(maximum, balance([n, -ledger['retained'], -ledger['returned'],
            -ledger['domain_exit'], -ledger['geometry_retired']], n, 'represented number'))
        maximum = max(maximum, balance([e, -ledger['retained_birth_j'], -ledger['returned_birth_j'],
            -ledger['domain_birth_j'], -ledger['retired_birth_j']], e, 'immutable birth energy'))
        maximum = max(maximum, balance([e, ledger['background_work_j'], ledger['wave_work_j'],
            -ledger['retained_kinetic_j'], -ledger['returned_event_j'], -ledger['domain_event_j'],
            -ledger['retired_event_j']], abs(e)+abs(ledger['background_work_j'])+abs(ledger['wave_work_j']),
            'current kinetic energy and signed external work'))
    require(ledgers[-1]['born'] > 0 and ledgers[-1]['returned'] > 0,
            'production source and absorbing returns were not exercised')

    sources = []
    for line in text.splitlines():
        if line.startswith('REDUCED_SOURCE event='):
            fields = dict(token.split('=', 1) for token in line.split()[1:])
            sources.append((int(fields['event']), float(fields['candidate']),
                            float(fields['no_focused_escape']), float(fields['committed'])))
    require([s[0] for s in sources] == list(range(STEPS)), 'source event history is incomplete')
    require(all(s[2] == 0 and s[1] == s[3] for s in sources), 'net source was filtered or rescaled')
    for i, ledger in enumerate(ledgers):
        balance([math.fsum(s[3] for s in sources[:i+1]), -ledger['born']], ledger['born'],
                'committed source and injected inventory')

    # Reconstruct each surviving segment's birth step from actual topology
    # commits, independently of its exported exposure denominator. Output is
    # sampled at the DEPARTURE of physical step 400, before geometry commit 400.
    initial = int(re.search(r'REDUCED_BACKGROUND epoch_s=0(?:\.0+e\+00)? .*?vertices=(\d+)', text).group(1))-1
    first, count = 0, initial
    births = dict.fromkeys(range(initial), 0)
    for i, event in enumerate(geometry):
        require(event['epoch_s'] == (i+1)*DT and event['generation'] == i+1,
                'material epoch/generation sequence changed')
        inserted, removed = int(event['inserted']), int(event['removed'])
        next_count = count+inserted-removed
        require(event['segments'] == next_count, 'material topology count does not close')
        if i == STEPS-1:
            break
        for label in range(first-inserted, first):
            births[label] = i+1
        first -= inserted
        count = next_count
    output = directory / 'PT/plots/amps.FieldLines.out=0.dat'
    exposure = [line.split() for line in output.with_name(output.name+'.sampling').read_text().splitlines()
                if line and not line.startswith('#')]
    require(len(exposure) == count, 'sampling sidecar has missing/duplicate segments')
    ages = []
    for index, fields in enumerate(exposure):
        line, label, epoch, generation, samples, seconds, fixed_dt, valid = map(float, fields)
        expected_samples = STEPS-births[first+index]
        require((line, label, epoch, generation) == (0, first+index, (STEPS-1)*DT, STEPS-1),
                'sampling sidecar has stale topology/epoch')
        require((samples, seconds, fixed_dt, valid) == (expected_samples, expected_samples*DT, DT, 1),
                'segment lifetime denominator or physical sampling dt changed')
        ages.append(samples)
    require(max(ages) == STEPS and min(ages) < STEPS, 'old/new material sampling histories absent')

    data = output.read_text()
    require('shock model=reduced-SSE/Parker;' in data, 'output names an inactive shock authority')
    require(float(re.search(r'TITLE="time=([^;]+); shock model=', data).group(1)) == (STEPS-1)*DT,
            'output timestamp is not its actual sampling departure epoch')
    names = re.findall(r'"([^"]*)"', next(line for line in data.splitlines() if line.startswith('VARIABLES=')))
    require(len(names) == len(set(names)) and ',,' not in data, 'malformed/duplicated output columns')
    numerical = [list(map(float, line.split())) for line in data.splitlines()
                 if line.strip() and line.lstrip()[0] in '-+.0123456789']
    require(len(numerical) == count+1 and all(len(row) == len(names) for row in numerical),
            'native output vertex/column count disagrees with material topology')
    primitives = ['Bx [nT]', 'By [nT]', 'Bz [nT]', 'Plasma Vx [m/s]', 'Plasma Vy [m/s]',
                  'Plasma Vz [m/s]', 'Plasma number density [1/m^3]', 'Plasma Temperature [K]', 'Plasma pressure [Pa]']
    for row in numerical:
        require(all(math.isfinite(row[names.index(name)]) for name in primitives), 'nonfinite active background primitive')
        require(math.hypot(*(row[names.index(name)] for name in primitives[:3])) > 0,
                'active magnetic field has an unresolved null')
        require(all(row[names.index(name)] > 0 for name in primitives[-3:]), 'nonpositive active plasma state')
    return {'physical_steps': STEPS, 'physical_dt_s': DT, 'final_epoch_s': ledgers[-1]['epoch_s'],
            'maximum_relative_budget_residual': maximum, 'sampling_segments': count,
            'minimum_segment_samples': min(ages), 'maximum_segment_samples': max(ages),
            'final_ledger': ledgers[-1], 'source_history': sources,
            'log': str(directory/'application.rank0.log'), 'process_receipt': str(directory/'process-result.json')}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign-dir', type=Path, required=True)
    args = parser.parse_args()
    root = args.campaign_dir.resolve()
    root.mkdir(parents=True, exist_ok=True)
    preflight_error = None
    try:
        plan = json.loads((root/'launch-plan.json').read_text())
        jobs = {job['mode']: job for job in plan['jobs']}
        hashes = json.loads((root/'build-sha256.json').read_text())
    except (OSError, KeyError, ValueError) as error:
        preflight_error = str(error)
    results = []
    for index, mode in enumerate(MODES, 1):
        row = {'test_id': 'SAPP%02d' % index, 'mode': mode, 'status': 'PASS',
               'seed': 0, 'acceptance_relative_budget_bound': TOLERANCE}
        try:
            if preflight_error is not None:
                raise ValueError(preflight_error)
            row.update(verify_mode(root, mode, jobs[mode], hashes))
        except GateFailure as error:
            row.update(status='FAIL', message=str(error))
        except (OSError, KeyError, ValueError, AttributeError, RuntimeError) as error:
            row.update(status='ERROR', message=str(error))
        results.append(row)
    if all(row['status'] == 'PASS' for row in results):
        if any(row['source_history'] != results[0]['source_history'] for row in results):
            for row in results:
                row.update(status='FAIL', message='mover selection changed the physical source history')
    for row in results:
        print('[%s] %s %s %s' % (row['test_id'], row['status'], row['mode'], row.get('message', 'complete fixed-step native campaign')))
    (root/'verified-results.json').write_text(json.dumps({'schema': 'srcsep-stage6-application-v1',
        'scope': 'synthetic prescribed reduced SSE/Parker full driver; not observational or moving-restart qualification',
        'results': results}, indent=2)+'\n')
    suite = ET.Element('testsuite', name='stage6-full-application', tests=str(len(results)))
    for row in results:
        case = ET.SubElement(suite, 'testcase', name=row['test_id'])
        if row['status'] != 'PASS':
            ET.SubElement(case, 'failure' if row['status'] == 'FAIL' else 'error', message=row['message'])
    ET.ElementTree(suite).write(root/'junit.xml', encoding='unicode')
    return int(any(row['status'] != 'PASS' for row in results))


if __name__ == '__main__':
    raise SystemExit(main())
