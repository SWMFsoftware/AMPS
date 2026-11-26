#!/usr/bin/env python3
"""Run and score a frozen family of the installed Moon chemistry development cases."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

from campaign_budget import CampaignBudget
import m3_qualification as q

TEST = Path(__file__).resolve().parent
ROOT = TEST.parents[1]
SCHEMA = 'srcMoon-M3-development-family-v1'
EXIT = {'PASS': 0, 'FAIL': 1, 'ERROR': 2, 'SKIPPED': 77}
ORDER = ('ERROR', 'FAIL', 'SKIPPED', 'PASS')
KINDS = {
    'i21-cohort': ('I21', ('plan',)),
    'i24-sodium-ensemble': ('I24', ('lock', 'mass_inventory')),
    'i24-sodium-pic': ('I24', ('lock', 'mass_inventory', 'producer')),
    # Approved 2026-11-25 controlled electron bath; the case JSON and frozen
    # admission files are pinned inside the runner, the profile/producer here.
    'i23-helium-ensemble': ('I23', ('profile',)),
    'i23-helium-pic': ('I23', ('profile', 'producer')),
}
RUNNERS = {
    'i21-cohort': TEST/'m3_qualification_executor.py',
    'i24-sodium-ensemble': TEST/'integration-tests/I24/run_ground_capture.py',
    'i24-sodium-pic': TEST/'integration-tests/I24/run_pic_ground_capture.py',
    'i23-helium-ensemble': TEST/'integration-tests/I23/run_helium_impact.py',
    'i23-helium-pic': TEST/'integration-tests/I23/run_pic_helium_impact.py',
}
# The underlying native runners separately bind their compilation and data
# dependencies. Retain the Python entry points and independent scorers as
# well, since their working-tree bytes can change without a new Git revision.
FAMILY_TOOLS = (
    Path(__file__).resolve(), *RUNNERS.values(), TEST/'integration_m3.py',
    TEST/'campaign_budget.py', TEST/'m3_qualification.py',
    TEST/'sodium_cohort/audit.py', TEST/'sodium_cohort/verify_numerical.py',
    TEST/'sodium_cohort/verify_ledger.py',
    TEST/'integration-tests/I23/prepare_helium_admission.py',
    TEST/'integration-tests/I23/reference/helium_case.json',
    TEST/'integration-tests/I23/reference/helium_admission_approval_20261125.json',
)


def write_json(path, value):
    # An interrupted later case leaves the preceding checkpoint readable.
    pending = path.with_suffix('.pending')
    pending.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    pending.replace(path)


def descriptor(path):
    return {'path': str(path), 'sha256': q.sha(path)}


def pin(record, base):
    if (not isinstance(record, dict) or set(record) != {'path', 'sha256'}
            or not isinstance(record['path'], str) or not record['path']
            or not q.digest_syntax(record['sha256'])):
        raise ValueError('each family input requires an explicit path and lowercase SHA-256')
    return {'path': str((base / record['path']).resolve()), 'sha256': record['sha256']}


def authenticate(inputs):
    """Check all available identities even when another file is missing."""
    missing = []
    for record in inputs.values():
        path = Path(record['path'])
        if not path.is_file():
            missing.append(str(path))
        elif q.sha(path) != record['sha256']:
            raise ValueError(f'pinned family input changed: {path}')
    if missing:
        raise q.MissingPrerequisite('pinned input unavailable: ' + ', '.join(missing))


def read_family(path, expected):
    """Freeze every case before inspecting outcomes or invoking any runner."""
    authenticate({'family': pin({'path': str(path), 'sha256': expected}, path.parent)})
    value = q.strict_object(path)
    q.approved_protocol()
    if (set(value) != {'schema', 'approved_addendum_sha256', 'cases'}
            or value['schema'] != SCHEMA or value['approved_addendum_sha256'] != q.APPROVED_SHA256):
        raise ValueError('family plan schema or approved protocol mismatch')
    if not isinstance(value['cases'], list) or not value['cases']:
        raise ValueError('family plan needs an explicit, nonempty ordered case list')
    names, selections, cases = set(), set(), []
    for case in value['cases']:
        if not isinstance(case, dict) or set(case) != {'name', 'kind', 'inputs'}:
            raise ValueError('each case needs name, kind and inputs')
        name, kind, inputs = case['name'], case['kind'], case['inputs']
        if not isinstance(name, str) or not re.fullmatch('[a-z0-9][a-z0-9-]{0,79}', name) or name in names:
            raise ValueError('use unique lowercase alphanumeric/hyphen case names (1 to 80 characters)')
        if not isinstance(kind, str) or kind not in KINDS:
            raise ValueError(f'unsupported development case kind: {kind}')
        identity, roles = KINDS[kind]
        if not isinstance(inputs, dict) or set(inputs) != set(roles):
            raise ValueError(f'{name} requires exactly these pinned input roles: {", ".join(roles)}')
        pins = {role: pin(inputs[role], path.parent) for role in roles}
        selection = (kind, tuple(pins[role]['sha256'] for role in roles))
        if selection in selections:
            raise ValueError(f'duplicate physical matrix selection: {name}')
        names.add(name)
        selections.add(selection)
        cases.append({'name': name, 'kind': kind, 'test_id': identity, 'inputs': pins})
    return cases


def case_command(case, target, predecessors, allowance):
    inputs = case['inputs']
    common = ['--output-dir', str(target)]
    runner = str(RUNNERS[case['kind']])
    if case['kind'] == 'i21-cohort':
        plan = inputs['plan']
        return [runner, '--case-plan', plan['path'],
                '--expected-plan-sha256', plan['sha256'], '--campaign-root', str(predecessors),
                '--development', '--wall-seconds', str(allowance), *common]
    if case['kind'] == 'i23-helium-ensemble':
        return [runner, '--test-runtime-seconds', str(allowance), '--mpi-ranks', '8',
                '--profile', inputs['profile']['path'], *common]
    if case['kind'] == 'i23-helium-pic':
        return [runner, '--producer-manifest', inputs['producer']['path'], '--refinement',
                '--wall-seconds', str(allowance), '--profile', inputs['profile']['path'], *common]
    source = ['--lock', inputs['lock']['path'], '--mass-inventory', inputs['mass_inventory']['path']]
    if case['kind'] == 'i24-sodium-ensemble':
        return [runner, '--test-runtime-seconds', str(allowance),
                '--mpi-ranks', '8', *source, *common]
    return [runner, '--producer-manifest', inputs['producer']['path'],
            '--refinement', '--wall-seconds', str(allowance), *source, *common]


def checks_from_case(kind, directory, result):
    """Index maintained independent checks, keeping the underlying artifacts."""
    if kind == 'i21-cohort' and result.get('execution') == 'actual-native-cohort-and-independent-audit':
        audit = q.strict_object(directory/'independent-audit/result.json')
        if (audit.get('status'), audit.get('checks')) != (result['status'], result.get('checks')):
            raise ValueError('I21 owner and independent audit disagree')
        checks = audit['checks']
    elif kind in ('i24-sodium-pic', 'i23-helium-pic'):
        layouts = result.get('checks')
        if not isinstance(layouts, dict) or any(
                not isinstance(values, dict) or any(type(passed) is not bool for passed in values.values())
                for values in layouts.values()):
            raise ValueError('invalid PIC layout or check status')
        checks = [{'name': f'{layout}/{name}', 'status': 'PASS' if passed else 'FAIL'}
                  for layout, values in layouts.items() for name, passed in values.items()]
    else:
        checks = result.get('checks', [])
    if not isinstance(checks, list) or any(not isinstance(item, dict) or
            not item.get('name') or item.get('status') not in EXIT for item in checks):
        raise ValueError('independent check list is invalid')
    return checks


def counts(rows):
    return {status: sum(row['status'] == status for row in rows) for status in EXIT}


def supervise(summary, command, directory, log, resources, budget, root, *, mpi_ranks=None):
    argv = [sys.executable, '-B', *command]
    summary.update(command=argv, runner=command[0], log=str(log))
    receipt = budget.run(argv, ROOT, log, output=resources, artifact_root=root,
                         allocated_cores=mpi_ranks or 1,
                         environment=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1'))
    summary['resources'] = receipt
    if receipt['reason']:
        raise RuntimeError(f"{receipt['reason']}. {receipt['remediation']}")
    path = directory/'result.json'
    result = q.strict_object(path)
    summary['result'] = descriptor(path)
    summary['run_manifest'] = descriptor(directory/'run_manifest.json')
    if EXIT.get(result.get('status')) != receipt['returncode'] or result.get('status') not in EXIT:
        raise ValueError('child result and exit disagree; see ' + str(log))
    summary.update(status=result['status'], reason=result.get('reason'))
    return result


def run(plan, expected, output, wall_seconds=3600, predecessors=None, amps=None):
    out = Path(output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    plan = Path(plan).resolve()
    predecessors = Path(predecessors).resolve() if predecessors is not None else out
    amps = Path(amps).resolve() if amps is not None else ROOT/'amps'
    (out/'logs').mkdir()
    result = {'schema': SCHEMA, 'status': 'ERROR', 'scope': 'linked M3 development verification',
              'cases': [], 'formal': [], 'reason': None}
    manifest = {'schema': SCHEMA, 'created_utc': datetime.now(timezone.utc).isoformat(),
                'source_revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                'wall_seconds': str(wall_seconds), 'mpi_ranks': [1, 8], 'openmp_threads': 1,
                'observational_scoring': False, 'predecessor_root': str(predecessors),
                'runner': descriptor(Path(__file__).resolve()), 'plan': str(plan), 'expected_plan_sha256': expected}

    def guarded(row, work):
        try:
            work()
        except q.MissingPrerequisite as exc:
            row.update(status='SKIPPED', reason=str(exc))
        except (OSError, ValueError, TypeError, KeyError, ArithmeticError, RuntimeError) as exc:
            row.update(status='ERROR', reason=str(exc))

    try:
        if not 0 < wall_seconds <= 3600:
            raise ValueError('development campaign wall allowance must be in (0, 3600] seconds')
        budget = CampaignBudget(wall_seconds, min(120, wall_seconds/10),
                                memory_bytes=q.approved_protocol()['execution']['aggregate_memory_bytes'])
        cases = read_family(plan, expected)
        frozen_plan = {'path': str(plan), 'sha256': expected}
        snapshot = out/'family_plan.json'
        snapshot.write_bytes(plan.read_bytes())
        if q.sha(snapshot) != expected:
            raise ValueError('family plan changed while taking its input snapshot')
        manifest['plan_snapshot'] = descriptor(snapshot)
        manifest['case_inputs'] = cases
        tool_pins = {path.name: descriptor(path) for path in FAMILY_TOOLS}
        manifest['tools'] = list(tool_pins.values())
        write_json(out/'run_manifest.json', manifest)

        def formal_case(row):
            identity = row['test_id']
            cmd = [str(TEST/'integration_m3.py'), identity, '--campaign-mode', 'qualification',
                   '--amps', str(amps), '--input', str(ROOT/'srcMoon/examples/m3_reaction_rates.in'),
                   '--predecessor-root', str(predecessors), '--output-dir', str(out/'formal')]
            supervise(row, cmd, out/'formal'/identity, out/'logs'/f'formal-{identity}.log',
                      out/'resources'/f'formal-{identity}', budget, out)

        def development_case(case, row):
            authenticate({'family': frozen_plan, **case['inputs']})
            directory = out/'cases'/case['name']
            log = out/'logs'/f"{case['name']}.log"
            command = case_command(case, directory, predecessors, budget.remaining())
            value = supervise(row, command, directory, log, out/'resources'/case['name'],
                              budget, out, mpi_ranks=8)
            authenticate({'family': frozen_plan, **case['inputs']})
            row['checks'] = checks_from_case(case['kind'], directory, value)
            if row['status'] == 'PASS' and (not row['checks'] or
                    any(check['status'] != 'PASS' for check in row['checks'])):
                raise ValueError('supporting PASS requires actual passing independent checks')

        # Validate pins for the entire case family before launching any case.
        # An unavailable member has its own status and does not erase others.
        for case in cases:
            row = {'name': case['name'], 'kind': case['kind'], 'test_id': case['test_id'],
                   'checks': [], 'status': None, 'runner': str(RUNNERS[case['kind']]),
                   'log': str(out/'logs'/f"{case['name']}.log")}
            guarded(row, lambda: authenticate(case['inputs']))
            if row['status'] is not None:
                Path(row['log']).write_text(f"{row['status']}: {row['reason']}\n")
            result['cases'].append(row)

        for identity in q.ENSEMBLE_LIMIT:
            row = {'test_id': identity}
            guarded(row, lambda: formal_case(row))
            result['formal'].append(row)

        write_json(out/'progress.json', dict(result, status='RUNNING'))
        for case, row in zip(cases, result['cases']):
            if row['status'] is None:
                print(f"RUNNING {row['test_id']} {row['name']}  log={row['log']}", flush=True)
                guarded(row, lambda: development_case(case, row))
                if not Path(row['log']).is_file():
                    Path(row['log']).write_text(f"{row['status']}: {row['reason']}\n")
            print(f"{row['status']} {row['test_id']} {row['name']}", flush=True)
            write_json(out/'progress.json', dict(result, status='RUNNING'))
        authenticate({'family': frozen_plan, **tool_pins})
        # Formal SKIPs describe incomplete chemistry and predecessors. An
        # infrastructure ERROR in the readiness assessment is actionable.
        result['status'] = next(s for s in ORDER if any(row['status'] == s for row in result['cases'])
                               or (s == 'ERROR' and any(row['status'] == s for row in result['formal'])))
        result['reason'] = 'see per-case independent checks and formal predecessor records'
    except q.MissingPrerequisite as exc:
        result.update(status='SKIPPED', reason=str(exc))
    except (OSError, ValueError, TypeError, KeyError, ArithmeticError, RuntimeError) as exc:
        result.update(status='ERROR', reason=str(exc))
    result['counts'] = counts(result['cases'])
    result['formal_counts'] = counts(result['formal'])
    result['independent_check_counts'] = counts([check for row in result['cases'] for check in row['checks']])
    result['wall_seconds'] = time.monotonic()-started
    with (out/'comparison.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(('test_id', 'case', 'check', 'status', 'reason', 'result_path'))
        for row in result['cases']:
            path = row.get('result', {}).get('path', '')
            writer.writerow((row['test_id'], row['name'], 'case outcome', row['status'], row.get('reason'), path))
            for check in row['checks']:
                writer.writerow((row['test_id'], row['name'], check['name'], check['status'], check.get('reason'), path))
    manifest['source_revision_at_end'] = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    if manifest['source_revision_at_end'] != manifest['source_revision']:
        result.update(status='ERROR', reason='AMPS revision changed during the family campaign')
    manifest['comparison'] = descriptor(out/'comparison.csv')
    write_json(out/'run_manifest.json', manifest)
    write_json(out/'result.json', result)
    # A monitoring client can follow this file until it reaches a final
    # status, including when preflight or the last child reports an error.
    write_json(out/'progress.json', result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--plan-sha256', required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--wall-seconds', type=float, default=3600)
    parser.add_argument('--predecessor-root', type=Path, help='existing output from the same qualification campaign')
    parser.add_argument('--amps', type=Path, help='executable identity for formal readiness records')
    args = parser.parse_args()
    result = run(args.plan, args.plan_sha256, args.output_dir, args.wall_seconds, args.predecessor_root, args.amps)
    print(json.dumps({key: result[key] for key in ('status', 'reason', 'counts', 'formal_counts',
                                                  'independent_check_counts', 'wall_seconds')}, indent=2))
    raise SystemExit(EXIT[result['status']])
