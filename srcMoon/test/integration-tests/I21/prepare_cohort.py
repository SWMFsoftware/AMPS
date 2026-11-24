#!/usr/bin/env python3
"""Freeze an I21 cohort from explicit physical inputs and a reference receipt."""
from __future__ import annotations

import argparse
from decimal import Decimal
import json
from pathlib import Path
import re
import shlex
import sys

TEST_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(TEST_ROOT))
import m3_qualification as qualification
from sodium_cohort.prepare_physical_case import (
    MissingConfiguration, prepare, validate_controls,
)
from sodium_cohort.verify_numerical import reference_receipt, source_rate


def descriptor(path):
    return {'path': str(path.resolve()), 'sha256': qualification.sha(path)}


def write_json(path, value):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def spice_sources(producer_manifest):
    """Bind the ordered kernel recipe and the producer's configured kernel root."""
    producer=qualification.strict_object(producer_manifest)
    configuration=shlex.split(producer['build_commands'][0])
    directories=[arg.partition('=')[2] for arg in configuration if arg.startswith('-spice-kernels=')]
    if len(directories)!=1 or not Path(directories[0]).is_absolute():
        raise ValueError('the source producer must bind one absolute SPICE kernel root')
    derived=producer['source_qualification_recipe']['derived_inputs']
    recipes=[Path(path) for path in derived if path.endswith('.input.spicekernels')]
    if len(recipes)!=1 or qualification.sha(recipes[0])!=derived[str(recipes[0])]:
        raise ValueError('the producer kernel recipe is missing or changed')
    contract_path=TEST_ROOT/'U04_rotating_frame/reference/acceptance.json'
    contract=qualification.strict_object(contract_path)
    names=re.findall(r'"([^"\n]+)"',recipes[0].read_text().split('SpiceKernels=',1)[1])
    if names!=contract['kernel_files']:
        raise ValueError('the cohort requires the reviewed ordered SPICE kernel set')
    checked=[descriptor(recipes[0]), descriptor(contract_path)]
    for name in names:
        record=descriptor(Path(directories[0])/name)
        if record['sha256']!=contract['data_hashes'][name]:
            raise ValueError('SPICE kernel differs from its reviewed identity: '+name)
        checked.append(record)
    return checked


def prepare_cohort(args):
    """Bind one cohort, varying only physical step count and MPI layout.

    The reference input receipt supplies independently frozen geometry and
    visibility. It is used to calculate an expected rate from D04 before the
    launches. The auditor independently recomputes that rate with the new
    receipts and requires the numerical geometry to match across refinements.
    No events, clocks or survival counts from the reference are read here.
    """
    if not qualification.digest_syntax(args.reference_sha256):
        raise ValueError('explicit lowercase reference receipt SHA-256 required')
    if qualification.sha(args.reference) != args.reference_sha256:
        raise ValueError('reference receipt SHA-256 mismatch')
    if args.branch_claims < 6:
        raise ValueError('branch family must include at least the six declared launches')
    case = qualification.strict_object(args.controls)
    _, controls, evidence = validate_controls(args.controls)
    reference = reference_receipt(qualification.strict_object(args.reference))
    varying = controls['case_kind'] == 'trajectory-integrated-hazard'
    if varying:
        # The numerical reference is an illuminated x-axis straight line.
        # Authenticated native receipts will additionally prove zero force,
        # the same velocity and finite-ray visibility over the whole horizon.
        velocity = [Decimal(controls[f'source_velocity_lso_ms_{axis}']) for axis in range(3)]
        if velocity[0] == 0 or any(velocity[1:]) or reference['point_sun_visible'] != 1:
            raise ValueError('varying exposure requires x-only nonzero velocity and an illuminated initial point')
    position = [Decimal(controls['source_position_lso_m_'+str(axis)]) for axis in range(3)]
    if (position != [Decimal(value) for value in reference['source_position_lso_m']]
            or controls['solar_endpoint'] != reference['solar_endpoint']
            or Decimal(controls['source_reference_distance_m']) != Decimal(reference['source_reference_distance_m'])
            or controls['source_reference_distance_identity'] != reference['source_reference_distance_identity']):
        raise ValueError('reference source position, endpoint or distance convention differs from case controls')
    data = {}
    rate, channel = source_rate(reference, data)
    kernels = spice_sources(args.producer)
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    # Controls copied into the run directory must retain the original resolved
    # policy evidence path. A relative path would be relative to the new case.
    case['controls']['transport_policy_evidence_file'] = evidence['path']
    plan = {'schema': 'srcMoon-M3-cohort-plan-v1', 'test_id': 'I21',
        'approved_addendum_sha256': qualification.APPROVED_SHA256,
        'cohort_id': case['cohort_id'], 'seed': int(controls['campaign_seed']),
        'original_parents': int(controls['original_parents']),
        'initial_represented_number': float(controls['represented_number_per_parent']),
        'horizon_s': float(controls['physical_horizon_s']),
        'test_kind': 'trajectory-integrated-hazard' if varying else
            'zero-hazard-invariant' if rate == 0 else 'homogeneous-fixed-position',
        'reference_total_rate_s': float(rate),
        'primary_estimator': 'integrated-exposure' if varying else 'zero-hazard' if rate == 0 else 'lifetime',
        'branch_claims': args.branch_claims,
        'source_channel_rates_s': [{'id': channel, 'rate_s': float(rate)}],
        'source_profile': descriptor(args.profile),
        'reference_binding': {'receipt': descriptor(args.reference), 'data_sha256': data,
            'spice_files': kernels},
        'native_adapter': {'schema': 'srcMoon-M3-normal-sodium-execution-v1',
            'producer_id': 'normal-source-bound-sodium-cohort-v1',
            'producer_manifest': descriptor(args.producer),
            'state_mapping': {}, 'runs': []}}
    if varying:
        plan['trajectory_reference'] = {
            'oracle': 'illuminated-zero-force-lso-x-axis-v1',
            'source_velocity_lso_ms': [str(value) for value in velocity]}
    horizon = Decimal(controls['physical_horizon_s'])
    for steps in (8, 16, 32):
        controls_path = output / f'controls-{steps}.json'
        case['controls'].update(physical_steps=str(steps), physical_step_s=str(horizon / steps))
        write_json(controls_path, case)
        prepared = output / f'steps-{steps}'
        prepare(args.profile, args.crosswalk, args.application_input, controls_path, prepared)
        preparation = qualification.strict_object(prepared / 'run_manifest.json')
        states = qualification.strict_object(prepared / 'case.json')['source_crosswalk']
        plan['native_adapter']['state_mapping'] = {
            key: states[key+'_state_id'] for key in ('parent','ion','electron')}
        root = descriptor(prepared / 'native.in')
        for ranks in (1, 8):
            plan['native_adapter']['runs'].append({'mpi_ranks': ranks, 'physical_steps': steps,
                'input': root, 'input_documents': [root, *preparation['application_documents']]})
    qualification.validate_plan(plan)
    plan_path = output / 'cohort_plan.json'
    write_json(plan_path, plan)
    return descriptor(plan_path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for role in ('profile','crosswalk','application-input','controls','reference','producer','output-dir'):
        parser.add_argument('--'+role,type=Path,required=True)
    parser.add_argument('--reference-sha256',required=True)
    parser.add_argument('--branch-claims',type=int,required=True)
    args = parser.parse_args()
    try:
        print(json.dumps(prepare_cohort(args), indent=2))
    except (FileNotFoundError, MissingConfiguration, qualification.MissingPrerequisite) as error:
        print(f'SKIPPED {error}', file=sys.stderr)
        return 77
    except (OSError, ValueError, ArithmeticError, KeyError, TypeError) as error:
        print(f'ERROR {error}', file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
