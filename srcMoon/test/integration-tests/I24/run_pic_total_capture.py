#!/usr/bin/env python3
"""Run one actual I24 total-capture transaction with each mover at 1 and 8 ranks.

Approved 2026-11-25 (TOTAL_CAPTURE_CONTRACT.md). The production owner
MoonSodiumTotalCaptureCohort is registered through normal startup. The audit
reconstructs the event time from the keyed source clock, the selected H(n)
shell from the keyed Kinematics draw, the zero-angle products, and the signed
proton-plus-radiation bath ledger with the independent Decimal oracle of
run_total_capture.py. It reuses the build, row, endpoint and plotting checks
of run_pic_ground_capture.py, whose controlled geometry it shares.
"""
import argparse
import csv
from decimal import Decimal as D, localcontext
import json
import math
import os
from pathlib import Path
import sys
import time
import traceback

from run_total_capture import (CASE, EXIT, HERE, ROOT, TEST, MissingPrerequisite, dump, keyed_uniform,
                               precise, prepare_inputs, sha)
from run_pic_ground_capture import (POSITION, TEMPLATE, WEIGHT, binary_spacing, checked_build, rows,
                                    same_endpoint, write_comparison)

sys.path.insert(0, str(TEST))
from build_m2_qualification import sha256
from campaign_budget import CampaignBudget
from common.moon_testlib import _probe_compile_command
from common.python_environment import describe as python_environment, require as require_python

SEED = 731993
COHORT = 'I24-PIC-total-capture-v1'
REFERENCE_WEIGHTS = {1: 56., 2: .875}   # Na+ and H, as installed by the probe
REGISTRATION_CONTROLS = ('non_lso_frame_refused', 'below_envelope_refused', 'above_envelope_refused',
                         'negative_density_refused', 'zero_root_effort_refused', 'wrong_elementary_charge_refused')


@precise
def audit(folder, ranks, hit, dt, oracle, energy_ev, density, steps=1, weight_scale=1):
    """Independently reconstruct the event, its shell and the continued packets."""
    live, initial, events = [rows(folder, table, ranks) for table in ('particles', 'products', 'events')]
    checks = {}
    checks['one_event_two_committed_and_two_live_particles'] = [len(events), len(initial), len(live)] == [1, 2, 2]
    if not checks['one_event_two_committed_and_two_live_particles']:
        return checks, {}, {}
    event, = events
    by_species = lambda records: {int(r['species']): r for r in records}
    births, final = map(by_species, (initial, live))
    checks['only_sodium_ion_and_neutral_hydrogen_remain'] = set(births) == set(final) == {1, 2}
    if not checks['only_sodium_ion_and_neutral_hydrogen_remain']:
        return checks, {}, {}
    with (folder / 'registration-controls.csv').open() as stream:
        controls, = csv.DictReader(stream)
    for name in REGISTRATION_CONTROLS:
        checks['native_' + name] = controls.get(name) == '1'
    origin, epoch = map(float, (event['origin_et'], event['epoch_et']))
    time_bound = 16 * max(binary_spacing(value) for value in (origin, epoch, origin + dt, dt))
    increment = dt / steps
    hit_step = math.ceil(hit / increment)
    hit_within_step = hit - (hit_step - 1) * increment
    absolute_epoch_error = abs((epoch - origin) - hit)
    checks['source_key_and_event_epoch'] = (int(event['parent']) == 0
        and abs(float(event['elapsed_s']) - hit_within_step) <= time_bound and absolute_epoch_error <= time_bound)
    by_step = {}
    for sample in rows(folder, 'steps', ranks):
        by_step.setdefault(int(sample['step']), []).append(sample)
    checks['committed_step_history'] = set(by_step) == set(range(1, steps + 1))
    for number, samples in by_step.items():
        before = number < hit_step
        counts = tuple(sum(int(row[key]) for row in samples) for key in
                       ('committed_events', 'neutral_parents', 'sodium_ions', 'hydrogen'))
        population = tuple(sum(float(row[key]) for row in samples) for key in
                           ('represented_Na', 'represented_Na_ion', 'represented_H'))
        checks['committed_step_history'] &= (len(samples) == ranks and
            all(float(row['relative_end_s']) == number * increment for row in samples) and
            counts == ((0, 1, 0, 0) if before else (1, 0, 1, 1)) and
            population == ((WEIGHT, 0, 0) if before else (0, WEIGHT, WEIGHT)))
    tol = 5e-12
    g = oracle.speed(energy_ev)
    rate = float(D(density) * oracle.sigma(D(energy_ev)) * g)
    fractions = oracle.fractions(D(energy_ev))
    uniform = keyed_uniform(SEED, COHORT, 0, 4, 0)
    shell = 1 if uniform < fractions[0] else 2 if uniform < fractions[0] + fractions[1] else 3
    ref = oracle.products(g, shell)
    checks['morgan_rate'] = abs(float(event['rate_s']) / rate - 1) < tol
    # The owner's keyed Kinematics draw must be the independent HMAC value.
    checks['keyed_shell_uniform_and_shell'] = float(event['shell_uniform']) == uniform and int(event['shell']) == shell
    checks['shell_Q'] = abs(float(event['Q_J']) / float(ref['q']) - 1) < tol
    excitation = float(oracle.excitation[shell - 1])
    checks['radiated_energy'] = (float(event['radiated_J']) == 0. if shell == 1 else
                                 abs(float(event['radiated_J']) / excitation - 1) < tol)
    velocity = {}
    max_position_error = 0.
    for species in (1, 2):
        wref = REFERENCE_WEIGHTS[species] * weight_scale
        begin, end = births[species], final[species]
        v = [float(end['v' + axis]) for axis in 'xyz']
        velocity[species] = v
        prefix = f'species_{species}_'
        checks[prefix + 'lineage'] = (int(begin['parent']) == int(end['parent']) == 0
                                      and int(begin['ordinal']) == int(end['child']))
        checks[prefix + 'weights'] = all(float(r['represented']) == WEIGHT and float(r['reference_weight']) == wref
                                         and float(r['correction']) * wref == WEIGHT for r in (begin, end))
        checks[prefix + 'velocity_after_continuation'] = all(float(begin['v' + a]) == v[i] for i, a in enumerate('xyz'))
        checks[prefix + 'clock'] = (abs(float(begin['epoch_et']) - epoch) <= time_bound
            and abs(float(begin['remaining_s']) - (hit_step * increment - hit)) <= time_bound
            and abs(float(end['epoch_et']) - (origin + dt)) <= time_bound)
        checks[prefix + 'continued_position'] = all(float(begin[a]) == POSITION[i]
            and abs(float(end[a]) - POSITION[i] - v[i] * (dt - hit)) <= abs(v[i]) * time_bound + 8 * binary_spacing(POSITION[i])
            for i, a in enumerate('xyz'))
        max_position_error = max(max_position_error, *(abs(float(end[a]) - POSITION[i] - v[i] * (dt - hit))
                                                        for i, a in enumerate('xyz')))
    mp, mi, mh, speed = map(float, (oracle.proton, oracle.ion, oracle.h, g))
    incident = [mp * speed * f for f in (.6, .8, 0.)]
    kinetic = .5 * mp * speed ** 2
    residual = math.dist([mi * a + mh * b for a, b in zip(velocity[1], velocity[2])], incident) / (mp * speed)
    energy_ratio = abs((mi * sum(a * a for a in velocity[1]) + mh * sum(a * a for a in velocity[2])) * .5
                       - kinetic - float(ref['q'])) / kinetic
    checks['independent_product_momentum_energy'] = max(residual, energy_ratio) < tol
    hydrogen_error = mh * math.dist(velocity[2], [float(x) for x in ref['hydrogen']]) / (mp * speed)
    ion_error = mi * math.dist(velocity[1], [float(x) for x in ref['ion']]) / float(oracle.electron * g)
    checks['zero_angle_hydrogen_vs_oracle'] = hydrogen_error < tol
    checks['sodium_ion_recoil_vs_oracle'] = ion_error < json.loads(CASE.read_text())['recoil_relative_tolerance']
    # Signed bath per represented atom: one consumed proton plus the radiated
    # excitation energy. Charge fields are per physical event (integers).
    bath_energy = float(ref['bath'])
    checks['signed_proton_and_radiation_moments'] = (float(event['bath_number']) == WEIGHT
        and int(event['bath_charge_e']) == -1 and int(event['tracked_charge_e']) == 1
        and float(event['bath_number']) * int(event['bath_charge_e'])
            + float(final[1]['represented']) * int(event['tracked_charge_e']) == 0
        and abs(float(event['bath_mass_kg']) / (-WEIGHT * mp) - 1) < tol
        and math.dist([float(event['bath_p' + a]) for a in 'xyz'], [-WEIGHT * v for v in incident]) / (WEIGHT * mp * speed) < tol
        and abs(float(event['bath_energy_J']) / (WEIGHT * bath_energy) - 1) < tol)
    checks['closed_committed_ledger'] = (event['integer_ok'] == event['kinematic_ok'] == '1'
        and abs(float(event['bath_energy_J']) + float(event['tracked_energy_J'])) < WEIGHT * kinetic * tol)
    metrics = dict(event_elapsed_s=float(event['elapsed_s']), continued_s=dt - hit, clock_bound_s=time_bound,
                   step_count=steps, step_s=increment, reference_weight_scale=weight_scale, horizon_s=dt,
                   physical_represented_number=WEIGHT, event_epoch_et_s=epoch, source_relative_hit_s=hit,
                   absolute_event_epoch_error_s=absolute_epoch_error, maximum_position_error_m=max_position_error,
                   relative_momentum_error=residual, relative_energy_error=energy_ratio, shell=shell,
                   shell_uniform=uniform, hydrogen_velocity_error=hydrogen_error, sodium_ion_recoil_error=ion_error,
                   radiated_energy_J=float(event['radiated_J']))
    physical = [tuple(row[k] for k in ('parent', 'child', 'species', 'represented', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'epoch_et'))
                for _, row in sorted(final.items())]
    return checks, metrics, physical


def main(args):
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    case = json.loads(CASE.read_text())
    result = {'schema': 'srcMoon-I24-native-PIC-total-capture-v1', 'status': 'ERROR', 'checks': {}, 'metrics': {},
              'formal_I24_status': 'SKIPPED / NOT QUALIFIED'}
    manifest = {'source_revision': __import__('subprocess').check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                'runner': str(Path(__file__).resolve()), 'source_cohort': COHORT, 'seed': SEED, 'mpi_layouts': [1, 8],
                'openmp_threads': 1, 'commands': [], 'model': case['model_id'], 'python': python_environment()}
    try:
        if not 0 < args.wall_seconds <= 3600:
            raise ValueError('development limit is one hour')
        budget = CampaignBudget(args.wall_seconds, min(60, args.wall_seconds / 10))

        def execute(command, label, cwd, cores, keep_independent_cases=False):
            folder = output / (label + '-resources')
            folder.mkdir()
            receipt = budget.run(command, cwd, output / (label + '.log'), output=folder, allocated_cores=cores,
                                 environment=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1'))
            manifest['commands'].append(receipt)
            if receipt.get('reason') or (receipt['returncode'] and not keep_independent_cases):
                raise RuntimeError(label + ' did not complete; see ' + str(output / (label + '.log')))
            return receipt
        producer = args.producer_manifest.resolve()
        build = checked_build(producer)
        manifest['producer'] = {'path': str(producer), 'sha256': sha(producer)}
        energy_ev = D(case['pic']['energy_kev_u']) * 1000
        density = case['ensemble']['proton_density_m3']
        with localcontext() as context:
            context.prec = 70
            argv, oracle, hashes = prepare_inputs(case, output)
            # After input authentication (a corrupt input stays ERROR) and before any
            # compilation or MPI launch: a missing package is SKIPPED, not a late ERROR.
            require_python(('matplotlib',), MissingPrerequisite)  # PNG/EPS figures
            g = oracle.speed(energy_ev)
            rate = D(density) * oracle.sigma(energy_ev) * g
            hit = -D.from_float(keyed_uniform(SEED, COHORT, 0, 1, 0)).ln() / rate
            dt = float(hit + D(case['pic']['continuation_s']))
            # Energy conservation bounds the fastest possible outgoing H for a
            # stationary target: all available energy of the most exothermic
            # shell (n = 1) given to H.
            h_speed = ((oracle.proton * g ** 2 + 2 * oracle.q1) / oracle.h).sqrt()
        manifest['sources'] = hashes
        manifest['test_files'] = {str(p): sha(p) for p in (Path(__file__).resolve(), HERE / 'pic_total_capture_probe.cpp',
            HERE / 'total_capture_setup.h', HERE / 'run_total_capture.py', HERE / 'run_pic_ground_capture.py', TEMPLATE, CASE)}
        manifest['analytical_first_event_s'] = str(hit)
        manifest['interval_s'] = dt
        executable = output / 'native-pic-total-capture'
        command = _probe_compile_command(build, executable)
        command.insert(1, '-DMPI_ON')
        command[command.index(str(TEST / 'common/production_kernel_probe.cpp'))] = str(HERE / 'pic_total_capture_probe.cpp')
        execute(command, 'compile', ROOT, 1)
        manifest['probe_sha256'] = sha(executable)
        refusal = [str(v) for v in case['ensemble']['refusal_energies_ev_per_u']]
        baseline = None
        layouts = [(order, ranks, 1, 1., f'{order}-{ranks}mpi') for order in ('first', 'second') for ranks in (1, 8)]
        if args.refinement:
            for steps, scale, level in ((2, 1., 'steps2'), (4, 1., 'steps4'), (1, .5, 'weight-half'), (1, .25, 'weight-quarter')):
                layouts.extend(('second', ranks, steps, scale, f'second-{level}-{ranks}mpi') for ranks in (1, 8))
        manifest['refinement_requested'] = args.refinement
        manifest['layouts'] = [dict(label=label, order=order, mpi_ranks=ranks, steps=steps, weight_scale=scale)
                               for order, ranks, steps, scale, label in layouts]
        for order, ranks, steps, weight_scale, label in layouts:
            folder = output / label
            folder.mkdir()
            (folder / 'PT').mkdir()
            input_file = folder / 'moon.in'
            text = TEMPLATE.read_text().replace('force-cut-surface-second-order', 'force-cut-surface-' + order + '-order')
            text = text.replace('#end application',
                f'hydrogen_characteristic_speed_m_s = {h_speed}\n'
                f'hydrogen_to_sodium_reference_weight_ratio = {case["pic"]["hydrogen_to_sodium_weight_ratio"]}\n#end application')
            input_file.write_text(text.replace('tracing_maximum_substep_s = 8192', 'tracing_maximum_substep_s = ' + repr(dt)))
            receipt = execute(['mpiexec', '-n', str(ranks), str(executable), *argv, str(input_file), str(folder),
                               repr(dt), str(energy_ev), *refusal, str(steps), str(weight_scale)],
                              label, folder, ranks, keep_independent_cases=args.refinement)
            if receipt['returncode']:
                result.setdefault('layout_errors', {})[label] = dict(returncode=receipt['returncode'], log=receipt['log'])
                result['checks'][label] = {'completed_native_launch': False}
                result['metrics'][label] = {}
                continue
            paths = (folder / 'spice-kernels.txt').read_text().splitlines()
            if not paths:
                raise ValueError('empty list of loaded SPICE kernels')
            manifest.setdefault('spice_sha256', {})[label] = {p: sha256(Path(p)) for p in paths}
            checks, metrics, physical = audit(folder, ranks, float(hit), dt, oracle, energy_ev, density, steps, weight_scale)
            if physical:
                write_comparison(folder, physical, dt - float(hit))
                if label == 'first-1mpi':
                    baseline = physical
                elif baseline is None:
                    checks['first_order_one_rank_anchor_available'] = False
                elif steps == 1 and weight_scale == 1:
                    checks['phase_space_matches_first_order_one_rank'] = physical == baseline
                else:
                    checks['common_endpoint_matches_unrefined_packets'] = same_endpoint(
                        baseline, physical, metrics['clock_bound_s'])
            result['checks'][label] = checks
            result['metrics'][label] = metrics
        if args.refinement:
            path = output / 'refinement.csv'
            columns = ('case', 'order', 'ranks', 'steps', 'weight_scale', 'status', 'shell', 'source_relative_hit_s',
                       'event_epoch_et_s', 'absolute_event_epoch_error_s', 'clock_bound_s', 'maximum_position_error_m',
                       'step_s', 'horizon_s', 'physical_represented_number', 'relative_energy_error',
                       'relative_momentum_error', 'hydrogen_velocity_error', 'sodium_ion_recoil_error')
            with path.open('w', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=columns, extrasaction='ignore')
                writer.writeheader()
                for order, ranks, steps, scale, label in layouts:
                    status = ('ERROR' if label in result.get('layout_errors', {}) else
                              'PASS' if all(result['checks'][label].values()) else 'FAIL')
                    writer.writerow(dict(result['metrics'][label], case=label, order=order, ranks=ranks, steps=steps,
                                         weight_scale=scale, status=status))
        result['status'] = 'PASS' if all(all(c.values()) for c in result['checks'].values()) else 'FAIL'
        result['check_count'] = sum(len(c) for c in result['checks'].values())
        if result.get('layout_errors'):
            result.update(status='ERROR', reason='native launch error in ' + ', '.join(result['layout_errors']))
    except MissingPrerequisite as exc:
        result.update(status='SKIPPED', reason=str(exc))
    except Exception as exc:
        (output / 'exception.log').write_text(traceback.format_exc())
        result.update(status='ERROR', reason=str(exc))
    result['wall_seconds'] = time.monotonic() - start
    result['counts'] = {s: int(result['status'] == s) for s in EXIT}
    manifest['result_files'] = {str(p.relative_to(output)): sha(p) for p in output.rglob('*')
        if p.is_file() and (p.suffix in ('.csv', '.in', '.log', '.png', '.eps', '.txt') or p.name == 'exception.log')}
    dump(output / 'run_manifest.json', manifest)
    dump(output / 'result.json', result)
    print(json.dumps({k: v for k, v in result.items() if k not in ('checks', 'metrics')}, indent=2))
    return EXIT[result['status']]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--producer-manifest', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--refinement', action='store_true',
                        help='add second-order fixed-horizon timestep and fixed-population weight controls at 1 and 8 ranks')
    parser.add_argument('--wall-seconds', type=float, default=1800)
    sys.exit(main(parser.parse_args()))
