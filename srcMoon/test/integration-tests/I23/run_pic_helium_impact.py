#!/usr/bin/env python3
"""Run one actual I23 electron-bath ionization through each mover at 1 and 8 ranks.

A single stationary 4He packet sits in the force-free frozen LSO frame inside
the approved controlled bath (kT_e = 10 eV, n_e = 1e7 m^-3). The normal Moon
startup, force tracer, reaction controller, allocator and MPI exchange run
the event; the reviewed production composer supplies the He+ packet and the
signed electron reservoir. The trial is conditioned on the keyed source draw
of parent 0 so the product is observable; lifetime statistics belong to the
separate 32,768-parent ensemble (run_helium_impact.py). A zero-density control
must keep the parent and commit nothing. Optional refinement varies the common
timestep (1/2/4) and reference weights (1, 1/2, 1/4) at a fixed horizon.
"""
import argparse
import csv
from decimal import Decimal as D, localcontext
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / 'I24'))
from run_helium_impact import (CASE, EXIT, ROOT, TEST, MissingPrerequisite, dump, prepare_inputs, sha,
                               source_uniform)
from run_pic_ground_capture import binary_spacing, checked_build, rows, same_endpoint
sys.path.insert(0, str(TEST))
from build_m2_qualification import sha256
from campaign_budget import CampaignBudget
from common.moon_testlib import _probe_compile_command

TEMPLATE = HERE / 'reference/pic_helium_impact.in'
MAPPING = {'Na': '23Na:ground', 'He': '4He:ground', 'He_PLUS': '4He+:ground'}
HE, HE_PLUS = 1, 2  # compiled indices fixed by the ordered mapping above


def audit(folder, ranks, hit, dt, reference, case, steps=1, weight_scale=1, zero_density=False):
    pic = case['pic']
    represented = float(pic['represented_atoms'])
    position = tuple(float(v) for v in pic['position_m'])
    live, initial, events = [rows(folder, table, ranks) for table in ('particles', 'products', 'events')]
    checks = {}
    # Production owner refusals evaluated natively before startup (rank 0).
    with (folder / 'registration-controls.csv').open() as stream:
        controls, = csv.DictReader(stream)
    checks['owner_refuses_invalid_registration_before_startup'] = all(v == '1' for v in controls.values()) and len(controls) == 6
    origin = None
    by_step = {}
    for sample in rows(folder, 'steps', ranks):
        by_step.setdefault(int(sample['step']), []).append(sample)
    increment = dt / steps
    if zero_density:
        # No exposure: the parent survives the whole horizon untouched and the
        # controller commits no event, product or bath entry.
        checks['no_event_product_or_bath_entry'] = (len(events), len(initial)) == (0, 0)
        checks['only_the_original_helium_parent_remains'] = (len(live) == 1 and int(live[0]['species']) == HE
            and float(live[0]['represented']) == represented
            and all(float(live[0][a]) == position[i] for i, a in enumerate('xyz')))
        checks['committed_step_history'] = set(by_step) == set(range(1, steps + 1)) and all(
            len(samples) == ranks and sum(int(r['committed_events']) for r in samples) == 0
            and sum(int(r['helium']) for r in samples) == 1 and sum(int(r['helium_ions']) for r in samples) == 0
            for samples in by_step.values())
        physical = [tuple(r[k] for k in ('parent', 'child', 'species', 'represented', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'epoch_et'))
                    for r in live]
        return checks, dict(horizon_s=dt, step_count=steps, events=0), physical
    checks['one_event_one_committed_and_one_live_ion'] = [len(events), len(initial), len(live)] == [1, 1, 1]
    if not checks['one_event_one_committed_and_one_live_ion']:
        return checks, {}, []
    event, = events
    birth, final = initial[0], live[0]
    checks['only_the_helium_ion_remains'] = int(birth['species']) == int(final['species']) == HE_PLUS
    origin, epoch = map(float, (event['origin_et'], event['epoch_et']))
    time_bound = 16 * max(binary_spacing(v) for v in (origin, epoch, origin + dt, dt))
    hit_step = math.ceil(hit / increment)
    hit_within_step = hit - (hit_step - 1) * increment
    absolute_epoch_error = abs((epoch - origin) - hit)
    checks['source_key_and_event_epoch'] = (int(event['parent']) == 0
        and abs(float(event['elapsed_s']) - hit_within_step) <= time_bound and absolute_epoch_error <= time_bound)
    checks['committed_step_history'] = set(by_step) == set(range(1, steps + 1))
    for number, samples in by_step.items():
        before = number < hit_step
        counts = tuple(sum(int(r[k]) for r in samples) for k in ('committed_events', 'helium', 'helium_ions'))
        population = tuple(sum(float(r[k]) for r in samples) for k in ('represented_He', 'represented_He_ion'))
        checks['committed_step_history'] &= (len(samples) == ranks
            and all(float(r['relative_end_s']) == number * increment for r in samples)
            and counts == ((0, 1, 0) if before else (1, 0, 1))
            and population == ((represented, 0) if before else (0, represented)))
    tol = case['numerical_bounds']['arithmetic_relative_tolerance']
    ev, me, I = (float(reference[k]) for k in ('ev', 'me', 'threshold'))
    rate = float(reference['coefficient']) * float(case['bath']['primary_density_m3'])
    checks['bath_rate_n_e_k_T_e'] = abs(float(event['rate_s']) / rate - 1) <= 1e-12
    checks['state_Q_equals_minus_physical_threshold'] = abs(float(event['Q_J']) / -I - 1) <= tol
    wref = float(pic['reference_weights']['He_PLUS']) * weight_scale
    checks['ion_lineage'] = int(birth['parent']) == int(final['parent']) == 0 and int(birth['ordinal']) == int(final['child'])
    checks['ion_weights'] = all(float(r['represented']) == represented and float(r['reference_weight']) == wref
                                and float(r['correction']) * wref == represented for r in (birth, final))
    # Inherited parent velocity is exactly zero, so the He+ must stay at the
    # parent position for the whole remaining flight.
    checks['ion_inherits_zero_velocity_and_stays_at_parent_position'] = all(
        float(r['v' + a]) == 0 for r in (birth, final) for a in 'xyz') and all(
        float(r[a]) == position[i] for r in (birth, final) for i, a in enumerate('xyz'))
    checks['ion_clock'] = (abs(float(birth['epoch_et']) - epoch) <= time_bound
        and abs(float(birth['remaining_s']) - (hit_step * increment - hit)) <= time_bound
        and abs(float(final['epoch_et']) - (origin + dt)) <= time_bound)
    # Signed electron reservoir: one incoming and two outgoing electrons per
    # physical event; ledger moments are represented-population totals.
    checks['signed_bath_electron_numbers'] = (float(event['bath_reactant_number']) == represented
        and float(event['bath_product_number']) == 2 * represented and float(event['bath_number_change']) == 1)
    checks['signed_bath_charge_closes_with_ion'] = (int(event['bath_charge_e']) == -1 and int(event['tracked_charge_e']) == 1
        and abs(float(event['bath_charge_change_C']) / -ev - 1) <= tol
        # removing N incoming electrons (-e each) adds +N e to the bath; adding
        # 2N outgoing electrons adds -2N e; the tracked He+ carries +N e.
        and float(event['bath_reactant_number']) - float(event['bath_product_number'])
            + float(final['represented']) == 0)
    checks['signed_bath_mass_energy_momentum'] = (abs(float(event['bath_mass_kg']) / (me * represented) - 1) <= tol
        and abs(float(event['bath_energy_J']) / (-I * represented) - 1) <= tol
        and abs(float(event['tracked_energy_J']) / (I * represented) - 1) <= tol
        and all(float(event['bath_p' + a]) == 0 and float(event['tracked_p' + a]) == 0 for a in 'xyz'))
    checks['incoming_and_outgoing_electron_energies'] = (abs(float(event['in_kinetic_J']) / float(reference['e_in']) - 1) <= tol
        and abs(float(event['out_kinetic_J']) / float(reference['e_out']) - 1) <= tol)
    checks['closed_committed_ledger'] = (event['integer_ok'] == event['kinematic_ok'] == '1'
        and abs(float(event['bath_energy_J']) + float(event['tracked_energy_J'])) <= I * represented * tol)
    checks['recoil_bound_within_approved_admission'] = (
        abs(float(event['heavy_velocity_bound_m_s']) / float(reference['upper_velocity']) - 1) <= 1e-9
        and float(event['heavy_velocity_bound_m_s']) <= float(reference['velocity_limit']))
    metrics = dict(event_elapsed_s=float(event['elapsed_s']), continued_s=dt - hit, clock_bound_s=time_bound,
                   step_count=steps, step_s=increment, reference_weight_scale=weight_scale, horizon_s=dt,
                   physical_represented_number=represented, event_epoch_et_s=epoch, source_relative_hit_s=hit,
                   absolute_event_epoch_error_s=absolute_epoch_error)
    physical = [tuple(final[k] for k in ('parent', 'child', 'species', 'represented', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'epoch_et'))]
    return checks, metrics, physical


def write_step_figure(output, histories):
    """Represented He and He+ after each committed step, every layout."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    with (output / 'step_populations.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(('case', 'step', 'relative_end_s', 'represented_He', 'represented_He_ion', 'committed_events'))
        figure, axis = plt.subplots(figsize=(6, 3.8))
        for label, samples, hit, dt in histories:
            ends, he, ion = [0.], [None], [None]
            for step in sorted(samples):
                group = samples[step]
                end = float(group[0]['relative_end_s'])
                rhe = sum(float(r['represented_He']) for r in group); rion = sum(float(r['represented_He_ion']) for r in group)
                writer.writerow((label, step, end, rhe, rion, sum(int(r['committed_events']) for r in group)))
                ends.append(end / dt); he.append(rhe); ion.append(rion)
            axis.plot(ends[1:], he[1:], 'o-', ms=3, color='#264b79')
            axis.plot(ends[1:], ion[1:], 's--', ms=3, color='#ac4615')
        axis.plot([], [], 'o-', color='#264b79', label='represented He')
        axis.plot([], [], 's--', color='#ac4615', label='represented He$^+$')
        axis.set(xlabel='Step end / horizon', ylabel='Represented atoms',
                 title='I23 linked PIC trial (controlled bath; not validation)')
        axis.legend(frameon=False)
        figure.tight_layout()
        for ext in ('png', 'eps'):
            figure.savefig(output / ('step_populations.' + ext), dpi=150)
        plt.close(figure)


def main(args):
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    case = json.loads(CASE.read_text())
    pic = case['pic']
    result = {'schema': 'srcMoon-I23-native-PIC-electron-bath-v1', 'status': 'ERROR', 'checks': {}, 'metrics': {},
              'formal_I23_status': 'SKIPPED / NOT QUALIFIED'}
    manifest = {'source_revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                'runner': str(Path(__file__).resolve()), 'source_cohort': pic['cohort'], 'seed': pic['seed'],
                'mpi_layouts': [1, 8], 'openmp_threads': 1, 'commands': []}
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
        build = checked_build(producer, MAPPING, 'Na, He, He+')
        manifest['producer'] = {'path': str(producer), 'sha256': sha(producer)}
        bundle, ref, hashes = prepare_inputs(case, args.profile, args.data_root, output)
        manifest['sources'] = hashes
        with localcontext() as context:
            context.prec = 70
            rate = D(case['bath']['primary_density_m3']) * ref['coefficient']
            hit = -D.from_float(source_uniform(pic['seed'], pic['cohort'], 0)).ln() / rate
            dt = float(hit + D(pic['continuation_after_event_s']))
        manifest['test_files'] = {str(p): sha(p) for p in (Path(__file__).resolve(), HERE / 'pic_helium_impact_probe.cpp',
            HERE / 'helium_impact_setup.h', HERE / 'run_helium_impact.py', TEMPLATE, CASE)}
        manifest['analytical_first_event_s'] = str(hit)
        manifest['interval_s'] = dt
        executable = output / 'native-pic-helium-impact'
        command = _probe_compile_command(build, executable)
        command.insert(1, '-DMPI_ON')
        command[command.index(str(TEST / 'common/production_kernel_probe.cpp'))] = str(HERE / 'pic_helium_impact_probe.cpp')
        execute(command, 'compile', ROOT, 1)
        manifest['probe_sha256'] = sha(executable)
        baseline = None
        layouts = [(order, ranks, 1, 1., '1e7', f'{order}-{ranks}mpi') for order in ('first', 'second') for ranks in (1, 8)]
        layouts += [('second', ranks, 1, 1., '0', f'second-zero-density-{ranks}mpi') for ranks in (1, 8)]
        if args.refinement:
            for steps, scale, level in ((2, 1., 'steps2'), (4, 1., 'steps4'), (1, .5, 'weight-half'), (1, .25, 'weight-quarter')):
                layouts.extend(('second', ranks, steps, scale, '1e7', f'second-{level}-{ranks}mpi') for ranks in (1, 8))
        manifest['refinement_requested'] = args.refinement
        manifest['layouts'] = [dict(label=label, order=order, mpi_ranks=ranks, steps=steps, weight_scale=scale, density_m3=n)
                               for order, ranks, steps, scale, n, label in layouts]
        weights = pic['reference_weights']
        histories = []
        for order, ranks, steps, weight_scale, density, label in layouts:
            folder = output / label
            folder.mkdir()
            (folder / 'PT').mkdir()
            input_file = folder / 'moon.in'
            text = TEMPLATE.read_text().replace('force-cut-surface-second-order', 'force-cut-surface-' + order + '-order')
            text = text.replace('#end application',
                f"helium_ion_characteristic_speed_m_s = {pic['helium_ion_characteristic_speed_m_s']}\n"
                f"helium_ion_to_sodium_reference_weight_ratio = {pic['helium_ion_to_sodium_reference_weight_ratio']}\n"
                f"helium_to_sodium_reference_weight_ratio = {pic['helium_to_sodium_reference_weight_ratio']}\n#end application")
            input_file.write_text(text.replace('tracing_maximum_substep_s = 8192', 'tracing_maximum_substep_s = ' + repr(dt)))
            argv = [*bundle, f'input={input_file}', f'output={folder}', f'horizon_s={dt!r}', f'steps={steps}',
                    f'weight_scale={weight_scale!r}', f'density_m3={density}', f"seed={pic['seed']}", f"cohort={pic['cohort']}",
                    f"x_m={pic['position_m'][0]}", f"y_m={pic['position_m'][1]}", f"z_m={pic['position_m'][2]}",
                    f"represented={pic['represented_atoms']}", f"weight_na={weights['Na']}", f"weight_he={weights['He']}",
                    f"weight_he_plus={weights['He_PLUS']}"]
            receipt = execute(['mpiexec', '-n', str(ranks), str(executable), *argv], label, folder, ranks,
                              keep_independent_cases=True)
            if receipt['returncode']:
                result.setdefault('layout_errors', {})[label] = dict(returncode=receipt['returncode'], log=receipt['log'])
                result['checks'][label] = {'completed_native_launch': False}
                result['metrics'][label] = {}
                continue
            paths = (folder / 'spice-kernels.txt').read_text().splitlines()
            if not paths:
                raise ValueError('empty list of loaded SPICE kernels')
            manifest.setdefault('spice_sha256', {})[label] = {p: sha256(Path(p)) for p in paths}
            checks, metrics, physical = audit(folder, ranks, float(hit), dt, {k: str(v) for k, v in ref.items()}, case,
                                              steps, weight_scale, zero_density=density == '0')
            samples = {}
            for sample in rows(folder, 'steps', ranks):
                samples.setdefault(int(sample['step']), []).append(sample)
            histories.append((label, samples, float(hit), dt))
            if physical and density != '0':
                if label == 'first-1mpi':
                    baseline = physical
                elif baseline is None:
                    checks['first_order_one_rank_anchor_available'] = False
                elif steps == 1 and weight_scale == 1:
                    checks['phase_space_matches_first_order_one_rank'] = physical == baseline
                else:
                    checks['common_endpoint_matches_unrefined_packets'] = same_endpoint(baseline, physical, metrics['clock_bound_s'])
            result['checks'][label] = checks
            result['metrics'][label] = metrics
        try:
            write_step_figure(output, histories)
        except ImportError as error:
            result['figure'] = 'SKIPPED: ' + str(error)
        path = output / 'refinement.csv'
        columns = ('case', 'order', 'ranks', 'steps', 'weight_scale', 'density_m3', 'status', 'source_relative_hit_s',
                   'event_epoch_et_s', 'absolute_event_epoch_error_s', 'clock_bound_s', 'step_s', 'horizon_s',
                   'physical_represented_number')
        with path.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=columns, extrasaction='ignore')
            writer.writeheader()
            for order, ranks, steps, scale, density, label in layouts:
                status = ('ERROR' if label in result.get('layout_errors', {}) else
                          'PASS' if all(result['checks'][label].values()) else 'FAIL')
                writer.writerow(dict(result['metrics'][label], case=label, order=order, ranks=ranks, steps=steps,
                                     weight_scale=scale, density_m3=density, status=status))
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
    parser.add_argument('--profile', type=Path, default=Path(json.loads(CASE.read_text())['profile_path']))
    parser.add_argument('--data-root', type=Path, default=Path('/data/vtenishe/moon_validation_data'))
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--refinement', action='store_true',
                        help='add second-order fixed-horizon timestep and fixed-population weight controls at 1 and 8 ranks')
    parser.add_argument('--wall-seconds', type=float, default=1800)
    sys.exit(main(parser.parse_args()))
