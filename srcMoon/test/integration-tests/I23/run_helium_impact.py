#!/usr/bin/env python3
"""I23 prescribed-electron-bath ensemble: production provider, independent audit.

Runs the compiled MPI adapter (helium_impact_ensemble_probe.cpp) that calls the
production HeliumElectronImpact::Evaluate/Compose and constant-hazard clock,
then scores every output against references computed here independently in
Decimal arithmetic from the authenticated sources:

* D05 Voronov coefficient and logarithmic slope (U16 reference expression);
* first-hit times -ln(U)/(n_e k) from the independently reproduced HMAC draws;
* endpoint-survival lifetime against the formal 3% I23 gate;
* every signed electron-reservoir moment (1 in, 2 out, net +1 e, -e, -I);
* the omitted-recoil bound at the upper end of the admitted slope error.

Status: PASS evaluated/met, FAIL evaluated/not met, ERROR infrastructure,
SKIPPED missing source/tool. A PASS here is supporting evidence only; formal
I23 stays SKIPPED / NOT QUALIFIED (predecessors and broader physics open).
"""
from __future__ import annotations
import argparse
import csv
from decimal import Decimal as D, localcontext
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
TEST = HERE.parents[1]
ROOT = TEST.parents[1]
CASE = HERE / 'reference/helium_case.json'
APPROVAL = HERE / 'reference/helium_admission_approval_20261125.json'
sys.path.insert(0, str(TEST / 'U16_electron_impact'))
sys.path.insert(0, str(TEST / 'sodium_cohort'))
sys.path.insert(0, str(HERE))
from voronov_verification import MissingPrerequisite, independent_rows, preflight, reference, reference_log_slope
from verify_numerical import source_uniform
sys.path.insert(0, str(TEST))  # srcMoon/test: shared common/ helpers
from common.python_environment import describe as python_environment, require as require_python
import prepare_helium_admission as admission_tool

EXIT = {'PASS': 0, 'FAIL': 1, 'ERROR': 2, 'SKIPPED': 77}
FRAME = 'LSO'
KINEMATIC_IDENTITY = 'I23-reviewed-helium-reservoir-kinematics-20261125'
SPECIES = ('HeliumElectronImpact', 'PhysicalReactionEvents', 'PhysicalStateProfile', 'ReducedKinematics',
           'IonCollector', 'ReactionExecution', 'PiecewiseHazard', 'QualificationRandom')


def sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def dump(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def profile_arguments(profile):
    args = []
    for name in ('isotopes', 'constants', 'states'):
        args += [f"{name}={profile[name]['path']}", f"{name}_sha256={profile[name]['sha256']}"]
    for source in profile['sources']:
        args.append('source=' + '|'.join(source[k] for k in ('id', 'path', 'sha256', 'role')))
    return args


def prepare_inputs(case, profile_path, data_root, output):
    """Authenticate sources, reproduce the frozen admission and build references."""
    profile_path = Path(profile_path)
    if not profile_path.is_file():
        raise MissingPrerequisite('helium state profile unavailable: ' + str(profile_path))
    if sha(profile_path) != case['profile_sha256']:
        raise ValueError('helium state profile SHA-256 mismatch')
    profile = json.loads(profile_path.read_text())
    hashes = {str(profile_path.resolve()): case['profile_sha256']}
    for descriptor in [profile[n] for n in ('isotopes', 'constants', 'states')] + profile['sources']:
        native = Path(descriptor['path'])
        if not native.is_file():
            raise MissingPrerequisite('profile source missing: ' + str(native))
        if sha(native) != descriptor['sha256']:
            raise ValueError('profile source changed: ' + str(native))
        hashes[str(native)] = descriptor['sha256']
    package, d05 = preflight(Path(data_root))
    table = package / 'raw/electron_impact/voronov/cfit.dat'
    if sha(table) != case['d05_cfit_sha256']:
        raise ValueError('D05 cfit.dat SHA-256 differs from the I23 case')
    hashes.update(d05)
    # Reproducibility: regenerate the admission and require byte equality
    # with the frozen, hash-pinned reference copies that Compose authenticates.
    regenerated = output / 'admission'
    result = admission_tool.run(profile_path, data_root, regenerated)
    if result['status'] != 'PASS':
        raise (MissingPrerequisite if result['status'] == 'SKIPPED' else ValueError)(
            'admission regeneration ' + result['status'] + ': ' + result.get('reason', ''))
    frozen = {}
    for key in ('crosswalk', 'model_error'):
        record = case['admission_files'][key]
        path = HERE / record['path']
        if sha(path) != record['sha256']:
            raise ValueError('frozen admission file changed: ' + str(path))
        if (regenerated / path.name).read_bytes() != path.read_bytes():
            raise ValueError('admission evidence is not reproduced byte-for-byte by the committed generator: ' + path.name)
        frozen[key] = (path, record['sha256'])
        hashes[str(path)] = record['sha256']
    hashes[str(APPROVAL)] = sha(APPROVAL)
    crosswalk = json.loads(frozen['crosswalk'][0].read_text())
    model = json.loads(frozen['model_error'][0].read_text())
    approval = json.loads(APPROVAL.read_text())
    states = {}
    with Path(profile['states']['path']).open(newline='') as stream:
        for row in csv.DictReader(stream):
            states[row['state_id']] = row
    constants = {}
    with Path(profile['constants']['path']).open(newline='') as stream:
        for row in csv.DictReader(stream):
            constants[row['name']] = D(row['value'])
    bounds = case['numerical_bounds']
    mech = model['mechanical']
    argv = [f'd05={table}', f"d05_sha256={case['d05_cfit_sha256']}", *profile_arguments(profile),
            f"thermal_energy_ev={case['bath']['thermal_energy_eV']}", f'frame={FRAME}',
            f"model_id={case['model_id']}", f"approval_identity={approval['schema']}:{approval['approval_date']}:{sha(APPROVAL)}",
            f"crosswalk={frozen['crosswalk'][0]}", f"crosswalk_sha256={frozen['crosswalk'][1]}",
            f"model_error={frozen['model_error'][0]}", f"model_error_sha256={frozen['model_error'][1]}",
            f"time_scale={case['time_scale']}", f"observable={model['observable']}",
            f"max_heavy_velocity_error_m_s={model['recoil']['admitted_maximum_heavy_velocity_error_m_s']}",
            f"max_heavy_energy_error_j={model['recoil']['admitted_maximum_heavy_kinetic_error_J']}",
            f"max_relative_mass_error={mech['relative_reference_inertia_error_max']}",
            f"max_classical_heavy_speed_m_s={mech['heavy_particle_speed_max_m_s']}",
            f"reservoir_mass_absolute_kg={bounds['reservoir_mass_absolute_kg']}",
            f"reservoir_mass_relative={bounds['reservoir_mass_relative']}",
            f"reservoir_time_projection_absolute_s={bounds['reservoir_time_projection_absolute_s']}",
            f"reservoir_momentum_absolute_kg_m_s={bounds['reservoir_momentum_absolute_kg_m_s']}",
            f"reservoir_momentum_relative={bounds['reservoir_momentum_relative']}",
            f"reservoir_energy_absolute_j={bounds['reservoir_energy_absolute_J']}",
            f"reservoir_energy_relative={bounds['reservoir_energy_relative']}",
            f"reservoir_unit_vector_tolerance={bounds['reservoir_unit_vector_tolerance']}",
            f"ledger_mass_absolute={bounds['ledger_mass'][0]}", f"ledger_mass_relative={bounds['ledger_mass'][1]}",
            f"ledger_momentum_absolute={bounds['ledger_momentum'][0]}", f"ledger_momentum_relative={bounds['ledger_momentum'][1]}",
            f"ledger_energy_absolute={bounds['ledger_energy'][0]}", f"ledger_energy_relative={bounds['ledger_energy'][1]}",
            f"kinematic_identity={KINEMATIC_IDENTITY}",
            f"threshold_crosswalk_max_difference_j={crosswalk['admitted_threshold_crosswalk_maximum_difference_J']}",
            f"max_abs_log_slope_error={model['conditional_mean']['admitted_maximum_absolute_log_slope_error']}",
            f"central_heavy_velocity_bound_m_s={model['recoil']['central']['rms_speed_bound_m_s']}"]
    with localcontext() as context:
        context.prec = 70
        fit = independent_rows(package)['He']
        temperature = D(case['bath']['thermal_energy_eV'])
        ev = constants['elementary charge']
        me = D(states['electron:free']['physical_rest_mass_kg'])
        m_ion = D(states['4He:ground']['physical_rest_mass_kg']) - me
        coefficient = reference(fit, temperature)
        slope = reference_log_slope(fit, temperature)
        threshold = D(states['4He+:ground']['formation_J'])
        e_in = temperature * ev * (D('1.5') + slope)
        upper_in = temperature * ev * (D('1.5') + slope + D(model['conditional_mean']['admitted_maximum_absolute_log_slope_error']))
        impulse = (2 * me * upper_in).sqrt() + (4 * me * (upper_in - threshold)).sqrt()
        ref = dict(coefficient=coefficient, slope=slope, ev=ev, me=me, m_ion=m_ion, threshold=threshold,
                   fit_threshold=fit['energy'], e_in=e_in, e_out=e_in - threshold,
                   upper_velocity=impulse / m_ion, upper_kinetic=impulse * impulse / (2 * m_ion),
                   velocity_limit=D(model['recoil']['admitted_maximum_heavy_velocity_error_m_s']),
                   kinetic_limit=D(model['recoil']['admitted_maximum_heavy_kinetic_error_J']),
                   mass_limit=D(mech['relative_reference_inertia_error_max']))
    return argv, ref, hashes


def relative(actual, wanted):
    actual, wanted = D(actual), D(wanted)
    return abs(actual / wanted - 1) if wanted else abs(actual)


def audit(directory, ranks, density, horizon, case, ref):
    ensemble = case['ensemble']
    count, tolerance = ensemble['parents'], D(repr(case['numerical_bounds']['arithmetic_relative_tolerance']))
    weight = D(ensemble['parent_weight'])
    checks = []

    def check(name, passed, actual, wanted):
        checks.append(dict(name=name, status='PASS' if passed else 'FAIL', actual=str(actual), reference=str(wanted)))
    with (directory / 'rate.csv').open() as stream:
        rate_row, = csv.DictReader(stream)
    with localcontext() as context:
        context.prec = 70
        coefficient = D(rate_row['coefficient_m3_s'])
        frequency = D(density) * ref['coefficient']
        check('coefficient vs independent Voronov (U16 1e-12)', relative(coefficient, ref['coefficient']) <= D('1e-12'),
              coefficient, ref['coefficient'])
        check('coefficient relative error <= 2% (formal I23 criterion)', relative(coefficient, ref['coefficient']) <= D('0.02'),
              relative(coefficient, ref['coefficient']), '0.02')
        rate = D(rate_row['rate_s'])
        check('rate n_e*k(T_e) [1/s]', (rate == 0) if density == 0 else relative(rate, frequency) <= D('1e-12'), rate, frequency)
        check('log slope d ln k / d ln T', relative(rate_row['log_slope'], ref['slope']) <= D('1e-12'), rate_row['log_slope'], ref['slope'])
        check('mean incident energy kT(3/2+slope) [J]', relative(rate_row['mean_incoming_energy_J'], ref['e_in']) <= D('1e-12'),
              rate_row['mean_incoming_energy_J'], ref['e_in'])
        check('physical threshold [J]', relative(rate_row['physical_threshold_J'], ref['threshold']) <= D('1e-15'),
              rate_row['physical_threshold_J'], ref['threshold'])
        for guard in ('negative_density_refused', 'nan_density_refused', 'nan_temperature_refused',
                      'negative_temperature_refused', 'unapproved_temperature_refused', 'drift_refused',
                      'moving_target_refused', 'missing_frame_refused'):
            check(guard, rate_row[guard] == '1', rate_row[guard], '1')
    rows = {}
    for rank in range(ranks):
        with (directory / f'parents-rank-{rank}.csv').open() as stream:
            for row in csv.DictReader(stream):
                if None in row or any(v is None for v in row.values()):
                    raise ValueError('malformed native parent CSV record')
                i = int(row['parent'])
                if not 0 <= i < count or i % ranks != rank or i in rows:
                    raise ValueError('duplicate or misplaced native parent ID')
                rows[i] = row
    if len(rows) != count:
        raise ValueError('native execution omitted a parent history')
    horizon_f = float(horizon)
    frequency_f = float(D(density) * ref['coefficient'])
    wrong_fates = 0
    worst_time = 0.
    events = []
    for i in range(count):
        row = rows[i]
        if row['fate'] not in ('ionized', 'survived'):
            raise ValueError('unknown native fate')
        hit = row['fate'] == 'ionized'
        t_first = -math.log(source_uniform(ensemble['seed'], ensemble['cohort'], i)) / frequency_f if frequency_f > 0 else math.inf
        wrong_fates += (t_first <= horizon_f) != hit
        worst_time = max(worst_time, abs(float(row['time_s']) - min(t_first, horizon_f)) / horizon_f)
        if hit:
            events.append(row)
    check('all parent first-hit classifications', wrong_fates == 0, wrong_fates, 0)
    check('event-time errors / horizon (arithmetic)', worst_time <= float(tolerance), worst_time, float(tolerance))
    survivors = count - len(events)
    if density == 0:
        check('zero density retains every parent', survivors == count and not events, survivors, count)
        return rows, checks, []
    fraction = survivors / count
    lifetime = -horizon_f / math.log(fraction) if 0 < fraction < 1 else None
    lifetime_error = abs(lifetime * frequency_f - 1) if lifetime else None
    check('endpoint-survival lifetime relative error <= 3% (formal I23 criterion)',
          lifetime_error is not None and lifetime_error <= case['acceptance']['ensemble_lifetime_relative_error_max'],
          lifetime_error, case['acceptance']['ensemble_lifetime_relative_error_max'])
    with (directory / 'event-controls.csv').open() as stream:
        guards, = csv.DictReader(stream)
    for name, value in guards.items():
        check(name, value == '1', value, '1')
    bad = {}

    def flag(name, ok):
        bad[name] = bad.get(name, 0) + (not ok)
    ev, me, I = ref['ev'], ref['me'], ref['threshold']
    with localcontext() as context:
        context.prec = 50
        for row in events:
            get = lambda k: D(row[k])
            flag('He+ product identity, zero inherited velocity and represented weight',
                 row['ion_species'] == '1' and all(get('ion_v' + a) == 0 for a in 'xyz')
                 and get('ion_reference_weight') == D(ensemble['ion_weight'])
                 and get('ion_represented') == weight and get('ion_correction') * get('ion_reference_weight') == weight)
            flag('bath terms: 1 incoming and 2 outgoing electrons per represented parent',
                 get('bath_reactant_represented') == weight and get('bath_product_represented') == 2 * weight)
            flag('integer topology and charge (+1 tracked, -1 bath)',
                 row['integer_ok'] == '1' and row['kinematic_ok'] == '1'
                 and row['tracked_charge_e'] == '1' and row['bath_charge_e'] == '-1')
            flag('incoming electron moment = kT(3/2+slope), zero momentum',
                 get('in_number') == 1 and relative(get('in_mass_kg'), me) <= tolerance
                 and relative(get('in_charge_C'), -ev) <= tolerance and relative(get('in_kinetic_J'), ref['e_in']) <= tolerance
                 and all(abs(get('in_p' + a)) <= D('1e-40') for a in 'xyz'))
            flag('outgoing aggregate: 2 electrons with E_in - I, zero momentum',
                 get('out_number') == 2 and relative(get('out_mass_kg'), 2 * me) <= tolerance
                 and relative(get('out_charge_C'), -2 * ev) <= tolerance and relative(get('out_kinetic_J'), ref['e_out']) <= tolerance
                 and all(abs(get('out_p' + a)) <= D('1e-40') for a in 'xyz'))
            flag('signed net bath change: +1 e, +m_e, -e, -I, zero momentum',
                 get('bath_number_change') == 1 and relative(get('bath_mass_change_kg'), me) <= tolerance
                 and relative(get('bath_charge_change_C'), -ev) <= tolerance
                 and relative(get('bath_energy_change_J'), -I) <= tolerance and relative(get('state_Q_J'), -I) <= tolerance)
            # Committed EventLedger moments are represented-population totals
            # (N physical events per packet), as in the reviewed I24 audit;
            # its integer stoichiometry/charge fields are per ONE event.
            flag('committed ledger closes with the reservoir (tracked +N I, bath -N I, +N m_e)',
                 relative(get('ledger_bath_energy_J'), -I * weight) <= tolerance
                 and relative(get('ledger_tracked_energy_J'), I * weight) <= tolerance
                 and relative(get('ledger_bath_mass_kg'), me * weight) <= tolerance
                 and all(abs(get('ledger_bath_p' + a)) <= D('1e-40') for a in 'xyz'))
            flag('represented bath charge/energy = N x per-event change',
                 relative(get('represented_bath_charge_C'), -ev * weight) <= tolerance
                 and relative(get('represented_bath_energy_J'), -I * weight) <= tolerance)
            flag('upper-end recoil bound reproduced and within the approved admission',
                 relative(get('heavy_velocity_bound_m_s'), ref['upper_velocity']) <= D('1e-9')
                 and relative(get('heavy_kinetic_bound_J'), ref['upper_kinetic']) <= D('1e-9')
                 and get('heavy_velocity_bound_m_s') <= ref['velocity_limit']
                 and get('heavy_kinetic_bound_J') <= ref['kinetic_limit']
                 and get('relative_mass_error') <= ref['mass_limit'])
    for name, failures in bad.items():
        check('every event: ' + name, failures == 0, failures, 0)
    check('ionization events present', len(events) > 0, len(events), '>0')
    times = sorted(float(r['time_s']) / horizon_f for r in events)
    return rows, checks, times


def plot(output, curves):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axis = plt.subplots(figsize=(5.6, 4), constrained_layout=True)
    grid = [j / 200 for j in range(201)]
    for (label, times, parents), color in zip(curves, ('#264b79', '#ac4615')):
        axis.step([0] + times + [1], [1] + [1 - (j + 1) / parents for j in range(len(times))] + [1 - len(times) / parents],
                  where='post', color=color, label='Native He, ' + label)
    axis.plot(grid, [math.exp(-v) for v in grid], 'k--', label='exp(-n_e k t)')
    axis.set(xlabel='Elapsed time / (n_e k)$^{-1}$', ylabel='Surviving He fraction',
             title='I23 controlled 10 eV bath (development, not validation)')
    axis.legend(framealpha=1)
    for suffix in ('png', 'eps'):
        fig.savefig(output / ('comparison.' + suffix))
    plt.close(fig)


def run(profile, data_root, output, ranks=8, wall_seconds=3600):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    case = json.loads(CASE.read_text())
    manifest = dict(schema='srcMoon-I23-helium-ensemble-manifest-v1',
                    scope='compiled MPI production provider, constant-hazard clock and reviewed reservoir composer',
                    source_revision=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                    seed=case['ensemble']['seed'], cohort=case['ensemble']['cohort'], parents=case['ensemble']['parents'],
                    mpi_layouts=list(dict.fromkeys((1, ranks))), openmp_threads=1, commands=[],
                    model=case['model_id'], formal_status='SKIPPED / NOT QUALIFIED',
                    python=python_environment())
    checks = []
    result = dict(status='ERROR', scope=manifest['scope'], formal_status=manifest['formal_status'], checks=checks)
    stdout, stderr = output / 'stdout.log', output / 'stderr.log'
    stdout.touch(); stderr.touch()
    try:
        if not 0 < wall_seconds <= 3600 or not math.isfinite(wall_seconds) or ranks < 1:
            raise ValueError('positive MPI ranks and a total wall allowance of at most 3600 s are required')

        def invoke(argv, directory):
            remaining = wall_seconds - (time.monotonic() - start)
            if remaining <= 0:
                raise TimeoutError('I23 total time allowance exhausted')
            argv = list(map(str, argv))
            manifest['commands'].append({'cwd': str(directory), 'argv': argv})
            with stdout.open('a') as out, stderr.open('a') as err:
                out.write(json.dumps(argv) + '\n'); out.flush()
                process = subprocess.run(argv, cwd=directory, stdout=out, stderr=err, timeout=remaining,
                                         env=dict(os.environ, OMP_NUM_THREADS='1', OMPI_MCA_hwloc_base_binding_policy='none'))
            manifest['commands'][-1]['returncode'] = process.returncode
            if process.returncode == 77:
                raise MissingPrerequisite('native provider reported an unavailable source; see ' + str(stderr))
            if process.returncode:
                raise RuntimeError('native command exited %d (see %s)' % (process.returncode, stderr))
        bundle, ref, hashes = prepare_inputs(case, profile, data_root, output)
        # After input authentication (a corrupt input stays ERROR) and before any
        # compilation or MPI launch: a missing package is SKIPPED, not a late ERROR.
        require_python(('matplotlib',), MissingPrerequisite)  # PNG/EPS figures
        manifest['input_sha256'] = hashes
        shutil.copyfile(CASE, output / 'case.json')
        dump(output / 'independent_reference.json', {k: str(v) for k, v in ref.items()})
        compiler, launcher = shutil.which('mpicxx'), shutil.which('mpirun')
        if not compiler or not launcher:
            raise MissingPrerequisite('mpicxx and mpirun required for the I23 MPI ensemble')
        manifest['compiler'] = subprocess.check_output([compiler, '--version'], text=True).splitlines()[0]
        manifest['mpi_launcher'] = subprocess.check_output([launcher, '--version'], text=True).splitlines()[0]
        species = ROOT / 'srcMoon/species'
        shared = ROOT / 'src/models/chemical_database'
        sources = [HERE / 'helium_impact_ensemble_probe.cpp', *(species / (n + '.cpp') for n in SPECIES),
                   ROOT / 'srcMoon/MoonTransportClock.cpp',
                   *(shared / (n + '.cpp') for n in ('DataIO', 'Photochemistry', 'ElectronImpact'))]
        flags = ['-std=c++17', '-O2', '-g', '-Wall', '-Wextra', '-pedantic', '-DOMPI_SKIP_MPICXX', '-I', str(ROOT)]
        dep = subprocess.run([compiler, *flags, '-MM', *map(str, sources)], cwd=ROOT, capture_output=True,
                             text=True, check=True, timeout=60)
        (output / 'compiler-dependencies.txt').write_text(dep.stdout)
        dependencies = set()
        for line in dep.stdout.replace('\\\n', ' ').splitlines():
            dependencies.update((ROOT / p).resolve() for p in line.split(':', 1)[1].split())
        dependencies.update((Path(__file__).resolve(), HERE / 'prepare_helium_admission.py', CASE, APPROVAL,
                             TEST / 'sodium_cohort/verify_numerical.py', TEST / 'U16_electron_impact/voronov_verification.py'))
        manifest['source_sha256'] = {str(p): sha(p) for p in sorted(dependencies)}
        binary = output / 'helium-impact-ensemble-probe'
        invoke([compiler, *flags, *sources, '-o', binary], ROOT)
        manifest['binary_sha256'] = sha(binary)
        ensemble = case['ensemble']
        primary = D(case['bath']['primary_density_m3'])
        curves = []
        for density_text in ('1e7', '1e6', '0'):
            density = D(density_text)
            with localcontext() as context:
                context.prec = 50
                # One optical depth of prescribed exposure; zero density reuses
                # the primary horizon so survival is tested over the same time.
                horizon = D(ensemble['optical_depth']) / ((density if density else primary) * ref['coefficient'])
            first = None
            for size in manifest['mpi_layouts']:
                directory = output / f'ne-{density_text}' / f'mpi-{size}'
                directory.mkdir(parents=True)
                invoke([launcher, '-np', size, binary, *bundle, f'density_m3={density_text}', f'horizon_s={horizon}',
                        f"seed={ensemble['seed']}", f"cohort={ensemble['cohort']}", f"parents={ensemble['parents']}",
                        f"parent_weight={ensemble['parent_weight']}", f"ion_weight={ensemble['ion_weight']}",
                        f'output={directory}'], directory)
                rows, metrics, times = audit(directory, size, float(density), horizon, case, ref)
                checks.extend(dict(c, density_m3=density_text, mpi_ranks=size) for c in metrics)
                if first is None:
                    first = rows
                    if density:
                        curves.append((f'n_e={density_text} m$^{{-3}}$', times, ensemble['parents']))
                else:
                    differences = sum(row != first[i] for i, row in rows.items())
                    checks.append(dict(name='per-parent MPI rank invariance', status='PASS' if not differences else 'FAIL',
                                       actual=str(differences), reference='0', density_m3=density_text, mpi_ranks=size))
        try:
            plot(output, curves)
        except ImportError as error:
            checks.append(dict(name='comparison figure', status='SKIPPED', actual=str(error), reference='matplotlib'))
        with (output / 'comparison.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=('density_m3', 'mpi_ranks', 'name', 'status', 'actual', 'reference'))
            writer.writeheader(); writer.writerows(checks)
        if any(sha(Path(p)) != digest for p, digest in manifest['source_sha256'].items()):
            raise RuntimeError('production provider or adapter changed during execution')
        dump(output / 'native_receipt.json', {'histories': {str(p.relative_to(output)): sha(p)
             for p in sorted(output.glob('ne-*/mpi-*/*.csv'))}, 'scope': manifest['scope']})
        failures = [c['name'] for c in checks if c['status'] == 'FAIL']
        result.update(status='FAIL' if failures else 'PASS',
                      reason=', '.join(sorted(set(failures))) if failures else 'all rate, clock, MPI and reservoir checks met')
    except MissingPrerequisite as error:
        result.update(status='SKIPPED', reason=str(error))
    except (OSError, ValueError, ArithmeticError, KeyError, TypeError, RuntimeError, ImportError,
            subprocess.SubprocessError) as error:
        result.update(status='ERROR', reason=str(error))
        with stderr.open('a') as stream:
            stream.write(str(error) + '\n')
    result['duration_s'] = time.monotonic() - start
    result['counts'] = {status: sum(c['status'] == status for c in checks) for status in EXIT}
    dump(output / 'run_manifest.json', manifest); dump(output / 'result.json', result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile', type=Path, default=Path(json.loads(CASE.read_text())['profile_path']))
    parser.add_argument('--data-root', type=Path, default=Path('/data/vtenishe/moon_validation_data'))
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--mpi-ranks', type=int, default=8)
    parser.add_argument('--test-runtime-seconds', type=float, default=3600)
    args = parser.parse_args()
    value = run(args.profile, args.data_root, args.output_dir, args.mpi_ranks, args.test_runtime_seconds)
    print(json.dumps({k: v for k, v in value.items() if k != 'checks'}, indent=2))
    sys.exit(EXIT[value['status']])
