#!/usr/bin/env python3
"""Compile the actual source-bound O I service; evaluate independent local oracles.

No network, application configuration, isotope selection or observations. A
local PASS cannot admit a lunar cascade provider; the result states that gate
separately. Source and output files are never repaired or overwritten.
"""
import argparse
import csv
from decimal import Decimal, localcontext
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import time

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DATABASE = Path('/data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261119-r3')
LEVEL_SHA = '8b0bfcbae54365ed668efa35f337e7b9ed6a60a55d0b4fbfe94693e22660c192'
LINK_SHA = 'bf4e17705e69cbaebbabd07f624e6034c5d96817b10aeb51ee8a77c96a21ed4f'
RELATIVE = 128 * sys.float_info.epsilon
N = 20000
ALPHA_PER_CLAIM = .01 / 112

class Missing(Exception):
    pass

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def strict_json(path):
    def pairs(items):
        d = {}
        for k, v in items:
            if k in d:
                raise ValueError('duplicate JSON key ' + k)
            d[k] = v
        return d
    return json.loads(Path(path).read_text(), object_pairs_hook=pairs,
                      parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))

def run(output, database, wall_seconds):
    started = time.monotonic()
    output.mkdir(parents=True, exist_ok=False)
    checks, identities, commands = [], [], []
    def check(name, passed, **metrics):
        checks.append(dict(name=name, status='PASS' if passed else 'FAIL', **metrics))
    def close(a, b):
        a, b = float(a), float(b)
        return math.isfinite(a) and math.isfinite(b) and (a == b or abs(a-b) <= RELATIVE*max(abs(a), abs(b)))
    def authenticate(path, expected):
        actual = sha(path)
        identities.append(dict(path=str(path), sha256=actual, expected_sha256=expected))
        if actual != expected:
            raise ValueError('source hash mismatch: ' + str(path))
    def execute(argv, log):
        remaining = wall_seconds-(time.monotonic()-started)
        if remaining <= 0:
            raise TimeoutError('local package wall budget exhausted')
        commands.append(dict(argv=list(map(str, argv)), cwd=str(REPO), log=str(log)))
        with log.open('w') as stream:
            p = subprocess.run(list(map(str, argv)), cwd=str(REPO), stdout=stream,
                               stderr=subprocess.STDOUT, timeout=remaining)
        if p.returncode:
            raise RuntimeError('native command returned {}: {}'.format(p.returncode, log))
    try:
        try:
            import numpy as np
            from scipy.linalg import expm
            from scipy.stats import binom
            import matplotlib
            matplotlib.use('Agg')
            from matplotlib import pyplot as plt
        except ImportError as e:
            raise Missing(str(e))
        compiler = shutil.which('g++')
        if not compiler:
            raise Missing('C++17 g++ compiler unavailable')
        if not (database/'database-lock.json').is_file():
            raise Missing('reviewed shared catalogue snapshot absent: ' + str(database))
        lock = strict_json(database/'database-lock.json')
        if lock.get('schema') != 'amps-chemical-database-lock-v1':
            raise ValueError('unsupported shared lock')
        export = database/lock['export_path']
        authenticate(export, lock['export_manifest_sha256'])
        # Detached export is a frozen self-contained scalar grammar. This
        # independent source reader interprets ONLY its literal index locators;
        # the maintained C++ parser authenticates the entire bundle separately.
        settings = {}
        for line in export.read_text().splitlines():
            if '=' in line:
                k, v = line.split('=', 1)
                settings[k.strip()] = v.strip()
        index = database/settings['table_index_path']
        authenticate(index, settings['table_index_sha256'])
        source_index = database/settings['source_index_path']
        authenticate(source_index, settings['source_index_sha256'])
        with source_index.open(newline='') as stream:
            originals = {row['source_id']: row for row in csv.DictReader(stream)}
        for sid, digest in [('chianti-oi:levels-native', '230bbdcd95c269fc0eea04a1fd37f534c1a42979d2e29c1b346e49dd71c688da'),
                            ('chianti-oi:radiative-native', 'f24efc2bd0fa5120baca6dce06ffb7831ae2a59e17cc186c233a43b5145866a8')]:
            if sid not in originals or originals[sid]['availability'] != 'available':
                raise Missing('complete permitted native source absent: '+sid)
            if originals[sid]['sha256'] != digest:
                raise ValueError('complete native source identity changed: '+sid)
            authenticate(database/originals[sid]['path'], digest)
        with index.open(newline='') as stream:
            tables = {row['table_id']: row for row in csv.DictReader(stream)}
        native = []
        for tid, digest in [('chianti-oi-levels-five', LEVEL_SHA),
                            ('chianti-oi-radiative-lowest-five', LINK_SHA)]:
            if tid not in tables:
                raise Missing('required permitted source table absent: ' + tid)
            descriptor = tables[tid]
            if descriptor['availability'] != 'available' or descriptor['evaluation_permission'] != 'allowed':
                raise Missing('required source table unavailable or permission-held: '+tid)
            if descriptor['sha256'] != digest:
                raise ValueError('frozen source identity changed: '+tid)
            path = database/descriptor['path']
            authenticate(path, digest)
            with path.open(newline='') as stream:
                native.append(list(csv.DictReader(stream)))
        levels, links = native
        if len(levels) != 5 or len(links) != 9:
            raise ValueError('selected graph is incomplete')
        sources = [REPO/'srcMoon/species/OxygenCascade.h', REPO/'srcMoon/species/OxygenCascade.cpp',
                   REPO/'src/models/chemical_database/Runtime.h', REPO/'src/models/chemical_database/Runtime.cpp',
                   REPO/'src/models/chemical_database/DataIO.h', REPO/'src/models/chemical_database/DataIO.cpp',
                   REPO/'src/models/chemical_database/ReducedModels.h', Path(__file__).resolve(),
                   REPO/'srcMoon/test/oxygen_cascade/probe.cpp', REPO/'srcMoon/test/oxygen_cascade/CONTRACT.md']
        identities.extend(dict(path=str(p), sha256=sha(p)) for p in sources)
        for p in [REPO/'srcMoon/development_and_validation.md',
                  Path('/data/vtenishe/moon_validation_data/development_and_validation.md'),
                  Path('/data/vtenishe/moon_validation_data/AGENTS.md')]:
            identities.append(dict(path=str(p), sha256=sha(p)))
        versions = dict(compiler=subprocess.check_output([compiler, '--version'], text=True, timeout=10).splitlines()[0],
                        python=sys.version, scipy=__import__('scipy').__version__, matplotlib=matplotlib.__version__)
        binary = output/'oxygen_cascade_probe'
        argv = [compiler, '-std=c++17', '-O2', '-Wall', '-Wextra', '-Werror',
                sources[1], sources[3], sources[5], sources[8], '-o', binary]
        execute(argv, output/'compile.log')
        identities.append(dict(path=str(binary), sha256=sha(binary)))
        execute([binary, database, lock['export_path'], lock['export_manifest_sha256'],
                 lock['release_id'], lock['release_manifest_sha256'], lock['profile_id'],
                 lock['profile_sha256']], output/'native.json')
        result = strict_json(output/'native.json')
        if (set(result['guards']) != {'level-domain-transaction', 'cdf-upper-exclusion',
                'stable-no-fabricated-decay', 'wait-open-support', 'missing-draws-transaction',
                'missing-source-preserves-network', 'source-hash-refusal', 'missing-constant-source'}
                or len(result['laws']) != 5 or len(result['links']) != 9 or len(result['recoil']) != 3
                or len(result['families']) != 8):
            raise ValueError('native receipt omitted a required predeclared control or source row')
        for name, value in result['guards'].items():
            check(name, value is True)
        with localcontext() as ctx:
            ctx.prec = 90
            h, c = Decimal('6.62607015e-34'), Decimal(299792458)
            rates = [Decimal(x['decay_rate_s_inverse_native']) for x in links]
            for i, law in enumerate(result['laws']):
                total = sum((rates[k] for k, row in enumerate(links) if int(row['upper_index']) == i+1), Decimal(0))
                check('source-total-A-{}'.format(i+1), close(law['rate'], total))
                check('source-lifetime-{}'.format(i+1), law['lifetime'] is None if total == 0 else close(law['lifetime'], 1/total))
            for k, (row, actual) in enumerate(zip(links, result['links'])):
                lo, hi = int(row['lower_index'])-1, int(row['upper_index'])-1
                total = sum((rates[j] for j, r in enumerate(links) if int(r['upper_index']) == hi+1), Decimal(0))
                check('source-link-A-{}'.format(k), close(actual['rate'], rates[k]))
                check('source-branch-{}'.format(k), close(actual['branch'], rates[k]/total))
                check('deterministic-positive-branch-{}'.format(k), actual['midpoint_selected'] == k)
                check('native-vacuum-wavelength-{}'.format(k), close(actual['vacuum_angstrom'], Decimal(row['signed_wavelength_angstrom_native'])))
                for basis in ['observed', 'theoretical']:
                    column = basis+'_energy_cm_inverse_native'
                    gap = (Decimal(levels[hi][column])-Decimal(levels[lo][column]))*100*h*c
                    check('source-{}-state-gap-{}'.format(basis, k), close(actual[basis+'_gap_j'], gap))
            for k, actual in enumerate(result['recoil']):
                mi, mf, c = (Decimal(str(actual[x])) for x in ['mi', 'mf', 'c'])
                gap = (mi-mf)*c*c
                photon = gap*(1-gap/(2*mi*c*c))
                p = photon/c
                mass_gamma = (mf*mf+(p/c)**2).sqrt()
                kinetic = p*p/(mass_gamma+mf)
                for name, expected in [('gap', gap), ('photon_j', photon), ('p', p), ('ke', kinetic), ('velocity', -p/mass_gamma)]:
                    check('independent-recoil-{}-{}'.format(k, name), close(actual[name], expected))
                check('independent-energy-conservation-{}'.format(k), abs(photon+kinetic-gap) <= Decimal('1e-80')*gap)
                required_guards = {'missing-isotope', 'zero-direction', 'mass-gap-mismatch',
                                   'constant-mismatch', 'negative-bound'}
                if k == 0:
                    required_guards.add('non-axis-zero-bound-refusal')
                if set(actual['guards']) != required_guards:
                    raise ValueError('native recoil receipt omitted a required control')
                for name, value in actual['guards'].items():
                    check('recoil-{}-{}'.format(k, name), value is True)
        # Independent forward equation over five populations and nine cumulative
        # edge counts. An acyclic downward history crosses any given edge at
        # most once, so every declared count has a binomial marginal.
        matrix = np.zeros((14,14))
        for row in links:
            u, l = int(row['upper_index'])-1, int(row['lower_index'])-1
            a = float(row['decay_rate_s_inverse_native'])
            matrix[u,l] += a
            matrix[u,u] -= a
        for k, row in enumerate(links):
            matrix[int(row['upper_index'])-1, 5+k] = float(row['decay_rate_s_inverse_native'])
        rows = []
        check('frozen-eight-families', len(result['families']) == 8)
        for family, actual in enumerate(result['families']):
            expected_initial = 4+family//4
            expected_time = [1,100,1e5,1e10][family%4]
            check('frozen-family-{}'.format(family), actual['initial'] == expected_initial and actual['duration_s'] == expected_time and actual['seed'] == 1846053+family and actual['n'] == N)
            initial = np.zeros(14); initial[expected_initial-1] = 1
            probabilities = initial.dot(expm(matrix*expected_time))
            check('matrix-population-normalization-{}'.format(family), abs(sum(probabilities[:5])-1) <= RELATIVE)
            if len(actual['population']) != 5 or len(actual['photons']) != 9:
                raise ValueError('native history count vector is incomplete')
            check('whole-cohort-population-{}'.format(family), sum(actual['population']) == N)
            counts = actual['population']+actual['photons']
            for j, (p, observed) in enumerate(zip(probabilities, counts)):
                if not math.isfinite(p) or p < -RELATIVE or p > 1+RELATIVE:
                    raise ValueError('independent matrix oracle exceeded declared arithmetic envelope')
                # Only the fixed numerical envelope handles matrix-rounding
                # endpoints. No source data or physical probabilities are fitted.
                p = max(0., min(1., float(p)))
                low = int(binom.ppf(ALPHA_PER_CLAIM/2, N, p))
                high = int(binom.ppf(1-ALPHA_PER_CLAIM/2, N, p))
                name = ('population-level-'+str(j+1)) if j < 5 else ('photon-link-'+str(j-5))
                check('family-{}-{}'.format(family,name), type(observed) is int and low <= observed <= high,
                      count=observed, lower=low, upper=high, independent_probability=p)
                rows.append(dict(family=family, initial_level=expected_initial, duration_s=expected_time,
                                 quantity=name, native_count=observed, independent_expected_count=N*p,
                                 lower_count=low, upper_count=high))
        with (output/'comparison.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        fig, axes = plt.subplots(2,1,figsize=(13,8))
        for ax, group in zip(axes, ['population','photon']):
            selected = [r for r in rows if r['quantity'].startswith(group)]
            x = range(len(selected))
            ax.plot(x, [r['native_count']/N for r in selected], '.', label='production supplied-uniform histories')
            ax.plot(x, [r['independent_expected_count']/N for r in selected], 'x', label='independent matrix exponential')
            ax.set_ylabel(group+' probability / per-history count')
            ax.set_xlabel('predeclared family/quantity index');ax.legend(loc='best')
        fig.suptitle('O I selected radiative-network mathematics — no observations or lunar admission')
        fig.tight_layout(rect=[0,0,1,.96])
        for suffix in ['png','eps']:
            fig.savefig(str(output/('comparison.'+suffix)), dpi=150)
        plt.close(fig)
        status = 'FAIL' if any(x['status']=='FAIL' for x in checks) else 'PASS'
        reason = 'evaluated local source/network/recoil contract'
    except Missing as e:
        status, reason = 'SKIPPED', str(e)
    except Exception as e:
        status, reason = 'ERROR', '{}: {}'.format(type(e).__name__, e)
    counts = {key:sum(x['status']==key for x in checks) for key in ['PASS','FAIL','SKIPPED','ERROR']}
    if status in ['SKIPPED','ERROR']:
        counts[status] += 1
    document = dict(schema='srcMoon-oxygen-cascade-local-v1', status=status, reason=reason,
                    formal_status='SKIPPED', physical_activation='DISABLED',
                    scope='source and mathematical network only; no admitted isotope/quenching/omission profile',
                    counts=counts, checks=checks, elapsed_s=time.monotonic()-started,
                    protocol=dict(N=N, seed=1846053, families=8, claims=112, familywise_alpha=.01,
                                  numeric_relative=RELATIVE), artifacts=[str(p) for p in sorted(output.iterdir())])
    (output/'result.json').write_text(json.dumps(document, indent=2, allow_nan=False)+'\n')
    (output/'run_manifest.json').write_text(json.dumps(dict(schema='srcMoon-oxygen-cascade-manifest-v1',
         commands=commands, identities=identities, tool_versions=locals().get('versions', {}), database_root=str(database), MPI_ranks=1,
         OpenMP_threads=1, source_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=str(REPO),text=True).strip()),indent=2)+'\n')
    print('{} oxygen cascade: {} {}'.format(status, counts, output/'result.json'))
    return dict(PASS=0, FAIL=1, ERROR=2, SKIPPED=77)[status]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--database-root', type=Path, default=DEFAULT_DATABASE)
    parser.add_argument('--wall-seconds', type=float, default=120)
    args = parser.parse_args()
    if not math.isfinite(args.wall_seconds) or not 0 < args.wall_seconds <= 3600:
        parser.error('--wall-seconds must be finite in (0,3600]')
    return run(args.output_dir.resolve(), args.database_root.resolve(), args.wall_seconds)

if __name__ == '__main__':
    sys.exit(main())
