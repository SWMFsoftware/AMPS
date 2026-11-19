#!/usr/bin/env python3
"""Audit actual sodium commissioning artifacts; no refinement/statistical PASS.

The single-level comparison uses the existing approved I21 timing allocation,
but never substitutes for three-level/1-and-8-rank convergence or a fixed-N
ensemble experiment. Source rates, input bounds, physical caps and observations
are neither fitted nor changed by this postprocessor.
"""
from __future__ import annotations
import argparse
import csv
from decimal import Decimal, localcontext
import json
from pathlib import Path
import subprocess
import sys
import traceback

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import audit
import verify_ledger as ledger
import verify_numerical as numerical
import verify_p02 as p02
import m3_qualification as q
from m2_initial_state import assignments

EXIT = {'PASS': 0, 'FAIL': 1, 'ERROR': 2, 'SKIPPED': 77}
BINDING_SCHEMA = 'srcMoon-M3-P02-owner-binding-v1'


def first_hit_comparisons(parents, rate, horizon, seed, cohort):
    """Independent HMAC/exponential oracle; no production clock is imported.

    Inputs are already validated original-parent rows and explicit SI scalars.
    Exact endpoint membership is event_time <= horizon, including a first hit
    at the endpoint. The unchanged 0.002 error allocation is supporting only;
    no convergence or ensemble claim is made by this function.
    """
    if not isinstance(rate, Decimal) or not rate.is_finite() or rate < 0:
        raise ValueError('finite nonnegative independent Decimal source rate required')
    if not isinstance(horizon, Decimal) or not horizon.is_finite() or horizon <= 0:
        raise ValueError('finite positive Decimal physical horizon required')
    limit = Decimal(str(q.NUMERICAL_LIMIT['I21']))
    comparisons = []
    with localcontext() as context:
        context.prec = 80
        for row in parents:
            parent = int(row['original_parent_id'])
            uniform = Decimal.from_float(numerical.source_uniform(seed, cohort, parent))
            event = None if rate == 0 else -uniform.ln() / rate
            expected = 'REACTED' if event is not None and event <= horizon else 'SURVIVED'
            same = row['status'] == expected
            error = None
            if same and expected == 'REACTED':
                actual = Decimal(row['relative_event_time_s'])
                if not actual.is_finite():
                    raise ValueError('nonfinite actual first-removal time')
                error = abs(actual - event) / horizon
            comparisons.append({'original_parent_id': parent,
                'expected_status': expected, 'actual_status': row['status'],
                'source_uniform': str(uniform),
                'reference_first_event_s': '' if event is None else str(event),
                'actual_first_event_s': row['relative_event_time_s'],
                'normalized_timing_error': '' if error is None else str(error),
                'timing_limit': str(limit),
                'membership_status': 'PASS' if same else 'FAIL',
                'single_level_timing_status': 'FAIL' if not same or
                    (error is not None and error > limit) else 'PASS'})
    return comparisons


def read_parents(path, receipt, channel):
    """Validate the existing native CSV WITHOUT manufacturing a formal plan.

    The qualification reader requires a statistically eligible nontrivial
    survival and a frozen family. A commissioning parser must also accept a
    valid small/zero-event cohort, so only its unchanged native columns and
    finite-number helpers are reused here. No chemistry law is implemented.
    """
    with Path(path).open(newline='') as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != q.PARENT_COLUMNS:
            raise ValueError('original-parent CSV schema mismatch')
        rows = list(reader)
    count = receipt['original_parents']
    if type(count) is not int or count <= 0 or len(rows) != count:
        raise ValueError('complete positive original-parent population required')
    horizon = q.finite(receipt['physical_horizon_s'], 'actual physical horizon', True)
    weight = q.finite(receipt['initial_represented_number'], 'actual represented number', True)
    seen = set()
    for row in rows:
        if None in row or any(v is None for v in row.values()): raise ValueError('ragged original-parent CSV')
        key = row['original_parent_id']
        if not key or not key.isascii() or not key.isdigit() or str(int(key)) != key or key in seen:
            raise ValueError('unique canonical original-parent IDs required')
        seen.add(key)
        if q.finite(row['initial_represented_number'], 'parent represented number', True) != weight:
            raise ValueError('parent represented number differs from actual receipt')
        if row['status'] == 'SURVIVED':
            if row['relative_event_time_s'] or row['channel_id']:
                raise ValueError('survivor cannot carry a removal event')
        elif row['status'] == 'REACTED':
            when = q.finite(row['relative_event_time_s'], 'first-removal time')
            if when < 0 or when > horizon or row['channel_id'] != channel:
                raise ValueError('actual removal outside declared horizon/source channel')
        else: raise ValueError('commissioning cannot silently censor an escaped/surface/unknown parent')
    if {int(key) for key in seen} != set(range(count)):
        raise ValueError('complete contiguous original-parent ID range required')
    return rows


def _complete_inputs(manifest, receipt, identities):
    """Bind the receipt's COMPLETE graph to the owner's ORIGINAL artifact map.

    The wrapper copied root input bytes before execution, but includes remain
    original source documents. Parse both owners with the maintained readers,
    including empty fragments, and require exact graph equality rather than
    accepting arbitrary entries supplied after the native run.
    """
    root = Path(manifest['input']['path']).resolve()
    if ledger.digest(root) != manifest['input']['sha256']:
        raise ValueError('original root input changed after completed owner')
    controls, native_docs = p02.native_controls(root)
    record = controls.get('application_input')
    if record is None:
        raise ValueError('actual native application input missing')
    raw = record['value']
    if len(raw) >= 2 and raw[0] == raw[-1] == '"': raw = raw[1:-1]
    application = Path(raw)
    application = (application if application.is_absolute() else record['file'].parent / application).resolve()
    # Application parser validates owner/duplicates/includes. The existing
    # preparation inventory preserves even includes without assignments.
    from prepare_physical_case import application_documents
    _, application_docs = application_documents(application)
    expected = {str(p.resolve()) for p in native_docs} | {r['path'] for r in application_docs}
    actual = audit.input_documents(receipt['input_documents'])
    if set(actual) != expected:
        raise ValueError('native receipt omits/adds an original input document')
    for path, digest in actual.items():
        if ledger.digest(path) != digest:
            raise ValueError('original native input document changed: ' + path)
        identities[path] = digest
    for key in ('profile_metadata', 'profile_native_manifest', 'isotopes',
                'constants', 'states', 'crosswalk', 'transport_policy_evidence'):
        file_record = controls.get(key + '_file'); digest_record = controls.get(key + '_sha256')
        if file_record is None or digest_record is None:
            raise ValueError('pinned native source descriptor absent: ' + key)
        name = file_record['value']
        if len(name) >= 2 and name[0] == name[-1] == '"': name = name[1:-1]
        source = Path(name)
        source = (source if source.is_absolute() else file_record['file'].parent / source).resolve()
        ledger.authenticate({'path': str(source), 'sha256': digest_record['value']}, root.parent, identities)
    return controls


def _write_comparisons(output, comparisons):
    if not comparisons: return
    with (output / 'comparison.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(comparisons[0]))
        writer.writeheader(); writer.writerows(comparisons)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(1, 2, figsize=(11, 4))
    for owner in sorted({r['owner'] for r in comparisons}):
        rows = [r for r in comparisons if r['owner'] == owner and r['reference_first_event_s']]
        reacted = [r for r in rows if r['actual_first_event_s']]
        axes[0].scatter([float(r['reference_first_event_s']) for r in reacted],
                        [float(r['actual_first_event_s']) for r in reacted], s=9, label=owner)
        timed = [r for r in rows if r['normalized_timing_error']]
        axes[1].scatter([r['original_parent_id'] for r in timed],
                        [float(r['normalized_timing_error']) for r in timed], s=9, label=owner)
    axes[0].set_xlabel('Independent HMAC/exponential first hit [s]')
    axes[0].set_ylabel('Actual native first hit [s]')
    axes[1].axhline(q.NUMERICAL_LIMIT['I21'], color='black', linestyle='--', label='existing I21 0.002 allocation')
    axes[1].set_xlabel('Original parent ID'); axes[1].set_ylabel('Single-level |timing error| / horizon')
    for axis in axes: axis.legend(framealpha=1)
    figure.suptitle('Commissioning support only — no refinement or ensemble qualification')
    figure.tight_layout()
    for suffix in ('png', 'eps'): figure.savefig(output / ('comparison.' + suffix))
    plt.close(figure)


def run(binding_path, output):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    identities = {}; checks = []; comparisons = []; diagnostics = []
    manifest = {'schema': 'srcMoon-M3-commissioning-audit-manifest-v1',
        'scope': 'completed-owner single-level independent commissioning audit only'}
    def check(name, status, **detail): checks.append({'name': name, 'status': status, **detail})
    try:
        binding_path = Path(binding_path).resolve()
        if not binding_path.is_file(): raise ledger.Unavailable('explicit commissioning binding unavailable')
        identities[str(binding_path)] = ledger.digest(binding_path)
        binding = ledger.strict_json(binding_path)
        if set(binding) != {'schema', 'profile', 'budget', 'owners'} or binding['schema'] != BINDING_SCHEMA:
            raise ValueError('exact P02 owner-binding schema required')
        # Bind the budget bytes, without interpreting/weakening caps here. The
        # separate verify_p02 owner performs those actual mechanical checks.
        try:
            ledger.authenticate(binding['budget'], binding_path.parent, identities)
            profile = ledger.authenticate(binding['profile'], binding_path.parent, identities)
        except FileNotFoundError as error:
            raise ledger.Unavailable('required source profile/budget file unavailable') from error
        owners = binding['owners']
        if not isinstance(owners, list): raise ValueError('explicit owners list required')
        if not owners: raise ledger.Unavailable('no completed native owner supplied')
        seen = set()
        for index, owned in enumerate(owners):
            if not isinstance(owned, dict) or set(owned) != {'owner_manifest', 'owner_result', 'mpi_ranks'} \
                    or type(owned['mpi_ranks']) is not int or owned['mpi_ranks'] not in (1, 8):
                raise ValueError('completed owner descriptors and actual 1/8 MPI layout required')
            try:
                mp = ledger.authenticate(owned['owner_manifest'], binding_path.parent, identities)
                rp = ledger.authenticate(owned['owner_result'], binding_path.parent, identities)
            except FileNotFoundError as error:
                raise ledger.Unavailable('completed native owner manifest/result unavailable') from error
            if mp in seen: raise ValueError('duplicate owner is not additional independent evidence')
            seen.add(mp); owner = ledger.strict_json(mp)
            launch = {'owner_manifest': {'path': str(mp), 'sha256': identities[str(mp)]},
                'owner_result': {'path': str(rp), 'sha256': identities[str(rp)]},
                'mpi_ranks': owned['mpi_ranks'], 'input': owner['input'],
                'output': str(mp.parent / 'native-output')}
            native = audit.authenticate_launch(launch)
            identities[str(native / 'input_receipt.json')] = ledger.digest(native / 'input_receipt.json')
            identities[str(native / 'parent_records.csv')] = ledger.digest(native / 'parent_records.csv')
            receipt = ledger.strict_json(native / 'input_receipt.json')
            if receipt.get('phase') != 'explicit-experimental-physical-input':
                raise ledger.Unavailable('completed physical native input receipt unavailable')
            if receipt.get('case_kind') != 'homogeneous-fixed-position':
                raise ValueError('variable hazard cannot use the homogeneous commissioning reference')
            _complete_inputs(owner, receipt, identities)
            mapping = receipt['mapping']
            if mapping.get('profile_sha256') != binding['profile']['sha256'] \
                    or receipt['profile']['sha256'] != binding['profile']['sha256']:
                raise ValueError('native profile differs from authenticated binding')
            reference = numerical.reference_receipt(receipt)
            rate, channel = numerical.source_rate(reference, identities)
            parents = read_parents(native / 'parent_records.csv', receipt, channel)
            paths = [native / ('event_metrics.rank%d.hex.csv' % rank) for rank in range(owned['mpi_ranks'])]
            events = ledger.read_events(paths)
            identities.update({str(path): ledger.digest(path) for path in paths})
            reacted = {row['original_parent_id'] for row in parents if row['status'] == 'REACTED'}
            if set(events) != reacted: raise ValueError('native parent and exact per-rank event sets differ')
            check('owner-%02d-original-input-parent-event-identity' % index, 'PASS')
            links = numerical.event_time_links(events, parents, receipt['transport_origin_et_s'])
            (output / ('event_time_links-%02d.json' % index)).write_text(json.dumps(links, indent=2) + '\n')
            check('owner-%02d-relative-and-absolute-clock-links' % index, links['status'])
            inverse = numerical.inverse_square_checks(events, rate)
            (output / ('inverse_square-%02d.json' % index)).write_text(json.dumps(inverse, indent=2) + '\n')
            check('owner-%02d-source-and-accepted-inverse-square' % index,
                inverse['status'] if events else 'SKIPPED',
                reason='accepted-event rate comparison' if events else 'no accepted event; empty-set identity is not rate verification')
            base = numerical.D04 / 'processed/phidrates'
            mapping = {**mapping, 'photon_source': {'reaction_id': channel,
                'solar_activity': '0' if reference['solar_endpoint'] == 'quiet' else '1',
                **{key: {'path': str(base / name), 'sha256': q.sha(base / name)}
                   for key, name in [('rates', 'photo_rates_1au.csv'), ('thresholds', 'thresholds.csv')]}}}
            bounds = {**receipt['bounds'], 'schema': 'srcMoon-M3-event-ledger-bounds-v1',
                      'scope': 'source-bound-experimental', 'event_frame': mapping['event_frame']}
            reconstructed = ledger.verify(paths, profile, mapping, bounds)
            identities.update(reconstructed.get('input_sha256', {}))
            check('owner-%02d-independent-committed-event-ledger' % index, reconstructed['status'])
            try:
                ledger.write_artifacts(reconstructed, output / ('ledger-%02d' % index))
            except ImportError as error:
                check('owner-%02d-ledger-visualization-dependency' % index, 'SKIPPED', reason=str(error))
            horizon = Decimal.from_float(float(receipt['physical_horizon_s']))
            rows = first_hit_comparisons(parents, rate, horizon, receipt['seed'], receipt['cohort_id'])
            for row in rows: row['owner'] = 'owner-%02d' % index
            comparisons.extend(rows)
            timing_status = 'PASS' if all(r['single_level_timing_status'] == 'PASS' for r in rows) else 'FAIL'
            check('owner-%02d-single-level-first-hit-error-allocation' % index, timing_status,
                relative_limit=str(q.NUMERICAL_LIMIT['I21']), refinement_status='NOT_EVALUATED')
            survivors = sum(r['status'] == 'SURVIVED' for r in parents)
            with localcontext() as context:
                context.prec = 80
                survival = (-rate * horizon).exp()
                diagnostics.append({'owner': 'owner-%02d' % index, 'original_parents': len(parents),
                    'survivors': survivors, 'removals': len(parents) - survivors,
                    'actual_survival_fraction': str(Decimal(survivors) / Decimal(len(parents))),
                    'independent_survival_probability': str(survival),
                    'independent_expected_survivors': str(survival * len(parents)),
                    'independent_rate_s': str(rate), 'horizon_s': str(horizon),
                    'statistical_qualification_status': 'NOT_EVALUATED',
                    'scope': 'counts/reference diagnostics only; no lifetime/branch acceptance'})
        (output / 'survival_diagnostics.json').write_text(json.dumps(diagnostics, indent=2) + '\n')
        _write_comparisons(output, comparisons)
    except (q.MissingPrerequisite, ledger.Unavailable) as error:
        check('required-evidence', 'SKIPPED', reason=str(error))
    except ImportError as error:
        check('visualization-dependency', 'SKIPPED', reason=str(error))
    except Exception as error:
        (output / 'exception.log').write_text(traceback.format_exc())
        check('owner-source-infrastructure', 'ERROR', reason=str(error))
    counts = {s: sum(c['status'] == s for c in checks) for s in EXIT}
    status = next(s for s in ('ERROR', 'FAIL', 'SKIPPED', 'PASS') if counts[s])
    result = {'schema': 'srcMoon-M3-independent-commissioning-audit-v1', 'status': status,
        'counts': counts, 'checks': checks, 'owners_evaluated': len(diagnostics),
        'parents_evaluated': len(comparisons), 'events_evaluated': sum(d['removals'] for d in diagnostics),
        'scope': 'single-level source-bound commissioning support only',
        'formal_M3_status': 'SKIPPED / NOT QUALIFIED', 'refinement_status': 'NOT_EVALUATED',
        'statistical_qualification_status': 'NOT_EVALUATED',
        'P02_cap_audit_status': 'SEPARATE_OWNER_NOT_EVALUATED_HERE', 'observational_scoring': False}
    manifest['input_sha256'] = identities
    manifest['source_revision'] = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=HERE.parents[2], text=True).strip()
    manifest['verifier_sha256'] = {str(p): q.sha(p) for p in
        (Path(__file__).resolve(), Path(audit.__file__), Path(ledger.__file__),
         Path(numerical.__file__), Path(p02.__file__), Path(q.__file__), HERE / 'prepare_physical_case.py')}
    (output / 'audit.log').write_text('\n'.join(c['status'] + ' ' + c['name'] +
        (': ' + c['reason'] if 'reason' in c else '') for c in checks) + '\n')
    (output / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    (output / 'run_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binding', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    try:
        result = run(args.binding, args.output_dir)
    except OSError as error:
        print(json.dumps({'status': 'ERROR', 'reason': str(error)})); raise SystemExit(2)
    print(json.dumps(result, indent=2)); raise SystemExit(EXIT[result['status']])
