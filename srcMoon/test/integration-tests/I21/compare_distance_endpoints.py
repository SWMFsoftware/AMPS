#!/usr/bin/env python3
"""I21 decision 1: the 1 AU / 2 AU native photoionization rate ratio (approved 2026-11-26).

Reads the actual native event metrics of the four prescribed-Sun cases
(sun-1au-quiet, sun-2au-quiet, sun-1au-active, sun-2au-active) as run by the
I21 owner. Each native event rate k(x_event) is mapped to the start point by
the exact inverse-square factor |S - x_event|^2 / |S - x_start|^2, using the
prescribed Sun and start recorded in the native receipt. For each solar
endpoint, mover and launch the ratio k_start(1 AU) / k_start(2 AU) must equal
4 within the existing I21 inverse-square tolerance of 0.1%. The active/quiet
ratio at each distance is also reported against the D04 table. No AMPS output
is used as a reference; the expected values are 4 and the D04 rates.

Usage: compare_distance_endpoints.py --campaign-dir DIR --output-dir NEW_DIR
where DIR holds <case>/run-<mover>/I21/trajectory/native-NN/ from the case driver.
"""
import argparse
import csv
from decimal import Decimal, localcontext
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
TEST = HERE.parents[1]
sys.path.insert(0, str(TEST / 'sodium_cohort'))
sys.path.insert(0, str(TEST))
import verify_ledger as ledger  # noqa: E402

LIMIT = Decimal('0.001')          # existing I21 inverse-square relative criterion
EXIT = {'PASS': 0, 'FAIL': 1, 'ERROR': 2, 'SKIPPED': 77}


def start_rates(native):
    """Native event rates mapped to the start point; one value per event."""
    receipt = json.loads((native / 'input_receipt.json').read_text())
    if receipt.get('photon_source_geometry') != 'prescribed-static-point-sun':
        raise ValueError('not a prescribed-Sun launch: ' + str(native))
    sun = [Decimal(v) for v in receipt['sun_position_lso_m']]
    start = [Decimal(v) for v in receipt['source_position_lso_m']]
    ranks = int(receipt.get('mpi_ranks', 0)) or len(list(native.glob('event_metrics.rank*.hex.csv')))
    events = ledger.read_events([native / f'event_metrics.rank{r}.hex.csv' for r in range(ranks)])
    rates = []
    with localcontext() as context:
        context.prec = 60
        d0 = sum((s - x) ** 2 for s, x in zip(sun, start))
        for event in events.values():
            rate = event['source_rate_event_s']
            position = [event[f'parent_event_x{i}'] for i in range(3)]
            d = sum((s - Decimal(x.numerator) / Decimal(x.denominator)) ** 2 for s, x in zip(sun, position))
            rates.append(Decimal(rate.numerator) / Decimal(rate.denominator) * d / d0)
    return receipt, rates


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--campaign-dir', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    checks, rows = [], []
    result = {'status': 'ERROR', 'checks': checks,
              'scope': 'I21 decision-1 1 AU / 2 AU native rate ratio; supporting only',
              'formal_I21_status': 'SKIPPED / NOT QUALIFIED'}
    try:
        table = {}
        for endpoint in ('quiet', 'active'):
            for distance in (1, 2):
                case = args.campaign_dir / f'sun-{distance}au-{endpoint}'
                for mover in ('first', 'second'):
                    trajectory = case / f'run-{mover}/I21/trajectory'
                    natives = sorted(trajectory.glob('native-*/native-output'))
                    if len(natives) != 6:
                        raise FileNotFoundError(f'expected six launches in {trajectory}')
                    for native in natives:
                        receipt, rates = start_rates(native)
                        if not rates:
                            raise ValueError('no native events in ' + str(native))
                        spread = (max(rates) - min(rates)) / min(rates)
                        table[(endpoint, distance, mover, native.parent.name)] = (receipt, rates)
                        rows.append({'endpoint': endpoint, 'distance_au': distance, 'mover': mover,
                                     'launch': native.parent.name, 'steps': receipt['physical_steps'],
                                     'events': len(rates), 'start_rate_s': str(sum(rates) / len(rates)),
                                     'relative_spread': str(spread)})
        with localcontext() as context:
            context.prec = 60
            d04 = None
            for endpoint in ('quiet', 'active'):
                for mover in ('first', 'second'):
                    for launch in sorted({k[3] for k in table}):
                        one, two = table[(endpoint, 1, mover, launch)], table[(endpoint, 2, mover, launch)]
                        k1, k2 = sum(one[1]) / len(one[1]), sum(two[1]) / len(two[1])
                        error = abs(k1 / k2 / 4 - 1)
                        checks.append({'name': f'rate ratio 1AU/2AU = 4 ({endpoint}, {mover}, {launch})',
                                       'status': 'PASS' if error <= LIMIT else 'FAIL', 'actual': str(k1 / k2),
                                       'relative_error': str(error), 'limit': str(LIMIT)})
            for distance in (1, 2):
                for mover in ('first', 'second'):
                    for launch in sorted({k[3] for k in table}):
                        quiet, active = table[('quiet', distance, mover, launch)], table[('active', distance, mover, launch)]
                        channels = {c['id']: c for c in quiet[0]['source_channels']}
                        channel, = channels.values()
                        expected = Decimal(channel['active_rate_1au_s']) / Decimal(channel['quiet_rate_1au_s'])
                        ratio = (sum(active[1]) / len(active[1])) / (sum(quiet[1]) / len(quiet[1]))
                        error = abs(ratio / expected - 1)
                        checks.append({'name': f'active/quiet = D04 ratio ({distance} AU, {mover}, {launch})',
                                       'status': 'PASS' if error <= LIMIT else 'FAIL', 'actual': str(ratio),
                                       'reference': str(expected), 'relative_error': str(error), 'limit': str(LIMIT)})
        with (output / 'rates.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader(); writer.writerows(rows)
        failures = [c['name'] for c in checks if c['status'] == 'FAIL']
        result.update(status='FAIL' if failures else 'PASS', reason=', '.join(failures) or 'all ratios within 0.1%')
    except FileNotFoundError as error:
        result.update(status='SKIPPED', reason=str(error))
    except (OSError, ValueError, KeyError, ArithmeticError) as error:
        result.update(status='ERROR', reason=str(error))
    result['counts'] = {s: sum(c['status'] == s for c in checks) for s in EXIT}
    (output / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'checks'}, indent=2))
    return EXIT[result['status']]


if __name__ == '__main__':
    sys.exit(main())
