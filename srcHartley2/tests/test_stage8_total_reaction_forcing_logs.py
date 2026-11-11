#!/usr/bin/env python3
"""Strictly validate Stage-8 mover-integrated forced-reaction logs.

Stage 8 is accepted only as a deterministic *two-branch* segment algorithm:
each incoming carrier produces an endpoint survivor with exp(-tau) of the
carrier and a reacted branch with 1-exp(-tau).  The reacted branch is placed
at a physical conditional reaction point and its products are advanced over
the remaining mover interval.  This parser independently recomputes the
branch, exponential-decay, piecewise-hazard, channel/packet, and transaction
identities from the durable log instead of trusting the producer's PASS token.

The runtime validator uses global deterministic event IDs.  Consequently the
homogeneous, position, decay, and composition records must agree at 1, 4, and
12 ranks.  Only the deliberately one-transaction-per-rank capacity counters
scale with MPI decomposition.
"""

from __future__ import annotations

import argparse
import math
import pathlib
import re
import sys


FLOAT = r"[+\-0-9.eE]+"
HEADER = re.compile(
    r"^Total-reaction forcing validation: events/fixture=(\d+) ranks=(\d+) "
    r"deterministic-global-events=on branch-mode=deterministic-two-branch$"
)
HOMOGENEOUS = re.compile(
    rf"^(PASS|FAIL) total-reaction-forcing homogeneous parent=(\d+) "
    rf"events=(\d+) steps=(\d+) hazard=({FLOAT}) duration=({FLOAT}) "
    rf"optical-depth=({FLOAT}) analytic-survival=({FLOAT}) "
    rf"represented-survival=({FLOAT}) analytic-reacted=({FLOAT}) "
    rf"represented-reacted=({FLOAT}) closure-error=({FLOAT})$"
)
POSITION = re.compile(
    rf"^(PASS|FAIL) total-reaction-forcing position "
    rf"fixture=(constant|piecewise) parent=(\d+) samples=(\d+) "
    rf"tau-first=({FLOAT}) tau-second=({FLOAT}) first-fraction=({FLOAT}) "
    rf"first-region-count=(\d+) first-region-expected=({FLOAT}) "
    rf"mean-fraction=({FLOAT}) expected-mean-fraction=({FLOAT}) "
    rf"cdf-distance=({FLOAT}) endpoint-violations=(\d+)$"
)
DECAY = re.compile(
    rf"^(PASS|FAIL) total-reaction-forcing decay parent=(\d+) events=(\d+) "
    rf"steps=(\d+) total-optical-depth=({FLOAT}) "
    rf"analytic-survival=({FLOAT}) represented-survival=({FLOAT}) "
    rf"analytic-reacted=({FLOAT}) represented-reacted=({FLOAT}) "
    rf"closure-error=({FLOAT})$"
)
COMPOSITION = re.compile(
    rf"^(PASS|FAIL) total-reaction-forcing composition reaction=(\S+) "
    rf"parent=(\d+) packets=(\d+) products=(\d+) "
    rf"occurrence-factor=({FLOAT}) channel-likelihood=({FLOAT}) "
    rf"expected-correction=({FLOAT}) packet-correction-sum=({FLOAT}) "
    rf"emitted=(\d+) momentum-error=({FLOAT}) energy-error=({FLOAT})$"
)
TRANSACTION = re.compile(
    rf"^(PASS|FAIL) total-reaction-forcing production-transaction ranks=(\d+) "
    rf"parents=(\d+) survivors=(\d+) reacted-branches=(\d+) emitted=(\d+) "
    rf"processed=(\d+) remaining-moves=(\d+) capacity-rollbacks=(\d+) "
    rf"byte-rollbacks=(\d+) rng-rollbacks=(\d+) path-rollbacks=(\d+) "
    rf"double-chemistry=(\d+) correction-error=({FLOAT}) "
    rf"position-error=({FLOAT})$"
)
MOVER_PRODUCTION = re.compile(
    rf"^(PASS|FAIL) total-reaction-forcing mover-production ranks=(\d+) "
    rf"mpi-face=(\d+) face-transfers=(\d+) initial-pending=(\d+) "
    rf"continuation-moves=(\d+) sweeps=(\d+) final-pending=(\d+) "
    rf"product-events=(\d+) sink-events=(\d+) wrong-owner-events=(\d+) "
    rf"final-particles=(\d+) survivors=(\d+) products=(\d+) "
    rf"uncleared-markers=(\d+) nonzero-remainders=(\d+) processed=(\d+) "
    rf"endpoint-checks=(\d+) moved-endpoints=(\d+) "
    rf"correction-error=({FLOAT}) event-position-error=({FLOAT}) "
    rf"remaining-fraction-error=({FLOAT}) endpoint-error=({FLOAT})$"
)
SUMMARY = re.compile(
    r"^Total-reaction forcing validation result: (PASS|FAIL) homogeneous=(\d+) "
    r"piecewise=(\d+) decay=(\d+) composition=(\d+) transactions=(\d+) "
    r"mover=(\d+) events=(\d+)$"
)


def finite(text: str, context: str) -> float:
    value = float(text)
    if not math.isfinite(value):
        raise ValueError(f"{context}: non-finite value {text!r}")
    return value


def close(left: float, right: float, relative: float = 3.0e-12) -> bool:
    return math.isclose(left, right, rel_tol=relative, abs_tol=1.0e-300)


def constant_expected_fraction(optical_depth: float) -> float:
    if optical_depth < 1.0e-4:
        return 0.5 - optical_depth / 12.0 + optical_depth**3 / 720.0
    return 1.0 / optical_depth - 1.0 / math.expm1(optical_depth)


def piecewise_expected_fraction(first: float, tau1: float, tau2: float) -> float:
    second = 1.0 - first
    rate1 = tau1 / first
    rate2 = tau2 / second
    e1 = math.exp(-tau1)
    e2 = math.exp(-tau2)
    probability = -math.expm1(-(tau1 + tau2))
    moment1 = 1.0 / rate1 - (first + 1.0 / rate1) * e1
    moment2 = e1 * (
        first * (1.0 - e2)
        + 1.0 / rate2
        - (second + 1.0 / rate2) * e2
    )
    return (moment1 + moment2) / probability


def binomial_allowed(samples: int, probability: float) -> float:
    sigma = math.sqrt(samples * probability * (1.0 - probability))
    roundoff = 4096.0 * sys.float_info.epsilon * max(
        1.0, samples * probability
    )
    return max(8.0 * sigma, roundoff)


def parse_log(path: pathlib.Path, ranks: int, expected_events: int) -> dict:
    text = path.read_text(encoding="utf-8", errors="strict")
    if any(line.startswith("FAIL total-reaction-forcing") for line in text.splitlines()):
        raise ValueError(f"{path}: contains a failed Stage-8 record")

    headers = []
    homogeneous = []
    positions: dict[str, dict] = {}
    decays = []
    compositions = []
    transactions = []
    mover_production = []
    summaries = []

    for line in text.splitlines():
        match = HEADER.match(line)
        if match:
            headers.append((int(match.group(1)), int(match.group(2))))
            continue

        match = HOMOGENEOUS.match(line)
        if match:
            (
                status,
                parent,
                events,
                steps,
                hazard,
                duration,
                optical_depth,
                analytic_survival,
                represented_survival,
                analytic_reacted,
                represented_reacted,
                closure_error,
            ) = match.groups()
            row = {
                "parent": int(parent),
                "events": int(events),
                "steps": int(steps),
                "hazard": finite(hazard, f"{path}:homogeneous:hazard"),
                "duration": finite(duration, f"{path}:homogeneous:duration"),
                "optical_depth": finite(
                    optical_depth, f"{path}:homogeneous:optical-depth"
                ),
                "analytic_survival": finite(
                    analytic_survival, f"{path}:homogeneous:analytic-survival"
                ),
                "represented_survival": finite(
                    represented_survival,
                    f"{path}:homogeneous:represented-survival",
                ),
                "analytic_reacted": finite(
                    analytic_reacted, f"{path}:homogeneous:analytic-reacted"
                ),
                "represented_reacted": finite(
                    represented_reacted,
                    f"{path}:homogeneous:represented-reacted",
                ),
                "closure_error": finite(
                    closure_error, f"{path}:homogeneous:closure-error"
                ),
            }
            if status != "PASS" or row["events"] != expected_events:
                raise ValueError(f"{path}: invalid homogeneous status/events")
            if row["parent"] < 0 or row["steps"] <= 0:
                raise ValueError(f"{path}: invalid homogeneous parent/steps")
            if row["hazard"] <= 0.0 or row["duration"] <= 0.0:
                raise ValueError(f"{path}: invalid homogeneous hazard/duration")
            expected_tau = row["hazard"] * row["duration"] * row["steps"]
            survival = math.exp(-expected_tau)
            reacted = -math.expm1(-expected_tau)
            if not close(row["optical_depth"], expected_tau):
                raise ValueError(f"{path}: homogeneous optical depth is wrong")
            for key, expected in (
                ("analytic_survival", survival),
                ("represented_survival", survival),
                ("analytic_reacted", reacted),
                ("represented_reacted", reacted),
            ):
                if not close(row[key], expected):
                    raise ValueError(f"{path}: homogeneous {key} is wrong")
            if row["closure_error"] < 0.0 or row["closure_error"] > 2.0e-12:
                raise ValueError(f"{path}: homogeneous carrier closure fails")
            homogeneous.append(row)
            continue

        match = POSITION.match(line)
        if match:
            (
                status,
                fixture,
                parent,
                samples,
                tau1,
                tau2,
                first_fraction,
                first_count,
                first_expected,
                mean_fraction,
                expected_mean,
                cdf_distance,
                endpoint_violations,
            ) = match.groups()
            if fixture in positions:
                raise ValueError(f"{path}: duplicate position fixture {fixture}")
            row = {
                "fixture": fixture,
                "parent": int(parent),
                "samples": int(samples),
                "tau1": finite(tau1, f"{path}:{fixture}:tau-first"),
                "tau2": finite(tau2, f"{path}:{fixture}:tau-second"),
                "first_fraction": finite(
                    first_fraction, f"{path}:{fixture}:first-fraction"
                ),
                "first_count": int(first_count),
                "first_expected": finite(
                    first_expected, f"{path}:{fixture}:first-region-expected"
                ),
                "mean_fraction": finite(
                    mean_fraction, f"{path}:{fixture}:mean-fraction"
                ),
                "expected_mean": finite(
                    expected_mean, f"{path}:{fixture}:expected-mean"
                ),
                "cdf_distance": finite(
                    cdf_distance, f"{path}:{fixture}:cdf-distance"
                ),
                "endpoint_violations": int(endpoint_violations),
            }
            if status != "PASS" or row["samples"] != expected_events:
                raise ValueError(f"{path}: invalid {fixture} position status/samples")
            if row["parent"] < 0 or row["endpoint_violations"] != 0:
                raise ValueError(f"{path}: invalid {fixture} position bounds")
            if fixture == "constant":
                if row["tau1"] <= 0.0 or row["tau2"] != 0.0:
                    raise ValueError(f"{path}: invalid constant optical depth")
                if row["first_fraction"] != 1.0:
                    raise ValueError(f"{path}: constant fixture does not span path")
                region_probability = 1.0
                expected = constant_expected_fraction(row["tau1"])
            else:
                if row["tau1"] <= 0.0 or row["tau2"] <= 0.0:
                    raise ValueError(f"{path}: invalid piecewise optical depth")
                if not 0.0 < row["first_fraction"] < 1.0:
                    raise ValueError(f"{path}: invalid piecewise split")
                region_probability = -math.expm1(-row["tau1"]) / -math.expm1(
                    -(row["tau1"] + row["tau2"])
                )
                expected = piecewise_expected_fraction(
                    row["first_fraction"], row["tau1"], row["tau2"]
                )
            if not close(row["first_expected"], region_probability):
                raise ValueError(f"{path}: {fixture} region probability is wrong")
            if abs(row["first_count"] - row["samples"] * region_probability) > (
                binomial_allowed(row["samples"], region_probability)
            ):
                raise ValueError(f"{path}: {fixture} region count is implausible")
            if not close(row["expected_mean"], expected):
                raise ValueError(f"{path}: {fixture} analytic mean is wrong")
            if not 0.0 <= row["mean_fraction"] <= 1.0:
                raise ValueError(f"{path}: {fixture} sample mean is outside path")
            mean_allowed = 8.0 / math.sqrt(row["samples"])
            if abs(row["mean_fraction"] - expected) > mean_allowed:
                raise ValueError(f"{path}: {fixture} mean is outside uncertainty")
            if not 0.0 <= row["cdf_distance"] <= 8.0 / math.sqrt(row["samples"]):
                raise ValueError(f"{path}: {fixture} conditional CDF fails")
            positions[fixture] = row
            continue

        match = DECAY.match(line)
        if match:
            (
                status,
                parent,
                events,
                steps,
                total_tau,
                analytic_survival,
                represented_survival,
                analytic_reacted,
                represented_reacted,
                closure_error,
            ) = match.groups()
            row = {
                "parent": int(parent),
                "events": int(events),
                "steps": int(steps),
                "total_tau": finite(total_tau, f"{path}:decay:tau"),
                "analytic_survival": finite(
                    analytic_survival, f"{path}:decay:analytic-survival"
                ),
                "represented_survival": finite(
                    represented_survival, f"{path}:decay:represented-survival"
                ),
                "analytic_reacted": finite(
                    analytic_reacted, f"{path}:decay:analytic-reacted"
                ),
                "represented_reacted": finite(
                    represented_reacted, f"{path}:decay:represented-reacted"
                ),
                "closure_error": finite(
                    closure_error, f"{path}:decay:closure-error"
                ),
            }
            if status != "PASS" or row["events"] != expected_events:
                raise ValueError(f"{path}: invalid decay status/events")
            if row["parent"] < 0 or row["steps"] < 2 or row["total_tau"] <= 0.0:
                raise ValueError(f"{path}: invalid decay fixture")
            survival = math.exp(-row["total_tau"])
            reacted = -math.expm1(-row["total_tau"])
            for key, expected in (
                ("analytic_survival", survival),
                ("represented_survival", survival),
                ("analytic_reacted", reacted),
                ("represented_reacted", reacted),
            ):
                if not close(row[key], expected):
                    raise ValueError(f"{path}: decay {key} is wrong")
            if row["closure_error"] < 0.0 or row["closure_error"] > 3.0e-12:
                raise ValueError(f"{path}: decay carrier closure fails")
            decays.append(row)
            continue

        match = COMPOSITION.match(line)
        if match:
            (
                status,
                reaction,
                parent,
                packets,
                products,
                occurrence,
                channel_likelihood,
                expected_correction,
                correction_sum,
                emitted,
                momentum_error,
                energy_error,
            ) = match.groups()
            row = {
                "reaction": reaction,
                "parent": int(parent),
                "packets": int(packets),
                "products": int(products),
                "occurrence": finite(occurrence, f"{path}:composition:P"),
                "channel_likelihood": finite(
                    channel_likelihood, f"{path}:composition:L"
                ),
                "expected_correction": finite(
                    expected_correction, f"{path}:composition:expected"
                ),
                "correction_sum": finite(
                    correction_sum, f"{path}:composition:sum"
                ),
                "emitted": int(emitted),
                "momentum_error": finite(
                    momentum_error, f"{path}:composition:momentum"
                ),
                "energy_error": finite(
                    energy_error, f"{path}:composition:energy"
                ),
            }
            if status != "PASS" or row["reaction"] != "C2_2C":
                raise ValueError(f"{path}: invalid composition status/reaction")
            if row["parent"] < 0 or row["packets"] != 8 or row["products"] != 2:
                raise ValueError(f"{path}: composition did not exercise K=2,N=8")
            if not 0.0 < row["occurrence"] < 1.0 or row["channel_likelihood"] <= 0.0:
                raise ValueError(f"{path}: invalid composition P or L")
            expected = row["occurrence"] * row["channel_likelihood"]
            if not close(row["expected_correction"], expected) or not close(
                row["correction_sum"], expected
            ):
                raise ValueError(f"{path}: Stage-8 P, Stage-7 L, Stage-6 /N do not compose")
            if row["emitted"] != 16 * ranks:
                raise ValueError(f"{path}: composition emitted count is wrong")
            if row["momentum_error"] > 1.0e-11 or row["energy_error"] > 1.0e-11:
                raise ValueError(f"{path}: composition conservation fails")
            compositions.append(row)
            continue

        match = TRANSACTION.match(line)
        if match:
            values = match.groups()
            status = values[0]
            row = {
                "ranks": int(values[1]),
                "parents": int(values[2]),
                "survivors": int(values[3]),
                "reacted": int(values[4]),
                "emitted": int(values[5]),
                "processed": int(values[6]),
                "remaining_moves": int(values[7]),
                "capacity_rollbacks": int(values[8]),
                "byte_rollbacks": int(values[9]),
                "rng_rollbacks": int(values[10]),
                "path_rollbacks": int(values[11]),
                "double_chemistry": int(values[12]),
                "correction_error": finite(
                    values[13], f"{path}:transaction:correction"
                ),
                "position_error": finite(
                    values[14], f"{path}:transaction:position"
                ),
            }
            if status != "PASS" or row["ranks"] != ranks:
                raise ValueError(f"{path}: invalid transaction status/ranks")
            for key, factor in (
                ("parents", 1),
                ("survivors", 1),
                ("reacted", 1),
                ("emitted", 16),
                ("processed", 17),
                ("remaining_moves", 16),
                ("capacity_rollbacks", 1),
                ("byte_rollbacks", 1),
                ("rng_rollbacks", 1),
                ("path_rollbacks", 1),
            ):
                if row[key] != factor * ranks:
                    raise ValueError(f"{path}: transaction {key} count is wrong")
            if row["double_chemistry"] != 0:
                raise ValueError(f"{path}: a Stage-8 branch ran chemistry twice")
            if row["correction_error"] > 2.0e-12 or row["position_error"] > 2.0e-12:
                raise ValueError(f"{path}: transaction correction/position fails")
            transactions.append(row)
            continue

        match = MOVER_PRODUCTION.match(line)
        if match:
            values = match.groups()
            status = values[0]
            integer_names = (
                "ranks",
                "mpi_face",
                "face_transfers",
                "initial_pending",
                "continuation_moves",
                "sweeps",
                "final_pending",
                "product_events",
                "sink_events",
                "wrong_owner_events",
                "final_particles",
                "survivors",
                "products",
                "uncleared_markers",
                "nonzero_remainders",
                "processed",
                "endpoint_checks",
                "moved_endpoints",
            )
            if len(values) != 23 or len(integer_names) != 18:
                raise ValueError(
                    f"{path}: internal mover-production schema mismatch"
                )
            row = {
                name: int(value)
                for name, value in zip(integer_names, values[1:19])
            }
            for name, value in zip(
                (
                    "correction_error",
                    "event_position_error",
                    "remaining_fraction_error",
                    "endpoint_error",
                ),
                values[19:23],
            ):
                row[name] = finite(value, f"{path}:mover-production:{name}")

            if status != "PASS" or row["ranks"] != ranks:
                raise ValueError(f"{path}: invalid mover-production status/ranks")
            expected_integers = {
                # These are exact fixture cardinalities, not merely positive
                # coverage counters.  Requiring the complete tuple prevents a
                # producer from claiming face coverage after moving only an
                # easy subset of the pending family.  The separate product,
                # sink, final-particle, and endpoint fields state the physical
                # outcomes explicitly.
                "face_transfers": 18,
                "initial_pending": 18,
                "continuation_moves": 18,
                "sweeps": 1,
                "final_pending": 0,
                "product_events": 16,
                "sink_events": 1,
                "wrong_owner_events": 0,
                "final_particles": 17,
                "survivors": 1,
                "products": 16,
                "uncleared_markers": 0,
                "nonzero_remainders": 0,
                "processed": 17,
                "endpoint_checks": 16,
                "moved_endpoints": 16,
            }
            for name, expected in expected_integers.items():
                if row[name] != expected:
                    raise ValueError(
                        f"{path}: mover-production {name}={row[name]}, "
                        f"expected {expected}"
                    )
            expected_mpi_face = 0 if ranks == 1 else 1
            if row["mpi_face"] != expected_mpi_face:
                raise ValueError(
                    f"{path}: mover-production MPI-face evidence is missing"
                )
            for name in (
                "correction_error",
                "event_position_error",
                "remaining_fraction_error",
                "endpoint_error",
            ):
                if row[name] < 0.0 or row[name] > 2.0e-12:
                    raise ValueError(
                        f"{path}: mover-production {name} exceeds tolerance"
                    )
            mover_production.append(row)
            continue

        match = SUMMARY.match(line)
        if match:
            summaries.append((match.group(1), *(int(item) for item in match.groups()[1:])))

    if headers != [(expected_events, ranks)]:
        raise ValueError(f"{path}: invalid header(s) {headers}")
    if len(homogeneous) != 1 or set(positions) != {"constant", "piecewise"}:
        raise ValueError(f"{path}: incomplete homogeneous/position coverage")
    if (
        len(decays) != 1
        or len(compositions) != 1
        or len(transactions) != 1
        or len(mover_production) != 1
    ):
        raise ValueError(
            f"{path}: incomplete decay/composition/transaction/mover coverage"
        )
    if summaries != [("PASS", 1, 1, 1, 1, ranks, 1, expected_events)]:
        raise ValueError(f"{path}: invalid summary {summaries}")

    return {
        "homogeneous": homogeneous[0],
        "positions": positions,
        "decay": decays[0],
        "composition": compositions[0],
        "transaction": transactions[0],
        "mover_production": mover_production[0],
    }


def compare_rows(left: dict, right: dict, context: str, ignored: set[str]) -> None:
    if set(left) != set(right):
        raise ValueError(f"{context}: fields differ")
    for key, left_value in left.items():
        if key in ignored:
            continue
        right_value = right[key]
        if isinstance(left_value, float):
            if not close(left_value, right_value):
                raise ValueError(f"{context}: {key} differs by MPI rank count")
        elif left_value != right_value:
            raise ValueError(f"{context}: {key} differs by MPI rank count")


def collect_paths(values: list[list[str]]) -> dict[int, pathlib.Path]:
    result: dict[int, pathlib.Path] = {}
    for rank_text, path_text in values:
        ranks = int(rank_text)
        path = pathlib.Path(path_text).resolve()
        if ranks in result:
            raise ValueError(f"duplicate --forcing-log for {ranks} ranks")
        if not path.is_file():
            raise ValueError(f"missing --forcing-log for {ranks} ranks: {path}")
        result[ranks] = path
    if set(result) != {1, 4, 12}:
        raise ValueError("--forcing-log is required for exactly 1, 4, and 12 ranks")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(
        description="validate Stage-8 deterministic forced-segment 1/4/12 logs"
    )
    parser.add_argument("--expected-events", type=int, required=True)
    parser.add_argument(
        "--forcing-log",
        nargs=2,
        action="append",
        required=True,
        metavar=("RANKS", "PATH"),
    )
    arguments = parser.parse_args()
    if arguments.expected_events < 1000:
        parser.error("--expected-events must be at least 1000")

    try:
        paths = collect_paths(arguments.forcing_log)
        records = {
            ranks: parse_log(path, ranks, arguments.expected_events)
            for ranks, path in sorted(paths.items())
        }
        baseline = records[1]
        for ranks in (4, 12):
            current = records[ranks]
            compare_rows(
                baseline["homogeneous"],
                current["homogeneous"],
                f"{paths[ranks]}:homogeneous",
                set(),
            )
            compare_rows(
                baseline["decay"], current["decay"], f"{paths[ranks]}:decay", set()
            )
            for fixture in ("constant", "piecewise"):
                compare_rows(
                    baseline["positions"][fixture],
                    current["positions"][fixture],
                    f"{paths[ranks]}:position:{fixture}",
                    set(),
                )
            # The composition fixture executes one real K=2,N=8 transaction
            # per rank.  Its physical/correction metadata is decomposition-
            # identical while its emitted counter deliberately scales.
            compare_rows(
                baseline["composition"],
                current["composition"],
                f"{paths[ranks]}:composition",
                {"emitted"},
            )
            # Transaction integer counts are independently checked as exact
            # rank multiples in parse_log().  Each rank deliberately runs its
            # transaction in its first local block with a rank-specific recoil
            # seed.  The error measures must each satisfy the strict physical
            # tolerance above, but their maxima need not be bitwise identical
            # across decompositions: MPI_MAX selects among different valid
            # local roundoff residues.  The deterministic global-event
            # analytic fixtures above remain exactly rank invariant.
            compare_rows(
                baseline["transaction"],
                current["transaction"],
                f"{paths[ranks]}:transaction",
                {
                    "ranks",
                    "parents",
                    "survivors",
                    "reacted",
                    "emitted",
                    "processed",
                    "remaining_moves",
                    "capacity_rollbacks",
                    "byte_rollbacks",
                    "rng_rollbacks",
                    "path_rollbacks",
                    "correction_error",
                    "position_error",
                },
            )
            # The fixture is one global family, not one family per rank.
            # Counts and physical outcomes therefore remain invariant.  Only
            # the rank metadata, the expected serial-vs-MPI face flag, and
            # independently bounded floating-point residues may differ.
            compare_rows(
                baseline["mover_production"],
                current["mover_production"],
                f"{paths[ranks]}:mover-production",
                {
                    "ranks",
                    "mpi_face",
                    "correction_error",
                    "event_position_error",
                    "remaining_fraction_error",
                    "endpoint_error",
                },
            )
    except (OSError, ValueError) as error:
        parser.error(str(error))

    print(
        "PASS Stage-8 total-reaction forcing integrated matrix "
        f"1rank=pass 4rank=pass 12rank=pass events={arguments.expected_events} "
        "homogeneous=1 piecewise=1 decay=1 composition=1 rollback=1 mover=1"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
