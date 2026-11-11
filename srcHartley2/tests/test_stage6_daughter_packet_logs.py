#!/usr/bin/env python3
"""Strictly validate the Stage-6 1/4/12-rank daughter-packet logs.

The C++ validator owns numerical acceptance.  This independent parser makes
the log schema an interface: it requires complete reaction/N coverage,
approved policy precedence, conservation bounds, variance reduction, and the
same deterministic global-event ensemble at every MPI decomposition.
"""

from __future__ import annotations

import argparse
import math
import pathlib
import re
import sys


ROW = re.compile(
    r"^(PASS|FAIL) daughter-packet reaction=(\S+) packets=(\d+) "
    r"products=(\d+) configured=(\d+) source=(\d+) selector=(\S+) "
    r"events=(\d+) daughters=(\d+) "
    r"velocity-mean=([+\-0-9.eE]+) velocity-variance=([+\-0-9.eE]+) "
    r"velocity-scaled=([+\-0-9.eE]+) density-mean=([+\-0-9.eE]+) "
    r"density-variance=([+\-0-9.eE]+) density-scaled=([+\-0-9.eE]+) "
    r"momentum-error=([+\-0-9.eE]+) energy-error=([+\-0-9.eE]+) "
    r"wall=([+\-0-9.eE]+)$"
)
HEADER = re.compile(
    r"^Daughter-packet validation: events/reaction=(\d+) ranks=(\d+) "
    r"packets=1,2,4,8 deterministic-global-events=on$"
)
SUMMARY = re.compile(
    r"^Daughter-packet validation result: (PASS|FAIL) "
    r"neutral-reactions=(\d+) sink-reactions=(\d+) packet-cases=(\d+)$"
)
PRODUCTION = re.compile(
    r"^(PASS|FAIL) daughter-packet production-transaction ranks=(\d+) "
    r"reaction=C2_2C packets=8 products=2 emitted=(\d+) reused=(\d+) "
    r"processed=(\d+) link-checks=(\d+) capacity-rollbacks=(\d+) "
    r"correction-error=([+\-0-9.eE]+) represented-momentum=([+\-0-9.eE]+) "
    r"represented-energy=([+\-0-9.eE]+) max-seed-attempts=(\d+)$"
)


def close(a: float, b: float) -> bool:
    return math.isclose(a, b, rel_tol=2.0e-12, abs_tol=1.0e-300)


def parse_log(path: pathlib.Path, ranks: int, expected_events: int):
    text = path.read_text(encoding="utf-8", errors="strict")
    if "FAIL daughter-packet" in text:
        raise ValueError(f"{path}: contains a failed daughter-packet row")

    rows = {}
    summaries = []
    headers = []
    production_records = []
    for line in text.splitlines():
        match = HEADER.match(line)
        if match:
            headers.append((int(match.group(1)), int(match.group(2))))
            continue
        match = ROW.match(line)
        if match:
            (
                status,
                reaction,
                packet_text,
                products_text,
                configured_text,
                source_text,
                selector,
                events_text,
                daughters_text,
                velocity_mean_text,
                velocity_variance_text,
                velocity_scaled_text,
                density_mean_text,
                density_variance_text,
                density_scaled_text,
                momentum_text,
                energy_text,
                wall_text,
            ) = match.groups()
            packet = int(packet_text)
            key = (reaction, packet)
            if key in rows:
                raise ValueError(f"{path}: duplicate row {reaction}/N={packet}")
            values = {
                "products": int(products_text),
                "configured": int(configured_text),
                "source": int(source_text),
                "selector": selector,
                "events": int(events_text),
                "daughters": int(daughters_text),
                "velocity_mean": float(velocity_mean_text),
                "velocity_variance": float(velocity_variance_text),
                "velocity_scaled": float(velocity_scaled_text),
                "density_mean": float(density_mean_text),
                "density_variance": float(density_variance_text),
                "density_scaled": float(density_scaled_text),
                "momentum_error": float(momentum_text),
                "energy_error": float(energy_text),
                "wall": float(wall_text),
            }
            if status != "PASS":
                raise ValueError(f"{path}: non-passing row {reaction}/N={packet}")
            if packet not in (1, 2, 4, 8):
                raise ValueError(f"{path}: unexpected packet count {packet}")
            if values["events"] != expected_events:
                raise ValueError(
                    f"{path}: {reaction}/N={packet} events="
                    f"{values['events']} expected={expected_events}"
                )
            if values["products"] <= 0 or values["daughters"] != (
                expected_events * packet * values["products"]
            ) or (
                values["daughters"] % (expected_events * packet) != 0
            ):
                raise ValueError(
                    f"{path}: invalid daughter cost for {reaction}/N={packet}"
                )
            if values["configured"] not in (2, 4, 8):
                raise ValueError(f"{path}: unexpected configured count")
            if values["source"] not in (2, 3, 4):
                raise ValueError(f"{path}: unexpected policy source")
            expected_selector = {2: "chemical-only", 3: "species", 4: "reaction"}[
                values["source"]
            ]
            if values["selector"] != expected_selector:
                raise ValueError(f"{path}: policy source/name mismatch")
            numeric = [value for value in values.values() if isinstance(value, float)]
            if not all(math.isfinite(value) for value in numeric):
                raise ValueError(f"{path}: non-finite row {reaction}/N={packet}")
            if values["velocity_variance"] < 0.0 or values["density_variance"] < 0.0:
                raise ValueError(f"{path}: negative variance")
            if values["momentum_error"] > 1.0e-11 or values["energy_error"] > 1.0e-11:
                raise ValueError(f"{path}: conservation tolerance exceeded")
            if values["wall"] <= 0.0:
                raise ValueError(f"{path}: non-positive timing")
            if packet == 1:
                if not close(values["velocity_scaled"], 1.0) or not close(
                    values["density_scaled"], 1.0
                ):
                    raise ValueError(f"{path}: N=1 scaled variance is not one")
            elif not (
                0.5 <= values["velocity_scaled"] <= 1.6
                and 0.5 <= values["density_scaled"] <= 1.6
            ):
                raise ValueError(f"{path}: variance reduction outside accepted band")
            rows[key] = values
            continue

        match = SUMMARY.match(line)
        if match:
            summaries.append((match.group(1), *(int(value) for value in match.groups()[1:])))
            continue
        match = PRODUCTION.match(line)
        if match:
            production_records.append(match.groups())

    if headers != [(expected_events, ranks)]:
        raise ValueError(
            f"{path}: expected one header events={expected_events} ranks={ranks}; "
            f"found {headers}"
        )
    if len(summaries) != 1:
        raise ValueError(f"{path}: expected one summary, found {len(summaries)}")
    status, neutral, sinks, cases = summaries[0]
    if status != "PASS" or neutral <= 0 or sinks <= 0 or cases != 4:
        raise ValueError(f"{path}: invalid summary {summaries[0]}")
    if len(rows) != neutral * cases:
        raise ValueError(
            f"{path}: row count {len(rows)} does not match {neutral}*{cases}"
        )
    if len(production_records) != 1:
        raise ValueError(
            f"{path}: expected one production transaction record, found "
            f"{len(production_records)}"
        )
    (
        production_status,
        production_ranks,
        emitted,
        reused,
        processed,
        link_checks,
        rollbacks,
        correction_error,
        momentum_error,
        energy_error,
        seed_attempts,
    ) = production_records[0]
    if production_status != "PASS" or int(production_ranks) != ranks:
        raise ValueError(f"{path}: invalid production transaction status/ranks")
    if (
        int(emitted) != 16 * ranks
        or int(reused) != ranks
        or int(processed) != 16 * ranks
        or int(link_checks) != ranks
        or int(rollbacks) != ranks
    ):
        raise ValueError(f"{path}: production transaction count mismatch")
    if not 1 <= int(seed_attempts) <= 256:
        raise ValueError(f"{path}: invalid production seed-search count")
    if (
        not math.isfinite(float(correction_error))
        or not math.isfinite(float(momentum_error))
        or not math.isfinite(float(energy_error))
        or float(correction_error) > 1.0e-12
        or float(momentum_error) > 1.0e-11
        or float(energy_error) > 1.0e-11
    ):
        raise ValueError(f"{path}: production represented-moment error")
    reactions = {reaction for reaction, _ in rows}
    for reaction in reactions:
        if {packet for name, packet in rows if name == reaction} != {1, 2, 4, 8}:
            raise ValueError(f"{path}: incomplete packet cases for {reaction}")
        policies = {
            (
                rows[(reaction, packet)]["products"],
                rows[(reaction, packet)]["configured"],
                rows[(reaction, packet)]["source"],
                rows[(reaction, packet)]["selector"],
            )
            for packet in (1, 2, 4, 8)
        }
        if len(policies) != 1:
            raise ValueError(f"{path}: policy changed across N for {reaction}")

    # Prove that the staged validation policy actually exercised group,
    # product-species, and reaction-ID precedence rather than merely compiling
    # those input records.
    sources = {values["source"] for values in rows.values()}
    if not {2, 3, 4}.issubset(sources):
        raise ValueError(f"{path}: missing policy precedence coverage {sources}")
    h2o = rows.get(("H2O_OH_H", 1))
    if (
        h2o is None
        or h2o["configured"] != 4
        or h2o["source"] != 4
        or h2o["selector"] != "reaction"
    ):
        raise ValueError(f"{path}: reaction-ID override was not resolved")
    if not any(
        values["configured"] == 8 and values["source"] == 3
        for values in rows.values()
    ):
        raise ValueError(f"{path}: C2 species override was not resolved")
    for reaction in ("C2H2_C2_H2", "C2H_C2_H"):
        record = rows.get((reaction, 1))
        if (
            record is None
            or record["configured"] != 8
            or record["source"] != 3
            or record["selector"] != "species"
        ):
            raise ValueError(f"{path}: named C2 species override missing for {reaction}")
    repeated = rows.get(("C2_2C", 1))
    if (
        repeated is None
        or repeated["products"] != 2
        or repeated["configured"] != 8
        or repeated["source"] != 4
    ):
        raise ValueError(f"{path}: repeated-product transaction override missing")

    return rows, (neutral, sinks, cases), production_records[0]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expected-events", type=int, required=True)
    parser.add_argument(
        "--packet-log", nargs=2, action="append", metavar=("RANKS", "PATH"), required=True
    )
    args = parser.parse_args()
    if args.expected_events < 1000:
        parser.error("--expected-events must be at least 1000")

    paths = {}
    for ranks_text, path_text in args.packet_log:
        ranks = int(ranks_text)
        if ranks in paths:
            parser.error(f"duplicate packet log for {ranks} ranks")
        paths[ranks] = pathlib.Path(path_text).resolve()
    if set(paths) != {1, 4, 12}:
        parser.error("packet logs must be supplied for exactly 1, 4, and 12 ranks")

    parsed = {
        ranks: parse_log(path, ranks, args.expected_events)
        for ranks, path in sorted(paths.items())
    }
    baseline_rows, baseline_summary, _ = parsed[1]
    for ranks in (4, 12):
        rows, summary, _ = parsed[ranks]
        if summary != baseline_summary or set(rows) != set(baseline_rows):
            raise ValueError(f"{paths[ranks]}: MPI coverage differs from one rank")
        for key, baseline in baseline_rows.items():
            current = rows[key]
            for field in (
                "products",
                "configured",
                "source",
                "selector",
                "events",
                "daughters",
            ):
                if current[field] != baseline[field]:
                    raise ValueError(
                        f"{paths[ranks]}: {key} field {field} differs by decomposition"
                    )
            for field in (
                "velocity_mean",
                "velocity_variance",
                "velocity_scaled",
                "density_mean",
                "density_variance",
                "density_scaled",
            ):
                if not close(current[field], baseline[field]):
                    raise ValueError(
                        f"{paths[ranks]}: {key} field {field} is not MPI invariant"
                    )

    neutral, sinks, cases = baseline_summary
    print(
        "PASS Stage-6 daughter-packet integrated matrix "
        f"1rank=pass 4rank=pass 12rank=pass reactions={neutral} "
        f"sinks={sinks} packet-cases={cases} events={args.expected_events}"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"FAIL Stage-6 log validation: {error}", file=sys.stderr)
        raise SystemExit(1)
