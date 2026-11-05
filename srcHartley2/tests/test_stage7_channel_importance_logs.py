#!/usr/bin/env python3
"""Strictly validate Stage-7 channel-importance logs and sidecars.

The C++ validator owns its numerical pass/fail decision.  This independent
checker treats both the text log and the diagnostic CSV/JSON as durable
interfaces.  It recomputes normalization, likelihood, weighted-yield,
parent-loss, and ESS identities instead of accepting a zero process status or
the producer's PASS token alone.  The deterministic global-event validator is
required at 1, 4, and 12 ranks, so all physical records must also agree across
those MPI decompositions.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import pathlib
import re
import sys
from collections import defaultdict


FLOAT = r"[+\-0-9.eE]+"
HEADER = re.compile(
    r"^Channel-importance validation: events/family=(\d+) ranks=(\d+) "
    r"deterministic-global-events=on$"
)
REACTION = re.compile(
    rf"^(PASS|FAIL) channel-importance reaction=(\S+) parent=(\d+) "
    rf"physical=({FLOAT}) biased=({FLOAT}) likelihood=({FLOAT}) "
    rf"analog=(\d+) selected=(\d+) weighted=({FLOAT}) expected=({FLOAT}) "
    rf"sigma=({FLOAT}) variance-ratio=({FLOAT}) sink=([01]) packets=(\d+)$"
)
FAMILY = re.compile(
    rf"^(PASS|FAIL) channel-importance family parent=(\d+) channels=(\d+) "
    rf"total-rate=({FLOAT}) trials=(\d+) reacted=(\d+) analog-loss=(\d+) "
    rf"biased-loss=(\d+) family-ess=({FLOAT}) ess-fraction=({FLOAT}) "
    r"configured=([01])$"
)
MOMENTS = re.compile(
    rf"^(PASS|FAIL) channel-importance moments parent=(\d+) "
    rf"analog-max-z=({FLOAT}) biased-max-z=({FLOAT})$"
)
SUMMARY = re.compile(
    r"^Channel-importance validation result: (PASS|FAIL) families=(\d+) "
    r"biased-families=(\d+) reactions=(\d+) sinks=(\d+) events=(\d+)$"
)
DIAGNOSTICS = re.compile(
    r"^(PASS|FAIL) channel-importance diagnostics output=(\d+) "
    r"csv=(\S+) json=(\S+)$"
)

SIDECAR_FIELDS = [
    "iteration",
    "output_index",
    "sampling_length",
    "mpi_ranks",
    "reaction_index",
    "reaction_id",
    "parent_species",
    "parent_symbol",
    "configured",
    "neutral_sink",
    "physical_probability",
    "biased_probability",
    "likelihood_ratio",
    "raw_selections",
    "likelihood_sum",
    "likelihood_square_sum",
    "weighted_channel_rate_s-1",
    "family_raw_selections",
    "family_likelihood_sum",
    "family_likelihood_square_sum",
    "family_effective_sample_size",
]


def close(left: float, right: float, *, rel: float = 3.0e-12) -> bool:
    """Compare printed binary64 values while retaining a tiny absolute floor."""

    return math.isclose(left, right, rel_tol=rel, abs_tol=1.0e-300)


def finite(value: float, context: str) -> float:
    if not math.isfinite(value):
        raise ValueError(f"{context}: non-finite value {value!r}")
    return value


def selection_allowed(samples: int, probability: float) -> float:
    sigma = math.sqrt(samples * probability * (1.0 - probability))
    floor = 4096.0 * sys.float_info.epsilon * max(1.0, samples * probability)
    return max(8.0 * sigma, floor)


def weighted_allowed(samples: int, physical: float, biased: float) -> float:
    variance = samples * physical * physical * (1.0 / biased - 1.0)
    sigma = math.sqrt(max(0.0, variance))
    floor = 4096.0 * sys.float_info.epsilon * max(1.0, samples * physical)
    return max(8.0 * sigma, floor)


def parse_log(path: pathlib.Path, ranks: int, expected_events: int) -> dict:
    text = path.read_text(encoding="utf-8", errors="strict")
    if any(line.lstrip().startswith("FAIL ") or line.strip() == "FAIL"
           for line in text.splitlines()):
        raise ValueError(f"{path}: contains a FAIL record")

    headers: list[tuple[int, int]] = []
    summaries: list[tuple] = []
    reactions: dict[str, dict] = {}
    families: dict[int, dict] = {}
    moments: dict[int, dict] = {}
    diagnostics: list[tuple[str, int, str, str]] = []

    for line in text.splitlines():
        match = HEADER.match(line)
        if match:
            headers.append((int(match.group(1)), int(match.group(2))))
            continue
        match = REACTION.match(line)
        if match:
            (
                status,
                reaction,
                parent,
                physical,
                biased,
                likelihood,
                analog,
                selected,
                weighted,
                expected,
                sigma,
                variance_ratio,
                sink,
                packets,
            ) = match.groups()
            if reaction in reactions:
                raise ValueError(f"{path}: duplicate reaction row {reaction}")
            row = {
                "status": status,
                "reaction": reaction,
                "parent": int(parent),
                "physical": finite(float(physical), f"{path}:{reaction}:q"),
                "biased": finite(float(biased), f"{path}:{reaction}:q*"),
                "likelihood": finite(float(likelihood), f"{path}:{reaction}:L"),
                "analog": int(analog),
                "selected": int(selected),
                "weighted": finite(float(weighted), f"{path}:{reaction}:weighted"),
                "expected": finite(float(expected), f"{path}:{reaction}:expected"),
                "sigma": finite(float(sigma), f"{path}:{reaction}:sigma"),
                "variance_ratio": finite(
                    float(variance_ratio), f"{path}:{reaction}:variance-ratio"
                ),
                "sink": int(sink),
                "packets": int(packets),
            }
            if status != "PASS":
                raise ValueError(f"{path}: non-passing reaction {reaction}")
            if not (0.0 < row["physical"] <= 1.0):
                raise ValueError(f"{path}: invalid physical q for {reaction}")
            if not (0.0 < row["biased"] <= 1.0):
                raise ValueError(f"{path}: invalid biased q* for {reaction}")
            if row["likelihood"] <= 0.0 or not close(
                row["likelihood"], row["physical"] / row["biased"]
            ):
                raise ValueError(f"{path}: q/q* != L for {reaction}")
            if row["packets"] < 1:
                raise ValueError(f"{path}: invalid packet count for {reaction}")
            reactions[reaction] = row
            continue

        match = FAMILY.match(line)
        if match:
            (
                status,
                parent,
                channels,
                total_rate,
                trials,
                reacted,
                analog_loss,
                biased_loss,
                ess,
                ess_fraction,
                configured,
            ) = match.groups()
            parent_index = int(parent)
            if parent_index in families:
                raise ValueError(f"{path}: duplicate family row parent={parent}")
            row = {
                "status": status,
                "parent": parent_index,
                "channels": int(channels),
                "total_rate": finite(float(total_rate), f"{path}:parent={parent}:rate"),
                "trials": int(trials),
                "reacted": int(reacted),
                "analog_loss": int(analog_loss),
                "biased_loss": int(biased_loss),
                "ess": finite(float(ess), f"{path}:parent={parent}:ESS"),
                "ess_fraction": finite(
                    float(ess_fraction), f"{path}:parent={parent}:ESS-fraction"
                ),
                "configured": int(configured),
            }
            if status != "PASS":
                raise ValueError(f"{path}: non-passing family parent={parent}")
            if row["channels"] <= 0 or row["total_rate"] <= 0.0:
                raise ValueError(f"{path}: invalid family parent={parent}")
            if row["trials"] != expected_events:
                raise ValueError(f"{path}: wrong trial count for parent={parent}")
            if not 0 < row["reacted"] <= row["trials"]:
                raise ValueError(f"{path}: invalid reacted count for parent={parent}")
            if row["analog_loss"] != row["reacted"] or row["biased_loss"] != row["reacted"]:
                raise ValueError(f"{path}: channel bias changed parent loss for parent={parent}")
            # The validator constructs dt so the analog reaction probability is
            # exactly 0.625.  Recheck that it did not silently validate only the
            # conditional selector while changing the total parent lifetime.
            reacted_allowed = selection_allowed(row["trials"], 0.625)
            if abs(row["reacted"] - 0.625 * row["trials"]) > reacted_allowed:
                raise ValueError(f"{path}: analog parent lifetime test fails for parent={parent}")
            families[parent_index] = row
            continue

        match = MOMENTS.match(line)
        if match:
            status, parent, analog_z, biased_z = match.groups()
            parent_index = int(parent)
            if parent_index in moments:
                raise ValueError(f"{path}: duplicate moment row parent={parent}")
            row = {
                "status": status,
                "analog_z": finite(float(analog_z), f"{path}:analog-z"),
                "biased_z": finite(float(biased_z), f"{path}:biased-z"),
            }
            if status != "PASS" or not (0.0 <= row["analog_z"] <= 1.0) or not (
                0.0 <= row["biased_z"] <= 1.0
            ):
                raise ValueError(f"{path}: invalid weighted moment row parent={parent}")
            moments[parent_index] = row
            continue

        match = SUMMARY.match(line)
        if match:
            summaries.append(
                (match.group(1), *(int(value) for value in match.groups()[1:]))
            )
            continue

        match = DIAGNOSTICS.match(line)
        if match:
            diagnostics.append(
                (match.group(1), int(match.group(2)), match.group(3), match.group(4))
            )

    if headers != [(expected_events, ranks)]:
        raise ValueError(f"{path}: invalid header(s) {headers}")
    if len(summaries) != 1:
        raise ValueError(f"{path}: expected one summary, found {len(summaries)}")
    expected_output = 7_000_000 + ranks
    expected_csv = (
        f"PT/diagnostics/comet_stage7_importance.out={expected_output}.csv"
    )
    expected_json = (
        f"PT/diagnostics/comet_stage7_importance.out={expected_output}.json"
    )
    if diagnostics != [("PASS", expected_output, expected_csv, expected_json)]:
        raise ValueError(f"{path}: invalid diagnostics record(s) {diagnostics}")
    status, family_count, biased_count, reaction_count, sink_count, events = summaries[0]
    if status != "PASS" or events != expected_events:
        raise ValueError(f"{path}: invalid summary {summaries[0]}")
    if family_count != len(families) or reaction_count != len(reactions):
        raise ValueError(f"{path}: summary coverage does not match parsed rows")
    if sink_count != sum(row["sink"] for row in reactions.values()):
        raise ValueError(f"{path}: summary sink count does not match rows")
    if set(families) != set(moments):
        raise ValueError(f"{path}: family and moment coverage differ")

    by_parent: dict[int, list[dict]] = defaultdict(list)
    for row in reactions.values():
        if row["parent"] not in families:
            raise ValueError(f"{path}: reaction {row['reaction']} has no family row")
        by_parent[row["parent"]].append(row)
    if set(by_parent) != set(families):
        raise ValueError(f"{path}: a family has no reaction rows")

    for parent, rows in by_parent.items():
        family = families[parent]
        if len(rows) != family["channels"]:
            raise ValueError(f"{path}: channel count mismatch for parent={parent}")
        if not close(sum(row["physical"] for row in rows), 1.0) or not close(
            sum(row["biased"] for row in rows), 1.0
        ):
            raise ValueError(f"{path}: q or q* is not normalized for parent={parent}")
        if sum(row["analog"] for row in rows) != family["reacted"] or sum(
            row["selected"] for row in rows
        ) != family["reacted"]:
            raise ValueError(f"{path}: selected channels do not close for parent={parent}")

        family_weight = 0.0
        family_weight2 = 0.0
        for row in rows:
            expected = family["reacted"] * row["physical"]
            weighted = row["selected"] * row["likelihood"]
            weighted_variance = family["reacted"] * row["physical"] ** 2 * (
                1.0 / row["biased"] - 1.0
            )
            weighted_sigma = math.sqrt(max(0.0, weighted_variance))
            analog_allowed = selection_allowed(family["reacted"], row["physical"])
            biased_allowed = selection_allowed(family["reacted"], row["biased"])
            if abs(row["analog"] - expected) > analog_allowed:
                raise ValueError(f"{path}: analog yield outside uncertainty for {row['reaction']}")
            if abs(row["selected"] - family["reacted"] * row["biased"]) > biased_allowed:
                raise ValueError(f"{path}: biased selection outside uncertainty for {row['reaction']}")
            if not close(row["weighted"], weighted) or not close(row["expected"], expected):
                raise ValueError(f"{path}: weighted-yield identity fails for {row['reaction']}")
            if not close(row["sigma"], weighted_sigma) or abs(weighted - expected) > weighted_allowed(
                family["reacted"], row["physical"], row["biased"]
            ):
                raise ValueError(f"{path}: weighted yield outside uncertainty for {row['reaction']}")
            analog_variance = row["physical"] * (1.0 - row["physical"])
            biased_variance = row["physical"] ** 2 * (1.0 / row["biased"] - 1.0)
            ratio = biased_variance / analog_variance if analog_variance > 0.0 else 1.0
            if not close(row["variance_ratio"], ratio):
                raise ValueError(f"{path}: variance ratio is wrong for {row['reaction']}")
            family_weight += weighted
            family_weight2 += row["selected"] * row["likelihood"] ** 2

        ess = family_weight * family_weight / family_weight2 if family_weight2 > 0.0 else 0.0
        if not close(family["ess"], ess) or not close(
            family["ess_fraction"], ess / family["reacted"]
        ):
            raise ValueError(f"{path}: ESS identity fails for parent={parent}")
        if not (1.0 <= ess <= family["reacted"] * (1.0 + 1.0e-12)):
            raise ValueError(f"{path}: invalid ESS bound for parent={parent}")
        if not family["configured"]:
            for row in rows:
                if not close(row["physical"], row["biased"]) or not close(
                    row["likelihood"], 1.0
                ):
                    raise ValueError(f"{path}: unconfigured family changed probabilities")

    if biased_count != sum(family["configured"] for family in families.values()):
        raise ValueError(f"{path}: biased-family summary count is wrong")
    if biased_count <= 0:
        raise ValueError(f"{path}: validation did not exercise a biased family")

    configured_rows = [
        row for row in reactions.values() if families[row["parent"]]["configured"]
    ]
    if not any(family["configured"] and family["channels"] == 2 for family in families.values()):
        raise ValueError(f"{path}: missing configured two-channel coverage")
    if not any(family["configured"] and family["channels"] >= 3 for family in families.values()):
        raise ValueError(f"{path}: missing configured many-channel coverage")
    if not any(row["physical"] > 0.5 for row in reactions.values()):
        raise ValueError(f"{path}: missing dominant-channel coverage")
    rare = [row for row in configured_rows if row["physical"] < 0.05]
    if not rare or not any(
        row["biased"] > row["physical"] and row["variance_ratio"] < 1.0
        for row in rare
    ):
        raise ValueError(f"{path}: missing variance-reduced rare-channel coverage")
    if not any(row["sink"] for row in configured_rows):
        raise ValueError(f"{path}: missing configured sink coverage")
    if not any(row["packets"] > 1 for row in configured_rows):
        raise ValueError(f"{path}: missing Stage-6 packet interaction coverage")

    return {
        "reactions": reactions,
        "families": families,
        "moments": moments,
        "summary": summaries[0],
        "diagnostics": diagnostics[0],
    }


def csv_bool(value: str, context: str) -> bool:
    if value == "true":
        return True
    if value == "false":
        return False
    raise ValueError(f"{context}: expected true/false, found {value!r}")


def normalize_sidecar_row(row: dict, context: str) -> dict:
    """Convert either a CSV or JSON record to one canonical typed schema."""

    def get(name: str):
        if name not in row:
            raise ValueError(f"{context}: missing field {name}")
        return row[name]

    csv_record = isinstance(get("configured"), str)
    configured = csv_bool(get("configured"), context) if csv_record else get("configured")
    sink = csv_bool(get("neutral_sink"), context) if csv_record else get("neutral_sink")
    if not isinstance(configured, bool) or not isinstance(sink, bool):
        raise ValueError(f"{context}: configured/sink fields are not booleans")
    typed = {
        "index": int(get("reaction_index") if csv_record else get("index")),
        "id": str(get("reaction_id") if csv_record else get("id")),
        "parent": int(get("parent_species")),
        "parent_symbol": str(get("parent_symbol")),
        "configured": configured,
        "sink": sink,
        "physical": finite(float(get("physical_probability")), f"{context}:q"),
        "biased": finite(float(get("biased_probability")), f"{context}:q*"),
        "likelihood": finite(float(get("likelihood_ratio")), f"{context}:L"),
        "raw": int(get("raw_selections")),
        "sum": finite(float(get("likelihood_sum")), f"{context}:sum-L"),
        "sum2": finite(float(get("likelihood_square_sum")), f"{context}:sum-L2"),
        "rate": finite(
            float(get("weighted_channel_rate_s-1")), f"{context}:weighted-rate"
        ),
        "family_raw": int(get("family_raw_selections")),
        "family_sum": finite(
            float(get("family_likelihood_sum")), f"{context}:family-sum"
        ),
        "family_sum2": finite(
            float(get("family_likelihood_square_sum")), f"{context}:family-sum2"
        ),
        "family_ess": finite(
            float(get("family_effective_sample_size")), f"{context}:family-ess"
        ),
    }
    if typed["index"] < 0 or typed["parent"] < 0 or typed["raw"] < 0 or typed[
        "family_raw"
    ] < 0:
        raise ValueError(f"{context}: negative index/count")
    if typed["rate"] < 0.0 or typed["sum"] < 0.0 or typed["sum2"] < 0.0:
        raise ValueError(f"{context}: negative additive statistic")
    return typed


def parse_sidecars(csv_path: pathlib.Path, json_path: pathlib.Path, ranks: int) -> dict:
    with csv_path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != SIDECAR_FIELDS:
            raise ValueError(f"{csv_path}: unexpected columns {reader.fieldnames}")
        csv_source = list(reader)
    document = json.loads(json_path.read_text(encoding="utf-8"))
    if document.get("schema") != "amps.comet.stage7-channel-importance.v1":
        raise ValueError(f"{json_path}: unexpected schema")
    if document.get("mpi_ranks") != ranks:
        raise ValueError(f"{json_path}: mpi_ranks does not match supplied rank")
    if not isinstance(document.get("iteration"), int) or document["iteration"] < 0:
        raise ValueError(f"{json_path}: invalid iteration")
    expected_output = 7_000_000 + ranks
    if document.get("output_index") != expected_output:
        raise ValueError(f"{json_path}: expected output_index={expected_output}")
    if not isinstance(document.get("sampling_length"), int) or document["sampling_length"] <= 0:
        raise ValueError(f"{json_path}: invalid sampling_length")
    if document.get("effective_sample_size_definition") != (
        "(sum L)^2/sum L^2 over every committed selection in one parent-channel family"
    ):
        raise ValueError(f"{json_path}: unexpected ESS definition")
    json_source = document.get("reactions")
    if not isinstance(json_source, list) or not json_source:
        raise ValueError(f"{json_path}: empty/missing reaction table")
    if len(csv_source) != len(json_source):
        raise ValueError(f"{csv_path}: CSV/JSON reaction counts differ")

    csv_rows: dict[str, dict] = {}
    json_rows: dict[str, dict] = {}
    for raw in csv_source:
        if int(raw["mpi_ranks"]) != ranks or int(raw["output_index"]) != expected_output:
            raise ValueError(f"{csv_path}: row provenance differs from file")
        if int(raw["iteration"]) != document["iteration"] or int(
            raw["sampling_length"]
        ) != document["sampling_length"]:
            raise ValueError(f"{csv_path}: CSV/JSON provenance differs")
        row = normalize_sidecar_row(raw, f"{csv_path}:{raw.get('reaction_id')}")
        if row["id"] in csv_rows:
            raise ValueError(f"{csv_path}: duplicate reaction {row['id']}")
        csv_rows[row["id"]] = row
    for raw in json_source:
        if not isinstance(raw, dict):
            raise ValueError(f"{json_path}: reaction row is not an object")
        row = normalize_sidecar_row(raw, f"{json_path}:{raw.get('id')}")
        if row["id"] in json_rows:
            raise ValueError(f"{json_path}: duplicate reaction {row['id']}")
        json_rows[row["id"]] = row
    if set(csv_rows) != set(json_rows):
        raise ValueError(f"{csv_path}: CSV/JSON reaction identifiers differ")

    for reaction, left in csv_rows.items():
        right = json_rows[reaction]
        for field in (
            "index",
            "id",
            "parent",
            "parent_symbol",
            "configured",
            "sink",
            "raw",
            "family_raw",
        ):
            if left[field] != right[field]:
                raise ValueError(f"{csv_path}: CSV/JSON {reaction}/{field} differs")
        for field in (
            "physical",
            "biased",
            "likelihood",
            "sum",
            "sum2",
            "rate",
            "family_sum",
            "family_sum2",
            "family_ess",
        ):
            if not close(left[field], right[field]):
                raise ValueError(f"{csv_path}: CSV/JSON {reaction}/{field} differs")

    indices = [row["index"] for row in csv_rows.values()]
    if sorted(indices) != list(range(len(indices))):
        raise ValueError(f"{csv_path}: reaction indices are not contiguous")

    by_parent: dict[int, list[dict]] = defaultdict(list)
    for row in csv_rows.values():
        by_parent[row["parent"]].append(row)
        if row["physical"] == 0.0:
            if row["biased"] != 0.0 or row["configured"] or row["raw"] != 0 or any(
                row[field] != 0.0 for field in ("sum", "sum2", "rate")
            ):
                raise ValueError(f"{csv_path}: inactive reaction {row['id']} has activity")
            if not close(row["likelihood"], 1.0):
                raise ValueError(f"{csv_path}: inactive reaction {row['id']} has L != 1")
        else:
            if not (0.0 < row["biased"] <= 1.0) or not close(
                row["likelihood"], row["physical"] / row["biased"]
            ):
                raise ValueError(f"{csv_path}: q/q* != L for {row['id']}")
            if not close(row["sum"], row["raw"] * row["likelihood"]) or not close(
                row["sum2"], row["raw"] * row["likelihood"] ** 2
            ):
                raise ValueError(f"{csv_path}: likelihood sums fail for {row['id']}")

    for parent, rows in by_parent.items():
        active = [row for row in rows if row["physical"] > 0.0]
        if not active:
            continue
        if not close(sum(row["physical"] for row in active), 1.0) or not close(
            sum(row["biased"] for row in active), 1.0
        ):
            raise ValueError(f"{csv_path}: sidecar family {parent} is not normalized")
        configured = {row["configured"] for row in active}
        if len(configured) != 1:
            raise ValueError(f"{csv_path}: sidecar family {parent} is partially configured")
        family_raw = sum(row["raw"] for row in rows)
        family_sum = sum(row["sum"] for row in rows)
        family_sum2 = sum(row["sum2"] for row in rows)
        family_ess = family_sum * family_sum / family_sum2 if family_sum2 > 0.0 else 0.0
        for row in rows:
            if row["family_raw"] != family_raw or not close(
                row["family_sum"], family_sum
            ) or not close(row["family_sum2"], family_sum2) or not close(
                row["family_ess"], family_ess
            ):
                raise ValueError(f"{csv_path}: repeated family metadata fails for parent={parent}")
        if not (0.0 <= family_ess <= family_raw * (1.0 + 1.0e-12)):
            raise ValueError(f"{csv_path}: invalid sidecar ESS for parent={parent}")

    if not any(row["configured"] for row in csv_rows.values()):
        raise ValueError(f"{csv_path}: no configured Stage-7 rows")
    if not any(row["configured"] and row["sink"] for row in csv_rows.values()):
        raise ValueError(f"{csv_path}: no configured sink in sidecar")

    return {
        "rows": csv_rows,
        "iteration": document["iteration"],
        "sampling_length": document["sampling_length"],
    }


def compare_runtime_and_sidecar(log: dict, sidecar: dict, path: pathlib.Path) -> None:
    for reaction, runtime in log["reactions"].items():
        if reaction not in sidecar["rows"]:
            raise ValueError(f"{path}: sidecar is missing runtime reaction {reaction}")
        row = sidecar["rows"][reaction]
        family_configured = bool(log["families"][runtime["parent"]]["configured"])
        if row["parent"] != runtime["parent"] or row["sink"] != bool(runtime["sink"]):
            raise ValueError(f"{path}: sidecar family/sink mismatch for {reaction}")
        if row["configured"] != family_configured:
            raise ValueError(f"{path}: sidecar configured flag mismatch for {reaction}")
        # The focused validator intentionally feeds exactly the same committed
        # selections into the production diagnostic transaction.  Therefore
        # its durable count and likelihood moments must agree with the log,
        # not merely have plausible signs and magnitudes.
        if row["raw"] != runtime["selected"]:
            raise ValueError(f"{path}: sidecar raw count differs for {reaction}")
        if not close(row["sum"], runtime["weighted"]):
            raise ValueError(f"{path}: sidecar likelihood sum differs for {reaction}")
        expected_rate = (
            log["families"][runtime["parent"]]["total_rate"]
            * row["sum"]
            / sidecar["sampling_length"]
        )
        # The writer accumulates parentRate*L once per committed carrier and
        # divides only after the sampling window.  Unlike the integer count
        # and the separately accumulated L moments, that repeated addition
        # may incur O(N*epsilon) binary64 rounding.  Bound it by eight ulps per
        # addition rather than weakening the global comparison tolerance used
        # for probabilities, likelihoods, ESS, and cross-rank identity.
        rate_relative_tolerance = max(
            3.0e-12, 8.0 * row["raw"] * sys.float_info.epsilon
        )
        if not math.isclose(
            row["rate"],
            expected_rate,
            rel_tol=rate_relative_tolerance,
            abs_tol=1.0e-300,
        ):
            raise ValueError(f"{path}: sidecar weighted rate differs for {reaction}")
        for sidecar_field, runtime_field in (
            ("physical", "physical"),
            ("biased", "biased"),
            ("likelihood", "likelihood"),
        ):
            if not close(row[sidecar_field], runtime[runtime_field]):
                raise ValueError(f"{path}: sidecar probability metadata differs for {reaction}")


def compare_rank_records(baseline: dict, current: dict, context: str) -> None:
    if set(baseline) != set(current):
        raise ValueError(f"{context}: record keys differ from one-rank baseline")
    for key, left in baseline.items():
        right = current[key]
        if set(left) != set(right):
            raise ValueError(f"{context}: fields differ for {key}")
        for field, left_value in left.items():
            right_value = right[field]
            if isinstance(left_value, float):
                # These fields are formed from repeated floating additions
                # and then an MPI reduction.  Different rank partitions
                # legitimately change that addition order, so use an explicit
                # N*epsilon roundoff bound.  Probabilities and likelihood
                # metadata retain the strict general tolerance; integers and
                # identifiers below must remain bit-for-bit identical.
                accumulated_fields = {
                    "sum",
                    "sum2",
                    "rate",
                    "family_sum",
                    "family_sum2",
                    "family_ess",
                }
                if field in accumulated_fields:
                    count_field = "family_raw" if field.startswith("family_") else "raw"
                    accumulation_relative_tolerance = max(
                        3.0e-12,
                        8.0
                        * max(left.get(count_field, 0), right.get(count_field, 0))
                        * sys.float_info.epsilon,
                    )
                    values_match = math.isclose(
                        left_value,
                        right_value,
                        rel_tol=accumulation_relative_tolerance,
                        abs_tol=1.0e-300,
                    )
                else:
                    values_match = close(left_value, right_value)
                if not values_match:
                    raise ValueError(f"{context}: {key}/{field} differs by MPI rank count")
            elif left_value != right_value:
                raise ValueError(f"{context}: {key}/{field} differs by MPI rank count")


def collect_rank_paths(
    values: list[list[str]], option: str, *, require_file: bool = True
) -> dict[int, pathlib.Path]:
    paths: dict[int, pathlib.Path] = {}
    for rank_text, path_text in values:
        ranks = int(rank_text)
        if ranks in paths:
            raise ValueError(f"duplicate {option} for {ranks} ranks")
        path = pathlib.Path(path_text).resolve()
        if require_file and not path.is_file():
            raise ValueError(f"missing {option} for {ranks} ranks: {path}")
        if not require_file and not path.is_dir():
            raise ValueError(f"missing {option} for {ranks} ranks: {path}")
        paths[ranks] = path
    if set(paths) != {1, 4, 12}:
        raise ValueError(f"{option} must be supplied for exactly 1, 4, and 12 ranks")
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(
        description="validate Stage-7 1/4/12-rank logs and importance sidecars"
    )
    parser.add_argument("--expected-events", type=int, required=True)
    parser.add_argument(
        "--importance-log",
        nargs=2,
        action="append",
        required=True,
        metavar=("RANKS", "PATH"),
    )
    parser.add_argument(
        "--diagnostics-root",
        nargs=2,
        action="append",
        required=True,
        metavar=("RANKS", "CASE_DIR"),
        help=(
            "repeat for ranks 1, 4, and 12; CASE_DIR contains PT/diagnostics/"
            "comet_stage7_importance.out=7000000+RANKS.{csv,json}"
        ),
    )
    arguments = parser.parse_args()
    if arguments.expected_events < 1000:
        parser.error("--expected-events must be at least 1000")

    try:
        logs = collect_rank_paths(arguments.importance_log, "--importance-log")
        roots = collect_rank_paths(
            arguments.diagnostics_root, "--diagnostics-root", require_file=False
        )
        csvs: dict[int, pathlib.Path] = {}
        jsons: dict[int, pathlib.Path] = {}
        for ranks in (1, 4, 12):
            diagnostic_directory = roots[ranks] / "PT" / "diagnostics"
            matching_csv = sorted(
                diagnostic_directory.glob("comet_stage7_importance.out=*.csv")
            )
            matching_json = sorted(
                diagnostic_directory.glob("comet_stage7_importance.out=*.json")
            )
            # Exactly one pair is intentional.  Accepting an older sidecar in
            # addition to the reserved validator output could silently mix
            # diagnostics from different executions in a reused case tree.
            if len(matching_csv) != 1 or len(matching_json) != 1:
                raise ValueError(
                    f"{diagnostic_directory}: expected exactly one Stage-7 "
                    f"CSV/JSON pair, found {len(matching_csv)}/{len(matching_json)}"
                )
            expected_stem = (
                f"comet_stage7_importance.out={7_000_000 + ranks}"
            )
            if matching_csv[0].stem != expected_stem or matching_json[0].stem != expected_stem:
                raise ValueError(
                    f"{diagnostic_directory}: Stage-7 sidecar output index is not "
                    f"{7_000_000 + ranks}"
                )
            csvs[ranks] = matching_csv[0]
            jsons[ranks] = matching_json[0]
        parsed_logs = {
            ranks: parse_log(path, ranks, arguments.expected_events)
            for ranks, path in sorted(logs.items())
        }
        parsed_sidecars = {
            ranks: parse_sidecars(csvs[ranks], jsons[ranks], ranks)
            for ranks in (1, 4, 12)
        }
        for ranks in (1, 4, 12):
            compare_runtime_and_sidecar(parsed_logs[ranks], parsed_sidecars[ranks], csvs[ranks])

        for ranks in (4, 12):
            compare_rank_records(
                parsed_logs[1]["reactions"],
                parsed_logs[ranks]["reactions"],
                str(logs[ranks]),
            )
            compare_rank_records(
                parsed_logs[1]["families"],
                parsed_logs[ranks]["families"],
                str(logs[ranks]),
            )
            compare_rank_records(
                parsed_logs[1]["moments"],
                parsed_logs[ranks]["moments"],
                str(logs[ranks]),
            )
            compare_rank_records(
                parsed_sidecars[1]["rows"],
                parsed_sidecars[ranks]["rows"],
                str(csvs[ranks]),
            )

        neutral_reactions=len(parsed_logs[1]["reactions"])
        families=len(parsed_logs[1]["families"])
        sidecar_rows=len(parsed_sidecars[1]["rows"])
        print(
            "PASS Stage-7 channel-importance integrated matrix "
            "1rank=pass 4rank=pass 12rank=pass "
            f"families={families} reactions={neutral_reactions} "
            f"sidecar-rows={sidecar_rows} events={arguments.expected_events}"
        )
        return 0
    except (OSError, ValueError, json.JSONDecodeError) as error:
        parser.error(str(error))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
