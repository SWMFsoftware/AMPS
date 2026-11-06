#!/usr/bin/env python3
"""Validate Stage-9 transported-daughter production-NTC evidence.

The application record is produced only after the exact C daughters created
by the Stage-8 ``C2_2C`` transaction have crossed a real block face, completed
their remaining local-timestep trajectory, and entered the general NTC kernel
with H2O.  This parser deliberately recomputes the record invariants rather
than trusting its leading ``PASS`` token.

Only fields that are mathematical/decomposition invariants are compared
across 1, 4, and 12 ranks.  The destination owner and collision trial counts
may differ because the physical destination block belongs to a different MPI
rank and the focused fixture uses the normal rank-specific random stream.
"""

from __future__ import annotations

import argparse
import math
import pathlib
import re
from typing import Dict, Iterable, Mapping


FLOAT = r"[+\-0-9.eE]+"
RECORD = re.compile(
    rf"^(PASS|FAIL) transported-daughter-collision "
    rf"ranks=(\d+) reaction=(\S+) product=(-?\d+) partner=(-?\d+) "
    rf"owner=(-?\d+) daughters=(\d+) partners=(\d+) "
    rf"daughter-velocity-changes=(\d+) partner-velocity-changes=(\d+) "
    rf"unchanged-records=(\d+) unfinished-trajectories=(\d+) "
    rf"wrong-owner=(\d+) cycles=(\d+) trials=(\d+) accepted=(\d+) "
    rf"updates=(\d+) probability-checks=(\d+) "
    rf"invalid-probabilities=(\d+) roundoff-clamps=(\d+) "
    rf"probability-min=({FLOAT}) probability-max=({FLOAT}) "
    rf"correction-error=({FLOAT}) represented-rate-error=({FLOAT}) "
    rf"momentum-error=({FLOAT})$"
)


def finite(text: str, label: str) -> float:
    value = float(text)
    if not math.isfinite(value):
        raise ValueError(f"{label}: non-finite value {text!r}")
    return value


def parse_log(path: pathlib.Path, expected_ranks: int) -> Dict[str, object]:
    """Parse one rank decomposition and enforce every C++ pass invariant."""

    text = path.read_text(encoding="utf-8", errors="strict")
    rows = []
    for line in text.splitlines():
        match = RECORD.match(line)
        if match is None:
            continue
        (
            status,
            ranks,
            reaction,
            product,
            partner,
            owner,
            daughters,
            partners,
            daughter_changes,
            partner_changes,
            unchanged_records,
            unfinished,
            wrong_owner,
            cycles,
            trials,
            accepted,
            updates,
            probability_checks,
            invalid_probabilities,
            roundoff_clamps,
            probability_minimum,
            probability_maximum,
            correction_error,
            represented_rate_error,
            momentum_error,
        ) = match.groups()
        rows.append(
            {
                "status": status,
                "ranks": int(ranks),
                "reaction": reaction,
                "product": int(product),
                "partner": int(partner),
                "owner": int(owner),
                "daughters": int(daughters),
                "partners": int(partners),
                "daughter_changes": int(daughter_changes),
                "partner_changes": int(partner_changes),
                "unchanged_records": int(unchanged_records),
                "unfinished": int(unfinished),
                "wrong_owner": int(wrong_owner),
                "cycles": int(cycles),
                "trials": int(trials),
                "accepted": int(accepted),
                "updates": int(updates),
                "probability_checks": int(probability_checks),
                "invalid_probabilities": int(invalid_probabilities),
                "roundoff_clamps": int(roundoff_clamps),
                "probability_minimum": finite(
                    probability_minimum, f"{path}:probability-min"
                ),
                "probability_maximum": finite(
                    probability_maximum, f"{path}:probability-max"
                ),
                "correction_error": finite(
                    correction_error, f"{path}:correction-error"
                ),
                "represented_rate_error": finite(
                    represented_rate_error, f"{path}:represented-rate-error"
                ),
                "momentum_error": finite(momentum_error, f"{path}:momentum-error"),
            }
        )

    if len(rows) != 1:
        raise ValueError(
            f"{path}: expected exactly one transported-daughter-collision record, "
            f"found {len(rows)}"
        )
    row = rows[0]
    if row["status"] != "PASS":
        raise ValueError(f"{path}: transported-daughter collision status is not PASS")
    if row["ranks"] != expected_ranks:
        raise ValueError(
            f"{path}: record ranks={row['ranks']} does not match {expected_ranks}"
        )
    if row["reaction"] != "C2_2C":
        raise ValueError(f"{path}: unexpected reaction {row['reaction']!r}")
    if row["product"] < 0 or row["partner"] < 0 or row["product"] == row["partner"]:
        raise ValueError(f"{path}: invalid product/partner species indices")
    if row["owner"] < 0 or row["owner"] >= expected_ranks:
        raise ValueError(f"{path}: destination owner is outside the MPI decomposition")

    # The Stage-8 policy creates K=2 products in N=8 packets.  H2O supplies
    # 64 equal-correction collision partners.  All 80 exact records must
    # survive the collision without a species, correction, marker, or
    # remainder mutation; only their velocities may change.
    if row["daughters"] != 16 or row["partners"] != 64:
        raise ValueError(f"{path}: expected 16 transported daughters and 64 partners")
    if row["unchanged_records"] != row["daughters"] + row["partners"]:
        raise ValueError(f"{path}: particle-record invariants were not retained")
    if not 8 <= row["daughter_changes"] <= row["daughters"]:
        raise ValueError(f"{path}: insufficient or impossible daughter velocity feedback")
    if not 8 <= row["partner_changes"] <= row["partners"]:
        raise ValueError(f"{path}: insufficient or impossible partner velocity feedback")
    if row["unfinished"] != 0 or row["wrong_owner"] != 0:
        raise ValueError(f"{path}: unfinished trajectory or wrong-owner collision")

    if not 1 <= row["cycles"] <= 16:
        raise ValueError(f"{path}: collision cycle count is outside the fixture bound")
    if row["trials"] <= 0 or row["accepted"] <= 0:
        raise ValueError(f"{path}: production NTC performed no collision work")
    if row["accepted"] > row["trials"]:
        raise ValueError(f"{path}: accepted collisions exceed NTC trials")
    if row["updates"] < 64:
        raise ValueError(f"{path}: too few production NTC particle updates")
    if row["probability_checks"] != 2 * row["accepted"]:
        raise ValueError(f"{path}: each accepted pair must check two update probabilities")
    if row["updates"] != row["probability_checks"]:
        raise ValueError(f"{path}: unit-probability collision did not update both particles")
    if row["invalid_probabilities"] != 0 or row["roundoff_clamps"] != 0:
        raise ValueError(f"{path}: invalid or clamped NTC update probability")
    if not (0.0 <= row["probability_minimum"] <= row["probability_maximum"] <= 1.0):
        raise ValueError(f"{path}: NTC probability range is outside [0,1]")
    if abs(row["probability_minimum"] - 1.0) > 2.0e-12 or abs(
        row["probability_maximum"] - 1.0
    ) > 2.0e-12:
        raise ValueError(f"{path}: equal-carrier update probabilities are not one")

    for key in ("correction_error", "represented_rate_error", "momentum_error"):
        if row[key] < 0.0:
            raise ValueError(f"{path}: {key} is negative")
    if row["correction_error"] != 0.0:
        raise ValueError(f"{path}: a collision changed an individual correction")
    if row["represented_rate_error"] != 0.0:
        raise ValueError(f"{path}: represented daughter rate changed")
    if row["momentum_error"] > 2.0e-12:
        raise ValueError(f"{path}: represented momentum does not close")
    return row


def collect_paths(entries: Iterable[Iterable[str]]) -> Mapping[int, pathlib.Path]:
    result: Dict[int, pathlib.Path] = {}
    for rank_text, path_text in entries:
        try:
            ranks = int(rank_text)
        except ValueError as exc:
            raise ValueError(f"invalid MPI rank count {rank_text!r}") from exc
        path = pathlib.Path(path_text).resolve()
        if ranks in result:
            raise ValueError(f"duplicate --collision-log for {ranks} ranks")
        if not path.is_file():
            raise ValueError(f"missing --collision-log for {ranks} ranks: {path}")
        result[ranks] = path
    if set(result) != {1, 4, 12}:
        raise ValueError(
            "--collision-log is required for exactly 1, 4, and 12 ranks"
        )
    return result


def compare_invariants(
    baseline: Mapping[str, object], current: Mapping[str, object], label: str
) -> None:
    # RNG-dependent collision work and destination ownership may differ.  The
    # chemical identity and represented-carrier invariants may not.
    fields = (
        "reaction",
        "product",
        "partner",
        "daughters",
        "partners",
        "unchanged_records",
        "unfinished",
        "wrong_owner",
        "correction_error",
        "represented_rate_error",
    )
    for field in fields:
        if baseline[field] != current[field]:
            raise ValueError(f"{label}: decomposition changed {field}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="validate Stage-9 transported-daughter NTC logs"
    )
    parser.add_argument(
        "--collision-log",
        nargs=2,
        action="append",
        required=True,
        metavar=("RANKS", "PATH"),
    )
    arguments = parser.parse_args()
    try:
        paths = collect_paths(arguments.collision_log)
        rows = {
            ranks: parse_log(path, ranks) for ranks, path in sorted(paths.items())
        }
        for ranks in (4, 12):
            compare_invariants(rows[1], rows[ranks], str(paths[ranks]))
    except (OSError, ValueError) as error:
        parser.error(str(error))

    print(
        "PASS Stage-9 transported-daughter collision matrix "
        "1rank=pass 4rank=pass 12rank=pass reaction=C2_2C "
        "daughters=16 partner-particles=64 production-ntc=1"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
