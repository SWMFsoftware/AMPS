#!/usr/bin/env python3
"""Validate Stage-1 startup policy tables and independent cadence logs."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path


POLICY = re.compile(
    r"population-policy species=(?P<species>\S+) "
    r"class=(?P<class>\S+) selector=(?P<selector>\S+) "
    r"locked=(?P<locked>yes|no) range=(?P<minimum>\d+):(?P<maximum>\d+) "
    r"surface-source-s-1=(?P<source>\S+) "
    r"surface-carrier-W-dt-s-1=(?P<carrier>\S+)"
)
CADENCE = re.compile(
    r"resolved radial population policies adaptation-cadence=(?P<mode>\S+)"
    r"(?: interval=(?P<interval>\d+))?"
)
ADAPTATION = re.compile(r"radial particle-weight adaptation iteration=(\d+)")


def parse_policy_log(path: Path) -> tuple[dict[str, dict[str, object]], tuple[str, int | None]]:
    text = path.read_text(encoding="utf-8", errors="replace")
    rows: dict[str, dict[str, object]] = {}
    for match in POLICY.finditer(text):
        species = match.group("species")
        if species in rows:
            raise AssertionError(f"duplicate startup policy for {species} in {path}")
        rows[species] = {
            "class": match.group("class"),
            "selector": match.group("selector"),
            "locked": match.group("locked") == "yes",
            "minimum": int(match.group("minimum")),
            "maximum": int(match.group("maximum")),
            "source": float(match.group("source")),
            "carrier": float(match.group("carrier")),
        }
    if not rows:
        raise AssertionError(f"no Stage-1 startup table in {path}")
    cadence_match = CADENCE.search(text)
    if cadence_match is None:
        raise AssertionError(f"no Stage-1 cadence record in {path}")
    interval = cadence_match.group("interval")
    return rows, (cadence_match.group("mode"), int(interval) if interval else None)


def validate_policy(path: Path, expected_chemical_range: tuple[int, int]) -> dict[str, dict[str, object]]:
    rows, cadence = parse_policy_log(path)
    if cadence != ("fixed-iteration", 10000):
        raise AssertionError(f"unexpected canonical cadence: {cadence}")
    classes = {row["class"] for row in rows.values()}
    if not {"direct-source", "chemical-only"} <= classes:
        raise AssertionError("startup table does not expose both derived groups")
    for species, row in rows.items():
        minimum, maximum = int(row["minimum"]), int(row["maximum"])
        if minimum <= 0 or maximum < 2 * minimum:
            raise AssertionError(f"invalid resolved range for {species}")
        if row["class"] == "direct-source":
            if row["selector"] != "DIRECT_SOURCE" or not row["locked"]:
                raise AssertionError(f"direct-source policy failed for {species}")
            if not float(row["source"]) > 0.0:
                raise AssertionError(f"direct source {species} has no physical rate")
        elif row["class"] == "chemical-only":
            if row["selector"] != "CHEMICAL_ONLY" or row["locked"]:
                raise AssertionError(f"chemical-only policy failed for {species}")
            if (minimum, maximum) != expected_chemical_range:
                raise AssertionError(f"wrong chemical-only range for {species}")
            if float(row["source"]) != 0.0:
                raise AssertionError(f"chemical-only species {species} has a direct rate")
    return rows


def compare_locked_sources(
    reference: dict[str, dict[str, object]],
    alternate: dict[str, dict[str, object]],
) -> None:
    if reference.keys() != alternate.keys():
        raise AssertionError("policy logs contain different species")
    for species, row in reference.items():
        if row["class"] != "direct-source":
            continue
        other = alternate[species]
        for field in ("source", "carrier"):
            left, right = float(row[field]), float(other[field])
            if not math.isclose(left, right, rel_tol=5.0e-14, abs_tol=0.0):
                raise AssertionError(
                    f"daughter target changed locked {species} {field}: {left} != {right}"
                )


def validate_cadence(path: Path, interval: int, required_iterations: list[int]) -> None:
    text = path.read_text(encoding="utf-8", errors="replace")
    match = CADENCE.search(text)
    if match is None or match.group("mode") != "fixed-iteration":
        raise AssertionError("cadence run did not select fixed-iteration adaptation")
    if int(match.group("interval") or -1) != interval:
        raise AssertionError("cadence run reports the wrong interval")
    observed = {int(value) for value in ADAPTATION.findall(text)}
    missing = set(required_iterations) - observed
    if missing:
        raise AssertionError(f"cadence run is missing adaptation iterations {sorted(missing)}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("policy_log", type=Path)
    parser.add_argument("--chemical-range", nargs=2, type=int, default=(50, 200))
    parser.add_argument("--alternate-log", type=Path)
    parser.add_argument("--alternate-chemical-range", nargs=2, type=int, default=(100, 400))
    parser.add_argument("--cadence-log", type=Path)
    parser.add_argument("--cadence-interval", type=int, default=50)
    parser.add_argument("--cadence-iterations", nargs="*", type=int, default=(0, 50, 100))
    args = parser.parse_args()

    reference = validate_policy(args.policy_log, tuple(args.chemical_range))
    if args.alternate_log is not None:
        alternate = validate_policy(
            args.alternate_log, tuple(args.alternate_chemical_range)
        )
        compare_locked_sources(reference, alternate)
    if args.cadence_log is not None:
        validate_cadence(
            args.cadence_log,
            args.cadence_interval,
            list(args.cadence_iterations),
        )
    print(f"PASS Stage-1 runtime population policies species={len(reference)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
