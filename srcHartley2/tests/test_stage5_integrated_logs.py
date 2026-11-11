#!/usr/bin/env python3
"""Strict checker for the Stage-5 1/4/12-rank runtime archive.

An MPI command returning zero is not sufficient evidence.  This checker
requires the focused controller record to prove that keep, split, roulette,
capacity rollback, restart, direct-source identity, and empty-list paths all
executed.  It reuses the Stage-4 chemistry protocol checker so every generated
reaction still has a passing exact-daughter ensemble record.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from test_stage4_integrated_logs import parse_fields, validate_chemistry


EXPECTED_RANKS = {1, 4, 12}
PASS_PREFIX = "PASS individual-correction window validation "
FAIL_PREFIX = "FAIL individual-correction window validation "
PARTICLES_PER_RANK = 65536


def integer(fields: dict[str, str], name: str, path: Path) -> int:
    try:
        return int(fields[name])
    except (KeyError, ValueError) as error:
        raise AssertionError(f"{path}: missing or invalid integer {name!r}") from error


def floating(fields: dict[str, str], name: str, path: Path) -> float:
    try:
        value = float(fields[name])
    except (KeyError, ValueError) as error:
        raise AssertionError(f"{path}: missing or invalid number {name!r}") from error
    if not math.isfinite(value):
        raise AssertionError(f"{path}: non-finite field {name!r}")
    return value


def validate_window(path: Path, ranks: int) -> tuple[float, float, float, int]:
    text = path.read_text(encoding="utf-8", errors="replace")
    records = [
        line.strip()
        for line in text.splitlines()
        if line.strip().startswith((PASS_PREFIX, FAIL_PREFIX))
    ]
    if len(records) != 1 or not records[0].startswith(PASS_PREFIX):
        raise AssertionError(f"{path}: expected one passing Stage-5 window record")
    fields = parse_fields(records[0], PASS_PREFIX, path)
    if fields.get("correction-mode") != "ON":
        raise AssertionError(f"{path}: correction-mode=ON evidence is absent")
    if integer(fields, "mpi-ranks", path) != ranks:
        raise AssertionError(f"{path}: MPI rank evidence is incorrect")
    if integer(fields, "owner-blocks", path) != ranks:
        raise AssertionError(f"{path}: owner-only traversal evidence is incorrect")

    lower = floating(fields, "lower", path)
    target = floating(fields, "target", path)
    upper = floating(fields, "upper", path)
    if not 0.0 < lower < target < upper:
        raise AssertionError(f"{path}: invalid correction-window ordering")
    maximum_splits = integer(fields, "max-split-particles", path)
    if maximum_splits < 2:
        raise AssertionError(f"{path}: split hard limit was not configured")

    for check in (
        "keep-checks",
        "split-checks",
        "roulette-checks",
        "capacity-rollback-checks",
        "restart-checks",
        "direct-source-checks",
        "empty-checks",
    ):
        if integer(fields, check, path) != 1:
            raise AssertionError(f"{path}: {check} did not execute exactly once")

    candidates = integer(fields, "roulette-candidates", path)
    survivors = integer(fields, "roulette-survivors", path)
    expected_candidates = PARTICLES_PER_RANK * ranks
    if candidates != expected_candidates:
        raise AssertionError(f"{path}: roulette candidate count is incomplete")
    probability = (0.5 * lower) / target
    expectation = candidates * probability
    sigma = math.sqrt(candidates * probability * (1.0 - probability))
    if abs(survivors - expectation) > 7.0 * sigma:
        raise AssertionError(f"{path}: roulette survivors violate the 7-sigma bound")
    return lower, target, upper, maximum_splits


def collect(records: list[list[str]], option: str) -> dict[int, Path]:
    result: dict[int, Path] = {}
    for rank_text, path_text in records:
        try:
            rank = int(rank_text)
        except ValueError as error:
            raise AssertionError(f"{option}: invalid MPI rank {rank_text!r}") from error
        path = Path(path_text)
        if rank in result or not path.is_file():
            raise AssertionError(f"{option}: duplicate rank or missing log for {rank}")
        result[rank] = path
    if set(result) != EXPECTED_RANKS:
        raise AssertionError(f"{option}: expected logs for ranks 1, 4, and 12")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--window-log", nargs=2, action="append", default=[])
    parser.add_argument("--chemistry-log", nargs=2, action="append", default=[])
    arguments = parser.parse_args()
    window_paths = collect(arguments.window_log, "--window-log")
    chemistry_paths = collect(arguments.chemistry_log, "--chemistry-log")

    policies = {
        rank: validate_window(window_paths[rank], rank) for rank in EXPECTED_RANKS
    }
    chemistry = {
        rank: validate_chemistry(chemistry_paths[rank]) for rank in EXPECTED_RANKS
    }
    if len(set(policies.values())) != 1:
        raise AssertionError("resolved correction-window policy differs across MPI ranks")
    if len(
        {
            (row["reactions"], row["unequal"], row["sinks"])
            for row in chemistry.values()
        }
    ) != 1:
        raise AssertionError("chemistry coverage differs across MPI ranks")

    print(
        "PASS Stage-5 correction-window integrated matrix "
        + " ".join(f"{rank}rank=pass" for rank in (1, 4, 12))
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
