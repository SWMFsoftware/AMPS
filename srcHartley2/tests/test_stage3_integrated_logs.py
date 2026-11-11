#!/usr/bin/env python3
"""Strictly validate the 1/4/12-rank Stage-3 runtime summaries."""

from __future__ import annotations

import argparse
import math
from pathlib import Path


PREFIX = "PASS systematic-residual validation "


def parse_record(path: Path) -> dict[str, str]:
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = [line.strip() for line in text.splitlines() if line.startswith(PREFIX)]
    if len(lines) != 1:
        raise AssertionError(
            f"{path}: expected exactly one Stage-3 PASS record, found {len(lines)}"
        )
    fields: dict[str, str] = {}
    for token in lines[0][len(PREFIX) :].split():
        if "=" not in token:
            raise AssertionError(f"{path}: malformed token {token!r}")
        name, value = token.split("=", 1)
        if not name or not value or name in fields:
            raise AssertionError(f"{path}: invalid or duplicate token {token!r}")
        fields[name] = value
    return fields


def integer(fields: dict[str, str], name: str, path: Path) -> int:
    try:
        return int(fields[name])
    except (KeyError, ValueError) as error:
        raise AssertionError(f"{path}: missing or invalid integer {name!r}") from error


def validate(path: Path, expected_ranks: int) -> dict[str, int | float]:
    fields = parse_record(path)
    if fields.get("mode") != "SYSTEMATIC":
        raise AssertionError(f"{path}: validator did not exercise SYSTEMATIC mode")
    if integer(fields, "mpi-ranks", path) != expected_ranks:
        raise AssertionError(f"{path}: MPI-rank count does not match invocation")

    owners = integer(fields, "owner-blocks", path)
    accumulators = integer(fields, "accumulators", path)
    if owners <= 0 or accumulators <= owners:
        raise AssertionError(f"{path}: owner accumulator state was not populated")
    if integer(fields, "same-carrier-packet-checks", path) <= 0:
        raise AssertionError(f"{path}: common-phase packet behavior was not checked")
    for name in ("plan-retry-checks", "carrier-epoch-reset-checks"):
        if integer(fields, name, path) != owners:
            raise AssertionError(f"{path}: {name} did not execute on every owner block")
    if integer(fields, "first-use-rollback-checks", path) != owners:
        raise AssertionError(
            f"{path}: first-use planning was not checked on every owner block"
        )
    for name in ("migration-checks", "restart-continuation-checks"):
        if integer(fields, name, path) != owners:
            raise AssertionError(f"{path}: {name} did not preserve every owner block")
    if integer(fields, "restart-particles", path) != expected_ranks:
        raise AssertionError(
            f"{path}: local-weight restart did not restore one particle per rank"
        )
    if integer(fields, "restart-particle-checksum", path) != 1:
        raise AssertionError(f"{path}: particle restart checksum did not continue")
    ghost_sync = integer(fields, "restart-ghost-sync-checks", path)
    if expected_ranks == 1 and ghost_sync != 0:
        raise AssertionError(f"{path}: a one-rank run cannot synchronize MPI ghosts")
    if expected_ranks > 1 and ghost_sync <= 0:
        raise AssertionError(f"{path}: restored owner state was not checked on ghosts")
    migrated = integer(fields, "migrated-nonzero-states", path)
    if expected_ranks == 1 and migrated != 0:
        raise AssertionError(f"{path}: a one-rank run cannot migrate MPI ownership")
    if expected_ranks > 1 and migrated <= 0:
        raise AssertionError(f"{path}: no nonzero state crossed an MPI ownership boundary")
    if integer(fields, "owner-only-state-checks", path) != owners:
        raise AssertionError(f"{path}: owner-only state guard was not exercised")
    if integer(fields, "partial-mutations", path) != 0:
        raise AssertionError(f"{path}: a rejected transaction mutated residual state")
    if integer(fields, "center-cell-output", path) != 1:
        raise AssertionError(f"{path}: center-cell diagnostic path was not exercised")
    if integer(fields, "sampling-report-output", path) != 1:
        raise AssertionError(f"{path}: systematic sampling report was not exercised")
    try:
        prefix_error = float(fields["max-prefix-error"])
    except (KeyError, ValueError) as error:
        raise AssertionError(f"{path}: invalid max-prefix-error") from error
    if not math.isfinite(prefix_error) or not (0.0 <= prefix_error < 1.0):
        raise AssertionError(f"{path}: systematic prefix error is not strictly below one")
    return {
        "owners": owners,
        "accumulators": accumulators,
        "packet_checks": integer(fields, "same-carrier-packet-checks", path),
        "rollback_checks": integer(fields, "plan-retry-checks", path),
        "first_use_checks": integer(fields, "first-use-rollback-checks", path),
        "carrier_checks": integer(fields, "carrier-epoch-reset-checks", path),
        "migration_checks": integer(fields, "migration-checks", path),
        "restart_checks": integer(fields, "restart-continuation-checks", path),
        "partial_mutations": integer(fields, "partial-mutations", path),
        "prefix_error": prefix_error,
        "migrated": migrated,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--log", action="append", nargs=2, required=True,
        metavar=("MPI_RANKS", "PATH")
    )
    arguments = parser.parse_args()
    supplied: dict[int, Path] = {}
    for rank_text, path_text in arguments.log:
        ranks = int(rank_text)
        if ranks in supplied:
            raise AssertionError(f"duplicate Stage-3 log for {ranks} ranks")
        supplied[ranks] = Path(path_text)
    if set(supplied) != {1, 4, 12}:
        raise AssertionError("Stage-3 validation requires exactly 1/4/12-rank logs")
    results = {rank: validate(path, rank) for rank, path in supplied.items()}
    invariant_names = (
        "owners",
        "accumulators",
        "packet_checks",
        "rollback_checks",
        "first_use_checks",
        "carrier_checks",
        "migration_checks",
        "restart_checks",
        "partial_mutations",
        "prefix_error",
    )
    for name in invariant_names:
        values = {results[rank][name] for rank in (1, 4, 12)}
        if len(values) != 1:
            raise AssertionError(
                f"Stage-3 decomposition-invariant field {name!r} changed "
                f"across MPI rank counts: {values}"
            )
    print(
        "PASS Stage-3 integrated systematic residual "
        + " ".join(
            f"{rank}rank-owners={results[rank]['owners']}"
            f"-migrated={results[rank]['migrated']}"
            for rank in (1, 4, 12)
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
