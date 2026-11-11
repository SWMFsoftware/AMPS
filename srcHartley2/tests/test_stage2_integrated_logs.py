#!/usr/bin/env python3
"""Validate the 1/4/12-rank Stage-2 component-carrier runtime records.

The C++ validator performs the physical/list checks because only it can read
authoritative AMPS owner lists, boundary copies, particle state, and buffer
capacity.  This script deliberately does not infer those properties from a
science output file.  Instead it treats the validator's single summary record
as a stable test protocol and rejects a run in which any required subtest was
omitted.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path


STAGE2_PREFIX = "PASS component-carrier validation "
LEGACY_PREFIX = "PASS particle-weight adaptation validation "


def parse_fields(line: str, prefix: str) -> dict[str, str]:
    """Parse a whitespace-separated ``name=value`` validation record."""

    if not line.startswith(prefix):
        raise AssertionError(f"validation line does not start with {prefix!r}")
    fields: dict[str, str] = {}
    for token in line[len(prefix) :].split():
        if "=" not in token:
            raise AssertionError(f"malformed validation token {token!r}")
        name, value = token.split("=", 1)
        if not name or not value or name in fields:
            raise AssertionError(f"invalid or duplicate validation field {token!r}")
        fields[name] = value
    return fields


def find_one_record(text: str, prefix: str, path: Path) -> dict[str, str]:
    matches = [line.strip() for line in text.splitlines() if line.strip().startswith(prefix)]
    if len(matches) != 1:
        raise AssertionError(
            f"expected exactly one {prefix.strip()!r} record in {path}, found {len(matches)}"
        )
    return parse_fields(matches[0], prefix)


def integer(fields: dict[str, str], name: str) -> int:
    try:
        return int(fields[name])
    except KeyError as error:
        raise AssertionError(f"missing Stage-2 field {name!r}") from error
    except ValueError as error:
        raise AssertionError(f"Stage-2 field {name!r} is not an integer") from error


def floating(fields: dict[str, str], name: str) -> float:
    try:
        value = float(fields[name])
    except KeyError as error:
        raise AssertionError(f"missing Stage-2 field {name!r}") from error
    except ValueError as error:
        raise AssertionError(f"Stage-2 field {name!r} is not a number") from error
    if not math.isfinite(value):
        raise AssertionError(f"Stage-2 field {name!r} is not finite")
    return value


def validate_log(path: Path, expected_ranks: int) -> dict[str, int | float]:
    text = path.read_text(encoding="utf-8", errors="replace")
    fields = find_one_record(text, STAGE2_PREFIX, path)
    # The older integrated validator remains a regression prerequisite.  A
    # Stage-2-only PASS may not hide a failure in emergency, face, chemistry,
    # or owner/ghost behavior that Stage 1 already validated.
    find_one_record(text, LEGACY_PREFIX, path)

    ranks = integer(fields, "mpi-ranks")
    owners = integer(fields, "owner-blocks")
    if ranks != expected_ranks:
        raise AssertionError(f"{path}: reports {ranks} ranks, expected {expected_ranks}")
    if owners != expected_ranks:
        raise AssertionError(
            f"{path}: Stage-2 validator must exercise one owner block on every rank"
        )

    # Each named check closes a distinct Stage-2 exit criterion.  Counts are
    # emitted instead of booleans so an MPI reduction bug cannot turn one
    # rank's success into an apparently global pass.
    per_owner = (
        "component-population-checks",
        "empty-block-checks",
        "two-ended-checks",
        "locked-source-checks",
        "reaction-allocation-checks",
        "capacity-rollback-checks",
    )
    for name in per_owner:
        if integer(fields, name) < owners:
            raise AssertionError(f"{path}: {name} did not execute on every owner rank")

    for name in ("split-conservation-checks", "roulette-conservation-checks"):
        if integer(fields, name) <= 0:
            raise AssertionError(f"{path}: {name} was not executed")

    ghosts = integer(fields, "ghost-weight-checks")
    if (expected_ranks == 1 and ghosts != 0) or (expected_ranks > 1 and ghosts <= 0):
        raise AssertionError(f"{path}: inconsistent ghost-weight check count {ghosts}")

    face_multiplicity = floating(fields, "max-face-multiplicity")
    if face_multiplicity > 1.0 + 256.0 * float.fromhex("0x1.0p-52"):
        raise AssertionError(f"{path}: face multiplicity exceeds one: {face_multiplicity}")

    reaction_allocation = integer(fields, "max-reaction-allocation")
    reaction_limit = integer(fields, "reaction-limit")
    if reaction_limit <= 0 or reaction_allocation > reaction_limit:
        raise AssertionError(
            f"{path}: reaction allocation {reaction_allocation} exceeds {reaction_limit}"
        )
    if integer(fields, "partial-mutations") != 0:
        raise AssertionError(f"{path}: capacity rejection partially mutated a transaction")

    return {
        "owners": owners,
        "ghosts": ghosts,
        "face": face_multiplicity,
        "reaction": reaction_allocation,
        "limit": reaction_limit,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--log",
        action="append",
        nargs=2,
        metavar=("MPI_RANKS", "PATH"),
        required=True,
        help="one Stage-2 validator log; supply the 1-, 4-, and 12-rank logs",
    )
    args = parser.parse_args()

    supplied: dict[int, Path] = {}
    for rank_text, path_text in args.log:
        ranks = int(rank_text)
        if ranks in supplied:
            raise AssertionError(f"duplicate log for {ranks} MPI ranks")
        supplied[ranks] = Path(path_text)
    if set(supplied) != {1, 4, 12}:
        raise AssertionError("Stage-2 integrated validation requires exactly 1/4/12 ranks")

    results = {ranks: validate_log(path, ranks) for ranks, path in supplied.items()}
    # The synthetic reaction setup is deterministic and independent of domain
    # decomposition.  A rank-dependent maximum would expose owner/ghost or
    # collective-planning drift even when each local value remains under cap.
    reaction_pairs = {(row["reaction"], row["limit"]) for row in results.values()}
    if len(reaction_pairs) != 1:
        raise AssertionError(f"reaction allocation result varies with MPI ranks: {reaction_pairs}")

    print(
        "PASS Stage-2 integrated component carriers "
        + " ".join(
            f"{ranks}rank-owners={results[ranks]['owners']}"
            f"-ghosts={results[ranks]['ghosts']}"
            for ranks in (1, 4, 12)
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
