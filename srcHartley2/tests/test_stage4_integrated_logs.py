#!/usr/bin/env python3
"""Strictly validate the correction-on Stage-4 1/4/12-rank matrix.

The executable's focused C++ validators own the authoritative particle and
block state. This checker treats their compact records as a test protocol: an
exit status of zero is insufficient unless the records prove that a wide
individual-correction range and every required MPI/restart/adaptation path
were actually exercised. Optional correction-off unit logs close the identity
limit by requiring both modes to pass independent physical tests with exactly
the same path-coverage signature at every MPI decomposition.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path


EXPECTED_RANKS = {1, 4, 12}
LOCAL_HEADER = "Local-timestep validation: "
LOCAL_RESULT = "Local-timestep validation result: "
LOCAL_AUDIT_PREFIX = "PASS individual-weight local-timestep validation "
ADAPT_PREFIX = "PASS particle-weight adaptation validation "
ADAPT_AUDIT_PREFIX = "PASS individual-weight adaptation validation "
COMPONENT_PREFIX = "PASS component-carrier validation "
CHEMISTRY_RESULT = "Chemistry ensemble validation result: "
CHEMISTRY_REACTION = "PASS chemistry reaction="
CORE_AUDIT_PREFIX = "PASS individual-weight core compatibility "
NTC_AUDIT_PREFIX = "PASS individual-weight NTC validation "
ROUND_OFF = 256.0 * float.fromhex("0x1.0p-52")
FLUX_CHECKSUM_LIMIT = 1.0e-9


def parse_fields(line: str, prefix: str, path: Path) -> dict[str, str]:
    """Parse one whitespace-separated name=value record without ambiguity."""

    if not line.startswith(prefix):
        raise AssertionError(f"{path}: record does not start with {prefix!r}")
    fields: dict[str, str] = {}
    for token in line[len(prefix) :].split():
        if "=" not in token:
            raise AssertionError(f"{path}: malformed validation token {token!r}")
        name, value = token.split("=", 1)
        if not name or not value or name in fields:
            raise AssertionError(f"{path}: invalid or duplicate field {token!r}")
        fields[name] = value
    return fields


def find_one(text: str, prefix: str, path: Path) -> dict[str, str]:
    lines = [
        line.strip()
        for line in text.splitlines()
        if line.strip().startswith(prefix)
    ]
    if len(lines) != 1:
        raise AssertionError(
            f"{path}: expected exactly one {prefix.strip()!r} record, "
            f"found {len(lines)}"
        )
    return parse_fields(lines[0], prefix, path)


def find_one_pass(text: str, pass_prefix: str, path: Path) -> dict[str, str]:
    """Require one status record and reject FAIL or mixed/duplicated status."""

    if not pass_prefix.startswith("PASS "):
        raise AssertionError("internal test error: PASS prefix required")
    fail_prefix = "FAIL " + pass_prefix[len("PASS ") :]
    lines = [
        line.strip()
        for line in text.splitlines()
        if line.strip().startswith((pass_prefix, fail_prefix))
    ]
    if len(lines) != 1 or not lines[0].startswith(pass_prefix):
        raise AssertionError(
            f"{path}: expected exactly one passing {pass_prefix.strip()!r} record"
        )
    return parse_fields(lines[0], pass_prefix, path)


def require_exact_line(text: str, prefix: str, expected: str, path: Path) -> None:
    lines = [
        line.strip()
        for line in text.splitlines()
        if line.strip().startswith(prefix)
    ]
    if lines != [prefix + expected]:
        raise AssertionError(
            f"{path}: expected exactly one {prefix + expected!r} final record"
        )


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


def positive(fields: dict[str, str], name: str, path: Path) -> int:
    value = integer(fields, name, path)
    if value <= 0:
        raise AssertionError(f"{path}: required check {name!r} did not execute")
    return value


def colon_range(
    fields: dict[str, str], name: str, path: Path, *, integral: bool = False,
    allow_equal: bool = False,
) -> tuple[float, float]:
    try:
        lower_text, upper_text = fields[name].split(":", 1)
        lower = float(lower_text)
        upper = float(upper_text)
    except (KeyError, ValueError) as error:
        raise AssertionError(f"{path}: missing or invalid range {name!r}") from error
    ordered = lower <= upper if allow_equal else lower < upper
    if (
        not math.isfinite(lower)
        or not math.isfinite(upper)
        or lower <= 0.0
        or not ordered
    ):
        raise AssertionError(f"{path}: invalid positive ordered range {name}={lower}:{upper}")
    if integral and (not lower.is_integer() or not upper.is_integer()):
        raise AssertionError(f"{path}: {name!r} must contain integer population bounds")
    return lower, upper


def require_correction_on(fields: dict[str, str], path: Path) -> None:
    if fields.get("individual-correction") != "on":
        raise AssertionError(f"{path}: individual-correction=on evidence is absent")


def require_correction_mode(fields: dict[str, str], path: Path) -> None:
    if fields.get("correction-mode") != "ON":
        raise AssertionError(f"{path}: correction-mode=ON evidence is absent")


def audit_correction_range(
    fields: dict[str, str], path: Path
) -> tuple[float, float, int]:
    minimum = floating(fields, "correction-min", path)
    maximum = floating(fields, "correction-max", path)
    values = integer(fields, "correction-values", path)
    if not (0.0 < minimum < 1.0 < maximum):
        raise AssertionError(f"{path}: correction extrema do not span unity")
    if maximum / minimum < 16.0:
        raise AssertionError(f"{path}: correction extrema are too narrow")
    if values < 5:
        raise AssertionError(f"{path}: fewer than five correction values were tested")
    error = floating(fields, "max-correction-error", path)
    if error < 0.0 or error > ROUND_OFF:
        raise AssertionError(f"{path}: individual correction was not preserved")
    if integer(fields, "partial-mutations", path) != 0:
        raise AssertionError(f"{path}: a rejected operation partially mutated state")
    return minimum, maximum, values


def validate_local(
    path: Path, expected_ranks: int, expected_pattern: str = "wide",
    expected_mode: str = "ON",
) -> dict[str, object]:
    text = path.read_text(encoding="utf-8", errors="replace")
    header = find_one(text, LOCAL_HEADER, path)
    audit = find_one_pass(text, LOCAL_AUDIT_PREFIX, path)
    require_exact_line(text, LOCAL_RESULT, "PASS", path)
    if integer(header, "ranks", path) != expected_ranks:
        raise AssertionError(f"{path}: local-timestep rank evidence is incorrect")
    if integer(header, "particles/case", path) <= 0:
        raise AssertionError(f"{path}: local-timestep particle count is not positive")
    if header.get("correction-pattern") != expected_pattern or audit.get(
        "correction-pattern"
    ) != expected_pattern:
        raise AssertionError(f"{path}: wrong correction pattern")
    header_range = colon_range(
        header, "correction-range", path, allow_equal=expected_pattern == "unit"
    )
    if header.get("deterministic-internal-energy") != "yes":
        raise AssertionError(f"{path}: deterministic internal-energy test is absent")
    if audit.get("correction-mode") != expected_mode:
        raise AssertionError(
            f"{path}: correction-mode={expected_mode} evidence is absent"
        )
    if integer(audit, "mpi-ranks", path) != expected_ranks:
        raise AssertionError(f"{path}: local correction audit reports wrong ranks")
    if expected_pattern == "wide":
        minimum, maximum, values = audit_correction_range(audit, path)
    elif expected_pattern == "unit":
        minimum = floating(audit, "correction-min", path)
        maximum = floating(audit, "correction-max", path)
        values = integer(audit, "correction-values", path)
        if (minimum, maximum, values) != (1.0, 1.0, 1):
            raise AssertionError(f"{path}: unit-correction identity case is absent")
        error = floating(audit, "max-correction-error", path)
        if error < 0.0 or error > ROUND_OFF:
            raise AssertionError(f"{path}: unit correction was not preserved")
        if integer(audit, "partial-mutations", path) != 0:
            raise AssertionError(f"{path}: a rejected operation partially mutated state")
    else:
        raise AssertionError(f"internal test error: unknown pattern {expected_pattern!r}")
    correction_range = (minimum, maximum)
    if correction_range != header_range:
        raise AssertionError(f"{path}: local header/audit correction ranges differ")
    if integer(header, "correction-bins", path) != values:
        raise AssertionError(f"{path}: local header/audit correction counts differ")
    flux_error = floating(audit, "max-flux-checksum-error", path)
    if flux_error < 0.0 or flux_error > FLUX_CHECKSUM_LIMIT:
        raise AssertionError(f"{path}: post-move physical flux checksum is inconsistent")
    for name in (
        "local-face-checks",
        "multiface-checks",
        "correction-preservation-checks",
        "clone-correction-checks",
        "internal-energy-checks",
    ):
        positive(audit, name, path)
    mpi_faces = integer(audit, "mpi-face-checks", path)
    if (expected_ranks == 1 and mpi_faces != 0) or (
        expected_ranks > 1 and mpi_faces <= 0
    ):
        raise AssertionError(f"{path}: inconsistent MPI-face check count")
    # Return the exact path-coverage signature as well as the correction range.
    # The optional correction-on/off identity audit below compares this
    # signature at each MPI decomposition.  Both arms still have independent
    # statistical acceptance tests in the C++ validator; equality here proves
    # that neither mode obtained an easier mover/MPI/internal-energy test.
    coverage_names = (
        "local-face-checks",
        "mpi-face-checks",
        "multiface-checks",
        "correction-preservation-checks",
        "clone-correction-checks",
        "internal-energy-checks",
    )
    return {
        "correction_range": correction_range,
        "coverage": tuple(integer(audit, name, path) for name in coverage_names),
        "flux_checksum_error": flux_error,
    }


def validate_component(fields: dict[str, str], path: Path, ranks: int) -> None:
    if integer(fields, "mpi-ranks", path) != ranks:
        raise AssertionError(f"{path}: component record reports the wrong rank count")
    owners = positive(fields, "owner-blocks", path)
    if owners != ranks:
        raise AssertionError(f"{path}: every MPI rank must own a validation block")
    for name in (
        "component-population-checks",
        "empty-block-checks",
        "two-ended-checks",
        "locked-source-checks",
        "reaction-allocation-checks",
        "split-conservation-checks",
        "roulette-conservation-checks",
        "capacity-rollback-checks",
    ):
        if integer(fields, name, path) < owners:
            raise AssertionError(f"{path}: {name} did not execute on every owner")
    ghosts = integer(fields, "ghost-weight-checks", path)
    if (ranks == 1 and ghosts != 0) or (ranks > 1 and ghosts <= 0):
        raise AssertionError(f"{path}: inconsistent ghost-weight check count")
    if floating(fields, "max-face-multiplicity", path) > 1.0 + ROUND_OFF:
        raise AssertionError(f"{path}: face multiplicity exceeds one")
    if integer(fields, "max-reaction-allocation", path) > integer(
        fields, "reaction-limit", path
    ):
        raise AssertionError(f"{path}: reaction allocation exceeds its hard limit")
    if integer(fields, "partial-mutations", path) != 0:
        raise AssertionError(f"{path}: rejected adaptation partially mutated state")


def validate_adaptation(path: Path, expected_ranks: int) -> dict[str, object]:
    text = path.read_text(encoding="utf-8", errors="replace")
    legacy = find_one_pass(text, ADAPT_PREFIX, path)
    audit = find_one_pass(text, ADAPT_AUDIT_PREFIX, path)
    component = find_one_pass(text, COMPONENT_PREFIX, path)
    core = find_one_pass(text, CORE_AUDIT_PREFIX, path)
    ntc = find_one_pass(text, NTC_AUDIT_PREFIX, path)
    require_correction_on(legacy, path)
    population_range = colon_range(legacy, "range", path, integral=True)
    validate_component(component, path, expected_ranks)

    # This compact correction-specific record is deliberately separate from
    # range=MIN:MAX above, which denotes the requested population per cell and
    # must never be misreported as an individual-correction range.
    require_correction_mode(audit, path)
    if integer(audit, "mpi-ranks", path) != expected_ranks:
        raise AssertionError(f"{path}: correction audit reports the wrong ranks")
    minimum, maximum, _ = audit_correction_range(audit, path)
    correction_range = (minimum, maximum)
    owners = integer(component, "owner-blocks", path)
    required_counts = {
        "split-correction-checks": owners,
        "roulette-correction-checks": owners,
        "emergency-correction-checks": 3 * owners,
        "corrected-moment-checks": 2 * owners,
        "capacity-rollback-checks": owners,
    }
    for name, expected in required_counts.items():
        if integer(audit, name, path) != expected:
            raise AssertionError(
                f"{path}: {name}={audit.get(name)!r}, expected {expected}"
            )
    ghosts = integer(audit, "ghost-weight-publication-checks", path)
    if (expected_ranks == 1 and ghosts != 0) or (
        expected_ranks > 1 and ghosts <= 0
    ):
        raise AssertionError(f"{path}: inconsistent correction ghost-publication count")
    require_correction_mode(core, path)
    if integer(core, "mpi-ranks", path) != expected_ranks:
        raise AssertionError(f"{path}: core compatibility record reports wrong ranks")
    for name in (
        "sampling-checks",
        "collision-checks",
        "boundary-checks",
        "restart-checks",
        "migration-checks",
        "collision-probability-checks",
    ):
        positive(core, name, path)
    probability_min = floating(core, "collision-probability-min", path)
    probability_max = floating(core, "collision-probability-max", path)
    if not 0.0 <= probability_min <= probability_max <= 1.0:
        raise AssertionError(f"{path}: collision update probabilities leave [0,1]")
    target_updates = positive(ntc, "target-updates", path)
    if integer(ntc, "updates", path) < target_updates:
        raise AssertionError(f"{path}: NTC work target was not reached")
    positive(ntc, "cycles", path)
    positive(ntc, "trials", path)
    positive(ntc, "accepted", path)
    if integer(ntc, "invalid-probabilities", path) != 0:
        raise AssertionError(f"{path}: NTC emitted an invalid update probability")
    ntc_probability_min = floating(ntc, "probability-min", path)
    ntc_probability_max = floating(ntc, "probability-max", path)
    if not 0.0 <= ntc_probability_min <= ntc_probability_max <= 1.0:
        raise AssertionError(f"{path}: focused NTC probability leaves [0,1]")
    initial_separation = floating(ntc, "initial-separation", path)
    final_separation = floating(ntc, "final-separation", path)
    if not 0.0 <= final_separation < initial_separation:
        raise AssertionError(f"{path}: cross-species NTC did not relax drift")
    return {
        "population_range": population_range,
        "correction_range": correction_range,
    }


def validate_chemistry(path: Path) -> dict[str, int]:
    text = path.read_text(encoding="utf-8", errors="replace")
    header = find_one(text, "Chemistry ensemble validation: ", path)
    require_correction_on(header, path)
    result_lines = [
        line.strip()
        for line in text.splitlines()
        if line.strip().startswith(CHEMISTRY_RESULT)
    ]
    pass_prefix = CHEMISTRY_RESULT + "PASS "
    if len(result_lines) != 1 or not result_lines[0].startswith(pass_prefix):
        raise AssertionError(f"{path}: chemistry ensemble final result did not pass")
    result = parse_fields(result_lines[0], pass_prefix, path)
    reactions = integer(result, "neutral-reactions", path)
    unequal = integer(result, "unequal-weight/timestep-reactions", path)
    sinks = integer(result, "sink-reactions", path)
    if reactions <= 0 or unequal <= 0 or sinks < 0:
        raise AssertionError(f"{path}: chemistry ensemble coverage is incomplete")
    reaction_lines = [
        line.strip()
        for line in text.splitlines()
        if line.strip().startswith((CHEMISTRY_REACTION, "FAIL chemistry reaction="))
    ]
    if len(reaction_lines) != reactions or any(
        not line.startswith(CHEMISTRY_REACTION) for line in reaction_lines
    ):
        raise AssertionError(
            f"{path}: expected {reactions} PASS reaction records, found "
            f"{len(reaction_lines)}"
        )
    for line in reaction_lines:
        fields = parse_fields(line, "PASS chemistry ", path)
        positive(fields, "events", path)
        positive(fields, "products", path)
    return {"reactions": reactions, "unequal": unequal, "sinks": sinks}


def collect_logs(records: list[list[str]], option: str) -> dict[int, Path]:
    supplied: dict[int, Path] = {}
    for rank_text, path_text in records:
        try:
            ranks = int(rank_text)
        except ValueError as error:
            raise AssertionError(f"{option}: invalid MPI rank {rank_text!r}") from error
        if ranks in supplied:
            raise AssertionError(f"{option}: duplicate log for {ranks} ranks")
        path = Path(path_text)
        if not path.is_file():
            raise AssertionError(f"{option}: log is not a regular file: {path}")
        supplied[ranks] = path
    if set(supplied) != EXPECTED_RANKS:
        raise AssertionError(f"{option}: requires exactly 1/4/12-rank logs")
    return supplied


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for option, destination, help_text in (
        ("--local-log", "local_log", "local-timestep validator log"),
        ("--unit-local-log", "unit_local_log", "unit-correction local-timestep log"),
        ("--adaptation-log", "adaptation_log", "weight-adaptation validator log"),
        ("--chemistry-log", "chemistry_log", "chemistry ensemble validator log"),
    ):
        parser.add_argument(
            option,
            dest=destination,
            action="append",
            nargs=2,
            required=True,
            metavar=("MPI_RANKS", "PATH"),
            help=help_text + "; supply exactly 1, 4, and 12 ranks",
        )
    parser.add_argument(
        "--off-unit-local-log",
        action="append",
        nargs=2,
        metavar=("MPI_RANKS", "PATH"),
        help=(
            "optional correction-off unit-correction reference log; when used, "
            "supply exactly 1, 4, and 12 ranks and compare it with the "
            "correction-on all-c=1 arm"
        ),
    )
    arguments = parser.parse_args()
    local_paths = collect_logs(arguments.local_log, "--local-log")
    unit_local_paths = collect_logs(arguments.unit_local_log, "--unit-local-log")
    adaptation_paths = collect_logs(arguments.adaptation_log, "--adaptation-log")
    chemistry_paths = collect_logs(arguments.chemistry_log, "--chemistry-log")
    off_unit_paths = None
    if arguments.off_unit_local_log:
        off_unit_paths = collect_logs(
            arguments.off_unit_local_log, "--off-unit-local-log"
        )

    local = {rank: validate_local(local_paths[rank], rank) for rank in EXPECTED_RANKS}
    unit_local = {
        rank: validate_local(unit_local_paths[rank], rank, "unit")
        for rank in EXPECTED_RANKS
    }
    adaptation = {
        rank: validate_adaptation(adaptation_paths[rank], rank)
        for rank in EXPECTED_RANKS
    }
    chemistry = {
        rank: validate_chemistry(chemistry_paths[rank]) for rank in EXPECTED_RANKS
    }
    off_unit_local = None
    if off_unit_paths is not None:
        off_unit_local = {
            rank: validate_local(off_unit_paths[rank], rank, "unit", "OFF")
            for rank in EXPECTED_RANKS
        }

    # All decompositions must exercise the same deliberate correction range,
    # population policy, and reaction network. A narrow/easier setup in one
    # subtest may not masquerade as complete Stage-4 coverage.
    correction_ranges = {
        local[rank]["correction_range"] for rank in EXPECTED_RANKS
    } | {adaptation[rank]["correction_range"] for rank in EXPECTED_RANKS}
    if len(correction_ranges) != 1:
        raise AssertionError(
            f"correction evidence differs across validators/ranks: {correction_ranges}"
        )
    if {
        unit_local[rank]["correction_range"] for rank in EXPECTED_RANKS
    } != {(1.0, 1.0)}:
        raise AssertionError("unit-correction identity evidence differs across MPI ranks")
    if len({adaptation[rank]["population_range"] for rank in EXPECTED_RANKS}) != 1:
        raise AssertionError("adaptation population range differs across MPI ranks")
    if len(
        {
            (row["reactions"], row["unequal"], row["sinks"])
            for row in chemistry.values()
        }
    ) != 1:
        raise AssertionError("chemistry ensemble coverage differs across MPI ranks")

    if off_unit_local is not None:
        for rank in EXPECTED_RANKS:
            if unit_local[rank]["coverage"] != off_unit_local[rank]["coverage"]:
                raise AssertionError(
                    "correction-on/off unit-correction path coverage differs at "
                    f"{rank} MPI ranks"
                )
        print(
            "PASS Stage-4 all-c=1 correction identity "
            + " ".join(f"{rank}rank=pass" for rank in (1, 4, 12))
        )

    print(
        "PASS Stage-4 correction-on integrated matrix "
        + " ".join(f"{rank}rank=pass" for rank in (1, 4, 12))
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
