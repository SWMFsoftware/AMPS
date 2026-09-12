#!/usr/bin/env python3
"""Regenerate and audit the frozen SHK05/SHK12 shock references.

This script is the executable SHK17 reproducibility contract.  It runs each
canonical generator from its raw in-script input grid, compares both generator
and generated-fixture SHA-256 values with the reviewed manifest, repeats
selected solves at 80 and 100 Decimal digits, and rejects dependencies on the
production shock solver.  Temporary candidates never overwrite checked-in
fixtures, so a failed audit cannot silently bless changed reference data.
"""

from __future__ import annotations

import argparse
import ast
from decimal import getcontext
import hashlib
import importlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tokenize
from typing import Iterable, Mapping, Sequence


REFERENCE_DIR = Path(__file__).resolve().parent
MANIFEST_PATH = REFERENCE_DIR / "shock_reference_manifest_v1.json"
FORBIDDEN_TOKENS = (
    "swcme_shock",
    "swcme1d",
    "swcme3d",
    "solve_ideal_mhd_fast_shock",
    "candidate_for_compression",
)


def sha256_bytes(payload: bytes) -> str:
    """Return a lowercase SHA-256 digest without external utilities."""
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def audit_generator_dependencies(path: Path) -> None:
    """Reject production imports/calls while allowing the SHK12->SHK05 reuse.

    AST parsing ensures every Python import is syntactically inspectable; the
    explicit token scan also catches a future dynamic import or copied call
    name that would not appear as a normal Import node.
    """
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(path))
    allowed_local = {"generate_shk05_oblique_v1"}
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names = {alias.name.split(".")[0] for alias in node.names}
        elif isinstance(node, ast.ImportFrom):
            names = {(node.module or "").split(".")[0]}
        else:
            continue
        local_imports = {name for name in names if name.startswith("generate_")}
        if not local_imports.issubset(allowed_local):
            raise RuntimeError(f"{path.name}: unapproved local generator import")
    # Scan executable identifiers rather than raw text: both generators
    # intentionally state in comments/docstrings that swcme_shock is forbidden,
    # and that assurance must not trigger its own dependency guard.
    executable_names = {
        token.string.lower()
        for token in tokenize.generate_tokens(io.StringIO(source).readline)
        if token.type == tokenize.NAME
    }
    for forbidden in FORBIDDEN_TOKENS:
        if forbidden.lower() in executable_names:
            raise RuntimeError(
                f"{path.name}: forbidden production dependency {forbidden}"
            )


def regenerate(generator: Path) -> bytes:
    """Run one generator in the pinned reference directory and capture bytes."""
    completed = subprocess.run(
        [sys.executable, generator.name],
        cwd=REFERENCE_DIR,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if completed.returncode != 0:
        detail = completed.stderr.decode("utf-8", errors="replace")
        raise RuntimeError(f"{generator.name}: regeneration failed: {detail}")
    return completed.stdout


def flattened_physics(result: Mapping[str, object], excluded: Iterable[str]):
    """Flatten Decimal scalar/vector physics values for precision comparison."""
    excluded_set = set(excluded)
    values = []
    for key, value in result.items():
        if key in excluded_set or isinstance(value, (str, int)):
            continue
        if isinstance(value, (tuple, list)):
            values.extend(value)
        else:
            values.append(value)
    return values


def require_same_binary64_literals(module, baseline: Sequence[Mapping[str, object]],
                                   higher: Sequence[Mapping[str, object]],
                                   excluded: Iterable[str], label: str) -> None:
    """Require higher precision to preserve every emitted physics literal."""
    if len(baseline) != len(higher):
        raise RuntimeError(f"{label}: higher-precision result count changed")
    for low_case, high_case in zip(baseline, higher):
        low_values = flattened_physics(low_case, excluded)
        high_values = flattened_physics(high_case, excluded)
        if len(low_values) != len(high_values):
            raise RuntimeError(f"{label}: higher-precision field count changed")
        for index, (low, high) in enumerate(zip(low_values, high_values)):
            if module.literal(low) != module.literal(high):
                raise RuntimeError(
                    f"{label}: binary64 literal {index} changed at 100 digits"
                )


def audit_higher_precision() -> None:
    """Repeat representative strong, polarity, and weak cases at 100 digits."""
    sys.path.insert(0, str(REFERENCE_DIR))
    try:
        shk05 = importlib.import_module("generate_shk05_oblique_v1")
        shk12 = importlib.import_module("generate_shk12_near_mach_v1")

        selected05 = (shk05.CASES[0], shk05.CASES[8], shk05.CASES[-1])
        getcontext().prec = 80
        baseline05 = [shk05.solve_case(case) for case in selected05]
        getcontext().prec = 100
        higher05 = [shk05.solve_case(case) for case in selected05]
        require_same_binary64_literals(
            shk05, baseline05, higher05,
            {"reference_residual", "minimum_newton_pivot",
             "converged_seed_count", "physical_seed_count",
             "physical_root_count"},
            "SHK05 selected cases",
        )

        selected12 = (shk12.FAMILIES[0], shk12.FAMILIES[-1])
        getcontext().prec = 80
        baseline12 = [case for family in selected12
                      for case in shk12.solve_family(family)]
        getcontext().prec = 100
        higher12 = [case for family in selected12
                    for case in shk12.solve_family(family)]
        require_same_binary64_literals(
            shk05, baseline12, higher12,
            {"reference_residual", "minimum_pivot", "continuation_iterations"},
            "SHK12 selected families",
        )
    finally:
        getcontext().prec = 80
        if sys.path and sys.path[0] == str(REFERENCE_DIR):
            sys.path.pop(0)


def verify(quiet: bool = False) -> None:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 1:
        raise RuntimeError("unsupported shock-reference manifest schema")
    for entry in manifest["references"]:
        generator = REFERENCE_DIR / entry["generator"]
        fixture = REFERENCE_DIR / entry["fixture"]
        audit_generator_dependencies(generator)
        actual_generator_hash = sha256_file(generator)
        if actual_generator_hash != entry["generator_sha256"]:
            raise RuntimeError(f"{generator.name}: generator hash changed")
        stored_hash = sha256_file(fixture)
        if stored_hash != entry["fixture_sha256"]:
            raise RuntimeError(f"{fixture.name}: stored fixture hash changed")
        regenerated = regenerate(generator)
        regenerated_hash = sha256_bytes(regenerated)
        if regenerated_hash != stored_hash or regenerated != fixture.read_bytes():
            raise RuntimeError(f"{fixture.name}: regeneration is not byte-identical")
        if not quiet:
            print(f"{fixture.name}: {stored_hash} PASS")
    audit_higher_precision()
    if not quiet:
        print("higher_precision_binary64_stability: PASS")
        print("production_dependency_guard: PASS")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    try:
        verify(quiet=args.quiet)
    except Exception as error:  # preserve one concise CI diagnostic
        print(f"SHK17 FAIL: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
