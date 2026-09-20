#!/usr/bin/env python3
"""Audit the source-level contracts introduced by Stage 3.

This test is intentionally dependency-light.  It does not claim that a
configured AMPS executable has been linked or that an observational campaign
has passed.  Instead, it protects the release-structure decisions that can be
verified from a clean source extraction:

* retired production sources and generated products are rejected by one
  AMPS-level manifest checker;
* every OV/EV case is registered with an explicit scientific role and input
  policy; and
* srcSEP3D exposes one immutable, proton-only AMPS species binding rather than
  silently copying one species record over an arbitrary AMPS species table.

Keeping these assertions in one small gate prevents documentation-only fixes:
the test inspects the actual runner, make targets, configuration contract, and
production source adapter that a release build consumes.
"""

from __future__ import annotations

import json
from pathlib import Path
import re
import sys


SRCSEP = Path(__file__).resolve().parents[1]
AMPS_ROOT = SRCSEP.parent
SRCSEP3D = AMPS_ROOT / "srcSEP3D"


class ContractFailure(RuntimeError):
    """Report one actionable Stage-3 source-contract violation."""


def require(condition: bool, message: str) -> None:
    """Raise a stable diagnostic instead of relying on disabled assertions."""
    if not condition:
        raise ContractFailure(message)


def read(path: Path) -> str:
    """Read UTF-8 source and retain the exact missing path in diagnostics."""
    try:
        return path.read_text(encoding="utf-8")
    except OSError as error:
        raise ContractFailure(f"cannot read required Stage-3 file {path}: {error}") from error


def check_release_hygiene() -> None:
    """Require the single release checker and physical retirement of mover.cpp."""
    checker = AMPS_ROOT / "tools" / "sep_package_hygiene.py"
    require(checker.is_file(), "missing AMPS-level tools/sep_package_hygiene.py")
    require(not (SRCSEP / "mover.cpp").exists(),
            "retired srcSEP/mover.cpp is still present")

    for relative in (
            "srcSEP/SOURCE_MANIFEST.json",
            "srcSEP3D/SOURCE_MANIFEST.json",
            "src/models/sep_common/SOURCE_MANIFEST.json",
            "src/models/swcme/SOURCE_MANIFEST.json"):
        manifest = AMPS_ROOT / relative
        require(manifest.is_file(), f"missing release allowlist {relative}")
        try:
            payload = json.loads(read(manifest))
        except json.JSONDecodeError as error:
            raise ContractFailure(f"invalid JSON manifest {manifest}: {error}") from error
        require(payload.get("schema") == "amps-sep-source-manifest-v1",
                f"unsupported manifest schema in {relative}")
        require("SOURCE_MANIFEST.json" in payload.get("support_globs", []),
                f"{relative} does not classify its own release manifest")


def check_validation_registration() -> None:
    """Require all seven Stage-3 scientific cases and their release semantics."""
    registry_path = SRCSEP / "validation" / "case_registry.json"
    try:
        registry = json.loads(read(registry_path))
    except json.JSONDecodeError as error:
        raise ContractFailure(f"invalid validation registry: {error}") from error
    cases = {str(item.get("id", "")).upper(): item
             for item in registry.get("cases", []) if isinstance(item, dict)}
    required_ids = {*(f"OV{i:02d}" for i in range(1, 6)), "EV01", "EV02"}
    require(required_ids <= cases.keys(),
            "validation registry is missing: " +
            ", ".join(sorted(required_ids - cases.keys())))

    for case_id in sorted(required_ids):
        descriptor = cases[case_id]
        require(descriptor.get("evidence_class") in {
                    "observational-validation", "campaign-evidence"},
                f"{case_id} has no machine-readable evidence_class")
        require(descriptor.get("release_status") in {
                    "release-gate", "diagnostic-only", "pilot-incomplete"},
                f"{case_id} has no machine-readable release_status")
        require(bool(descriptor.get("normalization_policy")),
                f"{case_id} has no declared normalization_policy")

    runner = read(SRCSEP / "test" / "run_tests.py")
    direct_runner = read(SRCSEP / "validation" / "run_case.py")
    for case_id in (f"OV{i:02d}" for i in range(1, 6)):
        require(case_id in runner and case_id in direct_runner,
                f"{case_id} is not protected by both fixed-input runner surfaces")

    makefile = read(SRCSEP / "makefile")
    for target in ("test-ov01-ov05-unit", "test-ev01-ev02-unit"):
        require(re.search(rf"(?m)^{re.escape(target)}\s*:", makefile) is not None,
                f"makefile does not expose {target}")


def check_species_contract() -> None:
    """Require one resolved species index and fail-closed AMPS binding."""
    header = read(SRCSEP3D / "runtime" / "run_configuration.h")
    implementation = read(SRCSEP3D / "runtime" / "run_configuration.cpp")
    parser = read(SRCSEP3D / "runtime" / "configuration_io.cpp")
    production = read(SRCSEP3D / "main_lib.cpp")

    require("ampsSpeciesIndex" in header,
            "SpeciesOptions has no immutable AMPS species index")
    require("ValidateSingleSpeciesBinding" in header and
            "ValidateSingleSpeciesBinding" in implementation,
            "AMPS-independent single-species validator is missing")
    require('field == "species.amps_index"' in parser,
            "configuration parser does not accept species.amps_index")
    require("ValidateSingleSpeciesBinding" in production,
            "amps_init does not validate the configured AMPS species table")
    require("source.species = 0" not in production,
            "shock injection still hard-codes AMPS species zero")
    require("source.speciesMassKg = PIC::MolecularData::GetMass(0)" not in production,
            "shock injection still hard-codes species-zero mass")


def main() -> int:
    try:
        check_release_hygiene()
        check_validation_registration()
        check_species_contract()
    except ContractFailure as error:
        print(f"FAIL STAGE3-CONTRACTS: {error}", file=sys.stderr)
        return 1
    print("PASS STAGE3-CONTRACTS: release hygiene, OV/EV registration, and "
          "single-species ownership are explicit")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
