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
* srcSEP3D enumerates and initializes the complete immutable AMPS species table
  rather than assuming a proton macro, slot zero, or runtime molecular data.

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
    """Require complete generated-table enumeration and fail-closed binding."""
    # The AMPS application deck is processed before compilation.  The
    # post-compile SEP input must not attempt to replace its species list.
    application_deck = read(AMPS_ROOT / "input" / "sep3d.input")
    species_assignments = re.findall(
        r"(?im)^\s*SpeciesList\s*=\s*([^!\r\n]+?)\s*(?:!.*)?$",
        application_deck)
    require(len(species_assignments) == 1 and species_assignments[0].split(),
            "input/sep3d.input must contain one non-empty SpeciesList")

    header = read(SRCSEP3D / "runtime" / "run_configuration.h")
    implementation = read(SRCSEP3D / "runtime" / "run_configuration.cpp")
    parser = read(SRCSEP3D / "runtime" / "configuration_io.cpp")
    production = read(SRCSEP3D / "main_lib.cpp")
    sep_initialization = read(SRCSEP / "main_lib.cpp")
    sep_prepopulation = read(SRCSEP / "shock_injection.cpp")
    sep_spherical_source = read(SRCSEP / "shock_analytical_model2D.cpp")

    require("CompiledSpeciesRecord" in header and
            "ValidateCompiledSpeciesBinding" in header and
            "ValidateCompiledSpeciesBinding" in implementation,
            "complete AMPS species-table binding contract is missing")
    for retired_key in ("species.amps_index", "species.name",
                        "species.mass_kg", "species.charge_c"):
        require(retired_key not in parser,
                f"runtime parser still accepts compiled species field {retired_key}")
    require("for (int speciesIndex = 0; speciesIndex < PIC::nTotalSpecies;" in
            production and
            "PIC::MolecularData::GetChemSymbol(speciesIndex)" in production,
            "production does not enumerate the complete generated ChemTable")
    require("PIC::MolecularData::SetMass(" not in production and
            "PIC::MolecularData::SetElectricCharge(" not in production,
            "production still overwrites immutable compiled molecular data")
    require("for (const auto& species : gCompiledSpecies)" in production,
            "particle numerics are not initialized for the full compiled table")
    require("for (const auto& compiled : gCompiledSpecies)" in production and
            "source.species = compiled.ampsIndex;" in production,
            "shock injection does not iterate every compiled species")
    require("ConfigureSpeciesSpectrum(" in production,
            "injection does not rebuild momentum bounds from each AMPS mass")

    # AMPS builds the final owner lists before allocating tree blocks, but its
    # cached BlockTable is refreshed separately.  All srcSEP3D initialization
    # passes use that cache, so it must be updated after allocation and before
    # the first background snapshot is collected.
    allocate_blocks = production.find("PIC::Mesh::mesh->AllocateTreeBlocks();")
    update_blocks = production.find(
        "PIC::DomainBlockDecomposition::UpdateBlockTable();", allocate_blocks)
    publish_background = production.find("FillAndPublishBackground();",
                                         update_blocks)
    require(allocate_blocks >= 0 and
            allocate_blocks < update_blocks < publish_background,
            "srcSEP3D does not refresh AMPS BlockTable before background fill")

    # Binding must occur after AMPS creates its base registries but before the
    # cell layout/mesh and Init_AfterParser.  Checking this order protects the
    # exact initialization boundary that first exposed the species mismatch.
    mesh_start = production.find("void amps_init_mesh()")
    pic_before = production.find("PIC::Init_BeforeParser();", mesh_start)
    bind_call = production.find("BindCompiledSpeciesTable();", pic_before)
    layout_freeze = production.find("PIC::Mesh::initCellSamplingDataBuffer();",
                                    bind_call)
    pic_after = production.find("PIC::Init_AfterParser();", mesh_start)
    require(mesh_start >= 0 and
            mesh_start < pic_before < bind_call < layout_freeze < pic_after,
            "AMPS species-table binding is not ordered before mesh/after-parser initialization")
    require("source.species = 0" not in production,
            "shock injection still hard-codes AMPS species zero")
    require("source.speciesMassKg = PIC::MolecularData::GetMass(0)" not in production,
            "shock injection still hard-codes species-zero mass")

    # srcSEP shares the same immutable AMPS SpeciesList contract.  Protect the
    # base numerical initialization and both legacy injection/prepopulation
    # paths against a regression to H_PLUS or slot-zero assumptions.
    require(re.search(
                r"for\s*\(int\s+species\s*=\s*0;\s*"
                r"species\s*<\s*PIC::nTotalSpecies;\s*\+\+species\)",
                sep_initialization) is not None,
            "srcSEP does not initialize every compiled AMPS species")
    require("for (int species=0;species<PIC::nTotalSpecies;++species)" in
            sep_prepopulation and "PopulateSegment(\n          species," in
            sep_prepopulation,
            "srcSEP field-line prepopulation does not enumerate every species")
    require("for (int spec=0;spec<PIC::nTotalSpecies;++spec)" in
            sep_spherical_source and
            "PIC::ParticleBuffer::SetI(spec,newParticleData)" in
            sep_spherical_source,
            "srcSEP spherical-shock injection does not enumerate every species")
    for source_name, source_text in (
            ("srcSEP initialization", sep_initialization),
            ("srcSEP prepopulation", sep_prepopulation),
            ("srcSEP spherical source", sep_spherical_source)):
        require("_H_PLUS_SPEC_" not in source_text,
                f"{source_name} still depends on _H_PLUS_SPEC_")


def main() -> int:
    try:
        check_release_hygiene()
        check_validation_registration()
        check_species_contract()
    except ContractFailure as error:
        print(f"FAIL STAGE3-CONTRACTS: {error}", file=sys.stderr)
        return 1
    print("PASS STAGE3-CONTRACTS: release hygiene, OV/EV registration, and "
          "complete compiled-species ownership are explicit")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
