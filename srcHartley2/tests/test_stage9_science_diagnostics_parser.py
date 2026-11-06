#!/usr/bin/env python3
"""Isolated parser tests for the Stage-9 science-diagnostic switch.

The application input names ``../hartley2.pl`` because configured run trees
carry a generated root copy.  These tests intentionally invoke the tracked
``input/hartley2.pl`` instead: a stale generated copy must not make a parser
test pass.  Each case receives a fresh source-only working tree, so no case
can inherit generated flags or policy tables from a previous configuration.
No AMPS compilation or MPI launch is performed.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


TEST_DIR = Path(__file__).resolve().parent
SOURCE_ROOT = TEST_DIR.parents[1]
PARSER = SOURCE_ROOT / "input" / "hartley2.pl"
INPUT = SOURCE_ROOT / "input" / "halley1d.input"

SPECIES = (
    "H", "H2", "C", "CH", "CH2", "N", "CH3", "NH", "CH4", "NH2",
    "O", "NH3", "OH", "H2O", "C2", "C2H", "C2H2", "CN", "HCN",
    "C2H4", "CO", "C2H5", "HCO", "C2H6", "H2CO", "CH3OH", "CO2",
)

ENHANCED_RECORDS = (
    "ParticleWeightCorrectionWindowMode=on",
    "ParticleWeightCorrectionWindow=1.0e-8 | 0.25 | 64.0",
    "ParticleWeightCorrectionWindowMaxSplits=512",
    "ChemistryDaughterPacketMode=on",
    "ChemistryDaughterPacketMaximum=8",
    "ChemistryDaughterPacketCount(GLOBAL)=1",
    "ChemistryDaughterPacketCount(REACTION:C2_2C)=8",
    "ChemistryChannelImportanceMode=on",
    "ChemistryChannelBiasedProbability(REACTION:C2_2C)=0.5",
    "ChemistryChannelBiasedProbability(REACTION:C2_ION)=0.5",
    "ChemistryTotalReactionForcingMode=on",
    "ChemistryTotalReactionForcing(SPECIES:C2)=on",
    "Stage9ScienceDiagnostics=on",
)


def input_block() -> list[str]:
    """Return only the application block from the maintained Halley input."""

    lines = INPUT.read_text(encoding="utf-8").splitlines()
    begin = next(index for index, line in enumerate(lines) if line.strip().lower().startswith("#block"))
    end = next(index for index in range(begin + 1, len(lines)) if lines[index].strip().lower() == "#endblock")
    # AssambleInputFile omits blank/comment-only records before invoking an
    # application parser.  Reproduce that contract here; feeding an empty
    # record directly to hartley2.pl would test an impossible block format.
    return [
        line
        for line in lines[begin + 1 : end]
        if line.strip() and not line.lstrip().startswith("!")
    ]


def replace_record(lines: list[str], key: str, value: str) -> None:
    pattern = re.compile(rf"^(\s*){re.escape(key)}\s*=.*$", re.IGNORECASE)
    matches = [index for index, line in enumerate(lines) if pattern.match(line)]
    if len(matches) != 1:
        raise AssertionError(f"expected one {key} record, found {len(matches)}")
    indentation = pattern.match(lines[matches[0]]).group(1)  # type: ignore[union-attr]
    lines[matches[0]] = f"{indentation}{key}={value}"


class ParserFixture:
    def __init__(self, root: Path):
        self.root = root
        self.work = root / "work"
        self.block = root / "case.Block"
        self.species = root / "case.Species"
        self.log = root / "parser.log"
        self._copy_sources()

    def _copy_sources(self) -> None:
        (self.work / "main").mkdir(parents=True)
        for pattern in ("*.cpp", "*.h", "*.dfn", "makefile"):
            for source in (SOURCE_ROOT / "srcHartley2").glob(pattern):
                shutil.copy2(source, self.work / "main" / source.name)

        required = {
            "pic": ("pic.h", "picGlobal.dfn"),
            "meshAMR": ("meshAMR_UserDefinitions.h", "meshAMRdef.h"),
            "models/dust": ("Dust.cpp", "Dust.dfn"),
            "models/exosphere": ("Exosphere_Chemistry.cpp",),
        }
        for relative, names in required.items():
            destination = self.work / relative
            destination.mkdir(parents=True)
            for name in names:
                shutil.copy2(SOURCE_ROOT / "src" / relative / name, destination / name)

        self.species.write_text(
            f"{len(SPECIES)}\n" + "\n".join(SPECIES) + "\n", encoding="utf-8"
        )

    def set_individual_corrections(self, enabled: bool) -> None:
        path = self.work / "pic" / "picGlobal.dfn"
        text = path.read_text(encoding="utf-8")
        replacement = (
            "_INDIVIDUAL_PARTICLE_WEIGHT_ON_"
            if enabled
            else "_INDIVIDUAL_PARTICLE_WEIGHT_OFF_"
        )
        text, count = re.subn(
            r"(?m)^(#define\s+_INDIVIDUAL_PARTICLE_WEIGHT_MODE_\s+)\S+.*$",
            rf"\g<1>{replacement}",
            text,
        )
        if count != 1:
            raise AssertionError("individual-correction macro was not unique")
        path.write_text(text + ("" if text.endswith("\n") else "\n"), encoding="utf-8")

    def write_block(self, lines: list[str]) -> None:
        # Application block files alternate provenance and source lines.
        assembled: list[str] = []
        for number, line in enumerate(lines, start=1):
            assembled.extend((f"{number} stage9-parser-test.input", line))
        self.block.write_text("\n".join(assembled) + "\n", encoding="utf-8")

    def run(self) -> subprocess.CompletedProcess[str]:
        result = subprocess.run(
            ["perl", str(PARSER), str(self.block), str(self.work), str(SOURCE_ROOT)],
            cwd=SOURCE_ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        self.log.write_text(result.stdout, encoding="utf-8")
        return result


def configured_lines(stage9: str | None = "on") -> list[str]:
    lines = input_block()
    replace_record(lines, "ChemistryResidualSampling", "bernoulli")
    insertion = next(index for index, line in enumerate(lines) if re.match(r"^\s*CometDiagnostics\s*=", line, re.I)) + 1
    records = list(ENHANCED_RECORDS)
    if stage9 is None:
        records = [record for record in records if not record.startswith("Stage9ScienceDiagnostics=")]
    elif stage9.lower() != "on":
        stage_index = next(
            index
            for index, record in enumerate(records)
            if record.startswith("Stage9ScienceDiagnostics=")
        )
        records[stage_index] = f"Stage9ScienceDiagnostics={stage9}"
    lines[insertion:insertion] = [f"  {record}" for record in records]
    return lines


def set_inserted_mode(lines: list[str], key: str, value: str) -> None:
    replace_record(lines, key, value)


def remove_records(lines: list[str], *keys: str) -> None:
    """Remove complete input records by their case-insensitive key."""

    patterns = [re.compile(rf"^\s*{re.escape(key)}(?:\([^)]*\))?\s*=", re.I) for key in keys]
    lines[:] = [line for line in lines if not any(pattern.match(line) for pattern in patterns)]


def disable_correction_windows(lines: list[str]) -> None:
    set_inserted_mode(lines, "ParticleWeightCorrectionWindowMode", "off")
    remove_records(
        lines,
        "ParticleWeightCorrectionWindow",
        "ParticleWeightCorrectionWindowMaxSplits",
    )


def disable_stage6(lines: list[str]) -> None:
    set_inserted_mode(lines, "ChemistryDaughterPacketMode", "off")
    remove_records(lines, "ChemistryDaughterPacketMaximum", "ChemistryDaughterPacketCount")


def disable_stage7(lines: list[str]) -> None:
    set_inserted_mode(lines, "ChemistryChannelImportanceMode", "off")
    remove_records(lines, "ChemistryChannelBiasedProbability")


def disable_stage8(lines: list[str]) -> None:
    set_inserted_mode(lines, "ChemistryTotalReactionForcingMode", "off")
    remove_records(lines, "ChemistryTotalReactionForcing")


class Stage9ScienceDiagnosticsParserTest(unittest.TestCase):
    def run_case(
        self,
        lines: list[str],
        *,
        corrections: bool = True,
        expected_error: str | None = None,
    ) -> tuple[ParserFixture, subprocess.CompletedProcess[str]]:
        temporary = tempfile.TemporaryDirectory(prefix="hartley2-stage9-parser.")
        self.addCleanup(temporary.cleanup)
        fixture = ParserFixture(Path(temporary.name))
        fixture.set_individual_corrections(corrections)
        fixture.write_block(lines)
        result = fixture.run()
        if expected_error is None:
            self.assertEqual(result.returncode, 0, result.stdout)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertIn(expected_error, result.stdout)
        return fixture, result

    def test_valid_on_sets_runtime_switch(self) -> None:
        fixture, _ = self.run_case(configured_lines())
        text = (fixture.work / "main" / "RadialSimulation.cpp").read_text(encoding="utf-8")
        self.assertIn("Stage9ScienceDiagnosticsEnabled=true;", text)

    def test_omitted_and_explicit_off_keep_legacy_default(self) -> None:
        for value in (None, "off"):
            with self.subTest(value=value):
                fixture, _ = self.run_case(configured_lines(value))
                text = (fixture.work / "main" / "RadialSimulation.cpp").read_text(encoding="utf-8")
                self.assertIn("Stage9ScienceDiagnosticsEnabled=false;", text)

    def test_invalid_token_is_rejected(self) -> None:
        self.run_case(
            configured_lines("sometimes"),
            expected_error="Stage9ScienceDiagnostics must be ON or OFF",
        )

    def test_dimension_prerequisite(self) -> None:
        lines = configured_lines()
        replace_record(lines, "Dimension", "3")
        replace_record(lines, "Geometry", "Triangulated")
        replace_record(lines, "ChemistryMode", "off")
        remove_records(lines, "ChemistryResidualSampling")
        replace_record(lines, "RadialSimulationMode", "off")
        disable_correction_windows(lines)
        disable_stage6(lines)
        disable_stage7(lines)
        disable_stage8(lines)
        self.run_case(lines, expected_error="Stage9ScienceDiagnostics=ON requires Dimension=1")

    def test_radial_prerequisite(self) -> None:
        lines = configured_lines()
        replace_record(lines, "RadialSimulationMode", "off")
        disable_correction_windows(lines)
        disable_stage6(lines)
        disable_stage7(lines)
        disable_stage8(lines)
        self.run_case(lines, expected_error="Stage9ScienceDiagnostics=ON requires RadialSimulationMode=ON")

    def test_chemistry_prerequisite(self) -> None:
        lines = configured_lines()
        replace_record(lines, "ChemistryMode", "off")
        remove_records(lines, "ChemistryResidualSampling")
        disable_stage6(lines)
        disable_stage7(lines)
        disable_stage8(lines)
        self.run_case(lines, expected_error="Stage9ScienceDiagnostics=ON requires ChemistryMode=ON")

    def test_comet_diagnostics_prerequisite(self) -> None:
        lines = configured_lines()
        replace_record(lines, "CometDiagnostics", "off")
        disable_stage7(lines)
        disable_stage8(lines)
        self.run_case(lines, expected_error="Stage9ScienceDiagnostics=ON requires CometDiagnostics=ON")

    def test_individual_correction_prerequisite(self) -> None:
        lines = configured_lines()
        disable_correction_windows(lines)
        disable_stage6(lines)
        disable_stage7(lines)
        disable_stage8(lines)
        self.run_case(
            lines,
            corrections=False,
            expected_error="Stage9ScienceDiagnostics=ON requires ParticleWeightCorrectionMode=ON",
        )

    def test_bernoulli_prerequisite(self) -> None:
        lines = configured_lines()
        replace_record(lines, "ChemistryResidualSampling", "systematic")
        disable_stage6(lines)
        disable_stage7(lines)
        disable_stage8(lines)
        self.run_case(
            lines,
            expected_error="Stage9ScienceDiagnostics=ON requires ChemistryResidualSampling=BERNOULLI",
        )

    def test_correction_window_prerequisite(self) -> None:
        lines = configured_lines()
        disable_correction_windows(lines)
        disable_stage6(lines)
        disable_stage7(lines)
        disable_stage8(lines)
        self.run_case(
            lines,
            expected_error="Stage9ScienceDiagnostics=ON requires ParticleWeightCorrectionWindowMode=ON",
        )

    def test_stage6_prerequisite(self) -> None:
        lines = configured_lines()
        disable_stage6(lines)
        disable_stage7(lines)
        disable_stage8(lines)
        self.run_case(lines, expected_error="Stage9ScienceDiagnostics=ON requires ChemistryDaughterPacketMode=ON")

    def test_stage7_prerequisite(self) -> None:
        lines = configured_lines()
        disable_stage7(lines)
        disable_stage8(lines)
        self.run_case(lines, expected_error="Stage9ScienceDiagnostics=ON requires ChemistryChannelImportanceMode=ON")

    def test_stage8_prerequisite(self) -> None:
        lines = configured_lines()
        disable_stage8(lines)
        self.run_case(lines, expected_error="Stage9ScienceDiagnostics=ON requires ChemistryTotalReactionForcingMode=ON")


if __name__ == "__main__":
    unittest.main(verbosity=2)
