#!/usr/bin/env python3
"""Integration checks for hartley2.pl brightness source generation.

The fixture is shared with the maintained Stage-9 parser tests so it copies a
fresh source-only AMPS work tree.  No configured repository build, I3 input,
or science output is touched.
"""

from __future__ import annotations

import re
import tempfile
import unittest
from pathlib import Path

from test_stage9_science_diagnostics_parser import (
    ParserFixture,
    input_block,
    replace_record,
)


def with_brightness() -> list[str]:
    lines = input_block()
    replace_record(lines, "BrightnessMode", "on")
    insertion = next(
        index
        for index, line in enumerate(lines)
        if re.match(r"^\s*CometDiagnostics\s*=", line, re.I)
    ) + 1
    lines[insertion:insertion] = [
        "BrightnessImpactParameterGrid=log | 6000 | 1.0e9 | 17",
        "BrightnessObservationGeometry=1|0|0|0|0|1|1|0|0",
        (
            "BrightnessEmission=H | H_LYMAN_ALPHA | 2.5e-3 | "
            "TEST_H_GFACTOR | thin | isotropic_4pi | line_center | "
            "1.0e-17 | 0.1 | TEST_H_SIGMA"
        ),
        (
            "BrightnessEmission=OH | OH_A2SIGMA_X2PI_00 | 2.0e-4 | "
            "TEST_OH_GFACTOR | thin | isotropic_4pi | none | 0 | 0 | NONE"
        ),
    ]
    return lines


def with_core_brightness_models() -> list[str]:
    lines = input_block()
    replace_record(lines, "BrightnessMode", "on")
    insertion = next(
        index
        for index, line in enumerate(lines)
        if re.match(r"^\s*CometDiagnostics\s*=", line, re.I)
    ) + 1
    lines[insertion:insertion] = [
        "BrightnessImpactParameterGrid=log | 6000 | 1.0e9 | 17",
        "BrightnessObservationGeometry=1|0|0|0|0|1|1|0|0",
        "BrightnessHeliocentricRadialVelocityKMS=-39.7",
        "BrightnessVelocityDomainPolicy=ABORT",
        "BrightnessSolarLymanAlphaFlux1AUPhotonCM2S=5.1040518742e11",
        (
            "BrightnessEmission=H | H_LYMAN_ALPHA | H_LYA_IKL2020 | "
            "CORE_IKL2020 | thin | isotropic_4pi | line_center | "
            "1.0e-17 | 0.1 | TEST_H_SIGMA | 0"
        ),
        (
            "BrightnessEmission=OH | OH_A2SIGMA_X2PI_00 | "
            "OH_AX00_SA1988_LOWELL_UNQUENCHED | CORE_LOWELL | thin | "
            "isotropic_4pi | none | 0 | 0 | NONE | 308.5"
        ),
    ]
    return lines


class BrightnessParserIntegrationTest(unittest.TestCase):
    def fixture(self, root: Path, lines: list[str]) -> ParserFixture:
        fixture = ParserFixture(root)
        fixture.write_block(lines)
        result = fixture.run()
        self.assertEqual(result.returncode, 0, result.stdout)
        return fixture

    def test_halley_explicit_off_keeps_brightness_disabled(self) -> None:
        with tempfile.TemporaryDirectory(prefix="hartley2-brightness-off.") as root:
            fixture = self.fixture(Path(root), input_block())
            source = (fixture.work / "main" / "CometBrightness.cpp").read_text(
                encoding="utf-8"
            )
            self.assertIn("bool Comet::Brightness::Enabled=false;", source)
            self.assertNotIn("SetGeneratedMetadata(0,", source)

    def test_enabled_case_generates_all_explicit_metadata(self) -> None:
        with tempfile.TemporaryDirectory(prefix="hartley2-brightness-on.") as root:
            fixture = self.fixture(Path(root), with_brightness())
            source = (fixture.work / "main" / "CometBrightness.cpp").read_text(
                encoding="utf-8"
            )
            makefile = (fixture.work / "main" / "makefile").read_text(
                encoding="utf-8"
            )
            self.assertIn("bool Comet::Brightness::Enabled=true;", source)
            self.assertIn("int Comet::Brightness::EmissionLineCount=2;", source)
            # Check the generated numerical arrays, not merely their
            # declarations. A missing regex escape in hartley2.pl once left
            # all five declarations at zero while metadata strings and the
            # enabled flag looked correct, which would make runtime
            # validation fail only after a complete AMPS build.
            self.assertIn(
                "int Comet::Brightness::EmissionSpecies"
                "[Comet::Brightness::MaxEmissionLines]={0,12};",
                source,
            )
            self.assertIn(
                "double Comet::Brightness::PhotonGFactor"
                "[Comet::Brightness::MaxEmissionLines]={0.0025,0.0002};",
                source,
            )
            self.assertIn(
                "int Comet::Brightness::TauDiagnosticMode"
                "[Comet::Brightness::MaxEmissionLines]={1,0};",
                source,
            )
            self.assertIn(
                "double Comet::Brightness::LineCenterCrossSection"
                "[Comet::Brightness::MaxEmissionLines]={1e-17,0};",
                source,
            )
            self.assertIn(
                "double Comet::Brightness::TauFlagThreshold"
                "[Comet::Brightness::MaxEmissionLines]={0.1,0};",
                source,
            )
            self.assertIn('SetGeneratedMetadata(0,"H_LYMAN_ALPHA"', source)
            self.assertIn('SetGeneratedMetadata(1,"OH_A2SIGMA_X2PI_00"', source)
            self.assertIn("TEST_H_GFACTOR", source)
            self.assertIn("TEST_OH_GFACTOR", source)
            self.assertIn("TEST_H_SIGMA", source)
            self.assertIn("PARTICLE_RESOLVED_DIM1_SHELL_GAUSS2", source)
            self.assertIn("pic.\"+symbol+\".brightness.s=", source)
            self.assertIn(".out=\"+std::to_string(outputIndex)+\".dat", source)
            self.assertIn('TITLE=\\"Particle-resolved comet brightness:', source)
            self.assertIn('VARIABLES=\\"Impact Parameter [m]\\"', source)
            self.assertIn("DATAPACKING=POINT", source)
            self.assertIn("DATASETAUXDATA CometHeliocentricVelocityMPS", source)
            self.assertIn("Line Center Tau", source)
            self.assertIn("tau_provenance=", source)
            # Center-node storage includes ghost cells, whereas the particle
            # list table contains only physical cells.  A former regression
            # reused the ghost-aware center-node index here and followed all
            # 20 million free-buffer entries on the first empty sampling
            # iteration.  Assert the two index spaces remain separate in the
            # configured source copied by the parser fixture.
            self.assertIn(
                "const long int particleCellNumber=localCell;", source
            )
            self.assertIn(
                "FirstCellParticleTable[particleCellNumber]", source
            )
            self.assertNotIn("FirstCellParticleTable[cellNumber]", source)
            self.assertRegex(makefile, r"(?m)^MAINLIBOBJ = .*\bCometBrightness\.o\b")

    def test_named_models_generate_conditions_and_model_arrays(self) -> None:
        with tempfile.TemporaryDirectory(prefix="hartley2-brightness-models.") as root:
            fixture = self.fixture(Path(root), with_core_brightness_models())
            source = (fixture.work / "main" / "CometBrightness.cpp").read_text(
                encoding="utf-8"
            )
            self.assertIn(
                "int Comet::Brightness::GFactorModel"
                "[Comet::Brightness::MaxEmissionLines]={1,2};",
                source,
            )
            self.assertIn(
                "double Comet::Brightness::PhotonGFactor"
                "[Comet::Brightness::MaxEmissionLines]={0,0};",
                source,
            )
            self.assertIn(
                "double Comet::Brightness::EffectivePhotonWavelength"
                "[Comet::Brightness::MaxEmissionLines]={0,3.085e-07};",
                source,
            )
            self.assertIn(
                "double Comet::Brightness::CometHeliocentricVelocity"
                "[3]={-39700,0,0};",
                source,
            )
            self.assertIn(
                "double Comet::Brightness::HeliocentricRadialDirection"
                "[3]={1,0,0};",
                source,
            )
            self.assertIn(
                "double Comet::Brightness::LineOfSightDirection"
                "[3]={0,0,1};",
                source,
            )
            self.assertIn(
                "bool Comet::Brightness::ObservationGeometryConfigured=true;",
                source,
            )
            self.assertIn(
                "int Comet::Brightness::OutOfRangeVelocityPolicy="
                "Comet::Brightness::VelocityDomainPolicy::Abort;",
                source,
            )
            self.assertEqual(
                source.count("Comet::Brightness::VelocityDomainPolicy::Abort;"),
                1,
                "the generated policy declaration must not leave an orphaned "
                "initializer continuation",
            )
            self.assertIn(
                "double Comet::Brightness::SolarLymanAlphaFlux1AU=5104051874200000;",
                source,
            )


if __name__ == "__main__":
    unittest.main()
