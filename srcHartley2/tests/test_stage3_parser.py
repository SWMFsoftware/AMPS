#!/usr/bin/env python3
"""Portable source-only tests for Stage-3 comet input records.

The historical shell test copied a configured ``build/`` tree and consumed
root-level ``*.Assembled.*`` files. That made the test depend on whichever
application had last been configured on the developer's machine. These tests
instead reuse the source-only parser fixture shared with Stage 9: each case
copies tracked source defaults into a temporary directory, constructs the
assembled application block from the maintained Halley input, invokes the
tracked parser, and discards the generated tree afterwards.
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


def remove_record(lines: list[str], key: str) -> None:
    """Remove exactly one simple ``KEY=VALUE`` input record."""

    pattern = re.compile(rf"^\s*{re.escape(key)}\s*=", re.IGNORECASE)
    matches = [index for index, line in enumerate(lines) if pattern.match(line)]
    if len(matches) != 1:
        raise AssertionError(f"expected one {key} record, found {len(matches)}")
    del lines[matches[0]]


class Stage3ParserTest(unittest.TestCase):
    def run_case(
        self, lines: list[str], expected_error: str | None = None
    ) -> ParserFixture:
        temporary = tempfile.TemporaryDirectory(prefix="hartley2-stage3-parser.")
        self.addCleanup(temporary.cleanup)
        fixture = ParserFixture(Path(temporary.name))
        fixture.write_block(lines)
        result = fixture.run()
        if expected_error is None:
            self.assertEqual(result.returncode, 0, result.stdout)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertIn(expected_error, result.stdout)
        return fixture

    def test_valid_systematic_configuration(self) -> None:
        fixture = self.run_case(input_block())
        model = (fixture.work / "main" / "CometModel.cpp").read_text(
            encoding="utf-8"
        )
        chemistry = (fixture.work / "main" / "CometChemistry.cpp").read_text(
            encoding="utf-8"
        )
        self.assertIn("std::uint64_t Comet::Model::RandomSeed=104729;", model)
        self.assertIn("bool Comet::Model::RandomSeedConfigured=true;", model)
        self.assertIn("ResidualSampling::RandomPhaseSystematic", chemistry)

    def test_omitted_records_retain_legacy_defaults(self) -> None:
        lines = input_block()
        remove_record(lines, "CometRandomSeed")
        remove_record(lines, "ChemistryResidualSampling")
        fixture = self.run_case(lines)
        model = (fixture.work / "main" / "CometModel.cpp").read_text(
            encoding="utf-8"
        )
        chemistry = (fixture.work / "main" / "CometChemistry.cpp").read_text(
            encoding="utf-8"
        )
        self.assertIn("bool Comet::Model::RandomSeedConfigured=false;", model)
        self.assertIn("ResidualSampling::IndependentBernoulli;", chemistry)

    def test_invalid_and_duplicate_seed_are_rejected(self) -> None:
        invalid = input_block()
        replace_record(invalid, "CometRandomSeed", "0")
        self.run_case(invalid, "CometRandomSeed must be an integer")

        duplicate = input_block()
        position = next(
            index
            for index, line in enumerate(duplicate)
            if re.match(r"^\s*CometRandomSeed\s*=", line, re.IGNORECASE)
        )
        duplicate.insert(position + 1, "  CometRandomSeed=5")
        self.run_case(duplicate, "duplicate CometRandomSeed")

    def test_invalid_and_duplicate_sampling_mode_are_rejected(self) -> None:
        invalid = input_block()
        replace_record(invalid, "ChemistryResidualSampling", "invalid")
        self.run_case(invalid, "must be BERNOULLI or SYSTEMATIC")

        duplicate = input_block()
        position = next(
            index
            for index, line in enumerate(duplicate)
            if re.match(
                r"^\s*ChemistryResidualSampling\s*=", line, re.IGNORECASE
            )
        )
        duplicate.insert(position + 1, "  ChemistryResidualSampling=bernoulli")
        self.run_case(duplicate, "duplicate ChemistryResidualSampling")

    def test_systematic_requires_seed(self) -> None:
        lines = input_block()
        remove_record(lines, "CometRandomSeed")
        self.run_case(lines, "SYSTEMATIC requires an explicit CometRandomSeed")

    def test_sampling_mode_requires_chemistry(self) -> None:
        lines = input_block()
        replace_record(lines, "ChemistryMode", "off")
        self.run_case(lines, "ChemistryResidualSampling requires ChemistryMode=ON")


if __name__ == "__main__":
    unittest.main(verbosity=2)
