#!/usr/bin/env python3
"""ARCHSCCM01 filesystem regressions followed by the real archive/ABI audit.

The temporary trees are adversarial fixtures, not substitute archives or MPI
qualification. Both this canonical entry point and make check-architecture
finish by running the actual built library's dependency/public-header checks.
"""
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import check_architecture as architecture


class ArchitectureSourceTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="sccm-source-audit-")
        self.root = Path(self.temporary.name) / "sep_coronal_cme"
        self.common = self.root.parent / "sep_common"
        for folder in (self.root/"include", self.root/"src", self.common):
            folder.mkdir(parents=True)

    def tearDown(self):
        self.temporary.cleanup()

    def test_binary_build_artifacts_and_source_named_directories_are_excluded(self):
        # Reproduce the reported invalid byte at offset 71. The old sep_*
        # glob attempted UTF-8 decoding of these adjacent compiled artifacts.
        binary = b"\x7fELF"+b"\x00"*67+b"\xc1\xff"
        for suffix in (".o", ".a", ".so", ".d"):
            (self.common/("sep_compiled"+suffix)).write_bytes(binary)
            (self.root/"src"/("compiled"+suffix)).write_bytes(binary)
        (self.common/"sep_directory.hpp").mkdir()
        (self.common/"sep_clean.h").write_text("#pragma once\n", encoding="utf-8")
        self.assertIsNone(architecture.source_dependency_error(self.root))

    def test_nested_and_alternate_source_suffixes_still_reject_forbidden_includes(self):
        nested = self.root/"include/nested"; nested.mkdir()
        for suffix, include in ((".hpp", "mpi.h"), (".CXX", "pic.h"), (".tpp", "swcme.h")):
            with self.subTest(suffix=suffix):
                source = nested/("dependency"+suffix)
                source.write_text('// fixture\n#include "'+include+'"\n', encoding="utf-8")
                error = architecture.source_dependency_error(self.root)
                self.assertIn("forbidden include", error)
                self.assertIn(str(source.relative_to(self.root))+":2", error)
                source.unlink()

    def test_neutral_sources_still_reject_upward_dependencies(self):
        nested = self.common/"nested"; nested.mkdir()
        for suffix in (".h", ".cpp", ".ipp"):
            with self.subTest(suffix=suffix):
                source = nested/("sep_dependency"+suffix)
                source.write_text('#include "sep_coronal_cme/constants.h"\n', encoding="utf-8")
                error = architecture.source_dependency_error(self.root)
                self.assertIn("neutral dependency points upward", error)
                self.assertIn("sep_common/nested/"+source.name, error)
                source.unlink()

    def test_invalid_utf8_in_real_sources_is_a_named_failure(self):
        for source in (self.common/"sep_invalid.h", self.root/"src/invalid.cpp"):
            with self.subTest(source=source.name):
                source.write_bytes(b"// invalid source \xc1\xff\n")
                error = architecture.source_dependency_error(self.root)
                self.assertIn("cannot read UTF-8 source", error)
                self.assertIn(source.name, error)
                source.unlink()

    def test_valid_unicode_source_comments_are_preserved(self):
        (self.common/"sep_unicode.h").write_text("// magnetic field δB, μ\n", encoding="utf-8")
        (self.root/"src/unicode.cpp").write_text("// plasma density ρ\n", encoding="utf-8")
        self.assertIsNone(architecture.source_dependency_error(self.root))


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(ArchitectureSourceTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():
        raise SystemExit(architecture.fail("source-audit regression failure"))
    raise SystemExit(architecture.main())
