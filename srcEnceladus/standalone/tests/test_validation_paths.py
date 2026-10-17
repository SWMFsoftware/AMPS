import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from validation_paths import (CO2_FIT_RELATIVE, JET_CATALOG_RELATIVE,  # noqa: E402
                              PLUME_FIT_RELATIVE, optional_root_path,
                              resolve_input_path, resolve_validation_root)


class ValidationPathTests(unittest.TestCase):
    def test_cli_root_overrides_environment_and_expands_standard_paths(self):
        with tempfile.TemporaryDirectory() as first, tempfile.TemporaryDirectory() as second:
            with patch.dict(os.environ, {"ENCELADUS_VALIDATION_DATA_ROOT": second}):
                root = resolve_validation_root(Path(first))
            self.assertEqual(root, Path(first).resolve())
            self.assertEqual(optional_root_path(None, root, JET_CATALOG_RELATIVE),
                             root/JET_CATALOG_RELATIVE)

    def test_environment_root_is_supported(self):
        with tempfile.TemporaryDirectory() as directory:
            with patch.dict(os.environ,
                            {"ENCELADUS_VALIDATION_DATA_ROOT": directory}):
                self.assertEqual(resolve_validation_root(None),
                                 Path(directory).resolve())

    def test_explicit_path_overrides_root_default(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            explicit = root/"custom"
            explicit.mkdir()
            self.assertEqual(
                resolve_input_path(explicit, root, PLUME_FIT_RELATIVE,
                                   "--data-root"), explicit.resolve())

    def test_missing_required_path_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "Required input does not exist"):
                resolve_input_path(None, Path(directory), CO2_FIT_RELATIVE,
                                   "--co2-data-root")


if __name__ == "__main__":
    unittest.main()
