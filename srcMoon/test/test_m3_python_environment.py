#!/usr/bin/env python3
"""Controls for the Python interpreter check of the I23/I24 supporting runners.

A missing third-party package is a missing dependency: each runner must report
SKIPPED before any compilation or MPI launch, name the interpreter, and record
the interpreter in its manifest. matplotlib is hidden in-process with
``sys.modules`` (restored afterwards). The check runs after input
authentication, so a corrupt input is still ERROR. No native code runs and
nothing here is physical evidence.
"""
from argparse import Namespace
import importlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

TEST = Path(__file__).resolve().parent
for folder in (TEST, TEST / 'integration-tests/I24', TEST / 'integration-tests/I23'):
    if str(folder) not in sys.path:
        sys.path.insert(0, str(folder))
from common import python_environment as environment  # noqa: E402

HIDDEN = {'matplotlib': None, 'matplotlib.pyplot': None}
# Import every available package once before any sys.modules patch, so the
# patch restore cannot drop a module first imported inside it (NumPy refuses
# to be imported twice).
environment.describe()


class Helper(unittest.TestCase):
    def test_absent_module_is_reported_with_the_interpreter(self):
        class Missing(Exception):
            pass
        name = 'amps_module_that_does_not_exist'
        self.assertEqual(environment.absent([name, 'json']), [name])
        with self.assertRaises(Missing) as caught:
            environment.require([name], Missing)
        self.assertIn(name, str(caught.exception))
        self.assertIn(sys.executable, str(caught.exception))
        self.assertIn(environment.MAINTAINED_INTERPRETER, str(caught.exception))

    def test_describe_records_versions_and_missing_packages(self):
        with patch.dict(sys.modules, HIDDEN):
            value = environment.describe(('matplotlib', 'json'))
        self.assertEqual(value['executable'], sys.executable)
        self.assertIsNone(value['packages']['matplotlib'])
        self.assertIsNotNone(value['packages']['json'])


class RunnersSkipBeforeNativeWork(unittest.TestCase):
    """The check follows input authentication and precedes all native work.

    Input preparation and the producer receipt check are replaced by stubs, so
    no real source, producer or build is needed. A stub that raises shows that
    a corrupt input remains ERROR when a package is also missing.
    """
    def run_case(self, module, call, corrupt=False):
        runner = importlib.import_module(module)
        def prepare(*args, **kwargs):
            if corrupt:
                raise ValueError('changed source identity: test stub')
            return [], {}, {}
        with tempfile.TemporaryDirectory() as work, patch.dict(sys.modules, HIDDEN), \
                patch.object(runner, 'prepare_inputs', side_effect=prepare), patch('sys.stdout'):
            producer = Path(work) / 'producer.json'
            producer.write_text('{"stub": "receipt check replaced"}')
            output = Path(work) / 'out'
            if hasattr(runner, 'checked_build'):
                with patch.object(runner, 'checked_build', return_value=Path(work)):
                    call(runner, output, producer)
            else:
                call(runner, output, producer)
            result = json.loads((output / 'result.json').read_text())
            manifest = json.loads((output / 'run_manifest.json').read_text())
        self.assertEqual(manifest['commands'], [])
        self.assertIsNone(manifest['python']['packages']['matplotlib'])
        self.assertEqual(manifest['python']['executable'], sys.executable)
        if corrupt:
            self.assertEqual(result['status'], 'ERROR', result)
            self.assertIn('changed source identity', result['reason'])
        else:
            self.assertEqual(result['status'], 'SKIPPED', result)
            self.assertIn('matplotlib', result['reason'])
            self.assertIn(sys.executable, result['reason'])

    CASES = {
        'run_ground_capture': lambda r, out, path: r.run(path, path, out, 8, 60),
        'run_pic_ground_capture': lambda r, out, path: r.main(Namespace(
            producer_manifest=path, lock=path, mass_inventory=path, output_dir=out,
            refinement=False, wall_seconds=60)),
        'run_total_capture': lambda r, out, path: r.run(out, 8, 60),
        'run_pic_total_capture': lambda r, out, path: r.main(Namespace(
            producer_manifest=path, output_dir=out, refinement=False, wall_seconds=60)),
        'run_helium_impact': lambda r, out, path: r.run(path, path, out, 8, 60),
        'run_pic_helium_impact': lambda r, out, path: r.main(Namespace(
            producer_manifest=path, profile=path, data_root=path, output_dir=out,
            refinement=False, wall_seconds=60)),
    }

    def test_missing_package_skips_after_authentication(self):
        for module, call in self.CASES.items():
            with self.subTest(runner=module):
                self.run_case(module, call)

    def test_corrupt_input_stays_error_without_the_package(self):
        for module, call in self.CASES.items():
            with self.subTest(runner=module):
                self.run_case(module, call, corrupt=True)


if __name__ == '__main__':
    unittest.main()
