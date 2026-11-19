#!/usr/bin/env python3
"""Manifest path/error controls ONLY: no native or observational evidence.

All temporary bytes below are deliberately labelled protocol fixtures. Native
execution is prohibited by a mock. The expected result is retained ERROR;
none of these metadata fixtures is submitted as a production PASS receipt.
"""
import hashlib
import json
from pathlib import Path
import tempfile
import types
import unittest
from unittest import mock

import run_m3_transport_clock as clock


class ClockManifestPaths(unittest.TestCase):
    def test_nonobject_build_receipt_retains_error_without_execution(self):
        with tempfile.TemporaryDirectory(prefix='amps-clock-build-shape-protocol-') as directory:
            out=Path(directory)/'expected-error'
            with mock.patch.object(clock,'_execute',side_effect=AssertionError('protocol must not execute')):
                result=clock.run_startup(out,30,build_evidence=[])
            self.assertEqual(result['status'],'ERROR')
            self.assertEqual(result['exception_type'],'ValueError')
            self.assertTrue((out/'infrastructure_exception.log').is_file())
            self.assertTrue((out/'result.json').is_file())
            self.assertEqual(result['commands'],[])

    def test_unexpected_local_attribute_error_is_retained_with_traceback(self):
        with tempfile.TemporaryDirectory(prefix='amps-clock-unexpected-protocol-') as directory:
            out=Path(directory)/'expected-error'
            with mock.patch.object(clock.shutil,'which',side_effect=AttributeError('PROTOCOL-ONLY unexpected dependency object')):
                result=clock.run_local(out,30)
            self.assertEqual(result['status'],'ERROR')
            self.assertEqual(result['exception_type'],'AttributeError')
            self.assertIn('AttributeError',(out/'infrastructure_exception.log').read_text())
            self.assertTrue((out/'run_manifest.json').is_file())

    def test_nonobject_producer_json_is_retained_as_error(self):
        with tempfile.TemporaryDirectory(prefix='amps-clock-producer-shape-protocol-') as directory:
            root=Path(directory);producer=root/'PROTOCOL-ONLY';producer.mkdir()
            manifest=producer/'run_manifest.json';manifest.write_text('[]\n')
            (producer/'result.json').write_text(json.dumps({'status':'PASS','scope':'PROTOCOL METADATA STUB ONLY; NEVER NATIVE EVIDENCE'}))
            out=root/'expected-error'
            with mock.patch.object(clock.shutil,'which',return_value='/not-executed/protocol-tool'), \
                 mock.patch.dict('sys.modules',{'matplotlib':types.SimpleNamespace()}), \
                 mock.patch.object(clock,'_execute',side_effect=AssertionError('protocol must never execute')):
                result=clock.run(out,manifest,30)
            self.assertEqual(result['status'],'ERROR')
            self.assertTrue((out/'infrastructure_exception.log').is_file())
            self.assertEqual(result['commands'],[])

    def test_json_string_paths_hash_identically_to_path_objects(self):
        with tempfile.TemporaryDirectory(prefix='amps-clock-path-protocol-') as directory:
            path=Path(directory)/'only-protocol-bytes.txt';path.write_bytes(b'NOT SCIENCE OR A NATIVE LIBRARY\n')
            expected=hashlib.sha256(path.read_bytes()).hexdigest()
            self.assertEqual(clock.sha256(str(path)),expected)
            self.assertEqual(clock.sha256(path),expected)

    def test_missing_string_generated_path_retains_error_before_any_launch(self):
        with tempfile.TemporaryDirectory(prefix='amps-clock-error-protocol-') as directory:
            root=Path(directory);producer=root/'PROTOCOL-ONLY-producer';producer.mkdir()
            library=root/'NOT-A-NATIVE-LIBRARY.a';library.write_text('PROTOCOL FIXTURE - NEVER EXECUTED\n')
            missing=producer/'generated/pic/missing-generated-header.dfn'
            document={'scope':'PROTOCOL FIXTURE ONLY, NOT NATIVE EVIDENCE','source_root':str(clock.REPO),
                'library':{'path':str(library),'sha256':clock.sha256(library)},
                'generated_snapshots':{str(missing):'0'*64},'production_source_hashes':{}}
            manifest=producer/'run_manifest.json';manifest.write_text(json.dumps(document))
            # Simulate prerequisite metadata solely to reach the JSON path
            # decoder. Expected ERROR is the only output this test accepts.
            (producer/'result.json').write_text(json.dumps({'status':'PASS','scope':'PROTOCOL-ONLY PREREQUISITE STUB; NEVER NATIVE EVIDENCE'}))
            out=root/'expected-error'
            with mock.patch.object(clock.shutil,'which',return_value='/not-executed/protocol-tool'), \
                 mock.patch.dict('sys.modules',{'matplotlib':types.SimpleNamespace()}), \
                 mock.patch.object(clock,'_execute',side_effect=AssertionError('protocol must never launch native work')):
                result=clock.run(out,manifest,30)
            self.assertEqual(result['status'],'ERROR')
            self.assertIn('missing-generated-header.dfn',result['reason'])
            self.assertTrue((out/'result.json').is_file())
            retained=json.loads((out/'run_manifest.json').read_text())
            self.assertEqual(retained['execution'],'not-run')
            self.assertEqual(result['commands'],[])


if __name__=='__main__':
    unittest.main()
