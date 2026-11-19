"""Authentication/ERROR controls for the single-workbook source entry point."""
from pathlib import Path
import contextlib
import hashlib
import importlib.util
import io
import json
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[2]
PATH=ROOT/"utility/chemical_database/new_surface_workbooks.py"
SPEC=importlib.util.spec_from_file_location("surface_workbooks",PATH)
M=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(M)
DECODER=hashlib.sha256(PATH.with_name("native_workbooks.py").read_bytes()).hexdigest()


class Protocol(unittest.TestCase):
    def test_corrupt_authenticated_archive_reports_error_without_output(self):
        """A matching hash authenticates bytes; it does not make them valid XLSX."""
        with tempfile.TemporaryDirectory()as tmp:
            source,out=Path(tmp)/"bad.xlsx",Path(tmp)/"out.csv";raw=b"not a ZIP workbook";source.write_bytes(raw)
            capture=io.StringIO()
            with contextlib.redirect_stdout(capture):
                code=M.main(["--source-file",str(source),"--source-sha256",hashlib.sha256(raw).hexdigest(),
                             "--decoder-sha256",DECODER,"--output-file",str(out)])
            self.assertEqual(code,2);self.assertEqual(json.loads(capture.getvalue())["status"],"ERROR")
            self.assertFalse(out.exists())

    def test_wrong_decoder_is_refused_before_parsing_source(self):
        with tempfile.TemporaryDirectory()as tmp:
            source,out=Path(tmp)/"native.xlsx",Path(tmp)/"out.csv";raw=b"administrative fixture";source.write_bytes(raw)
            with self.assertRaisesRegex(ValueError,"decoder SHA256"):
                M.project(source,out,expected_sha256=hashlib.sha256(raw).hexdigest(),decoder_sha256="0"*64)
            self.assertFalse(out.exists())

    def test_missing_native_prerequisite_is_skipped(self):
        with tempfile.TemporaryDirectory()as tmp:
            capture=io.StringIO()
            with contextlib.redirect_stdout(capture):
                code=M.main(["--source-file",str(Path(tmp)/"missing.xlsx"),"--source-sha256","0"*64,
                             "--decoder-sha256",DECODER,"--output-file",str(Path(tmp)/"out.csv")])
            self.assertEqual(code,77);self.assertEqual(json.loads(capture.getvalue())["status"],"SKIPPED")


if __name__=="__main__":unittest.main()
