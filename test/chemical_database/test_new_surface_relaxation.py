"""Source-parser fixtures only, not laboratory quenching observations."""
from pathlib import Path
import hashlib
import importlib.util
import io
import tempfile
import unittest
import zipfile

ROOT=Path(__file__).resolve().parents[2]
P=ROOT/"utility/chemical_database/new_surface_relaxation.py"
S=importlib.util.spec_from_file_location("new_relaxation",P)
M=importlib.util.module_from_spec(S); S.loader.exec_module(M)


class Protocol(unittest.TestCase):
    def test_duplicate_native_headers_are_preserved_by_position(self):
        rows=M.decode_table(b"N\tprobability\tN\tprobability\n1\t--\t1\t0\n", "fixture.txt","tab",4)
        self.assertEqual(rows[0], ["fixture.txt","1","1","N"])
        self.assertEqual(rows[2][-1], "N")
        self.assertEqual(rows[5][-1], "--")
        self.assertEqual(rows[7][-1], "0")

    def test_blank_field_is_not_zero_or_omitted(self):
        rows=M.decode_table(b"a\tb\n1\t\n", "fixture.txt","tab",2)
        self.assertEqual(rows[-1], ["fixture.txt","2","2",""])

    def test_signed_and_exponent_tokens_unchanged(self):
        rows=M.decode_table(b" -1.00 +2.0E-8 0 4 5 6\n", "fixture.txt","whitespace",6)
        self.assertEqual([r[-1]for r in rows][:3], ["-1.00","+2.0E-8","0"])

    def test_native_whitespace_block_separator_has_explicit_row_identity(self):
        """A source blank row is neither deleted nor six manufactured zeros."""
        rows=M.decode_table(b"1 2 3 4 5 6\n\n7 8 9 10 11 12\n", "fixture.txt","whitespace",6)
        self.assertEqual(len(rows),13)
        self.assertEqual(rows[6],["fixture.txt","2","0",""])
        self.assertEqual(rows[7],["fixture.txt","3","1","7"])

    def test_ragged_or_blank_row_refused(self):
        for raw in [b"1\t2\n3\n", b"1\t2\n\n", b""]:
            with self.subTest(raw=raw), self.assertRaises(ValueError):
                M.decode_table(raw,"fixture.txt","tab",2)

    def archive(self, extra=False):
        """Minimal administrative archive; its numbers are explicitly fixtures."""
        buffer=io.BytesIO()
        with zipfile.ZipFile(buffer,"w")as archive:
            archive.writestr("README.txt","definitions of the data columns; radians")
            for name in M.MEMBERS[7]:archive.writestr(name,"1 2 3 4 5 6\n")
            if extra:archive.writestr("unreviewed.txt","7")
        return buffer.getvalue()

    def test_complete_archive_and_source_preserved(self):
        with tempfile.TemporaryDirectory()as tmp:
            source=Path(tmp)/"source.zip"; out=Path(tmp)/"out"
            raw=self.archive();source.write_bytes(raw)
            value=M.project(source,out,archive_number=7,expected_sha256=hashlib.sha256(raw).hexdigest())
            self.assertEqual(len(value["tables"]),3)
            self.assertEqual(sum(t["cell_count"]for t in value["tables"]),18)
            self.assertEqual(source.read_bytes(),raw)
            with self.assertRaises(FileExistsError):
                M.project(source,out,archive_number=7,expected_sha256=hashlib.sha256(raw).hexdigest())

    def test_unreviewed_members_refused_before_any_output(self):
        with tempfile.TemporaryDirectory()as tmp:
            source=Path(tmp)/"source.zip";out=Path(tmp)/"out"
            raw=self.archive(extra=True);source.write_bytes(raw)
            with self.assertRaisesRegex(ValueError,"member set"):
                M.project(source,out,archive_number=7,expected_sha256=hashlib.sha256(raw).hexdigest())
            self.assertFalse(out.exists())

    def test_hash_mismatch_refused(self):
        with tempfile.TemporaryDirectory()as tmp:
            source=Path(tmp)/"source.zip";out=Path(tmp)/"out";source.write_bytes(self.archive())
            with self.assertRaisesRegex(ValueError,"SHA256"):
                M.project(source,out,archive_number=7,expected_sha256="0"*64)
            self.assertFalse(out.exists())


if __name__=="__main__":unittest.main()
