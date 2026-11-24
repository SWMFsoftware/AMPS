"""Test authenticated November 24 source views through the real catalogue."""

import json
import os
from pathlib import Path
from types import SimpleNamespace
import unittest
from zipfile import ZipFile

import test_protocol


class RealArchive(unittest.TestCase):
    requires_authoritative_source = True

    @classmethod
    def setUpClass(cls):
        supplied = os.environ.get("AMPS_M3_ARCHIVE_PACKET")
        if not supplied or test_protocol.CONTEXT is None:
            raise unittest.SkipTest("explicit real archive packet and shared runner required")
        cls.packet_path = Path(supplied)
        cls.packet = json.loads(cls.packet_path.read_text())
        cls.database = test_protocol.CONTEXT.python_database()

    def origin(self, identifier):
        row, = [source for source in self.packet["sources"] if source["id"] == identifier]
        path = self.packet_path.parent / row["path"]
        self.assertEqual(test_protocol.digest(path), row["sha256"])
        return path

    def test_whole_he_cx_view_matches_native_archive(self):
        identifier = "swaczyna-4555716:native-zip"
        origin = self.origin(identifier)
        self.assertEqual(test_protocol.digest(self.database.source_path(identifier)),
                         test_protocol.digest(origin))
        index = self.database.table("swaczyna-4555716:he-ct-index")
        angle_rows = self.database.table("swaczyna-4555716:he-ct-angles").rows
        self.assertEqual(index.row_count, 661)
        self.assertEqual(len(angle_rows), 507648)
        # Independently reconstruct each row from its saved member, preserving
        # the original angle order and spelling. Never use a projected CSV as
        # an oracle for another copy of itself.
        offset = 0
        with ZipFile(origin) as native:
            for row in index.dictionaries():
                member = row["native_member"]
                native_rows = [line.split() for line in native.read(member).decode("ascii").splitlines()
                               if not line.startswith("#")]
                count = int(row["angular_rows"])
                expected = [(row["log10_e_cm_x1000"], *cells) for cells in native_rows]
                self.assertEqual(list(angle_rows[offset:offset + count]), expected, member)
                offset += count
        self.assertEqual(offset, len(angle_rows))

    def test_native_cpp_and_python_index_agree(self):
        context = test_protocol.CONTEXT
        lock = json.loads(Path(context.lock).read_text())
        fixture = SimpleNamespace(root=self.database.root, fields=lock, export_path=lock["export_path"])
        identifier = "swaczyna-4555716:he-ct-index"
        native = context.native(fixture, "table", identifier, "table-access")
        self.assertTrue(native["query_ok"], native)
        self.assertEqual(native["rows"], [list(row) for row in self.database.table(identifier).rows])

    def test_held_sources_and_iaea_corroboration(self):
        live = self.origin("iaea-d119990-20261124:live-html")
        self.assertEqual(test_protocol.digest(self.database.source_path("iaea-d119990-20261124:live-html")),
                         test_protocol.digest(live))
        api = test_protocol.CONTEXT.python_api()
        held = (source for source in self.packet["sources"] if source["policy"] == "metadata-only")
        count = 0
        for source in held:
            self.origin(source["id"])
            with self.assertRaises(api.DatabaseError) as refused:
                self.database.source_path(source["id"])
            self.assertEqual(refused.exception.code, "PermissionRequired")
            count += 1
        self.assertEqual(count, 10)
