"""Protect the I22 decision-2 source packet (2026-11-26) in an installed release.

Checks are against facts read independently from the printed sources, not
against the packet's own output: the four OH+ appearance energies printed by
Bodi et al. (2014), page 3060, and the licence printed on that article. Every
other record must stay metadata-only, with no retrievable payload.
"""
from decimal import Decimal
import hashlib
import os
from types import SimpleNamespace
import unittest

import test_protocol

DATASET = "i22-decision2-sources-20261126"
BODI_SHA = "8882498b8d343b4603db90162a56eee4f0e424c94b2bf720e8d9461d8c1eaad2"
DOIS = {
    "tan1978": "10.1016/0301-0104(78)85080-0", "kronebusch1976": "10.1016/0020-7381(76)80088-5",
    "eland1975": "10.1016/0301-0104(75)80036-x", "powis1991": "10.1039/FT9918700921",
    "banna1986": "10.1063/1.450008", "truesdale1982": "10.1063/1.443058", "liu2009": "10.1063/1.3120443",
    "norwood1991": "10.1063/1.461334", "boyarkin2013": "10.1016/j.cplett.2013.03.007",
    "schorghofer2023": "10.3847/PSJ/acf19b", "schorghofer2025": "10.3847/PSJ/ade5b2",
    "joens2001": "10.1021/jp011833u", "sun2020": "10.1063/1674-0068/cjcp2003026",
    "ruscic2015": "10.1021/acs.jpca.5b01346", "bodi2014": "10.1039/C4SC00791C",
    "hultroos2018": "10.1063/1.5055593", "harper2019": "10.1063/1.5091966",
    "sage2010": "10.1080/00268971003596177", "suarez2018": "10.1039/C8CP03725F",
}


class I22Decision2Sources(unittest.TestCase):
    requires_authoritative_source = True

    @classmethod
    def setUpClass(cls):
        if not test_protocol.CONTEXT or os.environ.get("AMPS_CHECK_I22_DECISION2_RELEASE") != "1":
            raise unittest.SkipTest("I22 decision-2 source release was not selected")
        cls.db = test_protocol.CONTEXT.python_database()
        cls.files = {row["id"]: row for row in cls.db.describe(DATASET)["files"]}

    def test_printed_oh_plus_onsets_in_both_interfaces(self):
        # Read from the rendered page 4 (printed p. 3060), independently of the text layer.
        table = "bodi2014:oh-plus-appearance-energy"
        rows = self.db.table(table).rows
        self.assertEqual([(r[0], Decimal(r[1]), Decimal(r[2])) for r in rows], [
            ("Fig. 2a data set", Decimal("18.1177"), Decimal("0.0015")),
            ("Fig. 2b data set", Decimal("18.1183"), Decimal("0.0015")),
            ("Fig. 2c data set", Decimal("18.119"), Decimal("0.002")),
            ("reciprocal error-square weighted average", Decimal("18.1182"), Decimal("0.0009"))])
        identity = dict(self.db.identity)
        fixture = SimpleNamespace(root=self.db.root, fields=identity, export_path=identity["export_path"])
        native = test_protocol.CONTEXT.native(fixture, "table", table, "table-access")
        self.assertTrue(native["query_ok"], native)
        self.assertEqual(native["rows"], [list(row) for row in rows])

    def test_bodi_article_and_non_commercial_policy(self):
        pdf = self.db.source_path("bodi2014:article")
        self.assertEqual(hashlib.sha256(pdf.read_bytes()).hexdigest(), BODI_SHA)
        article = self.files["bodi2014:article"]
        self.assertEqual(article["policy"], "rsc-bodi-2014-cc-by-nc-3.0")
        self.assertEqual(article["rights"]["redistribution"], "non-commercial-with-attribution")
        self.assertIn("Non-commercial use only", article["rights"]["conditions"])
        self.assertEqual(article["quantity_metadata"]["licence"], "CC BY-NC 3.0 Unported")
        derived = self.files["bodi2014:oh-plus-appearance-energy:csv"]
        self.assertEqual(tuple(derived["derived_from"]), ("bodi2014:article",))
        self.assertEqual(derived["policy"], article["policy"])

    def test_every_other_record_is_held_and_inactive(self):
        for key, doi in DOIS.items():
            record = self.files["i22-decision2:" + key + ":crossref"]
            self.assertEqual(record["quantity_metadata"]["doi"], doi)
        allowed = {"bodi2014:article", "bodi2014:oh-plus-appearance-energy:csv"}
        self.assertEqual(len(self.files), 25)
        error = test_protocol.CONTEXT.python_api().DatabaseError
        for identifier, record in self.files.items():
            self.assertFalse(record["physical_activation"])
            if identifier in allowed:
                continue
            self.assertEqual(record["policy"], "metadata-only")
            with self.assertRaises(error) as held:
                self.db.source_path(identifier)
            self.assertEqual(held.exception.code, "PermissionRequired")
        # User-supplied publisher copies are identified by hash only.
        local = self.files["i22-decision2:tan1978:crossref"]["quantity_metadata"]["local_copy"]
        self.assertEqual(local["sha256"], "3e6d839dc3c8ae9fc3e2d6eb259216d28d70b1168a145538f4163a39db13efbd")
        self.assertTrue(local["storage"].startswith("held"))
