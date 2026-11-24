"""Check the r11 API projections against the published source tables and figure."""

import os
import statistics
import unittest
from types import SimpleNamespace

import test_protocol


class PhotonBenchmarks(unittest.TestCase):
    requires_authoritative_source = True

    @classmethod
    def setUpClass(cls):
        if os.environ.get("AMPS_CHECK_M3_PHOTON_BENCHMARKS") != "1" or test_protocol.CONTEXT is None:
            raise unittest.SkipTest("requires the selected solar and laboratory source packet")
        cls.db = test_protocol.CONTEXT.python_database()

    def test_charged_water_branches_have_independent_fits(self):
        # Independently read rows 13–16 in the printed Table 1. In
        # particular the O+ branch has H2 as the partner and a much
        # smaller prefactor than OH+. The source ± for A0 uses its exponent.
        rows = self.db.table("moore-mapaye2025:oh-water-fits").rows
        self.assertEqual([r[1] for r in rows], ["OH"] * 4 + ["H2O"] * 7)
        self.assertEqual(rows[-4:], (
            ("13", "H2O", "OH+ + H + e", "2.87", "0.814", "-9", "0.75", "0.035"),
            ("14", "H2O", "O+ + H2 + e", "4.00", "2.349", "-11", "1.15", "0.052"),
            ("15", "H2O", "H+ + OH + e", "3.32", "1.152", "-10", "0.89", "0.039"),
            ("16", "H2O", "H2O+ + e", "4.23", "0.771", "-8", "0.56", "0.025"),
        ))

    def test_posterior_center_separation_matches_published_splitting(self):
        # Only the mean is compared with the paper's rounded central
        # value. The marginal or paired sample dispersions and the
        # publication's uncertainty are distinct statistical quantities.
        table = self.db.table("laurell2025:argon-peak-centers")
        self.assertEqual(table.row_count, 800)
        self.assertEqual(table.rows[0], ("2", "14.90793419", "15.08822729"))
        self.assertEqual(table.rows[-1], ("801", "14.90396782", "15.06659597"))
        self.assertEqual(round(statistics.mean(1000 * (float(b) - float(a)) for _, a, b in table.rows)), 173)

    def test_native_reader_preserves_both_sets_of_tokens(self):
        fields = dict(self.db.identity)
        fixture = SimpleNamespace(root=self.db.root, fields=fields, export_path=fields["export_path"])
        for identifier in ("moore-mapaye2025:oh-water-fits", "laurell2025:argon-peak-centers"):
            with self.subTest(table=identifier):
                native = test_protocol.CONTEXT.native(fixture, "table", identifier, "table-access")
                self.assertTrue(native["query_ok"], native)
                self.assertEqual(native["rows"], [list(r) for r in self.db.table(identifier).rows])
