"""Publication plotting tests without SPICE, CUDA or archived mission data."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import cops_publication as pub
from cops_io import file_identity


class PublicationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.utc = ["2016-09-30T10:30:00", "2016-09-30T10:34:00", "2016-09-30T10:39:26"]
        self.table = dict(utc=self.utc, time=np.array([pub.parse_utc(u) for u in self.utc]),
                          iPoint=np.arange(3), altitude=np.array([500., 200., .1]),
                          inside=np.zeros(3), pNG=np.array([1., 2., 3.]), pRG=np.ones(3),
                          refNG=np.array([1., 2., 3.]), refRG=np.ones(3),
                          useNG=np.array([1, 0, 1]), useRG=np.ones(3))

    def tearDown(self):
        self.temp.cleanup()

    def test_window_controls_and_safe_ratios(self):
        window = pub.late_window(self.table)
        self.assertEqual(window["points"], [1, 2])
        with self.assertRaisesRegex(ValueError, "together"):
            pub.late_window(self.table, start=self.utc[0])
        with self.assertRaisesRegex(ValueError, "precede"):
            pub.late_window(self.table, start=self.utc[2], stop=self.utc[0])
        with np.errstate(all="raise"):
            result = pub.ratio([1, 2, np.inf], [0, 1, 1], self.table)
        np.testing.assert_equal(result, [np.nan, 2., np.nan])
        bundle = self.root / "bundle.json"
        bundle.write_text(json.dumps({"templates": {"late_ng_not_rg": {"interval_utc": self.utc[:2]}}}))
        self.assertEqual(pub.late_window(self.table, bundle)["points"], [0, 1])

    def test_layout_selection_and_rc_isolation(self):
        import matplotlib.pyplot as plt
        window = pub.late_window(self.table)
        seen = []

        def inspect(fig, prefix, overwrite=False):
            seen.append(fig)
            axes = fig.axes
            self.assertEqual(len(axes), 6 if str(prefix).endswith("_time") else 3)
            self.assertTrue(all(a.get_yscale() == "log" for a in axes))
            if str(prefix).endswith("_altitude"):
                self.assertEqual(axes[0].get_xscale(), "log")
                self.assertGreater(*axes[0].get_xlim())
            self.assertIsNotNone(axes[0].get_legend())
            # Ratio markers distinguish both-selected from at-least-one-excluded.
            ratio_axis = axes[4] if len(axes) == 6 else axes[2]
            self.assertEqual(len(ratio_axis.lines[0].get_ydata()), 2)
            self.assertEqual(len(ratio_axis.lines[1].get_ydata()), 1)
            return []

        original = matplotlib.rcParams["font.size"]
        with patch.object(pub, "save", side_effect=inspect):
            pub.pressure_figures(None, [("fit", self.table)], self.table, window, self.root / "plot")
        self.assertEqual(len(seen), 2)
        self.assertEqual(matplotlib.rcParams["font.size"], original)
        self.assertEqual(plt.get_fignums(), [])

    def test_misaligned_rows_and_overwrite_refused_before_writing(self):
        wrong = dict(self.table, iPoint=np.array([5, 6, 7]))
        window = pub.late_window(self.table)
        with self.assertRaisesRegex(ValueError, "matching iPoint/UTC"):
            pub.pressure_figures(None, [("bad", wrong)], self.table, window, self.root / "plot")
        existing = self.root / "plot_altitude.png"
        existing.write_bytes(b"preserve")
        with self.assertRaises(FileExistsError):
            pub.pressure_figures(None, [("fit", self.table)], self.table, window, self.root / "plot")
        self.assertEqual(list(self.root.iterdir()), [existing])
        self.assertEqual(existing.read_bytes(), b"preserve")

    def test_output_formats_and_nonfatal_short_run(self):
        products = pub.write_result_figures(self.root, self.table, self.table)
        self.assertEqual(len(products), 5)
        for name in products:
            self.assertTrue((self.root / name).is_file())
        self.assertTrue((self.root / "figures/pressures_time.eps").read_bytes().startswith(b"%!PS"))
        status = json.loads((self.root / products[-1]).read_text())
        self.assertEqual(status["status"], "complete")
        for name, identity in status["outputs"].items():
            self.assertEqual(identity["sha256"], file_identity(self.root / name)["sha256"])
        products = pub.write_result_figures(self.root / "short", self.table, self.table,
                                           interval=("2016-09-30T01:00:00", "2016-09-30T02:00:00"))
        self.assertEqual(products, ["figures/pressures.json"])
        self.assertEqual(json.loads((self.root / "short" / products[0]).read_text())["status"], "failed")

    def test_fine_table_hash_base_and_row_checks(self):
        base_path = self.root / "COPS_model.dat"
        base_path.write_text("synthetic base identity")
        base = dict(self.table, path=str(base_path))
        path = self.root / "pressure_comparison.txt"
        path.write_text("# iPoint UTC interval refNG refRG beforeNG beforeRG afterNG afterRG useNG useRG saturatedNG\n" +
                        "".join("%d %s 0 1 1 1 1 2 1 1 1 0\n" % (i, u) for i, u in enumerate(self.utc)))
        meta = dict(schema="liouville_cops.late_interval_fine_tuning.v1", state="complete",
                    optimizer={"success": True}, outputs={path.name: file_identity(path)},
                    definition={"table_sha256": file_identity(base_path)["sha256"]})
        (self.root / "fine_tuning.json").write_text(json.dumps(meta))
        before, after, _ = pub.read_fine_tuning(self.root, base)
        np.testing.assert_equal(before["pNG"], 1)
        np.testing.assert_equal(after["pNG"], 2)
        path.write_text(path.read_text() + " ")
        with self.assertRaisesRegex(ValueError, "does not match manifest"):
            pub.read_fine_tuning(self.root, base)


if __name__ == "__main__":
    unittest.main()
