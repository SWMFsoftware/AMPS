#!/usr/bin/env python3
"""Focused tests for VP05's independent DBM and pair-selection logic."""

import csv
from pathlib import Path
import tempfile
import unittest

from reference_solution import AU_KM, VSW_KM_S, arrival, dbm_state
from run_vp05 import select_pairs


class DBMReferenceTests(unittest.TestCase):
    def test_equal_wind_speed_is_exactly_ballistic(self) -> None:
        radius, speed = dbm_state(3600.0, 0.5 * AU_KM, VSW_KM_S)
        self.assertAlmostEqual(radius, 0.5 * AU_KM + VSW_KM_S * 3600.0)
        self.assertEqual(speed, VSW_KM_S)

    def test_fast_and_slow_speeds_relax_toward_wind(self) -> None:
        self.assertTrue(VSW_KM_S < dbm_state(86400.0, 0.5 * AU_KM, 800.0)[1] < 800.0)
        self.assertTrue(250.0 < dbm_state(86400.0, 0.5 * AU_KM, 250.0)[1] < VSW_KM_S)

    def test_arrival_hits_target(self) -> None:
        time_s, speed = arrival(0.5, 600.0, 1.0)
        radius, evaluated_speed = dbm_state(time_s, 0.5 * AU_KM, 600.0)
        self.assertAlmostEqual(radius / AU_KM, 1.0, places=13)
        self.assertAlmostEqual(speed, evaluated_speed, places=13)


class SelectionTests(unittest.TestCase):
    def test_outer_outcome_does_not_control_selection(self) -> None:
        fields = ["", "event", "lineupcat_id", "event_start_time", "spacecraft", "sc_heliodistance",
                  "sc_heeq_lon", "sc_heeq_lat", "cme_heeq_lon", "speed", "mo_bmax", "mo_bzmin", "catalog_id", "source"]
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "fixture.csv"
            with path.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader()
                base = {name: "" for name in fields}
                writer.writerow({**base, "lineupcat_id": "E", "event_start_time": "2024-01-01T00:00Z", "spacecraft": "PSP", "sc_heliodistance": "0.5", "sc_heeq_lon": "2", "speed": "500"})
                writer.writerow({**base, "lineupcat_id": "E", "event_start_time": "2024-01-03T00:00Z", "spacecraft": "Wind", "sc_heliodistance": "1.0", "sc_heeq_lon": "4", "speed": "300"})
            self.assertEqual(len(select_pairs(path)), 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
