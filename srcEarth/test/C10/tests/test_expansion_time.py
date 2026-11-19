"""Unit tests of the C10 maximum-expansion timing criterion.

The criterion is earth_test_runtime.evaluate_extremum_time(extremum="max"),
shared with C9.  The fixture is the 2006-12-14/16 low-rigidity expansion series
of the C10 ROUTINE run of 2026-11-19 (-O0 build): POES/MEPED observed expansion
with the standard error of each interval mean, and the AMPS modeled expansion.
The observed maximum (06:00, 5.635 +- 0.226 deg) is only 0.90 sigma above 00:00
(4.766 +- 0.941 deg), so the observation does not resolve the time of maximum.

Run from srcEarth/test/C10:  python3 -m unittest discover -s tests -v
"""
from __future__ import annotations

import sys
import unittest
from datetime import datetime, timezone
from pathlib import Path

TEST_ROOT = Path(__file__).resolve().parents[2]
if str(TEST_ROOT) not in sys.path:
    sys.path.insert(0, str(TEST_ROOT))

from earth_test_runtime import (  # noqa: E402
    EXTREMUM_TIME_FAIL,
    EXTREMUM_TIME_NOT_EVALUATED,
    EXTREMUM_TIME_PASS,
    evaluate_extremum_time,
)


def utc(text: str) -> datetime:
    return datetime.strptime(text, "%Y-%m-%dT%H:%M").replace(tzinfo=timezone.utc)


TIMES = [utc(t) for t in ("2006-12-14T06:00", "2006-12-14T14:00", "2006-12-14T18:00",
                          "2006-12-14T22:00", "2006-12-15T00:00", "2006-12-15T06:00",
                          "2006-12-15T12:00", "2006-12-15T18:00", "2006-12-16T04:00")]
OBSERVED = [(0.000, 0.244), (0.392, 0.254), (3.733, 0.374), (3.212, 0.665),
            (4.766, 0.941), (5.635, 0.226), (3.964, 0.438), (3.853, 0.297),
            (2.252, 0.787)]
MODELED = [0.000, 0.752, 3.596, 2.777, 5.639, 4.290, 3.389, 2.963, 0.802]


def observed(values=OBSERVED):
    return [(t, v, s, s) for t, (v, s) in zip(TIMES, values)]


def evaluate(obs, model, **kwargs):
    options = dict(extremum="max", tie_tolerance=0.25, min_separation_sigma=2.0,
                   max_error_minutes=180.0)
    options.update(kwargs)
    return evaluate_extremum_time(obs, list(zip(TIMES, model)), **options)


class ExpansionTimeCriterionTests(unittest.TestCase):
    def test_routine_storm_is_not_evaluated(self) -> None:
        result = evaluate(observed(), MODELED)
        self.assertEqual(result.status, EXTREMUM_TIME_NOT_EVALUATED)
        self.assertAlmostEqual(result.observed_extremum_separation_sigma,
                               (5.635 - 4.766) / (0.226 ** 2 + 0.941 ** 2) ** 0.5)
        # The old argmax comparison: model 00:00 vs observed 06:00.
        self.assertEqual(result.argextremum_time_error_minutes, 360.0)
        # The model maximum is clear (1.35 deg margin), so it is the only candidate.
        self.assertEqual(result.modeled_extremum_candidates, [TIMES[4]])

    def test_resolved_maximum_at_the_wrong_time_fails(self) -> None:
        sharp = list(OBSERVED)
        sharp[4] = (3.0, 0.3)          # 06:00 now a clear maximum
        result = evaluate(observed(sharp), MODELED)
        self.assertGreaterEqual(result.observed_extremum_separation_sigma, 2.0)
        self.assertEqual(result.status, EXTREMUM_TIME_FAIL)
        self.assertEqual(result.time_error_minutes, 360.0)

    def test_resolved_maximum_with_tied_model_passes(self) -> None:
        sharp = list(OBSERVED)
        sharp[4] = (3.0, 0.3)
        model = list(MODELED)
        model[5] = model[4] - 0.10     # 06:00 within the tie tolerance of the maximum
        result = evaluate(observed(sharp), model)
        self.assertEqual(result.modeled_extremum_candidates, [TIMES[4], TIMES[5]])
        self.assertEqual(result.status, EXTREMUM_TIME_PASS)
        self.assertEqual(result.time_error_minutes, 0.0)

    def test_maximum_and_minimum_are_mirror_images(self) -> None:
        sharp = list(OBSERVED)
        sharp[4] = (3.0, 0.3)
        as_max = evaluate(observed(sharp), MODELED)
        as_min = evaluate_extremum_time(
            [(t, -v, s, s) for t, (v, s) in zip(TIMES, sharp)],
            [(t, -v) for t, v in zip(TIMES, MODELED)],
            extremum="min", tie_tolerance=0.25, min_separation_sigma=2.0,
            max_error_minutes=180.0)
        self.assertEqual(as_max.status, as_min.status)
        self.assertEqual(as_max.time_error_minutes, as_min.time_error_minutes)
        self.assertAlmostEqual(as_max.observed_extremum_separation_sigma,
                               as_min.observed_extremum_separation_sigma)

    def test_invalid_extremum_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            evaluate(observed(), MODELED, extremum="median")


if __name__ == "__main__":
    unittest.main(verbosity=2)
