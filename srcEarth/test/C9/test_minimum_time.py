#!/usr/bin/env python3
"""Unit tests of the C9 storm-minimum timing criterion (evaluate_minimum_time()).

No AMPS run is needed.  The fixtures are the lowest-rigidity (0.4236 GV) rows of
the 2006-12-14/15 storm: the PAMELA values and uncertainties are those of
reference_C9_pamela_table_s1.csv, and the two model series are the AMPS outputs
of the -O0 and -O3 builds of 2026-11-19, whose argmin swapped between 23:55 and
03:03 although the two values differ by 0.13 and 0.004 deg.

Run:  python3 srcEarth/test/C9/test_minimum_time.py
"""

from __future__ import annotations

import importlib.util
import sys
import unittest
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load_c9():
    spec = importlib.util.spec_from_file_location("run_C9_under_test", HERE / "run_C9.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


C9 = load_c9()


def utc(text: str) -> datetime:
    return datetime.strptime(text, "%Y-%m-%dT%H:%M").replace(tzinfo=timezone.utc)


TIMES = [utc(t) for t in ("2006-12-14T09:49", "2006-12-14T14:31", "2006-12-14T16:05",
                          "2006-12-14T23:55", "2006-12-15T03:03", "2006-12-15T09:19",
                          "2006-12-16T04:08")]
# PAMELA lowest-rigidity cutoff [deg AACGM] and symmetric sigma [deg].
PAMELA = [(60.72, 0.48), (60.60, 0.51), (58.61, 0.60), (55.44, 0.73),
          (54.39, 0.61), (57.71, 0.53), (60.58, 0.49)]
MODEL_O0 = [63.340, 62.706, 60.401, 57.739, 57.865, 59.316, 61.810]
MODEL_O3 = [63.194, 62.772, 60.410, 57.931, 57.927, 59.625, 61.693]

TIE = 0.25     # deg: the default T50 grid step
SIGMA = 2.0
MAX_MIN = 100.0


def observed(values=PAMELA):
    return [(t, v, s, s) for t, (v, s) in zip(TIMES, values)]


def modeled(values):
    return list(zip(TIMES, values))


def evaluate(obs, model, **kwargs):
    options = dict(tie_tolerance_deg=TIE, min_separation_sigma=SIGMA,
                   max_error_minutes=MAX_MIN)
    options.update(kwargs)
    return C9.evaluate_minimum_time(obs, model, **options)


class MinimumTimeCriterionTests(unittest.TestCase):
    def test_rounding_swap_gives_the_same_status(self) -> None:
        # The -O0 and -O3 builds swap the argmin (188 min vs 0 min with the old
        # criterion).  Both minima are within the tie tolerance, so the status and
        # the tie-aware error no longer depend on the build.
        o0 = evaluate(observed(), modeled(MODEL_O0))
        o3 = evaluate(observed(), modeled(MODEL_O3))
        self.assertEqual(o0.argmin_time_error_minutes, 188.0)
        self.assertEqual(o3.argmin_time_error_minutes, 0.0)
        self.assertEqual(o0.status, o3.status)
        self.assertEqual(o0.time_error_minutes, 0.0)
        self.assertEqual(o3.time_error_minutes, 0.0)
        self.assertEqual(o0.modeled_minimum_candidates, [TIMES[3], TIMES[4]])
        self.assertEqual(o3.modeled_minimum_candidates, [TIMES[3], TIMES[4]])

    def test_storm_observation_does_not_resolve_the_minimum(self) -> None:
        # PAMELA: 54.39+-0.61 at 03:03 vs 55.44+-0.73 at 23:55, 1.1 sigma apart.
        result = evaluate(observed(), modeled(MODEL_O3))
        self.assertEqual(result.status, C9.MINIMUM_TIME_NOT_EVALUATED)
        self.assertAlmostEqual(result.observed_minimum_separation_sigma, 1.05 / (0.61**2 + 0.73**2) ** 0.5)
        self.assertNotEqual(result.status, C9.MINIMUM_TIME_PASS)

    def test_resolved_observation_and_tied_model_pass(self) -> None:
        sharp = list(PAMELA)
        sharp[4] = (52.00, 0.61)        # 03:03 now ~3.7 sigma below 23:55
        result = evaluate(observed(sharp), modeled(MODEL_O0))
        self.assertEqual(result.status, C9.MINIMUM_TIME_PASS)
        self.assertEqual(result.time_error_minutes, 0.0)

    def test_clear_model_minimum_at_the_wrong_time_fails(self) -> None:
        # A real timing error: the model minimum at 16:05 is 1 deg below every
        # other interval, so the observed time (03:03) is not a candidate.
        sharp = list(PAMELA)
        sharp[4] = (52.00, 0.61)
        wrong = list(MODEL_O3)
        wrong[2] = 56.90
        result = evaluate(observed(sharp), modeled(wrong))
        self.assertEqual(result.modeled_minimum_candidates, [TIMES[2]])
        self.assertEqual(result.status, C9.MINIMUM_TIME_FAIL)
        self.assertGreater(result.time_error_minutes, MAX_MIN)

    def test_model_minimum_beyond_tolerance_is_not_a_tie(self) -> None:
        sharp = list(PAMELA)
        sharp[4] = (52.00, 0.61)
        model = list(MODEL_O3)
        model[3] = model[4] - 0.30      # 23:55 lower by more than the tolerance
        result = evaluate(observed(sharp), modeled(model))
        self.assertEqual(result.modeled_minimum_candidates, [TIMES[3]])
        self.assertEqual(result.status, C9.MINIMUM_TIME_FAIL)
        self.assertEqual(result.time_error_minutes, 188.0)

    def test_missing_inputs(self) -> None:
        self.assertEqual(evaluate(observed(), []).status, C9.MINIMUM_TIME_FAIL)
        self.assertEqual(evaluate([], modeled(MODEL_O3)).status, C9.MINIMUM_TIME_NOT_EVALUATED)
        one = observed()[:1]
        self.assertEqual(evaluate(one, modeled(MODEL_O3)).status, C9.MINIMUM_TIME_NOT_EVALUATED)
        no_sigma = [(t, v, None, None) for t, v, _, _ in observed()]
        self.assertEqual(evaluate(no_sigma, modeled(MODEL_O3)).status,
                         C9.MINIMUM_TIME_NOT_EVALUATED)


if __name__ == "__main__":
    unittest.main(verbosity=2)
