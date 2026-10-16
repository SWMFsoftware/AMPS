#!/usr/bin/env python3
"""Deterministic U10 checks for the committed I17 RDR processor.

These fixtures are analytical software-verification inputs, not observations
and never count as I17 validation evidence.  They isolate decisions that the
large authoritative RDR run cannot conveniently diagnose after the fact:
locale-independent UTC decoding, native bit-mask selection, radiance-first
aggregation support, monotone native calibration inversion, and the published
Paige sensitivity substitution.  A failure is therefore a U10 kernel FAIL;
missing D03 files are irrelevant to this local test.
"""

from __future__ import annotations

import datetime as dt
import importlib.util
import math
from pathlib import Path
import sys


REPO = Path(__file__).resolve().parents[3]
PROCESSOR = REPO / "srcMoon/test/data_processing/prepare_d03_rdr_eclipse.py"


def load_processor():
    """Load the exact committed processor without copying its formulas."""
    specification = importlib.util.spec_from_file_location(
        "srcmoon_i17_rdr_processor", PROCESSOR)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load {PROCESSOR}")
    module = importlib.util.module_from_spec(specification)
    # Python 3.8 dataclasses resolve postponed annotations through
    # sys.modules while decorating the class, so register before execution.
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> int:
    processor = load_processor()

    parsed = processor.parse_utc("15-Jun-2011", "16:50:00.083")
    require(parsed == dt.datetime(
        2011, 6, 15, 16, 50, 0, 83000, tzinfo=dt.timezone.utc),
        "native SCET UTC was not decoded exactly")

    quality = {
        "activity_flag_allowed": [110],
        "maximum_emission_angle_degrees_exclusive": 15.0,
        "qca_allowed_mask": 0,
        "qge_allowed_mask": 12,
        "qmi_allowed_mask": 2,
    }
    begin = dt.datetime(2011, 6, 15, 16, 50, tzinfo=dt.timezone.utc)
    end = dt.datetime(2011, 6, 15, 23, 40, tzinfo=dt.timezone.utc)
    common = (begin, begin, end, 110, 0.25, 105.5, 2.0)
    # qge=12 retains both documented definitive-but-not-reprocessed bits and
    # qmi=2 retains the physical eclipse marker.  Calibration degradation,
    # predictive geometry, noise, and the exclusive emission limit must fail.
    require(processor.record_is_selected(*common, 0, 12, 2, quality),
            "documented eclipse/definitive-geometry record was rejected")
    require(not processor.record_is_selected(*common, 2, 12, 2, quality),
            "nearest-marker calibration record leaked through qca")
    require(not processor.record_is_selected(*common, 0, 16, 2, quality),
            "predictive pointing record leaked through qge")
    require(not processor.record_is_selected(*common, 0, 12, 32, quality),
            "noisy record leaked through qmi")
    emission_limit = (begin, begin, end, 110, 0.25, 105.5, 15.0)
    require(not processor.record_is_selected(
        *emission_limit, 0, 12, 2, quality),
        "exclusive emission-angle limit was treated as inclusive")
    require(processor.population_count((1 << 0) | (1 << 5) | (1 << 20)) == 3,
            "portable detector-mask population count changed")
    require(processor.mask_values((1 << 0) | (1 << 20), 1) == "1|21",
            "detector identities were not preserved as one-based values")

    # Archived radiance is the independent coordinate.  Duplicate detector
    # records at one rounded radiance must average native Tb before the
    # monotone inverse is formed; the result below is exactly 101 K.
    calibration = processor.NativeCalibration({
        1.0: [100.0 + 102.0, 2.0],
        2.0: [200.0, 1.0],
        **{float(index): [100.0 * index, 1.0]
           for index in range(3, 103)},
    })
    require(calibration.brightness_temperature(1.0) == 101.0,
            "duplicate-radiance native temperatures were not averaged")
    require(calibration.brightness_temperature(1.5) == 150.5,
            "native calibration inverse is not linear between stored nodes")

    bolometric, used = processor.tbol(
        [100.0, 110.0, 120.0, 80.0, 50.0, 30.0, 55.0])
    require(used == [55.0] * 7,
            "published next-longer-channel sensitivity chain changed")
    # 54.988557963459714... is the 50-decimal mpmath evaluation of the
    # dimensionless Planck integral over the seven published partitions.  It
    # is independent of the production C++ Simpson implementation.
    require(math.isclose(bolometric, 54.988557963459714,
                         rel_tol=1.0e-13, abs_tol=0.0),
            "independent Paige bolometric fixture changed")

    # Regression for a real cold-channel condition: a one-sigma C9 interval
    # can extend below zero even though its binned mean is physically valid.
    # Uncertainty propagation is a local Jacobian calculation and must neither
    # feed a non-positive temperature to the Planck operator nor clamp it to a
    # fabricated observation.  Large errors assigned to substituted C3--C8
    # channels have zero local effect by the published Paige rule.
    cold_temperatures = [100.0, 110.0, 120.0, 80.0, 50.0, 30.0, 1.0]
    cold_uncertainties = [1000.0] * 6 + [2.0]
    propagated = processor.tbol_uncertainty(
        cold_temperatures, cold_uncertainties)
    require(math.isfinite(propagated) and propagated > 0.0,
            "cold-C9 uncertainty propagation is not finite and positive")
    only_c9 = processor.tbol_uncertainty(
        cold_temperatures, [0.0] * 6 + [2.0])
    require(math.isclose(propagated, only_c9, rel_tol=1.0e-12),
            "substituted below-threshold channels affected local uncertainty")

    print("PASS I17 RDR processor contract checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
