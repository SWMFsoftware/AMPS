#!/usr/bin/env python3
"""Manufactured-reference tests for the Step-10 live/replay comparator."""

from __future__ import annotations

import tempfile
from pathlib import Path

from compare_cutoff_access import ComparisonFailure, compare_products, parse_product


HEADER = """TITLE=\"manufactured cutoff/access reference\"
AUXDATA SNAPSHOT_ID=\"field-v1-reference\"
AUXDATA SNAPSHOT_EPOCH_UTC=\"2024-05-10T12:00:00.000000000Z\"
AUXDATA SNAPSHOT_MESH_REVISION=\"swmf-mesh-v1-reference\"
AUXDATA SNAPSHOT_CONTENT_FINGERPRINT=\"swmf-state-v1-reference\"
AUXDATA OUTER_BOUNDARY_POLICY=\"BOX\"
VARIABLES=\"id\",\"Rc_GV\",\"access_state\",\"value_with_nan\"
ZONE T=\"points\" I=2 F=POINT
"""


def require_failure(callable_object, label: str) -> None:
    try:
        callable_object()
    except ComparisonFailure:
        return
    raise AssertionError(f"comparison unexpectedly passed: {label}")


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="step10_compare_") as directory:
        root = Path(directory)
        live = root / "live.dat"
        replay = root / "replay.dat"
        live.write_text(HEADER + "0 1.25 1 nan\n1 2.5 0 3.0\n", encoding="utf-8")
        replay.write_text(live.read_text(encoding="utf-8"), encoding="utf-8")

        summary = compare_products(parse_product(live), parse_product(replay),
                                   rtol=0.0, atol=0.0)
        assert summary["RESULT"] == "PASS"
        assert summary["max_abs_difference"] == 0.0
        print("PASS S10-U06 exact identical-snapshot artifact reference")

        replay.write_text(HEADER + "0 1.25 1 nan\n1 2.500001 0 3.0\n", encoding="utf-8")
        require_failure(
            lambda: compare_products(parse_product(live), parse_product(replay),
                                     rtol=0.0, atol=0.0),
            "changed cutoff under exact gate",
        )
        # An explicit tolerance is honored but never installed implicitly.
        compare_products(parse_product(live), parse_product(replay),
                         rtol=1.0e-6, atol=0.0)
        print("PASS S10-U07 exact default and explicit-only tolerance")

        replay.write_text(
            live.read_text(encoding="utf-8").replace(
                "field-v1-reference", "field-v1-different", 1
            ),
            encoding="utf-8",
        )
        require_failure(
            lambda: compare_products(parse_product(live), parse_product(replay),
                                     rtol=0.0, atol=0.0),
            "different snapshot identity",
        )

        replay.write_text(
            live.read_text(encoding="utf-8").replace(
                'AUXDATA SNAPSHOT_MESH_REVISION="swmf-mesh-v1-reference"\n', ""
            ),
            encoding="utf-8",
        )
        require_failure(
            lambda: compare_products(parse_product(live), parse_product(replay),
                                     rtol=0.0, atol=0.0),
            "missing replay provenance",
        )

        replay.write_text(
            live.read_text(encoding="utf-8").replace(
                'ZONE T="points" I=2 F=POINT',
                'ZONE T="points" I=1 J=2 F=POINT',
            ),
            encoding="utf-8",
        )
        require_failure(
            lambda: compare_products(parse_product(live), parse_product(replay),
                                     rtol=0.0, atol=0.0),
            "changed zone topology",
        )

        shue_header = HEADER.replace(
            'AUXDATA OUTER_BOUNDARY_POLICY="BOX"',
            'AUXDATA OUTER_BOUNDARY_POLICY="SHUE"\n'
            'AUXDATA SHUE_R0_RE="9.8"\n'
            'AUXDATA SHUE_ALPHA="0.62"\n'
            'AUXDATA SHUE_TAIL_CAP_X_M="-2.5e8"',
        )
        live.write_text(shue_header + "0 1.25 1 nan\n1 2.5 0 3.0\n", encoding="utf-8")
        replay.write_text(
            live.read_text(encoding="utf-8").replace(
                'AUXDATA SHUE_ALPHA="0.62"', 'AUXDATA SHUE_ALPHA="0.63"'
            ),
            encoding="utf-8",
        )
        require_failure(
            lambda: compare_products(parse_product(live), parse_product(replay),
                                     rtol=0.0, atol=0.0),
            "changed resolved Shue surface",
        )
        print("PASS S10-U08 stale/missing provenance fails closed")

    print("RESULT: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
