#!/usr/bin/env python3
"""Mutation tests for the strict Stage-8 runtime-log contract."""

from __future__ import annotations

import math
import pathlib
import tempfile
import unittest

import test_stage8_total_reaction_forcing_logs as stage8


EVENTS = 10_000


def make_log(ranks: int) -> str:
    homogeneous_tau = 0.5 * 0.8 * 2
    homogeneous_survival = math.exp(-homogeneous_tau)
    homogeneous_reacted = -math.expm1(-homogeneous_tau)
    constant_mean = stage8.constant_expected_fraction(homogeneous_tau)

    first_fraction = 0.35
    tau1 = 0.20
    tau2 = 1.40
    first_probability = -math.expm1(-tau1) / -math.expm1(-(tau1 + tau2))
    first_count = round(EVENTS * first_probability)
    piecewise_mean = stage8.piecewise_expected_fraction(
        first_fraction, tau1, tau2
    )

    decay_tau = 2.29
    decay_survival = math.exp(-decay_tau)
    decay_reacted = -math.expm1(-decay_tau)
    return "\n".join(
        (
            "Total-reaction forcing validation: "
            f"events/fixture={EVENTS} ranks={ranks} "
            "deterministic-global-events=on "
            "branch-mode=deterministic-two-branch",
            "PASS total-reaction-forcing homogeneous parent=0 "
            f"events={EVENTS} steps=2 hazard=5.00000000000000000e-1 "
            "duration=8.00000000000000000e-1 "
            f"optical-depth={homogeneous_tau:.17e} "
            f"analytic-survival={homogeneous_survival:.17e} "
            f"represented-survival={homogeneous_survival:.17e} "
            f"analytic-reacted={homogeneous_reacted:.17e} "
            f"represented-reacted={homogeneous_reacted:.17e} "
            "closure-error=0.00000000000000000e+0",
            "PASS total-reaction-forcing position fixture=constant parent=0 "
            f"samples={EVENTS} tau-first={homogeneous_tau:.17e} "
            "tau-second=0.00000000000000000e+0 "
            "first-fraction=1.00000000000000000e+0 "
            f"first-region-count={EVENTS} "
            "first-region-expected=1.00000000000000000e+0 "
            f"mean-fraction={constant_mean:.17e} "
            f"expected-mean-fraction={constant_mean:.17e} "
            "cdf-distance=0.00000000000000000e+0 endpoint-violations=0",
            "PASS total-reaction-forcing position fixture=piecewise parent=0 "
            f"samples={EVENTS} tau-first={tau1:.17e} tau-second={tau2:.17e} "
            f"first-fraction={first_fraction:.17e} "
            f"first-region-count={first_count} "
            f"first-region-expected={first_probability:.17e} "
            f"mean-fraction={piecewise_mean:.17e} "
            f"expected-mean-fraction={piecewise_mean:.17e} "
            "cdf-distance=0.00000000000000000e+0 endpoint-violations=0",
            "PASS total-reaction-forcing decay parent=0 "
            f"events={EVENTS} steps=5 total-optical-depth={decay_tau:.17e} "
            f"analytic-survival={decay_survival:.17e} "
            f"represented-survival={decay_survival:.17e} "
            f"analytic-reacted={decay_reacted:.17e} "
            f"represented-reacted={decay_reacted:.17e} "
            "closure-error=0.00000000000000000e+0",
            "PASS total-reaction-forcing composition reaction=C2_2C parent=12 "
            "packets=8 products=2 occurrence-factor=4.00000000000000000e-1 "
            "channel-likelihood=2.50000000000000000e-1 "
            "expected-correction=1.00000000000000000e-1 "
            "packet-correction-sum=1.00000000000000000e-1 "
            f"emitted={16 * ranks} momentum-error=0.00000000000000000e+0 "
            "energy-error=0.00000000000000000e+0",
            "PASS total-reaction-forcing production-transaction "
            f"ranks={ranks} parents={ranks} survivors={ranks} "
            f"reacted-branches={ranks} emitted={16 * ranks} "
            f"processed={17 * ranks} remaining-moves={16 * ranks} "
            f"capacity-rollbacks={ranks} byte-rollbacks={ranks} "
            f"rng-rollbacks={ranks} path-rollbacks={ranks} "
            "double-chemistry=0 correction-error=0.00000000000000000e+0 "
            "position-error=0.00000000000000000e+0",
            "PASS total-reaction-forcing mover-production "
            f"ranks={ranks} mpi-face={0 if ranks == 1 else 1} "
            "face-transfers=18 initial-pending=18 continuation-moves=18 "
            "sweeps=1 final-pending=0 product-events=16 sink-events=1 "
            "wrong-owner-events=0 final-particles=17 survivors=1 products=16 "
            "uncleared-markers=0 nonzero-remainders=0 processed=17 "
            "endpoint-checks=16 moved-endpoints=16 "
            "correction-error=0.00000000000000000e+0 "
            "event-position-error=0.00000000000000000e+0 "
            "remaining-fraction-error=0.00000000000000000e+0 "
            "endpoint-error=0.00000000000000000e+0",
            "Total-reaction forcing validation result: PASS "
            "homogeneous=1 piecewise=1 decay=1 composition=1 "
            f"transactions={ranks} mover=1 events={EVENTS}",
            "",
        )
    )


class Stage8ParserTest(unittest.TestCase):
    def parse(self, text: str, ranks: int = 1):
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory) / "run.log"
            path.write_text(text, encoding="utf-8")
            return stage8.parse_log(path, ranks, EVENTS)

    def test_complete_records_pass_at_all_rank_counts(self):
        for ranks in (1, 4, 12):
            with self.subTest(ranks=ranks):
                parsed = self.parse(make_log(ranks), ranks)
                self.assertEqual(set(parsed["positions"]), {"constant", "piecewise"})

    def test_wrong_survival_is_rejected(self):
        text = make_log(1).replace(
            "represented-survival=4.49328964117221563e-01",
            "represented-survival=5.49328964117221563e-01",
            1,
        )
        with self.assertRaisesRegex(ValueError, "represented_survival"):
            self.parse(text)

    def test_reaction_position_escape_is_rejected(self):
        text = make_log(1).replace("endpoint-violations=0", "endpoint-violations=1", 1)
        with self.assertRaisesRegex(ValueError, "position bounds"):
            self.parse(text)

    def test_wrong_occurrence_channel_packet_composition_is_rejected(self):
        text = make_log(1).replace(
            "packet-correction-sum=1.00000000000000000e-1",
            "packet-correction-sum=2.00000000000000000e-1",
        )
        with self.assertRaisesRegex(ValueError, "do not compose"):
            self.parse(text)

    def test_missing_rng_rollback_is_rejected(self):
        text = make_log(1).replace("rng-rollbacks=1", "rng-rollbacks=0")
        with self.assertRaisesRegex(ValueError, "rng_rollbacks"):
            self.parse(text)

    def test_double_chemistry_is_rejected(self):
        text = make_log(1).replace("double-chemistry=0", "double-chemistry=1")
        with self.assertRaisesRegex(ValueError, "twice"):
            self.parse(text)

    def test_transaction_roundoff_must_remain_bounded(self):
        # Cross-rank MPI_MAX residues need not be bitwise identical because
        # every rank uses a distinct local block/recoil seed, but every
        # residue is still subject to the production acceptance bound.
        text = make_log(1).replace(
            "correction-error=0.00000000000000000e+0",
            "correction-error=3.00000000000000000e-12",
        )
        with self.assertRaisesRegex(ValueError, "correction/position"):
            self.parse(text)

    def test_missing_mover_face_crossing_is_rejected(self):
        text = make_log(1).replace("face-transfers=18", "face-transfers=0")
        with self.assertRaisesRegex(ValueError, "face_transfers"):
            self.parse(text)

    def test_missing_mpi_face_crossing_is_rejected(self):
        text = make_log(4).replace("mpi-face=1", "mpi-face=0")
        with self.assertRaisesRegex(ValueError, "MPI-face evidence"):
            self.parse(text, 4)

    def test_uncleared_mover_marker_is_rejected(self):
        text = make_log(1).replace("uncleared-markers=0", "uncleared-markers=1")
        with self.assertRaisesRegex(ValueError, "uncleared_markers"):
            self.parse(text)

    def test_unfinished_product_remainder_is_rejected(self):
        text = make_log(1).replace("nonzero-remainders=0", "nonzero-remainders=1")
        with self.assertRaisesRegex(ValueError, "nonzero_remainders"):
            self.parse(text)

    def test_wrong_moved_product_endpoint_is_rejected(self):
        text = make_log(1).replace(
            "endpoint-error=0.00000000000000000e+0",
            "endpoint-error=3.00000000000000000e-12",
        )
        with self.assertRaisesRegex(ValueError, "endpoint_error exceeds tolerance"):
            self.parse(text)

    def test_wrong_mover_product_correction_is_rejected(self):
        # Replace only the mover record; the earlier transaction has an
        # independent correction-error field and a distinct failure message.
        old = (
            "moved-endpoints=16 correction-error=0.00000000000000000e+0 "
            "event-position-error="
        )
        new = (
            "moved-endpoints=16 correction-error=3.00000000000000000e-12 "
            "event-position-error="
        )
        text = make_log(1).replace(old, new)
        with self.assertRaisesRegex(ValueError, "correction_error exceeds tolerance"):
            self.parse(text)


if __name__ == "__main__":
    unittest.main(verbosity=2)
