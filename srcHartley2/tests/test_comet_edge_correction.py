#!/usr/bin/env python3
"""Tests of the spherical outer-boundary (edge) brightness correction.

Self-contained: every input is synthesized in a temporary run directory, so
the tests need neither an AMPS executable nor validation_runs/.
"""

from __future__ import annotations

import importlib.util
import math
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

VISUALIZATION = Path(__file__).resolve().parents[1] / "visualization"


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


EDGE = load_module("comet_edge_correction", VISUALIZATION / "comet_edge_correction.py")

R_EDGE = 1.0e8
CONF = """\
#main
  SpeciesList=H,OH,O,H2O
  ChemistryReaction=H2O_H_OH | H2O | H + OH | HELIOCENTRIC | 0 | 0 | 0 | 1 | 4.04
  ChemistryReaction=H2O_ION | H2O | NONE | HELIOCENTRIC | 0 | 0 | 0 | 1 | 0
  ChemistryReaction=OH_O_H | OH | O + H | HELIOCENTRIC | 0 | 0 | 0 | 1 | 1.27
  ChemistryReaction=OH_ION | OH | NONE | HELIOCENTRIC | 0 | 0 | 0 | 1 | 0
  ChemistryReaction=H_ION | H | NONE | HELIOCENTRIC | 0 | 0 | 0 | 1 | 0
#endmain
"""


def runtime_input(oh_loss: float, h_loss: float) -> str:
    # H2O decays within about 1 km of the 2.5-km nucleus, so H starts from
    # (nearly) the origin; with a negligible H loss rate the H density is
    # then ~1/r^2 out to R_BIG and the edge factor has a closed form.
    return f"""\
#begin application: comet
  nucleus.radius_m = 2500
  mesh.outer_radius_m = {R_EDGE:.6e}
  #begin case: 1
    heliocentric_distance_au = 1.0
    source.H2O.rate_s-1 = 1.0e28
    chemistry.rate.H2O_H_OH_1au_s-1 = 1.0
    chemistry.rate.H2O_ION_1au_s-1 = 0.0
    chemistry.rate.OH_O_H_1au_s-1 = 0.0
    chemistry.rate.OH_ION_1au_s-1 = {oh_loss:.6e}
    chemistry.rate.H_ION_1au_s-1 = {h_loss:.6e}
  #end case
#end application
"""


def make_run(root: Path, oh_loss: float = 1.0e-5, h_loss: float = 1.0e-13) -> Path:
    run = root / "case1_run"
    (run / "PT/plots").mkdir(parents=True)
    (run / "i3.in").write_text(runtime_input(oh_loss, h_loss), encoding="utf-8")
    (run / "i3.conf").write_text(CONF, encoding="utf-8")
    radii = np.geomspace(2500.0, 1.0e12, 400)
    rows = "\n".join(f"{r:.6e} 1.0e10 1000.0" for r in radii)
    (run / "PT/plots/pic.H2O.s=3.out=0.dat").write_text(
        'VARIABLES="X", "Number Density[1/m^3]", "Vx [m/s]"\n' + rows + "\n", encoding="utf-8"
    )
    return run


def inverse_square_factor(b: float, r_edge: float, r_big: float) -> float:
    """N_full/N_trunc for n ~ 1/r^2 (full column ends at r_big)."""

    return math.atan(math.sqrt(r_big**2 - b**2) / b) / math.atan(math.sqrt(r_edge**2 - b**2) / b)


class ParsingTest(unittest.TestCase):
    def test_products(self) -> None:
        self.assertEqual(EDGE.parse_products("H + OH"), ["H", "OH"])
        self.assertEqual(EDGE.parse_products("2*H + O"), ["H", "H", "O"])
        self.assertEqual(EDGE.parse_products("NONE"), [])
        with self.assertRaises(ValueError):
            EDGE.parse_products("H + ")

    def test_case_overlay_and_inputs(self) -> None:
        with tempfile.TemporaryDirectory(prefix="edge-correction-") as temporary:
            run = make_run(Path(temporary))
            inputs = EDGE.load_run_inputs(run / "i3.in", run / "i3.conf", 1)
            self.assertEqual(inputs.outer_radius_m, R_EDGE)
            self.assertEqual(inputs.nucleus_radius_m, 2500.0)
            self.assertEqual(inputs.sources, {"H2O": 1.0e28})
            self.assertAlmostEqual(inputs.loss_rate("OH"), 1.0e-5)
            self.assertEqual([c.products for c in inputs.channels["H2O"]], [["H", "OH"], []])
            with self.assertRaises(ValueError):
                EDGE.load_run_inputs(run / "i3.in", run / "i3.conf", 2)


class ChemistryVariantTest(unittest.TestCase):
    """Groups and table rates (srcHartley2/CHEMISTRY_VARIANTS.md, section 9)."""

    VARIANT_CONF = CONF.replace(
        "  ChemistryReaction=OH_O_H | OH | O + H | HELIOCENTRIC | 0 | 0 | 0 | 1 | 1.27\n",
        "  ChemistryReaction=OH_O_H | OH | O + H | HELIOCENTRIC | 0 | 0 | 0 | 1 | 1.27"
        " | GROUP=OH_LOSS:STATE_RESOLVED\n").replace(
        "  ChemistryReaction=OH_ION | OH | NONE | HELIOCENTRIC | 0 | 0 | 0 | 1 | 0\n",
        "  ChemistryReaction=OH_ION | OH | NONE | HELIOCENTRIC | 0 | 0 | 0 | 1 | 0"
        " | GROUP=OH_LOSS:STATE_RESOLVED\n"
        "  ChemistryGroup=OH_LOSS | DEFAULT=STATE_RESOLVED | VARIANTS=STATE_RESOLVED,SA1982_VH\n"
        "  ChemistryReaction=OH_O_H_SA1982 | OH | O + H | HELIOCENTRIC | 0 | 0 | 0 | 1 | 0.36"
        " | GROUP=OH_LOSS:SA1982_VH | RATE=TABLE:OH_TAU_SA1982:DISSOCIATION\n"
        "  ChemistryReaction=OH_ION_SA1982 | OH | NONE | HELIOCENTRIC | 0 | 0 | 0 | 1 | 0"
        " | GROUP=OH_LOSS:SA1982_VH | RATE=TABLE:OH_TAU_SA1982:IONIZATION\n"
        "  ChemistryRateTable=OH_TAU_SA1982 | MODEL=SA1982_OH_LIFETIME | SCALE_S=0.833e5"
        " | ION_FRACTION=0.05\n")

    def make_variant_run(self, root: Path, variant_line: str) -> Path:
        run = make_run(root)
        (run / "i3.conf").write_text(self.VARIANT_CONF, encoding="utf-8")
        (run / "oh.csv").write_text(
            "# test table\ntau_fuv_s,4.24e5\nv_h_km_s,ohtaus_1e5_s\n-2,2.66\n0,2.75\n2,2.80\n",
            encoding="utf-8")
        text = (run / "i3.in").read_text(encoding="utf-8").replace(
            "  #begin case: 1\n",
            f"{variant_line}  chemistry.table.OH_TAU_SA1982.file = oh.csv\n  #begin case: 1\n"
            "    brightness.heliocentric_radial_velocity_km_s = 0.4\n")
        (run / "i3.in").write_text(text, encoding="utf-8")
        return run

    def test_build_input_records(self) -> None:
        with tempfile.TemporaryDirectory(prefix="edge-variants-") as temporary:
            run = self.make_variant_run(Path(temporary), "")
            build = EDGE.parse_build_input(run / "i3.conf")
            self.assertEqual(build.groups["OH_LOSS"], ("STATE_RESOLVED", ["STATE_RESOLVED", "SA1982_VH"]))
            self.assertEqual(build.tables["OH_TAU_SA1982"], ("SA1982_OH_LIFETIME", 0.833e5, 0.05))
            # The excess energy is field 9, not the trailing GROUP/RATE field.
            self.assertEqual(build.reactions["OH_O_H"].energy_ev, 1.27)
            self.assertEqual(build.reactions["OH_O_H_SA1982"].energy_ev, 0.36)
            self.assertEqual(EDGE.parse_reactions(run / "i3.conf")["OH_ION"], ("OH", [], 0.0))

    def test_default_variant_uses_input_rates(self) -> None:
        with tempfile.TemporaryDirectory(prefix="edge-variants-") as temporary:
            run = self.make_variant_run(Path(temporary), "")
            inputs = EDGE.load_run_inputs(run / "i3.in", run / "i3.conf", 1)
            self.assertEqual([c.reaction_id for c in inputs.channels["OH"]], ["OH_O_H", "OH_ION"])
            self.assertAlmostEqual(inputs.loss_rate("OH"), 1.0e-5)

    def test_table_variant_from_input_and_override(self) -> None:
        k_d = 1.0 / 4.24e5 + 1.0 / (0.833e5 * (2.75 + 0.2 * 0.05))
        with tempfile.TemporaryDirectory(prefix="edge-variants-") as temporary:
            run = self.make_variant_run(Path(temporary), "  chemistry.variant.OH_LOSS = SA1982_VH\n")
            inputs = EDGE.load_run_inputs(run / "i3.in", run / "i3.conf", 1)
            channels = {c.reaction_id: c for c in inputs.channels["OH"]}
            self.assertEqual(sorted(channels), ["OH_ION_SA1982", "OH_O_H_SA1982"])
            self.assertAlmostEqual(channels["OH_O_H_SA1982"].rate_s / k_d, 1.0, delta=1e-12)
            self.assertAlmostEqual(channels["OH_ION_SA1982"].rate_s / (0.05 * k_d), 1.0, delta=1e-12)
            self.assertAlmostEqual(channels["OH_O_H_SA1982"].energy_j, 0.36 * EDGE.EV)
            # An explicit selection (e.g. from the output header) wins.
            inputs = EDGE.load_run_inputs(run / "i3.in", run / "i3.conf", 1,
                                          {"OH_LOSS": "STATE_RESOLVED"})
            self.assertEqual([c.reaction_id for c in inputs.channels["OH"]], ["OH_O_H", "OH_ION"])
            with self.assertRaises(ValueError):
                EDGE.load_run_inputs(run / "i3.in", run / "i3.conf", 1, {"OH_LOSS": "NOPE"})

    def test_variant_from_output_header(self) -> None:
        with tempfile.TemporaryDirectory(prefix="edge-variants-") as temporary:
            header = Path(temporary) / "pic.OH.brightness.s=3.out=0.dat"
            header.write_text(
                'TITLE="x"\nDATASETAUXDATA Chemistry.Variant.OH_LOSS="SA1982_VH; source=command line"\n'
                'DATASETAUXDATA Chemistry.Variant.OTHER="A; source=compiled default"\nZONE T="OH"\n',
                encoding="utf-8")
            self.assertEqual(EDGE.variants_from_output(header),
                             {"OH_LOSS": "SA1982_VH", "OTHER": "A"})


class ModelPiecesTest(unittest.TestCase):
    def test_recoil_conserves_momentum_and_energy(self) -> None:
        rng = np.random.default_rng(1)
        masses = np.array([1.008, 1.008, 15.999]) * EDGE.AMU
        energy = 0.697 * EDGE.EV
        u = EDGE.recoil(masses, energy, 1000, rng)
        momentum = (masses[None, :, None] * u).sum(axis=1)
        kinetic = 0.5 * (masses[None, :, None] * u**2).sum(axis=(1, 2))
        self.assertLess(np.abs(momentum).max(), 1.0e-30)
        np.testing.assert_allclose(kinetic, energy, rtol=1.0e-12)

    def test_shell_projection_of_inverse_square_density(self) -> None:
        edges = np.geomspace(1.0e3, 3.0e9, 3001)
        volume = 4.0 / 3.0 * np.pi * np.diff(edges**3)
        density = 4.0 * np.pi * np.diff(edges) / volume  # shell average of 1/r^2
        for b in (1.0e4, 1.0e6, 1.0e8, 2.0e9):
            exact = 2.0 / b * math.atan(math.sqrt(3.0e9**2 - b**2) / b)
            self.assertAlmostEqual(EDGE.project_shells(edges, density, [b])[0] / exact, 1.0, delta=2.0e-4)

    def test_factor_is_undefined_at_and_beyond_the_edge(self) -> None:
        b = np.array([1.0e7, R_EDGE, 2.0 * R_EDGE])
        factor = EDGE.correction_factor(np.ones(3), np.array([0.5, 0.5, 0.0]), b, R_EDGE)
        self.assertAlmostEqual(factor[0], 2.0)
        self.assertTrue(np.isnan(factor[1:]).all())


class EndToEndTest(unittest.TestCase):
    def test_inverse_square_limit_and_short_lived_species(self) -> None:
        b = np.array([1.0e6, 1.0e7, 5.0e7])
        with tempfile.TemporaryDirectory(prefix="edge-correction-") as temporary:
            run = make_run(Path(temporary))
            factors, inputs = EDGE.correction_for_run(run, 1, 0, b, parents=200_000, seed=7)
        for bb, got in zip(b, factors["H"]):
            expected = inverse_square_factor(bb, R_EDGE, EDGE.R_BIG_M)
            self.assertAlmostEqual(got / expected, 1.0, delta=0.03, msg=f"b={bb:g}")
        # OH lives ~1e5 s at ~1.5 km/s (~1.5e8 m scale before decay), so at
        # this small edge it is truncated too, but much less than H.
        self.assertTrue(np.all(factors["OH"] >= 1.0))
        self.assertTrue(np.all(factors["OH"] < factors["H"]))

    def test_large_edge_needs_no_correction(self) -> None:
        b = np.array([1.0e6, 1.0e7])
        with tempfile.TemporaryDirectory(prefix="edge-correction-") as temporary:
            run = make_run(Path(temporary), oh_loss=1.0e-3, h_loss=1.0e-3)
            factors, _ = EDGE.correction_for_run(run, 1, 0, b, parents=50_000, seed=3)
        # Lifetimes of 1e3 s give scale lengths of 1.5e6 m (OH) and 2.6e7 m
        # (H at 26 km/s), so the 1e8-m edge removes almost nothing.
        np.testing.assert_allclose(factors["OH"], 1.0, atol=1.0e-6)
        self.assertLess(abs(factors["H"][0] - 1.0), 0.01)


if __name__ == "__main__":
    unittest.main()
