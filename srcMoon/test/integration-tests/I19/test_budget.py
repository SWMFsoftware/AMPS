#!/usr/bin/env python3
"""Analytic transfer identities and corruption checks, never observations."""
import csv
from decimal import Decimal
import importlib.util
from pathlib import Path
import tempfile
import unittest

from analyze_budget import evaluate


class BudgetInfrastructure(unittest.TestCase):
    def fixture(self,path,change=None):
        """Known exact transfers: initial 10 gas+20 surface, +5 external,
        2 adsorbed, 3 released internally and 1 escaped =>15 gas+19 surface.
        Test values are bookkeeping fixtures, not simulated/observed science.
        """
        events=("emission","adsorption","reflection","photolytic_loss",
          "photolytic_product","escape","numerical_removal","reservoir_supply",
          "surface_deposit","surface_release")
        rows=[]
        for species,mass in (("Na","3.82e-26"),("Ar","6.6421562664e-26")):
            base={"completed_steps":0,"species":species,"mass_kg":mass,
              "initial_count":10,"gas_count":10,"initial_gas_atoms":10,"gas_atoms":10,
              "initial_surface_atoms":20,"surface_atoms":20,"count_residual":0,
              "gas_residual_atoms":0,"surface_residual_atoms":0,"combined_residual_atoms":0,
              "combined_residual_kg":0,"roundoff_bound_atoms":"1e-10",
              "arithmetic_epsilon":"1.08420217248550443401e-19"}
            for event in events:
                base[event+"_count"]=0;base[event+"_atoms"]=0
            rows.append(base)
            final=dict(base,completed_steps=1,gas_count=15,gas_atoms=15,surface_atoms=19)
            for event,n in (("emission",8),("adsorption",2),("surface_deposit",2),("surface_release",3),("escape",1)):
                final[event+"_count"]=n;final[event+"_atoms"]=n
            if change:change(final)
            rows.append(final)
        with path.open("w",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)

    def test_internal_transfers_cancel(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"ledger.csv";self.fixture(p)
            self.assertEqual(evaluate(p)["status"],"PASS")

    def test_false_integer_receipt_is_error(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"ledger.csv";self.fixture(p,lambda r:r.update(gas_count=16))
            with self.assertRaisesRegex(ValueError,"integer residual"):evaluate(p)

    def test_unexplained_gas_mass_is_failure(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"ledger.csv"
            self.fixture(p,lambda r:r.update(gas_atoms=16,gas_residual_atoms=1,
                combined_residual_atoms=1,combined_residual_kg=r["mass_kg"]))
            self.assertEqual(evaluate(p)["status"],"FAIL")

    def test_duplicated_replicated_surface_is_failure(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"ledger.csv"
            self.fixture(p,lambda r:r.update(surface_atoms=38,surface_residual_atoms=19,
                combined_residual_atoms=19,combined_residual_kg=str(19*Decimal(r["mass_kg"]))))
            self.assertEqual(evaluate(p)["status"],"FAIL")

    def test_m2_and_storage_build_physics_are_distinct(self):
        root=Path(__file__).resolve().parents[4]
        spec=importlib.util.spec_from_file_location("m2_builder",root/"srcMoon/test/build_m2_qualification.py")
        m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
        production=root/"input/moon.input"
        for purpose,chem,weight in (("m2","off","on"),("sampling","on","off")):
            text=m.derive_compile_input(production,"moon_test_fixture",purpose)
            self.assertIn("PhotolyticReactions="+chem,text)
            self.assertIn("ParticleWeightCorrectionMode="+weight,text)
        self.assertIn("SpeciesList=Na",production.read_text())


if __name__=="__main__":unittest.main()
