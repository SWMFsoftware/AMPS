"""Publication contracts use exact synthetic source vectors, not fitted truth."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"test"))
from case_fixtures import make_fixture
from prepare_forward_cases import prepare_cases,read_source


class ForwardCaseExportTests(unittest.TestCase):
    def test_all_measurements_and_bases_are_explicit_and_complete(self):
        with tempfile.TemporaryDirectory() as t:
            deck,data,selection,pub=make_fixture(Path(t))
            result=json.loads((pub/"publication.json").read_text())
            self.assertEqual(len(result["cases"]),28)
            for case,entry in result["cases"].items():
                block,species,_=read_source(pub/"cases"/case/"sources.in")
                self.assertEqual(block["fields"]["source_case"],case.split("--")[1])
                self.assertEqual(species,set(entry["species"]))
            with self.assertRaises(FileExistsError):prepare_cases(selection,pub)

    def test_matching_species_merge_without_background_or_rate_changes(self):
        with tempfile.TemporaryDirectory() as t:
            _,_,_,pub=make_fixture(Path(t),all_cases=False,species_pair=True,mode="shape")
            block,species,data=read_source(pub/"cases/E3--tenishev2014_hybrid_8/sources.in")
            self.assertEqual(species,{"H2O","CO2"})
            self.assertNotIn(b"background:",data)
            self.assertEqual(block["children"][0]["fields"]["strength.H2O"],"0")

    def test_incompatible_kinematics_and_duplicate_species_are_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);_,_,selection,_=make_fixture(root,all_cases=False)
            d=json.loads(selection.read_text());case=next(iter(d["cases"]))
            d["cases"][case]["source_files"]*=2
            selection.write_text(json.dumps(d))
            with self.assertRaisesRegex(ValueError,"same species"):prepare_cases(selection,root/"duplicate")
            self.assertFalse((root/"duplicate").exists())
            h2o=Path(d["cases"][case]["source_files"][0]);co2=next((root/"raw").glob("*CO2/*absolute.in"))
            co2.write_text(co2.read_text().replace("temperature_K=180","temperature_K=181",1))
            d["cases"][case]["source_files"]=[str(h2o),str(co2)]
            d["cases"][case]["physical_CO2_calibration"]="synthetic fixture only"
            selection.write_text(json.dumps(d))
            with self.assertRaisesRegex(ValueError,"different source geometry"):prepare_cases(selection,root/"kinematics")

    def test_no_physical_e3_co2_correction_is_invented(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);_,_,selection,_=make_fixture(root,all_cases=False)
            d=json.loads(selection.read_text());case=next(iter(d["cases"]))
            d["cases"][case]["source_files"]=[str(next((root/"raw").glob("*CO2/*absolute.in")))]
            selection.write_text(json.dumps(d))
            with self.assertRaisesRegex(ValueError,"deferred"):prepare_cases(selection,root/"physical")

if __name__=="__main__":unittest.main()
