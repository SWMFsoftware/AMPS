#!/usr/bin/env python3
"""REL3D01--03 protocol verification; production qualification stays separate.

Manufactured executables below prove argv/report ownership only. They are
explicitly synthetic and cannot establish MPI or observational qualification.
"""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"tools"))
from preprocessing.core import freeze_record, PreprocessingError
from release.qualification import (source_manifest, verify_sources, package_sources,
    verify_archive, qualify, run_applications, update_last_pass, file_digest)


def profile():
    return freeze_record(dict(schema="sccm-release-profile-v1", software_stage=13,
        stage14_is_release_dependency=False, required_tests=["REAL01"], event_profile=False))


def incomplete():
    return freeze_record(dict(test_results=[dict(id="REAL01", status="PASS")],
        evidence_kind="synthetic-protocol-verification", evidence_classes=["verification"]))


class REL3D01(unittest.TestCase):
    def test_manifest_ownership_hashes_and_no_false_qualification(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in ("src/kernel.cpp", "test/check.py", "examples/deck.in"):
                path = root/name; path.parent.mkdir(exist_ok=True); path.write_text("// fixture\n")
            registry = [dict(id="REAL01", owner="shared-model", full_suite_command=["runner", "--all"], implementation="test/check.py")]
            manifest = source_manifest(root, ["src/kernel.cpp", "test/check.py", "examples/deck.in"], registry,
                                       ["src/kernel.cpp"], ["examples/deck.in"])
            verify_sources(root, manifest)
            (root/"src/kernel.cpp").write_text("// mutation\n")
            with self.assertRaises(ValueError): verify_sources(root, manifest)
            with self.assertRaises(ValueError): source_manifest(root, ["test/check.py"], registry, ["missing.cpp"], [])
            result = qualify(profile(), incomplete())
            self.assertFalse(result["qualified"])
            self.assertEqual(result["status"], "INCOMPLETE")
            self.assertTrue(any("D10" in row for row in result["blockers"]))
            self.assertTrue(any("runtime" in row for row in result["blockers"]))
            with self.assertRaises(ValueError): update_last_pass(root/"last-pass.json", result, "0"*64)
            for state in ("FAIL", "SKIP", "ERROR"):
                data = incomplete(); data.pop("identity"); data["test_results"][0]["status"] = state
                self.assertTrue(any(state in row for row in qualify(profile(), freeze_record(data))["blockers"]))
        # Every baseline test remains selected by the owning full shared suite.
        listed = subprocess.run([sys.executable, str(ROOT/"test/run_tests.py"), "--list"], text=True, stdout=subprocess.PIPE, check=True).stdout
        for name in ("REL3D01", "REL3D02", "REL3D03", "ARCHSCCM01", "CAL3D01"):
            self.assertIn(name+"\t", listed)


class REL3D02(unittest.TestCase):
    def test_explicit_independent_applications_bundle_and_stale_report(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); bundle = root/"bundle.json"; bundle.write_text("synthetic immutable bundle\n")
            programs, commands = {}, {}
            for application in ("srcSEP", "srcSEP3D"):
                script = root/(application+".py")
                script.write_text('import json,hashlib,sys\nfrom pathlib import Path\n'
                    'data={"application":'+repr(application)+',"status":"PASS","evidence_kind":"synthetic-protocol-verification",'
                    '"bundle_sha256":hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest()}\n'
                    'Path(sys.argv[2]).write_text(json.dumps(data))\n')
                programs[application] = script
                commands[application] = [sys.executable, "{executable}", "{bundle}", "{report}"]
            result = run_applications(programs["srcSEP"], programs["srcSEP3D"], bundle, commands, root/"fresh")
            self.assertEqual([r["status"] for r in result["results"]], ["PASS", "PASS"])
            self.assertTrue(all(Path(r["log"]).is_file() for r in result["results"]))
            with self.assertRaises(ValueError): run_applications(programs["srcSEP"], programs["srcSEP"], bundle, commands, root/"duplicate")
            with self.assertRaises(FileExistsError): run_applications(programs["srcSEP"], programs["srcSEP3D"], bundle, commands, root/"fresh")
            programs["srcSEP"].write_text("raise SystemExit(0)\n")
            result = run_applications(programs["srcSEP"], programs["srcSEP3D"], bundle, commands, root/"later")
            self.assertEqual(result["results"][0]["status"], "ERROR")
            self.assertEqual(result["results"][1]["status"], "PASS")
            programs["srcSEP"].unlink()
            result=run_applications(programs["srcSEP"],programs["srcSEP3D"],bundle,commands,root/"missing-executable")
            self.assertEqual([r["status"] for r in result["results"]],["ERROR","PASS"])


class REL3D03(unittest.TestCase):
    def test_package_exact_members_hashes_and_generated_artifact_rejection(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); (root/"source.cpp").write_text("// synthetic source\n")
            rows = [dict(id="REAL01", owner="shared-model", full_suite_command=["runner", "--all"], implementation="source.cpp")]
            manifest = source_manifest(root, ["source.cpp"], rows, ["source.cpp"], [])
            archive = root/"source.tar.gz"; package_sources(root, manifest, archive); verify_archive(archive, manifest)
            embedded=root/"embedded.tar.gz"
            package_sources(root,manifest,embedded,"RELEASE_SOURCE_MANIFEST.json")
            verify_archive(embedded,manifest,"RELEASE_SOURCE_MANIFEST.json")
            for name in ("build/test.o", "test_output/results.json", "compiled.a", "../escape.cpp"):
                with self.assertRaises(ValueError): source_manifest(root, [name], rows, [], [])
            changed = copy.deepcopy(manifest); changed.pop("identity"); changed["files"]["source.cpp"] = "0"*64
            with self.assertRaises(ValueError): verify_archive(archive, freeze_record(changed))
            with tarfile.open(root/"extra.tar", "w") as tar:
                tar.add(root/"source.cpp", arcname="source.cpp"); tar.add(root/"source.cpp", arcname="extra.cpp")
            with self.assertRaises(ValueError): verify_archive(root/"extra.tar", manifest)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(); parser.add_argument("--test", choices=["REL3D01", "REL3D02", "REL3D03"])
    args = parser.parse_args()
    cases = [globals()[args.test]] if args.test else [REL3D01, REL3D02, REL3D03]
    suite = unittest.TestSuite(unittest.defaultTestLoader.loadTestsFromTestCase(case) for case in cases)
    raise SystemExit(0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1)
