#!/usr/bin/env python3
"""Exercise all catalog receipt modes through the actual C++ AMPS reader.

The caller supplies its configured dependency-light test binary, so no second
compiler/toolchain is guessed here. Synthetic nonuniform rates/geometries
test serialization and membership only; they are not observational fixtures.
"""
import argparse
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from standalone.tests.test_fit_source_export import fixture_config,fixture_rates
from fit_source_export import export_vectors
from source_contract import CASES,HYBRID,source_ids


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable",type=Path,required=True)
    args=parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="enceladus-source-roundtrip-") as tmp:
        for case in CASES:
            if case==HYBRID:
                continue  # Existing U03 checks every hybrid geometry reference.
            out=Path(tmp)/case
            index=export_vectors(fixture_config(case),[("test","E3","H2O",fixture_rates(case),"H2O")],out,{"fixture":"synthetic"})
            for files in index["artifacts"][0]["files"].values():
                path=out/files["source_file"]
                result=subprocess.run([str(args.executable.resolve()),"catalog",str(path),str(len(source_ids(case)))],
                                      text=True,capture_output=True)
                if result.returncode:
                    raise RuntimeError(result.stdout+result.stderr)
                # Finite subnormal strengths are legitimate fitted tails.
                # Replacing an exact zero changes neither a rounded shape
                # sum nor an absolute total; it must still preserve the tiny
                # nonzero input, rather than interpreting ERANGE as overflow.
                tiny=out/(path.stem+".subnormal.in")
                tiny.write_text(path.read_text().replace("strength.H2O=0\n",
                    "strength.H2O=3.7637234148938442e-317\n",1))
                accepted=subprocess.run([str(args.executable.resolve()),"catalog_subnormal",str(tiny),str(len(source_ids(case)))],
                                        text=True,capture_output=True)
                if accepted.returncode:
                    raise AssertionError("representable subnormal rejected: "+accepted.stderr)
                for invalid in ("1e-999","1e999"):
                    tiny.write_text(path.read_text().replace("strength.H2O=0\n",
                        "strength.H2O="+invalid+"\n",1))
                    rejected=subprocess.run([str(args.executable.resolve()),"catalog",str(tiny),str(len(source_ids(case)))],
                                            text=True,capture_output=True)
                    if rejected.returncode==0:
                        raise AssertionError("range-losing source strength accepted: "+invalid)
                # A valid field block with a noncanonical ID must still fail:
                # counting 41/98/99 syntactically valid points is insufficient.
                bad=out/(path.stem+".bad.in")
                bad.write_text(path.read_text().replace("#begin source: porco_001", "#begin source: unknown_001",1))
                rejected=subprocess.run([str(args.executable.resolve()),"catalog",str(bad),str(len(source_ids(case)))],
                                        text=True,capture_output=True)
                if rejected.returncode==0:
                    raise AssertionError("noncanonical catalog source was accepted")
    print("PASS: 98/41/99 shape and absolute receipts; malformed membership rejected")


if __name__=="__main__":
    main()
