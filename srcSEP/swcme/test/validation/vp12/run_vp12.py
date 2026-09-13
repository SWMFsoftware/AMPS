#!/usr/bin/env python3
"""Run VP12: SEP source records, identity, attribution, and source spectrum."""
from __future__ import annotations
import argparse
from pathlib import Path
import subprocess
import sys

CASE_DIR=Path(__file__).resolve().parent
sys.path.insert(0,str(CASE_DIR.parent))
from case_runner import run

def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download",action="store_true",help="materialize and verify the pinned benchmark manifest")
    parser.add_argument("--no-plots",action="store_true")
    parser.add_argument("--output-dir",type=Path,default=CASE_DIR/"output")
    args=parser.parse_args()
    if args.download:
        subprocess.run([sys.executable,str(CASE_DIR/"download_data.py")],check=True)
    return run("VP12",CASE_DIR,args.output_dir,args.no_plots)

if __name__=="__main__":
    try: raise SystemExit(main())
    except Exception as error:
        print(f"VP12 setup/analysis error: {error}",file=sys.stderr)
        raise SystemExit(2)

