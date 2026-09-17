#!/usr/bin/env python3
"""Materialize the immutable VP15 benchmark descriptor.

This case uses an analytical/literature/coupling fixture rather than a remote
bulk data product.  The acquisition command still creates a checksummed raw
input so offline runs and campaign manifests have the same provenance boundary
as observational cases with network downloads.
"""
from __future__ import annotations
import hashlib,json
from pathlib import Path
CASE="VP15"
PAYLOAD={"schema_version":1,"validation_id":CASE,"title":"SEP onset, profile, spectrum, peak intensity, and fluence","fixture_kind":"COUPLING_BENCHMARK","source":"case-owned reference_solution.py and public SWCME production API","network_required":False}
def main() -> int:
    target=Path(__file__).resolve().parent/"data/raw/benchmark.json"
    target.parent.mkdir(parents=True,exist_ok=True)
    data=(json.dumps(PAYLOAD,indent=2,sort_keys=True)+"\\n").encode()
    target.write_bytes(data)
    print(f"PASS {hashlib.sha256(data).hexdigest()}  {target}")
    return 0
if __name__=="__main__": raise SystemExit(main())

