#!/usr/bin/env python3
"""Relocate generated upstream regression-input provenance without physics edits.

The Enceladus exporter records absolute source paths in three comment lines.
Its existing reproducibility test compares the complete generated receipt, so
a fresh checkout must reproduce those path comments before running the test.
This preparation invokes the committed exporter, permits ONLY those three
path-comment differences, and preserves every other byte (including hashes,
fit provenance, geometry, source strengths and physical controls). It is an
infrastructure operation, never observation processing or validation evidence.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


PATH_COMMENTS=(b"! manifest=",b"! geometry_config=",b"! stripe_geometry=")
PAIRS=(("e3_equal_h2o_co2_startup_manifest.json","enceladus_e3_equal_h2o_co2_sources.in"),
       ("e5_equal_h2o_co2_startup_manifest.json","enceladus_e5_equal_h2o_co2_sources.in"))


def sha(path: Path) -> str:
    """Bind preparation evidence to exact existing source/receipt bytes."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def payload(data: bytes) -> bytes:
    """Remove exactly the exporter-defined host-path comments, nothing else.

    Require one occurrence of each key. Missing/duplicate metadata is corrupt
    input, not permission to strip arbitrary comments or adjust science rows.
    Hash comments remain part of the compared payload and cannot be changed.
    """
    lines=data.splitlines(keepends=True)
    for key in PATH_COMMENTS:
        if sum(line.startswith(key) for line in lines)!=1:
            raise ValueError("expected exactly one path comment: "+key.decode())
    return b"".join(line for line in lines if not line.startswith(PATH_COMMENTS))


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True,
                        help="new evidence directory; previous receipts are never overwritten")
    args=parser.parse_args();repo=args.repo_root.resolve();out=args.output_dir.resolve()
    if out.exists() and any(out.iterdir()):parser.error("output directory must be new or empty")
    if not (repo/"Config.pl").is_file():parser.error("repo-root is not an AMPS checkout")
    out.mkdir(parents=True,exist_ok=True)
    exporter=repo/"srcEnceladus/standalone/export_forward_sources.py"
    report={"schema":"srcMoon-core-regression-input-preparation-v1","status":"ERROR",
        "scope":"host-path metadata only; not a physics/observation test",
        "repo_root":str(repo),"exporter_sha256":sha(exporter),"files":[]}
    staged=[]
    try:
        # Prepare and verify BOTH receipts before replacing either template.
        # The real exporter uses exclusive creation in this new evidence root;
        # neither an earlier campaign receipt nor raw external data is edited.
        for manifest_name,input_name in PAIRS:
            manifest=repo/"srcEnceladus/standalone/results"/manifest_name
            source=repo/"input"/input_name;generated=out/input_name
            command=[sys.executable,str(exporter),"--manifest",str(manifest),
                     "--output",str(generated)]
            done=subprocess.run(command,cwd=repo,text=True,capture_output=True)
            (out/(input_name+".log")).write_text(done.stdout+done.stderr)
            if done.returncode:raise ValueError("exporter failed for "+input_name)
            if payload(source.read_bytes())!=payload(generated.read_bytes()):
                raise ValueError("refusing non-path/science change in "+input_name)
            previous=out/(input_name+".before");shutil.copy2(source,previous)
            report["files"].append({"path":str(source),"command":command,
                "manifest_sha256":sha(manifest),"before_sha256":sha(previous),
                "after_sha256":sha(generated),"non_path_payload_sha256":
                    hashlib.sha256(payload(generated.read_bytes())).hexdigest()})
            staged.append((generated,source))
        for generated,source in staged:shutil.copy2(generated,source)
        report["status"]="PASS"
    except (OSError,ValueError) as error:
        report["reason"]=str(error)
    (out/"result.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps(report,indent=2,sort_keys=True))
    return 0 if report["status"]=="PASS" else 2


if __name__=="__main__":raise SystemExit(main())
