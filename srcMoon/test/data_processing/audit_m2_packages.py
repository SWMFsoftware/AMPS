#!/usr/bin/env python3
"""Read-only D01/D02/D03 readiness and hash audit for M2 development.

Native metadata and existing processor receipts are evidence, not a reason to
invent package-level READY status. This audit never changes canonical/raw
files, reads model residuals, or exposes hold-out science rows. It hashes the
already frozen I17 history as bytes and verifies its source/processor chain.
"""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from integration_common import dataset_package_preflight,sha256,validation_data_root


def verify_file(path,expected):
    """A missing required file is unavailable; altered bytes are corruption."""
    if not path.is_file():raise FileNotFoundError(str(path))
    actual=sha256(path)
    if actual!=expected:raise ValueError(f"hash mismatch: {path}; expected {expected}, got {actual}")
    return {"path":str(path),"sha256":actual,"bytes":path.stat().st_size}


def inspect():
    root=validation_data_root();records=[]
    for dataset in ("D01","D02","D03"):
        record={"check":dataset+"-package","dataset":dataset,"scope":"data-readiness-only"}
        try:
            ready,reason,evidence=dataset_package_preflight(dataset)
            record.update(status="PASS" if ready else "SKIPPED",reason=reason,evidence=evidence)
            package=root/("diviner" if dataset=="D03" else "lola")
            provenance=package/"provenance.json"
            verified=[]
            if provenance.is_file():
                for source in json.loads(provenance.read_text()).get("raw_files",[]):
                    p=Path(source["path"])
                    if not p.is_absolute():p=package/p
                    verified.append(verify_file(p,source["sha256"]))
            record["verified_raw_files"]=verified
        except FileNotFoundError as error:record.update(status="SKIPPED",reason=f"required data unavailable: {error}")
        except (OSError,ValueError,KeyError) as error:record.update(status="ERROR",reason=str(error))
        records.append(record)
    qa_path=root/"diviner/processed/i17/d03_rdr_qa.json"
    record={"check":"D03-I17-integrity","dataset":"D03","scope":"byte-integrity-only-not-validation"}
    try:
        if not qa_path.is_file():raise FileNotFoundError(str(qa_path))
        qa=json.loads(qa_path.read_text())
        if qa["status"]!="READY":raise FileNotFoundError("I17 processor QA is not READY")
        chain=[verify_file(Path(qa["processed"]["path"]),qa["processed"]["sha256"]),
               verify_file(Path(qa["processor"]["path"]),qa["processor"]["sha256"]),
               verify_file(Path(qa["contract"]["snapshot_path"]),qa["contract"]["sha256"])]
        for raw in qa["raw_files"]:
            chain.append(verify_file(Path(raw["path"]),raw["sha256"]))
        record.update(status="PASS",verified_files=chain,
            reason="frozen processed/processor/selection/raw bytes match; package qualification and observation operators are separate")
    except FileNotFoundError as error:record.update(status="SKIPPED",reason=str(error))
    except (OSError,ValueError,KeyError) as error:record.update(status="ERROR",reason=str(error))
    records.append(record)
    counts={s:sum(r["status"]==s for r in records) for s in ("PASS","FAIL","SKIPPED","ERROR")}
    return {"schema":"srcMoon-M2-package-audit-v1","data_root":str(root),
            "raw_modified":False,"observational_scoring":False,"holdout_values_read":False,
            "counts":counts,"records":records}


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir",required=True,type=Path);a=p.parse_args()
    if a.output_dir.exists() and any(a.output_dir.iterdir()):p.error("output directory must be new or empty")
    a.output_dir.mkdir(parents=True,exist_ok=True)
    result=inspect()
    (a.output_dir/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"counts":result["counts"],"result":str(a.output_dir/"result.json")},indent=2))
    raise SystemExit(2 if result["counts"]["ERROR"] else 0)
