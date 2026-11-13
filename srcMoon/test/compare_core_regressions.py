#!/usr/bin/env python3
"""Audit current-core regression preservation by exact suite/case identity.

Reports are ordered explicitly, oldest first. A retry supplies the final
disposition of its selected cases, while EVERY actual earlier baseline PASS
remains protected. Missing individual rows never become PASS. Existing
upstream FAIL/ERROR/SKIP outcomes stay visible; this audit's PASS means only
that integration preserved the passing regressions, not that all tests or
scientific milestones passed. Protocol self-checks run before every audit.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import tempfile


STATUSES={"PASS","FAIL","ERROR","SKIP"}


def read_attempt(path: Path):
    """Read registry-selected cases and retain actual row/error provenance.

    An empty ERROR/FAIL suite has no measured per-case result; give its selected
    cases explicit infrastructure ERROR records instead of borrowing a prior
    PASS. A preflight SKIP remains a non-executed SKIP. Neither fallback can
    enter the protected baseline PASS set. The original report is untouched.
    """
    report=json.loads(path.read_text())
    selection=json.loads((path.parent/"selection.json").read_text())
    selected={r["id"]:r for r in selection["tests"]}
    rows={};expected=set()
    for suite in report["results"]:
        sid=suite["id"]
        if sid not in selected:raise ValueError("unselected suite in report: "+sid)
        declaration=selected[sid]
        names=declaration.get("selected_cases",declaration["expected_cases"])
        if not names or len(names)!=len(set(names)):raise ValueError("invalid selected case registry")
        declared={(sid,name) for name in names};expected.update(declared)
        for case in suite.get("cases",[]):
            key=(sid,case["id"]);status=case["status"]
            if key not in declared or key in rows or status not in STATUSES:
                raise ValueError("duplicate/unknown case or invalid status: "+str(key))
            rows[key]={"status":status,"actual_case_row":True,
                "message":case.get("message",""),"report":str(path),
                "runner":suite.get("runner"),"log":suite.get("log")}
        for key in declared-rows.keys():
            skipped=suite.get("status")=="SKIP"
            rows[key]={"status":"SKIP" if skipped else "ERROR",
                "actual_case_row":False,"report":str(path),
                "message":"no individual result; suite: "+suite.get("message",""),
                "runner":suite.get("runner"),"log":suite.get("log")}
    # A suite omitted from the final report is also missing evidence, even if
    # an earlier attempt happened to pass. Use the selection as the authority.
    for sid,declaration in selected.items():
        for name in declaration.get("selected_cases",declaration["expected_cases"]):
            key=(sid,name);expected.add(key)
            if key not in rows:rows[key]={"status":"ERROR","actual_case_row":False,
                "report":str(path),"runner":declaration.get("runner_path"),
                "log":None,"message":"selected suite/case absent from report"}
    return expected,rows


def combine(paths):
    """Keep the last explicit attempt, plus all real baseline PASS identities."""
    expected=set();latest={};ever_pass=set();sources=[]
    for path in paths:
        names,rows=read_attempt(path);expected.update(names);latest.update(rows)
        ever_pass.update(k for k,v in rows.items()
                         if v["actual_case_row"] and v["status"]=="PASS")
        sources.append({"path":str(path),"sha256":hashlib.sha256(path.read_bytes()).hexdigest()})
    return expected,latest,ever_pass,sources


def evaluate(baseline_paths,candidate_paths):
    before_expected,before,protected,before_sources=combine(baseline_paths)
    after_expected,after,_,after_sources=combine(candidate_paths)
    regressions=[];changes=[];missing=sorted(before_expected-after_expected)
    for key in sorted(before_expected|after_expected):
        a=before.get(key);b=after.get(key)
        record={"suite":key[0],"case":key[1],"baseline":a,"candidate":b}
        if (a or {}).get("status")!=(b or {}).get("status"):changes.append(record)
        if key in protected and (b or {}).get("status")!="PASS":
            regressions.append(record)
    newly_failing=[r for r in changes if (r["candidate"] or {}).get("status") in ("FAIL","ERROR")
                   and (r["baseline"] or {}).get("status") not in ("FAIL","ERROR")]
    return {"schema":"srcMoon-current-core-regression-comparison-v1",
        "scope":"individual regression preservation; not all-green execution or observational validation",
        "status":"FAIL" if regressions or missing or newly_failing else "PASS",
        "baseline_reports":before_sources,"candidate_reports":after_sources,
        "baseline_counts":dict(Counter(r["status"] for r in before.values())),
        "candidate_counts":dict(Counter(r["status"] for r in after.values())),
        "protected_baseline_passes":len(protected),
        "preserved_baseline_passes":sum((after.get(k) or {}).get("status")=="PASS" for k in protected),
        "regressions":regressions,"missing_case_ids":missing,
        "newly_failing_cases":newly_failing,"status_changes":changes},before,after


def self_checks():
    """Test false-green hazards using protocol fixtures, never science data."""
    with tempfile.TemporaryDirectory(prefix="amps-regression-protocol-") as name:
        root=Path(name)
        def fixture(label,status,cases):
            directory=root/label;directory.mkdir()
            (directory/"selection.json").write_text(json.dumps({"tests":[{
                "id":"SUITE","expected_cases":["A","B"],"selected_cases":["A","B"]}]}))
            path=directory/"report.json";path.write_text(json.dumps({"results":[{
                "id":"SUITE","status":status,"cases":[{"id":k,"status":v} for k,v in cases]}]}))
            return path
        good=fixture("good","SKIP",[("A","PASS"),("B","SKIP")])
        broken=fixture("broken","ERROR",[])
        skipped=fixture("skipped","SKIP",[("A","SKIP"),("B","SKIP")])
        result,_,_=evaluate([good],[good])
        assert result["status"]=="PASS" and result["protected_baseline_passes"]==1
        result,_,_=evaluate([good],[good,broken])
        assert result["status"]=="FAIL" and len(result["regressions"])==1
        result,_,_=evaluate([good,broken],[skipped])
        assert result["status"]=="FAIL" and result["protected_baseline_passes"]==1
    return 3


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline",type=Path,action="append",required=True,
                        help="repeat in explicit chronological/retry order")
    parser.add_argument("--candidate",type=Path,action="append",required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args();out=args.output_dir.resolve()
    if out.exists() and any(out.iterdir()):parser.error("output directory must be new or empty")
    out.mkdir(parents=True,exist_ok=True)
    try:
        checked=self_checks()
        report,before,after=evaluate([p.resolve() for p in args.baseline],
                                    [p.resolve() for p in args.candidate])
        report["infrastructure_protocol_checks_passed"]=checked
        with (out/"comparison.csv").open("w",newline="") as stream:
            writer=csv.writer(stream);writer.writerow(("suite","case","baseline","candidate","candidate_runner","candidate_log"))
            for key in sorted(before.keys()|after.keys()):
                a=before.get(key,{});b=after.get(key,{})
                writer.writerow((*key,a.get("status","MISSING"),b.get("status","MISSING"),b.get("runner"),b.get("log")))
    except (OSError,ValueError,KeyError,TypeError,AssertionError) as error:
        report={"status":"ERROR","scope":"regression audit infrastructure","reason":str(error)}
    (out/"result.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps(report,indent=2,sort_keys=True))
    return {"PASS":0,"FAIL":1,"ERROR":2}[report["status"]]


if __name__=="__main__":raise SystemExit(main())
