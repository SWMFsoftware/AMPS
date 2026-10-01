#!/usr/bin/env python3
"""Stage-13 release CLI; missing evidence is INCOMPLETE, never a release PASS."""
import argparse
import json
from pathlib import Path
from release.qualification import qualify, run_applications, update_last_pass
from preprocessing.core import load_json
from preprocessing.core import freeze_record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--srcsep-executable", type=Path)
    parser.add_argument("--srcsep3d-executable", type=Path)
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--application-output", type=Path)
    parser.add_argument("--last-pass", type=Path)
    args = parser.parse_args()
    profile, evidence = load_json(args.profile), load_json(args.evidence)
    selected = (args.srcsep_executable, args.srcsep3d_executable, args.bundle, args.application_output)
    run = None
    if any(selected):
        if not all(selected): parser.error("cross-application run requires both executables, bundle and fresh output directory")
        run = run_applications(*selected[:3], profile["application_commands"], selected[3])
        (selected[3]/"cross-application.json").write_text(json.dumps(run, indent=2)+"\n")
        # Running two applications is recorded separately. It never fabricates
        # missing physics/MPI/diagnostics into the caller's frozen evidence.
    result = qualify(profile, evidence)
    if run is not None:
        result.pop("identity")
        result["cross_application_run_identity"]=run["identity"]
        for row in run["results"]:
            if row["status"] != "PASS" or row.get("evidence_kind") != "production":
                result["blockers"].append("fresh cross-application "+row["application"]+" is "+row["status"]+" / "+row.get("evidence_kind","unqualified"))
        result["qualified"]=not result["blockers"]
        result["status"]="PASS" if result["qualified"] else "INCOMPLETE"
        result["observational_validation"]=bool(profile.get("event_profile") and result["qualified"])
        result=freeze_record(result)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+"\n")
    if args.last_pass:
        update_last_pass(args.last_pass, result, profile["reference_sha256"])
    print("release_qualification="+result["status"])
    for blocker in result["blockers"]: print("  "+blocker)
    return 0 if result["qualified"] else 1

if __name__ == "__main__":
    raise SystemExit(main())
