#!/usr/bin/env python3
"""Repeat the manuscript H2O/CO2 matrix without overwriting an old campaign.

Independent geometries run in separate processes with one BLAS thread each.
This is fit-process concurrency, not concurrent AMPS configuration/builds.
Every subprocess has its own immutable command, elapsed time and log receipt.
No model selection, calibration factor or source-rate averaging is implicit.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from source_contract import HYBRID, PORCO, SOLAR, CANDIDATE

ROOT = Path(__file__).resolve().parent
BASES = (HYBRID, PORCO, SOLAR)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation-data-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--paper", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=3)
    parser.add_argument("--resume", action="store_true",
                        help="retain successful steps; repair a completed fit's failed export without refitting")
    parser.add_argument("--resume-solver-repair", action="store_true",
                        help="record a fitting.py solver-only change and retry an incomplete failed fit, preserving its old directory/log")
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    destination = args.output_dir.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    # Existing provenance can have been prepared before execution. A campaign
    # receipt or scientific subdirectory, however, forbids an accidental rerun.
    for name in ("campaign.json", "execution", "h2o", "co2", "e17_search"):
        if (destination/name).exists() and not args.resume:
            parser.error("refusing to overwrite campaign: "+str(destination/name))
    execution = destination/"execution"
    execution.mkdir(exist_ok=args.resume)
    provenance = destination/"provenance"
    provenance.mkdir(exist_ok=True)
    if not args.resume:
        shutil.copy2(args.paper, provenance/"manuscript_requested.pdf")
    stamp=datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    snapshot = provenance/("standalone_software_resume_"+stamp if args.resume else "standalone_software")
    snapshot.mkdir()
    source_hashes = {}
    for path in sorted(ROOT.glob("*.py")):
        shutil.copy2(path, snapshot/path.name)
        source_hashes[path.name] = digest(path)
    shutil.copytree(ROOT/"config", snapshot/"config")
    config = ROOT/"config/tenishev2014_stage0.json"
    environment = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
                       MKL_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")

    def run(label, script, options):
        path = execution/(label+".json")
        if args.resume and path.is_file():
            previous=json.loads(path.read_text())
            if previous["status"] == "PASS":
                print("RETAIN "+label,flush=True)
                return
            # A fit may have completed all numerical products and failed only
            # in source serialization. Preserve both the failed log and the
            # fitted JSON. Repair only that handoff; never rerun/alter the
            # scientific objective merely to hide a transport-format error.
            if script.name in ("joint_fitting_study.py","co2_fitting_study.py"):
                study=Path(options[options.index("--output-dir")+1])
                fitted=study/"study_results.json"
                repair_path=execution/(label+"_source_export_repair.json")
                if repair_path.is_file():
                    print("RETAIN repaired "+label,flush=True)
                    return
                if fitted.is_file():
                    document=json.loads(fitted.read_text())
                    source=study/"amps_sources"
                    if source.exists():
                        source.rename(study/("amps_sources_failed_"+stamp))
                    from fit_source_export import export_study
                    species="H2O" if script.name.startswith("joint_") else "CO2"
                    export_study(document,document["resolved_geometry"],study,species)
                    repair_path.write_text(json.dumps({"status":"PASS", "operation":
                        "export completed fit without refitting", "fitted_result_sha256":digest(fitted),
                        "exporter_sha256":digest(ROOT/"fit_source_export.py"),
                        "original_failed_receipt":path.name},indent=2)+"\n")
                    print("PASS export repair "+label,flush=True)
                    return
            if args.resume_solver_repair and script.name in ("joint_fitting_study.py","co2_fitting_study.py"):
                if study.exists(): study.rename(study.with_name(study.name+"_failed_"+stamp))
                path.rename(execution/(label+"_failed_"+stamp+".json"))
                (execution/(label+".log")).rename(execution/(label+"_failed_"+stamp+".log"))
            else:
                raise RuntimeError("failed step has no complete fit to resume: "+label)
        command = [sys.executable, "-B", str(script), *map(str, options)]
        receipt = {"label": label, "command": command,
            "utc_start": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "log": label+".log", "status": "RUNNING"}
        path.write_text(json.dumps(receipt, indent=2)+"\n")
        print("START "+label, flush=True)
        start = time.monotonic()
        with (execution/(label+".log")).open("x") as stream:
            result = subprocess.run(command, env=environment, stdout=stream,
                                    stderr=subprocess.STDOUT, cwd=ROOT)
        receipt.update(seconds=time.monotonic()-start, returncode=result.returncode,
                       status="PASS" if result.returncode == 0 else "FAIL")
        path.write_text(json.dumps(receipt, indent=2)+"\n")
        print(receipt["status"]+" "+label+" {:.1f}s".format(receipt["seconds"]), flush=True)
        if result.returncode:
            raise RuntimeError("failed "+label+"; see "+str(execution/(label+".log")))

    def fit_pair(case, extra=None):
        base = PORCO if case == CANDIDATE else case
        common = ["--config", config, "--source-case", base,
                  "--validation-data-root", args.validation_data_root,
                  "--stripe-surface-mapping", "local_vertical_to_reference_sphere",
                  "--pooled-backgrounds", "per_observation"]
        if extra:
            common += ["--extra-jet", extra]
        h2o = destination/"h2o"/case
        co2 = destination/"co2"/case
        run("h2o_"+case, ROOT/"joint_fitting_study.py",
            [*common, "--quadrature-order", "192", "--output-dir", h2o])
        run("co2_"+case, ROOT/"co2_fitting_study.py",
            [*common, "--h2o-study-result", h2o/"study_results.json", "--output-dir", co2])

    campaign = {"status": "RUNNING", "paper_sha256": digest(args.paper),
        "software_sha256": source_hashes, "base_cases": BASES,
        "candidate_case": CANDIDATE, "expected_target_fits": 124,
        "decisions": {"stripe_mapping": "radial projection to 252100 m before fitting",
            "pooled_backgrounds": "per observation",
            "absolute_E3_E5_CO2": "DEFERRED until physical fragmentation correction selected"}}
    receipt_path = destination/"campaign.json"
    if args.resume:
        previous=json.loads(receipt_path.read_text())
        # Recovery permits exporter/report corrections, but refuses an
        # unrecorded change to the fitted physical/numerical operators.
        for name in ("fitting.py","plume_model.py","joint_fitting_study.py","co2_fitting_study.py"):
            if (previous["software_sha256"][name] != source_hashes[name]
                    and not (name=="fitting.py" and args.resume_solver_repair)):
                raise ValueError("scientific fit code changed; use a new campaign: "+name)
        previous["resume_software_sha256"]=source_hashes
        if args.resume_solver_repair:
            previous.setdefault("solver_repair_history",[]).append({"utc":stamp,
                "original_fitting_sha256":previous["software_sha256"]["fitting.py"],
                "current_fitting_sha256":source_hashes["fitting.py"],
                "scope":"same objective, active-set NNLS only after TRF failure; KKT <=1e-12"})
        campaign=previous
        campaign["status"]="RUNNING"
    receipt_path.write_text(json.dumps(campaign, indent=2)+"\n")
    try:
        with ThreadPoolExecutor(max_workers=min(3, args.jobs)) as executor:
            list(executor.map(fit_pair, BASES))
        run("h2o_geometry_comparison", ROOT/"compare_source_case_studies.py",
            ["--study-root", destination/"h2o"])
        search = destination/"e17_search"
        run("e17_search", ROOT/"extra_jet_e17_study.py",
            ["--config", config, "--validation-data-root", args.validation_data_root,
             "--co2-baseline-result", destination/"co2"/PORCO/"study_results.json",
             "--output-dir", search])
        selected = json.loads((search/"study_results.json").read_text())
        candidate = json.loads((ROOT/"e17_extra_jet_candidate.json").read_text())
        candidate["source"].update(selected["primary_candidate"])
        nearest=selected["nearest_catalog_jets"][0]
        candidate["selection"]={"target":"E17 mass-44/CO2 closest-approach feature",
            "time_window_s":[-selected["search_target"]["half_window_s"],selected["search_target"]["half_window_s"]],
            "constraint":"radial source sampled on modeled tiger-stripe polylines",
            "nearest_porco_source":nearest["id"],
            "nearest_porco_angular_distance_deg":nearest["angular_distance_deg"]}
        candidate["provenance"] = {
            "campaign_search_result": str(search/"study_results.json"),
            "study_result_sha256": digest(search/"study_results.json"),
            "interpretation": "Selected again on E17; not an independent vent detection"}
        candidate_path = destination/"selected_e17_candidate.json"
        candidate_path.write_text(json.dumps(candidate, indent=2)+"\n")
        fit_pair(CANDIDATE, candidate_path)
        run("candidate_comparison", ROOT/"compare_extra_jet_studies.py", [
            "--base-h2o-uvis", destination/"h2o"/PORCO/"study_results.json",
            "--augmented-h2o-uvis", destination/"h2o"/CANDIDATE/"study_results.json",
            "--base-co2", destination/"co2"/PORCO/"study_results.json",
            "--augmented-co2", destination/"co2"/CANDIDATE/"study_results.json",
            "--output-dir", destination/"candidate_comparison"])
        # These retained observation-side tools implement the manuscript's
        # phase proxies. They do not calculate geological stress or establish
        # a significant four-epoch correlation; their reports retain that limit.
        scripts = args.validation_data_root/"scripts"
        run("h2o_orbital_proxy", scripts/"plot_tidal_proxy_correlation.py", [
            "--fit-study-root", destination/"h2o", "--output-dir", destination/"h2o_orbital_proxy"])
        run("co2_orbital_proxy", scripts/"plot_co2_h2o_tidal_study.py", [
            "--h2o-study-root", destination/"h2o", "--co2-study-root", destination/"co2",
            "--output-dir", destination/"co2_orbital_proxy"])
        campaign["status"] = "FITS_COMPLETE_REPORT_PENDING"
    except BaseException as error:
        campaign.update(status="FAILED", error=str(error))
        raise
    finally:
        receipt_path.write_text(json.dumps(campaign, indent=2)+"\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
