#!/usr/bin/env python3
"""Assemble formal I19 from the existing, complete native numerical campaign.

This module contains no source, residence, force, or chemistry law. It joins
the already approved cluster executor and independent accounting/equivalence
oracles to the ordinary per-ID runner. A copied ``PASS`` summary is never
enough: all 288 selected native receipts, controls and bytes are checked again
and the independent numerical scorer is rerun. Missing scientific evidence is
SKIPPED; altered/contradictory bytes are ERROR; an evaluated miss is FAIL.

Long execution is a separate explicit opt-in. Development invocations and
ordinary calls use ``allow_execute=False`` even if an input names a reviewed
execution contract. Immutable raw/processed observations are never modified.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
PROTOCOL = Path(__file__).with_name("M2_I19_I20_PROTOCOL.json")
SCHEMA = "srcMoon-I19-formal-execution-v1"


class Unavailable(RuntimeError):
    """A real predecessor/contract/native artifact has not been supplied."""


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json(path):
    """Reject duplicate keys and non-finite JSON constants in frozen proofs."""
    def pairs(items):
        value = {}
        for key, item in items:
            if key in value:
                raise ValueError("duplicate frozen-contract key: " + key)
            value[key] = item
        return value
    def invalid(value):
        raise ValueError("non-finite JSON value: " + value)
    return json.loads(Path(path).read_text(), object_pairs_hook=pairs,
                      parse_constant=invalid)


def _write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True,
                                   allow_nan=False) + "\n")


def _file(record, base, label):
    """Bind a named file to reviewed bytes; relative paths use contract owner.

    There are no fallback locations and no unverified path-only records.
    Absence is different from corruption so operators can repair missing
    dependencies without an infrastructure fault being mistaken for physics.
    """
    if not isinstance(record, dict) or set(record) != {"path", "sha256"}:
        raise ValueError(label + " requires exactly path and sha256")
    path = Path(record["path"])
    path = (path if path.is_absolute() else base / path).resolve()
    expected = record["sha256"]
    if not isinstance(expected, str) or len(expected) != 64 or any(
            c not in "0123456789abcdef" for c in expected):
        raise ValueError(label + " has an invalid SHA-256")
    if not path.is_file():
        raise Unavailable(label + " unavailable: " + str(path))
    if _sha(path) != expected:
        raise ValueError(label + " SHA-256 mismatch: " + str(path))
    return path


def load_contract(path, executable, application_input, campaign_root, mpi_ranks):
    """Read explicit reviewed dispatch policy without opening model results.

    Review values are supplied by an investigator; this code never creates an
    approval. The canonical native-source signature includes dirty/untracked
    production files through the same builder helper that verifies binaries.
    A subsequent documentation-only Git revision need not change that native
    signature; exact source/build identities are still retained in receipts.
    """
    path = Path(path).resolve()
    if not path.is_file():
        raise Unavailable("reviewed I19 formal contract unavailable: " + str(path))
    value = _json(path)
    if value.get("schema") != SCHEMA:
        raise ValueError("unsupported I19 formal contract schema")
    review = value.get("review", {})
    if value.get("status") != "APPROVED" or review.get("approved") is not True:
        raise Unavailable("I19 execution contract is not investigator-approved")
    if not isinstance(review.get("reviewer"), str) or not review["reviewer"].strip():
        raise ValueError("I19 review must name its investigator")
    stamp = dt.datetime.fromisoformat(str(review["approved_utc"]).replace("Z", "+00:00"))
    if stamp.tzinfo is None or review.get("frozen_before_execution") is not True:
        raise ValueError("I19 review needs timezone-qualified approval and pre-execution freeze")
    if value.get("protocol_sha256") != _sha(PROTOCOL):
        raise ValueError("I19 approved numerical protocol identity differs")
    if value.get("action") not in ("import", "execute"):
        raise ValueError("I19 action must be import or execute")
    if value.get("qualification_root") is None:
        raise ValueError("I19 contract must identify the actual qualification campaign")
    qualification_root = Path(value["qualification_root"])
    if not qualification_root.is_absolute():
        qualification_root = path.parent / qualification_root
    if qualification_root.resolve() != Path(campaign_root).resolve():
        raise ValueError("I19 predecessor campaign identity differs")
    bound_input = _file(value["input"], path.parent, "frozen physical input")
    bound_executable = _file(value["executable"], path.parent, "native executable")
    if bound_input != Path(application_input).resolve() or bound_executable != Path(executable).resolve():
        raise ValueError("I19 contract is bound to a different input or executable")
    execution = value["execution"]
    positive = ("mpi_ranks", "workers", "allocated_cores", "core_hours",
                "memory_limit_bytes", "artifact_limit_bytes")
    for key in positive + ("host_free_floor_bytes", "disk_free_floor_bytes"):
        item = execution[key]
        if isinstance(item, bool) or not isinstance(item, (int, float)) or not math.isfinite(item):
            raise ValueError("non-numeric I19 execution control: " + key)
        if (key in positive and item <= 0) or (key not in positive and item < 0):
            raise ValueError("invalid I19 execution control: " + key)
    for key in ("mpi_ranks", "workers", "allocated_cores", "memory_limit_bytes",
                "artifact_limit_bytes", "host_free_floor_bytes", "disk_free_floor_bytes"):
        if int(execution[key]) != execution[key]:
            raise ValueError("I19 execution control must be an integer: " + key)
    if execution["mpi_ranks"] != mpi_ranks or not 1 <= mpi_ranks <= 1009:
        raise ValueError("I19 selected MPI layout differs from reviewed layout")
    if execution["workers"] not in (1, 2, 4) or execution["allocated_cores"] < mpi_ranks * execution["workers"]:
        raise ValueError("I19 allocation cannot supply its simultaneous MPI slots")
    files = {key: _file(value[key], path.parent, key) for key in
             ("sites_file", "forecast_file", "normalizations_file", "statistical_support_file")}
    # Pin the complete source map rather than Git HEAD alone. Existing lunar
    # implementation is intentionally uncommitted, so revision-only identity
    # would allow an older executable to be presented as current physics.
    from build_m2_qualification import production_source_identity
    signature = hashlib.sha256(json.dumps(production_source_identity(),
        sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    if value.get("production_source_signature_sha256") != signature:
        raise ValueError("I19 production source signature differs from reviewed native inputs")
    return value, files


def cluster_arguments(value, files, executable, application_input, campaign_root, mpi_ranks):
    """Adapt reviewed controls to the EXISTING cluster preflight API."""
    controls = value["execution"]
    return SimpleNamespace(execute=True, amps=Path(executable).resolve(),
        input=Path(application_input).resolve(), qualification_root=Path(campaign_root).resolve(),
        mpi_ranks=mpi_ranks, workers=int(controls["workers"]),
        allocated_cores=int(controls["allocated_cores"]), core_hours=float(controls["core_hours"]),
        **{key: int(controls[key]) for key in
           ("memory_limit_bytes", "artifact_limit_bytes", "host_free_floor_bytes", "disk_free_floor_bytes")},
        **files)


def _case_bytes(directory, case, executable, mpi_ranks, hashed, baseline=None,initial_qualification=None):
    """Validate authentic case provenance without rewriting its old evidence."""
    from integration_common import parse_receipt, moon_setting
    from m2_protocol import load_protocol, accounting_module, evaluate_cohorts
    required = ["result.json", "run_manifest.json", "receipt.txt", "application.in",
                "stdout.log", "stderr.log", "m2_budget.csv", "window_budget.csv"]
    if case["case"] == "cohort":
        required += ["m2_cohorts.csv", "m2_hazards.csv"]
    evidence = {}
    for name in required:
        file = directory / name
        if not file.is_file():
            raise Unavailable("selected I19 native artifact unavailable: " + str(file))
        evidence[name] = _sha(file)
    result = _json(directory / "result.json")
    native = parse_receipt(directory / "receipt.txt")
    identity = _json(directory / "run_manifest.json")
    if identity.get("execution_scope") != "cluster-numerical-design" or identity.get("case") != case:
        raise ValueError("I19 case is not its frozen full cluster design: " + case["key"])
    if identity.get("test_id") != "I19" or native.get("test_id") != "I19":
        raise ValueError("I19 native test identity differs: " + case["key"])
    if identity.get("executable", {}).get("sha256") != _sha(executable):
        raise ValueError("I19 native executable identity differs: " + case["key"])
    if identity.get("mpi_ranks") != mpi_ranks or native.get("mpi_ranks") != mpi_ranks or identity.get("openmp_threads") != 1:
        raise ValueError("I19 native MPI/OpenMP layout differs: " + case["key"])
    if identity.get("random_seed") != case["seed"] or native.get("seed") != case["seed"]:
        raise ValueError("I19 native stochastic seed differs: " + case["key"])
    if identity.get("protocol", {}).get("sha256") != _sha(PROTOCOL):
        raise ValueError("I19 native numerical protocol differs: " + case["key"])
    if moon_setting(directory / "application.in", "m2_verification_profile") != "cluster":
        raise ValueError("a short development case cannot become full I19 evidence")
    for key in ("dt_divisor", "weight_divisor"):
        if int(moon_setting(directory / "application.in", "m2_verification_" + key)) != case[key]:
            raise ValueError("I19 declared resolution control differs: " + case["key"])
    if not math.isclose(float(native["duration_s"]), float(load_protocol()["physical_duration_s"]),
                        rel_tol=1e-14, abs_tol=0):
        raise ValueError("I19 native physical duration differs from the approved long horizon")
    if not math.isclose(float(native["synchronized_timestep_s"]) * case["dt_divisor"],
                        float(native["coarse_timestep_s"]), rel_tol=1e-14, abs_tol=0):
        raise ValueError("I19 actual timestep does not implement its declared refinement")
    for species in ("na", "ar"):
        if not math.isclose(float(native[species + "_weight_atoms"]) * case["weight_divisor"],
                            float(native["reference_" + species + "_weight_atoms"]), rel_tol=1e-14):
            raise ValueError("I19 actual particle weight does not implement its declared refinement")
    if baseline is not None:
        for key in ("coarse_timestep_s", "reference_na_weight_atoms", "reference_ar_weight_atoms"):
            if not math.isclose(float(native[key]), float(baseline[key]), rel_tol=1e-14):
                raise ValueError("I19 physical baseline drifts across the numerical design: " + key)
    application_records = [(name, digest) for name, digest in identity.get("inputs", {}).items()
                           if Path(name).name == "application.in"]
    if len(application_records) != 1 or application_records[0][1] != evidence["application.in"]:
        raise ValueError("I19 case input snapshot identity differs: " + case["key"])
    if not identity.get("external_data"):
        raise ValueError("I19 native external-data identities are absent: " + case["key"])
    for name, expected in identity["external_data"].items():
        file = Path(name)
        if not file.is_absolute():
            from integration_common import SPICE_ROOT
            file = SPICE_ROOT / file
        key = (str(file.resolve()), expected)
        if key not in hashed:
            if not file.is_file():
                raise Unavailable("I19 consumed external data unavailable: " + str(file))
            if _sha(file) != expected:
                raise ValueError("I19 consumed external data changed: " + str(file))
            hashed.add(key)
    # A numerically plausible ledger cannot establish that the advertised
    # periodic state was actually applied. Bind this case's complete native
    # t0 state and physical staging before evaluating accounting or differences.
    from m2_initial_state_evidence import (Unavailable as StateUnavailable,
        verify_case_initial_state,STATE_SUFFIX)
    if initial_qualification is None:
        raise Unavailable("exact applied initial-state qualification unavailable for native import")
    try:
        applied=verify_case_initial_state(directory,initial_qualification)
    except StateUnavailable as error:
        raise Unavailable(str(error)) from error
    evidence['receipt.txt'+STATE_SUFFIX]=applied['sha256']
    if result.get("status") not in ("PASS", "FAIL"):
        raise Unavailable("I19 selected case has no evaluated scientific result: " + case["key"])
    # Recompute the signed atom AND kilogram identities from native rows.
    # Resolution equivalence alone would also accept equally wrong ledgers.
    accounting = accounting_module().evaluate(directory / "m2_budget.csv",
        require_coverage=False, allow_numerical_loss=True)
    failures = list(accounting["failures"])
    if result["status"] == "FAIL":
        failures.append({"case": case["key"], "reason": "retained native scientific result is FAIL"})
    if case["case"] == "cohort":
        cohort = evaluate_cohorts(directory / "m2_cohorts.csv", directory / "m2_hazards.csv",
                                  directory / "application.in")
        failures.extend(cohort["failures"])
    return evidence, native, failures


def evaluate_completed(directory, value, files, gate, executable, mpi_ranks,application_input=None):
    """Revalidate and independently score a complete, unchanged long campaign.

    The old numerical-score PASS is deliberately ignored. A failed native
    physics/accounting result remains a scientific FAIL, while absent native
    rows remain SKIPPED. Result files and cached summaries cannot substitute
    for actual signed atom/kg ledgers and passive-cohort records.
    """
    from m2_protocol import cases
    from m2_equivalence import score_campaign
    directory = Path(directory).resolve()
    summary_file = directory / "cluster_summary.json"
    plan_file = directory / "job_plan.json"
    preflight_file = directory / "preflight.json"
    if not all(p.is_file() for p in (summary_file, plan_file, preflight_file)):
        raise Unavailable("complete cluster summary/plan/preflight unavailable")
    summary, plan, previous_gate = (_json(p) for p in (summary_file, plan_file, preflight_file))
    if summary.get("execution") != "executed":
        raise Unavailable("cluster plan/preflight did not execute its native design")
    if plan.get("cases") != cases("cluster") or plan.get("protocol_sha256") != _sha(PROTOCOL):
        raise ValueError("completed cluster campaign differs from fixed 288-case design")
    if previous_gate.get("status") != "PASS":
        raise Unavailable("completed cluster campaign lacked passing scientific preflight")
    # Same-campaign means the ACTUAL retained predecessor and frozen-input
    # bytes agree, not merely that two unrelated directories both contain PASS.
    for key in ("I13", "I14", "I15", "I16", "I17", "I18",
                "sites_file", "forecast_file", "normalizations_file", "statistical_support_file",
                "initial_state_qualification", "initialization_cost", "normalization_physics"):
        if previous_gate.get("evidence", {}).get(key) != gate["evidence"].get(key):
            raise ValueError("completed I19 qualification evidence differs: " + key)
    native_root = directory / "cases"
    actual = {p.name for p in native_root.iterdir() if p.is_dir() and p.name != "monitors"} if native_root.is_dir() else set()
    expected = {case["key"] for case in cases("cluster")}
    if actual - expected:
        raise ValueError("unexpected native cases would change the frozen I19 family")
    if expected - actual:
        raise Unavailable("selected long I19 cases unavailable: " + ", ".join(sorted(expected - actual)))
    from m2_initial_state_evidence import Unavailable as StateUnavailable,load_qualification,bind_frozen_sites
    if application_input is None:
        raise Unavailable("explicit qualified physical template required to import the native initial-state family")
    try:
        qualification=load_qualification(_json(files['statistical_support_file']),application_input,executable,
            base=Path(files['statistical_support_file']).resolve().parent)
        qualification=bind_frozen_sites(qualification,_json(files['sites_file'])['sites'])
    except StateUnavailable as error:
        raise Unavailable(str(error)) from error
    evidence, hashed, baseline, failures = {}, set(), None, []
    for case in cases("cluster"):
        case_evidence, native, case_failures = _case_bytes(native_root / case["key"], case,
                                            executable, mpi_ranks, hashed, baseline,qualification)
        evidence[case["key"]] = case_evidence
        if baseline is None:
            baseline = native
        failures.extend({"case": case["key"], "failure": item} for item in case_failures)
    score = score_campaign(native_root, _json(files["normalizations_file"]),
                           _json(files["statistical_support_file"]),application_input=application_input,
                           executable=executable,support_base=Path(files['statistical_support_file']).resolve().parent,
                           sites=_json(files['sites_file'])['sites'])
    score.setdefault("failures", []).extend(failures)
    if failures and score.get('status') in ('PASS','FAIL'):
        score["status"] = "FAIL"
    return score, evidence


def write_comparisons(score, output, *, render=True):
    """Render evaluated normalized differences, not invented observations.

    Reference zero is the closed-form numerical-equivalence null. The CSV's
    uncertainty column carries the APPROVED simultaneous numerical confidence
    half-width, labelled as such in the title/score, not a measurement error.
    """
    rows = score.get("comparisons", [])
    if not rows:
        raise ValueError("evaluated I19 score has no numerical comparison rows")
    header = ["comparison_id", "x_value", "x_unit", "reference_value", "model_value",
              "y_unit", "uncertainty_value", "selection"]
    with (output / "comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=header, lineterminator="\n")
        writer.writeheader()
        for index, row in enumerate(rows):
            writer.writerow(dict(comparison_id="I19_normalized_equivalence",
                x_value=index, x_unit="frozen-comparison-index", reference_value=0,
                model_value=row["mean"], y_unit="normalized-difference",
                uncertainty_value=row["half_width"], selection="verification"))
    if render:
        sys.path.insert(0, str(Path(__file__).parent / "common"))
        from write_comparison_artifacts import render as draw
        draw(output / "comparison.csv", output / "comparison.png", output / "comparison.eps",
             "I19 numerical equivalence; bars are simultaneous numerical confidence intervals")


def run_formal_i19(executable, application_input, campaign_root, *, mpi_ranks=8,
                   allow_execute=False, wall_seconds=None):
    """Normal-runner entry: write canonical I19 artifacts and return its result.

    The caller must evaluate I13 first. This adapter repeats the full cluster
    predecessor/data/build/statistical/cost guards independently, so calling
    it directly cannot turn a supporting fixture into formal qualification.
    ``allow_execute`` must come from an explicit opt-in CLI, never development
    defaults. A reviewed import contract reads retained native results only.
    """
    from integration_common import git_value
    from m2_illumination import setting
    executable, application_input = Path(executable).resolve(), Path(application_input).resolve()
    campaign_root = Path(campaign_root).resolve()
    output = campaign_root / "I19"
    output.mkdir(parents=True, exist_ok=True)
    # Never overwrite scientific evidence on repeated dispatch. An empty/new
    # per-ID directory is needed; a previously preserved input/log is harmless
    # only if byte-identical and no completed result already exists.
    if (output / "result.json").exists():
        raise ValueError("refusing to overwrite an existing I19 result: " + str(output))
    snapshot = output / "application.in"
    if snapshot.exists() and _sha(snapshot) != _sha(application_input):
        raise ValueError("refusing to overwrite a different I19 input snapshot")
    if not snapshot.exists():
        shutil.copy2(application_input, snapshot)
    result = {"id": "I19", "status": "SKIPPED", "execution": "not-run",
              "scope": "formal-qualification", "metrics": {}, "artifacts": {}}
    identity = {"schema": "srcMoon-I19-formal-manifest-v1", "test_id": "I19",
        "execution_scope": "formal-qualification", "source_revision": git_value("rev-parse", "HEAD"),
        "executable": {"path": str(executable), "sha256": _sha(executable)} if executable.is_file() else None,
        "input": {"path": str(application_input), "sha256": _sha(application_input)},
        "protocol_sha256": _sha(PROTOCOL), "mpi_ranks": mpi_ranks, "openmp_threads": 1,
        "allow_long_execution": bool(allow_execute), "native_cases_executed_here": False}
    (output / "stdout.log").write_text("")
    (output / "stderr.log").write_text("")
    try:
        # Only actual absence is unavailable. Duplicate/malformed scalar
        # syntax is corruption and must remain ERROR, never a convenient SKIP.
        name = setting(application_input, "m2_formal_contract_file", "none")
        if not name or name.lower() == "none":
            raise Unavailable("explicit reviewed m2_formal_contract_file is unavailable")
        contract_path = Path(name)
        if not contract_path.is_absolute():
            contract_path = application_input.parent / contract_path
        from m2_contracts import resolve_contract
        contract_path, selection_evidence = resolve_contract(contract_path, "I19", unavailable=Unavailable)
        value, files = load_contract(contract_path, executable, application_input, campaign_root, mpi_ranks)
        identity["contract"] = {"path": str(contract_path.resolve()), "sha256": _sha(contract_path)}
        identity["contract_selection"] = selection_evidence
        shutil.copy2(contract_path, output / "execution_contract.json")
        from run_m2_cluster import preflight
        arguments = cluster_arguments(value, files, executable, application_input, campaign_root, mpi_ranks)
        gate = preflight(arguments)
        _write(output / "formal_preflight.json", gate)
        if gate["status"] == "ERROR":
            raise ValueError("; ".join(gate["errors"]))
        if gate["status"] != "PASS":
            raise Unavailable("; ".join(gate["gaps"]))
        identity["predecessor_and_build_evidence"] = gate["evidence"]
        if value["action"] == "execute":
            if not allow_execute:
                raise Unavailable("long I19 execution is not selected; explicit CLI opt-in required")
            command = [sys.executable, "-B", str(Path(__file__).with_name("run_m2_cluster.py")),
                       "--execute", "--output-dir", str(output / "long-campaign"),
                       "--amps", str(executable), "--input", str(application_input),
                       "--qualification-root", str(campaign_root)]
            for key in ("mpi_ranks", "workers", "allocated_cores", "core_hours", "memory_limit_bytes",
                        "artifact_limit_bytes", "host_free_floor_bytes", "disk_free_floor_bytes"):
                command += ["--" + key.replace("_", "-"), str(getattr(arguments, key))]
            for key, path in files.items():
                command += ["--" + key.replace("_", "-"), str(path)]
            identity["command"] = command
            # Reuse the existing process-tree supervisor and cluster worker.
            # The shorter enclosing caller deadline still wins over allocation.
            from campaign_budget import CampaignBudget
            seconds = arguments.core_hours / arguments.allocated_cores * 3600
            if wall_seconds is not None:
                seconds = min(seconds, float(wall_seconds))
            budget = CampaignBudget(seconds, max(.1, min(600, seconds * .05)),
                memory_bytes=arguments.memory_limit_bytes, artifact_bytes=arguments.artifact_limit_bytes,
                host_free_floor_bytes=arguments.host_free_floor_bytes,
                disk_free_floor_bytes=arguments.disk_free_floor_bytes)
            execution = budget.run(command, ROOT, output / "stdout.log", output=output,
                                   allocated_cores=arguments.allocated_cores)
            _write(output / "execution_budget.json", execution)
            if execution["returncode"] not in (0, 1):
                raise RuntimeError(execution.get("reason") or "selected cluster execution failed")
            directory = output / "long-campaign/campaign"
            identity["native_cases_executed_here"] = True
        else:
            # A contract names both the retained worker directory and its
            # summary bytes. Rescoring below still checks all native artifacts.
            summary = _file(value["completed_summary"], contract_path.parent, "completed cluster summary")
            directory = summary.parent
        score, evidence = evaluate_completed(directory, value, files, gate, executable, mpi_ranks,application_input)
        _write(output / "independent_numerical_score.json", score)
        _write(output / "native_evidence_manifest.json", evidence)
        if score["status"]=="ERROR":
            raise ValueError(score.get("reason","initial-state authentication failed"))
        if score["status"] not in ("PASS", "FAIL"):
            raise Unavailable(score.get("reason", "full I19 numerical evidence is not qualified"))
        write_comparisons(score, output)
        # A single synthetic native receipt would hide 288 distinct seeds and
        # resolutions. Bind the real native receipts directly instead. The
        # aggregate is explicitly orchestration provenance, not a kernel call.
        native_evidence = [{"path": str(directory / "cases" / key / "receipt.txt"),
                            "sha256": records["receipt.txt"]}
                           for key, records in evidence.items()]
        identity["native_evidence"] = native_evidence
        # Keep the receipt family exactly288 for the global receipt contract,
        # and expose ALL separately validated native/source/state bytes for
        # I20. Receipt-only membership cannot authenticate a cohort CSV or the
        # applied state from which it evolved.
        identity["native_source_evidence"]=[{"path":str(directory/"cases"/key/name),"sha256":digest}
            for key,records in evidence.items() for name,digest in records.items()]
        identity["applied_initial_state"]=score.get("applied_initial_state")
        identity["native_campaign_directory"] = str(directory)
        identity["native_evidence_manifest_sha256"] = _sha(output / "native_evidence_manifest.json")
        _write(output / "aggregate_native_receipt.json", {
            "schema": "srcMoon-I19-aggregate-native-evidence-v1", "test_id": "I19",
            "kind": "orchestration-receipt-not-native-kernel-output", "mpi_ranks": mpi_ranks,
            "native_case_count": len(native_evidence), "native_evidence": native_evidence,
            "native_source_evidence":identity["native_source_evidence"],
            "applied_initial_state":identity["applied_initial_state"]})
        result.update(status=score["status"], execution="executed" if identity["native_cases_executed_here"] else "imported-linked-results",
            reason="complete native I19 campaign revalidated and independently scored",
            metrics={"native_case_count": 288, "numerical_comparison_count": len(score["comparisons"]),
                     "failures": score.get("failures", []), "family_size_bound": score.get("family_size_bound")})
    except Unavailable as error:
        result.update(status="SKIPPED", reason=str(error))
    except (OSError, ValueError, KeyError, TypeError, RuntimeError, subprocess.TimeoutExpired) as error:
        result.update(status="ERROR", reason=str(error))
        (output / "stderr.log").write_text(str(error) + "\n")
    identity["formal_status"] = result["status"]
    identity["execution"] = result["execution"]
    _write(output / "run_manifest.json", identity)
    for path in output.iterdir():
        if path.is_file() and path.name != "result.json":
            result["artifacts"][path.name] = str(path)
    with (output / "stdout.log").open("a") as stream:
        stream.write(result["status"] + " I19: " + result["reason"] + "\n")
    _write(output / "result.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amps", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--qualification-root", type=Path, required=True)
    parser.add_argument("--mpi-ranks", type=int, default=8)
    parser.add_argument("--allow-long-execution", action="store_true")
    args = parser.parse_args()
    result = run_formal_i19(args.amps, args.input, args.qualification_root,
                            mpi_ranks=args.mpi_ranks, allow_execute=args.allow_long_execution)
    print(result["status"] + " I19: " + result["reason"])
    return {"PASS": 0, "FAIL": 1, "ERROR": 2, "SKIPPED": 77}[result["status"]]


if __name__ == "__main__":
    raise SystemExit(main())
