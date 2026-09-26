#!/usr/bin/env python3
"""Run a reproducible standalone SEP-in-geospace event campaign.

The runner owns orchestration and evidence only.  It renders one immutable input
for every ``(field model, epoch)``, executes the existing Step-7 standalone AMPS
path, and refuses to report PASS when an expected artifact, numerical gate, or
required observation comparison is absent.  It does not alter AMPS tolerances or
reinterpret a failed physics calculation.

Typical production invocation from the AMPS repository root::

  python3 srcEarth/standalone_campaign/run_campaign.py \
    --manifest srcEarth/examples/standalone_step8_campaign/campaign.json \
    --amps ./amps --output-dir test_output/step8 --restart
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

# Direct execution sets sys.path to this directory, not its parent.  Add only the
# package parent; no installed package or AMPS build-tree mutation is required.
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from standalone_campaign.adapters import (  # type: ignore
        adapt_observations,
        compare_with_observations,
        validate_instrument_response,
        write_canonical_csv,
        write_comparison_csv,
    )
    from standalone_campaign.campaign import (  # type: ignore
        CampaignError,
        artifact_inventory,
        atomic_write_json,
        cache_resources,
        campaign_epochs,
        epoch_token,
        evaluate_numerical_gates,
        format_utc,
        load_and_validate_manifest,
        parse_tecplot_table,
        parse_utc,
        read_prediction_csv,
        render_template,
        sha256_bytes,
        sha256_file,
        template_values,
        termination_evidence,
        verified_restart,
    )
else:
    from .adapters import (
        adapt_observations,
        compare_with_observations,
        validate_instrument_response,
        write_canonical_csv,
        write_comparison_csv,
    )
    from .campaign import (
        CampaignError,
        artifact_inventory,
        atomic_write_json,
        cache_resources,
        campaign_epochs,
        epoch_token,
        evaluate_numerical_gates,
        format_utc,
        load_and_validate_manifest,
        parse_tecplot_table,
        parse_utc,
        read_prediction_csv,
        render_template,
        sha256_bytes,
        sha256_file,
        template_values,
        termination_evidence,
        verified_restart,
    )


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path,
                        help="Frozen Step-8 campaign JSON manifest")
    parser.add_argument("--amps", type=Path, default=Path("./amps"),
                        help="Standalone AMPS executable")
    parser.add_argument("--output-dir", required=True, type=Path,
                        help="Campaign result directory")
    parser.add_argument("--restart", action="store_true",
                        help="Skip only hash-verified PASS runs with the same fingerprint")
    parser.add_argument("--validate-only", action="store_true",
                        help="Validate/cache all references and normalize observations; do not render/run AMPS")
    parser.add_argument("--dry-run", action="store_true",
                        help="Render inputs and commands but do not execute AMPS")
    parser.add_argument("--mpi-ranks", type=int, default=None,
                        help="Override the manifest MPI count without changing physics gates")
    parser.add_argument("--threads", type=int, default=None,
                        help="Override the manifest thread count without changing physics gates")
    return parser.parse_args(argv)


def _resource_by_id(manifest: Mapping[str, object]) -> Dict[str, Mapping[str, object]]:
    return {
        str(item["id"]): item
        for item in manifest["_verified_resources"]  # type: ignore[index]
    }


def _instrument_observations(
    manifest: Mapping[str, object], cached: Mapping[str, Path], output_root: Path
) -> Tuple[Dict[str, List[Dict[str, object]]], List[Dict[str, object]]]:
    """Normalize all observations before any expensive AMPS launch."""

    by_instrument: Dict[str, List[Dict[str, object]]] = {}
    records: List[Dict[str, object]] = []
    event = manifest["event"]  # type: ignore[index]
    event_start = parse_utc(str(event["start_utc"]))
    event_end = parse_utc(str(event["end_utc"]))
    for raw in manifest["instruments"]:  # type: ignore[index]
        instrument = raw  # type: ignore[assignment]
        identifier = str(instrument["id"])
        observation_path = cached[str(instrument["observation_resource"])]
        response_path = cached[str(instrument["response_resource"])]
        response_validation = validate_instrument_response(response_path, instrument)
        all_rows = adapt_observations(observation_path, instrument)
        # Reference archives commonly cover more than the selected event.  Use
        # the manifest's preregistered midpoint window; silently scoring adjacent
        # days would change sample counts and metrics when an archive is updated.
        rows = [
            row for row in all_rows
            if event_start <= parse_utc(str(row["midpoint_utc"])) <= event_end
        ]
        if not rows:
            raise CampaignError(
                "instrument %s has no quality-controlled rows in the event window"
                % identifier
            )
        by_instrument[identifier] = rows
        write_canonical_csv(output_root / "observations" / (identifier + ".csv"), rows)
        records.append({
            "instrument_id": identifier,
            "adapter": instrument["adapter"],
            "observation_resource": instrument["observation_resource"],
            "response_resource": instrument["response_resource"],
            "response_validation": response_validation,
            "row_count": len(rows),
            "archive_row_count": len(all_rows),
            "status": "PASS",
        })
    return by_instrument, records


def _expand_command(
    tokens: Sequence[str], values: Mapping[str, object]
) -> List[str]:
    result = []
    for token in tokens:
        try:
            result.append(token.format(**values))
        except KeyError as exc:
            raise CampaignError("execution.command uses unknown field: %s" % exc.args[0]) from exc
    return result


def _append_progress(path: Path, event: Mapping[str, object]) -> None:
    """Append one complete JSON line and flush it for live batch monitoring."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as stream:
        stream.write(json.dumps(dict(event), sort_keys=True) + "\n")
        stream.flush()
        os.fsync(stream.fileno())


def _run_process(
    command: Sequence[str], run_dir: Path, log_path: Path, environment: Mapping[str, str],
    prefix: str,
) -> int:
    """Stream merged stdout/stderr to both terminal and an unmodified run log."""

    with log_path.open("w", encoding="utf-8") as log:
        process = subprocess.Popen(
            list(command), cwd=str(run_dir), stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, text=True, bufsize=1, env=dict(environment),
        )
        assert process.stdout is not None
        with process.stdout:
            for line in process.stdout:
                log.write(line)
                log.flush()
                print("[%s] %s" % (prefix, line.rstrip()), flush=True)
        return process.wait()


def _gate_failure(gates: Sequence[Mapping[str, object]]) -> bool:
    return any(item.get("required", True) and item.get("status") != "PASS" for item in gates)


def _write_extracted_predictions(
    run_dir: Path, epoch, execution: Mapping[str, object]
) -> None:
    """Convert declared AMPS Tecplot columns into the canonical prediction CSV.

    An artifact pattern must resolve to exactly one file.  ``COLUMN`` copies one
    detector/product column; ``RATIO`` divides two columns (for example physical
    east/west rates).  The operation and row are frozen in the campaign manifest,
    so the runner never selects the most favorable location or detector after a
    comparison has been viewed.
    """

    extractors = execution.get("prediction_extractors", [])
    target_name = execution.get("prediction_artifact")
    if not extractors:
        return
    if not target_name:
        raise CampaignError("prediction_extractors require prediction_artifact")
    rows: List[Dict[str, object]] = []
    for index, raw in enumerate(extractors):
        extractor = raw  # type: ignore[assignment]
        matches = sorted(path for path in run_dir.glob(str(extractor["artifact"])) if path.is_file())
        if len(matches) != 1:
            raise CampaignError(
                "prediction extractor %d expected exactly one %s artifact, found %d"
                % (index, extractor["artifact"], len(matches))
            )
        variables, table = parse_tecplot_table(matches[0])
        expected_variables = [str(name) for name in extractor["expected_variables"]]
        if variables != expected_variables:
            raise CampaignError(
                "prediction extractor %d variable/order mismatch: expected %s, got %s"
                % (index, expected_variables, variables)
            )
        expected_rows = int(extractor["expected_row_count"])
        if len(table) != expected_rows:
            raise CampaignError(
                "prediction extractor %d expected %d rows, got %d"
                % (index, expected_rows, len(table))
            )
        columns = {name: position for position, name in enumerate(variables)}
        row_index = int(extractor.get("row_index", 0))
        if row_index < 0 or row_index >= len(table):
            raise CampaignError("prediction extractor row_index is outside the artifact")
        source = table[row_index]

        def column(name: object) -> float:
            label = str(name)
            if label not in columns:
                raise CampaignError("prediction artifact lacks column: %s" % label)
            return source[columns[label]]

        operation = str(extractor.get("operation", "COLUMN")).upper()
        if operation == "COLUMN":
            value = column(extractor["column"])
            lower = column(extractor["lower_column"]) if extractor.get("lower_column") else None
            upper = column(extractor["upper_column"]) if extractor.get("upper_column") else None
        else:
            numerator = column(extractor["numerator_column"])
            denominator = column(extractor["denominator_column"])
            if denominator == 0.0:
                raise CampaignError("prediction ratio denominator is zero")
            value = numerator / denominator
            lower = None
            upper = None
        if not math.isfinite(value):
            raise CampaignError("prediction extractor produced a nonfinite value")
        if lower is not None and upper is not None and not lower <= value <= upper:
            raise CampaignError(
                "prediction extractor produced inconsistent lower/value/upper bounds"
            )
        rows.append({
            "utc": format_utc(epoch), "instrument_id": extractor["instrument_id"],
            "platform": extractor["platform"], "channel": extractor["channel"],
            "direction": str(extractor["direction"]).upper(),
            "quantity": extractor["quantity"], "units": extractor["units"],
            "value": value,
            "lower": lower, "upper": upper,
        })
    target = run_dir / str(target_name)
    target.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "utc", "instrument_id", "platform", "channel", "direction",
        "quantity", "units", "value", "lower", "upper",
    ]
    with target.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _run_one(
    manifest: Mapping[str, object], cached: Mapping[str, Path], model: str, epoch,
    amps: Path, output_root: Path, mpi_ranks: int, threads: int, restart: bool,
    dry_run: bool, progress_path: Path, executable_digest: str,
) -> Dict[str, object]:
    execution = manifest["execution"]  # type: ignore[index]
    run_dir = output_root / "runs" / model / epoch_token(epoch)
    run_dir.mkdir(parents=True, exist_ok=True)
    input_path = run_dir / "amps.in"
    status_path = run_dir / "run_status.json"
    template_path = cached[str(execution["input_template_resource"])]  # type: ignore[index]
    template = template_path.read_text(encoding="utf-8")
    values = template_values(manifest, cached, model, epoch, run_dir)
    rendered = render_template(template, values)
    input_path.write_text(rendered, encoding="utf-8")

    command_values = {
        "amps": str(amps), "input": str(input_path), "output_dir": str(run_dir),
        "mpi_ranks": mpi_ranks, "threads": threads, "model": model,
        "epoch": format_utc(epoch), "epoch_token": epoch_token(epoch),
    }
    command = _expand_command(execution["command"], command_values)  # type: ignore[arg-type,index]
    fingerprint_document = {
        "campaign_manifest_sha256": manifest["_manifest_sha256"],
        "executable_sha256": executable_digest,
        "model": model,
        "epoch_utc": format_utc(epoch),
        "input_sha256": sha256_bytes(rendered.encode("utf-8")),
        "command": command,
        "mpi_ranks": mpi_ranks,
        "threads": threads,
    }
    fingerprint = sha256_bytes(
        json.dumps(fingerprint_document, sort_keys=True, separators=(",", ":")).encode("utf-8")
    )

    if restart and verified_restart(status_path, fingerprint):
        with status_path.open(encoding="utf-8") as stream:
            previous = json.load(stream)
        previous["restart_action"] = "SKIPPED_VERIFIED_PASS"
        _append_progress(progress_path, {
            "event": "SKIP", "model": model, "epoch_utc": format_utc(epoch),
            "reason": "verified PASS fingerprint and artifact hashes",
        })
        return previous

    if dry_run:
        status = {
            "status": "NOT_RUN", "restart_action": "DRY_RUN", "fingerprint": fingerprint,
            "model": model, "epoch_utc": format_utc(epoch), "command": command,
            "input_sha256": fingerprint_document["input_sha256"], "artifacts": [],
            "missing_artifacts": list(execution["expected_artifacts"]),  # type: ignore[index]
            "numerical_gates": [], "termination": {},
        }
        atomic_write_json(status_path, status)
        return status

    started = time.time()
    running = {
        "status": "RUNNING", "fingerprint": fingerprint, "model": model,
        "epoch_utc": format_utc(epoch), "command": command,
        "input_sha256": fingerprint_document["input_sha256"], "started_unix": started,
    }
    atomic_write_json(status_path, running)
    _append_progress(progress_path, {
        "event": "START", "model": model, "epoch_utc": format_utc(epoch),
        "run_dir": str(run_dir),
    })
    environment = dict(os.environ)
    environment["OMP_NUM_THREADS"] = str(threads)
    launch_error = None
    try:
        return_code = _run_process(
            command, run_dir, run_dir / "amps.log", environment,
            "%s %s" % (model, epoch_token(epoch)),
        )
    except OSError as exc:
        # A missing mpirun/wrapper or other launch failure still needs a durable
        # FAIL record.  Leaving only the earlier RUNNING file would make batch
        # diagnosis and restart ambiguous and would violate the Step-8 evidence
        # contract even though no particle calculation occurred.
        return_code = 127
        launch_error = "%s: %s" % (type(exc).__name__, exc)
        with (run_dir / "amps.log").open("a", encoding="utf-8") as log:
            log.write("STEP-8 LAUNCH ERROR: %s\n" % launch_error)

    extraction_error = None
    if return_code == 0:
        try:
            _write_extracted_predictions(run_dir, epoch, execution)
        except CampaignError as exc:
            extraction_error = str(exc)

    prediction_validation_error = None
    prediction_relative = execution.get("prediction_artifact")
    if return_code == 0 and extraction_error is None and prediction_relative:
        prediction_path = run_dir / str(prediction_relative)
        if prediction_path.is_file():
            try:
                # Validate externally generated canonical predictions at the
                # per-run boundary.  This keeps malformed CSV evidence in the
                # run status instead of aborting later during aggregation.
                read_prediction_csv(prediction_path)
            except CampaignError as exc:
                prediction_validation_error = str(exc)

    artifacts, missing = artifact_inventory(
        run_dir, execution["expected_artifacts"]  # type: ignore[arg-type,index]
    )
    try:
        termination = termination_evidence(run_dir)
    except CampaignError as exc:
        termination = {
            "files": [], "table_count": 0, "max_unresolved_fraction": None,
            "total_sampled": 0, "total_resolved": 0, "termination_counts": {},
            "error": str(exc),
        }
    gates = evaluate_numerical_gates(
        run_dir, termination,
        manifest["validation"]["numerical_gates"],  # type: ignore[index]
    )
    passed = (
        return_code == 0 and launch_error is None and extraction_error is None
        and prediction_validation_error is None and not missing
        and not _gate_failure(gates)
    )
    status = {
        "status": "PASS" if passed else "FAIL", "fingerprint": fingerprint,
        "model": model, "epoch_utc": format_utc(epoch), "command": command,
        "input_sha256": fingerprint_document["input_sha256"],
        "executable_sha256": executable_digest, "return_code": return_code,
        "elapsed_seconds": round(time.time() - started, 6), "artifacts": artifacts,
        "missing_artifacts": missing, "numerical_gates": gates,
        "termination": termination, "prediction_extraction_error": extraction_error,
        "prediction_validation_error": prediction_validation_error,
        "launch_error": launch_error,
        "restart_action": "EXECUTED",
    }
    atomic_write_json(status_path, status)
    _append_progress(progress_path, {
        "event": status["status"], "model": model, "epoch_utc": format_utc(epoch),
        "return_code": return_code, "missing_artifact_count": len(missing),
        "failed_gate_count": sum(item["status"] == "FAIL" for item in gates),
    })
    return status


def _load_predictions(
    manifest: Mapping[str, object], output_root: Path, run_statuses: Sequence[Mapping[str, object]]
) -> Dict[str, List[Dict[str, object]]]:
    execution = manifest["execution"]  # type: ignore[index]
    relative = execution.get("prediction_artifact")  # type: ignore[union-attr]
    by_model: Dict[str, List[Dict[str, object]]] = {}
    if not relative:
        return by_model
    for status in run_statuses:
        if status.get("status") != "PASS":
            continue
        model = str(status["model"])
        run_dir = output_root / "runs" / model / epoch_token_from_text(str(status["epoch_utc"]))
        path = run_dir / str(relative)
        if not path.is_file():
            continue
        rows = read_prediction_csv(path)
        for row in rows:
            row["model"] = model
            row["snapshot_epoch_utc"] = status["epoch_utc"]
        by_model.setdefault(model, []).extend(rows)
    return by_model


def epoch_token_from_text(value: str) -> str:
    # Keep this helper local to avoid exposing a text/token conversion in the
    # physics contract.  parse/format validation already occurred in the runner.
    from datetime import datetime, timezone
    text = value[:-1] + "+00:00" if value.endswith("Z") else value
    return datetime.fromisoformat(text).astimezone(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def _observation_gate_results(
    quantity: str, metrics: Mapping[str, object], required: bool
) -> List[Dict[str, object]]:
    """Apply the preregistered roadmap thresholds relevant to one quantity."""

    specifications = []  # type: List[Tuple[str, str, float, bool]]
    # A required comparison may not silently discard uncovered observation
    # windows.  Quality-control exclusions have already been applied by the
    # adapter and frozen in the manifest; everything that remains must match.
    specifications.append(("comparison_fraction", ">=", 1.0, False))
    if quantity in ("differential_flux", "integral_flux", "count_rate"):
        specifications.extend([
            ("median_abs_log10_ratio", "<=", 0.30, False),
            # The roadmap requires 70% pooled across an event and at least 60%
            # for each instrument.  This routine evaluates one instrument; the
            # stricter pooled check is reported separately by event analyses.
            ("factor_two_fraction", ">=", 0.60, False),
        ])
    elif quantity == "cutoff_latitude":
        specifications.extend([
            ("cutoff_mae_deg", "<=", 2.0, False),
            ("cutoff_bias_deg", "<=", 1.0, True),
        ])
    elif quantity == "directional_ratio":
        specifications.extend([
            ("directional_log10_rmse", "<=", 0.30, False),
            ("directional_sign_fraction", ">=", 0.80, False),
        ])
    if metrics.get("interval_coverage_fraction") is not None:
        specifications.extend([
            ("interval_coverage_fraction_min", ">=", 0.80, False),
            ("interval_coverage_fraction_max", "<=", 0.98, False),
        ])

    result: List[Dict[str, object]] = []
    for name, operator, threshold, absolute in specifications:
        metric_name = name
        if name.endswith("_min") or name.endswith("_max"):
            metric_name = "interval_coverage_fraction"
        value_object = metrics.get(metric_name)
        if value_object is None:
            result.append({
                "id": name, "status": "FAIL" if required else "NOT_EVALUATED",
                "value": None, "operator": operator, "threshold": threshold,
                "required": required, "error": "metric not produced",
            })
            continue
        value = abs(float(value_object)) if absolute else float(value_object)
        passed = value <= threshold if operator == "<=" else value >= threshold
        result.append({
            "id": name, "status": "PASS" if passed else "FAIL", "value": value,
            "operator": operator, "threshold": threshold, "required": required,
        })
    return result


def _write_flat_prediction_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    fields = [
        "model", "snapshot_epoch_utc", "utc", "instrument_id", "platform",
        "channel", "direction", "quantity", "units", "value", "lower", "upper",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _score_observations(
    manifest: Mapping[str, object], observations: Mapping[str, Sequence[Mapping[str, object]]],
    predictions_by_model: Mapping[str, Sequence[Mapping[str, object]]], output_root: Path,
) -> List[Dict[str, object]]:
    summaries: List[Dict[str, object]] = []
    instruments = {str(item["id"]): item for item in manifest["instruments"]}  # type: ignore[index]
    for model in manifest["models"]:  # type: ignore[index]
        model_text = str(model)
        predictions = predictions_by_model.get(model_text, [])
        pooled_flux_logs: List[float] = []
        pooled_flux_required = any(
            bool(item.get("comparison_required", True))
            and str(item.get("quantity", ""))
            in ("differential_flux", "integral_flux", "count_rate")
            for item in instruments.values()
        )
        for identifier, observation_rows in observations.items():
            instrument = instruments[identifier]
            required = bool(instrument.get("comparison_required", True))
            selected = [row for row in predictions if row["instrument_id"] == identifier]
            comparison, metrics = compare_with_observations(selected, observation_rows)
            quantity = str(instrument["quantity"])
            if quantity in ("differential_flux", "integral_flux", "count_rate"):
                pooled_flux_logs.extend(
                    float(row["log10_ratio"])
                    for row in comparison
                    if row.get("log10_ratio") is not None
                )
            gates = _observation_gate_results(quantity, metrics, required)
            if not comparison:
                gates.append({
                    "id": "comparison_produced", "status": "FAIL" if required else "NOT_EVALUATED",
                    "value": 0, "operator": ">", "threshold": 0, "required": required,
                    "error": "no model/observation comparison rows were produced",
                })
            else:
                gates.append({
                    "id": "comparison_produced", "status": "PASS", "value": len(comparison),
                    "operator": ">", "threshold": 0, "required": required,
                })
            destination = output_root / "comparisons" / model_text / identifier
            write_comparison_csv(destination / "comparison.csv", comparison)
            summary = {
                "model": model_text, "instrument_id": identifier, "quantity": quantity,
                "required": required, "metrics": metrics, "gates": gates,
                "status": "FAIL" if _gate_failure(gates) else "PASS",
                "comparison_path": str((destination / "comparison.csv").relative_to(output_root)),
            }
            atomic_write_json(destination / "summary.json", summary)
            summaries.append(summary)

        # The roadmap has two distinct amplitude-coverage gates: at least 60%
        # for each individual instrument (checked above) and at least 70%
        # pooled across the event.  Keeping the pooled decision as a separate
        # record prevents a large platform from hiding which per-instrument
        # comparison failed.
        if pooled_flux_required:
            pooled_fraction = (
                sum(abs(value) <= math.log10(2.0) for value in pooled_flux_logs)
                / len(pooled_flux_logs)
                if pooled_flux_logs else None
            )
            pooled_gate = {
                "id": "pooled_factor_two_fraction",
                "status": (
                    "PASS"
                    if pooled_fraction is not None and pooled_fraction >= 0.70
                    else "FAIL"
                ),
                "value": pooled_fraction,
                "operator": ">=",
                "threshold": 0.70,
                "required": True,
            }
            pooled_summary = {
                "model": model_text,
                "instrument_id": "POOLED_EVENT",
                "quantity": "flux_amplitude",
                "required": True,
                "metrics": {
                    "log_comparison_count": len(pooled_flux_logs),
                    "factor_two_fraction": pooled_fraction,
                },
                "gates": [pooled_gate],
                "status": pooled_gate["status"],
                "comparison_path": None,
            }
            destination = output_root / "comparisons" / model_text / "POOLED_EVENT"
            atomic_write_json(destination / "summary.json", pooled_summary)
            summaries.append(pooled_summary)
    return summaries


def _aggregate_termination(run_statuses: Sequence[Mapping[str, object]]) -> Dict[str, object]:
    maximum = None
    sampled = 0
    resolved = 0
    counts: Dict[str, int] = {}
    for status in run_statuses:
        evidence = status.get("termination", {})
        value = evidence.get("max_unresolved_fraction") if isinstance(evidence, Mapping) else None
        if value is not None:
            maximum = float(value) if maximum is None else max(maximum, float(value))
        if isinstance(evidence, Mapping):
            sampled += int(evidence.get("total_sampled", 0))
            resolved += int(evidence.get("total_resolved", 0))
            for name, count in evidence.get("termination_counts", {}).items():
                counts[str(name)] = counts.get(str(name), 0) + int(count)
    return {
        "run_count": len(run_statuses), "max_unresolved_fraction": maximum,
        "total_sampled": sampled, "total_resolved": resolved,
        "termination_counts": counts,
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.validate_only and args.dry_run:
        raise SystemExit("--validate-only and --dry-run are mutually exclusive")
    manifest = load_and_validate_manifest(args.manifest)
    output_root = args.output_dir.expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    cached = cache_resources(manifest["_verified_resources"], output_root / "cache")  # type: ignore[arg-type,index]
    observations, adapter_records = _instrument_observations(manifest, cached, output_root)

    execution = manifest["execution"]  # type: ignore[index]
    mpi_ranks = args.mpi_ranks if args.mpi_ranks is not None else int(execution.get("mpi_ranks", 1))
    threads = args.threads if args.threads is not None else int(execution.get("threads", 1))
    if mpi_ranks <= 0 or threads <= 0:
        raise CampaignError("MPI ranks and threads must be positive")

    amps = args.amps.expanduser().resolve()
    if not (args.validate_only or args.dry_run):
        if not amps.is_file() or not os.access(str(amps), os.X_OK):
            raise CampaignError("AMPS executable is missing or not executable: %s" % amps)
        executable_digest = sha256_file(amps)
    else:
        executable_digest = sha256_file(amps) if amps.is_file() else "NOT_EVALUATED"

    preflight = {
        "schema_version": manifest["schema_version"],
        "campaign_id": manifest["campaign_id"], "event_id": manifest["event_id"],
        "role": manifest["role"],
        "frozen": manifest["frozen"], "manifest_path": manifest["_manifest_path"],
        "manifest_sha256": manifest["_manifest_sha256"],
        "normalization_policy": manifest["normalization_policy"],
        "resource_count": len(cached), "resources": [
            {
                "id": item["id"], "kind": item["kind"], "sha256": item["sha256"],
                "size_bytes": item["size_bytes"],
                "provenance": item["provenance"],
                "cached_path": str(cached[str(item["id"])]),
            }
            for item in manifest["_verified_resources"]  # type: ignore[index]
        ],
        "observation_adapters": adapter_records,
        "mpi_ranks": mpi_ranks, "threads": threads,
        "status": "PASS",
    }
    atomic_write_json(output_root / "campaign_preflight.json", preflight)
    if args.validate_only:
        summary = {
            "campaign_id": manifest["campaign_id"], "event_id": manifest["event_id"],
            "status": "PREFLIGHT_PASS",
            "runs_planned": len(manifest["models"]) * len(campaign_epochs(manifest["event"])),  # type: ignore[arg-type,index]
            "runs_executed": 0, "comparisons_produced": 0,
        }
        atomic_write_json(output_root / "campaign_summary.json", summary)
        print("STEP-8 CAMPAIGN PREFLIGHT PASS: %s" % manifest["campaign_id"])
        return 0

    progress_path = output_root / "progress.jsonl"
    statuses: List[Dict[str, object]] = []
    for model in manifest["models"]:  # type: ignore[index]
        for epoch in campaign_epochs(manifest["event"]):  # type: ignore[arg-type,index]
            statuses.append(_run_one(
                manifest, cached, str(model), epoch, amps, output_root, mpi_ranks,
                threads, args.restart, args.dry_run, progress_path, executable_digest,
            ))

    termination = _aggregate_termination(statuses)
    atomic_write_json(output_root / "termination_summary.json", termination)
    convergence = {
        "runs": [
            {
                "model": item.get("model"), "epoch_utc": item.get("epoch_utc"),
                "status": item.get("status"), "numerical_gates": item.get("numerical_gates", []),
            }
            for item in statuses
        ]
    }
    atomic_write_json(output_root / "convergence_summary.json", convergence)

    observation_summaries: List[Dict[str, object]] = []
    predictions_by_model: Dict[str, List[Dict[str, object]]] = {}
    if not args.dry_run:
        predictions_by_model = _load_predictions(manifest, output_root, statuses)
        all_predictions = [row for rows in predictions_by_model.values() for row in rows]
        _write_flat_prediction_csv(output_root / "campaign_predictions.csv", all_predictions)
        observation_summaries = _score_observations(
            manifest, observations, predictions_by_model, output_root
        )
        atomic_write_json(output_root / "observation_summary.json", {
            "comparisons": observation_summaries,
            "status": "FAIL" if any(item["status"] == "FAIL" for item in observation_summaries) else "PASS",
        })

    run_fail = any(item.get("status") == "FAIL" for item in statuses)
    comparison_fail = any(item.get("status") == "FAIL" for item in observation_summaries)
    if args.dry_run:
        overall = "NOT_RUN"
    else:
        overall = "FAIL" if run_fail or comparison_fail else "PASS"
    summary = {
        "campaign_id": manifest["campaign_id"], "event_id": manifest["event_id"],
        "role": manifest["role"], "frozen": manifest["frozen"],
        "status": overall,
        "manifest_sha256": manifest["_manifest_sha256"],
        "executable_sha256": executable_digest,
        "run_count": len(statuses),
        "run_pass_count": sum(item.get("status") == "PASS" for item in statuses),
        "run_fail_count": sum(item.get("status") == "FAIL" for item in statuses),
        "run_not_executed_count": sum(item.get("status") == "NOT_RUN" for item in statuses),
        "restart_skip_count": sum(item.get("restart_action") == "SKIPPED_VERIFIED_PASS" for item in statuses),
        "observation_comparison_count": len(observation_summaries),
        "observation_comparison_fail_count": sum(item.get("status") == "FAIL" for item in observation_summaries),
        "observation_comparison_row_count": sum(
            int(item.get("metrics", {}).get("comparison_count", 0))
            for item in observation_summaries
            if item.get("instrument_id") != "POOLED_EVENT"
        ),
        "required_observation_comparison_count": sum(
            bool(item.get("required"))
            for item in observation_summaries
            if item.get("instrument_id") != "POOLED_EVENT"
        ),
        "required_observation_comparison_fail_count": sum(
            bool(item.get("required")) and item.get("status") == "FAIL"
            for item in observation_summaries
            if item.get("instrument_id") != "POOLED_EVENT"
        ),
        "products": manifest["products"], "models": manifest["models"],
        "normalization_policy": manifest["normalization_policy"],
        "resource_count": len(cached),
        "all_required_numerical_gates_pass": (
            not args.dry_run
            and bool(statuses)
            and all(
                bool(item.get("numerical_gates"))
                and not _gate_failure(item.get("numerical_gates", []))
                for item in statuses
            )
        ),
    }
    atomic_write_json(output_root / "campaign_summary.json", summary)
    print("STEP-8 CAMPAIGN %s: %s" % (overall, manifest["campaign_id"]), flush=True)
    return 0 if overall in ("PASS", "NOT_RUN") else 2


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except CampaignError as exc:
        print("STEP-8 CAMPAIGN ERROR: %s" % exc, file=sys.stderr)
        raise SystemExit(2)
