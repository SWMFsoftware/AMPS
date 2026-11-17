"""Reviewed, immutable I17 hold-out release and acceptance contracts.

This module never approves a threshold and never derives a tolerance from a
December residual.  A reviewer supplies a separate frozen JSON artifact after
the June numerical/training evidence exists.  Missing approval is a scientific
prerequisite miss; changed bytes or contradictory declarations are errors.
Temporary records used by test_i17_holdout.py exercise this protocol only and
are never observations or evidence that the production thermal model agrees
with Diviner.
"""
from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
from typing import Any


class HoldoutUnavailable(RuntimeError):
    """The investigator has not supplied a complete reviewed release."""


# These fields change the selected epoch, output ownership or observation
# transport, rather than the material/geometry/forcing model.  The complete
# reviewed hold-out input is ALSO byte-bound below, so exclusion here does not
# authorize a caller to change its independent spin-up interval silently.
_EPOCH_OR_TRANSPORT_KEYS = frozenset((
    "simulation_start_utc", "thermal_spinup_cycle_start_utc",
    "thermal_spinup_cycle_end_utc", "thermal_time_step_s",
    "diviner_eclipse_selection", "diviner_eclipse_freeze_file",
    "diviner_eclipse_authorization_file",
    "diviner_eclipse_history_file", "diviner_eclipse_selection_contract_file",
    "surface_cea_file", "surface_tecplot_file",
    "surface_inventory_checkpoint_file", "m2_formal_contract_file",
))


def sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _without_comment(line: str) -> str:
    """Remove an input ! comment while preserving a quoted literal !.

    This is a fail-closed identity projection for already native-validated
    flat Moon receipts, not an alternative production parser. Unsupported
    constructs raise an error instead of omitting physical input from a hash.
    """
    quote = None
    escaped = False
    for index, char in enumerate(line):
        if escaped:
            escaped = False
        elif char == "\\" and quote:
            escaped = True
        elif char in ("'", '"'):
            quote = None if quote == char else char if quote is None else quote
        elif char == "!" and quote is None:
            return line[:index].strip()
    if quote is not None:
        raise ValueError("I17 identity requires closed quoted input values")
    return line.strip()


def settings(path: Path) -> dict[str, str]:
    """Read the canonical flat assignment spelling for identity checks only."""
    result: dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = _without_comment(raw)
        if not line or line.startswith("#"):
            continue
        match = re.fullmatch(r"([A-Za-z][A-Za-z0-9_]*)\s*=\s*(.+)", line)
        if match is None:
            raise ValueError("I17 frozen identity requires flat one-line assignments")
        key, value = match.group(1).lower(), match.group(2).strip()
        if key in result:
            raise ValueError("duplicate I17 input identity key: " + key)
        result[key] = value
    return result


def unquote(value: str) -> str:
    if len(value) >= 2 and value[0] == value[-1] and value[0] in ("'", '"'):
        return value[1:-1]
    return value


def selection(path: Path) -> str:
    """Legacy receipts remain June-only; a new selection must be explicit."""
    value = unquote(settings(path).get("diviner_eclipse_selection", "training")).lower()
    if value not in ("training", "holdout"):
        raise ValueError("I17 selection must be training or holdout")
    return value


def freeze_path(template: Path) -> Path:
    value = unquote(settings(template).get("diviner_eclipse_freeze_file", "none"))
    if not value or value.lower() == "none":
        raise HoldoutUnavailable("I17 December requires a reviewed frozen threshold artifact")
    path = Path(value)
    return (path if path.is_absolute() else template.parent / path).resolve()


def physics_identity(path: Path) -> str:
    """Hash every supplied non-transport setting without numerical rounding."""
    declared = {key: value for key, value in settings(path).items()
                if key not in _EPOCH_OR_TRANSPORT_KEYS}
    for key in ("spice_path", "lola_image_file", "lola_label_file"):
        if key in declared:
            value = Path(unquote(declared[key]))
            declared[key] = str((value if value.is_absolute() else path.parent/value).resolve())
    return hashlib.sha256(json.dumps(declared, sort_keys=True,
        separators=(",", ":")).encode("utf-8")).hexdigest()


def _record(base: Path, record: Any, label: str) -> Path:
    if not isinstance(record, dict) or not isinstance(record.get("path"), str):
        raise ValueError("invalid frozen I17 file record: " + label)
    frozen = record.get("sha256", "")
    if not isinstance(frozen, str) or not re.fullmatch(r"[0-9a-f]{64}", frozen):
        raise ValueError("invalid frozen I17 SHA-256: " + label)
    path = Path(record["path"])
    path = (path if path.is_absolute() else base / path).resolve()
    if not path.is_file():
        raise HoldoutUnavailable("reviewed I17 evidence unavailable: " + str(path))
    if sha256(path) != frozen:
        raise ValueError("reviewed I17 evidence changed: " + label)
    return path


def _finite_limit(mapping: dict, key: str, *, fraction: bool = False) -> float:
    value = mapping.get(key)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("I17 reviewer must explicitly supply numerical " + key)
    value = float(value)
    if not math.isfinite(value) or value < 0 or (fraction and value > 1):
        raise ValueError("invalid frozen I17 acceptance limit: " + key)
    return value


def reviewed_git_identity(path: Path) -> dict[str, str]:
    """Require the separately committed frozen artifact promised by I17.

    This is read-only Git access. It neither creates a commit nor treats an
    untracked/edited candidate as approved. Stage copies are validated later
    by the hash recorded here, so private campaign directories need no Git.
    """
    # The review can remain in its original committed checkout while native
    # verification runs in a private copy. Its OWN enclosing Git repository
    # controls this immutable artifact; model identity is checked separately.
    owner = subprocess.run(["git", "-C", str(path.parent), "rev-parse", "--show-toplevel"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
    if owner.returncode:
        raise HoldoutUnavailable("I17 reviewed freeze must belong to a committed Git repository")
    repo = Path(owner.stdout.strip()).resolve()
    relative = path.resolve().relative_to(repo).as_posix()
    frozen = subprocess.run(["git", "show", "HEAD:"+relative], cwd=repo,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if frozen.returncode:
        raise HoldoutUnavailable("I17 reviewed threshold is not committed at HEAD")
    if hashlib.sha256(frozen.stdout).hexdigest() != sha256(path):
        raise ValueError("I17 committed threshold artifact has unreviewed local edits")
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip()
    blob = subprocess.check_output(["git", "rev-parse", "HEAD:"+relative], cwd=repo, text=True).strip()
    return {"repository": str(repo), "source_revision": revision,
            "blob_oid": blob, "relative_path": relative}


def validate_release(template: Path, executable: Path, *, history: Path,
                     contract: Path, processor: Path, campaign: Path) -> dict:
    """Validate approval/identity BEFORE staging any December science rows.

    Reading file hashes is metadata qualification, not residual inspection.
    The caller supplies the actual executable and exact canonical evidence.
    The JSON schema deliberately has no approval-writing helper: an aggregate
    runner cannot review its own candidate.  The campaign identifier binds one
    approved release to one planned output directory.  Repeated execution in
    that directory is refused separately by claim_release().
    """
    path = freeze_path(template)
    if not path.is_file():
        raise HoldoutUnavailable("reviewed I17 threshold file unavailable: " + str(path))
    doc = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(doc, dict):
        raise ValueError("I17 threshold release must be a JSON object")
    if doc.get("schema") != "srcMoon-I17-reviewed-holdout-v1":
        raise ValueError("unsupported I17 threshold-release schema")
    if doc.get("status") != "APPROVED_FROZEN":
        raise HoldoutUnavailable("I17 threshold candidate has not been approved/frozen")
    git_identity = reviewed_git_identity(path)
    if (doc.get("frozen_before_holdout_residuals") is not True or
            not isinstance(doc.get("approved_by"), str) or
            not doc["approved_by"].strip()):
        raise ValueError("I17 approval requires named reviewer and pre-holdout freeze")
    try:
        approved = dt.datetime.fromisoformat(doc["approved_utc"].replace("Z", "+00:00"))
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("I17 approval requires an explicit UTC timestamp") from error
    if approved.utcoffset() != dt.timedelta(0):
        raise ValueError("I17 approval timestamp must explicitly use UTC")
    if approved > dt.datetime.now(dt.timezone.utc):
        raise ValueError("I17 approval timestamp is in the future")
    if doc.get("holdout_campaign_id") != campaign.name:
        raise ValueError("I17 review does not authorize this planned campaign ID")
    if doc.get("initialization") != "independent-observation-free-periodic-spinup":
        raise ValueError("I17 December must use the approved independent spin-up method")
    for name, actual in (("processed_history_sha256", history),
                         ("selection_contract_sha256", contract),
                         ("processor_sha256", processor),
                         ("executable_sha256", executable)):
        if doc.get(name) != sha256(actual):
            raise ValueError("I17 frozen identity differs: " + name)
    from build_m2_qualification import production_source_identity
    native = production_source_identity()
    identity = hashlib.sha256(json.dumps(native, sort_keys=True,
        separators=(",", ":")).encode("utf-8")).hexdigest()
    if doc.get("production_source_identity_sha256") != identity:
        raise ValueError("I17 source/build provenance differs from reviewed June model")
    records = doc.get("training_evidence", {})
    if not isinstance(records, dict):
        raise ValueError("I17 training evidence must be a file-record object")
    verified = {key: _record(path.parent, records.get(key), key) for key in
                ("input", "candidate", "result", "manifest", "model_csv", "convergence_csv")}
    holdout_input = _record(path.parent, doc.get("reviewed_holdout_input"), "holdout input")
    if sha256(holdout_input) != sha256(template):
        raise ValueError("I17 hold-out receipt differs from the reviewed explicit input")
    if physics_identity(verified["input"]) != physics_identity(template):
        raise ValueError("June/December material, geometry or forcing settings differ")
    values = settings(template)
    for name in ("thermal_spinup_cycle_start_utc", "thermal_spinup_cycle_end_utc"):
        if name not in values or unquote(values[name]).lower() in ("", "none"):
            raise HoldoutUnavailable("December independently requires explicit " + name)
    if any(not unquote(values[name]).startswith("2011-") for name in
           ("thermal_spinup_cycle_start_utc", "thermal_spinup_cycle_end_utc")):
        raise ValueError("I17 spin-up interval must belong to the frozen 2011 case")
    # Do not silently reuse June's end state at the December epoch. CSPICE in
    # RunI17 validates full UTC ordering and that observations follow the end.
    if unquote(values["thermal_spinup_cycle_end_utc"])[:10] != "2011-12-10":
        raise ValueError("December spin-up must explicitly end on its frozen hold-out date")
    if unquote(values.get("simulation_start_utc", ""))[:10] != "2011-12-10":
        raise ValueError("December initialization must explicitly use its frozen hold-out date")
    candidate = json.loads(verified["candidate"].read_text())
    result = json.loads(verified["result"].read_text())
    manifest = json.loads(verified["manifest"].read_text())
    if not all(isinstance(value, dict) for value in (candidate, result, manifest)):
        raise ValueError("I17 training evidence must contain JSON objects")
    if not isinstance(manifest.get("executable"), dict):
        raise ValueError("June manifest lacks executable identity object")
    if (candidate.get("selection") != "training" or
            candidate.get("status") != "REVIEW_REQUIRED_NOT_FROZEN" or
            manifest.get("executable", {}).get("sha256") != sha256(executable)):
        raise ValueError("frozen I17 evidence is not the matching June training run")
    if not isinstance(result.get("metrics"), dict):
        raise ValueError("June result lacks numerical metrics object")
    convergence = result["metrics"].get("thermal_convergence", {})
    if not isinstance(convergence, dict):
        raise ValueError("June numerical gate is not an object")
    updates = convergence.get("maximum_update_changes_K")
    count = convergence.get("compared_bins")
    if (type(convergence.get("failed_bins")) is not int or convergence["failed_bins"] != 0 or
            type(count) is not int or count < 1 or
            not isinstance(updates, list) or len(updates) != 3 or
            any(isinstance(x, bool) or not isinstance(x, (int, float)) or
                not math.isfinite(x) or not 0 <= x < 50 for x in updates)):
        raise ValueError("June evidence does not establish the unchanged thermal numerical gate")
    with verified["model_csv"].open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows or any(row.get("selection", "training") != "training" for row in rows):
        raise ValueError("reviewed June model table is empty or contains another selection")
    if len(rows) != count or candidate.get("comparison_count") != count:
        raise ValueError("June candidate, numerical gate and model-row counts differ")
    limits = doc.get("acceptance", {})
    if not isinstance(limits, dict):
        raise ValueError("I17 acceptance must be an explicit metric-limit object")
    for key in ("weighted_rmse_k_max", "absolute_weighted_bias_k_max",
                "normalized_residual_rms_max", "coverage_sigma"):
        _finite_limit(limits, key)
    _finite_limit(limits, "minimum_coverage_fraction", fraction=True)
    # Phase/amplitude are deliberately defined on reviewed, frozen site/bin
    # groups; the scorer never searches for a group that happens to agree.
    if limits.get("phase_amplitude_operator") != "sampled-extrema-on-frozen-bin-groups":
        raise HoldoutUnavailable("reviewed I17 phase/amplitude operator is unavailable")
    for key in ("sampled_phase_error_s_max", "sampled_amplitude_error_k_max"):
        _finite_limit(limits, key)
    groups = doc.get("phase_amplitude_groups")
    if not isinstance(groups, list) or not groups:
        raise HoldoutUnavailable("reviewed pre-holdout phase/amplitude groups unavailable")
    seen = set()
    for group in groups:
        if not isinstance(group, dict):
            raise ValueError("I17 phase/amplitude group must be an object")
        ids = group.get("site_time_bin_ids", [])
        if (not isinstance(group.get("id"), str) or group["id"] in seen or
                not isinstance(ids, list) or len(ids) < 3 or
                not all(isinstance(value, str) and value for value in ids)):
            raise ValueError("invalid frozen I17 phase/amplitude group")
        if len(ids) != len(set(ids)):
            raise ValueError("duplicate frozen I17 phase/amplitude bin")
        seen.add(group["id"])
    # A tolerance freeze cannot by itself establish that omitted atmospheric
    # light, terrain-disc sampling or the grey observation operator are valid.
    # Preserve separately reviewed evidence of model-approximation scope.
    scope = doc.get("model_scope_qualification", {})
    if (not isinstance(scope, dict) or scope.get("status") != "READY" or
            not isinstance(scope.get("evidence"), list) or not scope["evidence"]):
        raise HoldoutUnavailable("I17 model-approximation/operator scope is not qualified")
    for index, record in enumerate(scope["evidence"]):
        _record(path.parent, record, "model scope " + str(index))
    return {"path": str(path), "sha256": sha256(path), "document": doc,
            "training_evidence": {key: str(value) for key, value in verified.items()},
            "physics_identity_sha256": physics_identity(template),
            "reviewed_git_identity": git_identity,
            "source_identity_sha256": identity}


def claim_release(directory: Path, release: dict) -> None:
    """Claim one approved hold-out execution before opening its science CSV.

    An interrupted/repeated claim cannot silently become a second inspection.
    Investigator review must authorize a new campaign/release if a restart is
    needed. Refinement children reuse the already staged selection; they do not
    claim additional independent observational experiments.
    """
    path = directory / "holdout_release_claim.json"
    try:
        with path.open("x", encoding="utf-8") as stream:
            json.dump({"schema": "srcMoon-I17-holdout-release-claim-v1",
                       "freeze_sha256": release["sha256"],
                       "freeze_path": release["path"],
                       "scope": "one approved hold-out campaign, including dt refinements"},
                      stream, indent=2)
            stream.write("\n")
    except FileExistsError as error:
        raise HoldoutUnavailable("I17 hold-out release already claimed; preserve prior evidence and request review before another execution") from error


_AUTHORIZATION_KEYS = frozenset(("schema", "selection", "review_status",
    "history_file", "history_sha256", "freeze_file", "freeze_sha256",
    "validation_file", "validation_sha256", "source_identity_sha256",
    "physics_identity_sha256", "review_token"))


def write_native_authorization(path: Path, *, history: Path, freeze: Path,
                               validation: Path, release: dict) -> None:
    """Write a small native guard ONLY after Python has validated the review.

    Native AMPS does not have a JSON/SHA library. The fixed key/value sidecar
    prevents accidental direct hold-out use with an arbitrary existing JSON
    file: C++ checks its schema/status/selection/paths/token structure before
    opening observations. Python remains the authority for actual byte hashes,
    reviewer approval, Git identity, source/model scope and one-time release.
    This is an audit/accident guard, not a signature/security mechanism.
    """
    fields = {"schema": "srcMoon-I17-native-authorization-v1", "selection": "holdout",
        "review_status": "APPROVED_FROZEN", "history_file": str(history.resolve()),
        "history_sha256": sha256(history), "freeze_file": str(freeze.resolve()),
        "freeze_sha256": sha256(freeze), "validation_file": str(validation.resolve()),
        "validation_sha256": sha256(validation),
        "source_identity_sha256": release["source_identity_sha256"],
        "physics_identity_sha256": release["physics_identity_sha256"],
        "review_token": release["sha256"]+"/"+release["reviewed_git_identity"]["blob_oid"]}
    if any("\n" in value or "\r" in value for value in fields.values()):
        raise ValueError("I17 native authorization cannot encode multiline paths/tokens")
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write("".join(key+"="+value+"\n" for key, value in fields.items()))


def authorization(application_input: Path) -> dict[str, str]:
    """Verify staged native guard and all referenced hashes before scoring."""
    values = settings(application_input)
    name = unquote(values.get("diviner_eclipse_authorization_file", "none"))
    if not name or name.lower() == "none":
        raise HoldoutUnavailable("I17 hold-out native authorization is unavailable")
    path = Path(name)
    path = (path if path.is_absolute() else application_input.parent/path).resolve()
    if not path.is_file():
        raise HoldoutUnavailable("I17 native authorization sidecar is missing")
    if path.parent != freeze_path(application_input).parent:
        raise ValueError("I17 authorization/review do not share the staged campaign context")
    fields = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        key, separator, value = line.partition("=")
        if not separator or key in fields or not value:
            raise ValueError("invalid/duplicate I17 native authorization field")
        fields[key] = value
    if (set(fields) != _AUTHORIZATION_KEYS or
            fields["schema"] != "srcMoon-I17-native-authorization-v1" or
            fields["selection"] != "holdout" or fields["review_status"] != "APPROVED_FROZEN"):
        raise ValueError("I17 native authorization schema/status/selection differs")
    expected_history = Path(unquote(values["diviner_eclipse_history_file"]))
    if not expected_history.is_absolute(): expected_history = application_input.parent/expected_history
    if fields["history_file"] != str(expected_history.resolve()):
        raise ValueError("I17 native authorization names another history")
    if fields["freeze_file"] != str(freeze_path(application_input)):
        raise ValueError("I17 native authorization names another review artifact")
    if fields["validation_file"] != str(path.parent/"holdout_release_validation.json"):
        raise ValueError("I17 native authorization names another validation receipt")
    for key in ("history", "freeze", "validation"):
        target = Path(fields[key+"_file"])
        if not target.is_absolute() or not target.is_file():
            raise HoldoutUnavailable("I17 authorized "+key+" file unavailable")
        if fields[key+"_sha256"] != sha256(target):
            raise ValueError("I17 authorized "+key+" bytes changed")
    return fields


def staged_release(application_input: Path) -> dict:
    """Load the byte-verified release staged BEFORE the linked execution."""
    auth = authorization(application_input)
    path = freeze_path(application_input)
    record = path.parent/"holdout_release_validation.json"
    if not record.is_file():
        raise HoldoutUnavailable("I17 pre-execution release validation receipt unavailable")
    release = json.loads(record.read_text())
    if release.get("sha256") != sha256(path):
        raise ValueError("staged I17 review changed after pre-execution validation")
    if release.get("document") != json.loads(path.read_text()):
        raise ValueError("staged I17 validation and reviewed document disagree")
    if release.get("physics_identity_sha256") != physics_identity(application_input):
        raise ValueError("I17 launched physics differs from pre-execution review")
    claim = json.loads((path.parent/"holdout_release_claim.json").read_text())
    if claim.get("freeze_sha256") != release["sha256"]:
        raise ValueError("I17 hold-out launch claim differs from frozen review")
    if (auth["source_identity_sha256"] != release["source_identity_sha256"] or
            auth["physics_identity_sha256"] != release["physics_identity_sha256"] or
            auth["review_token"] != release["sha256"]+"/"+release["reviewed_git_identity"]["blob_oid"]):
        raise ValueError("I17 native guard differs from reviewed source/physics/release identity")
    return release


def score_holdout(rows: list[dict[str, str]], release: dict) -> tuple[str, dict, str]:
    """Apply only predeclared limits and groups; never tune or interpolate."""
    doc = release["document"]
    if (doc.get("status") != "APPROVED_FROZEN" or
            doc.get("frozen_before_holdout_residuals") is not True):
        raise HoldoutUnavailable("I17 scorer requires an approved pre-holdout release")
    limits = doc["acceptance"]
    if not rows or any(row.get("selection") != "holdout" for row in rows):
        raise ValueError("I17 formal score requires nonempty hold-out-only model rows")
    identities = [row["site_time_bin_id"] for row in rows]
    if len(set(identities)) != len(identities):
        raise ValueError("I17 hold-out model has duplicate bin identities")
    residuals, normalized, weights = [], [], []
    for row in rows:
        reference = float(row["observed_bolometric_brightness_temperature_K"])
        model = float(row["model_bolometric_brightness_temperature_K"])
        sigma = float(row["observed_uncertainty_K"])
        if not all(math.isfinite(value) and value > 0 for value in (reference, model, sigma)):
            raise ValueError("invalid I17 hold-out temperature or documented uncertainty")
        residuals.append(model-reference)
        normalized.append((model-reference)/sigma)
        weights.append(1/sigma**2)
    total_weight = sum(weights)
    metrics = {
        "comparison_count": len(rows), "selection": "holdout",
        "weighted_rmse_k": math.sqrt(sum(w*r*r for w, r in zip(weights, residuals))/total_weight),
        "weighted_bias_k": sum(w*r for w, r in zip(weights, residuals))/total_weight,
        "normalized_residual_rms": math.sqrt(sum(r*r for r in normalized)/len(rows)),
        "uncertainty_coverage_fraction": sum(abs(r) <= limits["coverage_sigma"] for r in normalized)/len(rows),
        "threshold_freeze_sha256": release["sha256"],
        "phase_amplitude": [],
    }
    by_id = dict(zip(identities, rows))
    for group in doc["phase_amplitude_groups"]:
        if any(key not in by_id for key in group["site_time_bin_ids"]):
            raise ValueError("hold-out omits a frozen phase/amplitude bin; no gap filling allowed")
        subset = [by_id[key] for key in group["site_time_bin_ids"]]
        if len({(r["latitude_deg"], r["longitude_east_deg"]) for r in subset}) != 1:
            raise ValueError("phase/amplitude group combines distinct native model sites")
        times = [float(r["ephemeris_time_tdb_s"]) for r in subset]
        if len(set(times)) != len(times) or not all(math.isfinite(value) for value in times):
            raise ValueError("phase/amplitude requires distinct finite native sample epochs")
        observed = [float(r["observed_bolometric_brightness_temperature_K"]) for r in subset]
        predicted = [float(r["model_bolometric_brightness_temperature_K"]) for r in subset]
        if observed.count(min(observed)) != 1:
            return "SKIPPED", metrics, "frozen observed sampled minimum is ambiguous; phase is not identifiable"
        if predicted.count(min(predicted)) != 1:
            return "FAIL", metrics, "model sampled minimum is ambiguous in a declared phase group"
        # These are sampled extrema on exactly the declared bins. Real gaps
        # remain gaps: the code claims neither a continuously resolved minimum
        # nor a reconstructed eclipse curve between native samples.
        metrics["phase_amplitude"].append({"id": group["id"],
            "sampled_phase_error_s": abs(times[predicted.index(min(predicted))]-times[observed.index(min(observed))]),
            "sampled_amplitude_error_k": abs((max(predicted)-min(predicted))-(max(observed)-min(observed))),
            "sample_count": len(subset)})
    checks = {
        "weighted_rmse": metrics["weighted_rmse_k"] <= limits["weighted_rmse_k_max"],
        "weighted_bias": abs(metrics["weighted_bias_k"]) <= limits["absolute_weighted_bias_k_max"],
        "normalized_rms": metrics["normalized_residual_rms"] <= limits["normalized_residual_rms_max"],
        "uncertainty_coverage": metrics["uncertainty_coverage_fraction"] >= limits["minimum_coverage_fraction"],
        "sampled_phase": all(r["sampled_phase_error_s"] <= limits["sampled_phase_error_s_max"] for r in metrics["phase_amplitude"]),
        "sampled_amplitude": all(r["sampled_amplitude_error_k"] <= limits["sampled_amplitude_error_k_max"] for r in metrics["phase_amplitude"]),
    }
    metrics["acceptance_checks"] = checks
    return ("PASS" if all(checks.values()) else "FAIL", metrics,
            "reviewed I17 hold-out criteria evaluated; no parameter or threshold fitting")
