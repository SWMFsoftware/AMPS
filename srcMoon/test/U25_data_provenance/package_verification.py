#!/usr/bin/env python3
"""U25: validation-data package provenance, immutability and reproduction.

U25 is the shared data-package preflight of every data-dependent srcMoon test.
It judges only the *package*: whether the selected Dxx package is complete,
unaltered, reproducible from its immutable raw files by its own committed
processing script, QA-ready, and (for observation datasets) has a frozen
selection.  It never compares AMPS with observations, so a PASS is data-package
QA only and never observational validation (acceptance.json, status_semantics).

Dataset selection
-----------------
The contract says the dataset list is "declared by the requesting test".  The
requesting tests are the U contracts that list U25 among their
``predecessor_gates``; their ``required_data`` entries name the packages
(e.g. U15 -> "D04", U26 -> "D07-D09/D12-D15").  ``declared_datasets()`` parses
those committed contracts at run time, so this module holds no second copy of
the mapping that could silently drift.  An untargeted run evaluates the union
and reports one aggregate status (user decision, 2026-11-23): PASS only when
every declared dataset passes, otherwise the most serious of ERROR, FAIL and
SKIPPED.  ``--dataset`` / ``--for-test`` select a specific target.

Order of the gates for one dataset (each gate stops the evaluation):

1. registry entry and package directory present          -> else SKIPPED
2. QA readiness value in acceptance.accepted_readiness_values
   (``READY``/``PASS``); corrupt or schema-inconsistent QA -> ERROR,
   any other native value                                  -> SKIPPED
3. observation datasets: a registered selection file whose status is one of
   that package's registered frozen values; none are registered while no
   package has frozen its split, so this gate fails closed -> SKIPPED
4. integrity: required artifacts (missing -> SKIPPED), hash manifests,
   manifest completeness of raw/ and processed/, raw inventory, provenance
   output/script hashes, no symlinks                        -> ERROR on mismatch;
   provenance declaring a prohibited method (manual edit, synthetic data,
   digitization, OCR) as used                              -> FAIL
5. reproduction: the package's documented processing command is run on a
   private copy (never on the canonical package), with network proxies
   pointed at a discard port.  Nonzero exit/timeout        -> ERROR;
   canonical package changed during the run                -> ERROR;
   any reproduced file different                           -> FAIL.

A gate that stops the evaluation leaves the later gates unevaluated; the
result says so explicitly (``integrity_evaluated`` / ``reproduction``).

Timestamp rule (user decision, 2026-11-23)
------------------------------------------
Regenerated provenance/QA files carry a new generation time.  A reproduced file
may differ from the canonical one only on lines that carry a field named in
acceptance.json ``timestamp_normalization`` *and* registered for that file in
packages.json, and only when both values are ISO-8601 UTC strings.  For JSON
files the two documents must additionally be equal after deleting exactly the
registered key paths, so an identically named field elsewhere cannot hide a
difference.  The package's own SHA256SUMS.txt then necessarily differs in the
hash entries of exactly those files; every other entry, and the listed path
set, must be identical.  Everything else must be byte-identical.

Metadata completeness (user decision, 2026-11-23)
-------------------------------------------------
The roadmap's C6 minimum provenance schema is reported field by field
(present/absent) without mapping differently named equivalents and without
gating: no current package uses it literally.  Absence of a declared manual
edit is established directly by byte-identical reproduction instead.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple

HERE = Path(__file__).resolve().parent
TEST_ROOT = HERE.parent
REPO_ROOT = TEST_ROOT.parents[1]
ACCEPTANCE_PATH = HERE / "reference" / "acceptance.json"
REGISTRY_PATH = HERE / "reference" / "packages.json"
DEFAULT_DATA_ROOT = Path("/data/vtenishe/moon_validation_data")
DATA_ROOT_ENVIRONMENT = "MOON_VALIDATION_DATA_ROOT"

# Most serious first.  The aggregate and every per-dataset result use the
# AGENTS.md status semantics: ERROR (infrastructure/data corruption), FAIL
# (evaluated criterion not met), SKIPPED (data/capability unavailable), PASS.
STATUS_ORDER = ("ERROR", "FAIL", "SKIPPED", "PASS")

# "D04", "D10 global ...", or a range such as "D07-D09".  Free text such as
# "target-specific Dxx package" deliberately does not match.
DATASET_TOKEN = re.compile(r"\bD(\d{2})(?:\s*-\s*D(\d{2}))?\b")

# Generation times as written by the packages (datetime.isoformat() + "Z").
ISO_UTC = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d+)?Z$")

# The C6 minimum provenance schema of development_and_validation.md, reported
# for metadata completeness only.
C6_FIELDS = (
    "dataset", "mission", "instrument", "source_product_ids",
    "source_product_versions", "source_urls_or_dois", "raw_files",
    "processing_script", "processing_script_sha256", "transformations",
    "quality_filters", "interpolation", "smoothing", "manual_edits",
    "synthetic_data_used", "output_files", "notes",
)

# Lines of a sha256sum manifest: 64 lowercase hex digits, two spaces, path.
MANIFEST_LINE = re.compile(r"^([0-9a-f]{64})  (\S.*)$")

# Proxy variables honoured by urllib/requests; pointing them at the discard
# port makes an undeclared network access fail instead of silently fetching.
OFFLINE_PROXY = "http://127.0.0.1:9"


class Stop(Exception):
    """End the evaluation of one dataset with a status and a reason."""

    def __init__(self, status: str, reason: str, **details: Any) -> None:
        super().__init__(reason)
        self.status = status
        self.reason = reason
        self.details = details


# ---------------------------------------------------------------------------
# Small helpers


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def sha256_file(path: Path) -> str:
    """SHA-256 of a file, streamed so that large raw products fit in memory."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path, what: str) -> Any:
    """Read package JSON; a corrupt file is data corruption, i.e. ERROR."""
    try:
        with path.open(encoding="utf-8") as stream:
            return json.load(stream)
    except (OSError, ValueError) as error:
        raise Stop("ERROR", f"cannot read {what} {path}: {error}")


def lookup(document: Any, key_path: List[str], what: str) -> Any:
    """Follow a registered key path; a missing key is a schema mismatch."""
    value = document
    for key in key_path:
        if not isinstance(value, dict) or key not in value:
            raise Stop("ERROR", f"{what}: registered key path {key_path} is absent")
        value = value[key]
    return value


def ignored(relative: Path, names: List[str]) -> bool:
    return any(part in names for part in relative.parts)


def package_files(root: Path, ignored_names: List[str]) -> Dict[str, Path]:
    """All regular files of a package (POSIX relative path -> path).

    Symlinks are refused: a link could make the private reproduction copy
    write through to the canonical package, and it hides the real bytes from
    the hash manifests.
    """
    files: Dict[str, Path] = {}
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root)
        if ignored(relative, ignored_names):
            continue
        if path.is_symlink():
            raise Stop("ERROR", f"symbolic link in package: {relative.as_posix()}")
        if path.is_file():
            files[relative.as_posix()] = path
    return files


def inside(root: Path, relative: str, what: str) -> Path:
    """Resolve a package-relative path and refuse anything outside the root."""
    candidate = Path(relative)
    if candidate.is_absolute() or ".." in candidate.parts:
        raise Stop("ERROR", f"{what}: path {relative!r} is not package-relative")
    return root / candidate


def worst(statuses: List[str]) -> str:
    for status in STATUS_ORDER:
        if status in statuses:
            return status
    return "PASS"


# ---------------------------------------------------------------------------
# Contracts and registry


def load_contract() -> Dict[str, Any]:
    with ACCEPTANCE_PATH.open(encoding="utf-8") as stream:
        return json.load(stream)


def load_registry() -> Dict[str, Any]:
    with REGISTRY_PATH.open(encoding="utf-8") as stream:
        return json.load(stream)


def expand_dataset_tokens(texts: List[str]) -> List[str]:
    """Return the Dxx identifiers named in ``required_data`` strings."""
    found: List[str] = []
    for text in texts:
        for match in DATASET_TOKEN.finditer(str(text)):
            first = int(match.group(1))
            last = int(match.group(2)) if match.group(2) else first
            if last < first:
                raise ValueError(f"descending dataset range in {text!r}")
            for number in range(first, last + 1):
                name = f"D{number:02d}"
                if name not in found:
                    found.append(name)
    return found


def declared_datasets(test_root: Path = TEST_ROOT) -> Dict[str, List[str]]:
    """Map each U test that lists U25 as predecessor to its declared Dxx.

    The mapping is read from the committed acceptance contracts every run.  A
    requesting test without any parseable Dxx would make the untargeted U25
    silently ignore it, so that is an infrastructure ERROR.
    """
    result: Dict[str, List[str]] = {}
    for contract_path in sorted(test_root.glob("U[0-9][0-9]_*/reference/acceptance.json")):
        with contract_path.open(encoding="utf-8") as stream:
            contract = json.load(stream)
        if contract.get("id") == "U25" or "U25" not in contract.get("predecessor_gates", []):
            continue
        datasets = expand_dataset_tokens(contract.get("required_data", []))
        if not datasets:
            raise Stop("ERROR", f"{contract.get('id')} lists U25 as predecessor but declares no Dxx dataset")
        result[str(contract["id"])] = datasets
    return result


def check_registry(contract: Dict[str, Any], registry: Dict[str, Any]) -> None:
    """Registered timestamp fields may only use the contract's frozen names."""
    frozen = contract["timestamp_normalization"]
    allowed_keys = set(frozen["json_keys"])
    allowed_prefixes = set(frozen["text_line_prefixes"])
    for name, spec in registry["packages"].items():
        fields = (spec.get("reproduction") or {}).get("timestamp_fields", {})
        for file_name, rule in fields.items():
            for key_path in rule.get("json_keys", []):
                if not key_path or key_path[-1] not in allowed_keys:
                    raise Stop("ERROR", f"packages.json {name}/{file_name}: timestamp key {key_path} is not frozen in acceptance.json")
            for prefix in rule.get("line_prefixes", []):
                if prefix not in allowed_prefixes:
                    raise Stop("ERROR", f"packages.json {name}/{file_name}: line prefix {prefix!r} is not frozen in acceptance.json")


# ---------------------------------------------------------------------------
# Gates 2 and 3: readiness and frozen selections


def check_readiness(entry: Dict[str, Any], package: Path, accepted: List[str]) -> Dict[str, Any]:
    rule = entry["readiness"]
    qa_path = inside(package, rule["file"], "readiness file")
    if not qa_path.is_file():
        raise Stop("SKIPPED", f"{rule['file']} is absent")
    value = lookup(read_json(qa_path, "QA report"), rule["key"], rule["file"])
    record = {"file": rule["file"], "key": rule["key"], "value": value,
              "accepted_values": accepted}
    if value not in accepted:
        raise Stop("SKIPPED", f"{rule['file']} {'.'.join(rule['key'])}={value}",
                   readiness=record)
    return record


def check_selection(entry: Dict[str, Any], package: Path) -> List[Dict[str, Any]]:
    """Observation datasets need a frozen calibration/hold-out or selection.

    Roadmap C8 requires a committed split before final scoring.  The frozen
    state is accepted only through values registered per package after
    review; no package has frozen its split yet, so none are registered and
    this gate stays closed rather than guessing what "frozen" looks like.
    """
    if not entry.get("observation_dataset"):
        return []
    files = entry.get("selection_files", [])
    if not files:
        raise Stop("SKIPPED", "observation dataset without a calibration/hold-out or selection file")
    records = []
    for rule in files:
        path = inside(package, rule["file"], "selection file")
        if not path.is_file():
            raise Stop("SKIPPED", f"{rule['file']} is absent", selections=records)
        value = lookup(read_json(path, "selection file"), rule["status_key"], rule["file"])
        frozen = rule.get("frozen_status_values", [])
        records.append({"file": rule["file"], "key": rule["status_key"], "value": value,
                        "registered_frozen_values": frozen,
                        "sha256": sha256_file(path)})
        if value not in frozen:
            raise Stop("SKIPPED", f"{rule['file']} {'.'.join(rule['status_key'])}={value} (selection not frozen)",
                       selections=records)
    return records


# ---------------------------------------------------------------------------
# Gate 4: integrity


def parse_manifest(path: Path) -> List[Tuple[str, str]]:
    entries = []
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        match = MANIFEST_LINE.match(line)
        if not match:
            raise Stop("ERROR", f"{path.name} line {number} is not '<sha256>  <path>'")
        entries.append((match.group(1), match.group(2)))
    return entries


def verify_integrity(spec: Dict[str, Any], package: Path) -> Dict[str, Any]:
    """Recompute every recorded hash and check raw immutability.

    Returns the verified hash table for verified_hashes.json.  Any recorded
    hash that disagrees with the bytes on disk is data corruption (ERROR).
    """
    ignored_names = spec.get("ignored_names", [])
    missing = [name for name in spec["required_artifacts"] if not (package / name).is_file()]
    if missing:
        raise Stop("SKIPPED", "required package artifacts absent: " + ", ".join(missing),
                   missing_artifacts=missing)
    files = package_files(package, ignored_names)
    hashes = {relative: sha256_file(path) for relative, path in files.items()}
    checks: Dict[str, Any] = {}

    # Hash manifests (SHA256SUMS.txt at the package root and/or raw/).
    listed: Dict[str, str] = {}
    for manifest in spec["hash_manifests"]:
        manifest_path = inside(package, manifest["file"], "hash manifest")
        base = inside(package, manifest.get("base", "."), "manifest base")
        count = 0
        for expected, relative in parse_manifest(manifest_path):
            target = inside(base, relative, manifest["file"])
            key = target.relative_to(package).as_posix()
            if key not in hashes:
                raise Stop("ERROR", f"{manifest['file']} lists a missing file: {relative}")
            if hashes[key] != expected:
                raise Stop("ERROR", f"{manifest['file']}: hash mismatch for {key}",
                           expected_sha256=expected, actual_sha256=hashes[key])
            listed[key] = expected
            count += 1
        checks[manifest["file"]] = {"entries": count, "status": "PASS"}

    # Completeness: nothing in raw/ or processed/ escapes the manifests.
    for directory in spec.get("complete_manifest_directories", []):
        prefix = directory.rstrip("/") + "/"
        unlisted = sorted(key for key in hashes if key.startswith(prefix) and key not in listed)
        if unlisted:
            raise Stop("ERROR", f"files under {directory}/ not in any hash manifest: " + ", ".join(unlisted[:5]),
                       unlisted_files=unlisted)
        checks[f"manifest_completeness:{directory}"] = "PASS"

    # Raw inventory: per-file hash, size, immutability flag, and exact set.
    rule = spec["raw_inventory"]
    inventory = read_json(inside(package, rule["file"], "raw inventory"), "raw inventory")
    entries = lookup(inventory, rule["list_key"], rule["file"])
    inventoried = set()
    for item in entries:
        relative = str(item[rule["path_key"]])
        key = inside(package, relative, rule["file"]).relative_to(package).as_posix()
        if key not in hashes:
            raise Stop("ERROR", f"{rule['file']} lists a missing raw file: {relative}")
        if hashes[key] != item[rule["sha256_key"]]:
            raise Stop("ERROR", f"raw file changed since inventory: {key}",
                       expected_sha256=item[rule["sha256_key"]], actual_sha256=hashes[key])
        if rule.get("bytes_key") and files[key].stat().st_size != item[rule["bytes_key"]]:
            raise Stop("ERROR", f"raw file size differs from inventory: {key}")
        if rule.get("immutable_key") and item.get(rule["immutable_key"]) is not True:
            raise Stop("ERROR", f"raw inventory does not declare {key} immutable")
        inventoried.add(key)
    raw_prefix = spec["raw_directory"].rstrip("/") + "/"
    on_disk = {key for key in hashes if key.startswith(raw_prefix)}
    if on_disk != inventoried:
        raise Stop("ERROR", "raw/ differs from the raw inventory",
                   not_inventoried=sorted(on_disk - inventoried),
                   inventoried_but_absent=sorted(inventoried - on_disk))
    checks["raw_inventory"] = {"entries": len(inventoried), "status": "PASS"}

    # Provenance: output and script hashes, self-references, prohibited flags.
    prov_rule = spec["provenance"]
    provenance = read_json(inside(package, prov_rule["file"], "provenance"), "provenance")
    for group in ("outputs", "scripts"):
        group_rule = prov_rule[group]
        recorded = set()
        for item in lookup(provenance, group_rule["list_key"], prov_rule["file"]):
            relative = str(item[group_rule["path_key"]])
            key = inside(package, relative, prov_rule["file"]).relative_to(package).as_posix()
            if key not in hashes:
                raise Stop("ERROR", f"provenance {group} entry is missing on disk: {relative}")
            if hashes[key] != item[group_rule["sha256_key"]]:
                raise Stop("ERROR", f"provenance {group} hash mismatch: {key}")
            if group_rule.get("bytes_key") and files[key].stat().st_size != item[group_rule["bytes_key"]]:
                raise Stop("ERROR", f"provenance {group} size mismatch: {key}")
            recorded.add(key)
        cover = group_rule.get("covers_directory")
        if cover:
            prefix = cover.rstrip("/") + "/"
            present = {key for key in hashes if key.startswith(prefix)}
            pattern = group_rule.get("covers_suffix")
            if pattern:
                present = {key for key in present if key.endswith(pattern)}
            if present != recorded:
                raise Stop("ERROR", f"provenance {group} do not cover {cover}/ exactly",
                           unrecorded=sorted(present - recorded), recorded_but_absent=sorted(recorded - present))
        checks[f"provenance_{group}"] = {"entries": len(recorded), "status": "PASS"}
    for reference in prov_rule.get("self_hashes", []):
        relative = lookup(provenance, reference["path_key"], prov_rule["file"])
        expected = lookup(provenance, reference["sha256_key"], prov_rule["file"])
        key = inside(package, str(relative), prov_rule["file"]).relative_to(package).as_posix()
        if hashes.get(key) != expected:
            raise Stop("ERROR", f"provenance reference hash mismatch: {key}")
    declared_true = []
    for key_path in prov_rule.get("must_be_false", []):
        value = lookup(provenance, key_path, prov_rule["file"])
        if value is not False:
            declared_true.append({"key": key_path, "value": value})
    for key in ("manual_edits", "synthetic_data_used"):
        if key in provenance and provenance[key] is not False:
            declared_true.append({"key": [key], "value": provenance[key]})
    if declared_true:
        raise Stop("FAIL", "provenance declares a prohibited method as used: "
                   + ", ".join(".".join(item["key"]) for item in declared_true),
                   declared_prohibited_methods=declared_true)
    checks["prohibited_methods"] = "PASS"

    completeness = {field: (field in provenance) for field in C6_FIELDS}
    return {"files": hashes, "checks": checks,
            "metadata_completeness_c6": {
                "gating": False,
                "present": [field for field, present in completeness.items() if present],
                "absent": [field for field, present in completeness.items() if not present],
                "note": "literal C6 field names only; differently named equivalents are not mapped"}}


# ---------------------------------------------------------------------------
# Gate 5: reproduction


def timestamp_only_lines(canonical: str, reproduced: str, json_keys: List[str],
                         prefixes: List[str]) -> Optional[str]:
    """Return None when the texts differ only on registered timestamp lines.

    Otherwise return a short description of the first disallowed difference.
    Line counts must agree; each differing line pair must carry the same
    registered field (same indentation and trailing comma for JSON) with an
    ISO-8601 UTC value on both sides.
    """
    left = canonical.splitlines(keepends=True)
    right = reproduced.splitlines(keepends=True)
    if len(left) != len(right):
        return f"line count {len(left)} -> {len(right)}"
    key_pattern = re.compile(r'^(\s*)"([^"]+)": "([^"]*)"(,?)(\r?\n)?$')
    for number, (a, b) in enumerate(zip(left, right), 1):
        if a == b:
            continue
        ma, mb = key_pattern.match(a), key_pattern.match(b)
        if (ma and mb and ma.group(2) in json_keys and ma.group(2) == mb.group(2)
                and ma.group(1) == mb.group(1) and ma.group(4) == mb.group(4)
                and ma.group(5) == mb.group(5)
                and ISO_UTC.match(ma.group(3)) and ISO_UTC.match(mb.group(3))):
            continue
        matched = False
        for prefix in prefixes:
            if a.startswith(prefix) and b.startswith(prefix):
                va, vb = a[len(prefix):].rstrip("\r\n"), b[len(prefix):].rstrip("\r\n")
                if (ISO_UTC.match(va) and ISO_UTC.match(vb)
                        and a[len(prefix) + len(va):] == b[len(prefix) + len(vb):]):
                    matched = True
                    break
        if not matched:
            return f"line {number} differs"
    return None


def delete_key_path(document: Any, key_path: List[str]) -> None:
    value = document
    for key in key_path[:-1]:
        if not isinstance(value, dict) or key not in value:
            return
        value = value[key]
    if isinstance(value, dict):
        value.pop(key_path[-1], None)


def compare_normalized(canonical: Path, reproduced: Path, rule: Dict[str, Any]) -> Optional[str]:
    """Apply the frozen timestamp rule to one registered metadata file."""
    try:
        left = canonical.read_text(encoding="utf-8")
        right = reproduced.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as error:
        return f"not comparable as text: {error}"
    key_paths = rule.get("json_keys", [])
    problem = timestamp_only_lines(left, right, [path[-1] for path in key_paths],
                                   rule.get("line_prefixes", []))
    if problem:
        return problem
    if key_paths:
        try:
            a, b = json.loads(left), json.loads(right)
        except ValueError as error:
            return f"not valid JSON: {error}"
        for key_path in key_paths:
            delete_key_path(a, key_path)
            delete_key_path(b, key_path)
        if a != b:
            return "differs outside the registered timestamp key paths"
    return None


def reproduce(spec: Dict[str, Any], package: Path, work: Path,
              canonical_hashes: Dict[str, str]) -> Dict[str, Any]:
    """Rerun the package's documented processing on a private copy."""
    rule = spec.get("reproduction")
    if not rule:
        raise Stop("SKIPPED", "no package-specific reproduction command is registered",
                   reproduction={"evaluated": False})
    interpreter = Path(rule["interpreter"])
    if not interpreter.is_absolute():
        interpreter = package / interpreter
    if not interpreter.exists():
        raise Stop("SKIPPED", f"documented interpreter is unavailable: {interpreter}",
                   reproduction={"evaluated": False})
    ignored_names = spec.get("ignored_names", [])
    copy = work / "copy"
    if copy.exists():
        raise Stop("ERROR", f"reproduction directory already exists: {copy}")
    work.mkdir(parents=True, exist_ok=True)
    shutil.copytree(package, copy, symlinks=False,
                    ignore=shutil.ignore_patterns(*rule.get("copy_exclude", [])))
    copied = package_files(copy, ignored_names)
    if set(copied) != set(canonical_hashes):
        raise Stop("ERROR", "private copy does not contain the canonical file set")

    command = [part.replace("{python}", str(interpreter)).replace("{copy}", str(copy))
               for part in rule["command"]]
    environment = os.environ.copy()
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    if not rule.get("network", False):
        for name in ("http_proxy", "https_proxy", "HTTP_PROXY", "HTTPS_PROXY", "ALL_PROXY", "all_proxy"):
            environment[name] = OFFLINE_PROXY
        environment["no_proxy"] = environment["NO_PROXY"] = ""
    version = subprocess.run([str(interpreter), "--version"], capture_output=True, text=True, check=False)
    report: Dict[str, Any] = {
        "evaluated": True, "command": command, "cwd": str(work),
        "interpreter": str(interpreter),
        "interpreter_version": (version.stdout or version.stderr).strip(),
        "network_blocked_by_proxy": not rule.get("network", False),
        "documented_in": rule.get("documented_in"),
        "log": str(work / "reproduction.log"),
    }
    started = time.monotonic()
    with (work / "reproduction.log").open("w", encoding="utf-8") as log:
        log.write("COMMAND: " + " ".join(command) + "\n")
        log.flush()
        try:
            completed = subprocess.run(command, cwd=work, env=environment, stdout=log,
                                       stderr=subprocess.STDOUT, check=False,
                                       timeout=rule["timeout_seconds"])
            report["returncode"] = completed.returncode
        except subprocess.TimeoutExpired:
            report["returncode"] = None
            report["elapsed_seconds"] = time.monotonic() - started
            raise Stop("ERROR", f"reproduction exceeded {rule['timeout_seconds']} s", reproduction=report)
    report["elapsed_seconds"] = time.monotonic() - started

    # The canonical package must be byte-for-byte untouched by the run.
    after = {key: sha256_file(path) for key, path in package_files(package, ignored_names).items()}
    if after != canonical_hashes:
        changed = sorted(key for key in set(after) | set(canonical_hashes)
                         if after.get(key) != canonical_hashes.get(key))
        report["canonical_changed"] = changed
        raise Stop("ERROR", "canonical package changed during reproduction: " + ", ".join(changed[:5]),
                   reproduction=report)
    if report["returncode"] != 0:
        raise Stop("ERROR", f"reproduction command exited {report['returncode']}", reproduction=report)

    # Compare every file of the canonical package with its reproduction.
    reproduced = package_files(copy, ignored_names)
    timestamp_fields = rule.get("timestamp_fields", {})
    comparison: Dict[str, str] = {}
    differences: List[Dict[str, str]] = []
    normalized_ok = set()
    for key in sorted(set(canonical_hashes) | set(reproduced)):
        if key not in reproduced:
            differences.append({"path": key, "difference": "absent after reproduction"})
            continue
        if key not in canonical_hashes:
            differences.append({"path": key, "difference": "created by reproduction, absent from package"})
            continue
        if sha256_file(reproduced[key]) == canonical_hashes[key]:
            comparison[key] = "identical"
            continue
        if key in timestamp_fields:
            problem = compare_normalized(package / key, reproduced[key], timestamp_fields[key])
            if problem is None:
                comparison[key] = "identical after registered timestamp fields"
                normalized_ok.add(key)
                continue
            differences.append({"path": key, "difference": problem})
            continue
        comparison[key] = "pending manifest rule"
    # Hash manifests may differ only in the entries of normalized files.
    manifest_files = {inside(package, item["file"], "hash manifest").relative_to(package).as_posix(): item
                      for item in spec["hash_manifests"]}
    for key, state in list(comparison.items()):
        if state != "pending manifest rule":
            continue
        if key not in manifest_files:
            differences.append({"path": key, "difference": "bytes differ"})
            del comparison[key]
            continue
        base_rel = manifest_files[key].get("base", ".")
        left = {inside(package / base_rel, path, key).relative_to(package).as_posix(): digest
                for digest, path in parse_manifest(package / key)}
        right = {inside(copy / base_rel, path, key).relative_to(copy).as_posix(): digest
                 for digest, path in parse_manifest(reproduced[key])}
        bad = sorted(path for path in set(left) | set(right)
                     if left.get(path) != right.get(path) and not (
                         path in left and path in right and path in normalized_ok))
        if bad:
            differences.append({"path": key, "difference": "manifest entries differ: " + ", ".join(bad[:5])})
            del comparison[key]
        else:
            comparison[key] = "identical except hashes of timestamp-normalized files"
    report["files_compared"] = len(set(canonical_hashes) | set(reproduced))
    report["comparison_counts"] = {state: list(comparison.values()).count(state)
                                   for state in sorted(set(comparison.values()))}
    report["normalized_files"] = sorted(normalized_ok)
    differences.sort(key=lambda item: item["path"])
    report["differences"] = differences
    if differences:
        raise Stop("FAIL", f"{len(differences)} reproduced file(s) differ: "
                   + ", ".join(item["path"] for item in differences[:5]), reproduction=report)
    return report


# ---------------------------------------------------------------------------
# One dataset, many datasets


def evaluate_package(spec: Dict[str, Any], package: Path, work: Path) -> Dict[str, Any]:
    """Integrity then reproduction of one package (shared by its datasets)."""
    result: Dict[str, Any] = {"package": str(package), "integrity_evaluated": False,
                              "reproduction": {"evaluated": False}}
    try:
        integrity = verify_integrity(spec, package)
        result["integrity_evaluated"] = True
        result["integrity"] = {key: value for key, value in integrity.items() if key != "files"}
        result["verified_hashes"] = integrity["files"]
        result["reproduction"] = reproduce(spec, package, work, integrity["files"])
        result.update(status="PASS", reason="package complete, unaltered and reproducible")
    except Stop as stop:
        result.update(status=stop.status, reason=stop.reason)
        result.update(stop.details)
    return result


def evaluate_dataset(dataset: str, registry: Dict[str, Any], contract: Dict[str, Any],
                     data_root: Path, work_root: Path, cache: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    entry = registry["datasets"].get(dataset)
    result: Dict[str, Any] = {"dataset": dataset}
    try:
        if entry is None:
            raise Stop("SKIPPED", "no U25 recipe is registered for this dataset")
        package = data_root / entry["package_directory"]
        result["package"] = str(package)
        result["description"] = entry.get("description")
        if not package.is_dir():
            raise Stop("SKIPPED", f"package directory is absent: {package}")
        result["readiness"] = check_readiness(entry, package, contract["accepted_readiness_values"])
        result["selections"] = check_selection(entry, package)
        spec_name = entry.get("verification")
        if not spec_name:
            raise Stop("SKIPPED", entry.get("verification_unregistered_reason",
                                            "no integrity/reproduction recipe is registered"))
        spec = registry["packages"][spec_name]
        if spec_name not in cache:
            cache[spec_name] = evaluate_package(spec, package, work_root / spec_name)
        shared = cache[spec_name]
        result["package_evaluation"] = spec_name
        result.update({key: value for key, value in shared.items() if key not in ("package",)})
    except Stop as stop:
        result.update(status=stop.status, reason=stop.reason)
        result.update(stop.details)
    result.setdefault("integrity_evaluated", False)
    result.setdefault("reproduction", {"evaluated": False})
    return result


def run_controls(output: Path) -> Dict[str, Any]:
    """Run the verifier's own planted-defect controls (fixtures, not data)."""
    runner = HERE / "test_package_verification.py"
    command = [sys.executable, "-B", str(runner)]
    process = subprocess.run(command, cwd=REPO_ROOT, capture_output=True, text=True, check=False)
    log = output / "verifier-controls.log"
    log.write_text(process.stdout + process.stderr, encoding="utf-8")
    return {"status": "PASS" if process.returncode == 0 else "ERROR", "command": command,
            "returncode": process.returncode, "log": str(log),
            "scope": "planted-defect fixture packages; they exercise U25 itself and are never validation data"}


def summary_reason(results: List[Dict[str, Any]]) -> str:
    """One line for the nightly table, which truncates long reasons.

    Counts and the passing datasets come first, then ERROR/FAIL details, then
    the (usually long and least urgent) SKIPPED details.
    """
    groups: Dict[str, List[Dict[str, Any]]] = {}
    for item in results:
        groups.setdefault(item["status"], []).append(item)
    counts = ", ".join(f"{len(groups[status])} {status}" for status in STATUS_ORDER if status in groups)
    parts = [f"{len(results)} datasets: {counts}"]
    if "PASS" in groups:
        parts.append("PASS: " + ", ".join(item["dataset"] for item in groups["PASS"]))
    for status in ("ERROR", "FAIL", "SKIPPED"):
        if status in groups:
            parts.append(f"{status}: " + ", ".join(f"{item['dataset']} ({item['reason']})"
                                                  for item in groups[status]))
    return ". ".join(parts)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def run_u25(output_root: Path, options: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Evaluate the selected (default: all declared) datasets."""
    options = options or {}
    output = output_root / "U25"
    output.mkdir(parents=True, exist_ok=True)
    try:
        contract = load_contract()
        registry = load_registry()
        check_registry(contract, registry)
        declared = declared_datasets()
    except Stop as stop:
        return {"status": stop.status, "reason": stop.reason}
    except (OSError, ValueError, KeyError) as error:
        return {"status": "ERROR", "reason": f"U25 contract/registry unreadable: {error}"}

    controls = run_controls(output)
    if controls["status"] != "PASS":
        return {"status": "ERROR", "reason": "U25 verifier controls did not PASS",
                "verifier_controls": controls}

    targets: List[str] = []
    for test_id in options.get("for_tests", []):
        if test_id not in declared:
            return {"status": "ERROR", "reason": f"{test_id} does not list U25 as a predecessor",
                    "declared_datasets": declared}
        targets.extend(dataset for dataset in declared[test_id] if dataset not in targets)
    targets.extend(dataset for dataset in options.get("datasets", []) if dataset not in targets)
    selection = "explicit" if targets else "declared by requesting tests"
    if not targets:
        targets = sorted({dataset for datasets in declared.values() for dataset in datasets})

    data_root = Path(options.get("data_root") or os.environ.get(DATA_ROOT_ENVIRONMENT, DEFAULT_DATA_ROOT))
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    work_root = output / "packages" / stamp
    cache: Dict[str, Dict[str, Any]] = {}
    results = []
    for dataset in targets:
        item = evaluate_dataset(dataset, registry, contract, data_root, work_root, cache)
        directory = output / dataset
        write_json(directory / "verified_hashes.json",
                   {"dataset": dataset, "package": item.get("package"),
                    "integrity_evaluated": item["integrity_evaluated"],
                    "files": item.get("verified_hashes", {}),
                    "selections": item.get("selections", [])})
        write_json(directory / "reproduction_report.json",
                   {"dataset": dataset, "reproduction": item["reproduction"]})
        brief = {key: value for key, value in item.items() if key != "verified_hashes"}
        write_json(directory / "result.json", brief)
        results.append(brief)

    status = worst([item["status"] for item in results])
    aggregate = {
        "status": status,
        "reason": summary_reason(results),
        "selection": selection,
        "targets": targets,
        "declared_datasets": declared,
        "data_root": str(data_root),
        "datasets": {item["dataset"]: {"status": item["status"], "reason": item["reason"]}
                     for item in results},
        "verifier_controls": controls,
        "validation_status": "NOT VALIDATED (data-package QA only)",
    }
    write_json(output / "verified_hashes.json",
               {item["dataset"]: str(output / item["dataset"] / "verified_hashes.json") for item in results})
    write_json(output / "reproduction_report.json",
               {item["dataset"]: item["reproduction"] for item in results})
    return aggregate


def main() -> int:
    """CLI used by U25_data_provenance/test.py."""
    sys.path.insert(0, str(TEST_ROOT / "common"))
    from moon_testlib import DEFAULT_BUILD, DEFAULT_OUTPUT, STATUS_EXIT, run_one

    parser = argparse.ArgumentParser(description="Run srcMoon U25 (validation-data package preflight)")
    parser.add_argument("--dataset", action="append", default=[], metavar="Dxx",
                        help="evaluate this dataset (repeatable)")
    parser.add_argument("--for-test", action="append", default=[], metavar="UXX",
                        help="evaluate the datasets declared by this requesting test (repeatable)")
    parser.add_argument("--data-root", type=Path, default=None,
                        help=f"validation-data root (default ${DATA_ROOT_ENVIRONMENT} or {DEFAULT_DATA_ROOT})")
    parser.add_argument("--build-dir", type=Path, default=DEFAULT_BUILD, help="unused; accepted for runner symmetry")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    for dataset in args.dataset:
        if not re.fullmatch(r"D\d{2}", dataset):
            parser.error(f"not a dataset identifier: {dataset}")
    options = {"datasets": args.dataset, "for_tests": args.for_test,
               "data_root": str(args.data_root) if args.data_root else None}
    result = run_one("U25", args.build_dir.resolve(), args.output_dir.resolve(), options=options)
    print(json.dumps(result, indent=2, sort_keys=True))
    return STATUS_EXIT[result["status"]]


if __name__ == "__main__":
    raise SystemExit(main())
