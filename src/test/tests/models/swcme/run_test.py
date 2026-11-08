#!/usr/bin/env python3
"""Import the complete *available* SWCME test surface without overstating it.

SWCME currently has three materially different local registries:

* the four-case R1 runner is executable and produces JSON/JUnit evidence;
* a 127-case C++ registry is preserved in source, but six translation units
  named by its Makefile are absent; and
* the VP01--VP16 observational registry is preserved in ``campaign.json``,
  but its authoritative campaign runner is absent.

The global adapter therefore executes only R1.  It emits explicit SKIP records
for the dormant extended registry and SKIP/BLOCKED records for the
observational registry.  In particular, a four-case R1 PASS can never be
mistaken for full SWCME or observational qualification.

Every registry is validated in full before ``selected_cases`` is honored.
That ordering is intentional: selecting one case must not hide a changed or
truncated local registry from the global control plane.
"""

from __future__ import annotations

import ast
import hashlib
import json
import re
import shutil
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


TEST_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(TEST_ROOT))

from lib.runner_support import (  # noqa: E402
    RunnerError, case, command, error_result, load_context, main_guard,
    run_command, write_result)


# These definitions mirror the executable tuple in the authoritative local
# Python runner.  Groups and titles are pinned as well as IDs so a semantic
# registry edit cannot pass unnoticed while the global manifest remains stale.
EXPECTED_R1: Tuple[Tuple[str, str, str], ...] = (
    ("R1CFG01", "COMMON", "Canonical constant identity"),
    ("R1D01", "1D", "One-dimensional prepared-state smoke"),
    ("R3D01", "3D", "Three-dimensional prepared-state smoke"),
    ("R13D01", "1D3D", "Common ambient-state identity"),
)


# Declaration order is execution order in core/test_swcme.cpp.  The registry
# is deliberately retained here even though it is not executable: this gives
# every dormant local case a stable global identity and makes additions,
# removals, or reordering an explicit manifest/adapter update rather than a
# silently narrowed qualification claim.
EXPECTED_EXTENDED_IDS: Tuple[str, ...] = tuple("""
PST01 PST02 PST03 PST06 PST04 OUT02 OUT03 OUT05 OUT04 OUT06 OUT01 OUT07
PST07 PST05 PST08 OUT08 CFG03 CFG04 KIN09 CON09 CON10
SHK05 SHK12 SHK16 SHK15 SHK06 SHK07 SHK08 SHK09 SHK10 SHK11 SHK17
DEN03 DEN04 DEN02 DEN05 PAR04 PAR05 PAR06 PAR02 PAR03 PAR01 CROSS01
CFG01 CFG02 DEF01 DEF02 DEF03 DEF04
SEP01 SEP02 SEP03 SEP04 SEP05 SEP06 DEN01 GEO01 MESH01 REG01 CON11
SAN01 THR01 COV01 PERF01 REP01 EVT01 V1 V2 V3 V4 V5 V6
GEO02 GEO03 GEO04 GEO05 GEO06 GEO07 GEO08
MSH01 MSH02 MSH03 MSH04 MSH05
KIN01 KIN02 KIN03 KIN04 KIN05 KIN06 KIN07 KIN08
SHK01 SHK02 SHK03 SHK04 SHK13 SHK14
CON01 CON02 CON03 CON04 CON05 CON06 CON07 CON08
REG02 REG03 REG04 REG05 ACC01 ACC02 ACC03 ACC04 ACC05
ERR01 ERR02 ERR03 ERR04 ERR05 DIV01 DIV02 DIV03
1D_AMBIENT_01 1D3D01 1D3D02 1D3D03
""".split())


EXPECTED_MISSING_EXTENDED_SOURCES: Tuple[str, ...] = (
    "3d/test_parker.cpp",
    "3d/test_geometry.cpp",
    "3d/test_mesh.cpp",
    "3d/test_connectivity.cpp",
    "3d/test_shock_surface.cpp",
    "test_1d3d.cpp",
)


EXPECTED_VP_DEPENDENCIES: Mapping[str, Tuple[str, ...]] = {
    "VP01": (),
    "VP02": (),
    "VP03": (),
    "VP04": ("VP01", "VP02", "VP03"),
    "VP05": ("VP01", "VP02", "VP03"),
    "VP06": ("VP01", "VP02", "VP03", "VP04", "VP05"),
    "VP07": ("VP01", "VP02", "VP03", "VP04", "VP05"),
    "VP08": ("VP01", "VP02", "VP03", "VP04", "VP05", "VP06", "VP07"),
    "VP09": ("VP01", "VP02", "VP03", "VP04", "VP05", "VP06", "VP07"),
    "VP10": ("VP01", "VP02", "VP03", "VP04", "VP05", "VP06", "VP07", "VP08", "VP09"),
    "VP11": ("VP01", "VP02", "VP03", "VP04", "VP05", "VP06", "VP07", "VP08", "VP09"),
    "VP12": ("VP01", "VP02", "VP03", "VP04", "VP05", "VP06", "VP07", "VP08", "VP09", "VP10", "VP11"),
    "VP13": ("VP01", "VP02", "VP03", "VP04", "VP05", "VP06", "VP07", "VP08", "VP09", "VP10", "VP11"),
    "VP14": ("VP01", "VP02", "VP03", "VP04", "VP05", "VP06", "VP07", "VP08", "VP09", "VP10", "VP11"),
    "VP15": ("VP01", "VP02", "VP03", "VP04", "VP05", "VP06", "VP07", "VP08", "VP09", "VP10", "VP11"),
    "VP16": ("VP01", "VP02", "VP03", "VP04", "VP05", "VP06", "VP07", "VP08", "VP09", "VP10", "VP11", "VP12", "VP13", "VP14", "VP15"),
}
EXPECTED_VP_IDS: Tuple[str, ...] = tuple(EXPECTED_VP_DEPENDENCIES)


def r1_global(local_id: str) -> str:
    """Keep the short established R1 spelling while adding model ownership."""

    return "SWCME-" + local_id


def extended_global(local_id: str) -> str:
    """Convert the one legacy underscore ID to the global ID grammar."""

    return "SWCME-EXT-" + local_id.replace("_", "-")


def vp_global(local_id: str) -> str:
    return "SWCME-" + local_id


R1_GLOBAL_IDS: Tuple[str, ...] = tuple(r1_global(item[0]) for item in EXPECTED_R1)
EXTENDED_GLOBAL_IDS: Tuple[str, ...] = tuple(
    extended_global(item) for item in EXPECTED_EXTENDED_IDS)
VP_GLOBAL_IDS: Tuple[str, ...] = tuple(vp_global(item) for item in EXPECTED_VP_IDS)
ALL_GLOBAL_CASES: Tuple[str, ...] = (
    R1_GLOBAL_IDS + EXTENDED_GLOBAL_IDS + VP_GLOBAL_IDS)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_r1_registry(path: Path) -> Tuple[Tuple[str, str, str], ...]:
    """Read ``TESTS`` through Python's AST without importing executable code."""

    try:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    except (OSError, SyntaxError) as error:
        raise RunnerError("cannot parse SWCME R1 registry {}: {}".format(path, error))
    for node in tree.body:
        # The local declaration is annotated today, but accepting an ordinary
        # assignment keeps this parser focused on the literal registry rather
        # than on one harmless Python spelling of that declaration.
        if isinstance(node, ast.AnnAssign):
            is_tests = isinstance(node.target, ast.Name) and node.target.id == "TESTS"
            value = node.value
        elif isinstance(node, ast.Assign):
            is_tests = any(isinstance(target, ast.Name) and target.id == "TESTS"
                           for target in node.targets)
            value = node.value
        else:
            continue
        if not is_tests:
            continue
        if not isinstance(value, (ast.Tuple, ast.List)):
            break
        records: List[Tuple[str, str, str]] = []
        for entry in value.elts:
            if (not isinstance(entry, ast.Call) or
                    not isinstance(entry.func, ast.Name) or
                    entry.func.id != "Definition" or len(entry.args) != 3):
                raise RunnerError("SWCME R1 TESTS contains a non-Definition entry")
            values = []
            for argument in entry.args:
                if not isinstance(argument, ast.Constant) or not isinstance(argument.value, str):
                    raise RunnerError("SWCME R1 Definition fields must be string literals")
                values.append(argument.value)
            records.append((values[0], values[1], values[2]))
        return tuple(records)
    raise RunnerError("SWCME R1 runner does not define a literal TESTS registry")


def parse_extended_registry(path: Path) -> Tuple[Tuple[str, str, str], ...]:
    """Extract the exact C++ ``kTests`` array, including wrapped entries."""

    try:
        text = path.read_text(encoding="utf-8")
        start = text.index("const swcme_test::TestCase kTests[] = {")
        stop = text.index("constexpr std::size_t kTestCount", start)
    except (OSError, ValueError) as error:
        raise RunnerError("cannot locate extended SWCME registry: {}".format(error))
    pattern = re.compile(
        r'\{\s*"([^"]+)"\s*,\s*"([^"]+)"\s*,\s*"([^"]+)"\s*,'
        r'\s*[A-Za-z_][A-Za-z0-9_]*\s*\}', re.MULTILINE)
    return tuple(pattern.findall(text[start:stop]))


def missing_extended_sources(makefile: Path, test_root: Path) -> Tuple[str, ...]:
    """Resolve the Makefile's own source inventory instead of guessing files."""

    try:
        text = makefile.read_text(encoding="utf-8")
        start = text.index("SOURCES :=") + len("SOURCES :=")
        stop = text.index("OBJECTS :=", start)
    except (OSError, ValueError) as error:
        raise RunnerError("cannot parse SWCME Makefile SOURCES: {}".format(error))
    tokens = text[start:stop].replace("\\\n", " ").split()
    sources = [token for token in tokens if token.endswith(".cpp")]
    return tuple(source for source in sources if not (test_root / source).is_file())


def load_vp_registry(path: Path, validation_root: Path) -> Tuple[Dict, ...]:
    """Validate campaign ordering, state, dependencies, and case contracts."""

    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise RunnerError("cannot load SWCME VP registry {}: {}".format(path, error))
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        raise RunnerError("SWCME VP registry must be a schema_version 1 object")
    raw_cases = payload.get("cases")
    if not isinstance(raw_cases, list):
        raise RunnerError("SWCME VP registry cases must be a list")
    records: List[Dict] = []
    for index, raw in enumerate(raw_cases, start=1):
        if not isinstance(raw, dict):
            raise RunnerError("SWCME VP registry entry {} is not an object".format(index))
        local_id = raw.get("id")
        expected_id = "VP{:02d}".format(index)
        if local_id != expected_id or raw.get("priority") != index:
            raise RunnerError("SWCME VP priority/ID drift at entry {}".format(index))
        if raw.get("state") != "implemented":
            raise RunnerError("SWCME {} state drifted from implemented".format(local_id))
        if raw.get("case_dir") != local_id.lower():
            raise RunnerError("SWCME {} case_dir is not canonical".format(local_id))
        dependencies = tuple(raw.get("dependencies", ()))
        if dependencies != EXPECTED_VP_DEPENDENCIES[local_id]:
            raise RunnerError("SWCME {} dependency registry drift".format(local_id))

        case_dir = validation_root / raw["case_dir"]
        contract_path = case_dir / "case.json"
        try:
            contract = json.loads(contract_path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as error:
            raise RunnerError("cannot load {} contract: {}".format(local_id, error))
        expected_runner = [
            "{python}", "{case_dir}/run_" + local_id.lower() + ".py",
            "--output-dir", "{output_dir}",
        ]
        if (contract.get("schema_version") != 1 or
                contract.get("id") != local_id or
                contract.get("runner") != expected_runner or
                contract.get("result_file") != "{}_result.json".format(local_id.lower())):
            raise RunnerError("SWCME {} case contract drift".format(local_id))
        record = dict(raw)
        # VP01's case contract currently names an absent case runner.  Retain
        # that additional local limitation in evidence/messages; it is not a
        # reason to delete VP01 from the authoritative campaign registry.
        record["case_runner_present"] = (
            case_dir / "run_{}.py".format(local_id.lower())).is_file()
        records.append(record)
    if tuple(record["id"] for record in records) != EXPECTED_VP_IDS:
        raise RunnerError("SWCME VP registry membership/order drift")
    return tuple(records)


def audit_registries(repo: Path, artifacts: Path) -> Tuple[
        Tuple[Tuple[str, str, str], ...], Tuple[Dict, ...], Path]:
    """Validate all local authorities and retain their hashes as evidence."""

    model_root = repo / "src" / "models" / "swcme"
    test_root = model_root / "test"
    r1_path = test_root / "run_tests.py"
    extended_path = test_root / "core" / "test_swcme.cpp"
    makefile = test_root / "Makefile"
    validation_root = test_root / "validation"
    campaign_path = validation_root / "campaign.json"

    r1 = parse_r1_registry(r1_path)
    if r1 != EXPECTED_R1:
        raise RunnerError("SWCME R1 registry drift: expected {} observed {}".format(
            EXPECTED_R1, r1))
    extended = parse_extended_registry(extended_path)
    observed_extended_ids = tuple(record[0] for record in extended)
    if observed_extended_ids != EXPECTED_EXTENDED_IDS:
        raise RunnerError(
            "SWCME extended registry drift: expected {} IDs, observed {}".format(
                len(EXPECTED_EXTENDED_IDS), len(observed_extended_ids)))
    missing = missing_extended_sources(makefile, test_root)
    if missing != EXPECTED_MISSING_EXTENDED_SOURCES:
        raise RunnerError("SWCME extended-source blocker drift: expected {} observed {}".format(
            EXPECTED_MISSING_EXTENDED_SOURCES, missing))
    vp = load_vp_registry(campaign_path, validation_root)

    # A restored campaign runner changes the correct behavior from typed SKIP
    # to execution.  Refuse to keep reporting the old limitation in that case.
    campaign_runner = validation_root / "run_validation.py"
    if campaign_runner.exists():
        raise RunnerError(
            "SWCME observational runner now exists; replace the typed-SKIP adapter")

    artifacts.mkdir(parents=True, exist_ok=True)
    audit_path = artifacts / "swcme-registry-audit.json"
    audit_path.write_text(json.dumps({
        "schema_version": 1,
        "r1": {"count": len(r1), "source": str(r1_path),
               "sha256": sha256(r1_path), "ids": [item[0] for item in r1]},
        "extended": {"count": len(extended), "source": str(extended_path),
                     "sha256": sha256(extended_path),
                     "ids": list(observed_extended_ids),
                     "missing_sources": list(missing)},
        "observational": {"count": len(vp), "source": str(campaign_path),
                          "sha256": sha256(campaign_path),
                          "ids": [item["id"] for item in vp],
                          "missing_case_runners": [item["id"] for item in vp
                                                   if not item["case_runner_present"]],
                          "authoritative_runner": "absent"},
        "global_case_count": len(ALL_GLOBAL_CASES),
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return extended, vp, audit_path


def copy_model_isolated(source: Path, destination: Path) -> None:
    """Copy sources without carrying mutable build products into the run."""

    if destination.exists():
        raise RunnerError("isolated SWCME work directory already exists: {}".format(
            destination))
    ignore = shutil.ignore_patterns(
        "output", "build", "__pycache__", "*.a", "*.o", "*.d", "*.pyc")
    shutil.copytree(str(source), str(destination), symlinks=True, ignore=ignore)


def execute_r1(context: Mapping, repo: Path, work: Path, artifacts: Path,
               timeout: int) -> Tuple[Dict[str, Dict], List[str]]:
    """Run all four local cases once, then return globally namespaced records."""

    isolated_model = work / "swcme-isolated"
    copy_model_isolated(repo / "src" / "models" / "swcme", isolated_model)
    local_output = artifacts / "r1-local"
    jobs = int(context["execution"].get("build_jobs", 1))
    if jobs < 1:
        raise RunnerError("execution.build_jobs must be positive")
    argv = command(context, "python") + [
        str(isolated_model / "test" / "run_tests.py"),
        "--all", "--jobs", str(jobs), "--verbose",
        "--output-dir", str(local_output),
    ]
    completed = run_command(argv, isolated_model / "test", timeout=timeout)
    transcript = artifacts / "swcme-r1-transcript.txt"
    transcript.write_text(completed.stdout or "", encoding="utf-8")
    report = local_output / "swcme-r1-tests.json"
    evidence = [str(transcript)]
    if not report.is_file():
        message = "authoritative R1 runner exited {} without its JSON report".format(
            completed.returncode)
        return ({r1_global(local): case(r1_global(local), "ERROR", message)
                 for local, _, _ in EXPECTED_R1}, evidence)

    evidence.extend([str(report), str(local_output / "swcme-r1-tests.xml")])
    try:
        payload = json.loads(report.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise RunnerError("invalid SWCME R1 report: {}".format(error))
    if payload.get("schema") != "swcme-r1-test-summary-v1":
        raise RunnerError("unexpected SWCME R1 report schema")
    results = payload.get("results")
    if not isinstance(results, list):
        raise RunnerError("SWCME R1 report results must be a list")
    observed = [item.get("test_id") for item in results if isinstance(item, dict)]
    expected_local = [item[0] for item in EXPECTED_R1]
    if observed != expected_local or len(set(observed)) != len(observed):
        raise RunnerError("SWCME R1 report membership/order drift: {}".format(observed))
    allowed = {"PASS", "FAIL", "SKIP", "ERROR"}
    totals = {status: sum(item.get("status") == status for item in results)
              for status in allowed}
    if payload.get("totals") != {status: totals[status]
                                 for status in ("PASS", "FAIL", "SKIP", "ERROR")}:
        raise RunnerError("SWCME R1 report totals disagree with records")
    expected_exit = 2 if totals["ERROR"] else (1 if totals["FAIL"] else 0)
    if completed.returncode != expected_exit:
        raise RunnerError("SWCME R1 exit/report mismatch: exit {} expected {}".format(
            completed.returncode, expected_exit))

    by_global: Dict[str, Dict] = {}
    definitions = {item[0]: item for item in EXPECTED_R1}
    for item in results:
        local_id = item["test_id"]
        status = item.get("status")
        if status not in allowed:
            raise RunnerError("SWCME {} has invalid status {}".format(local_id, status))
        _, group, title = definitions[local_id]
        message = "local {} [{}] {}: {}".format(
            local_id, group, title, item.get("message", ""))
        global_id = r1_global(local_id)
        by_global[global_id] = case(global_id, status, message)
    return by_global, evidence


def skipped_extended(records: Sequence[Tuple[str, str, str]]) -> Dict[str, Dict]:
    missing = ", ".join(EXPECTED_MISSING_EXTENDED_SOURCES)
    return {
        extended_global(local_id): case(
            extended_global(local_id), "SKIP",
            "local {} [{}] {}: extended registry is not executable; "
            "missing Makefile sources: {}".format(local_id, group, title, missing))
        for local_id, group, title in records
    }


def skipped_vp(records: Sequence[Mapping]) -> Dict[str, Dict]:
    """Preserve implemented state and distinguish roots from blocked cases."""

    results: Dict[str, Dict] = {}
    for record in records:
        local_id = str(record["id"])
        dependencies = tuple(record["dependencies"])
        if dependencies:
            typed = "BLOCKED"
            detail = "dependencies {} cannot be orchestrated".format(
                ",".join(dependencies))
        else:
            typed = "SKIP"
            detail = "root validation case cannot be orchestrated"
        global_id = vp_global(local_id)
        case_runner_note = ""
        if not record.get("case_runner_present", False):
            case_runner_note = "; declared case runner is also absent"
        results[global_id] = case(
            global_id, "SKIP",
            "local {} status={}; state={}; layer={}; title={}; {} because authoritative "
            "validation/run_validation.py is absent{}; limitation: no observational "
            "qualification was executed".format(
                local_id, typed, record["state"], record["layer"],
                record["title"], detail, case_runner_note))
    return results


def aggregate_status(records: Iterable[Mapping]) -> str:
    statuses = {str(record["status"]) for record in records}
    for status in ("ERROR", "FAIL", "SKIP", "PASS"):
        if status in statuses:
            return status
    raise RunnerError("cannot aggregate an empty SWCME result set")


def main() -> int:
    context = load_context("Import authoritative SWCME model-local registries")
    repo = Path(context["repo_root"])
    work = Path(context["work_dir"])
    artifacts = Path(context["artifact_dir"])
    timeout = int(context["test"].get("timeout_seconds", 1800))

    try:
        extended, vp, audit_path = audit_registries(repo, artifacts)

        # The central manifest is part of the qualification boundary.  Refuse
        # to run if it omits any locally registered case or adds an identity
        # that this adapter cannot source from a local authority.
        configured = tuple(context["test"].get("expected_cases", ()))
        if configured != ALL_GLOBAL_CASES:
            raise RunnerError(
                "SWCME manifest case drift: expected {} ordered cases, observed {}".format(
                    len(ALL_GLOBAL_CASES), len(configured)))
        requested = tuple(context["test"].get("selected_cases", configured))
        if (not requested or len(set(requested)) != len(requested) or
                any(identifier not in ALL_GLOBAL_CASES for identifier in requested)):
            raise RunnerError("invalid SWCME selected_cases: {}".format(requested))

        all_results: Dict[str, Dict] = {}
        evidence = [str(audit_path)]
        if any(identifier in R1_GLOBAL_IDS for identifier in requested):
            r1_results, r1_evidence = execute_r1(
                context, repo, work, artifacts, timeout)
            all_results.update(r1_results)
            evidence.extend(r1_evidence)
        all_results.update(skipped_extended(extended))
        all_results.update(skipped_vp(vp))

        selected_results = [all_results[identifier] for identifier in requested]
        status = aggregate_status(selected_results)
        counts = {value: sum(item["status"] == value for item in selected_results)
                  for value in ("PASS", "FAIL", "SKIP", "ERROR")}
        message = ("{} selected SWCME cases: PASS={} FAIL={} SKIP={} ERROR={}; "
                   "registry inventory R1=4 extended=127 observational=16").format(
                       len(selected_results), counts["PASS"], counts["FAIL"],
                       counts["SKIP"], counts["ERROR"])
        return write_result(
            context, status, message, cases=selected_results,
            evidence=evidence,
            metrics={"registered_cases": len(ALL_GLOBAL_CASES),
                     "r1_cases": len(R1_GLOBAL_IDS),
                     "extended_cases": len(EXTENDED_GLOBAL_IDS),
                     "observational_cases": len(VP_GLOBAL_IDS)})
    except Exception as error:
        return error_result(context, error)


if __name__ == "__main__":
    main_guard(main)
