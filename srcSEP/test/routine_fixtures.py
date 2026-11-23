#!/usr/bin/env python3
"""Exact Routine selection and reviewed fixture routing for the maintained runner.

The maintained Python runner and Make gate use this descriptor-aware dispatch.
It imports the maintained runner's execution/report helpers and preserves its native callback
and validation-case protocols. Unimplemented native fixtures are ERROR, never
an excuse to omit a required ID or fabricate material ownership/thermodynamics.
The companion unit tests use explicitly labelled software protocol fixtures;
their PASS results make no claim about SEP physics or native model readiness.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "srcsep-routine-fixtures-v1"
CHECKPOINT_IDS = frozenset(("CP05C01", "CP05C02", "CP05C03", "CP05C04", "CP05C05"))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def framework_module():
    """Reuse the existing runner by code location, not an ambient import path."""
    spec = importlib.util.spec_from_file_location("sep_existing_runner", ROOT / "test/run_tests.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_catalog(text):
    """Read actual class/initialization metadata; never infer it from an ID.

    The native list is a documented six-column interface. The description may
    contain a pipe, so split only the five schema separators. Duplicate IDs or
    unknown classes/initialization spellings reject instead of being silently
    deduplicated or filtered out of the requested routine catalog.
    """
    rows, seen, current = [], set(), None
    for line in text.splitlines():
        detail = line.strip()
        # Native descriptors print these policy lines after their six-column
        # row. Preserve the ACTUAL seed/repeat/isolation/tags metadata as well
        # as the raw catalog. A changed capability is a manifest mismatch, not
        # permission to run a different seed schedule or skip a required role.
        if detail.startswith(("seed:", "state/isolation:", "seed-capability:", "tags:")):
            require(current is not None, "descriptor policy appears before its catalog row")
            if detail.startswith("seed-capability:"):
                match = re.fullmatch(r"seed-capability: (context|deterministic|fixed-or-unconverted); repeat-supported: (yes|no)", detail)
                require(match is not None and "seed_capability" not in current, "malformed/duplicate seed and repeat policy")
                current["seed_capability"], current["repeat_supported"] = match.group(1), match.group(2) == "yes"
            else:
                prefix, value = detail.split(":", 1)
                key = {"seed": "seed_policy", "state/isolation": "state_isolation", "tags": "tags"}[prefix]
                require(key not in current, "duplicate descriptor policy " + key)
                current[key] = value.strip().split() if key == "tags" else value.strip()
            continue
        if "|" not in line:
            continue
        fields = [field.strip() for field in line.split("|", 5)]
        if fields[0] == "ID":
            require(fields[:4] == ["ID", "group", "class", "initialization"], "foreign native catalog header")
            continue
        require(len(fields) == 6 and re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]*[0-9][A-Za-z0-9_-]*", fields[0]),
                "malformed native descriptor row")
        require(fields[0] not in seen, "duplicate native descriptor: " + fields[0])
        require(fields[2] in ("routine", "extended"), "unknown native runtime class")
        require(fields[3] in ("none", "field-line-model", "native-storage"), "unknown native initialization class")
        seen.add(fields[0])
        current = dict(zip(("id", "group", "class", "initialization", "build_modes", "description"), fields))
        rows.append(current)
    require(rows, "empty native catalog")
    return sorted((row for row in rows if row["class"] == "routine"), key=lambda row: row["id"])


def inventory_digest(rows):
    return hashlib.sha256(json.dumps(rows, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def validate_manifest(manifest, rows):
    """Bind exact routine membership and reviewed routes before execution.

    A manifest is a fixture/configuration declaration, not expected-PASS data.
    No entry may authorize a different native callback or change acceptance
    metrics. Unknown/missing routes fail closed with the full required catalog.
    """
    require(manifest.get("schema") == SCHEMA, "unknown routine fixture schema")
    ids = [row["id"] for row in rows]
    require(manifest.get("expected_ids") == ids, "routine descriptor membership/order differs from reviewed manifest")
    require(manifest.get("routine_inventory_sha256") == inventory_digest(rows), "routine descriptor metadata changed")
    require(set(manifest.get("fixtures", {})) == set(ids), "missing/foreign fixture mapping")
    for row in rows:
        route = manifest["fixtures"][row["id"]]
        kind = route.get("kind")
        require(kind in ("native", "native-storage", "validation", "material-checkpoint-v2", "birth-operator", "unavailable"), "unknown fixture kind")
        # Only the exact NativeStorage contract below owns startup arguments.
        # Other routes cannot smuggle a second input/selector/report after the
        # canonical hashed request, even if its child happened to emit one PASS.
        require(kind == "native-storage" or "startup_args" not in route,
                "startup_args are forbidden outside an exact native-storage contract")
        if kind == "native":
            require(row["initialization"] == "none", "native storage/particle fixture cannot be assumed from allocation metadata")
        if kind == "native-storage":
            require(row["initialization"] == "native-storage", "storage route requires explicit allocation-only metadata")
            approved = (["--shock-model", "analytical"] if row["id"] in ("FTE01", "PARKER01") else
                        ["--shock-model", "swcme1d", "--turbulence-model", "wave-number-resolved"] if row["id"] in ("PNW01", "PNW02", "PNW03", "PNW04") else [])
            require(route.get("startup_args", []) == approved, "storage startup differs from declared fixture contract")
        if kind == "material-checkpoint-v2":
            require(row["id"] in CHECKPOINT_IDS, "checkpoint protocol cannot be assigned to another test ID")
        if kind == "birth-operator":
            require(row["id"] in {"BIRTHN01", "BIRTHN02", "BIRTHN03", "BIRTHN04", "BIRTHN05", "BIRTHN06"}, "birth protocol cannot be assigned to another test ID")
        if kind == "unavailable":
            require(isinstance(route.get("reason"), str) and route["reason"], "unavailable fixture needs an explicit reason")


def quote_protocol(value):
    """Match C++ std::quoted without guessing nested file/path semantics."""
    require(not any(character in value for character in "\0\r\n"), "invalid multiline protocol path")
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'


def checkpoint_request(route, directory, source_root):
    """Stage only the existing approved CP05 protocol and a fresh output path.

    The absolute startup/event paths keep their existing parser-owned relative
    asset bases. Their hashes are required; no input value, physical dt, source
    spectrum, epoch or normalization is synthesized here. The native protocol
    implementation remains the sole owner of its already-qualified physics.
    """
    inputs = {}
    for name in ("startup_input", "event"):
        record = route[name]
        path = Path(record["path"])
        path = (path if path.is_absolute() else source_root / path).resolve()
        require(path.is_file() and digest(path) == record["sha256"], "changed/missing checkpoint fixture " + name)
        inputs[name] = path
    request = directory / "checkpoint-request.args"
    prefix = directory / "checkpoint"
    text = "SEP_CHECKPOINT_FRESH_PROCESS 2\nrole \"routine\"\nprefix " + quote_protocol(str(prefix))
    text += "\nstartup_input " + quote_protocol(str(inputs["startup_input"]))
    text += "\nevent " + quote_protocol(str(inputs["event"])) + "\n"
    request.write_text(text)
    return request, inputs


def birth_campaign_arguments(arguments):
    """Translate only existing campaign controls into complete fresh role groups.

    Model/physics flags cannot be silently ignored by a delegated protocol.
    Native CLI last-occurrence-wins semantics are preserved for the two
    supported controls. The role runner validates ranges and invokes the same
    native DeriveTestCaseSeed for each explicit replicate index.
    """
    mapping = {"--test-seed": "--seed", "--test-repeat": "--repeat-count"}
    values = {}
    iterator = iter(arguments)
    for argument in iterator:
        option, separator, value = argument.partition("=")
        require(option in mapping, "birth protocol does not support forwarded native option " + option)
        if not separator:
            try:
                value = next(iterator)
            except StopIteration:
                raise ValueError("missing native campaign-control value")
        values[mapping[option]] = value
    return [value for option in mapping.values() if option in values for value in (option, values[option])]


def load_result(framework, path, identity, code):
    """Preserve native status/metrics/seed; reject foreign or contradictory data."""
    report = framework._load_report(path)
    rows = report.get("results", [])
    require(len(rows) == 1 and rows[0].get("id") == identity, "child returned missing/duplicate/foreign case records")
    result = dict(rows[0])
    status = result.get("status")
    require(status in ("PASS", "FAIL", "SKIP", "ERROR"), "invalid child status")
    require(not (status in ("PASS", "SKIP") and code != 0), "child success status contradicts its process exit")
    return result


def run(args, framework=None):
    """Execute all reviewed routine IDs, retaining missing fixture errors.

    This dispatcher has no permissive mode or physics-guard override. Model
    arguments are forwarded unchanged; selectors/report/input arguments are
    reserved so they cannot silently replace the case or its reviewed fixture.
    All callbacks retain their actual compiled seed/repeat behavior. Required
    SKIP is preserved as evidence but cannot close this dispatcher's full gate.
    """
    framework = framework or framework_module()
    output = args.output.resolve()
    require(not output.exists(), "choose a fresh routine evidence directory")
    output.mkdir(parents=True)
    executable = args.amps.resolve()
    before = digest(executable)
    manifest_bytes = args.manifest.read_bytes()
    manifest_hash = hashlib.sha256(manifest_bytes).hexdigest()
    manifest = json.loads(manifest_bytes)
    (output / "parsed-fixture-manifest.json").write_bytes(manifest_bytes)
    expected = manifest.get("expected_ids", [])
    require(isinstance(expected, list) and expected and
            all(isinstance(identity, str) and re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]*[0-9][A-Za-z0-9_-]*", identity) for identity in expected) and
            len(set(expected)) == len(expected), "invalid required ID catalog")
    reserved = {"--test", "--tests", "--input", "-input", "--test-group", "--all-tests", "--list-tests", "--test-input", "--test-output-dir", "--test-output", "--test-json", "--test-junit"}
    require(not any(value.partition("=")[0] in reserved for value in args.native_args), "native arguments override reserved case/fixture selection")
    log = output / "driver.log"
    log.touch()
    discovery = output / "discovery"
    discovery.mkdir()
    reports, commands, rows, preflight = [], [], [], None
    try:
        command = [str(executable), "--list-tests"]
        status = framework._run_streaming(command, discovery, log, dict(os.environ), args.timeout)
        require(status == 0, "native discovery process failed")
        catalog = (output / "catalog.txt")
        # Discovery has no physics input; retain the exact recorded list while
        # using the same helper's complete log rather than a truncated tail.
        text = log.read_text()
        catalog.write_text(text)
        rows = parse_catalog(text)
        validate_manifest(manifest, rows)
        requested=getattr(args,"requested_ids",None)
        if requested is not None:
            # An explicit caller selection narrows execution, never the full
            # catalog admission audit. --routine leaves this unset and still
            # requires the complete catalog. --all reuses every Routine result before adding
            # every Extended descriptor, and single-ID adapters reuse exactly
            # their reviewed startup without fabricating another callback.
            require(requested and len(set(requested))==len(requested) and set(requested)<=set(expected),
                    "requested fixture IDs are absent/duplicate/non-Routine")
            expected=sorted(requested)
    except (ValueError, OSError, framework.RunnerError) as error:
        preflight = str(error)
    if preflight and rows:
        # A stale/missing manifest entry is a gate ERROR, not authorization to
        # drop a compiled Routine test. Preserve the union for error evidence;
        # no callbacks execute when descriptor membership or metadata differs.
        expected = sorted(set(expected) | {row["id"] for row in rows})
    by_id = {row["id"]: row for row in rows}
    for identity in expected:
        started = time.monotonic()
        directory = output / identity
        directory.mkdir()
        native_json, native_xml = directory / "native.json", directory / "native.xml"
        normalized, junit = directory / "normalized.json", directory / "normalized.xml"
        case_log = directory / "case.log"
        case_log.touch()
        route = manifest.get("fixtures", {}).get(identity, {})
        command, status, error, inputs = [], None, preflight, {}
        try:
            if error:
                raise ValueError(error)
            kind = route["kind"]
            if kind == "unavailable":
                raise ValueError("required native fixture is not implemented: " + route["reason"])
            if kind == "birth-operator":
                # The maintained role driver owns the fixed physical startup,
                # complete independent operator controls, and final native
                # comparisons. Repeating a role in one process would reuse its
                # geometry owner, so repeat the complete experiment instead.
                command = [sys.executable, str(framework.ROOT / "test/run_birth_operator_stage_tests.py"),
                           "--amps", str(executable), "--test-id", identity,
                           "--mpi-ranks", str(getattr(args, "mpi_np", None) or 2), "--mpi-exec", getattr(args, "mpiexec", "mpiexec"), "--output", str(directory / "birth")]
                command += birth_campaign_arguments(args.native_args)
                native_json = directory / "birth/srcsep-tests.json"
            elif kind == "validation":
                require(not args.native_args, "maintained validation protocol owns seed/repeat/model arguments; forwarded native overrides are unsupported")
                require(identity in framework._validation_case_ids(), "foreign validation-case route")
                command = [sys.executable, str(framework.ROOT / "validation/run_case.py"),
                           "--amps", str(executable), "--case", identity, "--output-dir", str(directory / "validation")]
                if args.timeout is not None:
                    command += ["--timeout", str(args.timeout)]
                native_json = directory / "validation/srcsep-tests.json"
            else:
                require(identity not in framework._validation_case_ids(), "validation case requires its maintained input/reference workflow")
                command = [str(executable), "--test", identity, "--test-json", str(native_json), "--test-junit", str(native_xml)]
                if kind == "material-checkpoint-v2":
                    request, inputs = checkpoint_request(route, directory, framework.ROOT)
                    command += ["--test-input", str(request), "--test-output-dir", str(directory / "artifacts")]
                else:
                    command += ["--test-output", str(directory / "artifacts")]
                command += route.get("startup_args", [])
                command += args.native_args
                if getattr(args, "mpi_np", None):
                    command = [getattr(args, "mpiexec", "mpiexec"), "-n", str(args.mpi_np), *command]
            commands.append(command)
            # The maintained validation runner owns its several native and
            # independent-reference commands; do not give the enclosing case
            # an artificially shorter aggregate deadline.
            timeout = None if kind in ("validation", "birth-operator") else args.timeout
            status = framework._run_streaming(command, directory, case_log, dict(os.environ), timeout)
            result = load_result(framework, native_json, identity, status)
            for name, path in inputs.items():
                require(digest(path) == route[name]["sha256"], "fixture input changed during execution")
            require(digest(executable) == before, "executable changed during execution")
            require(digest(args.manifest) == manifest_hash, "parsed fixture manifest changed during execution")
        except (ValueError, OSError, KeyError, framework.RunnerError) as failure:
            error = str(failure)
            result = framework._synthetic_error(identity, error, time.monotonic() - started, case_log)
            if digest(args.manifest) != manifest_hash:
                preflight = "parsed fixture manifest changed during execution; dependent callbacks blocked"
        result.setdefault("configuration", []).append("routine_fixture_route=" + route.get("kind", "unknown"))
        framework._write_json(normalized, {"schema": "srcsep-component-tests-v1", "exit_code": 2 if result["status"] in ("FAIL", "ERROR") else 0,
                                          "totals": framework._report_totals([result]), "results": [result]})
        framework._write_junit(junit, {"totals": framework._report_totals([result]), "results": [result]})
        framework._write_json(directory / "process.json", {"argv": command, "cwd": str(directory), "actual_exit_code": status,
                                                            "launch_error": error, "native_report": str(native_json), "metadata": by_id.get(identity)})
        reports.append(normalized)
    # Even a campaign consisting entirely of declared unavailable fixtures
    # must bind the exact bytes parsed at startup. Never relabel old execution
    # with a later on-disk manifest hash. A changed declaration cannot qualify.
    final_manifest_hash = digest(args.manifest)
    if final_manifest_hash != manifest_hash:
        # Keep all required IDs as ERROR and retain every actual native/raw
        # report. Previously completed PASS rows cannot qualify a campaign
        # whose controlling declaration changed partway through execution.
        for path in reports:
            row = framework._load_report(path)["results"][0]
            invalid = framework._synthetic_error(row["id"],
                "parsed fixture manifest changed during campaign; original declaration preserved in parsed-fixture-manifest.json", 0.0, log)
            framework._write_json(path, {"results": [invalid], "totals": framework._report_totals([invalid]), "exit_code": 2})
    aggregate = framework._merge_reports(reports, output / "srcsep-tests.json")
    require([row["id"] for row in aggregate["results"]] == expected, "aggregate lost/reordered required routine IDs")
    incomplete = any(row["status"] != "PASS" for row in aggregate["results"])
    # Keep native SKIP verbatim while making missing required evidence visible
    # through the authoritative aggregate process status; never turn it into a
    # false PASS or discard its record from the mandatory catalog.
    aggregate["exit_code"] = max(int(aggregate["exit_code"]), 2 if incomplete and aggregate["totals"].get("skipped") else 0)
    aggregate["scope"] = "exact routine routing; unavailable fixtures remain required ERROR"
    framework._write_json(output / "srcsep-tests.json", aggregate)
    framework._write_junit(output / "srcsep-tests.xml", aggregate)
    framework._write_json(output / "run_manifest.json", {"schema": "srcsep-routine-isolation-v1", "expected_ids": expected,
        "catalog": rows, "fixture_manifest": str(args.manifest), "fixture_manifest_sha256": manifest_hash, "fixture_manifest_final_sha256": final_manifest_hash,
        "executable": str(executable), "executable_sha256": before, "commands": commands, "actual_exit_code": aggregate["exit_code"],
        "qualification": "required native fixtures and full validation workflows; no missing fixture waiver"})
    return int(aggregate["exit_code"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amps", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--timeout", type=float)
    parser.add_argument("--mpi-np", type=int)
    parser.add_argument("--mpiexec", default="mpiexec")
    parser.add_argument("native_args", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    require(args.mpi_np is None or args.mpi_np > 0, "MPI rank count must be positive")
    if args.native_args and args.native_args[0] == "--":
        args.native_args = args.native_args[1:]
    return run(args)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, OSError) as error:
        print("ERROR: " + str(error), file=sys.stderr)
        sys.exit(2)
