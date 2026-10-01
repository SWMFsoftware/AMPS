"""Checksummed clean-source packaging and explicit cross-application evidence.

The source-only protocol tests can qualify this machinery, not an AMPS event.
Production release requires independently supplied application runs, all ten
diagnostics, convergence/MPI evidence and the declared validity/limitations.
Stage-14 research/campaign evidence never becomes a baseline release dependency.
"""
from __future__ import annotations
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import subprocess
import tarfile
from datetime import datetime, timezone
from preprocessing.core import require, digest, freeze_record, verify_frozen, load_json, finite, canonical

PRODUCTS = {"D"+str(n) for n in range(1, 11)}
EVIDENCE_CLASSES = {"verification", "integration", "mpi", "convergence", "cross-model"}
FORBIDDEN_PARTS = {"build", "test_output", "__pycache__", ".git"}
FORBIDDEN_SUFFIXES = {".o", ".a", ".so", ".pyc", ".pyo"}


def file_digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def safe_source(name):
    path = PurePosixPath(name)
    return (not path.is_absolute() and bool(path.parts) and ".." not in path.parts
            and not FORBIDDEN_PARTS.intersection(path.parts)
            and path.suffix not in FORBIDDEN_SUFFIXES)


def source_manifest(root, paths, registry, build_sources, examples):
    """Freeze explicitly enumerated sources; do not search sibling applications.

    Registries declare their owning suite and command. Canonical IDs belong to
    their scientific authority, while manifests bind the actual implementation
    and examples by content. A duplicate local public header is a packaging error.
    """
    root = Path(root).resolve()
    require(len(paths) == len(set(paths)) and bool(paths), "empty/duplicate source manifest")
    files, public, modes = {}, set(), {}
    canonical_public={Path(name).name:name for name in paths if "/include/sep_coronal_cme/" in name}
    for name in sorted(paths):
        require(safe_source(name), "unsafe/generated package path: "+name)
        path = root/name
        require(path.is_file() and not path.is_symlink(), "missing/nonregular source: "+name)
        require(not path.read_bytes().startswith(b"\x7fELF"), "compiled executable in source package")
        require(Path(name).name not in canonical_public or name==canonical_public[Path(name).name],
                "duplicate local public header: "+name)
        if "/include/sep_coronal_cme/" in name:
            basename = Path(name).name
            require(basename not in public, "duplicate local public header")
            public.add(basename)
        files[name] = file_digest(path)
        modes[name] = path.stat().st_mode & 0o777
    require(set(build_sources)|set(examples) <= set(files), "build/example missing from package manifest")
    ids = [row["id"] for row in registry]
    require(ids and len(ids) == len(set(ids)), "empty/duplicate owning registry")
    for row in registry:
        require(row["owner"] in {"shared-model", "srcSEP", "srcSEP3D", "core"}
                and row["full_suite_command"] and row["implementation"] in files,
                "unowned/unselected release test")
    return freeze_record(dict(schema="sccm-source-manifest-v1", files=files, file_modes=modes,
        registry=registry, build_sources=sorted(build_sources), examples=sorted(examples)))


def verify_sources(root, manifest):
    verify_frozen(manifest)
    require(manifest["schema"] == "sccm-source-manifest-v1", "wrong source manifest")
    for name, checksum in manifest["files"].items():
        require(safe_source(name) and file_digest(Path(root)/name) == checksum,
                "missing/changed source: "+name)


def package_sources(root, manifest, archive, embedded_manifest_name=None):
    """Copy only frozen source members and verify the bytes before publication."""
    verify_sources(root, manifest)
    with tarfile.open(archive, "w:gz") as tar:
        for name in sorted(manifest["files"]):
            data = (Path(root)/name).read_bytes()
            member = tarfile.TarInfo(name); member.size = len(data); member.mode = manifest.get("file_modes",{}).get(name,0o644)
            # Preserve source timestamps so overlay extraction does not make
            # changed code appear older than installed object files to make.
            member.mtime=int((Path(root)/name).stat().st_mtime)
            tar.addfile(member, io.BytesIO(data))
        if embedded_manifest_name:
            require(safe_source(embedded_manifest_name) and embedded_manifest_name not in manifest["files"],"unsafe/self-referential embedded manifest")
            data=canonical(manifest);member=tarfile.TarInfo(embedded_manifest_name);member.size=len(data);member.mode=0o644
            tar.addfile(member,io.BytesIO(data))
    verify_archive(archive, manifest, embedded_manifest_name)


def verify_archive(archive, manifest, embedded_manifest_name=None):
    verify_frozen(manifest)
    with tarfile.open(archive) as tar:
        members = tar.getmembers()
        names = [member.name for member in members]
        expected=set(manifest["files"])|({embedded_manifest_name} if embedded_manifest_name else set())
        require(len(names) == len(set(names)) and set(names) == expected,
                "archive omits/adds/duplicates manifest members")
        for member in members:
            require(member.isfile() and safe_source(member.name), "archive contains generated/nonregular member")
            data=tar.extractfile(member).read()
            if member.name==embedded_manifest_name:
                require(data==canonical(manifest),"embedded source manifest drift")
                continue
            require(member.mode==manifest.get("file_modes",{}).get(member.name,0o644),"archive source-mode mismatch")
            require(hashlib.sha256(data).hexdigest() == manifest["files"][member.name],
                    "archive/reference checksum mismatch")


def qualify(profile, evidence):
    """Retain every failing/missing required item; never infer from source presence.

    Observation tolerances are active only after a frozen event registration.
    Protocol/synthetic checks and generic Parker host initialization cannot be
    relabelled observational validation or a coronal runtime qualification.
    """
    verify_frozen(profile); verify_frozen(evidence)
    require(profile["schema"] == "sccm-release-profile-v1", "wrong release profile")
    require(profile["software_stage"] == 13 and profile["stage14_is_release_dependency"] is False,
            "Stage-14 research cannot redefine the baseline release")
    required = profile["required_tests"]
    require(required and len(required) == len(set(required)), "invalid mandatory release registry")
    rows = evidence["test_results"]
    require(len(rows) == len({r["id"] for r in rows}), "duplicate release result")
    known = {r["id"]: r for r in rows}
    tests = [{"id": name, "status": known.get(name, {}).get("status", "MISSING")} for name in required]
    blockers = ["required test "+r["id"]+" is "+r["status"] for r in tests if r["status"] != "PASS"]
    diagnostics = evidence.get("diagnostics", {})
    def verified_asset(record, path_key="path", checksum_key="content_sha256"):
        try:
            checksum = record[checksum_key]
            return (isinstance(checksum, str) and len(checksum) == 64 and
                    file_digest(record[path_key]) == checksum)
        except (KeyError, OSError, TypeError):
            return False
    for name in sorted(PRODUCTS):
        record = diagnostics.get(name)
        if (not record or not verified_asset(record) or
                not verified_asset(record, "provenance_path", "provenance_sha256") or record.get("finite") is not True):
            blockers.append("missing/nonfinite/provenance-incomplete "+name)
    if not EVIDENCE_CLASSES <= set(evidence.get("evidence_classes", [])):
        blockers.append("missing mandatory integration/MPI/convergence/cross-model evidence")
    if evidence.get("clean_tree_builds") != ["sep_common", "sep_coronal_cme", "srcSEP3D", "srcSEP"]:
        blockers.append("four clean-tree builds not established")
    # Booleans or source-presence claims cannot replace machine reports. Each
    # build, application and numerical-evidence class must have immutable
    # bytes; missing hashes/paths leave the release blocked even if all labels
    # claim PASS. The profile freezes the required report identities.
    required_reports = profile.get("required_machine_reports", {})
    mandatory_owners = {"sep_common", "sep_coronal_cme", "srcSEP3D", "srcSEP"}|EVIDENCE_CLASSES
    if not mandatory_owners <= set(required_reports):
        blockers.append("profile lacks mandatory machine-report ownership")
    supplied = evidence.get("machine_reports", {})
    for owner, checksum in required_reports.items():
        record = supplied.get(owner, {})
        if not verified_asset(record) or record.get("content_sha256") != checksum:
            blockers.append("missing/changed machine report "+owner)
        else:
            data = load_json(record["path"])
            if data.get("status") != "PASS" or data.get("evidence_kind") != "production":
                blockers.append("unqualified machine report "+owner)
    if evidence.get("evidence_kind") != "production" or evidence.get("coronal_runtime_qualified") is not True:
        blockers.append("production coronal runtime qualification not established")
    if not evidence.get("domain_of_validity") or not evidence.get("known_limitations"):
        blockers.append("domain of validity/limitations not frozen")
    if profile.get("event_profile"):
        event = evidence.get("event_registration", {})
        if not event.get("frozen_before_observations") or not event.get("dataset_checksums") or not event.get("metrics_and_tolerances"):
            blockers.append("observational acceptance not preregistered/checksummed")
        if profile.get("campaign_c") and (event.get("complete_candidate_table") is not True or event.get("all_n4_caps_pass") is not True):
            blockers.append("campaign-C candidate product/N4 event-grade gate incomplete")
    return freeze_record(dict(schema="sccm-release-evidence-v1", profile_identity=profile["identity"],
        evidence_identity=evidence["identity"], qualified=not blockers, status="PASS" if not blockers else "INCOMPLETE",
        blockers=blockers, tests=tests, classifications=evidence.get("evidence_classes", []),
        observational_validation=bool(profile.get("event_profile") and not blockers),
        domain_of_validity=evidence.get("domain_of_validity"), known_limitations=evidence.get("known_limitations")))


def run_applications(srcsep, srcsep3d, bundle, commands, output, timeout=1800):
    """Use only explicit independent executable/bundle paths and bounded argv.

    Commands are profile-owned argument lists with placeholders, not shell
    strings. Each invocation has a fresh report path; failed/missing reports
    are retained as ERROR. The caller must still qualify the returned evidence.
    """
    paths = {"srcSEP": Path(srcsep).resolve(), "srcSEP3D": Path(srcsep3d).resolve()}
    require(paths["srcSEP"] != paths["srcSEP3D"], "applications need independent executable paths")
    bundle = Path(bundle).resolve()
    require(bundle.is_file(), "explicit immutable bundle is absent")
    finite([timeout]); require(timeout > 0, "invalid application deadline")
    output = Path(output); output.mkdir(parents=True, exist_ok=False)
    results = []
    for application, executable in paths.items():
        report = output/(application+".json"); log = output/(application+".log")
        values = dict(executable=str(executable), bundle=str(bundle), report=str(report))
        require(application in commands and isinstance(commands[application], list), "missing explicit application argv")
        argv = [part.format(**values) for part in commands[application]]
        require(argv and str(executable) in argv and str(bundle) in argv and str(report) in argv,
                "argv omits explicit executable, bundle or fresh report")
        try:
            before = file_digest(executable)
            with log.open("wb") as stream:
                process = subprocess.run(argv, stdout=stream, stderr=subprocess.STDOUT, timeout=timeout)
            data = load_json(report)
            require(file_digest(executable) == before and data["application"] == application,
                    "application identity/executable changed")
            require(data["bundle_sha256"] == file_digest(bundle), "application consumed another bundle")
            require(data["status"] in {"PASS", "FAIL", "ERROR", "SKIP"} and
                    process.returncode == {"PASS": 0, "SKIP": 0, "FAIL": 1, "ERROR": 2}[data["status"]],
                    "application report/exit mismatch")
            results.append(dict(data, log=str(log), executable_sha256=before))
        except (OSError, ValueError, KeyError, subprocess.TimeoutExpired) as error:
            # Missing executables/reports are local phase failures. Retain their
            # diagnosis and continue the independent application's invocation.
            with log.open("ab") as stream: stream.write((str(error)+"\n").encode("utf-8"))
            results.append(dict(application=application, status="ERROR", message=str(error), log=str(log)))
    return freeze_record(dict(schema="sccm-cross-application-run-v1", results=results,
                              bundle_sha256=file_digest(bundle), independently_selected=True))


def update_last_pass(path, result, reference_sha256):
    """Record qualified evidence only; a synthetic machinery PASS is insufficient."""
    verify_frozen(result)
    require(result.get("qualified") is True, "unqualified evidence cannot update last-pass")
    require(isinstance(reference_sha256, str) and len(reference_sha256) == 64 and
            all(c in "0123456789abcdef" for c in reference_sha256), "invalid reference identity")
    record = freeze_record(dict(schema="sccm-last-pass-v1", release_identity=result["identity"],
        reference_sha256=reference_sha256, utc=datetime.now(timezone.utc).isoformat()))
    Path(path).write_text(json.dumps(record, indent=2)+"\n")
    return record
