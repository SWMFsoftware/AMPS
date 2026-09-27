#!/usr/bin/env python3
"""Audit and package the AMPS SEP source components from explicit allowlists.

The checker is intentionally located at the AMPS level because ``srcSEP`` and
``srcSEP3D`` share code from ``src/models/sep_common`` and
``src/models/swcme``.  Auditing either application in isolation could accept a
stale shared archive or omit a required shared source file.  Each component
owns a ``SOURCE_MANIFEST.json``; this program interprets all four manifests in
one transaction and refuses to package any unclassified, generated, retired,
or unsafe path.

The default operation is read-only.  ``--create-archive`` writes a deterministic
gzip-compressed tar archive only after the complete audit passes, then reopens
the archive and verifies its member set against the same allowlists.  Runtime
test evidence is deliberately not part of that archive; it belongs in a
separate artifact tied to the archive SHA-256.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import gzip
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import sys
import tarfile
import tempfile
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple


SCHEMA = "amps-sep-source-manifest-v1"
COMPONENTS: Tuple[str, ...] = (
    "srcSEP",
    "srcSEP3D",
    "src/models/sep_common",
    "src/models/swcme",
)
ALLOWLIST_FIELDS: Tuple[str, ...] = (
    "production_source_globs",
    "test_source_globs",
    "documentation_globs",
    "retained_input_globs",
    "support_globs",
)

# These AMPS-level files cannot be named by an application-local manifest.
# Keep the list deliberately small: adding a new top-level deliverable must be
# a reviewed release-policy change rather than an accidental archive side
# effect.
TOP_LEVEL_ALLOWED: Tuple[str, ...] = (
    # This application deck is a required source input to AMPS' compile-time
    # generator.  Omitting it lets an installation retain the historical
    # ELECTRON selection even when all srcSEP3D runtime files are updated.
    "input/sep3d.input",
    # These core guides document the application-facing associated-data and
    # output callback ABI used by the packaged SEP applications.  They live
    # outside every application-owned manifest, so their inclusion must remain
    # an explicit AMPS-level release-policy decision.
    "src/README.md",
    "src/pic/README.md",
    # srcSEP3D's population controller deliberately disables AMPS' automatic
    # splitter, but it still shares the public linked-list split/merge helpers.
    # Ship the guarded implementation with the application sources so an
    # overlay cannot retain an older core file containing empty/singleton
    # undefined-behaviour paths while the integration gate expects the fix.
    "src/pic/pic_particle_spliting.cpp",
    "tools/sep_package_hygiene.py",
    "STAGE1_STAGE2_IMPLEMENTATION.md",
    "STAGE3_IMPLEMENTATION.md",
)

# An allowlist and a required-member list answer different release questions.
# ``TOP_LEVEL_ALLOWED`` says which AMPS-level files may enter an SEP source
# archive when they exist.  ``TOP_LEVEL_REQUIRED`` says which of those files
# must exist for the archive to be usable.  Keep the guarded split/merge source
# in both sets: srcSEP3D calls that AMPS-core implementation, and accepting an
# archive without it would leave an older, unsafe destination copy active after
# an overlay installation.
TOP_LEVEL_REQUIRED: Tuple[str, ...] = (
    "src/pic/pic_particle_spliting.cpp",
)

if not set(TOP_LEVEL_REQUIRED).issubset(TOP_LEVEL_ALLOWED):
    raise RuntimeError(
        "every required AMPS-level release member must also be allowlisted")


@dataclass(frozen=True)
class Issue:
    """One stable audit finding suitable for terminal and CI output."""

    category: str
    path: str
    detail: str


@dataclass(frozen=True)
class AuditResult:
    """Complete source audit and the exact file set eligible for packaging."""

    allowed_files: Tuple[str, ...]
    issues: Tuple[Issue, ...]


class HygieneError(RuntimeError):
    """Describe malformed policy or an unsafe archive operation."""


def _read_manifest(root: Path, component: str) -> Mapping[str, object]:
    """Load and validate one component policy before inspecting its files."""
    path = root / component / "SOURCE_MANIFEST.json"
    try:
        with path.open("r", encoding="utf-8") as stream:
            value = json.load(stream)
    except (OSError, json.JSONDecodeError) as error:
        raise HygieneError(f"cannot read source manifest {path}: {error}") from error
    if not isinstance(value, dict) or value.get("schema") != SCHEMA:
        raise HygieneError(f"unsupported or missing manifest schema in {path}")
    if value.get("component") != component:
        raise HygieneError(
            f"manifest component mismatch in {path}: expected {component!r}, "
            f"found {value.get('component')!r}")
    for field in (*ALLOWLIST_FIELDS, "generated_globs", "retired_sources"):
        entries = value.get(field)
        if not isinstance(entries, list) or any(not isinstance(item, str) or not item
                                                for item in entries):
            raise HygieneError(f"{path}: {field} must be an array of nonempty strings")
    return value


def _matches(path: str, patterns: Iterable[str]) -> bool:
    """Match component-relative POSIX paths with manifest glob semantics.

    ``PurePosixPath.match`` handles the ``**`` patterns used by the manifests
    but treats a trailing ``/**`` as requiring at least one child.  The prefix
    branch below makes that intent explicit for all descendants while still
    keeping ordinary ``*.cpp`` patterns directory-local.
    """
    for pattern in patterns:
        normalized = pattern.strip("/")
        if normalized.endswith("/**"):
            prefix = normalized[:-3].rstrip("/")
            if path == prefix or path.startswith(prefix + "/"):
                return True
        # ``pathlib.Path.match`` historically treats ``**`` differently
        # across supported Python releases.  Translate the small manifest glob
        # language directly so ``*`` never crosses a directory boundary while
        # ``**`` always can.  This keeps ``*.cpp`` local to the component root
        # and makes ``test/**/*.py`` cover both direct and nested test files.
        expression = ""
        index = 0
        while index < len(normalized):
            character = normalized[index]
            if character == "*" and index + 1 < len(normalized) and normalized[index + 1] == "*":
                index += 2
                if index < len(normalized) and normalized[index] == "/":
                    expression += "(?:.*/)?"
                    index += 1
                else:
                    expression += ".*"
                continue
            if character == "*":
                expression += "[^/]*"
            elif character == "?":
                expression += "[^/]"
            else:
                expression += re.escape(character)
            index += 1
        if re.fullmatch(expression, path) is not None:
            return True
    return False


def _component_files(component_root: Path) -> Iterable[Path]:
    """Yield files and symlinks without following directory symlink escapes."""
    for directory, subdirectories, filenames in os.walk(component_root,
                                                          followlinks=False):
        directory_path = Path(directory)
        # Directory symlinks are package members too.  Remove them from the
        # traversal so an external target can never enlarge the audited tree.
        retained_subdirectories: List[str] = []
        for name in subdirectories:
            child = directory_path / name
            if child.is_symlink():
                yield child
            else:
                retained_subdirectories.append(name)
        subdirectories[:] = retained_subdirectories
        for name in filenames:
            yield directory_path / name


def _safe_symlink(root: Path, path: Path) -> bool:
    """Accept only relative symlinks whose resolved target stays below root."""
    if not path.is_symlink():
        return True
    target = os.readlink(path)
    if os.path.isabs(target):
        return False
    try:
        (path.parent / target).resolve().relative_to(root.resolve())
    except (OSError, ValueError):
        return False
    return True


def audit(root: Path) -> AuditResult:
    """Classify every file in all four SEP components.

    Generated and retired patterns take precedence over allowlists.  This is
    important when a broad documentation or support pattern would otherwise
    hide a stale product.  An unclassified file is also a release failure: the
    author must either add it deliberately to a manifest or remove it.
    """
    root = root.resolve()
    allowed: Set[str] = set()
    issues: List[Issue] = []
    for component in COMPONENTS:
        component_root = root / component
        if not component_root.is_dir():
            issues.append(Issue("missing-component", component,
                                "required SEP source component is absent"))
            continue
        manifest = _read_manifest(root, component)
        allow_patterns = [str(item) for field in ALLOWLIST_FIELDS
                          for item in manifest[field]]
        generated_patterns = [str(item) for item in manifest["generated_globs"]]
        retired = {str(item).strip("/") for item in manifest["retired_sources"]}

        for path in _component_files(component_root):
            relative_component = path.relative_to(component_root).as_posix()
            relative_root = path.relative_to(root).as_posix()
            if not _safe_symlink(root, path):
                issues.append(Issue("unsafe-symlink", relative_root,
                                    "symlink target leaves the audited AMPS root"))
            elif relative_component in retired:
                issues.append(Issue("retired", relative_root,
                                    "retired source is physically present"))
            elif _matches(relative_component, generated_patterns):
                issues.append(Issue("generated", relative_root,
                                    "generated product is forbidden in a source release"))
            elif _matches(relative_component, allow_patterns):
                allowed.add(relative_root)
            else:
                issues.append(Issue("unclassified", relative_root,
                                    "path matches no release allowlist category"))

    for relative in TOP_LEVEL_ALLOWED:
        path = root / relative
        if path.exists() or path.is_symlink():
            if _safe_symlink(root, path):
                allowed.add(relative)
            else:
                issues.append(Issue("unsafe-symlink", relative,
                                    "top-level symlink target leaves AMPS root"))
        elif relative in TOP_LEVEL_REQUIRED:
            # Absence is a release error, not permission to silently omit the
            # file.  This distinction is essential for overlay archives: a
            # missing source member does not delete a stale destination file.
            issues.append(Issue(
                "missing-required", relative,
                "required AMPS-level release member is absent"))

    return AuditResult(tuple(sorted(allowed)),
                       tuple(sorted(issues, key=lambda item: (item.path,
                                                              item.category))))


def _tar_filter(info: tarfile.TarInfo) -> tarfile.TarInfo:
    """Remove host identity and timestamps from a release archive member."""
    info.uid = 0
    info.gid = 0
    info.uname = ""
    info.gname = ""
    info.mtime = 0
    # Group/other write bits are never useful in immutable source releases.
    info.mode &= ~0o022
    return info


def create_archive(root: Path, output: Path, files: Sequence[str]) -> str:
    """Write a deterministic archive and return its SHA-256 digest."""
    root = root.resolve()
    output = output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    try:
        # Supplying an empty gzip filename avoids embedding a workstation path
        # in the header; mtime=0 makes identical sources byte reproducible.
        with temporary.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
                with tarfile.open(fileobj=compressed, mode="w", format=tarfile.PAX_FORMAT) as archive:
                    for relative in files:
                        archive.add(root / relative, arcname=relative,
                                    recursive=False, filter=_tar_filter)
        os.replace(temporary, output)
    finally:
        if temporary.exists():
            temporary.unlink()
    verify_archive(output, set(files))
    digest = hashlib.sha256()
    with output.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_archive(path: Path, expected_files: Optional[Set[str]] = None) -> None:
    """Reject traversal, duplicate, special, or unexpected archive members."""
    try:
        with tarfile.open(path, "r:gz") as archive:
            members = archive.getmembers()
    except (OSError, tarfile.TarError) as error:
        raise HygieneError(f"cannot read candidate archive {path}: {error}") from error
    names: List[str] = []
    for member in members:
        name = PurePosixPath(member.name).as_posix()
        if name.startswith("/") or ".." in PurePosixPath(name).parts:
            raise HygieneError(f"unsafe archive member path: {member.name}")
        if not (member.isfile() or member.issym()):
            raise HygieneError(f"unsupported archive member type: {member.name}")
        names.append(name)
    if len(names) != len(set(names)):
        raise HygieneError("candidate archive contains duplicate member names")
    if expected_files is not None and set(names) != expected_files:
        missing = sorted(expected_files - set(names))
        extra = sorted(set(names) - expected_files)
        raise HygieneError(f"archive member mismatch: missing={missing}, extra={extra}")


def _write_fixture_manifest(path: Path, component: str) -> None:
    """Create a minimal valid manifest used only in an isolated self-test."""
    value = {
        "schema": SCHEMA,
        "component": component,
        "production_source_globs": ["*.cpp"],
        "test_source_globs": [],
        "documentation_globs": [],
        "retained_input_globs": [],
        "support_globs": ["SOURCE_MANIFEST.json"],
        "generated_globs": ["*.o", "test_output/**", "**/__pycache__/**",
                            "**/*.pyc"],
        "retired_sources": ["mover.cpp"],
    }
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def run_self_test() -> None:
    """Exercise positive packaging and every required negative control."""
    with tempfile.TemporaryDirectory(prefix="sep-package-hygiene-") as temporary:
        root = Path(temporary) / "AMPS"
        for component in COMPONENTS:
            directory = root / component
            directory.mkdir(parents=True)
            _write_fixture_manifest(directory / "SOURCE_MANIFEST.json", component)
            (directory / "main.cpp").write_text("// reviewed source fixture\n",
                                                encoding="utf-8")
        tool = root / "tools" / "sep_package_hygiene.py"
        tool.parent.mkdir(parents=True)
        tool.write_text("# fixture tool\n", encoding="utf-8")

        # The production policy requires the guarded AMPS split/merge source.
        # Model that external member explicitly so the positive fixture and the
        # missing-member negative control exercise the same release contract.
        required_core = root / "src" / "pic" / "pic_particle_spliting.cpp"
        required_core.parent.mkdir(parents=True)
        required_core.write_text("// guarded core fixture\n", encoding="utf-8")

        baseline = audit(root)
        if baseline.issues:
            raise HygieneError(f"positive self-test failed: {baseline.issues}")

        # Prove that a required member cannot regress to the historical
        # "allowed when present" behavior.  Removing the core source must make
        # the audit fail before any archive is constructed.
        required_core.unlink()
        required_findings = audit(root).issues
        if not any(issue.category == "missing-required" and
                   issue.path == "src/pic/pic_particle_spliting.cpp"
                   for issue in required_findings):
            raise HygieneError(
                "missing required AMPS core source did not fail the audit")
        required_core.write_text("// guarded core fixture\n", encoding="utf-8")

        controls = (
            ("stale.o", "generated"),
            ("mover.cpp", "retired"),
            ("__pycache__/cached.pyc", "generated"),
            ("test_output/result.json", "generated"),
            ("unreviewed.dat", "unclassified"),
        )
        component = root / "srcSEP"
        for relative, category in controls:
            injected = component / relative
            injected.parent.mkdir(parents=True, exist_ok=True)
            injected.write_bytes(b"negative-control\n")
            findings = audit(root).issues
            if not any(issue.category == category and
                       issue.path == f"srcSEP/{relative}" for issue in findings):
                raise HygieneError(
                    f"negative control {relative} did not produce {category}")
            injected.unlink()
            parent = injected.parent
            while parent != component and not any(parent.iterdir()):
                parent.rmdir()
                parent = parent.parent

        archive = Path(temporary) / "fixture.tar.gz"
        digest = create_archive(root, archive, baseline.allowed_files)
        if len(digest) != 64:
            raise HygieneError("archive self-test did not produce a SHA-256 digest")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Audit/package srcSEP, srcSEP3D, sep_common, and SWCME source.")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1],
                        help="AMPS root containing the four manifest-owned components")
    parser.add_argument("--self-test", action="store_true",
                        help="also run isolated positive and negative controls")
    parser.add_argument("--list-files", action="store_true",
                        help="print every allowlisted archive member")
    parser.add_argument("--create-archive", type=Path,
                        help="write and verify a deterministic source tar.gz")
    parser.add_argument("--verify-archive", type=Path,
                        help="verify archive paths/types without creating an archive")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    try:
        result = audit(args.root)
        if result.issues:
            for issue in result.issues:
                print(f"ERROR [{issue.category}] {issue.path}: {issue.detail}",
                      file=sys.stderr)
            print(f"FAIL SEP-PACKAGE-HYGIENE: {len(result.issues)} issue(s)",
                  file=sys.stderr)
            return 1
        if args.list_files:
            for relative in result.allowed_files:
                print(relative)
        if args.self_test:
            run_self_test()
        if args.verify_archive:
            verify_archive(args.verify_archive)
        if args.create_archive:
            digest = create_archive(args.root, args.create_archive,
                                    result.allowed_files)
            print(f"archive={args.create_archive.expanduser().resolve()}")
            print(f"sha256={digest}")
        print(f"PASS SEP-PACKAGE-HYGIENE: {len(result.allowed_files)} "
              "allowlisted file(s); zero generated, retired, or unclassified paths")
        return 0
    except HygieneError as error:
        print(f"ERROR SEP-PACKAGE-HYGIENE: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
