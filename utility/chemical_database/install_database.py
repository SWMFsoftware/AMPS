#!/usr/bin/env python3
"""Install a verified chemical publication into a new/empty explicit directory.

The installer never replaces an existing installation or removes unknown user
files. It verifies staged data/software first, copies with exclusive file
creation, verifies the destination, and publishes the exact lock last. Source
archives and the staged snapshot remain untouched. Runtime uses only the
installed relative paths; no symlink to the development checkout is created.
"""

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent / "python"))
from amps_chemical_database import Database, DatabaseError
from amps_chemical_database.api import member, sha256, strict_json


def _verify_software(root):
    manifest = strict_json(root / "software_manifest.json")
    if manifest.get("schema") != "amps-chemical-software-manifest-v1":
        raise DatabaseError("UnsupportedSchema", "unrecognized installed software manifest")
    for relative, digest in manifest["files"].items():
        if sha256(member(root, relative)) != digest:
            raise DatabaseError("CorruptSource", "installed software/documentation changed: " + relative)
    return len(manifest["files"])


def _exclusive_copy(source, target):
    target.parent.mkdir(parents=True, exist_ok=True)
    with source.open("rb") as incoming, target.open("xb") as outgoing:
        shutil.copyfileobj(incoming, outgoing, 1024 * 1024)
    if sha256(source) != sha256(target):
        raise DatabaseError("CorruptSource", "destination copy failed its exact hash: " + str(target))
    # Published source/data files are read-only by default. The CLI launcher
    # remains executable; runtime disables bytecode writes before imports.
    target.chmod(0o555 if source.name == "chemdb" and source.parent.name == "bin" else 0o444)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--staged-root", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True,
                        help="new completion receipt outside the immutable installed archive")
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args(argv)
    source = args.staged_root.resolve()
    # Do not resolve away a destination symlink before explicitly refusing it.
    if args.destination.is_symlink():
        raise DatabaseError("InvalidInput", "installation destination must not be a symlink")
    target = args.destination.resolve()
    receipt_path = args.receipt.resolve()
    marker = target / ".installation-in-progress.json"
    lock_written = False
    receipt = {"schema": "amps-chemical-installation-v1", "started_utc": datetime.now(timezone.utc).isoformat(),
               "staged_root": str(source), "destination": str(target), "status": "ERROR",
               "source_archives_modified": False, "physical_activation": False}
    try:
        if target == source or target in source.parents or source in target.parents or "raw" in target.parts:
            raise DatabaseError("InvalidInput", "installation target overlaps source or immutable raw storage")
        if receipt_path == target or target in receipt_path.parents:
            raise DatabaseError("InvalidInput", "installation receipt must be outside immutable destination")
        if receipt_path.exists():
            raise DatabaseError("InvalidInput", "refusing to overwrite an existing receipt")
        if target.exists() and (not target.is_dir() or any(target.iterdir())):
            raise DatabaseError("InvalidInput", "installation requires a new/empty destination; existing contents are preserved")
        if strict_json(source / "publication.json").get("publication_status") != "READY":
            raise DatabaseError("Unavailable", "staged publication is incomplete")
        if any(path.is_symlink() for path in source.rglob("*")):
            raise DatabaseError("InvalidInput", "staged publication contains a nonportable symlink")
        database = Database.from_lock(source / "database-lock.json")
        receipt["staged_integrity"] = database.verify(deep=True)
        receipt["software_files_checked"] = _verify_software(source)
        receipt["lock_sha256"] = sha256(source / "database-lock.json")
        if args.check_only:
            receipt.update(status="PASS", scope="preinstallation verification only; no destination write")
        else:
            target.mkdir(parents=True, exist_ok=True)
            with marker.open("x", encoding="utf-8") as output:
                json.dump(receipt, output, indent=2)
            files = sorted(path for path in source.rglob("*") if path.is_file())
            for path in files:
                relative = path.relative_to(source)
                if relative.as_posix() == "database-lock.json":
                    continue
                _exclusive_copy(path, target / relative)
            pins = strict_json(source / "database-lock.json")
            installed = Database.open(target, pins["release_id"], pins["release_manifest_sha256"],
                                      pins["profile_id"], pins["profile_sha256"],
                                      pins["export_path"], pins["export_manifest_sha256"])
            receipt["installed_integrity"] = installed.verify(deep=True)
            _verify_software(target)
            _exclusive_copy(source / "database-lock.json", target / "database-lock.json")
            lock_written = True
            command = [sys.executable, "-B", str(target / "bin/chemdb"),
                       "--lock", str(target / "database-lock.json"), "verify"]
            result = subprocess.run(command, capture_output=True, text=True)
            receipt["installed_cli"] = {"command": command, "returncode": result.returncode,
                                        "stdout": result.stdout, "stderr": result.stderr}
            if result.returncode:
                raise DatabaseError("NumericalFailure", "fresh installed CLI failed integrity verification")
            marker.unlink()
            receipt.update(status="PASS", scope="complete byte-verified installation and fresh installed CLI execution",
                           installed_file_count=len(files))
    except (DatabaseError, OSError, ValueError, KeyError) as error:
        if lock_written:
            # Only this invocation's newly published lock is removed. The
            # installer refuses pre-existing contents before any mutation.
            (target / "database-lock.json").unlink()
        receipt["error"] = {"code": getattr(error, "code", "InfrastructureError"), "message": str(error)}
    receipt["finished_utc"] = datetime.now(timezone.utc).isoformat()
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    with receipt_path.open("x", encoding="utf-8") as output:
        json.dump(receipt, output, indent=2)
        output.write("\n")
    print(json.dumps(receipt, indent=2))
    return 0 if receipt["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
