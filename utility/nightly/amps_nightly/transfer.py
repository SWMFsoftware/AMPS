"""Moving results from a test machine to the results server: one tarball per run (one column, one night).

On the test machine, after `amps-nightly run` (script 2: nightly_pack.sh = `amps-nightly pack`):

    {workdir}/nights/<night>/runs/<column>/     run.json results.jsonl runners.json logs/... DONE
        -> {workdir}/outbox/<run_id>.tar.gz     everything the server needs for the pages and the archive
        -> the machine file's "upload": a local directory (the server itself) or [user@]host:directory (rsync)

On the server (script 3: nightly_update.sh = `amps-nightly update`), store.ingest unpacks every
incoming/<run_id>.tar.gz with `unpack`, which checks every member first, and archives the run.

A run whose orchestrator died (reboot, out of memory, killed) has no DONE. `finalize` completes it from
the progress the orchestrator writes after every runner and marks it incomplete, so that night still
reaches the report.
"""

from __future__ import annotations

import datetime as dt
import fcntl
import json
import os
import re
import shutil
import subprocess
import tarfile
import zlib
from pathlib import Path, PurePosixPath
from typing import Callable, Dict, List, Optional, Tuple

from .records import DONE_FILE, RecordError, read_run_dir, unique_rows, utc_text, write_run_dir
from .registry import UNASSIGNED

SUFFIX = ".tar.gz"
PACKED_FILE = "PACKED"                  # in a run directory: it is in the outbox
PARTIAL_RUN_FILE = "run.partial.json"   # run.json as known when the run started (for finalize)
PROGRESS_FILE = "progress.jsonl"        # rows so far, rewritten after every runner
RUNNING_FILE = "running"                # in the workdir: name of the night `amps-nightly run` is executing
LOCAL_ONLY = (PACKED_FILE, PARTIAL_RUN_FILE, PROGRESS_FILE)
KEEP_SENT_DAYS = 14
MAX_MEMBERS = 100000
MAX_UNPACKED_BYTES = 8 * 2 ** 30
_TOP = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._+-]{0,127}$")


class TransferError(Exception):
    pass


# ----------------------------------------------------------------------------- test machine
def running_night(workdir: Path) -> Optional[str]:
    """The night an `amps-nightly run` is executing now ("?" if unknown), or None if none is running."""
    lock = Path(workdir) / "lock"
    if not lock.exists():
        return None
    with open(str(lock), "a") as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            try:
                return (Path(workdir) / RUNNING_FILE).read_text().strip() or "?"
            except OSError:
                return "?"
        fcntl.flock(stream, fcntl.LOCK_UN)
    return None


def runs_to_pack(workdir: Path) -> Tuple[List[Path], List[Path]]:
    """(finished run directories not packed yet, unfinished ones of nights that no run is executing)."""
    busy = running_night(workdir)
    finished, abandoned = [], []
    for rundir in sorted((Path(workdir) / "nights").glob("*/runs/*")):
        if not rundir.is_dir() or (rundir / PACKED_FILE).exists():
            continue
        if (rundir / DONE_FILE).is_file():
            finished.append(rundir)
        elif busy is None or (busy != "?" and rundir.parent.parent.name != busy):
            abandoned.append(rundir)
    return finished, abandoned


def _interrupted(night_log: Path, column: str, finished: set) -> str:
    """The runner the orchestrator was running for this column when it stopped (from orchestrator.log)."""
    try:
        lines = night_log.read_text(errors="replace").splitlines()
    except OSError:
        return ""
    current, started = "", ""
    for line in lines:
        text = line.split(" ", 1)[1] if " " in line else line
        if text.startswith("--- "):
            current = text[4:].split(" ", 1)[0]
        elif current == column:
            match = re.match(r"^([A-Za-z0-9_.-]+): start \(limit", text)
            if match:
                started = match.group(1)
    return started if started and started not in finished else ""


def finalize(rundir: Path, family_of: Callable[[str], str] = lambda runner: "") -> Dict:
    """Complete a run directory without DONE (its orchestrator died) and mark the run incomplete."""
    rundir = Path(rundir)
    try:
        run = json.loads((rundir / PARTIAL_RUN_FILE).read_text())
    except (OSError, ValueError) as error:
        raise TransferError("{} has no DONE and no readable {} ({}); it cannot be completed".format(
            rundir, PARTIAL_RUN_FILE, error))
    progress = rundir / PROGRESS_FILE
    try:
        rows = [json.loads(l) for l in progress.read_text().splitlines() if l.strip()] if progress.is_file() else []
        runners = json.loads((rundir / "runners.json").read_text()) if (rundir / "runners.json").is_file() else []
    except (OSError, ValueError) as error:
        raise TransferError("cannot read the progress of {}: {}".format(rundir, error))
    night_log = rundir.parent.parent / "orchestrator.log"
    (rundir / "logs").mkdir(parents=True, exist_ok=True)
    if night_log.is_file():
        shutil.copyfile(str(night_log), str(rundir / "logs" / "orchestrator.log"))
    why = "the orchestrator stopped before the night finished (the machine rebooted or its process was killed)"
    runner = _interrupted(night_log, rundir.name, {r.get("runner") for r in runners})
    if runner:
        log = "logs/{}/runner.log".format(runner)
        rows.append({"family": family_of(runner) or UNASSIGNED, "test": runner, "status": "FAILED", "stage": "run",
                     "reason": "{}: {} was running; its results are missing".format(why, runner), "source": runner,
                     "log": log if (rundir / log).is_file() else ""})
        runners.append({"runner": runner, "state": "interrupted", "problem": why, "log": log})
        (rundir / "runners.json").write_text(json.dumps(runners, indent=2) + "\n")
    stamps = [p.stat().st_mtime for p in (progress, night_log, rundir / PARTIAL_RUN_FILE) if p.is_file()]
    run.update(complete=False, finished=utc_text(dt.datetime.fromtimestamp(max(stamps), dt.timezone.utc)),
               notes="; ".join(x for x in (run.get("notes", ""), why) if x))
    write_run_dir(rundir, run, unique_rows(rows), done=True)
    return run


def pack(rundir: Path, outbox: Path) -> Path:
    """Write <outbox>/<run_id>.tar.gz from a finished run directory (DONE goes last) and mark it PACKED."""
    rundir, outbox = Path(rundir), Path(outbox)
    try:
        run, _ = read_run_dir(rundir)              # never send what the server would reject
    except RecordError as error:
        raise TransferError("{}: {}".format(rundir, error))
    outbox.mkdir(parents=True, exist_ok=True)
    name = run["run_id"] + SUFFIX
    final, part = outbox / name, outbox / ("." + name + ".part")
    files = [p for p in rundir.rglob("*") if p.is_file() and not p.is_symlink()
             and not (p.parent == rundir and p.name in LOCAL_ONLY)]
    files.sort(key=lambda p: (p == rundir / DONE_FILE, str(p)))
    with tarfile.open(str(part), "w:gz", compresslevel=6) as tar:
        for path in files:
            info = tar.gettarinfo(str(path), arcname="{}/{}".format(run["run_id"], path.relative_to(rundir).as_posix()))
            info.uid = info.gid = 0
            info.uname = info.gname = ""
            info.mode = 0o644
            with open(str(path), "rb") as stream:
                tar.addfile(info, stream)
    os.replace(str(part), str(final))
    sent = outbox / (name + ".sent")
    if sent.exists():
        sent.unlink()
    (rundir / PACKED_FILE).write_text(str(final) + "\n")
    return final


def is_remote(target: str) -> bool:
    return ":" in target.split("/")[0] and not target.startswith(("/", "~"))


def send(tarball: Path, target: str, timeout_s: int = 3600) -> None:
    """Deliver a tarball to the server's incoming directory; the file appears there only when complete."""
    tarball = Path(tarball)
    if not target:
        raise TransferError("the machine file has no \"upload\" destination")
    try:
        if is_remote(target):
            # rsync writes into a hidden temporary file and renames it at the end
            subprocess.run(["rsync", "-t", "--chmod=F644", str(tarball), target.rstrip("/") + "/"], check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout_s)
        else:
            incoming = Path(target).expanduser()
            incoming.mkdir(parents=True, exist_ok=True)
            part = incoming / ("." + tarball.name + ".part")
            shutil.copyfile(str(tarball), str(part))
            os.replace(str(part), str(incoming / tarball.name))
    except subprocess.CalledProcessError as error:
        raise TransferError("rsync to {} failed: {}".format(target, (error.stderr or b"").decode(errors="replace").strip()))
    except (OSError, subprocess.SubprocessError) as error:
        raise TransferError("sending to {} failed: {}".format(target, error))
    (tarball.parent / (tarball.name + ".sent")).write_text("{} {}\n".format(
        utc_text(dt.datetime.now(dt.timezone.utc)), target))


def unsent(outbox: Path) -> List[Path]:
    outbox = Path(outbox)
    if not outbox.is_dir():
        return []
    return sorted(p for p in outbox.glob("*" + SUFFIX)
                  if not p.name.startswith(".") and not (outbox / (p.name + ".sent")).exists())


def clean_outbox(outbox: Path, keep_days: int = KEEP_SENT_DAYS) -> int:
    """Delete tarballs that were sent more than keep_days ago; unsent ones stay until they are sent."""
    outbox, removed = Path(outbox), 0
    limit = dt.datetime.now().timestamp() - keep_days * 86400
    for marker in (outbox.glob("*" + SUFFIX + ".sent") if outbox.is_dir() else []):
        if marker.stat().st_mtime < limit:
            tarball = outbox / marker.name[:-len(".sent")]
            if tarball.exists():
                tarball.unlink()
            marker.unlink()
            removed += 1
    return removed


# ----------------------------------------------------------------------------- server
def unpack(tarball: Path, staging: Path) -> Path:
    """Extract a run tarball into staging after checking every member; return the run directory.

    Only regular files and directories below one top-level directory are accepted: no absolute paths,
    no "..", no links or devices, at most MAX_MEMBERS members and MAX_UNPACKED_BYTES in all.
    """
    tarball, staging = Path(tarball), Path(staging)
    shutil.rmtree(str(staging), ignore_errors=True)
    staging.mkdir(parents=True)
    try:
        with tarfile.open(str(tarball), "r:gz") as tar:
            members, top, total = [], None, 0
            for member in tar:
                members.append(member)
                if len(members) > MAX_MEMBERS:
                    raise TransferError("more than {} members".format(MAX_MEMBERS))
                path = PurePosixPath(member.name)
                parts = path.parts
                if not parts or path.is_absolute() or ".." in parts:
                    raise TransferError("unsafe member name {!r}".format(member.name))
                if not (member.isfile() or member.isdir()):
                    raise TransferError("{} is not a regular file or directory".format(member.name))
                if member.isfile() and len(parts) < 2:
                    raise TransferError("{} is not inside the run directory".format(member.name))
                top = top or parts[0]
                if parts[0] != top:
                    raise TransferError("more than one top-level directory ({}, {})".format(top, parts[0]))
                total += member.size
                if total > MAX_UNPACKED_BYTES:
                    raise TransferError("more than {} GB unpacked".format(MAX_UNPACKED_BYTES // 2 ** 30))
            if top is None or not _TOP.match(top):
                raise TransferError("no run directory with a valid name inside")
            for member in members:
                target = staging.joinpath(*PurePosixPath(member.name).parts)
                if member.isdir():
                    target.mkdir(parents=True, exist_ok=True)
                    continue
                target.parent.mkdir(parents=True, exist_ok=True)
                source = tar.extractfile(member)
                with open(str(target), "wb") as stream:
                    shutil.copyfileobj(source, stream)
                os.utime(str(target), (member.mtime, member.mtime))
    except (tarfile.TarError, EOFError, zlib.error, OSError) as error:
        raise TransferError("cannot unpack {}: {}".format(tarball.name, error))
    return staging / top
