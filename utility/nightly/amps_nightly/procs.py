"""Running one command for the orchestrator: own session, time limit, log file, and a clean finish.

Every process started for a runner carries the environment variable AMPS_NIGHTLY_TAG=<tag>. Runners
start their own children in new sessions (src/test suites, the Earth runner, mpiexec), so killing the
runner's process group is not enough; `sweep` finds every process of this user whose environment carries
the tag (or a tag below it, "<tag>/...") and kills it.
"""

from __future__ import annotations

import os
import resource
import signal
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence

TAG = "AMPS_NIGHTLY_TAG"
GRACE_S = 60          # after SIGINT at the time limit, before SIGKILL
SAMPLE_S = 30         # memory sampling interval


@dataclass
class Outcome:
    exit_code: Optional[int]     # None: could not be started
    timed_out: bool
    elapsed_s: float
    cpu_s: float
    peak_rss_mb: float
    error: str = ""              # why it could not be started


def _tagged(tag: str) -> List[int]:
    """PIDs of this user's processes whose AMPS_NIGHTLY_TAG is tag or below it."""
    want = (TAG + "=" + tag).encode()
    me, uid, found = os.getpid(), os.getuid(), []
    for entry in os.listdir("/proc"):
        if not entry.isdigit() or int(entry) == me:
            continue
        try:
            if os.stat("/proc/" + entry).st_uid != uid:
                continue
            with open("/proc/{}/environ".format(entry), "rb") as stream:
                variables = stream.read().split(b"\0")
        except OSError:
            continue
        if any(v == want or v.startswith(want + b"/") for v in variables):
            found.append(int(entry))
    return found


def _rss_mb(pids: Sequence[int]) -> float:
    total = 0
    for pid in pids:
        try:
            with open("/proc/{}/statm".format(pid)) as stream:
                total += int(stream.read().split()[1])
        except (OSError, ValueError, IndexError):
            pass
    return total * os.sysconf("SC_PAGE_SIZE") / 2.0 ** 20


def sweep(tag: str, grace_s: float = 5) -> int:
    """Terminate, then kill, every leftover process carrying tag; returns how many there were."""
    pids = _tagged(tag)
    count = len(pids)
    for sig in (signal.SIGTERM, signal.SIGKILL):
        for pid in pids:
            try:
                os.kill(pid, sig)
            except OSError:
                pass
        deadline = time.time() + grace_s
        while sig == signal.SIGTERM and time.time() < deadline and _tagged(tag):
            time.sleep(0.2)
        pids = _tagged(tag)
        if not pids:
            break
    return count


def run(argv: Sequence[str], cwd: Path, env: Dict[str, str], log: Path, timeout_s: float, tag: str,
        header: str = "") -> Outcome:
    """Run argv with stdout and stderr appended to log; stop it at timeout_s; leave no process behind."""
    log.parent.mkdir(parents=True, exist_ok=True)
    env = dict(env)
    env[TAG] = tag
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.time()
    with open(str(log), "ab") as stream:
        stream.write("\n### {} {}\n### cwd {}\n### $ {}\n".format(
            time.strftime("%Y-%m-%d %H:%M:%S"), header, cwd, " ".join(argv)).encode())
        stream.flush()
        try:
            process = subprocess.Popen(list(argv), cwd=str(cwd), env=env, stdin=subprocess.DEVNULL, stdout=stream,
                                       stderr=subprocess.STDOUT, start_new_session=True)
        except OSError as error:
            stream.write("### could not start: {}\n".format(error).encode())
            return Outcome(None, False, 0.0, 0.0, 0.0, str(error))
        timed_out, peak, next_sample = False, 0.0, started
        while process.poll() is None:
            now = time.time()
            if now >= next_sample:
                peak = max(peak, _rss_mb(_tagged(tag)))
                next_sample = now + SAMPLE_S
            if now - started > timeout_s and not timed_out:
                timed_out = True
                stream.write("\n### time limit of {:.0f} s reached: SIGINT, then SIGKILL after {} s\n".format(
                    timeout_s, GRACE_S).encode())
                stream.flush()
                try:
                    os.killpg(process.pid, signal.SIGINT)
                except OSError:
                    pass
                try:
                    process.wait(GRACE_S)
                except subprocess.TimeoutExpired:
                    pass
                sweep(tag, grace_s=0)
            time.sleep(0.5)
        leftovers = sweep(tag)
        elapsed = time.time() - started
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        cpu = (after.ru_utime - before.ru_utime) + (after.ru_stime - before.ru_stime)
        stream.write("### exit {} after {:.0f} s{}{}\n".format(
            process.returncode, elapsed, " (time limit)" if timed_out else "",
            ", {} leftover processes killed".format(leftovers) if leftovers else "").encode())
    return Outcome(process.returncode, timed_out, elapsed, cpu, peak)


def shell(command: str) -> List[str]:
    return ["/bin/bash", "-c", command]
