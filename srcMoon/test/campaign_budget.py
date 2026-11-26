"""One monotonic campaign budget, including nested MPI/build process groups.

This module controls execution only. A timeout is ERROR and cannot change a
science acceptance criterion. Linux /proc identities protect against PID reuse;
all process groups observed below the launched child are retained even if a
launcher exits before its MPI ranks. No unrelated host process is signalled.
"""
from __future__ import annotations

import json
import ctypes
import math
import os
from pathlib import Path
import signal
import shutil
import subprocess
import time


def process_table():
    """Return pid -> (parent, group, start ticks, state), never command/env data."""
    result = {}
    for path in Path('/proc').iterdir():
        if not path.name.isdecimal():
            continue
        try:
            fields = (path/'stat').read_text().rsplit(')', 1)[1].split()
            result[int(path.name)] = (int(fields[1]), int(fields[2]),
                                     int(fields[19]), fields[0])
        except (OSError, ValueError, IndexError):
            continue
    return result


def descendants(table, root, known=None):
    """Retain verified descendant identities across reparenting and new sessions."""
    found = {} if known is None else dict(known)
    if root in table:
        found[root] = table[root][2]
    changed = True
    while changed:
        changed = False
        for pid, (parent, _, identity, state) in table.items():
            if pid not in found and parent in found and table.get(parent, (0,0,-1))[2] == found[parent]:
                found[pid] = identity
                changed = True
    return {pid: birth for pid, birth in found.items()
            if pid in table and table[pid][2] == birth and table[pid][3] != 'Z'}


def terminate_tree(process, known, grace_s=5.0):
    """Reap nested sessions, then the direct child, within a finite grace period."""
    def send(sig):
        table = process_table()
        live = descendants(table, process.pid, known)
        # Signal individual verified descendants as well as their groups: an
        # mpiexec/make may create a session outside the worker's original group.
        for pid, birth in live.items():
            if pid == os.getpid():
                continue
            try:
                os.kill(pid, sig)
            except ProcessLookupError:
                pass
        return live
    live = send(signal.SIGTERM)
    stop = time.monotonic() + grace_s
    while live and time.monotonic() < stop:
        time.sleep(0.05)
        live = descendants(process_table(), process.pid, live)
    if live:
        send(signal.SIGKILL)
    try:
        process.wait(timeout=max(0.1, grace_s))
    except subprocess.TimeoutExpired:
        process.kill()
        process.wait(timeout=max(0.1, grace_s))


class CampaignBudget:
    """A single clock; every child consumes the same remaining work allowance.

    Values are wall seconds, bytes and allocated CPU-core hours. The enclosing
    supervisor also covers hashing, Python plotting and otherwise unbounded
    child restoration. The reserved portion belongs to termination/reporting.
    """
    def __init__(self, wall_s=3600.0, reserve_s=600.0, memory_bytes=200_000_000_000,
                 artifact_bytes=50_000_000_000, core_hours=None,
                 host_free_floor_bytes=8*1024**3,disk_free_floor_bytes=100*1024**3):
        if not math.isfinite(wall_s) or wall_s <= 0 or not 0 <= reserve_s < wall_s:
            raise ValueError('wall budget must be positive and reserve smaller than total')
        if core_hours is not None and (not math.isfinite(core_hours) or core_hours <= 0):
            raise ValueError('allocated core-hour budget must be finite and positive')
        if memory_bytes<=0 or artifact_bytes<=0:
            raise ValueError('memory/artifact byte limits must be positive')
        self.started = time.monotonic()
        self.wall_s, self.reserve_s = wall_s, reserve_s
        self.memory_bytes, self.artifact_bytes = memory_bytes, artifact_bytes
        self.core_hours = core_hours
        self.host_free_floor_bytes=host_free_floor_bytes
        self.disk_free_floor_bytes=disk_free_floor_bytes

    def remaining(self):
        return max(0.0, self.wall_s-self.reserve_s-(time.monotonic()-self.started))

    def run(self, command, cwd, log, *, output, allocated_cores=0, environment=None,
            artifact_root=None):
        """Run a real command; record limits, resource history and exact argv.

        ``allocated_cores`` applies to a reserved MPI job, not the observer.
        CPU-hour usage includes setup/failed attempts for that allocation. A
        caller running several jobs must retain one allocation ledger; this
        class never grants a fresh campaign allowance to each seed.
        Use artifact_root for a family of sequential children to enforce the
        disk limit on the whole family while keeping individual resource logs.
        """
        output.mkdir(parents=True, exist_ok=True)
        log.parent.mkdir(parents=True, exist_ok=True)
        artifact_root = output if artifact_root is None else Path(artifact_root)
        # Become a Linux child subreaper before launching. A daemonized MPI
        # child is then adopted here even if its original parent exits between
        # two /proc polls. Pre-existing children are explicitly excluded.
        libc=ctypes.CDLL(None,use_errno=True)
        if libc.prctl(36,1,0,0,0)!=0:  # PR_SET_CHILD_SUBREAPER
            raise OSError(ctypes.get_errno(),"cannot supervise orphaned campaign descendants")
        baseline_children={pid for pid,record in process_table().items() if record[0]==os.getpid()}
        begin = time.monotonic()
        known = {}; reason = None; peak = 0; last_scan = -10.0; artifact_size = 0
        # Per-process user+system CPU counters are cumulative. Retain the
        # maximum for each (pid,start-ticks) identity after that process exits;
        # summing each poll would count the same work repeatedly. Processes
        # born and reaped between polls may be missed, so this is explicitly
        # a lower bound, distinct from charged allocation or scheduler CPU.
        cpu_ticks={};ticks_per_second=os.sysconf('SC_CLK_TCK')
        proc = None
        with log.open('a') as stream, (output/'budget_resources.csv').open('w') as trace:
            stream.write('COMMAND: '+json.dumps([str(x) for x in command])+'\n'); stream.flush()
            trace.write('wall_s,process_tree_rss_bytes,artifact_bytes,allocated_core_hours,sampled_cpu_core_hours_lower_bound\n')
            if self.remaining() <= 0:
                reason = 'aggregate campaign work deadline exhausted before launch'
            else:
                proc = subprocess.Popen(command, cwd=cwd, env=environment,
                    stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
            while proc is not None:
                now = time.monotonic(); table = process_table()
                adopted={pid:record[2] for pid,record in table.items()
                         if record[0]==os.getpid() and pid not in baseline_children}
                known.update(adopted)
                known = descendants(table, proc.pid, known)
                rss = 0
                for pid in known:
                    try:
                        rss += int((Path('/proc')/str(pid)/'statm').read_text().split()[1])*os.sysconf('SC_PAGE_SIZE')
                        stat=(Path('/proc')/str(pid)/'stat').read_text().rsplit(')',1)[1].split()
                        if int(stat[19])==known[pid]:
                            key=(pid,known[pid]);cpu_ticks[key]=max(cpu_ticks.get(key,0),int(stat[11])+int(stat[12]))
                    except (OSError, ValueError, IndexError):
                        pass
                peak = max(peak, rss)
                code = proc.poll()
                # A quick final child can add output between scans. Check the
                # complete shared directory even when that child already exited.
                if now-last_scan >= 5 or (code is not None and artifact_root != output):
                    artifact_size=0
                    for path in artifact_root.rglob('*'):
                        try:
                            if path.is_file() and not path.is_symlink():artifact_size+=path.stat().st_size
                        except FileNotFoundError:
                            pass # A generated build receipt may be atomically renamed.

                    last_scan = now
                # Include preflight and every earlier phase charged to this
                # allocation; a new child must never reset the core-hour clock.
                used = allocated_cores*(now-self.started)/3600
                sampled_cpu=sum(cpu_ticks.values())/ticks_per_second/3600
                trace.write(f'{now-self.started:.6f},{rss},{artifact_size},{used:.9f},{sampled_cpu:.9f}\n');trace.flush()
                if code is not None and not known:
                    # /proc enumeration is not atomic. A fast launcher can
                    # fork after the preceding table's PID list was captured,
                    # then become a zombie before its own stat was read. The
                    # table consequently contains neither a live launcher nor
                    # its newly adopted child. poll() now reaps that launcher;
                    # take a NEW identity/adoption snapshot before declaring
                    # its work complete. Keep the same baseline exclusions and
                    # birth-tick checks; never broaden ownership to unrelated
                    # children or grant the orphan a fresh execution budget.
                    table = process_table()
                    known.update({pid: record[2] for pid, record in table.items()
                        if record[0] == os.getpid() and pid not in baseline_children})
                    known = descendants(table, proc.pid, known)
                    if not known and artifact_size < self.artifact_bytes and self.remaining() > 0:
                        break
                if self.remaining() <= 0:
                    reason = 'aggregate campaign work deadline exhausted (cleanup reserve protected)'
                elif rss >= self.memory_bytes:
                    reason = 'aggregate campaign process RSS reached byte limit'
                elif artifact_size >= self.artifact_bytes:
                    reason = 'campaign artifact byte limit reached'
                elif self.core_hours is not None and used >= self.core_hours:
                    reason = 'aggregate allocated CPU-core-hour allowance exhausted'
                available=-1
                for line in Path('/proc/meminfo').read_text().splitlines():
                    if line.startswith('MemAvailable:'):available=int(line.split()[1])*1024
                if 0<=available<=self.host_free_floor_bytes:
                    reason='host available memory reached configured reserve'
                elif shutil.disk_usage(artifact_root).free<=self.disk_free_floor_bytes:
                    reason='filesystem free space reached configured reserve'
                if reason:
                    terminate_tree(proc, known, min(5.0, max(.1, self.reserve_s/2)))
                    break
                # A direct child that left live descendants is incomplete. Keep
                # supervising them instead of turning daemonized work into PASS.
                time.sleep(.1)
        reap_deadline=time.monotonic()+min(1.0,max(.1,self.reserve_s/4))
        adopted=set(known)-({proc.pid} if proc is not None else set())
        while adopted:
            for pid in tuple(adopted):
                try:
                    done,_=os.waitpid(pid,os.WNOHANG)
                    if done:adopted.remove(pid)
                except ChildProcessError:adopted.remove(pid)
            if not adopted or time.monotonic()>=reap_deadline:break
            time.sleep(.01)
        elapsed = time.monotonic()-self.started
        code = 124 if reason else proc.returncode
        result = {'command':[str(x) for x in command], 'returncode':code,
            'status':'ERROR' if reason or code not in (0,1,2,77) else ('PASS' if code==0 else 'CHILD_RESULT'),
            'reason':reason, 'wall_seconds':elapsed,'whole_wall_limit_seconds':self.wall_s,
            'cleanup_reserve_seconds':self.reserve_s,'peak_process_tree_rss_bytes':peak,
            'allocated_cores':allocated_cores,'allocated_core_hours':allocated_cores*elapsed/3600,
            'sampled_cpu_core_hours_lower_bound':sum(cpu_ticks.values())/ticks_per_second/3600,
            'actual_scheduler_cpu_core_hours':None,
            'scheduler_cpu_reason':'scheduler-neutral local executor; retain actual batch accounting separately for cluster qualification',
            'log':str(log), 'remediation':'Inspect the retained phase log. Reduce selected development work or defer long qualification to the explicit cluster tier; do not relax science gates.' if reason else None}
        (output/'budget_result.json').write_text(json.dumps(result,indent=2)+'\n')
        return result
