"""Workspaces: a private, configured copy of the night's commit for each runner.

Runners reconfigure and rebuild AMPS in the tree they run in and do not restore it, so every runner gets
its own git worktree of the same commit, made from local mirrors (updated once per night):
    {workdir}/mirrors/<name>.git          bare mirrors of AMPS and the companion repositories
    {night}/ws/<config>/<runner>/         a worktree of AMPS at the night's commit, with companions
                                          (worktrees of their mirrors, or links to machine data)
"""

from __future__ import annotations

import os
import re
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from . import procs
from .machines import Config, Machine
from .runners import Checkout


class WorkspaceError(RuntimeError):
    pass


def _git(*args: str, cwd: Optional[Path] = None, timeout: int = 1800) -> str:
    done = subprocess.run(["git"] + list(args), cwd=str(cwd) if cwd else None, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, universal_newlines=True, timeout=timeout)
    if done.returncode != 0:
        raise WorkspaceError("git {} failed: {}".format(" ".join(args), (done.stderr.strip().splitlines() or ["?"])[-1]))
    return done.stdout.strip()


def _name(repository: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]", "_", repository.rstrip("/").split("/")[-1].split(":")[-1]).replace(".git", "")


class Mirrors:
    """Bare mirrors of the tested repository and its companions."""

    def __init__(self, root: Path, checkout: Checkout):
        self.root = Path(root)
        self.checkout = checkout
        self.commits: Dict[str, str] = {}     # path ("" = AMPS itself) -> commit

    def path(self, repository: str) -> Path:
        return self.root / (_name(repository) + ".git")

    def update(self, log) -> None:
        """Clone missing mirrors, fetch the others, and pin each repository's branch head for tonight."""
        self.root.mkdir(parents=True, exist_ok=True)
        repos = [("", self.checkout.repository, self.checkout.branch)] + \
                [(c.path, c.repository, c.branch) for c in self.checkout.companions if c.repository]
        for path, repository, branch in repos:
            mirror = self.path(repository)
            if mirror.is_dir():
                log("fetch " + repository)
                _git("fetch", "--prune", "--quiet", "origin", cwd=mirror)
            else:
                log("clone --mirror " + repository)
                _git("clone", "--mirror", "--quiet", repository, str(mirror))
            self.commits[path] = _git("rev-parse", "refs/heads/" + branch, cwd=mirror)
            log("  {} {} = {}".format(path or "AMPS", branch, self.commits[path]))

    @property
    def commit(self) -> str:
        return self.commits[""]


def create(target: Path, mirrors: Mirrors, machine: Machine, values: Dict[str, str], log) -> None:
    """A worktree of the night's commit at target, with the companions in place."""
    if target.exists():
        remove(target, mirrors)
    target.parent.mkdir(parents=True, exist_ok=True)
    _git("worktree", "add", "--detach", "--quiet", str(target), mirrors.commit, cwd=mirrors.path(mirrors.checkout.repository))
    for companion in mirrors.checkout.companions:
        where = target / companion.path
        if companion.repository:
            if where.exists():
                shutil.rmtree(str(where))
            _git("worktree", "add", "--detach", "--quiet", str(where), mirrors.commits[companion.path],
                 cwd=mirrors.path(companion.repository))
        else:
            source = values.get(companion.link_var, "")
            if source and Path(source).expanduser().exists():
                if where.is_symlink() or where.exists():
                    where.unlink() if where.is_symlink() else shutil.rmtree(str(where))
                where.symlink_to(Path(source).expanduser())
            else:
                log("note: {} not linked: machine variable {} is not set or does not exist".format(
                    companion.path, companion.link_var))


def configure(target: Path, machine: Machine, config: Config, env: Dict[str, str], log_file: Path, tag: str,
              timeout_s: float) -> Tuple[bool, str]:
    """The configuration's Config.pl calls and after_install lines; (ok, failing command)."""
    env = dict(env, AMPS_TABLE_HOST=machine.table_host)
    steps = [["./Config.pl"] + shlex.split(call) for call in config.config_pl] + \
            [procs.shell(line) for line in config.after_install]
    for argv in steps:
        outcome = procs.run(with_setup(argv, config.setup), target, env, log_file, timeout_s, tag, header="configure")
        if outcome.exit_code != 0:
            return False, " ".join(argv)
    return True, ""


def with_setup(argv: List[str], setup: List[str]) -> List[str]:
    """Run argv after the configuration's setup lines (module loads) in one bash."""
    if not setup:
        return argv
    return procs.shell("; ".join(setup) + " && exec " + " ".join(shlex.quote(a) for a in argv))


def remove(target: Path, mirrors: Optional[Mirrors] = None) -> None:
    """Delete a workspace and forget its worktrees."""
    if target.exists():
        shutil.rmtree(str(target), ignore_errors=True)
    if mirrors:
        for repository in [mirrors.checkout.repository] + [c.repository for c in mirrors.checkout.companions if c.repository]:
            mirror = mirrors.path(repository)
            if mirror.is_dir():
                subprocess.run(["git", "worktree", "prune"], cwd=str(mirror), stdout=subprocess.DEVNULL,
                               stderr=subprocess.DEVNULL)


def size_mb(path: Path) -> float:
    total = 0
    for root, _dirs, files in os.walk(str(path)):
        for name in files:
            try:
                total += os.lstat(os.path.join(root, name)).st_size
            except OSError:
                pass
    return total / 2.0 ** 20
