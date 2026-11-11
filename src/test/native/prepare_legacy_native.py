#!/usr/bin/env python3
"""Preflight a legacy native build before configuration rewrites production hooks.

Legacy make targets may be invoked with unlimited ``-j``. Their enclosing make
processes are expected ancestors, but another build, AMPS executable or test
runner in the same checkout must not overlap removal/configuration. This hook
does not select an application: the immediately following Config.pl recipe
remains the established source of the generated configuration.
"""
from pathlib import Path
import os
import subprocess
import sys


def main():
    root = Path(__file__).resolve().parents[3]
    if Path.cwd().resolve() != root or not (root / "Config.pl").is_file():
        raise RuntimeError("native preflight requires the AMPS repository root")
    ancestors = {os.getpid()}
    pid = os.getppid()
    while pid > 1 and pid not in ancestors:
        ancestors.add(pid)
        try:
            pid = int(Path(f"/proc/{pid}/stat").read_text().rsplit(")", 1)[1].split()[1])
        except (OSError, ValueError, IndexError):
            break
    active = []
    for process in Path("/proc").iterdir():
        if not process.name.isdigit() or int(process.name) in ancestors:
            continue
        try:
            cwd = (process / "cwd").resolve()
            cwd.relative_to(root)
            argv = (process / "cmdline").read_bytes().split(b"\0")
            name = Path(os.fsdecode(argv[0])).name
            arguments = b" ".join(argv).decode(errors="replace")
        except (OSError, ValueError):
            continue
        if name in {"make", "gmake", "amps", "Config.pl", "ampsConfig.pl"} or any(
                token in arguments for token in ["run_tests.py", "run_all_stage_tests.sh"]):
            active.append(f"pid={process.name} {arguments}")
    if active:
        raise RuntimeError("native preflight found active work: " + "; ".join(active))
    build = root / "build"
    if build.is_symlink():
        raise RuntimeError("refusing a symlink at the generated build path")
    print(f"[preflight] root={root} active_native_processes=0 removal=build/", flush=True)
    # The repository's required native sequence names exactly this generated
    # directory. Do not use a configurable recursive target or touch site data.
    subprocess.run(["rm", "-rf", "--", "build"], cwd=root, check=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"native preflight: {error}", file=sys.stderr)
        raise SystemExit(2)
