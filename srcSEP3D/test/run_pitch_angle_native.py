#!/usr/bin/env python3
"""srcSEP3D native smoke test of the pitch-angle library binding (PAN01-PAN05).

Runs the built ./amps on the reduced-front injection deck
(examples/application-input/amps-injection-smoke.in) switched to the
focused-diffusion mover with the shared pitch-angle library:

    transport = focused-diffusion
    pitch_angle_diffusion_model = pitch-angle-library
    #begin physics: pitch_angle_diffusion   (model and rate-fraction nu_max*dt
                                             of the first Milstein row of the
                                             library certification matrix,
                                             pitch_angle_certified_pairs.inc)

for --steps time steps on one and on four MPI ranks (each in its own working
directory), and checks

  PAN01.rankR  the run completes and the startup log names the installed model
               (pitch_angle_diffusion_model=<certified model>) and its 64-hex
               fingerprint;
  PAN02.rankR  the global particle ledger of the last publication
               (shocks.csv: per-step injected/escaped/absorbed/failed) shows
               injected particles and zero mover failures (a library or step
               rule failure inside the mover would be counted as 'failed');
  PAN03        the per-step ledger rows are identical on one and four ranks
               (keyed per-particle random streams make the transport rank
               independent);
  PAN04        a PA-QFORM deck stops at start-up as an incompatible
               coefficient-kernel pair (unbounded drift at mu = 0);
  PAN05        the certified deck with reflecting-euler-maruyama stops at
               start-up (Euler-Maruyama is not certified);
  PAN06        a PA-EPS (eps = 0.048) Milstein deck stops at start-up as an
               uncertified pair (only while the matrix has no PA-EPS/Milstein
               row; the case reports SKIP otherwise).

This is a native integration smoke test of the binding, not a physics
validation of focused transport (the library's acceptance campaign and the
srcSEP3D component tests PAD3D01-03 carry the physics evidence).

Usage (from the AMPS root, after the native srcSEP3D build):
  python3 srcSEP3D/test/run_pitch_angle_native.py --amps ./amps \\
      --output-dir test_output/pitch-angle-native
Output lines '[ID] STATUS message' (read_bracket_text); exit 1 on any failure.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SMOKE = ROOT / "srcSEP3D" / "examples" / "application-input"
DECK = SMOKE / "amps-injection-smoke.in"
LIBRARY = ROOT / "src" / "models" / "sep_common" / "pitch_angle_diffusion"


MODEL_NAMES = {"Isotropic": "PA-ISO", "Eprem": "PA-EPREM", "EpsilonForm": "PA-EPS"}


def certified_rows() -> list[list[str]]:
    """Milstein rows of the library matrix, as field lists
    {ModelId::X, Kernel::Milstein, nu_max*dt, eps_min, eps_max, "evidence"}."""
    rows = []
    text = (LIBRARY / "pitch_angle_certified_pairs.inc").read_text()
    for line in text.splitlines():
        if line.strip().startswith("//") or "Kernel::Milstein" not in line:
            continue
        rows.append([field.strip() for field in line.strip(" {},").split(",")])
    return rows


def certified_deck_lines(rows: list[list[str]]) -> tuple[str, str, str]:
    """(model name, model lines, nu_max*dt) of the first Milstein row; a PA-EPS
    row uses the middle of its certified epsilon range."""
    if not rows:
        raise SystemExit("no Milstein row in pitch_angle_certified_pairs.inc "
                         "(src/models/sep_common/pitch_angle_diffusion/"
                         "IMPLEMENTATION_STATUS.md)")
    fields = rows[0]
    name = MODEL_NAMES[fields[0].split("::")[-1]]
    lines = f"model = {name}\n"
    if name == "PA-EPS":
        lines += f"epsilon = {0.5 * (float(fields[3]) + float(fields[4])):.17g}\n"
    return name, lines, fields[2]


def block(model_lines: str, fraction: str) -> str:
    return ("\n#begin physics: pitch_angle_diffusion\n" + model_lines +
            "lambda_reference_m = 4.48793612e10\n"
            "lambda_reference_rigidity_v = 1e9\n"
            "lambda_rigidity_exponent = 1/3\n"
            "lambda_reference_radius_m = 1.495978707e11\n"
            "lambda_radial_exponent = 2/3\n"
            "step_rule = rate-fraction\n"
            f"rate_fraction = {fraction}\n"
            "#end physics\n")


def write_deck(work: Path, name: str, steps: int, model_lines: str, scheme: str,
               fraction: str) -> Path:
    """Copy the smoke deck: sep3d.in edited (mover, selector, scheme) and
    referenced by absolute path; the library block appended."""
    sep3d = (SMOKE / "sep3d.in").read_text()
    for old, new in [("transport = parker", "transport = focused-diffusion"),
                     ("pitch_angle_diffusion_model = jokipii-1966",
                      "pitch_angle_diffusion_model = pitch-angle-library"),
                     ("pitch_angle_scheme = reflecting-milstein",
                      f"pitch_angle_scheme = {scheme}")]:
        if old not in sep3d:
            raise SystemExit(f"sep3d.in anchor not found: {old}")
        sep3d = sep3d.replace(old, new, 1)
    (work / f"{name}-sep3d.in").write_text(sep3d)
    text = DECK.read_text()
    text, budget = re.subn(r"maximum_time_steps = \d+", f"maximum_time_steps = {steps}",
                           text, count=1)
    if budget != 1:
        raise SystemExit(f"{DECK}: step budget anchor not found")
    text = text.replace('#include "sep3d.in"',
                        f'#include "{(work / (name + "-sep3d.in")).as_posix()}"')
    text = re.sub(r'#include "([^"/]+)"',
                  lambda m: f'#include "{(SMOKE / m.group(1)).as_posix()}"', text)
    deck = work / f"{name}.in"
    deck.write_text(text + block(model_lines, fraction))
    return deck


def emit(test_id: str, status: str, message: str, results: list) -> None:
    results.append(status)
    print(f"[{test_id}] {status} {message}", flush=True)


def run(args, deck: Path, rundir: Path, ranks: int, extra: list[str]) -> tuple[int, str]:
    shutil.rmtree(rundir, ignore_errors=True)
    rundir.mkdir(parents=True)
    log = rundir / "run.log"
    command = [args.mpiexec, "-n", str(ranks), str(Path(args.amps).resolve()),
               "-input", str(deck), *extra]
    try:
        with log.open("w") as stream:
            code = subprocess.run(command, cwd=rundir, stdout=stream,
                                  stderr=subprocess.STDOUT, timeout=args.timeout).returncode
    except subprocess.TimeoutExpired:
        return -1, f"timeout; log {log}"
    return code, log.read_text(errors="replace")


def ledger(rundir: Path) -> list[tuple[int, ...]] | None:
    publications = sorted((rundir / "publications").glob("sep3d.*"))
    if not publications:
        return None
    rows = []
    with (publications[-1] / "shocks.csv").open() as stream:
        for row in csv.DictReader(stream):
            rows.append(tuple(int(row[key]) for key in (
                "step", "species", "injected_count", "escaped_count",
                "absorbed_count", "failed_count")))
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--amps", default=str(ROOT / "amps"))
    parser.add_argument("--mpiexec", default="mpiexec")
    parser.add_argument("--output-dir", default=str(ROOT / "test_output/pitch-angle-native"))
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--timeout", type=int, default=3600)
    args = parser.parse_args()
    work = Path(args.output_dir).resolve()
    work.mkdir(parents=True, exist_ok=True)
    results: list = []
    certified = certified_rows()
    model_name, model_lines, fraction = certified_deck_lines(certified)
    valid = write_deck(work, "valid", args.steps, model_lines, "reflecting-milstein",
                       fraction)
    ledgers = {}
    for ranks in (1, 4):
        rundir = work / f"rank{ranks}"
        code, text = run(args, valid, rundir, ranks, ["--output-dir", str(rundir / "publications")])
        fingerprint = re.search(r"pitch_angle_diffusion_configuration_fingerprint=([0-9a-f]{64})",
                                text)
        ok = (code == 0 and f"pitch_angle_diffusion_model={model_name}" in text
              and fingerprint and "End of the run" in text)
        emit(f"PAN01.rank{ranks}", "PASS" if ok else "FAIL",
             (f"{model_name} installed (fingerprint {fingerprint.group(1)[:16]}...), {args.steps} steps "
              f"completed" if ok else f"exit {code}; log {rundir / 'run.log'}"), results)
        rows = ledger(rundir) if ok else None
        if rows is None:
            emit(f"PAN02.rank{ranks}", "ERROR", "no ledger (run unavailable)", results)
            continue
        injected = sum(r[2] for r in rows)
        failed = sum(r[5] for r in rows)
        escaped = sum(r[3] for r in rows)
        absorbed = sum(r[4] for r in rows)
        good = injected > 0 and failed == 0
        emit(f"PAN02.rank{ranks}", "PASS" if good else "FAIL",
             f"ledger over {len(rows)} steps: injected {injected}, escaped {escaped}, "
             f"absorbed {absorbed}, failed {failed}", results)
        ledgers[ranks] = rows
    if 1 in ledgers and 4 in ledgers:
        same = ledgers[1] == ledgers[4]
        emit("PAN03", "PASS" if same else "FAIL",
             "per-step ledger rows identical on 1 and 4 ranks" if same
             else "ledger rows differ between 1 and 4 ranks", results)
    else:
        emit("PAN03", "ERROR", "a run is unavailable", results)
    refusals = [
        ("PAN04", "qform", "model = PA-QFORM\nq = 5/3\ngap_h = 0.05\n", "reflecting-milstein",
         "unbounded"),
        ("PAN05", "euler", model_lines, "reflecting-euler-maruyama", "euler-maruyama"),
    ]
    if any(MODEL_NAMES[r[0].split("::")[-1]] == "PA-EPS" for r in certified):
        emit("PAN06", "SKIP", "the matrix certifies a PA-EPS/Milstein pair", results)
    else:
        refusals.append(("PAN06", "eps", "model = PA-EPS\nepsilon = 0.048\n",
                         "reflecting-milstein", "no certification entry"))
    for test_id, name, lines, scheme, expected in refusals:
        deck = write_deck(work, name, args.steps, lines, scheme, fraction)
        rundir = work / name
        code, text = run(args, deck, rundir, 1,
                         ["--initialization-only", "--initialization-output-dir", str(rundir / "init")])
        ok = code != 0 and expected in text
        emit(test_id, "PASS" if ok else "FAIL",
             f"refused at start-up (exit {code}) with '{expected}'" if ok
             else f"exit {code}, expected refusal containing '{expected}'; log {rundir / 'run.log'}",
             results)
    bad = sum(status in ("FAIL", "ERROR") for status in results)
    skipped = sum(status == "SKIP" for status in results)
    print(f"pitch-angle native summary: {len(results) - bad - skipped} passed, "
          f"{bad} failed/errored, {skipped} skipped")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
