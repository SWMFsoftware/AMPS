#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/restart_grain_start.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/restart_grain_start.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 8bcd4bc84a4ff93fbff013aaed43afb537aab4514b34588dc486985283a7c6da
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: copies a grain fit's restart file and sets its grain start (A, B) to a fine tuning's tuned values.
# USED FOR: every 12M fit/fine-tuning cycle and the late-weights fit (start = previous fit + fine-tuned grain).
# RUN: .venv/bin/python <this> --restart <fit dir>/fit_restart.npz --finetune <its grain fine-tuning dir> --out <new npz>
# OUT: a new restart file (never overwrites) with a history entry; checks that the fine tuning was made from it.
# ==== end of header ====
"""A copy of a grain fit's restart file whose grain start is the fine-tuned grain (user decision 2026-11-05, for the
12M fit/fine-tuning cycles: "restart after the fine tuning ... the fine-tuned values" for the grain).

    .venv/bin/python <this> --restart <fit dir>/fit_restart.npz --finetune <grain fine-tuning dir> --out <new file>

* Only x[grain_A_uPa] and x[grain_B_uPa_per_100s] change, to A_Pa / 1e-6 and B_Pa_per_s / 1e-8 of
  <finetune>/grain_finetune.json (grain_crossing.py units). gpu_fit takes every start value by name from the
  restart's (names, x), so the next fit starts there; the stored predictions, activity and the other parameters are
  the fit's, unchanged. A history entry records the change; the file is new (restart.write never replaces).
* The fine tuning must have been made from this restart (checked: its grain_term.json source is this fit dir).
"""

import argparse
import datetime as dt
import json
import sys
from pathlib import Path

APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
sys.path[:0] = [str(APP), str(APP / "tools")]

A_UNIT, B_UNIT = 1e-6, 1e-6 / 100.0      # grain_crossing.py


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--restart", required=True, type=Path)
    ap.add_argument("--finetune", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()
    from cops_cuda import provenance, restart as rs
    g = json.loads((a.finetune / "grain_finetune.json").read_text())
    if Path(g["source_grain_term"]).parent.resolve() != a.restart.parent.resolve():
        sys.exit("restart_grain_start: %s was not tuned from %s" % (a.finetune, a.restart))
    r = rs.read(a.restart).copy()
    names = list(r.meta["names"])
    iA, iB = names.index("grain_A_uPa"), names.index("grain_B_uPa_per_100s")
    old = (float(r.x[iA]), float(r.x[iB]))
    r.arrays["x"][iA] = g["A_Pa"] / A_UNIT
    r.arrays["x"][iB] = g["B_Pa_per_s"] / B_UNIT
    r.meta.setdefault("history", []).append({
        "utc": dt.datetime.now(dt.timezone.utc).isoformat(), "action": "grain start from the fine tuning",
        "finetune": str(a.finetune), "grain_finetune_json": provenance.file_identity(a.finetune / "grain_finetune.json"),
        "source_restart": provenance.file_identity(a.restart),
        "grain_A_uPa": [old[0], float(r.x[iA])], "grain_B_uPa_per_100s": [old[1], float(r.x[iB])]})
    ident = rs.write(a.out, r)
    print(json.dumps({"out": str(a.out), "sha256": ident["sha256"], "grain_A_uPa": [old[0], float(r.x[iA])],
                      "grain_B_uPa_per_100s": [old[1], float(r.x[iB])]}, indent=1), flush=True)


if __name__ == "__main__":
    main()
