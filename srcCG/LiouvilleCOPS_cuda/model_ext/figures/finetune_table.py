#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/finetune_table.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/finetune_table.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 63d3a149bf7772d2c18364b6e343f159ba5907c982ada54c11470f9d309ccbe8
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: writes a COPS_model.dat-style table of a fine-tuned model, so model_obs_figures.py can plot it.
# RUN: .venv/bin/python <this> <fine-tuning dir> <fit dir it tuned> <output file>
# ==== end of header ====
"""A COPS_model.dat-style table of a fine-tuned model, for model_obs_figures.py.

    finetune_table.py <fine-tuning dir> <fit dir it tuned> <output file>

fine_tune_gpu.py (and grain_finetune.py) write pressure_comparison.txt (iPoint, refNG/RG, before/after NG/RG,
useNG/RG) but no COPS_model.dat. This takes UTC, et and altitude from the tuned fit's COPS_model.dat and the
"after" pressures (the fine-tuned model; with the grain, NG includes the tuned grain) from pressure_comparison.txt.
The products are only read; the table is written next to the figures.
"""
import sys
from pathlib import Path

ft, fit, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])


def table(path, header):
    names, rows = None, {}
    for line in path.read_text().splitlines():
        if line.startswith(header):
            names = line[1:].split()
            continue
        if line.startswith("#") or not line.strip() or names is None:
            continue
        v = line.split()
        rows[int(v[0])] = dict(zip(names, v))
    return rows


base = table(fit / "COPS_model.dat", "#iPoint")
tuned = table(ft / "pressure_comparison.txt", "# iPoint")
if set(base) != set(tuned):
    sys.exit("point sets differ: %s vs %s" % (fit, ft))
lines = ["# fine-tuned model of %s (pressures 'after' from %s/pressure_comparison.txt; geometry from %s)" % (ft.name, ft, fit),
         "#iPoint UTC et altitude pNG pRG refNG refRG useNG useRG"]
for i in sorted(base):
    b, t = base[i], tuned[i]
    lines.append(" ".join([str(i), b["UTC"], b["et"], b["altitude"], t["afterNG"], t["afterRG"], t["refNG"], t["refRG"],
                           t["useNG"], t["useRG"]]))
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text("\n".join(lines) + "\n")
