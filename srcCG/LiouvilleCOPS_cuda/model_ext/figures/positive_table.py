#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-07, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/positive_table.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/positive_table.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: bc7b5670b859c8c24e55964e590ab5ee03d8abc8700f8cf5524cc0b91f88e3a5
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: figure helper. Observation rows of a positive fit's COPS_model.dat (the 401 model-only grid rows removed,
#   checked against late_window_model.dat), for model_obs_figures.py. 2026-11-07.
# USED BY: figures/figures_12m_positive_scan.sh.
# RUN: .venv/bin/python positive_table.py <fit segment dir> <output file>
# ==== end of header ====
"""Observation rows of a positive-activity fit's COPS_model.dat, for model_obs_figures.py (user 2026-11-07).

    positive_table.py <fit segment dir> <output file>

The positive fits (model_ext/positive_fit.py) add 401 model-only rows of the 1-s late grid to COPS_model.dat
(kind "grid" in late_window_model.dat; reference 0, use 0; pNG = instantaneous surface + grain): they are not
observations and must not be drawn (their header line PositiveFitGridRows says so). This copies COPS_model.dat with
those rows removed; every other line, header included, is kept verbatim. Checks: the removed rows are exactly the
grid rows of late_window_model.dat and all have reference 0 and use 0. Read only; the table goes next to the figures.
"""
import sys
from pathlib import Path

seg, out = Path(sys.argv[1]), Path(sys.argv[2])
grid = set()
for line in (seg / "late_window_model.dat").read_text().splitlines():
    if not line.startswith("#"):
        v = line.split()
        if v[3] == "grid":
            grid.add(int(v[0]))
names, keep, dropped = None, [], 0
for line in (seg / "COPS_model.dat").read_text().splitlines():
    if line.startswith("#iPoint"):
        names = line[1:].split()
    if line.startswith("#") or not line.strip():
        keep.append(line)
        continue
    v = line.split()
    if int(v[0]) in grid:
        c = {n: i for i, n in enumerate(names)}
        if float(v[c["refNG"]]) > 0 or float(v[c["refRG"]]) > 0 or v[c["useNG"]] != "0" or v[c["useRG"]] != "0":
            sys.exit("grid row %s has data or a use flag: %s" % (v[0], seg))
        dropped += 1
        continue
    keep.append(line)
if dropped != len(grid):
    sys.exit("removed %d rows, late_window_model.dat lists %d grid rows: %s" % (dropped, len(grid), seg))
out.parent.mkdir(parents=True, exist_ok=True)
keep.insert(1, "# positive_table.py: the %d grid rows (model-only 1-s points) of %s/COPS_model.dat removed; observation rows only"
            % (dropped, seg))
out.write_text("\n".join(keep) + "\n")
print("%s: %d observation rows (%d grid rows removed)" % (out, sum(1 for l in keep if l and not l.startswith("#")), dropped))
