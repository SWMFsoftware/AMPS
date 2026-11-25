# ==== version-controlled copy: header added 2026-11-23, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/12M/full_set_q002_C_f10_allheld/checks/switched_off_1b.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/12M/full_set_q002_C_f10_allheld/checks/switched_off_1b.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: e8136efc5911733f6e7118cfb9d35f9ef2ebdda0bf7385f4720976fdf26127c8
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: decision-1b switched-off share of the all-held candidate's final state (2026-11-23)
# RUN: srun on normal128cores (checks/switched_off_1b.txt)
# ==== end of header ====
# 2026-11-23: switched-off share of the seen area of this state by user decision 1b (2026-11-21), the definition of
# figures_cluster/scripts/q002_active_area.py (T from q002_decisions_20261121/active_area.json; option C has no floor).
# Run: srun -p normal128cores -A general -N1 -n1 -c 2 --mem=24G -t 00:10:00 .venv/bin/python -I switched_off_1b.py
import json, numpy as np
from pathlib import Path
R = Path("/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5")
st = R / "nonneg/localnn/slopetest_c_q002_f10_pm0.01_tploff_20261122T195607Z/chain/c5/tools_state_ft"
T = json.load(open("/data001/highresolu_vt/cops_output/figures_cluster/12M/q002_decisions_20261121/active_area.json"))["T"]
with np.load(st / "finetune_state.npz") as z:
    A = z["A_final"].astype(float)
with np.load(st / "fit_state/nonneg_facets.npz") as z:
    seen = z["seen"].astype(bool)
n = len(A)
d = np.memmap(R / "base_q002/elements.bin", dtype="<f8", mode="r", offset=12, shape=(n, 34))
area = np.array(d[:, 6]); orig = np.array(d[:, 32])
print("T", T, "T recomputed", 0.01 * float(np.median(orig[seen])))
aseen = area[seen].sum()
for f, tag in ((1, "T"), (0.1, "T x0.1"), (10, "T x10")):
    print("off %-7s %.4f" % (tag, area[seen & (A <= f * T)].sum() / aseen))
print("A<=0 seen area", area[seen & (A <= 0)].sum() / aseen, "A_min", A.min())
