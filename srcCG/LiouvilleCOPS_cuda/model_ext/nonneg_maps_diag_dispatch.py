#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-12, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/nonneg_maps_diag_dispatch.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/nonneg_maps_diag_dispatch.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: c61c3d5dcfdbc68894fbac6c5ec0fdef12b730abb2341504614533e874ed21c2
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: diagnostic. Shows numpy's CPU dispatch of exp / log1p (AVX-512 on EPYC 9554 vs AVX2 on EPYC 7513): explains the <= 4.4e-16
#   differences of the recomputed softplus activity; the maps job sets NPY_DISABLE_CPU_FEATURES for the bitwise check.
# USED BY: job 124320 (diagnostic).
# ==== end of header ====
"""Diagnostic for nonneg_maps.py (2026-11-12): does nonneg_mult.fields reproduce the fit's activity A (nonneg_facets.npz)
bitwise on this CPU? The fits ran on gpu001/gpu002 (AMD EPYC 7513, AVX2); the CPU nodes of normal128cores are AMD EPYC
9554 (AVX-512), where numpy 1.24 dispatches other SIMD kernels of exp / log1p. Run once with the default dispatch and
once with NPY_DISABLE_CPU_FEATURES set to the AVX-512 features; prints the mismatch counts and where they are.
usage: nonneg_maps_diag_dispatch.py <fit dir> <ft dir> <base dir>
"""
import json
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nonneg_maps as M  # noqa: E402

fit, ft, base = (Path(a) for a in sys.argv[1:4])
from numpy.core._multiarray_umath import __cpu_features__ as F  # noqa: E402
print("NPY_DISABLE_CPU_FEATURES=%r; AVX512F %s AVX512_SKX %s AVX2 %s FMA3 %s" % (
    os.environ.get("NPY_DISABLE_CPU_FEATURES"), F.get("AVX512F"), F.get("AVX512_SKX"), F.get("AVX2"), F.get("FMA3")))
ftj = json.load(open(ft / "fine_tuning.json"))
nm = M.load_module(ftj["modules"]["nonneg_mult.py"]["path"], "nm_diag", ftj["modules"]["nonneg_mult.py"]["sha256"])
man = json.load(open(fit / "fit_manifest.json"))
cfg = man["model_extension"]["nonneg"]["config"]
model, s_ratio, scale = cfg["model"], cfg["s_ratio"], float(man["parameters"]["enhancement_scale"])
el = M.read_elements(base / "elements.bin")
with np.load(fit / "nonneg_facets.npz") as z:
    u, TS, A, bins, held = (np.asarray(z[k]) for k in ("u", "TS", "A", "bin", "held"))
with np.load(fit / "fit_restart.npz") as z:
    m = np.asarray(z["template_values"], float)[:7]
s = s_ratio * el["AH2O"] if model == "softplus" else None
gh = np.exp(held)
A_re = nm.fields(model, u, TS, s, m[bins.astype(np.intp)], gh, scale, False)[0]
bad = A_re != A
print("model %s: %d of %d facets differ; max rel %.3e; of them held != 0: %d; A > 0: %d" % (
    model, int(bad.sum()), len(A), float(np.max(np.abs(A_re[bad] / A[bad] - 1))) if bad.any() else 0.0,
    int(np.sum(bad & (held != 0))), int(np.sum(bad & (A > 0)))))
q0 = nm.fields(model, u, TS, s, m[bins.astype(np.intp)], None, scale, False)[0]
print("  without exp(h): differs from A at %d facets with held == 0" % int(np.sum((q0 != A) & (held == 0))))
