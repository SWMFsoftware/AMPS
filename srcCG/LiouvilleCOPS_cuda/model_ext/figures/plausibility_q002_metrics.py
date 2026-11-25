#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/plausibility_q002_metrics.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/plausibility_q002_metrics.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: f4df6102ea8082dacd68bc7663dc38be27603c5a49d567e86d35e01b03de7c65
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: activity-spread measures (zero area, graded active area, effective area, f50/f90, Gini) and the energy-limit
#   test (share of Q_H2O with F_H2O > S0 / (r_h^2 dH)) of the ten new-law fits, for the plausibility ranking (2026-11-21)
# RUN: CPU job jobs/plausibility_q002_metrics.sbatch; output figures_cluster/12M/q002_plausibility/
# ==== end of header ====
"""How much of the surface each new-law (q_min = 0.02 q_max) fit keeps active (2026-11-21).

For the ten runs of the q002 round (D2_q002_warm and the nine regularized stage-2 chains) this reads the
final H2O activity per facet, A = finetune_state.npz["A_final"] (the field every map of the round is drawn
from; pub_maps_figures.load_fields checks it against the maps' *_static.vtu), the seen mask of the fit
(nonneg_facets.npz["seen"]) and the facet areas of the base (elements.bin), and writes area-weighted
measures of how the activity is spread over the surface, for the whole surface and for the seen part:

  zero_area        area fraction with A == 0 (exactly; facts.json "A_zero" gives the same number)
  ge_1pct, ge_10pct, ge_50pct
                   area fraction with A >= 0.01, 0.1, 0.5 x the area-weighted mean activity of that region
  effective_area   participation ratio (sum a A)^2 / (sum a * sum a A^2): the area fraction a uniform
                   activity would need to give the same first and second moments (1 = uniform)
  f50, f90         smallest area fraction that carries 50 % / 90 % of sum a A (uniform: 0.5 / 0.9)
  top1pct_share    share of sum a A carried by the most active 1 % of the area
  gini             area-weighted Gini coefficient of A (0 uniform, 1 all on one facet)
  ge_1pct_trim, ge_10pct_trim
                   as ge_1pct / ge_10pct, but relative to the trimmed mean (the mean without the most active
                   1 % of the area), so that a few extreme facets do not set the reference level

Energy limit (added 2026-11-21): no facet can sublimate more water than the sunlight it absorbs pays for.
With the solar constant S0 = 1361 W m^-2 (TSI), zero albedo, normal incidence, no thermal emission and no
heat conduction, the flux can never exceed F_E = S0 / (r_h[AU]^2 * dH), dH = 51.0 kJ/mol / N_A (latent heat
of sublimation of water ice, about 2.83 MJ/kg). This is an upper bound, not an estimate (real surfaces
reradiate most of the absorbed energy at 3.8 AU). From each run's H2O flux map at 10:33:10 (the *_flux_*.vtu
the figures use, read with the maps' frozen nonneg_maps.read_vtu) the tool gives the maximum F_H2O, the
area and the share of Q_H2O with F_H2O > F_E, and (from maps_manifest.json "times" statistics) the largest
F_H2O / F_E over all map times.

Inputs come from each run's maps_manifest.json (base, fine-tuning and fit folders), so the numbers belong to
exactly the state the figures show. Nothing is modified.

  python3 plausibility_q002_metrics.py --out-dir <dir>      (writes metrics.json, metrics.csv; never overwrites)
"""
import argparse
import csv
import hashlib
import importlib.util
import json
import os
import struct
import time
from pathlib import Path

import numpy as np

S0_W_M2 = 1361.0                          # total solar irradiance at 1 AU
DH_J_MOL = 51.0e3                         # latent heat of sublimation of water ice
N_A = 6.02214076e23
MAP_TIME = "2016-09-30T10:33:10"

ROOT = Path("/data001/highresolu_vt/cops_output")
RUNS_JSON = ROOT / "figures_cluster/scripts/data/source_rate_q002_runs_20261118.json"
ORDER = ["D2_q002_warm", "softplus_q002_reg", "A_q002_reg", "B_q002_reg", "D_q002_reg", "C_f1_pm1em4_q002_reg",
         "C_f10_pm1em4_q002_reg", "C_f1_pm0p01_q002_reg", "C_f10_pm0p01_q002_reg", "C_f10_pm0p01_tploff_q002_reg"]
T0 = time.time()


def say(t):
    print("plausibility_q002_metrics [%5.0f s]: %s" % (time.time() - T0, t), flush=True)


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


def read_area(path):
    """Facet areas of elements.bin version 1 (column 6 of 34 doubles per facet), read column-wise."""
    size = os.path.getsize(path)
    with open(path, "rb") as f:
        version, n, nb = struct.unpack("<iii", f.read(12))
    if version != 1 or nb != 25 or size != 12 + n * 34 * 8:
        raise SystemExit("unsupported element file %s" % path)
    mm = np.memmap(path, dtype="<f8", mode="r", offset=12, shape=(n, 34))
    area = np.array(mm[:, 6], dtype=np.float64)
    del mm
    return area


def measures(A, a):
    """Area-weighted spread measures of activity A over facets of area a."""
    tot_a = a.sum()
    w = a * A
    tot_w = w.sum()
    mean = tot_w / tot_a
    out = {"area_m2": float(tot_a), "mean_activity": float(mean), "zero_area": float(a[A == 0].sum() / tot_a)}
    for name, t in (("ge_1pct", 0.01), ("ge_10pct", 0.1), ("ge_50pct", 0.5)):
        out[name] = float(a[A >= t * mean].sum() / tot_a)
    out["effective_area"] = float(tot_w ** 2 / (tot_a * np.sum(a * A * A)))
    order = np.argsort(-A, kind="stable")
    ca = np.cumsum(a[order]) / tot_a                 # area fraction of the most active facets
    cw = np.cumsum(w[order]) / tot_w                 # their share of sum a A
    out["f50"] = float(ca[np.searchsorted(cw, 0.5)])
    out["f90"] = float(ca[np.searchsorted(cw, 0.9)])
    out["top1pct_share"] = float(cw[np.searchsorted(ca, 0.01)])
    # Gini of A with area weights: 1 - 2 * area under the Lorenz curve (ascending order)
    asc = order[::-1]
    la = np.concatenate(([0.0], np.cumsum(a[asc]) / tot_a))
    lw = np.concatenate(([0.0], np.cumsum(w[asc]) / tot_w))
    out["gini"] = float(1.0 - np.sum((la[1:] - la[:-1]) * (lw[1:] + lw[:-1])))
    out["max_over_mean"] = float(A.max() / mean)
    # reference without the most active 1 % of the area
    keep = order[np.searchsorted(ca, 0.01) + 1:]
    trim = float(w[keep].sum() / a[keep].sum())
    out["trimmed_mean_activity"] = trim
    for name, t in (("ge_1pct_trim", 0.01), ("ge_10pct_trim", 0.1)):
        out[name] = float(a[A >= t * trim].sum() / tot_a)
    return out


def frozen_maps_module(code_sha):
    for d in sorted((ROOT / "model_ext_frozen").glob("nonneg_maps_*")):
        p = d / "nonneg_maps.py"
        if p.exists() and sha256(p) == code_sha:
            spec = importlib.util.spec_from_file_location("nonneg_maps_%s" % code_sha[:12], p)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            return mod, str(p)
    raise SystemExit("no model_ext_frozen/nonneg_maps_* with sha256 %s" % code_sha)


def energy_numbers(maps, man, area):
    """F_H2O against the energy limit at MAP_TIME (flux map) and at every map time (manifest statistics)."""
    tm = [t for t in man["times"] if t["utc"] == MAP_TIME]
    if len(tm) != 1:
        raise SystemExit("%s: no map at %s" % (maps, MAP_TIME))
    tm = tm[0]
    nm, code = frozen_maps_module(man["code_sha256"]["nonneg_maps.py"])
    rd, fd, _ = nm.read_vtu(maps / tm["vtu"])
    F = np.asarray(rd["F_H2O_per_m2_s"], np.float64)
    if len(F) != len(area):
        raise SystemExit("%s: flux map has %d facets, base %d" % (maps, len(F), len(area)))
    Q = float(np.sum(F * area))
    if abs(Q - tm["Q_H2O_s-1"]) > 1e-9 * Q:
        raise SystemExit("%s: Q from the flux map %.6e differs from the manifest %.6e" % (maps, Q, tm["Q_H2O_s-1"]))
    rh = float(tm["rh_AU"])
    FE = S0_W_M2 / rh ** 2 / (DH_J_MOL / N_A)
    over = F > FE
    over10 = F > 0.1 * FE
    worst = max(t["stats"]["F_H2O"]["max"] * (t["rh_AU"] ** 2) for t in man["times"]) / (S0_W_M2 / (DH_J_MOL / N_A))
    return {"utc": MAP_TIME, "rh_AU": rh, "F_E_m-2_s-1": FE, "Q_H2O_s-1": Q, "F_H2O_max": float(F.max()),
            "F_max_over_F_E": float(F.max() / FE), "area_fraction_over_F_E": float(area[over].sum() / area.sum()),
            "Q_share_over_F_E": float(np.sum(F[over] * area[over]) / Q),
            "Q_share_over_0p1_F_E": float(np.sum(F[over10] * area[over10]) / Q),
            "facets_over_F_E": int(over.sum()), "max_F_over_F_E_all_map_times": float(worst),
            "map_times": len(man["times"]), "flux_vtu": str(maps / tm["vtu"]), "nonneg_maps_code": code}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    out = Path(a.out_dir)
    finals = [out / "metrics.json", out / "metrics.csv"]
    for p in finals:
        if p.exists():
            raise SystemExit("%s exists (never overwritten; move it aside first)" % p)
    out.mkdir(parents=True, exist_ok=True)
    runs = json.load(open(RUNS_JSON))["runs"]
    area_cache, res = {}, {"schema": "plausibility_q002_metrics_v1", "runs_json": str(RUNS_JSON),
                           "runs_json_sha256": sha256(RUNS_JSON), "runs": {}}
    for key in ORDER:
        maps = Path(runs[key]["maps"])
        man = json.load(open(maps / "maps_manifest.json"))
        base, ft, fit = Path(man["base"]), Path(man["state"]["ft_dir_resolved"]), Path(man["state"]["fit_dir_resolved"])
        if str(base) not in area_cache:
            say("areas of %s" % base)
            area_cache[str(base)] = read_area(base / "elements.bin")
        area = area_cache[str(base)]
        with np.load(ft / "finetune_state.npz", allow_pickle=False) as z:
            A = np.asarray(z["A_final"], np.float64)
        with np.load(fit / "nonneg_facets.npz", allow_pickle=False) as z:
            seen = np.asarray(z["seen"]).astype(bool)
        if not (len(A) == len(area) == len(seen)):
            raise SystemExit("%s: facet counts differ (A %d, area %d, seen %d)" % (key, len(A), len(area), len(seen)))
        if np.any(A < 0) or not np.all(np.isfinite(A)):
            raise SystemExit("%s: A has negative or non-finite values" % key)
        r = {"maps": str(maps), "base": str(base), "ft_dir": str(ft), "fit_dir": str(fit),
             "finetune_state_sha256": sha256(ft / "finetune_state.npz"),
             "whole": measures(A, area), "seen": measures(A[seen], area[seen]),
             "seen_area_fraction": float(area[seen].sum() / area.sum()),
             "energy": energy_numbers(maps, man, area)}
        res["runs"][key] = r
        w, s, e = r["whole"], r["seen"], r["energy"]
        say("%-30s F_max %.3e = %.3g F_E (F_E %.3e at %.4f AU); Q share above F_E %.4f, above 0.1 F_E %.4f; "
            "area above F_E %.2e; worst over map times %.3g F_E; trimmed >=1%% %.3f/%.3f >=10%% %.3f/%.3f"
            % (key, e["F_H2O_max"], e["F_max_over_F_E"], e["F_E_m-2_s-1"], e["rh_AU"], e["Q_share_over_F_E"],
               e["Q_share_over_0p1_F_E"], e["area_fraction_over_F_E"], e["max_F_over_F_E_all_map_times"],
               w["ge_1pct_trim"], s["ge_1pct_trim"], w["ge_10pct_trim"], s["ge_10pct_trim"]))
        say("%-30s zero %.3f/%.3f  >=1%% %.3f/%.3f  >=10%% %.3f/%.3f  eff %.4f/%.4f  f50 %.4f/%.4f  f90 %.4f/%.4f  "
            "gini %.3f/%.3f" % (key, w["zero_area"], s["zero_area"], w["ge_1pct"], s["ge_1pct"], w["ge_10pct"],
                                s["ge_10pct"], w["effective_area"], s["effective_area"], w["f50"], s["f50"], w["f90"],
                                s["f90"], w["gini"], s["gini"]))
    cols = ["zero_area", "ge_1pct", "ge_10pct", "ge_50pct", "ge_1pct_trim", "ge_10pct_trim", "effective_area", "f50",
            "f90", "top1pct_share", "gini", "max_over_mean", "mean_activity", "trimmed_mean_activity"]
    ecols = ["F_H2O_max", "F_max_over_F_E", "Q_share_over_F_E", "Q_share_over_0p1_F_E", "area_fraction_over_F_E",
             "max_F_over_F_E_all_map_times"]
    res["energy_limit"] = {"S0_W_m-2": S0_W_M2, "dH_J_mol-1": DH_J_MOL, "N_A": N_A,
                           "definition": "F_E = S0 / (r_h[AU]^2 * dH / N_A): zero albedo, normal incidence, no "
                                         "thermal emission, no conduction (an upper bound)"}
    (out / "metrics.json").write_text(json.dumps(res, indent=1))
    with open(out / "metrics.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["run"] + ["whole_" + c for c in cols] + ["seen_" + c for c in cols] + ["seen_area_fraction"]
                    + ["energy_" + c for c in ecols])
        for key in ORDER:
            r = res["runs"][key]
            wr.writerow([key] + ["%.6g" % r["whole"][c] for c in cols] + ["%.6g" % r["seen"][c] for c in cols]
                        + ["%.6g" % r["seen_area_fraction"]] + ["%.6g" % r["energy"][c] for c in ecols])
    say("wrote %s" % ", ".join(str(p) for p in finals))


if __name__ == "__main__":
    main()
