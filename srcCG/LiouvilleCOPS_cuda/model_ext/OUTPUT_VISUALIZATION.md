# OUTPUT_VISUALIZATION.md — how every figure of the 12M fits is made, and how to redo them

Kept current with every new or changed figure (user request 2026-11-12: "save instructions on how to plot and how
is plotting in OUTPUT_VISUALIZATION.md -- it should contain all necessary information to redo all the plots you
have done today and do later in case the model is rerun again; keep the file current so it can be used directly to
redo all plots next time"). Change log at the end. Paths are as on Voyager (`COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output`);
the scripts run from the copies in `$COPS_OUTPUT_ROOT` (`figures_cluster/scripts/`, `jobs/`, `model_ext/`); their
version-controlled copies are in this folder (`figures/`, `jobs/`, top level; see README.md, `RUN_COPIES.sha256`).

---

## 1. Conventions (apply to every figure)

* **Format:** `.eps` and `.png` at 300 dpi (`cops_cuda/figures.py` `style()` / `save()`; the tools refuse to overwrite
  existing files unless `--overwrite`).
* **Content:** model against the COPS observations only (no cluster-vs-passion plots). Labels in words, no internal
  model or stage codes; labels must not contain `=` where a script splits `LABEL=PATH`.
* **Late window:** shaded and labelled with its UTC and altitude range (`figures.LATE_FILL`) in the model/data
  figures; NOT in the publication NG/RG altitude figure (user 2026-11-12).
* **Ground tracks:** Rosetta impact only (no Philae).
* **Every figure folder** gets `HOW-TO-RERUN.md` (the exact commands with all arguments, the AMPS commit, the commit of
  the version-controlled code copies, job ids), `code/` (the scripts that ran + `SHA256SUMS`), and per figure a
  `_caption.txt` where the figure type has one.
* **Never overwrite a product folder:** regenerate into a new folder, or move the old files aside first
  (`<name>.superseded_<UTC>/` or `superseded_<UTC>/` with a `SUPERSEDED*.txt` that says why). Never delete.
* **Where to run:** model/data figures, tables and the publication figure on the login node (no GPU). Anything
  that rasterizes the 12M mesh or ray-casts (templates, ground tracks, map tool) inside a Slurm job (GPU node;
  `--gpus-per-node=1` is enough; a shared node is fine with `oom_score_adj 1000`, CLAUDE.md rule 7). When the user
  asks for speed, spread independent jobs over every free GPU node.
* **Environment** (login node or job):
  ```bash
  export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
  cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda       # the figure scripts import cops_cuda from here
  PY=".venv/bin/python -W ignore"; S=$COPS_OUTPUT_ROOT/figures_cluster/scripts
  RUN=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5
  ```

## 2. Inputs: where a fit's results are

* A fit stage: `<chain>/c<k>/fit/final/` (`fit_result.txt`, `COPS_model.dat`, `late_window_model.dat`,
  `fit_restart.npz`, `nonneg_facets.npz`, `nonneg_bins.txt`, `NONNEG_MODEL.txt`); a fine tuning:
  `<chain>/c<k>/ft/final/` (`finetune_result.txt`, `finetune_state.npz` with `A_final`, `COPS_model.dat`,
  `late_window_model.dat`). The final state of a chain = `c<N>/ft/final`.
* Current best fit (2026-11-13): softplus x slope multipliers >= 1 after 5 cycles, data cost 4.2861e-3,
  `$RUN/nonneg/softplus/chain_5cycles_smm1_20261113T041215Z/` (job 124664; maps 124667, geometry
  `$RUN/maps_geometry/job_124667/geometry`; contributions 124668). Before: softplus x slope multipliers >= 0 after 5 cycles
  (4.3668e-3), `$RUN/nonneg/softplus/chain_5cycles_20261112T012553Z/` (job 124310); squared: `.../nonneg/squared/chain_5cycles_20261112T012553Z/` (124313).
  Reference: clipped refit D2 `$RUN/clip_refit/D2/final`; best positive fit `$RUN/positive_r1e-5/lam1e-4/c3/ft`.
* **Observation table of a fit** (every model/data figure starts from it; drops the 401 model-only 1-s grid rows
  and takes the late-window values W10x5 from `late_window_model.dat`):
  ```bash
  $PY $S/positive_table.py <fit or ft dir> <out>.dat
  ```
* Template bundle (late window, NG-only/RG-only facets, cliff templates): `$RUN/bundle/restart_templates.npz`;
  observation-only base for the ground track/impact as in the original figures: `$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier/base`
  with `.../bundle/restart_templates.npz`.

## 3. Figure catalogue

### 3.1 Model vs observations (time, altitude, NG, RG, NG+RG with crossings)

Engine: `figures_cluster/scripts/model_obs_figures.py` (writes `<prefix>_time`, `_altitude`, `_altitude_NG`,
`_altitude_RG`, `_altitude_NG_RG`; panels NG, RG, NG+RG; full descent left, last 25 min right; late window shaded;
NG = RG crossings marked; at most 8 curves: the palette has 8 colours).
```bash
$PY $S/model_obs_figures.py --bundle $RUN/bundle/restart_templates.npz --overwrite --out <dir>/<prefix> \
    --ref <table.dat> --title "<title>" \
    --fit "<label>=<table.dat>:solid" --fit "<label>=<table.dat>:dashed" ...     # styles solid|dashed|dotted|dashdot
```
Wrapper scripts (each writes its folder under `figures_cluster/12M/`, with tables/, summary, HOW-TO-RERUN.md, code/;
login node; the chain paths are variables at the top of each script — change them for a rerun):

| script | folder | content |
|---|---|---|
| `figures_12m_nonneg9_smm1.sh` (NC=9, JOBB, JOBB2, MAPJB, CONJB; OUT) | `model_vs_observations_softplus_smm1_9cycles/` | the same chain extended to 9 cycles (jobs 124664 c1-c5, 124772 c6-c9): c9 vs c5 vs multipliers >= 0 vs D2 (`12M_softplus_smm1_9c_final_*`), the cycles split `_cycles_1-5_*` and `_cycles_5-9_*`; chain_stages.txt with the cost change per cycle, summary.txt (cost ratios, 3-column multiplier table, A = 0 fractions), MAPS.txt; refuses unless the chain's last event is c9ft rc 0 (2026-11-13; run `MAPJB=124775 CONJB=124776 figures_12m_nonneg9_smm1.sh`) |
| `figures_12m_nonneg5_smm1.sh` (TAGB, JOBB, NC; OUT) | `model_vs_observations_softplus_smm1/` | softplus with multipliers >= 1 after 5 cycles vs multipliers >= 0 vs D2; the bounded chain's cycle-1 fit and the state after each cycle; chain_stages.txt and summary.txt with termination / nfev, the parameters at a bound per stage (from fit_restart.npz x, lo, hi) and the multipliers at the bound 1; bin tables of both runs; MAPS.txt (2026-11-13) |
| `figures_12m_nonneg5.sh` (TAG, NC) | `model_vs_observations_nonneg_5cycles/` | softplus and squared after 5 cycles vs D2 and the best positive fit; per model the cycle-1 fit and the state after each cycle; chain_stages.txt, summary.txt, bin tables, MAPS.txt |
| `figures_12m_nonneg.sh` (TAG) | `model_vs_observations_nonneg/` | the 3-cycle runs; every stage; cycle 1 vs final; `12M_all_families_*` (best fit of every family) |
| `figures_12m_model_families.sh` | `model_vs_observations_model_families/` | best fit of each family (D2, positive, slope multipliers, bumps) |
| `figures_12m_slope_mult.sh`, `figures_12m_slope_bumps.sh` | `model_vs_observations_slope_mult/`, `_slope_bumps/` | slope-multiplier and bump models |
| `figures_12m_positive_cycles.sh`, `figures_12m_clip_vs_positive.sh` | `model_vs_observations_positive_cycles*/`, `_clip_vs_positive/` | positive-model cycles; clipped refits |

Also `tools/pub_figures.py pressures` (tracked, unchanged): `--base <base> --bundle <bundle> --fit "<label>=<dir>:solid" [--ref <dir>] [--no-base] --out <prefix>`
(writes `_time`, `_altitude`; keep labels short: its legend has 3 columns).

### 3.2 Publication figure: COPS and model NG/RG vs altitude

`figures_cluster/scripts/pub_ng_rg_altitude.py` — one panel, legend on the right, "Model NG/RG" labels, the gap in
the NG data hatched (consecutive NG altitudes differing by > x2 with RG data inside: 700–297 m, 10:24:38–10:33:10 UTC),
no late-window shading, no NG = RG point (user 2026-11-12; the crossing altitude is only printed).
```bash
O=$COPS_OUTPUT_ROOT/figures_cluster/12M/publication_softplus
$PY $S/positive_table.py $RUN/nonneg/softplus/chain_5cycles_20261112T012553Z/c5/ft/final $O/softplus_c5ft.dat
$PY $S/pub_ng_rg_altitude.py --table $O/softplus_c5ft.dat --out $O/12M_softplus_NG_RG_altitude [--model-label Model]
```
Bounded refit (2026-11-13): `figures_cluster/12M/publication_softplus_smm1/12M_softplus_smm1_NG_RG_altitude` (the same two
commands with `chain_5cycles_smm1_20261113T041215Z/c5/ft/final` and `softplus_smm1_c5ft.dat`; model NG = RG at 475.5 m);
after 9 cycles: `figures_cluster/12M/publication_softplus_smm1_9cycles/12M_softplus_smm1_9c_NG_RG_altitude` (`c9/ft/final`,
`softplus_smm1_c9ft.dat`; NG = RG at 472.6 m).

### 3.3 Activity and source-rate maps for ParaView / Tecplot (data files, not figures)

`jobs/nonneg_maps.sbatch` -> `model_ext/nonneg_maps.py` (+ `nonneg_maps_geom.cpp`, `nonneg_maps_shadow.cpp`): for the
last cycle of a chain, per facet: H2O activity A_final, kappa x ACO2, seen flag, slope bin, area; at 00:00 … 10:00
hourly, 10:33:10, 10:39:26 UTC: cos SZA, illumination (day / night / cast shadow), H2O and CO2 flux [m^-2 s^-1] with
the fitter's law (q = qmin at night or in cast shadow, else qmin + (qmax − qmin) cos SZA; r_h; SPICE sun; shadow ray
for every facet). Output `<chain>/maps/`: `<model>_flux.pvd` + 13 `.vtu`, `<model>_static.vtu`, `<model>_maps.dat`
(Tecplot 360 only; ParaView's Tecplot reader mis-reads shared zones), README.txt, maps_manifest.json; geometry and
build in `$RUN/maps_geometry/job_<id>/`. CPU job (normal128cores), ~6 min; validates the geometry against the base's
kernel records first.
```bash
sbatch --export=ALL,CHAINS=<chain dir>,MODEL_EXT=$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_maps_20261112T021623Z $COPS_OUTPUT_ROOT/jobs/nonneg_maps.sbatch
# after a chain: add --dependency=afterany:<chain job ids>; one chain per job (jobs 124324 softplus, 124325 squared)
```
Output names (2026-11-13, for a later state of the SAME chain; defaults unchanged): `MAPS_NAME` (default `maps`) and, for
`jobs/nonneg_contrib.sbatch`, `CONTRIB_NAME` (default `contributions`) with `CYCLE`; `EXPECT_CYCLES=<N>` refuses unless the
chain's last event is c<N>ft rc 0. The c9 state of the bounded chain: `MAPS_NAME=maps_c9 EXPECT_CYCLES=9` (job 124775,
geometry `maps_geometry/job_124775/`), `CONTRIB_NAME=contributions_c9 CYCLE=9` (job 124776); the c5 products stay in
`<chain>/maps` and `<chain>/contributions`.

### 3.4 Full figure set: templates, gravity (effective slope), ground tracks, model/data

`jobs/pub_maps_figures.sbatch` (GPU, 1 GPU, ~15 min) -> folder `figures_cluster/12M/full_set_softplus_5cycles/`
(set by OUT). Steps (STEPS, set in the environment, NOT in `--export=` because the comma would split it):

| step | folder | how |
|---|---|---|
| set1 | `visibility/` | `tools/pub_figures.py templates` and `groundtrack` (tracked, UNCHANGED) on the observation-only base `runs/shap8_012m_flat_earlier/base` + its bundle: rows late-visibility NG-only/RG-only, effective gravitational slope (gravity cache named in the bundle), cliff templates ≥ 45/≥ 60°; global + impact zoom; ground track and impact (10:39:33.8 UTC, lon 12.164°, lat 31.685°). Pixel-identical to the original figures. |
| maps (SETS=activity,production) | `activity/`, `production_rate/<UTC>/`, `checks/` | `figures_cluster/scripts/pub_maps_figures.py`: same layout with the visibility row replaced by the activity A (set 2) or by H2O, CO2 and total production rate at 10:33:10 and 10:39:26 (set 3); ground tracks coloured accordingly; fields bitwise-checked against the maps' .vtu; track/impact from the observation-only base. Colour ranges (user 2026-11-12): log10 A 13..23, F_H2O and total 12..22, F_CO2 −28..−18 m^-2 s^-1; non-zero below range green, zero grey, above range = top colour. |
| set4, pressures | `model_data/` | tables (positive_table.py), model_obs_figures.py (fit + D2), tools/pub_figures.py pressures, copy of the publication figure |

```bash
cd $COPS_OUTPUT_ROOT/jobs
OUT=$COPS_OUTPUT_ROOT/figures_cluster/12M/<new folder> STEPS=set1,set4,pressures,maps sbatch --nodelist=<gpu node> --export=ALL pub_maps_figures.sbatch
```
For another fit (2026-11-13: inputs from the environment, defaults = the first fit): `CHAIN`, `NCYC` (or `FT`, `FIT`),
`MAPS`, `GEOM` (`$RUN/maps_geometry/job_<map job>/geometry`), `PUB`, `PUBSTEM`, `PUBTABLE`, `TAG` (file names
`12M_<TAG>_vs_COPS_*`, `_pressures_*`, `tables/<TAG>_final.dat`), `MODEL_LABEL`, `MODEL_SHORT`, `TITLE_PRESSURES` (labels
must not contain "=": use "≥"). The `maps` step writes `activity/` and `production_rate/` (fixed ranges) and the set-1
ground-track captions; it collides with the version-2 folders, so for a complete set in ONE folder run
`STEPS=set1,set4,pressures` plus the version-2 jobs (3.5), write the model-data captions with
`figures_cluster/scripts/model_data_captions.py` and copy the fit-independent set-1 ground-track captions after checking that
the set-1 PNGs are pixel-identical (done for `full_set_softplus_smm1_5cycles/`, job 124722, CAPTIONS_COPIED.txt).
The c9 state of the bounded chain (2026-11-13): `full_set_softplus_smm1_9cycles/` (3.4 job 124843 with NCYC=9,
MAPS=<chain>/maps_c9, GEOM=job_124775; version-2 jobs 124777-124783 with MAPS/CONTRIB of c9, OUT=the same folder,
cross-check `<chain>/v2_total_crosscheck_c9`; `comparison_vs_c5/` 124784 and `comparison_vs_unbounded/` 124785).

### 3.5 Version 2 of the activity / source-rate maps, contribution maps, cliff plot (IN PROGRESS 2026-11-12)

User request 2026-11-12: final activity; discrete log scale with 15–20 intervals over the entire non-zero range on the
seen facets within the shown window; colour only on the facets seen during the fit (the fit's seen flag), grey
elsewhere; the same for the TOTAL production rate (no H2O/CO2 split); template figures without the top text block and
each panel also as its own file; per post-gap observation (95, 10:33:10–10:39:26) the smallest highest-contribution
facet set supplying 50/60/70/80/90 % of the model's RG (and NG) surface signal, union maps and per-facet count maps,
NG and RG side by side and separately; all captions also in one separate file; a plot showing the importance of the
cliffs. Folder: `figures_cluster/12M/full_set_softplus_5cycles_v2/`; contribution data: `<chain>/contributions/`.
**Figures (done 2026-11-12; jobs 124565-124569, gpu002, 1 GPU and 40 GB each, ~3-4 min each):** 49 figures, each
`.eps` + `.png` + `_caption.txt`; all captions collected in `CAPTIONS.md`; colour-class tables in `checks/`.

| folder | figures | script |
|---|---|---|
| `activity/` | `12M_templates_activity` (rows A, slope, cliffs; no top text block), `panels/12M_templates_activity_{a..f}_*` (each panel alone), `12M_groundtrack_activity_{latlon,3d}` | `figures_cluster/scripts/pub_maps_v2_figures.py --quantity activity` |
| `production_rate/` | the same for the total rate F = F_H2O + F_CO2 at `103310UTC` and `103926UTC` | `pub_maps_v2_figures.py --quantity total --time HH:MM:SS` |
| `source_areas/` | `12M_source_area_union_{NG,RG,NG_RG}` (unions S_50 ⊂ … ⊂ S_90), `12M_source_area_occurrence_P{50..90}_{NG,RG,NG_RG}` (N_P classes 1, 2–4, 5–9, 10–19, 20–39, 40–59, 60–79, 80–94, 95); global map + close-up lon 11.7–18.8°, lat 30.3–34.0° | `figures_cluster/scripts/pub_source_area_figures.py` |
| `cliffs/` | `12M_cliff_slope_bins` (per slope bin: area and seen fractions, multiplier, production shares), `12M_cliff_pressure_fraction_{NG,RG,NG_RG}` (stacked share of the model surface pressure per slope-bin group vs gauge altitude; NG gap hatched; late window boxed) | `figures_cluster/scripts/cliff_importance_figures.py` (CPU) |

Colours (user request: colour only on seen facets): seen > 0 'plasma' at the class centres (discrete log10 classes,
round boundaries, 15–20 classes over the entire non-zero range of the seen facets in each displayed window, scanned
per panel; each window has its own colour bar); seen and exactly 0: #A1D76A (relief-shaded); not seen: grey #CCCCCC
(relief-shaded). Activity: 18 classes of 20 decades, 1e-320 … 1e40 (the non-zero seen activity spans 3.5e-315 …
1.9e21: softplus underflow x slope multipliers 1e-22…1e-31; denormals below 2.2e-308 stated in the captions); total
rate: 20 (global) / 18 (impact area) classes of 2.5 decades, 1e-27.5 … 1e22.5. Track and impact from the
observation-only base (as set 1 of 3.4). `checks/` verifies: fields bitwise = the maps' .vtu; every value inside its
class and class counts complete; unions = the contribution files and nested; cliff shares = the verified fields.

Commands (from `<AMPS>/srcCG/LiouvilleCOPS_cuda`; exact ones with every argument in the folder's `HOW-TO-RERUN.md`):
```bash
P=".venv/bin/python -W ignore"; S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; R=$RUN
CH=$R/nonneg/softplus/chain_5cycles_20261112T012553Z; TB=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier
MX=$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_maps_20261112T021623Z; OUT=$COPS_OUTPUT_ROOT/figures_cluster/12M/full_set_softplus_5cycles_v2
$P $S/pub_maps_v2_figures.py --base $R/base --track-base $TB/base --bundle $TB/bundle/restart_templates.npz \
   --ft-dir $CH/c5/ft/final --fit-dir $CH/c5/fit/final --maps-dir $CH/maps --geometry $R/maps_geometry/job_124324/geometry \
   --maps-code $MX/nonneg_maps.py --out-dir $OUT --quantity activity --classes 15,20 --device 0          # GPU
#   ... the same with --quantity total --time 10:33:10   and   --quantity total --time 10:39:26
$P $S/pub_source_area_figures.py --base $R/base --track-base $TB/base --bundle $TB/bundle/restart_templates.npz \
   --contrib-dir $CH/contributions --maps-code $MX/nonneg_maps.py --out-dir $OUT --device 0             # GPU
$P $S/cliff_importance_figures.py --contrib-dir $CH/contributions --base $R/base \
   --bundle $TB/bundle/restart_templates.npz --v2-dir $OUT --out-dir $OUT/cliffs                         # after the 'total' runs
$P $S/collect_captions.py $OUT --title "..."                                                             # CAPTIONS.md
```
As Slurm jobs (set every variable in the submitting shell, not inside `--export=`: a comma there splits the value):
```bash
cd $COPS_OUTPUT_ROOT/jobs
QUANTITY=activity                sbatch --nodelist=<gpu> --mem=40G --export=ALL pub_maps_v2.sbatch
QUANTITY=total TIME=10:33:10     sbatch --nodelist=<gpu> --mem=40G --export=ALL pub_maps_v2.sbatch
QUANTITY=total TIME=10:39:26     sbatch --nodelist=<gpu> --mem=40G --export=ALL pub_maps_v2.sbatch
STEPS=sources,cliffs             sbatch --nodelist=<gpu> --mem=40G --export=ALL pub_source_areas.sbatch   # after the three above
```
For another fit set (environment; defaults = today's values; each job writes them to `provenance-<job>/inputs.txt`):
`RUN`, `CHAIN` + `NCYC` (or `FT`, `FIT`), `MAPS`, `GEOM` (geometry of the chain's maps job), `MX`, `CONTRIB`
(`<chain>/contributions` from the contribution job above), `TBASE`, `TBUNDLE`, `OUT`, `QUANTITY`, `TIME` (a map time
of the maps job), `CLASSES`, `STEPS`.
Order after a rerun: maps job (3.3) -> contribution job -> the three `pub_maps_v2` jobs (in parallel, any free GPUs)
-> `pub_source_areas` -> `collect_captions.py`.

Colour-range and relief options (user decisions 2026-11-12; CURRENT figures: jobs 124589 activity, 124594/124595 H2O rate,
124576 cliffs, 124565 source areas): `pub_maps_v2_figures.py --low-cut-area F` (per displayed window the classes run from the
area-weighted F-quantile of the seen non-zero values to their maximum; 0 = entire range, 0.001 = the 99.9 % rule used now;
values below take the lowest class colour, triangle on the colour bar, counts in the caption; class steps: round 0.1 … 2.5,
5, 10, 20, 25, 50, 100 decades, else 3, 4, 6, 8, 12, 15, 30, 40) and `--template-relief S` (template value row and its
panels: class colour darkened/lightened by the hillshade, unchanged on level ground; 0 = flat, 0.7 used for the total rate;
activity flat). Job variables: `LOWCUT` (default 0.001), `TEMPLATE_RELIEF` (default 0.7 for QUANTITY=total, 0 for activity).
Cliff descent figure: RG shares shown through the NG data gap, gap edges dashed and labelled (built in, no option).
H2O production rate (user 2026-11-12, replaces the total in the map figures): `--quantity h2o --time HH:MM:SS` (F_H2O =
q_H2O A; file stems `h2o_production_rate_<hhmmss>UTC`; zero-activity seen facets exactly 0 = green). Job defaults:
`LOWCUT` 0.02 for activity and h2o/total (user 2026-11-12); `TEMPLATE_RELIEF` 0.7 for h2o/total, 0 for activity. Preview with another
cut: `QUANTITY=h2o TIME=... LOWCUT=0.02 OUT=<preview folder> sbatch ... pub_maps_v2.sbatch`. Captions: `collect_captions.py
$OUT --skip preview_ --title "..."`. The cliff job's cross-check reads `checks/slope_bins_total_*.json` of a
`QUANTITY=total` run (set `V2` to a folder with them, or run QUANTITY=total first).
```bash
cd $COPS_OUTPUT_ROOT/jobs; export OUT=<folder>
QUANTITY=activity                sbatch --nodelist=<gpu> --mem=32G --export=ALL pub_maps_v2.sbatch     # --low-cut-area 0.02
QUANTITY=h2o TIME=10:33:10       sbatch --nodelist=<gpu> --mem=32G --export=ALL pub_maps_v2.sbatch     # --low-cut-area 0.02 --template-relief 0.7
QUANTITY=h2o TIME=10:39:26       sbatch --nodelist=<gpu> --mem=32G --export=ALL pub_maps_v2.sbatch
STEPS=sources,cliffs             sbatch --nodelist=<gpu> --mem=40G --export=ALL pub_source_areas.sbatch  # after the two total jobs
```

**Bounded refit (2026-11-13), complete set in ONE folder `figures_cluster/12M/full_set_softplus_smm1_5cycles/`:** set 1, set 4
and pressures (3.4, job 124722), activity (124677), H2O rate (124678, 124679), source areas (124682), cliffs (124683, groups
`all`: one colour per slope bin since no bin can be switched off), CAPTIONS.md (59 captions), README.txt with this fit's
facts, `comparison_vs_unbounded/` (`figures_cluster/scripts/source_area_compare.py` via `jobs/source_area_compare.sbatch`:
unions, area-weighted locations, distance bands from the last gauge position, Jaccard overlaps and slope-bin pressure shares
of two contribution folders; job 124727). The total-rate cross-check of the cliff job: QUANTITY=total into the non-figure
folder `<chain>/v2_total_crosscheck` (124680, 124681), then `STEPS=cliffs V2=<chain>/v2_total_crosscheck`.
New options (defaults keep the earlier captions): `MODEL_NAME` / `--model-name` (v2 maps, source areas, cliffs: the name of
the fitted state in the captions); the multiplier-bound note (from fit_restart.npz names in meta_json, x, lo), the CO2
sentence (maps_manifest Q_CO2 / Q_H2O) and the contributions' job id and verification (contrib_manifest.json) are computed;
cliffs `--fit-dir`, `--maps-dir`, `--groups default|all` (job variables `NCYC`/`FIT`, `MAPS`, `CLIFF_GROUPS`).
```bash
cd $COPS_OUTPUT_ROOT/jobs; export CHAIN=<chain> GEOM=$RUN/maps_geometry/job_<map job>/geometry OUT=<folder> MODEL_NAME="<name>"
DEP="--dependency=afterok:<contrib job> --kill-on-invalid-dep=yes --exclude=gpu003"
QUANTITY=activity sbatch $DEP --mem=24G --export=ALL pub_maps_v2.sbatch;  QUANTITY=h2o TIME=10:33:10 sbatch ...; QUANTITY=h2o TIME=10:39:26 sbatch ...
QUANTITY=total TIME=10:33:10 OUT=$CHAIN/v2_total_crosscheck sbatch ...;  QUANTITY=total TIME=10:39:26 OUT=$CHAIN/v2_total_crosscheck sbatch ...
STEPS=sources sbatch $DEP --mem=32G --export=ALL pub_source_areas.sbatch
STEPS=cliffs V2=$CHAIN/v2_total_crosscheck CLIFF_GROUPS=all sbatch --dependency=afterok:<total jobs> --mem=16G --export=ALL pub_source_areas.sbatch
```

**Contribution data (done 2026-11-12, CPU job 124549, 3 min 20 s):** `jobs/nonneg_contrib.sbatch` ->
`model_ext/nonneg_contrib.py` (frozen `model_ext_frozen/nonneg_contrib_20261112T210012Z`). For the final state of a
chain: each kernel record's NG and RG contribution as the fitted model computes it (the fit's FP32 record fields via
`cops_cuda.device_store.fp32_record_fields`, FP64 arithmetic, F = q A with the night/shadow rule, the fit's T and
kappa, the fit_cops gauge constants, the chain's W10x5 for the 95 late observations); NG grain term not attributed
to facets. Verified: per-row sums = the model tables to <= 4.7e-10 (print precision), late sums = late_window_model.dat
to 2.2e-10; all contributions >= 0. Names (proposed 2026-11-12, after Schmid 1994): **P% source area** S_P(k, g) =
the smallest set of facets, in decreasing order of contribution (ties: lower facet index), supplying P % of
observation k's modelled surface signal of gauge g; **source-area occurrence** N_P,g(e) = the number of the 95
post-gap observations whose S_P contains facet e. Products `<chain>/contributions/`: `source_areas_NG.npz`,
`source_areas_RG.npz` (per observation totals, grain, per P counts/areas/thresholds, the full ranking, occurrence
5 x nfaces, union masks), `slope_bins.npz`/`.txt` (per base row the 7 slope-bin sums and fractions of NG and RG; per
bin multiplier, areas, production shares at 10:33:10 and 10:39:26), `source_areas_per_observation.txt`, README.txt.
```bash
sbatch --export=ALL,CHAIN=<chain dir>,CYCLE=<N>,GEOM=$RUN/maps_geometry/job_<map job>/geometry,MODEL_EXT=$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_contrib_20261112T210012Z \
    $COPS_OUTPUT_ROOT/jobs/nonneg_contrib.sbatch          # CPU (normal128cores); GEOM: the validated geometry of the chain's maps job (3.3)
```

### 3.7 Movies: the nucleus as seen from RG and NG, with the time series (2026-11-12)

User request 2026-11-12 (1-min movies: the nucleus seen from RG, from NG, and a combined movie with the model/COPS time
series and a current-time line; resolution for projectors and online). Answers: whole descent 00:00:04–10:39:26 UTC;
camera-like sunlit view with cast shadows; the same view (camera at the gauge = Rosetta + GaugeOffsetSC, looking along the
RG boresight, image up = NG boresight) for both gauges, tinted by each gauge's acceptance weight as used in the fit (RG
-l.b_RG, 0..1; NG 1 - l.b_NG in front, the shield fraction behind, 0..2); frame times uniform in log altitude (16.7 km ->
5.5 m; the late window is the second half); NG and RG time series stacked with a current-time line. Perspective while the
nucleus is within 60 deg of the RG axis (frames 0-687), azimuthal-equidistant fisheye after (labelled on the frame).
Products `figures_cluster/12M/movies_softplus_5cycles/`: `movie{1_RG_view,2_NG_view,3_combined}_1920x1080.mp4` (H.264,
30 fps, CRF 18) and `_1280x720_web.mp4` (CRF 26, faststart), 60 s, 1800 frames; `stills/`, `validation/`, README.txt,
CAPTIONS.md, HOW-TO-RERUN.md, code/. Scripts: `figures_cluster/scripts/cops_movies.py`, `cops_movies_geom.cpp`; jobs
`jobs/cops_movies_{prep,render,encode,qc}.sbatch`; ffmpeg: user-space `cops_output/tools/imageio_ffmpeg` (outside the
.venv). Work folder `<chain>/movies_work/` (frames, composed frames, plan). Render: 8 GPUs ~1:14 wall (0.29 s/frame);
the result is independent of the GPU split (checked bitwise with 4 GPUs).
```bash
cd $COPS_OUTPUT_ROOT/jobs
export WORK=<chain>/movies_work_<new> OUT=$COPS_OUTPUT_ROOT/figures_cluster/12M/<new movie folder>
STEPS=prep,plan,validate sbatch --nodelist=<gpu> --export=ALL cops_movies_prep.sbatch
NGPU=8 WORKER0=0 NWORKERS=8 sbatch --nodelist=<gpu> --gpus-per-node=8 --cpus-per-task=32 --mem=96G --export=ALL cops_movies_render.sbatch
STEPS=compose COMPOSED=$WORK/composed_samples FRAMES=0,688,1799 sbatch --export=ALL cops_movies_encode.sbatch   # sample stills
STEPS=compose,encode PROCS=32 sbatch --export=ALL cops_movies_encode.sbatch
sbatch --export=ALL cops_movies_qc.sbatch
```
**Observable-part movies with the spacecraft shadow (2026-11-13; folder `figures_cluster/12M/movies_softplus_smm1_5cycles/`,
work `<chain smm1>/movies_work`; user decisions in section 5):** colour only inside each gauge's IK field of view (RG ±45° x ±45°,
NG ±77.5° x ±77.5°, getfov, tested as Geometry.cpp InsideFOV), tinted by the fit's acceptance weight (colour bars 0..1 /
0..2 kept), sunlit grey outside, the aperture boundary as a white line; the spacecraft (bus box + two solar-array wings,
ROS_STRUCT_V8 / ROS_SA_2016_V0042, SA orientation held after 10:39:19.851 UTC) shades the nucleus with the Sun as a disc
(64 samples), movies only. New steps: `fov` (fov.npz), `validate_fov`, `sc` (sc.npz; `--sa-hold-utc`, `--origin-offset`:
the structure is relative to ROS_SPACECRAFT's origin -226000, the launch-vehicle interface point; no kernel links it to
ROSETTA -226 (centre of gravity), offset 0 used), `validate_sc`, `frametable`, `diagnose` (dark-pixel classification with
the exact shadow ray; `jobs/cops_movies_diag.sbatch`); render options `--mask-rg/--mask-ng none|fit|hemisphere|ik_fov|cone:<deg>`,
`--tint-rg/--tint-ng weight|uniform`, `--sc-shadow --sc-samples N`, `--no-boundary`; compose `--ts-label`. Defaults (none)
reproduce the first movies bit for bit (job 124707). Share of the modelled signal from outside the IK FOVs:
`figures_cluster/scripts/fov_share.py` via `jobs/fov_share.sbatch` (CPU, the frozen contribution code; needs
NPY_DISABLE_CPU_FEATURES as nonneg_contrib.sbatch).
```bash
cd $COPS_OUTPUT_ROOT/jobs; export RUN=... CHAIN=<chain> GEOM=$RUN/maps_geometry/job_<map job>/geometry WORK=$CHAIN/movies_work
STEPS=prep,fov,plan,validate,validate_fov sbatch --exclude=gpu003 --export=ALL cops_movies_prep.sbatch
STEPS=sc,validate_sc sbatch --dependency=afterok:<previous> --exclude=gpu003 --export=ALL cops_movies_prep.sbatch
MASK_RG=ik_fov MASK_NG=ik_fov SC_SHADOW=1 SC_SAMPLES=64 NGPU=8 WORKER0=0 NWORKERS=8 sbatch --nodelist=<gpu> --gpus-per-node=8 --cpus-per-task=32 --mem=120G --export=ALL cops_movies_render.sbatch
export OUT=<movie folder> TABLE=<the fit's positive_table .dat> TS_LABEL='COPS data and the\nsoftplus activity model,\nmultipliers ≥ 1\n(final state, 5 cycles)'
STEPS=frametable,compose,encode PROCS=32 sbatch --export=ALL cops_movies_encode.sbatch; sbatch --export=ALL cops_movies_qc.sbatch
OUTP=$OUT/validation/fov_share/<prefix> sbatch --export=ALL fov_share.sbatch
WORK=<first movies' work> GEOM=... OUTD=<folder> BOXES="1500,700,830,330,550" FRAMES=1500,1799 sbatch --exclude=gpu003 --export=ALL cops_movies_diag.sbatch
```
Variables: all scripts RUN, CHAIN, WORK; prep GEOM, STEPS (prep, fov, plan, validate, validate_fov, sc, validate_sc),
SA_HOLD_UTC; render GEOM, NGPU, WORKER0, NWORKERS, FRAMES_DIR, ONLY, MASK_RG, MASK_NG, TINT_RG, TINT_NG, SC_SHADOW,
SC_SAMPLES, NO_BOUNDARY; encode also TS_LABEL, FRAMETABLE (STEPS frametable);
encode STEPS (compare, compose, encode), COMPOSED, FRAMES, PROCS, OUT, TABLE (the model table for the time series),
CMP_A, CMP_B; qc STILLS. For another fit: CHAIN (camera, illumination are fit-independent; the time-series TABLE and the
maps GEOM change). Split check: render again with FRAMES_DIR=$WORK/frames_check NGPU=4 and STEPS=compare.

**c9 movies, fixed camera and the touchdown flag (2026-11-13; user decisions in section 5).** Two new sets of the bounded
chain's c9 state: `figures_cluster/12M/movies_softplus_smm1_9cycles/` (gauge view as the c5 observable-part movies;
work `<chain>/movies_work_c9`, jobs 124786-124788, touchdown 124863, compose/encode 124866, qc 124867, fov_share 124789)
and `figures_cluster/12M/movies_softplus_smm1_9cycles_fixedcam/` (work `<chain>/movies_work_c9_fixedcam`, jobs 124864
prep/plan/validate_fixed/touchdown/validate_sc, 124865 render, 124868 compose/encode, 124869 qc). Both with the waving
touchdown flag. New in `cops_movies.py` (defaults unchanged; the earlier movies reproduce bit for bit: default render, plan
and compose paths checked in `tests/movies_fixedcam_touchdown_20261113/`, jobs 124819, 124844):
* `plan --camera fixed --camera-altitude 100 --fov-deg 90 --base --kernels --geom-exe`: the camera = COPS position and
  RG / NG boresights at the time the COPS altitude is 100 m (10:37:31.80 UTC; interpolated in log altitude as the frame
  times; SPICE as prep), frozen in 67P/C-G_CK; perspective, 90° horizontal; `camera_checks.json`.
* render with a fixed-camera plan (masks `ik_fov` required): per frame and gauge, from the TRUE gauge position and
  attitude, a camera pixel is tinted when its surface point is inside the gauge's IK FOV, its facet faces the gauge and the
  segment to the gauge crosses no facet (perspective depth map from the gauge over the window of the selected directions,
  exact segment tests against the front-most facets of the 3 x 3 texels and the near-plane facets, then exact first-hit
  rays along every edge of the unseen parts, iterated: `Renderer.gauge_los`, `refine_los`); sunlit relief elsewhere; the
  FOV boundary on the surface as the white line; "+" the sub-spacecraft point, a ring the gauge position (both depth-tested
  from the camera). Validation `validate_fixed` (FOV vs the instrument-frame rectangle, line of sight vs exact rays at
  random and boundary samples, diagnostic images); `frametable` writes `frames.tsv.camera_check.json` (camera and face
  buffer identical in all frames).
* `touchdown --track-base runs/shap8_012m_flat_earlier/base --impact-json <ground-track _impact.json> [--kernels
  --geom-exe]`: the ground-track figures' impact point (figures.impact_point, compared: 0.0 m), per frame its pixel, the
  pole direction and the visibility from the camera; SPICE check (RORB trajectory crossing the 12M surface) ->
  `touchdown.npz/.json`. `compose --touchdown <touchdown.npz>`: the flag and its legend.
* job variables: prep `FROM_WORK` + STEPS `copy`, `CAMERA=fixed CAMERA_ALT FOV_DEG`, STEPS `validate_fixed` (MASK_RG,
  MASK_NG, SAMPLES), `touchdown` (TRACK_BASE, IMPACT_JSON); encode `TOUCHDOWN`; all three `MOVIES_PY` (the cops_movies.py
  to run; the c9 jobs ran the frozen copy `model_ext_frozen/cops_movies_fixedcam_20261113T215921Z/`, identical to the
  installed script).
```bash
cd $COPS_OUTPUT_ROOT/jobs; export RUN=... CHAIN=<chain> GEOM=$RUN/maps_geometry/job_<map job>/geometry
export IMPACT_JSON=<full set>/visibility/12M_groundtrack_visibility_impact.json
WORK=$CHAIN/movies_work_<gauge> STEPS=touchdown sbatch --exclude=gpu003 --export=ALL cops_movies_prep.sbatch       # after its plan
WORK=$CHAIN/movies_work_<fixed> FROM_WORK=$CHAIN/movies_work_<gauge> CAMERA=fixed CAMERA_ALT=100 FOV_DEG=90 \
  STEPS=copy,plan,validate_fixed,touchdown,validate_sc sbatch --exclude=gpu003 --mem=64G --export=ALL cops_movies_prep.sbatch
WORK=$CHAIN/movies_work_<fixed> MASK_RG=ik_fov MASK_NG=ik_fov SC_SHADOW=1 SC_SAMPLES=64 NGPU=8 WORKER0=0 NWORKERS=8 \
  sbatch --nodelist=<gpu> --gpus-per-node=8 --cpus-per-task=32 --mem=160G --export=ALL cops_movies_render.sbatch
WORK=... OUT=<movie folder> TABLE=<the fit's table> TS_LABEL='...' TOUCHDOWN=$WORK/touchdown.npz STEPS=frametable,compose,encode \
  PROCS=32 sbatch --export=ALL cops_movies_encode.sbatch; WORK=... OUT=... sbatch --export=ALL cops_movies_qc.sbatch
```

**Detach-200 m set (2026-11-13/14; user decisions in section 5):** `figures_cluster/12M/movies_softplus_smm1_9cycles_detach200/`
(work `<chain>/movies_work_c9_detach200`; jobs 124894 contact/prep/fov/sc/plan/validate_fixed/touchdown/validate_sc,
124895 render, 124896 frametable/compose/encode, 124897 qc; fov_share with circles 124885). Circular fields of view
(`--mask-rg cone:45 --mask-ng cone:77.5`); the camera at COPS (gauge view) until the first frame at or below 200 m
(frame 992, 10:35:12.2 UTC), then frozen there (perspective 90°); the spacecraft model (bus, wings as five panels each,
HGA disc; Lambert; near clipping 1 m), its shadow incl. the HGA, the translucent FOV cones (24 generators) and the
footprints; frames after the last observation uniform in time to the first contact of the model with the surface
(10:39:30.624 UTC, -Y wing). New in `cops_movies.py` (defaults unchanged, tested bit for bit, jobs 124887/124892 + CPU
compose): `contact`, `prep --extend-to-utc`, `sc --hga-ck --hga-hold-utc`, `plan --camera detach --detach-altitude`,
the detach render; `cops_movies_geom hgageom`; `fov_share.py --cone-rg/--cone-ng` (job variables CONE_RG/CONE_NG).
HGA C-kernel: `cops_output/kernels_extra/20261113T232248Z/ROS_HGA_2016_V0043.BC` (ESA SPICE archive; SOURCE.txt).
```bash
export WORK=$CHAIN/movies_work_<new> HGA_CK=<ROS_HGA_2016_V0043.BC> HGA_HOLD_UTC=2016-09-30T10:39:19.851 IMPACT_JSON=<...>
STEPS=contact,prep,fov,sc,plan,validate_fixed,touchdown,validate_sc EXTEND_TO_CONTACT=1 CAMERA=detach DETACH_ALT=200 FOV_DEG=90 \
  MASK_RG=cone:45 MASK_NG=cone:77.5 TEST_FRAMES=<frames after the detach> sbatch --exclude=gpu003 --mem=80G --time=02:00:00 --export=ALL cops_movies_prep.sbatch
MASK_RG=cone:45 MASK_NG=cone:77.5 SC_SHADOW=1 SC_SAMPLES=64 NGPU=8 WORKER0=0 NWORKERS=8 sbatch --nodelist=<gpu> --gpus-per-node=8 \
  --cpus-per-task=32 --mem=200G --export=ALL cops_movies_render.sbatch
OUT=... TABLE=... TS_LABEL=... TOUCHDOWN=$WORK/touchdown.npz STEPS=frametable,compose,encode PROCS=32 sbatch --export=ALL cops_movies_encode.sbatch
OUTP=<folder>/validation/fov_share/<prefix> CONE_RG=45 CONE_NG=77.5 FOV=$WORK/fov.npz NCYC=9 sbatch --export=ALL fov_share.sbatch
```
**Detach-300 m set (2026-11-14 UTC; user decisions in section 5):** `figures_cluster/12M/movies_softplus_smm1_9cycles_detach300/`
(work `<chain>/movies_work_c9_detach300`, prep/fov/sc/contact copied from the detach-200 m work; jobs 124909 copy/plan/
validate_fixed/touchdown/validate_sc with `DETACH_ALT=300`, 124910 render, 124911 compose/encode with `FLAG_ON_TOP=1`,
124912 qc). As the detach-200 m set but the camera detaches at the first frame at or below 300 m (frame 901, 10:33:05.7
UTC, 299.98 m), and the touchdown flag is drawn on top (`compose --flag-on-top`: no depth test against the nucleus or the
spacecraft, wherever its point projects into the image; tested: without the option the earlier frames reproduce bit for
bit). Distances on the surface at the last frame (10:39:30.624): sub-bus point 0.16 m from the flag, the bus's first
contact (10:39:30.753) 1.29 m, the wing's first contact 16.3 m (`tests/movies_detach300_20261113/flag_distances.json`).

**Panorama set (2026-11-14 UTC; user decisions in section 5):** `figures_cluster/12M/movies_softplus_smm1_9cycles_detach300_pan/`
(work `<chain>/movies_work_c9_detach300_pan`, prep/fov/sc/contact copied from the detach-300 m work; jobs 124930 copy/plan/
validate_fixed/touchdown/validate_sc with `PAN_TO=analysis/poses/pan_B.json PAN_END=1800`, 124931 render (8 A100, 9:45),
124936 frametable/compose/encode with `FLAG_ON_TOP=1 NO_SUB_MARK=1`, 124937 qc). As the detach-300 m set, but after the
detach (frame 901, 300 m) the camera moves to pose B (120 m north-north-west of the touchdown point, 220 m above it,
looking toward azimuth 99°, 57° down, level horizon, field of view 105.5°): `plan --camera detach --pan-to <pose.json>
[--pan-start-frame K] [--pan-end-frame K]` = smootherstep in the frame index, position on an arc around the pivot (the
touchdown point; direction slerp, distance geometric), orientation by quaternion slerp, focal length geometric; the final
pose reached at frame 1800 (10:39:26.03 UTC) and held to the end; checks in `camera_checks.json` "pan" (clearance, pivot
and COPS in the image, steps) and the pan check of `frametable` (each frame the plan's pose; one face buffer per pose).
`compose --no-sub-mark`: no sub-spacecraft "+". Analysis behind the pose (`analysis/`): `ng_source_map.py` (per-facet NG
source sums, late and all observations; job 124913), `ng_source_regions.py` (lobes by the neck plane, shares by lobe,
distance and azimuth, skyline and relief, map figure), `ng_lobe_share_vs_altitude.py` (big-lobe share per late
observation), `pan_candidates.py` (grid of final poses with target framing, visible shares of the big-lobe and relief NG
sources; jobs 124915, 124918), candidate stills A/B/C (`cops_movies_pan_stills.sbatch`, jobs 124919, 124920). Defaults of
`cops_movies.py` unchanged (compose, plan and render bit for bit: `validation/regression/`).
```bash
STEPS=copy,plan,validate_fixed,touchdown,validate_sc CAMERA=detach DETACH_ALT=300 FOV_DEG=90 PAN_TO=<pose.json> PAN_END=1800 \
  MASK_RG=cone:45 MASK_NG=cone:77.5 TEST_FRAMES=<frames along the path> FROM_WORK=$CHAIN/movies_work_c9_detach300 \
  sbatch --nodelist=gpu002 --mem=80G --export=ALL cops_movies_prep.sbatch
# render as the detach-300 m set; then
FLAG_ON_TOP=1 NO_SUB_MARK=1 OUT=... TABLE=... TS_LABEL=... TOUCHDOWN=$WORK/touchdown.npz STEPS=frametable,compose,encode \
  PROCS=32 sbatch --export=ALL cops_movies_encode.sbatch
```

### 3.6 Other analysis figures and tables

* `figures_cluster/scripts/positive_table.py`, `positive_cycles_table.py`, `finetune_table.py`, `scan_summary.py`:
  tables used by the wrappers.
* `analysis/shadow_zero_activity_12m/shadow_classes.py` (CPU job): illumination of the seen surface with and without
  activity (no figure).

## 4. Redo all figures after a model rerun (recipe)

1. Identify the new chain(s): `<chain>/c<N>/ft/final` complete (`finetune_result.txt`), and the chain job id(s).
2. **Maps** (3.3): `sbatch --export=ALL,CHAINS=<chain>,MODEL_EXT=<frozen nonneg_maps> jobs/nonneg_maps.sbatch`; note the
   job id (its geometry folder `$RUN/maps_geometry/job_<id>/geometry` is an input of 3.4).
3. **Model/data figures** (3.1): copy the closest wrapper (e.g. `figures_12m_nonneg5.sh`) to a new name, set `TAG`/`NC`
   (and the job ids in its `JOB`/`MAPJ` arrays) and the output folder; run it on the login node.
4. **Publication figure** (3.2): the two commands with the new `c<N>/ft/final`, into a new folder (or move the old
   figure to `superseded_<UTC>/`).
5. **Full set** (3.4) and **version 2** (3.5): edit the input block of `jobs/pub_maps_figures.sbatch` (or its v2
   successor), new OUT folder, submit on a free GPU node.
6. **Check** every figure visually (overlaps, empty panels, legends); for maps compare the plotted fields with the
   maps' .vtu (the jobs do this: `checks/`).
7. **Record:** HOW-TO-RERUN.md + code/ in every folder (the scripts write them), an entry in
   `figures_cluster/12M/STATUS.md` and in `/data001/highresolu_vt/AMPS/handoff.md` §10, version-controlled copies of
   new or changed scripts here (header + `RUN_COPIES.sha256` row + README row; `./check_copies.sh`), and update THIS file.

## 5. User decisions on figures (dated)

* 2026-10-27: model against observations only; a panel with NG and RG together (to see when they flip); no ratio panel.
* 2026-11-03/06: HOW-TO-RERUN.md with exact commands, AMPS and code commits, and code/ next to every comparison plot set.
* 2026-11-11: production-rate maps = instantaneous H2O and CO2 flux as the model computes it, at 00:00…10:00 hourly,
  10:33:10, 10:39:26 UTC; final state only; VTK and Tecplot.
* 2026-11-12: full set for the softplus 5-cycle fit; production-rate times 10:33:10 and 10:39:26; log colour scale;
  10-decade ranges (A 1e13–1e23, F 1e12–1e22, CO2 1e-28–1e-18), below range green, zero grey; set 1 on the
  observation-only base (95 points, impact 10:39:33.8 UTC), the same track/impact in sets 2–3.
* 2026-11-12: publication NG/RG figure: legend right; "Model NG/RG"; no late-window shading; the NG data gap marked
  and named; no Model NG = RG point.
* 2026-11-12: version 2 (3.5): seen-only colour, discrete 15–20 log intervals over the entire range in the shown
  window, total rate only, no top text, panels as separate files, contribution maps NG/RG, cliff plot, captions file.
* 2026-11-12: activity AND total production rate: the 15–20 log classes cover the range holding 99.9 % of the seen
  non-zero area in each window (no separate below-range class: values below take the lowest class colour);
  cliff plot: the RG slope-bin shares filled in through the NG data gap, the gap's edges marked.
* 2026-11-12: total production-rate template figures: class colours modulated by the relief shading so the nucleus
  topography stays visible (colour bar = unshaded class colours); the other panels and figures unchanged.
* 2026-11-12: activity colour range: 2 % area cut (--low-cut-area 0.02); production-rate maps show the H2O rate
  F_H2O = q A only (no CO2; zero-activity facets exactly 0), with the 99.9 % rule.
* 2026-11-12: production rate also with the 2 % area cut (--low-cut-area 0.02).
* 2026-11-12/13: figures of the softplus refit with slope multipliers >= 1 (user "yes, make maps, contributions and
  comparison plots"; 2026-11-13 "when the fitting is done regenerate all plots and the movies"): the complete set (3.4 set 1,
  set 4, pressures and the version-2 maps, source areas, cliffs) in one folder; captions computed from this fit's data.
* 2026-11-13 (movies; AskUserQuestion answers): "colored should be only that part of the nucleus that can be observed
  (separately for NG and RG -- use appropriate); use actual aperture of RG": RG observable part = the IK rectangle ±45° x
  ±45° (getfov -226332), NG = the IK rectangle ±77.5° x ±77.5° (getfov -226331), tested as Geometry.cpp InsideFOV; inside:
  tint by the fit's acceptance weight as before (colour bars 0..1 / 0..2 kept), outside untinted sunlit grey, boundary a thin
  line; captions state the weight range inside each FOV (RG 0.577..1, NG 1.155..2), that the fit uses the RG hemisphere and the
  NG srccg_shield acceptance, and the share of the modelled signal from outside the IK FOVs; movie 3 with the fit's table;
  the old movies stay.
* 2026-11-13 (movies; AskUserQuestion answers): spacecraft shadow = bus + the two solar-array wings only; the last
  solar-array orientation held after the SA C-kernel ends (10:39:19.851 UTC); the Sun as a disc (lit fraction = unblocked
  fraction of the disc, sampled); movies only (not the fit, no impact analysis); caption sentence that the fitted model has
  no spacecraft shadow.
* 2026-11-13 (movies; user "agree"): the offset between the spacecraft structure origin (-226000, used by
  ROS_STRUCT_V8.BSP) and the trajectory point ROSETTA (-226) is taken as 0 (no loaded kernel gives it; shadow position
  uncertain by ~1 m; `cops_movies.py --origin-offset` stays at its default 0).
* 2026-11-13 (user "yes, keep per-bin colours and regenerate everything after c9"): the cliff stacked figures with one
  colour per slope bin from now on (`cliff_importance_figures.py --groups all`, job variable `CLIFF_GROUPS=all`; the tool
  default `default` stays only to reproduce the figures of 2026-11-12); the bounded softplus chain extended to 9 cycles and
  every figure and movie regenerated for its c9 state into new folders (the c5 folders stay).
* 2026-11-13 (movies; user "modify movies script: freeze the location of the camera at about 100 m altitude in the frame
  moving with the nucleus ..."; AskUserQuestion answers): an EXTRA movie set `movies_softplus_smm1_9cycles_fixedcam/` (the
  gauge-view c9 set stays): camera = the COPS position (ROSETTA + GaugeOffsetSC) at the time the COPS altitude is 100 m
  (time interpolated in log altitude as the frame plan), frozen in 67P/C-G_CK; axis = the RG boresight of that time, image
  up = the NG boresight of that time projected (frozen); perspective, 90° horizontal field of view, no fisheye, no
  re-framing; the whole descent with the same timeline (1800 frames uniform in log altitude, 60 s, 30 fps); everything
  else as the c9 gauge-view movies (IK FOVs from the TRUE gauge position and attitude of each frame, lighting and the
  spacecraft shadow of the frame time); the moving sub-spacecraft point (+) and the gauge position (when in the view)
  marked; the FOV outline where its boundary crosses the visible surface.
* 2026-11-13 (movies, fixed camera; user via the lead): "the part each gauge can see" includes the line of sight: a surface
  point is tinted for a gauge only if (1) its direction from the true gauge position is inside the gauge's IK FOV, (2) its
  facet faces the gauge (n·(gauge − point) > 0), (3) the segment point → gauge crosses no facet (a depth map from the gauge
  position plus exact segment tests, validated against exact rays on samples); everything else untinted sunlit relief.
* 2026-11-13 (movies; user "also, mark on the nucleus surface the location of the touchdown", then "make the touchdown
  location a Formula 1 waving finish flag in the movies"): both new c9 sets (gauge view and fixed camera), not the old
  sets: the touchdown point of the ground-track figures (Rosetta impact 10:39:33.8 UTC, straight-line extrapolation of the
  last two observation rows, `cops_cuda/figures.impact_point` as `tools/pub_figures.py groundtrack`, on the 12M facet hit)
  drawn as a generic chequered finish flag (no logo or trademark; 6 × 4 squares, dark outline) on a thin pole planted at
  the projected point (pole along the projected surface normal when it is ≥ 20° from the line of sight, else image up),
  waving with a travelling sine wave whose phase is frame / 30 s in all three movies, fixed size in output pixels,
  anti-aliased, shown only in frames where the point is visible from the camera (depth-tested at the base point); legend
  "touchdown (impact) 10:39:33.8 UTC" with a small static flag icon; check against the SPICE trajectory (RORB) crossing
  the 12M surface.

* 2026-11-13 (movies; user "redo the movie: 1) is the field of view of RG a square? I thought it is a circle -- verify that
  and if it is a circle then adjust the rendering; 2) during the descent -- detach the camera from the spacecraft at the
  altitude of 200 m and freeze its location and orientation, and then render the descending spacecraft, its shadows and the
  field of view of each gauge ..."; AskUserQuestion answers, relayed by the lead): a NEW set
  `movies_softplus_smm1_9cycles_detach200/`: circular fields of view RG 45°, NG 77.5° half-angle around the boresights
  (the RG entrance is a round 6 mm opening into a 60 mm spherical cavity, Pestoni et al. 2021; the IK rectangle is only
  SPICE's encoding), same inclusion rule as the fixed camera (cone from the true gauge position, facing, line of sight),
  tint = the fit's weight; camera as the gauge-view movies until COPS 200 m, then detached and frozen (position and axis
  in 67P/C-G_CK) with a fixed 90° perspective (the switch labelled); the spacecraft as a 3D model (bus box, both
  solar-array wings, the 2.2 m HGA disc, Sun-lit Lambert; no self- or nucleus shadow on it), its shadow on the surface
  (bus, wings, HGA), the gauges' fields of view as the tinted footprint AND translucent 3D cone lines from the gauge
  positions; COPS booms NOT drawn (no documented geometry; the fit's gauge point marked "COPS (fit)"); HGA orientation held
  after its C-kernel ends (10:39:19.851 UTC), as the solar arrays; frames: the 1800 log-altitude frames to 5.5 m, then
  uniform in time at the last step to the first contact of the model with the 12M surface (computed from the trajectory
  and the model; no penetration shown); the time series ends at 10:39:26 UTC (said on the frames); the touchdown flag at
  the ground-track impact point as before; fov_share recomputed for the circles.
* 2026-11-13 (movies; user ~19:25 CST "draw the flag on top so it is always visible -- make Rosetta touchdown near the
  flag; detach the camera at the altitude of 300 m"; lead's specification): a NEW set `movies_softplus_smm1_9cycles_detach300/`
  (the detach-200 m set stays): the flag drawn on top in every frame where its point projects into the image, also in the
  attached part, legend without "not visible", captions say it is drawn on top; the flag stays at the ground-track impact
  point (centre-of-mass crossing 10:39:33.874) and the movie still ends at the first contact (10:39:30.624, -Y wing),
  because at the last frame the sub-bus point (bus centre along the local normal) is 0.16 m from the flag (limit given:
  ~5 m); camera detached at COPS 300 m with the same rules (frozen in 67P/C-G_CK, fixed 90° perspective, labelled, near
  handling as before); everything else as the detach-200 m set.
* 2026-11-13 (movies; user ~20:05 CST "after detaching the camera from the spacecraft, slowly move the camera and change
  its orientation such that in the last frame the camera [shows] the fields of view of both gauges as well as a panoramic
  view on the mountains and the other comet's lobe that contributes most to the NG data; remove the white cross"; lead's
  specification; user choices ~20:55 CST): a NEW set `movies_softplus_smm1_9cycles_detach300_pan/` (the detach-300 m set
  stays); the white cross (sub-spacecraft "+") removed in every frame (`compose --no-sub-mark`); "the other lobe" = the big
  lobe (body), from the c9 NG source analysis (big lobe 76 % of the all-observation NG surface signal, 46 % of the late;
  FINDINGS.md A10); "the mountains" = the rise 77-82 m east of the impact point (37-39 m high, 45 % of the late NG signal)
  with the big lobe beyond the head's edge; the camera moves smoothly (ease-in/ease-out, arc around the touchdown point,
  orientation slerp, field of view 90° -> 105.5°) from the detach pose to the final pose B chosen by the user from three
  candidates (A west 300 m up, B north-north-west 220 m up, C south-south-west 300 m up; stills in
  `stills_candidates/`); the final pose is reached at frame 1800 (10:39:26, last COPS point) and held during the touchdown
  (frames 1800-1958); the RG field of view is accepted to be invisible in the last frames (a ~2.8 m disc under the
  spacecraft at contact; visible earlier during the motion; said in the captions); flag on top, cones, spacecraft model and
  timeline as the detach-300 m set.

## 6. Change log of this file

* 2026-11-12: created (sections 1–6; 3.5 pending).
* 2026-11-12: 3.5 contribution data (P% source area, source-area occurrence, slope-bin fractions; job 124549).
* 2026-11-12: 3.5 version-2 figures (activity, total rate, source areas, cliffs; jobs 124565-124569), commands and variables.
* 2026-11-12: 3.5 options --low-cut-area (99.9 % rule) and --template-relief (relief under the rate colours); RG gap fill in the cliff figure.
* 2026-11-12: 3.5 --quantity h2o (H2O rate only), activity cut 0.02, preview of the rate with cut 0.02.
* 2026-11-12: 3.5 final: LOWCUT 0.02 default for every quantity; current figures jobs 124589 (activity), 124594/124595 (H2O rate), 124576 (cliffs), 124565 (source areas).
* 2026-11-12: 3.7 movies (jobs 124615-124634).
* 2026-11-13: bounded softplus refit (multipliers >= 1, chain 124664): 3.1 wrapper `figures_12m_nonneg5_smm1.sh`
  (`model_vs_observations_softplus_smm1/`), 3.2 `publication_softplus_smm1/`, 3.4/3.5 complete set in
  `full_set_softplus_smm1_5cycles/` (jobs 124677-124683, 124722, 124727), data-driven captions (`--model-name`, bounds, CO2,
  contribution job), `pub_maps_figures.sbatch` inputs from the environment, `model_data_captions.py`,
  `source_area_compare.py`; 3.7 observable-part movies with the spacecraft shadow (`movies_softplus_smm1_5cycles/`, jobs
  124708-124735), `fov_share.py`, the dark-surface diagnosis of the first movies (`movies_softplus_5cycles/validation/
  artifact_check/`, jobs 124710, 124723); section 5 decisions of 2026-11-13.
* 2026-11-13 (c9 of the bounded chain, chain job 124772): 3.1 `figures_12m_nonneg9_smm1.sh`, 3.2 the c9 publication figure,
  3.3 `MAPS_NAME` / `CONTRIB_NAME` options (c9 maps 124775, contributions 124776), 3.4/3.5 `full_set_softplus_smm1_9cycles/`
  (cliffs with one colour per slope bin, section 5), 3.7 fixed camera (`--camera fixed`), gauge line of sight,
  `validate_fixed`, `touchdown`, `compose --touchdown` (waving finish flag), job variables `MOVIES_PY`, `FROM_WORK`,
  `CAMERA`, `TOUCHDOWN`; the c9 movie sets `movies_softplus_smm1_9cycles/` and `movies_softplus_smm1_9cycles_fixedcam/`;
  section 5 decisions (fixed camera, line of sight, touchdown flag).
* 2026-11-14 (UTC): 3.7 the detach-200 m set (circular FOVs, detach camera, spacecraft model with HGA, FOV cones, frames
  to the first contact; `contact`, `prep --extend-to-utc`, `sc --hga-ck`, `plan --camera detach`, `hgageom`,
  `fov_share.py --cone-*`), section 5 decisions of 2026-11-13 (detach set).
* 2026-11-14 (UTC): 3.7 the detach-300 m set and `compose --flag-on-top` (job variable FLAG_ON_TOP), the copy step also
  copies contact.json; section 5 decisions (flag on top, detach at 300 m).
* 2026-11-14 (UTC): 3.7 the panorama set (`plan --pan-to`, `--pan-start-frame`, `--pan-end-frame`, `compose --no-sub-mark`;
  job variables PAN_TO, PAN_START, PAN_END, NO_SUB_MARK; `ng_source_map.py`, `ng_source_regions.py`,
  `ng_lobe_share_vs_altitude.py`, `pan_candidates.py`, `cops_movies_pan_stills.sbatch`); section 5 decisions of 2026-11-13
  (panorama, white cross removed, pose B, hold from frame 1800, RG FOV invisible at the end accepted).
