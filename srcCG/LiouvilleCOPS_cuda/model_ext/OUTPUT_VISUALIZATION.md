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

**H2O night-side fix (2026-11-18; user decision "(b) -- use 0.02; x10 was not intended"):** every map, production-rate
figure, contribution, cliff and illumination product made BEFORE 2026-11-18 used H2O q_min = 0.2 q_max (the law of all
fits so far); the fixed law is q_min = 0.02 q_max (as srcCG). `model_ext/nonneg_maps.py` now holds the ratio as
`H2O_QMIN_RATIO` (0.02; snapshots `model_ext_frozen/nonneg_maps_20261118T072729Z`, `nonneg_contrib_20261118T072729Z`).
Make maps of fits done with the new law with these snapshots. The older snapshots refuse with the fixed
`fit_cops.py` (law_guard: "the reference source law changed"); maps of the old (0.2) fits can only be remade with the
pre-fix fit_cops (git history, commit 4b1e0c2 or earlier) and their own snapshots.
**Maps on base_q002 (2026-11-18 13:30 CST):** base_q002 has no run_manifest.json by design, so the geometry validation
needs the EXPLICIT shape manifest: `sbatch --export=ALL,CHAINS=<...>,BASE=$R/base_q002,SHAPE_MANIFEST=$R/base_q002/cpp_run_manifest.json,SHAPE_REFERENCE=$R/base/run_manifest.json,MODEL_EXT=model_ext_frozen/nonneg_maps_20261118T192126Z jobs/nonneg_maps.sbatch`
(`validate-geometry --shape-manifest ... --shape-manifest-must-equal ...`); a validated geometry exists:
`$R/maps_geometry/job_136691/geometry` (use it as GEOM for maps / contributions / the localnn and d2 tools, with
`--contrib-code model_ext_frozen/nonneg_contrib_20261118T192126Z`). D2 on base_q002: `jobs/d2_tools.sbatch` with
`D2DIR=$R/clip_refit/D2_q002_warm D2_BASE=$R/base_q002` (job 136692: `D2_q002_warm/maps`).
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

**Half-space set (2026-11-14 UTC; user decisions in section 5):** `figures_cluster/12M/movies_softplus_smm1_9cycles_detach300_pan_hemi/`
(work `<chain>/movies_work_c9_detach300_pan_hemi`, products copied from the panorama work; jobs 124942 copy/scrows/plan/
validate_fixed/touchdown, 124943 render, 124944 `fov_share.py --hemi-sc` (replaces fov_share for this set), encode and qc;
regression 124945, tests 124940). As the panorama set (pose B, pan to frame 1800 and hold, flag on top, no "+"), but each
gauge's field of view = the fitted model's acceptance region: the front half-space of its boresight (`--mask-rg hemisphere
--mask-ng hemisphere`; RG tint -l·b_RG 0–1, NG tint the front factor 1 - l·b_NG 1–2), minus the directions in which the
spacecraft model (bus box, ten solar-array panels, HGA disc) hides the surface from the deck point (`--sc-block-from
rosina-deck`: SPICE ROS_ROSINA (0.8928, -0.7447, 2.6058) m in ROS_SPACECRAFT moved onto the bus top face z = 2.655 m + 1 mm;
exact segment test `sc_segments_blocked`, checked against render.cu ray_first_hit on the drawn triangles: 0 mismatches);
direction, weight, facing and line of sight from the fit's gauge point as before; line of sight for a half-space by cube
faces (`gauge_los_cube`); no 3D cone lines (outline on the surface only); NG's shield part behind is not shown (captions).
`scrows` writes the spacecraft geometry at the observation rows for `fov_share.py --hemi-sc` (share not coloured = NG back
part + spacecraft-blocked front part; RG spacecraft-blocked part; blocked solid-angle fraction of each half-space).
```bash
STEPS=copy,scrows,plan,validate_fixed,touchdown HGA_CK=<ROS_HGA_2016_V0043.BC> HGA_HOLD_UTC=... SA_HOLD_UTC=... CAMERA=detach DETACH_ALT=300 \
  FOV_DEG=90 PAN_TO=<pan_B.json> PAN_END=1800 MASK_RG=hemisphere MASK_NG=hemisphere SC_BLOCK_FROM=rosina-deck \
  FROM_WORK=$CHAIN/movies_work_c9_detach300_pan sbatch --nodelist=gpu002 --mem=100G --export=ALL cops_movies_prep.sbatch
MASK_RG=hemisphere MASK_NG=hemisphere SC_SHADOW=1 SC_BLOCK_FROM=rosina-deck NGPU=8 ... cops_movies_render.sbatch
NCYC=9 FOV=$WORK/fov.npz OUTP=<folder>/validation/hemi_share/hemi_share_c9 HEMI_SC=$WORK/sc_rows.npz HPREP=$WORK/prep.npz \
  HMOVIES=<cops_movies.py> sbatch --export=ALL fov_share.sbatch
FLAG_ON_TOP=1 NO_SUB_MARK=1 ... cops_movies_encode.sbatch
```

**Contribution set (2026-11-14 UTC; user decisions in section 5):** `figures_cluster/12M/movies_softplus_smm1_9cycles_detach300_pan_contrib/`
(work `<chain>/movies_work_c9_detach300_pan_contrib`; frames base `runs/shap8_012m_lategrid_dz1.5/frames_base_pan1959` (frames
0-1808, job 124961, stopped at the end of the RG temperature table) + `frames_base_pan1959_tail` (frames 1808-1958, job 124969,
TemperatureRGMode constant 253.438 K); contributions `<chain>/frame_contrib_c9_pan1959` (job 124972); jobs 124966 prep, 124975
render, 124976 compose/encode, 124977 qc). As the half-space set's camera, spacecraft and flag, but the colour is the fitted
model's contribution per unit surface area at the frame time, exactly as the fit computes it: the fitting base's executable
and boundary.in run at the 1959 frame times (`frames_base.sbatch`), per record the frozen nonneg_contrib arithmetic with the c9
state (`frame_contrib.py`): NG density F (√π β/2) K_n / A [m⁻⁵] (the fit's NG weight incl. the shield side inside K_n), RG
flux F (K_Φ + V·b √π β/2 K_V) / A [m⁻⁴ s⁻¹]; no spacecraft blocking; facets without a kernel record uncoloured. Fixed log
range per gauge (`contrib_range.py`): bottom = the per-area level above which the facets carry 90 % of the surface signal at
frame 0, top = the 99.99th percentile of the per-area values of frames 1790-1799 (the last 10 frames with COPS data; the
literal rule on frames 1949-1958 fails for RG: the model's RG flux falls to 1.4e-18 m⁻² s⁻¹ at contact). No field-of-view
outlines, mask boundary lines or cone lines; log colour bars in physical units. `compose --hold-seconds 15`: 450 frames
holding the last frame with the flag waving ("touchdown — final frame held"). Verification (`verification.json`): the
decomposition equals record_contributions to 1.6e-15; with the base's own source the sums equal the frames base's
COPS_model.dat n / Φ columns to 4.9e-11 (print precision); frame 1799 = base row 1549 record for record.
```bash
FB=$RUN/frames_base_pan1959 OBSF=<frame-time observation file> sbatch --nodelist=gpu002 --export=ALL frames_base.sbatch
FB=$RUN/frames_base_pan1959_tail OBSF=<same> FIRST=1808 TRG_CONST=253.438 sbatch --nodelist=gpu002 --export=ALL frames_base.sbatch
CHAIN=... FB=... FB_TAIL=... OUTC=$CHAIN/frame_contrib_c9_pan1959 sbatch --export=ALL frame_contrib.sbatch
contrib_range.py $CHAIN/frame_contrib_c9_pan1959 $CHAIN/frame_contrib_c9_pan1959/color_range_last_obs.json
STEPS=copy,plan,touchdown CAMERA=detach DETACH_ALT=300 FOV_DEG=90 PAN_TO=<pan_B.json> PAN_END=1800 FROM_WORK=<half-space work> ... cops_movies_prep.sbatch
SC_SHADOW=1 COLOR_CONTRIB=<frame_contrib dir> COLOR_RANGE=<color_range_last_obs.json> NGPU=8 ... cops_movies_render.sbatch
FLAG_ON_TOP=1 NO_SUB_MARK=1 HOLD_SECONDS=15 ... cops_movies_encode.sbatch
```

### 3.8 Contribution density maps and region shares (2026-11-14)

In the full set of a state (`figures_cluster/12M/full_set_softplus_smm1_9cycles/`, first made for c9):
* `contribution_density/12M_contrib_density_<HHMMSS|pre_gap|gap|late>_{NG_RG,NG,RG}`: lon/lat maps (global + a close-up fitted
  to the maps' 90 % source areas near the impact) of the contribution per unit surface area exactly as the fit computes it
  (NG density F (√π β/2) K_n / A [m⁻³ per m²], the shield side included; RG flux F (K_Φ + V·b √π β/2 K_V) / A [m⁻² s⁻¹ per
  m²]; no spacecraft blocking): snapshots at base rows (10:33:10, 10:37:31.282 = 100 m, 10:39:09.282 = 20 m, 10:39:26 =
  5.5 m; instantaneous, no time response) and phase time means (NG and RG before the NG gap, RG in the gap, the late window
  as W10x5 of the grid rows). Log colour per map: top = 99.99th percentile of the values > 0, bottom = the level above which
  99 % of the map's signal lies or 3 decades below the top if lower; white contours = the 50 % (solid) / 90 % (dashed)
  source areas of the map (smoothed pixel mask). Captions give the share of each map's signal from the facets drawn in the
  lon/lat maps (overhang facets are not drawn).
* `regions/12M_region_share_{lobe,distance}_{NG,RG,NG_RG}`: stacked share of the modelled surface pressure by lobe (neck
  plane of FINDINGS A10) and by distance shell from the touchdown point (0-30, 30-100, 100-300, 300-1000, 1000-3000, > 3000 m)
  along the descent, in the style of the cliff figures (NG gap hatched in NG, dashed in RG, late window boxed).
```bash
CHAIN=<chain> OUTD=<chain>/contrib_phase_c<N> sbatch --export=ALL contrib_phase_data.sbatch                # CPU
DATA=<chain>/contrib_phase_c<N> OUTF=<full set> sbatch --nodelist=gpu002 --export=ALL pub_contrib_density.sbatch   # 1 GPU, ~5 min
regions_figures.py --data <chain>/contrib_phase_c<N> --base <fit base> --bundle <bundle> --out-dir <full set>/regions   # CPU
```

### 3.9 Stage-1 non-negative activity options vs the data and the softplus c9 (2026-11-14)

`figures_cluster/12M/model_vs_observations_nonneg_options_stage1/` (the four options of FINDINGS / model_ext README
"Alternative non-negative activity models": A sum of squares rank 3, B background bg 1e-3/1e-2/1e-1, C local non-negative
basis f 0.1/1/10, D log form with spread penalty f 0.1/1/10; single fits from D2):
* `12M_options_<A|B|C|D>_{time,altitude,altitude_NG,altitude_RG,altitude_NG_RG}`: every setting of one option (dashed /
  dash-dot / dotted) vs the COPS data and the bounded softplus c9 (solid), model_obs_figures style.
* `12M_options_best_*`: the lowest-data-cost setting of each option + c9 + the clipped D2 (dotted).
* `summary.txt` / `summary.csv` / `best.json`: per fit data cost, penalty, nfev / termination, rms ln NG/RG (all, before the
  10:33:10 gap, late), multipliers, parameters at bounds, Q_H2O at 10:33:10, active seen-area fraction (A above 1 % of the
  srcCG median over the seen facets), seen area holding 50/90 % of the production, NG/RG shares of the active area; read
  from each run's `summary.txt.json` (A, B, D) or `analysis/report.json` (C), nothing recomputed.
* `activity_maps/12M_options_activity_best`: A of the best setting of each option and of c9 (rows), global + impact
  close-up, continuous log colour above the 1 % threshold on one scale, light green seen below it, grey unseen.
* Only runs whose summary / report and `fit/final/COPS_model.dat` exist are included; rerun into a new folder when the
  others finish (the wrapper refuses an existing OUT).
* v2 (`..._stage1_v2/`, 2026-11-14 evening): the same fits with `DIAG=figures_cluster/scripts/data/stage1_convergence_20261114.json`
  (summary convergence column + section: verdict of the CPU diagnosis checks/stage1_diag_20261114 with its scaled gradient,
  Gauss-Newton predicted decrease and best trial step; also in the captions), `RESID=1`
  (`12M_options_best_residuals`: ln(model/obs) vs altitude, (a) NG (b) RG, best of each option + c9 + D2, every fit solid in
  its colour of the best figure, the fits' own residuals incl. W10x5 late, rms recomputed and checked against the summary)
  and `ACT_FROM=<first folder>` (activity maps copied after checking that the best settings, roots and inputs are
  unchanged; COPIED_FROM.txt). Option-C captions state the GPU cross-check from gpu_eval/gpu_check.json.
* final (`..._stage1_final/`, 2026-11-15): `RUNS=figures_cluster/scripts/data/stage1_final_runs_20261115.json` (the latest
  completed fit per setting: B/A continuations stage1c_*, D refits with the fixed-level penalty, C f 0.1/1/10 and the
  multiplier-penalty variant f 10 P_m 0.01), `TEXTS=.../stage1_final_texts_20261115.json` (caption texts), `CONV_TERM=1`
  (convergence from each fit's termination: ftol / not converged), `RESID=1`; summary.txt/.csv carry the stage-1 costs each
  continuation started from (prev_*) and the C multiplier penalty (penalty_mult); activity maps by a new GPU job
  (`options_activity_maps.sbatch`, job 131968), since the best C, A, B, D settings changed.
```bash
cd $COPS_OUTPUT_ROOT
# v2 variables: DIAG=figures_cluster/scripts/data/stage1_convergence_20261114.json RESID=1 ACT_FROM=<earlier folder>
# final variables: RUNS=figures_cluster/scripts/data/stage1_final_runs_20261115.json TEXTS=figures_cluster/scripts/data/stage1_final_texts_20261115.json CONV_TERM=1 RESID=1
OUT=<folder> srun -p normal128cores -A general -N1 -n1 -c 8 --mem=64G -t 00:50:00 bash figures_cluster/scripts/figures_12m_nonneg_options_stage1.sh
sbatch --export=ALL,OPT_DIR=<folder> jobs/options_activity_maps.sbatch          # 1 GPU, short
.venv/bin/python figures_cluster/scripts/options_howto.py --folder <folder> --wrapper-log <log> --maps-job <job>
.venv/bin/python figures_cluster/scripts/collect_captions.py <folder> --title "..." --order ""   # CAPTIONS.md incl. the maps
```

### 3.10 Final stage-2 comparison (2026-11-17)

`figures_cluster/12M/model_vs_observations_nonneg_options_stage2_final/` (`stage2_figures.py`, CPU, about 1 min; spec
`figures_cluster/scripts/data/stage2_final_spec_20261117.json`): the final states (last cycle's fine tuning, c5/ft/final) of the
regularized stage-2 chains (stage2r_*: fine tuning with P_delta, f_delta 10, cap ln 100, templates on) and, as context, the
unregularized ones (stage2_*; softplus chain_5cycles_smm1 after 5 and 9 cycles) and D2. Figures (model_obs_figures.py style,
late window shaded and labelled): `12M_stage2_reg_best_*` (softplus reference, D2, A, B, best C, D), `12M_stage2_C_reg_*`,
`12M_stage2_C_unreg_*`, `12M_stage2_<A|B|D>_*` (regularized vs unregularized), `12M_stage2_softplus_*`; residuals
`12M_stage2_{reg,unreg}_best_residuals`; `summary.txt/.csv` (final data cost, P_delta, fine-tuning and last-fit terminations
and nfev (segments), fit penalty, rms all / late, multipliers, parameters at bounds, data cost after each cycle; pending
chains). A chain enters when its c5 fine tuning is complete (C: also c5/analysis/report.json). No maps; the script stops
below 25 GB free. To add chains: rerun into a new folder (or move the folder aside).
State 2026-11-17 14:45 CST: all nine regularized chains (incl. C f 10 / f 1 with P_m 1e-4, cycles 1-3 and part of 4 on
CPU, the rest on the GPU) and the unregularized context; `kappa_check/kappa_table.txt` (`kappa_table.py`): kappa of
every fit (first and last fit of each chain, stage 1, softplus chains, D2) with the at-bound test.
```bash
srun -p normal128cores -A general -N1 -n1 -c 8 --mem=64G -t 00:50:00 bash -c ".venv/bin/python -W ignore figures_cluster/scripts/stage2_figures.py --spec figures_cluster/scripts/data/stage2_final_spec_20261117.json --out <new folder>"
.venv/bin/python figures_cluster/scripts/stage2_howto.py --folder <folder> --spec <spec> --log <build log> [--note ...]
.venv/bin/python figures_cluster/scripts/kappa_table.py --spec figures_cluster/scripts/data/stage2_final_spec_20261117.json --out-dir <folder>/kappa_check
```

### 3.11 Full figure sets of the stage-2 states and D2 (2026-11-17)

`figures_cluster/12M/full_set_stage2_{softplus_reg,A_reg,B_reg,D_reg,C_f1_pm1em4_reg}/` and `full_set_D2/`: the six cases of
12M_stage2_reg_best, each at its final state, with every folder of `full_set_softplus_smm1_9cycles/` except the c9-specific
comparison_vs_* (possible comparisons listed in each README) and no movies: visibility/ (set 1, fit-independent, PNGs
pixel-identical to the c9 set's; ground track with the Rosetta impact only), model_data/ (incl. the publication NG/RG figure
of `publication_stage2_<key>/` / `publication_D2/`), activity/ (LOWCUT 0.02, flat), production_rate/ (H2O, 10:33:10 and
10:39:26 UTC, LOWCUT 0.02, relief 0.7), source_areas/ (P50-90 of the 95 post-gap observations), cliffs/ (one colour per slope
bin), contribution_density/, regions/, checks/ (incl. facts.json / facts.txt), CAPTIONS.md, README.txt (Facts computed by
`full_set_facts.py`), HOW-TO-RERUN.md, code/, provenance. Same tools, styles and user decisions as the c9 set.
State products next to each state: maps_c5 / contributions_c5 / contrib_phase_c5 / v2_total_crosscheck_c5 (chains), C
`c5/{maps_ft,contributions_ft,contrib_phase_ft,v2_total_crosscheck_ft,tools_state_ft}` (agent C's localnn_tools), D2
`clip_refit/D2/{tools_state,maps,contributions,contrib_phase,v2_total_crosscheck,tools_state_verify.json}` (d2_tools.py).
Code (2026-11-17, lead approval; defaults unchanged, tested): A / B / D maps and contributions with the new frozen
`nonneg_maps_20261117T224340Z` / `nonneg_contrib_20261117T224340Z` (state_activity re-derives the alternative models bitwise);
C and D2 through adapters (localnn_tools.py, d2_tools.py) and the figure-script fallbacks (bin edges from the bin names,
multiplier bounds from the export's manifest, contrib_phase_data --state-adapter, model_data captions without fine tuning).
```bash
cd /data001/highresolu_vt/cops_output/jobs
# maps + contributions of the state (nonneg chains; C: localnn_tools.sbatch STEPS=maps; D2: d2_tools.sbatch), then
KEY=<key> CHAIN=<chain> MAPS=.. GEOM=.. CONTRIB=.. PUB=<publication folder> PUBSTEM=.. PUBTABLE=.. MODEL_NAME=.. MODEL_LABEL=.. \
  MODEL_SHORT=.. TITLE_PRESSURES=.. [FT= FIT= XC= PHASE= MX= CX= ADAPTER= MAPS_JOB= CONTRIB_JOB=] bash full_set_submit.sh
# CPU afterwards: visibility captions, regions/, model_data captions, facts, CAPTIONS.md
KEY=<key> FT=.. FIT=.. MAPS=.. CONTRIB=.. PHASE=.. [CHAIN= NCYC=] JOB_SET4=<set1/4 job> MODEL_NAME=.. [CAVEATS=..] bash finish_full_set.sh
.venv/bin/python figures_cluster/scripts/full_set_docs.py --spec figures_cluster/scripts/data/full_sets_stage2_20261117.json --key <key>
```

### 3.12 Word document comparing the six full sets plot by plot (2026-11-17)

`figures_cluster/12M/full_sets_stage2_comparison/12M_stage2_full_sets_comparison.docx` (user: "prepare word doc that compares
all the cases: group the same plot for all cases and provide summary and that for every plot"); `fullsets_comparison_docx.py`
(version-controlled copy `figures/`). For each of the 85 PNGs of a full set (3.11) one figure group: the six cases in one grid
(order softplus, A, B, D, C, D2; portrait or landscape, 3 x 2 or 2 x 3, chosen for the largest images), each labelled with its
case, caption "Figure N: <plot> — all six cases", then a Summary (what the plot shows from its caption; a table of the case
numbers; how the cases compare; observations from the images, marked as such). Plots whose six PNGs are byte-identical
(visibility/, the geometry panels (c)-(f)) are shown once. Front matter: cases table (data cost, terminations, rms ln all /
late, Q_H2O total / seen, kappa, multipliers, README caveats) and an overall summary (FINDINGS A12 / A14 / A15); Appendix A:
the source PNG paths (SHA256 prefixes). Every printed number is read from a file and listed in `numbers_used.json` (literal
numbers of hand-written sentences checked against the captions / READMEs / FINDINGS / this file; the build refuses to write
the document otherwise). Images: downscaled copies in `images/` (maps JPEG q78, line plots 256-colour PNG, 1000-1400 px).
Interpreter: separate venv `$COPS_OUTPUT_ROOT/tools/docx_venv` (python-docx, Pillow, defusedxml; not the project `.venv`).
No PDF: LibreOffice is not installed on Voyager (login and compute nodes); `validate` checks counts, page fit and zip
entries and draws layout previews; the docx skill's `validate.py` checks the OOXML schema. CPU job, ~20 s; refuses an existing
document (move the folder aside as `.superseded_<UTC>`).
```bash
OUT=/data001/highresolu_vt/cops_output/figures_cluster/12M/full_sets_stage2_comparison; PYD=/data001/highresolu_vt/cops_output/tools/docx_venv/bin/python
srun -p normal128cores -c 8 --mem=32G -t 00:25:00 bash -c "$PYD -I figures_cluster/scripts/fullsets_comparison_docx.py build --out $OUT --jobs 8 && $PYD -I figures_cluster/scripts/fullsets_comparison_docx.py validate --out $OUT"
```
To compare other sets: edit `CASES` (set folder, file stem, summary.csv key, residual table names) and the catalogue in the
script (it asserts 85 plot types).

### 3.13 Slope at equal illumination: does the source rate depend on the slope angle? (2026-11-17)

`figures_cluster/12M/cliff_effect_stage2/` (user: "test several definitions: >40, >60, >80, as well as 0-40,40-50,50-60,
60-70,70-80, 80-90; for the six cases individually ... find whether at the same illuminations condition the source rate
depends on the slope angle ... word doc ... pack the same plot for all six cases in a cluster and add a summary").
Scripts (version-controlled copies `figures/`, `jobs/`): `cliff_effect_data.py` (one CPU job per case,
`jobs/cliff_effect_data.sbatch`, 3-5 min), `cliff_effect_figures.py` (`jobs/cliff_effect_figures.sbatch`, all cases, ~3 min;
labels and caveats from `figures/data/cliff_effect_cases_20261117.json`), `cliff_effect_docx.py` (docx_venv; imports
`fullsets_comparison_docx.py`). Per case `<key>/`: P1_NG / P1_RG (F, A, ratios to 0-40, active fraction, lit / shadowed /
night strip, per illumination class, gauge-seen surface, four snapshots), P1_whole (whole surface, maps' .vtu, 10:33:10 /
10:39:26), P2_NG / P2_RG (g = K/dS = 1/r^2 x cos theta x W x visible fraction at the centroid), P3_NG / P3_RG (matched cells
SZA class x 10 % distance x 0.1 cos theta: contribution-weighted distribution of the ratio per m^2), P4_NG / P4_RG (signal
share / area share = activity x illumination x geometry x covariance, snapshots and phase means); answer.csv / .json (the
SZA-controlled ratio, lead's definition; additions: under-sampled-reference flag (grey columns, dagger), variant over
well-sampled classes, ratio of sums), p2_summary.json, p3.csv, p4.csv, P1_whole_tables.csv. Slope categories = the fit's
bins and their unions (>40, >60, >80 include 90-180; 90-180 alone as an addition). Document
`12M_stage2_cliff_effect.docx`: question, model facts, definitions, checks, overall answer per case (tables, caveats), nine
clusters (the same plot of the six cases, labelled (a)-(f)) each followed by a summary with numbers and the answer for that
plot; `numbers_used.json`; validation (python-docx reopen, page fit, layout previews, OOXML schema). Checks in each
`data/<key>/cliff_data.json` (contributions = the fit's bitwise; F vs the maps <= 4e-8; row totals = the full sets'
snapshots). Rerun: README.txt and HOW-TO-RERUN.md of the folder.

### 3.14 Illumination maps at the production-rate map times (2026-11-18)

`figures_cluster/12M/illumination_maps_103310_103926/` (user: "check the production rate plots -- i dont see the effect of the
shadow or night side decrease of the productions -- is it there? generate illumination maps for the time intervals that
corresponds to that on the production rate maps"). Script `figures_cluster/scripts/illumination_maps_figures.py` (copy
`figures/`), jobs `jobs/illumination_numbers.sbatch` (CPU, ~2 min) and `jobs/illumination_maps.sbatch` (1 GPU, `TIME=`, ~7 min);
then `illumination_maps_figures.py docs` on the login node. Quantities (fixed classes, the same in every panel):
illumination class (day #E69F00, night #0072B2, cast shadow #CC79A7), cos SZA of the day facets (10 classes of 0.1,
'YlOrBr'; night / shadow in their class colours), q/q_max = 0.2 (night, shadow; colour-bar triangle) or 0.2 + 0.8 cos SZA
(8 classes of 0.1, 'YlOrBr' reversed) and log10 q/q_max (14 classes of 0.05 decade); masks `whole_surface` and `seen_only`
(unseen grey). Layouts = every layout of `full_set_*/production_rate/` (templates with slope and cliff rows, panels a / b,
lon/lat and 3-D ground tracks; relief 0.7 in templates and panels, flat colours in the ground tracks, Rosetta impact only);
`side_by_side/`: production rate (classes checked equal to the full set's manifest) | q/q_max on the same facets | q/q_max
on the whole surface. Illumination is fit-independent (cos_SZA, illumination and seen identical bitwise in the six cases'
maps), so the folder serves every full set. Numbers and the answer: README.txt, checks/illumination_numbers.txt. Rerun:
HOW-TO-RERUN.md of the folder (figure jobs ran the script version sha256 1980967985..., the docs step 829697321f...).
**2026-11-18 (H2O night-side fix):** this folder shows the 0.2 law (q/q_max = 0.2 at night, 0.2 + 0.8 cos SZA by day) of
all fits so far. The script no longer has the literal 0.2: `QRATIO` = the maps' own `nonneg_maps.H2O_QMIN_RATIO`
(`set_qratio`; the numbers json records it as `qratio_h2o`), and the factor classes, ticks and texts derive from it
(identical to the figures above for 0.2). For the fixed law (0.02) they become: f = 0.02 at night / in shadow,
0.02 + 0.98 cos SZA by day, 10 classes of 0.1 from 0 to 1 (triangle 0.02), log10 f 34 classes of 0.05 from -1.70
(ticks every 0.2 / 0.4 decade); a factor-50 (1.70-decade) range. That layout is derived, not a user decision: confirm it
when the maps are redone after a refit. A rerun needs maps made with `nonneg_maps_20261118T072729Z` (an older maps
folder has no `H2O_QMIN_RATIO` and is refused).

### 3.15 Source-rate maps of the runs fitted with the new H2O law (2026-11-18)

`figures_cluster/12M/source_rate_q002/<run>/` (user: "prepare the surface source rate maps for completed runs for inspection
of the source rate"; the H2O night-side fix of 2026-11-18: q_min = 0.02 q_max, srcCG's value; before, every fit and figure
used 0.2). Per completed run (`figures/source_rate_q002_submit.py --run <key>`, GPU job `jobs/source_rate_q002.sbatch`):
activity/ and production_rate/ exactly as the full sets (pub_maps_v2_figures.py: templates + panels, lon/lat and 3-D ground
tracks; seen-only colouring, area cut 0.02, relief 0.7 for the rate; captions state the law and that the old figures used 0.2;
the activity model described in the model's own terms via `--model-text`) and inspection/ (source_rate_q002_inspect.py, per
map time: (a, b) new-law F_H2O, (c, d) the old-law counterpart's F_H2O on the SAME log classes, (e, f) the new q/q_max on the
whole surface; checks/inspect_<hhmmss>.{json,txt}: Q_H2O total and seen, night + cast-shadow shares, active area, new vs old).
README.txt / CAPTIONS.md / HOW-TO-RERUN.md / code/ by source_rate_q002_docs.py. Runs, their maps, old-law counterparts and
caption texts: `figures/data/source_rate_q002_runs_20261118.json` (D2_q002_warm and the nine stage-2 chains of the new law;
C f 10 / f_m 1e-4 has no old-law maps: no inspection figure). `cluster/` (source_rate_q002_cluster.py, CPU): the production-rate
panels of all completed runs, one row per run, per time (each panel keeps its own classes). GPU needed (cops_cuda.figures
rasterizer); gpu003 is not schedulable for us (2026-11-18).

### 3.16 Comparison set of the new-law stage-2 round (2026-11-19)

`figures_cluster/12M/model_vs_observations_nonneg_options_stage2_q002/` (user: "prepare a full set of comparison plots for the
last round of fitting"; the round = the refitted D2 and the nine regularized stage-2 chains with H2O q_min = 0.02 q_max on
base_q002). The set of 3.10 for the new-law states (`stage2_q002_figures.py`, CPU job `jobs/stage2_q002_figures.sbatch`, ~3
min; docs `stage2_q002_docs.py`): 12M_q002_best_* (softplus, D2, A, B, best C, D vs COPS), 12M_q002_C_* (the five C chains),
12M_q002_vs_old_<run>_* (new law solid vs old law dashed for softplus, A, B, D, D2 and each C setting), residuals
(12M_q002_best_residuals, 12M_q002_C_residuals, 12M_q002_vs_old_<run>_residuals for softplus, A, B, D, D2, best C),
12M_q002_convergence_<chain> (data cost after each cycle's fit and fine tuning, new vs old) and 12M_q002_convergence_all;
summary.txt / .csv (new and old rows: costs, P_delta, terminations / nfev, rms all / late, multipliers, bounds, kappa,
Q_H2O total and seen at 10:33:10, old-law cost and change). Every legend label carries "new law q_min 0.02 q_max" or
"old law q_min 0.2 q_max"; every caption states both laws. 78 figures.

### 3.17 Word document clustering the plots of the ten new-law runs (2026-11-19)

`figures_cluster/12M/stage2_q002_comparison/12M_stage2_q002_comparison.docx` (user: "similarly to what you have done
priorly, generate word doc clastering the same plot from each of the runs and summarizing them"); `q002_comparison_docx.py`
(version-controlled copy `figures/`; derived from 3.12's generator, same rules). Inputs: 3.15 (`source_rate_q002/<run>/`) and
3.16. Front matter: purpose, the law correction, a runs table from summary.csv (new / old data cost and change,
terminations, rms all / late, multipliers, kappa, Q_H2O total / seen), a summary of the round (FINDINGS A18) with caveats.
Section 3: the multi-run figures of 3.16 (12M_q002_best_*, 12M_q002_C_*, convergence_all), one image each. Sections 4-7: one
cluster per per-run plot type (new-vs-old model/data figures, residuals, convergence; activity; production rate at 10:33:10 /
10:39:26; inspection), the runs in the fixed order D2, softplus, A, B, D, C f1/1e-4, C f10/1e-4, C f1/0.01, C f10/0.01,
C f10/0.01 off, all on one page (2-5 columns), a placeholder with the README's reason where a run lacks the plot; plot types
byte-identical in all runs shown once. Every figure caption states the law (validate checks it). Summaries: number table per
run (summary.csv, inspection json, captions, manifests, slope_bins json, tables/*.dat), comparison, observations marked as
such; `numbers_used.json`. No PDF (no LibreOffice on Voyager); docx skill `validate.py` for the schema. Notes in the document:
the full-set categories (contribution density, source areas, cliffs, regions) were not produced for the new-law runs; the
tall `source_rate_q002/cluster/` page is not inserted (it repeats the production-rate panel clusters).
```bash
OUT=/data001/highresolu_vt/cops_output/figures_cluster/12M/stage2_q002_comparison; PYD=/data001/highresolu_vt/cops_output/tools/docx_venv/bin/python
srun -p normal128cores -c 8 --mem=32G -t 00:25:00 bash -c "$PYD -I figures_cluster/scripts/q002_comparison_docx.py build --out $OUT --jobs 8 && $PYD -I figures_cluster/scripts/q002_comparison_docx.py validate --out $OUT"
```

### 3.18 Full figure sets of the ten new-law states (2026-11-20)

`figures_cluster/12M/full_set_q002_<key>/` for D2 (refitted), softplus, A, B, D and the five option-C chains (user: "prepare the
full set of plots for the last round of fits with the fixed surface source distribution"): exactly the 3.11 sets (visibility/,
model_data/ with the publication NG/RG figure, activity/, production_rate/, source_areas/, cliffs/, contribution_density/,
regions/, checks/, CAPTIONS.md, README.txt with Facts, HOW-TO-RERUN.md, code/, provenance; no movies) on base_q002 with the
new-law snapshots (nonneg_maps / nonneg_contrib_20261118T192126Z) and each state's maps / contributions; plus inspection/
(new vs old law, 3.15). activity/, production_rate/ and inspection/ are sha256-checked COPIES of `source_rate_q002/<run>/`
(same tool, settings, inputs); everything else was made for the set (`figures/full_set_q002.py`, `jobs/full_set_submit.sh` with
the 2026-11-20 overrides, `jobs/full_set_q002_finish.sbatch`). Every caption carries the law via the model name (the visibility/
captions are the fit-independent copies of 3.4). Model-data reference: the refitted D2 (same law). The publication NG/RG figures:
`12M/publication_q002_<key>/`.

### 3.19 ONE Word document summarizing all fits of the new-law round (2026-11-20)

`figures_cluster/12M/q002_all_fits_summary/12M_q002_all_fits_summary.docx` (user: "prepare the full set of plots for the last
round of fits with the fixed surface source distribution; also prepare a single word doc summarizing all the fit");
`q002_all_fits_summary_docx.py` (version-controlled copy `figures/`; from 3.12 / 3.17's generators, same rules). For this round it
supersedes 3.17's document. Inputs: 3.16 (comparison set), 3.18 (the ten full sets `full_set_q002_<key>/`, discovered by name or
`--sets-spec`), 3.15 (inspection). Structure: runs table new vs old (summary.csv) + the full sets' facts; summary of the round
(A18: fit quality, source rates, multiplier shift new [old], full-set facts incl. the overhang share of NG at the last
observation, caveats); the 13 multi-run figures; 7 per-run new-vs-old clusters; one cluster per full-set plot type (85; ten runs in
the fixed order; byte-identical types shown once; activity/ and production_rate/ of the full sets are copies of 3.15 and appear
once); 2 inspection clusters; Appendix A old-law results, B sources with SHA-256, C build. Every caption states the law
(validate checks it). Observations from the images: data file `figures/data/q002_all_fits_observations_20261120.json`, written
after viewing contact sheets (`q002_all_fits_sheets.py`); a cluster without an entry gets none. No PDF (no LibreOffice on
Voyager); docx skill `validate.py` for the schema.
```bash
P=/data001/highresolu_vt/cops_output/tools/docx_venv/bin/python; S=figures_cluster/scripts/q002_all_fits_summary_docx.py
OBS=figures_cluster/scripts/data/q002_all_fits_observations_20261120.json; OUT=figures_cluster/12M/q002_all_fits_summary
srun -p normal128cores -c 4 --mem=16G -t 00:10:00 $P -I $S check --observations $OBS
srun -p normal128cores -c 8 --mem=32G -t 00:30:00 bash -c "$P -I $S build --out $OUT --observations $OBS --jobs 8 && $P -I $S validate --out $OUT"
```

### 3.6 Other analysis figures and tables

* `figures_cluster/scripts/positive_table.py`, `positive_cycles_table.py`, `finetune_table.py`, `scan_summary.py`:
  tables used by the wrappers.
* `analysis/shadow_zero_activity_12m/shadow_classes.py` (CPU job): illumination of the seen surface with and without
  activity (no figure).

## 4. Redo all figures after a model rerun (recipe)

1. Identify the new chain(s): `<chain>/c<N>/ft/final` complete (`finetune_result.txt`), and the chain job id(s).
2. **Maps** (3.3): `sbatch --export=ALL,CHAINS=<chain>,MODEL_EXT=<frozen nonneg_maps> jobs/nonneg_maps.sbatch`; note the
   job id (its geometry folder `$RUN/maps_geometry/job_<id>/geometry` is an input of 3.4). For fits made after the
   2026-11-18 H2O night-side fix use `MODEL_EXT=model_ext_frozen/nonneg_maps_20261118T072729Z` (and
   `nonneg_contrib_20261118T072729Z` for the contributions); older snapshots refuse with the fixed fit_cops.
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
* 2026-11-13 (movies; user ~22:20 CST "modify the movie: 1) the field of view of NG is 2pi half space (everything that is
  not blocked by the spacecraft body); 2) verify the field of view of RG: is that actually 45 degrees? what was used in the
  fitting model"; the lead verified that the fit uses RGAcceptance = hemisphere (front half-space, weight -l·b_RG, no 45°
  cut; the 45° is only the IK's nominal width); user choices): a NEW set `movies_softplus_smm1_9cycles_detach300_pan_hemi/`
  (the panorama set stays); RG = the fit's front half-space (2π) with tint -l·b_RG (0–1) minus directions blocked by the
  spacecraft body; NG = the front half-space of the NG boresight (2π) with tint the fit's front factor 1 - l·b_NG (1–2)
  minus directions blocked by the spacecraft body; the fit's srccg_shield part behind NG is NOT shown (said in the
  captions); blocking tested from the SPICE ROSINA point moved onto the +Z deck (it lies 4.9 cm inside the SPICE bus box;
  z = 2.655 m + 1 mm), against the spacecraft model (bus, both wings, HGA); a surface point is coloured for a gauge when it
  is in the gauge's front half-space (direction from the fit's gauge point, as the fit computes l and the weight), faces
  the gauge, has a free line of sight to the fit's gauge point and the segment deck point -> surface point misses the
  spacecraft; raw factors and colour bars kept (RG 0–1, NG 0–2); no translucent cone lines (a half-space has no cone),
  only the outline on the surface; the shares of the modelled signal the movie does not colour replace fov_share for
  this set; everything else as the panorama set (pose B, pan to frame 1800 then hold, flag on top, no white cross,
  Rosetta model, end at the first contact).
* 2026-11-13/14 (movies; user ~23:10 CST "prepare recolored movies; use the fixed log range; add 15 sec pause time after
  the last frame with the waving flag; closer to the end of the movies an isoline appeared on the surface -- remove that; at
  some point during descent, a circle is shown presenting the limited field of view of the RG. that circle is shown only for
  a few frames -- remove that circle at all"; lead's specification, user decision on the colour): a NEW set
  `movies_softplus_smm1_9cycles_detach300_pan_contrib/` (the half-space set stays): colour = the contribution per unit
  surface area exactly as the fitted model computes it (no spacecraft blocking; NG including its back side through the
  shield; RG its full front hemisphere; facets the gauge cannot see contribute 0 as in the kernels): NG density per m² of
  surface [m⁻³ m⁻²], RG flux into RG per m² [m⁻² s⁻¹ m⁻²]; computed per frame from fitting-base kernels at the frame
  times; fixed log range per gauge over the whole movie (bottom: 90 % of the signal at the first frame; top: 99.99th
  percentile at the end - taken over the last 10 frames with COPS data, because the model's RG surface flux falls
  from 5e15 m⁻² s⁻¹ at 10:39:28.9 to 1.4e-18 from 10:39:29.8 UTC, COPS 2.1 m: RG then sees only facets at the activity floor); below the bottom uncoloured relief, above the top the top colour; ALL field-of-view outlines and mask
  boundary lines removed (the "isoline" was the white half-space edge line on the surface; the RG "circle" was the 90°
  half-space edge drawn across the sky of the attached fisheye frames 688-900, and the 45° cone outline in the earlier
  cone sets); no cone lines; a 15 s hold of the last frame with the flag waving; everything else as the panorama /
  half-space sets (pose B, pan to frame 1800 then hold, flag on top, no white cross, Rosetta model, end at first contact).
* 2026-11-14 (user: "the range rule is OK"): the colour-range top is the 99.99th percentile of the per-area contribution
  over the last 10 frames WITH COPS data (frames 1790-1799, up to 10:39:26 UTC), not over the last frames of the movie;
  bottom unchanged (level enclosing 90 % of the first frame's signal). Ranges of the c9 set: NG 1.185e7-7.48e14
  m⁻³ m⁻², RG 4.81e9-2.27e15 m⁻² s⁻¹ m⁻² (`validation/contributions/color_range*.json`). Use this rule for later reruns.

* 2026-11-14 (user: "make comparison plot of the results with the data and the current softplus fit"; lead's
  specification): a new folder for the stage-1 option fits: per option all settings vs the data and c9; the best (lowest
  data cost) setting of each option + c9 + D2; a table per fit read from the runs' summaries; activity maps of the best
  settings and c9 side by side; CPU where possible; update when the remaining fits (C f1) finish.

* 2026-11-14 (~20:00 CST, user: "prepare the comparison plots for stage-1 fits completed to this point"; lead): v2 folder
  with the completed fits (C f 1 still running, left out and said so), convergence verdicts of the CPU diagnosis in the
  table and captions, a residual figure ln(model/obs) vs altitude (NG, RG; best of each + c9 + D2), the activity maps
  copied from the first folder (inputs unchanged); the first folder kept.

* 2026-11-15 (11:30 CST, lead: all stage-1 fits complete): the FINAL stage-1 comparison in a new folder (v1 / v2 kept) with
  the latest completed fit per setting (B/A continuations, D fixed-level penalty, C incl. f 1 and the P_m variant), the
  pre-continuation stage-1 costs in the table, convergence from the terminations (ftol = converged in the optimizer's
  sense, max_nfev = not converged), the residual figure, and new activity maps (best settings changed).

* 2026-11-17 (12:30 CST, lead): the final stage-2 comparison: regularized chains (and the unregularized ones as context),
  same figure style as the final stage-1 set (NG/RG vs data, late window, residuals, cost table with final cost /
  termination / active bounds), no maps (shared disk, small products); the two running C chains (P_m 1e-4) added when posted.

* 2026-11-17 (15:35 CST, lead; user: "yes, free the 44M store and generate the full sets", "I need the full set" - surface
  activity, production rate, contribution, ground track): full figure sets, exactly like the c9 set, for the six cases of
  12M_stage2_reg_best (softplus, A, B, D regularized, best C = f 1 / P_m 1e-4, D2), no movies, no comparison_vs_* folders;
  the needed adapters approved (15:50 CST) with the bitwise conditions.

* 2026-11-17 (20:30 CST, lead; user verbatim in 3.13): slope at equal illumination for the six stage-2 cases
  individually, definitions >40, >60, >80, 0-40 ... 80-90; Word document packing the same plot of the six cases in a cluster
  with a summary. Defaults set by the lead (not changed by the user): illumination = 10-deg SZA classes + shadowed + night;
  times = four snapshots + phase means; surface = gauge-seen (main) + whole surface (secondary); P3 cells 10 deg x 10 % x 0.1;
  answer = SZA-controlled ratio (area-weighted over the classes, weights = the category's area). 90-180 shown separately as an
  addition. My additions, labelled in the document: the under-sampled-reference flag and the variant ratio (sensitivity
  check; at 20 m / 5.5 m the specified ratio can exceed 1e100 because a few m^2 of 0-40 deg facets carry the reference).
* 2026-11-18 (00:15 CST, lead; user question verbatim in 3.14): illumination maps at 10:33:10 and 10:39:26 UTC: class, cos SZA, q/q_max (linear and log), whole surface and seen facets, in every production-rate layout with its relief and track; one side-by-side figure per time for the softplus reg case; one set for all cases after the cross-case identity check.
* 2026-11-18 (user, relayed by the lead; a model decision that changes every production-rate figure): "(b) -- use 0.02;
  x10 was not intended; use the same law but with decreased activity": H2O q_min = 0.02 q_max (was 0.2). Every figure and
  map made before 2026-11-18 shows the 0.2 law. Captions now take q_min / q_max from the maps' law (pub_maps_figures.py,
  pub_maps_v2_figures.py, illumination_maps_figures.py); nothing was redrawn (no refit yet).

* 2026-11-18 (14:00 CST, lead; user: "prepare the surface source rate maps for completed runs for inspection of the source
  rate"): per completed new-law run the source-rate set as in the full sets' production_rate/ and activity/ (LOWCUT 0.02,
  relief 0.7, seen-only) plus an inspection figure per time (new-law rate | old-law rate | new q/q_max); captions state the
  law and that the old figures used 0.2; a cluster page of all completed runs. Lead 14:12 CST: the old 0.2-law sets stay as
  they are (superseded; their map captions describe every model with the softplus formula), the new sets use each model's own
  text (`--model-text`). Mine (stated in the captions): common log classes for the new and old rate in the inspection figure.

* 2026-11-19 (17:05 CST, lead; user: "prepare a full set of comparison plots for the last round of fitting"): the set of
  3.10 for the new-law round: best of each option and all C chains vs COPS, per option new vs old law (law in every legend
  and caption), residuals for the new best set and new vs old, cost per cycle per chain, summary with kappa, Q_H2O and the
  cost change. Mine: new-vs-old residual figures only for softplus, A, B, D, D2 and the best C; the C f10 / 1e-4 Q_H2O total
  from its final maps and the seen part from its report (no inspection figure: no old-law maps).

* 2026-11-20 (14:45 CST, lead; user: "prepare the full set of plots for the last round of fits with the fixed surface source
  distribution; also prepare a single word doc summarizing all the fit"): full sets of the ten new-law states like 3.11 (the
  Word document by the docx agent). Mine: activity / production_rate / inspection copied from the source-rate sets (identical
  inputs) instead of regenerated; the maps' own geometry folders used (content-identical to job_136691); the refitted D2 as the
  model-data reference.

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
* 2026-11-14 (UTC): 3.7 the half-space set (`render/validate_fixed --sc-block-from`, `hemisphere` in the detach render with
  the cube-face line of sight, `scrows`, `fov_share.py --hemi-sc`; job variables SC_BLOCK_FROM, step scrows, HEMI_SC,
  HPREP, HMOVIES, DECK); section 5 decisions of 2026-11-13 (half-spaces minus the spacecraft, RG verified as the fit's
  hemisphere, shield part not shown, deck point, no cone lines).
* 2026-11-14 (UTC): 3.7 the contribution set (`render --color-contrib/--color-range`, `compose --hold-seconds`,
  `frames_base.sbatch`, `frame_contrib.py`, `contrib_range.py`; job variables COLOR_CONTRIB, COLOR_RANGE, HOLD_SECONDS,
  FIRST, TRG_CONST); section 5 decisions of 2026-11-13/14 (contribution colour, fixed log range, no outlines, no RG
  circle, 15 s hold).
* 2026-11-14: §5 user confirmation of the colour-range rule of the contribution movies (top from the last frames with COPS data).
* 2026-11-14 (UTC): 3.8 contribution density maps and region shares (`contrib_phase_data.py`, `pub_contrib_density_figures.py`,
  `regions_figures.py`, jobs `contrib_phase_data.sbatch`, `pub_contrib_density.sbatch`) in the c9 full set; the c9 full set,
  model_vs_observations and publication sets checked complete and current (nothing regenerated).
* 2026-11-14 (UTC): 3.9 stage-1 option comparison (`figures_12m_nonneg_options_stage1.sh`, `options_table.py`,
  `options_captions.py`, `options_howto.py`, `pub_options_activity_maps.py`, job `options_activity_maps.sbatch`); §5 decision.
* 2026-11-14 (UTC 2026-11-15): 3.9 v2 (`DIAG`, `RESID`, `ACT_FROM`; `options_residual_figure.py`, `data/stage1_convergence_20261114.json`; caption wording "completed so far", "stopped on ftol"; table column width); §5 decision.
* 2026-11-15 (UTC): 3.9 final (`RUNS`, `TEXTS`, `CONV_TERM`; `options_table.py` prev_* / penalty_mult / --conv-from-termination / setting width, `options_captions.py --texts`, `options_howto.py` ftol line; data files `stage1_final_runs_20261115.json`, `stage1_final_texts_20261115.json`; maps job 131968); §5 decision.
* 2026-11-17 (UTC): 3.10 final stage-2 comparison (`stage2_figures.py`, `stage2_howto.py`, spec `stage2_final_spec_20261117.json`; `options_residual_figure.py` one legend column for long labels); §5 decision.
* 2026-11-17 (UTC): 3.10 rebuilt with all nine regularized chains; `kappa_table.py` (kappa_check/).
* 2026-11-17 (UTC): 3.11 full figure sets of the stage-2 states and D2 (`full_set_submit.sh`, `finish_full_set.sh`, `full_set_facts.py`, `full_set_docs.py`, `d2_tools.py`, the alternative-model state_activity, figure-script fallbacks); §5 decision.
* 2026-11-17 (CST; UTC 2026-11-18): 3.12 Word document comparing the six stage-2 full sets plot by plot (`fullsets_comparison_docx.py`, separate venv `tools/docx_venv`; `full_sets_stage2_comparison/`).
* 2026-11-17 (CST; UTC 2026-11-18): 3.13 slope at equal illumination (`cliff_effect_data.py`, `cliff_effect_figures.py`, `cliff_effect_docx.py`, jobs `cliff_effect_{data,figures}.sbatch`, spec `cliff_effect_cases_20261117.json`; `cliff_effect_stage2/`); §5 decision.
* 2026-11-18 (CST): 3.14 illumination maps at the production-rate map times (`illumination_maps_figures.py`, jobs `illumination_{numbers,maps}.sbatch`; `illumination_maps_103310_103926/`); §5 decision.
* 2026-11-18 (01:45 CST): H2O night-side fix (q_min/q_max 0.2 -> 0.02): notes in 3.3, 3.14, 4 (step 2) and 5: every
  production-rate / illumination / map product so far used 0.2; the new snapshots nonneg_maps / nonneg_contrib
  20261118T072729Z; illumination_maps_figures.py QRATIO from the maps' code with derived classes (0.02: 10 linear classes
  from 0, 34 log classes from -1.70; to confirm); captions of pub_maps_figures.py / pub_maps_v2_figures.py data-driven.
  No figure redrawn. Diagnostic of the old states with the new law: tests/nightfix_diagnostic_20261118/ (FINDINGS A18).
* 2026-11-18 (13:35 CST): 3.3: maps on base_q002 (explicit --shape-manifest; SHAPE_MANIFEST / SHAPE_REFERENCE of
  jobs/nonneg_maps.sbatch; snapshots nonneg_maps / nonneg_contrib 20261118T192126Z; geometry job_136691; D2_q002 maps).
* 2026-11-18 (CST): 3.15 source-rate maps of the new-law runs (`source_rate_q002_{inspect,submit,docs,cluster}.py`, `jobs/source_rate_q002.sbatch`, spec `source_rate_q002_runs_20261118.json`; `pub_maps_v2_figures.py --model-text`); §5 decision.
* 2026-11-19 (CST): 3.16 comparison set of the new-law stage-2 round (`stage2_q002_figures.py`, `stage2_q002_docs.py`, `jobs/stage2_q002_figures.sbatch`, spec `stage2_q002_spec_20261119.json`); §5 decision.
* 2026-11-19 (CST): 3.17 Word document clustering the plots of the ten new-law runs (`q002_comparison_docx.py`, separate venv `tools/docx_venv`; `stage2_q002_comparison/`).
* 2026-11-20 (CST): 3.18 full figure sets of the new-law states (`full_set_q002.py`, `full_set_q002_finish.sbatch`, job-script overrides, `full_set_docs.py` spec fields); §5 decision.
* 2026-11-20 (CST): 3.19 ONE Word document summarizing all fits of the new-law round (`q002_all_fits_summary_docx.py`, `q002_all_fits_sheets.py`, `figures/data/q002_all_fits_observations_20261120.json`; `q002_all_fits_summary/`).
