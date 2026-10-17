# Enceladus Stage-0 neutral-plume runbook

## 1. Case description

This standalone application fits a steady, collisionless, gravity-free neutral plume to Cassini observations without building AMPS. It supports the `tenishev2014_hybrid_8` case (eight discrete jets, seven tiger-stripe line components, and two empirical backgrounds), the `porco2014_98` case, and the `portyankina2022_solar_41` case. The default source law is the production-normalized Tenishev (2010) point-source law. H2O is fitted to E3/E5 INMS density and the 2007/2010 UVIS occultations; CO2 is fitted to E3, E5, E14, E17, and E18. The local `e17_extra_jet_candidate.json` is an optional 99th-source sensitivity case and is never enabled unless `--extra-jet` is explicitly present.

All observation products and generated results remain outside the code tree. `ENCELADUS_VALIDATION_DATA_ROOT` must point to the machine-local Cassini validation-data directory containing `data/processed/plume_fit_v1`, `data/processed/co2_fit_v1`, `data/derived/porco2014_jet_catalog.json`, and `studies`. Every fitting output is opened exclusively; select a new `RUN_TAG` or remove/move an unwanted prior run deliberately before repeating these commands.

## 2. Complete commands

Each command below is one complete shell line. Run them in order. Replace only the two `/absolute/path/...` values and the run tag.

### Configure the machine and Python environment

```bash
export ENCELADUS_STANDALONE_ROOT=/absolute/path/to/AMPS/srcEnceladus/standalone
export ENCELADUS_VALIDATION_DATA_ROOT=/absolute/path/to/cassini_validation_data
export RUN_TAG=enceladus_stage0_001
cd "$ENCELADUS_STANDALONE_ROOT"
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
test -f "$ENCELADUS_VALIDATION_DATA_ROOT/data/processed/plume_fit_v1/release_manifest.json" -a -f "$ENCELADUS_VALIDATION_DATA_ROOT/data/processed/co2_fit_v1/release_manifest.json" -a -f "$ENCELADUS_VALIDATION_DATA_ROOT/data/derived/porco2014_jet_catalog.json"
.venv/bin/python -m unittest discover -s tests -v
mkdir -p "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/single_fit"
```

### Fit the default hybrid model independently to E3 and E5 H2O

```bash
.venv/bin/python plume_fit.py --config config/tenishev2014_stage0.json --source-case tenishev2014_hybrid_8 --validation-data-root "$ENCELADUS_VALIDATION_DATA_ROOT" --flyby e3 --output "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/single_fit/e3_tenishev2014_hybrid_8.json"
.venv/bin/python plume_fit.py --config config/tenishev2014_stage0.json --source-case tenishev2014_hybrid_8 --validation-data-root "$ENCELADUS_VALIDATION_DATA_ROOT" --flyby e5 --output "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/single_fit/e5_tenishev2014_hybrid_8.json"
```

### Run the complete H2O INMS/UVIS study for all three source cases

```bash
.venv/bin/python joint_fitting_study.py --config config/tenishev2014_stage0.json --source-case tenishev2014_hybrid_8 --validation-data-root "$ENCELADUS_VALIDATION_DATA_ROOT" --output-dir "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/h2o_uvis/tenishev2014_hybrid_8" --quadrature-order 192
.venv/bin/python joint_fitting_study.py --config config/tenishev2014_stage0.json --source-case porco2014_98 --validation-data-root "$ENCELADUS_VALIDATION_DATA_ROOT" --output-dir "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/h2o_uvis/porco2014_98" --quadrature-order 192
.venv/bin/python joint_fitting_study.py --config config/tenishev2014_stage0.json --source-case portyankina2022_solar_41 --validation-data-root "$ENCELADUS_VALIDATION_DATA_ROOT" --output-dir "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/h2o_uvis/portyankina2022_solar_41" --quadrature-order 192
.venv/bin/python compare_source_case_studies.py --study-root "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/h2o_uvis"
```

### Run the five-flyby CO2 study for all three source cases

```bash
.venv/bin/python co2_fitting_study.py --config config/tenishev2014_stage0.json --source-case tenishev2014_hybrid_8 --validation-data-root "$ENCELADUS_VALIDATION_DATA_ROOT" --h2o-study-result "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/h2o_uvis/tenishev2014_hybrid_8/study_results.json" --output-dir "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/co2/tenishev2014_hybrid_8"
.venv/bin/python co2_fitting_study.py --config config/tenishev2014_stage0.json --source-case porco2014_98 --validation-data-root "$ENCELADUS_VALIDATION_DATA_ROOT" --h2o-study-result "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/h2o_uvis/porco2014_98/study_results.json" --output-dir "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/co2/porco2014_98"
.venv/bin/python co2_fitting_study.py --config config/tenishev2014_stage0.json --source-case portyankina2022_solar_41 --validation-data-root "$ENCELADUS_VALIDATION_DATA_ROOT" --h2o-study-result "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/h2o_uvis/portyankina2022_solar_41/study_results.json" --output-dir "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/co2/portyankina2022_solar_41"
.venv/bin/python "$ENCELADUS_VALIDATION_DATA_ROOT/scripts/plot_co2_h2o_tidal_study.py" --co2-study-root "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/co2" --h2o-study-root "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/h2o_uvis" --output-dir "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/co2/figures"
```

### Optional: reproduce the E17 location search

```bash
.venv/bin/python extra_jet_e17_study.py --config config/tenishev2014_stage0.json --validation-data-root "$ENCELADUS_VALIDATION_DATA_ROOT" --co2-baseline-result "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/co2/porco2014_98/study_results.json" --output-dir "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/e17_extra_jet_search"
```

### Optional: repeat the complete Porco workflow with the fixed candidate as source 99

```bash
.venv/bin/python joint_fitting_study.py --config config/tenishev2014_stage0.json --source-case porco2014_98 --validation-data-root "$ENCELADUS_VALIDATION_DATA_ROOT" --extra-jet e17_extra_jet_candidate.json --output-dir "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/porco99/h2o_uvis" --quadrature-order 192
.venv/bin/python co2_fitting_study.py --config config/tenishev2014_stage0.json --source-case porco2014_98 --validation-data-root "$ENCELADUS_VALIDATION_DATA_ROOT" --extra-jet e17_extra_jet_candidate.json --h2o-study-result "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/porco99/h2o_uvis/study_results.json" --output-dir "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/porco99/co2"
.venv/bin/python compare_extra_jet_studies.py --base-h2o-uvis "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/h2o_uvis/porco2014_98/study_results.json" --augmented-h2o-uvis "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/porco99/h2o_uvis/study_results.json" --base-co2 "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/co2/porco2014_98/study_results.json" --augmented-co2 "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/porco99/co2/study_results.json" --output-dir "$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG/porco99/comparison_98_vs_99"
```

## 3. Expected results

The complete run is written below `$ENCELADUS_VALIDATION_DATA_ROOT/studies/$RUN_TAG`. No observation or source file is modified.

- `single_fit/` contains the two independent default-hybrid JSON results. Each records the resolved configuration, component rates, modeled density, residuals, rank and conditioning diagnostics, metadata, and input hashes.
- `h2o_uvis/<source-case>/` contains `study_results.json`, `fitting_study.md`, component and source-rate CSV tables, cross-prediction and ablation tables, and `fit_profiles.png`/`.pdf`. `h2o_uvis/` also receives the cross-case comparison products from `compare_source_case_studies.py`.
- `co2/<source-case>/` contains five observed/model CO2 profiles, all fitted-rate strategies, H2O/CO2 ratio and H2O-proxy tables, diagnostics, and `co2_fit_profiles.png`/`.pdf`. `co2/figures/` contains CO2 and direct-plus-proxy H2O plots versus mean anomaly, second-focus angle, and crack-opening proxy.
- `e17_extra_jet_search/` contains the surface and tiger-stripe searches, time-window and direction sensitivity tables, cross-flyby tests, JSON provenance, and PNG/PDF figures. The selected location is a diagnostic candidate and not an independent vent detection.
- `porco99/h2o_uvis/` and `porco99/co2/` contain complete 99-source refits. `porco99/comparison_98_vs_99/` contains the 31 direct 98-versus-99 comparisons, strategy summaries, heatmaps, and `comparison_report.md`.

A successful run does not imply that every fitted coefficient is physically identifiable. Inspect each generated `fitting_study.md`, condition number, effective rank, and cross-prediction metric. In particular, the Porco 98/99 and Portyankina 41 inversions are ill-conditioned, the CO2-to-H2O conversion is a conditional composition proxy, and the extra source should remain disabled unless the explicit 99-source sensitivity case is intended.
