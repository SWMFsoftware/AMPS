# CME model-local registry adapters

These adapters execute model-owned registries and translate only their typed
result records into the global schema. They do not duplicate physics expected
values. Every local registry is validated in full before a global `--case`
selection is applied, so selecting one case cannot conceal local registry
drift.

The staged coronal-CME adapter imports 222 rows as **software verification**;
its three campaign-evidence cases remain typed `SKIP`, and the local
`production_qualified=false` and `observational_campaign_qualified=false`
claims are retained. The full-volume Corona--SWCME adapter imports 15 bounded
code markers, but BG3D-4 and G1 remain unqualified and deferred. The reduced
adapter imports 97 scope-qualified rows. Routine global execution uses the
model runner's `--skip-native` mode, producing 62 executed portable rows and
35 explicit native `SKIP` rows rather than silently dropping native scope.

All Makefile-driven sources are copied below the suite work directory before
building. This excludes stale checkout objects and supports read-only source
mounts. The authoritative local JSON/JUnit/text products are retained below
the global artifact directory.

Run the registered suites with:

```sh
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-MODEL-SEP-CORONAL-CME-001
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-MODEL-SEP-CORONA-SWCME-FULL-001
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-MODEL-SEP-CORONA-SWCME-REDUCED-001
```

A complete native reduced-front campaign remains a separate model-local run:

```sh
python3 srcSEP3D/test/run_reduced_shock_front.py --output-dir /absolute/evidence/path
```

That campaign is zero-particle and reduced-front-plus-ambient only. Neither a
geometric observer arrival nor its PASS summary asserts an accepted fast shock
or a spatial sheath/contact/ejecta volume.

