# I16 — Diviner-driven surface-temperature boundary

Purpose: Score temperatures returned through the production surface callback against quality-filtered D03 GCP hold-out samples without gap filling.

Run this ID through the aggregate runner so its I13 predecessor is evaluated
in the same immutable campaign directory:

```bash
python3 srcMoon/test/run_all_test.py --campaign-mode qualification --only I11 --only I12 --only I13 --only I16 --mpi-ranks 8
```

The runner always writes `stdout.log`, `stderr.log`, `application.in`,
`run_manifest.json`, and `result.json`. If I13 is not `PASS`, this test is
`SKIPPED` with `execution=not-run`; a U-series pass is not promoted to linked
evidence. Required external evidence: D03 Diviner GCP. Raw archive files are immutable,
all units/frames come from native labels, and missing package provenance or QA
is `SKIPPED`, never `PASS`.

The runner reproducibly streams all 18 immutable GCP latitude-band products,
selects only native `tbol` records, preserves `-9998/-9999` as gaps, and
freezes a checkerboard training/hold-out split before model residuals are
examined. It derives an `application.in` selecting `diviner-map`; the linked
C++ adapter queries the ordinary production surface-temperature callback at
the SPICE epoch's physical local time. The scorer independently performs only
the adjacent native quarter-hour interpolation when both required hold-out
bins exist. It never fills a missing bin or substitutes a zero uncertainty.

Selection uses the actual body-fixed Sun state at the callback's SPICE epoch
with the same input `illumination_aberration_correction` as the map operator.
Eligibility is tested at the physical query local time, rather than only the
nearest archive-bin centre: a supported node does not authorize interpolation
through a missing adjacent corner. The frozen nearest-bin window and all
observational membership/quality criteria are preserved.

The staged input now explicitly records
`diviner_map_kinetic_interpretation = isothermal-grey`. Under the user-approved
isothermal grey resolved-facet approximation, native D03 bolometric brightness
`T_bol` becomes the production kinetic callback value
`T_kin=T_bol/thermal_infrared_emissivity^(1/4)`. This has no implicit fallback
for old map inputs and does not resolve anisothermality or spectral emissivity.
The staging receipt `diviner_boundary_interpretation.json` records this policy
and input/data hashes. Real table gaps remain gaps; unavailable temperatures
are not filled with an initial temperature.

`receipt.txt.I16_model.csv` retains kinetic temperature [K], forward-model
bolometric brightness [K] and dimensionless emissivity separately. The scorer
independently checks the forward identity and compares only the brightness
column with D03 brightness. A stale kinetic-only receipt or inconsistent
conversion is ERROR. This changes neither frozen hold-out membership nor the
15 K RMSE criterion. The local protocol fixtures run with
`/usr/bin/python3 -B -m unittest discover -s srcMoon/test -p test_i16_grey_boundary.py`;
they read no observations and cannot make I16 PASS. Native operator/cache
checks require the newly rebuilt U10 production probe.

The individual `I16/test.py` runner always executes those local receipt checks
before delegating to the shared scientific runner. Thus `run_all_test.py`
includes them even when I13 or D03 makes scientific I16 `SKIPPED`. Their logs
and infrastructure-only summary are under `I16/protocol_checks/`; success
does not override scientific status. A subordinate check failure writes the
per-ID `result.json` as `ERROR` and prevents a linked model launch.

The generated `d03_gcp_qa.json` qualifies only that reproducible processing
step. At present the canonical external `diviner/` package has no package-level
`README.md`, `provenance.json`, or `qa_report.json`; therefore I16 stops as
`SKIPPED` before it reads the generated table. Adding downloaded TAB/XML files
alone cannot change that status. When the three package records exist, their
D03 status must agree and be exactly `READY`, and their hashes are added to the
linked run manifest.

`PASS` requires the processed QA status and hashes to be READY/correct, every
evaluated row to retain `GOOD`/`hold-out` selection, at least one eligible
comparison, and non-PSR RMSE no greater than the roadmap's frozen 15 K gate.
Mean bias and maximum absolute residual are reported but are not hidden
secondary tuning criteria.

An executed comparison must additionally preserve the exact row-level
`comparison.csv` and render those unchanged rows as `comparison.png` and
`comparison.eps`. A skipped/non-executed test creates none of these three
files; empty or synthetic comparison graphics are prohibited.
