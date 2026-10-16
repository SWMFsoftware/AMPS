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
