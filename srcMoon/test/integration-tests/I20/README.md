# I20 — Permanent-shadow and polar cold-trapping physics

Purpose: Collocate independent D02 PSR and D03 polar temperatures with production residence statistics, warming release, and inventory sensitivity.

Run this ID through the aggregate runner so its I13 predecessor is evaluated
in the same immutable campaign directory:

```bash
python3 srcMoon/test/run_all_test.py --campaign-mode qualification --only I11 --only I12 --only I13 --only I20 --mpi-ranks 8
```

The runner always writes `stdout.log`, `stderr.log`, `application.in`,
`run_manifest.json`, and `result.json`. If I13 is not `PASS`, this test is
`SKIPPED` with `execution=not-run`; a U-series pass is not promoted to linked
evidence. Required external evidence: D02 LPSR/AVGVISIB and D03 PRP. Raw archive files are immutable,
all units/frames come from native labels, and missing package provenance or QA
is `SKIPPED`, never `PASS`.

The runner requires both package-level statuses, not merely generated CSVs.
D02 remains `PARTIAL`, while the canonical D03 package currently lacks its
required `README.md`, `provenance.json`, and `qa_report.json`. I20 also waits
for linked I15, I17, I18, and I19 PASS evidence, so it remains
`SKIPPED / NOT VALIDATED` until every earlier layer is satisfied.

An executed comparison must additionally preserve the exact row-level
`comparison.csv` and render those unchanged rows as `comparison.png` and
`comparison.eps`. A skipped/non-executed test creates none of these three
files; empty or synthetic comparison graphics are prohibited.
