# I15 — Topographic illumination and shadowing

Purpose: Compare the production terrain-horizon and Earth-eclipse service with independent D02 polar visibility/PSR evidence and published temporal sampling.

Run this ID through the aggregate runner so its I13 predecessor is evaluated
in the same immutable campaign directory:

```bash
python3 srcMoon/test/run_all_test.py --campaign-mode qualification --only I11 --only I12 --only I13 --only I15 --mpi-ranks 8
```

The runner always writes `stdout.log`, `stderr.log`, `application.in`,
`run_manifest.json`, and `result.json`. If I13 is not `PASS`, this test is
`SKIPPED` with `execution=not-run`; a U-series pass is not promoted to linked
evidence. Required external evidence: D01 LOLA topography; D02 AVGVISIB/LPSR. Raw archive files are immutable,
all units/frames come from native labels, and missing package provenance or QA
is `SKIPPED`, never `PASS`.

An executed comparison must additionally preserve the exact row-level
`comparison.csv` and render those unchanged rows as `comparison.png` and
`comparison.eps`. A skipped/non-executed test creates none of these three
files; empty or synthetic comparison graphics are prohibited.
