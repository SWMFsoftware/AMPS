# I14 — LOLA topography integration into production surface boundary

Purpose: Prove that production particle segments intersect the time-dependent body-fixed LOLA triangulation and meet the frozen D01 control-point, topology, and ray-hit gates.

Run this ID through the aggregate runner so its I13 predecessor is evaluated
in the same immutable campaign directory:

```bash
python3 srcMoon/test/run_all_test.py --campaign-mode qualification --only I11 --only I12 --only I13 --only I14 --mpi-ranks 8
```

The runner always writes `stdout.log`, `stderr.log`, `application.in`,
`run_manifest.json`, and `result.json`. If I13 is not `PASS`, this test is
`SKIPPED` with `execution=not-run`; a U-series pass is not promoted to linked
evidence. Required external evidence: D01 LOLA LDEM_4. Raw archive files are immutable,
all units/frames come from native labels, and missing package provenance or QA
is `SKIPPED`, never `PASS`.

The linked adapter reports the exact production triangulation size, all-face
orientation/degeneracy counts, 512 approximately uniform vertex samples,
body-centred radial BVH intersections, and a real second-order-mover surface
crossing on every MPI rank that owns a suitable facet. The Python scorer
independently decodes the signed little-endian LDEM_4 raster using the detached
label's 0.5 m scale, nearest native latitude row and periodic longitude column.
It never reads the production triangulation writer as its reference.

`PASS` requires maximum vertex-versus-D01 radius error at most 50 m, radial
ray-versus-vertex error at most the input collision-distance tolerance, zero
degenerate or inward faces, a hit for every radial query, and a production
mover collision on every rank that ran the fixture. D01 must be READY in both
its QA and provenance records and every frozen raw SHA-256 must match.

An executed comparison must additionally preserve the exact row-level
`comparison.csv` and render those unchanged rows as `comparison.png` and
`comparison.eps`. A skipped/non-executed test creates none of these three
files; empty or synthetic comparison graphics are prohibited.
