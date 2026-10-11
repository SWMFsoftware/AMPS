# U05 — Sodium radiation-pressure geometry and shadow kernels

The complete executable contract—including purpose, production entry points,
configuration, units, frames/time, oracle, controls, procedure, metrics,
acceptance criteria, failure modes, data/hashes, predecessor gates, artifacts,
and status semantics—is in [reference/acceptance.json](reference/acceptance.json).

Run from the AMPS repository root:

```bash
python3 srcMoon/test/U05_na_radiation_shadow/test.py
```

Exit codes are `0=PASS`, `1=FAIL`, `2=ERROR`, and `77=SKIPPED`. A skipped
capability is intentionally not converted into a passing source-preflight result.

## Independent radiation-pressure reference

The original machine-readable digitization used to populate
`src/species/Na.cpp` is unavailable. U05 therefore does not freeze the current
array and call it truth. The independent absolute check uses the upper panel of
Figure 7 in Combi, DiSanti, and Fink, *Icarus* 130 (1997), DOI
`10.1006/icar.1997.5832`, as reproduced on printed page 13/PDF page 25 of NASA
NTRS report `19990024948`.

The reference directory records two pixel-coordinate traces made from 150 dpi
and 300 dpi rasterizations. Points hidden by the publication's arrows or
markers are excluded solely by the documented legibility rule. The preparation
script performs the axis conversions, checks that both traces agree within
`0.5 cm s^-2`, averages them, and writes the frozen CSV, provenance, QA report,
and checksums. No AMPS output or value from `Na.cpp` is read while constructing
the reference.

To reproduce the reference after downloading the authoritative report:

```sh
python3 srcMoon/test/U05_na_radiation_shadow/reference/prepare_reference.py \
  --source-pdf /path/to/NTRS-19990024948.pdf
```

The expected PDF SHA-256 is
`af3becc93b3c3f1bd1a7d0b208a80b2fd7aa4d8399c31e91fc9ac9e121dc30d9`.
The script rejects a different file. U05 verifies the hashes of every committed
reference artifact before executing the production kernel.

The `0.5 cm s^-2` acceptance band is a graphical digitization uncertainty, not
an adjustable model tolerance. U05 PASS qualifies the published curve at that
resolution; it does not reconstruct the unavailable historical table and does
not replace I07's full-mover shadow-gating test.
