# U04 — Third-body gravity and rotating-frame kernels

The complete executable contract—including purpose, production entry points,
configuration, units, frames/time, oracle, controls, procedure, metrics,
acceptance criteria, failure modes, data/hashes, predecessor gates, artifacts,
and status semantics—is in [reference/acceptance.json](reference/acceptance.json).

Run from the AMPS repository root:

```bash
./Config.pl -application=moon \
  -spice-path=/home/vtenishe/SPICE/cspice \
  -spice-kernels=/home/vtenishe/SPICE/Kernels
rm -rf build
make -j
python3 srcMoon/test/U04_rotating_frame/test.py \
  --output-dir test_output/srcMoon/U04
```

The fixed epoch is `2009-01-24T00:00:00 UTC`. The declared coordinate
contract is J2000 inertial, Moon-centred LSO rotating, and MOON_ME_DE421
body-fixed, with geometric aberration correction `NONE`. The exact seven-file
DE421 kernel list and SHA-256 values are in `reference/acceptance.json`.

The production side of the comparison calls the helpers used by
`Moon::TotalParticleAcceleration` and obtains LSO angular velocity through
`xf2rav`. The independent oracle evaluates the differential point-mass and
cross-product equations directly, and reconstructs angular velocity from the
lower-left derivative block of the SPICE state transformation using
`W=(dR/dt)R^T`. A non-axial fixture makes sign or component-order mistakes
observable. This is local linked-kernel verification, not the full AMPS mover
trajectory required by I06/I10.

Artifacts are written below `<output-dir>/U04/`:

- `kernel_hashes.json` records paths, sizes, expected and actual hashes;
- `term_vectors.json` records every scalar component used by the acceptance
  decision;
- `result.json` records compilation, stdout/stderr, revision, and final status.

Kernel availability is checked before CSPICE is called. A missing declared
kernel is `SKIPPED`; a hash mismatch or infrastructure/parser failure is
`ERROR`; evaluated numerical disagreement is `FAIL`.

Exit codes are `0=PASS`, `1=FAIL`, `2=ERROR`, and `77=SKIPPED`. A skipped
capability is intentionally not converted into a passing source-preflight result.
