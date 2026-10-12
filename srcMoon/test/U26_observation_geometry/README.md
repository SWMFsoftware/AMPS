# U26 — Observation-operator and instrument geometry

The complete executable contract—including purpose, production entry points,
configuration, units, frames/time, oracle, controls, procedure, metrics,
acceptance criteria, failure modes, data/hashes, predecessor gates, artifacts,
and status semantics—is in [reference/acceptance.json](reference/acceptance.json).

Run from the AMPS repository root:

```bash
python3 srcMoon/test/U26_observation_geometry/test.py
```

Exit codes are `0=PASS`, `1=FAIL`, `2=ERROR`, and `77=SKIPPED`. A skipped
capability is intentionally not converted into a passing source-preflight result.

The generic CLI `application: moon` block now carries the agreed adjustable observation-geometry
receipt (run UTC, common altitude grid, LSO SZA list, body-fixed pole and shell
flags, solid-angle shell weighting, projected-Sun limb branch selection,
anti-lunar map width, J2000/ECLIPJ2000 coordinate selection, display
orientation, circular-field diameter, exposure duration, unresolved beam
response, and exposure-average requirement). This closes configuration
ambiguity only. The current production samplers still lack the complete
pole/shell operators and true finite-beam/time convolution, and the
instrument-specific authoritative geometry cases remain incomplete; U26
therefore remains `SKIPPED`, not `PASS`.

The D14 package supports the committed default values of a 1 degree circular
WHAM field and 120 s exposure, and supplies both J2000 RA/Dec and ECLIPJ2000
longitude/latitude for the four published brightest beams. It does not supply
the full numerical raster, per-beam timestamps, spectra, or a radial beam
response. The parser consequently requires `tail_beam_response=unresolved`;
no assumed Gaussian or uniform-disk response is accepted as validation input.
