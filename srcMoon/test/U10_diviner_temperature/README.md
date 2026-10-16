# U10 — Diviner ingest and interpolation

The complete executable contract—including purpose, production entry points,
configuration, units, frames/time, oracle, controls, procedure, metrics,
acceptance criteria, failure modes, data/hashes, predecessor gates, artifacts,
and status semantics—is in [reference/acceptance.json](reference/acceptance.json).

Run from the AMPS repository root:

```bash
python3 srcMoon/test/U10_diviner_temperature/test.py
```

Exit codes are `0=PASS`, `1=FAIL`, `2=ERROR`, and `77=SKIPPED`. A skipped
capability is intentionally not converted into a passing source-preflight result.

## Assumptions, physics, and implementation decisions

U10 is local verification, not a comparison with Diviner observations.  Its
2-by-2-by-2 temperature table is an affine analytical fixture used solely to
prove the production CSV loader, complete-cell rule, periodic longitude/local-
time axes, and uncertainty propagation.  It is never stored in D03 and is
never described as measured lunar temperature.

The bolometric fixture follows the published Paige et al. C3--C9 construction:
the wavelength partitions end at 8.075, 8.40, 13, 25, 50, 100, and 1000
micrometres, and C3--C8 sensitivity limits are 190, 190, 180, 95, 60, and 40 K.
A channel below its limit inherits the next longer usable channel; C9 is
mandatory.  The independent expected value deliberately remains slightly
below 55 K because the published final partition stops at 1000 micrometres.

The model fixture distinguishes kinetic temperature from a blackbody-
equivalent observation.  A grey surface at kinetic temperature `T` and
bolometric emissivity `epsilon` emits `epsilon sigma T^4`, so its comparable
bolometric brightness temperature is `epsilon^(1/4) T`.  Directly comparing
`T` with Diviner `T_bol` would silently assume unit emissivity and is forbidden.

`test_processor.py` additionally verifies the exact frozen I17 software
decisions without opening observational data: locale-independent SCET UTC
parsing; native qca/qge/qmi bit-mask semantics; the exclusive emission-angle
cut; averaging paired native temperatures at a rounded radiance before forming
the monotone inverse; the Paige substitution chain; and branch-preserving
first-order uncertainty propagation.  The last fixture deliberately gives C9
a one-sigma uncertainty larger than its positive mean: the implementation must
differentiate locally on the valid positive-temperature branch, not evaluate a
fabricated negative Planck temperature.  C3--C8 inputs already removed by the
published low-sensitivity substitution have zero local derivative.  The test
also verifies that this remains true when those unused channels carry large
uncertainties.  These deterministic fixtures can verify processing logic but
cannot make D03 READY or make I17 PASS.  The large RDR processor still has to
preserve native labels, hashes, real gaps, actual channel
footprints/times/detectors, and the untouched December hold-out under I17.
