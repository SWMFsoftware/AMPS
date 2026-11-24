# Production trajectory exposure controls

Run from the AMPS root:

```bash
python3 -B srcMoon/test/test_m3_trajectory_exposure.py --output-dir /tmp/moon-trajectory-NEW
```

The output must be new. The runner compiles the actual `TrajectoryExposure`
and `PiecewiseHazard` production sources with strict warnings, preserves each
query's stdout/stderr, and writes `result.json`, `run_manifest.json` and
`comparison.csv`, `comparison.png` and `comparison.eps`. Exit codes are
0/1/2/77 for PASS/FAIL/ERROR/SKIPPED. These are
analytical numerical fixtures, not observations or full I21–I24 qualification.
The global Moon supporting runner registers this group as `trajectory-exposure`.

The prescribed radial refinement writes separate `radial_refinement.csv`,
PNG and EPS files. It compares the production clock against an analytical
first hit at 0.8 s at 8, 16 and 32 steps and checks a surviving threshold
on the same 2 s horizon.
[I21 varying-exposure contract](../integration-tests/I21/VARYING_EXPOSURE_CONTRACT.md)

Linear and accelerated radial arcs have independent elementary inverse-square
integrals. Finite-distance opaque-sphere shadow entry/exit is derived from the
line/sphere tangency equation; the short-shadow case has two visible endpoints
and proves that endpoint-only sampling is insufficient. The terrain fixture
uses a declared rectangular occulting screen and an analytical visibility
callback to isolate the integrator. Actual production terrain/BVH and controller
wiring require separate native controls; this fixture does not claim those.

The smooth approximation bound follows `|integral(f-linear(f))| <= h³ sup|f''|/12`
on each interval. The point-Sun silhouette-root time envelope is reported
separately. The tests compare against the sum of these bounds plus an explicit
binary64 comparison envelope, never a changed physical acceptance threshold.
All rates, distances, geometry and tolerances in the fixtures are declared SI
inputs. Invalid input, missing geometry, excessive requested root error and
exhausted segment capacity must refuse without publishing a partial history.
The full test specification is [M3_TRAJECTORY_EXPOSURE_CONTRACT.md](../M3_TRAJECTORY_EXPOSURE_CONTRACT.md).

For actual PIC/controller wiring, reuse a fresh successful all-species native
producer BEFORE its ordinary-Moon restore, in the same disposable checkout:

```bash
python3 -B srcMoon/test/run_m3_trajectory_exposure.py \
  --producer-manifest /absolute/path/native-particle-backend/run_manifest.json \
  --output-dir /absolute/path/native-trajectory-NEW --wall-seconds 600
```

The runner verifies the actual library, generated ABI and maintained production
source identities, then runs both mover orders with one and eight MPI ranks.
It writes52 actual controls with per-rank event CSVs, logs and manifests. Its
spatial-rate and zero-Q packet are explicitly mathematical fixtures; the test
does not activate a physical Na model. `run_m3_development.py --isolated-checkout`
includes this runner in its single supervised native sequence and restores
ordinary Moon afterward. That option must be used in a disposable checkout;
it acknowledges ownership of derived configuration, rather than creating a
temporary checkout automatically.

Independent review of the first implementation found two counterexamples:
nearby polynomial roots could merge and understate a short shadow, and a far
Sun could cancel an edge normal computed from two large Sun-centred vectors.
The maintained implementation now constructs outward interval coefficients,
uses the edge form with a nearby plane origin, and proves a strict polynomial
sign on every complement of the unresolved root brackets. It charges the
union width of those actual brackets to the exposure error budget; it uses no
root-count-times-epsilon estimate. `maximumGeometryIntervals` bounds that work
explicitly. Cases that cannot meet the declared depth budget refuse instead
of publishing an underestimated bound. Fast-math builds are rejected for the
interval service. The far diagonal Sun and close-root examples are permanent
controls. Interval geometry and BVH/intersection accuracy remain distinct
numerical layers; formal terrain gates are not waived.

Native controls also exercise a7e-6 s^-1 analytical hazard through many accepted
flights, two caller intervals and MPI continuation. The exact compensated
clock keeps both components; it does not force them through a single long
double. Provider approximation and native clock-storage error bounds persist
with particle histories and accepted-event ledgers. Error diagnostics are
numerical bounds, not observational uncertainties.

A translated sub-ULP transverse-flight fixture also checks the distinction
between the ideal quadratic arc and its binary64 callback point. On every
root-free cell the service certifies identical strict signs for all silhouette
constraints at the ideal midpoint and represented callback point. Any cell
whose classification cannot be certified is charged in full, together with
the root-bracket union. The wide-budget fixture must cover the independently
derived approximately 0.8 optical-depth discrepancy; the restrictive-budget
fixture must refuse. This position-projection bound does not cover unknown
BVH/ray arithmetic or external frame-transform errors; those remain separate
geometry qualification gates. No threshold has been relaxed.
