# Stage 12: offline observation inference and provenance

This implementation guide maps Stage 12 and `PROV3D01--03`/`CAL3D01` in
`model/testing_validation.md` to executable code. The canonical r4 model
remains the specification. No event, spacecraft coordinates, observation
epochs, or measured calibration constants are embedded in C++.

## Code ownership and execution

`tools/preprocess_observations.py` is the offline CLI.
`tools/preprocessing/core.py` owns deterministic serialization, covariance
algebra, source ownership, immutable publication and verification.
`inference.py` owns parameter recovery and uncertainty propagation;
`campaign.py` owns the preregistered five-factor candidate product and
likelihood; `protocols.py` owns reviewed field-line requests and frozen
second-event/analytic-versus-MHD manifests. All use the Python standard
library. Python is never called by the AMPS cell or particle loops.

## Inference products and stated uncertainty

Every operator receives an explicit source JSON, parameters and metadata. It
produces a neutral JSON asset retaining the source checksum and the complete
inference declarations. Covariance matrices must be finite, symmetric and
positive definite; singular fits and unsupported extrapolation fail. General
least squares retains correlations rather than replacing covariance by
independent error bars. Nonlinear transforms propagate their Jacobians; the
positive-profile log transform is a stated first-order approximation and
rejects large relative input errors.

| Operator | Prepared product and uncertainty authority |
| --- | --- |
| `magnetogram-harmonics` | Real orthonormal Condon--Shortley coefficients from generalized least squares, including a fitted monopole diagnostic and coefficient covariance; the supplied removal limit gates the monopole before use. |
| `image-ellipse` | Projected conic/ellipse parameters and propagated covariance, with explicit pixel or physical scale. One image is a projected constraint, not a uniquely reconstructed 3-D CME. |
| `reconstructed-ellipsoid` | Three positive axes fitted to independently reconstructed 3-D points, conditional on the declared center and SO(3) attitude; fixed parameters are explicitly retained. |
| `front-kinematic-fit` | Linear or quadratic center/axis histories with joint fit covariance, covered time support and positive-axis checks, conditional on the declared fixed attitude. |
| `ephemeris-transform` | Position, velocity and full 6-by-6 covariance under an explicit proper rotation/translation; source/target frames and transformation epoch are stored. |
| `power-law-profile` | Positive density, temperature, speed or pressure profiles with reference radius, exponent, full log-parameter covariance and bounded support. |
| `plasma-sheet-width` | Gaussian angular width/contrast fitted above an independently declared baseline, with fit covariance. This is an inference asset for base-state construction, not a cell density overwrite. |
| `shock-compression` | Density ratio and uncertainty using the joint upstream/downstream covariance. |
| `wsa-speed-fit` | WSA-like target speeds under preregistered shape parameters, with fitted slow/fast speeds and covariance. This is an empirical target, not a substitute wind dynamics solution. |
| `critical-mach-table` | Ordered beta/obliquity table, full value covariance, declared Mach convention/EOS and bounded bilinear interpolation. |
| `instrument-response` | Response-folded product and full covariance; the response matrix is checksummed and its physical definitions retained. |
| `declared-table` | Closure, boundary, kinematic, normalization or comparison values with explicit SI column units/definitions, support, and full covariance or named ensemble members. |
| `formation-constraint` | Time-resolved radio frequency/covariance and explicit harmonic hypotheses, or a joint height/density inversion with the density-model checksum and full joint covariance. |

The synthetic recovery tests compare against independently constructed maps,
geometries and transforms. They do not estimate their expected answers by
calling the same fitting function.

## Immutable assets and data-use roles

A job has schema `sep-observation-preprocess-job-v1`, a processing version,
reference UTC epoch, output frame and asset requests. Each request gives an
operator, parameters, source file and metadata. Required metadata include
asset ID/kind, source URI/SHA-256, units, frame, epoch, processing version,
support, independence ID and exactly one data-use role:

- `construction`: determines a background, front, closure or trace seed;
- `qualification`: independently gates or weights candidates;
- `withheld-validation`: may be evaluated after selection is frozen.

An input file's raw bytes must reproduce the declared SHA-256. Reusing those
bytes, an independence group or an instrument response across roles fails
even after renaming the asset. Mixed processing versions/reference epochs,
unknown units/frames, duplicated IDs and improper frame rotations fail.
The job explicitly states when an ephemeris is transformed into another
frame; an implicit conversion is never inferred from its label.

Publication atomically creates a new directory containing canonical assets,
the original source bytes and a checksum-owned `manifest.json`. Members are
sorted, JSON keys are sorted, nonfinite numbers/duplicate keys are rejected,
and no timestamp or host path is added. Processing identical bytes and
declarations gives identical output bytes. The verifier checks the complete
manifest, each member, the retained originals and the role/epoch contracts.
An existing output is never silently replaced. Frozen records carry an
identity over every field; editing a weight or a declaration invalidates it.

## Formation-height campaigns

Preregister all members and priors of the exact product
`magnetogram x field_scale x coupling_radii x wind_density x front`. The
coupling-radii asset represents the complete `(R_b,R_i,R_scs)` authority.
Freeze the independent constraint, hard D6/topology/D1 limits, the weighting
rule, and the density construction procedure/equations checksum before
evaluating the realized product. The frozen record binds all asset content.

Each realized row must contain its five-member tuple exactly once, independent
qualification assets for D6/topology/D1/D2, all diagnostic results, separate
first-fast/first-supercritical heights and times, sampled front radii, and a
`density_construction` record. That record binds the same tuple, frozen
procedure/equations and original source checksums for all five axes.
Missing, duplicated or extra tuples fail. Rejected members remain in the
published selection with zero weight and their diagnostics intact.

The CLI verifies and scores a fully realized product; the background builder
must rebuild each tuple using the frozen construction equations. It cannot
copy one density or scale an Alfvén speed across changed field/radii members.
Supply either a covered positive power-law reconstruction with covariance,
or `sep-front-sampled-density-v1` containing the exact front radii/densities,
full covariance and an explicit geometry uncertainty authority. The latter
allows fully prepared wind/background outputs to enter the likelihood without
fitting them to a universal density approximation.

For the preferred radio case, `f_pe=(1/2pi)*sqrt(n*e^2/(epsilon0*m_e))` is
recomputed from each tuple's density. Fundamental and harmonic predictions
have preregistered probabilities. Their Gaussian likelihoods include observed
frequency covariance, propagated density covariance and any explicitly
declared density/frequency cross covariance; normalization includes the
covariance determinant. Branch likelihoods are combined before multiplying
by the five-factor prior and normalizing over gated candidates. No universal
height window or SEP output enters this calculation.

For a pre-inferred-height constraint, require the density-model checksum and
`joint-density-inference` declaration. The residual combines all inferred
heights with the inferred log reference-density coordinate and uses the
complete joint covariance. It is one shared inversion, not independent copies
for candidates sharing that density authority. This supplied branch needs an
explicit density normalization coordinate; exact front-sampled profiles are
therefore accepted only by the radio branch. Unsupported height inversions
fail rather than inventing a density-independent error.

Construction data cannot qualify a candidate, and withheld SEP cannot gate or
weight it. `load_withheld_after_freeze` verifies a qualified frozen selection
before returning the withheld product for downstream validation. Role labels
and checksum ownership are auditable contracts; scientific independence still
requires reviewed observation provenance.

## Reviewed multi-spacecraft requests

The `field-line-requests` command consumes a construction ephemeris and
explicit observer/line IDs, times, trace controls, solar/outer radii and flux
measures. Every seed is derived from a covered ephemeris sample. An explicit
review authority, frame transformation, epoch, covariance and ephemeris
identity are stored. IDs are sorted deterministically. The neutral request
record is ready for an application adapter to map into the Stage-10 tracing
API; it does not embed coordinates in a test executable or silently activate
an AMPS active corridor.

## Frozen transfer and cross-model comparison

`freeze-transfer` freezes equations, source-efficiency/spectrum/scattering
forms, inference procedures, reference-surface/return policies and background
methodology before withheld SEP is examined. Transferable constants and
permitted event-specific inputs are disjoint. Observer-path records must
declare ICME screening, continuous background/SEP coverage and independent
front reconstruction. An unrepresented prior cloud forces a stress-test
classification; the 2012 May 17 event is a stress test unless its preceding
cloud background is qualified. `bind-transfer` accepts only the preregistered
event inputs, each bound to a content checksum and frozen inference procedure,
and computes a fresh per-run calibration fingerprint.

`compare-mhd` consumes two independently prepared, frozen
`sep-offline-background-samples-v1` records. Use `core.freeze_record` and
`core.canonical` in the offline producer to freeze those records; a hash alone
does not qualify the source model. Each record retains model/version/run,
source location/checksum, magnetogram/preprocessing, magnetic/open-flux
normalizations, epoch/frame, variable units/definitions, composition/EOS,
cadence, support, interpolation, masks and common front/history identity.
Coordinates include time and must match exactly within declared support.
Magnetogram, preprocessing, front/history and source bindings must be actual
64-hex SHA-256 digests. Absolute magnetic/open-flux normalization records each
contain a positive `value`, explicit `units` (`T`/`Wb`) and a physical
`definition`; a shared arbitrary label does not constitute matched normalization.

Every required density, pressure, temperature, vector magnetic/velocity,
Alfvén/fast-speed, signed/unsigned open-flux and same-front D2 metric is
reported with full covariance, masked support and finite differences.
Topology/connectivity agreement is reported where definitions are comparable.
Frame/epoch/support/units/mask/history mismatches fail. Different inner-boundary
or normalization authorities require `--allow-combined-discrepancy` and force
the explicit combined boundary-plus-model label. Residual tuning is forbidden.
A complete comparison claims neither truth validation nor qualification of a
runtime imported-MHD provider.

`EVT3D01` and `XMD3D01` remain reserved event-campaign gates in the
specification. The Stage-12 software prepares/validates their manifests and
tests the contracts synthetically; it does not register future event evidence
as a completed generic release test.

## Running and verification

See [the synthetic walkthrough](../examples/stage12/README.md) for fully
specified source files, a two-member campaign and executable CLI commands.

```sh
make -j8 test-stage11  # cumulative 202 canonical gates
make -j8 test-stage12  # cumulative 206 canonical gates
make -j8 test          # all currently registered shared-model gates
python3 test/run_tests.py --test CAL3D01
python3 test/individual/PROV3D01/test.py
```

`test/test_stage12.py` supplies the four new canonical callbacks. `PROV3D01`
checks synthetic parameter/frame/covariance recovery; `PROV3D02` checks
byte-deterministic reprocessing, version/source identity and CLI publication;
`PROV3D03` checks provenance, frame/unit/epoch/covariance and role-leakage
failures, with synthetic transfer/MHD manifest checks. `CAL3D01` checks the
complete manufactured candidate product, an independently evaluated radio
likelihood and covariance-correct joint height likelihood, and prohibited
selectors/mutations. Individual launchers and cumulative runs share the same
registry. JSON/JUnit evidence is generated only under `build/`.

## Host and production boundary

The shared C++ Stage-11 providers and these offline Stage-12 assets are
independently buildable. Schema 5 and field-line bundle schema 3 are unchanged.
The neutral observation JSON is not an added schema-5 `.ini` selector: a host
adapter must prepare the existing typed background/closure/trace records and
pass their original physical qualification gates. An arbitrary curved sheath,
finite PFSS/SCS transition, general imported-MHD provider and real-event
production calibration are not qualified by the manufactured tests. Full
AMPS/MPI mover integration must be verified in the complete application tree.
