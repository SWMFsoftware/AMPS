# SWCME-to-SEP source interface and campaign automation

## Purpose

This remediation closes the boundary between the validated SWCME background /
shock model and an energetic-particle transport consumer.  Earlier revisions
exposed a local `ShockAccelerationState`, but deliberately stopped short of
assigning transport-facing spectrum units, surface-patch bookkeeping, or an
observer-cobpoint record.  Validation also required manually invoking the C++
test executable one test at a time.

The present change adds two layers without duplicating production physics:

1. `swcme_sep_source.hpp` defines the dimensional, serializable SEP source
   contract shared by 1-D and 3-D.
2. `swcme_sep_interface.hpp` adapts prepared SWCME states to background/source
   queries suitable for AMPS or another focused-transport solver.
3. `test/run_tests.py` orchestrates the single C++ validation executable,
   reproducible manifests, parameter sweeps, convergence checks, event/reference
   comparison and report generation.

`tools/sep_reference.cpp` is an intentionally small example consumer of the same
public interface.  It does not contain independent Parker, shock, connectivity,
or source physics.

## Source-state contract

`swcme::sep::SEPSourceState` contains all information needed to decide whether a
source exists and to construct the configured spectrum:

```text
status
active
connection_evaluated / connected
acceleration mode
time_s / source_id
position_m[3] / normal[3]
patch_area_m2 / active_surface_area_m2 / area_fraction
relative_patch_weight
compression / theta_Bn_rad / fast_mach / normal_speed_m_s
upstream_density_m3 / upstream_B_T
upstream_pressure_Pa
focusing_length_m / field_line_path_length_m
q_phase_space
momentum_intensity_index
nonrel_energy_intensity_index
SpectrumConfig
relative_source_weight_per_area
```

There are no ambiguous zero-valued sentinels for connectivity.  `NO_CONNECTION`
and `SOURCE_INACTIVE` are explicit expected `ModelStatus` outcomes; numerical,
configuration and shock-solver errors remain failures.

### Background adapter units

`BackgroundState` uses SI throughout:

```text
density                  m^-3
pressure                 Pa
velocity                 m s^-1
magnetic field           T
magnetic-field magnitude T
div(V)                   s^-1
focusing length          m
position                 m
```

The 3-D single-point adapter calls the production checked Cartesian evaluator
with `N=1`.  The 1-D adapter similarly calls the production checked radial
path.  No new density, Parker, region, or divergence equation exists in the SEP
adapter.

`BackgroundState` also copies the prepared model identity and configuration
digest.  Parker path length and focusing are evaluated in the common
`swcme_solarwind.hpp` production component.  Directional/surface sources carry
local focusing and an unavailable (`NaN`, serialized `NA`) observer path;
connected cobpoint sources carry the selected production connectivity length.

## Spectrum convention

The SOURCE representation supplies the test-particle DSA phase-space slope

```text
f(p) proportional to p^(-q)
q = 3 r_c / (r_c - 1).
```

For an isotropic population, `J(E)=p^2 f(p)` because `dE/dp=v`.  The adapter
therefore uses the exact relativistic normalized shape

```text
J(E)/J(E_ref) = [p(E)/p(E_ref)]^(2-q),
p c = sqrt[K (K + 2 m c^2)].
```

This is preferable to hard-wiring a non-relativistic energy power law across a
source energy interval that may extend into relativistic energies.  The
non-relativistic energy index `(q-2)/2` is retained only as a diagnostic.

The public kinetic-energy interval and reference energy use MeV.  Momentum is
computed in SI.  `rigidity_GV_from_kinetic_MeV()` provides an explicit rigidity
diagnostic, and differential-intensity conversion helpers map between

```text
particles cm^-2 s^-1 sr^-1 MeV^-1
```

and

```text
particles m^-2 s^-1 sr^-1 J^-1.
```

## Normalization modes

`SpectrumConfig::normalization` has two explicit choices.

### `RelativeOnly`

This is the controlled/default mode.  The source shape and spatial patch
weights are dimensionless.  SWCME makes no assertion about an absolute
injection rate or measured flux normalization.  This is the appropriate mode
for isolating connectivity, perpendicular diffusion, and transport effects.

### `ReferenceDifferentialIntensity`

This mode requires a positive `reference_differential_intensity_SI` at
`reference_energy_MeV`.  The returned physical spectrum is

```text
J(E) = J(E_ref) [p(E)/p(E_ref)]^(2-q)
```

with the declared SI units.  Validation rejects an absent/non-positive physical
reference normalization.

The explicit enum prevents a dimensionless relative source from being silently
published as a physical differential intensity.

## Spatial source representation

### 1-D

`Interface1D::source_at_shock()` obtains the already validated 1-D
`ShockAccelerationState` and converts it through the common
`make_source_state()` helper.  No 1-D-specific spectrum algebra exists.

### 3-D directional source

`Interface3D::source_at_direction()` calls
`Model::shock_acceleration_state_checked()`.  This checked wrapper uses the
canonical surface-owned local shock state introduced by the earlier shock-state
remediation.  A finite surface outside the SSE cap is distinguished from a
solver failure.

### Shock-surface source

`Interface3D::build_shock_surface_source()` builds the corrected triangular
shock mesh and its physical triangle metrics.  Each cell receives one source
record sampled at the triangle-centroid direction through the same directional
source path.  For active fast-shock cells,

```text
area_fraction_i = A_i / sum_active A_j
relative_patch_weight_i = source_weight_per_area * area_fraction_i.
```

Sub-fast cells remain in the returned surface for auditing but have
`active=false` and zero source weight.  Consequently a uniform source is
weighted by actual surface area rather than by the number or indexing of mesh
cells.

### Observer cobpoint

`Interface3D::source_at_observer_cobpoint()` first calls the production
`observer_connectivity()` solver.  If no Parker-field-line/shock intersection
exists, the adapter returns `NO_CONNECTION`.  If a root exists, its selected
direction is routed through the same checked directional source path.  There is
no second cobpoint equation in the adapter.

## Source serialization and manifests

`source_csv_header()` and `serialize_source_csv()` produce a deterministic audit
record.  Non-applicable quantities such as an inactive DSA slope are serialized
as `NA`, not as zero.

`Interface1D::resolved_manifest()` and `Interface3D::resolved_manifest()` append
the complete SEP spectrum/normalization contract to the existing deterministic
SWCME resolved-configuration manifest.  Campaign manifests archive this text
and its SHA-256.

## Prepared-state and parallel-use convention

The intended transport call sequence is:

1. construct and validate one `Interface1D` or `Interface3D` from immutable run
   parameters;
2. call `prepare(time_s)` once per SWCME background update;
3. distribute/read the resulting prepared step for particle/background/source
   queries during that interval;
4. prepare a new step only when the background update time advances.

The interface itself owns no MPI communicator, OpenMP region, random-number
state, or mutable global source cache.  Shock-surface source sampling likewise
does not own an RNG.  AMPS remains responsible for its existing MPI/thread
ownership and may use `area_fraction` / `relative_patch_weight` with its native
parallel random-number machinery.

## Campaign manager

`test/run_tests.py` is the only higher-level validation manager.  It always
executes tests through `test/output/test_swcme`, which keeps the C++ registry as
the authoritative list of deterministic/statistical validations.

### Profiles

```text
SMOKE    short development gate
ROUTINE  broad deterministic gate; excludes MSH05 and CON05
FULL     every registered C++ validation
EVENT    FULL + configured event analysis
```

Profiles are plain version-controlled files under `test/profiles/` and are
validated against the live C++ `--list` output before execution.

### Campaign artifacts

Every run writes:

```text
manifest.json
summary.json
summary.csv
logs/<test>.log
```

FULL/EVENT can additionally write `sep_reference.csv`.  Event analysis may add
`sweeps/`, `convergence/`, `comparisons/`, and `plots/` artifacts.

The manifest records:

- UTC creation time and selected profile/tests;
- random seed;
- git commit, branch and dirty flag;
- compiler command/path/version/flags;
- host, platform, machine, Python version and CPU count;
- visible MPI and OpenMP environment variables;
- deterministic source-tree SHA-256;
- complete resolved SWCME+SEP parameters, their SHA-256 and provenance;
- event JSON plus SHA-256 when EVENT is used.

### Parameter sweeps

A sweep contains a JSON parameter dictionary and command-template array.  The
runner expands the Cartesian product and substitutes parameter names plus
`{seed}`, `{root}`, `{test_dir}` and `{output_dir}`.  Optional `metric_regex`
extracts one scalar from stdout and optional min/max bounds make it an
acceptance criterion.  Every case receives its own log and one combined CSV is
written.

`sep_reference --probe-time-hours` deliberately emits simple `KEY=value` data
so event sweeps can capture quantities such as `SHOCK_RADIUS_AU`, `FAST_MACH`,
or `COMPRESSION` without parsing a full CSV history.

### Convergence

A convergence analysis accepts explicit positive `(x,error)` arrays or data
from a named sweep.  It performs a least-squares fit of

```text
log(error) = p log(x) + constant
```

and reports the observed order `p`, with optional minimum/maximum acceptance
bounds.

### Event/reference comparison

CSV comparison supports row-aligned data or key matching, with optional numeric
key tolerance.  Each configured data column has independent absolute and/or
relative tolerances.  Failures are preserved in a machine-readable JSON report.

### Plotting

Plots are intentionally post-processing only.  They are made from campaign CSV
data using Matplotlib when available.  A missing Matplotlib installation causes
an optional plot to SKIP rather than changing a physics result; `required=true`
turns it into a campaign failure.

### Exit codes

```text
0 = all required work passed
1 = validation/reference/event failure
2 = CLI or input-configuration error
3 = build failure
```

These codes are stable so shell/CI/HPC batch scripts do not have to parse text
output.

## Reference exporter

`tools/sep_reference.cpp` demonstrates the public 3-D SEP adapter.  It supports:

```text
--print-manifest
--output FILE
--start-hours / --end-hours / --step-hours
--energies E1,E2,...
--probe-time-hours
--v0-kms
--gamma-km-inv
--half-width-deg
--observer-lon-deg / --observer-lat-deg
```

Time-history output is CSV.  Probe output is exclusively `KEY=value` text so
campaign metric extraction has an unambiguous wire format.

## Validation

Six new C++ tests protect the integration contract:

- `SEP01`: AMPS-facing background adapter equals the direct production query;
- `SEP02`: spectrum units, relativistic rigidity and normalization;
- `SEP03`: byte-identical complete 1-D/3-D source records;
- `SEP04`: finite-SSE source-surface active-area weighting;
- `SEP05`: observer source reuses the production cobpoint;
- `SEP06`: resolved compression exposes no prescribed SEP source and manifest
  serialization is deterministic.

The campaign manager has separate Python regression tests for profile expansion,
convergence fits, parameter sweeps, CSV/event comparison, manifest/report
creation, and deterministic command/configuration errors.

This remediation does not claim an empirical injection-efficiency law.  The
source contract deliberately separates validated shock/geometry quantities from
optional absolute spectral normalization so future event calibration can be
added without changing the transport-facing state semantics.
