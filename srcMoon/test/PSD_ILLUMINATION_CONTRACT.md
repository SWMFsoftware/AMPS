# M2 PSD illumination production-wiring contract — 2026-11-16

This extension completes the R2.2 source coupling identified during the M2
code audit. It changes illumination dispatch for explicitly enabled M2 only;
it does not introduce a PSD coefficient, spectrum, yield, source population,
observation or permanent U/I ID. The minimal shared-core interface was approved
after the inventory below and is now implemented. Native campaign results,
including regression preservation, must still establish its executed scope;
implementation or registration alone is not a numerical PASS.

## Inventory and source timing

* The production per-element kernel is
  `Exosphere::SourceProcesses::PhotonStimulatedDesorption::GetSurfaceElementProductionRate(int,int,void*)`
  in `Exosphere.h`. Its law is the existing population [atoms] times effective
  cross section [m2] times 1-AU photon flux [m^-2 s^-1], multiplied by
  directional incidence and inverse-square distance.
* `Exosphere_SourceProcesses.cpp` supplies the two-argument callback used by
  the ordinary discrete surface CDF. It delegates to that same kernel.
* `Exosphere_Parallel.cpp::ExchangeSurfaceAreaDensity` calls the kernel and
  updates total `SourceRate`, per-area maxima and per-element diagnostics at
  the end of `PIC::TimeStep`. Those caches alone can be stale at the next
  source interval or initial particle-weight selection.
* `Exosphere.cpp::totalProductionRate` and its ordinary source selector read
  the PSD total; the ordinary CDF is initialized before per-event departure
  scopes. The runtime total must therefore be refreshed from the same
  current per-element rates when the optional geometry hook is installed.
* `Exosphere_FiniteInventory.cpp` derives PSD competing first-order hazard
  from this exact production rate divided by actual population before packet
  departure times are drawn. The already approved finite operator holds its
  coefficients fixed over the particle/source interval. Its departure scope
  is for placement and motion, not a different rate normalization.
* `main_lib.cpp::amps_time_step` records `etBegin`, advances thermal state
  over the interval and changes global orbit state before `PIC::TimeStep`.
  Consequently the PSD adapter must not infer source-coefficient time from
  the updated global end epoch or from a fractional departure scope. The
  application explicitly supplies source-interval-begin TDB, initialized to
  the configured starting epoch before the first source/weight query. Frozen
  mode retains its declared fixed external-state epoch.

## Test specification

| Field | Contract |
|---|---|
| IDs/name | Existing U09 analytical illumination and I19 linked production source/inventory checks; preserve all frozen IDs. |
| Purpose | Apply one production terrain/Earth/direct-photon illumination to PSD and verify default dispatch, actual formula and cached-source normalization. |
| Production capability | Existing PSD kernel, source CDF, total selector and integrated finite inventory hazard. |
| Production entries | Inline PSD per-element/total APIs, ordinary source callbacks, `InjectFiniteSurfaceInventory`, Moon production surface illumination. |
| Configuration | Optional NULL-default geometry-factor hook; registration only after explicitly enabled M2/logical-surface initialization; M0/M1 and other applications leave it NULL. |
| Inputs | Actual logical stocks and immutable source coefficients; explicit actual point/normal/Sun/Earth geometry and source-begin epoch; deterministic analytical local fixtures. |
| Independent reference | Existing law's exact zero/full/half geometric scaling, unchanged legacy expression, summed per-element/total identity, signed atom/kg budgets; no observational comparison. |
| Units | Hook returns dimensionless direct photon flux divided by 1-AU photon flux: jointly visible projected incidence divided by actual point-to-Sun distance in AU squared. Output is atoms s^-1. No thermal Bond-albedo absorption enters PSD. |
| Frame/time | Body-fixed production terrain geometry; explicit source-begin TDB; actual departure geometry/mover scopes remain separate. |
| Deterministic/stochastic controls | No geometry/kernel/CDF normalization RNG draws; original stochastic departure/energy/angle streams remain unchanged in NULL mode. |
| Procedure | Preserve old source expression in NULL branch; install/read production geometry factor only in M2; refresh runtime totals from production rates; test analytical full/half/zero masks through actual PSD kernel and default invariance; run M0 nightly/core regressions and source-budget development. |
| Metrics | Per-element rate, exact scaling, recomputed total/current population, cache state, independent masks versus joint mask, exact nightly diff and unchanged core outcomes. |
| Acceptance | Existing physical law and coefficients unchanged; deterministic geometry/formula checks within declared floating arithmetic tolerances; M0 baseline exact; prior passing regression identities retained; existing I19 accounting/coverage/numerical gates unchanged. |
| Expected PASS | Evaluated local/source wiring and ledger gates; not a Na observational or M4 source-strength validation. |
| Failure modes | Albedo multiplied into PSD; stale source cache; wrong normal/distance; ignored terrain/eclipse; mismatched CDF/hazard epoch; invalid population/factor; unexpected NULL-mode changes; corrupted layout/source evidence. |
| Required data / hashes | None for analytical fixtures; linked Moon build uses declared LOLA/SPICE and source/build/input identities, retained in normal manifests. |
| Predecessors | M0/local geometry for support; same-campaign I13 and original I15/I19 prerequisites for formal qualification. |
| Artifacts | Actual source probes, source/total cache metrics, test logs/results, source/build/executable/data hashes, full code comments/README, nightly zero-byte diff, signed ledgers and preserved attempts. |
| Status semantics | Local/supporting PASS stays scoped. Missing external inputs/predecessors SKIPPED, malformed/invalid runtime geometry ERROR, evaluated assertion miss FAIL. No source-only result becomes formal validation. |

The parent agent approved the minimal interface after this inventory and
preserved the existing dirty core files before editing. The implemented
shared-core API is `PhotonStimulatedDesorption::RuntimePhotonFluxFactor`, an
optional `double(int surfaceElement, void* logicalSphere)` callback, initially
NULL. The application owns an immutable serially prepared source-begin-epoch
cache; the callback only reads it, without SPICE, random draws or stock writes.
`RefreshRuntimeProductionRates` updates `SourceRate`, `maxLocalSourceRate` and
per-area diagnostics from the same element kernel without invoking the
inventory-mutating `ExchangeSurfaceAreaDensity` routine. NULL mode retains
the original element expression and cached-total getter. Negative/nonfinite
physical input, a negative/nonfinite geometric factor, overflow, invalid
element indices or missing allocations terminate with the AMPS diagnostic.
Exactly zero stock/cross section/photon flux skips geometry after validating
the physical quantities; unsupported species never acquire a source.

Existing U09 now includes the compiled `psd-illumination` probe and frozen
analytical fixture contract. Its separately invoked
`psd-illumination-invalid` mode tests the fatal nonfinite-factor diagnostic.
These tests use no observations and no modified species coefficients. See
[U09 README](U09_terrain_shadow/README.md) for commands and scoped meanings.
Native builds, M0 regression and linked source/budget verification remain
the parent campaign's responsibility; an unexecuted assertion is not a PASS.
No conditional or independently multiplied terrain/Earth surrogate is
permitted. The same joint production geometry service supplies direct
irradiance and PSD's unabsorbed direct photon factor.

## Finite-disc and fixed-site geometry

`Moon::SurfaceEnvironment::EvaluateSolarDisc` samples a uniformly bright
apparent disc at explicit radial/azimuth resolution. For one source element,
the dimensionless direct-photon factor is the quadrature sum of the joint
terrain/opaque-Earth visibility mask times each positive incidence cosine,
including the projected-disc-to-solid-angle Jacobian, diluted by the actual
point-to-Sun distance in AU squared. The source then evaluates

```text
rate_e [atoms/s] = population_e [atoms] * cross_section [m2]
                  * photon_flux_1AU [1/(m2 s)] * geometric_factor_e
```

The thermal entry point additionally integrates directional absorptivity
`1 - BondAlbedo(local incidence angle)` ray by ray before advancing the
existing column solver. PSD does not include that absorptivity. Multiplying
separately averaged terrain and Earth fractions is incorrect when their masks
overlap; replacing incident photon flux by any-disc occurrence or disc area
also changes the source law. Both source and thermal diagnostics retain those
different quantities explicitly.

The optional `.in` selection is `illumination_terrain_solar_model=uniform-disc`
with mandatory positive `illumination_disc_radial_points` and azimuth points
at least four and divisible by four. Legacy point-Sun receipts retain their
declared mode; point-Sun mode refuses unused finite-disc grid controls. A ray
grid is a numerical approximation and cannot establish the absence of an
arbitrarily small exposed sliver without qualification. No grid value is a
substitute for missing D02 product timing, observer height or temporal sampling.

`ApparentSolarSystemGeometryAtPoint` calls NAIF `spkcpo_c` for the actual fixed
observer in MOON_ME_DE421 relative to MOON, with observer-epoch output axes and
explicit `NONE`, `LT`, `LT+S`, `CN` or `CN+S` reception corrections. Observer
metres are converted to SPICE kilometres and both returned target vectors are
converted back to SI metres. This includes the observer's rotating-site
geometry for light-time and stellar corrections; correcting a body-centre
vector and then subtracting a site position is not equivalent. Failure leaves
both outputs unchanged and reports the actual SPICE error; it does not select
an uncorrected fallback. Calls belong to the serialized application SPICE
phase, not particle worker threads. The helper preserves the caller's SPICE
error action, reports a pre-existing error without erasing it, and restores
the action after its own calls.

## SOURCE-BEGIN cache and baseline boundary

Initialization converts the explicit `simulation_start_utc` to TDB before the
first particle-weight/source query. `BeginProductionParticleStep(etBegin)`
refreshes `RefreshProductionPhotonFluxCache` serially for that source-BEGIN
epoch. The ordinary PSD element kernel, total refresh, surface CDF and
finite-inventory competing hazards subsequently read the same immutable
factors through `ProductionPhotonFluxFactor`. No callback performs SPICE,
draws a random number, modifies stock or changes the source epoch.

The cache is bound to the actual logical-surface pointer, element count and
source-BEGIN epoch. Missing/stale cache, wrong surface, invalid index,
nonfinite geometry or overflow causes the existing detailed AMPS fatal
diagnostic; no unit factor or stale value is silently substituted. Frozen
environment mode retains its declared external-state epoch. Evolving mode
supplies the beginning of each interval even after the global orbital state
advances. Fractional packet departure scopes remain placement/mover geometry,
not a second source-rate clock. Finite-stock coefficients remain constant
within the approved split interval, requiring the unchanged time/weight
convergence gates for formal I19.

The shared-core callback starts NULL and is registered only after explicit
M2 surface initialization. NULL keeps the original element formula and cached
total behavior, with no geometry work or extra random draws. This preserves
the M0/M1 and other-application baseline modes; their actual regression result
must still be checked in the retained campaign. Exactly zero physical stock,
cross section or photon flux returns zero after physical-input validation,
without tracing empty elements. Unsupported species do not acquire a source.

## Execution and evidence

Use the interpreter selected by this machine's profile, `/usr/bin/python3`.
Check `numpy`, `scipy`, `matplotlib` and `mpmath` in that interpreter before a
campaign; packages installed only in another shell-selected environment do
not establish its readiness. The ordinary complete M2 development command in
[RUN.md](../RUN.md) includes native controls and the real operator exporters.
For a dedicated matching current generated build, parser/refusal checks are:

```sh
/usr/bin/python3 -B srcMoon/test/test_m2_native_controls.py \
  --build-dir build --input srcMoon/examples/m2_development.in \
  --output-dir test_output/srcMoon/M2-native-controls-NEW
```

The runner links a thin adapter against actual production objects. Its
negative PSD case succeeds only when the real core rejects the nonfinite
factor with its declared diagnostic; an arbitrary crash is not accepted.
It preserves each case's input, command, streams, expected/actual code and
result. Pure operator-runner control checks are separately executable:

```sh
/usr/bin/python3 -B srcMoon/test/test_m2_operator_checks.py
```

Those temporary schema/reporting fixtures do not launch MPI, solve a thermal
column or supply observations. The actual `run_m2_operator_checks.py` exporter
run remains independently necessary and retains actual native CSV, SPICE/data
identities, receipts and resource logs. A deliberately millisecond thermal
replay checks the production state exporter; it does not qualify a synodic
initial condition or annual polar temperature. Current execution outcomes
belong in STATUS.md and the campaign report, not this predeclared contract.
M1 remains unresolved, D02 remains PARTIAL and formal M2 is NOT QUALIFIED.
