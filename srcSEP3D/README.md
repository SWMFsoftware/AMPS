# srcSEP3D

## V01–V05 controlled extensions and evidence gates

See `V01_V05_IMPLEMENTATION.md` for the algorithms, configuration, tests, and
release rules. V01 adds controlled perpendicular diffusion and relativistic
gradient-B/curvature drift. V02 uses distinct compiled 1-D and 3-D producers;
V03–V05 add native MPI, scientific-campaign, and release-evidence gates. R8 is
deferred, so live-SWMF validation remains explicitly blocked.

`srcSEP3D` is the AMPS application for three-dimensional solar-energetic-
particle and energetic-electron transport in the heliosphere. Its target
physics is the Parker or focused transport equation in an analytic Parker or
coupled SWMF/AWSoM background, with Alfvén-turbulence scattering and SWCME
shock/source parameters.

## Complete standalone initialization input (schema 3)

New standalone shock-injection runs use the strict schema-3 deck:

```sh
./amps --input srcSEP3D/examples/sep3d_analytic_parker.in --dry-run
./amps --input srcSEP3D/examples/sep3d_analytic_parker.in
```

The first command parses, canonically resolves, validates, fingerprints, and
resource-preflights the request without allocating AMPS. The production command
repeats the same parser on initialization; no cached parse result or generated
configuration source is used. Unknown/duplicate sections or keys, missing
fields, invalid enumerations, non-finite values, and inconsistent physical
descriptions fail before mesh allocation.

The file is INI syntax. Section/key names are case-insensitive, `#` begins a
comment, and each assignment is `key = value`. Application dimensional keys
are bare SI values with the unit in the key. Values inside `[swcme]` retain
their unit token and are resolved by the canonical model-owned parser. Schema
3 requires even mode-inactive application fields so a later mode edit cannot
silently revive a C++ default. The only deliberately excluded model values are
deprecated compatibility parameters that no longer affect SWCME physics.

### Mesh, line, run, and resource syntax

| Section | Required keys | Contract |
|---|---|---|
| `[run]` | `schema_version`, `intent`, `transport`, `time_step_s`, `maximum_time_steps`, `campaign_seed`, `background_cadence_steps`, `injection_cadence_steps` | Schema is `3`; intent is `shock-injection`; transport is `parker3d` or `focused3d`. Time step is positive SI seconds, seed and step counts are nonzero, and injection cadence must be exactly one because `samples_per_step` is a count for every simulation step. |
| `[domain]` | `preset`, `inner_radius_m`, `inner_boundary`, `outer_radius_mode`, `outer_radius_m`, `outer_boundary`, `coordinate_frame`, `origin_x_m`, `origin_y_m`, `origin_z_m` | Presets are `solar`, `one-au`, or `mars`; outer mode is `preset` or `explicit`. The current heliocentric implementation requires the declared origin `(0,0,0)`, absorbing inner boundary, and a domain containing fixed observers and mesh references. |
| `[parker_spiral]` | `origin_x_m`, `origin_y_m`, `origin_z_m`, `initial_x_m`, `initial_y_m`, `initial_z_m`, `length_m`, `point_count` | Finite diagnostic/refinement centreline. The origin equals the domain origin; the initial point lies on the inner sphere and matches the tube longitude/colatitude. Length is positive arc length and count includes both endpoints. |
| `[mesh]` | `global_cell_size_m`, `minimum_cell_size_m`, `cells_per_block_edge`, `maximum_level`, `memory_budget_bytes`, `block_overhead_bytes` | Global/floor resolution, AMPS block shape, realizable AMR depth, and pre-allocation resource ceiling. |
| `[mesh.solar]` | `enabled`, `surface_cell_size_m`, `transition_outer_radius_m`, `profile`, `exponent` | Resolution at the Sun and its radial degradation to the global value. Profiles are `linear`, `power-law`, or `smoothstep`; exponent is positive. |
| `[mesh.tube]` | `enabled`, `source_longitude_rad`, `source_colatitude_rad`, `reference_radius_m`, `radius_at_reference_m`, `radius_mode`, `center_cell_size_m`, `transverse_profile`, `transverse_exponent` | Parker-centreline location, physical/angular tube radius, centre resolution, and degradation in the perpendicular plane. Radius mode is `physical-constant` or `constant-angular-width`; profile choices match `[mesh.solar]`. |
| `[memory]` | `base_cell_bytes`, `base_node_bytes`, `block_structure_bytes`, `communication_bytes_per_block`, `particle_bytes`, `particles_per_cell`, `halo_fraction`, `safety_margin_fraction` | Explicit build-dependent coefficients used by the allocation-free memory preflight; fractions are finite and nonnegative. |

The near-Sun interpolation is

\[
h_\mathrm{sun}=h_\mathrm{surface}+P(s)
(h_\mathrm{global}-h_\mathrm{surface}),\qquad
s=\frac{r-r_\mathrm{in}}{r_\mathrm{transition}-r_\mathrm{in}},
\]

with `s` clamped to `[0,1]`. The tube uses the same profile with normalized
transverse arc distance. In constant-angular-width mode its radius is
`radius_at_reference_m*r/reference_radius_m`. The finer of radial and tube
requests wins, then the result is clamped between the declared minimum and
global cell sizes. The same function drives dry-run planning, the standalone
octree tests, and AMPS `localResolution()`.

### Background, transport, source, and particle syntax

| Section | Required keys | Contract |
|---|---|---|
| `[background]` | `provider`, `external_script` | Standalone schema 3 requires `analytic-parker` and `external_script = false`. A coupled SWMF host uses the parser-free typed interface instead. |
| `[background.parker]` | `reference_radius_m`, `radial_field_at_reference_t`, `solar_rotation_rate_rad_per_s`, `solar_wind_speed_m_per_s`, `magnetic_polarity`, `number_density_at_reference_m3`, `temperature_k`, `validity_cadence_s` | Complete analytic Parker state in SI. Wind, rotation, source radius, field normalization, density, temperature, and polarity must agree with `[swcme]`. |
| `[turbulence]` | `authority`, `delta_b_over_b`, `k_min_per_m`, `k_max_per_m`, `spectral_index`, `correlation_length_m`, `missing_data`, `resonance_range`, `self_consistent_3d` | Schema 3 standalone authority is `prescribed`; self-consistent 3-D is false. Missing-data policy is `fail` or `ballistic`; resonance policy is `reject` or `power-law-extension`. |
| `[transport]` | `cell_crossing_fraction`, `diffusion_fraction`, `focusing_fraction`, `cooling_fraction`, `field_variation_fraction`, `shock_crossing_fraction`, `minimum_substep_s`, `maximum_substeps`, `pitch_angle_scheme`, `perpendicular_diffusion`, `constant_kappa_perpendicular_m2_per_s`, `kappa_perpendicular_to_parallel_ratio`, `drifts` | Positive timestep limiters. Pitch scheme is `reflecting-milstein` or `reflecting-euler-maruyama`. Perpendicular mode is `none`, `constant`, or `constant-ratio`; drift is `none`, `gradient-b`, `curvature`, or `gradient-curvature`. Selected extensions require their positive coefficient/storage. |
| `[shock]` | `authority` | Must be `swcme`. Schema 3 rejects the retired constant-radius/speed/compression surrogate fields. |
| `[source]` | `enabled`, `physical_particle_rate_per_s`, `injection_efficiency`, `minimum_energy_j`, `maximum_energy_j`, `samples_per_step` | Source is enabled. Rate is the physical seed rate before efficiency and patch partition; energies, efficiency, and exact global computational count are positive. No free spectral index is accepted: each patch's canonical compression ratio determines its DSA slope. |
| `[species]` | `amps_index`, `name`, `mass_kg`, `charge_c`, `macroparticle_weight` | Current reviewed contract is exactly one index-zero proton. Mass and charge must match both AMPS and `[swcme]`; base weight is derived by the conservation identity below. |
| `[storage]` | `magnetic_gradient`, `velocity_gradient`, `sampling_bytes_per_cell` | Explicit associated-data layout. Required transport choices may force a gradient on before the layout fingerprint freezes. |

For schema 3 the base AMPS weight must satisfy

\[
W_0=\frac{\dot N_\mathrm{seed}\,\epsilon\,
w_\mathrm{surface}\,\Delta t}{N_\mathrm{macro}},
\]

where `physical_particle_rate_per_s` is \(\dot N_\mathrm{seed}\), the application
and canonical SWCME injection efficiencies are identical, canonical
`shock.relative_source_weight_per_area` is \(w_\mathrm{surface}\),
`run.time_step_s` is \(\Delta t\), and `source.samples_per_step` is
\(N_\mathrm{macro}\). Injection cadence is one, so no hidden cadence multiplier
exists. The parser rejects a mismatched `species.macroparticle_weight`.

At an active boundary, the runtime assigns exactly `samples_per_step` over all
active SWCME surface patches. It first reserves one weighted representative for
each nonzero patch, then apportions the remainder by physical patch weight with
a deterministic largest-remainder rule and stable source-ID tie-break. If the
requested count is smaller than the active patch count, initialization fails;
it never drops a physical patch. Each patch's individual AMPS weight correction
is its represented physical population divided by its exact assigned count, so
integer apportionment does not alter total physical number.

### Observer and output syntax

Every `[observer.ID]` requires all of these fields:

| Fields | Meaning |
|---|---|
| `kind`, `normalization` | Kind is `fixed-cartesian`, `fixed-heliographic`, `moving-cartesian`, `spherical-shell`, or `field-connected`; normalization is `represented-particles` or `differential-intensity`. |
| `position_x_m`, `position_y_m`, `position_z_m`, `follows_trajectory` | Initial observer location and explicit trajectory choice. A fixed observer must lie inside the heliocentric shell. |
| `velocity_x_m_per_s`, `velocity_y_m_per_s`, `velocity_z_m_per_s` | Cartesian trajectory velocity; explicit even for a fixed observer. |
| `collection_radius_m`, `shell_radius_m` | Positive spatial acceptance geometry. |
| `cadence_s`, `energy_bins`, `pitch_angle_bins` | Positive cadence/counts; cadence must be an integer multiple of the simulation step. |
| `minimum_energy_j`, `maximum_energy_j`, `minimum_mu`, `maximum_mu` | Energy and pitch-cosine acceptance with ordered bounds and `-1 <= mu <= 1`. |
| `species`, `products` | Comma-separated AMPS species indices (currently only `0`) and requested product names. |

`[output]` requires `cadence_steps`, `checkpoint_cadence_steps`, `directory`,
`prefix`, `initialization_mesh_tecplot_file`, and
`initialization_parker_line_tecplot_file`. `[restart]` requires `input_path` and
`output_path`; the literal `none` selects a fresh run. Output/restart relocation
does not change the physics fingerprint, while every mesh, observer, source,
species, and canonical SWCME input does.

After final AMR refinement and decomposition, the initialization stage writes
the actual distributed AMPS tree to `initialization_mesh_tecplot_file` and a
rank-zero ordered Parker-centreline zone to
`initialization_parker_line_tecplot_file`. The line contains arc length,
Cartesian position, heliocentric radius, and requested resolution, all with SI
unit-bearing variable names. Both paths are mandatory and write failures are
fatal.

### Complete canonical `[swcme]` syntax

The common model keys are the same unit-aware fields used by srcSEP:

| Group | Required keys |
|---|---|
| Scenario/ambient | `preset`; `ambient.wind_speed`, `ambient.density_1au`, `ambient.magnetic_field_1au`, `ambient.proton_temperature`, `ambient.adiabatic_index`, `ambient.alpha_to_proton_ratio`, `ambient.electron_temperature`, `ambient.alpha_temperature`, `ambient.thermodynamic_closure` |
| Parker/CME | `parker.radial_polarity`, `parker.sin_theta`, `parker.source_radius`, `parker.solar_rotation_rate_rad_per_s`; `cme.kinematics`, `cme.launch_radius`, `cme.launch_speed`, `cme.drag_coefficient`, `cme.extrapolation` |
| Shock/regions | `shock.region_mode`, `shock.acceleration_mode`, `shock.relative_source_weight_per_area`; `geometry.sheath_thickness_1au`, `geometry.ejecta_thickness_1au`; all three `smoothing.*_width_1au` keys; `sheath.ramp_power`, `sheath.leading_edge_speed_factor`, `ejecta.density_factor`, `ejecta.speed_factor` |
| Event/source | `event.launch_epoch`, `event.valid_from`, `event.valid_until`; `source.particle_mass`, `source.charge_number`, `source.energy_min`, `source.energy_max`, `source.reference_energy`, `source.injection_efficiency`, `source.normalization`, `source.reference_intensity_si` |

The genuinely three-dimensional keys are all required:

```ini
geometry.shape = sphere
geometry.axis_ratio_y = 1
geometry.axis_ratio_z = 1
geometry.half_width_rad = 1.5707963267948966
geometry.cme_direction_x = 1
geometry.cme_direction_y = 0
geometry.cme_direction_z = 0
geometry.solar_rotation_axis_x = 0
geometry.solar_rotation_axis_y = 0
geometry.solar_rotation_axis_z = 1
surface.theta_intervals = 13
surface.phi_points = 24
```

`cme.kinematics = data_driven` additionally requires both `cme.data_times` and
`cme.data_radii`; the pair is forbidden otherwise. The canonical resolver
validates all units and model combinations, then emits a normalized manifest
and fingerprint. The standalone AMPS crossing operator is currently spherical,
so schema 3 deliberately accepts only `geometry.shape = sphere`,
`shock.region_mode = shock_only`, and `shock.acceleration_mode = source`.
Ellipsoid or SSE input is rejected rather than approximated by a sphere.
`source.normalization` must be `relative_only` because the application-owned
physical rate supplies absolute population normalization.

The resolved SWCME and application descriptions must also agree on Parker wind,
rotation, source radius, reference-latitude field normalization, +Z rotation
axis, density, temperature, polarity, particle mass/charge, source energies,
and injection efficiency. Before AMPS allocates its mesh, the provider evaluates
the complete surface at `event.valid_from`, solves every canonical MHD shock
state, validates every source patch, and proves the requested exact particle
count can represent all active patches. A delayed valid epoch is allowed and
the provider remains inactive before it; an invalid surface or empty source is
not.

Here SWCME's `ambient.magnetic_field_1au` is total magnitude at one AU and
`parker.sin_theta`, whereas the analytic provider accepts radial field at the
application's arbitrary `background.parker.reference_radius_m = r_ref`. The
required value is therefore

\[
B_r(r_\mathrm{ref})=
\frac{|B|(1\,\mathrm{AU})}
{\sqrt{1+[\Omega(1\,\mathrm{AU}-r_0)\sin\theta/V_\mathrm{sw}]^2}}
\left(\frac{1\,\mathrm{AU}}{r_\mathrm{ref}}\right)^2.
\]

The polarity is checked separately. This conversion permits any physically
valid declared reference radius; it does not silently treat a total field at
one AU as a radial field somewhere else.

The authoritative, fully commented input is
[`examples/sep3d_analytic_parker.in`](examples/sep3d_analytic_parker.in).
Schemas 1 and 2 remain available only for existing transport campaigns; schema
3 is the no-hidden-default initialization path described here.

## B01 source-distribution baseline

The 3-D application is independently classified by
[`SOURCE_MANIFEST.json`](SOURCE_MANIFEST.json). Production code, standalone
tests, documentation, retained validation inputs, and generated paths are
declared separately. No manifest or runner inside `srcSEP3D` inspects
`srcSEP`; the AMPS-level `tools/sep_package_hygiene.py` release gate checks the
two peer applications and their shared model directories without creating a
runtime or build dependency between them.

Release archives contain no standalone `test/stage1` binary, object/archive
files, Python cache, or `test_output` evidence. These are regenerated by the
documented runner after extraction. The local `.gitignore` is deliberately
narrow so frozen kernel records and reviewed validation manifests cannot be
hidden accidentally.

## Implemented scope

The production tree implements the rebaseline and shared foundations (R0–R2),
**Phase M Mesh and Storage**, **Phase B Background Providers and Snapshots**,
**Phase T Turbulence and Scattering Inputs**, **Phase P Transport Cores**,
**Phase A AMPS Mover and Source Adapters**, and **Phase O Sampling, Output, and
Restart**, **Phase V Integration and Scientific Validation**, and production
runtime improvements **R01–R07**.

`amps_time_step()` now enters the typed Runtime particle phase, calls the AMPS
step through the installed srcSEP3D mover, closes a global conservation ledger,
and executes due snapshot, source, observer, and checkpoint transactions at the
joined boundary. `make prepare-production` installs the configured mover hook;
`amps_init()` installs the immutable local resolver automatically. No legacy
mover, source, sampler, or fallback physics is substituted.

### R0–R2 foundation

- The obsolete `SEP3D.cpp`, axisymmetric mover, legacy global sampler, wedge
  mesh, Maxwellian prepopulation, and placeholder output paths are absent.
- `src/models/sep_common/sep_common.a` is the single owner of general SEP
  transport/coefficient kernels used by both SEP applications.
- `src/models/swcme/swcme.a` is the single compiled SWCME implementation.
- srcSEP3D never inspects or requires the independent `srcSEP` application
  directory.
- `RunConfiguration3D` is immutable after validation and fingerprints every
  trajectory-relevant option. Output-only options remain outside that physics
  fingerprint.
- `Runtime` owns the lifecycle, active background generation, output cadence,
  checkpoint sequence, and restart counters.
- The standalone and SWMF adapters enter the same `Runtime` acquisition and
  publication path.
- AMPS mover integers are translated only in `amps/amps_mover_status.h`, where
  `static_assert` binds the adapter to the actual `pic.h` ABI.

### Improvements R01–R07: production runtime closure

- R01 installs one idempotent generated AMPS mover hook, one immutable mover
  context, and collectively closed per-step/species particle ledgers.
- R02 consumes the full requested AMPS interval through accepted substeps,
  re-resolving cells and coefficients before every substep.
- R03 stages background and turbulence together and publishes their immutable
  generation only after every MPI rank accepts the complete candidate.
- R04 derives physical time and all recurring events from one persisted integer
  tick and verifies PIC, snapshot, shock, and Runtime clocks before motion.
- R05 reconstructs analytic or coupled shock states, converts physical patch
  rates to weighted AMPS particles with semantic stochastic rounding, and
  records number/energy/momentum source ledgers.
- R06 makes observers first-class configuration, gathers complete records in
  stable-ID order, reports statistical uncertainty, and resets a sampling
  window only after atomic publication succeeds.
- R07 writes and transactionally restores a versioned complete state image,
  including identities/layout, clocks/events, provider generations, shock,
  particles/RNG tuples, ledgers, and pending sampling state.

See [RUNTIME_INTEGRATION.md](RUNTIME_INTEGRATION.md) for the algorithms,
failure semantics, MPI ownership rules, and restart ordering.

### Improvements C01–C05: production configuration and preflight

- The standalone host reads one versioned, sectioned input file through
  `runtime/configuration_io.cpp`. Every dimensional key declares its SI unit;
  duplicate, unknown, unitless, malformed, or missing required entries fail
  before AMPS initializes.
- `RunConfiguration3DOptions` is the single typed contract for file-driven and
  coupled SWMF construction. It contains complete Parker, turbulence,
  transport, shock, source, species, observer, mesh, storage, memory, output,
  and restart groups.
- Domain radii use explicit `preset` or `explicit` modes rather than a numeric
  sentinel. Solar, one-AU, and Mars presets are normalized before bounds,
  observers, shock extent, mesh levels, and fingerprints are validated.
- `core/parker_geometry.cpp` owns one polarity-independent Parker curve and
  tangent for both background evaluation and tube refinement. Magnetic
  polarity reverses the field only; it cannot move the refined tube.
- Composite near-Sun and transverse-tube profiles meet the global resolution
  continuously. Their overlap requests the finer size. The tube radius is
  defined at an explicit reference distance and may retain physical width or
  constant angular width.
- `--dry-run` reports normalized physics identity, resolution extrema,
  estimated blocks by level, and full resident/particle/halo/sampling/safety
  memory without allocating the AMPS mesh.

For a standalone preflight and run:

```bash
./amps --input srcSEP3D/examples/sep3d_analytic_parker.in --dry-run
./amps --input srcSEP3D/examples/sep3d_analytic_parker.in \
  --output-dir output/production
```

See [CONFIGURATION.md](CONFIGURATION.md) for the schema, algorithms, formulas,
coupled-host contract, and C01–C05 acceptance evidence.

### Phase M: mesh and storage

- Earth and Mars heliospheric domain presets are represented as Cartesian
  cubes containing a physical inner sphere and requested outer sphere.
- One resolution law provides named near-Sun degradation and an optional
  finite-width Parker-spiral tube with a continuous transverse profile.
- The same AMPS-independent law drives the standalone octree verifier and the
  production `localResolution()` callback.
- The standalone octree enforces 2:1 face balance, deterministic global leaf
  and cell IDs, owner assignment, leaf histograms, and a pre-allocation memory
  estimate.
- The complete static/sampling byte layout is frozen in `RunConfiguration3D`
  before AMPS initializes its cell buffer. AMPS requests those exact bytes
  through its model allocation callbacks.
- Production cell population iterates only owner-local blocks. The physical
  cells inside the inner boundary remain allocated and zero-initialized; they
  are not presented as valid heliospheric background cells.
- Least-squares scalar/vector gradients support mixed coarse/fine neighbor
  distances and reject rank-deficient stencils.

See [MESH_STORAGE.md](MESH_STORAGE.md) for formulas, storage order, production
initialization order, and invariants.

### Phase B: backgrounds and immutable snapshots

- `BackgroundProvider` is active and independent of AMPS/MPI.
- `AnalyticParkerProvider` supplies Cartesian magnetic field, analytic
  magnetic/velocity gradients, focusing length, curvature, density,
  temperature, pressure, Alfvén speed, and validity metadata. Its Cartesian
  form has finite polar limits.
- `SwmfAwsomProvider` imports read-only SI or documented AWSoM coupling units,
  validates frame/epoch/ownership/completeness, and commits a generation only
  after the entire candidate succeeds.
- `BackgroundSnapshotBuilder` evaluates into temporary storage and publishes
  an immutable snapshot only after all points and fields validate.
- `SnapshotBuffer` supports current/next generations and linear interpolation
  only inside a compatible epoch bracket; extrapolation is rejected.
- Batch evaluation returns one status per point and never overwrites output for
  a failed point.
- The AMPS boundary stores every complete background field at its frozen
  offset, then publishes the same immutable object through either the
  standalone or SWMF `Runtime` adapter.

See [BACKGROUND_FIELD.md](BACKGROUND_FIELD.md) for units, field completeness,
coupling ownership, and atomicity rules.

### Phase T: turbulence and scattering inputs

- `TurbulenceProvider` is a separate authority from `BackgroundProvider`.
- The prescribed provider produces a normalized finite-band Kolmogorov
  spectrum with explicit amplitude, wave-number bounds, radial scaling, and
  correlation length.
- The AWSoM adapter accepts `w+` propagating along `+B` and `w-` propagating
  against `+B`, both in J/m³. It uses `δB²=μ₀w` and resolves outward/inward
  labels from the sign of `B·r`, so polarity reversals cannot silently swap
  coupling semantics.
- Out-of-band resonances use an explicit reject or power-law-extension policy;
  extension evaluates the physical resonant wave number and never clamps it to
  a band edge.
- Missing waves fail by default. Zero-rate ballistic transport exists only as
  an explicit configuration and is represented by a typed state.
- `CoefficientBridge` calls `sep_common` directly for Dμμ, mean-free-path,
  and parallel-diffusion conversions. srcSEP3D contains no copied coefficient
  formula.
- The provider API in `turbulence_models.h` is intentionally independent of
  `sep_common`. Coefficient users opt into `coefficient_bridge.h`; this keeps
  the AMPS-facing `main_lib.cpp` compilable under historic `Makefile.conf`
  rules that do not propagate application-local include flags.
- Self-consistent 3-D turbulence remains a hard reserved feature until a mesh
  wave-energy equation and conservation tests exist.

See [TURBULENCE_SCATTERING.md](TURBULENCE_SCATTERING.md) for conventions,
normalization, policies, and the shared-kernel boundary.

### Phase P: transport cores

- The Parker core advances the gyrotropic tensor
  `kappa_perpendicular I + (kappa_parallel-kappa_perpendicular) b b` with the complete Itô drift, including the
  field-aligned coefficient gradient, field-line curvature, and `div(b)`.
- The focused core advances full gyrotropic focusing and flow coefficients
  with a symmetric split, reflecting pitch boundaries, and a declared
  Milstein or Euler–Maruyama stochastic scheme.
- Cell crossing, diffusion, focusing, cooling, background variation, shock
  crossing, and snapshot validity are separate named timestep limits.
- Counter-based random streams are keyed by campaign, particle, step,
  substep, and physical purpose, making histories independent of iteration
  order and worker ownership.
- Controlled constant/constant-ratio perpendicular diffusion and selectable
  gradient-B/curvature drifts are implemented by V01. Current-sheet drift is
  still excluded because the required sheet geometry is unspecified.
- The AMPS local-state resolver now supplies the Parker core's required
  `b·grad(kappa_parallel)` term. It reevaluates the canonical local coefficient
  chain one cell crossing in both field-aligned directions, uses a centered
  difference when both samples exist, and uses an explicit one-sided difference
  at a boundary. If neither neighbor is usable it returns a typed failure; it
  never substitutes the former unconditional zero drift.

See [TRANSPORT_CORES.md](TRANSPORT_CORES.md) for the equations, splitting
algorithm, reproducibility contract, and Phase-P acceptance tests.

### Phase A: AMPS mover and SWCME source adapters

- One AMPS particle-buffer entry point validates and dispatches exactly the
  tensor Parker or split focused core selected by immutable configuration.
- A packed particle extension persists stable ID, stochastic step/substep,
  shock generation, momentum, pitch cosine, and gyrophase through migration
  and AMPS checkpointing.
- The adapter performs deterministic gyrotropic-to-Cartesian velocity
  reconstruction, exact destination-list insertion, and explicit terminal
  deletion/return-code mapping.
- Inner absorption, outer escape, invalid background, and failed transport are
  distinct semantic outcomes. Per-step/species integer ledgers require exact
  closure of active, injected, escaped, absorbed, and failed counts.
- Moving spherical shock crossings use the first analytic segment/surface root
  and are de-duplicated by shock generation.
- The SWCME adapter consumes the canonical common `SEPSourceState`, maps the
  DSA law to the shared `sep_common` injection sampler, and uses independent
  semantic random streams for momentum, pitch, gyrophase, and stable identity.

See [AMPS_ADAPTERS.md](AMPS_ADAPTERS.md) for buffer layout, dispatch, shock
geometry, DSA spectrum mapping, conservation, and host configuration.

### Phase O: sampling, output, and restart

- Read-only particle observations are sorted by stable ID and reduced with a
  specified compensated sum into cell moments, virtual-spacecraft spectra and
  anisotropy, field-line projections, and closed-ledger shock diagnostics.
- Output uses SI unit-bearing CSV schemas, per-artifact hashes, and a manifest
  containing configuration, code, and snapshot identities. A staging-directory
  rename publishes the complete sequence atomically.
- The independent parser verifies manifest keys, exact schemas, and hashes;
  corrupted or partial products never replace caller state.
- The restart codec writes a canonical versioned little-endian image rather
  than C++ object memory. It includes Runtime cadence/checkpoint counters,
  stochastic identity, all active particle state, snapshot/turbulence/source
  generations, sampling state, next stable ID, and closed ledger rows.
- Snapshot mismatch has an explicit reject or bounded-wait policy. Failed
  checkpoint writes roll Runtime back to `SnapshotReady` without incrementing
  the checkpoint sequence.

See [SAMPLING_OUTPUT_RESTART.md](SAMPLING_OUTPUT_RESTART.md) for algorithms,
file schemas, atomicity, restart contents, and lifecycle rules.

### Phase V: integration and scientific validation

- Rank-local observations are merged in stable-ID order, independently of
  rank and receive order. Duplicate physical identities are rejected globally.
- Closed integer particle ledgers are reduced by `(step,species)` with overflow
  checks and exact global conservation. Load imbalance, wall time, and memory
  are evaluated against explicit budgets.
- The positive-series comparator performs log-linear interpolation only within
  model coverage and reports log-space errors, correlation, onset, peak, and
  fluence metrics.
- Absolute, one-global-amplitude, and unit-peak normalization are distinct
  policies; shape-only results cannot be reported as absolute-flux validation.
- Controlled validation calls the production Parker, focused, and SWCME source
  kernels against independent analytical distributions and characteristics.
- Linked `NAT3D`/`MPI3D`, cross-model `XM3D`, and observational `OV3D` cases
  share the public CLI. Missing prerequisites are `SKIP`; malformed or
  checksum-invalid evidence is `ERROR`.
- Scientific cross-model results are immutable exported evidence. The bounded
  V2D01 development gate deliberately compiles distinct cores from both source
  trees; the scientific validation runner itself never searches for srcSEP.

See
[INTEGRATION_SCIENTIFIC_VALIDATION.md](INTEGRATION_SCIENTIFIC_VALIDATION.md)
and [validation/README.md](validation/README.md) for algorithms, metrics,
schemas, case roles, commands, and physical limitations.

## Current limitations

The following are intentionally not enabled:

- current-sheet drift and arbitrary tensor-valued perpendicular closures;
- self-consistent 3-D turbulence evolution;
- external-script background providers;
- unconfigured direct access to mutable SWMF state from mover workers;
- bundled linked/MPI and observational evidence. Those Phase-V gates require
  the configured target executable and independently reviewed evidence bytes.

For coupled operation, the host configures SWMF authority and publishes a
complete initial background, selected turbulence provider, and optional shock
provider before `amps_init()`. Later candidates may be staged only while the
Runtime is joined at `SnapshotReady`; srcSEP3D performs the collective commit,
mover-context update, global observation gather, and checkpoint coordination.
A standalone analytic run constructs the same provider interfaces from the
immutable input configuration.

## Source layout

```text
srcSEP3D/
├── core/                         semantic types and test registry
├── mesh/                         Phase-M resolution, octree, storage, gradients
├── background/                   Phase-B providers and immutable snapshots
├── turbulence/                   Phase-T providers, spectra, coefficient bridge
├── transport/                    Phase-P Parker/focused cores, timestep, RNG
├── adapters/                     Phase-A neutral dispatch, ledger, SWCME source
├── output/                       Phase-O sampling, publication, restart
├── validation/                   Phase-V audits, metrics, registry, runner
├── runtime/                      immutable configuration and lifecycle
├── examples/                     annotated versioned production input
├── amps/                         AMPS-only ABI adapters
├── MESH_STORAGE.md
├── CONFIGURATION.md
├── BACKGROUND_FIELD.md
├── TURBULENCE_SCATTERING.md
├── TRANSPORT_CORES.md
├── AMPS_ADAPTERS.md
├── SAMPLING_OUTPUT_RESTART.md
├── INTEGRATION_SCIENTIFIC_VALIDATION.md
├── RUNTIME_INTEGRATION.md          R01–R07 algorithms and invariants
├── MIGRATION_MANIFEST.md
├── SEP3D.h                       production/coupling interface
├── main_lib.cpp                  AMPS mesh/storage/provider boundary
├── main.cpp                      standalone typed host and driver
├── makefile
└── test/
    ├── run_tests.py              srcSEP-style unified runner
    ├── stage1.cpp                AMPS/MPI-free C++ test executable
    ├── individual-test/          component acceptance callbacks
    └── frozen/                   reviewed byte-exact references
```

Generated objects, archives, binaries, reports, and `test_output/` are not
source deliverables. Overlaying this package on an older checkout cannot
delete stale files; remove an old `AMPS/srcSEP3D/SEP3D.cpp` explicitly if
`BLDL3D02` reports it.

## Layering contract

| Layer | Location | AMPS/MPI allowed? | Responsibility |
|---|---|---:|---|
| L0/L1 | `core/`, `mesh/` | No | types, resolution, standalone octree/storage, gradients |
| L1 | `background/` | No | analytic/imported ambient state and snapshots |
| L1 | `turbulence/` | No | scattering authority, spectra, AWSoM mapping, coefficient bridge |
| L1/L2 | `runtime/`, `transport/`, `adapters/`, `output/`, `validation/` | No | lifecycle, transport, neutral coupling, diagnostics/restart, validation metrics |
| L2 | `amps/` | Yes | model-to-AMPS ABI translation |
| L3 | `SEP3D.h`, `main_lib.cpp`, `main.cpp` | Yes | AMPS allocation, owner-local filling, host entry points |

Every directory except `amps/` and L3 must not include AMPS/MPI headers or
refer to the AMPS namespace. `LAY01`, `LAY02`, and `BLD01` enforce this with a source
scan, negative control, AMPS-free link, and symbol-table inspection.

Because generated `pic.h` exposes `SEP3D.h` to generic AMPS translation units,
the umbrella uses forward declarations for source, restart, and provider
objects. Concrete headers that require `src/models/sep_common` or SWCME are
included only by srcSEP3D implementation files compiled with those paths.
Historic `Makefile.conf` recipes do not consistently expand `CPPFLAGS`,
`CXXFLAGS`, or `INCLUDE`, so the srcSEP3D makefile also exports the canonical
model roots through `CPLUS_INCLUDE_PATH` **only for `MAINLIBOBJ` and
`MAINOBJ`**. This target-scoped environment reaches `mpicxx` even for a fixed
generic recipe without exposing model headers to unrelated AMPS objects.
`BLDL3D06` enforces both the transitive-header boundary and this real compiler
search-path contract.

## Test runner

Run from `AMPS/srcSEP3D`:

```bash
# Discover all IDs, groups, and suites without building.
test/run_tests.py --list

# Fast development gate.
test/run_tests.py --routine --amps-source .. --rebuild

# One implemented phase.
test/run_tests.py --suite phase-m --rebuild
test/run_tests.py --suite phase-b --rebuild
test/run_tests.py --suite phase-t --rebuild
test/run_tests.py --suite phase-p --rebuild
test/run_tests.py --suite phase-a --rebuild
test/run_tests.py --suite phase-o --rebuild
test/run_tests.py --suite phase-v --rebuild
test/run_tests.py --suite improvements-c --rebuild
test/run_tests.py --suite improvements-r --rebuild

# Configured Phase-V executable and independently owned evidence.
test/run_tests.py --suite phase-v --amps ../amps \
  --validation-data /path/to/evidence \
  --validation-launch-prefix "mpiexec -n 8" \
  --output-dir test_output/phase-v

# Complete source and configured-production evidence.
env MAKEFLAGS="-j16" test/run_tests.py --all \
  --amps-source .. --make-config ../Makefile.conf \
  --output-dir test_output/all --rebuild
```

`MAKEFLAGS` controls recursive GNU Make compilation, including the enclosing
AMPS build. Tests themselves remain intentionally sequential. Every run writes
JSON and JUnit summaries. A missing real `Makefile.conf` makes `BLDL3D01`
**SKIP**, never PASS.

For a detached application directory, provide the canonical shared dependency
explicitly:

```bash
test/run_tests.py --routine \
  --sep-common-dir /path/to/AMPS/src/models/sep_common \
  --sep-common-archive /path/to/AMPS/src/models/sep_common/sep_common.a \
  --amps-source /path/to/AMPS
```

See [test/README.md](test/README.md) for the complete evidence catalog,
selection rules, exit codes, and troubleshooting.

## Production build gate

Within a configured AMPS tree:

```bash
make -C srcSEP3D prepare-production
make -C srcSEP3D strict-production
```

`prepare-production` installs the idempotent mover declaration/macro in the
already configured `build/pic/picGlobal.dfn`. The strict target then delegates
to the enclosing `make amps` workflow and
audits `AMPS/build/main/mainlib.a` and `main.a`. This is required because AMPS
copies `srcSEP3D` to `build/main`; a direct source-directory compile lacks the
generated include/definition set and is not production evidence. All paths are
resolved from the active makefile, so source and copied locations find the
same `AMPS/Makefile.conf`, `src/models/sep_common`, and `src/models/swcme`.

Production application objects are intentionally rebuilt on every application
archive invocation. Deterministic release archives use normalized timestamps;
without this guard, overlaying a new package can retain an older
`mesh_model.o` whose unchanged symbols resolve but whose C03/C05 ABI does not.
The archive step also verifies every required member and the current
normalized-domain/preflight definitions before the final Fortran-driver link.
`BLDL3D07` enforces this freshness contract.

## Implemented acceptance groups

| Group | Scope |
|---|---|
| `BLDL3D`, `ARCH3D`, `SWCME3D` | production routing, retired-symbol audit, canonical archives |
| `HARN`, `RUNNER`, `LAY`, `BLD`, `UTIL` | runner, layering, binary boundary, frozen common kernels |
| `LIFE3D01–04` | immutable configuration and complete lifecycle transition matrix |
| `R3D01–07` | mover hook, subcycling, transactional snapshots, clock/events, source, observers, complete restart |
| `CFG3D01–07` | input/CLI, typed contracts, domains, shared Parker geometry, mesh/memory preflight, finite-line schema, proton-only AMPS binding |
| `MSH3D01–10` | resolution bounds/laws, tube geometry, balance, octrees, memory, ownership, presets, gradients, finite-line/origin identities |
| `BGP3D01–06` | analytic Parker identities, component laws, focusing, wind derivatives, polar limits |
| `SNAP3D01–08` | completeness, finite values, units, epochs, atomicity, interpolation, batch status, frame |
| `TUR3D01–04` | spectrum normalization, AWSoM mapping, resonance range, missing-data policy |
| `COEF3D01–02`, `COEF3D06` | six-decade conversions, bitwise shared-kernel identity, and nonzero field-aligned kappa-gradient stencils |
| `COEF3D03–05`, `PRK3D01–08` | tensor assembly/Itô drift and Parker transport behavior |
| `FTE3D01–07`, `RNG3D01–03` | focused transport, pitch boundaries, strong-scattering limit, keyed reproducibility |
| `ADP3D01`, `NAT3D04–05/08`, `SHK3D01–04` | mover dispatch, boundaries, ledger, moving shocks, common SWCME source |
| `NAT3D06–07`, `RST3D01–03` | sampling isolation, transactional output/schema, complete restart |
| `INT3D01–03`, `VFY3D01–05` | deterministic rank audit, scientific metrics, analytical Parker/focused/source validation |
| `NAT3D01–03/09–12`, `MPI3D01–02` | registered configured-host integration and multi-rank gates |
| `XM3D01–06`, `OV3D01–04` | checksum-owned cross-model and observational campaign gates |

## Stage 3 single-species ownership

The current application is explicitly proton-only. `[species] amps_index=0`
and `name=proton` are required configuration, and the index participates in
the physics fingerprint. Every observer species filter must select that same
index. After PIC has initialized its species table, but before any time step,
weight, or injected-particle field is assigned, the application requires one
AMPS species at index zero and compares its SI mass and signed charge with the
immutable configuration to relative tolerance `1e-12`.

The validated index is then reused for PIC time-step/weight assignment,
observer filtering, and source injection; those paths do not maintain separate
species literals. `CFG3D07` exercises the successful binding and independent
negative controls for count, index, name, mass, charge, observer selection,
and fingerprint identity. Multi-species operation remains unsupported because
it requires species-specific transport, source, storage, ledger, output, and
restart semantics rather than removal of a single guard.

## Remaining release evidence

Phase-V algorithms, case registration, evidence validation, and controlled
physics gates are implemented. Closing the production release still requires
running the registered linked cases on the configured AMPS host and supplying
the reviewed cross-model/observational bundles. The R01 hook, pinned resolver,
integer schedule, global observation gather, and restart coordinator are now
implemented; a configured `BLDL3D01` run on the target checkout remains
required after every production-boundary change.

## Schema-version-2 Parker initialization

`examples/sep3d_analytic_parker.in` now uses schema version 2 and adds a
required `[parker_spiral]` section.  It supplies the origin, initial point,
physical arc length, and total number of points in SI units.  Version-1 files
remain accepted and normalize to a deterministic line derived from their
existing domain/tube settings; this preserves archived campaigns and typed
SWMF-host construction.

During standalone initialization the input is parsed before the runtime
lifecycle enters mesh setup.  The finite line is materialized with a
second-order midpoint tangent integration and appears in the dry-run summary.
The production AMPS `localResolution()` callback continues to use the same
analytic Parker geometry, so line sampling density cannot imprint artificial
facets on the refined tube.  The origin is carried through the mesh law and all
radial/tube distances are origin-relative.  The current analytic/SWMF physics
contract still requires a heliocentric zero origin; a nonzero production
origin fails validation rather than being only partially honored.

`CFG3D06` enforces the complete version-2 input and source consistency;
`MSH3D10` enforces point count, arc length, and origin-relative AMR invariance.
`CFG3D07` enforces the AMPS/configuration proton binding. All three are part of
the normal `test/run_tests.py --all` manifest.
