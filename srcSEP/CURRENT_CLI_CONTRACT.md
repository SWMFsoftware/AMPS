# Current srcSEP CLI contract through WP41

The authoritative parser remains `SEP::Util::CLI::ParseCommandLine` in
`util/sep_cli.cpp`.  Step 1 is additive.

## Preserved behavior

- No arguments retain production defaults: coupling/cascade/reflection on,
  integrated turbulence, the canonical `fte-dmumu` mover (the historical FTE
  implementation), 300 injected particles per
  iteration, TestManager off, and spectrum output interval 100.
- Option values are case-insensitive.  Existing Boolean values and aliases,
  including `--cascase` and `--no-cascase`, remain accepted.
- Value-taking options accept both `--option value` and `--option=value`.
- Repeated legacy options retain last-occurrence-wins behavior.
- `-h`/`--help` exits before model initialization; malformed/unknown options
  fail and cannot fall through to production.
- Parsed production options are applied after post-compile input.
- Legacy TestManager controls retain their prior continue-to-production
  behavior and do not imply component-test mode.
- SWMF library entry points do not parse arguments or execute the registry.

## Additive Step 1 behavior

- `--list-tests` lists sorted registry metadata and exits before initialization.
- `--test ID`/`--test=ID` are repeatable.
- `--test-group GROUP`/`--test-group=GROUP` are repeatable.
- `--all-tests` selects the bounded routine set only.
- IDs/groups are case-insensitive; overlaps are de-duplicated and run once in
  stable-ID order.
- Listing cannot be combined with execution.  `--all-tests` cannot be combined
  with explicit ID/group selectors.  Missing/unknown selectors fail early.
- Any execution selector is test-only and returns before TestManager or the
  production timestep loop.  FAIL/ERROR return nonzero; SKIP remains explicit.

## Step 4 mover behavior

- `--particle-mover`, `--mover`, and `--sep-mover` select only `parker`,
  `fte-dmumu`, or `fte-mfp`.
- `--list-movers` lists those three descriptors and exits before initialization.
- Unambiguous legacy names map for one transition release with a warning.
- Ambiguous, direct-wave, and 3-D historical mover names fail before
  initialization.
- Startup metadata prints the canonical name, capability contract, and active
  coefficient provider.

## Step 5 scope behavior

- The CLI continues to expose only `parker`, `fte-dmumu`, and `fte-mfp`; all
  transferred Cartesian, Boris, drift, and Parker3D names remain hard errors.
- The linked application must be compiled with field-line mode and
  field-line-segment particle attachment. There is no runtime option that can
  switch srcSEP back to Cartesian particle transport.
- `--run-test-manager` remains functional, but it no longer needs a runtime
  fallback for builds without field lines because such builds are rejected at
  compile time.

## Step 10 coefficient behavior

- `--coefficient-source` accepts `prescribed`, `self-consistent`, or `swmf`.
- `--spatial-diffusion-provider` accepts `from-dmumu` or `from-mfp`.
- `--pitch-angle-diffusion-provider` accepts `configured`, `constant`,
  `jokipii-1966`, or `florinskiy`. Coupled turbulence sources reject
  `configured` because that compatibility callback is not source-bound.
- `--mean-free-path-provider` accepts `qlt`, `qlt1`, `tenishev-2005`,
  `chen-2024`, or `from-spatial`.
- `--invalid-coefficient-policy` accepts `fail` or `ballistic`; ballistic is
  the exact `lambda=+infinity` zero-event-rate state, not a finite clamp.
- Defaults preserve prior behavior: prescribed source, spatial-from-Dmumu,
  configured Dmumu, Tenishev-2005 MFP, and fail policy.
- Invalid values, conversion cycles, and incompatible source/model pairs fail
  before initialization. Startup metadata prints every canonical selection.

## WP11--WP20 coefficient and error-control behavior

- `--constant-dmumu <1/s>` explicitly overrides the constant provider with a
  finite non-negative SI rate. If absent, a post-compile input value is retained
  and validated.
- `--resonance-gap-policy` accepts `reject` or `ballistic`.
- `--turbulence-amplitude-policy` accepts `reject` or
  `limit-to-mean-field`; limiting is reported by a diagnostic counter.
- `--prescribed-delta-b-over-b`, `--coefficient-correlation-length`,
  `--coefficient-reference-radius`, `--coefficient-k-min`,
  `--coefficient-k-max`, and both radial-exponent options define the named
  prescribed spectrum and correlation scales.
- `--coefficient-quadrature-absolute` and
  `--coefficient-quadrature-relative` control adaptive Dmumu-to-kappa
  integration.
- `--transport-geometry-fraction`, `--transport-deterministic-tolerance`,
  `--transport-stochastic-mu-rms`, `--transport-cooling-log-change`,
  `--transport-focusing-mu-change`, `--transport-shock-fraction`, and
  `--transport-min-step` define one validated mover error budget.
- Parker preflight rejects any selected policy that can deliver an infinite
  kappa to its finite diffusion operator. `fte-mfp` accepts typed ballistic
  lambda and advances with zero scattering rate.
- Per-species source abundance, signed charge, mass, nucleon count, injection
  efficiency, spectral index, and energy convention are currently configured
  through the C++ `SEP::Transport::SpeciesSource` initialization API. Installing
  an explicit table requires an entry for every injected species; no CLI syntax
  for that table is introduced in this work package.

## Step 11 turbulence behavior

- `--turbulence-source` accepts `prescribed`,
  `self-consistent-integrated`, `self-consistent-spectral`, `swmf-read-only`,
  or `swmf-initial-then-local`.
- `--turbulence-model` selects integrated or wave-number-resolved authority;
  self-consistent source and representation remain synchronized.
- `--turbulence-coupling-policy` accepts `disabled` or
  `streaming-energy-exchange`; existing `--coupling` forms remain aliases.
- Inner/outer boundaries independently select `specified-incoming-energy`,
  `specified-incoming-flux`, `transparent-outflow`, or `fixed-reservoir`, with
  SI values supplied by the matching `--turbulence-*-value` option.
- `--turbulence-advection`, `--shock-injection`, `--reflection`, and
  `--cascade` control the documented operator phases.
- Reflection/cascade coefficients, CFL safety, perpendicular correlation
  length [m], spectral k range [1/m], bin count, output cadence, and
  conservation tolerance are validated before initialization.
- Startup metadata prints the complete turbulence configuration.

Step 12 adds no decomposition selector: reproducibility keys intentionally
exclude thread count and MPI rank.

## WP30 authoritative run controls

- `--total-iterations <N>` sets the positive standalone iteration count.
- `--shock-model <analytical|swcme1d>` and `--cme-scenario <fast|slow>` select
  shock/background behavior before initialization.
- `--field-line-seed-area <m2>` supplies the positive area used once to form
  each line's conserved magnetic-flux record.
- `--shock-turbulence-efficiency <0..1>` and
  `--shock-turbulence-plus-fraction <0..1>` configure the WP25 source formula.
- `--merge-minimum <N>` and `--merge-maximum <N>` configure the coupled
  merge/split population range; maximum must not be below minimum.

These values are copied into a validated immutable `RunConfiguration` after
input-file and CLI resolution. Startup prints its fingerprint and source layer.
A restart whose fingerprint differs is rejected rather than silently adopting
the checkpoint or current defaults.

## D01-D02 SWCME background and source controls

- `--swcme-failure-policy` accepts `strict`, `clamp-radius`, or
  `diagnostic-fallback`. `strict` is the production default. Recovery choices
  and their use counters are never hidden.
- `--swcme-fallback-density <m^-3>`, `--swcme-fallback-speed <m/s>`, and
  `--swcme-fallback-divergence <s^-1>` define the validated diagnostic sample.
  They have no physical effect unless `diagnostic-fallback` is selected.
- `--swcme-override key=value` is repeatable and accepts both separated and
  equals-sign CLI forms. Shell quoting is required when the value contains a
  space before its unit, for example
  `--swcme-override 'cme.launch_speed=1400 km/s'`.
- The same canonical keys can be placed in a `SWCME1D on` block in the
  post-compile input. Resolution is preset, then input block, then command-line
  overrides. Duplicate keys within one layer, unknown keys, missing/wrong
  units, and incompatible physics choices fail before mesh initialization.
- `--cme-scenario fast|slow` remains the convenient preset selector and maps to
  the same canonical `preset` key. Existing option spelling and precedence are
  preserved.

The authoritative key/unit schema is `src/models/swcme/swcme1d_input.hpp`, not
the application parser. It covers ambient plasma, Parker normalization, CME
kinematics (including data-driven tables), shock representation, sheath/ejecta
geometry and smoothing, launch/validity epochs, and the one-species source
record. `sheath.compression_floor` is explicitly rejected because it has no
effect in the canonical RH model. Calibrated reference-intensity normalization
is canonical SWCME functionality but is explicitly rejected by the current
srcSEP swept-volume adapter; `relative_only` is the supported 1-D source mode.

After unit conversion and validation, startup prints the normalized manifest
and its 16-hex fingerprint. That canonical fingerprint is embedded in the
frozen run configuration, background snapshot metadata, serialized restart
record, and compatibility check. SWMF/coupled hosts call
`SW1DAdapter::Configure(ConfigurationRequest)` directly; no process arguments,
environment variables, or PARAM parsing occur in the coupled API.

The linked native registry contains extended IDs `D01`, `D02`, and `D03PRE`.
They are discoverable through `--list-tests` and individually selectable with
`--test`. Python `test/run_tests.py --all` includes them because it runs every
discovered ID in a separate process. Native `--all-tests` and Python
`--routine` select only routine-class tests and therefore intentionally exclude
these process-global SWCME reconfiguration checks.

## D03 test-runner controls

D03 does not add production-model options. Its build and campaign controls
belong only to `test/run_tests.py`. `--native-manifest`, `--no-build`, and
`--release` require `--suite d03-native-integration`. `--rebuild` is also
accepted with native `--all`: it runs the enclosing AMPS clean and relocation-
safe `strict-production` build before asking the newly linked executable for
`--list-tests`. `--amps-source` and `--make-config` identify that enclosing
build in either mode. This separation prevents a test orchestrator path from
becoming a second production input parser.

The gate appends the selected canonical `--particle-mover` value to each site
case and rejects a duplicate mover option in manifest arguments. All other
physical configuration continues through the production CLI and the canonical
D02 SWCME schema. `--release` requires `--rebuild`, rejects `--no-build`, and
makes missing configured-build/campaign prerequisites nonzero `INCOMPLETE`
results instead of development `SKIP` results.

## WP31 turbulence operator controls

- `--turbulence-operator-safety <value>` sets the dimensionless local-rate
  accuracy safety in `(0,1]`.
- `--turbulence-max-source-fraction <value>` limits the diagnosed source change
  relative to current wave energy.
- `--turbulence-max-cascade-fraction <value>` limits cascade transfer per stage.
- `--turbulence-min-substep <s>` sets the positive minimum accepted stage.
- `--turbulence-max-substeps <N>` sets the positive stage-count guard.

All options support both `--name value` and `--name=value`. The scheme is
currently the fingerprinted `LieFirstOrder` implementation; no CLI option
silently selects a different temporal order.
