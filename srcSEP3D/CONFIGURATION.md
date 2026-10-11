# Configuration, Geometry, Preflight, and Compiled Species Binding

## Compile-time species selection

The AMPS application deck and the srcSEP3D runtime deck have different roles.
`input/sep3d.input` is processed while AMPS is configured and compiled. Its
`SpeciesList` is the only authority for species count, order, chemical symbol,
mass, and signed charge. The runtime file passed to `./amps -input FILE`
(synonym `--input FILE`; default `./amps.in`) is parsed only after compilation
and cannot select, relabel, add, remove, or mutate an AMPS species. Its
`species` block therefore contains only the common numerical
`macroparticle_weight` (`species.macroparticle_weight`), and that block is
optional and the key forbidden in the derived reduced-front injection mode
described below.

An executable generated with `SpeciesList=ELECTRON` is valid and injects an
electron. A mixed `SpeciesList=H_PLUS ELECTRON` build initializes and injects
both entries in generated index order. Runtime never searches for a proton
macro and never overwrites either molecular-table entry.

For an AMPS checkout already configured as `sep3d`, changing the source deck
does not necessarily refresh the root-level working copy. Use
`cp input/sep3d.input sep3d.input`, run
`./ampsConfig.pl -input sep3d.input -no-compile`, and then perform the normal
clean site build. Before execution, inspect the generated `nTotalSpecies` and
`ChemTable` declarations in `build/pic/pic.h`; startup prints the same table
with mass and signed charge after binding.

## V01 perpendicular diffusion and guiding-centre drift

In the `application: sep` / `transport` block
(`application.sep.transport.*`), `perpendicular_diffusion` accepts `none`, `constant`, or
`constant-ratio`. The selected closure requires its positive SI coefficient:
`constant_kappa_perpendicular_m2_per_s` or
`kappa_perpendicular_to_parallel_ratio`. `drifts` accepts `none`,
`gradient-b`, `curvature`, or `gradient-curvature`. These choices enter the
physics fingerprint. Drift uses the signed charge from the compiled AMPS table and forces
magnetic-gradient storage before the mesh layout freezes. Current-sheet drift
is unsupported because its geometry is unspecified.

This document describes the production input boundary introduced by the
configuration/preflight improvements and extended by Stage 3 species
ownership. The implementation is split between
`runtime/configuration_io.{h,cpp}`, `runtime/run_configuration.{h,cpp}`,
`core/parker_geometry.{h,cpp}`, and `mesh/mesh_model.{h,cpp}`. None of these
files depends on AMPS or MPI.

## C01: input file and standalone CLI

### One generic runtime input format

srcSEP3D reads exactly one runtime file format: the generic AMPS block input
specified normatively in [`../src/INPUT_FILE.md`](../src/INPUT_FILE.md)
(sections 3-5), with the **shared SEP content schema 6** that srcSEP3D and
srcSEP both read (section 5; the block/key table with per-application support
is [`../src/models/sep_common/sep_input_layout.h`](../src/models/sep_common/sep_input_layout.h)).
The retired srcSEP3D INI decks (`[section]` headers, `#` comments), the
retired one-dash shared container (`#section begin: sep3d`,
`#subsection begin: ...`) and the generic-format `#begin application: sep3d`
files with content schemas 1-5 no longer parse: `#section`, `#subsection` and
`[section]` lines are fatal errors, and `application.sep.schema_version`
other than 6 is a typed "retired" error that names the converter
(`tools/convert_amps_input.py sep6-from-srcsep3d`).  The typed
`RunConfiguration3DOptions::inputSchemaVersion` keeps its contract levels
1-5 for parser-free typed hosts and fixtures; every text file installs level
6, which satisfies every "level >= N" completeness rule.

Grammar summary (the specification is authoritative):

- `#begin KIND` or `#begin KIND: NAME` ... `#end KIND` blocks, nestable; every
  `KEY = VALUE` assignment is inside a block and a key may appear once per
  block;
- `!` starts a comment (there is no `#` comment);
- a trailing `\` continues a logical line;
- `#include "PATH"` inserts a block-balanced file; a relative path resolves
  from the directory of the including file, not from the working directory.

A key's **qualified name** is the dot-joined block path, where
`#begin kind: name` contributes `kind.name`: for example
`mesh.domain.inner_radius_m`, `application.sep.schema_version`,
`physics.background.parker.temperature_k`, `source.shock.samples_per_step`,
`observer.earth.cadence_s`.  Diagnostics, tests and this document use these
names.  Fingerprints are computed from normalized values, so comments,
layout and include structure do not change a run's identity.

Shared-schema decisions (normative text: `src/INPUT_FILE.md` 5.1):

- **D1** one application block `#begin application: sep` with
  `schema_version = 6` in both applications;
- **D2** the turbulence description is `application: sep` / `turbulence`;
  it is the only block whose content differs (srcSEP3D keeps its prescribed
  spectrum keys, srcSEP its self-consistent Alfven-wave keys);
- **D3** a setting only srcSEP implements (for example
  `observer.ID.heliocentric_radius_m`, observer `kind = field-line-radius`,
  a srcSEP turbulence key) stops with the typed error `... is a srcSEP setting
  that is not supported by srcSEP3D (shared SEP input schema 6,
  src/INPUT_FILE.md section 5)`; srcSEP rejects the srcSEP3D-only settings
  with the reverse message;
- **D4** every physical quantity is written once: with `physics: swcme` the
  SWCME ambient owns the solar wind and `source: shock` owns the particle
  source (see "SWCME shock-injection contract" below);
- **D5** observers share `kind`, `minimum_energy_j`, `maximum_energy_j`,
  `energy_bins`, `energy_spacing`, `pitch_angle_bins`; srcSEP3D implements the
  spatial kinds.

File selection and CLI:

```bash
./amps -input srcSEP3D/examples/sep3d_analytic_parker.in --dry-run
./amps --input srcSEP3D/examples/sep3d_analytic_parker.in --dry-run  # exact synonym
./amps                       # reads ./amps.in in the working directory
mpiexec -n 4 ./amps -input srcSEP3D/examples/sep3d_analytic_parker.in \
  --initialization-only --initialization-output-dir mesh-preview
./amps -input run.in --output-dir products
./amps -input run.in --restart restart/sep3d.chk --log-level verbose
```

`-input FILE` and `--input FILE` are exact synonyms.  Without either, the
standalone executable reads `./amps.in`; a coupled SWMF run has no command
line and reads `./amps.in` in `amps_init_mesh` when no host installed a typed
configuration.  Native-test selection modes do not fall back to `./amps.in`;
they require their explicit reviewed deck (`--test-input`, unchanged).
`--dry-run` and `--restart` apply to every input file (previously they were
unavailable in the one-dash mode).

`--output-dir` and `--restart` override only their corresponding file fields.
`--dry-run` parses, normalizes, validates, freezes, fingerprints, samples the
resolution law, and estimates memory; it does not allocate the AMPS mesh.
`--initialization-only` instead follows the production path through AMPS mesh
construction and complete application initialization, writes all three declared
Tecplot products, synchronizes all MPI ranks, and exits before the first
particle step. `--initialization-output-dir` requires that mode and changes
only the parent directory of the three initialization filenames.
Unknown options and conflicting test selectors are usage errors.

The public test interface remains `test/run_tests.py`; it provides `--list`,
`--test`, `--group`, `--routine`, `--all`, and `--suite`. The production
executable recognizes test-selector spellings only to reject an accidental
simulation launch with a diagnostic directing the user to the test runner.

### Block layout (shared SEP schema 6)

Owner: `runtime/configuration_io.cpp`, which walks the document against the
shared table `SEP::InputLayout::kBlocks`/`kKeys`
(`src/models/sep_common/sep_input_layout.h`).  A block or key that is not in
the table, at that nesting position, is a fatal error; one that the table
marks srcSEP-only is a typed "not supported by srcSEP3D" error (D3).

| Block path | Keys (qualified prefix) | Former spelling (retired) |
|---|---|---|
| `run` | `run.time_step_s`, `run.time_step_mode`, `run.maximum_particle_speed_m_s`, `run.time_step_margin_factor`, `run.maximum_time_steps`, `run.campaign_seed`, `run.stop_shock_radius_m` | INI `[run]` |
| `run` / `restart` | `run.restart.input_path`, `run.restart.output_path` | INI `[restart]` |
| `application: sep` | `application.sep.schema_version` (= 6), `.transport`, `.intent`, `.background_cadence_steps`, `.injection_cadence_steps` | `application: sep3d` (schemas 1-5); `.particle_mover` retired (the mover family follows from `transport`) |
| `application: sep` / `transport`, `population_control`, `storage` | `application.sep.transport.*`, `application.sep.population_control.*`, `application.sep.storage.*` | `application.sep3d.*` |
| `application: sep` / `turbulence` | `application.sep.turbulence.*` (prescribed spectrum, listed below) | `physics: turbulence` |
| `mesh` | `mesh.*` (global mesh keys) | INI `[mesh]` |
| `mesh` / `domain` | `mesh.domain.*` | INI `[domain]` |
| `mesh` / `field_line` | `mesh.field_line.*` | INI `[parker_spiral]` |
| `mesh` / `refinement: solar`, `refinement: tube` | `mesh.refinement.solar.*`, `mesh.refinement.tube.*` | INI `[mesh.solar]`, `[mesh.tube]` |
| `mesh` / `active_region`, `memory` | `mesh.active_region.*`, `mesh.memory.*` | INI `[mesh.active_region]`, `[memory]` |
| `physics: background` (+ nested `parker`) | `physics.background.provider`, `.model_id`; `physics.background.parker.validity_cadence_s` and, only without `physics: swcme`, the seven ambient keys (D4) | `.external_script` and `parker.number_density_at_reference_m3` retired |
| `physics: shock` | `physics.shock.authority` (`none` or `swcme`) | the constant-speed surrogate keys `active_from_s`, `active_until_s`, `initial_radius_m`, `maximum_radius_m`, `speed_m_per_s`, `compression_ratio` retired |
| `physics: swcme` | canonical SWCME3D vocabulary (units in the value); no `source.*` keys (D4) | SWCME `source.*` moved to `source: shock` |
| `physics: parallel_diffusion` | shared parallel-diffusion library vocabulary | -- |
| `physics: shock_front` | reduced shock-front event (dotted keys), normally `#include`d from its `.event` file | INI `[background] model_asset = FILE` |
| `source: shock` | `source.shock.*`, including the SWCME reference spectrum (D4) | `source.shock.spectral_index` retired |
| `species` | `species.macroparticle_weight` | -- |
| `observer: ID` | `observer.ID.*`; `energy_bins` | `energy_channels` |
| `output` | `output.*`, `output.initialization_field_line_tecplot_file` | `output.initialization_parker_line_tecplot_file` |

The srcSEP3D turbulence block (`application.sep.turbulence.*`) keeps the
former `physics: turbulence` keys unchanged: `authority`, `model`,
`amplitude_model`, `delta_b_over_b`,
`wave_energy_density_at_reference_j_per_m3`,
`wave_energy_density_radial_exponent`, `normalized_cross_helicity`,
`reference_radius_m`, `k_min_per_m`, `k_max_per_m`, `k_min_radial_exponent`,
`k_max_radial_exponent`, `spectral_index`, `correlation_length_m`,
`correlation_length_radial_exponent`, `validity_cadence_s`, `missing_data`,
`resonance_range`, `self_consistent_3d`.  srcSEP's 28 self-consistent
turbulence keys (`source`, `representation`, `coupling`, ...) are typed
"not supported by srcSEP3D" errors.

`application.sep.transport` is required and accepts `parker`,
`focused-diffusion` or `focused-scattering`; the aliases `parker3d` and
`focused3d` and the legacy `false` spellings of
`application.sep.transport.perpendicular_diffusion`/`drifts` are retired
(write `none`).

**Completeness.**  Every file contains the blocks `run`, `run` / `restart`,
`application: sep` with its `transport`, `population_control`, `storage` and
`turbulence` children, `mesh`, `mesh` / `domain`, `mesh` / `field_line`,
both `refinement` blocks, `mesh` / `active_region`, `mesh` / `memory`,
`physics: background` with `parker`, `physics: shock`, `source: shock`,
`output`, and `species` (except in the derived reduced-front injection mode).
`physics: swcme` is required unless `physics.background` selects the reduced
shock front (`provider = runtime-model`, `model_id =
sep-corona-swcme-shock-front-v1`, plus a `physics: shock_front` block).  At
least one `observer: ID` block is required except for a `shock-propagation`
run, and every observer block contains all 21 srcSEP3D observer keys.  Every
srcSEP3D key of a required block is required, even when the selected mode
makes it inactive; the only exceptions are the mode-dependent keys checked by
their own rules: the derived-numerics keys (below), the optional
`run.stop_shock_radius_m`, the corner-cube group (`mesh.domain.box_geometry`,
`corner_direction_*`, `corner_margin_m`, `mesh.field_line.end_radius_m`,
`mesh.refinement.solar.anchor`, `mesh.active_region.solar_sphere_radius_m`,
required together for a corner box), `physics.background.model_id`
(runtime-model only), and the D4 keys (Parker ambient keys without SWCME
only; reference-spectrum keys with SWCME only).  A disabled feature still has
its block, making the decision visible in code review.  Duplicate keys,
unknown keys, invalid enumerations, missing blocks, non-finite values, and
incompatible choices fail in the parser or immutable factory, with a
`file:line:` diagnostic naming the actual (possibly included) file.

Old files are converted mechanically, keeping comments (except those attached
to a removed key) and values (unit conversions are exact):

```bash
python3 tools/convert_amps_input.py sep6-from-srcsep3d old.in  -o new.in   # generic-format srcSEP3D schemas 1-5
python3 tools/convert_amps_input.py srcsep3d           old.ini -o new.in   # INI srcSEP3D (chains into schema 6)
python3 tools/convert_amps_input.py event              old.event -o new.event
```

The reduced shock-front event file now contains a
`#begin physics: shock_front` ... `#end physics` block (its own dotted keys,
such as `ambient.reference_radius_m`, are unchanged).  It replaces the retired
`[background] model_asset = X.event` reference and is included at the top
level of the deck:

```text
#begin physics: background
  provider = runtime-model
  model_id = sep-corona-swcme-shock-front-v1
  ...
#end physics
#include "positive_1au.event"     ! supplies the physics: shock_front block
```

A `physics: shock_front` block is accepted only together with
`physics.background.provider = runtime-model` and
`physics.background.model_id = sep-corona-swcme-shock-front-v1`.  Relative
asset paths inside the block (for example `assets.harmonics_file`) resolve from
the directory of the file that holds the block.

Every dimensional application key includes its SI unit, for example
`inner_radius_m`, `time_step_s`, `radial_field_at_reference_t`, and
`minimum_energy_j`. A key without the declared suffix is unknown. The parser
does not guess whether a number is AU, solar radii, km/s, nT, or MeV.

### Reduced-front injection: derived particle numerics

The keys below, formerly reachable only through the retired one-dash
container, are now ordinary `run` and `source: shock` keys.  The governing
equations are in [`../src/INPUT_FILE.md`](../src/INPUT_FILE.md), section 5.

| Qualified key | Values | Meaning |
|---|---|---|
| `run.time_step_mode` | `fixed` (default) or `derived-from-mesh` | `derived-from-mesh`: after mesh allocation, `dt = time_step_margin_factor * h_min / maximum_particle_speed_m_s`, with `h_min` the global smallest allocated characteristic cell length |
| `run.maximum_particle_speed_m_s` | `(0, c]` | `v_max` of the derived step |
| `run.time_step_margin_factor` | `(0, 1]` | dimensionless margin `f` |
| `source.shock.normalization` | `configured-rate` (default) or `accepted-shock-incident-flux` | incident flux: `Ndot_s = sum over accepted curved faces of n_s [V_n - U_1 . n]_+ A_face`, `W_s = Ndot_s dt / samples_per_step` for every compiled species |
| `source.shock.normalization_radius_m` | positive, no default | heliocentric apex radius at which `Ndot_s` (hence `W_s`) is evaluated |
| `source.shock.weight_model` | `constant-statistical-weight` or `log-uniform-momentum-importance` | the importance mode is a typed **reserved** mode that stops startup until its individual-weight correction is implemented and qualified |
| `source.shock.maximum_events_per_species_per_step` | positive integer | fatal Poisson runaway guard; never a cap or renormalization |

Constraints enforced by the parser:

- `run.time_step_mode = derived-from-mesh` and
  `source.shock.normalization = accepted-shock-incident-flux` must be selected
  together;
- in that mode `run.time_step_s`, `source.shock.physical_particle_rate_per_s`,
  `source.shock.injection_efficiency` and `species.macroparticle_weight` must
  be absent (the `species` block is optional), and all five derivation inputs
  above must be present; in the configured modes those five inputs are
  rejected;
- the mode requires `application.sep.intent = shock-injection`,
  `source.shock.enabled = true`, `physics.shock.authority = none`,
  `physics.background.provider = runtime-model` with
  `model_id = sep-corona-swcme-shock-front-v1`, and a `physics: shock_front`
  block; consequently there is no `physics: swcme`, the SWCME
  reference-spectrum keys of `source: shock` must be absent, and the seven
  Parker ambient keys must be present (D4).  The mover family follows from the
  required `application.sep.transport`; there is no separate mover key.

Spellings of the retired one-dash container map as follows:

| Retired one-dash key | Current spelling |
|---|---|
| `shock_model`, `background_plasma_model`, `source_model` | `physics.background.provider`/`model_id`, a `physics: shock_front` block, and `source.shock.normalization` |
| `particles_per_iteration` | `source.shock.samples_per_step` |
| `statistical_weight_model` | `source.shock.weight_model` |
| `phase_space_power_model = compression-ratio` | `source.shock.spectrum_model = local-compression-dsa` |
| `phase_space_power_model = constant` | `source.shock.spectrum_model = fixed-phase-space-power-law` with `source.shock.phase_space_power_index` |
| `source_normalization_radius_m` | `source.shock.normalization_radius_m` |
| `maximum_events_per_species_per_step` | `source.shock.maximum_events_per_species_per_step` |
| `maximum_particle_speed_m_s`, `time_step_margin_factor` | `run.maximum_particle_speed_m_s`, `run.time_step_margin_factor` |

The one-dash mode took every other value from compiled
`RunConfiguration3DOptions` defaults; the current example
[`examples/application-input/sep3d.in`](examples/application-input/sep3d.in)
writes those values explicitly, marked `[former default]`.  Status: parsing,
the constraints above and the CLI spellings are covered by `CFG3D13`, and the
schema-6 spellings of every example deck by `CFG3D18`.
Native evidence after the generic-format migration (2026-10-10, recorded
with the schema-5 spellings; the schema-6 rerun of 2026-10-11 in
`test_output/run_all_test/srcsep3d-schema6` reproduces the same counts, `sep3d` build whose only
compiled species is `H_PLUS`): the migrated reduced-front decks pass the
reduced shock-front runner (97/97) and the positive 1-AU qualification
(28/28) in `run_all_test.sh`; the reduced event identity is unchanged.  The
derived-numerics deck `examples/application-input/amps-injection-smoke.in`
injects `H_PLUS` with an identical receipt on one and four ranks (weight
4.27644841713363414e32, 993 then 1064 Poisson candidates, all injected, zero
disconnected; configuration fingerprint `b882ccfb611c1b9a`).  The compiled
`H_PLUS` (AMPS mass 1 u) is identified as the ambient proton population by
charge number and mass number; see
[`examples/application-input/README.md`](examples/application-input/README.md).

The complete commented schema-6 example is
[`examples/sep3d_analytic_parker.in`](examples/sep3d_analytic_parker.in),
exercised directly by `CFG3D01`.

`mesh.domain.inner_radius_m` does **not** set the radius of the Sun. It controls the
Parker/CME source surface, injection/transport cutoff, and beginning of the
configured finite field line. srcSEP3D separately registers the physical
photosphere in AMPS at the fixed reviewed constant
`Core::Const::R_sun = 6.957e8 m`, centered on the validated domain origin.
There is deliberately no duplicate photospheric-radius input key. The source
shell may equal `R_sun` but may not lie below it; the shipped standalone
examples retain their physically separate 20-solar-radius source shell.

A flat assignment outside any block (for example `scattering = ...` at the
top level) is a grammar error rather than being guessed into
`application: sep` / `turbulence` or `application: sep` / `transport`.  This is
intentional: those two blocks represent different physical contracts.

### Schema 6 complete contract: SWCME shock injection

Every schema-6 file contains the required `mesh` / `field_line` block (see
"Required finite field-line block" below), explicit initialization Tecplot
paths, complete observer acceptance fields, and an exact particle source.  A
standalone SWCME shock-injection run additionally requires:

- `application.sep.intent = shock-injection`,
  `physics.shock.authority = swcme`, and `source.shock.enabled = true`
  (or `intent = shock-propagation` with the source disabled);
- `application.sep.injection_cadence_steps = 1`, so
  `source.shock.samples_per_step` means exactly that many computational particles for
  every compiled AMPS species on every active simulation step;
- `physics.background.provider = analytic-parker` or `swcme` and prescribed
  turbulence (`application.sep.turbulence.authority = prescribed`) selected as
  `kolmogorov`, `kraichnan`, or explicit `power-law`. The recognized
  `python-interpolator` value is reserved and fails before allocation; SWMF/AWSoM
  coupling remains available through the parser-free typed host interface;
- a canonical shock of shape `sphere` or `sse` with `shock_only` region
  behavior and source-mode acceleration (or the source-free
  `full_icme`/`resolved_compression` mesh background), and
  `source.shock.spectrum_normalization = relative-only`.  An ellipsoid is
  rejected because no exact ellipsoid crossing operator exists;
- a complete `physics: swcme` assignment layer: the common SWCME1D fields plus
  shape, axis ratios, half width, CME direction, solar rotation axis/rate, and
  surface theta/phi resolution, but **no** `source.*` key. Data-driven
  kinematics additionally owns its two knot lists;
- the SWCME Parker source radius equal to `mesh.domain.inner_radius_m`, the
  SWCME `parker.sin_theta` equal to the sine of
  `mesh.refinement.tube.source_colatitude_rad`, and the rotation axis +Z; and
- nonempty `output.initialization_mesh_tecplot_file`,
  `output.initialization_field_line_tecplot_file`, and
  `output.initialization_data_tecplot_file`.

**D4: one owner per physical quantity.**  With `physics: swcme` the SWCME
ambient is the only written description of the solar wind, so the Parker
keys `physics.background.parker.reference_radius_m`,
`radial_field_at_reference_t`, `solar_rotation_rate_rad_per_s`,
`solar_wind_speed_m_per_s`, `magnetic_polarity`,
`number_density_at_one_au_m3` and `temperature_k` are illegal (only
`validity_cadence_s` remains).  The parser copies the resolved SWCME ambient
into the typed Parker provider: wind speed, rotation rate and axis, polarity,
density at 1 AU, proton temperature, adiabatic index, composition, electron
and alpha temperatures, thermodynamic closure and reference sin(colatitude).
The radial field is expressed at the 1-AU reference radius.  SWCME supplies
the *total* `|B|` at one AU and its reference latitude, so the parser first
removes the one-AU Parker winding,

```text
w     = Omega * (1 AU - r_source) / V_sw * sin(theta_ref)
B_r(1 AU) = |B|_1AU / sqrt(1 + w^2) ,
```

and then applies `B_r proportional to r^-2` (identity at the 1-AU
reference).  `CFG3D08` checks this derivation against the former
hand-computed values `3.585667175612118e-9 T` (20 R_sun source) and
`3.6305743892405979e-9 T` (25 R_sun source).  Species identity remains
exclusively in AMPS' compiled table.  Without `physics: swcme` (the reduced
shock front) the seven ambient keys are required and the Parker record is
the ambient description used by the mesh and diagnostics.

`source: shock` owns the particle source, so the SWCME `source.*` keys are
illegal in `physics: swcme`.  The reference spectrum is written once in
`source: shock` and is required exactly when `physics: swcme` is present:
`spectrum_reference_energy_j`, `spectrum_particle_mass_kg`,
`spectrum_charge_number`, `spectrum_normalization` (`relative-only` or
`reference-differential-intensity`; srcSEP3D requires `relative-only`) and
`spectrum_reference_intensity_per_m2_s_sr_j`.  Together with
`minimum_energy_j`, `maximum_energy_j` and `injection_efficiency` the parser
synthesizes the SWCME `source.particle_mass`, `source.charge_number`,
`source.energy_min`, `source.energy_max`, `source.reference_energy`,
`source.injection_efficiency`, `source.normalization` and
`source.reference_intensity_si` assignments in exact round-trip decimal text
with units J, kg and si (`SEP::InputLayout::SynthesizeSwcmeSource`).  Each
synthesized assignment carries the line of the `source: shock` key that
supplied it, so a canonical-resolver diagnostic points at an editable line.

Complete prescribed turbulence input additionally declares `model`,
`normalized_cross_helicity`, `reference_radius_m`, independent radial exponents
for both wave-number bounds, a correlation-length radial exponent, and provider
validity cadence. Kolmogorov requires `q=5/3`, Kraichnan requires `q=3/2`, and
the general power-law choice requires an explicit finite `q>1`; a contradictory
named model/index pair fails rather than being normalized silently.

All application fields and every field of every `observer: ID` block are
explicit, including values inactive under the chosen mode. The
`physics: swcme` resolver is the
only code that interprets its unit-bearing values. Its normalized manifest and
fingerprint are frozen into the application identity; raw spelling is not used
as a substitute for resolved physics.

Before AMPS allocation, the standalone provider evaluates the canonical model
at `event.valid_from`, builds and validates the complete shock surface, solves
the MHD jump/source state on every patch, and verifies that
the per-species `source.shock.samples_per_step` is at least the active patch count.
This distinguishes
a legitimately delayed event (inactive before `valid_from`) from an invalid or
empty source.

The required particle-weight identity is

```text
species.macroparticle_weight =
  source.shock.physical_particle_rate_per_s
  * source.shock.injection_efficiency
  * swcme.shock.relative_source_weight_per_area
  * run.time_step_s
  / source.shock.samples_per_step
```

Each species' count is apportioned with a deterministic largest-remainder
allocation after reserving one representative per active physical patch.
Per-patch individual weight corrections preserve represented physical number
exactly. The global and every local AMPS timestep/weight are then set for every
compiled species from the frozen configuration.

The full field-by-field tables and `physics: swcme` key list are in the top-level
[`README.md`](README.md); the annotated file is executable acceptance input,
not pseudocode.

### Transport-domain and statistical controls

Schema 6 requires three independent decisions that must not be inferred from
one another:

1. `mesh` / `refinement: tube` defines where AMR is refined.
2. `mesh` / `active_region` defines which finalized leaf blocks participate
   in the calculation.
3. `application: sep` / `population_control` defines how many computational particles represent
   the physical distribution in each occupied cell/species.

The active-region syntax is:

```text
#begin mesh
  #begin active_region
    mode = full-domain                    ! or parker-tube
    reference_radius_m = 1.495978707e11
    radius_at_reference_m = 0             ! positive for parker-tube
    radius_mode = constant-angular-width  ! or physical-constant
    buffer_blocks = 0                     ! at least 1 for parker-tube
  #end active_region
#end mesh
```

`parker-tube` uses the same polarity-independent centreline as refinement and
the analytic Parker field. Its radius is an independent, normally wider,
physical boundary. Configuration checks both endpoints and rejects a radius
narrower than the refinement tube anywhere on the finite active line, including
mixed constant-width/angular-width selections. The physical seed is the
finite `mesh.field_line` curve, clipped at the physical outer sphere, and uses
conservative curve-capsule/leaf-box intersection. `buffer_blocks` is an exact
count of coarse/fine-aware touching-neighbour layers; it is not a distance in
units of each leaf's diagonal. Bounded inactive cavities are filled, and the
final active leaves must form one face-connected source-to-endpoint component.
A fixed observer is accepted when its collection sphere intersects the finite
physical active tube or its end cap. Because AMPS freezes block activation before allocation,
a moving or coupled field-connected observer is incompatible with a static
active tube until a reviewed dynamic reactivation/initialization protocol
exists.

Transport selection uses `application.sep.transport` plus explicit
coefficient choices in `application.sep.transport.*`:

| Mover | Required local coefficient | Selectors |
|---|---|---|
| `parker` | spatial `kappa_parallel` and its field-aligned gradient | `spatial_diffusion_model = mean-free-path` uses the selected mean-free-path model and `kappa=v lambda/3`; `pitch-angle-integral` evaluates the selected Dμμ model with the declared adaptive-quadrature tolerances. The old `correlation-mean-free-path` spelling is an input-only alias. |
| `focused-diffusion` | `D_mumu` and `dD_mumu/dmu` | `pitch_angle_diffusion_model = jokipii-1966`, `florinskiy`, or `constant` |
| `focused-scattering` | `lambda_parallel` | `mean_free_path_model = correlation`, `constant`, or `radial-rigidity-power-law`; frame is `plasma-frame-isotropic` or `alfven-wave-frame-isotropic` |

`constant` Dμμ requires positive `constant_dmumu_per_s`; constant mean free
path requires positive `constant_mean_free_path_m`. Each inactive constant is
zero in complete input. Discrete focused scattering currently rejects
perpendicular diffusion because no event-partition-invariant transverse
operator has passed validation. The resolver evaluates only the quantities the
selected mover consumes.

Source spectral shape is selected independently of the mover:

```text
#begin source: shock
  spectrum_model = local-compression-dsa
  phase_space_power_index = 0
  ...
#end source
```

`local-compression-dsa` uses the canonical SWCME value
`q=3*r/(r-1)` independently for every active shock patch and requires the
fixed-index field to be the explicit zero sentinel. For an externally
prescribed isotropic seed law, use:

```text
#begin source: shock
  spectrum_model = fixed-phase-space-power-law
  phase_space_power_index = 5
  ...
#end source
```

Here the input is the positive \(q\) in \(f(p)\propto p^{-q}\), not the signed
textual exponent `-5` and not the sampled `dN/dp` exponent. The adapter applies
the spherical momentum-space Jacobian and samples
`dN/dp proportional to p^(-(q-2))`. Both selectors retain the declared total
kinetic-energy bounds and convert those bounds separately with every compiled
AMPS species mass. The selector does not change source normalization:
`source.shock.physical_particle_rate_per_s`,
`source.shock.injection_efficiency`, patch weight, cadence, and the exact
sample count retain their documented meanings (in the reduced-front derived
mode the normalization is instead the accepted incident flux).

The radial–rigidity closure is
`lambda=lambda_ref*(r/r_ref)^a*(R/R_ref)^b`, with
`R=p*c/abs(q)` in volts. Its five keys are
`mean_free_path_reference_m`, `mean_free_path_reference_radius_m`,
`mean_free_path_reference_rigidity_v`,
`mean_free_path_radial_exponent`, and
`mean_free_path_rigidity_exponent`. The three reference values must be
positive and all values finite when that model is selected. For either other
MFP model all five are mandatory zero sentinels in schema 6. This rule makes
the selected physics reviewable and prevents hidden reference scales from
becoming active after a one-line model change.

Population control syntax is:

```text
#begin application: sep
  #begin population_control
    mode = off                          ! or split-merge
    minimum_particles_per_cell_per_species = 0
    target_particles_per_cell_per_species = 0
    maximum_particles_per_cell_per_species = 0
    cadence_steps = 1
  #end population_control
#end application
```

For `split-merge`, limits satisfy
`2 <= minimum <= target <= maximum`; `off` requires all three limits to be
zero. Counts are per allocated AMR cell and per compiled species. At a control
boundary, a nonempty count below minimum is split to target and a count above
maximum is merged to target; counts inside the band and empty cells are left
unchanged. The controller runs after that boundary's shock injection, so the
upper bound describes the representation entering the next transport step.
The dry-run summary prints all selectors, active-tube values, and population
limits alongside the resolved fingerprint.

## C02: typed run contract and fingerprints

The parser produces `RunConfiguration3DOptions`, the same SI-only record a
coupled SWMF host constructs directly. A host that installs this typed record
does not parse a file; when no host installed one, the coupled run reads
`./amps.in` (see C01). There is no second configuration representation. The successful
`RunConfiguration3D::Create` call returns an immutable normalized object.

The typed contract includes:

- run intent, transport core, time step, maximum steps, random campaign, and
  integer background/injection/sampling/checkpoint cadences;
- explicit domain, boundary, coordinate-frame, refinement, active-block
  corridor, and storage choices;
- complete SWCME-backed analytic Parker/solar-wind parameters and a typed
  reserved Python-interpolator authority;
- turbulence model, amplitude, cross helicity, radial scalings, spectrum,
  cadence, and out-of-range/missing-data policies;
- one of three transport movers, independent spatial/pitch/MFP coefficient
  selectors, scattering frame/event budget, and numerical quadrature values;
- per-cell/per-species split-merge hysteresis and cadence;
- the canonical SWCME3D manifest/fingerprint (the typed record still carries
  the constant-speed surrogate shock interval/radius/speed/compression for
  parser-free hosts at contract levels 1-2; no text input can set them, the
  former `physics.shock.*` surrogate keys being retired);
- source efficiency, physical particle rate, energy interval, spectrum, and
  maximum samples per injection event;
- species AMPS index, name, mass, signed charge, and macroparticle weight;
- named observer geometry/trajectory, collection or shell radius, cadence,
  species/pitch acceptance, energy range/bin counts, normalization, and products;
- output/restart paths and cadence.

The run intent prevents ambiguous source behavior. `transport-only` rejects an
enabled source. `shock-injection` requires both SWCME shock authority and an
enabled source, except in the reduced-front injection mode, which
requires `physics.shock.authority = none` and the reduced shock front (C01). SWMF turbulence requires SWMF background authority, and an
`imported-coverage` outer boundary likewise requires SWMF background.

The physics fingerprint includes all trajectory, source, observer, mesh,
layout, and resource-model inputs. Output directory, output prefix, restart
paths, and output cadence are included in the resolved manifest but excluded
from the physics fingerprint. Consequently a cosmetic output relocation can
resume the same physics, while a species weight, observer, mesh, shock, or
transport change cannot masquerade as the same run. Restart loading compares
the frozen physics/layout identities transactionally.

The photospheric radius is not an input option, but the resolved manifest and
physics fingerprint explicitly record
`solar_boundary=amps-absorbing-sphere-v1` and the fixed SI
`solar_radius_m`. Introducing the AMPS cut-cell boundary changes mesh/transport
semantics, so the fingerprint domain is `sep3d-physics-v8`; checkpoints
produced before this boundary contract cannot be resumed as though their
geometry were identical.

### Complete compiled AMPS binding

Immediately after `PIC::Init_BeforeParser()`, production enumerates
`0..PIC::nTotalSpecies-1` and copies `GetChemSymbol`, `GetMass`, and
`GetElectricCharge` into AMPS-independent `CompiledSpeciesRecord` values.
`ValidateCompiledSpeciesBinding` requires the copied count to equal the
generated count, indices to be contiguous, symbols to be unique and nonempty,
masses to be finite and positive, and charges to be finite and nonzero. A
neutral entry fails because the selected focused SEP scattering operator
requires charge; treating it as charged would not be physically correct.

The adapter does not call the AMPS molecular-data setters. It then installs the
explicit runtime timestep and common base statistical weight in every compiled
global slot and every allocated local block. Observer comma-separated numeric
indices are checked against the generated count during this binding, when that
count is authoritative and available. `species = all` selects the entire
compiled table without embedding a build-specific count in the runtime file;
use numeric indices only for an intentional observer subset.

Source controls have per-compiled-species semantics. Every species receives
`source.shock.samples_per_step` computational particles and the declared per-species
physical rate. The kinetic-energy bounds are converted to momentum separately
with that species' AMPS mass; SWCME shock geometry and compression remain
common. This is complete multi-species initialization/injection without a
second runtime species-definition authority.

## C03: domain and boundary contract

`outer_radius_mode=preset` resolves a named radius during normalization:

| Preset | Resolved outer radius |
|---|---:|
| `solar` | 0.30 AU |
| `one-au` (`earth` input alias) | 1.00 AU |
| `mars` | 1.666 AU |

`outer_radius_mode=explicit` uses `outer_radius_m`; no numeric sentinel means
“automatic.” The current geometry is heliocentric, so the coordinate origin
must be `(0,0,0)` and the named frame must agree with the Parker or supplied
SWMF frame. Observers, shock initial radius, shock maximum radius, Parker
reference radius, and tube reference radius are validated against the
normalized domain before mesh initialization.

The configured Parker/CME transport source shell is absorbing. An outward
segment that crosses the outer sphere returns `DomainExit`; the message
distinguishes ordinary escape from leaving imported SWMF coverage. Direction
matters: a segment entering the outer sphere is not an escape, and a segment
moving outward across the source shell is not inward absorption. This
transport classifier is separate from the fixed AMPS photosphere described
above.

## C04: one Parker geometry

`ParkerSpiralGeometry` is the single geometric definition used by both the
mesh tube and analytic background. For source position \(\mathbf{x}_0\),
rotation axis \(\hat{a}\), source radius \(r_0\), rotation rate \(\Omega\), and
wind speed \(V\), the curve rotates the source direction by

\[
\Delta\phi(r)=-\frac{\Omega}{V}
\left[(r-r_0)-r_0\ln\left(\frac r{r_0}\right)\right],\qquad r\ge r_0.
\]

This is the exact integral curve of the initialized SWCME field ratio
\(B_\phi/B_r=-\Omega(r-r_0)\sin\theta/V\). Its arc-length derivative is

\[
\frac{ds}{dr}=\sqrt{1+
\left[\frac{\Omega\sin\theta}{V}(r-r_0)\right]^2}.
\]

Initialization inverts the corresponding analytic arc-length expression, so
`length_m` and the finite active mask share one geometry and do not depend on
the visualization `point_count`.

Its local outward tangent is

\[
\hat{t}=\operatorname{normalize}
\left[\hat{r}-\frac{\Omega\max(0,r-r_0)}{V}
(\hat{a}\times\hat{r})\right].
\]

Magnetic polarity is deliberately absent from this geometry. The analytic
provider applies polarity only after computing the unsigned tangent, reversing
`B` and pitch orientation without moving the mesh tube. Tube distance is the
spherical transverse arc
`r*atan2(|r_hat cross t_hat|, r_hat dot t_hat)` at the point's radius.

## C05: composite refinement and whole-run memory

The global target is `global_cell_size_m`. Near the Sun, normalized coordinate

\[
s=(r-r_{in})/(r_{transition}-r_{in})
\]

is clamped to `[0,1]` and mapped by `linear`, `power-law`, or `smoothstep`.
The result interpolates from `surface_cell_size_m` to the global target.

The tube has a declared physical radius at a declared heliocentric reference
distance. `physical-constant` keeps that radius fixed;
`constant-angular-width` scales it as `radius_at_reference_m*r/reference_r`.
The transverse profile interpolates from `center_cell_size_m` on the Parker
centerline to the global target at the tube boundary. Where near-Sun and tube
regions overlap, the finer request wins. The final value is clamped between
the declared minimum and global sizes.

The dry-run preflight reports sampled minimum/maximum requests and locations,
estimated blocks by AMR level, tube radius at the reference distance, and a
conservative whole-run memory model. Its categories are:

- base cell objects;
- frozen associated background/turbulence bytes;
- sampling bytes;
- mesh node objects;
- block structures and application block overhead;
- resident particles;
- communication buffers and halo fraction;
- explicit safety margin.

All byte/object coefficients and both fractions are input parameters because
the exact native sizes depend on the AMPS build. The estimate is a planning
gate; native integration tests compare it with actual AMPS counts and peak
memory. Configuration also proves that `maximum_level` can realize the finest
requested cell size. An impossible level cap or an estimate above the memory
budget fails before full allocation.

## R04–R06 runtime controls

The production clock and recurring services are configured with integer
cadences:

| Key | Meaning |
|---|---|
| `application.sep.background_cadence_steps` | request the next background/turbulence generation |
| `application.sep.injection_cadence_steps` | evaluate shock patches and inject a source event |
| `output.cadence_steps` | gather and publish observer products |
| `output.checkpoint_cadence_steps` | write a complete R07 restart; zero disables periodic checkpoints |
| `application.sep.transport.maximum_substeps` | hard accepted-substep limit for one requested AMPS interval |
| `application.sep.transport.maximum_scattering_events_per_substep` | hard event budget inside one focused-scattering substep |
| `application.sep.population_control.cadence_steps` | apply split/merge after injection at this integer boundary cadence |
| `source.shock.physical_particle_rate_per_s` | physical source rate normalized over the injection interval (absent in the derived reduced-front mode) |

Observer blocks (`observer: ID`) additionally accept `kind`, Cartesian velocity, collection
radius, shell radius, minimum/maximum energy, minimum/maximum pitch cosine,
`species = all` or comma-separated subset IDs, and `normalization`. Observer cadence in seconds
must be an exact integer multiple of `run.time_step_s`; normalization therefore
cannot acquire a drifting fractional-tick window. In the derived reduced-front
mode the step is known only after mesh allocation, and observer cadences are
moved up to the first exact multiple of the derived `dt`.

All trajectory, cadence, source, and observer acceptance fields participate in
the physics fingerprint. Checkpoint and output paths remain relocation-only
manifest fields. `R3D04–06` cover clock/event restoration, source cadence
identity/conservation, and resolved observer geometry plus commit-only reset.

## Acceptance evidence

| ID | Contract |
|---|---|
| `CFG3D01` | schema 6, CLI, early errors, layout- and comment-independent configuration identity, D1-D4 rejections (retired schema, srcSEP-only settings, SWCME-owned Parker keys, SWCME `source.*` keys), and dry-run |
| `CFG3D02` | typed groups, field classification, compatibility, and SWMF/analytic factory parity |
| `CFG3D03` | presets, explicit override, photosphere/source-shell ordering, containment, and directional transport-boundary status |
| `CFG3D04` | shared tangent and polarity-independent tube geometry |
| `CFG3D05` | monotone composite profiles, tube scaling, AMR levels, memory categories, and level rejection |
| `CFG3D06` | required complete finite field-line block and fail-closed source consistency |
| `CFG3D07` | complete generated table, mixed signed charges, count/index/symbol/mass/charge failures, observer bounds, and fingerprint identity |
| `CFG3D08` | complete SWCME shock-injection input, D4 Parker derivation (B_r at 1 AU against hand-computed values), missing-field, weight, and per-step-cadence rejection |
| `CFG3D13` | unified generic input: includes, continuation, derived reduced-front numerics keys, `-input`/`--input`/default `amps.in` CLI spellings, immutable commit, and `file:line` diagnostics |
| `CFG3D11` | active corridor/observer, population hysteresis, fixed/local source spectra, mover/coefficient compatibility, and dry-run fields |
| `MSH3D12` | conservative finite Parker-capsule/block intersection |
| `MSH3D13` | exact Parker curve/tangent and arc-length agreement |
| `MSH3D14` | hole-free whole-octree mask and exact AMR-neighbour halo layers |
| `MSH3D15` | fixed photospheric geometry, source-shell separation, and conservative solid-cell classification |
| `MSH3D11` | deterministic unit-labeled initialization Parker Tecplot output |
| `BLDL3D10` | AMPS sphere registration, absorption callback, active-mask composition, and cell-measure ordering |
| `CFG3D18` | shared SEP vocabulary: the schema-6 table's SWCME 3-D key list equals `swcme::input3d::detail::geometry_keys()`, an srcSEP deck stops with a typed "not supported by srcSEP3D" error, every srcSEP3D `run`/`output`/observer table key is understood by the parser, and every shipped example deck parses |
| `R3D08` | canonical source-surface preflight and exact per-species, per-step particle allocation |

Status (2026-10-11): `./test/stage1 --test-group CFG3D` and the full stage-1
suite (154/154) pass after the schema-6 migration; the full native
regressions with the converted decks are pending.  Note that `CFG3D18` is
registered in the stage-1 binary (`test/individual-test/test_configuration.cpp`);
`test/run_tests.py` does not list it by ID, so select it through the
`CFG3D` group of `./test/stage1`.

## Required finite field-line block

Every schema-6 file requires `mesh` / `field_line` (`mesh.field_line.*`) with
these keys (the first eight are shared with srcSEP; `start_mode` is a
srcSEP3D key, and `end_radius_m` is required only by the corner-cube
geometry):

```text
#begin mesh
  #begin field_line
    start_mode = explicit        ! or cme-launch-point (requires physics: swcme)
    origin_x_m = 0
    origin_y_m = 0
    origin_z_m = 0
    initial_x_m = 1.3914e10
    initial_y_m = 0
    initial_z_m = 0
    length_m = 2.0e11
    point_count = 4001
  #end field_line
#end mesh
```

`point_count` includes both endpoints and must be at least two.  `length_m` is
arc length, not final heliocentric radius.  The origin must equal the domain
origin, and the initial point relative to that origin must lie on
`mesh.domain.inner_radius_m` in the direction declared by the tube source
longitude and colatitude.  With `start_mode = cme-launch-point` the initial
point must also equal the SWCME launch apex (`cme.launch_radius` times the
normalized CME direction), and that launch radius must equal
`mesh.domain.inner_radius_m`.  Solar wind speed and rotation are owned by
`physics: swcme` (or, without SWCME, by `physics.background.parker`); the line
block does not duplicate them.  There is no fallback: omitting a required key
is an error (`CFG3D06`).  The finite line is written by
`--initialization-only` to `output.initialization_field_line_tecplot_file`,
and the run summary prints `initialization_field_line=`.

Run only these gates with:

```bash
test/run_tests.py --suite improvements-c --rebuild \
  --output-dir test_output/improvements-c
```
