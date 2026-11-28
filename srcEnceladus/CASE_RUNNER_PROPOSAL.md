# Enceladus measurement/source-basis cases: implementation proposal

Date: 2026-11-13. Status: **implemented; native regression evidence is recorded
in `STATUS.md`.** The numbered sections retain the accepted design contract.
Runtime domain/refinement were added by direct request; see `RUNTIME_INPUT.md`
for the complete key/geometry syntax.

## 1. Intended user interface and the source-data contract

Configure the location of the model data once, compile the Enceladus
application, and select a measurement/source-basis combination at runtime:

```bash
./Config.pl -application=enceladus \
  -model-data-path=/data/vtenishe/cassini_validation_data/studies

# After the required clean/configure/build sequence described in section 16:
mpiexec -n 1 ./amps -i input/enceladus_cases.in \
  -case E3--tenishev2014_hybrid_8

mpiexec -n 1 ./amps -i input/enceladus_cases.in \
  -case uvis_2010_e10_solar--porco2014_98
```

The selector has exactly the requested form:

```text
MEASUREMENT_ID--SOURCE_BASIS_ID
```

**The Monte Carlo model must not select, reconstruct, or interpret a fitting
strategy.** It consumes a complete published source prescription. Independent,
pooled, regularized, or future fitting methods belong to the standalone
application. Their names may occur in provenance, but are not Monte Carlo
configuration keys, selection rules, or physical algorithms.

In particular, there will be no `fit_strategy` CLI option or runtime key.
The standalone publication step associates one chosen source prescription
with each case ID. If several alternatives exist, that association must be
made upstream; the forward model must not guess which one to use.

An absolute source prescription must include every physical source's ID,
geometry, direction, temperature, bulk speed and species production rates.
The existing schema-2 `.in` source format already supplies those quantities.
Source fitting does not replace the forward model's numerical controls,
reaction selection, gravity, sampling, or observational response operator.

The previously authorized shape-file mode remains supported: its fractions
do not define an absolute production rate, so the user must supply each
injected species' total rate in the runtime input. That is an explicit
source-normalization input, not knowledge of the fitting strategy.

## 2. Existing implementation to reuse

| Existing component | Verified behavior | Proposed use |
| --- | --- | --- |
| `input/i3.in` | Common application settings and direct named `case` children | Same input organization for Enceladus |
| `srcHartley2/RadialMain.cpp` | `-i`/`--input`, `-case`/`--case`, MPI runtime validation | Match the familiar CLI spellings |
| `srcHartley2/CometRuntimeInput.cpp` | Calls generic `ResolveCase`; checks canonical input across MPI ranks before applying it | Reuse the pattern without importing comet physics |
| `src/general/input/amps_input.h` | Named overlays, includes, source locations, canonical text | Use the existing grammar; no parallel case parser |
| `Config.pl` | `-model-data-path=PATH` records `MODELINPUTDATA` through the established configuration workflow | No new configuration option needed |
| `ampsConfig.pl` | Generates `PIC::UserModelInputDataPath` from `MODELINPUTDATA` | Supply the configured root to the application resolver |
| `srcEnceladus/ForwardModel.cpp` | Dependency-light source/configuration parsing and numerical planning | Retain physical validators and add a case-aware resolution layer |
| `srcEnceladus/ForwardModelAMPS.cpp` | Native installation, observations, source injection, loss and output | Install an already resolved case without changing source physics |
| Standalone `fit_source_export.py` | Complete absolute/shape `.in` exports with hashes and geometry snapshots | Inputs to an upstream publication step |

The baseline CLI supported explicit `-input`/`--input` and local
`--validate-input`, with schema 2 observation children. The implemented CLI
now also accepts the named-case commands below and schema 3 catalogs.
Case-free schema 2 inputs retain their existing interpretation.

## 3. Case identities

Use existing fitting target IDs wherever possible, avoiding an extra naming
translation:

| Measurement ID | Instrument / event | Present standalone source exports |
| --- | --- | --- |
| `E3` | INMS E3 | H2O absolute/shape; CO2 conditional shape, physical absolute deferred |
| `E5` | INMS E5 | H2O absolute/shape; CO2 conditional shape, physical absolute deferred |
| `E14` | INMS E14 | CO2 absolute/shape |
| `E17` | INMS E17 | CO2 absolute/shape |
| `E18` | INMS E18 | CO2 absolute/shape |
| `uvis_2007_zeta_orionis` | UVIS stellar occultation | H2O absolute/shape |
| `uvis_2010_e10_solar` | UVIS solar occultation | H2O absolute/shape |

The instrument is an explicit catalog field. It is not inferred from a
filename, capitalization, or the presence of `UVIS` in a user string.

Supported basis IDs remain exactly:

| Basis ID | Geometry |
| --- | --- |
| `tenishev2014_hybrid_8` | Eight jets plus seven stripe components |
| `porco2014_98` | 98 catalog jets |
| `portyankina2022_solar_41` | The canonical 41-jet selection with local-normal axes |
| `porco2014_98_plus_e17_candidate` | 98 catalog jets plus the explicitly selected E17 candidate |

The current export inventory allows **28 potential measurement/basis
identities**: seven measurements times four bases. They become runnable only
when a publication supplies the requested physical source data and observation
assets. This is not 124 cases: the 124 fitting results include alternatives
that must be resolved on the standalone side before publication.

Case names are case-sensitive, matching the shared input grammar. Each
component must be nonempty, contain only letters, digits, single hyphens or
underscores, and must not itself contain `--`. Require exactly one `--`.
Reject missing/extra separators and mismatches between the case name, catalog
measurement/basis fields and the source file's `source_case`.

E7/E13/E21 and other observations must not silently reuse E3 source rates.
Having an observation product does not mean a fitted source prescription for
it exists. Such cases require an explicit upstream publication. Any future
transfer experiment must record the source-training/normalization epoch and
the target observation separately; it must not fabricate a target-specific fit.

## 4. Data ownership and publication

The Monte Carlo application reads a **published source set**, not a directory
of competing optimization outputs.

```text
standalone fitting
    -> existing detailed results and alternative source exports
    -> explicit standalone publication selection
    -> versioned, complete source/observation bundle
    -> Monte Carlo case resolution
    -> native source/transport/reaction/observation calculation
```

Add a standalone tool, implemented as `standalone/prepare_forward_cases.py`
(provisionally named `publish_forward_cases.py` in this proposal). It will:

1. Read a publication-selection file identifying the complete exported source
   artifact(s) to publish for each measurement/basis. Selection occurs here,
   outside the Monte Carlo input and executable.
2. Verify export hashes, canonical membership, SI units, frame, surface
   projection, species, source law and source kinematics.
3. Copy or produce one canonical `sources.in` for each case, without refitting,
   replacing coefficients, adding empirical backgrounds, discarding exact zeros
   or choosing an unrequested CO2 fragmentation correction.
4. Package the measurement's qualified observation tables, required trajectory
   or ray tables, column constraints and their existing metadata.
5. Write the case catalog, source/asset hashes, machine-readable provenance and
   a human-readable publication receipt. Preserve upstream method information
   only as provenance.
6. Publish to a new versioned directory after all validation succeeds. Never
   overwrite the raw fitting campaign or modify a publication while it is used.

For initial E3 H2O development, the already prepared explicit-file run uses
the independent E3 export and `Q=3.1424817803383656e28 molecules/s`. That is
an existing example, not a proposed rule for selecting other cases. The
publication tool must require an explicit mapping instead of automatically
choosing independent, all-data, or any other alternative for every case.

### Multiple species

The existing native `Source` representation shares geometry, temperature,
bulk speed and source law, with per-species strengths. Publication can merge
H2O and CO2 files only if these common physical parameters agree source by
source. Compare a canonical physical-geometry/kinematics fingerprint, not
whole JSON snapshot hashes: snapshots legitimately contain different molecular
species metadata even when common source parameters are identical.

If the parameters differ, publication must fail clearly. A future per-species
kinematics extension would require an explicit source-schema/native-model
change; the publisher must not silently use one species' temperature for both.
Likewise, absolute and shape strengths cannot be mixed under the current
file-wide `mode`. A publication must resolve its mode explicitly upstream.

At present, an absolute E3/E5 publication can contain H2O. It cannot acquire
physical CO2 rates by treating apparent mass-44 rates as physical. CO2 remains
unavailable in that absolute bundle until a correction is explicitly supplied
upstream. An explicitly published shape bundle plus user-defined physical
totals remains the authorized alternative. UVIS data do not infer CO2 rates.

## 5. Configured root and portable layout

Use the existing generated `PIC::UserModelInputDataPath` as the anchor.
Recommend configuring a stable parent directory:

```text
Config.pl -model-data-path=/data/vtenishe/cassini_validation_data/studies
```

The runtime input selects a publication **by directory**, not fitting method:

```text
source_set=fitting_compain_111126/published/run_sources_v1
```

This allows new source publications under the same configured root without
recompiling. It is also valid to configure the campaign itself as `PATH` and
use `source_set=published/run_sources_v1`. The path semantics must be identical
and documented on every machine. Changing the compiled root requires the
normal configuration/rebuild; changing the selected case, source-set directory
or runtime controls does not.

Proposed tree, relative to the configured root:

```text
fitting_compain_111126/
  h2o/ ...                         existing fitting results; unchanged
  co2/ ...                         existing fitting results; unchanged
  published/
    run_sources_v1/
      enceladus_cases.in           catalog, schema 1
      publication.json            upstream provenance; not a physical parser
      cases/
        E3--tenishev2014_hybrid_8/
          sources.in              full schema-2 source prescription
          provenance.json
        uvis_2010_e10_solar--porco2014_98/
          sources.in
          provenance.json
      observations/
        E3/e3_h2o_sticking_corrected.csv
        uvis_2010_e10_solar/optical_depth.csv
        uvis_2010_e10_solar/observer_geometry.csv
      trajectories/
        E3/e3_spice_ground_track.csv
```

All operational asset references are bundle-relative. Resolve `source_set` from the
configured root, and catalog assets from the selected publication directory.
Normalize/check paths and reject traversal or symlink escape outside that
publication. If existing raw exports reference external files, the publisher
must package them deliberately rather than relying on the launch directory.
Record both configured and resolved locations in provenance. Historical
absolute paths may remain in copied source-file comments or upstream
provenance; they are descriptive strings, not runtime dependencies. The
source `.in` already embeds its physical geometry/kinematics completely.

Generic runtime `#include` remains unchanged and relative to its including
file. Do not invent `#include model:...` syntax or expand environment variables
inside source filenames. Legacy schema-2 explicit-file path rules stay as they
are. Proposed catalog-relative paths are a separate, unambiguous context.

## 6. Published catalog syntax

Use the existing AMPS grammar for native-readable catalog data. This avoids
making the C++ application parse the standalone optimizer's JSON reports.
The following is a structural example; `PUBLISHER_SHA256` is a placeholder for
the actual 64-hex-digit hash, not a runnable production value.

```text
#begin case_catalog: enceladus
schema_version=1
publication_id=fitting_compain_111126_run_sources_v1
frame=IAU_ENCELADUS
reference_radius_m=252100

#begin case: E3--tenishev2014_hybrid_8
measurement_id=E3
instrument=INMS
source_basis=tenishev2014_hybrid_8
source_distribution_file=cases/E3--tenishev2014_hybrid_8/sources.in
source_sha256=PUBLISHER_SHA256
available_species=H2O
unavailable_source.CO2=Physical absolute rate not published; calibration decision deferred
provenance_file=cases/E3--tenishev2014_hybrid_8/provenance.json
provenance_sha256=PUBLISHER_SHA256
flyby_trajectory_file=trajectories/E3/e3_spice_ground_track.csv
flyby_trajectory_sha256=PUBLISHER_SHA256

#begin observation: inms_e3_h2o
instrument=INMS
event_id=E3
data_file=observations/E3/e3_h2o_sticking_corrected.csv
data_sha256=PUBLISHER_SHA256
time_column=time_ca_s
x_column=x_m
y_column=y_m
z_column=z_m
comparison_species=H2O
value_column=density_m3
uncertainty_column=sigma_m3
#end observation
#end case

#begin case: uvis_2010_e10_solar--porco2014_98
measurement_id=uvis_2010_e10_solar
instrument=UVIS
source_basis=porco2014_98
source_distribution_file=cases/uvis_2010_e10_solar--porco2014_98/sources.in
source_sha256=PUBLISHER_SHA256
available_species=H2O
provenance_file=cases/uvis_2010_e10_solar--porco2014_98/provenance.json
provenance_sha256=PUBLISHER_SHA256

#begin observation: uvis_2010_e10_solar
instrument=UVIS
event_id=uvis_2010_e10_solar
data_file=observations/uvis_2010_e10_solar/optical_depth.csv
data_sha256=PUBLISHER_SHA256
ray_geometry=finite_observer_to_sun
observer_geometry_file=observations/uvis_2010_e10_solar/observer_geometry.csv
observer_geometry_sha256=PUBLISHER_SHA256
published_h2o_column_m-2=9.0e19
published_h2o_column_sigma_m-2=2.3e19
#end observation
#end case
#end case_catalog
```

The source file remains the authority for source coordinates, kinematics,
source law, mode and rates. Catalog identity/species/frame/radius fields are
cross-checks. A contradictory catalog is an error, not an override.
An optional `unavailable_source.SPECIES` reason is availability metadata, not
a correction algorithm. It lets listing/validation explain why H2O is ready
while physical CO2 is not; Monte Carlo cannot resolve that gap itself.
Trajectory fields are optional assets; E3 has a prepared dense trajectory.
Do not manufacture such a trajectory for an event with only observation-point
positions. INMS observation sampling can operate without the legacy dense
flyby output.

Preserve source geometry schema 2. Do not embed particle counts, reaction
switches, mesh settings, or new emission laws in a fitted-source file.

## 7. Runtime input with I3-style cases

Add `input/enceladus_cases.in` using application schema 3. Keep the existing
schema-2 `input/enceladus.in` and explicit E3 H2O example working unchanged.
This example illustrates complete per-case species controls. Its publication
paths are proposed; the source set does not exist merely because it appears
in this document. Numerical values retain the earlier user controls and are
not a new convergence recommendation for either instrument.

```text
#begin application: enceladus
schema_version=3
source_set=fitting_compain_111126/published/run_sources_v1
mesh.domain_min_m=-1260500,-1260500,-1260500
mesh.domain_max_m=1260500,1260500,1260500
mesh.background_cell_size_m=600000
mesh.surface_cell_size_north_m=300000
mesh.surface_cell_size_south_m=150000
reaction_products=discard
solar_activity=active
plasma_background_mode=none
sample_length_steps=1000
total_iterations=10000
volume_field_output=false
particle_splitting_merging=off
random_seed=184467

#begin case: E3--tenishev2014_hybrid_8
simulated_species=H2O
injected_species=H2O
chemistry_species=H2O
selected_reactions=burger2007.photo.h2o_to_oh_h,burger2007.photo.h2o_to_h2_o,burger2007.photo.h2o_to_h_h_o,burger2007.photo.h2o_to_h2op_e,burger2007.photo.h2o_to_h_ohp_e,burger2007.photo.h2o_to_h2_op_e,burger2007.photo.h2o_to_oh_hp_e
primary_macroparticles_per_step.H2O=100
characteristic_speed.H2O_m_s-1=300
courant.H2O=0.5
flyby_density_output=true
flyby_mesh_refinement=true
flyby_refinement_target_cell_size_m=20000
flyby_refinement_radius_m=150000
flyby_refinement_time_min_s=-120
flyby_refinement_time_max_s=100
#end case

#begin case: uvis_2010_e10_solar--porco2014_98
simulated_species=H2O
injected_species=H2O
chemistry_species=H2O
selected_reactions=burger2007.photo.h2o_to_oh_h,burger2007.photo.h2o_to_h2_o,burger2007.photo.h2o_to_h_h_o,burger2007.photo.h2o_to_h2op_e,burger2007.photo.h2o_to_h_ohp_e,burger2007.photo.h2o_to_h2_op_e,burger2007.photo.h2o_to_oh_hp_e
primary_macroparticles_per_step.H2O=100
characteristic_speed.H2O_m_s-1=300
courant.H2O=0.5
flyby_density_output=false
flyby_mesh_refinement=false
#end case

#begin case: E14--tenishev2014_hybrid_8
simulated_species=CO2
injected_species=CO2
chemistry_species=none
selected_reactions=none
primary_macroparticles_per_step.CO2=100
characteristic_speed.CO2_m_s-1=300
courant.CO2=0.5
flyby_density_output=false
flyby_mesh_refinement=false
#end case
#end application
```

Common assignments are overlaid by the selected direct case using the generic
`ResolveCase` implementation. Do not put species-specific numerical controls
or measurement observation blocks in the common section when cases select
different species/measurements: generic overlays replace matching entries,
but do not delete unrelated entries. Otherwise a CO2 case could inherit
inactive H2O keys, or an E3 case could retain an E5 observer.

The publisher may additionally generate a runtime-case template from a
user-supplied numerical-control template. New observations then need a
publication entry and a regenerated/edited case block, not C++ changes.
Generation must not guess particle counts, characteristic speeds, chemistry
or missing trajectory/refinement parameters from the fitting results.

Common `schema_version` and `source_set` are application-level identities;
do not permit a case to change them. Runtime controls may vary by case. A
later source-set version is selected by editing the common `source_set`, not
by teaching Monte Carlo about fitting methods.

In schema-3 case mode, source files and primary observation data come from
the catalog. Reject a competing runtime `source_distribution_file`,
`observation_data_root`, trajectory pathname or observation block. Legacy
explicit-file mode retains those capabilities. Additional diagnostic
observations can be designed later with explicit ownership; do not silently
combine two competing sources of scientific case data in the first version.

If the published source mode is `shape`, the relevant case must explicitly
add every `total_source_rate.SPECIES_s-1`. For `absolute`, reject those keys
exactly as today. The runtime does not choose another mode or publication if
the requested species/rate information is unavailable.

## 8. CLI specification

| CLI | Proposed behavior |
| --- | --- |
| `-i FILE`, `-input FILE`, `--input FILE`, `--input=FILE` | One runtime input path; retain old spellings and add I3's short alias |
| `-case ID`, `--case ID`, `--case=ID` | Select one direct application case, exactly `measurement--basis` |
| `--list-cases` | Local inventory of declared cases and catalog availability, without allocating MPI/mesh/particles |
| `--describe-case -case ID` | Resolve/validate the selected data and print its full receipt without starting the simulation |
| `--validate-input -case ID` | Existing local pre-MPI validation extended through case/catalog/source/observation validation |
| `--validate-runtime-input -case ID` | Collective MPI agreement/installation validation, then exit before mesh/particle allocation |
| `--help` | Explain the two input layers, case syntax and model-data root |

Duplicate selectors/options, missing values, malformed IDs and unknown options
are errors. Case blocks require an explicit selector; never select the first
case. A case selector applied to a legacy case-free deck is an error.
`--list-cases` is exclusive with case/describe/validation/run actions.

For listing, distinguish the declared runtime cases from publications and
show `READY`, `MISSING_PUBLICATION`, `MISSING_DATA`, `DEFERRED_SOURCE`, or
`INVALID` with reasons. `READY` means input-ready, not physically qualified.
List available species and source mode. Do not turn a missing publication into
an equal-rate startup case. Unselected missing data may be reported without
preventing a selected, complete case from running.

`--describe-case` reports Q, source counts, kinematic ranges, active species,
reaction IDs, observation operators, numerical inputs and asset hashes.
Actual timestep/particle weights remain **pending AMR realization** until the
minimum allocated cell size exists; do not substitute the requested mesh
width for the realized cell size to print a falsely final weight.

## 9. Resolution and validation algorithm

1. Parse CLI syntax and runtime input with the generic AMPS parser. Require
   exactly the owned Enceladus application block and reject unknown structure.
2. Check all case declarations structurally, including duplicate/unknown keys
   and nested-case rejection. Resolve the selected case with `ResolveCase`.
   Validate selected numerical values and active-species controls; do not load
   every unselected case's external datasets.
3. Obtain the configured data root from the thin PIC adapter and resolve
   `source_set/enceladus_cases.in`. Reject empty/missing/invalid roots in new
   case mode; legacy explicit-file mode does not acquire this requirement.
4. Find exactly the matching catalog case. Cross-check the two ID components
   and declared measurement/basis; reject duplicate or contradictory entries.
5. Read selected asset bytes, verify their hashes, and retain those exact
   bytes for parsing and snapshotting. Do not verify one file version and then
   reopen a mutable pathname for a different version.
6. Parse the existing complete source file and run its unchanged membership,
   geometry, species, frame, law, surface and rate validators. Preserve exact
   zeros and representable subnormal strengths.
7. Select observation requests for the declared measurement and active species.
   Reject a requested comparison lacking qualified data; preserve explicitly
   labelled context-only observations as context, never score them as density.
8. Resolve optional dense trajectory assets only when dense flyby output is
   requested. Translate the prepared case into existing `Configuration`,
   `SourceDistribution` and observation structures, without fitting logic.
9. Run existing Enceladus physical/numerical checks, including compiled species,
   source mode, inactive controls, chemistry and surface-radius agreement.
10. For an MPI run, establish all-rank success and content agreement before
    installing state/callbacks. Then construct the mesh and evaluate the
    actual timestep/weights through existing planning functions.

Selected unknown keys, missing assets, incompatible species, malformed
observations, source/basis mismatches and bad hashes are errors with file,
line, selected case and logical asset path. Do not catch them and substitute
another physical case. Error messages should list available case IDs or
species where useful.

## 10. MPI startup and state ownership

Preserve a dependency-free resolver. Pass a plain data-root string into it;
do not include PIC/MPI headers in `ForwardModel.h` or a new catalog library.

Split native installation into preparation and commit. Preparation returns
an immutable candidate or an error without setting global installed state.
Local inspection modes remain genuinely pre-MPI. During an actual case run,
every rank must participate in the same validation collectives even if its
local data read failed. No rank may exit into a different startup phase while
others enter mesh collectives.

`PIC::InitMPI()` establishes ranks/communicators without constructing the AMR
mesh or particle buffer. The case-aware path can initialize that communicator,
reduce preparation status and compare content fingerprints, then commit the
resolved state and continue the existing initialization sequence. Keep the
legacy explicit-file behavior intact unless a separate, tested migration is
approved.

Fingerprint the canonical effective configuration, selected source bytes,
observation/trajectory bytes and publication identity. Comparing only the
runtime input text misses differing files behind identical paths. Keep local
absolute paths as provenance, but do not make machine-specific mount names
part of the physical-content identity. A collective validation mode must test
both missing-file and rank-divergent-file failures without hanging.

For SHA-256, there is an existing dependency-free implementation in
`src/models/sep_coronal_cme/src/sha256.*`, currently in an internal CME
namespace. Avoid linking Enceladus physics to CME internals. Introduce a
small generic input-digest helper using that reviewed algorithm with recorded
provenance and independent standard vectors; leave existing CME callers and
files unchanged in this task. No runtime Python/shell hash subprocess or new
mandatory OpenSSL dependency is needed.

## 11. Physical and numerical behavior to preserve

Case resolution must not change the source law, emitted thermal/drift
distribution, source geometry, gravity, surface absorption, photolytic rates,
splitting/merging algorithms, or generic PIC behavior.

For each active species s, preserve the current planning equations:

```text
dt = min_s [ Courant_s * realized_minimum_cell_size / characteristic_speed_s ]
Q_s = sum_i Q_i,s                        (absolute source mode)
Q_s = user_total_source_rate_s           (shape mode)
W_s = Q_s * dt / requested_primaries_per_step_s
```

Units are seconds, metres, m/s, molecules/s and molecules per model particle.
The requested primary count is the exact global integer count per step, not
a separate count on every MPI rank. Characteristic speed is a numerical
timestep input and does not replace the fitted bulk/thermal velocities.

Retain neutral H2O/CO2 transport and selected H2O parent-loss photolysis with
discarded products. Do not add plasma, daughter production or new CO2
photochemistry as part of this case-selection feature.

Empirical fitting backgrounds are not physical surface injection. They stay
in upstream provenance/reference results. A computed Monte Carlo density
therefore need not equal a fitted curve that includes an empirical background.
The steady, gravity-free standalone model also differs from the native model's
transport/loss physics; profile agreement is a separate scientific validation.

## 12. INMS and UVIS observation handling

For INMS, retain qualified number-density tables with their exact species,
coordinates, time reference, uncertainty columns and calibration metadata.
Dense flyby output is optional and requires a real frozen trajectory; it is
not fabricated from a source case name. Report domain masks/cell sizes and
exclude out-of-domain placeholders from comparison statistics.

For UVIS, use the published ray/observer tables and the declared operator.
Preserve the existing 2007 tangent-track-normal chord and 2010 finite
observer-to-Sun contracts and their stated validity limits. Output raw H2O
column density in m^-2. Optical-depth shape/response scaling and a published
absolute column constraint remain distinct comparisons; do not relabel an
optical-depth proxy as a separately calibrated column in every row.

Preserve the inherited standalone 2007/2010 geometry approximations in
provenance. In particular, native and fitted-reference ray operators may
differ for the documented 308 outward 2010 integrations. A case selector
cannot make those operators equivalent or close a physical validation gate.

The first implementation reuses the existing AMR and column-integration
operators. INMS E3 can retain its existing trajectory tube. UVIS cases must
not be silently routed through an E3 tube or assigned a nonexistent flyby.
Explicit runtime boxes, balls, shells, hemisphere surface requests and a
south-polar cone are implemented alongside the existing E3 trajectory tube.
Their conservative block-intersection criterion and noncubic aspect conversion
are tested independently and natively. Automatic extraction of refinement
tubes from arbitrary INMS/UVIS rays remains a future feature. Mesh/timestep/
particle convergence is open; a case name alone implies no resolved density
or column profile.

A future measurement-group ID must identify one physically consistent source
state. Pooled fitting is not permission to switch the total rate separately
inside each observation callback. Batch execution of separate cases is often
the appropriate operation for different source normalizations.

## 13. Output and provenance

Use a fresh run directory for each invocation. Retain existing Tecplot and
PNG comparison filenames so current postprocessors continue working. Add
case/publication IDs to file headers and a new `run_manifest.json` containing:

- Original CLI, runtime input/includes and selected case ID.
- Measurement, instrument, source basis and publication identity.
- Configured/resolved data paths and hashes of all consumed assets.
- An exact source snapshot and canonical resolved runtime input.
- Source totals, available/injected species and source-mode normalization.
- Seed, reaction IDs, numerical controls, realized mesh minimum, timestep,
  weights, rank count and executable/source build identity.
- Observation operator/calibration metadata, domain masks and comparison scope.
- Completion status, per-species particle budgets and output paths.

Method provenance from the standalone publisher can be retained verbatim or
referenced by hash. It must not drive native configuration. Record the chosen
bundle, not a misleading inferred fitting-strategy label.

## 14. Concrete implementation units

| Location | Responsibility |
| --- | --- |
| `srcEnceladus/standalone/prepare_forward_cases.py` (new; proposed as `publish_forward_cases.py`) | Explicit upstream publication from existing full exports and qualified observation assets |
| `srcEnceladus/CaseSelection.h/.cpp` (new) | Case identity, generic overlay integration and schema-3 selection, no PIC/MPI |
| `srcEnceladus/SourceCatalog.h/.cpp` (new) | Catalog/asset resolution and immutable source-data preparation, no optimizer knowledge |
| `src/general/input/amps_input_digest.h` (proposed new helper) | Portable content digests, independently tested; no physical-model dependency |
| `srcEnceladus/ForwardModel.h/.cpp` | Preserve legacy entry point; add explicit case-resolution context/result and translate to existing validators |
| `srcEnceladus/ForwardModelAMPS.h/.cpp` | Supply generated data-root string, prepare/commit native state and MPI agreement |
| `srcEnceladus/main.cpp` | CLI aliases, case/list/describe/validation modes and safe startup ordering |
| `input/enceladus_cases.in` (new) | Common controls and generated/maintained direct named cases |
| `srcEnceladus/test/` and standalone tests | Independent fixtures, source publication, CLI, native and MPI receipts |
| Global test registries/catalog/manual snapshot | Register actual new test identities and documentation once implemented |

`Config.pl` and `ampsConfig.pl` already support the required option; do not
rewrite their option handling unnecessarily. Verify the generated path and
its C++ string-length/escaping behavior, including a path containing spaces.
Any required quoting/length fix must be a small separately tested generic
change. Never edit `.amps.conf` directly.

Do not change `SourceDistribution::sourceCase` into a fitting-method selector.
It continues to identify the physical source basis; the new runtime case ID
combines that basis with the measurement.

## 15. Implementation order and acceptance gates

| Stage | Work | Acceptance evidence |
| --- | --- | --- |
| C0: freeze interfaces | Review this proposal, source publication ownership and exact catalog/schema fields | Examples and contracts agree; no fitting-strategy key in forward configuration |
| C1: upstream publication | Publisher, selection file, complete merged sources, catalog and hash inventory | Nonuniform synthetic rates/zeros preserved; no background injection or guessed physical CO2; deterministic repeat output |
| C2: dependency-free resolution | ID parser, generic overlays, root/source-set paths, selected catalog/assets | Independent parser/path/species/mode/hash/error fixtures; legacy parser cases pass unchanged |
| C3: native CLI | Aliases, list/describe/local validation, prepare/commit and receipts | Native executable resolves both INMS and UVIS cases; no allocation in inspection modes |
| C4: MPI agreement | Collective failure handling and content identity | One/four-rank validation; divergent/missing rank-local assets fail promptly without deadlock |
| C5: integration campaign | H2O/CO2 source bases, INMS sampling, both UVIS operators and isolated output | Existing gates remain; bounded native runs close particle budgets; independent observation fixtures pass |
| C6: documentation/release | READMEs, source publication guide, catalog/table updates and migration examples | Exact commands/known limitations/evidence paths; proposed behavior no longer advertised as unimplemented only after evidence exists |

Source publication for all current cases requires the upstream scientific
selection to be supplied. A missing choice is a data-publication gap, not a
reason to encode a default fitting strategy in Monte Carlo.

## 16. Build, validation and runtime commands after implementation

The new input parser/CLI and source catalog code require one native rebuild.
Follow the repository preflight exactly; the following is a documented future
sequence, not an instruction to rebuild during this proposal task:

```bash
cd /home/vtenishe/Mars2/AMPS
# Confirm this root and no active build/test process before proceeding.
rm -rf -- build
./Config.pl -application=enceladus \
  -model-data-path=/data/vtenishe/cassini_validation_data/studies
./ampsConfig.pl -input enceladus.input -no-compile
env MAKEFLAGS="-j16" make amps
```

Preserve site configuration/libraries and old outputs. The Enceladus source
tree has a large generated `test_output` directory; use the documented safe
source-generation procedure rather than recursively copying it into a native
build or deleting it. No Git-state-changing operation is part of this work.

After a source set has actually been published:

```bash
./amps -i input/enceladus_cases.in --list-cases
./amps -i input/enceladus_cases.in \
  -case E3--tenishev2014_hybrid_8 --describe-case
./amps -i input/enceladus_cases.in \
  -case E3--tenishev2014_hybrid_8 --validate-input
mpiexec -n 4 ./amps -i input/enceladus_cases.in \
  -case E3--tenishev2014_hybrid_8 --validate-runtime-input
```

For a real simulation, first create a fresh output directory and run the
absolute executable/input paths there. Changing only `-case` selects another
published measurement/basis with the same executable. Publishing a new
source-set directory under the configured root also requires no rebuild.

## 17. Required tests and regression boundaries

Use synthetic independently specified sources/observations for correctness;
real fitting campaigns are data-compatibility evidence, not expected-value
generators. Implemented checks and remaining scientific boundaries:

| Group | Required checks |
| --- | --- |
| CLI/identity | Every supported spelling, duplicate/missing flags, both instruments, all bases, malformed separator/name, missing/unknown selector, case-free legacy rejection |
| Overlay ownership | Common/selected precedence, source locations, includes, no nested cases, no inherited inactive-species controls or unrelated observations, unknown keys rejected |
| Config/root portability | Generated root honored from an unrelated working directory; two machine roots with identical contents give identical physical receipts; missing/overlong/quoted paths diagnosed |
| Publication | Explicit mapping required, exact geometry/kinematics and rates, zero/subnormal preservation, species merge compatibility, no background terms, no inferred fragmentation factor, atomic new publication |
| Catalog/integrity | Wrong case/basis/frame/radius/species, duplicate entries, modified source/observer bytes, traversal/symlink escape and incomplete source vectors fail |
| Source modes | Absolute Q and weights, shape plus explicit totals, mixed-mode conflicts, missing/extra total controls, unavailable physical CO2 and unsupported UVIS CO2 requests |
| Source physics | Existing moment/quantile, geometry and neutral parent-loss tests remain unchanged; case selection cannot alter their reference values |
| INMS | Independent constant/analytic field at fixed points, exact time/units mapping, valid domain masks and current plotting compatibility |
| UVIS | Independent column of a constant field over known clipped chord/finite ray, correct orientation and path length, occultation/domain masks and distinct column/optical-depth meanings |
| MPI | Same contents agree on one/four ranks; one rank with missing or changed assets fails collectively; inspection mode allocates no particles; global injection count independent of rank count |
| Native budget | Bounded runs for all four bases with H2O/CO2 where supported; exact injected/remaining/lost/escaped accounting and output cadence |
| Legacy | Existing 43-result application campaign, six startup gates, explicit source paths and existing default deck retain behavior; do not count previous smoke evidence as new-case evidence |

Retain per-case logs and references, and print exact test IDs. Register new
actual tests in `src/test/case_registries.json` and traceability as appropriate;
regenerate `src/test/TEST_CATALOG.md` and `doc/TestSuiteCatalog.tex` in the same
implementation change. Extend core-input coverage for the digest helper if
introduced. Do not invent a fixed test count before implementation or weaken
existing tolerances to accommodate changed case selection.

Scientific validation remains separate: study cell size, particle number,
seed, sampling duration, losses and instrument response. A successful catalog
load, native smoke, fitted residual or source-file integrity check alone is
not proof of converged agreement with Cassini.

## 18. Decisions fixed versus still upstream

Fixed by the user and this interface proposal:

- CLI identity is `measurement--base`.
- The model-data anchor is configured through `Config.pl -model-data-path`.
- Forward physics consumes complete published sources and does not know the
  fitting strategy. Runtime particle/reaction/numerical controls stay explicit.
- Preserve existing physical constraints and the prior E3/E5 CO2 deferral.

Still an upstream scientific/publication choice, not a missing Monte Carlo
algorithm: which already fitted source vector is published for each case and
which physical CO2 correction, if any, is approved later. The publisher must
require those choices explicitly. Additional observers/bases/species require
qualified data or corresponding model support; do not populate gaps with
equal source rates, assumed compositions or guessed trajectories.

This proposal changes no application physics, native code, fitted products,
build configuration or Git metadata. Implementation should begin only after
the user requests it.

## 19. Source references for implementation review

- [I3 case guide](../srcHartley2/i3/README.md) and
  [I3 input](../input/i3.in).
- [Generic runtime-input specification](../src/INPUT_FILE.md) and
  [generic parser/ResolveCase](../src/general/input/amps_input.h).
- [Configuration option](../Config.pl) and
  [native configuration generator](../ampsConfig.pl).
- [Current Enceladus CLI](main.cpp), [dependency-light model](ForwardModel.cpp)
  and [native installation](ForwardModelAMPS.cpp).
- [Source handoff format](SOURCE_DISTRIBUTION_FILE.md),
  [runtime documentation](RUNTIME_INPUT.md) and
  [standalone exporter](standalone/fit_source_export.py).


## 20. Implemented evidence (2026-11-13)

The accepted C0–C6 interface/implementation sequence is implemented. This is
an engineering acceptance record; scientific convergence remains separate.
Runtime mesh additions were explicitly authorized after the original proposal.

| Stage | Evidence |
|---|---|
| C0 | Complete runtime/catalog/source ownership contracts; no fitting-strategy key; schema 2 compatibility retained. |
| C1 | Four independent publisher test methods preserve all 28 synthetic IDs/rate vectors, zero entries, compatible species merge and CO2 deferral. |
| C2 | Nine actual C++ resolver-driven test methods plus geometry/AMR-limit fixtures and original reference checks. |
| C3 | Native list/describe/local validation of all 28 synthetic cases; final real inventory correctly reports only one READY prescription. |
| C4 | One/four-rank MPI content agreement with zero particle buffer; inconsistent or missing rank inputs and occupied output directories fail without deadlock/overwrite. |
| C5 | Nine bounded native transport runs across all four bases and both species; both UVIS operators; noncubic runtime box/ball/shell refinement; exact source weights and global particle ledgers. |
| C6 | Updated runtime/source/publisher/test READMEs, detailed STATUS ledger, case registry/traceability, generated catalog/manual snapshot, and real-root executable restoration. |

Final combined focused validation: **77 PASS / 0 FAIL / 0 SKIP / 0 ERROR**
across four suites. Full standalone discovery separately passes **66 tests**.
Exact logs and superseded failed attempts are in `STATUS.md`; the clean
combined report is:

`test_output/enceladus-case-runtime-20261113/global-results/run-20261113T092040Z-2290774/amps-core-tests.json`

Only the previously explicit independent E3 H2O hybrid prescription has been
packaged as real initial data. The other real case declarations deliberately
remain unavailable until their complete upstream source prescriptions are
selected. The 28-case synthetic bundle is never a fallback for real runs.
No new long-run Cassini agreement/convergence or absolute E3/E5 physical CO2
qualification is claimed.
