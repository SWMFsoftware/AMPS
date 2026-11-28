# Enceladus forward-model runtime input

The `measurement--basis` CLI, portable source-data catalog and runtime mesh
controls are implemented. The design/evidence contract is maintained in
[CASE_RUNNER_PROPOSAL.md](CASE_RUNNER_PROPOSAL.md). Schema 2 explicit-file
runs remain supported; schema 3 adds named cases. Neither interface asks the
forward model to select or interpret a fitting strategy.

## Case input and commands

`input/enceladus_cases.in` declares all 28 identities: seven measurements
(E3, E5, E14, E17, E18, UVIS 2007 and UVIS 2010) and four source bases.
A declaration is not a claim that its fitted source prescription is ready.
Currently the local real-data bundle contains only the explicitly selected
independent E3 H2O hybrid fit. Other identities report `MISSING_PUBLICATION`;
synthetic tests exercise all 28. E3/E5 absolute physical CO2 remains deferred.

The data root is compiled through the existing configuration option:

```bash
# At repository root, with no other build/test running:
pwd
rm -rf -- build
./Config.pl -application=enceladus -model-data-path=/data/vtenishe/cassini_validation_data/studies
./ampsConfig.pl -input enceladus.input -no-compile
env MAKEFLAGS=-j16 make amps
./amps -i input/enceladus_cases.in --list-cases
./amps -i input/enceladus_cases.in -case E3--tenishev2014_hybrid_8 --describe-case
./amps -i input/enceladus_cases.in -case E3--tenishev2014_hybrid_8 --validate-input
mpiexec -n 4 ./amps -i input/enceladus_cases.in -case E3--tenishev2014_hybrid_8 --validate-runtime-input
```

`Config.pl` copies the application tree: keep generated test output outside
that tree during configuration if it is large. The native test owner holds
and restores `srcEnceladus/test_output` automatically. Do not edit `.amps.conf`.
Local inspection modes do not initialize MPI or build a mesh. Collective
validation initializes MPI, verifies identical input/asset/executable hashes,
installs the source/observation readers, and exits without a particle buffer.
Only an actual run computes final timestep and particle weights from allocated
cells. `--describe-case` labels those values `pending_realized_AMR`.

Start a production case in a fresh output directory:

```bash
run_dir=$(mktemp -d "$PWD/test_output/e3_case.XXXXXX")
mkdir -p "$run_dir/PT/plots"
cd "$run_dir"
mpiexec -n 1 /home/vtenishe/Mars2/AMPS/amps \
  -i /home/vtenishe/Mars2/AMPS/input/enceladus_cases.in \
  -case E3--tenishev2014_hybrid_8
```

The example retains 10,000 iterations, 1,000-step windows, 100 primary
macroparticles/step, 300 m/s characteristic speed and Courant 0.5. These
numerical settings do not establish convergence. Fitted gas temperatures,
bulk velocities, axes and rates come from the complete source file.

Aliases `-i`, `-input`, `--input`, `--input=FILE` and `-case`, `--case`,
`--case=ID` are accepted. Duplicate options, malformed IDs, unknown arguments,
missing selectors and unavailable data are errors. No unselected case or
legacy source table is substituted. List/describe/local validate/MPI validate
are mutually exclusive. A case run refuses an existing run manifest, mesh
file, or populated `PT/plots` directory. Legacy flat runs retain their old
startup validation ordering and source-path rules.

## Runtime domain and refinement (schema 2 or 3)

All coordinates are SI metres in body-fixed `IAU_ENCELADUS`; the box must
strictly enclose the complete 252100 m reference sphere on every axis.
A cell-size target is its **maximum edge length**, not its diagonal. AMR
bisects blocks, so actual edges are normally smaller than the request.
Explicit mesh controls bypass the old `TESTMODE` mesh multiplier; absence
of every new key preserves the legacy resolution/domain law. Schema 3
requires the three mandatory keys below. Schema 2 can opt in to the same law.

```text
#begin application: enceladus
schema_version=3
source_set=enceladus_sources_20261113
mesh.domain_min_m=-1260500,-1260500,-1260500
mesh.domain_max_m=1260500,1260500,1260500
mesh.background_cell_size_m=600000
! Optional surface pair; north is z>0, south z<=0.
mesh.surface_cell_size_north_m=300000
mesh.surface_cell_size_south_m=150000
! Optional infinite south-polar cone pair: angle from the -z axis.
mesh.south_polar_half_angle_deg=45
mesh.south_polar_cell_size_m=300000
! Common physical/numerical keys omitted here; see the complete input file.
#begin refinement: instrument_region
geometry=box
min_m=350000,-30000,-230000
max_m=450000,30000,-170000
cell_size_m=20000
#end refinement
#begin case: E3--tenishev2014_hybrid_8
! Scalars override common keys; identically named refinements overlay by name.
! Case-specific controls are shown in input/enceladus_cases.in.
#end case
#end application
```

Repeat named `refinement` blocks at common or case level. Supported geometries:

| `geometry` | Required geometry keys | Meaning |
|---|---|---|
| `box` | `min_m`, `max_m` (three components each) | Closed axis-aligned volume |
| `sphere` | `center_m`, `radius_m` | Closed ball |
| `shell` | `center_m`, `inner_radius_m`, `outer_radius_m` | Closed spherical shell; inner radius may be zero |

Every region also requires positive finite `cell_size_m`. Radius/extent
ordering, unknown keys, incomplete optional pairs, disjoint regions and
nonfinite values are rejected. Overlaps take the finest request, including
background, surface, polar cone, explicit regions and the existing E3 flyby
tube. Regions may extend outside the domain; only intersecting AMR blocks
exist. A conservative block/region intersection test catches small off-centre
features even when none of the old centre/corner probes lies inside them.
The tube uses exact segment-to-box distance. Noncubic domains use their actual
cell aspect ratio when converting the core's diagonal resolution metric.

The compiled block dimensions, ghost cells and maximum AMR level remain build
capabilities. Input inspection rejects targets below that level's attainable
maximum edge; no silent target clipping is accepted. Final
`dt=min_species(Courant * smallest_allocated_cell_edge / characteristic_speed)`
and `W_species=Q_species*dt/N_species` use the actual mesh and active primaries.
A smaller refinement can therefore reduce every species' global timestep and
statistical weight. Surface tessellation and domain changes can increase memory;
these controls do not qualify spatial or statistical convergence.

## Local source-data bundle

The runtime `source_set` is a relative directory below the configured data
root. It contains `enceladus_cases.in`, complete source files, observations,
optional observer/trajectory tables, SHA-256 inventory and opaque provenance.
Absolute/traversing asset paths and escaping symlinks are rejected. Asset bytes
are checked before use and retained, so native readers consume the verified
snapshot rather than reopen changing data. A portable content fingerprint
excludes installation-specific absolute paths; actual MPI runs add the
executable hash and require all ranks to agree before mesh allocation.

Prepare a new bundle from an explicit upstream selection; existing directories
are never overwritten:

```bash
python3 srcEnceladus/standalone/prepare_forward_cases.py \
  --selection /path/to/selection.json \
  --output-dir /configured/model-data-root/enceladus_sources_VERSION
```

The selection syntax and complete catalog contract are documented in the
[standalone README](standalone/README.md#prepare-forward-case-data).
`mode=absolute` carries every rate and forbids runtime total-rate overrides;
`mode=shape` requires explicit `total_source_rate.SPECIES_s-1` in the run input.
Multiple species can share a case only with identical geometry, temperature,
bulk speed, emission law and source mode. No fitted nuisance background is
injected. New observation IDs require explicit catalog data and native-qualified
observation/ray schemas; a name alone creates no trajectory or calibration.

Every actual run writes `run_manifest.json` plus `case_inputs/`: selected
canonical input, source bytes, checksummed asset copies, case/publication
identity, binary hash, rank count, realized minimum cell edge, limiting species,
global timestep, weights, completed iterations and phase. A short synthetic
transport PASS validates these interfaces, not observational agreement.

## Run E3 H2O from the independent hybrid fit

`input/enceladus_e3_h2o_fitted.in` is a ready-to-run, pre-MPI-validated
single-species example for the independent E3 H2O fit in
`fitting_compain_111126`. It reads the hybrid 8+7 absolute source file directly:
`Q=3.1424817803383656e28 molecules/s`. It retains the existing 10,000 steps,
1,000-step sampling, 100 global injections/step, active-Sun parent-loss
photolysis, no plasma/daughters/splitting, and E3 trajectory refinement and
observed-density comparison. Its 300 m/s characteristic speed and Courant
0.5 are retained user numerical inputs; source gas velocities still come
from the fitted file. Neither the speed control nor this example claims
long-run numerical convergence. Other fitting strategies have separate
files and must be selected explicitly with `source_distribution_file`.

With the current configured Enceladus executable, no rebuild is needed:

```bash
cd /home/vtenishe/Mars2/AMPS
./amps --validate-input --input input/enceladus_e3_h2o_fitted.in
run_dir=$(mktemp -d "$PWD/test_output/e3_h2o_fitted.XXXXXX")
mkdir -p "$run_dir/PT/plots"
cd "$run_dir"
mpiexec -n 1 /home/vtenishe/Mars2/AMPS/amps \
  --input /home/vtenishe/Mars2/AMPS/input/enceladus_e3_h2o_fitted.in \
  > run.log 2>&1
```

Each invocation uses a fresh output directory. The current working directory
owns `run.log` and `PT/plots/`, including the dense E3 flyby Tecplot file and
the named `inms_e3_h2o` observation product. The application-global E3 test
wrapper assumes both H2O and CO2; use this direct command for the H2O-only
example. Its full 10,000-step simulation has not been executed merely by
validating the new input deck.

`input/enceladus.input` is read only by the AMPS build/configuration workflow.
It selects compiled capabilities and the union of available species; it is not
a run deck. Every executable invocation reads `input/enceladus.in` by default,
before MPI initialization. The exact `srcSEP`-compatible overrides are
`-input FILE`, `--input FILE`, and `--input=FILE`.

The checked-in `input/enceladus.in` is the requested longer E3 run deck.
It selects the temporary equal-rate H2O+CO2 receipt, active-Sun H2O neutral
parent loss, no plasma, 100 primaries per species per step, 300 m/s,
Courant 0.5, ten 1,000-step sampling windows, seed 184467, and no split/merge.
The E5 source receipt remains an explicit one-line alternative, but changing
flyby or any numerical control creates a distinct run configuration.

New standalone fits export separate source receipts automatically for hybrid
8+7, Porco 98, Portyankina 41 and explicit Porco 99 geometries. Select the
desired target/strategy `.in` from its `amps_sources/index.json` with
`source_distribution_file`; no fit changes this default run deck. Absolute
mode supplies physical rates directly. Shape mode additionally requires the
per-species `total_source_rate.SPECIES_s-1` inputs. E3/E5 apparent CO2 has no
automatic physical absolute file until a fragmentation correction is chosen.
The exact handoff and geometry rules are in
[`SOURCE_DISTRIBUTION_FILE.md`](SOURCE_DISTRIBUTION_FILE.md).

The forward model reads one explicit `application: enceladus` block. It uses
the repository-wide grammar in `src/INPUT_FILE.md`; unknown or duplicated keys
are fatal.

```text
#begin application: enceladus
schema_version=2
source_distribution_file=enceladus_e3_equal_h2o_co2_sources.in
! For E5: source_distribution_file=enceladus_e5_equal_h2o_co2_sources.in

simulated_species=H2O,CO2
injected_species=H2O,CO2
chemistry_species=H2O
reaction_products=discard

selected_reactions=burger2007.photo.h2o_to_oh_h,burger2007.photo.h2o_to_h2_o,burger2007.photo.h2o_to_h_h_o,burger2007.photo.h2o_to_h2op_e,burger2007.photo.h2o_to_h_ohp_e,burger2007.photo.h2o_to_h2_op_e,burger2007.photo.h2o_to_oh_hp_e
solar_activity=active
plasma_background_mode=none
primary_macroparticles_per_step.H2O=100
primary_macroparticles_per_step.CO2=100

characteristic_speed.H2O_m_s-1=300
characteristic_speed.CO2_m_s-1=300
courant.H2O=0.5
courant.CO2=0.5
...one speed and Courant number for every simulated species...

sample_length_steps=1000
total_iterations=10000
volume_field_output=false
flyby_density_output=true
flyby_id=E3
flyby_trajectory_file=../srcEnceladus/standalone/figure_data/e3_spice_ground_track.csv
flyby_mesh_refinement=true
flyby_refinement_target_cell_size_m=20000
flyby_refinement_radius_m=150000
flyby_refinement_time_min_s=-120
flyby_refinement_time_max_s=100
random_seed=184467
particle_splitting_merging=off
#end application
```

`simulated_species` is the runtime-active subset of the compiled union.
Schema 2 permits only primary neutral H2O and/or CO2 in both
`simulated_species` and `injected_species`; the checked-in stage selects both.
`chemistry_species` lists the neutral parents on which the selected reactions
operate. `reaction_products=discard` is mandatory: a reaction event removes
the parent and creates no daughter particle, including for tabulated
photoionization channels. Every selected reaction must be photolytic and its
parent must be a selected primary. Charge exchange and daughter tracking are
typed startup errors in this schema. Use the literal `none` for both
`chemistry_species` and `selected_reactions` in a transport-only run; listing
chemistry species with no reactions is rejected.

The global step is computed after AMR allocation:

```text
dt = min_species(courant.SPECIES * minimum_allocated_cell_size
                 / characteristic_speed.SPECIES)
```

Every simulated species requires its own finite positive speed and Courant
number in `(0,1]`. The limiting species and cell size are printed in the run
receipt. For primary species, the fixed statistical weight is
`Q * dt / primary_macroparticles_per_step`.

Run-control keys are exact: source totals and primary counts for a species not
in `injected_species`, or speed/Courant keys for a species not in
`simulated_species`, are rejected instead of being silently ignored.

The injector creates exactly that requested integer number of primary model
particles globally on each iteration. All MPI ranks replay the same seeded
event sequence, native AMR ownership selects the creating rank, and each
particle receives an independent uniform birth time within the global step.
The requested count is not a Poisson mean.

For a reacting parent, the runtime sums the selected Burger photo frequencies
and applies the exact constant-rate step probability
`1-exp(-dt*sum(nu))`. On an event the parent is deleted. No daughter species,
statistical weight, momentum partition, or plasma feedback is created or
inferred.

`sample_length_steps` controls the output/sampling cadence;
`total_iterations` is the positive number of native timesteps to execute.
They are intentionally separate. The checked-in run uses ten complete
sampling windows.

`volume_field_output` defaults to `true`. Setting it to `false` suppresses
the generic full-volume H2O/CO2 AMR dumps while retaining completed density
sampling, surface diagnostics, conservation ledgers, flyby/observation
Tecplot products, and sampling-window resets. This is the checked-in mode for
the multi-observation run: ten generic active-species dumps would otherwise
duplicate roughly 14 GiB of fields that the observation operators do not
consume. It is an output-storage control, not a physics or sampling switch.

`reproducible_parallel_run` (optional, `true` or `false`, default `false`;
schema 2 and 3) makes a multi-rank run with a fixed `random_seed` and MPI rank
count reproducible regardless of machine load. Without it two such runs can
differ (2026-11-28, 8 ranks, 5000-step UVIS 2010: 1,132 vs 1,135 resident
macroparticles; columns up to a few percent apart), because two inputs to the
particle histories depend on timing:

- the AMPS core rebalances the domain (emergency rebalancing triggered by the
  measured MPI latency, and rebalancing at every sampling break) by measured
  wall-clock execution time, which changes block ownership and each rank's
  random-number sequence; and
- particles received in the exchange enter the cell lists in MPI
  message-arrival order, and each particle draws its random numbers in list
  order.

`true` sets the core switches `PIC::Parallel::DynamicLoadRebalancingEnabled=false`
(the decomposition made at mesh initialization is kept) and
`PIC::Parallel::DeterministicParticleExchangeOrder=true` (every cell's
particle list is sorted by particle signature after each exchange, the
mechanism the nightly test mode uses). Physics, sampling and output are
unchanged; the run time can grow on an uneven plume load, and which particle
receives which random number differs from a run without it. The end-to-end
regression decks (`test/end_to_end/`) set it. The receipt reports
`reproducible_parallel_run=on|off`.

## Flyby density observation operator

`flyby_density_output=true` requires both `flyby_id` and
`flyby_trajectory_file`; setting it to `false` forbids both dependent keys.
The trajectory path is resolved relative to the run deck. Its CSV must use the
frozen SPICE trajectory schema documented in `standalone/README.md`, every row
must carry the selected flyby ID, all SI values must be finite, altitude must
be nonnegative, and times must increase strictly. Malformed or inaccessible
trajectories fail before MPI initialization.

At every ordinary AMPS output, the adapter evaluates the completed AMPS
number-density sampling buffer at each spacecraft position. The current build
uses the native cell-centered constant stencil, so refinement is required to
reduce containing-cell plateaus; the code does not claim linear interpolation.
The result is a sampling-window mean in m^-3, not an instantaneous particle
count. MPI ownership is reduced to rank zero and all output zones are appended to
`PT/plots/flyby_density.FLYBY.dat`. Each Tecplot row includes time from the
declared flyby reference, body-fixed spacecraft position, SPICE ellipsoid
altitude, `in_domain`, the containing leaf's maximum Cartesian edge
`cell_max_edge_m`, the active primary densities, and their sum.

`flyby_mesh_refinement=true` requires all four refinement controls shown in
the example and requires flyby output. The selected trajectory interval is
treated as a body-fixed polyline. AMR requests the target cell width throughout
the union of radius-`flyby_refinement_radius_m` capsules around its segments;
outside that tube the legacy Enceladus resolution law is unchanged. The target
is a maximum requested edge, not a promise of an exact edge: AMR is dyadic and
may select the next smaller realizable width. The checked-in 20 km request
produces 15.75625 km cells along the in-domain E3 path, compared with the prior
126.05 km cells. Setting refinement to `false` forbids all four dependent keys.

The refinement interval and tube radius are numerical-domain choices, not
instrument smoothing parameters. A resolved physical claim still requires a
second mesh level, particle-count/seed convergence, and an observational
response model. Reducing cell volume without increasing particle statistics
can increase Monte Carlo variance; the output cell-width column makes that
tradeoff explicit in retained evidence.

The complete SPICE trajectory is retained even when it extends beyond the AMR
box. Such rows have `in_domain=0` and zero placeholder density and must be
excluded from physical analysis. The supplied plotting utility does this
automatically:

```text
python3 test/plot_flyby_density.py \
  PT/plots/flyby_density.E3.dat \
  --output PT/plots/flyby_density.E3.png
```

The plotter can also overlay the released E3 INMS observables. The INMS tables
use time from closest approach, while the Tecplot product uses time from the
declared reference UTC. Therefore the conversion offset is explicit and
required; for the frozen E3 geometry, the independently resolved closest
approach is at -0.5 s on the reference-UTC axis:

```text
python3 test/plot_flyby_density.py \
  test_output/srcEnceladus/e3_10000_run/PT/plots/flyby_density.E3.dat \
  --output test_output/srcEnceladus/e3_10000_run/PT/plots/flyby_density.E3.inms_comparison.png \
  --inms-h2o /data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1/e3_h2o_sticking_corrected.csv \
  --inms-co2 /data/vtenishe/cassini_validation_data/data/processed/co2_fit_v1/e3_co2_density.csv \
  --closest-approach-offset-s -0.5
```

The H2O input is the released, sticking-corrected, INMS-derived outbound
profile. Its displayed bars are the independent digitization uncertainty; its
23% common scale uncertainty is stated on the plot but is not folded into each
bar. The E3 mass-44 data are deliberately labeled `apparent CO2`: at this
encounter speed, impact fragmentation can enhance mass-44 counts. Their bars
show Poisson counting uncertainty; the 15% common calibration uncertainty and
the published factor-2-to-5 fragmentation range are stated but not plotted as
independent point errors. A zero-count mass-44 sample is placed at the released
one-count sensitivity with an open downward triangle, rather than being
converted to a positive observation on the logarithmic axis.

The reader fails on non-E3 rows, missing/non-finite values, non-increasing
times, a mass channel other than 44, a CO2 interpretation other than
`apparent_CO2`, or inconsistent zero counts. Observations and model values are
displayed only where the final Tecplot zone is inside the current AMR domain;
this is a direct density overlay, not an INMS response-function convolution or
a numerical-convergence qualification.

## Repeatable INMS and UVIS comparison requests

The `flyby_*` keys above retain the dense frozen E3 trajectory and its AMR
refinement contract. General comparison output is selected independently with
repeatable named `observation` child blocks. The executable contains no list
of allowed Cassini event names: adding a future event means preparing a CSV in
one of the documented canonical forms and adding another uniquely named block.
Comment every line from `#begin observation` through `#end observation` to
omit that data set from a run. Unknown keys, duplicate block names, missing
columns, non-finite values, non-increasing times, or a point on/inside the
252.1 km reference sphere fail before MPI.

One root makes a relocated observation archive a one-line input change.
Ordinary relative `data_file` and `observer_geometry_file` values resolve
against it; absolute paths are also accepted. A `data_file=input:PATH` value
instead resolves `PATH` against the runtime deck directory. The latter keeps
repository-released products portable without coupling them to the external
validation archive:

```text
observation_data_root=/data/vtenishe/cassini_validation_data/data/processed

#begin observation: inms_e3_h2o
instrument=INMS
event_id=E3
data_file=plume_fit_v1/e3_h2o_sticking_corrected.csv
time_column=time_ca_s
x_column=x_m
y_column=y_m
z_column=z_m
comparison_species=H2O
value_column=density_m3
uncertainty_column=sigma_m3
#end observation

#begin observation: uvis_2010_e10_solar
instrument=UVIS
event_id=uvis_2010_e10_solar
data_file=plume_fit_v1/uvis_2010_e10_solar_optical_depth.csv
ray_geometry=finite_observer_to_sun
observer_geometry_file=plume_fit_v1/uvis_2010_e10_solar_observer_range.csv
published_h2o_column_m-2=9.0e19
published_h2o_column_sigma_m-2=2.3e19
#end observation
```

### INMS contract

`instrument=INMS` samples the completed local AMPS number-density buffer at
the body-fixed spacecraft positions named by the four column-mapping keys.
`comparison_species` is `H2O`, `CO2`, or `none`. A calibrated comparison
requires value and positive-uncertainty columns. With `none`, omitting all
observable keys creates an explicitly model-only profile. Alternatively, the
following three keys retain a native instrument observable as a labelled
context panel without comparing it with model density:

```text
comparison_species=none
context_observable=INMS_OSNB_H2O_counts_per_integration
context_value_column=counts_per_integration
context_uncertainty_column=poisson_sigma_counts
```

Context values must be finite and nonnegative and uncertainties finite and
positive. A request cannot declare both a density comparison and context. The
Tecplot schema uses the unit-neutral names `observed_context_value` and
`observed_context_sigma`; the plot/JSON label it
`display_only_no_density_or_count_response` so a count-to-density fit cannot
be inferred accidentally.

Each request produces `PT/plots/observation.REQUEST.dat`. Every zone records
time, position, spherical altitude, AMR ownership, cell edge, all active model
densities, and—only when declared—the observed density and uncertainty. Rows
outside the AMR box have `in_domain=0` and zero model placeholders.

The checked-in deck requests all eight currently represented INMS plume
flybys. Qualified H2O density exists for E3/E5; qualified CO2 or
apparent-CO2 density exists for E3/E5/E7/E14/E17/E18. E7's 43-row mass-44
density vector is exported from SwRI `Flybys2.mat` and matched count-for-count
to official PDS geometry by `standalone/export_forward_observations.py`.
E13 remains model-only because its mass-13/CH4 table supplies geometry but no
H2O/CO2 density. E21 carries the released open-source neutral-beaming H2O and
CO2 counts as native context; they are not scored until the AMPS phase-space
distribution is accumulated through the supplied velocity-, voltage-, and
attitude-dependent SIMION response. Raw counts or CH4 proxies must never be
relabelled as density.

### UVIS contract

UVIS canonical optical-depth CSVs contain tangent position, source time,
time from minimum tangent altitude, raw optical-depth proxy, and its Poisson
uncertainty. Two ray modes are implemented:

- `infinite_chord_tangent_track_normal` reconstructs the stellar line of
  sight as the continuous unit vector parallel to
  `r_tangent x d(r_tangent)/dt`; a nonuniform three-point derivative is used,
  and AMPS integrates both directions from the tangent point.
- `finite_observer_to_sun` reads a row-aligned auxiliary table containing the
  signed observer coordinate and body-fixed ray direction. For a negative
  coordinate the observer is `r_tangent+s_observer*l` and integration proceeds
  toward the Sun; the direction is projected perpendicular to the tangent
  radius exactly as in the qualified standalone geometry. Positive-coordinate
  rows are integrated outward from the tangent point.

The native AMPS column integrator uses a trapezoidal step equal to 0.3 of the
local AMR characteristic cell size and reduces owner contributions across
MPI. Output columns have units m^-2 and cover only the finite AMPS box; each
row records its integrated in-domain path length. The current UVIS operator
integrates H2O and every other active neutral but compares only H2O.

The broadband optical-depth proxy is not treated as an absolute column at
every sample. `test/plot_observation_comparisons.py` fits one nonnegative
column-to-optical-depth coefficient plus one signed baseline for shape, while
the two explicit published column values provide separate absolute peak-scale
tests. The current deck requests the two fit-ready products used by the
standalone model: the 2007 zeta Orionis stellar occultation and the 2010 E10
solar occultation. Published inventories identify additional 2005, 2011,
2016, and 2017 events, but the current local validation release contains no
qualified reduced spectra plus explicit row-aligned ray tables for them.
They are therefore not input requests: adding them first requires a
reproducible spectral/background reduction and geometry release, not an
assumed profile.

The present ray products also retain two documented geometry limitations.
The 2007 file has tangent positions but no row-aligned SPICE star-direction
table, so `infinite_chord_tangent_track_normal` is a tangent-track
reconstruction rather than a SPICE-qualified stellar ray. In the 2010
auxiliary file, 308 of 1,133 rows have a nonnegative observer coordinate but
do not retain the full three-dimensional Cassini position; those outward
branches begin at the tangent point. The other 825 rows use the finite lower
limit on the tangent-line parameter. Both products are diagnostic until those
geometry tables are completed, and their output must not be described as
publication-qualified occultation columns.

To validate ten zones and generate one PNG per request:

```text
python3 test/plot_observation_comparisons.py \
  test_output/srcEnceladus/e3_10000_extended_observations/PT/plots/observation.*.dat \
  --output-dir test_output/srcEnceladus/e3_10000_extended_observations/PT/plots/comparisons \
  --expected-zones 10
```

The JSON summary reports in-domain counts, density-qualified INMS chi-square
and median positive model/observation ratios, display-only native-context
counts, UVIS weighted shape fit, and model peak relative to the published
column. These are diagnostics, not a refit or a convergence qualification.

## Plasma-background authority

The neutral-only mode used by the current E3 run is:

```text
plasma_background_mode=none
```

It forbids cache and ICES keys, launches no external process, reads no plasma
file, and writes no `ices.data.dat`. After AMR allocation the adapter sets the
entire allocated coupled DATAFILE center-node background region to zero. This
explicit initialization is required because the compiled Boris implementation
interpolates E and B even for a neutral particle whose charge makes the
Lorentz acceleration zero. Enceladus gravity and surface/outer boundaries
remain active in this mode.

The later coupled modes remain available but were not used in the current run.
They are one-way: particles read plasma density, bulk velocity, pressure,
temperature, electric field, and magnetic field but never modify the plasma.
Select an explicit background authority:

```text
plasma_background_mode=cache_only
plasma_cache_directory=plasma_cache
```

`plasma_cache_directory` is resolved relative to the runtime deck and must
already be an accessible directory. The generic AMPS data-file reader names
the file from the exact AMR mesh signature plus the stable identity
`ENCELADUS-BATSRUS`; therefore a cache produced for another mesh cannot be
silently reused. In `cache_only` mode, a missing mesh-matched cache is fatal
and no external process is launched.

The existing ICES extraction path is available only when explicitly enabled:

```text
plasma_background_mode=cache_or_ices
plasma_cache_directory=plasma_cache
ices_installation_directory=/authoritative/path/to/ICES/Models
ices_case=AUTHORITATIVE_CASE_NAME
```

`cache_or_ices` first loads the same mesh-matched cache when present. Only
when it is absent does it invoke the compiled ICES coordinate/extraction
workflow, read the resulting SWMF data, and save the cache. Both ICES values
are mandatory in this mode and forbidden in `cache_only`; the legacy
developer-specific ICES path and case are never used as implicit defaults.
Because the inherited ICES routine composes an unquoted shell command, the
installation path is restricted to letters, digits, `/`, `_`, `-`, and `.`,
and the case must be one non-`.`/`..` component using the same set without
`/`. Unsupported whitespace or shell metacharacters are startup errors.
The build must also contain the established ICES coordinate-list capability,
otherwise the fallback terminates with a typed error.

## Burger 2007 reaction identifiers

Photo rates are the quiet/active endpoints at 9.5 AU from Burger et al.,
“Understanding the escape of water from Enceladus,” JGR 112 A06219 (2007),
doi:10.1029/2006JA012086, Table 2:

- `burger2007.photo.h2o_to_oh_h`
- `burger2007.photo.h2o_to_h2_o`
- `burger2007.photo.h2o_to_h_h_o`
- `burger2007.photo.h2o_to_h2op_e`
- `burger2007.photo.h2o_to_h_ohp_e`
- `burger2007.photo.h2o_to_h2_op_e`
- `burger2007.photo.h2o_to_oh_hp_e`

The library also retains Burger Table 1 ion-neutral cross sections under
stable IDs for later model stages and independent constant checks:

- `burger2007.cx.hplus_h2o_to_h_h2op`
- `burger2007.cx.oplus_h2o_to_o_h2op`
- `burger2007.cx.h2op_h2o_resonant`
- `burger2007.cx.ohplus_h2o_to_oh_h2op`

They are not selectable in application schema 2. Any `burger2007.cx.*`
identifier, charged simulated species, ion fraction, or charge-exchange weight
reference is rejected. The H2O+ + H2O -> OH + H3O+ channel remains absent from
the table entirely, and H3O+ remains outside the model.

The loss processor is the configured AMPS user-particle pass, executed once
before the ordinary mover. Existing particles receive the exact competing
photo-channel probability `1-exp(-sum(nu)*dt)` for a full global step.
Fresh primaries are injected and moved through a uniformly sampled remaining
fraction of that step into AMPS temporary lists, so their first chemistry pass
occurs at the start of the following step. Their reaction exposure during that
initial partial step is therefore omitted by the current operator splitting.
This known `O(nu*dt)` birth-time error must be measured in timestep convergence
tests; no claim of exact continuous-time chemistry is made.

Both flyby receipts contain CO2 as a neutral transport primary, while
`chemistry_species=H2O`. No Burger-2007 CO2 photo channel exists, so a CO2
reaction identifier is not offered and startup rejects invented CO2
chemistry. The E3/E5 receipt values are the standalone report's apparent-CO2
mass-44 rates, with a documented factor-2--5 fragmentation bracket. Select a
physical value within that bracket only when an external scientific choice
justifies the enhancement factor; the startup fixture makes no such choice.

## Sampling and particle population control

`sample_length_steps` sets the fixed AMPS sampling window. Population control
is selected with `particle_splitting_merging=off|velocity_shift|scatter`.
Enabled modes additionally require
`min_particles_per_cell < max_particles_per_cell`. `velocity_shift` also
requires positive dimensionless
`split_position_shift_cell_diagonal_fraction` and
`split_velocity_shift_thermal_speed_fraction`. The generic controller scales
the first by the local cell diagonal and the second by the measured local
thermal speed; they are not metres and m/s. `nonuniform_split_position_shift`
is `true` or `false` and selects the generic controller's nonuniform radial
position-shift proposal. These values map directly to the existing generic
PIC split/merge controller; no Enceladus-specific algorithm replaces it.
`particle_splitting_merging=off` forbids every split/merge parameter key, and
`scatter` forbids the velocity-shift fractions and nonuniform-position flag,
so inactive controls cannot remain silently in a run receipt.

## Transport mover

The compile deck selects generic PIC pointer-mode mover dispatch. After the
runtime deck and mesh are validated, the adapter installs the existing Boris
mover. The schema-2 active populations are neutral, so their trajectory uses
Enceladus point-mass gravity and generic internal-sphere handling; no Lorentz
force or charged daughter is present. In `none` mode the coupled-background
storage is zeroed. In either cache mode the established one-way background is
loaded, but it receives no particle feedback.

At normal termination, the adapter reports per-species injection,
photolytic-parent loss, surface absorption, true outer-domain escape, and
resident counts. It traverses native cell particle lists independently and
prints `ledger_residual = injected - photolytic - surface - escape - resident`.
A qualified discrete run requires zero residual for every injected species.
