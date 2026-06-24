# AMPS Earth Energetic-Particle Model

This directory contains the Earth/geospace energetic-particle tools used by AMPS in standalone runs and in SWMF-coupled PT runs. The code supports three main execution paths:

1. `-mode gridless` — backward trajectory tracing with direct magnetic-field evaluation.
2. `-mode 3d` — backward trajectory tracing with magnetic fields stored/interpolated on the AMPS AMR mesh.
3. `-mode 3d_forward` — forward Monte Carlo transport with boundary injection and volumetric AMR sampling.

In SWMF-coupled builds, the same mesh-backed Mode3D backward-product logic can be called from the coupled PT time step using the magnetic field imported from SWMF.

The most important distinction is this:

- `gridless` and mesh-backed `3d` cutoff/density/flux are **backward-access / transmissivity** calculations.
- `3d_forward` is a **forward particle transport / residence-time sampling** calculation.

They answer related but not identical questions.

---

## 1. Source layout

Important source files:

```text
srcEarth/main.cpp
    Standalone CLI dispatch for -mode gridless, -mode 3d, and -mode 3d_forward.

srcEarth/util/amps_param_parser.h
srcEarth/util/amps_param_parser.cpp
    Strict AMPS_PARAM.in parser used by the gridless, Mode3D, and 3d_forward paths.

srcEarth/util/cutoff_cli.h
srcEarth/util/cutoff_cli.cpp
    Command-line option parser and help text.

srcEarth/gridless/CutoffRigidityGridless.cpp
srcEarth/gridless/DensityGridless.cpp
srcEarth/gridless/GridlessParticleMovers.cpp
srcEarth/gridless/AnisotropicSpectrum.cpp
    Direct-field, gridless backward tracing, cutoff, density, spectrum, flux,
    and boundary anisotropy support.

srcEarth/3d/Mode3D.cpp
srcEarth/3d/CutoffRigidityMode3D.cpp
srcEarth/3d/DensityMode3D.cpp
srcEarth/3d/ElectricField.cpp
srcEarth/3d/GlobalMagneticField.cpp
    Mesh-backed 3-D backward products: cutoff, directional cutoff maps,
    density/flux, standalone time snapshots, field initialization, and global
    replicated B-field materialization.

srcEarth/3d_forward/Mode3DForward.cpp
srcEarth/3d_forward/ForwardParticleMovers.cpp
srcEarth/3d_forward/Density3D.cpp
srcEarth/3d_forward/SphereFlux3D.cpp
    Standalone forward 3-D particle injection, propagation, cell density sampling,
    and inner-sphere impact flux sampling.

srcEarth/3d_forward_swmf/Mode3DForwardSWMF.cpp
srcEarth/3d_forward_swmf/Mode3DForwardSWMF.h
    SWMF-coupled bridge. In coupled builds it supports both historical forward
    injection and the newer Mode3D backward products driven by live SWMF fields.

parallel_affinity.h
parallel_affinity.cpp
    Main AMPS/PIC affinity helpers under namespace PIC::Parallel. These files
    are not srcEarth-local files; they are intended to be compiled by the main
    AMPS build system and made available through the AMPS include path. The
    Mode3D direct-thread cutoff/density backends call these helpers to widen
    the MPI-rank CPU affinity mask before std::thread workers are created.
```

---

## 2. Execution modes

### 2.1 `-mode gridless`

Run command:

```bash
mpirun -np 8 ./amps -mode gridless -i AMPS_PARAM.in
```

The gridless mode does not build or sample a 3-D magnetic-field mesh. Each trajectory step evaluates the selected background magnetic field directly. This mode is useful for fast access studies, analytic/dipole validation, cutoff calculations, and transmissivity-based density/flux calculations.

Supported main products:

```text
CALC_TARGET CUTOFF_RIGIDITY
CALC_TARGET DENSITY_SPECTRUM
```

The gridless density/flux calculation is a backward-tracing calculation. It does not inject particles forward in time. For each observation point and energy, the code samples arrival directions, backtraces each direction, computes the allowed fraction, multiplies the boundary spectrum by that transmissivity, and integrates the result into density and flux.

Typical output files:

```text
cutoff_gridless_points.dat
cutoff_gridless_shells.dat
cutoff_gridless_dir_map_point_####.dat

gridless_points_density.dat
gridless_points_spectrum.dat
gridless_points_flux.dat
gridless_shell_<ALT>km_density_channels.dat
```

### 2.2 `-mode 3d`

Run command:

```bash
mpirun -np 8 ./amps -mode 3d -i AMPS_PARAM.in
```

Mode3D uses the AMPS AMR mesh as the field backend. The standalone code initializes the requested background magnetic field on mesh cell centers, then materializes a replicated read-only B-field snapshot on every MPI rank so any rank can backtrace through the whole domain.

This mode now supports three product selections from the same prepared mesh-field snapshot:

```text
CALC_TARGET CUTOFF_RIGIDITY
CALC_TARGET DENSITY_SPECTRUM
CALC_TARGET CUTOFF_RIGIDITY+DENSITY_SPECTRUM
```

The aliases below are also accepted by the Mode3D product selector:

```text
CALC_TARGET ALL
CALC_TARGET BOTH
```

Mode3D products:

- cutoff rigidity;
- optional directional cutoff sky maps;
- gridless-style density, local spectrum, and integral flux, but using the meshed magnetic field;
- optional energy-channel integral fluxes;
- standalone single-snapshot runs;
- standalone time-series runs over multiple Tsyganenko-driver snapshots.

The expensive Mode3D backward products can use one of three intra-rank
shared-memory backends, selected by the `#NUMERICAL` keywords
`DENSITY_PARALLEL` and `DENSITY_THREADS` or by the standalone CLI options
`-density-parallel` and `-density-threads`. The names are historical: the
same settings now apply to both `CUTOFF_RIGIDITY` and `DENSITY_SPECTRUM`.

Typical output files for a single snapshot:

```text
cutoff_3d_points.dat
cutoff_3d_shells.dat
cutoff_3d_dir_map_loc_000000.dat
cutoff_3d_dir_map_loc_000001.dat

mode3d_points_density.dat
mode3d_points_spectrum.dat
mode3d_points_flux.dat
mode3d_shell_<ALT>km_density_flux.dat
```

For time-series standalone runs, a snapshot suffix is inserted before `.dat`, for example:

```text
cutoff_3d_points_snapshot_000000_2024_05_10T00_00_00_000.dat
mode3d_points_density_snapshot_000000_2024_05_10T00_00_00_000.dat
mode3d_points_spectrum_snapshot_000000_2024_05_10T00_00_00_000.dat
mode3d_points_flux_snapshot_000000_2024_05_10T00_00_00_000.dat
cutoff_3d_dir_map_loc_000000_snapshot_000000_2024_05_10T00_00_00_000.dat
```

### 2.3 `-mode 3d_forward`

Run command:

```bash
mpirun -np 8 ./amps -mode 3d_forward -i AMPS_PARAM.in
```

The forward 3-D mode is a particle-in-cell style Monte Carlo calculation. Particles are injected at the outer rectangular boundary, propagated forward in time through the configured B/E fields, and sampled in AMR cells and on the inner absorbing sphere.

This path is not the same as the gridless/Mode3D transmissivity density/flux calculation.

Main products:

- AMR volumetric sampled density in energy bins;
- inner-sphere incident flux maps;
- optional initialized field mesh dump;
- optional particle trajectory records when AMPS is compiled with particle tracking.

Typical output files include AMPS mesh data files plus:

```text
sphere_flux3d.out=####.dat
sphere_flux3d_total.out=####.dat
amps_3dforward_initialized.data.dat     # only with -mode3d-output-initialized
```

The cell-density variables are written into the standard AMPS mesh output by `Density3D.cpp` as variables like:

```text
n_dens[Elo-Ehi_MeV]_m3
n_dens_total_m3
```

---

## 3. SWMF-coupled behavior

In a SWMF-coupled build, the PT component does not normally enter through `./amps -mode 3d`. Instead, the SWMF/PT execution path calls the hooks in:

```text
srcEarth/3d_forward_swmf/Mode3DForwardSWMF.cpp
```

The coupled bridge has two behaviors.

### 3.1 Historical coupled forward mode

If `CALC_TARGET` is not explicitly set to a backward product, the coupled path remains the historical `3d_forward` style mode: particles are injected and propagated forward using fields received from SWMF.

This behavior is preserved intentionally for backward compatibility.

### 3.2 Coupled Mode3D backward-product mode

If `#CALCULATION_MODE / CALC_TARGET` explicitly requests cutoff and/or density/flux, the coupled bridge switches to mesh-backed backward products:

```text
CALC_TARGET CUTOFF_RIGIDITY
CALC_TARGET DENSITY_SPECTRUM
CALC_TARGET CUTOFF_RIGIDITY+DENSITY_SPECTRUM
CALC_TARGET ALL
```

The coupled backward-product path uses the same Mode3D shared-memory backend
as standalone Mode3D. In particular, `DENSITY_PARALLEL THREADS` enables
direct `std::thread` workers for cutoff and density calculations inside each
MPI rank, and `DENSITY_THREADS` sets the worker count.

At every accepted coupled snapshot, the bridge:

1. uses the current SWMF/PT magnetic field in the AMPS coupler buffers;
2. materializes the global cell-centered B field on every MPI rank;
3. runs the requested Mode3D backward products;
4. writes files with a SWMF simulation-time suffix.

Example coupled output names:

```text
cutoff_3d_points.swmf_n000000_t000000000.000s.dat
mode3d_points_density.swmf_n000000_t000000000.000s.dat
mode3d_points_spectrum.swmf_n000000_t000000000.000s.dat
mode3d_points_flux.swmf_n000000_t000000000.000s.dat
amps_coupled_data.swmf_n000000_t000000000.000s.dat
```

The coupled cadence is controlled by:

```text
#TEMPORAL
FIELD_UPDATE_DT  <minutes, floating point allowed>
```

In the coupled case, `FIELD_UPDATE_DT` means: run the expensive backward products approximately every requested number of minutes of SWMF/PT simulation time, using `PIC::SimulationTime::TimeCounter` for the file stamp. Fractional minutes are accepted, for example `FIELD_UPDATE_DT 0.5` requests a 30-second cadence.

---

## 4. Algorithms

### 4.1 Backward trajectory classifier

Both gridless and Mode3D backward products use the same physical classification:

```text
Launch time-reversed particle from observation point x0.

If the trajectory exits the outer domain box first:
    direction is ALLOWED.

If the trajectory hits the inner absorbing sphere first:
    direction is FORBIDDEN.

If the trajectory exceeds max time, max steps, or max trace distance:
    direction is treated as FORBIDDEN conservatively.
```

The outer domain is defined by `#DOMAIN_BOUNDARY`. The inner loss sphere is defined by `R_INNER`.

### 4.2 Cutoff rigidity

For each observation point and arrival direction, the solver searches for the transition between forbidden and allowed access over the rigidity/energy range specified in `#CUTOFF_RIGIDITY`.

For `CUTOFF_SAMPLING VERTICAL`, one local vertical direction is tested.

For `CUTOFF_SAMPLING ISOTROPIC`, the code samples multiple directions and reports an effective/minimum cutoff for the location.

Optional directional sky maps are controlled independently by:

```text
DIRECTIONAL_MAP T
DIRMAP_LON_RES  <deg>
DIRMAP_LAT_RES  <deg>
```

Directional maps write one Tecplot file per observation location. Direction labels follow the gridless convention: SM lon/lat when SPICE is available for SM-to-GSM rotation; GSM fallback when SPICE is unavailable.

### 4.3 Density, local spectrum, and integral flux

The density/flux calculation in `gridless` and mesh-backed `3d` is transmissivity-based.

For each observation point and energy:

```text
T(E; x0) = allowed_direction_weight / total_direction_count
J_local(E; x0) = T(E; x0) * J_boundary(E)
```

For isotropic boundary mode, each allowed direction has unit weight.

For anisotropic boundary mode, each allowed trajectory is weighted using the boundary pitch-angle and spatial factors evaluated at the outer-boundary exit point.

The omnidirectional number density and flux are then integrated as:

```text
n(x0) = 4*pi * integral[ J_local(E; x0) / v(E) dE ]
F(x0) = 4*pi * integral[ J_local(E; x0) dE ]
```

Additional user-defined energy channels can be integrated using `#ENERGY_CHANNELS`.

### 4.4 3-D forward Monte Carlo transport

The forward mode uses a different algorithm:

```text
For each forward iteration:
    inject particles through the outer box boundary;
    sample energy from the configured injection distribution;
    sample inward velocity direction;
    assign statistical weight;
    move particles forward through B/E fields;
    accumulate cell residence density and inner-sphere impact flux.
```

This is useful for explicit forward propagation and residence-time diagnostics. It is not a direct replacement for the transmissivity-based density/flux calculation.

---

## 5. Magnetic and electric fields

### 5.1 Gridless magnetic fields

The gridless branch evaluates magnetic fields directly. The source supports:

```text
DIPOLE
T96
T05
```

### 5.2 Standalone Mode3D magnetic fields

Standalone Mode3D initializes the magnetic field on the AMPS AMR mesh. Supported field models include:

```text
DIPOLE
T96
T05
TA16
```

The field is first written into owner-rank DATAFILE cell buffers. Then `GlobalMagneticField::MaterializeCellCenteredMagneticFieldForCutoff()` assigns dense leaf IDs, allocates missing leaf blocks on all MPI ranks, gathers owner-cell B values, and fills all replicated blocks/ghost cells. After that step, the normal AMPS interpolation stencil can be used from any MPI rank.

### 5.3 SWMF-coupled Mode3D magnetic fields

In coupled mode, the field comes from the SWMF coupler data buffer, not from standalone Tsyganenko/DIPOLE initialization.

The same global materialization helper is used, but the source offset is:

```cpp
PIC::CPLR::SWMF::MagneticFieldOffset
```

### 5.4 Electric field options

The standalone 3-D field initializer can also populate cell-centered electric fields. Recognized model names include:

```text
NONE
COROTATION
VOLLAND_STERN
COROTATION_VOLLAND_STERN
```

The backward cutoff/density access problem is primarily a magnetic-field access calculation. The forward transport path can use the configured electric field in particle motion.

---

## 6. MPI model

### 6.1 Gridless

Gridless mode does not need a field mesh. Work is distributed over observation locations, energies, and/or directions. The magnetic field is evaluated directly along each trajectory.

### 6.2 Standalone Mode3D

Standalone Mode3D no longer uses independent private MPI domains for cutoff calculations. It uses the normal distributed AMPS mesh initialization and then builds a replicated read-only magnetic-field snapshot for tracing.

The intended sequence is:

```text
PIC::InitMPI()
amps_init_mesh()
amps_init()
InitMeshFields()
GlobalMagneticField::MaterializeCellCenteredMagneticFieldForCutoff()
RunCutoffRigidity() and/or RunDensityAndFlux()
```

This gives all ranks access to the global B field while preserving the normal AMPS MPI domain decomposition during initialization.

### 6.3 SWMF coupled

In coupled mode, the global B-field materialization is repeated for every accepted SWMF/PT snapshot before the backward products are computed.

### 6.4 Intra-rank shared-memory backends for Mode3D backward products

Mode3D cutoff and density/flux products can use three shared-memory backends
inside each MPI rank:

```text
DENSITY_PARALLEL SERIAL
DENSITY_PARALLEL OPENMP
DENSITY_PARALLEL THREADS
```

Despite the `DENSITY_` prefix, these settings control both the cutoff and
density/flux backward products. The prefix is kept for compatibility with
older input files and command-line options.

`SERIAL` disables intra-rank parallelism. MPI decomposition over locations is
still used.

`OPENMP` preserves the OpenMP implementation. For density it parallelizes the
energy/direction work where enabled. For cutoff it uses OpenMP over the local
observation-location loop. This backend is useful when the OpenMP runtime and
MPI launcher are known to place threads correctly.

`THREADS` uses direct `std::thread` workers over local observation locations.
It intentionally suppresses nested OpenMP inside those workers to avoid
oversubscription. This backend is often preferable when OpenMP interpolation
stencil state or MPI/OpenMP placement is problematic.

The worker count is controlled by:

```text
DENSITY_THREADS <N>
```

or by the standalone CLI option:

```bash
-density-threads <N>
```

`N > 1` requests that many workers per MPI rank. `N = 0` asks the code to
choose automatically from environment/runtime information. For production
runs, the number of MPI ranks per node times `DENSITY_THREADS` should not
exceed the number of CPU cores allocated on that node.

For the direct-thread backend, the code calls:

```cpp
PIC::Parallel::SetWideAffinityForScheduler();
```

before the cutoff or density worker pool is created. This reproduces the
manual operation:

```bash
taskset -apc <CPUSET> <rank_pid>
```

for each MPI rank. It does not pin individual worker threads to fixed cores;
it widens the allowed CPU mask of the rank so the Linux scheduler can move
those threads among the allowed CPUs. By default the CPU set is read from:

```text
/sys/devices/system/cpu/online
```

The automatic CPU set can be overridden with either environment variable:

```bash
export AMPS_MODE3D_DENSITY_CPUSET=0-39
export PIC_PARALLEL_CPUSET=0-39
```

If the batch scheduler or MPI runtime restricts a rank to one CPU with a
cgroup/cpuset, affinity widening cannot escape that restriction. The actual
mask can be checked at runtime from the job log through the diagnostic printed
by `PIC::Parallel::PrintCurrentAffinity()`, or manually with:

```bash
grep Cpus_allowed_list /proc/<pid>/status
ps -L -p <pid> -o pid,tid,psr,pcpu,state,comm
```

---

## 7. Time dependence and snapshots

### 7.1 Single-snapshot standalone run

By default, standalone Mode3D uses one magnetic-field snapshot at:

```text
#BACKGROUND_FIELD
EPOCH <UTC>
```

If a driver file is present but `TEMPORAL_MODE` remains `SNAPSHOT`, the table is sampled once at the background-field epoch.

### 7.2 Standalone time-series run

Standalone Mode3D can loop through multiple background-field snapshots using a Tsyganenko driver table. The table can be specified in either of two ways.

Option A:

```text
#TEMPORAL
TS_INPUT_MODE FILE
TS_INPUT_FILE  ts_driving_file.txt
```

Option B:

```text
#BACKGROUND_FIELD
DRIVER_FILE    ts_driving_file.txt
```

Then request a time series:

```text
#TEMPORAL
TEMPORAL_MODE   TIME_SERIES
EVENT_START     2024-05-10T00:00:00
EVENT_END       2024-05-10T12:00:00
FIELD_UPDATE_DT 15
```

For each snapshot, the code interpolates the driver table to the requested epoch, updates the background-field parameters, rebuilds/materializes the mesh field, and runs the requested products.

SPICE must be available for time-series operation because UTC strings and driver-table timestamps are converted to ephemeris time.

### 7.3 SWMF-coupled time dependence

In coupled mode, the magnetic-field snapshots come from SWMF itself. `FIELD_UPDATE_DT` is used as the calculation cadence, not as a request to load a Tsyganenko file.

---

## 8. AMPS_PARAM.in sections

The parser is strict: unknown sections or unknown keywords are fatal. This is intentional so typos do not silently produce wrong physics.

### 8.1 `#RUN_INFO`

```text
#RUN_INFO
RUN_ID  my_run_name
```

`RUN_ID` is a user label. It does not by itself alter the calculation.

### 8.2 `#CALCULATION_MODE`

```text
#CALCULATION_MODE
CALC_TARGET        CUTOFF_RIGIDITY
FIELD_EVAL_METHOD  GRIDLESS
```

Supported `CALC_TARGET` values:

```text
CUTOFF_RIGIDITY
DENSITY_SPECTRUM
CUTOFF_RIGIDITY+DENSITY_SPECTRUM
ALL
BOTH
DENSITY_3D      # forward mode validation path
```

Recommended `FIELD_EVAL_METHOD` values:

```text
GRIDLESS   # direct field evaluation
GRID_3D    # mesh-backed standalone 3D
SWMF       # set internally in SWMF-coupled mode
```

Notes:

- `-mode gridless` with density requires `FIELD_EVAL_METHOD GRIDLESS`.
- `-mode 3d` uses the mesh-backed Mode3D driver even if the field method label is not the only selector.
- In SWMF-coupled backward-product mode, `CALC_TARGET` must be explicit; otherwise the bridge preserves historical forward-injection behavior.

### 8.3 `#PARTICLE_SPECIES`

```text
#PARTICLE_SPECIES
SPECIES_NAME       PROTON
SPECIES_CHARGE     1
SPECIES_MASS_AMU   1.0073
```

Aliases are also accepted:

```text
SPECIES
CHARGE
MASS_AMU
```

### 8.4 `#BACKGROUND_FIELD`

```text
#BACKGROUND_FIELD
FIELD_MODEL    T05
EPOCH          2024-05-10T00:00:00
PDYN           2.0
DST            -20.0
IMF_BX         0.0
IMF_BY         0.0
IMF_BZ         -5.0
SW_VX          -400.0
SW_N           5.0
T05_W1         0.0
T05_W2         0.0
T05_W3         0.0
T05_W4         0.0
T05_W5         0.0
T05_W6         0.0
```

Common model-specific keys:

```text
DIPOLE_MOMENT
DIPOLE_TILT or DIPOLE_TILT_DEG
G1 G2 G3
W1..W6 or T05_W1..T05_W6 or TS05_W1..TS05_W6
BZ1..BZ6 or TA15_BZ1..TA15_BZ6 or TA16_BZ1..TA16_BZ6
XIND or TA15_XIND
TA16_COEFF_FILE
DRIVER_FILE or MAGNETIC_DRIVER_FILE
```

Model aliases like `TS05` are normalized to the canonical names used internally.

### 8.5 `#ELECTRIC_FIELD`

```text
#ELECTRIC_FIELD
EFIELD_MODEL       NONE
COROTATION_SCALE   1.0
VS_POTENTIAL_KV    60.0
VS_GAMMA           2.0
VS_REFERENCE_L     10.0
VS_SCALE           1.0
EFIELD_RMIN        20.0
EFIELD_LMIN        1.0e-3
```

`EFIELD_MODEL` values:

```text
NONE
COROTATION
VOLLAND_STERN
COROTATION_VOLLAND_STERN
```

### 8.6 `#DOMAIN_BOUNDARY`

```text
#DOMAIN_BOUNDARY
DOMAIN_XMIN   -35 Re
DOMAIN_XMAX    35 Re
DOMAIN_YMIN   -35 Re
DOMAIN_YMAX    35 Re
DOMAIN_ZMIN   -35 Re
DOMAIN_ZMAX    35 Re
R_INNER        1 Re
```

Lengths can be provided using inline units where supported by the parser, for example `Re` or `km`. The solver internally uses SI meters.

The outer box defines escape. `R_INNER` defines the absorbing loss sphere.

### 8.7 `#OUTPUT_DOMAIN`

POINTS example:

```text
#OUTPUT_DOMAIN
OUTPUT_MODE  POINTS
COORDS       GSM
POINTS_BEGIN
  POINT 1.1 Re 0.0 0.0
  POINT 2.0 Re 0.0 0.0
POINTS_END
```

SHELLS example:

```text
#OUTPUT_DOMAIN
OUTPUT_MODE      SHELLS
SHELL_ALTS_KM    500 9000
SHELL_RES_DEG    10
```

TRAJECTORY example:

```text
#OUTPUT_DOMAIN
OUTPUT_MODE  TRAJECTORY
TRAJ_FRAME   GSM
TRAJ_FILE    spacecraft_path.txt
FLUX_DT      1.0
```

Supported output modes depend on the solver:

```text
gridless cutoff:        POINTS, TRAJECTORY, SHELLS
gridless density/flux:  POINTS, TRAJECTORY, SHELLS
Mode3D cutoff:          POINTS, TRAJECTORY, SHELLS
Mode3D density/flux:    POINTS, TRAJECTORY, SHELLS
3d_forward:             AMR mesh sampling + inner sphere; not POINTS/SHELLS transmissivity
```

Coordinate note: the trajectory and newer parser paths can convert some trajectory frames to GSM when SPICE is available, but the physical solvers operate in GSM. Verify coordinate conventions before using non-GSM inputs.

### 8.8 `#NUMERICAL`

```text
#NUMERICAL
DT_TRACE             1.0
MAX_STEPS            300000
MAX_TRACE_TIME       7200.0
MAX_TRACE_DISTANCE   0.0
DENSITY_PARALLEL     THREADS
DENSITY_THREADS      8
```

Meanings:

```text
DT_TRACE            initial trajectory time step [s]
MAX_STEPS           maximum integration steps per trajectory
MAX_TRACE_TIME      maximum integration time per trajectory [s]
MAX_TRACE_DISTANCE  cumulative path-length cap [Re]; <=0 disables it
DENSITY_PARALLEL    SERIAL, OPENMP, or THREADS for Mode3D cutoff/density products
DENSITY_THREADS     shared-memory worker count per MPI rank; 0 means automatic
```

`DENSITY_PARALLEL` and `DENSITY_THREADS` apply to both Mode3D cutoff and
density/flux calculations. The `DENSITY_` prefix is historical. In SWMF-coupled
backward-product runs these keywords can be used in the PT `AMPS_PARAM.in` to
select direct-thread or OpenMP execution without standalone CLI options.

Forward-mode controls also accepted here:

```text
FORWARD_N_PARTICLES
FORWARD_N_ITERATIONS
FORWARD_INJECTION_ENERGY
FORWARD_INJECTION_ENERGY_DISTRIBUTION
FORWARD_ENERGY_SAMPLING
FORWARD_INJECTION_SCHEME
```

The parser recognizes forward mover keywords for future use, but the active `3d_forward` mover is currently selected by the command line or falls back to `BORIS`.

### 8.9 `#CUTOFF_RIGIDITY`

```text
#CUTOFF_RIGIDITY
CUTOFF_EMIN            1.0
CUTOFF_EMAX            20000.0
CUTOFF_NENERGY         50
CUTOFF_MAX_PARTICLES   500
CUTOFF_MAX_TRAJ_TIME   60.0
CUTOFF_SAMPLING        ISOTROPIC
DIRECTIONAL_MAP        T
DIRMAP_LON_RES         30
DIRMAP_LAT_RES         30
```

Meanings:

```text
CUTOFF_EMIN / CUTOFF_EMAX     energy/range bounds [MeV/n]
CUTOFF_NENERGY                number of samples used by the cutoff search
CUTOFF_MAX_PARTICLES          direction/sample cap per location
CUTOFF_MAX_TRAJ_TIME          cutoff-specific trajectory time cap [s]
CUTOFF_SAMPLING               VERTICAL or ISOTROPIC
DIRECTIONAL_MAP               write directional sky-map products
DIRMAP_LON_RES/LAT_RES        sky-map angular resolution [deg]
```

### 8.10 `#DENSITY_SPECTRUM`

```text
#DENSITY_SPECTRUM
DS_EMIN            1.0
DS_EMAX            20000.0
DS_NINTERVALS      50
DS_ENERGY_SPACING  LOG
DS_MAX_PARTICLES   500
DS_MAX_TRAJ_TIME   60.0
DS_BOUNDARY_MODE   ISOTROPIC
```

Meanings:

```text
DS_EMIN / DS_EMAX       density/flux energy bounds [MeV/n]
DS_NINTERVALS           number of energy intervals; number of grid points is +1
DS_ENERGY_SPACING       LOG or LINEAR
DS_MAX_PARTICLES        total trajectory cap per observation point; 0 means no cap
DS_MAX_TRAJ_TIME        density-specific trajectory time cap [s]
DS_BOUNDARY_MODE        ISOTROPIC or ANISOTROPIC
```

For Mode3D density/flux, the direction set follows the gridless convention: a 24 x 48 direction grid is built, then optionally reduced deterministically if `DS_MAX_PARTICLES` is positive.

### 8.11 `#SPECTRUM`

The `#SPECTRUM` section defines the outer-boundary differential intensity used by density/flux and forward injection.

The parser stores the raw key/value section and also builds a typed spectrum view. The runtime spectrum evaluator is in:

```text
srcEarth/boundary/spectrum.cpp
srcEarth/boundary/spectrum.h
```

Typical spectrum inputs include a spectrum type, normalization, energy range, and optional file/time-dependent spectrum settings. Consult `boundary/spectrum.h` and `util/amps_param_parser.cpp` for the complete recognized key list.

### 8.12 `#BOUNDARY_ANISOTROPY`

Required when:

```text
DS_BOUNDARY_MODE ANISOTROPIC
```

Example:

```text
#BOUNDARY_ANISOTROPY
BA_PAD_MODEL          COSALPHA_N
BA_PAD_EXPONENT       2.0
BA_SPATIAL_MODEL      DAYSIDE_NIGHTSIDE
BA_DAYSIDE_FACTOR     1.0
BA_NIGHTSIDE_FACTOR   0.2
```

Recognized PAD models:

```text
ISOTROPIC
SINALPHA_N
COSALPHA_N
BIDIRECTIONAL
```

Recognized spatial models:

```text
UNIFORM
DAYSIDE_NIGHTSIDE
```

Anisotropic mode evaluates weights at the boundary exit point of allowed backtraced trajectories.

### 8.13 `#ENERGY_CHANNELS`

Optional. Adds integral-flux channels to density/flux output.

```text
#ENERGY_CHANNELS
CH_BEGIN
  P10_100      10.0     100.0
  P100_1000   100.0    1000.0
  P1G_10G    1000.0   10000.0
CH_END
```

Channels may overlap. Channels are clipped to the available density/flux energy grid.

### 8.14 `#DENSITY_3D`

Used by `-mode 3d_forward`, not by the transmissivity-based Mode3D density/flux solver.

```text
#DENSITY_3D
DENS_EMIN             1.0
DENS_EMAX             20000.0
DENS_NENERGY          30
DENS_ENERGY_SPACING   LOG
```

This section controls the energy bins used by forward-mode AMR cell-density sampling and inner-sphere flux sampling.

### 8.15 `#PARTICLE_TRAJECTORY`

Used by `-mode 3d_forward` particle trajectory record initialization.

```text
#PARTICLE_TRAJECTORY
INITIALIZE_TRAJECTORIES  T
N_TRAJECTORIES           1000
```

Trajectory output also requires the AMPS particle tracker to be enabled at compile time.

### 8.16 `#TEMPORAL`

Standalone Mode3D time-series example:

```text
#TEMPORAL
TEMPORAL_MODE    TIME_SERIES
EVENT_START      2024-05-10T00:00:00
EVENT_END        2024-05-10T12:00:00
FIELD_UPDATE_DT  15
TS_INPUT_MODE    FILE
TS_INPUT_FILE    ts05_driving_2024_05_10.txt
```

Coupled SWMF usage:

```text
#TEMPORAL
FIELD_UPDATE_DT  15
```

In standalone Mode3D, `FIELD_UPDATE_DT` is the magnetic driver snapshot spacing in minutes. In SWMF-coupled Mode3D products, it is the calculation cadence in minutes of simulation time.

---

## 9. Command-line options

Main usage:

```bash
./amps -h
./amps -mode gridless   -i AMPS_PARAM.in
./amps -mode 3d         -i AMPS_PARAM.in
./amps -mode 3d_forward -i AMPS_PARAM.in
```

Common options:

```text
-h, --help
    Print help and exit.

-mode <gridless|3d|3d_forward>
    Select standalone execution path.

-i <file>
    AMPS_PARAM input file.

-mover <BORIS|RK2|RK4|RK6|GC2|GC4|GC6|HYBRID>
    Select backward tracing mover for gridless/Mode3D, or the supported subset for 3d_forward.

-density-mode <ISOTROPIC|ANISOTROPIC>
    Override DS_BOUNDARY_MODE for density/flux calculations.

-max-trace-distance <Re>
    Override MAX_TRACE_DISTANCE.

-mode3d-output-initialized
    Write initialized 3-D field mesh diagnostic.

-mode3d-field-eval <INTERPOLATION|ANALYTIC>
    For standalone Mode3D, choose whether tracing uses mesh interpolation or direct analytic/background evaluation.

-density-parallel <OPENMP|THREADS|SERIAL>
    Select the shared-memory backend for Mode3D backward products within each MPI process. Despite the name, this option now applies to both cutoff and density/flux calculations. OPENMP preserves the OpenMP loops. THREADS uses direct std::thread workers over observation locations and suppresses nested OpenMP inside those workers. SERIAL disables intra-rank shared-memory parallelism.

-density-threads <N>
    Number of shared-memory workers per MPI process for Mode3D cutoff and density/flux products. For OPENMP this sets omp_set_num_threads(N). For THREADS this sets the number of std::thread workers. N=0 means automatic.
```

Forward-mode options:

```text
-forward-niter <int>
-forward-nparticles <int>
-forward-injection-energy <SPECTRUM|LOG_UNIFORM>
-forward-injection-emin <MeV/n>
-forward-injection-emax <MeV/n>
-forward-injection-energy-range <Emin> <Emax>
-forward-boundary-dist <ISOTROPIC>
-forward-track-trajectories
-forward-no-track-trajectories
-forward-n-trajectories <int>
```

---

## 10. Example inputs

### 10.1 Cutoff only, standalone Mode3D

```text
#CALCULATION_MODE
CALC_TARGET        CUTOFF_RIGIDITY
FIELD_EVAL_METHOD  GRID_3D

#CUTOFF_RIGIDITY
CUTOFF_EMIN            1.0
CUTOFF_EMAX            20000.0
CUTOFF_NENERGY         50
CUTOFF_MAX_PARTICLES   500
CUTOFF_SAMPLING        ISOTROPIC
DIRECTIONAL_MAP        F
```

Run:

```bash
mpirun -np 8 ./amps -mode 3d -i AMPS_PARAM.in
```

### 10.2 Cutoff plus density/flux from one standalone mesh snapshot

```text
#CALCULATION_MODE
CALC_TARGET        CUTOFF_RIGIDITY+DENSITY_SPECTRUM
FIELD_EVAL_METHOD  GRID_3D

#CUTOFF_RIGIDITY
CUTOFF_EMIN            1.0
CUTOFF_EMAX            20000.0
CUTOFF_NENERGY         50
CUTOFF_MAX_PARTICLES   500
CUTOFF_SAMPLING        ISOTROPIC
DIRECTIONAL_MAP        T
DIRMAP_LON_RES         30
DIRMAP_LAT_RES         30

#DENSITY_SPECTRUM
DS_EMIN            1.0
DS_EMAX            20000.0
DS_NINTERVALS      50
DS_ENERGY_SPACING  LOG
DS_MAX_PARTICLES   500
DS_BOUNDARY_MODE   ISOTROPIC
```

This produces cutoff, directional maps, density, local spectrum, and flux files from the same mesh field.

To use the direct-thread backend for both cutoff and density/flux in this run,
add to `#NUMERICAL`:

```text
#NUMERICAL
DENSITY_PARALLEL THREADS
DENSITY_THREADS  8
```

or use the standalone CLI override:

```bash
mpirun -np 4 ./amps -mode 3d -i AMPS_PARAM.in -density-parallel THREADS -density-threads 8
```

For the `THREADS` backend, each MPI rank attempts to widen its CPU affinity
mask before creating worker threads. If needed, override the detected node CPU
set, for example:

```bash
export AMPS_MODE3D_DENSITY_CPUSET=0-39
```

### 10.3 Density/flux only

```text
#CALCULATION_MODE
CALC_TARGET        DENSITY_SPECTRUM
FIELD_EVAL_METHOD  GRID_3D

#DENSITY_SPECTRUM
DS_EMIN            1.0
DS_EMAX            20000.0
DS_NINTERVALS      50
DS_ENERGY_SPACING  LOG
DS_MAX_PARTICLES   500
DS_BOUNDARY_MODE   ISOTROPIC
```

### 10.4 Standalone Tsyganenko time series

```text
#BACKGROUND_FIELD
FIELD_MODEL     T05
DRIVER_FILE     ts05_driving_2024_05_10.txt
EPOCH           2024-05-10T00:00:00

#TEMPORAL
TEMPORAL_MODE    TIME_SERIES
EVENT_START      2024-05-10T00:00:00
EVENT_END        2024-05-10T12:00:00
FIELD_UPDATE_DT  15
```

The requested products are repeated for every snapshot in the event window.

### 10.5 SWMF-coupled cutoff plus density/flux

In the AMPS/PT `AMPS_PARAM.in` used by the coupled run:

```text
#CALCULATION_MODE
CALC_TARGET        CUTOFF_RIGIDITY+DENSITY_SPECTRUM
FIELD_EVAL_METHOD  SWMF

#TEMPORAL
FIELD_UPDATE_DT    15

#NUMERICAL
DENSITY_PARALLEL   THREADS
DENSITY_THREADS    8
```

The coupled bridge reads `AMPS_PARAM.in`, waits for SWMF field updates, and runs the requested products at the configured cadence. The `DENSITY_PARALLEL` and `DENSITY_THREADS` settings are used by both the cutoff and density portions of the coupled backward-product calculation.

---

## 11. Output-file summary

### Gridless cutoff

```text
cutoff_gridless_points.dat
cutoff_gridless_shells.dat
cutoff_gridless_dir_map_point_####.dat
```

### Gridless density/flux

```text
gridless_points_density.dat
gridless_points_spectrum.dat
gridless_points_flux.dat
gridless_shell_<ALT>km_density_channels.dat
```

### Mode3D cutoff

```text
cutoff_3d_points.dat
cutoff_3d_shells.dat
cutoff_3d_dir_map_loc_000000.dat
cutoff_3d_points_dipole_compare.dat
cutoff_3d_shells_dipole_compare.dat
```

The dipole comparison files are diagnostic/benchmark products for DIPOLE runs.

### Mode3D density/flux

```text
mode3d_points_density.dat
mode3d_points_spectrum.dat
mode3d_points_flux.dat
mode3d_shell_<ALT>km_density_flux.dat
```

### Mode3D standalone time series

```text
*_snapshot_000000_YYYY_MM_DDTHH_MM_SS_000.dat
*_snapshot_000001_YYYY_MM_DDTHH_MM_SS_000.dat
```

### Mode3D SWMF-coupled snapshots

```text
*.swmf_n000000_t000000000.000s.dat
*.swmf_n000001_t000000900.000s.dat
```

### 3d_forward

```text
amps_3dforward_initialized.data.dat
sphere_flux3d.out=####.dat
sphere_flux3d_total.out=####.dat
```

and standard AMPS mesh outputs containing `Density3D` variables.

---

## 12. Validation and development recommendations

1. Start with `DIPOLE` and `CUTOFF_SAMPLING VERTICAL` when validating a new build.
2. Compare gridless and Mode3D with `-mode3d-field-eval ANALYTIC` to separate interpolation issues from trajectory-integration issues.
3. Then switch to `-mode3d-field-eval INTERPOLATION` to validate the AMR field materialization and interpolation path.
4. For time-series Tsyganenko runs, first run one or two snapshots before launching a long event window.
5. Use `DIRECTIONAL_MAP T` only when needed; it can multiply runtime by the number of sky-map direction cells.
6. In SWMF-coupled runs, verify that output suffixes correspond to the intended SWMF simulation times.
7. Treat `3d_forward` density/flux outputs separately from transmissivity-based outputs; they are produced by different statistics.

---

## 13. Known scope boundaries

The current Earth energetic-particle model is primarily a magnetic-access and test-particle transport tool. The backward products do not yet represent a full radiation-belt transport solver.

Not included in the backward transmissivity products:

- pitch-angle diffusion;
- radial diffusion;
- wave-particle scattering;
- Coulomb collisions;
- atmospheric energy loss;
- self-consistent particle feedback on fields;
- time-varying fields during a single traced trajectory.

The inner boundary is an absorbing sphere, not a detailed atmosphere/ionosphere model.

The anisotropic boundary model is factored and simplified:

```text
J_boundary(E, direction, x_exit) = J_iso(E) * f_PAD(direction, B) * f_spatial(x_exit)
```

It is useful for controlled SEP/radiation access studies but is not a full 5-D boundary distribution.

---

## 14. Troubleshooting

### `DIRECTIONAL_MAP T` does not write files

Check that the run is requesting cutoff:

```text
CALC_TARGET CUTOFF_RIGIDITY
```

or

```text
CALC_TARGET CUTOFF_RIGIDITY+DENSITY_SPECTRUM
```

Also check that both angular resolutions are positive:

```text
DIRMAP_LON_RES > 0
DIRMAP_LAT_RES > 0
```

### Standalone time series does not loop over snapshots

Check:

```text
TEMPORAL_MODE TIME_SERIES
EVENT_START   ...
EVENT_END     ...
FIELD_UPDATE_DT > 0
```

and verify that a driver file is loaded through `TS_INPUT_FILE` or `DRIVER_FILE`.

### Coupled SWMF run still behaves like forward injection

In coupled mode, backward products require explicit `CALC_TARGET`. If `CALC_TARGET` is absent, the coupled bridge preserves the historical forward-injection path.

### SPICE-related failures

Standalone time-series operation and frame transforms require SPICE. Rebuild without `_NO_SPICE_CALLS_` if the run needs UTC-to-ET conversion, driver-file interpolation, SM/GSM directional-map rotation, or trajectory frame conversion.

### Direct-thread backend creates workers but all run on one CPU

If `DENSITY_PARALLEL THREADS` creates the requested worker count but `top -H`
or `ps -L` shows that all worker threads run on the same `PSR` CPU, check the
rank affinity mask:

```bash
grep Cpus_allowed_list /proc/<pid>/status
ps -L -p <pid> -o pid,tid,psr,pcpu,state,comm
```

The direct-thread backend calls `PIC::Parallel::SetWideAffinityForScheduler()`
before creating cutoff/density workers. This is equivalent to:

```bash
taskset -apc <CPUSET> <pid>
```

and should widen the rank affinity mask when the scheduler allows it. The CPU
set is auto-detected from `/sys/devices/system/cpu/online`, but it can be
overridden with:

```bash
export AMPS_MODE3D_DENSITY_CPUSET=0-39
```

or:

```bash
export PIC_PARALLEL_CPUSET=0-39
```

If `Cpus_allowed_list` remains a single CPU, the batch system or MPI runtime is
still restricting the rank through a cpuset/cgroup, and the code cannot move
threads outside that allowed set.

### No useful density/flux values

Check:

- the outer domain is large enough for allowed trajectories to escape;
- `R_INNER` is reasonable;
- `MAX_TRACE_TIME`, `MAX_STEPS`, and `MAX_TRACE_DISTANCE` are not too restrictive;
- the boundary spectrum is initialized and covers the requested energy range;
- `DS_MAX_PARTICLES` is not so small that each energy gets too few directions.

---

## 15. Minimal command examples

Gridless cutoff:

```bash
mpirun -np 8 ./amps -mode gridless -i cutoff_gridless.in
```

Gridless density/flux:

```bash
mpirun -np 8 ./amps -mode gridless -i density_gridless.in -density-mode ISOTROPIC
```

Standalone mesh-backed cutoff plus density/flux:

```bash
mpirun -np 8 ./amps -mode 3d -i mode3d_products.in
```

Standalone mesh-backed cutoff/density products using direct threads inside each MPI rank:

```bash
mpirun -np 4 ./amps -mode 3d -i mode3d_products.in -density-parallel THREADS -density-threads 8
```

This backend applies to `CUTOFF_RIGIDITY`, `DENSITY_SPECTRUM`, and `CUTOFF_RIGIDITY+DENSITY_SPECTRUM`. Choose the MPI-rank count and thread count so that `ranks_per_node * density_threads` does not exceed the cores allocated per node.

In SWMF-coupled PT runs the standalone CLI is usually not used; the same backend can be selected from the input file with `DENSITY_PARALLEL THREADS` and `DENSITY_THREADS 8` in `#NUMERICAL` or `#DENSITY_SPECTRUM`, or with environment variables `AMPS_MODE3D_DENSITY_PARALLEL=THREADS` and `AMPS_MODE3D_DENSITY_THREADS=8`.

Standalone mesh-backed direct-field cross-check:

```bash
mpirun -np 8 ./amps -mode 3d -i mode3d_products.in -mode3d-field-eval ANALYTIC
```

Forward Monte Carlo 3-D run:

```bash
mpirun -np 8 ./amps -mode 3d_forward -i forward3d.in -forward-niter 1000 -forward-nparticles 1000
```

