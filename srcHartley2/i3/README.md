# 3I/ATLAS six-case radial campaign

Related documentation: [collisionless reference](../COLLISIONLESS_REFERENCE.md),
[edge-correction method](../EDGE_CORRECTION.md),
[MAVEN/IUVS model comparison](COMPARISON.md), and
[how to compile and run a given chemistry, case, and production or test build](HOW_TO.md). The comparison tables describe
the COMET1D campaign; production archives are separate from this source tree.

This directory owns the observation-derived data and deterministic runtime
input generation for the six MAVEN/IUVS 3I/ATLAS cases supplied in the
source-controlled test artifact
`srcHartley2/tests/data/i3/Archive.zip`. It does not put a 3I constant into
the generic comet C++ application.

Run:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 srcHartley2/i3/prepare_i3_cases.py
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest \
  srcHartley2/i3/tests/test_prepare_i3_cases.py
```

The first command reproduces the single loadable `input/i3.in` file from the
reviewed CSV files under `reference/`. Common model settings appear once;
six direct `case: 1` through `case: 6` sections contain only the values that
vary with the observation. Byte-exact regeneration is required for review.

The executable is compiled once from `input/i3.conf`:

```sh
./Config.pl -application=i3
rm -rf -- build
make -j8 EXTRACOMPILEROPTIONS=-O3
```

`Config.pl` first looks for `input/i3.input`; because this application has no
such file, it selects `input/i3.conf`. The `.conf` file fixes species,
dimension, compiled model/reaction/collision topology, block layout, particle
storage and allocation maxima. It contains no observation case. The runtime
file supplies physical sources, case reaction rates, mesh targets, run
controls, fluorescence geometry and solar conditions:

```sh
mpiexec -n 12 ./amps -i input/i3.in -case 1
```

`-case` is mandatory because `i3.in` contains case sections. Omitting it or
requesting a name outside `1`--`6` terminates before mesh or particle
allocation. All six cases use the same executable; no rebuild is required
when only the selected case or runtime values change. Use a distinct run
directory (or archive its `PT` tree) for every launch.

For a fast executable-level receipt check, add `--validate-runtime-input`.
It resolves the selected case, verifies that all MPI ranks obtained the same
canonical configuration, runs the complete comet/chemistry/collision/cooling/
brightness/radial validation chain, and exits before mesh or particle
allocation:

```sh
mpiexec -n 2 ./amps -i input/i3.in -case 1 --validate-runtime-input
```

## Consistency with MAVEN/IUVS: campaign `qfit_20261123`

Added 2026-11-26 at the user's request. The I3 model is consistent with the MAVEN/IUVS observations of OH
(0-0), H Lyman-alpha and O I 130.4 nm for all six cases. This holds with the OH chemistry that was
implemented later, `OH_LOSS = SA1982_VH` (2026-11-21, section "OH loss variants"), and Q(H2O) fitted to the
IUVS OH.

The table gives edge-corrected model / observation for the converged runs (70000-100000 iterations), as the
geometric mean over the MAVEN/IUVS points:

| Case | 1 | 2 | 3 | 4 | 5 | 6 |
|---|--:|--:|--:|--:|--:|--:|
| Q(H2O) (s^-1) | 1.853e28 | 3.775e28 | 8.893e28 | 2.102e29 | 1.527e29 | 1.483e29 |
| OH (0-0) | 1.03 | 1.03 | 1.05 | 1.00 | 1.03 | 1.02 |
| H Lyman-alpha | 0.62 | 0.79 | 1.15 | 1.08 | 0.78 | 0.84 |
| O I 130.4 nm at 1e5 km (R) | 0.161 | 0.504 | 4.45 | 92.0 | 56.4 | 68.9 |

- **OH:** Q(H2O) was fitted to OH, so OH agrees by construction.
- **H Lyman-alpha:** individual points range from 0.47 to 1.91 of the observations.
- **O I 130.4 nm:** case 5 is the only case with an IUVS O profile. IUVS / model is 0.60-1.23 per bin, with
  geometric mean 0.87. The O values of the other cases are model predictions.
- **O model:** direct (KVB2022) + Ly-beta Bowen, optically thin; see "Chemical model and g-factors in
  detail".

Everything needed to rebuild and rerun this campaign from git on another machine is in
`campaigns/qfit_20261123/`. Its `README.md` gives the steps, resources and reference results. The folder
holds:

- the deck, which carries a comment with this statement;
- the build input, `BlockCells=12`;
- scripts to build, validate, run and compare.

The default deck `input/i3.in` keeps the earlier `STATE_RESOLVED` OH loss with the generator's source rates.
Both variants remain runnable for every case, with either deck: select one with
`--chemistry-variant OH_LOSS=...` (section "OH loss variants").

## Nightly regression test

The nightly test `i3` (block `Name=i3` in `MakefileTest/Table`, family
`comet-dust`) runs this application end to end and compares its output with
references:

```sh
make test_i3                            # compile, run directory, run, check (4 MPI processes)
make test_i3 TESTMPIRUN4="mpirun -np 8" # the same with the second reference layout, N = 8
cat 'test_i3[5].diff'                   # empty: the test passed (test_i3.diff without a tag)
```

Since 2026-11-23 the runs use `$(or ${TESTMPIRUN4},mpirun -np 4)`, like the
other Table tests that need a fixed layout, and no longer `${MPIRUN}`. With
`${MPIRUN}` a plain `make test_i3` ran with `NP=2` (`Makefile.def`), and every
result was reported as `reference ... ref_np2.gz is missing`. References exist
only for N = 4 and N = 8. The nightly scripts set `TESTMPIRUN4="mpirun -np 4"`
(`utility/TestScripts/AMPS-gpu/RunGNU.sh`), so the nightly test uses N = 4.

- **Build.** `input/i3.conf` built in test mode (`-amps-test=on`), like the
  other nightly tests, with a removed `build/` directory. Test mode switches
  on the debugger checks (particle lists, particle positions, block
  crossings) and sorts the particle lists after every particle exchange, so
  its results differ from those of a production build. Until 2026-11-16 test
  mode aborted every DIM1 run: the particle-list check
  (`src/pic/pic_pbuffer.cpp`) compared three particle coordinates with
  `node->xmin`/`node->xmax`, which have `_MESH_DIMENSION_` elements, and the
  block-crossing check of the radial mover (`src/pic/pic_mover.cpp`) compared
  the Cartesian x coordinate, not the radius, with the radial block limits.
  Until 2026-11-23 the DIM1 build also gave 47
  `-Waggressive-loop-optimizations` warnings in `src/pic`. They came from
  3D-only code, which indexes `xmin`/`xmax`/`downNode` past their
  `_MESH_DIMENSION_` size, and none of it is reached in i3. Since then:
  - the child-node loops use `(1<<_MESH_DIMENSION_)`;
  - the bodies of 25 3D-only functions are compiled only for
    `_MESH_DIMENSION_ == 3`, and otherwise stop with `exit()`;
  - `Split::Scatter`'s even-block check is done by the preprocessor.

  The build is free of warnings, and the i3 outputs were byte-identical before
  and after this change. Record:
  `test_output/i3/rerun_2026-11-23_baseline/compile_warnings_proposal.md` and
  `test_output/i3/warning_fix_2026-11-23/`.
- **Input.** `input/i3-test.in` holds cases 1 and 5 of `input/i3.in`
  unchanged; the reduced nightly entry selects case 1 in two OH-loss variants.
  Case 5 remains available in the deck for explicit runs. Only three run
  controls of the common block differ:
  `simulation.reference_injection_particles = 8` (59) and
  `simulation.max_iterations = 10000` (150000), with
  `simulation.sampling_length = 10000` as before.
- **Run length (test-mode build).** Like the other applications (e.g. the 3D
  driver `srcHartley2/main.cpp`, `srcCG/main.cpp`), the test-mode build sets
  the run length in the code, not in the input. With
  `_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_` (`-amps-test=on`),
  `Comet::RadialSimulation::ApplyNightlyTestRunControls()` (called in
  `RadialMain.cpp` before `RadialSimulation::Validate()`) replaces the
  runtime-input values with the named constants of `RadialSimulation.h`:
  `NightlyTestSamplingLength = 2000`, `NightlyTestMaxIterations = 2000` and,
  when the input interval is positive, `NightlyTestWeightAdaptationInterval =
  2000`. Each case thus runs one sampling window, writes output `out=0` at
  iteration 2000, and then makes one regular weight adaptation (interval =
  window, as in `input/i3.in`). The log line `nightly test mode
  (_PIC_NIGHTLY_TEST_MODE_=on): run controls replaced: ...` gives the
  executed and the runtime-input values. A production build is unaffected.
  So the 10000-iteration controls of `input/i3-test.in` apply only when the
  deck is run with a production build. Until 2026-11-23 the test ran 10000
  iterations. The references were then regenerated for 2000 (see
  "References" below).
  `prepare_i3_cases.py` does not write this file; it is derived from the
  generated `input/i3.in`. It predates the O I 130.4-nm brightness records
  added to `input/i3.in` on 2026-11-19 (three emission records and two solar
  fluxes per case, see "O I 130.4-nm brightness" below) and does not contain
  them; its references therefore hold the H and OH brightness only. The O
  records are optional for the runtime parser, so the deck remains valid.
  Since 2026-11-21 it also has the two OH-loss keys of `input/i3.in`
  (`chemistry.variant.OH_LOSS = STATE_RESOLVED` and
  `chemistry.table.OH_TAU_SA1982.file`, see
  [CHEMISTRY_VARIANTS.md](../CHEMISTRY_VARIANTS.md)). The table path
  `../srcHartley2/i3/reference/oh_lifetime_sa1982.csv` is relative to the
  directory of the deck, and is valid both in `input/` and in the copy in
  `run_test_i3/`; the output headers record it in that relative form, so the
  references do not depend on where the tree is checked out.
- **Run.** The two runs both select case 1. `run_test_i3/case1` uses the
  default OH loss `STATE_RESOLVED`; `case1sa` uses the alternative selected on the
  command line, `--chemistry-variant OH_LOSS=SA1982_VH`, in
  `run_test_i3/case1sa`. Each run's log is in `amps.log` there, and the
  results are copied unchanged to `run_test_i3/res/<run>.<species>.dat` and
  `run_test_i3/res/<run>.<species>.brightness.dat`. That makes 18 files: the
  seven radial profiles and the H and OH brightness profiles of each run.
- **Check.** Each result is compared whole (`share/Scripts/DiffNum.pl`,
  default tolerances) with
  `output/test_i3/<case>.<species>[.brightness]<tag>.ref_np<N>.gz`. These
  references live in SWMF_data under `PT/AMPS/data/output/test_i3`. They are
  the unchanged output files, and N is the number of MPI processes in the
  case's `amps.log`.
  - References stay compressed in SWMF_data. The check uses `gunzip -c` to
    write each unpacked copy to `run_test_i3/res/<output>.ref`, preserving
    the original `.gz` file. A missing or corrupt gzip reference is recorded
    in the corresponding `.diff`; corrupt data are not passed to `DiffNum.pl`.
  - No column is removed. In test mode every column is reproducible: two
    runs per N gave identical files (for N = 4 with two builds), including
    column 3, `Temp_ID` (the allocation counter of the mesh nodes), and
    columns 12 and 13.
  - Columns 12 and 13 are the parallel and tangential translational
    temperatures. Their sampling direction (`src/pic/pic.cpp`) points from the
    origin to the cell centre, whose inactive Cartesian components are zero in
    DIM1. The DIM1 spherical movers keep every particle at x = (r, 0, 0) with
    v[0] the radial velocity, so the parallel temperature is the radial one.
    The sampled cell values satisfy T = (T_parallel + 2 T_tangential)/3; the
    output values, interpolated to the cell corners, agree with it where the
    cells are well sampled, and far from the nucleus T_tangential goes to zero
    (radial free flow).
- **References.** The results depend on N (each rank has its own random
  stream) and on the build, and a run is bit-for-bit reproducible for a given
  executable and N. The references are therefore made per machine and
  compiler. Those of passion (GCC 9.4, Open MPI 4.0.3) carry the tag `[5]`:
  the block's `Ref=5{gcc}@@(0+passion)` sets `$(TESTi3-REF)=[5]` there, as for
  the other tests on passion. Every other machine uses the untagged names, so
  the references of a second machine can be added as untagged files without
  touching those of passion. As of 2026-11-23, however, no untagged `test_i3`
  reference has ever been committed to SWMF_data (none is in its history). On
  any machine other than passion, every result is therefore reported as
  `reference ... is missing` and is not compared.

  **2026-11-23 regeneration (2000-iteration nightly run length).** The 36
  passion files used by the test (`case1`, `case1sa`; N = 4 and 8) were
  regenerated after the nightly run length moved into the code
  (`_PIC_NIGHTLY_TEST_MODE_`, see "Run length" above). Evidence, in this order:
  1. Before the change, a fresh `make test_i3` with N = 8 reproduced the
     previous 10000-iteration references byte for byte (18/18).
  2. The new build (one executable, SHA-256 `5ab92edc…bbe8a51`) gave
     byte-identical results in two independent runs per layout (18/18 for
     N = 4 and for N = 8).
  3. The pre-change executable, run with only the three run controls of a
     copy of `input/i3-test.in` set to 2000 (no line moved), reproduced the
     new N = 4 results byte for byte (18/18). The code change therefore only
     replaces those three run controls.
  4. After installation, `make test_i3_check` passed on all four runs,
     including the two independent repeats.

  No rates, seeds, model inputs or comparison tolerances were changed. The 18
  `case5` files were not touched. The previous files, their SHA-256 sums, the
  logs and the scripts are kept in
  `test_output/i3/refgen_nightly_macro_2026-11-23/` (not committed).
  At iteration 2000 the regular weight adaptation runs but changes no weight
  (adaptation log at iteration 2000: `particles=47544->47544 split=0
  roulette=0`, N = 4). The 10000-iteration
  test exercised actual splits (482 at iteration 10000). The shorter nightly
  run therefore no longer covers a weight change.

  **2026-11-23, second regeneration (comment update of `input/i3-test.in`).**
  The brightness headers record the input line of
  `chemistry.variant.OH_LOSS` (`source=runtime input ../i3-test.in:NN`). The
  comment-only update of the deck moved that key from line 54 to line 71, so
  `case1.H.brightness` and `case1.OH.brightness` changed in that header line.
  The 36 references were therefore made again, with the build that includes
  the DIM<3 warning fix. Evidence:
  1. With the unchanged deck, that build reproduced the installed references
     byte for byte (N = 4 and 8, 18/18 each).
  2. Two independent runs per layout with the updated deck were byte-identical.
  3. Against the installed references, the only difference was `:54` → `:71`
     in those two files. The other 16 outputs per layout were byte-identical.

  Record: `test_output/i3/refgen_comment_update_2026-11-23/`. Adding or
  removing lines above that key in `input/i3-test.in` requires new references.

  The passion references were previously regenerated on
  2026-11-22 from `b4980503e35c498d42a689188e0f6fc2e4913d94`, including
  the DIM1 domain-boundary callback fix `848e9c5`, for all 27 outputs and
  N = 4 and 8 (54 gzip files: 36 replacements and 18 new `case1sa` files).
  A clean test-mode build and the unchanged `input/i3-test.in` were used.
  Independent fresh repeats on both layouts reproduced all files byte for
  byte and passed the standard `make test_i3_check` comparisons, including
  every column and all brightness chemistry-provenance headers. No rates,
  seeds, sampling controls or comparison tolerances were changed. Previous
  references, inputs, executable identity and logs are retained under
  `test_output/srcHartley2/test_i3_references_20261123T004840Z/`.
  The reduced nightly suite uses only `case1` and `case1sa`, so it requires
  18 of those files per layout (36 for N = 4 and 8 together). Case5 references
  from the full generation remain archived and are not selected by the test.
  The untagged references (Voyager: GCC 11.5.0, Open MPI 4.1.7rc1, AMD EPYC
  7763) were remade on 2026-11-21 with the test-mode build for all 27 outputs
  and N = 4 and 8: two runs per N gave byte-identical results, a third check
  run per N passed `make test_i3_check`, and the radial profiles of cases 1
  and 5 are byte-identical to those of the code before the chemistry variants
  (the brightness files differ only in the `Chemistry.*` provenance lines).
  The earlier untagged references, made with the production build, were moved
  aside (not deleted). Record of that work:
  `validation_runs/test_i3_sa1982_20261121/` (not committed).
  That untagged generation predates `848e9c5`; it needs its own regeneration
  on the corresponding machine/compiler. The 2026-11-22 passion results do
  not replace or qualify untagged references.
  To regenerate the references of one machine and one N after a `make test_i3`
  run with that N (`T` is the machine's tag, `[5]` on passion, empty
  elsewhere; N = 4 is `make test_i3`, N = 8 is
  `make test_i3 TESTMPIRUN4="mpirun -np 8"`):

  ```sh
  N=4; T='[5]'
  for f in run_test_i3/res/*.dat; do
    gzip -n -9 -c $f > "data/output/test_i3/$(basename $f .dat)$T.ref_np$N.gz"
  done
  ```

- **Timing of the 2000-iteration test (passion, 2026-11-23).** The complete
  `make test_i3` (N = 4) took 2.6 min. The compilation (8 jobs) took about
  1.9 min, and the two runs 42 s together. With N = 8 the two runs took about
  24 s. The adaptation log line at iteration 2000 reports `particles=47544`
  (N = 4; 240796 at iteration 10000 in the N = 8 baseline).
- **Historical timing measurements (10000 iterations).** These include older
  case5 workloads; the current Table runs only case1 and case1sa. On passion (Intel i9-10980XE,
  18 cores, shared with other jobs;
  2026-11-16) the compilation takes 1.6 min, and the run with 4 MPI
  processes 20 to 22 min (case 1 8 min, case 5 12 to 14 min). With 8
  processes it took 19 and 30 min (case 1 6.4 and 10.3 min) at a load
  average of 40 to 50. On an AMD EPYC 7763 the production build took 3.8 +
  6.7 min with 4 processes and 8 min in all with 8. On Voyager (AMD EPYC 7763,
  2026-11-21, test-mode build, four test runs at once on one node) the
  compilation took 1 min with 32 jobs; with 4 processes case 1, case 5 and
  `case1sa` took 5.4, 8.6 and 4.9 min (19 min in all), with 8 processes 3.8,
  5.0 and 3.7 min (12 min in all). There are about 1.9e5
  model particles at iteration 10000. Each MPI process uses about 1.5 GB of
  memory on passion, whatever the particle number, because `input/i3.conf`
  requests a buffer of 2e7 particles.

## Archived observation table

The supplied archive has SHA-256
`344d33e69c913a5e76d57ddb7c92c1a8d1e52a92a03c06f2f5849f70e45597be`.
Its canonical repository location is
`srcHartley2/tests/data/i3/Archive.zip`; the portable orchestrator verifies
that hash before regenerating the numerical brightness reference.
Its embedded `Model_Table_6cases.pdf` has SHA-256
`a88eadce0b9dbdd3f4fda0974ff00b6a791d516ff6be497e732b57bdb7f34e0c`.
`reference/i3_cases.csv` transcribes its exposure times, heliocentric distance,
radial velocity, spatial resolution, and table production rates. Rates are
interpreted as molecules per second. H2O values inferred from the raster red
curves are retained in the separate `plot_fit_q_h2o_s-1` column for future
comparison and are not used as source rates.

The observations are in 2025. This is independently consistent with the NASA
mission report; the archive file creation date is not an observation date.

## Solar-activity correction and reaction rates

The rate calculation implements the user-selected Moore and Mapaye (2025)
per-channel law

```text
j_1AU = A0 * F10.7P^alpha,
F10.7P = (daily adjusted F10.7 + centered 81-day mean)/2.
```

For every exposure endpoint, Sun-centered J2000 ecliptic vectors for Earth
and 3I were obtained from the NASA/JPL Horizons API (3I packed designation
`K25N010`). The same solar face was found with Carrington sidereal rotation
`360/25.38 deg/day`; the resulting Earth-equivalent interval was averaged
using official NRC adjusted 20-UT daily flux. `solar_activity_audit.csv`
retains endpoint shifts, unshifted controls, final shifted values, and the
method. The source snapshots used during derivation had these hashes:

- NRC daily flux table: `5962ab4b3dda01beab622fcfd75634908aca00b716bfdb6c2f22d5d10119330a`;
- Horizons 3I vectors: `82c17649f9a06df8c6542cc7f741ead0673531d088b113a58cce18a660710aab`;
- Horizons Earth vectors: `d809d1ab766343839afbde35d2ff7a5fe8a1f589a5dd9ea742e15191e3febd93`.

The authoritative services are
`https://www.spaceweather.gc.ca/solar_flux_data/daily_flux_values/fluxtable.txt`
and `https://ssd.jpl.nasa.gov/api/horizons.api`. No runtime network access is
needed. The generated runtime chemistry records contain the calculated 1-AU value;
the existing `HELIOCENTRIC` law then applies the case's inverse-square factor
once.

`photochemistry_channels.csv` contains the full neutral H2O descendant
network selected for this study: three H2O dissociation channels, three
state-labeled OH channels, two H2 dissociation channels, and neutral-system
photoionization/dissociative-ionization sinks for H2O, OH, H2, H, and O. O
electronic products are collapsed into the compiled O species, but reaction
IDs, rates, and excess energies remain distinct. H2O/H2 energies are from
Huebner and Mukherjee (2015) as tabulated by Carberry Mogan et al. (2022).
OH excess energies are state-labeled published values: 1.27 eV for O(3P)
from Bryans and Pesnell (2012), and 7.73/10.0 eV for O(1D)/O(1S) from
Bisikalo et al. (2015), Table 1. Bryans and Pesnell list the close 7.90/9.80
eV alternatives; the generator deliberately uses one cited table rather than
silently averaging the sources.

CO and CO2 are direct collision partners only. No CO/CO2 photochemistry is
invented.

## OH loss variants

The OH loss is a chemistry group with two variants, both compiled into the I3
executable (design, decisions, and tests: `srcHartley2/CHEMISTRY_VARIANTS.md`):

- `STATE_RESOLVED` (default): `OH_O3P_H`, `OH_O1S_H`, `OH_O1D_H`, and `OH_ION`
  with the per-case F10.7P-scaled rates above. All campaigns before
  2026-11-21 used it.
- `SA1982_VH`: one OH -> O + H channel (`OH_O_H_SA1982`, 0.36 eV) with
  k_d = 1/tau_FUV + 1/(S * OHTAUS(v_h)) at 1 AU and OH ionization
  (`OH_ION_SA1982`) at 0.05 k_d, both scaled by (1 AU / r_h)^2 and not by
  solar activity. OHTAUS and tau_FUV = 4.24e5 s come from
  `reference/oh_lifetime_sa1982.csv` (the legacy `BLOCKDATA OHPHOT` table,
  -60..+60 km/s in steps of 2 km/s, linear interpolation, clamped outside);
  S = 0.833e5 s and the fraction 0.05 are fixed in `input/i3.conf`. v_h is the
  case's `brightness.heliocentric_radial_velocity_km_s`. For the six cases the
  OH lifetime at 1 AU is 1.7-2.5 times that of `STATE_RESOLVED`.

`prepare_i3_cases.py` writes `chemistry.variant.OH_LOSS` (option `--oh-loss`,
per case `--oh-loss-case CASE=VARIANT`) and the table path
(`--oh-lifetime-table`, written relative to the output directory) into
`input/i3.in`, and always writes the state-resolved rates, so either variant
can also be chosen at run time with `--chemistry-variant OH_LOSS=...`. Another
solar-activity condition is represented by another table file of the same
format named in the input. `input/i3-test.in` has no variant keys and runs the
default. The edge correction (`srcHartley2/visualization/comet_edge_correction.py`)
takes the variant recorded in the run's brightness header.

## Source distribution, mesh, and collisions

All parents use a uniform 175-K, zero-drift surface source. This does not
sample an ordinary volume Maxwellian. The AMPS boundary injector samples

```text
p(v_vector) proportional to v_normal f_M(v_vector), v_normal > 0.
```

Consequently the speed marginal is `v` times the ordinary Maxwell speed
distribution and the angular marginal is the cosine law `p(mu)=2*mu`.
Residence weighting recovers the outgoing half-volume-Maxwellian used below.

For parent `j`, the surface number density is

```text
n_j = Q_j / [4*pi*R_n^2 * sqrt(2*k*T/(pi*m_j))].
```

For each direct parent `i`, the generator evaluates

```text
lambda_i = 1 / [sqrt(2) * sum_j(n_j*sigma_ij)]
```

with the exact collision records selected by the compiled capability at
175 K. The sampled-temperature floor and occupancy threshold are runtime
model settings in `i3.in`. The requested first-cell width is the minimum of
the H2O, CO, and CO2 mixture paths (CO2 is limiting in every case). It is
floored 1% above the finest cell the compiled AMR tree allows, which is the
radial extent / (BlockCells x 2^MaxMeshRefinmentLevel) = 0.349 m; AMPS
aborts on a finer request. At the current rates, case 1 is resolved at
0.444 m and cases 2--6 (mean free paths 0.076--0.27 m) are floored at
0.353 m, and the generated `i3.in` records each mean free path. The domain's
cometocentric outer radius is 3e9 m (3x10^6 km), with 2,000 ideal radial
cells and a 2.5-km spherical nucleus. The first campaign used 1e9 m, which
removed 14--37% of the H column at the median observation point.

### Production rates

`reference/i3_cases.csv` keeps the archived table rates unchanged. The
generator multiplies all three parent rates of each case by
`parent_rate_factor` in `reference/i3_rate_scaling.csv`. The current factors
(3.95, 4.81, 2.97, 1.74, 1.25, and 1.53 for cases 1--6) match the MAVEN/IUVS
OH(0-0) profiles to the converged AMPS OH, keeping the archive CO/H2O and
CO2/H2O ratios. The derivation is in
`validation_runs/i3_campaign_20261108/analysis_q_scaling/README.md`.

These factors are provisional. They were fitted to runs made before the DIM1
mover fix of 2026-11-22 (`src/pic/pic_mover.cpp`). In those runs the mover
deleted the molecules that returned to the nucleus instead of re-emitting
them thermally, so only 77--80% of the input H2O rate reached the coma; this
was a defect, not a physical return flux. With the fix, every returning
molecule is re-emitted and the steady net outflow is the input rate. The
same factors therefore supply about 1.25--1.3 times as much H2O to the coma
as in the fitted runs, and they must be refitted to runs with the fix. The
evidence is in `validation_runs/inner_boundary_probe_20261121/README.md`.
For the `SA1982_VH` OH loss, the rates refitted to runs with the fix are in the
deck of `campaigns/qfit_20261123/` (section "Consistency with MAVEN/IUVS").

### Time step

Each species' characteristic speed is a conservative upper bound on its
particle speed: 2,500 m/s for the direct parents (the first campaign reached
bulk + 3 sigma = 2.2 km/s), and parent bound plus maximum recoil for each
product. The Courant number of 0.8 therefore limits the fastest particle of
each species to 0.8 cells per step. Compared with the first campaign (400 m/s
bounds at Courant 0.25), the time step is 3.1 times longer for H, 2.9 for H2,
2.2 for O, and 1.6 for OH, which shortens the slow transit of collisionally
slowed H; the parent time step is 0.51 times as long.

The collision hierarchy is explicit. Three CO2-parent pairs retain Rubin et
al. (2011) Table 3. The selected completion adds Combi (1996) H2O/OH/H2/H/O/CO
values, followed only where needed by CHEMKIN-diameter/Lorentz hard-sphere
estimates. The unsupported OH--CO2, O--CO2, H2--CO2, and H--CO2 pairs remain
labeled estimates in the generated compiled capability table. Collision-pair
topology is compile-time because it defines the complete callback database;
no case-specific collision value is hidden in the executable.

The direct-source injection reference is derived for at least 800 particles
per first-cell parent population and is 59 events per parent time step. Regular block-local adaptation keeps the validated 800--3,200
range for directly represented populations; chemically produced species use
the separate 50--200 range, and direct-source weights remain locked outside
emergency common-factor coarsening. These are statistical controls and do not
alter the physical production rates above.

## Brightness scope and limitations

Each case requests 200 logarithmic impact parameters from the nucleus
surface to the finite domain edge (3x10^6 km since campaign 2). The input selects the reusable core
models `OH_AX00_SA1988_LOWELL_UNQUENCHED` and `H_LYA_IKL2020`, and supplies
the observation's heliocentric radial velocity, explicit Sun/comet/observer
geometry, and 1-AU LASP V4 Ly-alpha photon flux. The Lowell source quantity
is an integrated energy rate; each case explicitly selects 308.5 nm to
reproduce the archived photon conversion.

Brightness is sampled from individual model particles during every AMPS
sampling iteration. At each line-of-sight quadrature point, the DIM1 local
three-component velocity is rotated into the physical radial direction and
combined with the comet velocity before evaluating the heliocentric radial
Doppler speed. Thus H Ly-alpha uses the simulated velocity distribution; the
archived table's scalar comet radial velocity only supplies the bulk vector
`V_comet=Vr*e_r`. The velocity-dependent core model is evaluated per particle
without clamping or extrapolation. A conservative line-center cross section
`9.5e-18 m2` with threshold `tau=0.1` flags where optically thin output needs
full frequency-dependent transfer. It does not attenuate the result.

The geometry vectors in `reference/i3_cases.csv` were derived at each
exposure midpoint from geometric J2000-ecliptic JPL Horizons states for the
comet and Mars center: `e_r` points from the Sun to the comet, `LOS` points
from the comet to Mars, and the positive impact direction is the projection
of `e_r` onto the observer sky plane. The vectors are explicit case data, not
compiled constants. They are unit normalized and `LOS dot IMPACT=0` is
enforced by both the generator and runtime parser. The archived scalar radial
velocities remain authoritative; current-solution Horizons velocity was not
substituted for those values.

The frozen LASP source subset used for derivation had SHA-256
`672c9577fad0019d8ee7a6a8198604afd4671a5895ef66ee3e0cc544aee9a5db`.
The authoritative dataset is the LASP/LISIRD Version-4 daily composite
(Machol et al. 2019). OH factors follow Schleicher and A'Hearn (1988), using
the Lowell Observatory service. A Sun-directed H radiation-pressure force is
not applied in exact spherical symmetry; replacing it with a cometocentric
radial force would be physically incorrect. Particle-resolved brightness is
written as ordinary Tecplot data using
`PT/plots/pic.<species>.brightness.s=<species-index>.out=<output-index>.dat`.
The integral stops at the finite domain and performs no unreported outer-tail
continuation. The far-coma (edge) correction used in the comparison plots is
applied only in post-processing; see "What does 'corrected for the 3e9-m edge
(collisionless reference)' mean?" below.

### O I 130.4-nm brightness

Since 2026-11-19 `input/i3.in` also samples the atomic-oxygen 130.4-nm
triplet (130.217 + 130.486 + 130.603 nm), as three additive components of
the same particle-resolved sampler. They are written as three column groups of
`PT/plots/pic.O.brightness.s=2.out=<k>.dat`, each with its own column
density (identical), effective g factor, brightness, and optical-depth flag:

| record | transition label | core model | process |
|---|---|---|---|
| `brightness.emission.O_DIRECT_KVB2022` | `OI_1304_TRIPLET_DIRECT_KVB2022` | `O_I_1304_KVB2022_DIRECT_J2` | direct solar pumping at 130.217 nm, Gladstone (1992) line shape |
| `brightness.emission.O_DIRECT_SUMER2001` | `OI_1304_TRIPLET_DIRECT_SUMER2001` | `O_I_1304_SUMER2001_DIRECT_J2` | the same, SOHO/SUMER atlas line shape |
| `brightness.emission.O_LYB_BOWEN` | `OI_1304_TRIPLET_LYB_BOWEN_SUMER2001` | `O_I_1304_SUMER2001_LYB_BOWEN_J2` | solar Lyman-beta pumping at 102.576 nm and cascade (Bowen fluorescence) |

The total 130.4-nm brightness is one direct component plus the Bowen
component. The two direct components bracket the solar line wings, which no
available full-disk measurement resolves. At the comet's radial velocity, the
1-AU g factors (photons per second per O atom in 2p4 3P2) are:

| case | rdot (km/s) | direct, KVB2022 | direct, SUMER2001 | Lyman-beta (Bowen) |
|---|---|---|---|---|
| 1 | -39.7 | 3.3e-8 | 1.7e-6 | 3.9e-7 |
| 2 | -36.2 | 1.4e-7 | 2.7e-6 | 4.1e-7 |
| 3 | -29.6 | 1.3e-6 | 4.7e-6 | 4.4e-7 |
| 4 | -5.8 | 1.3e-5 | 1.3e-5 | 4.2e-7 |
| 5 | +0.4 | 1.0e-5 | 1.3e-5 | 4.2e-7 |
| 6 | +7.1 | 1.4e-5 | 1.1e-5 | 4.2e-7 |

Near rest (cases 4--6) the two line shapes agree within 30% and the Bowen term
is 3--4% of the total. In cases 1--3 the comet samples the far wing of the
solar lines: the case-1 total is 0.42--2.1e-6 and the Bowen term is a large or
dominant part of it, so the O brightness of these cases carries a factor of
2--5 line-shape uncertainty. See `src/models/resonance_fluorescence/README.md`
for the sources, the 2p4 3P2 ground-level assumption, and the exclusions
(electron impact, prompt emission, optical depth).

Two solar inputs per case come from `reference/solar_fuv_lines.csv`, written by

```sh
PYTHONDONTWRITEBYTECODE=1 python3 srcHartley2/i3/derive_solar_fuv_lines.py          # from the frozen snapshot
PYTHONDONTWRITEBYTECODE=1 python3 srcHartley2/i3/derive_solar_fuv_lines.py --check  # is the CSV current?
```

The source is the LASP/LISIRD FISM2 daily spectrum (`fism_daily_hr`, 0.1-nm
bins), averaged over the same Earth-equivalent interval as F10.7P and Lyman
alpha. `brightness.solar_oi1304_triplet_flux_1au_photon_cm-2_s-1` is the
FISM2 130.1--130.7 nm flux relative to 2009-03-01 (the SORCE SOLSTICE day of
the Killen et al. table), 1.21--1.27, times the triplet flux that table
implies (8.41e9 photons cm^-2 s^-1); the ratio removes the FISM2/SOLSTICE
calibration offset. `brightness.solar_lyb_flux_1au_photon_cm-2_s-1` is the
absolute FISM2 flux in 102.4--102.8 nm, 4.82--4.96e9 photons cm^-2 s^-1.
The frozen LISIRD responses (retrieved 2026-11-19) are in `reference/fism2/`:

- `fism2_daily_hr_1300_1308A_20251001_20251031.csv`: `57d238957a16a6e589038e4ef0aee1b3df987934c477e998e657ef1de4f97634`;
- `fism2_daily_hr_1300_1308A_20090301.csv`: `48d69005b640880277fd33b3b4d926f731e88c381515a4b5eb9c1dabea0eeafb`;
- `fism2_daily_hr_1024_1028A_20251001_20251031.csv`: `291fbd8a1f2a7d3ba555347fc12e09001806665202d8da2c6fd6d069e73ada1c`;
- `fism2_daily_hr_1024_1028A_20090301.csv`: `8eaabfa43b1b3d2dac738392dd0e78223e2b581dd4f75bd11768ca842edc3953`.

The O line-centre optical-depth diagnostic uses the 130.217-nm line
(f = 0.0520) with a 1 km/s Doppler parameter, 1.0e-17 m^2, and flags
tau > 0.1. An O column of 1e17 m^-2 already reaches tau = 1, so the flag is
expected near the nucleus of the high-production cases; as for H, the
output stays optically thin and the flag marks where it is not valid.

### Shared velocity-dependent fluorescence data

The authoritative H, OH, and O reference data and evaluators now live in the
application-independent core module
`src/models/resonance_fluorescence/`. Its manifest records source URLs,
citations, units, valid velocity intervals, retrieval snapshot, and SHA-256
hashes. Runtime values come from the immutable compiled representation: there
is no network or data-file access in a sampling loop.

The H model evaluates the irradiance-dependent IKL equations directly, so a
case-specific H velocity table is neither stored nor needed. The OH model
interpolates all 1,201 reviewed, non-extrapolated Lowell samples from -60 to
+60 km/s. Both models reject velocities outside their supported intervals.
Core tests reproduce the six archived scalar case values (within the stated
rounding tolerances), verify inverse-square scaling, and audit every source
hash. Run them with:

```sh
make -C src/models/resonance_fluorescence test
```

## Chemical model and g-factors in detail

Added 2026-11-26. The three observed emissions depend on two parts of the model:

- **Chemistry and transport.** These set how many OH molecules, H atoms and O atoms are where.
- **The g-factor.** This is the number of photons that one molecule or atom emits per second.

This section describes both as implemented. The sections above give the inputs, files and commands.
The numbers are for the six cases. Q(H2O) does not enter the rates or g-factors.

### Species

`input/i3.conf` compiles seven neutral species: H, H2, O, OH, H2O, CO and CO2.

- **Parents.** H2O, CO and CO2 are injected at the nucleus surface (175 K).
- **Daughters.** OH, H, H2 and O are produced only by the chemistry below.
- **CO and CO2** are collision partners only. They have no photochemistry.
- **Ions** are not tracked. Ionization, and dissociative ionization, remove the neutral; no neutral product
  is created.
- **No UV shielding.** The solar UV is not attenuated in the coma, so every rate below applies everywhere.

### Photolysis network and rates

Each reaction rate is computed at 1 AU from the per-channel solar-activity law of Moore and Mapaye (2025):

```text
k_1AU = A0 * F10.7P^alpha,   F10.7P = (daily adjusted F10.7 + centred 81-day mean) / 2
k(r_h) = k_1AU * (1 AU / r_h)^2
```

F10.7P is averaged over each case's Earth-equivalent interval (section "Solar-activity correction and
reaction rates"). It is 136.2, 142.9, 139.7, 143.6, 150.4 and 138.1 sfu for cases 1--6.

The generator writes `k_1AU` per case as `chemistry.rate.<ID>_1au_s-1`. The `HELIOCENTRIC` law applies
`(1 AU / r_h)^2` once. In the 1-D model, r_h is constant within a case, so each rate is a constant of the
case.

A0, alpha and the excess energies are in `reference/photochemistry_channels.csv`. In the table, "1-AU rate"
gives the range over the six cases. Within a channel the cases differ by 1--12%; the channels with larger
alpha differ more. "Product speeds" are the recoil speeds in the
parent's rest frame (explained after the table).

| ID | Reaction | A0 (s^-1) | alpha | 1-AU rate (s^-1) | Excess energy (eV) | Product speeds (km/s) |
|---|---|--:|--:|---|--:|---|
| `H2O_H_OH` | H2O -> H + OH | 3.60e-6 | 0.24 | 1.17--1.20e-5 | 4.04 | H 27.0, OH 1.60 |
| `H2O_H2_O` | H2O -> H2 + O | 1.80e-7 | 0.36 | 1.06--1.09e-6 | 3.94 | H2 18.3, O 2.31 |
| `H2O_2H_O` | H2O -> 2H + O | 2.40e-7 | 0.36 | 1.41--1.46e-6 | 0.697 | three-body, random split |
| `H2O_OHPLUS_H_ION`, `H2O_OPLUS_H2_ION`, `H2O_HPLUS_OH_ION`, `H2O_ION` | H2O -> sink | see CSV | see CSV | 8.1--8.7e-7 together | 0 | none |
| `OH_O3P_H` | OH -> O(3P) + H | 3.97e-6 | 0.059 | 5.31--5.34e-6 | 1.27 | H 15.1, O 0.95 |
| `OH_O1D_H` | OH -> O(1D) + H | 1.38e-6 | 0.35 | 7.71--7.98e-6 | 7.73 | H 37.3, O 2.35 |
| `OH_O1S_H` | OH -> O(1S) + H | 1.74e-7 | 0.36 | 1.02--1.06e-6 | 10.0 | H 42.4, O 2.67 |
| `OH_ION` | OH -> sink | 1.53e-8 | 0.68 | 4.32--4.63e-7 | 0 | none |
| `H2_2H_EXCITED` | H2 -> 2H | 1.03e-8 | 0.36 | 6.0--6.3e-8 | 0.488 | each H 6.8 |
| `H2_2H` | H2 -> 2H | 6.84e-9 | 0.50 | 8.0--8.4e-8 | 8.22 | each H 28.1 |
| `H2_ION`, `H2_DISS_ION` | H2 -> sink | see CSV | see CSV | 1.4--1.5e-7 together | 0 | none |
| `H_ION` | H -> sink | 1.20e-8 | 0.45 | 1.10--1.15e-7 | 0 | none |
| `O_ION` | O -> sink | 1.81e-8 | 0.67 | 4.87--5.21e-7 | 0 | none |

The energy sources are as follows:

- **H2O and H2 channels:** Huebner and Mukherjee (2015), as tabulated by Carberry Mogan et al. (2022).
- **OH channels:** Bryans and Pesnell (2012) for O(3P), and Bisikalo et al. (2015) for O(1D) and O(1S).
- **Ionization rates:** Moore and Mapaye (2025).

O(1D) and O(1S) are the electronic states at birth. The model has only one O species, so all of them become
ordinary O atoms: there is no O(1D) 630-nm or O(1S) 557.7-nm emission and no quenching. The channels keep
separate IDs, rates and excess energies.

At 1 AU, the lifetimes and branching are:

- **H2O:** total lifetime 6.5--6.7e4 s. The branches are H + OH 78%, 2H + O 9.4%, H2 + O 7.1% and ionization
  5.5%. At the comet this lifetime is 1.2e5 s (cases 4--6) to 2.0e5 s (case 1).
- **OH:** in the default (`STATE_RESOLVED`) OH loss, the lifetime is 6.7--6.9e4 s. The branches are
  O(1D) + H 53--54%, O(3P) + H 36--37%, O(1S) + H 7.1% and ionization 3.0%.
- **Other species:** H lives 8.7--9.1e6 s, O 1.9--2.1e6 s and H2 3.4--3.6e6 s.

**How a reaction is carried out** (`srcHartley2/CometChemistry.cpp`):

1. A parent particle reacts at the total rate of its species, using its species-local time step.
2. One channel is chosen with probability k_i / k_total.
3. All products of that event are created together. In I3, reactions occur stochastically in the ordinary
   way: no forced reactions and no channel biasing.
4. The product velocities are the parent velocity plus an isotropic recoil. The recoil conserves momentum,
   and its kinetic energy in the parent frame equals the channel's full excess energy.

All of the excess energy goes into translation; no internal excitation of OH or H2 is kept. The speeds in
the table follow from this:

- **Two-body channels:** each product's speed is fixed (momentum balance).
- **H2 -> 2H:** the two H atoms share the energy equally.
- **H2O -> 2H + O:** the energy is split at random among the three products.

The parent's own velocity is then added.

### OH loss variants used in the comparisons

The OH loss is the `OH_LOSS` chemistry group (section "OH loss variants" and
`srcHartley2/CHEMISTRY_VARIANTS.md`):

- **`STATE_RESOLVED`** uses the four OH channels above.
- **`SA1982_VH`** replaces them with one lumped channel, `OH_O_H_SA1982` (OH -> O + H, 0.36 eV), plus
  ionization `OH_ION_SA1982`:

```text
k_d(1 AU)   = 1/tau_FUV + 1/(S * OHTAUS(v_h)),   tau_FUV = 4.24e5 s,  S = 0.833e5 s
k_ion(1 AU) = 0.05 * k_d
both * (1 AU / r_h)^2, no solar-activity scaling
```

The terms in this formula are:

- `OHTAUS` is the heliocentric-velocity-dependent part of the OH photodissociation lifetime, from the legacy
  Fortran `BLOCKDATA OHPHOT` (which cites Schleicher and A'Hearn 1982, ApJ 258, 864). It is tabulated in
  `reference/oh_lifetime_sa1982.csv` (-60..+60 km/s, 2 km/s steps, linear interpolation, clamped outside).
- `tau_FUV` is the velocity-independent FUV/EUV part (mostly Lyman-alpha).
- `v_h` is the case's `brightness.heliocentric_radial_velocity_km_s`.

The OH lifetime uses this one bulk value per case. The g-factors below are evaluated per particle.

| Case | v_h (km/s) | OHTAUS (1e5 s) | k_d, 1 AU (s^-1) | k_ion, 1 AU (s^-1) | OH lifetime, 1 AU (s) | `STATE_RESOLVED` lifetime, 1 AU (s) |
|---|--:|--:|--:|--:|--:|--:|
| 1 | -39.7 | 3.755 | 5.56e-6 | 2.78e-7 | 1.71e5 | 6.91e4 |
| 2 | -36.2 | 3.320 | 5.97e-6 | 2.99e-7 | 1.59e5 | 6.83e4 |
| 3 | -29.6 | 2.860 | 6.56e-6 | 3.28e-7 | 1.45e5 | 6.87e4 |
| 4 | -5.8 | 2.035 | 8.26e-6 | 4.13e-7 | 1.15e5 | 6.82e4 |
| 5 | +0.4 | 2.760 | 6.71e-6 | 3.35e-7 | 1.42e5 | 6.74e4 |
| 6 | +7.1 | 3.420 | 5.87e-6 | 2.93e-7 | 1.62e5 | 6.89e4 |

`SA1982_VH` OH lives 1.7--2.5 times longer. Its H atoms from OH leave at 8.0 km/s and its O atoms at
0.51 km/s; the state-resolved channels give 15--42 km/s and 0.95--2.7 km/s. The H2O network is the same in
both variants.

The run log and every brightness header record:

- the active variant (`Chemistry.Variant.OH_LOSS`);
- the table rates (`Chemistry.TableRate.*`);
- the SHA-256 of the table file.

The comparison with MAVEN/IUVS in `validation_runs/i3_qfit_campaign_20261123` uses `SA1982_VH`. Its deck,
build input and scripts are in `campaigns/qfit_20261123/`. Q(H2O) is fitted to the IUVS OH: 1.85e28, 3.78e28, 8.89e28, 2.10e29, 1.53e29 and 1.48e29 s^-1 for cases 1--6. CO and CO2
are as in the generated `input/i3.in`.

### g-factors: what AMPS evaluates

All three emissions are treated as optically thin, isotropic (4 pi) resonance fluorescence. A molecule or
atom moving with heliocentric velocity `V` sees the solar line Doppler-shifted by

```text
v_r = (V_comet + v_particle) . e_r,     V_comet = rdot * e_r
```

where `rdot` is the case's `brightness.heliocentric_radial_velocity_km_s` and `e_r` points from the Sun to
the comet.

The brightness sampler (`srcHartley2/CometBrightness.cpp`) handles each line-of-sight point as follows:

1. It rotates the particle's DIM1 local velocity into the physical frame.
2. It evaluates the g-factor of that particle at its own `v_r`, which includes the Swings effect particle by
   particle.
3. It adds `column x g / 1e10` to the brightness in rayleigh.

No velocity is clamped or extrapolated. A particle outside a model's velocity range is either excluded and
reported, or stops the run, according to `brightness.velocity_domain_policy`.

The output tables `PT/plots/pic.<species>.brightness.s=<i>.out=<k>.dat` record, per impact parameter:

- the column;
- the brightness;
- the "Effective Photon g Factor" (brightness x 1e10 / column), i.e. the column-weighted mean g that was
  applied.

At start-up AMPS also evaluates every line at the bulk `rdot`, as a check of the inputs. All g-factors
scale as `(1 AU / r_h)^2`.

The models live in `src/models/resonance_fluorescence/`:

- its manifest records the sources, URLs, units, velocity ranges and SHA-256;
- its tests check the evaluators against the source data.

The case inputs are in `reference/i3_cases.csv` (OH, H) and `reference/solar_fuv_lines.csv` (O).

#### OH A-X (0-0) band, 308.5 nm (`OH_AX00_SA1988_LOWELL_UNQUENCHED`)

Source: the fluorescence efficiency of Schleicher and A'Hearn (1988, ApJ 331:1058), from the Lowell
Observatory comet g-factor service. The service gives the (0-0) band's energy g-factor at 1 AU
(erg s^-1 molecule^-1) as a function of heliocentric radial velocity. The data were retrieved from -60 to +60
km/s in 0.1 km/s steps; all 1,201 values are reviewed, interpolated linearly and never extrapolated.

- **Unquenched.** These are the collision-free values. Collisions in the inner coma change the OH
  ground-state populations, and with them the g-factor; that effect is not included.
- **Photon conversion.** Energy is converted to photons at the effective band wavelength of 308.5 nm, which
  each case selects explicitly to reproduce the archived conversion.
- **No solar-activity scaling.** The values are for one fixed solar spectrum.

The solar spectrum near 309 nm is full of absorption lines, so g depends strongly on velocity (Swings
effect). Within +-20 km/s it is lowest at -0.3 km/s, 1.49e-15 erg s^-1 at 1 AU. It reaches 2.68e-15 at
+6 km/s and 2.47e-15 at -10 km/s, and its largest value is 5.93e-15 at +15 km/s.

| Case | rdot (km/s) | energy g, 1 AU (erg s^-1) | photon g, 1 AU (s^-1) | photon g at r_h (s^-1) | applied / bulk at 1e5 km | at 1e6 km |
|---|--:|--:|--:|--:|--:|--:|
| 1 | -39.7 | 3.55e-15 | 5.52e-4 | 1.80e-4 | 0.99 | 0.99 |
| 2 | -36.2 | 3.33e-15 | 5.18e-4 | 1.88e-4 | 1.02 | 1.03 |
| 3 | -29.6 | 3.06e-15 | 4.75e-4 | 2.00e-4 | 0.98 | 0.94 |
| 4 | -5.8 | 1.71e-15 | 2.66e-4 | 1.44e-4 | 1.01 | 1.01 |
| 5 | +0.4 | 1.50e-15 | 2.33e-4 | 1.26e-4 | 1.09 | 1.21 |
| 6 | +7.1 | 2.70e-15 | 4.20e-4 | 2.24e-4 | 1.02 | 1.05 |

"Applied / bulk" is the Effective Photon g Factor of the comparison runs (`i3_qfit_campaign_20261123`,
output 4) divided by the g at `rdot`. In case 5 the comet sits at the minimum of the g curve, so the spread
of OH velocities raises the applied g-factor by 9% at 1e5 km and 21% at 1e6 km.

#### H Lyman-alpha, 121.6 nm (`H_LYA_IKL2020`)

Source: the solar Lyman-alpha line model of Kowalska-Leszczynska et al. (2020, ApJS 247:62), the authors'
Table 1 (`table1v4.txt`). The model gives the radiation-pressure factor `mu` (Lyman-alpha radiation pressure
relative to solar gravity) as an analytic function of radial velocity:

- a kappa function for the line;
- minus a Gaussian for its self-reversed core;
- plus a linear background.

Each of the nine coefficients depends linearly on the total solar Lyman-alpha flux. The model covers
-200..+200 km/s and is evaluated directly; no H table is stored.

The flux is the LASP/LISIRD composite Lyman-alpha, version 4 (Machol et al. 2019), averaged over the
case's Earth-equivalent interval: 5.10--5.45e11 photons cm^-2 s^-1 at 1 AU (`h_lya_photon_flux_1au` in
`reference/i3_cases.csv`).

`mu` is converted into a photon scattering rate by momentum balance. Each scattered photon transfers
`h nu / c` of radial momentum, so

```text
g(1 AU) = mu * (G M_sun m_H / AU^2) / (h nu / c),     g(r_h) = g(1 AU) * (1 AU / r_h)^2
```

| Case | rdot (km/s) | Ly-alpha flux, 1 AU (cm^-2 s^-1) | g, 1 AU (s^-1) | g at r_h (s^-1) | applied / bulk at 1e5 km | at 1e6 km |
|---|--:|--:|--:|--:|--:|--:|
| 1 | -39.7 | 5.10e11 | 2.96e-3 | 9.66e-4 | 0.96 | 0.95 |
| 2 | -36.2 | 5.28e11 | 2.95e-3 | 1.07e-3 | 0.97 | 0.95 |
| 3 | -29.6 | 5.44e11 | 2.85e-3 | 1.20e-3 | 0.97 | 0.95 |
| 4 | -5.8 | 5.24e11 | 2.42e-3 | 1.31e-3 | 1.01 | 1.01 |
| 5 | +0.4 | 5.45e11 | 2.54e-3 | 1.37e-3 | 1.01 | 1.02 |
| 6 | +7.1 | 5.45e11 | 2.55e-3 | 1.36e-3 | 1.03 | 1.03 |

At 1 AU the g-factor is highest in cases 1--3, despite their lower flux. At -30 to -40 km/s the comet sees
the brighter part of the solar line instead of its self-reversed centre.

The H output stays optically thin. A conservative line-centre cross section (9.5e-18 m^2: b = 8 km/s,
f = 0.4164) flags where tau > 0.1, but the result is not attenuated.

The radiation-pressure acceleration that defines `mu` is not applied to the H trajectories: in exact
spherical symmetry a Sun-directed force cannot be represented.

#### O I 130.4-nm triplet (`O_I_1304_KVB2022_DIRECT_J2` + `O_I_1304_SUMER2001_LYB_BOWEN_J2`)

The g-factor counts photons of the whole triplet (130.217 + 130.486 + 130.603 nm), which IUVS does not
resolve. It is given per O atom in the 2p4 3P2 ground level; the 3P1 and 3P0 levels relax radiatively within
3.1 h and 16 h. The comparison uses the sum of two terms (section "O I 130.4-nm brightness").

**Direct scattering at 130.217 nm (KVB2022).** The source is the data behind Fig. 14 of Killen, Vervack and
Burger (2022, ApJS 263:37):

- **The source table.** It gives the g-factors of the three lines from -50 to +50 km/s in 0.1 km/s steps, at
  0.352 au. They are computed from Gladstone (1992) two-Gaussian solar line shapes scaled to the SORCE
  SOLSTICE spectrum of 2009-03-01.
- **What we use.** We take the 130.217-nm (3P2 - 3S1) table and multiply it by 0.352^2 to get 1 AU. We use it
  within +-48 km/s, where the generated profiles end.
- **Scaling to the case.** The table is multiplied by the case triplet flux divided by the 8.41e9 photons
  cm^-2 s^-1 that the table implies.
- **The case triplet flux.** It is that reference times the FISM2 130.1--130.7-nm flux of the case's
  Earth-equivalent interval divided by the FISM2 flux of 2009-03-01. The ratio is 1.21--1.27 and cancels the
  FISM2/SOLSTICE calibration offset.

**Lyman-beta (Bowen) fluorescence.** Solar Lyman-beta pumps 3P2 to 3d 3D at 102.576 nm, and the 1128.7-nm and
844.6-nm cascade reaches 3s 3S1 and gives one triplet photon.

- **Line shape.** The solar line shape is from the SOHO/SUMER quiet-Sun disk-centre atlas (Curdt et al. 2001,
  from BASS2000). It is normalized to the case's FISM2 102.4--102.8-nm flux, 4.82--4.96e9 photons cm^-2 s^-1
  at 1 AU.
- **Atomic data.** Oscillator strengths and the branching to 1128.7 nm (0.287--0.289) are from NIST ASD.
- **Size.** The term is about 4.2e-7 s^-1 at 1 AU and almost independent of velocity.

| Case | rdot (km/s) | direct, 1 AU (s^-1) | Bowen, 1 AU (s^-1) | total at r_h (s^-1) | Bowen share | applied / bulk at 1e5 km |
|---|--:|--:|--:|--:|--:|--:|
| 1 | -39.7 | 3.30e-8 | 3.91e-7 | 1.38e-7 | 92% | 1.02 |
| 2 | -36.2 | 1.44e-7 | 4.09e-7 | 2.01e-7 | 74% | 1.05 |
| 3 | -29.6 | 1.29e-6 | 4.42e-7 | 7.30e-7 | 25% | 1.14 |
| 4 | -5.8 | 1.27e-5 | 4.18e-7 | 7.11e-6 | 3% | 0.97 |
| 5 | +0.4 | 1.04e-5 | 4.25e-7 | 5.85e-6 | 4% | 1.01 |
| 6 | +7.1 | 1.40e-5 | 4.17e-7 | 7.66e-6 | 3% | 1.02 |

**Velocity spread.** In cases 1--3, `rdot` lies on the steep wing of the solar line, so the velocity spread
of the O atoms raises the applied g-factor. In case 3 it is 37% higher at 2/3 of the domain radius.

**Line-shape uncertainty.** The direct term of cases 1--3 depends on the far wing of the assumed line shape.
The SUMER-atlas shape (`O_I_1304_SUMER2001_DIRECT_J2`, also written to the output) gives 53, 18 and 3.6 times
the KVB2022 value for cases 1, 2 and 3. Within |v| < 10 km/s the two shapes agree within 30%.

**Not included in the O model:**

- electron impact;
- prompt emission from photodissociation of H2O, CO or CO2;
- optical depth. The 130.217-nm line-centre tau exceeds 0.1 in all IUVS bins of case 5; it is flagged but
  not applied.

### Where the numbers come from

| Quantity | File | Written by |
|---|---|---|
| A0, alpha, excess energies | `reference/photochemistry_channels.csv` | reviewed transcription |
| F10.7P per case | `reference/i3_cases.csv`, `reference/solar_activity_audit.csv` | derivation recorded in this README |
| per-case rates `chemistry.rate.*` | `input/i3.in` | `prepare_i3_cases.py` |
| `SA1982_VH` lifetime table | `reference/oh_lifetime_sa1982.csv` | legacy `BLOCKDATA OHPHOT` |
| OH, H and O g-factor data and evaluators | `src/models/resonance_fluorescence/` | `tools/generate_reference_tables.py` |
| Ly-alpha flux and bulk OH/H g at r_h | `reference/i3_cases.csv` | `prepare_i3_cases.py` inputs |
| O triplet and Ly-beta fluxes | `reference/solar_fuv_lines.csv` | `derive_solar_fuv_lines.py` |
| applied g-factor | "Effective Photon g Factor" columns of the brightness tables | AMPS |

The "applied / bulk" ratios and the O table come from the comparison runs; that figure and its data are
in `validation_runs/i3_qfit_campaign_20261123/report_20261124/` (`plots/o1304_gfactor_out4.png`), which is
not in git.

## What does "corrected for the 3e9-m edge (collisionless reference)" mean?

### Short answer

AMPS deletes every particle that crosses the outer boundary at 3e9 m (3
million km). H atoms travel about 5e11 m before they are ionized, so every
line of sight misses the part of the H coma that lies beyond the boundary.

We estimate the missing part with a separate, much simpler model of the same
coma. It uses the same reaction rates, product energies, and AMPS gas flow,
but has no collisions after photodissociation. That model computes each
column twice: once for an unbounded coma, and once counting only atoms that
stayed inside 3e9 m, which is what AMPS can hold. The corrected curve is the
AMPS brightness multiplied by the ratio of the two.

- For H, the correction is +1–9% at the case-1 points and up to +30–73% at
  the outermost points of cases 4–6.
- OH is unchanged: it never reaches the boundary.

### Detailed answer

**What is missing, and why.** The outer boundary removes two kinds of atoms:

- atoms that would be beyond 3e9 m on the line of sight;
- atoms that would be inside, but whose parent (H2O, OH, or H2) left the
  domain before decaying, so in AMPS they are never born.

How much this matters depends on how far each species travels before it is
destroyed:

| Species | Speed | Lifetime at r_h = 1.36–1.75 au | Distance travelled | Compared with 3e9 m |
| --- | --- | --- | --- | --- |
| H | about 27 km/s | 1.6–2.8e7 s | 4–8e11 m | 150–250 times farther |
| OH | about 1.6 km/s | 1.3–2.1e5 s | about 3e8 m | well inside |

Inside the domain, the H density therefore falls almost exactly as 1/r^2.
For a 1/r^2 coma, the fraction of the column lost at impact parameter b is
1 - (2/pi) arctan(sqrt(R^2 - b^2)/b). That is 2% at b = 1e8 m, 9% at
4e8 m, and 38% at 1.7e9 m.

**The collisionless reference.** It is an independent steady-state Monte
Carlo model built from the run's own files:

- **Inputs:** case rates, nucleus radius, and `mesh.outer_radius_m` from
  `i3.in`; reaction topology and recoil energies from `i3.conf`; the AMPS H2O
  radial speed profile at the same output.
- **Parents:** water molecules leave the nucleus along the AMPS bulk flow
  and decay at sampled radii.
- **Products:** they receive the same recoil law as in AMPS (momentum
  conserved, full excess energy as translation). They fly in straight lines;
  OH -> O + H and H2 -> 2H are followed through the chain.
- **Steady-state density:** obtained by placing each H or OH at a
  log-uniformly sampled age with weight age * ln(1e9) * exp(-beta * age).
  The density is binned in 1,000 spherical shells and projected exactly onto
  each impact parameter.
- **Two columns from the same particles:**
  - N_full(b): every atom;
  - N_trunc(b): only atoms whose whole history (parent decay point, every
    intermediate decay point, and present position) stayed within 3e9 m.

  For straight-line motion this equals the population AMPS contains.

**The factor.** C(b) = N_full(b) / N_trunc(b), and corrected(b) =
AMPS(b) * C(b).

- Only the ratio enters, so the reference model's absolute normalization
  cancels.
- Because both columns share the same random particles, the noise in C is
  small: 0.14% between random seeds with 3e6 parents.
- C is shown only up to 2e9 m, two-thirds of the edge. Near the edge
  N_trunc goes to zero and C diverges.

**Size of the correction.** Median (range) of C at the H observation points,
campaign 2, latest outputs:

| Case | C for H | Model/obs before | Model/obs after |
| ---: | --- | ---: | ---: |
| 1 | 1.051 (1.014–1.094) | 0.55 | 0.60 |
| 2 | 1.046 (1.012–1.085) | 0.72 | 0.75 |
| 3 | 1.086 (1.017–1.175) | 0.94 | 1.02 (not converged) |
| 4 | 1.144 (1.036–1.296) | 1.01 | 1.16 |
| 5 | 1.193 (1.039–1.437) | 0.73 | 0.89 |
| 6 | 1.272 (1.045–1.727) | 0.91 | 1.19 |

C = 1.000 for OH in every case. Cases 4–6 have the largest corrections
because their observations reach 1.1–1.9e9 m.

**Why it can be trusted.**

- **AMPS against the truncated reference.** OH agrees within 2%. H agrees
  within 1–8%, and that excess is the collisional effect AMPS includes and
  the reference omits.
- **Analytic limit.** In the 1/r^2 limit the model reproduces the closed form
  within 0.6% (`srcHartley2/tests/test_comet_edge_correction.py`).
- **Independent implementation.** A separate implementation written for the
  campaign analysis agrees within 0.7%.

**What it is not.**

- It is not AMPS output: the raw AMPS curve is always plotted beside it.
- It does not model collisions.
- It ignores solar radiation pressure. The effect is about 1% at
  b <= 5e8 m and a few percent at the outermost points of cases 5–6.
- It assumes that the missing atoms have the same effective g-factor as the
  included atoms on each ray. That g varies by 1–5% across a case.
- It is a steady state.
- It does not fix differences between the model and the observations inside
  the domain. Case 1 is still 0.60 after correction.

**More detail and reproduction.**

- Method, I3 example, and the 2D/3D extension: `srcHartley2/EDGE_CORRECTION.md`.
- Complete description of the reference model and both implementations:
  `srcHartley2/COLLISIONLESS_REFERENCE.md`.
- Code: `srcHartley2/visualization/comet_edge_correction.py`.

Plot one case with the corrected curve:

```sh
python3 srcHartley2/visualization/plot_i3_brightness_comparison.py --no-show \
  --edge-correction --xmax 30 RUN/PT/plots/pic.H.brightness.s=0.out=K.dat --case N
```

Tabulate C(b):

```sh
python3 srcHartley2/visualization/comet_edge_correction.py \
  RUN/PT/plots/pic.H.brightness.s=0.out=K.dat --case N --csv OUT.csv
```

## Validation boundary

The focused generator test verifies all six shifted solar indices, reaction
count and reachability, minimum mixture paths, injection count, recoil-derived
characteristic speeds, brightness records, six named case sections, and
byte-for-byte regeneration of the tracked `input/i3.in`. The generic input
test independently verifies the overlay algorithm and fail-closed selector
rules. Case 1 also has a clean C++17/MPI build from `input/i3.conf`, a runtime
load through `-i input/i3.in -case 1`, a 200,000-sample-per-parent test of the
production crossing-flux injector, and a
20,000-event-per-neutral-reaction compiled chemistry test. All eight neutral
product channels passed mass, momentum, excess-energy, and ensemble checks;
all nine neutral-system ionization channels were recognized as sinks. The
case-1 AMR mesh built 376 blocks/3,008 dyadic cells and resolved the requested
surface width.

The pure brightness tests verify the inverse DIM1 velocity rotation, exact
shell chords, nucleus occultation, Rayleigh units, tau flagging, strict input
policy, velocity-domain policy, geometry, and exact source generation.  The
repaired 12-rank case-1 run has now produced finite, 200-ray H and OH Tecplot
tables at its first completed 10,000-iteration sampling boundary.  The
companion `visualization/plot_i3_brightness_comparison.py` overlays those
profiles on the reviewed numerical raster digitization in
`visualization/references/i3_brightness_MAVEN/`.  The digitization preserves
source-pixel anchors, hashes and limitations and can be rebuilt with
`visualization/extract_i3_brightness_maven.py`; ordinary plotting does not
read the archive. None of this evidence establishes convergence,
optically thin H validity, or agreement with the observations. Those require
later sampling windows, original quantitative uncertainties beyond what can
be resolved in the supplied raster, and review of every H tau flag.
