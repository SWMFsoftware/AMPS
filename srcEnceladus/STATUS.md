# Enceladus forward-model status

Last updated: 2026-11-01 (America/Chicago)

## Current outcome

The requested longer E3 forward run is complete: 10,000 iterations in ten
1,000-step sampling windows. All eleven run-level checks pass. At each AMPS
output the completed H2O/CO2 number-density buffer was evaluated along the
frozen 3,201-point Cassini E3 SPICE trajectory, appended as one Tecplot zone,
and the ten zones were rendered to PNG. A separate comparison PNG now overlays
the released E3 INMS-derived H2O and mass-44 observations without changing the
retained model-only artifact.

This remains a one-rank neutral-transport and H2O parent-loss photolysis run,
not a qualified observational validation or a converged production
calculation. The direct data overlay has no instrument-response convolution.
The flyby profile is cellwise constant on the current coarse mesh. The
near-closest-approach depression is controlled by changes of the containing
126.05 km cell and must not be interpreted as a resolved physical depletion or
as evidence of agreement with INMS.

The qualified evidence is:

- receipt: `test_output/srcEnceladus/e3_first_run_qualified/test_receipt.json`;
- native log: `test_output/srcEnceladus/e3_first_run_qualified/run.log`;
- AMPS fields and surface samples:
  `test_output/srcEnceladus/e3_first_run_qualified/PT/plots/`;
- application-independent receipt:
  `test_output/srcEnceladus/test_receipt.json`;
- native startup receipt:
  `test_output/srcEnceladus/native_startup/test_receipt.json`.

The new long-run evidence is:

- receipt: `test_output/srcEnceladus/e3_10000_run/test_receipt.json`;
- native log: `test_output/srcEnceladus/e3_10000_run/run.log`;
- all AMPS fields and surface samples:
  `test_output/srcEnceladus/e3_10000_run/PT/plots/`;
- flyby Tecplot product:
  `test_output/srcEnceladus/e3_10000_run/PT/plots/flyby_density.E3.dat`;
- flyby PNG:
  `test_output/srcEnceladus/e3_10000_run/PT/plots/flyby_density.E3.png`;
- AMPS/INMS comparison PNG:
  `test_output/srcEnceladus/e3_10000_run/PT/plots/flyby_density.E3.inms_comparison.png`.

## 10,000-iteration E3 run (2026-11-01)

The checked-in `input/enceladus.in` now selects 10,000 total iterations,
1,000-step sampling, and the E3 flyby density observation operator. All other
physical/numerical controls are unchanged from the first qualified E3 run:
active Sun, neutral H2O+CO2, seven Burger (2007) H2O loss-only channels,
plasma off, 100 injected macroparticles/species/iteration, characteristic
speed 300 m/s, Courant 0.5, random seed 184467, and split/merge off.

Command:

```text
cd /home/vtenishe/Mars2/AMPS/srcEnceladus
python3 test/run_e3_simulation.py \
  --executable /home/vtenishe/Mars2/AMPS/amps \
  --launcher mpiexec \
  --expected-iterations 10000 \
  --expected-sample-length 1000 \
  --output test_output/srcEnceladus/e3_10000_run
```

Result: **PASS, 11 checks, 0 failures**, elapsed wall time 473.5 s. The
process returned zero through ordinary MPI finalization, produced no plasma
output, wrote outputs 0--9 with exactly 3,201 flyby rows per zone, and retained
exactly 186 in-domain positions in every zone.

| species | injected | photolytic loss | surface absorption | outer escape | resident | ledger residual |
|---|---:|---:|---:|---:|---:|---:|
| H2O | 1,000,000 | 425 | 3,738 | 994,987 | 850 | 0 |
| CO2 | 1,000,000 | 0 | 7,774 | 991,242 | 984 | 0 |

The particle-buffer resident total was 1,834, exactly the sum of the two
independently traversed species lists. The global timestep remained
210.08333333333334 s; each sampling window therefore spans 210,083.3333 s
(58.3565 h) of model time and the complete run spans 2,100,833.3333 s
(24.3152 days). This is model evolution time, not Cassini trajectory time;
every AMPS window is separately evaluated against the same E3 geometry.

### Flyby product and coverage

The trajectory contains E3 times -1600 through +1600 s at 1 s cadence. The
current Cartesian AMR root box covers only the 186 spacecraft locations from
-100 through +85 s. The Tecplot product retains all rows and marks coverage
with `in_domain`; outside-domain rows contain explicit zero placeholders and
are excluded from the PNG and physical analysis.

For the final 9,001--10,000 sampling window, the raw cell-sampled maxima
within that covered interval are:

| species | maximum density, m^-3 | E3 reference time, s |
|---|---:|---:|
| H2O | 2.3369410627079513e14 | +22 |
| CO2 | 9.020153492626351e11 | +22 |

These raw values are not convergence-qualified. The profile has abrupt cell
transitions and near-zero cell values around closest approach; with a
126.05 km minimum cell width and finite Monte Carlo population, the plot
diagnoses the present numerical resolution rather than validating an
instrument-resolved flyby profile. Mesh refinement, an expanded domain,
particle/seed studies, and an observational response operator remain required
before a quantitative physical validation against Cassini data.

### Closest-approach depression diagnosis (2026-11-01)

The large depression in both species is not a fitted physical feature. In the
final sampling window the decisive transitions are:

| E3 reference time, s | altitude, km | H2O, m^-3 | apparent CO2, m^-3 |
|---:|---:|---:|---:|
| -15 | 119.71 | 2.291543e12 | 6.966671e9 |
| -14 | 111.50 | 3.184311e-4 | 1.283412e-6 |
| 0 | 52.52 | 3.062876e-4 | 1.588486e-6 |
| +10 | 88.97 | 5.059817e-4 | 1.798881e-6 |
| +11 | 95.76 | 2.008343e13 | 4.651938e10 |

The effectively empty interval is -14 through +10 s. The entry at -14 s
coincides with the trajectory crossing the `x=-252.10 km` and `y=-252.10 km`
cell planes; the exit at +11 s coincides with crossing the `z=-252.10 km`
plane. The 16-order changes occur at cell boundaries, and the reported values
remain exactly constant for several consecutive spacecraft positions within a
cell.

The implementation comment in `ForwardModelAMPS.cpp` currently says "linear
center-node interpolation", but that statement is incorrect for this build.
`PIC::Mesh::GetCenterNodesInterpolationCoefficients` dispatches through the
configured coupler stencil, and `src/pic/picGlobal.dfn` selects
`_PIC_COUPLER__INTERPOLATION_MODE__CELL_CENTERED_CONSTANT_`. The selected
stencil assigns weight one to the containing cell. This code/comment and
documentation discrepancy was discovered after the retained run; no results
were regenerated or changed during the diagnosis.

The south-polar source geometry and solid-body occultation can produce a real
density reduction along this part of E3, but the current calculation does not
resolve its depth or timing. The temporary forward source file also contains
only the 15 surface components and omits the standalone H2O fit's constant and
radial background terms, so a shadowed coarse cell can become effectively
empty. An explicit spatial observation operator plus mesh/particle convergence
is required before assigning physical meaning to the depression.

### INMS comparison overlay (2026-11-01)

The model-only PNG and its run receipt remain unchanged. The new comparison
PNG was rendered from the same ten-zone Tecplot file with two locally released
E3 observational products:

| observable | released CSV | CSV SHA-256 | loaded | plotted in -100 to +85 s AMR overlap |
|---|---|---|---:|---:|
| sticking-corrected, INMS-derived H2O | `/data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1/e3_h2o_sticking_corrected.csv` | `5233ac2f3bf6677839089fae2ac03edf5a9cf95bb34407192c848051f719fd05` | 79 | 16 |
| INMS mass 44, apparent CO2 | `/data/vtenishe/cassini_validation_data/data/processed/co2_fit_v1/e3_co2_density.csv` | `6c9ea9fd234d5ac5a169bf650bc36e1e6aa7f60f094894fdd505811a58955a8a` | 44 | 44 |

The corresponding metadata hashes are
`9f4bf56f92be2550e32852595f5ae24f0c69a7781eae1f9833cdca37ab0457aa`
for H2O and
`702dd6368ed03368fb9682270becb55832178c58389f3c966949aa25c7205e4c`
for apparent CO2. The 79-point H2O release spans +10 to +400 s from closest
approach; 16 points lie inside the present model coverage. All 44 CO2 rows lie
inside it: 37 have positive counts and seven are zero-count samples.

The absence of inbound H2O points in this PNG is therefore a data-product
selection, not an INMS non-detection. The selected fit-authoritative H2O
release is deliberately outbound-only because the raw mass-18 counts require
the time-dependent inlet/surface-sticking response; a scalar count-to-density
conversion was not substituted. By contrast, the selected closed-source mass
44 product covers approximately +/-50 s and consequently includes inbound
points. A separate digitized full-flyby published H2O trace is retained in
`standalone/figure_data/figure3_published_curves.csv` over approximately
-400 to +395 s. It is not merged with the sticking-corrected fit table because
it is a digitized published curve with different provenance and uncertainty.

The observational tables use time from closest approach. The plot applies
`t_reference = t_CA - 0.5 s`, because the independent 0.1 s SPICE search
resolved E3 closest approach at -0.5 s relative to reference UTC
`2008-03-12T19:06:12`. Thus the plotted H2O observations occupy +9.5 through
+84.5 s, and the mass-44 data occupy -49.975 through +48.743 s. This offset is
an explicit command-line argument and is not inferred from the data.

H2O error bars show the released independent digitization uncertainty; the
23% common density-scale uncertainty is stated below the plot but is not added
independently to every point. E3 mass 44 is labeled **apparent CO2**, not
composition-resolved CO2. Its error bars show Poisson counting uncertainty;
the 15% common calibration uncertainty and published factor-2-to-5 impact
fragmentation enhancement range are stated but not converted into pointwise
Gaussian bars. On the logarithmic plot, seven zero-count samples appear as
open downward triangles at the one-count sensitivity. Lower Poisson bars that
terminate at zero are omitted rather than extended to an arbitrary positive
density.

Exact rendering command:

```text
cd /home/vtenishe/Mars2/AMPS/srcEnceladus
python3 test/plot_flyby_density.py \
  test_output/srcEnceladus/e3_10000_run/PT/plots/flyby_density.E3.dat \
  --output test_output/srcEnceladus/e3_10000_run/PT/plots/flyby_density.E3.inms_comparison.png \
  --inms-h2o /data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1/e3_h2o_sticking_corrected.csv \
  --inms-co2 /data/vtenishe/cassini_validation_data/data/processed/co2_fit_v1/e3_co2_density.csv \
  --closest-approach-offset-s -0.5
```

Result: **PASS**, ten model zones, 186 in-domain trajectory points, and both
observation series parsed and plotted. The comparison PNG is 2160x1512 RGBA,
292,394 bytes, SHA-256
`65702676de31f1d222886c4c22cf30de344e4c94d9e6bcb5551a0bdc7f995e43`.
The plotting utility is SHA-256
`fb9d0761e10abfc513f5105a94ed0526e262e16b8c87a6c84340232ef1c1db79`.
Python bytecode validation, the observational rendering, and a model-only
backward-compatibility rendering all passed. A negative CLI check omitted the
required closest-approach offset and correctly exited 2 before creating an
output. No AMPS rebuild or native rerun was needed because this change reads
the retained Tecplot and observational products without altering the
simulation.

This overlay is not an independent fit result. The forward file assigns the
same 1/15 relative strength to every one of the eight jets and seven stripe
components for both H2O and apparent CO2. Only the species totals came from
the E3 standalone studies: H2O from the independent E3 fit and apparent CO2
from the all-five-flyby shared-shape study. The strongly unequal fitted H2O
component fractions and the fitted CO2 shared shape were not transferred to
the forward run. Apparent visual normalization agreement, particularly for
outbound H2O, is therefore partly expected from reuse of the same observation
family and does not validate the equal-strength spatial source fixture.

Artifact identities and sizes:

| artifact | size | SHA-256 |
|---|---:|---|
| Tecplot flyby density | 6,981,187 bytes | `22199dfd18438e37442177737dcdd559670614490c59a935ea9ce01a6d50d876` |
| PNG (1800x1368 RGBA) | 183,956 bytes | `031959cba97ea9756991e2dde5a462f4a46fa99e6c7d3e9bf9a41d4ffa405d93` |
| INMS comparison PNG (2160x1512 RGBA) | 292,394 bytes | `65702676de31f1d222886c4c22cf30de344e4c94d9e6bcb5551a0bdc7f995e43` |
| native log | 2,256,569 bytes | `e8347db0fdd8c62fc329e9f865aa1c66e6fac89794f8f197ec9a2032481e694d` |
| run receipt | 4,488 bytes | `e770d3562ef917378f8e6227f5de6044a6c6a12625351b7428615bda4862a023` |

The complete run directory is 2.6 GiB because generic AMPS also writes field
and surface files for the compiled species union.

### SPICE provenance and independent reproduction

The suitable E3 kernels already present under `~/SPICE/Kernels` are:

- LSK `NAIF/naif0012.tls`, SHA-256
  `678e32bdb5a744117a467cd9601cd6b373f0e9bc9bbde1371d5eee39600a039b`;
- event PCK `CASSINI/cpck10Mar2008.tpc`, SHA-256
  `e1a4e70401b0c25abd1f9d7efa038c162596083f51c5fba8d736521b2169094d`;
- reconstructed SPK `CASSINI/080228AP_SCPSE_08050_08076.bsp`, SHA-256
  `b3eee97309e6ae59e0e00292c9d60aeb5b90c2040ba982b01157378475b1859e`.

The local CSPICE N0067 static library hash is
`53d74623fcfd22dcd5f9d01d1ebc4889f7b1034a5bd6738a681c1361b4d60ee1`.
`standalone/tools/spice_ground_tracks.cpp` was rebuilt and run with observer
CASSINI, target ENCELADUS, frame IAU_ENCELADUS, aberration `NONE`, reference
UTC `2008-03-12T19:06:12`, and -1600:1:+1600 s. The independently regenerated
CSV and checked-in trajectory both hash to
`e5707e456fb03586a36958524007394c5238b5e48fbfa2db3bf41030ace5e663`.
AMPS itself remains built without SPICE; it consumes this frozen, independently
reproducible trajectory and therefore does not silently select kernels at run
time.

### Implementation and verification changes

- Runtime keys `flyby_density_output`, `flyby_id`, and
  `flyby_trajectory_file` select the observation operator. Dependent controls,
  unsafe IDs, malformed headers/rows, mixed flyby IDs, non-finite values,
  negative altitude, and non-monotonic times fail before MPI.
- The trajectory file is parsed transactionally during runtime installation.
  At each ordinary AMPS output a registered external sampling callback reads
  the completed number-density buffer, uses the build-selected native
  cell-centered constant stencil, enforces one AMR owner per point, reduces
  across MPI, and appends a Tecplot zone with SI values and coverage flags.
- `test/plot_flyby_density.py` strictly parses the zone schema and renders only
  in-domain values. Its optional INMS readers validate E3 identity, finite and
  strictly increasing values, the mass-44/apparent-CO2 interpretation, and
  count/density consistency before plotting. `test/run_e3_simulation.py`
  checks requested iteration, sample, and particle totals plus all Tecplot
  zones and the model-only PNG artifact.
- Application-independent tests remain **PASS, 39/39**; receipt SHA-256 is
  `b55599e7604d7f7da8ea2b9defa94a7ad3f3770018c2b9008e609ee6c11078b5`.
- Native startup tests remain **PASS, 6/6**; receipt SHA-256 is
  `c5a9d4051198bb38e5693461ab2b67c5550639c2d57d4deb7b488f7b9ed1b1a7`.

The first clean rebuild exposed and stopped on an incorrect namespace for the
generic interpolation-stencil capacity. After correcting that compile-time
API name, the full mandatory preflight was repeated from the repository root:
process check, `rm -rf -- build`, `./Config.pl -application=enceladus`,
`./ampsConfig.pl -input enceladus.input -no-compile`, and
`env MAKEFLAGS="-j16" make amps`. The second clean build passed. `.amps.conf`
was never edited directly. The resulting executable is 31,591,216 bytes,
SHA-256 `5b7643cebb05129058cec10906373c00bd493309e1bc9cb01841a5f6cd8cafc3`;
the long-run deck hash is
`733f6fbaa0e307e15c19453365fe4ef37263db958486d75c137afe431c7da0bc`.

## Qualified E3 run configuration

The following section records the retained, earlier 800-step qualification
run and should not be confused with the current 10,000-step result.

The run deck used for that retained qualification was complete rather than a
fail-closed template. It used:

| control | value |
|---|---|
| source case | `tenishev2014_hybrid_8`, E3 absolute receipt |
| primary species | neutral H2O and neutral CO2 |
| chemistry | seven Burger (2007) H2O photo channels |
| reaction products | discarded; parent loss only |
| solar activity | active |
| plasma | `none`; allocated coupled-background arrays explicitly zeroed |
| characteristic speed | 300 m/s for H2O and CO2 |
| Courant number | 0.5 for H2O and CO2 |
| primary particles | exactly 100/species/iteration globally |
| sample length | 800 steps |
| total length | 800 steps, one sampling window |
| random seed | 184467, folded with MPI rank for native RNG |
| split/merge | off |

After AMR allocation, the global rule
`dt=min(Courant*minimum_cell_width/characteristic_speed)` produced
`dt=210.08333333333334 s`, limited by H2O (the two configured limits are
equal). The implied minimum allocated cell width is 126,050 m. The 800-step
window represents 168,066.6667 s (46.6852 h).

The fixed primary statistical weights printed by the executable are:

| species | absolute receipt rate, molecules/s | macroparticle weight |
|---|---:|---:|
| H2O | 2.0609099506734099e28 | 4.3296283213730557e28 |
| apparent CO2 | 5.356000000000001e25 | 1.1252063333333338e26 |

The CO2 rate is the fitted mass-44 apparent rate. The retained factor-2--5
impact-fragmentation interpretation brackets physical E3 CO2 at
`1.0712e25--2.678e25 molecules/s`; no value inside that interval was selected.
Consequently this startup run must not be cited as a physical-CO2 inference.

## Run results and conservation

Command:

```text
cd /home/vtenishe/Mars2/AMPS/srcEnceladus
python3 test/run_e3_simulation.py \
  --executable /home/vtenishe/Mars2/AMPS/amps \
  --launcher mpiexec \
  --output test_output/srcEnceladus/e3_first_run_qualified
```

Result: **PASS, 9 checks, 0 failures**, elapsed wall time 40.0 s. The native
loop completed all 800 steps, every rank reported 800 injection calls, and the
executable terminated through the ordinary MPI-finalization path.

| species | injected | photolytic loss | surface absorption | outer escape | resident | ledger residual |
|---|---:|---:|---:|---:|---:|---:|
| H2O | 80,000 | 28 | 288 | 78,896 | 788 | 0 |
| CO2 | 80,000 | 0 | 568 | 78,513 | 919 | 0 |

For each row, `injected = photolytic + surface + escape + resident` exactly.
CO2 is transport-only and recorded zero photolytic losses. No
`ices.data.dat` was produced. The combined resident count of 1,707 agrees
with the particle-buffer count and the sum of per-species linked-list counts.

The run-level receipt SHA-256 is
`4f9570d54c460a4ad25b0b5290e03ae79c9a192fb31e36e31c218d133adb53eb`;
the native log SHA-256 is
`690c4aad886fe06eecb1bfd90a5a941bcf2346751e43a57cc824fcb6d718dc3d`.

## Implemented forward-model behavior

### Source handoff

- `tenishev2014_hybrid_8` is the only runtime source distribution currently
  implemented. The schema is extensible through the source-case identifier.
- The canonical source set has eight point jets and seven tiger-stripe
  polylines. Every published stripe direction is projected radially to the
  252.100 km Enceladus sphere. All 102 generated stripe nodes pass a `1e-6 m`
  compiled-surface check.
- The E3 and E5 startup receipts distribute each species total equally over
  the 15 sources only to make the forward model runnable. This is not a claim
  that the physical sources have equal activity. The standalone fit will
  replace these fixtures with a shape receipt or individually fitted rates.
- E3 totals are H2O `2.060909950673409e28` and apparent CO2 `5.356e25`
  molecules/s. E5 totals are H2O `4.520015885781454e28` and apparent CO2
  `1.013e26` molecules/s.
- `SOURCE_DISTRIBUTION_FILE.md` is the normative standalone/forward handoff
  syntax for both absolute component rates and normalized shape files.

### Runtime parser and physics adapter

- `input/enceladus.input` remains the infrequently changed build-generator
  deck. `input/enceladus.in` is parsed for every run. Unknown, duplicate,
  inactive, missing, or inconsistent controls fail before the simulation.
- `--validate-input` validates the complete per-run/source pair without
  initializing MPI or allocating the AMR mesh.
- Runtime controls select simulated/injected/chemistry species, reactions,
  solar activity, plasma authority, particle counts, speeds, Courant numbers,
  sampling and total lengths, seed, and generic split/merge mode.
- `plasma_background_mode=none` is a typed mode. The adapter zeroes every
  allocated DATAFILE background variable because the compiled Boris mover
  still interpolates E/B even when the active particles are neutral.
- The production injector creates the exact requested global count on each
  step. A fixed species weight is `Q*dt/N`; event time is uniform within the
  step; deterministic event replay and native AMR ownership choose one owner.
- Source velocities use the production-normalized Tenishev flux law. Gravity,
  the internal spherical boundary, and outer-domain deletion remain active.
- H2O loss uses the exact per-step constant-rate probability
  `1-exp[-dt*sum(nu_r)]`. No daughter, charge, recoil, or plasma feedback is
  created. CO2 has no invented Burger-2007 reaction and is transport-only.
- The final receipt counts injection, photolysis, surface absorption, true
  outer-domain escape, and resident particles separately. A thread-local
  synchronous marker distinguishes internal-sphere deletion from Boris's
  shared `_PARTICLE_LEFT_THE_DOMAIN_` return; per-species conservation is
  then checked against independently traversed native particle lists.

## Build and validation evidence

The final rebuild followed the required repository workflow. `.amps.conf` was
not edited directly:

```text
cd /home/vtenishe/Mars2/AMPS
# confirmed this root and no active build/test process
rm -rf -- build
./Config.pl -application=enceladus
./ampsConfig.pl -input enceladus.input -no-compile
env MAKEFLAGS="-j16" make amps
```

Build result: **PASS**. Only inherited warnings in legacy dust, legacy
Enceladus, and legacy filename-formatting code appeared; `ForwardModel.cpp`
and `ForwardModelAMPS.cpp` compiled without a reported warning or error.

- executable: `/home/vtenishe/Mars2/AMPS/amps`;
- size: 31,522,360 bytes;
- timestamp: 2026-11-01 01:16:57.464491013 -0500;
- SHA-256:
  `ac398e461ab52576915cf67b8028869cc1800518a88ae1df4b85f779baba3de3`;
- run deck SHA-256:
  `cef51eb2283cb616495fca83acc3b0577e8b7b4b78c7cc3d62439b8a854e2bb9`;
- `.amps.conf` SHA-256:
  `affef8bbc6d754159a42a4298faa26c43ae8de78c23a5cd89511c6d5ca0b029d`.

Application-independent suite:

```text
python3 test/run_tests.py
```

Result: **PASS, 39 checks, 0 failures**. It includes independent quadrature,
source normalization, parser/reference, reaction constants, timestep/weight,
projected-stripe, and 15 typed negative checks. Receipt SHA-256:
`12fd71b7c07692f90541ccc3b683eb6cf18b89fdaab66ca2004390a63e638e8e`.
The adjacent `reference_comparison.svg` visualizes calculated/reference error.

Native startup suite:

```text
python3 test/run_native_startup.py \
  --executable /home/vtenishe/Mars2/AMPS/amps
```

Result: **PASS, 6 checks, 0 failures**. Receipt SHA-256:
`a88b36a5e74268c956c401c0ddebac961f730a8f4076d4e51707b1bb9644163c`.

### Retained superseded attempts

Two first attempts are retained under `e3_first_run/` and
`e3_first_run_singleton/`. Both failed before `MPI_Init` because the restricted
execution sandbox exposed no network interface for Open MPI. They are
environmental failures, not model runs. The successful command was therefore
run with host-level local-MPI socket access.

An earlier successful native run is retained under `e3_first_run_native/`.
Its dynamics and total particle count completed, but its diagnostic wrapper
classified sphere absorptions both as surface losses and domain escapes
because Boris shares one return code for those deletions. That receipt is
superseded by `e3_first_run_qualified/`, whose separated sink ledger closes
exactly.

## Limits and next validation gates

- The current result is one seed, one rank, one coarse mesh, one timestep and
  one particle count over ten 1,000-step windows. The longer duration checks
  execution, repeated output and exact discrete bookkeeping; it is not a
  mesh/timestep/particle convergence study.
- Fresh particles receive chemistry beginning with the following full step;
  their initial partial-birth-step exposure is omitted. This operator-split
  error is `O(nu*dt)` and requires a timestep study or an exact birth-to-end
  loss treatment before chemistry accuracy is qualified.
- Only H2O parent decay is modeled. No daughters, CO2 photolysis, charged
  species, charge exchange, H3O+, plasma, or back-coupling is active.
- The apparent CO2 rate and equal 15-way fixture are not final physical source
  products. A fitting-model receipt with justified CO2 interpretation remains
  required.
- Four-rank ownership/rank-equivalence, multiple seeds, expanded domain,
  mesh/timestep/particle convergence, an explicit spatial observation
  operator, analytic profile comparison, and observation-operator comparison
  remain open. The current cell-centered constant flyby sampler is not
  adequate for interpreting the closest-approach structure.
- The many zero-valued field files for compiled-but-runtime-inactive species
  are generic AMPS output behavior; they do not mean those species were
  simulated or produced.
- `PT/amps_error.log` contains the generic pre-allocation warning that no mesh
  blocks yet existed and that connectivity-list mode was adjusted. The same
  run subsequently allocated and advanced 64 blocks and returned process code
  zero; the warning is retained rather than suppressed.

No observational agreement, production readiness, or unique source recovery
is claimed by the current E3 run.
