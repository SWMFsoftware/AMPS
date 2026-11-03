# Enceladus forward-model status

Last updated: 2026-11-01 (America/Chicago)

## Current outcome — E7 density and E21 native-count observations

The runtime observation mechanism now distinguishes three contracts without
mixing units:

- response-qualified local number density (`comparison_species=H2O|CO2`);
- model-only geometry (`comparison_species=none`); and
- native instrument context (`comparison_species=none` plus the three
  `context_*` keys), which is displayed but never scored against density.

`data_file=input:PATH` resolves a repository-released table relative to the
runtime deck. Ordinary relative observation paths continue to resolve through
`observation_data_root`, so future external products require only a new named
`observation` block and no recompilation. Duplicate names, unknown keys,
missing columns, a mixed density/context contract, invalid uncertainties, and
nonfinite values fail before MPI. `RUNTIME_INPUT.md` contains the complete
syntax and units.

The checked-in `/home/vtenishe/Mars2/AMPS/input/enceladus.in` now requests 13
products:

- quantitative INMS H2O density for E3 and E5;
- quantitative INMS apparent-CO2 density for E3, E5, E7, E14, E17, and E18;
- model-only H2O/CO2 profiles on E13 geometry;
- E21 H2O and CO2 open-source neutral-beaming counts as separate native-count
  context products; and
- the 2007 zeta-Orionis and 2010 E10 UVIS occultation products.

### New observation releases and provenance

`standalone/export_forward_observations.py` reproducibly exports the added
tables from `/data/vtenishe/cassini_validation_data` and writes per-product
metadata plus `release_manifest.json`.

- E7: 43 mass-44 density rows from `Flybys2.mat`, matched count-for-count to
  official PDS Level-1A geometry. The optimized time alignment is 0.321 s and
  its maximum residual is 0.006 s. The point uncertainty propagates Poisson
  counting error through the local density/count factor; a separate 15%
  common calibration uncertainty is recorded but not folded into each point.
  CSV SHA-256:
  `1df70acf4549488aac187b0d783da01edd53c562be3181ec286b82cd0c461acb`.
- E21 H2O: 704 native count rows with interpolated supplied trajectory
  geometry. The 322 earlier count rows outside the supplied trajectory range
  are excluded rather than extrapolated. CSV SHA-256:
  `7b4507bf58b3015a23987815f68cc14b448cc4b57a3006f613ec97211ddbc292`.
- E21 CO2: 481 native count rows, all within the supplied trajectory range.
  CSV SHA-256:
  `7a8e7867bef97a16df6a98f969875b0fdf854334ac718a36531741847f44f2a5`.

The E21 values remain counts per integration with Poisson count uncertainty.
A quantitative AMPS comparison needs the model phase-space distribution
accumulated through the released velocity-, voltage-, and attitude-dependent
SIMION response. No scalar count-to-density factor was invented.

The exporter now forces LF-only records. The native validator caught the
initial CRLF final-header token before any run; the regenerated release passes
the strict C++ table reader on all three products.

### Output control

The new `volume_field_output=false` setting suppresses only the generic
full-volume active-species AMR dumps. Completed density buffers, window resets,
surface output, conservation ledgers, dense E3 flyby output, and all 13
observation operators remain active. This avoids approximately 14 GiB of
redundant H2O/CO2 field files across ten windows. A two-window native smoke run
proved that the collecting buffer is flushed between observation-only outputs;
memory also stayed constant at 2.889238 GiB throughout the long run.

### Native 10,000-iteration run

The current equal-per-source E3 distribution was transferred without
refitting. Its 15 vertically projected sources sum to
H2O=`2.060909950673409e28 s^-1` and apparent-CO2=`5.356e25 s^-1`.
The one-rank neutral run used 100 injected macroparticles per species per
iteration, the 26.260416666666668 s global timestep, H2O parent-loss-only
photolysis, no daughter particles, no plasma, no splitting/merging, 10,000
iterations, and ten 1,000-step sampling windows.

```text
cd /home/vtenishe/Mars2/AMPS/srcEnceladus
python3 test/run_e3_simulation.py \
  --executable ../amps \
  --input ../input/enceladus.in \
  --output test_output/srcEnceladus/e3_10000_extended_observations \
  --launcher direct --timeout-s 10800 \
  --expected-iterations 10000 --expected-sample-length 1000 \
  --expect-volume-fields off
```

The process completed in **2,521.6 s** and passed all **13/13** native gates.
Exactly 1,000,000 particles of each species were injected. Both ledgers close
with zero residual:

| Species | escaped | surface loss | photolytic loss | resident |
| --- | ---: | ---: | ---: | ---: |
| H2O | 990,189 | 3,234 | 368 | 6,209 |
| CO2 | 986,230 | 6,531 | 0 | 7,239 |

All 13 observation products contain ten zones. The dense E3 output contains
3,201 rows/zone and 186 in-domain rows/zone. Every in-domain E3 sample uses an
actual 15,756.25 m cell edge, below the requested 20 km maximum.

The complete 938 MiB evidence directory is:

```text
/home/vtenishe/Mars2/AMPS/srcEnceladus/test_output/srcEnceladus/e3_10000_extended_observations
```

Important retained artifacts and hashes are:

- `test_receipt.json`:
  `bef3784396a587b24bea65fa4ab50470c95a2c78f166479d413627c5dea5c2e3`;
- `enceladus.in`:
  `c48fa7c5ecd6b5e9cb312098a5e0bca3446c1aeb7bd730cc01121c084a7c423f`;
- `enceladus_e3_equal_h2o_co2_sources.in`:
  `eb30bcbcf5852b60ebdd072faf740e5d4f776f01561a8b80865dded11c70986d`;
- `forward_observation_release_manifest.json`:
  `0a15046d633381823f2e7519f1484c9e597193e665bcf03fc77cae0ff6377357`;
- `PT/plots/observation.<request>.dat` for each request;
- `PT/plots/comparisons/observation.<request>.comparison.png` for each
  request; and
- `PT/plots/comparisons/observation_comparison_summary.json`:
  `c7b1c21542f68dfbb34a8e4f1fd3b751474a4968ebfb3c466f859e4af2a43e32`.

### Comparison results

The plotter passed all **13/13** schema/zone/PNG gates. PASS means that a
consistent product was generated; it is not a claim that the unrefitted E3
source distribution agrees with every event. Final-window density metrics are:

| INMS product | in-domain rows | median model/observed | chi-square/row |
| --- | ---: | ---: | ---: |
| E3 H2O | 16 | 1.177 | 300.54 |
| E5 H2O | 13 | 0.499 | 89.97 |
| E3 apparent CO2 | 44 | 3.141 | 284.15 |
| E5 apparent CO2 | 44 | 0.390 | 226.83 |
| E7 apparent CO2 | 43 | 2.505 | 617.98 |
| E14 apparent CO2 | 371 | 1.665 | 721.00 |
| E17 apparent CO2 | 377 | 1.516 | 809.76 |
| E18 apparent CO2 | 366 | 1.266 | 588.71 |

E13 has 88 in-domain model-only rows. E21 H2O contains 704 count samples, of
which 530 lie in the current AMPS box; all 481 E21 CO2 count samples lie in the
box. E21 reports `display_only_no_density_or_count_response` and no goodness
of fit.

| UVIS event | samples | shape chi-square/sample | AMPS peak H2O column (m^-2) | published column (m^-2) | ratio |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2007 zeta Orionis | 97 | 9.213 | 3.366e20 | 1.400e20 | 2.404 |
| 2010 E10 solar | 1,133 | 7.557 | 3.625e20 | 9.000e19 | 4.028 |

The mismatch is real evidence that one equal-source E3 distribution is not an
all-event fit. No observation was used to tune this rerun.

### Build, smoke, and regression evidence

The required native preflight was followed from
`/home/vtenishe/Mars2/AMPS`: the working directory and inactive process state
were checked, only `build/` was removed, `./Config.pl
-application=enceladus` was used without editing `.amps.conf`, production
hooks were regenerated with `./ampsConfig.pl -input enceladus.input
-no-compile`, and the top-level `env MAKEFLAGS=-j16 make amps` build passed.

The first hook-generation attempt exposed an existing packaging behavior:
the configured application source copy recursively included the prior 18 GiB
untracked `test_output` tree and filled the partial build. The process was
stopped, only the failed `build/` was removed, `test_output` was moved to a
sibling hold path during source-tree generation, and it was restored intact
before compilation. The repeated clean workflow passed.

- checked-in native input validation: **PASS**;
- `python3 test/run_tests.py`: **PASS, 42/42**;
- `python3 test/run_native_startup.py --executable ../amps`: **PASS, 6/6**;
- Python syntax checks for exporter, runner, and plotter: **PASS**;
- 200-step/two-window native smoke run: **PASS, 13/13**;
- smoke comparison generation: **PASS, 13/13**;
- native 10,000-step integration: **PASS, 13/13**;
- ten-zone comparison generation: **PASS, 13/13**;
- `git diff --check`: **PASS** after the final documentation audit.

### Current scientific limitations

- Only the E3 trajectory has a dedicated 20 km refinement tube. This run does
  not establish mesh convergence for E5/E7/E13/E14/E17/E18/E21.
- The finite box includes only 16/79 E3 and 13/79 E5 qualified H2O rows.
- E21 has no quantitative model-count response yet; a detector convolution
  must use phase space, voltage, attitude, and SIMION acceptance.
- E13 has no qualified H2O/CO2 measurement and remains model-only.
- UVIS columns include only the ray segment inside the finite AMPS domain.
  The 2007 ray is reconstructed from the tangent track, and 308 outward 2010
  branches lack a retained full Cassini vector.
- Published inventories identify additional 2005, 2011, 2016, and 2017 UVIS
  events, but the current local release lacks the qualified reduced spectra
  and explicit ray tables needed to add them without inventing data.
- The one-rank, 100-particle/species/step run has not passed particle-count,
  source-rate, mesh, timestep, box-size, or multi-rank convergence studies.

## Previous refined E3-only outcome

The requested refined E3 run is complete: 10,000 iterations in ten 1,000-step
sampling windows. All 13 run-level checks pass. A trajectory-following AMR
refinement capsule reduces the Cartesian cell edge at every in-domain E3
spacecraft sample from 126.05 km to 15.75625 km, an eightfold linear and
512-fold volume refinement. At each AMPS output the completed H2O/CO2
number-density buffer was evaluated along the frozen 3,201-point Cassini E3
SPICE trajectory and appended as one Tecplot zone. Model-only and INMS-overlay
PNGs were regenerated from the ten refined zones.

The coarse-mesh, effectively empty plateau from -14 through +10 s is removed:
the refined profile changes on one-to-two-second scales and no longer contains
the old 25-second, approximately 15-order-of-magnitude cellwise plateau. The
geometry is now resolved substantially better, but the density estimate is not
particle-converged. With the user-requested 100 injected
macroparticles/species/iteration and no splitting, 106 of the 186 in-domain
final-window H2O samples and 99 of 186 CO2 samples are exactly zero. In
particular H2O is empty from -10 through +7 s in the final window and CO2 is
empty from -8 through +4 s. These gaps are finite-particle occupancy, not a
resolved physical zero. Increasing spatial resolution alone therefore fixes
the cut-cell artifact but exposes Monte Carlo sparsity.

This remains a one-rank neutral-transport and H2O parent-loss photolysis run,
not a qualified observational validation or a mesh/timestep/particle-converged
production calculation. The direct data overlay has no instrument-response
convolution. The flyby sampler uses the compiled cell-centered-constant AMPS
stencil; it does not interpolate across neighboring cells.

The qualified evidence is:

- receipt: `test_output/srcEnceladus/e3_first_run_qualified/test_receipt.json`;
- native log: `test_output/srcEnceladus/e3_first_run_qualified/run.log`;
- AMPS fields and surface samples:
  `test_output/srcEnceladus/e3_first_run_qualified/PT/plots/`;
- application-independent receipt:
  `test_output/srcEnceladus/test_receipt.json`;
- native startup receipt:
  `test_output/srcEnceladus/native_startup/test_receipt.json`.

The retained superseded coarse long-run evidence is:

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

The current refined evidence is retained at the following absolute directory:

```text
/home/vtenishe/Mars2/AMPS/srcEnceladus/test_output/srcEnceladus/e3_10000_refined
```

- receipt: `test_output/srcEnceladus/e3_10000_refined/test_receipt.json`;
- native log: `test_output/srcEnceladus/e3_10000_refined/run.log`;
- exact runtime deck and equal-source startup receipt:
  `test_output/srcEnceladus/e3_10000_refined/enceladus.in` and
  `test_output/srcEnceladus/e3_10000_refined/enceladus_e3_equal_h2o_co2_sources.in`;
- refined mesh: `test_output/srcEnceladus/e3_10000_refined/mesh.dat` and
  `test_output/srcEnceladus/e3_10000_refined/mesh.msh`;
- flyby Tecplot product:
  `test_output/srcEnceladus/e3_10000_refined/PT/plots/flyby_density.E3.dat`;
- model-only PNG:
  `test_output/srcEnceladus/e3_10000_refined/PT/plots/flyby_density.E3.png`;
- AMPS/INMS comparison PNG:
  `test_output/srcEnceladus/e3_10000_refined/PT/plots/flyby_density.E3.inms_comparison.png`.

The complete 15 GiB AMPS field-output directory remains at
`/tmp/enceladus_e3_10000_refined_20261101`; the durable 171 MiB evidence set
above intentionally retains the mesh, log, receipt, inputs, Tecplot profile,
and PNGs rather than duplicating ten 742 MB volume files per active species in
the nearly full repository filesystem.

## Refined closest-approach E3 run (2026-11-01)

### Refinement and output implementation

Five fail-closed runtime controls were added to `input/enceladus.in`:

```text
flyby_mesh_refinement=true
flyby_refinement_target_cell_size_m=20000
flyby_refinement_radius_m=150000
flyby_refinement_time_min_s=-120
flyby_refinement_time_max_s=100
```

The mesh-resolution callback loads the selected segment of the same frozen E3
trajectory used by the observation operator and applies the target resolution
inside a 150 km-radius polyline capsule. The AMPS mesh builder tests a cell
characteristic diagonal, so the callback converts the requested maximum
Cartesian edge to its corresponding diagonal. Dyadic AMR then realizes an
actual 15,756.25 m edge throughout the in-domain track. The legacy south-polar
source refinement remains active outside this superposed capsule.

The runtime parser requires all four dependent refinement values when the mode
is on, rejects them when it is off, and validates positive lengths plus an
increasing time interval. The flyby Tecplot schema now records
`cell_max_edge_m` for each row. The integration runner requires every in-domain
value to be positive and at most the configured 20 km limit.

The refined mesh contains 3,165 leaf blocks and reaches level 5. Its minimum
edge is 15,756.25 m, compared with 126,050 m at level 2 in the superseded
coarse run. The global rule
`dt=min(Courant*minimum_cell_width/characteristic_speed)` consequently gives
`dt=26.260416666666668 s` for Courant 0.5 and characteristic speed 300 m/s.
Each 1,000-step sampling window spans 26,260.4167 s (7.29456 h), and the full
run spans 262,604.1667 s (3.03940 days) of model time.

The compiled species union remains unchanged, but generic volume-field output
is now enabled only for runtime-simulated species. This prevents inactive
neutral, ion, electron, and dust fields from consuming hundreds of MB apiece
on the refined mesh. The legacy zero-dust `mesh.data.dat` snapshot is also
suppressed only when the forward runtime is installed. Active H2O/CO2 mover,
sampling, chemistry, storage, and surface behavior are unchanged. Small
generic surface diagnostic files for the compiled species union remain.

### Build and staged run evidence

An initial clean-build invocation was not accepted because the configuration
generator was still running when compilation was started. No binary from that
overlap was qualified. The complete required preflight was repeated from
`/home/vtenishe/Mars2/AMPS`: the root and absence of active build/test processes
were confirmed, only `build/` was removed, and the application was regenerated
and built with the established workflow:

```text
cd /home/vtenishe/Mars2/AMPS
rm -rf -- build
./Config.pl -application=enceladus
./ampsConfig.pl -input enceladus.input -no-compile
env MAKEFLAGS="-j16" make amps
```

The final build passed. `.amps.conf` was not edited directly. Only inherited
warnings in legacy dust, Enceladus, and filename-formatting code appeared; the
modified forward-model sources compiled without a reported warning. Binary:

- `/home/vtenishe/Mars2/AMPS/amps`;
- 31,656,784 bytes, timestamp `2026-11-01 01:42:30.966129988 -0600`;
- SHA-256
  `1f3110e4447d77742f7563c683da17f67034cc6f8396308bb2f265e972499a02`.

Application-independent tests pass **41/41**. They now include parsed
refinement target/radius references and 17 typed negative configuration checks.
Receipt SHA-256:
`6c756838f7e607964d960d47d46678b2ae9f41a149e2a004d6327792eefbaa1b`.
Native startup tests pass **6/6**; refined-startup receipt SHA-256:
`84434d49d142258017a55df77e8596892e45d974219911230044dc3704163bcd`.

The first 100-step refinement pilot demonstrated the correct 15,756.25 m edge,
3,165-block mesh, and 26.2604167 s timestep, but it was deliberately stopped
during output after generic AMPS began writing the full refined field for each
compiled inactive species. That was an output-volume issue, not a solver
failure. After restricting generic volume output to active H2O/CO2, the fresh
pilot at `/tmp/enceladus_e3_refinement_smoke_active_20261101` passed all
**13/13** gates in 246.2 s. It injected 10,000 particles of each species,
closed both ledgers exactly, wrote only active volume fields, and verified all
186 in-domain trajectory points at 15,756.25 m.

The full command was:

```text
cd /home/vtenishe/Mars2/AMPS/srcEnceladus
python3 test/run_e3_simulation.py \
  --executable /home/vtenishe/Mars2/AMPS/amps \
  --launcher mpiexec \
  --expected-iterations 10000 \
  --expected-sample-length 1000 \
  --maximum-flyby-cell-edge-m 20000 \
  --timeout-s 14400 \
  --output /tmp/enceladus_e3_10000_refined_20261101
```

Result: **PASS, 13 checks, 0 failures**, elapsed wall time 5,235.4 s. The
process returned zero through ordinary MPI finalization, wrote outputs 0--9
with exactly 3,201 flyby rows and 186 in-domain positions in every zone, wrote
only H2O/CO2 generic volume fields, produced no plasma file, and verified
15,756.25 m as both the minimum and maximum in-domain flyby cell edge.

| species | injected | photolytic loss | surface absorption | outer escape | resident | ledger residual |
|---|---:|---:|---:|---:|---:|---:|
| H2O | 1,000,000 | 368 | 3,234 | 990,189 | 6,209 | 0 |
| CO2 | 1,000,000 | 0 | 6,531 | 986,230 | 7,239 | 0 |

The final native particle-buffer total was 13,448, exactly the sum of the two
independently traversed species lists. CO2 remained transport-only. H2O alone
underwent the selected Burger (2007) parent-loss channels, with daughters
discarded as configured.

### Refined flyby result and closest approach

The final 9,001--10,000 window has the following raw maxima inside the AMR
overlap:

| species | maximum density, m^-3 | E3 reference time, s | positive/186 cells |
|---|---:|---:|---:|
| H2O | 1.6602874525068922e14 | +27 | 80 |
| CO2 | 6.148649007991746e11 | +33 | 87 |

The nearest 1 s trajectory sample is at `t=0 s`, altitude 52.51995134 km; both
species have zero sampled particles in its containing 15.75625 km cell in the
final window. Selected points demonstrate the difference from the old coarse
plateau:

| E3 time, s | altitude, km | H2O, m^-3 | CO2, m^-3 |
|---:|---:|---:|---:|
| -14 | 111.497 | 0 | 1.438280e10 |
| -12 | 96.321 | 9.685010e12 | 0 |
| -11 | 89.421 | 1.245216e13 | 0 |
| -10 | 83.021 | 0 | 1.078710e10 |
| 0 | 52.520 | 0 | 0 |
| +5 | 62.937 | 0 | 3.595701e9 |
| +7 | 71.653 | 0 | 7.191402e9 |
| +8 | 76.890 | 5.534292e12 | 1.438280e10 |
| +11 | 95.762 | 2.628788e13 | 4.314841e10 |

The old values of order `1e-4 m^-3`/`1e-6 m^-3` from -14 through +10 s no
longer occur. Instead the finer mesh distinguishes plume-bearing and empty
cells at roughly one-second trajectory cadence. The remaining exact zeros must
not be interpreted physically: reducing the cell edge by eight reduces cell
volume by 512, while the fixed-particle timestep reduction increases the
resident macroparticle population only about eightfold. The expected local
occupancy is therefore much lower. A statistically resolved closest-approach
profile requires authorization to increase injected particles, enable a
qualified splitting scheme, lengthen/ensemble the sampling, or implement a
separately validated spatial observation estimator. None was assumed here.

The refined INMS overlay was regenerated with:

```text
python3 test/plot_flyby_density.py \
  test_output/srcEnceladus/e3_10000_refined/PT/plots/flyby_density.E3.dat \
  --output test_output/srcEnceladus/e3_10000_refined/PT/plots/flyby_density.E3.inms_comparison.png \
  --inms-h2o /data/vtenishe/cassini_validation_data/data/processed/plume_fit_v1/e3_h2o_sticking_corrected.csv \
  --inms-co2 /data/vtenishe/cassini_validation_data/data/processed/co2_fit_v1/e3_co2_density.csv \
  --closest-approach-offset-s -0.5
```

Result: **PASS**, ten zones, 186 in-domain positions, and 123 observational
rows loaded. The H2O observations remain outbound-only by product design, and
the mass-44 observable remains labeled apparent CO2. No source distribution
was fitted during this run: both species still use the equal 1/15 startup
fixture with E3 total rates inherited from the standalone observational
studies.

Artifact identities:

| artifact | bytes | SHA-256 |
|---|---:|---|
| run log | 2,254,848 | `cd986a979d8c3d9a26e19fa466ee8ee47699836d98e53cfb4000690bcf63a442` |
| run receipt | 4,775 | `81518157ed1bb0ad139b00325682c5e0d9be51cd4f2017fb09dcf8cd64a6f89b` |
| refined `mesh.dat` | 166,871,507 | `82b49cfa8693d0258d84f534c1cfacfd48a9331b08dbab4dbb67b51e537ff423` |
| refined `mesh.msh` | 716,945 | `c54a48b70ba5b4481341f7da8b34d58ce249cfb880e9ecb5550564711017637b` |
| runtime deck | 4,628 | `0b0e7c61579370ce6a081ef6cb5c49255b536b6d56bccee8d5be23d1dd84caca` |
| equal-source receipt | 15,874 | `eb30bcbcf5852b60ebdd072faf740e5d4f776f01561a8b80865dded11c70986d` |
| flyby Tecplot | 7,749,445 | `4a0770c56f6a441116f5d240efa9766584afc451e3a0d3029a3b93a6ad8f18e5` |
| model PNG, 2160x1512 | 389,007 | `a4a5753640769b53b9d48caf09fd91ebfd4d964ebb8237eaad7d115ac861069a` |
| INMS PNG, 2160x1512 | 465,977 | `104547f6743b6f51c4202971028351810526767369cedea40e2172e796d4cabb` |

## Superseded coarse 10,000-iteration E3 run (2026-11-01)

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
- Application-independent tests now pass **41/41**; receipt SHA-256 is
  `6c756838f7e607964d960d47d46678b2ae9f41a149e2a004d6327792eefbaa1b`.
- Native startup tests pass **6/6**; refined-startup receipt SHA-256 is
  `84434d49d142258017a55df77e8596892e45d974219911230044dc3704163bcd`.

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
- size: 31,656,784 bytes;
- timestamp: 2026-11-01 01:42:30.966129988 -0600;
- SHA-256:
  `1f3110e4447d77742f7563c683da17f67034cc6f8396308bb2f265e972499a02`;
- run deck SHA-256:
  `0b0e7c61579370ce6a081ef6cb5c49255b536b6d56bccee8d5be23d1dd84caca`;
- `.amps.conf` SHA-256:
  `affef8bbc6d754159a42a4298faa26c43ae8de78c23a5cd89511c6d5ca0b029d`.

Application-independent suite:

```text
python3 test/run_tests.py
```

Result: **PASS, 41 checks, 0 failures**. It includes independent quadrature,
source normalization, parser/reference, reaction constants, timestep/weight,
projected-stripe, flyby-refinement configuration references, and 17 typed
negative checks. Receipt SHA-256:
`6c756838f7e607964d960d47d46678b2ae9f41a149e2a004d6327792eefbaa1b`.
The adjacent `reference_comparison.svg` visualizes calculated/reference error.

Native startup suite:

```text
python3 test/run_native_startup.py \
  --executable /home/vtenishe/Mars2/AMPS/amps \
  --output test_output/srcEnceladus/native_startup_refined_active_output
```

Result: **PASS, 6 checks, 0 failures**. Receipt SHA-256:
`84434d49d142258017a55df77e8596892e45d974219911230044dc3704163bcd`.

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

- The current result is one seed, one rank, one refined mesh, one timestep and
  one particle count over ten 1,000-step windows. It checks execution, repeated
  output, requested trajectory resolution and exact discrete bookkeeping; it
  is not a mesh/timestep/particle convergence study.
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
  remain open. The 15.75625 km cell-centered-constant flyby sampler resolves
  the old cut-cell plateau but remains too sparsely populated at 100 injected
  particles/species/iteration to interpret exact closest-approach zeros.
- Generic refined volume fields are now limited to runtime-active H2O and CO2.
  Small surface files for compiled inactive species remain generic diagnostics;
  they do not mean those species were simulated or produced.
- `PT/amps_error.log` contains the generic pre-allocation warning that no mesh
  blocks yet existed and that connectivity-list mode was adjusted. The same
  run subsequently allocated and advanced 3,165 blocks and returned process
  code zero; the warning is retained rather than suppressed.

No observational agreement, production readiness, or unique source recovery
is claimed by the current E3 run.
