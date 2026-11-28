# Enceladus end-to-end Table tests: E3 INMS and UVIS 2010

Four `MakefileTest/Table` tests run the production Enceladus executable
from a fitted plume source distribution and compare the model with Cassini
observations:

| Table test | Tier | Deck | Steps / windows | Observation compared |
|---|---|---|---|---|
| `Enceladus--E3-INMS` | nightly (`Family=enceladus-e2e`) | `e3_inms.nightly.in` | 5000 / 5 | INMS E3 H2O density along the spacecraft track |
| `Enceladus--UVIS-2010` | nightly (`Family=enceladus-e2e`) | `uvis_2010.nightly.in` | 5000 / 5 | UVIS 2010 E10 solar occultation H2O column |
| `Long--Enceladus--E3-INMS` | long, `Ranks=8` (`Family=enceladus-long`) | `e3_inms.long.in` | 10000 / 10 | as above |
| `Long--Enceladus--UVIS-2010` | long, `Ranks=8` (`Family=enceladus-long`) | `uvis_2010.long.in` | 10000 / 10 | as above |

Long-tier names may not begin with a nightly test's name (the nightly
collects result files as `test_<name>*.diff`, and the legacy nightly test is
called `Enceladus`), hence the `Long--` prefix.

## Inputs

All inputs are SWMF_data files (`data/` in the AMPS root links to
`SWMF_data/PT/AMPS/data`), each directory with a `README.md` and
`SHA256SUMS`:

- `data/input/Enceladus/fitted_plume/`: the `tenishev2014_hybrid_8`
  (8 jets + 7 stripes) absolute H2O source distributions fitted to each
  observation alone in the 2026-11-11 standalone campaign
  (`independent_each_observation`): E3, 3.142482e28 molecules/s; UVIS 2010,
  4.576082e27 molecules/s.
- `data/input/Enceladus/observation_files/`: the processed E3 INMS H2O
  profile and the UVIS 2010 optical depth with its row-aligned ray table.
- `srcEnceladus/standalone/figure_data/e3_spice_ground_track.csv` (in AMPS):
  the dense E3 track of the flyby density product and refinement tube.

The decks copy the physics and numerics of the production cases
`E3--tenishev2014_hybrid_8` and `uvis_2010_e10_solar--tenishev2014_hybrid_8`
of `input/enceladus_cases.in`: H2O only, active-Sun Burger (2007) parent-loss
photolysis with discarded products, no plasma, 100 primaries per step,
300 m/s and Courant 0.5, seed 184467, no split/merge, the explicit runtime
mesh (box +/-1260.5 km; 600 km background; 300/150 km north/south surface;
300 km south-polar cone) and, for E3 only, the 150 km tube of <=20 km cells
around the E3 track. The two decks of one observation differ only in
`total_iterations`. Every path is relative to the deck.

The decks also set `reproducible_parallel_run=true`
(`srcEnceladus/RUNTIME_INPUT.md`). Without it, two 8-rank 5000-step UVIS runs
with the same seed differed on 2026-11-28 (1,132 vs 1,135 resident
macroparticles; columns up to a few percent apart): AMPS rebalances the domain
by measured wall-clock time, and exchanged particles enter the cell lists in
MPI message-arrival order. The key keeps the initial decomposition and sorts
the cell lists after every exchange, so the same seed and rank count
reproduce exactly, which the regression requires.

## Stages

- **Compile**: `./Config.pl -application=enceladus -spice-path=nospice
  -spice-kernels=nospice -amps-test=off` and `make amps`: the production
  executable, not the nightly test mode. The explicit mesh makes the mesh
  independent of the test mode anyway.
- **Run**: `${MPIRUN} ./amps -input …nightly.in` (nightly; 8 ranks on
  passion) or `$(TESTMPILAUNCH) 8 ./amps -input …long.in` (long), output in
  `run_test_<NAME>/run.log` and `run_test_<NAME>/PT/plots/`.
- **Check** (`check.sh`, gating; the `.diff` file must be empty):
  1. the run finished (`final_receipt` and `End of the run:` in `run.log`);
  2. it used the references' rank count (`PT/thread*`);
  3. the observation file has one zone per sampling window (5 or 10);
  4. the particle ledger closes (injected = expected, `ledger_residual=0`);
  5. its products equal the references under `DiffNum.pl`'s default
     tolerances (relative 1e-12):
     `output/test_Enceladus/<NAME>--observation.ref_np<N>.gz` (the
     observation Tecplot file without its `#` comment lines, which hold
     absolute input paths) and `…--ledger.ref_np<N>.gz` (final receipt,
     species ledger and global time step from `run.log`).
- **Compare** (`compare.sh`, never gating): plots the last window of the
  model against the observation into `run_test_<NAME>/comparison/`
  (`observation.<request>.comparison.png`, for E3 also
  `flyby_density.E3.png`) with `srcEnceladus/test/plot_observation_comparisons.py`
  and prints its metrics in the test log: INMS chi-square and median
  model/observation ratio; UVIS column-to-optical-depth shape fit and the
  model peak column relative to the published 9.0e19 m^-2. It needs a python3
  with matplotlib (`ENCELADUS_PYTHON`, else `python3`, else
  `/usr/bin/python3`); without one it prints why it skipped. The nightly
  keeps logs, not run directories, so the plots are in the local or retained
  run directory only.

With `reproducible_parallel_run=true`, repeated 8-rank runs of each deck are
bit-identical (verified 2026-11-28; `STATUS.md`), so the regression has no
Monte Carlo tolerance: a difference means a change of the code, the inputs,
the build or the rank count.

## Run locally

From the AMPS root, after the native preflight of `AGENTS.md` (this
reconfigures and rebuilds the shared build), with the rank count of the
references:

```sh
make test_Enceladus--E3-INMS MPIRUN="mpirun -np 8" TESTCOMPILEJOBS=16
make test_Enceladus--UVIS-2010 MPIRUN="mpirun -np 8" TESTCOMPILEJOBS=16
make test_Long--Enceladus--E3-INMS TESTMPILAUNCH="mpirun -np" TESTCOMPILEJOBS=16
make test_Long--Enceladus--UVIS-2010 TESTMPILAUNCH="mpirun -np" TESTCOMPILEJOBS=16
ls -l test_Enceladus--E3-INMS.diff     # empty = PASS
```

`Makefile.test` and `Makefile.test.long` are generated from
`MakefileTest/Table` by `./Config.pl -install` (`utility/TestScripts/table_parser.cpp`);
regenerate them after editing the Table. The rundir stage moves `./amps`
into the run directory, so rebuild the production executable afterwards.

## Regenerate the references

Only after a reviewed change of the model, the inputs or the decks, from a
run on the reference rank count, in the same change:

```sh
for L in observation ledger; do
  gzip -n -c run_test_<NAME>/products/$L.dat > output/test_Enceladus/<NAME>--$L.ref_np8.gz
done
```

The references live in SWMF_data (`PT/AMPS/data/output/test_Enceladus/`)
and must be committed there for the nightly, which tests a fresh SWMF_data
clone.

## Limitations

These are regression and comparison tests, not convergence or validation.
The fitted sources come from a steady, collisionless, gravity-free model, so
AMPS (with gravity and photolysis) is not expected to reproduce the fitted
profiles; the 5000-step run is not shown to be steady. The UVIS run has no
refinement along the rays (as in the production case), the UVIS optical depth
is a shape proxy with one published peak column for the absolute scale, and
308 of its 1,133 rays start at the tangent point
(`srcEnceladus/RUNTIME_INPUT.md`).
