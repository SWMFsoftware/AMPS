# How to compile and run I3 for a given chemistry, case, and production or test mode

Written 2026-11-26 at the user's request. It covers three things:

- which inputs belong together;
- how to build production and test-mode executables;
- how to run one case with a given OH chemistry.

Background:

- [README.md](README.md): model, inputs, chemistry and g-factors, nightly test;
- [campaigns/qfit_20261123/README.md](campaigns/qfit_20261123/README.md): the campaign consistent with
  MAVEN/IUVS;
- [README_I3_RUN.md](README_I3_RUN.md): the general build-and-run guide;
- [../CHEMISTRY_VARIANTS.md](../CHEMISTRY_VARIANTS.md): the OH-loss variants.

Commands are given from the AMPS repository root unless stated otherwise. `C` stands for the campaign
folder:

```sh
C=srcHartley2/i3/campaigns/qfit_20261123
```

## 1. Choices

| Choice | When | How |
|---|---|---|
| Build input (`.conf`) | compile time | the `i3.conf` in `input/` of the build tree |
| Production or test mode | compile time | `Config.pl -amps-test=off` (production) or `-amps-test=on` (test) |
| Runtime deck (`.in`) | run time | `-i <deck>` |
| Case | run time | `-case N`, required, 1-6 |
| OH chemistry | run time | the deck's `chemistry.variant.OH_LOSS`, or `--chemistry-variant OH_LOSS=SA1982_VH` or `STATE_RESOLVED` on the command line |

One executable runs every case with both OH chemistries. Only a change of `.conf` or mode needs a rebuild.

For the OH chemistry, the precedence is:

1. the command line;
2. the case block of the deck;
3. the common block of the deck;
4. the compiled default, `STATE_RESOLVED`.

The OH chemistries are:

- **`STATE_RESOLVED`:** the four state-resolved OH channels with F10.7-scaled rates. This was the original
  default.
- **`SA1982_VH`:** one lumped OH -> O + H channel with the heliocentric-velocity-dependent lifetime of
  Schleicher and A'Hearn (1982). It was implemented later (2026-11-21) and is used by the campaign
  consistent with MAVEN/IUVS.

## 2. The three standard setups

| Setup | `.conf` | Mode | Deck | Deck's OH chemistry |
|---|---|---|---|---|
| MAVEN/IUVS campaign `qfit_20261123` | `$C/i3.conf` (`BlockCells=12`) | production | `$C/i3.in`: six cases, Q(H2O) fitted to OH, 4e9-m domain, O lines | `SA1982_VH` |
| Default deck | `input/i3.conf` (`BlockCells=8`) | production | `input/i3.in`: six cases, archive Q(H2O) x `reference/i3_rate_scaling.csv`, 3e9-m domain, O lines | `STATE_RESOLVED` |
| Nightly regression test `i3` | `input/i3.conf` | test | `input/i3-test.in`: cases 1 and 5 (the test runs case 1) | `STATE_RESOLVED`; `SA1982_VH` from the command line |

**The `.conf` files** differ only in `BlockCells` (12 or 8). Species, reactions, both OH-loss variants,
collisions and buffers are the same. The nightly references were made with 8; the campaign ran with 12.

**The decks** are written by the same generator (`prepare_i3_cases.py`), so chemistry rates, seeds, collision
settings and brightness geometry agree. The campaign deck differs in:

- Q(H2O);
- the domain radius;
- the surface-resolution floor;
- `EXCLUDE_AND_REPORT` instead of `ABORT`;
- its default OH chemistry.

`input/i3-test.in` uses 8 instead of 59 injection particles per step and has no O lines.

## 3. Compile

### 3a. Campaign build (separate directory, checkout untouched)

```sh
$C/scripts/compile_amps.sh --dest /path/build_prod --jobs 32                 # production, commit 848e9c5
$C/scripts/compile_amps.sh --dest /path/build_test --jobs 32 --test-mode     # test mode
$C/scripts/compile_amps.sh --dest /path/build_8053739 --commit 8053739        # another commit
```

The script does the following:

1. clones the repository at the commit into `--dest`;
2. clones SWMF `share` (at `c5f8beef`) and `util`;
3. copies `$C/i3.conf` to `input/i3.conf` of that tree;
4. runs `Config.pl -install -compiler=gfortran,gcc_mpicc` and
   `Config.pl -application=i3 -spice-path=nospice -spice-kernels=nospice -amps-test=off|on`;
5. runs `rm -rf -- build` and `make EXTRACOMPILEROPTIONS=-O3`;
6. checks `_BLOCK_CELLS_X_ 12` and `_MAX_REFINMENT_LEVEL_ 30`;
7. writes `build_provenance.txt`;
8. compares the build with the campaign build's fingerprints.

`--help` lists the options. On Voyager, run it inside a Slurm job (for example `--ntasks=32 --nice=1000`);
it took about 1 min with 32 jobs.

The campaign code is commit `848e9c5`. Master `8053739` gave byte-identical test-mode results with the
campaign deck (case 1, `SA1982_VH`, 2000 iterations). Later commits are untested; see the campaign README
for how to check one.

### 3b. Default build in your checkout (uses `input/i3.conf`)

```sh
./Config.pl -install -compiler=gfortran,gcc_mpicc          # once per fresh checkout (clones share and util if missing)
./Config.pl -application=i3 -spice-path=nospice -spice-kernels=nospice -amps-test=off   # -amps-test=on for test mode
rm -rf -- build                                            # always before a rebuild
make -j8 EXTRACOMPILEROPTIONS=-O3                          # writes ./amps
```

Do not copy the campaign `i3.conf` into `input/` of your checkout: the nightly test builds `input/i3.conf`
and its references need `BlockCells=8`. Use 3a for the campaign.

### Production or test mode

Changing between production and test mode means running `Config.pl` with the other `-amps-test` value,
then `rm -rf -- build` and `make` again.

- **Production (`-amps-test=off`).** Use it for every science run.
- **Test mode (`-amps-test=on`).** It turns on the debugger checks (particle lists, particle positions, block
  crossings) and sorts the particle lists after every MPI exchange, so a run is deterministic for a fixed
  rank count. Since commit `1379c68` it also forces the run length in the code: 2000 iterations, one
  2000-iteration sampling window and weight adaptation at 2000. The run log reports the replaced deck values.
  Test mode is therefore for regression checks and code comparisons only. (The test-mode build of `848e9c5`
  does not force the length; it uses the deck's values.)

## 4. Before running: check the inputs

`--validate-runtime-input` resolves the case and the OH chemistry and checks the whole input chain (comet,
chemistry, collisions, cooling, brightness, radial). It exits before the mesh is built, in seconds:

```sh
mpiexec -n 2 /path/amps -i <deck> -case N [--chemistry-variant OH_LOSS=...] --validate-runtime-input
```

All six cases with both OH chemistries of the campaign deck:

```sh
$C/scripts/validate_inputs.sh --amps /path/build_prod/amps --work /path/validate     # 12 checks
```

## 5. Run

Run each case in its own directory. AMPS writes `PT/`, the mesh files and the outputs into the current
directory.

**The `SA1982_VH` table.** `SA1982_VH` reads `srcHartley2/i3/reference/oh_lifetime_sa1982.csv`. The deck gives
its path relative to the deck's own directory:

- `input/i3.in` points to `../srcHartley2/i3/reference/oh_lifetime_sa1982.csv`, so it works where it is.
- The campaign deck expects `oh_lifetime_sa1982.csv` next to itself. If you run it in place, `SA1982_VH` stops
  at start-up with "cannot open the rate table file", while `STATE_RESOLVED` runs (checked 2026-11-26).
  Run the campaign deck through `run_case.sh`, which copies the deck, the table, `$C/i3.conf` and the
  executable into the run directory and writes `SHA256SUMS` and `provenance.txt`.

**Recorded variant.** Every run records the OH chemistry it used and where it was selected (deck or command
line) in its log (`AMPS: chemistry Chemistry.Variant.OH_LOSS=...; source=...`) and in the headers of the
brightness files.

### Campaign deck

```sh
# case 3, SA1982_VH, production, 128 ranks, runs until converged (1.7-3.1 days in the campaign)
$C/scripts/run_case.sh --amps /path/build_prod/amps --case 3 --out /path/qfit_rerun
#   -> /path/qfit_rerun/sa1982_vh/case3_np128/production

# the same case with STATE_RESOLVED (same deck and Q(H2O); passes --chemistry-variant OH_LOSS=STATE_RESOLVED)
$C/scripts/run_case.sh --amps /path/build_prod/amps --case 3 --out /path/qfit_rerun --chemistry state_resolved
#   -> /path/qfit_rerun/state_resolved/case3_np128/production

# short test: 2000 iterations, 4 ranks (max_iterations, sampling_length and
# weight_adaptation_interval set to 2000 in the run's copy of the deck)
MPIRUN="mpiexec --bind-to none" $C/scripts/run_case.sh --amps /path/build_test/amps --case 1 \
    --out /path/tests --np 4 --short 2000

# Slurm, all six cases (adapt partition and time limit)
for N in 1 2 3 4 5 6; do
  sbatch --partition=<partition> --time=7-00:00:00 --ntasks=128 --job-name=i3qfit-case$N \
    --export=ALL,CASE=$N,AMPS=/path/build_prod/amps,OUT=/path/qfit_rerun $C/scripts/production.sbatch
done
```

Options of `run_case.sh`:

- `--np NP`, default 128;
- `--setup-only`, to prepare the run directory without running;
- `$MPIRUN`, the launcher, default `mpiexec`.

The Q(H2O) of the campaign deck were fitted for `SA1982_VH`. With `STATE_RESOLVED`, whose OH lifetime is
1.7-2.5 times shorter, the deck does not fit OH.

Afterwards, make the figures and compare with the campaign's reference results:

```sh
python3 srcHartley2/i3/tools/make_qfit_figures.py --output final --campaign /path/qfit_rerun --report-dir /path/qfit_rerun/report
python3 $C/scripts/compare_results.py /path/qfit_rerun/report/data/qfit_numbers_final.json
```

### Default deck

The default deck uses the build of section 3b:

```sh
mkdir -p /path/run_case5 && cd /path/run_case5
# STATE_RESOLVED (the deck's default)
mpiexec -n 128 /path/AMPS/amps -i /path/AMPS/input/i3.in -case 5 > case5.run.log 2>&1
# SA1982_VH
mpiexec -n 128 /path/AMPS/amps -i /path/AMPS/input/i3.in -case 5 --chemistry-variant OH_LOSS=SA1982_VH > case5.run.log 2>&1
```

For a record that survives later changes to the checkout, copy the executable and the inputs into the run
directory first, as section 7 of [README_I3_RUN.md](README_I3_RUN.md) shows. That section also adjusts
the table path in the copied deck.

### Nightly regression test

From the AMPS root:

```sh
make test_i3                                  # 4 ranks
make test_i3 TESTMPIRUN4="mpirun -np 8"       # 8 ranks
```

The test does the following:

1. reconfigures the checkout in test mode and builds `input/i3.conf`;
2. moves `amps` into `run_test_i3/`;
3. runs case 1 of `input/i3-test.in` twice: `case1` with `STATE_RESOLVED`, and `case1sa` with
   `--chemistry-variant OH_LOSS=SA1982_VH`;
4. compares the 18 output tables with the references in SWMF_data.

References exist for passion only (tag `[5]`, 4 and 8 ranks). On other machines each result reports
"reference ... is missing".

After `make test_i3`, the checkout is configured in test mode and has no `./amps`. For a production build run
3b again with `-amps-test=off`.

## 6. Quick reference

| I want to run | Build | Command |
|---|---|---|
| the campaign, case N | `compile_amps.sh` (production) | `run_case.sh --amps B/amps --case N --out O` |
| the campaign deck, case N, `STATE_RESOLVED` | `compile_amps.sh` (production) | `run_case.sh ... --chemistry state_resolved` |
| the default deck, case N, `STATE_RESOLVED` | 3b, `-amps-test=off` | `mpiexec -n NP amps -i input/i3.in -case N` |
| the default deck, case N, `SA1982_VH` | 3b, `-amps-test=off` | `... --chemistry-variant OH_LOSS=SA1982_VH` |
| a deterministic short comparison | `compile_amps.sh --test-mode` or 3b with `-amps-test=on` | `run_case.sh ... --short 2000 --np 4` |
| the nightly test | (done by the test) | `make test_i3` |
| only check the inputs | any | `... --validate-runtime-input` or `validate_inputs.sh` |
