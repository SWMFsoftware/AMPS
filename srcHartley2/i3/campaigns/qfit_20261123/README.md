# I3 campaign `qfit_20261123`: SA1982_VH OH loss, Q(H2O) fitted to MAVEN/IUVS OH

This folder holds everything needed to rebuild and rerun, from git on another machine, the six-case I3
(3I/ATLAS) campaign that is consistent with the MAVEN/IUVS observations. Only the run output stays outside git.

The campaign uses the OH chemistry that was implemented later (`OH_LOSS = SA1982_VH`, 2026-11-21). The
earlier `STATE_RESOLVED` OH loss remains the default of `input/i3.in` and stays runnable for every case,
with this deck too (step 6).

The original runs are `validation_runs/i3_qfit_campaign_20261123` (not in git):

- jobs 138167-138172 on Voyager, 128 ranks each;
- submitted 2026-11-23, all converged by 2026-11-26.

Created 2026-11-26 at the user's request:

> i want the compain we ran here be capable to run out of git on another machine

> keep thecase with the state resolved OG chamistry -- it should still remain runnable with cases. add a
> comment in the .in and in i3's README that the comde is consistent wil MAVEN/IUVIS observations for all
> three OH,H, and O 6 cases with the OH chemistry that was implemented later

## Result: consistency with MAVEN/IUVS

With this deck and build, the model is consistent with the MAVEN/IUVS observations of OH (0-0), H
Lyman-alpha and O I 130.4 nm for the six cases. The numbers below are for the last (converged) output of
each run. The model is edge-corrected and shown as run, not scaled.

| Case | Q(H2O) (s^-1) | output (iterations) | OH model/obs | H Ly-alpha model/obs | O I 130.4 nm at 1e4 / 1e5 km (R) |
|---|--:|---|---|---|---|
| 1 | 1.853e28 | 6 (70000) | 1.03 (0.86-1.22, 4 points) | 0.62 (0.47-0.89, 11) | 0.388 / 0.161 |
| 2 | 3.775e28 | 6 (70000) | 1.03 (0.82-1.29, 8) | 0.79 (0.54-1.33, 14) | 1.19 / 0.504 |
| 3 | 8.893e28 | 6 (70000) | 1.05 (0.91-1.21, 5) | 1.15 (0.73-1.91, 14) | 9.41 / 4.45 |
| 4 | 2.102e29 | 8 (90000) | 1.00 (0.82-1.14, 3) | 1.08 (1.04-1.22, 7) | 205 / 92.0 |
| 5 | 1.527e29 | 6 (70000) | 1.03 (0.99-1.08, 2) | 0.78 (0.73-0.91, 8) | 127 / 56.4 |
| 6 | 1.483e29 | 9 (100000) | 1.02 (0.97-1.08, 2) | 0.84 (0.76-0.94, 9) | 156 / 68.9 |

How to read the table:

- **Model/obs columns:** geometric mean over the MAVEN/IUVS points, followed by the range and the number
  of points.
- **OH:** Q(H2O) was fitted to OH, so these ratios are close to 1 by construction.
- **O I 130.4 nm:** case 5 is the only case with an IUVS O profile (digitized,
  `srcHartley2/i3/tools/inputs/`). IUVS / model per bin is 0.60, 1.08, 1.23, 0.93, 0.80, 0.73 and 0.85,
  with geometric mean 0.87 (the fit predicted 0.95). For cases 1-4 and 6 the O values are model
  predictions.
- **O model:** direct (KVB2022) + Ly-beta Bowen, optically thin. The line-centre optical depth exceeds 0.1
  in all case-5 IUVS bins.
- **What I3 leaves out:** O from CO/CO2 and electron impact.

The g-factors and the chemistry are described in `srcHartley2/i3/README.md`, section "Chemical model and
g-factors in detail".

The figures `qfit_OH_H_O_vs_maven_final.png` and `qfit_o_case5_vs_iuvs_final.png` are made by
`srcHartley2/i3/tools/make_qfit_figures.py --output final` (step 5). All numbers are in
`expected_results.json`.

## Files

| File | Content |
|---|---|
| `i3.conf` | Build input: `input/i3.conf` of `3f71104` (both OH-loss variants) with `BlockCells=12,1,1`; byte-identical to the campaign's |
| `i3.in` | Runtime deck of the six cases: the campaign's deck plus a 17-line comment at the top (consistency, `STATE_RESOLVED`). After line 17 it is byte-identical to the deck the runs used |
| `SHA256SUMS` | Checksums of `i3.conf` and `i3.in`; the scripts check them |
| `expected_results.json` | Reference results (the table above, the case-5 bins, the SHA-256 of every brightness file) |
| `scripts/compile_amps.sh` | Clone, configure and compile AMPS at a given commit with `i3.conf`, in a separate directory |
| `scripts/validate_inputs.sh` | `--validate-runtime-input`, six cases x both OH-loss variants |
| `scripts/run_case.sh` | Set up a run directory and run one case (also a short test run, and `STATE_RESOLVED`) |
| `scripts/production.sbatch` | SLURM example for one case |
| `scripts/compare_results.py` | Compare a rerun with `expected_results.json` |

`srcHartley2/i3/reference/oh_lifetime_sa1982.csv` is the `SA1982_VH` table. The deck names it relative to
itself, and `run_case.sh` copies it into the run directory.

## How the inputs were made

**Code.** Commit `848e9c5`. The campaign executable was built from `c42b555` plus the DIM1
domain-boundary fix patch. `848e9c5` is that fix: its `src/pic/pic_mover.cpp` is byte-identical to the
patched file, and its other changes are documentation only.

The executable was built with:

```text
Config.pl -install -compiler=gfortran,gcc_mpicc
Config.pl -application=i3 -spice-path=nospice -spice-kernels=nospice -amps-test=off
rm -rf -- build
make -j32 EXTRACOMPILEROPTIONS=-O3
```

The build environment was GCC 11.5.0 and Open MPI 4.1.7rc1. Build fingerprints, which `compile_amps.sh`
compares:

- `build/pic/picGlobal.dfn` without its two build-date comment lines: `d82c0e8a...`;
- `build/pic/pic_mover.cpp`: `c9e0b745...`.

The campaign's `amps` had SHA-256 `5c3534f9...`; another build has another one.

**Build input `i3.conf`.** The repository's `input/i3.conf` keeps `BlockCells=8`, which the nightly test
references require. This campaign uses 12 cells per block: the finest cell of the 4e9-m domain is then
0.310 m instead of 0.466 m.

**Deck `i3.in`.** It is `prepare_i3_cases.py --oh-loss sa1982_vh` output with these changes:

| Setting | `input/i3.in` (generator) | this deck |
|---|---|---|
| `chemistry.variant.OH_LOSS` | `STATE_RESOLVED` | `SA1982_VH` (both rate sets are present) |
| `mesh.outer_radius_m`, brightness grid end | 3e9 m | 4e9 m |
| `mesh.surface_resolution_m`, cases 2-6 | 0.3527 m (floor for `BlockCells=8`, 3e9 m) | 0.3135 m (floor for `BlockCells=12`, 4e9 m) |
| `brightness.velocity_domain_policy` | `ABORT` | `EXCLUDE_AND_REPORT` (a rare O atom beyond the +-48 km/s of KVB2022 is excluded from that line and counted, instead of stopping the run) |
| `chemistry.table.OH_TAU_SA1982.file` | `../srcHartley2/i3/reference/oh_lifetime_sa1982.csv` | `oh_lifetime_sa1982.csv` (copied into the run directory) |
| `source.H2O.rate_s-1` | archive rate x `reference/i3_rate_scaling.csv` | fitted to MAVEN/IUVS OH (below) |

The 4e9-m domain, the 0.3135-m floor and `EXCLUDE_AND_REPORT` date from the O I 130.4-nm campaign of
2026-11-19.

**Fitted Q(H2O).** Each case's rate is the previous rate times a fit factor:

- The fit factor is the geometric mean of MAVEN/IUVS OH / edge-corrected AMPS OH.
- It was computed at output 4 of the `SA1982_VH` runs of `validation_runs/i3_fixed_campaign_20261122`;
  those runs were cancelled once the fit was made.
- The factors, by case, are 0.5865, 0.5606, 0.5650, 0.4647, 0.4362 and 0.3728.
- The new rate is rounded to 4 significant digits.

The fit numbers are in git: `srcHartley2/i3/tools/inputs/ohfit_numbers_output4.json`.

CO and CO2 rates, seeds, chemistry rates, mesh and brightness settings are those of the generator.

## Rebuilding and rerunning on another machine

**Requirements:**

- git access to the AMPS repository and to `github.com:SWMFsoftware/share` and `util`. The script clones
  them; share is checked out at `c5f8beef`, as in the campaign build.
- gfortran, gcc and MPI (`mpicxx`, `mpiexec`).
- Python 3 with numpy and matplotlib, for the figures.
- No SPICE, no network access at run time.

**Resources per case:**

- 128 MPI ranks, about 1.8 GB per rank (about 225 GB per node);
- 1.7-3.1 days of wall time on the 128-core nodes of the campaign (convergence after 70000-100000 iterations);
- about 50-60 MB of output per run directory.

Commands are given from the AMPS root of the checkout. `C=srcHartley2/i3/campaigns/qfit_20261123`.

1. **Build** (1-2 min with 32 jobs in the test). The build tree is a separate clone, so the checkout is not
   changed:

   ```sh
   $C/scripts/compile_amps.sh --dest /path/to/build_qfit --jobs 32
   ```

   By default this builds commit `848e9c5`, the exact campaign code. `--commit <rev>` builds another
   commit. Master `8053739` gave byte-identical test-mode results ("Reproduction test"). The script stops if the compiled
   block layout is not `_BLOCK_CELLS_X_ 12` and `_MAX_REFINMENT_LEVEL_ 30`. It writes
   `build_provenance.txt` and compares the build fingerprints with the campaign build.

2. **Check the inputs** (12 checks, under 1 min):

   ```sh
   $C/scripts/validate_inputs.sh --amps /path/to/build_qfit/amps --work /path/to/validate
   ```

3. **Optional short test** (1000 iterations, 8 ranks, about 1 min):

   ```sh
   MPIRUN="mpiexec --bind-to none" $C/scripts/run_case.sh --amps /path/to/build_qfit/amps --case 1 \
       --out /path/to/test --np 8 --short 1000
   ```

4. **Production**, six runs of 128 ranks. Either run each case directly:

   ```sh
   $C/scripts/run_case.sh --amps /path/to/build_qfit/amps --case N --out /path/to/qfit_rerun
   ```

   or submit them with SLURM (adapt the partition and time limit to the machine):

   ```sh
   for N in 1 2 3 4 5 6; do
     sbatch --partition=<partition> --time=7-00:00:00 --ntasks=128 --job-name=i3qfit-case$N \
       --export=ALL,CASE=$N,AMPS=/path/to/build_qfit/amps,OUT=/path/to/qfit_rerun $C/scripts/production.sbatch
   done
   ```

   Each run directory is `<out>/sa1982_vh/caseN_np128/production`. It holds the executable, deck, build
   input, OH table, `SHA256SUMS` and `provenance.txt`, and the run log `caseN.run.log`. A run stops by
   itself when it has converged: the relative L1 change stays below 0.05 for 3 consecutive 10000-iteration
   windows.

5. **Figures and comparison:**

   ```sh
   python3 srcHartley2/i3/tools/make_qfit_figures.py --output final --campaign /path/to/qfit_rerun \
       --report-dir /path/to/qfit_rerun/report
   python3 srcHartley2/i3/tools/make_o_gfactor_figure.py --output final --campaign /path/to/qfit_rerun \
       --report-dir /path/to/qfit_rerun/report
   python3 $C/scripts/compare_results.py /path/to/qfit_rerun/report/data/qfit_numbers_final.json
   ```

   A rerun is stochastic: another machine, compiler or rank count gives other random sequences, so it is
   not bit-identical. `compare_results.py` lists rerun / reference for OH, H, O and the case-5 bins, and
   flags differences larger than 5% (`--tolerance`).

6. **The earlier OH chemistry (`STATE_RESOLVED`)** is runnable for every case. With this deck, add
   `--chemistry state_resolved` to `run_case.sh`; this passes
   `--chemistry-variant OH_LOSS=STATE_RESOLVED` and writes to `<out>/state_resolved/...`. The Q(H2O) of this
   deck were fitted for `SA1982_VH`, so `STATE_RESOLVED` does not fit OH with them: its OH lifetime is 1.7-2.5
   times shorter.

   The repository's `input/i3.in` runs `STATE_RESOLVED` by default, with the generator's source rates:

   ```sh
   ./amps -i input/i3.in -case N
   ```

   Build it from `input/i3.conf` (`BlockCells=8`), as described in `srcHartley2/i3/README.md`.

## Relation to the nightly test

The nightly test `i3` (`MakefileTest/Table`) does not run this campaign. It is a regression check of the
code output:

- it builds `input/i3.conf` (`BlockCells=8`) in test mode;
- it runs case 1 of `input/i3-test.in`, with both OH-loss variants, for 2000 iterations on 4 ranks;
- it has no O lines and a 3e9-m domain.

The physics code is the same, but the campaign's block size, domain, O emission records and Q(H2O) are
covered only by the checks in this folder.

## Reproduction test (2026-11-26)

The test was run in `validation_runs/i3_qfit_repro_test_20261126` (not in git), jobs 139743 and 139748, on
one 128-core Voyager node (GCC 11.5.0, Open MPI 4.1.7rc1). It used a fresh clone of master `8053739` with
the files of this folder, the restored `input/i3.conf` and `input/i3.in`, and the updated tools copied in.
Only the scripts of this folder were used:

- **Build:** `compile_amps.sh` (default commit `848e9c5`) took 1 min. Both fingerprints equal those of the
  campaign build.
- **Input checks:** `validate_inputs.sh` passed 12 of 12 (six cases x `SA1982_VH` and `STATE_RESOLVED`).
- **Short runs:** `run_case.sh --short 1000 --np 8` ran case 1 and case 5 with `SA1982_VH`, and case 1 with
  `STATE_RESOLVED`. All finished, wrote the 10 output tables and logged the selected OH-loss variant.
- **Later code:** the test-mode builds of `848e9c5` and of master `8053739` gave byte-identical results in
  all 10 output tables. The run was case 1 with this deck, `SA1982_VH`, 2000 iterations and 4 ranks. Test
  mode is deterministic (the nightly references rely on that), so it allows a byte comparison; production
  builds were not compared byte for byte. So `--commit 8053739` builds the same I3 code as the campaign
  for this path.
- **Not checked:** `STATE_RESOLVED` between the two commits, the other cases, and full-length runs.
- **Commits after `8053739`:** not checked. To use one, repeat the comparison:
  1. `compile_amps.sh --test-mode`, once with `--commit 848e9c5` and once with `--commit <rev>`;
  2. `run_case.sh --short 2000 --np 4 --case 1` with each build;
  3. compare the `PT/plots/*.dat` files byte for byte.
- **Figures from the clone:** `make_qfit_figures.py --output final` and `make_o_gfactor_figure.py --output
  final`, run from the clone on the original runs, gave figures and CSV files byte-identical to those made
  in the working tree. `compare_results.py` gave 0.000.
