# I3 comparison figures: `make_qfit_figures.py`, `make_o_gfactor_figure.py`

Draws the comparison of the I3 (3I/ATLAS) AMPS runs with MAVEN/IUVS for all six cases:

- `qfit_OH_H_O_vs_maven_out<K>.png`: 3 x 6 panels. Rows: OH (0-0) and H Lyman-alpha against the
  MAVEN/IUVS points, and O I 130.4 nm (direct KVB2022 + Ly-beta Bowen). The case 5 O panel also shows the
  digitized MAVEN/IUVS O I 130.4-nm profile.
- `qfit_o_case5_vs_iuvs_out<K>.png`: case 5 O I 130.4 nm against the IUVS profile and the IUVS team
  model of the same slide, with observed / model per IUVS bin.

The runs are those of the SA1982_VH OH-loss campaign with Q(H2O) fitted to the MAVEN/IUVS OH
(`validation_runs/i3_qfit_campaign_20261123`); they are shown as run, not scaled. The deck, build input and
scripts to rerun that campaign on another machine are in `../campaigns/qfit_20261123/`.

`make_o_gfactor_figure.py` draws the O I 130.4-nm g-factor of the same runs (section "O I 130.4-nm g-factor"
below).

## Requirements

- Python 3 with numpy and matplotlib (tested with Python 3.12.9, numpy 2.3.2, matplotlib 3.10.5).
- This folder and `srcHartley2/visualization` of the same checkout (edge correction; MAVEN/IUVS OH and H
  points in `srcHartley2/visualization/references/i3_brightness_MAVEN`).
- The AMPS run output (not in git). The campaign directory must hold, for each case N = 1-6,
  `sa1982_vh/caseN_np128/production/` with
  - `i3.in` (runtime input: Q(H2O), r_h, r-dot),
  - `i3.conf` (compile-time input of the edge correction; if missing, the repository's `input/i3.conf`),
  - `caseN.run.log` (convergence state),
  - `PT/plots/pic.{H,O,OH}.brightness.s={0,2,3}.out=K.dat` for the chosen output K, and the density
    profiles `PT/plots/pic.<species>.s=*.out=K.dat` that the edge correction reads.

## How to run

From the AMPS repository root (any directory works; paths are found from the script's location):

```
python3 srcHartley2/i3/tools/make_qfit_figures.py
```

This uses the newest output written by all six runs of the default campaign and writes to
`validation_runs/i3_qfit_campaign_20261123/report_20261124/`. It takes about 1 minute on the login node.
For converged runs use `--output final`: the last output of each run, i.e. the window in which it converged
(the runs converge after different numbers of iterations).

Options:

| Option | Default | Meaning |
|---|---|---|
| `--output K` or `--output final` | newest output written by all six runs | output index; output K covers iterations K x 10000 + 1 to (K + 1) x 10000. `final`: the last output of each run (file names end in `_final`) |
| `--campaign DIR` | `validation_runs/i3_qfit_campaign_20261123` | campaign directory (layout above) |
| `--np NP` | 128 | MPI ranks in the run-directory names `caseN_np<NP>` |
| `--report-dir DIR` | `<campaign>/report_20261124` | where `plots/` and `data/` are written (created if missing) |
| `--predicted FILE` | `inputs/ohfit_numbers_output4.json` | fit numbers of the cancelled runs ("predicted by the fit" values, Q(H2O) check) |

Examples:

```
# a given output
python3 srcHartley2/i3/tools/make_qfit_figures.py --output 4

# the converged (last) output of each run
python3 srcHartley2/i3/tools/make_qfit_figures.py --output final

# write somewhere else, leaving the report folder untouched
python3 srcHartley2/i3/tools/make_qfit_figures.py --output 4 --report-dir /tmp/i3_figures

# list the options
python3 srcHartley2/i3/tools/make_qfit_figures.py --help
```

Files with the same output index are overwritten; other outputs are kept (names end in `_out<K>`).

## Output

In `<report-dir>`:

| File | Content |
|---|---|
| `plots/qfit_OH_H_O_vs_maven_out<K>.png` | the 3 x 6 figure (300 dpi) |
| `plots/qfit_o_case5_vs_iuvs_out<K>.png` | the case 5 O figure (300 dpi) |
| `data/qfit_numbers_out<K>.json` | per case: run, output, iterations, Q(H2O), SHA-256 of every brightness file, OH and H model/obs (geometric mean, median, min, max), O at 1e3-1e6 km, O optical-depth range, predicted values; case 5 IUVS / model per bin |
| `data/qfit_o_case5_vs_iuvs_out<K>.csv` | case 5: IUVS radiance, model bin mean and their ratio per bin |

The script prints, per case, model/obs for OH and H and the O brightness at 1e5 km, and the case 5
IUVS / model per bin.

## What the script checks (it stops with a message otherwise)

- output K exists for all six runs (`output K of the OH-fitted run is incomplete`);
- the output headers record the SA1982_VH OH loss (`the run is not SA1982_VH`);
- each run's input Q(H2O) is the fitted Q(H2O) of `--predicted` within 0.1 % (`input Q(H2O) ... is not
  the fitted ...`). For runs with other source rates, give a matching `--predicted` file;
- all IUVS bins lie inside 2/3 of the domain radius.

## Method

- All brightness is edge-corrected (`srcHartley2/visualization/comet_edge_correction.py`) and shown to
  2/3 of the domain radius.
- OH and H model/obs: model interpolated (log-log) to the MAVEN/IUVS points inside the tabulated range.
- Case 5 O: model averaged over each IUVS bin along the slit (1-D), then observed / model.
- O is optically thin brightness; the panel text gives the impact parameter inside which the line-centre
  optical depth exceeds 0.1. Model O comes from H2O only (no CO/CO2 photolysis or electron impact in I3).
- Runs that have not converged are labelled "(running)" and "NOT converged".

## Files in this folder

| File | Content |
|---|---|
| `make_qfit_figures.py` | the script (copy of `validation_runs/i3_qfit_campaign_20261123/report_20261124/scripts/make_qfit_figures.py`, made to run from git; `--output final` and `--np` added 2026-11-26) |
| `make_o_gfactor_figure.py` | O I 130.4-nm g-factor figure (copy of the script of the same name in that report, made to run from git, 2026-11-26) |
| `i3_figure_helpers.py` | table reader, run loader, edge-correction call, bin average, palette and axis style (copied from the I3 report scripts in `validation_runs/`) |
| `inputs/` | case 5 IUVS O I 130.4-nm profile, the slide's model curve, the OH fit numbers; sources and SHA-256 in `inputs/README.md` |

The input folder is called `inputs/`, not `data/`: the repository's `.gitignore` excludes every directory
named `data`.

## O I 130.4-nm g-factor: `make_o_gfactor_figure.py`

```
python3 srcHartley2/i3/tools/make_o_gfactor_figure.py [--output K|final] [--campaign DIR] [--report-dir DIR] [--np NP]
```

Same options and defaults as above, except that `--output` defaults to 4. The figure has three panels:

- (a) g at 1 AU against heliocentric radial velocity: direct (KVB2022) and Ly-beta Bowen, with the case-5
  solar fluxes, and the six cases at their r-dot;
- (b) g of each case at r-dot and r_h;
- (c) the g-factor applied in the run (the column-weighted "Effective Photon g Factor" of the brightness
  tables) divided by the g of (b).

The g evaluation repeats `src/models/resonance_fluorescence/resonance_fluorescence.cpp` on the compiled tables
(`generated/reference_tables.inc`).

Output, with tag `out<K>` or `final`:

- `plots/o1304_gfactor_<tag>.png`;
- `data/o1304_gfactor_cases_<tag>.csv`;
- `data/o1304_gfactor_<tag>.json`;
- `data/o1304_gfactor_1au_case5_fluxes.csv`.

Sources and method: `srcHartley2/i3/README.md`, section "Chemical model and g-factors in detail".
