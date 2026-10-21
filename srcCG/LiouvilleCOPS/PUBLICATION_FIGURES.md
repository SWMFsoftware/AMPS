# CPU result figures matching the CUDA pressure layout

`tools/cops_publication.py` implements the pressure-result layout from
`LiouvilleCOPS_cuda/cops_cuda/figures.py`, locally and without importing the
CUDA branch, CuPy, SPICE, or a GPU runtime. This is a presentation change:
no optimizer, selected observation, pressure prediction, kernel, or physical
parameter changes. CUDA **ground-track, impact, and mesh/template rendering**
are separate geometry products and are not ported by this change. Existing CPU
activity-map and residual diagnostics remain available unchanged.

## What is drawn

- Time: three rows (NG, RG, NG/RG), full descent and last-25-minute columns.
  The zoom ends at the last observation, even if the highlighted window ends
  earlier. `--zoom-minutes` changes only the view.
- Altitude: the same rows, logarithmic altitude reversed so descent proceeds
  left to right. All series use the **reference table's altitude**. For different
  meshes this is a common trajectory coordinate, not a claim that their
  nearest-surface distances are identical.
- Pressures are in Pa; pressure and ratio axes are logarithmic. Invalid,
  nonpositive, and inside-nucleus predictions break model lines. No smoothing,
  interpolation, or new pressure evaluation is performed. NG/RG is diagnostic,
  not a third objective. Its marker is filled only if both gauges were selected.
- Selected observations are filled black; excluded ones are hollow. Saturated
  readings retain their actual selection flags. Existing A5 tables/plots
  separately mark saturation and retain selection reasons.
- The optional unfitted base is dashed gray; fits use the Okabe–Ito palette.
  Dashed/dotted curves are above solid curves so coincident solutions remain
  distinguishable. Shading is opaque for EPS compatibility.
- Each figure is PNG (300 dpi) and EPS. Styling is scoped to these plots and
  does not affect other diagnostics.

## Automatic outputs

Ordinary `fit_cops.py` fits (including profile-runner fits and each lambda
directory) and `fine_tune_activity.py` add:

```text
figures/pressures_time.png
figures/pressures_time.eps
figures/pressures_altitude.png
figures/pressures_altitude.eps
figures/pressures.json
```

Files are indexed by the fit/fine-tuning manifest; ordinary fits also index
them in the A5 report. Ordinary fits use the selected `late_ng_not_rg` template
definition for shading, or default to 2016-09-30 10:33:10–10:39:26 UTC, matching
the GPU publication default. This default is only a visualization window:
it does not select data or imply that a visibility template was fitted.

Fine tuning overlays its frozen fit (dashed) and corrected fit (solid), and
highlights its actual `end_utc - minutes` through `end_utc` interval. Zero
eligible groups legitimately yield coincident curves. Its older
`pressure_comparison.png` is retained.

Plotting failures are recorded in `figures/pressures.json` without changing
optimizer success. A short synthetic trajectory outside the mission window,
for example, cannot supply that window's altitude range; use explicit bounds
with the standalone command. When imported by the CUDA compute driver, the CPU
fitter defers automatic publication output to that driver, avoiding filename
collisions and preserving its `--no-figures` behavior. No CUDA source was
changed for this implementation.

## Plot existing output, without rerunning fitting

Run from LiouvilleCOPS. Relative input/output names follow `COPS_OUTPUT_ROOT`
(default `<AMPS>/cops_output`); outputs inside LiouvilleCOPS are refused.

```bash
python3 tools/pub_figures.py pressures \
  --base output_fit_base --fit 'Fitted model=completed_fit:solid' \
  --late-start 2016-09-30T10:33:10 --late-stop 2016-09-30T10:39:26 \
  --zoom-minutes 25 --title 'COPS pressure reconstruction' \
  --out publication/pressures
```

The completed external 199k campaign can be plotted with:

```bash
export COPS_OUTPUT_ROOT=/data/vtenishe/rosetta_data_end_of_mission_cops/cops_output
python3 tools/pub_figures.py pressures \
  --base cpu_199k_full_20261020T214018Z/base \
  --fit 'Joint fit=cpu_199k_full_20261020T214018Z/joint_fit:dashed' \
  --fit 'Fine tuning=cpu_199k_full_20261020T214018Z/fine_tuning:solid' \
  --zoom-minutes 25 --title 'SHAP5 199k: joint fit and late fine tuning' \
  --out cpu_199k_full_20261020T214018Z/publication/pressures
```

This creates new publication products without modifying inputs, restart
archives, or historical reports.

| Option | Meaning |
|---|---|
| `--base RUN` | Required original base directory or `COPS_model.dat`; unfitted curve and fine-tuning geometry. |
| `--fit LABEL=DIR[:style]` | Required, repeatable fit/forward directory, model table, or completed fine-tuning directory. Styles: `solid`, `dashed`, `dotted`, `dashdot`. |
| `--ref RUN` | Observation values/selection and altitude coordinate; default last `--fit`. Does not change predictions. |
| `--no-base` | Hide the unfitted curve. |
| `--bundle FILE` | `.npz`/`.json` metadata with `late_ng_not_rg.interval_utc` defining the shading window. |
| `--late-start UTC --late-stop UTC` | Paired explicit bounds, overriding bundle/default. Must include a positive-altitude outside-surface row. |
| `--zoom-minutes N` | Positive finite late-time panel width; default 25. |
| `--out PREFIX` | Required prefix; creates `_time.{png,eps}`, `_altitude.{png,eps}`, `_figures.json`. |
| `--title TEXT` | Optional title. |
| `--overwrite` | Permit replacing this prefix's figure products. Otherwise existing targets cause refusal before figures are written. |

When the reference is fine tuning and neither bundle nor explicit bounds is
given, its recorded fine-tuning window is used. Otherwise the GPU late-window
default applies. Bundle/explicit windows annotate a view; they do not recompute
visibility masks. Exact ordered `iPoint`/UTC correspondence is required across
all series. Different-length/misaligned trajectories are rejected, not paired
by position or interpolated.

The CLI retains existing CPU table/report/manifest checks. Fine-tuning tables
must match their manifest output hash and the original base-table hash.
Historical unmanifested forward tables remain readable with legacy provenance
in the figure JSON. The JSON records input/output hashes, series provenance,
reference choice, and window. A good-looking plot is not evidence of optimizer
convergence or scientific validation; consult the fit manifest.

## Validation

```bash
python3 -m unittest discover -s tests -p 'test_publication_figures.py' -v
python3 -m unittest discover -s tests -p 'test_fine_tuning.py' -v
```

Tests cover layout, markers, ratio safety, window controls, row correspondence,
overwrite refusal, style isolation, PNG/EPS output, hash inventory, nonfatal
plot failure, and fine-tuning provenance rejection.
