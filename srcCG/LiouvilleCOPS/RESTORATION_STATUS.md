# Restored baseline and physical-template development

Development started 2026-10-19 in three concurrent, separately owned tracks:
gravity repair, joint fitting, and baseline campaign preparation. This record
distinguishes implemented code from completed scientific production results.

## Current deliverables

- `GRAVITY_REPAIR.md`: independent mesh audit, repaired DISABLE gravity paths,
  cache invalidation, regression tests, and pending full-field validation.
- `JOINT_TEMPLATE_FIT.md`: jointly fitted degree-0–4 activity, temperatures,
  CO2 multiplier, and nonnegative physical-template amplitudes; combined
  regularization and exact per-facet C++ export.
- `RESTORED_BASELINE.md`: read-only campaign preparation and explicit commands
  for baseline B, NG-only N, cliff-only C, and combined NC comparisons.

The fitter checks repaired gravity provenance when selecting cliff columns.
Renaming an old bundle cannot make its cliff columns valid. The geometry-only
`late_ng_not_rg` column is unaffected by gravity and remains selectable for N.
The baseline B and NG-only N calculations have no physical dependency on the
corrected cliff products.

## Integrated verification

The completed development suite ran 104 tests successfully (97 passed, seven
opt-in tests skipped). The skipped set comprises six mission/solver integration
tests and the external C00 gravity fixture test. That C00 test was exercised
separately on eight targets and passed its 1e-12 vector-relative gate. The
suite includes full synthetic joint result/plot writing and a compiled C++
facet-reader replay of the exported joint field; this is not a substitute for
the full mission SPICE replay required for production acceptance.

The preparation tool was run against the actual archived and 788k manifests:
930 NG and 1146 RG selected observations are common and unchanged, and point
IDs 1146–1149 are newly selected for each gauge. Non-shape recorded dependency
hashes match. B and N have no preparation blockers when the quarantined bundle
is selected for its unaffected NG-only column; C and NC await corrected gravity
products. Preparation is read-only and the fitter still performs the full
input/cache validation on launch.

## Quarantined evidence

The following original names under the current external `cops_output` root
were renamed by appending `_heal_defect` (after the existing extension):

- `physical_templates_shap5_788k.gravity.npz`
- `physical_templates_shap5_788k.gravity_100k_rejected.npz`
- `physical_templates_shap5_788k.npz`
- `physical_templates_shap5_788k.json`
- `physical_templates_shap5_788k.png`

No bytes were deleted. Existing harmonic-only results are unaffected. The
earlier 100k rejection must be reassessed using corrected gravity rather than
being interpreted as a validated resolution study.

## Production gates still outstanding

1. Regenerate and validate the corrected full 788k gravity/slopes/cliff bundle,
   or import an independently validated exact-mesh field with complete provenance.
2. Run the restored 29-parameter 788k baseline and inspect convergence, bounds,
   clipping, and differences on common archived observations.
3. Run the 30-parameter N and C and 31-parameter NC comparisons, and repeat
   scientifically important comparisons with the transition-cut profile.
4. Validate final selected pressures using the C++ forward model and compare
   the full mission and final-descent intervals separately.

The archived fit stopped at its evaluation limit and hit unusual temperature
bounds. Restoring its parameter freedom is a controlled reproduction of its
numerics, not certification of the inferred temperatures or source map.

Warm-start support, resumable optimizer checkpoints, and eager output from each
lambda of a scan are not included in this first joint-physics change. Existing
initialization remains the original run coefficients/temperatures and zero
template amplitudes. Additional starts must not redefine the original prior.

Use 36 threads for appropriately parallel CPU calculations, subject to actual
backend control (polyhedral gravity uses TBB). One SciPy lambda fit remains
one optimizer process; avoid 36-fold kernel replication. This development
does not install or modify CUDA software and does not claim a GPU speedup.

## CPU production launch: 2026-10-19 18:40 UTC

The requested independent B29, N30, and corrected-gravity jobs were launched
under a detached supervisor. This is a launch record, **not** a statement of
completed optimization or accepted cliff products. Consult the live status:

```text
/data/vtenishe/rosetta_data_end_of_mission_cops/cops_output/
  restored_cpu_campaign_20261019T184004Z/
    status.json
    B.log
    N.log
    gravity.log
    supervisor.log
    preparation.json
    launch.json
    launch_recipe.py
    source_hashes.json
    code_snapshot/
```

`status.json` records job PIDs and, after termination, exit codes. Exit zero
means the command finished, not necessarily optimizer convergence: inspect
the fit's stopping reason, optimality and bounds. The logs may be quiet during
cache validation, loading, or optimization; the fitter prints final optimizer
status rather than every iteration. The operational recipe and source snapshot
are preserved outside the application for reproducibility. Do not rerun the
recipe into occupied output paths.

Outputs beneath the same external root are:

- `fit_shap5_788k_restored_B`: 29 free parameters, no physical templates.
- `fit_shap5_788k_restored_N`: 30 free parameters, only `late_ng_not_rg` selected
  from the quarantined mixed bundle. No corrupted gravity/cliff column is used.
- `physical_templates_shap5_788k_corrected.npz` and companion JSON/PNG: intended
  corrected bundle, with both 45- and 60-degree cliff hypotheses and NG mask.
- `physical_templates_shap5_788k_corrected.gravity.npz`: intended reusable
  audited-DISABLE CPU gravity cache on the exact 788k transport surface.

Both fits use `all_available`, lambda zero, degree/enhancement-degree 4,
independent 30–400 K temperature bounds, 1000 maximum evaluations, and 1e-10
ftol/xtol/gtol. All baseline coefficients, temperatures and kappa remain free.
Each optimizer has one worker and one BLAS/OpenMP thread; the two fits are
independent concurrent processes, not parallel iterations of one optimizer.

The gravity process is restricted with CPU affinity to CPUs 0–35 and runs at
nice priority 10 so the two fits can progress. Its package uses TBB; the
affinity is the actual resource restriction, not an assumption that OpenMP
variables control TBB. Density is 533 kg/m³, outward offset 0.1 m, chunk size
2048, and the offset-halving slope gate is 0.25 degrees. The existing builder
performs the independent topology audit and provenance checks. No CUDA code
or imported gravity field is involved.

Before launch, four joint-fit tests, two cliff-provenance tests, and seven
gravity/topology tests passed; the external C00 test was skipped in that quick
run. Full-field comparison against all C00 reference targets and final mission
replay remain production acceptance work after products are available.
