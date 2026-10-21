# CPU/CUDA-compatible fitter restarts

`tools/fit_cops.py` now saves and reads the same version-1 restart format as
`../LiouvilleCOPS_cuda/cops_cuda/restart.py` and its `tools/gpu_fit.py` driver.
The CUDA directory was inspected, not modified. The CPU implementation is
`tools/fit_restart.py`; it requires neither the sibling directory nor CuPy.

## Run and resume

Ordinary optimization saves `fit_restart.npz` automatically in its result
directory. Choose a new output directory on every invocation:

```bash
python3 tools/run_fit_profile.py \
  --run output_fit_shap5_788k_base --profile all_available \
  --degree 4 --enhancement-degree 4 --lambda 0 \
  --jobs 1 --threads-per-process 1 --max-nfev 1000 \
  --out fit_cpu_initial

python3 tools/run_fit_profile.py \
  --run output_fit_shap5_788k_base --profile all_available \
  --degree 4 --enhancement-degree 4 --lambda 0 \
  --jobs 1 --threads-per-process 1 --max-nfev 1000 \
  --restart-from fit_cpu_initial/fit_restart.npz \
  --out fit_cpu_resumed
```

The restart path can instead identify a CUDA-produced `fit_restart.npz`.
Relative product paths follow `COPS_OUTPUT_ROOT` using the normal path helpers;
absolute paths are accepted. The specified fitting base must exist, have
compatible static kernels, and match the restart's exact mesh/facet ordering.
The examples use one worker because one lambda is one optimization task.

The direct fitter supports the same options:

```bash
python3 tools/fit_cops.py output_fit_shap5_788k_base \
  --profile all_available --degree 4 --lambda 0 \
  --restart-from /absolute/path/to/fit_restart.npz \
  --max-nfev 1000 --jobs 1 --out fit_cpu_from_restart
```

`--no-save-restart` disables writing the new restart, without disabling loading
an existing one. `--data-only` and `--check-jacobian` do not save fitted states.
With a restart, the Jacobian diagnostic uses the mapped restart start.

For a multi-lambda scan, every worker starts independently from the same mapped
restart (or from the run prior without one). Every `lambda=.../` result receives
its own `fit_restart.npz`; the scan root has no ambiguous single restart.
Selecting several profiles similarly starts each profile independently.

## Meaning and limitations

This is a **physical-parameter warm start**, not recovery of SciPy's internal
trust-region state. Iteration counters and stopping tolerances start afresh.
There is no periodic checkpoint during optimization. A fit that returns at
`--max-nfev` still saves its final parameters and records nonconvergence; a
process interrupted before result writing does not obtain a new restart.

The new command controls lambda, degree, bounds, fixed/free parameters,
selection, and optimizer settings. Old settings are recorded in the file,
not automatically replayed. Repeat the original options when continuing the
same experiment. A restart never changes the new run's regularization prior
(`P.x0`), its underlying CO2 map, or its selected observations.

Mapping is by parameter name, not position:

- Legacy H2O coefficients are converted to the new optimizer scaling.
- Harmonic modes and template amplitudes are matched by mode/template name.
- Physical temperatures map to free temperatures or ordered-temperature
  differences. Bounds may clip the start; every clipping decision is reported.
- Kappa is restored only if free. Fixed parameters keep the new base's values.
- New parameters absent from the file keep the ordinary problem start;
  parameters no longer free are not reintroduced.

With identical parameterization/scales, optimizer coordinates are preserved
exactly. If the enhancement normalization changes, the CPU loader multiplies
harmonic/template amplitudes by `old_A_scale / new_A_scale` to preserve their
physical contribution. This is an explicit safety improvement over the inspected
CUDA loader, which directly maps these amplitudes. For cross-backend continuation
with changed priors/scales, review the mapping rather than assuming the old CUDA
loader performs this rescaling. Same-base CPU/CUDA round trips are tested.

When no `--activity-template` is provided, the loader restores the saved list
of template source bundles and verifies their file hashes. The embedded fields
do not bypass current bundle, mesh, or repaired-cliff provenance checks.
If bundles moved, pass their new paths explicitly with `--activity-template`.
Explicit templates take precedence, but any reused template name must have
exactly the same per-facet field. A new name starts as a new parameter.

A different mesh is rejected, even if it has the same number of facets.
Cross-mesh transfer remains the separately validated CUDA
`tools/restart_tool.py translate` workflow; the CPU reader does not silently
interpolate or remap meshes.

## File format and provenance

The schema is intentionally still `liouville_cops_cuda.fit_restart`, version 1.
The compressed, pickle-free NPZ contains:

- `a` (25 H2O coefficients), `coef_prior` (3×25), `T`, and `kappa`;
- harmonic and template amplitudes, template columns, and unclipped activity;
- facet centroids and areas defining the mesh identity;
- point IDs, UTCs, altitudes, reference/predicted pressures in Pa, selection,
  and saturation flags;
- optimizer coordinates `x`, prior `x0`, and bounds `lo`/`hi`;
- `meta_json` with basis/template names, source identities, settings, optimizer
  outcome, and edit/fit history.

Numeric array SHA256 uses the CUDA definition: dtype string, shape tuple, and
contiguous C-order bytes. Reads verify recorded hashes and mesh identity and
reject malformed/nonfinite physical-state arrays. Existing restart files are
never overwritten. A large mesh produces a substantial file because geometry
and per-facet fields are included; compression adds time and memory to final
output writing.

Each CPU result also saves immutable `fit_restart_manifest.json`, an input to
the restart's provenance. The live `fit_manifest.json` indexes both products.
The snapshot avoids a circular hash dependency between that manifest and the
restart it describes. Resumed files retain earlier history and identify their
input restart. The manifest's `restart.mapping` records restored/defaulted/
clipped parameters, scale ratio, and fixed-parameter information.

The existing CUDA driver imports the CPU fitter and owns its own restart
writer. Its `compute_backend.restart` record makes the CPU writer defer to it,
preventing duplicate output or interference with its `--no-save-restart` flag.
This does not change or broaden the CUDA driver's own lambda-scan behavior.

## Tests

```bash
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 \
  python3 -m unittest discover -s tests -p test_fit_restart.py -v
```

Tests cover round trips, corruption, mesh mismatch, temperature mapping,
missing names, clipping, scaling, template-source checks, prior preservation,
writer ownership, history, opt-out, and profile-runner forwarding. When the
read-only CUDA sibling is present, tests also exercise its real reader/writer
and a complete tiny CPU CLI fit/resume/two-worker scan using its binary-format
fixture generator. Generated files remain in temporary directories. No GPU
or mission fit is required; these are numerical/format tests, not a new
scientific validation of mission results.
