# Joint baseline and physical-template fitting

This implementation restores free degree-0–4 H2O activity coefficients while
allowing nonnegative coefficients for declared per-facet physical templates.
With two templates, free analytical temperatures, and a free CO2 scale, the
parameter vector contains 31 entries: 25 activity coefficients, two template
amplitudes, three temperatures, and kappa. No higher-degree correction is needed.

The field is `max(0, Y a + A_scale B c)`. There is exactly one clipping
operation, applied after adding all terms. Template masks remain fixed and
are evaluated throughout the trajectory. The analytic Jacobian differentiates
this same expression; its activity derivative is zero at and below the clip.
Temperature and kappa derivatives include the current, template-modified H2O
field. Template amplitudes have lower bound zero; baseline coefficients are
unbounded, as in the archived fitting scheme.

`A_scale` is the area-weighted RMS of the original run's *unclipped* baseline
activity. It does not change while fitting. Baseline optimizer coordinates are
scaled by the largest absolute original H2O coefficient. Template coordinates
are dimensionless amplitudes multiplying `A_scale`. Independent temperature
bounds remain the existing defaults unless overridden; ordered temperature is
an explicit alternative and is not silently enabled.

## Regularization and numerical interpretation

For positive lambda the joint penalty is the surface integral of the squared
TOTAL unclipped change from the original run map, divided by the original
map's squared surface norm. It includes baseline/template and template/template
cross terms. An eigenvalue square root of the symmetric Gram matrix supports
dependent predictors without adding a fictitious ridge penalty. Negative
roundoff eigenvalues are clipped to zero; significantly negative eigenvalues
raise an error. Data Jacobian ranks must still be inspected: a valid map norm
does not guarantee unique fitted coefficients. With lambda zero the penalty
does not enter optimization.

Existing legacy-only and fixed-baseline harmonic fits retain their existing
regularizers. Higher-degree harmonics still require all baseline, temperature,
and kappa parameters fixed; that distinct experimental configuration has not
been broadened by this change.

## Running and outputs

After validating or rebuilding the physical-template bundle, run from this
application directory (the bundle name below must identify repaired products):

```bash
export COPS_OUTPUT_ROOT=/data/vtenishe/rosetta_data_end_of_mission_cops/cops_output
python3 tools/run_fit_profile.py \
  --run output_fit_shap5_788k_base --profile all_available \
  --degree 4 --enhancement-degree 4 \
  --activity-template validated_templates.npz late_ng_not_rg \
  --activity-template validated_templates.npz cliff_ge_60deg \
  --lambda 0 --jobs 1 --threads-per-process 1 \
  --max-nfev 1000 --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 \
  --check-tolerance 1e-6 --out fit_788k_joint_templates
```

Do not pass `--fix-legacy-activity`, `--fix-temperature`, or `--fix-kappa`
for the 31-parameter case. Remove template arguments for the 29-parameter
baseline; select one template for either 30-parameter comparison. Ordinary
SciPy iterations remain serial. Independent fits can use separate workers
within measured memory limits. C++ verification may use 36 OpenMP threads.

`ActivityFit.txt` contains the fitted baseline coefficients (and scaled CO2).
`ActivityEnhancementFit.txt` contains the explicit residual per facet and binds
it to that fitted baseline file's SHA256 identity. `cops_fit.in` selects both.
Activity diagnostic maps show the total change from the original run prior,
including baseline changes and templates. Temperature limits, clipping, rank,
and nonnegative-amplitude bounds require scientific review. Optimizer stopping
alone does not validate the physical interpretation.

## Tests and remaining work

`tests/test_joint_template_fit.py` checks all 31 joint analytic derivative
columns, zero-template equivalence, and the area-integral regularization
including dependent columns. It runs the complete synthetic result writer,
checks fitted baseline hashes and report fields, produces comparison plots,
and (when g++ is available) compiles the independent C++ facet reader to replay
a facet exported from the same fitted state. Existing facet tests additionally
check rejection behavior. This
development does not itself execute a production fit or validate cliff gravity.
Warm starts, resumable optimizer checkpoints, and per-lambda eager output
remain separate work; the initial state is still the original run parameters
and zero template amplitudes. The prior must remain distinct from initialization
when warm-start support is added.
