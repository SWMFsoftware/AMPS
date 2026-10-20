# Restored baseline and physical-template campaign

This workflow prepares the archived successful fitting scheme on the 787,754
facet SHAP5 surface. It does not promote the archived temperatures or surface
map to physical truth. It restores the adjustable parameters so mesh and
template effects can be assessed without simultaneously freezing the baseline.

## What was verified

The reference is `output_fit_full_current16` and its sibling base in
`/home/vtenishe/ROS-END-OF-MISSION/RES-COPS-101826-6pm-7459fd3c350733bb38c87799df62a3cdcbdf9b34`.
The new base is
`/data/vtenishe/rosetta_data_end_of_mission_cops/cops_output/output_fit_shap5_788k_base`.
The read-only audit of these manifests and tables found:

- Effective settings differ only in `ShapeModel`, `OutputDirectory`, and
  `TemperatureRGTable` (the latter is a pathname change).
- Recorded hashes agree for the observation input, all 14 SPICE inputs, and
  the HK table. Shape content differs, as intended.
- All 930 archived selected NG pressures and 1,146 archived selected RG
  pressures are present and unchanged among the newly selected points.
- New-mesh all-available selection adds points 1146–1149 for both gauges,
  giving 934 NG and 1,150 RG. Consequently, compute comparison statistics on
  the common points as well as all available points; the additional final
  descent points must not silently change the meaning of the comparison.

This lightweight preparation checks existence of current dependencies and
rehashes current dependency files up to 1 MiB. It compares recorded hashes for
larger files; it does **not** claim to have revalidated their current contents
or the multi-gigabyte transport cache. The fitter's full audit remains required.

## Read-only campaign preparation

From LiouvilleCOPS, run:

```bash
python3 tools/prepare_restored_campaign.py \
  --archive-fit /home/vtenishe/ROS-END-OF-MISSION/RES-COPS-101826-6pm-7459fd3c350733bb38c87799df62a3cdcbdf9b34/output_fit_full_current16 \
  --archive-base /home/vtenishe/ROS-END-OF-MISSION/RES-COPS-101826-6pm-7459fd3c350733bb38c87799df62a3cdcbdf9b34/output_fit_full_current16_base \
  --run /data/vtenishe/rosetta_data_end_of_mission_cops/cops_output/output_fit_shap5_788k_base \
  --output-root /data/vtenishe/rosetta_data_end_of_mission_cops/cops_output \
  --ng-bundle /data/vtenishe/rosetta_data_end_of_mission_cops/cops_output/physical_templates_shap5_788k.npz_heal_defect
```

This prints JSON with commands, blockers, setting differences, and common/new
observation IDs. It never executes a command or writes a result. Supply
`--template-bundle /absolute/path/to/corrected_templates.npz` after the corrected
gravity products pass validation. An absent bundle blocks C and NC preparation.
Never execute a command with nonempty `precheck_blockers`; passing the precheck
still requires the shared fitter's authoritative cache and template checks.
Without a bundle, blocked commands contain `None` as a placeholder and are not
runnable recipes.

`--ng-bundle` is optional and applies only to N. Selecting its `late_ng_not_rg`
column is permissible from the quarantined mixed bundle because visibility
was computed independently of gravity. The old cliff columns are invalid.
The helper refuses cliff selections unless the gravity definition records
the corrected algorithm and audited DISABLE orientation policy. The fitter
also validates the actual selected template against the exact facet ordering.

## Baseline command and parameter meaning

```bash
export COPS_OUTPUT_ROOT=/data/vtenishe/rosetta_data_end_of_mission_cops/cops_output
python3 tools/run_fit_profile.py \
  --run output_fit_shap5_788k_base \
  --profile all_available \
  --degree 4 --enhancement-degree 4 \
  --lambda 0 --jobs 1 --threads-per-process 1 \
  --tbounds 30 400 --max-nfev 1000 \
  --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 \
  --check-tolerance 1e-6 \
  --out fit_shap5_788k_restored_B
```

No `--fix-*` flags are present: 25 degree-0–4 activity coefficients, three
independent analytical temperatures, and the nonnegative CO2 scale are free.
Degree 4 for the enhancement disables additional degree-5–8 harmonics. Lambda
zero reproduces the unregularized objective; temperatures retain independent
30–400 K bounds. `all_available` includes saturated NG and late data excluded
by the base's original `ExcludeNG`, matching the archived selection policy.

The initial baseline uses the base's original parameters. A separate archived
solution start should use supported initialization controls and a distinct
output directory; never replace the base prior or geometry to simulate a warm
start. The archive's 1,000-evaluation limit was reached, so the same limit is a
reproduction setting, not evidence of convergence. Assess optimizer termination,
bound hits, and sensitivity to initialization before comparing interpretations.

## Enhancement comparisons and execution dependencies

| Label | Free parameters | Additional templates | Prerequisite |
| --- | ---: | --- | --- |
| B | 29 | None | Audited 788k transport base |
| N | 30 | `late_ng_not_rg` | Joint fitter and mesh-verified visibility column |
| C | 30 | `cliff_ge_60deg` | Joint fitter and corrected validated gravity |
| NC | 31 | Both | Both preceding template prerequisites |

Preparation supplies commands with distinct outputs for these cases. The
combined correction enters before one final nonnegative activity clip, with
nonnegative template amplitudes. Keep masks fixed, and evaluate the NG-only
surface throughout the trajectory rather than turning it on at the jump.
Compare earlier pressures and final NG/RG separately; improved aggregate RMS
does not identify a causal cliff production law. Follow with the transition-cut
profile and a separate 45-degree slope sensitivity case, each in new outputs.

The gravity repair and baseline preparation can proceed concurrently. Joint
fitter implementation can use synthetic templates until corrected gravity
exists. Production C/NC fits depend on those products and validated derivatives.
These instructions do not start any expensive fits or gravity calculations.

A single lambda gives one optimizer task: `--jobs 36` does not parallelize its
iterations. Use one worker here. C++ replay may use `OMP_NUM_THREADS=36`; the
gravity backend's own threading must be measured because TBB does not
necessarily honor OpenMP controls. Overlapping jobs require an explicit shared
CPU/RAM budget. Archive the code version, input hashes, commands, full optimizer
diagnostics, per-gauge common-point statistics, and selected C++ replay checks
with the pressure plots in each output directory.
