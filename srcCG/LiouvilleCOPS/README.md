# LiouvilleCOPS: forward model, COPS data, and surface-source fitting

LiouvilleCOPS is a standalone application for predicting ROSINA/COPS nude-gauge
(NG, G0) and ram-gauge (RG, G1) pressures during the 2016-09-30 Rosetta descent,
and for fitting a surface source model to those observations. Its C++ solver
integrates collisionless molecular transport from a triangulated nucleus.
Its Python fitter reuses saved geometric kernels to adjust an H2O activity map,
three surface-temperature parameters, and a global CO2 multiplier.

This description was checked on **2026-10-17** against files inside
LiouvilleCOPS, including the current C++ source, Python scripts, the fitter,
data sources, and the results that were present during the audit. References to legacy names such
as `srcCG`, `shap5`, and `RosinaMeasurements.h` describe identifiers and
provenance retained in these files. The Makefile builds only this application's
`src/*.cpp`; it does not compile or link surrounding application code.

The application has its own implementation, but the supplied configuration
still points to external observation, shape, and SPICE files. Those external
files were not inspected in this review. Their locations below are configuration
values, not confirmation that they are installed or scientifically validated.

Two distinctions are essential:

- **Measurements selected for fitting** are observations and selection flags.
  The multi-source comparison script displays these; it does not fit a model.
- **Fitted/reconstructed pressures** are predictions produced by the fitter or
  a C++ rerun with fitted coefficients. They are displayed by the model plotter.

[RUN.md](RUN.md) covers building and running the application and its tools;
[quick_guide.md](quick_guide.md) is the fitting quick start.

The CPU fitter now saves `fit_restart.npz` by default and accepts
`--restart-from` files from either CPU or CUDA fits. The wire format matches
the CUDA C12 version-1 format, without a CUDA runtime dependency. See
[FIT_RESTART.md](FIT_RESTART.md) for parameter mapping, mesh/template checks,
multi-lambda outputs, and the distinction between a warm start and an
iteration checkpoint. The new fit's prior and selection remain unchanged.

An optional second stage, `tools/fine_tune_activity.py`, freezes a completed
fit and adds nonnegative activity on facets first contributing in successive
late-mission intervals. Each interval has up to three amplitudes: NG-dominated,
RG-dominated, and all newly contributing facets. No harmonic degree is raised.
Earlier cached pressures are protected by exact structural-zero checks, not a
small-response cutoff. Existing temperature, CO2, activity, and data-selection
settings remain frozen. [FINE_TUNING_STEP.md](FINE_TUNING_STEP.md) documents the
algorithm, reusable response cache, overlap/identifiability rules, commands,
limitations, validation, and the array contract for a later GPU implementation.
This is a separate refinement workflow, not a change to ordinary fitter defaults.

For external SHAP8 PLY/BDF products, large node-ID support, and the separate
44M collinear-triangle repair workflow, see [PLY_CONVERSION.md](PLY_CONVERSION.md).
Forward-surface validation and gravity/cliff admission are separate checks;
the repair does not relax the gravity topology policy.
[ROADMAP.md](ROADMAP.md) gives a staged implementation and validation plan
for future work; its proposed features are not claims of current capability.
[amendement_development_plan.md](amendement_development_plan.md) tracks the
localized-activity A-stages and their individual implementation status.
[SCRIPTS.md](SCRIPTS.md) is a script overview. Where older descriptions or
script comments conflict with the implementation and local-file checks described
here, use the qualifications in this README.
[`HARMONIC_BASIS.md`](HARMONIC_BASIS.md) and
[`ACTIVITY_ENHANCEMENT_FILE.md`](ACTIVITY_ENHANCEMENT_FILE.md) are the
normative mathematical and coefficient-file specifications for the optional
higher-degree H2O forward-model component. [`A6_VALIDATION.md`](A6_VALIDATION.md)
defines its implementation-validation scope, tests, fixed thresholds, and
claim limits.

Contents:

- [Directory and script inventory](#1-directory-contents-and-implementation-map)
- [Workflow](#2-how-information-flows-through-the-application)
- [Data and their limitations](#3-data-paths-units-and-limitations)
- [Forward-model physics](#4-forward-model-physics-and-numerical-implementation)
- [Rotation-driven variability and its fitting controls](#47-rotation-driven-source-variability-and-controls-in-the-fit)
- [Fitting procedure](#5-what-is-fitted-and-how)
- [Numerical workflow and diagnostic interpretation](#57-numerical-workflow-and-diagnostic-interpretation)
- [Commands and options](#6-running-the-application)
- [Outputs and file formats](#7-results-and-file-formats)
- [Interpretation of historical fits](#8-what-the-historical-fits-do-and-do-not-establish)
- [Publication and interpretation checklist](#81-publication-and-interpretation-checklist)
- [Verification](#9-verification-performed-for-this-documentation)
- [Limitations and improvement priorities](#10-limitations-and-priorities-for-improvement)

## Generated-output location

New scientific products go to **`<AMPS>/cops_output/`**, where `<AMPS>` is
the AMPS directory containing `srcCG/LiouvilleCOPS`. No machine-specific
absolute output path is built into the application.
Relative output arguments (`--out`, `--save`, `--plots`, `--fit-data`,
digitizer `--out-dir`, the HK output filename, and C++ `OutputDirectory`)
resolve under that root. Absolute output paths are honored only when they are
outside LiouvilleCOPS; a path that resolves inside the application is rejected.
Parent directories are created automatically. Set `COPS_OUTPUT_ROOT` to an absolute directory outside LiouvilleCOPS to
override the root for both Python and C++; rebuild the executable after updating
the sources. The C++ default derives from the executable's application directory;
keep the executable in LiouvilleCOPS. Python derives it from the tools directory.

Generated-product inputs (run directories, fit tables, HK tables, activity
files, and digitized CSVs) prefer the output root, then existing files in the
working/application directory. An absolute path or explicit `./` / `../` input
selects that location directly. Historical local results have now been removed.
Source images, L4 data, configurations, and external scientific
inputs stay in their existing locations; build objects, the executable, and
Python bytecode are still local build/runtime artifacts. Shell redirections
and commands such as `cp` do not use the application's path resolver.

External mission inputs default to
`/data/vtenishe/rosetta_data_end_of_mission_cops`: raw BG files are under
`data_COPS_L2/mtp35_cops_bg`, while the L4 table and source images are under
`sources/`, and the supplied SHAP5 mesh is `SHAP5_Volume.bdf` at the root.
`COPS_DATA_ROOT` controls the Python raw/L4/image readers; the C++ mesh path is
set explicitly by `ShapeModel` in the input file. Override the Python common root with `COPS_DATA_ROOT`, or use
`COPS_L2_ROOT`, `COPS_SOURCES_ROOT`, digitizer `--dir`, and collector
`--sources-dir` for separate locations.

Output names below are relative to this root unless explicitly described as
historical application-local files. See [RUN.md](RUN.md#output-location) for
examples. `tools/cops_paths.py` and `src/Paths.h` implement these conventions;
`tests/test_output_paths.py` checks routing and opt-in end-to-end workflows.

### Provenance, integrity, and result states (W01)

Every new observation-mode C++ run writes `run_manifest.json` after reading
its completed table and binary caches back from disk. The version-1 JSON
schema is `liouville_cops.run_manifest`. It records the exact input file,
effective settings including defaults, executable identity, external input
paths/sizes/SHA-256 hashes, output hashes, kernel format/count/record count,
and independent process/cache/scientific-review states. SHA-256 is streamed
in bounded memory; large meshes and kernel collections are not loaded merely
to hash them. `input.snapshot.in` preserves the exact run input inside the
output directory.

Every new Python fitting operation writes `fit_manifest.json`, schema
`liouville_cops.fit_manifest` version 1. It records the parent run, strict
cache audit, exact selection hash/counts, prior identity, parameter names and
bounds, objective and optimizer settings, Python/library/code identities,
optimizer result, declared output hashes, and these separate states. Stage A5
also records the profile source/reason hashes and counts, legacy coefficient
source, enhancement basis/degree/order/scale/clipping, fixed/free component
status, objective decomposition, data and augmented Jacobian ranks/singular
values, coefficient-file hashes, and a per-result report identity.

| State | Meaning; what it does **not** mean |
|---|---|
| `process` | The command completed or failed; completion does not imply optimizer convergence. |
| `optimizer` | `not_requested`, `converged`, `not_converged`, or `mixed` for a scan; convergence does not establish uniqueness or physical validity. |
| `cache` | The table, elements, and kernels passed structural and correspondence checks. |
| `base_python_cpp_agreement` | The cached Python model reproduced the base C++ aggregate pressures at the configured tolerance. |
| `cpp_prediction_agreement` | A separate C++ rerun reproduced the fitted aggregate pressures; initially `not_run`. |
| `scientific_review` | Always `not_reviewed` unless a later, explicit scientific review process is defined. Passing software tests never changes it. |

New C++ runs and fits refuse an existing output path. This prevents merging a
new experiment with stale kernels or historical results. The boundary rebuild
wrapper is the one controlled exception: it creates a new directory itself,
then explicitly permits the solver to use it while the solver still refuses
any pre-existing scientific product in that directory.

Before fitting, `tools/cops_io.py` rejects malformed rows, duplicate point IDs,
missing/extra kernels, unsupported versions, truncated or appended
binary data, table/kernel point or geometry disagreement, invalid record
counts, duplicate/out-of-range face indices, non-finite cache values, and
manifest/output hash disagreement. A base created before W01 has no manifest
and is refused by default. `--allow-legacy-base` is an explicit compatibility
choice: the same structural audit is performed, and the missing provenance is
recorded in the fit manifest rather than reconstructed by assumption.

Use the read-only cache validator with:

```bash
python3 tools/validate_cops_run.py output_fit_base
```

After running a fit's generated C++ input, register numerical reproduction:

```bash
./liouville_cops output_fit_new/cops_fit.in
python3 tools/validate_cops_run.py \
  --fit output_fit_new --verification output_fit_new_verify
```

The second command verifies declared fit outputs, the verification run
manifest, parent dependency identities when the parent is manifested, exact
point/UTC correspondence, and aggregate NG/RG pressures. It changes only
`cpp_prediction_agreement`; it does not certify calibration, model adequacy,
convergence, uniqueness, uncertainty, or publication readiness.

## 1. Directory contents and implementation map

### Top-level files and directories

| Path | Purpose and current contents |
|---|---|
| `Makefile` | Builds `liouville_cops` from local C++ sources using CSPICE and OpenMP. |
| `liouville_cops` | Existing compiled executable. Rebuild after changing source or build settings. |
| `build/` | Compiled object files; not scientific input. |
| `cops.in` | Example full forward run, currently writes `output/`. |
| `fit_base.in` | Forward run that writes `output_fit_base/` with kernels and surface elements for fitting. |
| `m004.in` | Geometry-only M004 DSK run that writes `output_M004/COPS_model.dat` for the collector's exact nearest-surface altitude; gas transport and fitting kernels are disabled. |
| `selftest.in` | Analytic/numerical verification using a generated sphere. |
| `HARMONIC_BASIS.md`, `ACTIVITY_ENHANCEMENT_FILE.md` | Normative Stage-A2 basis and Stage-A3 coefficient-file/model contracts. |
| `A6_VALIDATION.md` | Stage-A6 validation scope, requirement-to-test map, predeclared tolerances, commands, and limits on what passing tests establish. |
| `examples/activity_enhancement_zero_l5.txt` | Complete all-zero degree-5 forward-model regression example for `builtin:shap5` in `67P/C-G_CK`. |
| `src/` | Entire C++ implementation and embedded numerical tables. |
| `tools/` | Data preparation, fitting support, and plotting scripts; inventory below. |
| `tests/test_rebuild_cops_boundary.py` | Standard-library tests for boundary rebuilding, with an opt-in local C++/mission-data integration check. |
| `tests/test_a2_activity_basis.py`, `tests/test_a3_activity_enhancement.py`, `tests/test_a4_enhancement_fit.py`, `tests/test_a5_reporting.py`, `tests/test_a6_validation.py` | Cross-language harmonic tests, strict optional-forward-model tests, fitter parameter/Jacobian/regularization tests, report identity/table/plot tests, and deterministic Stage-A6 synthetic recovery/observability tests. |
| `data/` | Input-data directory; generated tables and plots have been removed. |
| External `<COPS_DATA_ROOT>/SHAP5_Volume.bdf` | SHAP5 NASTRAN surface mesh used by `cops.in` and `fit_base.in`; the supplied inputs currently use its absolute path. |
| External `<COPS_DATA_ROOT>/sources/` | L4 table and labels plus three source images; default root is `/data/vtenishe/rosetta_data_end_of_mission_cops`. |
| `README.md`, `RUN.md`, `SCRIPTS.md`, `quick_guide.md` | Detailed description, operating guide, script overview, and fitting quick start. |
| `ROADMAP.md` | Future implementation work packages, dependencies, tests, acceptance criteria, and scientific decision gates. Not an implemented feature list. |
| `<AMPS>/cops_output/` | Destination for newly generated scientific products; outside this application. |

The following historical inventory records what was present during the model
audit. These directories and logs were removed from LiouvilleCOPS after their
documented properties were recorded; they are not currently available files:

| Removed output directory/group | What the audited files contained |
|---|---|
| `output_fit_base/` | Base `COPS_model.dat`, `elements.bin`, and 1,150 point-kernel files; about 1.3 GB on disk. |
| `output_fit_base_fit/` | Historical single-fit parameters, activity file, generated rerun input, fitted pressures, and data-selection plots. |
| `output_fit_base_fit_verify/` | C++ verification run, including elements and kernels; about 1.3 GB. |
| `output_fit_base_lcurve/` | `lcurve.txt`, plot, selection files, and 15 `lambda=.../` fit directories. |
| `output_fit_base_selection/` | Data-only selection table and plots, with 930 NG and 1,146 RG points. |
| `output_M001/`, `output_M004/` | Retained pressure/geometry tables for the DSK shape models. |
| `output_check_M001/`, `output_check_M004/` | Additional shape-check output tables. |
| `output_rot_none/`, `output_rot_exact/` | Saved static/flight-time-rotation comparison tables. |
| `output_esa_ratt01_T6/`, `output_esa_ratt04/` | Saved attitude-kernel comparison tables. |
| `output_ng_full/`, `output_ng_hemisphere/`, `output_ng_ik_fov/`, `output_ng_srccg_shield/` | NG acceptance comparison tables. |
| `output_flat/`, `output_rough/`, `output_flat_offset/`, `output_test/` | Other boundary-condition, geometry, or test runs; inspect their headers before interpreting them. |

Some removed historical headers referenced temporary input files that were not
included. Their complete configuration was not recoverable from directory names
alone; removed generated fit inputs also contained obsolete absolute paths.

### C++ modules

| Local source | Responsibility |
|---|---|
| `main.cpp` | Driver, verification tests, data selection, per-time evaluation, text and binary output. |
| `Input.cpp/.h` | Parses `Key = Value`, comments, repeated keys, options, and typed values. |
| `Observations.cpp/.h` | Parses observation times and NG/RG arrays from the configured text header at runtime. |
| `Mesh.cpp/.h` | NASTRAN/DSK readers, triangle geometry, topology/orientation checks, distances, test sphere. |
| `BVH.cpp/.h` | Bounding-volume hierarchy for visibility and shadow ray intersections; brute-force reference traversal. |
| `Geometry.cpp/.h` | SPICE loading, spacecraft/Sun states, body rotation, gauge pointing and field of view. |
| `Source.cpp/.h` | Activity basis and coefficients, illumination law, shadow state, temperatures, `elements.bin`. |
| `ActivityEnhancementBasis.cpp/.h` | Versioned normalized real spherical-harmonic evaluator through degree 8 specified in `HARMONIC_BASIS.md`. |
| `ActivityEnhancement.cpp/.h` | Strict version-1 harmonic and exact per-facet replay parsers, mesh/frame/baseline compatibility checks, residual evaluation, and final clipping for the optional H2O forward-model enhancement. |
| `SrcCGTables.h` | Locally embedded activity coefficient arrays and tabulated temperatures. |
| `Gauges.cpp/.h` | NG angular response, RG acceptance, species sensitivity, pressure conversion, housekeeping temperature interpolation. |
| `Liouville.cpp/.h` | Static geometric kernels and speed-integrated flight-time-rotation solver. |
| `Provenance.cpp/.h` | Streaming SHA-256, run manifests, atomic JSON output, and post-write table/cache validation. |
| `Constants.h`, `Vec3.h` | Species/constants and vector operations. |

### Python scripts

| Script | Purpose | Main inputs | Main outputs |
|---|---|---|---|
| `tools/cops_hk_table.py` | Extracts pressures and housekeeping from raw L2 BG files. | Directory containing `bg_*.tab`. | User-named housekeeping table. |
| `tools/rebuild_cops_boundary.py` | Regenerates a complete static fitting base for a selected shape and `flat` or `rough` surface emission; see section 6.6. | Original observation input, existing compiled solver, selected mesh, configured external data. | New input snapshot, table, elements, kernels, log, and completion metadata in a fresh directory. |
| `tools/collect_cops_sources.py` | Combines local measurement sources, compares them, and constructs a candidate fit-selection table. | HK table, L4 table, digitized CSVs, SHAP5/M004 runs, base fit run. | Combined JSON, comparison plots, fit-selection table and plot, optional actual-fit profile plot. |
| `tools/digitize_cops_figures.py` | Extracts curves from the three supplied source images. | JPG files in `<COPS_DATA_ROOT>/sources/` or `--dir`. | Five digitized CSVs; optional diagnostic PNGs. |
| `tools/plot_cops_model.py` | Compares predicted NG/RG pressures and reference observations; strictly verifies manifested fit/report identities and labels legacy inputs. | `COPS_model.dat` files. | Interactive four-panel plot, or `--save` image whose labels identify profile/degree/lambda. |
| `tools/plot_cops_model_altitude.py` | Wrapper around the preceding script with altitude as the default horizontal axis. | Same. | Same. |
| `tools/plot_lcurve.py` | Displays fit quality versus regularization, temperatures, and activity diagnostics. | `lcurve.txt`. | Interactive plot or `--save` image. |
| `tools/fit_cops.py` | Fits either the legacy degree-0--4 source/temperature/CO2 parameters or the controlled degree-5--8 H2O residual using static kernels. | Base run, kernels, elements, local coefficient tables, original input, HK table, optional replacement fit data or named selection profile. | Parameters, fitted activity, fitted pressures, rerun input, exact data selection, self-identifying A5 tables/plots; enhancement exchange file when selected; optionally an L-curve. |
| `tools/run_fit_profile.py` | Validates and runs one or all checked-in selection profiles with safe, distinct output directories. | Static fitting base and `fit_profiles/*.json`. | The normal fitter products for each requested profile. |
| `tools/run_a6_validation.py` | Runs the complete unit/synthetic gate, optionally enables local C++/mission-data integration tests, and runs the C++ self-test. | Compiled executable; integration mode also uses configured external shape, observation, SPICE, and source inputs. | Terminal evidence and an optional immutable JSON report; test run products stay in temporary directories. |
| `tools/benchmark_fit_performance.py` | Measures the Stage-A7 production fitter startup, basis, residual/Jacobian, memory, independent-fit scaling, and result-package size. | A current manifested static base, matching HK table, profile/degrees/job candidates, and a declared campaign size. | Terminal resource tables plus optional new JSON/Markdown reports; temporary output probes are deleted. |
| `tools/convert_dsk_to_bdf.py` | Converts one-segment type-2 PDS DSKs into the metre-unit NASTRAN surface format used by the original SHAP5 file, with strict topology repair/audit. | DSKs, mission frame kernel, and NAIF `dskbrief`/`dskexp`. | New `.bdf` files plus per-file conversion manifests; existing destinations are refused. |
| `tools/cops_selection.py` | Shared selection/profile parser used by the fitter and profile runner; not a standalone command. | Model rows, base intervals, saturation mask, and a profile JSON. | Selected flags, reason codes, counts, and a canonical reason hash. |
| `tools/activity_basis.py` | Reference Python implementation of the versioned normalized real spherical-harmonic basis used by tests and the Stage-A4 fitter; not a standalone command. | Cartesian directions or angular coordinates supplied by callers. | Basis arrays in the exact ordering defined by `HARMONIC_BASIS.md`; no files. |
| `tools/cops_fit_reporting.py` | Stage-A5 report writer imported by the fitter; not a standalone command. | Final facet activity, pressures, selection/reason masks, and immutable report-identity payload. | Per-facet/degree/region/residual tables, one activity-map figure, three pressure/residual figures, and `fit_report.json`. |
| `tools/benchmark_activity_basis.py` | Measures evaluation time and memory for the Stage-A2 basis on centroids from an actual `elements.bin` without changing that cache. | A version-1 `elements.bin`, degree range, repeat count. | Console report and optional JSON report. |
| `tools/build_physical_activity_templates.py` | Constructs mesh-bound late-NG-only visibility and effective-slope cliff hypotheses without fitting their amplitudes. | Manifested base, explicit UTC interval, homogeneous-gravity surface/density, optional resolution-check surface, and SPICE helper. | Strict `.npz` template bundle, reusable per-facet gravity cache, JSON metadata, and body-frame diagnostic map. |
| `tools/activity_templates.py` | Strict template-bundle writer/reader imported by the builder and fitter; not a standalone command. | Per-facet fields plus the exact active `elements.bin` centroid/area sequence. | Versioned bundle or validated in-memory columns tied to one mesh identity. |
| `tools/gravity_cache.py` | Strict reusable gravity-cache I/O imported by the physical-template builder; not a standalone command. | Exact modeled centroid/normal/area ordering, gravity shape/hash, density, offset, self gravity, and reference SPICE rotation. | Validated self, centrifugal/effective-gravity, and slope arrays; epoch-dependent fields are cheaply re-derived when rotation differs. |
| `tools/spice_angular_velocity.cpp` | Small compiled helper using the same `Geometry.cpp` and SPICE configuration as the forward solver to obtain body-frame angular velocity. | Observation input and one or more UTC epochs. | Machine-readable angular-velocity rows on stdout. |
| `tools/cops_io.py` | Shared strict product parser, cache validator, streaming hashes, atomic manifests, and selection identities; imported by other tools. | Model tables, manifests, elements and static kernels. | Validated in-memory metadata; no standalone output. |
| `tools/validate_cops_run.py` | Validates a static base or records fitted-pressure agreement with a manifested C++ verification run. | Base run, or fit plus verification run. | Diagnostic report; in verification mode, an atomic update of the fit manifest's C++-agreement state. |

Both the collector and fitter are present as ordinary working files. The
fitter's numerical implementation originated in former RCS revision 1.9;
the current file additionally routes products through `cops_paths.py`. The
RCS directories were removed only after all archived files were confirmed to
have working counterparts.

## 2. How information flows through the application

```text
Raw COPS L2 BG files
  -> cops_hk_table.py
  -> data/COPS_HK_20160930.txt
       -> forward solver: RG sensor temperature
       -> fitter: NG saturation identification
       -> source collector: L2 pressure/housekeeping dataset

ObservationFile + shape + SPICE + source/gauge settings
  -> liouville_cops fit_base.in
  -> output_fit_base/COPS_model.dat + elements.bin + kernels/
       -> fit_cops.py + optional --data selection table
       -> ActivityFit.txt + fitted COPS_model.dat + cops_fit.in
            -> liouville_cops cops_fit.in: independent C++ evaluation

HK + L4 + digitized figures + saved geometry/reference tables
  -> collect_cops_sources.py
  -> comparison JSON/plots + COPS_fit_data.txt
       -> optional fit_cops.py --data COPS_fit_data.txt
```

The collector does **not** merge all plotted sources into one fitting dataset.
Its fit table takes pressures from the chosen base run's `refNG/refRG`
columns. L4 and digitized curves are comparison datasets unless the user
explicitly constructs a different compatible `--data` table.

Changing the raw L2 directory, rebuilding the HK table, or regenerating the
collector plots does not automatically replace the C++ observation pressure
arrays or an existing base run's reference pressures.

## 3. Data, paths, units, and limitations

### 3.1 Which files are used

All example commands below assume the current directory is the LiouvilleCOPS root.

| Input | Location or setting | Actual use |
|---|---|---|
| Raw L2 BG files | `/data/vtenishe/rosetta_data_end_of_mission_cops/data_COPS_L2/mtp35_cops_bg` by default. | Input to `cops_hk_table.py`; no raw BG reader in the C++ solver. |
| Processed L2 table | `<AMPS>/cops_output/data/COPS_HK_20160930.txt`. | Generate with `cops_hk_table.py`; contains pressures, temperature, currents, ranges, flags, and quality. |
| Observation table | `ObservationFile` in `cops.in` and `fit_base.in`, currently `/home/vtenishe/ROS-END-OF-MISSION/AMPS/srcCG/RosinaMeasurements.h`. | Runtime data input. The file is not compiled into this application and is not bundled inside this directory. |
| Base observations for Python fitting | `<AMPS>/cops_output/output_fit_base/COPS_model.dat`: `UTC`, `refNG`, `refRG`, `T_RG`, `inside`, `useNG/useRG`. | Generate with `liouville_cops fit_base.in`; these saved columns are read by the fitter. |
| Optional edited fit data | `<AMPS>/cops_output/data/COPS_fit_data.txt` or `fit_cops.py --data FILE`. | Generate with the collector; overrides pressures and selection at existing base-run points. |
| L4 density | `<COPS_DATA_ROOT>/sources/COPS_L4_MTP35.ASC`; accompanying `.LBL` and `COPS_TS_TABLE.FMT`. | NG density and uncertainty for comparison. |
| Digitized source figures | Three JPG inputs in `<COPS_DATA_ROOT>/sources/`; generated CSVs go under `<AMPS>/cops_output/data/sources/`. | Generate with `digitize_cops_figures.py` for comparisons to ESA-blog NG and EGU-slide NG/RG curves. |
| SHAP5 shape mesh | `/data/vtenishe/rosetta_data_end_of_mission_cops/SHAP5_Volume.bdf`, set by `ShapeModel` in `cops.in` and `fit_base.in`. | NASTRAN surface used for the standard forward and fitting-base runs; it is external and not included in this Git directory. |
| SPICE and alternate shapes | `KernelDirectory`, repeated `Kernel` keys, and optionally `ShapeModel` in the input files. | Forward geometry; paths remain external configuration. |

The raw-data location is specified when regenerating the table:

```bash
python3 tools/cops_hk_table.py \
  /data/vtenishe/rosetta_data_end_of_mission_cops/data_COPS_L2/mtp35_cops_bg \
  data/COPS_HK_20160930.txt
```

This writes the processed table under `<AMPS>/cops_output`. The script reads every matching
`bg_*.tab` in the supplied directory; it does not select a date range or
download missing files. Use a directory containing the intended intervals.

### 3.2 Observation format and timestamps

`ReadObservations()` expects a text file containing three C/C++ array
initializers named:

```text
ObservationTime
NudeGaugeReferenceData
RamGaugeReferenceData
```

It strips comments, extracts quoted UTC strings and comma-separated numerical
values, and requires equal array lengths. A generic CSV or the HK table cannot
be substituted directly without changing this reader or converting the input.

The audited, subsequently removed base run had 1,150 points, from `2016-09-30T0:00:04` through
`2016-09-30T10:39:26`. Some hour fields are not zero-padded. The Python
tools account for this in their timestamp parsing.

The HK table has 1,155 rows. Its row time is the **midpoint** of the product's
`START_TIME..STOP_TIME` interval; it also preserves both endpoints. The
reference times mostly correspond to rounded interval **endpoints**. Do not
compare the two tables by row number or by their displayed midpoint times.

A local-file audit matched each reference UTC to the nearest HK stop time,
accepting a match only when the difference was less than one second:

- 1,148 of the 1,150 reference times matched.
- `09:13:04` and `09:14:04` had no match within that tolerance.
- Among matched, positive values, maximum relative differences from
  `100 * p_HK_mbar` were about **0.478% NG** and **0.397% RG**.
  These small differences are consistent with the lower precision visible
  in the saved reference pressures, but the conversion history is not established here.
- At `01:30:04`, the reference NG value is zero whereas the matched HK row
  ending `01:30:03.989` contains `7.33790e-10 mbar`.
- Matched RG rows had no disagreement in positive-versus-missing status.

Thus the local reference and HK tables are closely related but **not identical**.
The origin of the two unmatched times and the NG missing-value discrepancy
remains unresolved. Older comments claiming exact identity are too strong.
This comparison used the local saved reference columns, not the external
observation header or raw BG archive.

The local HK table has one midpoint dated September 29 and 1,154 dated
September 30. Older descriptions of “three previous-day files” were not
established by this local table and should not be used as an inventory claim.

### 3.3 Units and gauge conventions

- Raw HK pressures and digitized CSV pressures are in **mbar**.
- Forward output and the fit table use **Pa**: `1 mbar = 100 Pa`.
- `ReferenceDataUnit = mbar` converts the observation table on reading;
  the C++ fallback if that key is absent is `Pa`.
- The tools label pressures as N2-equivalent instrument pressures at
  approximately 20 C. They are not direct species partial pressures.
- Zero/nonpositive reference pressures are treated as missing observations,
  not physical measurements of zero gas density.
- In the combined JSON, missing values become `null`.
- In `fit_data_used.txt`, `-1` means that gauge's point was not selected.

The collector calculates NG equivalent density as `p/(k * 293 K)`.
Its RG equivalent flux uses
`p / sqrt(2*pi*m_H2O*k*T_RG)`, with the H2O mass even in the column
named `PhiRG_N2eq_m2s`. Its `_H2O` density/flux columns multiply by
0.893. These are specified conversion conventions, not a measured
H2O/CO2 decomposition. The forward model instead sums species responses
with their own masses and sensitivity factors.

The sensor temperature comes from `ROSINA_COPS_SENSOR_TEMP`, converted
from degrees Celsius by adding 273.15. The model uses it as `T_RG`.
The available implementation and processed files do not establish that
this particular telemetry sensor measures precisely the thermodynamic
temperature required by the ram-gauge pressure relation. Treat that mapping
as a modeling assumption requiring instrument documentation/calibration.

### 3.4 Missing and saturated NG data

In the base reference table, NG is nonpositive at 216 points:

| UTC interval in the reference table | Number of missing NG points |
|---|---:|
| 01:12:04–01:30:04 | 19 |
| 03:12:04–03:29:04 | 18 |
| 05:12:04–05:29:04 | 18 |
| 07:12:04–07:28:04 | 17 |
| 09:13:04–09:29:04 | 17 |
| 10:24:42–10:33:06 | 127 |

The long final gap ends before the first subsequent positive point at
10:33:14. The input's `ExcludeNG` interval is
`10:24:40..10:40:00`, covering the recorded HIGH-to-LOW filament range
change and the remaining descent.

The saturation detector marks positive readings equal to the maximum NG
pressure in the HK table's HIGH filament range. It matches those records
to reference times through the stop time with a one-second tolerance.
For the present data this identifies 13 reference points from
10:23:50 through 10:24:38. The corresponding BG intervals span
10:23:46.168 through 10:24:38.251, with pressure
`2.26330e-8 mbar` and zero recorded filament current.

This is a dataset-specific plateau detector. It does not use the zero
current as its decision criterion, fit censored observations, or implement
a general detector-saturation calibration. Both the fitter and collector
include those points by default. See the explicit exclusion options below.

Neither fitter nor collector automatically applies a general
`DATA_QUALITY_ID` filter or uncertainty weighting from housekeeping.
The quality field is retained by preprocessing.

### 3.5 Other sources and their evidential status

The combined JSON has seven datasets:

| Dataset ID | Source and interpretation |
|---|---|
| `PSA_L2_BG` | Processed BG housekeeping, with measured NG/RG pressures and sensor information. |
| `srcCG_table` | Legacy dataset ID for the reference columns from the selected SHAP5 run; the collector does not read an external source-code tree. |
| `fit_data` | Reference pressures plus candidate selection flags derived from the fit run and collector options. |
| `PSA_L4` | NG density and uncertainty from the local L4 ASCII table. |
| `ESA_blog_2016` | Digitized NG interpretation of the supplied ESA-blog figure. |
| `EGU2019_s02` | Digitized NG/RG versus the altitude shown in slide 2. |
| `EGU2019_s16` | Digitized NG/RG versus the altitude shown in slide 16. |

The local L4 file contains 901 rows dated September 30, from 00:00:58 to
10:36:41. Its density uncertainties are approximately 20%, with numerical
rounding. The label's overall `STOP_TIME` is 10:39:26, so that metadata
must not be mistaken for the time of the last retained density row.
The collector plots L4 as `n*k*293 K`; this is a density-derived pressure
equivalent, not the original NG gauge reading. The source metadata describe
composition correction using DFMS. L4 uncertainties do not enter the fitter's
default objective.

The digitizer uses fixed pixel-axis calibrations, color masks/median curve
positions for the EGU plots, and black-dot/component extraction for the blog
figure. It is specific to these images. EGU curves have altitude coordinates,
not observation times. Their altitude definition is not established by the
files here, and the collector preserves their plotted altitude rather than
silently converting it to SHAP5 or M004.

The collector's comments and JSON contain historical comparison statements:
a roughly 230-second blog time offset, different low-altitude EGU curves,
and several numerical agreement estimates. These statements are embedded
text, not tests or fitted corrections performed on each execution. They
were not independently rederived in this review. In particular, the script
does not shift the blog curve by 230 seconds or fill the fit's archive gap
with digitized points.

The same distinction applies to the hardcoded touchdown time and
`not_included` literature/archive notes in the JSON. They record prior
interpretations/searches; they are not evidence that the present program
has recomputed touchdown or checked current archive availability.

## 4. Forward-model physics and numerical implementation

### 4.1 Scope

The solver treats emitted neutral molecules as collisionless. It includes
H2O and CO2 only (`Const::nSpecies = 2`). CO coefficients are retained
in activity files for format compatibility but CO does not contribute to
the calculated pressures.

There is no intermolecular collision operator, gravity, photochemical loss,
distributed coma source, or gas energy equation. Surface temperature is a
prescribed parametrization. Agreement with pressure measurements alone does
not demonstrate that the collisionless approximation is valid at every
altitude; no mean-free-path/Knudsen-number test is implemented.

### 4.2 Shape, visibility, and SPICE

`ShapeModelFormat = nastran` reads GRID/GRID* nodes and CTRIA3 surface
triangles; CTETRA volume cells are ignored. NASTRAN coordinates are used
numerically as supplied, with no unit conversion: supply meters.
`dsk` accepts a single-segment SPICE type-2 plate model, converts its
kilometer coordinates to meters, and checks its frame against `BodyFrame`.

The mesh checks directed-edge consistency, closure by default, degenerate
triangles, and signed-volume orientation. It flips all triangle windings
if required for positive volume. Visibility uses the first triangle hit by
a ray from the gauge. A ray-count parity test marks positions inside the body.

At every observation time, SPICE supplies:

- Geometric spacecraft state relative to the nucleus, in `BodyFrame`,
  using aberration option `NONE`.
- Apparent Sun direction and heliocentric distance using `LT+S`.
- Gauge boresights and FOV corners, obtained from the instrument kernel
  and rotated into the body frame.
- Body angular velocity from the body-to-J2000 state transformation.

The supplied frame is `67P/C-G_CK`; the configured instrument IDs are
-226331 (NG) and -226332 (RG). The input explicitly lists its kernels
rather than loading a meta-kernel automatically.

Both gauges use one position: spacecraft reference position plus
`GaugeOffsetSC` rotated from the spacecraft frame. The supplied offset is
zero. The code does not separately locate the two sensor heads or add the
time derivative of this positional offset to spacecraft velocity.

Altitude in the C++ table is the **unsigned minimum distance to the
triangulated surface**, with `inside` as a separate flag. The collector
and fit-data exporter apply a negative sign for inside points. Altitude
is not distance from the nucleus center.

The removed, audited 100k SHAP5 base had four inside points, indices
1146–1149. The 788k SHAP5 base generated for the resolution experiment has no
inside points and puts the final gauge position 0.340 m outside. The audited
M001 and M004 tables had no inside points; their final altitudes were about
5.989 m and 5.928 m. Inside status is therefore a property of each saved
geometry, not a universal data flag or shape-model correction. Activity
coefficients are evaluated at each new shape's face centroids; they are not
independently recalibrated when the shape is changed.

### 4.3 Emitted source and temperature

For face e and species s, the implemented source is

```text
A_s(e) = max(0, sum_j a_s,j Y_j(e))
F_s(e,t) = q_s(e,t) A_s(e)                         [m^-2 s^-1]
```

The 25 basis functions are explicitly defined in `Source.cpp` and
evaluated at the face centroid's colatitude and longitude. Use these exact
functions and their ordering when interpreting or replacing coefficients.
They are not supplied as a normalized general-purpose harmonics library;
indices 11 and 13, for example, explicitly contain
`5*(cos(colatitude)^2 - 1)*sin(colatitude)*sin/cos(longitude)`.
Do not substitute a library basis and reuse the coefficient arrays without
an explicit basis transformation.

#### Optional independently normalized degree-5--8 H2O residual (Stages A2--A5)

The application contains a second harmonic evaluator for the localized H2O
activity amendment. Its versioned identifier is
`liouville_cops.real_orthonormal_cs_v1`; the normative equations are in
[`HARMONIC_BASIS.md`](HARMONIC_BASIS.md). The authoritative implementations
are `tools/activity_basis.py` and
`src/ActivityEnhancementBasis.{h,cpp}`.

This is a fully orthonormal **real** basis with the Condon--Shortley phase.
For `m>0`, positive signed order is the cosine branch and negative signed
order is the sine branch. Coefficients are ordered by increasing degree and
then `m=-l..+l`. The supported enhancement range starts at degree 5 and contains
11, 24, 39, or 56 coefficients through maximum degree 5, 6, 7, or 8.

The evaluator uses mesh/body-frame Cartesian centroid directions:
`theta=acos(z/r)` and `phi=atan2(y,x)`. It defines `phi=0` at the poles, where
every nonzero-order mode vanishes. The mesh axes are interpreted as the input
`BodyFrame`; the supplied manifested base records `67P/C-G_CK`, while the
NASTRAN mesh itself has no independent frame metadata. Nonfinite coordinates,
the zero vector, invalid modes, and degrees outside 0--8 are rejected.

The two implementations use the same explicit associated-Legendre recurrence
and log-gamma normalization rather than relying on differing library defaults.
Tests cover normalization, phase, ordering, poles, axes, near-pole and generic
points, scale invariance, invalid inputs, and every mode through degree 8. On
the deterministic test set, the measured maximum Python/C++ absolute
difference is `1.33e-15`; numerical quadrature's maximum deviation from the
81-mode identity Gram matrix is `1.47e-14`. The enforced cross-language limit
is `2e-14`, far tighter than the later pressure-agreement target of `1e-6`.

An actual-cache benchmark on the current 100,000-face manifested SHAP5 base
constructed the 100,000-by-56 degree-5--8 table in a median `0.0987 s` over
five samples. The returned table is exactly 44,800,000 bytes (42.725 MiB); the
measured process maximum-RSS increase was 115,688 KiB, including temporary
arrays and allocator effects. These values describe that machine/run, not a
portable performance guarantee.

Stage A3 connects this evaluator to `cSourceModel` as a default-off, H2O-only
residual. For a mesh centroid `x`, the solver now supports

```text
DeltaA(x) = ActivityEnhancementScale * sum b_lm Y_lm(x)
A_H2O(x) = max(0, A0_linear(x) + DeltaA(x)).
```

The legacy field is deliberately kept linear until the residual has been
added, and clipping occurs exactly once afterward. CO2 and CO cannot receive
enhancement coefficients. `ActivityEnhancementModel = none` with degree 4 is
the default. `real_spherical_harmonics_v1` reads the strict format in
[`ACTIVITY_ENHANCEMENT_FILE.md`](ACTIVITY_ENHANCEMENT_FILE.md); it checks the
basis ID, H2O species, body frame, complete degree/order set, scale, clipping
rule, and identity of the active legacy map before propagation.

The legacy 25 basis values and version-1 `elements.bin` layout remain
unchanged. When enhancement is enabled, the existing `activity_H2O` value in
that cache is the final applied (clipped) activity, and the separate
`ActivityEnhancementFacets.dat` preserves legacy-linear, residual,
total-linear, and final values for interpretation. When disabled, that new
diagnostic is not written and old table/binary schemas are unchanged.

Stage A4 connects the same basis to `tools/fit_cops.py`. It is intentionally a
controlled residual-only mode: requesting degree 5--8 currently requires all
of `--fix-legacy-activity`, `--fix-temperature`, and `--fix-kappa`. Thus the
new free vector contains only `b_lm`; the 25 legacy H2O coefficients, the three
analytical-temperature values, and `kappa=1` are copied exactly from the base.
Degree 4 is an explicit zero-parameter control when the same three fix options
are supplied. The counts through degrees 5, 6, 7, and 8 are respectively 11,
24, 39, and 56.

The fitter constructs the normalized facet table once, before optimization,
and shares it across lambda workers. The conditioning scale is not free and
cannot be chosen to improve a result. It is computed from the base mesh as

```text
A_scale = sqrt(sum_e area_e A0_linear(e)^2 / sum_e area_e).
```

`--enhancement-scale VALUE` is only an assertion: the supplied value must
equal that derived value to relative tolerance `1e-12`. Coefficients start at
zero and are unbounded. For coefficient regularization, the implemented
dimensionless deviation is

```text
D^2 = sum_e area_e DeltaA(e)^2 / sum_e area_e A0_linear(e)^2,
DeltaA(e) = A_scale sum_(l=5..L,m) b_lm Y_lm(e).
```

This penalty acts on the **unclipped residual** and pulls it toward zero over
all facets, including facets not seen by the fit trajectory. The data term
uses the once-clipped total activity. Its analytic derivative follows the
declared subgradient: it is active only where
`A0_linear + DeltaA > 0`, and is zero at and below the clipping boundary.
Central finite differences are valid only while neither perturbation crosses
that boundary; one-sided checks must be interpreted with the same convention
at a kink.

An enhanced fit writes the unchanged three-by-25 `ActivityFit.txt`, plus the
strict C++ exchange file `ActivityEnhancementFit.txt` and the Python-side
`ActivityEnhancementFacets.dat` with centroid, area, visibility, legacy,
residual, total-linear, and clipped fields. The exchange file hashes the exact
generated `ActivityFit.txt` as its legacy identity. `cops_fit.in` selects both
files so `./liouville_cops <fit>/cops_fit.in` independently recomputes the
result with Stage A3. The manifest records the basis/specification identities,
scale definition/value, degrees, bounds, clipping derivative, regularization,
and warns when `lambda=0` leaves activity unregularized. A base run that
already has an enhancement enabled is rejected; fit from an unenhanced,
manifested static-kernel base.

#### Self-identifying diagnostics and plot integration (Stage A5)

For each optimizer result, the fitter constructs one canonical report identity
from the effective profile and its file hash, selected-point and reason hashes,
legacy-activity hash, complete fitted parameter-vector and activity hashes,
final temperature and CO2 scale, activity basis/degrees, derived enhancement
scale, regularizer identifier, and lambda. This
is schema `liouville_cops.fit_report` version 1. The full payload and identity
are stored in `fit_report.json`; the identity is repeated in the
`COPS_model.dat`/diagnostic-table headers, PNG text metadata, and matching
`fit_manifest.json` optimizer record. The manifest separately hashes every
generated file, including the exact `ActivityFit.txt` and optional
`ActivityEnhancementFit.txt`. Consequently two results with different data
selection, coefficients, basis, degree, or lambda cannot share a report
identity.

The following A5 products are written for both legacy and residual fits:

| Product | Interpretation |
|---|---|
| `activity_facets.dat` | Facet centroid/area, whether the facet contributes to any selected-data kernel, legacy linear activity, fitted residual, total linear field, and once-clipped applied activity. |
| `activity_degree_power.dat` | Coefficient and area-weighted residual power by degree. `spherical_field_power` is supplied only for the normalized enhancement basis; it is `nan` for the unnormalized legacy basis. |
| `activity_region_diagnostics.dat` | Area fractions and area-weighted baseline/residual/total moments for all, seen, and unseen facets, including the fraction where the total linear field is clipped. Active fraction is one minus clipped fraction. |
| `activity_maps.png` | Body-frame longitude/latitude bin averages of baseline, residual, and clipped total. These are area-weighted mesh-bin summaries, not an equal-area or shape-unwrapped map. |
| `pressure_residuals.dat` | Every NG/RG observation and prediction, selected/saturated flags, reason code, `ln(model/observed)`, and its equal-gauge objective value divided by `sqrt(N_gauge)`; residual fields exist only where the datum entered the objective. |
| `pressure_comparison_time.png`, `pressure_comparison_altitude.png` | Model and observations with selected points filled, excluded points hollow, saturated NG crossed, and configured instrument-transition intervals shaded. |
| `pressure_residuals.png` | Selected-data log residuals versus UTC and altitude, with saturated NG identified. |

“Seen” is a numerical contribution flag, not a proof of spatial
identifiability: it means that the facet occurs with a nonzero response in at
least one kernel loaded for the union of selected NG/RG points. A seen facet
can still be weakly constrained or correlated with many others. Likewise,
degree power measures fitted field structure; it is not evidence that a
localized physical source is unique. Interpret it together with data and
augmented Jacobian rank/singular values, neighboring lambdas/degrees/profiles,
unseen-region behavior, clipping, and independent C++ reproduction.

`tools/plot_cops_model.py` treats schema-v2 Python fits strictly. Before
plotting, it verifies the table hash declared by `fit_manifest.json`, the
`fit_report.json` hash, and agreement of report, profile, selection, and basis
identities across all three files. A missing or incompatible fit manifest is
an error. Check-box labels and the reference-panel title show the effective
profile, enhancement degree, lambda, and short report identity, so comparisons
do not silently mix configurations. Historical unmanifested forward tables
remain readable but are visibly marked `legacy`; that label is a provenance
warning, not verification.

Manifested C++ forward/verification tables receive the same table-byte check
against `run_manifest.json`. Their labels show basis, maximum degree, lambda
when applicable, and a short enhancement-coefficient or legacy-activity
identity. A Python fit and its C++ replay therefore do not appear as anonymous
curves with potentially different activity inputs.

#### Unit, synthetic, and integration validation (Stage A6)

Stage A6 is a reproducible implementation-consistency gate, driven by
`tools/run_a6_validation.py`. The normal gate discovers the complete Python
suite and then runs `liouville_cops selftest.in`. The opt-in integration gate
also enables tests that execute small real C++ bases, zero-enhancement
regression, a Python degree-5 fit, its generated C++ replay, manifest
validation, plotting, selection profiles, boundary rebuilding, and external
source handling:

```bash
python3 tools/run_a6_validation.py --threads 2
python3 tools/run_a6_validation.py --integration --threads 2 \
  --report a6_validation/validation_YYYYMMDD.json
```

The deterministic synthetic tests use the production fitter residual and
Jacobian. They check an exactly recoverable degree-5 mode, all degree-5--6
Jacobian columns, several noisy modes, rank loss for a deliberately unseen
mode, and weak/strong regularization limits. A separate round trip sends all
56 degree-5--8 coefficients through the production Python writer and strict
C++ reader. Existing tests cover Python/C++ basis equality, profile counts and
reason codes, and the legacy workflows. The integration replay requires the
maximum relative NG and RG pressure error to be no larger than `1e-6` and
records the achieved values rather than merely assuming agreement.

Passing A6 demonstrates consistency among these implementations and fixtures.
It is not an independent validation of coma physics, COPS calibration, a
full-descent optimizer result, or the uniqueness of a localized source. The
complete requirements, numerical thresholds, current measured diagnostics,
and failure interpretation are in [`A6_VALIDATION.md`](A6_VALIDATION.md).

#### Performance characterization (Stage A7)

Stage A7 measures the production fitter rather than extrapolating from the
isolated harmonic-basis benchmark.  `tools/benchmark_fit_performance.py`
loads a real manifested base in a fresh process per degree, constructs the
same problem used by `fit_cops.py`, times repeated residual and analytic-
Jacobian calls, records peak RSS, measures equal-work scaling across
independent processes, and makes then deletes one-evaluation Stage-A5 output
packages to measure their size.  Those output probes are not fitted results.

On the measured 100,000-facet, 1,150-point base, degrees 5--8 required about
89--90 s to reach a prepared problem, about 1.90--1.93 s per residual, and
5.87--8.19 s per Jacobian.  Prepared-process peak RSS was about 11.0 GiB and
each complete output probe was 32.4 MiB.  Independent-job efficiencies for
1, 2, 4, 8, and 14 jobs were 1.000, 0.931, 0.774, 0.667, and 0.453.  The
declared 60% efficiency and conservative 70%-of-available-memory gates
therefore recommend `--jobs 8 --threads-per-process 1` on that host when at
least eight lambdas are available.  This measured recommendation is supplied
explicitly; it does not silently change the profile runner's portable default.

For the declared planning case of three profiles, degrees 5--8, five lambdas,
and 1,000 optimizer evaluations, the measured-rate estimate is 147.53 serial
CPU-hours, 27.90 hours at the recommended throughput, and 1.9 GiB of output.
It is an upper capacity scenario, not a convergence forecast: actual SciPy
call counts, retries, C++ verification, and I/O contention are not predicted.

The exact command, tables, metric definitions, safe output-name convention,
and interpretation limits are in
[`A7_PERFORMANCE.md`](A7_PERFORMANCE.md).  Repeat it on the target host before
a campaign or after changing the shape, selection, kernels, libraries, or
fitter implementation.

The same document inventories higher-resolution 67P models. The current
external SHAP5 BDF has 100,000 faces. The downloaded NASA PDS SPC SHAP5
products at nominal 199k, 384k, 788k, 1M, and 3M resolution have been
converted to metre-unit `GRID*`/`CTRIA3` BDFs. The raw lower four DSKs contain
small numbers of reversed exact duplicate triangle pairs that violate the
strict two-manifold check; the converter removes only those zero-thickness
fin pairs, audits every remaining edge, and records exact removals and hashes
in a manifest. Start a controlled geometry-resolution test with the converted
199k BDF and one fitter job. Larger SHAP7/SHAP8
PLY products are not directly supported, and no 1M-or-larger fitting-base
memory/runtime claim has been measured.  Increasing facet count refines
geometry, shadowing, and visibility; it does not increase the degree of the
separate harmonic activity parameterization.

#### Mesh-bound late-visibility and effective-slope hypotheses

`tools/build_physical_activity_templates.py` adds two deliberately narrow,
auditable hypotheses to the controlled residual fit. It does not infer a
region from the pressure residual and does not fit anything itself.

The `late_ng_not_rg` field is binary. For every cached observation in the
explicit `--late-start`/`--late-stop` interval, a facet enters the union when
its NG kernel has `Kn > 0` while that same observation has exactly
`KPhi = KV = 0` for RG. Thus “seen” means nonzero response in the implemented
instrument/transport kernel, not merely geometric line of sight. The union is
evaluated observation by observation because rotation can move one facet into
the RG aperture later. Integrated NG/RG sensitivities are saved as diagnostics
but never used to tune or threshold this region. The interval must be declared
from instrument evidence; the tool will not discover a favorable interval by
looking at the fit residuals.

The cliff fields are binary effective-slope thresholds. With outward facet
normal `n`, body-frame position `r`, and SPICE angular velocity `omega`, the
implemented acceleration and slope are

```text
g_eff(r) = g_polyhedron(r; rho) - omega x (omega x r)
slope(r) = acos[(-n . g_eff) / |g_eff|].
```

`g_polyhedron` is the exact homogeneous-polyhedron field evaluated by the
pinned ESA/TU Munich `polyhedral-gravity` implementation. Density is an
explicit fixed input, not fitted. The angular velocity is evaluated at the
start, midpoint, and stop with `spice_angular_velocity`, which compiles the
same `Geometry.cpp` and reads the same SPICE input as LiouvilleCOPS. The
reported median vector is used only after the tool records the maximum sampled
variation. Evaluation is offset a declared distance along the outward normal;
a deterministic surface-wide sample must pass an offset-halving slope test.
An optional second gravity shape must pass a sampled resolution check. Winding
audit, package version/backend, shape hashes, density, rotation vector,
offset checks, slope range, and mask areas are all stored in the companion
JSON.

The production 788k experiment uses the same closed
`cg_spc_shap5_788k_cart.bdf` for the transport surface and homogeneous gravity
body. The older 100k `SHAP5_Volume.bdf` is not treated as an interchangeable
reduction: its minimum-Z extent differs by about 70 m, its signed-volume
magnitude differs by about 1.8%. Its earlier 199k gravity comparison used
the defective package HEAL operation and is inconclusive. The comparison
must be repeated with audited input winding and HEAL disabled. That old
cache and the original 788k gravity/template products are quarantined with
`_heal_defect` suffixes and must not supply cliff masks. A
cross-shape comparison at a 0.1 m offset is meaningful only when both surfaces
place the evaluation points on the same side of their boundaries. The
corrected campaign uses the active 788k surface as the exact gravity boundary.
The offset-halving audit alone cannot detect a corrupted polyhedron: both
evaluations can use the same wrong face orientations. Independent topology,
orientation and reference-gravity checks are required. See
[`HEAL-corruption.md`](HEAL-corruption.md) and
[`GRAVITY_REPAIR.md`](GRAVITY_REPAIR.md).

The expensive per-facet calculation is also persisted in
`<bundle-stem>.gravity.npz`. It contains homogeneous self gravity and the
reference SPICE-dependent centrifugal/effective-gravity and slope arrays.
`tools/gravity_cache.py` hashes the exact centroid, normal, and area sequence;
the gravity-shape bytes, density, offset, and algorithm are part of the cache
definition. A subsequent builder invocation can supply this file through
`--gravity-cache`. Exact self gravity is reused, while centrifugal/effective
gravity and slope are recomputed inexpensively if the requested SPICE angular
velocity differs. A mesh, ordering, shape-file, density, offset, checksum, or
algorithm mismatch is fatal. This cache is diagnostic/derived data and is not
an additional fitted field.

The default `cliff_ge_60deg` threshold follows the published statistical
cliff definition. A `cliff_ge_45deg` field can be generated as a sensitivity
case because 45--90 degree terrain has also been classified morphologically
as cliffs. These are categorical predictors, not a claimed continuous
production-versus-slope law: the available cited observations motivate
testing enhanced activity at cliffs, but do not supply a calibrated gas-source
function for this end-of-mission state.

The fitter selects a bundle column with repeatable
`--activity-template BUNDLE NAME`. Each selected field receives one
nonnegative coefficient `c_name`. In the controlled mode,

```text
DeltaA(e) = A_scale [sum b_lm Y_lm(e) + sum c_k T_k(e)],  c_k >= 0
A_H2O(e) = max(0, A0_linear(e) + DeltaA(e)).
```

Legacy activity, the three analytical-temperature parameters, and CO2
`kappa` must be fixed. Harmonics and templates share one area-weighted
regularizer, including cross terms, so an overlapping template is not granted
an independent penalty. Exact linear dependence is rejected. The `.npz`
bundle contains a hash of the complete ordered centroid/area table and cannot
be loaded with another mesh. A fitted mixed field is exported as the exact
per-facet format documented in `ACTIVITY_ENHANCEMENT_FILE.md`; C++ verifies
every facet centroid and area before replaying it.

These additions improve hypothesis testing, not spatial uniqueness. NG-only
and cliff fields can overlap each other and the degree-5--8 harmonics, and a
nonnegative best-fit amplitude can still be weakly constrained or pinned at
zero. A late-NG amplitude that disappears when the instrument-transition data
are excluded is selection sensitivity, not evidence for a surface source.
Neither template adds thermal inertia, volatile depletion, cliff age, dust
mantling, a sublimation energy balance, gas collisions, or a time-varying
activity map. Those missing mechanisms must not be attributed to the fitted
coefficient.

Scientific inputs for this implementation are the globally homogeneous
gravity interpretation of Pätzold et al.
([Nature 530, 63--65](https://doi.org/10.1038/nature16535)); the rotating,
uniform-density effective-slope method and 45--90 degree terrain class of
Groussin et al. ([A&A 583, A32](https://doi.org/10.1051/0004-6361/201525660));
the explicit slope-above-60-degree statistical cliff mask of Vincent et al.
([MNRAS 469, S329--S338](https://doi.org/10.1093/mnras/stx1691)); and the
observational active-cliff hypothesis of Vincent et al.
([A&A 587, A14](https://doi.org/10.1051/0004-6361/201527159)). The numerical
gravity evaluator is the validated ESA/TU Munich implementation described by
Wittmann et al. ([JOSS 9, 6384](https://doi.org/10.21105/joss.06384)). These
references motivate the declared test; none supplies a calibrated mapping
from local slope to the H2O amplitude fitted here.

`ActivityModel = shap5` or `shap5_1` selects local embedded coefficients.
`file` reads three rows of 25 coefficients, for H2O, CO2, CO.
`ActivityFileHeliocentricModel` then controls the H2O distance law.
The coefficients carry the source amplitude in this convention; the
activity map is not constrained to [0,1] and is not an active-area probability.

With heliocentric distance r_h in AU:

```text
H2O: q_max = 600 r_h^-p,  q_min = 0.2 q_max
     p = 4.2143229 for shap5, 4.8 for shap5_1
CO2: q_max = 600 r_h^-2,  q_min = 0.1 q_max

sunlit, unshadowed: q = q_min + (q_max-q_min) cos(SZA)
night or shadow:   q = q_min
```

These q factors must be used together with the amplitude convention of the
coefficient table. The total production at a time is
`Q_s(t) = sum_e F_s(e,t) area_e`, in molecules/s.

`SZAMode = srccg_table` replaces the production-law SZA by its truncated
integer degree (limited to 89 on the illuminated branch);
`continuous` uses the actual cosine. This option does not quantize the
analytical temperature law.

Cast shadows are tested from face centroids toward the Sun. Shadowed faces
use the night-side production and temperature. Subdivision for observer
integration does not refine the activity map or the static face shadow state.

For `TemperatureModel = analytical`:

```text
a_T = (Tmax-T75)/(1-1/cos(75 deg))
b_T = Tmax-a_T

T_raw = max(Tmin, a_T/cos(SZA)+b_T)   if SZA <= 75 deg and unshadowed
T_raw = Tmin                         otherwise
T_surface = SurfaceTemperatureCorrectionFactor * T_raw
```

The supplied analytical parameters are 182.09922, 136.11149, and 133 K.
The C++ solver also supports `constant` and `bjorn` temperatures.
The latter uses the local temperature table with a configured distance
index 0–5 (default 4), not a continuous lookup at the current heliocentric
distance. The fitter supports only the analytical law.

The analytical formula can jump at 75 degrees when T75 exceeds Tmin.
Without ordered fitting, parameter names do not guarantee
`Tmax >= T75 >= Tmin`, and clipping can make nominal temperature
parameters ineffective.

### 4.4 Static collisionless transport

For surface temperature T and molecular mass m:

```text
beta = sqrt(m/(2 k T))
h(v) = 2 beta^4 v^3 exp(-beta^2 v^2)
integral_0^infinity h(v) dv = 1
integral_0^infinity h(v)/v dv = sqrt(pi) beta/2
```

The emitted angular flux law is:

| `SurfaceBC` | Emitted flux per surface area and solid angle |
|---|---|
| `flat` | `F cos(theta)/pi`, outward hemisphere. |
| `rough` | `F/(2*pi)`, outward hemisphere; an assumed angular law, not explicitly resolved surface roughness. |

For each triangle, the static solver subdivides into 4^level congruent
triangles and samples their centroids. It increases level until the maximum
edge is no greater than `SubdivisionAccuracy * distance_to_triangle`,
or reaches `MaxSubdivisionLevel`. The supplied values are 0.05 and 8.
The maximum level can stop refinement before the requested geometric
criterion is met; these controls are not a certified pressure-error bound.

Let l point from surface to gauge, b be the RG boresight, W_NG the NG
response, and W_RG the RG acceptance. Summing over visible sub-elements
with angular/geometric weight w gives face kernels:

```text
Kn   = sum w W_NG
KPhi = sum w W_RG (-l.b)
KV   = sum w W_RG
```

The predictions for each species are

```text
n_NG,s = sum_e F_s,e sqrt(pi) beta_s,e/2 Kn_e
Phi_RG,s = sum_e F_s,e [KPhi_e + (V.b) sqrt(pi) beta_s,e/2 KV_e]
```

`n_NG,s` is acceptance-weighted density. It equals an unweighted local
density only with the appropriate full-acceptance/no-head-on settings.
`SpacecraftVelocity = off` sets the V.b term to zero.

In the static solver, acceptances are evaluated using l; spacecraft motion
enters through the additive RG flux term. Thus this is not a full
velocity-dependent integration of spacecraft-frame acceptance. That
distinction matters if velocity or narrow angular cutoffs become important.

### 4.5 Gauge response and pressure conversion

NG acceptance options:

| Option | Implemented behavior |
|---|---|
| `full` | All incident directions accepted. |
| `hemisphere` | Front hemisphere only. |
| `ik_fov` | Directions inside the instrument-kernel FOV. |
| `srccg_shield` | Full front acceptance; partial rear transmission through a modeled tube, and only for directions arriving from the RG front. |

When `NGHeadOnWeight = on`, front-side NG acceptance is additionally
multiplied by `1-l.b_NG`, even when another NG acceptance option is used.
The shielding table is generated with a fixed random seed, by default
2,001 angular bins and 20,000 samples/bin. It is an instrument-response
approximation, not a spacecraft CAD occlusion model.

RG supports `hemisphere` or `ik_fov`. FOV clipping currently supports
RECTANGLE/POLYGON shapes only.

The pressures are

```text
p_NG = sum_s n_NG,s k T_NG / S_s
p_RG = sum_s Phi_RG,s C_s(T_RG) / S_s
```

where S_H2O = 0.893 and S_CO2 = 0.704 in the supplied configuration;
`SpeciesSensitivity = off` uses S_s = 1.

| `RGPressureModel` | C_s(T_RG) |
|---|---|
| `equilibrium_sphere` | `sqrt(2*pi*m_s*k*T_RG)` |
| `calibration_srccg` | `sqrt(pi*m_s*k*T_RG/2)` |
| `dsmc_srccg` | `(pi/2)*k*T_RG*sqrt(m_s/(2*k*T_RG))` |

The supplied input files explicitly choose `equilibrium_sphere`, while
the C++ fallback when the key is omitted is `calibration_srccg`.
The former factor is exactly twice the latter. Changing this convention
can strongly change inferred source rates and mixture.

Comments associate the equilibrium-sphere relation with the ROSINA
instrument documentation. The formula is confirmed in the code; the cited
instrument manuals themselves are not bundled in the inspected directory,
so their applicability and precise calibration interpretation were not
independently verified here.

`T_NG` defaults to 293 K. `T_RG` is constant or linearly interpolated
between HK midpoint samples. Within the outermost measurement intervals
the nearest endpoint temperature is used; a time outside the overall
interval coverage is an error. Gas/surface temperature and gauge temperature
are distinct quantities.

### 4.6 Flight-time rotation

Nucleus orientation at each observation time is included in **both**
`FlightTimeRotation = none` and `exact`. The switch concerns rotation
while an individual molecule travels.

The `exact` branch integrates speed numerically, rotates the emitting
point/normal back to an estimated emission time, includes emitting-surface
motion and a retarded-time Jacobian, and checks the corresponding curved
path in body coordinates. Source illumination is rotated back in time.
Gauge acceptance uses the relative velocity computed in that branch.

The name `exact` does not remove numerical/model approximations:

- Angular velocity is held constant during each flight.
- Flight time uses six fixed iterations, without a convergence test.
- Gauss–Legendre quadrature is truncated at `beta(T_reception)*v = 5.5`
  with 16 nodes by default.
- The path approximation has a maximum of 256 segments.
- Cast-shadow flags remain those calculated at reception, rather than
  being ray-traced afresh at emission.
- A heuristic prefilter discards faces sufficiently far behind their planes.
- Heliocentric distance is held at reception.

The use of body-frame spacecraft velocity in this branch, along with
surface rotation and positional offsets, deserves an independent frame
consistency test before extending it to substantially different conditions.
The existing self-tests use zero spacecraft velocity for the rotation tests.

Version-2 exact-rotation kernels depend on temperature/species. The present
fitter rejects them and uses static version-1 kernels. An exact-rotation
rerun is a subsequent check, not a temperature-consistent exact-rotation fit.

### 4.7 Rotation-driven source variability and controls in the fit

There are two distinct rotation effects. The nucleus changes orientation
over the observation interval, and it can also rotate while an individual
molecule travels to the spacecraft. **The current fit includes the first
effect, but not the second.** `FlightTimeRotation = none` does not freeze the
nucleus for the approximately ten-hour dataset: it freezes transport geometry
separately at each observation time.

#### Observation-time geometry and illumination

For every requested observation timestamp, `main.cpp` calls
`cSpiceGeometry::Evaluate()` and then `cSourceModel::Evaluate()`. SPICE supplies
spacecraft position/velocity, gauge pointing, and the Sun direction in the
configured nucleus-fixed `BodyFrame` (supplied: `67P/C-G_CK`). The frame and
loaded orientation kernels determine the rotation; the fitter does not infer
the spin period, pole, rotational phase, or an illumination time offset.

For facet e, the local normal stays fixed in body coordinates while the
body-frame Sun direction changes with time:

```text
mu_e(t) = cos(SZA_e(t)) = normal_e . SunDirection_body(t)
```

When `CastShadows = on`, the forward solver ray-traces from each illuminated
facet's centroid toward the Sun at that observation time. A facet is treated
as night/shadow when `mu <= 0` or the Sun ray is blocked by the nucleus.
Visibility from the spacecraft and gauge acceptance are separate geometric
tests, also recomputed for each observation. Thus pressure variability
contains both changing illumination/source strength and changing observer
geometry; it is not a direct measurement of the global production rate.

#### Temporal source rate: fixed activity map, changing illumination

The spatial activity factor is attached to the nucleus:

```text
A_s,e = max(0, sum_j a_s,j Y_j(e))

F_s,e(t) = A_s,e [q_min,s + (q_max,s-q_min,s) mu_e(t)]  sunlit, unshadowed
F_s,e(t) = A_s,e q_min,s                               night or shadow
Q_s(t)   = sum_e F_s,e(t) area_e
```

These expressions use the supplied `SZAMode = continuous`. The q factors
also vary with heliocentric distance through the prescribed laws in section
4.3. Night/shadow production factors are 0.2 of q_max for H2O and 0.1 for
CO2. `SZAMode = srccg_table` instead truncates the illuminated production-law
angle to integer degrees, limited to 89 degrees; it does not quantize the
analytical temperature law.

One set of H2O activity coefficients and one CO2 multiplier are fitted across
the selected interval. They are not independent functions of time. The
facet fluxes and total production can nevertheless vary as active areas
rotate into/out of sunlight and shadows. The current CLI has no free
per-observation source amplitude, time-dependent activity map, or fitted
night/day production ratio.

#### Temporal surface temperature: an instantaneous illumination law

For the analytical law used by the fitter:

```text
a = (Tmax-T75) / (1-1/cos(75 degrees))
b = Tmax-a

T_raw,e(t) = max(Tmin, a/mu_e(t)+b)  if unshadowed and SZA <= 75 degrees
T_raw,e(t) = Tmin                    otherwise
T_surface,e(t) = SurfaceTemperatureCorrectionFactor * T_raw,e(t)
```

`Tmax`, `T75`, and `Tmin` are shared parameters across the fitted interval,
not a separate temperature for each observation. Each facet's temperature
changes as its illumination changes. Beyond 75 degrees or in shadow the
analytical model immediately uses Tmin, so it can jump at the threshold
when T75 exceeds Tmin. The temperature correction factor must be 1 for the
current Python/C++ fitting agreement.

Production and temperature are separate prescribed functions of illumination:
the emitted flux is not derived from the temperature through a sublimation
law. There is no solved surface energy balance, thermal inertia, subsurface
heat conduction, delayed heating/cooling, or memory of previous illumination.
Morning and evening facets with the same activity and instantaneous
illumination/shadow conditions follow the same prescriptions. The fit cannot
estimate an afternoon thermal lag with its current parameters.

Surface temperature changes molecular speeds, hence the NG density and the
spacecraft-velocity contribution to RG flux (sections 4.4 and 6.6). The NG/RG
instrument temperatures used for pressure conversion are separate quantities;
the fitted surface-temperature law does not replace the measured HK sensor
temperature.

#### How the optimizer retains rotation-driven variability

The static base run stores a separate point-kernel file for each observation.
Each contributing face record includes `Kn`, `KPhi`, `KV`, `cos(SZA)`, and
the shadow flag; point headers include the Sun direction, heliocentric
distance, and spacecraft velocity projection. These are the time-dependent
geometric/illumination states of the SPICE base run, not a single averaged
orientation or SZA for the interval.

At each parameter trial the fitter uses these cached states to recompute
activity-scaled source fluxes, analytical temperatures, and thermal speed
factors, then sums the contributions for each observation. It does not call
SPICE or retrace shadows inside the optimizer. Changing the activity or
analytical temperature parameters does not require new geometric kernels;
changing the orientation geometry, observation times, shadow treatment, or
other base assumptions requires a new consistent forward base. Editing the
old input/header alone does not update the cached illumination.

#### Controls and their actual scope

| Control | What it does / does not do |
|---|---|
| `BodyFrame`, `KernelDirectory`, repeated `Kernel` entries | Choose the SPICE orientation/geometry used at each timestamp. There is no independent rotation-period or phase-fitting switch. Changes require a new base; preserve exact kernel identities for comparisons. |
| `FlightTimeRotation = none` | Required for the current fitter. Includes observation-to-observation rotation, but neglects rotation during molecular flight. |
| `FlightTimeRotation = exact` | C++ forward-only flight-time correction; not supported by the current fitter or boundary-rebuild driver. See section 4.6 for approximations. |
| `CastShadows = on/off` | Enable/disable topographic sunlight obstruction. `off` does not disable the day/night test or nucleus rotation. Rebuild before fitting a changed setting. |
| `SZAMode = continuous/srccg_table` | Select the source-rate illumination prescription described above, not a rotation model. Rebuild to keep the base table and fit model consistent. |
| `TemperatureModel = analytical` | Required by the fitter. `constant`/`bjorn` are C++ forward alternatives, not additional fitter modes. |
| `AnalyticalTmax`, `AnalyticalT75`, `AnalyticalTmin` | Initial surface-temperature-law parameters from the base input; fitted unless fixed. |
| `--fix-temperature` | Freeze those three parameters at base values, not the time-dependent facet temperatures. Illumination-driven temperature variability remains. |
| `--ordered-temperature`, `--tbounds LO HI` | Control parameter ordering/bounds, not thermal lag. In ordered mode only Tmin receives the stated upper bound; see section 5.1. |
| `SpacecraftVelocity = on/off` | Controls a velocity contribution to the gauge calculation; it is not an on/off switch for nucleus illumination changes. |
| `FirstPoint`, `LastPoint`, `PointStep` | Choose which observation timestamps receive forward evaluations/kernels. They do not specify the rotation period or a physical time-integration step. |
| `FitInterval`, exclusion intervals, fitter selection switches | Select observations for the misfit; do not change their SPICE orientation or physical source prescription. |

For a deliberately constant-surface-temperature experiment compatible with
the fitter, use an analytical base with all three temperature parameters
equal, and fit with `--fix-temperature`. This keeps T constant but still
allows rotation-driven source-rate and visibility changes. Setting only
`--fix-temperature` on the ordinary unequal-parameter base does not do that.
There is no standard input switch to hold the entire nucleus orientation
fixed across the observation window.

#### Assess the neglected rotation during molecular flight

After fitting with the supported static base, make a separate copy of its
generated `cops_fit.in`. Keep the fitted activity and temperature parameters,
set `FlightTimeRotation = exact`, and choose a fresh `OutputDirectory`.
Run the C++ executable directly and compare it with the static C++ verification
run (section 6.4). `SpeedNodes`, `SpeedRangeMax`, and `PathTolerance` control
the numerical flight-time calculation, not the long-duration spin evolution.

The exact branch rotates source illumination back to the estimated emission
time and includes surface motion, but holds angular velocity constant during
each flight and retains reception-time cast-shadow flags. It adds no thermal
memory. Its species/temperature-dependent version-2 kernels cannot be passed
to the present optimizer. This comparison measures the forward correction
at fixed fitted parameters; it is not a refit with that correction included.
Do not assume it is negligible or important solely because the observations
span approximately ten hours: the flight-time effect concerns the much
shorter travel time of each contributing molecule and the evolving geometry.

## 5. What is fitted and how

Implementation: current `tools/fit_cops.py` (whose numerical implementation
originated in former RCS revision 1.9), especially
`Model`, `Problem`, `fit_one()`, and `main()`.

### 5.1 Free and fixed quantities

| Quantity | Treatment |
|---|---|
| H2O activity | 25 coefficients by default, in the local basis through nominal degree 4. `--degree 0..4` frees only the corresponding lower-degree groups. Higher coefficients retain prior values; the total map is not truncated. |
| Surface temperature | Tmax, T75, Tmin fitted unless `--fix-temperature`. |
| CO2 | One nonnegative multiplier kappa on the base CO2 activity map. `--fix-kappa` keeps it at 1. |
| Shape, trajectory, pointing, gauge responses, cast shadows, angular emission law | Fixed by the base-run geometry/kernels and header. |
| T_NG, T_RG, sensitivities, distance law, night/day ratios | Fixed, not fitted. |
| Observational noise or systematic calibration | No free noise, timing-offset, or calibration-nuisance parameters. |
| CO | Not simulated or fitted. |

The usual parameter count is 25 + 3 + 1 = 29.
H2O coefficients are unbounded; nonnegative emission is enforced by
clipping the resulting map at zero. kappa is bounded below by zero.
Coefficients are internally divided by the maximum absolute prior H2O
coefficient for scaling.

Ordinary temperature fitting bounds each of Tmax, T75, Tmin independently
to `--tbounds`, default 30–400 K. It does not impose their ordering.

`--ordered-temperature` instead fits
`Tmin, dT75=T75-Tmin, dTmax=Tmax-T75`, with nonnegative increments.
**In this mode the upper bound is applied only to Tmin; the increments
have no upper bound, so Tmax and T75 can exceed 400 K.**
The base temperature parameters must already be ordered.

`SurfaceTemperatureCorrectionFactor` is implemented in C++ but is not
read/applied by the current Python model. For the intended fitting
workflow keep it at 1. A different effective temperature law requires a
matching Python implementation, not a relaxed consistency tolerance.

#### Explicit parameter vector and parameter counts

For the default degree-4 fit without freezing options, the physical parameter
vector is:

```text
x_physical = [a_H2O,0 ... a_H2O,24, Tmax, T75, Tmin, kappa_CO2]
```

The optimizer uses scaled H2O coefficients internally. With ordered
temperature enabled it replaces the last temperature triple by
`[Tmin, T75-Tmin, Tmax-T75]`; this still contributes three free parameters.
The same coefficients, temperature-law parameters, and CO2 multiplier apply
to the entire selected observation interval. There are no independent
per-observation source amplitudes or temperatures.

For the supported nominal degree d from 0 to 4:

```text
N_parameters = (d+1)^2
               + 3 if temperature parameters are fitted
               + 1 if the CO2 multiplier is fitted
```

| `--degree` | Free H2O coefficients | Total with temperatures and CO2 scale free | Total with both `--fix-temperature --fix-kappa` |
|---|---:|---:|---:|
| 0 | 1 | 5 | 1 |
| 1 | 4 | 8 | 4 |
| 2 | 9 | 13 | 9 |
| 3 | 16 | 20 | 16 |
| 4 | 25 | 29 | 25 |

At degree 4, `--fix-temperature` alone leaves 26 parameters, while
`--fix-kappa` alone leaves 28. There is no dedicated CLI switch to freeze
all H2O coefficients: even `--degree 0` fits the degree-zero coefficient.
Reducing degree leaves all higher-degree coefficients at their base values;
it does not zero them or turn the full map into a degree-d truncation.

The H2O coefficients may be signed and have no finite box bounds. Only the
resulting activity `max(0, Y a)` is constrained to nonnegative emission.
There is no independent fitted total-production parameter or upper bound
on total production. CO2 retains the base spatial pattern and only its
nonnegative scalar multiplier changes. Temperature parameters describe one
common surface-temperature law, not separate H2O/CO2 temperature fields.

These are nominal optimizer dimensions, not numbers of independently
measured physical quantities. Activity clipping, the temperature floor,
parameter degeneracy, and regularization can make some directions weakly
constrained or locally ineffective. For example, a nominal temperature
parameter can have no effect where the floor always selects Tmin. Inspect
the fit's Jacobian/conditioning and bound diagnostics rather than equating
29 free parameters with 29 independently determined values.

#### Physics prescribed rather than optimized

Here "fixed" means the optimizer cannot change the setting or its prescribed
time history. It does not necessarily mean constant with time.

| Quantity / assumption | What is fixed during the fit |
|---|---|
| Shape and spatial basis | Mesh, facet geometry/areas/order, and the local basis definition. H2O coefficients change, but the basis and shape do not. |
| Rotation and observation geometry | SPICE orientation, spacecraft position/velocity, observation timestamps, gauge position offset and pointing. No spin-period, pole, phase, trajectory, or timing-offset fit. |
| Angular surface boundary condition | `flat` or `rough` selected when the base was generated; not an optimized discrete choice or a fitted mixture. |
| Illumination and shadows | Each observation's cached solar angle/shadow state and the chosen `CastShadows`/`SZAMode` prescription. |
| Source illumination/distance law | H2O night/day ratio 0.2, CO2 ratio 0.1, H2O distance exponent 4.2143229 (`shap5`) or 4.8 (`shap5_1`), CO2 exponent 2, and the implemented normalization convention. These are not fitted separately from activity. |
| CO2 spatial shape | Base map; only kappa may vary. `--fix-kappa` holds kappa at 1 relative to that base map. |
| Surface-temperature functional form | Analytical solar-angle law, floor and 75-degree transition. Only its three temperature parameters can vary; the transition angle is not fitted. |
| Temperature correction | `SurfaceTemperatureCorrectionFactor = 1` is required for this workflow, not a fitted extra scaling parameter. |
| Gauge angular response | NG shielding/head-on weight, RG acceptance, instrument geometry, and other acceptance parameters encoded in the base. |
| Pressure calibration | RG pressure-conversion choice, species sensitivities and NG/RG sensor temperatures. The supplied base uses sensitivity factors 0.893 (H2O), 0.704 (CO2), `TemperatureNG = 293 K`, and an HK-derived RG temperature. |
| Transport treatment | Collisionless model, the selected spacecraft-velocity contribution, static flight geometry, and numerical integration settings of the base. |
| Species set | H2O and CO2 only; CO coefficients are retained for file-format compatibility but CO contributes no modeled pressure. |

In particular, rotation, solar angle, shadowing, visibility, heliocentric
distance, source flux, surface temperature, and measured RG temperature can
still vary across the dataset. Their geometry or governing prescriptions
are fixed; the optimizer changes the allowed source/temperature parameters
against those time-dependent states. `--fix-temperature` fixes the three
analytical-law parameters, not the temperature of every facet at every time.
See section 4.7 for the temporal treatment.

Changing cached physical assumptions requires a new consistent forward base,
then a fresh fit. The boundary-change driver in section 6.6 handles this for
`SurfaceBC`. Analytical-temperature and activity changes made as optimizer
trials already use the static cache correctly; do not edit an old base's
headers/input merely to make incompatible cached results appear current.

#### User choices, derived outputs, and missing physics

`--lambda`, degree, parameter-freezing switches, bounds, data-selection
options, and convergence tolerances are user-controlled problem/algorithm
settings, not physical quantities inferred by the optimizer. A lambda scan
runs independent fits at prescribed weights; its curvature diagnostic is
post-processing, not a fitted physical parameter. The H2O activity prior is
fixed, and its regularization penalty does not directly constrain temperature
or the CO2 multiplier.

NG and RG have equal aggregate misfit weight, with equal point weights within
each selected gauge dataset. The model does not fit observation uncertainties,
correlations, calibration offsets/gains, or extra noise parameters. Selection
flags, measured reference pressures, and missing-data handling are not free
optimization variables.

Total production rates, facet fluxes, predicted/reconstructed pressures, and
facet temperatures are derived from the fitted parameters. They are not
additional independent fitting parameters. Obtain fully recomputed production
and species diagnostics from the C++ verification run; the Python fitted
table retains some base-run fields (section 7.3).

Thermal inertia/lag, subsurface conduction, a coupled sublimation/energy-balance
solution, time-dependent activity coefficients, collisions, gravity, and
photochemical losses are not adjustable frozen parameters here: they are
absent from the implemented fitting model. Similarly, the implemented
exact-flight-time-rotation forward solver is not available as a fit option.
Adding these effects requires model/optimizer development, not enabling a
hidden free parameter.

### 5.2 Data selection, including defaults that change the base flags

The C++ run marks a gauge usable when its reference pressure is positive,
the gauge is outside the nucleus, the time is inside `FitInterval` if
supplied, and it is outside that gauge's exclusion intervals.
Interval boundaries are inclusive. This only flags observations; kernels
are still generated at the requested points.

The Python fitter then applies additional rules:

1. Load base pressures/flags, or apply `--data FILE`.
2. Read exclusion intervals from the base run's original input file.
3. By default, **without `--data`**, re-enable otherwise valid NG points
   in `ExcludeNG`, while honoring `FitInterval`.
4. With `--data`, normally use the table's flags as provided. The
   default does not automatically re-enable NG points removed in that table.
5. `--apply-exclude-ng` forces NG exclusion in the input's `ExcludeNG`
   intervals, including when `--data` is supplied.
6. `--exclude-ng-saturated` removes saturation-flagged NG points.
   Without it those points are not removed.
7. Build positive-pressure log residuals for each gauge's selected points.

For new profile comparisons, prefer `--profile NAME` or
`tools/run_fit_profile.py`. A named profile replaces steps 3--6 with one
versioned policy implemented by `tools/cops_selection.py`; it cannot be mixed
with `--apply-exclude-ng` or `--exclude-ng-saturated`. The checked-in profiles
are:

| Profile | NG/RG, 100k reference | NG/RG, 788k base | Policy on the audited 1,150 rows |
|---|---:|---:|---|
| `legacy_base_mask` | 839 / 1,146 | 839 / 1,150 | Positive outside-shape data, base FitInterval, and base `ExcludeNG`/`ExcludeRG` intervals. The final four NG points are already excluded by the run interval. |
| `all_available` | 930 / 1,146 | 934 / 1,150 | Every positive outside-shape datum in the base FitInterval, including the 13 saturated and post-transition NG points. |
| `transition_cut_100m` | 869 / 1,146 | 873 / 1,150 | Same candidates, but NG is excluded on the closed interval `10:24:42`--`10:37:10`; it resumes at `10:37:14` (99.48 m on the 100k reference, 113.75 m on the 788k shape). |

The profiles derive flags from pressure, shape status, FitInterval, and their
declared interval policy; they do not trust stale source `useNG/useRG` flags.
All interval endpoints are inclusive. Each checked-in profile stores both the
100k reference counts and shape-independent counts after all data/time/profile
rules but before the `inside` test. The fitter permits final counts to change
only by the number of those otherwise eligible points that the active mesh
classifies inside. A changed pressure, timestamp, FitInterval, exclusion, or
saturation decision remains fatal. `fit_data_used.txt` records a reason code
for both gauges at every point. `selection_profile.json` preserves the effective policy, its
source hash, reason counts, and reason hash; `fit_manifest.json` records the
same identities. When combined with `--data`, a named profile requires the
replacement table to contain every base-run point; it replaces pressures but
the profile, not the table's use flags, determines selection.

The HK table is read for saturation identification even when
`--exclude-ng-saturated` is absent. The current fitter therefore needs it
for a normal run and for `--data-only`.

For the removed, audited SHAP5 base, the selections were:

| Selection | NG | RG |
|---|---:|---:|
| C++ base flags; also Python `--apply-exclude-ng` without custom data | 839 | 1,146 |
| Current Python default without `--data` | 930 | 1,146 |
| Default plus `--exclude-ng-saturated` | 917 | 1,146 |
| Both exclusion switches | 826 | 1,146 |

These counts refer to the inspected data and input intervals. They are not
constants for other tables or shapes.

`--data` requires columns
`iPoint UTC pNG_Pa pRG_Pa useNG useRG T_RG_K`. The supplied collector
table includes additional columns for altitude, saturation, and reasons.
The fitter checks point membership, exact UTC-string equality, binary
use flags, positive selected pressures, outside-shape status, and agreement
of T_RG within relative 1e-6. Base-run points missing from the replacement
table are disabled. A different UTC spelling, even for the same physical
instant, fails the string check.

The replacement table can alter selection beyond the original base flags.
Apart from the explicit NG exclusion switch, the fitter does not generally
reapply all original FitInterval/ExcludeRG constraints to a custom table.
The descriptive `reason_*` and `ng_saturated` columns are not authoritative
fit controls: the fitter uses the flags and recomputes saturation from HK.
Duplicate replacement-table point IDs are rejected. Validate edited tables.

The collector's selection implementation is similar but not identical:
its default re-enables points classified as `ExcludeNG`. Because exclusion
reasons are assigned before testing FitInterval, overlapping intervals can
cause it to re-enable a point outside FitInterval. Inspect the actual fit's
`fit_data_used.txt` after any selection change.

### 5.3 Cached model evaluation

The fitter reads `elements.bin` for areas, basis values, and prior species
activity, and reads each version-1 point file for face indices, Kn/KPhi/KV,
SZA, shadow status, heliocentric distance, and V.b.

It recomputes H2O activity, per-contribution temperature, and species fluxes
for each parameter trial, then sums with NumPy and sparse matrices. No
SPICE queries or new ray tracing occur in the optimizer.

The static NG prediction scales locally as `F/sqrt(T)`. The main static
RG flux term scales with F and does not depend on gas temperature;
gas temperature enters its spacecraft-velocity term. The separately fixed
sensor temperature enters the pressure conversion. This explains why combining NG and RG can help separate
source amplitude and temperature, while composition and angular-response
assumptions still introduce degeneracy.

### 5.4 Objective and regularization

For each gauge g, selected point i, and N_g selected observations:

```text
r_g,i = [ln(p_model,g,i) - ln(p_observed,g,i)] / sqrt(N_g)

cost = 1/2 [mean_NG[(ln pressure ratio)^2]
            + mean_RG[(ln pressure ratio)^2]
            + lambda D^2]
```

Both gauges have equal aggregate weight. Individual observations within
one gauge have equal weight regardless of cadence, interval duration,
quality flags, or uncertainty. Consequently the densely sampled final
descent has more influence than an equal-duration sparsely sampled interval.
There is no robust loss, censored-data likelihood, or temporal covariance
model. Both gauges need nonempty selected data; the code does not provide
a guarded single-gauge fitting mode.

For the **unclipped** linear H2O map A_lin = Y a:

```text
D^2 = sum_e area_e [Y_e(a-a0)]^2
      / sum_e area_e [Y_e a0]^2
```

The regularization is a mesh-area-weighted distance from the prior map,
not a penalty on total production, temperature, CO2, or the clipped physical
activity. A Gram matrix `Y_fit^T diag(area) Y_fit` and Cholesky factor
construct the residual operator. Only varied coefficients enter its update.

A positive lambda penalizes changes over the whole surface, including
unobserved regions. It **does not guarantee that unseen regions remain
exactly at the prior**: global basis coefficients couple seen and unseen
areas. Neither the prior nor the data uniquely establish the global source
outside observed regions. The CLI default is **lambda = 0**, i.e. no
regularization; specify lambda explicitly for a regularized fit.

### 5.5 Optimizer and derivative checks

The code calls SciPy `least_squares` with trust-region reflective method,
analytic Jacobian, bound constraints, and `x_scale="jac"`.
Defaults are `max_nfev=1000` and ftol/xtol/gtol = 1e-10.

Derivatives include H2O map coefficients, all three temperature parameters,
CO2 scale, and the regularization residuals. Activity clipping and the
temperature floor create nondifferentiable boundaries; the code chooses
branches using strict positivity/temperature comparisons.
Model pressures are floored at 1e-300 before logarithms. Nonpositive
predictions are therefore not a physically satisfactory fit merely because
the log calculation remains finite.

Before optimization the fitter:

1. Reconstructs prior H2O activity and checks agreement with `elements.bin`
   (relative scale tolerance 1e-10).
2. Evaluates the Python model at the base parameters and compares it with
   saved C++ pressures (default `--check-tolerance 1e-6`).
3. Applies selection and writes the actual data table/plots.
4. Optionally stops with `--data-only`, or compares the Jacobian with
   central finite differences using `--check-jacobian`.

`--check-jacobian` prints errors and returns; it has no pass/fail error
threshold. It checks at the prior parameters, not all possible clipping
boundaries or the final solution. `--data-only` still loads kernels and
performs the forward consistency checks; it is not a lightweight text export.

### 5.6 L-curve and uncertainties

Several `--lambda` values launch independent fits from the same prior.
The implementation uses `multiprocessing.get_context("fork")` and a
worker pool (`--jobs`, default min(number of lambdas, 12)). It is designed
around a platform with fork support. Kernel arrays and per-evaluation
temporaries can use substantially more RAM than the binary files on disk.

The L-curve table records rho = sqrt(msNG+msRG), D, temperatures, kappa,
activity diagnostics, evaluations, termination status, and optimality.
With at least four positive-lambda, positive-D results, a numerical signed
curvature in (ln rho, ln D), parameterized by ln lambda, selects an interior
point. The selector does **not** reject fits that reached their evaluation
limit. A marked corner is a diagnostic, not a validated optimum.

Let J_d be the weighted data Jacobian and R the regularization Jacobian
including sqrt(lambda). The reported local uncertainty calculation is:

```text
H = J_d^T J_d + R^T R
Hplus = pseudoinverse(H), cutoff rcond=1e-14
p_eff = trace(Hplus J_d^T J_d)
v = sum(data residuals^2) / max(N_data-p_eff, 1)
covariance approximation = v Hplus
sigma = sqrt(diagonal(covariance approximation))
```

This is a local, regularized Gauss–Newton approximation in the chosen
parameter coordinates. It is not a full posterior sample or an uncertainty
budget for calibration, shape, temperature telemetry, time correlation,
or source-law errors. Bounds and clipping complicate its interpretation.
A null direction discarded by the pseudoinverse can yield zero reported
sigma; that does **not** mean a physically exact parameter estimate.
With ordered temperatures the reported parameter sigmas concern Tmin and
increments; marginal sigmas for reconstructed Tmax/T75 are not separately
propagated in the printed parameter rows.

Activity diagnostics use the clipped map, area weighted. Their
`activity_ratio` is not a total production ratio: it omits illumination
and distance factors. “Seen” marks faces contributing to either response at
any selected point, not a quantitative resolution or identifiability test
for the gauge actually selected at every such point.

### 5.7 Numerical workflow and diagnostic interpretation

This section describes the current implementation, not the proposed changes
in [ROADMAP.md](ROADMAP.md). The numerical code derives from former RCS
revision 1.9; the current working fitter additionally has shared-output path handling
described above.

#### From the forward model to one optimization

1. **Generate the forward base.** C++ evaluates SPICE geometry and source
   illumination at each observation time. It integrates visible sub-triangle
   contributions using adaptive geometric subdivision and saves static
   `Kn/KPhi/KV` kernels, solar angles, shadows, and point geometry. The
   subdivision size criterion is not an a posteriori pressure-error estimate.

2. **Load and verify the cached problem.** Python loads the element/basis
   arrays and version-1 kernels. It verifies the H2O base activity against
   its coefficients and the initial pressure predictions against C++.
   Failure means the fitting base or assumptions are inconsistent; increasing
   `--check-tolerance` is not a substitute for finding the cause.

3. **Construct the data mask.** The fitter reads base or replacement
   observations, applies its NG exclusion/saturation rules, and exports
   `fit_data_used.txt`. The effective dataset is this exported selection,
   not necessarily the original C++ `useNG/useRG` flags.

4. **Construct the parameter vector.** Free H2O coefficients are scaled by
   the maximum absolute base H2O coefficient. Surface temperatures use
   either direct or ordered-increment coordinates; kappa is a separate
   nonnegative parameter. Higher-degree frozen coefficients remain in the
   evaluated map. See section 5.1 for exact dimensions and bounds.

5. **Evaluate a trial.** Form and clip the H2O activity map, scale CO2,
   recompute analytical temperatures and thermal speed factors at the cached
   illumination states, and sum pressures with array/sparse operations.
   No new SPICE query or visibility/shadow ray tracing occurs in a trial.

6. **Build residuals and derivatives.** Use the log-pressure residuals
   below, optionally append the activity-prior residuals, and evaluate the
   analytic Jacobian. Activity clipping and temperature floors make this
   a piecewise-smooth nonlinear problem; parameter names alone do not ensure
   all derivative directions are active or identifiable.

7. **Optimize.** Call SciPy `least_squares` with trust-region reflective
   method (`method="trf"`), bound constraints, analytic Jacobian, and
   `x_scale="jac"`. The latter adds Jacobian-based optimizer scaling on top
   of the explicit coefficient scaling. The code uses the ordinary quadratic
   least-squares loss; it does not specify a robust loss. Defaults are
   `max_nfev=1000`, `ftol=xtol=gtol=1e-10`.

8. **Report, then verify independently.** Compute local diagnostics and
   write the fit products even when termination was the evaluation limit.
   Inspect those diagnostics, then recompute the fitted physical model with
   the generated C++ input. Successful output writing alone is not an
   optimizer-convergence or scientific-validation flag.

The objective is, for each selected gauge g with N_g points:

```text
r_g,i = ln(p_model,g,i / p_observed,g,i) / sqrt(N_g)
cost = 0.5 * [sum_i r_NG,i^2 + sum_i r_RG,i^2 + lambda D^2]

D^2 = sum_e area_e [Y_e (a-a_prior)]^2
      / sum_e area_e [Y_e a_prior]^2
```

The regularization uses the linear, unclipped activity map and acts only on
H2O coefficient departures from the base prior. A mesh-area Gram matrix and
Cholesky factor implement that penalty. It is not a direct constraint on
total production, temperature, or kappa. A positive lambda does not freeze
unobserved regions at the prior because basis coefficients are global.

Both gauges receive equal aggregate data weight; within a gauge each selected
point receives equal weight. This is not a calibrated measurement-uncertainty
likelihood. Changing cadence therefore changes the relative influence of
different time intervals. Saturation, time correlation, reported L4 errors,
and general housekeeping-quality fields are not automatically incorporated
into this residual model. Both gauges must have selected data; empty-gauge
handling is not a supported single-gauge inference mode.

#### What the numerical diagnostics mean

| Diagnostic | Interpretation and caveat |
|---|---|
| RMS log residual, `exp(RMS)` | Describes multiplicative pressure agreement. Roughly 6–7% agreement is not a 6–7% uncertainty on source activity, temperature, or production. |
| `status`, termination message, `nfev` | Status 0 is an evaluation-limit stop, not convergence. Positive status indicates a stopping criterion was met, not uniqueness or physical validity. |
| `ftol`, `xtol`, `gtol` | Stopping tests for objective change, scaled step, and first-order optimality. Tight values do not compensate for model mismatch or poor identifiability. |
| Singular values / condition number | Sensitive to the parameter coordinates, scaling, and whether regularization rows are included. The current fitter reports the stacked data-plus-regularization Jacobian; historical reports may differ. |
| Bound hits / temperature-floor dominance | Can reveal that parameters are constrained by imposed limits or inactive branches instead of the observations. Inspect alongside residuals. |
| `p_eff` | A local regularized effective-dimension diagnostic, not a count of independently measured global physical quantities. |
| `sigma` | Local approximation described in section 5.6; not a complete calibration/model-error budget. Zero can reflect a discarded null direction. |
| L-curve corner | A post-processing curvature diagnostic of the sampled fits, not a fitted parameter or an automatic scientific choice of lambda. Current selection can choose a nonconverged fit. |

The status and solver-option meanings are documented by
[SciPy's least-squares reference](https://docs.scipy.org/doc/scipy-1.15.1/reference/generated/scipy.optimize.least_squares.html).
Record the actually installed SciPy version for each experiment rather than
assuming this linked documentation specifies the local runtime version.

The printed covariance approximation forms `H = J_d^T J_d + R^T R` and
uses its pseudoinverse with a numerical cutoff. This is local to the final
parameters and depends on scaling, rank decisions, regularization, and the
chosen residual normalization. It is not posterior sampling. Bounds,
clipping, temperature floors, and unmodeled correlations/systematics limit
the interpretation. In ordered mode, physical Tmax/T75 uncertainties require
propagation of the full covariance of Tmin and increments, including their
cross-covariances; the printed increment sigmas alone are insufficient.

#### Minimum checks for a new fitting experiment

- Use a complete base with unchanged paired inputs/table/elements/kernels;
  save the command, source revision, data selection, and dependency versions.
- Inspect a data-only export before launching an expensive fit. Keep the
  zero/missing, saturation, range-change, and inside-nucleus policies explicit.
- Check initial Python/C++ consistency and the analytic Jacobian. The current
  Jacobian check only prints differences at the prior; it is not a universal
  derivative certificate or an automated pass/fail check across clipping.
- Check optimizer termination, active bounds, temperature branch usage, and
  singular directions. Compare alternatives when different parameterizations
  give nearly identical pressure residuals.
- Verify through C++ and compare pressure predictions using matching point
  identities. Use recomputed C++ production/species outputs, not inherited
  production fields in the Python fitted table.
- Vary numerical resolution and selected physical assumptions before treating
  fine structure in a map or near-surface pressure as resolved.

## 6. Running the application

### 6.1 Requirements and paths

C++ requirements: a C++17 compiler, OpenMP, CSPICE headers and static library.
The Makefile's CSPICE default is `/home/vtenishe/SPICE/cspice`.

Python requirements: NumPy and Matplotlib for data/plotting; SciPy for fitting
and image digitization. Use a Python environment containing those packages.
There is no pinned dependency environment in this directory; older notes
about packages in a particular user's `~/.local` are machine-specific.

```bash
make
# If needed, override the toolkit location:
make CSPICE=/path/to/cspice

OMP_NUM_THREADS=4 ./liouville_cops selftest.in
```

Before a new observation run, review these input settings:

The supplied SHAP5 shape path is
`/data/vtenishe/rosetta_data_end_of_mission_cops/SHAP5_Volume.bdf`, and the SPICE
root is `/home/vtenishe/SPICE/Kernels/ROSETTA/esa_2022`. Commented DSK
alternatives point to `dsk/ROS_CG_M001_OSPGDLR_N_V1.BDS` and
`dsk/ROS_CG_M004_OSPGDLR_N_V1.BDS` under that SPICE root. These resources
are not bundled within LiouvilleCOPS; configure their locations for the
machine used to run the application.

| Setting | Meaning / supplied behavior |
|---|---|
| `ObservationFile`, `ReferenceDataUnit` | Runtime pressure/time data file and its units. |
| `ShapeModel`, `ShapeModelFormat`, `RequireClosedSurface` | Surface file, format, topology requirement. |
| `KernelDirectory`, repeated `Kernel` | SPICE root and explicitly loaded kernel files. |
| `BodyFrame`, `BodyName`, `SpacecraftName/Frame`, instrument IDs | Geometry and pointing conventions. |
| `GaugeOffsetSC` | Common gauge position offset in meters, spacecraft coordinates. |
| `TemperatureRGMode/Table/TableColumn` | Constant or HK gauge temperature; column is 1-based, default 6. |
| `OutputDirectory`, `WriteKernels` | Output destination; fit requires kernels and elements. |
| `AllowExistingOutputDirectory` | `off` by default. Reserved for a trusted wrapper that already created a fresh directory; even when enabled, existing table/elements/kernels/manifest are refused. Do not use it to resume or overwrite. |
| `GeometryOnly` | `off` by default. With `on`, write exact SPICE position, BVH nearest-surface altitude, inside flag, references, and selection flags, but skip source/transport calculations and write zero model-pressure columns. Requires `WriteKernels = off`; this output is not a fitting base. |
| `FirstPoint/LastPoint/PointStep` | Subset of observation indices to evaluate. |
| `FitInterval`, repeated `ExcludeNG/ExcludeRG` | C++ selection flags; Python NG override rules still apply. |
| `ActivityModel/File`, temperature and SZA keys | Surface source described above. |
| `ActivityEnhancementModel/Degree/File/Scale` | Optional H2O-only residual. Default `none` requires degree 4. Enabled model is `real_spherical_harmonics_v1`; file degree/scale must agree with duplicated input values. See `ACTIVITY_ENHANCEMENT_FILE.md`. |
| `SurfaceBC`, subdivision, velocity and rotation keys | Transport approximation and numerical accuracy. |
| Gauge acceptance, pressure model, sensitivity keys | Instrument-response conventions. |
| `Threads` | Optional OpenMP thread count; otherwise the OpenMP environment applies. |

C++ ordinary source-data paths resolve from the process working directory,
not from the input file's directory. `OutputDirectory` uses the shared output
root; `TemperatureRGTable`, `ActivityFile`, and `ActivityEnhancementFile` use generated-product lookup.
Relative `Kernel` paths are appended to `KernelDirectory`. Unknown/unused keys in observation
mode cause an error, and conditionally irrelevant options can therefore be
rejected. A single-valued repeated key uses its last value.

For a new manifested base, the fitter reads selection intervals and generates
its verification input from the archived `input.snapshot.in`, not a mutable
original pathname in the table header. Pre-W01 legacy products have no
snapshot and retain the old fallback to the recorded input path; this is one
reason legacy use requires explicit acknowledgement.

### 6.2 Prepare a fit base and inspect the observations

After configuring the external inputs and preparing the HK table:

```bash
./liouville_cops fit_base.in
python3 tools/plot_cops_model.py output_fit_base
python3 tools/plot_cops_model_altitude.py output_fit_base
```

The supplied `fit_base.in` writes `<AMPS>/cops_output/output_fit_base/`.
The solver refuses that destination if it already exists; use a distinct
`OutputDirectory` for every experiment. A successful new run contains
`input.snapshot.in` and `run_manifest.json` in addition to its scientific
products. Validate a fitting base independently with
`python3 tools/validate_cops_run.py output_fit_base`.

Required fitting settings are `FlightTimeRotation = none`,
`WriteKernels = on`, `TemperatureModel = analytical`, and, for Python/C++
temperature agreement, `SurfaceTemperatureCorrectionFactor = 1`.
Changing geometry, acceptance, surface emission law, shadows, species
sensitivities, or fixed gauge temperatures requires a compatible new base.

A table alone is insufficient for fitting: point kernels are loaded for all
base-table rows, including rows later excluded from the misfit.

#### Apply an optional Stage-A3 H2O enhancement

Copy an observation input to a new filename, choose a fresh
`OutputDirectory`, and replace its disabled enhancement settings with:

```ini
ActivityEnhancementModel = real_spherical_harmonics_v1
ActivityEnhancementDegree = 5
ActivityEnhancementFile = examples/activity_enhancement_zero_l5.txt
ActivityEnhancementScale = 1.0
```

Then run the same executable, for example with 14 OpenMP threads:

```bash
OMP_NUM_THREADS=14 ./liouville_cops my_enhanced_run.in
```

The checked-in file is an all-zero degree-5 regression example for
`ActivityModel = shap5`; it must reproduce the baseline pressures. Replace it
with a scientifically obtained file conforming to
`ACTIVITY_ENHANCEMENT_FILE.md` to apply a nonzero residual. Degree and scale
can be omitted from the input when enabled because the file supplies them; if
duplicated, they must match exactly. A degree-4 file is the explicit
zero-length control. With `ActivityEnhancementModel = none`, degree must be 4
and file/scale settings are rejected rather than silently ignored.

An enabled run writes `ActivityEnhancementFacets.dat`, records the complete
configuration and coefficient-file hash in `run_manifest.json`, and adds an
enhancement identity header to `COPS_model.dat`. It does not change the table
columns or binary cache layout. Stage A4's Python fitter creates this file in
the controlled residual-only mode described in section 4.3.

### 6.3 Multi-source pressure comparison

The L4 table and images are external inputs, but the HK table, digitized CSVs,
and geometry tables are generated products. From a clean checkout or after
removing prior outputs, build all prerequisites before running the collector:

```bash
python3 tools/cops_hk_table.py
python3 tools/digitize_cops_figures.py
env OMP_NUM_THREADS=4 ./liouville_cops fit_base.in
env OMP_NUM_THREADS=4 ./liouville_cops m004.in
python3 tools/collect_cops_sources.py
python3 tools/collect_cops_sources.py --show
```

The final two commands are alternative invocations; both write plots. The
latter also opens windows. `m004.in` uses the ESA M004 DSK with
`GeometryOnly = on` and `WriteKernels = off` because the collector reads only
its SPICE position, inside flag, cometocentric distance, and exact
nearest-surface altitude. Its model-pressure columns are intentionally zero,
and it must not be used as a forward-pressure result or fitting base. The collector
performs a complete prerequisite check before parsing or writing anything and
reports regeneration commands when products are missing. Existing default
output files are overwritten once all inputs are available.

| Collector option | Meaning / default |
|---|---|
| `--shap5-run DIR` | SHAP5 geometry/reference table; default `output_fit_base`. |
| `--m004-run DIR` | M004 geometry table; default `output_M004`. |
| `--fit-run DIR` | Base run supplying the candidate fit pressures/flags; default `output_fit_base`. |
| `--fit-data FILE` | Output selection table; default `data/COPS_fit_data.txt`. Its plot uses the same stem plus `.png`. |
| `--fit-used FILE` | Actual fitter selection for the separate profile plot; default `output_fit_base_selection/fit_data_used.txt`. |
| `--apply-exclude-ng` | Retain base NG exclusions instead of re-enabling them. |
| `--exclude-ng-saturated` | Remove plateau-flagged NG points from the candidate table. |
| `--out FILE` | Combined JSON, default `data/COPS_EoM_sources.json`. |
| `--plots PREFIX` | Comparison PNG prefix, default `data/COPS_EoM_sources`. |
| `--show` | Also display windows. |

The collector reads HK and digitized CSVs from the output root first. It reads
the L4 source from `--sources-dir`, defaulting to `<COPS_DATA_ROOT>/sources`.
It has no raw-data-directory option; another raw location is supplied to the HK
preprocessing step or through `COPS_L2_ROOT`.

Its geometry is interpolated in time from saved run tables, with out-of-range
altitudes set to missing. It does not rerun SPICE. File/axis labels assume
SHAP5 and M004; replacing those runs with another shape does not make the
labels automatically truthful.

To regenerate the digitized CSVs from the supplied images:

```bash
python3 tools/digitize_cops_figures.py --check
```

`--dir DIR` changes the input image directory (default
`<COPS_DATA_ROOT>/sources`); `--out-dir DIR` selects the
output directory (default `data/sources` under the shared root). `--check` additionally
writes `digitized_check.png` and
`EGU2019_slide02_on_slide16_axes.png`. Without it the CSVs are still written.

### 6.4 Fit, inspect, and verify

The working fitter is included. Preview a particular
selection without running the optimizer:

```bash
python3 tools/fit_cops.py output_fit_base \
  --data-only --exclude-ng-saturated --out output_selection_new
```

Check derivatives separately:

```bash
python3 tools/fit_cops.py output_fit_base \
  --check-jacobian --lambda 1e-3 --out output_jacobian_check
```

Example fit with an explicit regularization weight and selection:

```bash
python3 tools/fit_cops.py output_fit_base \
  --lambda 3e-3 --exclude-ng-saturated --ordered-temperature \
  --out output_fit_new
```

The lambda here is illustrative, not a scientifically established choice.
Add `--apply-exclude-ng` if the late NG interval should also be blocked.
Review temperature bounds in ordered mode as explained above.

To use a selection table:

```bash
python3 tools/fit_cops.py output_fit_base \
  --data data/COPS_fit_data.txt --lambda 3e-3 \
  --out output_fit_table_new
```

For multiple regularization weights:

```bash
python3 tools/fit_cops.py output_fit_base \
  --lambda 1e-4 1e-3 3e-3 1e-2 1e-1 \
  --exclude-ng-saturated --jobs 2 --out output_lcurve_new

python3 tools/plot_lcurve.py output_lcurve_new \
  --save output_lcurve_new/lcurve.png
```

Before judging the L-curve, inspect each fit's termination status and
temperature/CO2 bounds. Increase resources/evaluations or change the model
only on the basis of those diagnostics.

Verify a newly generated fit through C++, then compare:

```bash
./liouville_cops output_fit_new/cops_fit.in

python3 tools/validate_cops_run.py \
  --fit output_fit_new --verification output_fit_new_verify

python3 tools/plot_cops_model.py \
  output_fit_base output_fit_new output_fit_new_verify \
  --ref output_fit_new --save COPS_fit_comparison.png
```

This C++ rerun recomputes the physics and species outputs. Its selection
flags still follow the copied base input, not necessarily the Python
selection used in the fit. Use `fit_data_used.txt` to identify the actual
misfit data.

The plotting command deliberately selects the Python fit with `--ref`: that
table carries the fitter's final `useNG`/`useRG` masks. Selecting the base or
C++ verification instead uses the copied input exclusions for hollow/filled
markers and can make re-enabled observations look unused.

For an exact-rotation check, make a separate copy of the newly generated
input, set `FlightTimeRotation = exact` and a distinct output directory,
and rerun. Do not fit its version-2 kernels with the current fitter.

The former historical `cops_fit.in` files were removed with their result
directories. Newly generated fit inputs use the shared output-path handling.

### 6.5 Fitter option reference

| Option | Default | Effect |
|---|---|---|
| Positional `run` | Required | Base directory with model table, elements, and kernels. |
| `--out DIR` | `<run-name>_fit` | New output directory under the shared root. Any existing path is refused; there is no overwrite/resume mode. |
| `--degree N` | 4 | Frees basis groups through nominal degree N, intended range 0–4. |
| `--enhancement-degree N` | 4 | Maximum normalized residual degree, 4–8; degree 4 is the zero-length control. Degrees 5–8 require all three fixed-physics switches below. |
| `--fix-legacy-activity` | off | Holds all 25 legacy H2O coefficients fixed and enables the controlled residual/control parameterization. |
| `--enhancement-scale VALUE` | derived | Optional assertion of the required area-weighted legacy-linear RMS scale; it cannot override that scale. |
| `--activity-template BUNDLE NAME` | none | Adds one mesh-bound, nonnegative H2O template amplitude; repeatable. At enhancement degree 4 it supports jointly free legacy activity, temperatures, and `kappa`; degrees 5–8 retain their fixed-physics contract. Cliff columns require repaired gravity provenance. |
| `--lambda L [L ...]` | 0 | Regularization weight(s); multiple values produce an L-curve. |
| `--jobs N` | min(number of lambdas, 12) | Concurrent L-curve fits. |
| `--fix-temperature` | off | Keeps all three base temperature parameters. |
| `--ordered-temperature` | off | Uses nonnegative temperature increments; see bounds caveat. |
| `--fix-kappa` | off | Keeps CO2 multiplier at 1. |
| `--tbounds LO HI` | 30 400 | Independent T bounds, or Tmin bounds in ordered mode. |
| `--max-nfev N` | 1000 | Maximum function evaluations. |
| `--ftol`, `--xtol`, `--gtol` | 1e-10 each | Cost, step, and optimality stopping tolerances. |
| `--check-tolerance` | 1e-6 | Allowed prior-model discrepancy from the C++ table. |
| `--check-jacobian` | off | Prints analytic/finite-difference comparison at the prior, then stops. |
| `--data FILE` | Base table | Replaces observed pressures and selection flags. |
| `--apply-exclude-ng` | off | Applies the base input's NG exclusion intervals. |
| `--exclude-ng-saturated` | off | Removes HK plateau-flagged NG points. |
| `--hk-table FILE` | `data/COPS_HK_20160930.txt`, output-root first | HK used to identify saturation; not a replacement for saved base-run T_RG. |
| `--data-only` | off | Writes selected observations/plots after model checks, then stops. |
| `--allow-legacy-base` | off | Explicitly accept a base without `run_manifest.json` after strict cache validation; the provenance limitation is recorded. Prefer regenerating the base. |

### 6.6 Change the surface boundary condition or shape and rebuild the fitting base

`tools/rebuild_cops_boundary.py` automates the forward recalculation required
to change `SurfaceBC` or the surface mesh. It is a Python 3.8+ standard-library driver for the
existing C++ executable, not another transport solver or a binary-cache
converter. It does not compile the executable or run the optimizer. Build
with the Makefile first if necessary, and configure the external paths in
the source input as described above.

#### What must be recalculated?

At each visible surface sub-element, the static geometric weight is:

```text
flat:  w = (dS / r^2) cos(theta) / pi
rough: w = (dS / r^2) / (2 pi)
```

`theta` is the emission angle relative to the outward local surface normal.
The same total emitted flux is distributed differently in angle. Since the
angle varies across sub-elements and observations, there is no single
flat-to-rough scale factor for a saved face kernel or pressure curve. The
version-1 cache contains sums, not all sub-element angles/weights needed to
reweight those sums exactly. Reevaluate the geometry/integration instead.

| Quantity | Effect of changing only `SurfaceBC` |
|---|---|
| `Kn`, `KPhi`, `KV` for each observation/face | Recompute: all contain the boundary-dependent angular weight. These are saved in the version-1 point files. |
| `KnBelowDeck`, `NGBelowDeckFraction` | Recompute the NG directional diagnostic. `KnBelowDeck` is used during the forward run but is not saved in version-1 point records. |
| NG densities, RG fluxes, species and total pressures | Recompute using the new weights and the retained source/temperature parameters. |
| Mesh geometry, basis, base activity, trajectory, pointing, Sun geometry, shadow state, observed pressures, gauge temperatures, selection intervals | Do not intrinsically change when only the angular emission law changes. The driver still reruns the complete solver to keep the products consistent. |
| `elements.bin` | Not intrinsically boundary-dependent; regenerated alongside the new kernels so face ordering/source data cannot be accidentally mixed between runs. |
| Total production `Q_H2O`, `Q_CO2` | Unchanged for the same source parameters: both angular laws integrate to the same prescribed surface flux. A subsequent refit can change production. |
| Best-fit activity coefficients, temperatures, CO2 multiplier, residuals, regularization diagnostics | Not transferable as a new optimum. Refit with the rebuilt base if the goal is a fitted model for the new boundary condition. |

The driver regenerates all observation indices selected by
`FirstPoint`/`LastPoint`/`PointStep`, including observations excluded from the
misfit. This is necessary because the fitter checks the full base model.

#### How temperature enters

For static transport (`FlightTimeRotation = none`), `Kn`, `KPhi`, and `KV`
contain geometry, angular emission, and gauge acceptance but no surface
temperature. For species s and face e:

```text
beta_s,e = sqrt(m_s / (2 k T_surface,e))
n_NG,s   = sum_e F_s,e (sqrt(pi)/2) beta_s,e Kn_e
Phi_RG,s = sum_e F_s,e [KPhi_e + (V.b) (sqrt(pi)/2) beta_s,e KV_e]
```

Surface temperature sets the emitted speed distribution. At fixed emission
flux, increasing it reduces the NG density contribution as `1/sqrt(T)`;
it also changes the spacecraft-velocity part of RG flux. The leading static
RG flux term has no surface-temperature dependence. The implemented source
law `F_s,e` is prescribed from activity and illumination, not computed from
a temperature-dependent sublimation/energy-balance solution.

The fitter can therefore vary the analytical surface-temperature parameters
without rerunning SPICE or ray tracing: it recomputes the temperature and
speed factors on the cached geometry. Changing `SurfaceBC`, in contrast,
requires the new geometric kernels. The rebuild tool preserves all supplied
temperature settings; it does not derive or fit a new temperature merely
because the angular emission law changed. Use the new base in a fresh fit
to determine whether different fitted temperatures are preferred.

Instrument temperatures are separate from surface temperature. NG pressure
conversion uses `TemperatureNG`; RG conversion uses `TemperatureRG` or the
HK-derived `T_RG`. In `equilibrium_sphere` mode, before any sensitivity
correction, `p_RG = Phi_RG sqrt(2 pi m k T_RG)`. The other implemented RG
conversions also scale with `sqrt(T_RG)` but have different constants. These
gauge temperatures remain fixed during fitting. Changing their input values
requires a new consistent base table for the fitter's initial check, even
though they do not change static angular kernels.

In `FlightTimeRotation = exact`, speed and temperature also affect emission
time, moving-surface geometry, and path obstruction. Its version-2 kernels
are temperature-dependent. The current fitter cannot read them, so this
rebuild tool explicitly rejects exact-rotation inputs rather than silently
converting them to static transport. Use a separate direct C++ forward run
for exact-rotation comparisons (section 6.4).

#### Commands and options

Preview the generated input and command without creating files or running
the solver:

```bash
python3 tools/rebuild_cops_boundary.py \
  --input fit_base.in --surface-bc rough \
  --out output_fit_base_rough --threads 4 --dry-run
```

Then run the same command without `--dry-run`:

```bash
python3 tools/rebuild_cops_boundary.py \
  --input fit_base.in --surface-bc rough \
  --out output_fit_base_rough --threads 4
```

Select the converted 199k SHAP5 surface as an explicit simulation parameter
while retaining the flat boundary law:

```bash
python3 tools/rebuild_cops_boundary.py \
  --input fit_base.in --surface-bc flat \
  --shape-model /data/vtenishe/rosetta_data_end_of_mission_cops/shape_models/SHAP5_PDS_V2/cg_spc_shap5_199k_cart.bdf \
  --shape-model-format nastran \
  --out output_fit_base_shap5_199k --threads 16 --dry-run
```

After reviewing the generated input, repeat without `--dry-run`. A shape
change recomputes the facet table, BVH, signed altitude and inside/outside
classification, illumination and shadows, visible facets, pressures, and
every point kernel. The resulting directory is the base passed to the fitter.
Shape is a controlled forward-model input, not a quantity optimized by
`fit_cops.py`.

| Option | Default | Meaning |
|---|---|---|
| `--input FILE` | `fit_base.in` | Source observation input. Use the actual input for the geometry/source you want, not just an old output directory. |
| `--surface-bc flat\|rough` | Required | Target angular emission law. |
| `--out DIR` | Required | New destination; any existing path, including an empty directory or a dangling symlink, is refused. No overwrite/resume switch. |
| `--threads N` | Preserve input/environment | Positive OpenMP thread count inserted into the generated input; overrides any original `Threads`. |
| `--shape-model FILE` | Preserve source `ShapeModel` | Surface mesh for this simulation. Relative paths resolve from the application root and are stored as absolute paths in the generated input. |
| `--shape-model-format nastran\|dsk` | Preserve source value or `nastran` default | Select the matching C++ mesh reader. The driver deliberately does not infer format from a filename suffix. |
| `--executable FILE` | Application-root `liouville_cops` | Existing compiled solver to execute; no automatic rebuild. |
| `--dry-run` | off | Validate the source configuration and print the proposed input/command; create nothing and skip executable/external-file availability checks. |

Relative `--out` paths are interpreted under the shared output root. Other
relative CLI paths and ordinary source-data paths use the application root,
even when launched elsewhere. A generated input can also be found in the
output root; HK and activity files use generated-product lookup. Relative `Kernel` entries
are interpreted under the effective `KernelDirectory`. Generated data paths
are absolute; there is no shell expansion inside input values. Whitespace,
quotes, `#`, `$`, and `~` in generated input paths are rejected for compatibility
with the C++ input and Python fit-header readers.

The source must be an observation-mode, static, analytical-temperature input
with `SurfaceTemperatureCorrectionFactor = 1`. Missing keys use the same
defaults as the C++ reader for those checks. Unsupported modes are rejected,
not silently changed. The driver preserves the source file, comments,
repeated kernels/exclusion intervals, source parameters, and gauge settings.
It overrides `SurfaceBC`, `OutputDirectory`, `WriteKernels = on`, and
optionally `Threads`, `ShapeModel`, and `ShapeModelFormat`, apart from
normalizing paths. Point-range settings must
be nonnegative integral indices with a positive integral step.

Real execution checks the executable and explicitly referenced data files
before creating the destination. SPICE coverage, file contents, mesh validity,
and unknown/unused input keys are still checked by the C++ solver; a dry run
does not establish that these are valid. A full rebuild costs approximately
another full forward run and its disk storage, not an inexpensive file edit.

#### Products and failure handling

All products go under `--out`:

| File | Meaning |
|---|---|
| `source.in` | Snapshot of the original input text. |
| `boundary.in` | Generated input with absolute data/output paths; the new model table records this filename. Preserve it for fitting. |
| `COPS_model.dat`, `elements.bin`, `kernels/point=*.bin` | Newly computed static forward base. |
| `rebuild.log` | Combined solver stdout/stderr, also streamed to the terminal. |
| `rebuild.json` | Command, working directory, source-input SHA-256, start/end times, status, and on success the verified observation count. |

Success requires a zero solver exit status, the correct boundary/static table
header, the full expected sequence of observation indices (allowing the C++
last-index clamp), and matching version-1 element/kernel headers and file
lengths. These are completeness/format checks, not proof of numerical
convergence or scientific validity. The source-input hash is not a full
provenance hash of the mesh, executable, or mission kernels.

On solver failure, interruption, or output-check failure, the script returns
nonzero and retains diagnostic/partial products with `status = failed` in
`rebuild.json` once that metadata file has been initialized. Do not fit a
failed directory. Retry with a new destination after fixing the problem.
Existing runs and source inputs are never overwritten by the driver. Do not
manually rerun `boundary.in` into its existing destination unless you intend
to bypass that protection. Moving the new directory also requires updating
the absolute input/output references; keep it in place for fitting.

After successful rebuilding, use the included `tools/fit_cops.py`
(section 1), inspect selection, then fit:

```bash
python3 tools/fit_cops.py output_fit_base_rough \
  --apply-exclude-ng --exclude-ng-saturated --data-only \
  --out output_rough_selection

python3 tools/fit_cops.py output_fit_base_rough \
  --apply-exclude-ng --exclude-ng-saturated --ordered-temperature \
  --lambda 3e-3 --out output_rough_fit
```

Use `--surface-bc flat` with a different `--out` to rebuild the flat case.
Compare the cases with identical data-selection and fitting options. The
example regularization weight is illustrative; the rebuild tool does not
select it or transfer any old best-fit parameters automatically.

#### Tests

Standard-library tests exercise input preservation, repeated/path settings,
dry-run behavior, refusal to overwrite, unsupported configurations, process
failure, and incomplete/corrupt cache detection using a fixture executable:

```bash
python3 -m unittest discover -s tests -p 'test_rebuild_cops_boundary.py' -v
```

An optional integration test uses the local compiled executable and the
external resources referenced by `fit_base.in`. It performs three one-point
SHAP5 runs in temporary directories: flat, rough, and flat with doubled
analytical surface temperatures. It disables spacecraft velocity and uses
full unweighted NG acceptance to isolate the dependencies above. It checks
that changing the boundary changes the kernels, doubling surface temperature
does not, NG density/pressure scales by `1/sqrt(2)`, and RG flux/pressure and
production remain unchanged in this controlled temperature test. Temporary
test products are removed automatically; scientific products outside the
temporary directories are not modified.

```bash
env COPS_RUN_INTEGRATION=1 python3 -m unittest discover \
  -s tests -p 'test_rebuild_cops_boundary.py' -v
```

## 7. Results and file formats

### 7.1 Forward `COPS_model.dat`

The header names the columns; do not rely on fixed column positions across
historical versions.

| Columns | Meaning |
|---|---|
| `iPoint UTC et` | Observation index, UTC string, SPICE ephemeris time. |
| `x y z altitude inside rh_AU` | Gauge coordinates (m), unsigned surface distance (m), inside flag, heliocentric AU. |
| `VdotB_RG nVisible` | RG velocity projection (m/s) and number of contributing faces. |
| `Q_H2O Q_CO2` | Whole-surface production, molecules/s. |
| `n_H2O n_CO2` | NG-response-weighted densities, m^-3. |
| `PhiRG_H2O PhiRG_CO2` | Accepted incoming RG fluxes, m^-2 s^-1. |
| `pNG_H2O pNG_CO2 pNG pRG_H2O pRG_CO2 pRG` | Species contributions and summed pressures, Pa. |
| `refNG refRG` | Observation pressures converted to Pa. |
| `NGBelowDeckFraction` | Fraction of weighted NG density from the code's geometrical below-deck direction test. |
| `T_RG useNG useRG` | Sensor temperature used (K), C++ data-selection flags. |

Values at inside-body points are retained in the table but excluded from
fitting. They are not valid external-coma model predictions.

When Stage A3 is enabled, `ActivityEnhancementFacets.dat` contains one row per
mesh face with

```text
face x y z area legacy_linear_h2o enhancement_h2o
total_linear_h2o activity_h2o clipped
```

Coordinates and area are in meters and square meters; activity columns are
dimensionless. `activity_h2o=max(0,total_linear_h2o)`. The diagnostic is not
created when the enhancement is disabled.

### 7.2 Binary caches

C++ writes native `int` and `double` values with separate writes, not
a serialized padded struct. The current Python reader assumes 32-bit
integers, 64-bit floats and little-endian record fields. Treat these as
implementation-specific caches, not a portable archival format.

`elements.bin`:

```text
int32 version = 1, nFaces, nBasis = 25
for each face:
  float64 centroid[3], normal[3], area, basis[25], activity_H2O, activity_CO2
```

Each face record is 34 doubles = 272 bytes. Its face order must match every
kernel file. The 25 saved basis values always remain the legacy basis. With an
enabled Stage-A3 residual, `activity_H2O` is the final applied activity; use
`ActivityEnhancementFacets.dat` to distinguish its linear legacy and residual
components.

Version-1 `kernels/point=NNNN.bin`:

```text
int32 version = 1, iPoint
float64 et, VdotB_RG, xGauge[3], SunDirection[3], r_h_m, altitude_m
int32 nElements
for each contributing face:
  int32 face
  float64 Kn, KPhi, KV, cosSZA
  int32 shadow
```

The header is 92 bytes and each record 40 bytes on the expected platform.
Only contributing faces are stored. Face indices refer to the base mesh;
the number of records can be zero.

Version 2 adds Omega[3] to the header and an `nSpecies` integer before
`nElements`. Each record contains face, Gn[nSpecies], GPhi[nSpecies],
cosSZA, temperature, and shadow. These temperature-dependent response
kernels cannot be used by the current static fitter.

The binary headers retain their original compact format. New runs add
`run_manifest.json`, which hashes inputs and outputs and identifies the kernel
set, and the fitter independently validates every table/cache relationship.
Do not mix kernels, elements, and tables manually: the manifest and validator
will reject a detected mismatch. Pre-W01 runs remain explicitly identifiable
as unmanifested legacy products.

### 7.3 Fit products

Every fitter invocation that reaches selection export writes:

| File under `--out` | Meaning |
|---|---|
| `fit_data_used.txt` | Actual selected pressures from the fit problem; unused gauge entries are -1. |
| `fit_data_used.png` | Selected points and other reference data versus UTC/altitude. |
| `fit_profiles_used.png` | Selected profiles only; lines break at unused points or gaps over 90 seconds. |

Each completed optimization writes:

| File | Meaning |
|---|---|
| `fit_result.txt` | Selection/options, objective, convergence diagnostics, parameters, approximate sigmas, activity diagnostics. |
| `ActivityFit.txt` | H2O fitted coefficients, scaled CO2 coefficients, unchanged CO row. |
| `ActivityEnhancementFit.txt` | Strict normalized-basis residual coefficients and baseline binding; enhancement mode only. |
| `cops_fit.in` | Input for a C++ rerun into `<fit directory>_verify`. |
| `COPS_model.dat` | Schema-v2 fitted aggregate pressure predictions; unrecomputed diagnostics are `nan`. |
| `fit_manifest.json` | Versioned input/cache/selection/prior/settings/output identities and independent process, optimizer, C++-agreement, and review states. |
| `fit_report.json` | Canonical A5 report identity and the exact scientific-choice payload used to form it. |
| `activity_facets.dat`, `activity_degree_power.dat`, `activity_region_diagnostics.dat` | Universal facet fields, residual power by degree, and all/seen/unseen area-weighted diagnostics. |
| `activity_maps.png` | Separate legacy baseline, fitted residual, and clipped-total maps. |
| `pressure_residuals.dat` | Per-gauge selection/reason/saturation flags, observations, predictions, and selected-only log residuals. |
| `pressure_comparison_time.png`, `pressure_comparison_altitude.png`, `pressure_residuals.png` | Self-identifying fit comparisons and residual diagnostics. |

**The Python-fitted model table is not a complete recomputed forward run.**
Its aggregate `pNG`/`pRG` are recomputed outside the nucleus. Values that the
cached Python model does not recompute—`Q_H2O/Q_CO2`, species densities,
species fluxes, species pressures, and `NGBelowDeckFraction`—are explicitly
written as `nan`, never copied from the base. Geometry, observations, gauge
temperature, and selection fields retain their documented input meanings.
Use the C++ verification run for complete physical diagnostics, then register
the pressure comparison with `tools/validate_cops_run.py`.

For multiple lambdas, each result goes into `lambda=<value>/`;
`lcurve.txt` and the selection exports are in the parent directory.
The generated C++ input preserves base exclusions and does not encode an
arbitrary replacement pressure table or all Python selection switches.

### 7.4 Collector products

| Default output | Content |
|---|---|
| `data/COPS_EoM_sources.json` | Seven datasets, conversions, source metadata, geometry references, and historical event notes. |
| `data/COPS_EoM_sources_time.png` | NG/RG source comparison versus UTC, with candidate fit points marked red. |
| `data/COPS_EoM_sources_altitude.png` | Source comparison versus altitude; M004 for time-tagged sources, original axes for EGU. |
| `data/COPS_fit_data.txt` | Candidate fit-selection table, suitable for `fit_cops.py --data`. |
| `data/COPS_fit_data.png` | Selected/excluded observations and exclusion reasons. |
| `data/COPS_fit_profiles.png` | Actual fitter-selected profiles from `--fit-used`, if that file exists. |

The actual-fit profile plot is independent of the collector's newly constructed
selection. Changing collector options alone does not update `--fit-used`.
Point it to the matching fitter output to avoid displaying two different
selections as if they were the same.

With custom output paths, the first two PNGs are `PREFIX_time.png` and
`PREFIX_altitude.png`; the selection PNG follows `--fit-data`'s stem.
`COPS_fit_profiles.png` is written beside the plot prefix. Missing
`--fit-used` skips that plot, but does not delete an older existing PNG.
Parent directories for custom collector outputs are created automatically.

### 7.5 Model and L-curve plotters

```bash
python3 tools/plot_cops_model.py output_fit_base output_fit_base_fit
python3 tools/plot_cops_model.py output_M004 --x altitude --save m004.png
python3 tools/plot_lcurve.py output_fit_base_lcurve
```

The model plotter accepts run directories or explicit table filenames,
`--x time|altitude`, `--save FILE`, `--theme light|dark`,
and `--ref NAME`. NAME is a loaded run's directory basename (or its
suffix after `output_`), not an arbitrary unloaded directory path.
With no runs specified it loads `output*/COPS_model.dat` one level below
the output root; if none are found, it scans the application root for legacy
runs. It does not merge those scans. Nested lambda directories must be specified
explicitly.

For a schema-v2 Python fit, the plotter requires and verifies the matching
fit manifest, declared table/report hashes, and the report/profile/selection/
basis identities. It rejects mismatches. Its labels include profile, degree,
lambda, and report ID; unmanifested historical tables are explicitly labelled
legacy. Manifested C++ tables are also checked against their run-manifest
table hash and labelled with their activity identity. Saved comparison PNG
metadata lists all loaded fit report identities (or a legacy marker).

Its panels show combined NG/RG, NG, RG, and NG/RG. Reference points with
use=0 are hollow. Without `--ref`, it chooses the most recently modified
loaded table with selection flags. Set it explicitly for reproducible
geometry/selection comparisons.

The L-curve plotter accepts a directory or `lcurve.txt` and optional
`--save FILE`. Open circles identify evaluation-limit results; a marked
curvature point may still be unconverged.

## 8. What the historical fits do and do not establish

A single-fit report inspected before output cleanup used 839 NG and 1,146 RG
points and reported:

- RMS log residuals 0.0605 (NG) and 0.0637 (RG).
- Termination at the 200-evaluation limit, not optimizer convergence.
- Tmax about 107.93 K, T75 about 75.04 K, Tmin about 115.14 K,
  and kappa about 0.5043.
- A weighted-Jacobian condition number about 1.1e16.

Because Tmin exceeds both other temperature parameters, the implemented
temperature floor makes this reported solution effectively constant at
Tmin. Zero sigmas printed for Tmax/T75 do not demonstrate that those two
parameters were determined. This older report does not contain all fields
written by the current fitter; do not retroactively assign it
current default options.

The removed historical L-curve marked lambda = 3.162e-3 as its maximum-curvature point.
That fit reports status 0 at 1,000 evaluations, RMS log residuals about
0.0668/0.0676, D = 0.95015, and p_eff about 15.876 of 29 parameters.
Its T75 is at the 30 K lower bound. About 64% of its clipped activity
integral lies on faces classified as unseen, which cover about 28% of the
surface area. These diagnostics show the importance of prior/selection
and model assumptions; they do not establish a unique global production map.

The removed fits' selection exports contained 839/1,146 points. A later
selection table inspected before cleanup contained 930/1,146. Re-running with
today's defaults therefore does not reproduce those historical fits' selection.

These recorded historical diagnostics do not by themselves validate detector
calibration, missing-data reconstruction, the physical source prior, or
collisionless conditions. The underlying result files are no longer local.

### 8.1 Publication and interpretation checklist

The checks below are project-specific recommendations, not a claim that the
historical results passed them or a statement of any journal's policy.
They separate a useful methods/exploratory result from a strong physical
inference. Implementation tasks and decision gates are in
[ROADMAP.md](ROADMAP.md).

1. **Resolve convergence before interpreting a preferred fit.** The historical
   single fit stopped at 200 evaluations, with a condition number near
   1e16 and effectively constant temperature caused by its Tmin floor.
   The marked L-curve solution stopped at 1,000 evaluations and has T75
   at its lower bound. These are numerical experiments, not demonstrated
   converged estimates. More evaluations alone may not fix inactive or
   degenerate parameters; compare lower degree, fixed/free temperature and
   CO2, and multiple initial guesses with the same prior and observations.

2. **Separate pressure agreement from map uniqueness.** Different activity,
   temperature, and CO2 combinations can predict similar pressures.
   Approximately 64% of the selected historical solution's clipped activity
   integral lies on faces classified as unseen; that is not a gas-production
   fraction. Report sensitivity of global maps/production to degree, prior,
   regularization, and physical assumptions. A visible-face mask is not a
   quantitative spatial-resolution assessment.

3. **Define the observation product and selection policy.** Retain original
   units, calibration conventions, interval endpoints/midpoints, missing
   values, saturation and range-change flags, and the final selection table.
   Resolve or explicitly document the unmatched reference/HK times and
   missing-versus-positive discrepancy in section 3.2. Do not silently
   replace missing observations with model predictions. Label digitized,
   derived, and reconstructed sources separately; different representations
   of the same underlying measurements are not independent validation.

4. **State the actual statistical objective.** Equal-gauge normalized log
   least squares is not a reported-instrument-error likelihood. Higher-cadence
   intervals have more influence per unit time; serial correlations and
   saturation are not modeled. Motivate any alternative weighting/noise model
   with evidence and show its effect. Do not equate residual scatter with
   parameter uncertainty or introduce arbitrary errors merely to obtain a
   preferred reduced statistic.

5. **Report uncertainty at the level justified by the analysis.** Printed
   local sigmas do not include calibration, shape, pointing, telemetry-to-gauge
   temperature mapping, or model discrepancy. Bounds and null directions need
   explicit treatment. Profile analyses, sensitivity ensembles, blocked
   resampling, or posterior methods are possible follow-up work; they are
   not current validated outputs and require a defensible observation model.

6. **Distinguish numerical verification from physical applicability.**
   Passing sphere/cache tests does not establish the collisionless assumption
   near the surface or validate NG/RG calibration. Assess mesh/integration
   convergence, acceptance and sensitivity choices, pressure conversion,
   head positions, sensor temperature, and static versus flight-time rotation.
   Thermal inertia and a coupled sublimation/energy balance are absent;
   fitted temperatures may compensate for errors in prescribed physics.
   Any collisionality estimate must identify its cross sections, density,
   length scales, and uncertainty rather than assume collisionless validity.

7. **Use scientifically appropriate products.** C++ verification supplies
   recomputed production and species diagnostics. Python fitted-table Q
   columns remain inherited from the base; activity ratios omit illumination
   factors and are not production ratios. Preserve `fit_data_used.txt`
   separately because the C++ verification run need not reproduce Python's
   selection. Publish modeled values as modeled, not newly measured pressures.

8. **Test prediction and recoverability.** Consider contiguous held-out time
   intervals, with training-only parameter/lambda selection, and synthetic
   recovery experiments with known sources. Avoid assuming adjacent randomly
   split observations are independent. Same-model synthetic recovery tests
   check implementation/identifiability under that model, not validity of its
   physical assumptions. Include deliberate model mismatch to probe bias.

9. **Archive a reproducible experiment.** Retain source and executable identity,
   compiler/flags, dependency versions, inputs, kernel/mesh/data hashes,
   selection exports, bounds, lambda, prior, optimizer reports, and verified
   physical outputs. Document the exact implemented activity basis and its
   coefficient ordering, not just "degree-4 spherical harmonics". Provide
   references and access instructions for external resources when they cannot
   be redistributed. Keep rejected/nonconverged alternatives in the audit
   trail rather than reporting only a preferred curve.

The recommended progression is: data/provenance integrity, converged and
diagnostically acceptable fits, C++ verification, identifiability/sensitivity
analysis, uncertainty and physical-applicability assessment, then claims
limited to what survives those tests. A methods paper can transparently
describe current approximations and exploratory results without claiming a
unique global activity, temperature, or composition solution.

## 9. Verification performed for this documentation

The initial documentation audit used only files within this application:

- Inspected local C++ implementation, current scripts, input files, and the
  fitter source (then available as RCS revision 1.9).
- Inspected local L4 label/format and compared the saved reference pressures
  with the processed HK table using stop-time matching.
- Counted current fit selections and inside-body points from actual tables.
- Ran `OMP_NUM_THREADS=4 ./liouville_cops selftest.in`:
  **ALL PASSED, 0 failures**. This tested the existing executable;
  `make -n` reported it up to date, but no fresh rebuild was performed.
- Loaded the Python fitter in memory and evaluated base points
  0, 300, 600, 900, 1050, 1145 using their local kernels. Maximum relative
  differences from saved C++ pNG/pRG were approximately 1.93e-11 and
  5.14e-10.
- Checked all 29 analytic Jacobian columns on those six points at the prior,
  including lambda=1e-3 regularization, against central finite differences.
  Worst relative column discrepancy was about 5.87e-5, for a11.
  This is a sampled check, not a full-domain proof.
- For Stage A3, compiled and ran strict parser/application tests covering the
  degree-4 control, zero and single-mode degree-5--8 fields, H2O-only scope,
  frame and legacy-map identity, completeness, nonfinite values, and final
  clipping. A fresh one-point, 100,000-facet SHAP5 run with the all-zero
  degree-5 example matched the disabled run and the previously manifested
  baseline data row exactly. Its facet diagnostic had all 100,000 rows and
  its JSON manifest recorded the expected 11 modes, basis ID, scale, file
  identity, and output hash.

The built-in tests cover BVH versus brute force, flat/rough uniform-sphere
density and flux, zero-rotation consistency, a rotating uniform sphere on
the axis, the static spacecraft-velocity term, and convex-body shadows.
Off-axis rotation changes are reported rather than checked against an
independent analytic reference. The largest static sphere relative error
observed was about 4.18e-4; the configured sphere pass threshold is 0.01.

No new full descent, full optimization, raw-archive comparison, or independent
instrument-calibration validation was performed for this documentation.
Historical full-run consistency errors are recorded in
`output_fit_base_fit.log`; they are historical evidence, not fresh results.

The subsequent shared-output-path update was checked with a fresh
`make -j4` build and `OMP_NUM_THREADS=4 ./liouville_cops selftest.in`
(all passed, zero failures). The following suite passed all 28 tests:

```bash
env COPS_RUN_INTEGRATION=1 python3 -m unittest discover -s tests -v
```

These include temporary-directory digitization, collection, plotting, a
two-point C++ base, data selection, a one-evaluation fit, its C++ verification,
and three one-point boundary/temperature checks. Solver integration uses the
configured external shape/observation/SPICE inputs; it requires them to be
installed. The collector test also checks that external source files are unchanged.
Without `COPS_RUN_INTEGRATION=1`, three integration tests are skipped.
These are routing/regression checks, not a new full-descent fit or calibration.

Stage A6 was subsequently run with the current executable and all integration
tests enabled:

```bash
python3 tools/run_a6_validation.py --integration --threads 2 \
  --report /tmp/liouville_cops_a6_validation_20261018.json
```

All 66 Python unit, synthetic, and integration tests passed with no skips, and
the following C++ self-test passed with zero failures. The enhanced Python fit
and C++ replay differed by at most `5.360e-8` for NG and `5.368e-8` for RG,
inside the predeclared `1e-6` threshold. Synthetic recovery, Jacobian,
regularization, and unsupported-mode measurements are recorded in
`A6_VALIDATION.md`. These results establish the tested implementation
contracts, not a full-descent fit or independent scientific validation.

## 10. Limitations and priorities for improvement

These are implications of the inspected implementation and data, not changes
already implemented.

1. **Make each experiment reproducible.** Keep the precise input, dependency
   versions, source revision, selection table, mesh/kernel hashes, and
   original data provenance with each run. Remove stale absolute paths
   from relocated rerun inputs. The working fitter is now restored; version
   and preserve its path-handling updates along with the helper module.
   Reject duplicate table indices and mismatched point/file identities.

2. **Resolve the observation/HK discrepancies.** Investigate the two unmatched
   times and the missing-versus-positive NG point before calling the tables
   interchangeable. Retain original precision and explicitly distinguish
   midpoint, endpoint, and averaging-interval conventions.

3. **Choose and record a defensible data-selection policy.** Test the effects
   of saturated-point exclusion, the final NG range change, and the changing
   cadence. A plateau observation may require a censored likelihood rather
   than an ordinary log residual. Uncertainty/quality or time-interval
   weighting would require a stated noise model. Harmonize collector/fitter
   selection logic, especially overlapping exclusion and fit intervals.

4. **Separate identifiable temperature parameters from inactive ones.**
   Inspect the fraction of contributing elements on each temperature branch,
   active bounds, singular directions, and sensitivity to starting values.
   Ordered temperatures help enforce a chosen physical ordering, but their
   present upper-bound implementation needs care. Consider smoother,
   physically motivated temperature/source laws only with validation.

5. **Assess nonuniqueness before interpreting a source map.** Try lower
   degrees, several starts, fixed versus free temperature/CO2, alternative
   lambdas, and held-out time intervals. Global basis functions couple
   unobserved faces to the inferred map. Compare actual production from
   C++ reruns, not just activity integrals or inherited fit-table Q columns.
   Document/verify the exact angular basis before changing its convention.

6. **Finish optimizer diagnostics before selecting an L-curve corner.**
   Revisit evaluation-limit solutions and bound hits. Curvature selection
   should account for convergence and sensitivity to lambda sampling.
   A low pressure misfit is insufficient evidence of a stable parameter fit.

7. **Validate gauge physics and frames.** Establish the sensor-temperature
   mapping, species sensitivities, RG conversion, NG head-on factor,
   angular acceptance, and physical head locations from instrument evidence.
   Test nonzero spacecraft velocity with rotation/offsets in consistent
   frames. The existing tests do not establish those calibration choices.

8. **Check numerical and physical applicability near the surface.**
   Vary shape resolution, subdivision threshold and maximum level; assess
   centroid-based shadows and near-surface positions. For rotation, vary
   speed nodes/range, iteration/path settings, and evaluate retarded shadows.
   Estimate collisionality independently before extending a collisionless
   model into conditions where collisions may matter.

9. **Improve uncertainty reporting.** Distinguish null directions from small
   uncertainty, account for bounds/clipping, and test sensitivity to
   correlated errors and systematic response/geometry errors. Bootstrap,
   profile-likelihood, or posterior sampling would need a defensible
   observation likelihood and additional computation.

10. **Improve output integrity and generality.** Recompute or explicitly mark
    all fitted-table fields and headers; avoid presenting base production as
    fitted production. Make raw/processed input paths configurable in the
    collector, validate date/shape assumptions and empty saturation sets,
    and ensure its actual-fit plot uses the intended selection. Reduce
    repeated kernel loading and measure memory requirements before scaling
    to finer meshes or many workers.

The collisionless transport equations, local fitting objective, and file
handling can be inspected and tested here. The physical correctness of
instrument calibration, the origin of unresolved data discrepancies, and
the global uniqueness of inferred activity cannot be supplied by assumptions
where this directory lacks the necessary evidence.
