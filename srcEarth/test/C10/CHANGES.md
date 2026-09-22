# Changes from the attached C10 scaffold

## Gzip reference update

- The production comparison file is now
  `reference_C10_poes_meped_boundary.csv.gz`.
- `build_poes_reference.py` requires a `.csv.gz` reference output and records
  gzip compression in the summary.
- `poes_sem2.py` writes deterministic gzip bytes so repeated identical builds
  have stable SHA-256 values.
- `run_C10.py` reads the compressed CSV directly while retaining explicit
  legacy plain-CSV support for developer fixtures.
- The uncompressed checked-in placeholder was replaced by a compressed
  placeholder, and commands/documentation were updated accordingly.


## Modified existing files

- `AMPS_PARAM_C10_gridless.in`
  - replaces synthetic rigidities with P6–P9 nominal-threshold rigidities;
  - uses and documents the `Rc_effective` half-transmission analogue;
  - adds detailed scientific and numerical comments.
- `AMPS_PARAM_C10_mode3d.in`
  - same observable changes as GRIDLESS;
  - documents Mode3D mesh and interpolation sensitivity.
- `README.md`
  - completely rewritten as the scientific, data-acquisition, installation,
    execution, plotting, validation, and reproducibility guide.
- `build_poes_reference.py`
  - replaces the unimplemented NCEI stub with a real archive-driven reference
    builder.
- `reference_C10_poes_meped_boundary.csv.gz`
  - removes the synthetic numerical values and becomes a safe header-only
    placeholder that instructs the user to build the real reference.
- `requirements.txt`
  - adds the numerical, plotting, AACGM, and optional CDF dependencies.
- `run_C10.py`
  - requires archive-derived reference rows for scientific runs;
  - compares POES measurements with AMPS in UTC/rigidity/hemisphere/MLT;
  - adds C9-style `FULL_SCAN` and GRIDDED-only `DIRECT_ACCESS` products;
  - uses common exact-state `ACCESS_T50` as the default observable;
  - retains `Rc_lower`, `Rc_effective`, and `Rc_upper` as FULL_SCAN diagnostics;
  - optionally verifies raw access-state agreement between the products;
  - adds C9-style time/residual, scatter, and MLT comparison plots;
  - adds detailed provenance/coverage fields to comparison output;
  - prevents silent low-rigidity baseline substitution;
  - retains GRIDLESS and GRIDDED branches and C9-style result JSON.

## New files

- `poes_sem2.py` — documented archive readers and boundary extraction library.
- `download_poes_sem2.py` — official NCEI downloader and SHA-256 manifest writer.
- `REFERENCES.md` — ten core dataset/instrument/method references and an explicit
  assessment of why papers cannot replace the archive.
- `data/README.md` — local source-data and TS05-driver layout.
- `tests/test_reference_pipeline.py` — parser, rigidity, crossing, and aggregation tests.
- `tests/test_runner_plots.py` — smoke tests for all three plot products.
- `tests/test_cutoff_modes.py` — access parser, T50, command-selection, and
  FULL_SCAN/DIRECT_ACCESS consistency tests.
- `legacy/` — original synthetic CSV and prototype validation document, clearly
  marked as non-scientific and excluded from defaults.
- `VALIDATION.md` — implementation-validation record and remaining archive test.

## No AMPS C++ source changes

The current Earth cutoff implementation already contains the capabilities added
for C9: Mode3D `RIGIDITY_LIST`, FULL_SCAN companion exact-rigidity access files,
and the common three-state access contract.  C10 reuses those products, so this
update changes only test code and documentation.  The executable must have been
built from a source tree containing the C9 common-T50 implementation.
