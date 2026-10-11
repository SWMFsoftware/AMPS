# SA01: standing instructions for Codex

## Scope and installation

Place this file and its companions in:

`/data/vtenishe/GEO/AMPS/srcEarth/test/SA01/`

These instructions govern the SA01 runners, comparison tools, figures, documentation,
and related work requested for this test. Read applicable ancestor `AGENTS.md` files
before modifying the repository. Preserve unrelated user changes. Do not overwrite
an existing instructions file without reading and merging its relevant requirements.

## Required reading and continuity

At the start of this task and every resumed session:

1. Read this file.
2. Read `SA01_TASK.md` for the complete execution specification.
3. Read `SA01_PROGRESS.md` for completed work, evidence, blockers, and the next action.
4. Verify the recorded repository revision, outputs, and process state before resuming.

Maintain `SA01_PROGRESS.md` after each meaningful implementation stage, build,
verification, physical run, comparison, or scientific decision, and before stopping.
Record failed attempts as well as successes. A progress file is an explicit work
record; its existence does not establish that a run or comparison occurred.

Keep completed history append-only. Correct a mistaken entry with a dated correction,
not by erasing the original. Update the current-state section and next action in place.
Do not restart completed expensive work unless artifacts are missing, changed, or
invalid for the current revision. Explain the reason for any rerun.

## Authoritative locations

| Item | Location |
| --- | --- |
| Live AMPS repository | `/data/vtenishe/GEO/AMPS` |
| Application | `srcEarth` under the live repository |
| New test implementation | `srcEarth/test/SA01` |
| Extracted data package | `/data/vtenishe/GEO/sep-in-geospace-sa01-sa04-standalone-20260926` |
| Processed SA01 data | `validation_data/processed/SA01` under that package |
| Raw source manifest | `validation_data/raw/SA01/source_manifest.json` under that package |

The data path is an extracted directory, not a tar file. Resolve it through the
runner CLI. Do not hard-code it in scientific routines. Render absolute driver,
boundary-spectrum, and trajectory paths; packaged `../../../` paths depend on the
launch directory. Do not copy the older packaged source over the live repository.

## Mandatory rebuild procedure

Whenever recompilation is needed, remove only the verified live repository's
`build` directory, then invoke `make -j` from the repository root:

```bash
cd /data/vtenishe/GEO/AMPS
rm -rf build
make -j
```

In Python, validate the resolved repository path before deletion. Refuse a `build`
symlink or an unexpected deletion target. Capture the build log and return code.
Do not rebuild while another build or physical run uses this repository's products.
Check observed process state first. Do not terminate another session's work.
Use the user's existing build configuration; do not invent `Makefile.conf` contents.
Use this procedure when building, not an alternate build command.

## Scientific contract

- Event: `2014-01-06T09:00:00Z <= t < 2014-01-10T00:00:00Z`.
- Eligible spacecraft: RBSP-A and RBSP-B; inspect the packaged quality/provenance rules.
- Species: protons; standalone T05 field with Qin-Denton drivers.
- Upstream spectrum: SOHO/ERNE, independent of the target REPT observations.
- Comparison domain: seven REPT differential channels, 18.5–116 MeV; L >= 3.2.
- Baseline boundary: isotropic and spatially uniform; no electric field.
- Products: differential spectrum, seven channel integrals, total integral intensity,
  and energetic-proton density over the same finite energy interval.
- REPT density is isotropic-equivalent and derived from spectra; do not label it a
  direct density measurement or count it as an independent observation dataset.
- Successive static field snapshots do not establish dynamic electric-field transport.
- Keep all supplied observations immutable. Do not redownload, reprocess, or change
  event selections merely to improve agreement. Record justified processing changes
  as separate versions; preserve the original reference and unchanged baseline.

## Fixed acceptance requirements

For spectrum, channel flux, total flux, and density separately:

| Metric | Requirement |
| --- | --- |
| Median absolute log10(model/reference) | <= 0.30 |
| Fraction within factor two, each spacecraft | >= 0.60 |
| Fraction within factor two, pooled | >= 0.70 |
| Maximum unresolved fraction | <= 0.01 |
| Matched epochs per eligible spacecraft | >= 12; also report counts by product |
| Integrated-product numerical convergence | <= 0.02 fractional change on a common comparison population |

Preserve stricter applicable repository gates. Never relax thresholds, omit failed
channels, fit spacecraft-specific scales, or turn missing comparisons into PASS.
Physical zero predictions against positive references are failed predictions, not
rows to remove. Numerical noncompletion is unresolved, not magnetic exclusion.

## Implementation and evidence

- Keep new runners and SA01-specific support code in this directory.
- Reuse correct existing study modules after reviewing their scientific behavior.
- Add detailed comments/docstrings to all new code and update relevant READMEs.
- Provide useful CLI help and examples; keep generated results in a selected output directory.
- Isolate RBSPA/RBSPB outputs and run the branches sequentially by default.
- Record code/executable/input/data fingerprints, commands, dependencies, time stamps,
  resource layout, selection counts, return codes, metrics, and artifact paths.
- Provide standalone comparison/plot modes that do not rerun AMPS.
- Distinguish synthetic pipeline checks, data-integrity checks, physical execution,
  observational agreement, and numerical convergence.
- Do not invent completion, runtime, uncertainty, or scientific performance numbers.
- Do not commit, push, or publish repository changes unless the user requests it.

## Completion and handoff

Follow the completion checklist in `SA01_TASK.md`. If a required step fails or is
blocked, document the exact evidence and next action rather than declaring success.
Before ending, update `SA01_PROGRESS.md`, list actual changes and executed checks,
and identify the report/figure locations. Do not leave a detached physical run without
recording how to identify it and whether it is still running.
