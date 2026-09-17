# VP10: Parker mapping, field-line path length, and focusing

## Purpose

This case implements validation-plan item **VP10** in the extensible case-owned campaign layout. The numerical target is analytical. It validates production physics and the comparison path; any observational-release limitation is recorded in the result JSON.

## What is compared

The runner obtains candidate values through the public SWCME production API (or, for transport-only quantities, through the explicitly labelled controlled coupling benchmark) and compares them row by row with `reference_solution.py`. The raw model/probe output, comparison CSV, result JSON, artifact manifest, and both PNG and EPS figures are preserved.

## Why it is tested

Parker mapping, field-line path length, and focusing is a distinct validation contract: a passing lower-level unit test does not by itself prove that the case selection, units, independent calculation, acceptance threshold, plotting, and campaign reporting remain mutually consistent.

## Data and provenance

Run `python3 download_data.py` to materialize `data/raw/benchmark.json`. It is intentionally offline and immutable; `data/PROVENANCE.json` identifies the evidence class and scope. The numerical target is analytical. It validates production physics and the comparison path; any observational-release limitation is recorded in the result JSON.

## How to run

```sh
python3 download_data.py
python3 test_reference_solution.py
python3 run_vp10.py
# or from swcme/test:
make validation-case CASE=VP10
```

The last command can add `VALIDATION_ARGS=--download` on a fresh checkout. Outputs are written beneath this case's `output/` for direct execution and beneath the timestamped global campaign directory when orchestrated.

## Expected result

The case passes only when every criterion in `vp10_result.json` is true. A PASS is scoped to the evidence class recorded there; limitations are machine-readable and must be carried into any science-release claim.

