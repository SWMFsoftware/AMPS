# VP12: SEP source records, identity, attribution, and source spectrum

## Purpose

This case implements validation-plan item **VP12** in the extensible case-owned campaign layout. This is a controlled coupling benchmark. It does not claim an external AMPS science result because AMPS is not distributed with SWCME.

## What is compared

The runner obtains candidate values through the public SWCME production API (or, for transport-only quantities, through the explicitly labelled controlled coupling benchmark) and compares them row by row with `reference_solution.py`. The raw model/probe output, comparison CSV, result JSON, artifact manifest, and both PNG and EPS figures are preserved.

## Why it is tested

SEP source records, identity, attribution, and source spectrum is a distinct validation contract: a passing lower-level unit test does not by itself prove that the case selection, units, independent calculation, acceptance threshold, plotting, and campaign reporting remain mutually consistent.

## Data and provenance

Run `python3 download_data.py` to materialize `data/raw/benchmark.json`. It is intentionally offline and immutable; `data/PROVENANCE.json` identifies the evidence class and scope. This is a controlled coupling benchmark. It does not claim an external AMPS science result because AMPS is not distributed with SWCME.

## How to run

```sh
python3 download_data.py
python3 test_reference_solution.py
python3 run_vp12.py
# or from swcme/test:
make validation-case CASE=VP12
```

The last command can add `VALIDATION_ARGS=--download` on a fresh checkout. Outputs are written beneath this case's `output/` for direct execution and beneath the timestamped global campaign directory when orchestrated.

## Expected result

The case passes only when every criterion in `vp12_result.json` is true. A PASS is scoped to the evidence class recorded there; limitations are machine-readable and must be carried into any science-release claim.

