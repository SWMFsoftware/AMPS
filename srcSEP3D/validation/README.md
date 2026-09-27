# Phase-V validation runner

## V03–V05 qualification

Run the following command for a real MPI matrix:

```bash
python3 validation/run_native_matrix.py \
  --amps /path/to/native-test-amps \
  --test-input /path/to/reviewed-sep3d.in \
  --profile small \
  --launcher "mpiexec -n {ranks}" \
  --output-dir test_output/native-small
```

Profiles are declared in `native_profiles.json`; source-only execution is
never MPI evidence. The literal `{ranks}` field is mandatory so the runner can
prove that each profile row requests its declared process count. Before the
matrix starts, the runner invokes `--list-tests` with one rank and requires the
linked host to advertise every profile case. An ordinary production `amps`
driver, the AMPS-independent `test/stage1`, or a stale linked host therefore
fails before any evidence is accepted.

Each matrix invocation receives the same explicit `--test-input`. The summary
hashes that deck, the executable, the profile registry, every native JSON
report, and every stdout log. Existing per-case `native.json` files are removed
before launch, so an interrupted prior campaign cannot supply a stale PASS.
`v04_campaign.json` defines the ordered scientific ladder. Its live-SWMF rung
and `SWMF3D01` remain blocked by deferred R8. Release profiles and evidence
generation are documented in `../release/README.md`.

`run_validation.py` is the external-evidence half of the srcSEP3D test system.
The public entry point remains `test/run_tests.py`; this lower-level runner is
also available for batch campaigns that need only linked, cross-model, or
observational cases.

## Evidence ownership

- srcSEP3D owns the registry, thresholds, parsers, metrics, and reports.
- A linked AMPS executable owns native `NAT3D`/`MPI3D` callbacks.
- The producer of a cross-model reference owns its exported bytes and
  provenance. The runner never searches for another application checkout.
- The observational-data preparer owns digitization, uncertainty, instrument
  response, and publication attribution.

This boundary prevents a passing controlled calculation from being presented
as observational validation and prevents a mutable sibling source tree from
silently changing a comparison.

`cases/2013-04-11/` contains the reviewed OV3D01 setup blueprint. It records
publication-backed CME/connectivity/source/transport values separately from
the SWCME surrogate and numerical parameters that still require event data or
convergence. The registry validates this blueprint on every invocation, but
the blueprint is not observational evidence and is intentionally not an
executable input deck. See
[`cases/2013-04-11/README.md`](cases/2013-04-11/README.md) before preparing an
OV3D01 campaign.

That case includes `event_fit.template.json` and `prepare_case.py`. The
preparer refuses incomplete event fits, derives coordinate-dependent input
values from a reviewed rotation, optionally runs the rendered deck through the
linked executable's `--dry-run`, and creates checksum-owned scientific bundles
from already calibrated two-column products. It deliberately does not perform
instrument calibration or fit unresolved SWCME parameters.

## CLI

```bash
validation/run_validation.py --list
validation/run_validation.py --case XM3D01 \
  --evidence-root /path/to/evidence --output-dir test_output/XM3D01
validation/run_validation.py --case NAT3D01 --amps ../native-test-amps \
  --test-input /path/to/reviewed-sep3d.in \
  --launch-prefix "mpiexec -n 8" --output-dir test_output/NAT3D01
validation/run_validation.py --all --amps ../native-test-amps \
  --test-input /path/to/reviewed-sep3d.in \
  --evidence-root /path/to/evidence --output-dir test_output/phase-v
```

Selectors are mutually exclusive. Repeated `--case` values are
case-insensitive, de-duplicated, and executed in stable ID order. The optional
launch prefix is parsed into direct process arguments; it is not evaluated by
a shell. Supplying `--amps` for a linked case also requires `--test-input`;
neither the validation runner nor a native callback may choose hidden mesh,
transport, source, or population-control defaults.

## Outputs

Each case writes `CASE_ID/result.json`. The run also writes
`validation-summary.json` and `validation-summary.xml` transactionally. Reports
preserve evidence class and role alongside status, message, metrics, artifacts,
and elapsed time.

Exit codes match the common test convention:

| Code | Meaning |
|---:|---|
| 0 | every selected case passed or explicitly skipped |
| 1 | at least one valid case missed an acceptance threshold |
| 2 | evidence, executable, schema, checksum, or runner error |

## Bundle preparation

Copy the appropriate file from `templates/` to
`EVIDENCE_ROOT/CASE_ID/manifest.json`, add the model/reference products, compute
SHA-256 for every declared file, and replace all instructional text with
specific provenance. The runner rejects path traversal, missing/empty
provenance, unordered coordinates, nonpositive values, checksum changes, and
insufficient overlap before interpreting scientific metrics.

See [../INTEGRATION_SCIENTIFIC_VALIDATION.md](../INTEGRATION_SCIENTIFIC_VALIDATION.md)
for equations, algorithms, case roles, physical limitations, and release-gate
interpretation.
