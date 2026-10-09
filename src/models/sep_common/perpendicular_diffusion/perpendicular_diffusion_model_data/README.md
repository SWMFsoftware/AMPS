# Perpendicular diffusion specification 2.1: verification assets

This archive accompanies `PERPENDICULAR_DIFFUSION_COEFFICIENT_MODEL.md`.
It provides reproducible equation checks and source-transcription evidence.
It contains no production AMPS implementation or particle-simulation results.
The final section of the specification gives the implementation roadmap.

Unzip with its directory structure intact. The specification is at the archive
root and these assets are in `perpendicular_diffusion_model_data/`.

## Reproduce the checks

Use Python with NumPy, SciPy and mpmath. The successful check report records the exact
versions used; these are execution provenance, not mandatory physics inputs.

From the data directory run:

```bash
python3 reference_verification.py
python3 review_additions_verification.py
python3 document_fixtures.py
python3 document_audit.py
```

The first command verifies the checksums, recomputes the equation fixtures,
compares them with the stored JSON, and checks the paired NLGCE backends.
The added review script independently computes the short-FLPD roots and branch
parameters, nonanalytic-tail behavior, ordinary/backtracking time transforms,
source integral lengths, closed-form and paired bounds, and estimator/angular
identities. It also checks the stated rounding of the printed spectral constants
and GCD coefficients. mpmath uses 50 decimal digits for these calculations; this
is a numerical control, not physical precision.
The table script compares all 12 printed numeric tables with the fixture JSON.
An alternative document location can be provided using `--document`.
These normal commands read the delivered assets without modifying them.

For an independent source-table check, download the primary Qin–Zhang PDF
identified in `NLGCE_source_audit.json`, extract it with
`pdftotext -layout source.pdf source.txt`, and run:

```bash
python3 source_table_verification.py --primary-text source.txt --primary-pdf source.pdf
```

`--generate-fixtures`, `--write-results`, and the table script's `--write` are
maintenance options for a deliberate specification revision. Changing assets
requires regenerating the checksum manifest; use `--skip-checksums` only during
that maintenance step, then verify the final archive without it.

## Files

| File | Purpose |
|---|---|
| `reference_verification.py` | Mathematical and numerical checks for the perpendicular equations |
| `review_additions_verification.py` | Independent revision 2.1 algebraic/numerical checks and printed-constant audit |
| `review_additions_fixtures.json` | Computed review-addition values with explicit synthetic-input notes |
| `review_additions_results.json` | Actual successful results, dependency versions and engineering tolerances |
| `benchmark_points.json` | Dimensionless fixtures, input conventions, and numerical controls |
| `verification_results.json` | Results of the actual successful equation/fixture verification |
| `document_fixtures.py` | Exact comparison of generated, rounded Markdown tables with fixture data |
| `paired_reference_verification.py` | Direct NLGCE-N and polynomial NLGCE-F pair checks; independent Decimal polynomial evaluation |
| `paired_benchmark_points.json` | Seven corrected paired states reused from parallel specification 1.3 |
| `NLGCE_F_2014_parallel.csv` | All 288 parallel polynomial coefficients with explicit indices |
| `NLGCE_F_2014_perpendicular.csv` | All 288 perpendicular polynomial coefficients with explicit indices |
| `NLGCE_source_audit.json` | All 576 coefficients compared by index with primary Tables 3/4 |
| `source_table_verification.py` | Repeat the indexed primary-table comparison on locally supplied source text |
| `source_audit.json` | Source versions, equation roles, inspection scope, and unresolved source gates |
| `document_audit.py` | Repeat the document structure and table checks |
| `document_audit.json` | Structural checks, defined equation/reference identifiers, and table checks |
| `SHA256SUMS` | SHA-256 digest for each delivered data/code/audit asset, excluding the manifest itself |

The Markdown document is also copied into this data directory so the checksum
manifest covers its exact text. Its top-level archive copy is identical.

## Evidence and limitations

All unit-scale states, sample Parker inputs, nonunit restoration checks, and
artificial finite-time fit parameters are **defined test inputs**. They are not
observations or recommended physical defaults. The fixtures report numerical
precision; they do not claim a closure has that physical accuracy.

The checks include reduced/independent integrals, spectrum normalization and
moments, scalar residuals, refined quadrature, Gaussian moments, running-model
chain rules, tensor covariance/eigenvalues, and Cartesian tensor divergence.
The paired checks use two polynomial arithmetic methods and compare nonlinear
solutions under quadrature refinement and a changed starting scale. A mismatch
between NLGCE-F and NLGCE-N is retained as a surrogate discrepancy, not hidden
with a multiplier.

The GCD source's printed exact prefactor and its nearby approximate value are
inconsistent. The specification uses the inspected conference formula and an
independent Gaussian integral, and explicitly leaves a journal-version
comparison unresolved. The Kuhlen perpendicular calibration parameters absent
from the inspected table remain required inputs. Unresolved historical presets
remain unavailable under the source-gate policy. In particular, no Corti angular
width is inferred; missing Hall calibration-domain metadata and non-Kolmogorov
Snodin ordered-field definitions remain gated. Frozen-line sharing and coordinate
choices must be explicit in the future transport adapter. The model roadmap
states these implementation requirements without assuming missing physics.

Source PDF hashes identify the primary documents inspected. The source PDFs
themselves are not included. Checksums demonstrate byte integrity, not physical
correctness. The future implementation must perform the coefficient, derivative,
transport-operator, boundary, memory, and stochastic-convergence acceptance
tests specified in the roadmap before a scientific release.
