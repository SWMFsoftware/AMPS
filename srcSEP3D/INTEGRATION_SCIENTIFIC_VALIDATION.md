# Phase V — Integration and Scientific Validation

Phase V converts the numerical model assembled in Phases M/B/T/P/A/O into an
auditable scientific application. It does not equate component verification
with coupled or observational validation. Four evidence classes remain
separate in every JSON/JUnit report:

| Evidence class | What it establishes | What it cannot establish |
|---|---|---|
| controlled (`INT3D`, `VFY3D`) | deterministic reductions, metrics, analytical limits, and source statistics | correctness of an AMPS/MPI build or a heliospheric event |
| linked (`NAT3D`, `MPI3D`) | real mesh/storage/list/MPI/restart behavior | agreement with another model or observations |
| cross-model (`XM3D`) | agreement with an independently exported calculation | observational fidelity |
| observational (`OV3D`) | agreement with reviewed spacecraft evidence under declared metrics | physics outside the configured model scope |

A missing linked executable or evidence bundle is `SKIP`, not `PASS`. A
present but malformed/checksum-invalid bundle is `ERROR`. Valid evidence that
misses a tolerance is `FAIL`.

## 1. Deterministic integration audit

`validation/validation_metrics.{h,cpp}` receives one `RankPartition` per host
rank. Each partition carries immutable particle observations, closed integer
ledger rows, elapsed time, and peak resident memory. The audit is deliberately
MPI-independent: the host may collect records with MPI, but the acceptance
algorithm never uses an order-dependent floating collective.

The merge is:

1. Validate unique rank labels and finite nonnegative resource measurements.
2. Concatenate observations and sort by persistent `stableId`.
3. Reject stable ID zero or a duplicate anywhere in the global set. This
   catches duplicated migration/restart records even when each rank is locally
   valid.
4. Sum closed ledger rows by `(step,species)` with overflow checks.
5. Require, in exact integer arithmetic,

   \[
   N_\mathrm{start}+N_\mathrm{inject}=
   N_\mathrm{end}+N_\mathrm{escape}+N_\mathrm{absorb}+N_\mathrm{fail},
   \qquad N_\mathrm{advanced}=N_\mathrm{end}.
   \]

6. Report load imbalance as `max(N_rank)/mean(N_rank)`, maximum rank wall
   time, and summed peak resident bytes.

Structural validity and budget acceptance are separate. A resource miss leaves
the measured record usable and sets `withinBudget=false`; corrupted identity
or conservation state returns a typed error and no budget claim.

`INT3D01–03` exercise repartition/order invariance, global conservation and
duplicate detection, and load/wall/memory thresholds. The external linked
cases `NAT3D01–03/09–12` and `MPI3D01–02` must still run through a configured
AMPS executable.

## 2. Scientific curve comparison

`ComparePositiveSeries()` accepts strictly increasing coordinates and positive
finite values. The manifest supplies the physical names and units; the C++
routine is dimension-agnostic and can therefore compare time profiles, radial
profiles, or energy spectra without guessing a column meaning.

The model is interpolated in `ln(value)` at covered reference coordinates. No
extrapolation or artificial intensity floor is permitted. This choice is
appropriate for SEP intensities spanning orders of magnitude and prevents a
short model interval from silently omitting onset or decay data.

Three normalization policies are explicit:

- `absolute`: no fitted factor;
- `one-global-amplitude`: one factor
  \(A=\exp\langle\ln(J_\mathrm{ref}/J_\mathrm{model})\rangle\) for an
  explicitly uncalibrated source area;
- `unit-peak`: shape-only diagnostic; it is never reported as absolute-flux
  validation.

The common metrics are coverage, log10 RMSE, median absolute log10 error,
log-space correlation, onset-coordinate error at one percent of peak,
peak-coordinate error, peak ratio, and trapezoidal fluence ratio. One fitted
amplitude cannot alter timing, profile shape, anisotropy, or spectral slope.
Each threshold is stored in `validation/case_registry.json` and copied into the
case report.

## 3. Controlled physics validations

The AMPS-free `phase-v` suite validates the machinery using production kernels:

| ID | Algorithm and physical reference |
|---|---|
| `VFY3D01` | recovers an identical profile under exactly one declared global amplitude and verifies every shape metric |
| `VFY3D02` | rejects incomplete coverage scientifically and malformed coordinates structurally |
| `VFY3D03` | projects the 3-D Parker SDE along three field orientations and compares empirical CDFs with the 1-D Gaussian Green function |
| `VFY3D04` | demonstrates second-order convergence of focused transport to the exact focusing characteristic `mu=tanh(atanh(mu0)+at)` |
| `VFY3D05` | compares sampled SWCME/DSA momenta with the independent truncated-power-law CDF and closes total event weight |

These are scientific verification prerequisites. They do not satisfy an
`XM3D` or `OV3D` release gate by themselves.

## 4. External evidence contract

`validation/run_validation.py` never locates the sibling `srcSEP` tree. A
cross-model calculation is run independently and exports its result into a
Phase-V bundle. This preserves application independence and records exactly
which bytes were compared.

For a series case, place below the supplied evidence root:

```text
EVIDENCE_ROOT/
└── XM3D01/
    ├── manifest.json
    ├── model.csv
    └── reference.csv
```

Both CSVs have the exact grammar:

```csv
coordinate_si,value_si
0.0,1.0
1.0,2.0
```

`manifest.json` follows
`validation/templates/scientific_evidence.template.json`. It records case ID,
quantity names and units, both SHA-256 values, and nonempty model/reference
provenance. Paths must remain inside the case directory. `XM3D04` additionally
requires at least one `exact_pairs` entry whose independently hashed source
artifacts are byte-identical.

`XM3D05` uses the convergence template. At least three positive
`(resolution,error)` pairs are ordered coarse-to-fine and evaluated with

\[
p_i=\frac{\ln(e_i/e_{i+1})}{\ln(h_i/h_{i+1})}.
\]

The minimum observed order and the declared default-to-finest difference are
both release criteria.

## 5. Registered scientific campaign

| Cases | Role |
|---|---|
| `XM3D01` | Parker 3-D versus independently exported field-line profiles; release gate |
| `XM3D02` | focused 3-D versus focused field-line onset/peak/profile; release gate |
| `XM3D03` | disconnected-longitude diagnostic; remains SKIP until perpendicular diffusion is released |
| `XM3D04` | exact common SWCME source plus transported spectrum; release gate |
| `XM3D05` | mesh/timestep convergence; release gate |
| `XM3D06` | independent 3-D finite-volume reference; extended |
| `OV3D01` | 2013-04-11 ACE/GOES/SOHO comparison; release gate |
| `OV3D02` | 2020-05-29 PSP/STEREO-A comparison; release gate |
| `OV3D03` | 2014-01-06 PAMELA connection-sensitive diagnostic |
| `OV3D04` | near-relativistic multi-spacecraft electron diagnostic |

The observational manifests must name the publication/dataset, instrument,
processing, uncertainty treatment, coordinate epoch/frame, and model
configuration fingerprint. The generic template contains the minimum fields;
event-specific provenance may add more fields without changing the parser.

## 6. Physical scope of interpretation

Release 1 has three-dimensional geometry and local background evaluation, but
parallel-only particle transport. It should reproduce well-connected events
within the limits of the background, source, and scattering model. It is
expected to under-predict poorly connected observers because
`kappa_perpendicular=0` and gradient/curvature drifts are disabled. Therefore
`OV3D03`, `OV3D04`, and `XM3D03` quantify a known physics gap and are diagnostic,
not evidence that missing perpendicular physics has been validated.

## 7. Commands

From `AMPS/srcSEP3D`:

```bash
# Controlled Phase-V prerequisites plus explicit external SKIPs.
test/run_tests.py --suite phase-v --rebuild \
  --output-dir test_output/phase-v

# Cross-model or observational evidence already exported by its owner.
test/run_tests.py --suite phase-v --validation-data /path/to/evidence \
  --output-dir test_output/phase-v-evidence

# Linked native cases. The prefix is argv, not a shell command.
test/run_tests.py --suite phase-v --amps ../amps \
  --validation-launch-prefix "mpiexec -n 8" \
  --output-dir test_output/phase-v-linked

# Run only the external campaign interface.
validation/run_validation.py --all --amps ../amps \
  --evidence-root /path/to/evidence \
  --output-dir test_output/phase-v-campaign
```

The linked executable must advertise a requested case through `--list-tests`
and emit the standard test JSON. An executable from an older source tree is an
`ERROR`, not a skipped or passing integration result.
