# Approved I19/I20 numerical verification and execution

The user approved the v3 design on 2026-11-13. The active, hashed contract is
[M2_I19_I20_PROTOCOL.json](M2_I19_I20_PROTOCOL.json); the preceding analysis
is in [M2_APPROVED_IMPLEMENTATION_CONTRACT.md](M2_APPROVED_IMPLEMENTATION_CONTRACT.md).
No new permanent test ID is assigned. These are extensions of U12/I19/I20 and
the existing global M2 supporting phase. Current outcomes, source revision and
all attempts are recorded in [STATUS.md](../STATUS.md).

## Development execution

From an isolated, fully installed AMPS checkout with no concurrent build:

```sh
python3 -B srcMoon/test/test_runner.py --phase m2-development --mpi-ranks 8 \
  --m2-input srcMoon/examples/m2_development.in --stream \
  --campaign-wall-seconds 3600 --cleanup-reserve-seconds 600 \
  --output-dir test_output/srcMoon/M2-development-NEW
```

The enclosing supervisor shares one monotonic deadline across configuration,
clean native builds, baseline/local checks, data QA, linked cases, plotting
and restoration. Linux process identities and a child subreaper retain orphan
MPI/build descendants even when they create another session. Only descendants
of this campaign are signalled. The work deadline is 3,000 seconds, with a
600-second cleanup reserve within the whole hour. Completion before the limit
must be measured; termination alone is not a successful campaign. Limits are
200,000,000,000 bytes of aggregate process RSS, 50,000,000,000 bytes of retained
artifacts, 8 GiB host available-memory reserve and 100 GiB free disk reserve.
A limit breach is ERROR with resource traces and completed case evidence.

The default run performs configuration-separated M0 requalification, then
builds the dedicated Na/Ar binary with chemistry OFF, individual weights ON,
and the existing AMPS `RND=MERSENNETWISTER` option. Each native build starts by
removing its generated `build/` tree; configuration uses explicit
`-amps-test=off`, the declared SPICE installation, and `make -j16 'OPT3=-O3 -g'`.
Canonical M0/M1 RNG settings, physics and nightly reference bytes are preserved.
The coordinator owns final restoration; the M0 child omits an intermediate
normal build that the immediately following clean M2 build would discard.
At completion, the checkout is restored to the normal Na-only Moon/SPICE build.

Supplying `--m2-amps` skips those clean builds/M0 and restoration. Reuse requires
matching binary, current native source, generated chemistry/layout and RNG
hashes. It is a focused run and cannot prove timing of the complete clean
command. Use a new output directory for every attempt. No prior scientifically
meaningful result or validation-data file is overwritten.

## Physical cases and numerical controls

Both cases use the ordinary production `amps_time_step`, sources, moving LOLA
boundary, force mover, thermal/environment update and surface callbacks. The
development volume mesh is the existing bounded custom regression fixture;
the surface still uses the requested LOLA DEM resolution. Formal cluster
cases preserve the reviewed template's complete volume AMR and physical
configuration instead of applying the I18 coarsening transformation. Short
results do not establish spatial or thermal convergence. Species are neutral Na and Ar, with chemistry
explicitly OFF; the broader chemistry network belongs to M3.

The source-fed case starts with zero gas and zero surface stock and retains
the declared external sources. The cohort case also starts with zero gas,
turns external MIV/Ar sources OFF, and places a declared physical stock on each
selected logical surface element. It keeps finite-stock thermal/PSD release,
real transport, return, adsorption and subsequent release. These are stated
initial-value verification fixtures, not measured lunar abundances, equilibria,
or substitutes for qualified thermal initial conditions.

The fixed short design has two cases, a full 3 by 3 grid and two base seeds:
36 native cases. Common timestep is bounded by the minimum final initialized
Na/Ar timestep, including the existing Ar 0.1 copy. The coarse reference is
frozen at full precision from the first actual receipt. Refinements are
`dt`, `dt/2`, `dt/4`, independently crossed with `W`, `W/2`, `W/4`; reference
weights are frozen before applying either numerical axis. Coarse 24 steps
become 48/96 with exactly the same physical duration. There are eight reporting
intervals; epochs use integer-indexed elapsed time plus the initial SPICE TDB
epoch. No unsupported fractional step or interpolated inventory is invented.

The original positive-coverage I18 fixture retains its frozen base 731993
and its declared rank offset; numerical-template seeds cannot replace it.
Its manifest records the actual input seed. This is a contract restoration,
not selection of a seed by statistical agreement.

Numerical base seeds are `2000003 + 1009*j`, j=0,1; actual rank seed is base plus rank.
Seeding occurs after production initialization. Two cases may share the paired
seed label without having identical subsequent random draws. Diagnostics take
no random draws. Explicit Mersenne Twister selection avoids the small 31-bit
state space of the default generator for long studies, but does not itself
prove independence or justify a confidence interval. Long statistical support
must still be declared and reviewed. OpenMP threads are fixed to one.

The complete optional `.in` control group is shown in
[examples/m2_development.in](../examples/m2_development.in):

```text
#begin application: moon
! Include the full production surface/thermal/species settings from the example.
m2_verification_case = source-fed
m2_verification_profile = development
m2_verification_coarse_steps = 24
m2_verification_dt_divisor = 1
m2_verification_weight_divisor = 1
m2_verification_horizon_s = short
m2_verification_coarse_timestep_s = auto
m2_verification_na_weight_atoms = auto
m2_verification_ar_weight_atoms = auto
m2_verification_cohort_packets = 32
m2_verification_surface_elements = 3000
m2_verification_source_integral_tolerance = 1e-8
#end application
```

This excerpt is not a standalone physical input. `cohort` and `normalization`
are the other active cases. `cluster` requires an explicit positive SI horizon
and aligned coarse step count; `normalization` integrates source rates only.
Explicit positive timestep/weight references in the template are honored in
the first case. Auto references are resolved by native initialization, then replaced with
explicit full-precision numbers in subsequent per-case receipts. Divisors must
be 1, 2 or 4; coarse steps are a positive multiple of eight. Active controls
are accepted only with the I19 budget adapter, not in ordinary science runs.
Omitting the group leaves legacy behavior unchanged. A partial/inactive group,
unsupported configuration, nonfinite value or excessive counter fails explicitly.

Surface elements are logical reservoir indices, not independent DEM triangles
or independent observations. Up to 48 site tags are supported. Several native
sites may map to the same model support; their stocks are additive and tags
are distinct, but their physical histories are correlated. This does not create
new spatial resolution or independent evidence. Site metadata/collocation must
be frozen before a formal run. The single development element 3000 is an
explicit computational fixture without an observational or PSR claim.

## Cohort implementation and independent oracles

`MoonM2Cohorts` registers passive atom-fraction fields in actual particle records
before allocation. The normal MPI particle migration therefore carries tags.
Replicated surface tag stocks are counted once; owner-only deposits/losses are
reduced once at completed steps. Optional Exosphere finite-release hooks observe
the actual stock, competing thermal/PSD rates, exact remaining stock and emitted
packet weight. The actual surface-credit and particle-loss paths update tags
before the particle is destroyed. Null hooks preserve ordinary behavior and
consume no RNG draws. Detailed comments specify collective/owner semantics.

The initial cohort's first-departure history is separate from surface occupancy.
A returned atom is part of current surface stock but never becomes a new
first resident. The diagnostics preserve never-departed survival S(t),
first-release fraction 1-S(t), restricted mean first residence
`integral_0^H S(t) dt` (including censored atoms), and current surface/gas/
escape/numerical-loss fractions. Gross repeated release is a separate ledger.

`m2_protocol.py` independently integrates every recorded constant-coefficient
hazard interval. Zero release gives S=1 and restricted residence H; a constant
rate k gives S=exp(-kH) and `-expm1(-kH)/k`. It independently evaluates the
input-selected Na/Ar thermal laws at the actual recorded temperature, checks
thermal plus PSD hazards sum, and checks tag mass using a forward arithmetic
bound computed from counted operations and MPI reductions. It does not fit an
arithmetic allowance to residuals. These checks are conditional on the actual
temperature history and cannot qualify that thermal history.

U12 additionally imposes linear warming spanning H/tau=0.01 to 1000, using
production residence and finite-depletion kernels. These temperatures are
computational stress-test forcing, not lunar data. A separate SciPy quadrature
integrates temperature-dependent hazards and censored survival, with error
budget below 1% of the numerical margin. Production midpoint steps are 0.5,
0.25 and 0.125 seconds over H=60 s. The finest survival error must be <=0.01,
restricted-residence error <=0.01 H, and refinement differences decrease.
This is a release-kernel check, not a linked thermal or observational PASS.

The original I19 positive-coverage fixture remains separate. Exact gas model-
particle count and independent signed gas/surface/combined atom and kilogram
identities retain their existing criteria. Continuous surface stock has no
model-particle count. External supply excludes internal repeated release;
source-fed imbalance remains below the existing 0.5% supply limit. Numerical
box/mover losses remain distinct from positive-energy lunar escape. Short
development cases retain the original no-numerical-loss gate.

## Artifacts and status

See [RUN.md](../RUN.md#m2-supporting-verification-before-m1-closure) for the
full output tree. Every native case retains input/config/source/build/executable/
SPICE identities, actual seed/layout, stdout/stderr, receipt, result, resources,
cumulative and completed-window ledger CSVs. Completed stages produce
`comparison.csv/.png/.eps`; cohorts also produce hazard/occupancy CSVs and
`cohorts.png/.eps`. The comparison is explicitly against conservation and
analytical kinetic references, never fabricated observations. Skipped formal
tests produce no invented comparison values or blank success plots.

The global runner automatically registers the infrastructure checks and every
36-case short design result in a separate supporting namespace. Supporting
FAIL/ERROR makes the command fail even when formal results are SKIPPED.
PASS means the declared supporting criterion ran and was met. Missing
capability/data/predecessor is SKIPPED. Build/parser/runtime/corruption/resource
or deadline failures are ERROR. Two seeds and a short physical interval do not
establish formal stochastic precision, convergence, M1, M2 or validation.

## Cluster qualification

`test_runner.py --phase m2-cluster` exports a plan only. The fixed complete
design is two cases by nine grid cells by sixteen paired seeds (288 runs),
with candidate H=2,550,240 SI seconds. It is not an exact lunation, D02 averaging
interval or I17 spin-up cycle. Integer alignment may imply millions of steps;
no feasibility is claimed without measured cost. Aggregate planned allocation
is 300 CPU-core hours including coordinator/preflight, retries, monitoring and
plots, not 300 per case. The executable submits no remote batch job.

Execution additionally requires all of these explicit global arguments:

- `--cluster-execute --m2-amps <retained binary> --qualification-root <campaign>`;
- `--sites-file <frozen native mapping.json>`;
- `--forecast-file <measured allocation forecast.json>`;
- `--normalizations-file <frozen physical atom scales.json>`;
- `--statistical-support-file <reviewed interval and RNG evidence.json>`;
- `--cluster-core-hours 300 --cluster-allocated-cores <reserved slots>`;
- `--cluster-memory-limit-bytes`, `--cluster-artifact-limit-bytes`,
  `--cluster-host-free-floor-bytes`, `--cluster-disk-free-floor-bytes`.

Use `--help` for argument syntax; resource limits belong to the actual cluster
allocation and are never guessed from workstation limits. Required proof
fields and exact identity checks are implemented in `run_m2_cluster.preflight`.
Supporting/preflight status cannot satisfy a formal predecessor, and the
native test ID/MPI layout must agree with its manifest. Missing evidence
yields SKIPPED with no native launch; corrupt evidence yields
ERROR. Current missing I13/D02 evidence deliberately blocks formal execution.

The deterministic model-input normalization can be prepared separately:

```sh
python3 -B srcMoon/test/prepare_m2_normalization.py \
  --tier cluster \
  --record-initial-state \
  --amps /absolute/path/to/amps-moon-m2-na-ar \
  --input srcMoon/examples/m2_development.in \
  --sites-file /absolute/path/to/frozen-sites.json \
  --horizon-s 2550240 --mpi-ranks 8 \
  --output-dir test_output/srcMoon/M2-normalization-NEW
```

This integrates actual external Na MIV/Ar production rates with SPICE
heliocentric radius using a converged composite Simpson rule, without advancing
particles or drawing RNG. Explicitly selected production periodic thermal
initialization still occurs before the source-only adapter. Independent
constant Ar Q*H is checked. The
positive per-species scales are fixed before grid residuals: expected external
supply for zero-stock source-fed startup; initial declared stock for cohorts.
It is a reproducible model-input reference, not observational evidence. Default
single-site preparation without a frozen site file is developmental only.

Long numerical scoring uses all predeclared adjacent and finest-grid contrasts
(18), eight common checkpoints, six global metrics/species and three endpoint
cohort metrics/site, at most 48 sites. Bonferroni bounds the family by 8,640
comparisons. Paired Student-t intervals over sixteen seeds are approximate and
require justified marginal/stream assumptions. Mean plus half-width must fit
the 1% scale/fraction/residence margin; half-width must be <=half the margin.
Zero empirical variance alone is insufficient without independently justified
deterministic-zero support. Resolved increasing timestep sensitivity fails.
The long per-run numerical-domain loss ceiling is 0.1% of the fixed atom scale,
separate from roundoff accounting and physical escape. No physical coefficient
or observational threshold is changed by these numerical allocations.

Even complete numerical PASS cannot supply independent PSR/thermal evidence.
I20 remains SKIPPED without its required data/operator/site comparison and
same-campaign predecessors. December hold-out residuals are not inspected.

The resource traces distinguish allocated CPU-core hours from sampled local
user/system CPU. Allocation includes coordinator/preflight and all earlier
phases on the shared clock. The sampled local CPU counter is explicitly a
lower bound because a short-lived process can be reaped between polls; actual
scheduler CPU/charging is unavailable to this scheduler-neutral executor and
must be retained separately for cluster work. Process-tree RSS is a local-host
measurement; it is not a measurement of memory on unobserved remote nodes.

The active JSON retains dated review/planning fields such as implementation
gap descriptions. They are the approved design's historical requirements,
not current execution evidence. STATUS.md and retained campaign receipts
control implementation/verification status. No proposal field or READY
model-input normalization constitutes formal or observational PASS.

The statistical-support document must also contain
`actual_epoch_thermal_qualification`: READY status, the exact input-template
SHA-256, `simulation_start_utc`, and a nonempty `evidence_files` list of
`path`/`sha256` records. Those files are verified before long execution.
A general I17 PASS at another epoch cannot supply this initial-state proof;
missing proof remains SKIPPED. The current development initialization is a
declared transient fixture, not a qualified lunar thermal initial state.

The formal staging correction is specified in
[FORMAL_STAGING_CONTRACT.md](FORMAL_STAGING_CONTRACT.md). Cluster input now
requires Ar already enabled, an evolving environment and zero initial surface
inventory. It never enables a missing physical prerequisite or replaces AMR,
source strengths, thermal controls, initialization mode or epoch. Missing
required configuration is SKIPPED. Development retains the bounded I18
transformation. Formal case staging may change only the approved seed and
numerical controls, output ownership and checkpoint/output suppression.

`m2_initial_state.physics_identity_sha256(path)` hashes canonical physical
settings from the actual staged `.in`. Known input paths resolve from their
assigning file, including assignments in `#include` fragments; case and literal
backslashes are preserved. Comments and continuations follow the native
grammar, and shell/environment expansion is forbidden. Unknown fields remain
identity-bearing, including unknown `m2_verification_*` fields; exclusions use
an exact named list. Completion-only approval and diagnostic-query paths
(`m2_formal_contract_file`, `m2_illumination_queries_file`,
`m2_polar_queries_file`) are independently bound by their formal contracts and
may vary without changing initial physical coefficients. The initial thermal mode and every spin-up/material/time
control are retained. The native omitted-mode default is represented as
`declared-uniform`, without treating that default as qualified. File-content
hashes remain the responsibility of the existing data/state evidence gates.

Normalization preparation also has an explicit staging tier. Its default
`development` mode preserves the existing bounded regression reference;
`--tier cluster` preserves the reviewed physical template and is required for
formal use. Both use the declared positive physical horizon for deterministic
source quadrature, without injecting particles. READY normalization records
the actual staged physics hash and input/manifest/receipt descriptors, plus
the native coarse initial stock and the declared cohort-packet count from
which its ordinary-initialization reference weights are derived. This binds
cohort normalization to the same AMR-dependent initial weights as production.
The original template hash alone cannot qualify a bounded reference for a
different formal volume mesh. Older references without this evidence remain
SKIPPED until regenerated; no reference is automatically promoted or frozen.

`--record-initial-state` (Python API `record_initial_state=True`) also retains
the actual I19 initial full-column export from this cheap source-only job.
Cluster staging already requests that diagnostic; development can opt in
without a transport run. Requested missing/incorrect CSV, receipt counts,
method or actual SPICE/DEM identity is ERROR. The result records an
`candidate_initial_state` with `qualified: false`; model-input READY never
means the state is qualified, reviewed, observational truth or an I19 PASS.
The global development command requests this export in its normal source
normalization stage so the real I19 dispatcher remains covered without an
extra transport case.

Every formal case requests
`m2_verification_record_initial_thermal_state = true`, records its staged
physics hash and a UTC prelaunch `started_utc` in `run_manifest.json`, and
retains the actual native full-column state. A reviewed spin-up document alone
cannot prove that a process used that state. The new initialization-evidence
validator separately checks the selected native initialization algorithm,
actual start epoch and state bytes. No spin-up interval, physical initial
temperature or qualification criterion is selected by this staging code.

The compact cohort PNG/EPS previews show cohort 0. Complete per-site/species
statistics and hazard histories remain in the CSVs; all tags are evaluated.

The corrected clean development command was measured at **1849.678 s** on
2026-11-13, with six local PASS, seven formal SKIPPED and supporting 44 PASS /
two SKIPPED, zero FAIL/ERROR. All 36 numerical cases pass; this demonstrates
the development execution tier, not formal M2 qualification. Source snapshots,
failed attempts and exact commands are in STATUS.md section 5.13.
