# M3 approved qualification infrastructure

[CONTRACT.md](CONTRACT.md) records the implementation inventory and test
contract. [The approved addendum](../M3_QUALIFICATION_ADDENDUM_20261118.json)
preserves the original frozen protocol/criteria and records the user's
2026-11-18 choices. It changes no previous result or observational membership.

The production `species/QualificationRandom` is a stateless HMAC-SHA256
counter stream, using the maintained shared DataIO hash and the construction
defined in [RFC2104](https://www.rfc-editor.org/rfc/rfc2104.html). Its public
seed is a reproducibility key, not an authentication secret. Explicit cohort,
original-parent, genealogy, purpose and draw-index fields give the same history
under rank/reordering changes; physical daughters get distinct genealogy.
An exact52-bit midpoint gives strict open uniforms. This is a finite-grid
construction/statistical model; finite optical-depth tails and maximum draw
family remain explicit numerical/statistical evidence requirements. No mutable
generator, clock advance, PIC record/address or allocation-order seed is used.

`m3_qualification.py` implements independent endpoint-survival rate/lifetime
planning and scoring, exact conditional branch intervals, rare-branch precision
planning and numerical-bias checks. I22 uses total-removal `-log(S/N)/T`; branch
frequencies/partial rates and the event-time `D/exposure` estimator remain
separately named diagnostics. It never changes estimator after outputs.
Every complete source branch is retained, including source-zero branches.
Original parents have equal initial represented number within each binomial
cohort; daughters/splits/repeated samples are not extra trials. Unsupported
censoring or unequal initial weights are refused. Exact Clopper–Pearson/Binomial
planning requires SciPy; absence is SKIPPED, never a Gaussian substitute.

The readiness boundary is explicit:

```sh
/usr/bin/python3 -B srcMoon/test/m3_qualification_executor.py \
  --case-plan /path/to/frozen-cohort-plan.json \
  --campaign-root /path/to/same-campaign-predecessors \
  --output-dir test_output/srcMoon/m3-qualification-NEW
```

The plan must bind its approved protocol, complete source-channel list,
actual represented horizon, original cohort/seed and source-state profile.
The runner verifies actual artifacts/hashes and predecessor statuses, writes
result/manifest/input-plan/stdout/stderr and reports SKIPPED while complete
physical cohort providers/executors/numerical/ledger proofs are unavailable.
It invents no native CLI/receipt/figure and accepts no Boolean shortcut. Even
all-green infrastructure predecessor fixtures cannot create a native provider.
I21–I24 and full U24 remain unqualified when their required scope is missing.

Local supporting controls run independently of an AMPS rebuild:

```sh
/usr/bin/python3 -B srcMoon/test/test_m3_qualification.py \
  --output-dir test_output/srcMoon/m3-qualification-controls-NEW
```

They compile the actual RNG/DataIO, compare independent Python HMAC vectors,
check strict-open endpoints, semantic key/shard identities and fixed statistical
controls, then verify the maintained independent scorer/planner on explicitly
analytical fixtures. These fixtures are not observations or source profiles.
CSV/PNG/EPS comparisons, compiler/source/executable identities and all outcomes
remain in the owned new output directory. Native application campaigns require
coordinated source freeze and one3600s/200GB package (larger cohorts need an
explicit cluster tier); these supporting controls do not configure AMPS.

The source-installed I21 adapter now connects this scorer to the actual normal
cohort owner. It still requires complete authoritative state/Q inputs, a fresh
source-matched compiled producer and frozen normal `.in` include identities.
`m3_qualification_executor.native_request` validates the concrete adapter schema;
`sodium_cohort/audit.py` reconciles every removal with per-rank committed event
records before independent ledger and first-hit audits. Three increasing common
step levels at both1/8MPI are mandatory. The whole actual attempt uses one
original resource clock and retains attempted-launch journals on failure.
The partial current Na profile remains SKIPPED; I22–I24 have no admitted complete
provider. Full matrix/milestone qualification is distinct from one cohort result.
