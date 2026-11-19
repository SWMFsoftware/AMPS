# Independent P02 mechanical audit

`verify_p02.py` independently checks the approved, explicitly pinned P02 budget
against authenticated source masses and actual committed native Na event
velocities. It is supporting verification of the controlled sodium model,
not full I21/M3 qualification or observational validation. The free electron
is an explicitly modelled reservoir product here, not a PIC electron particle.
The preimplementation contract is [P02_CONTRACT.md](P02_CONTRACT.md).

## Run after actual native execution

Run the maintained independent verifier from an AMPS checkout, supplying the
original completed source-matched native owner's pinned artifacts. The verifier
does not need a generated PIC build and does not reconfigure its checkout:

```sh
cd /absolute/path/to/AMPS
/usr/bin/python3 -B srcMoon/test/sodium_cohort/verify_p02.py \
  --binding /absolute/path/to/frozen-p02-owner-binding.json \
  --output-dir /absolute/path/to/NEW-p02-independent-audit
```

The new binding contains source/approval identities and original completed
owner identities, rather than arbitrary native scalar/reference pairs:

```json
{
  "schema": "srcMoon-M3-P02-owner-binding-v1",
  "profile": {"path": "/absolute/profile.json", "sha256": "PINNED_SHA256"},
  "budget": {
    "path": "/absolute/M3_P02_NUMERICAL_BUDGET.json",
    "sha256": "PINNED_SHA256"
  },
  "owners": [{
    "owner_manifest": {"path": "/absolute/native-owner/run_manifest.json", "sha256": "PINNED_SHA256"},
    "owner_result": {"path": "/absolute/native-owner/result.json", "sha256": "PINNED_SHA256"},
    "mpi_ranks": 8
  }]
}
```

Replace the marked values with the exact original producer/profile/budget
paths and digests. Relative paths belong to this binding file. The caller
must freeze it before audit; this component cannot authenticate a caller's
external approval process merely by hashing a document at scoring time.
Multiple actual 1/8-MPI owners may be supplied; duplicate owners are errors.
There is no native executable invocation or build inside this verifier.

The original owner result must be PASS with a native-run, mode-run manifest,
the declared MPI count and OMP1. Original wrapper artifact hashes authenticate
`native-output/input_receipt.json`, every original parent and every rank's
exact `.hex.csv` events. The complete native/application include graph must
still match its receipt. A source-only preflight is never an actual event.
P02's event set must exactly equal all REACTED parents in a complete,
nonduplicated original-parent roster. No survivor can become an extra event.

The budget file is `M3_P02_NUMERICAL_BUDGET.json`, schema
`srcMoon-M3-P02-approved-mechanical-budget-v1`. Its four approved decimal values
are parsed from the pinned file; there are no numeric cap defaults in code.
Every native parent/ion/electron limit and its native speed guard must equal
the approved values. A loose native limit cannot approve itself.

## Independent quantities and frame scope

The exact source reference mass is reconstructed from neutral isotope masses
and integer charge times electron mass. The physical rest mass is read from
the authenticated selected state, whose source identity/role is checked. The
scope requires 23Na neutral ground, 23Na+ ground and free-electron states, and
the explicitly named source-effective-photon convention. No missing ion mass
is replaced by a reference mass. This check uses the supplied central source
values; source uncertainty propagation and broad model/data completeness
remain their own gates.

The relative inertia denominator is the physical mass; the relative
acceleration denominator is the reference mass. Kinetic error means the
relative difference between classical and relativistic kinetic energy with
the SAME reference mass, divided by relativistic kinetic energy. Reference-
versus-physical mass error is checked separately and not counted twice.
The exact zero-speed limit is0. The implementation uses an exact rational
squared inequality for its acceptance decision, so tiny positive speeds
cannot underflow to a false zero. Decimal output is diagnostic only.

Native parent/electron inertial vectors remain exact hexadecimal rationals.
The actual committed ion vector is recorded in LSO; its declared inertial
vector is recovered by exact inversion of the recorded affine R/offset.
It does not assume a rounded rotation's transpose is an exact inverse. The
proper-rotation bounds are checked separately, and an inadmissible transform
cannot produce a passing inertial-speed claim. The actual original input
must explicitly name SUN-origin inertial axes and common-time transport.
This is a binding to the native frame operator; it does not independently
verify SPICE ephemerides or frame classification.

The speed cap is checked in both the declared inertial frame and actual LSO
event coordinates, with separate labels. Relativistic discrepancy is assessed
only in the declared inertial frame. The worst-case declared maximum speed
is checked independently too, matching production's conservative admission.
Actual compiled parent/ion PIC masses additionally receive the same physical
mass error checks. No observation, guessed velocity, fitted residual or
native verified boolean supplies an oracle.

## Results and supporting tests

A NEW output directory contains `result.json`, `comparison.csv`, `checks.log`
and `run_manifest.json` with exact owner/source/budget/verifier identities.
These are scalar guard checks, so no physical comparison plot is implied.
PASS/0 means the declared criteria were evaluated and met; FAIL/1 means an
evaluated cap or equality failed; ERROR/2 means corrupt/malformed/contradictory
source/owner infrastructure; SKIPPED/77 means necessary evidence was unavailable.
A later missing input preserves an earlier evaluated FAIL. No actual event
means source/guard checks alone cannot establish actual event-speed verification.

Run the maintained software/source controls without an AMPS build:

```sh
/usr/bin/python3 -B srcMoon/test/sodium_cohort/test_p02.py
```

Analytical arithmetic and affine fixtures are labelled software controls.
The real source-dependent class declares `requires_authoritative_source=True`.
It uses the retained approved Na profile and leaves unavailable ion entries
unavailable. Tests never fabricate a successful native physical run or use
synthetic observations.

For the separate complete homogeneous-cohort conservation/statistics/timing
audit, retain its frozen case plan and original launch registry, then run:

```sh
/usr/bin/python3 -B srcMoon/test/sodium_cohort/audit.py \
  --plan /absolute/frozen-cohort-plan.json \
  --runs /absolute/native_runs.json \
  --output-dir /absolute/NEW-independent-cohort-audit
```

That audit additionally checks original-parent membership, source ledger,
source-keyed first-hit times, endpoint statistics and matched three-level
1/8-MPI refinement. Its formal predecessors, unchanged ensemble/numerical
thresholds and source/operator gates remain independent of this P02 audit.

### Correct dimensionless rotation scale

The recorded proper-rotation criterion is the column Gram relation `R^T R=I`.
Every residual uses the product of unit-vector norms, hence dimensionless
scale1, including off-diagonal entries whose target is0. The determinant
comparison also uses scale1. Scaling a tiny off-diagonal dot product by itself
would incorrectly demand an error near1e−32 when the original native bound is
0+1e−15. The first43-event audit exposed precisely that verifier error and is
retained. Correcting the scale does not alter the native bound or coefficients.
Maintained controls include the actual rounded matrix, exact1e−15 boundary,
a rational value above that boundary and an improper determinant.
