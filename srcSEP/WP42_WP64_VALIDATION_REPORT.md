# WP42--WP64 validation report

## Executed in the delivered source environment

The following commands were run from the srcSEP directory:

```sh
./test/run_wp42_wp64_tests.sh
./test/run_wp31_wp41_tests.sh
make test-sanitizer
./test/run_step7_tests.sh
./test/run_step8_tests.sh
./test/run_step9_tests.sh
./test/run_step14_tests.sh
```

All returned exit status zero. The WP42--WP64 executable was built with C++11,
`-Wall -Wextra -Wpedantic -Werror`, AddressSanitizer, and
UndefinedBehaviorSanitizer. It reported one PASS for every work package.

Notable measured results were:

| Work package | Metric | Observed |
|---|---|---:|
| WP53 | Brownian parent-minus-children residual | `0` |
| WP53 | accepted / rejected / maximum depth | `3 / 2 / 2` |
| WP54 | coarse L1 error | `0.00385425` |
| WP54 | fine L1 error | `0.00122031` |
| WP54 | observed refinement order | `1.6592` |
| WP31 regression | observed combined-operator order | `0.969891` |

The source-integration assertions also confirmed that production contains a
persistent `TurbulenceRuntimeStore`, has no per-step local `State state;`, has
no legacy `WaveParticleCouplingManager` call in the production turbulence
adapter, drains coupling through the transactional API, uses the campaign seed
and stable source identity, evaluates distinct velocity derivatives, normalizes
shock processing with upstream relative speed, and copies per-event MFP
wave-frame momentum/branch metadata rather than shell-wide endpoints.
The same gate now rejects any reintroduction of the sampling manager's
`goto end`; root-only output is expressed as a structured conditional so native
C++17 compilers do not diagnose a jump across initialized local objects.

## What the controlled tests establish

- Runtime turbulence ownership and restart round trips preserve complete state.
- Particle--wave batches reject duplicates/invalid records without partial
  mutation and close particle-plus-wave energy for valid records.
- Curved-field velocity derivatives and continuity residuals are distinct.
- Shock validation, upstream processed flux, and Mach compression are bounded.
- Source IDs and random-purpose streams survive scheduler restart.
- Parker measure conversion, dynamic resonance, 90-degree closure, wave action,
  conservative cascade, and versioned profiles satisfy analytical invariants.
- Brownian bridges preserve paths; adaptive rejection is transactional.
- MUSCL/SSPRK2 advection improves under refinement and conserves the integral.
- Remap, mover failure, lineage statistics, detector response, native-evidence
  level, decomposition signatures, refinement fits, power plans, external
  manifests, and scaling gates accept valid fixtures and reject negative ones.
- The preceding WP31--WP41 suite still passes after the production queue-owner
  migration.
- The complete dependency-light sanitizer aggregate, Parker, both focused
  movers, and documentation/static-analysis gates remain green.

## Explicitly blocked evidence

The runner reported these conditions as `BLOCKED`, not `PASS`:

- **WP59--WP60 native:** no linked AMPS executable/matrix/restart command was
  supplied through `SRCSEP_NATIVE_GATE`.
- **WP63 external:** the source archive contains no authenticated real SWMF
  replay or held-out spacecraft manifests.
- **WP64 scaling:** no frozen multi-node hardware/environment campaign was
  supplied.

These are release gates. A source-only analytical PASS must not be promoted to
native, coupled, observational, or scaling evidence.

## Required completion commands

```sh
make test-wp59-wp64-native \
  SRCSEP_NATIVE_GATE='/reviewed/site/runner --matrix --restart --mpi'
make test-swmf-validation SWMF_MANIFEST=/evidence/swmf-replay.json
make test-observational-validation \
  OBSERVATIONAL_MANIFESTS='/evidence/event-1.json /evidence/event-2.json'
```

The native command must archive its executable checksum, compiler/flags,
configuration fingerprints, callback traces, decomposition/restart signatures,
and frozen performance environment. External commands must retain input bytes,
checksums, roles, preprocessing, exclusions, uncertainty, and metric results.
