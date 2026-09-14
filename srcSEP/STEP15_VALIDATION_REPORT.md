# Step 15 implementation and validation report

Date: 2026-09-13

## Implemented evidence

| Class | Case | Result in this environment |
|---|---|---:|
| Numerical verification | `VAL01` Parker manufactured solution | PASS |
| Cross-mover verification | `VAL02` matched Dmumu/MFP campaign | PASS |
| Cross-model verification | `VAL03` independent finite-volume solver | PASS |
| Coupled integration | `VAL04-SWCME` actual SWCME-to-srcSEP replay | PASS |
| Coupled integration | `VAL04-SWMF` real native SWMF replay | INCOMPLETE |
| Observational validation | held-out spacecraft events | INCOMPLETE |

The overall Step 15 release status is therefore **INCOMPLETE**. This is the
scientifically correct status: the source-only numerical successes do not
stand in for real SWMF output or spacecraft observations.

## Source-only command and observed metrics

`./test/run_step15_tests.sh` completed successfully under strict warnings,
AddressSanitizer, and UndefinedBehaviorSanitizer. Its internal `--release`
negative test correctly returned nonzero. Representative deterministic metrics
were:

- `VAL01`: mean displacement 0.568 standard errors from the exact mean,
  variance relative error `6.28e-4`, momentum relative error zero;
- `VAL02`: absolute mean-pitch differences `0.00625` and `0.00607`, Dmumu/MFP
  mean-square-displacement relative errors `0.123`/`0.00200`, onset/peak/fluence
  relative differences `0.0714`/`0.327`/`0.0459`;
- `VAL03`: pitch-density L1 error `0.0343`, mean-pitch error `4.23e-4`, and
  second-moment error `0.00234`;
- `VAL04-SWCME`: zero characteristic arc/momentum discrepancy to reported
  precision, positive background density, and one active physical source.

The authoritative machine-readable result of a retained execution is produced
at runtime because it records that execution's input checksums, source-tree
digest, compiler version, platform, complete configurations, seeds, output
schemas, and metrics.

## Required remaining execution

A release campaign must supply a completed
`validation/manifests/swmf_replay.template.json` derivative produced by a real
native SWMF-to-srcSEP run. It must also supply completed held-out event manifests
whose genuine spacecraft products jointly evaluate onset, anisotropy, spectra,
fluence, decay, and multi-spacecraft longitude. The manifest validator reopens
the local data bytes and verifies every SHA-256 before accepting either claim.

The enclosing AMPS `Makefile.conf`, generated PIC types, MPI executable, SWMF
run products, and spacecraft archives are absent from this environment. Native
build/MPI/coupling and observational closure are consequently not reported as
passed.
