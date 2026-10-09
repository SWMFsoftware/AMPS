# Mean-free-path development ledger

Specification: `MEAN_FREE_PATH_MODEL.md` version 1.2. Companion dependencies:
PARALLEL version 1.4 and PERPENDICULAR version 2.1. Companion integrity is
defined by `MEAN_FREE_PATH_MODEL_DATA/SHA256SUMS`.

| Stage | Status | Evidence and remaining boundary |
| --- | --- | --- |
| MF00 | PASS | Audited srcSEP/srcSEP3D MFP, kappa and D_mumu consumers. Confirmed `Chen2024AA` source equation and recorded `Tenishev2005AIAA` as code-derived legacy behavior. |
| MF01 | PASS | Read-only JSON source loader, raw scalar parser and manifest identity implemented. `reference_verification.py` passes 5,542 checks over all 34 fixture groups; `independent_physics_checks.py` passes 38/38. |
| MF02 | PASS | Relativistic proton/electron/alpha-capable kinematics and typed SI boundaries; F-KIN-01 representative checks pass. |
| MF03 | PASS | Lambda tags, Eq. (2), Eq. (3), round trips and singular-angle rejection implemented and tested. |
| MF04 | PASS | Explicit Eq. (9) evaluator and source identities implemented; U-11 requires input. Convenience construction of every companion-data row remains future work. |
| MF05 | PASS | q, epsilon, isotropic, printed variants, EPREM HalfD, Dröge nominal/effective distinction and stable He-Wan evaluation implemented; algebraic fixtures pass. Host SDE convergence is MF12 work. |
| MF06 | PARTIAL | Typed static local snapshot inputs and provenance implemented. Resolved wave-spectrum/provider delegation remains blocked until host integration. |
| MF07 | PASS | TS2003 ion and Zank forms implemented with named turbulence/variance conventions and stable Zank auxiliary evaluation. |
| MF08 | PASS | Electron RS and three distinct DT printings implemented with no default variant. |
| MF09 | PASS_WITH_BLOCKER | Bohm, Afanasiev and M-FLAMPA implemented/tested; SHOCK-PARASOL remains U-6 blocked. |
| MF10 | PASS_WITH_BLOCKERS | Fully specified GCR functional forms implemented with explicit VALUE normalizations and U-13/D-21/D-33 choices. Source-incomplete/time-series parameter rows remain typed external input or source gates. |
| MF11 | PASS_WITH_BLOCKERS | Confirmed Chen mapping and species-gated derived lambda implemented. Minoshima, Bobik, Jiang, Luo and ambiguous AMPS table entries remain typed gates. |
| MF12 | NOT_STARTED | User limited this task to a standalone library. `INTEGRATION_PLAN.md` defines the future srcSEP3D binding. |
| MF13 | NOT_STARTED | Requires host/science comparison campaign. |
| MF14 | NOT_STARTED | Requires a time-dependent host provider. |
| MF15 | PARTIAL | Standalone C++ regression, both companion scripts and ASan/UBSan pass. End-to-end SEP/GCR transport remains outside the standalone boundary. |
| MF16 | PARTIAL | API, schemas, tests, blockers and integration plan documented; host handoff awaits MF12--MF15. |

## Commands and current results

From `src/models/sep_common/mean_free_path`:

```sh
make clean && make test
```

Current results:

- C++ model/parser suite: 55 PASS, 0 FAIL;
- companion reference verification: 5,542 PASS, 0 FAIL across 34 groups;
- independent physics checks: 38 PASS, 0 FAIL;
- archive PIC/MPI dependency boundary: PASS; and
- ASan/UBSan C++ suite: 55 PASS, 0 FAIL with leak detection disabled because
  LeakSanitizer cannot operate under the environment's ptrace wrapper. The
  first leak-enabled attempt reported that environment limitation, not a leak.

For this run `mpmath` 1.3.0 was installed only under
`/tmp/amps-mfp-python`; it did not modify repository or system Python files.

No AMPS native rebuild or srcSEP3D run is part of this standalone stage.
