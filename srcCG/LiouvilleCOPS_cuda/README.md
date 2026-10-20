# LiouvilleCOPS CUDA + CuPy backend (development)

This tree ports the gravity, masking, pressure and analytic-Jacobian work of
`../LiouvilleCOPS` to two RTX A6000 GPUs. The CPU application stays the reference. The plan
is in `ROADMAP.md` (revision 2) and progress is in `STATUS.md`.

The primary fitting target is the free 29-parameter baseline (B29) and its NG-only and cliff
template extensions (N30/C30/NC31). The fixed-baseline higher-harmonic modes are
secondary.

Stages completed:
* C00: CPU reference and the approved tolerance policy;
* C01: environment;
* C02: backend contract, device layout, memory policy, provenance;
* C03: exact FP64 polyhedral gravity on one GPU.

| File | Purpose |
|---|---|
| `ROADMAP.md` | development plan, revision 2 (incorporates `AMMENDED_ROADMAP.md`); `ROADMAP_v1.md` is the original |
| `AMMENDED_ROADMAP.md`, `AMENDMENT_REVIEW.md` | the amendment (revised fitting scheme, gravity admission) and the review comments on it |
| `quick_guide.md` | **start here**: exact commands to set up, run GPU fits (production and CPU comparison), compare, replay and plot |
| `RUN.md` | the two interpreters, environment rebuild, doctor, C00 capture, C02 and C03 evidence commands |
| `STATUS.md` | per-stage evidence, measurements, limitations |
| `VALIDATION.md`, `validation/tolerances.json` | numerical acceptance policy (approved 2026-10-19; Jacobian and end-to-end fit items provisional) |
| `GRAVITY_REFERENCE.md` | polyhedral-gravity 3.3 as implemented: equations, singular cases, summation order, the HEAL defect, GPLv3 obligations, and the CUDA translation (§8) |
| `cops_cuda/` | the backend package; module guide and C02 design in `cops_cuda/README.md` |
| `cops_cuda/kernels/` | CUDA C sources compiled by NVRTC; conventions in `cops_cuda/kernels/README.md` |
| `tools/c00_capture_reference.py` | C00 reference capture (CPU reference interpreter) |
| `tools/cuda_doctor.py` | C01 device diagnostics (CUDA environment) |
| `tools/c02_device_load.py` | C02 evidence on the campaign run: subset/full load, cycles, capacity refusal, second GPU |
| `tools/derive_outward_surface.py` | exact global winding reversal of a consistently inward closed surface (for `SHAP5_Volume.bdf`), with a provenance record |
| `tools/c04_gravity_equivalence.py` | C04: full-field CPU-builder vs GPU gravity on the same body and targets, plus the corrected 100k-vs-788k slope comparison |
| `tools/c03_gravity_check.py` | C03 evidence: GPU gravity against the C00 fixtures, package state, difficult samples, determinism, full 788k field |
| `tests/` | `pytest` suite; `tests/synthetic_run.py` writes small valid base runs for the tests |
| `third_party/` | the pinned polyhedral-gravity 3.3 sdist and wheel, with their digests |

`cops_cuda/gravity.py` and `cops_cuda/kernels/gravity.cu` are translated from GPLv3 source and
are GPLv3; see `GRAVITY_REFERENCE.md` §7–8.
