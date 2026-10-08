# Parallel diffusion implementation status

Updated: 2026-10-08. Working tree: uncommitted.

## Input contract

- Scientific specification: `PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md`,
  revision 1.3.
- Required companion bundle: not present in the checkout, `/home/vtenishe/Mars2`,
  `/tmp`, or the preserved `sources--100626--3am.tar` archive.
- Expected coefficient digests from the specification:
  `7cdc5ab9cddda0c7295a25bfcaff4ba3422002da12ac1b4660a98dcc64eaa762`
  (parallel) and
  `7371b6557f40c5efa7a48460a5f1b2971374a2c43d0eb3e209f0a8e62670356a`
  (perpendicular). They have not been verified against absent files.
- Approved scientific departures: none.

## Repository decisions

- Reusable core/API/tests: `src/models/sep_common/parallel_diffusion/`.
- Language: repository-supported C++17.
- Standalone archive: `libparallel_diffusion.a`; verification target:
  `make -f makefile verify`.
- Dependencies: C++ standard library only; no PIC, MPI, background, or
  application headers.
- Host adapters: planned in existing `srcSEP` and `srcSEP3D` adapter/provider
  seams; not implemented in this stage.

## Stage state

| Stage | State | Evidence or blocker |
| --- | --- | --- |
| PD00 | blocked | Repository/build/consumer map recorded, but the Section 22 companion bundle and its `SHA256SUMS` are absent, so the mandated reference-data import/digest gate cannot pass. |
| PD01 | in progress | API, SI kinematics, statuses, complete registry, provenance, fingerprints, parser bridge, and active function pointer implemented and locally tested. Formal acceptance remains downstream of PD00 and lacks the full-precision `benchmark_points.json` fixture. |
| PD02 | in progress | Five explicit prescriptions implemented with analytical rigidity slopes and focused tests. Formal stage acceptance remains downstream of incomplete PD00/PD01. |
| PD03 | pending | Requires accepted PD01. |
| PD04 | pending | Requires accepted PD02 and PD03. |
| PD05 | pending | Requires accepted PD01/PD03 and the absent exact coefficient/reference assets. |
| PD06 | pending | Requires accepted PD03/PD04 and independent backend fixtures. |
| PD07 | pending | Requires accepted PD03/PD04. |
| PD08 | pending | Requires accepted PD02--PD07 and supplied tables/assets. |
| PD09 | pending | Requires accepted PD02/PD04--PD08. |
| PD10 | pending | Requires accepted PD09. |
| PD11 | pending | Host plan is in `INTEGRATION_PLAN.md`; no mover binding claimed. |
| PD12 | pending | Release closure requires PD11. |

## Supported inventory

Implemented and selectable: `constant_lambda`, `constant_kappa`,
`power_law_lambda`, `broken_rigidity_kappa`, and `bohm`.

Registered and explicitly unavailable: `qlt_slab_spectrum`,
`qlt_slab_inertial`, `prescribed_lambda_mu_shape`, `broadened_slab`,
`nlpa_given_perp`, `nlgc_e`, `nlgce_n`, `nlgce_f_2014`,
`turbulence_adapter`, `wave_spectrum_adapter`, and `tabulated_parallel`.

The implemented explicit models return the scalar pair and analytical
log-rigidity slopes. Only constant models return an explicit zero spatial
gradient. Perpendicular values, `D_mumu`, non-constant spatial gradients,
batches, bounds, and fallbacks are unavailable.

## Acceptance evidence

| Command | Outcome | Evidence/notes |
| --- | --- | --- |
| `cd src/models/sep_common/parallel_diffusion && make -f makefile verify` | PASS, 24/24 | Warning-clean C++17 build; readable records on stdout and machine-readable `build/test-report.json`. |
| `g++ -std=c++17 -O1 -g -Wall -Wextra -Wpedantic -Werror -fsanitize=address,undefined -fno-omit-frame-pointer parallel_diffusion.cpp test_parallel_diffusion.cpp -o /tmp/parallel_diffusion_sanitize` followed by `ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 /tmp/parallel_diffusion_sanitize --json /tmp/parallel_diffusion_sanitize.json` | PASS, 24/24 | Strict compilation plus AddressSanitizer and UndefinedBehaviorSanitizer completed with exit status zero. |
| Same sanitizer binary with leak detection enabled | ERROR after 21/21 assertions | LeakSanitizer reports that it cannot operate under this environment's ptrace layer; this is not counted as a pass. |
| `cd src/models/sep_common && make verify` | PASS, 1/1 | Existing `SEP_COMMON01` dependency/member/unique-symbol boundary remains clean. |
| `python3 -m json.tool src/models/sep_common/SOURCE_MANIFEST.json` | PASS | Updated nested source/test/documentation manifest is valid JSON. |
| `python3 srcSEP/test/check_stage3_contracts.py` | FAIL before this module is examined | Existing clean tracked `srcSEP/mover.cpp` is present while the repository probe declares it retired. The file is outside this work and was not modified or removed. |

Primary repeat command:

```sh
cd src/models/sep_common/parallel_diffusion
make -f makefile verify
```

The generated report is not source-controlled. These results qualify the
implemented arithmetic and failure contracts only; they do not close PD00,
the missing full-precision fixture comparison, later models, or host transport
integration.

## Known gaps and next action

The first unfinished prerequisite is PD00 reference-data import. Supply the
exact `parallel_diffusion_model_data/` bundle described in Section 22, verify
its `SHA256SUMS`, run its independent `reference_verification.py`, and record
the environment/results here. Then close PD01/PD02 against the full-precision
fixtures before starting PD03. Do not reconstruct the absent JSON/CSV audit
bundle from rounded Markdown tables.
