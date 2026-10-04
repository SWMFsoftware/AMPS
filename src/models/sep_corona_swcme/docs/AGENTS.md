# AGENTS.md — AMPS CME background program

## Applicability and precedence

These instructions apply to the entire repository. More deeply nested
`AGENTS.md` files add application-specific rules for their subtrees. Direct
user instructions take precedence. Preserve unrelated baseline, generic PIC,
application, CLI, particle, and legacy SWCME rules; this file changes the CME
program's scope and ordering rather than replacing those contracts.

The normative execution and physics documents are:

- `CODEX_CME_PLAN.md` (active execution/evidence ledger);
- `src/models/sep_coronal_cme/docs/CME_SURFACE_LAUNCH_ROADMAP.md`, revision
  1.3 (ordered gates and acceptance scope);
- the maintained modules under `src/models/sep_coronal_cme/model/` (coronal
  ambient, EOS, geometry, local shock, configuration, architecture, tests);
- `src/models/sep_corona_swcme/model.md` (continuous Corona–SWCME handoff and
  regional spatial CME closure); and
- the maintained public contracts under `src/models/swcme` and
  `src/models/sep_common`.

## Workspace and preservation boundary

- Work only in `/home/vtenishe/Mars2/AMPS`.
- Treat `/home/vtenishe/Mars1/AMPS` as read-only reference material. Never
  edit, clean, build, commit, or discard its working-tree changes.
- Preserve all user-owned dirty-tree changes. Do not commit or push.
- Transfer Mars1 code only after classifying it as background-only or safely
  separating mixed-file background edits from particle edits. Revalidate every
  transferred algorithm independently in Mars2.
- Keep baseline particle-only files and particle sections of mixed files.
  Preserve generic PIC behavior, existing particle tests/schemas/CLI, and
  legacy SWCME.

## Mandatory program order and authorization

The only active gate is the earliest unqualified gate below:

1. **G1 / BG3D-0 through BG3D-9:** complete and qualify the full spatial CME
   plasma/IMF background in srcSEP3D, with zero allocated and zero injected
   particles, from a prescribed low-coronal launch through a continuous
   SWCME handoff and actual propagation to 1 AU.
2. **G2 / BGSEP-0 through BGSEP-3:** only after G1 passes, couple and qualify
   the same shared background in srcSEP, still with zero particles.
3. **G3 / P3D-0 through P3D-4:** only after G2 passes, implement and qualify
   srcSEP3D particle modeling.
4. **G4 / PSEP-0 through PSEP-3:** only after G3 passes, implement and qualify
   srcSEP particle modeling.

Follow stages within each gate sequentially. Do not begin srcSEP background
work or any new particle/source/feedback implementation during G1. Historical
L0–L13 and particle assertions retain their identities and remain open until
their assigned later gate passes.

The current task terminates at G1. Even after every BG3D-0 through BG3D-9
acceptance gate passes, update the plan/READMEs, report changes, exact tests,
evidence paths and limitations, and stop. G2, G3 and G4 require new explicit
user authorization. Any required skip or unresolved failure leaves G1 open.

## G1 physical completion contract

- A front, local jump, or quiet ambient behind a front is not a full CME
  background. Implement declared, sampled spatial ambient, shock, sheath,
  contact, ejecta, wake/transition regions and their plasma/IMF state.
- Report shock position, normal, normal speed, upstream-relative inflow, fast
  Mach number, obliquity, compression, both-sided states, admissibility, and
  conservation residuals. Subfast geometry remains typed geometry; never add a
  compression floor or convert solver failure into ambient state.
- Select and document analytical regional closures and independently check
  mass, magnetic flux, induction, thermodynamics, interfaces, material-map
  Jacobians, and prescribed force/heating/work budgets. Separate numerical
  residuals from physical/model residuals.
- This is a prescribed analytical/semi-empirical model, not a global MHD
  solver or an eruption-instability prediction. A spheromak is required only
  if the selected ejecta closure itself requires one.
- Preserve one composite authority, event identity, geometry/material state,
  and generation lineage. The handoff must not reset radius, speed, shape,
  ambient normalization, mass, or flux.
- Background epochs and readiness must not depend on source planning,
  particle allocation/initialization, or particle-wave feedback.
- Shared physics belongs under `src/models`; application adapters remain thin.
  Public physics libraries must remain usable without PIC, MPI, or application
  headers.

## Native build preflight — required before every rebuild

Before **every** native AMPS rebuild, perform and record this sequence exactly:

1. Confirm the working directory is `/home/vtenishe/Mars2/AMPS` and verify no
   build or test process is running.
2. Remove only the generated directory with `rm -rf -- build` from that
   confirmed root. Do not remove any source, site configuration, library,
   Mars1 path, or broader directory.
3. Regenerate the selected application's configuration and its production
   hooks through the repository's established configuration workflow.
4. Compile through the established top-level workflow with `-j16` (for
   example `env MAKEFLAGS="-j16" make amps`).

Preserve `Makefile.conf`, site libraries, application inputs, and local site
settings. Inspect real scripts/decks and eliminate hard-coded Mars1 build or
output paths. Never treat a direct source-only compile as a native build.

## Validation and evidence

- Add detailed physical and numerical reasoning comments to every new or
  modified CME/background implementation, including earlier code touched by
  this task. Comments must state equations/closures, SI units, frames and sign
  conventions, assumptions/approximations and validity limits, algorithms,
  tolerances, numerical methods and typed failure behavior; do not merely
  paraphrase statements. Cover regional classification/shock diagnostics,
  launch/handoff/outer propagation, state ownership/epochs and, when reached,
  native storage and MPI synchronization.
- Maintain relevant READMEs alongside the code. They must distinguish
  implemented from independently validated behavior and document input
  parameters, commented examples, exact build/run/test commands, validation
  criteria, outputs and known limitations. Planned native/MPI behavior must
  remain explicitly unimplemented until its stage supplies real evidence.
  Apply this documentation requirement to every subsequent change while
  preserving baseline particle and unrelated documentation.
- Tests must use independent equations or fixtures rather than
  production-generated expected values.
- G1 requires actual evolving one-rank and four-rank native runs with zero
  allocated/injected particles. Retain multi-epoch region, shock, plasma,
  owner-readback, received-ghost, temporal-slot, derivative, generation,
  mesh/layout, input, closure, and binary evidence.
- Keep the short handoff-smoke deck separate from the full low-corona-to-1-AU
  deck. Initialization, unit tests, final-only receipts, and short outer-only
  propagation cannot close G1.
- Record exact commands and PASS/FAIL/SKIP/ERROR outcomes. Missing data/tools,
  incomplete implementation, and physical rejection are distinct. Never
  weaken a gate or tolerance, invent a result, or mark a particle assertion
  passed with source-free evidence.
- Update `CODEX_CME_PLAN.md` after every meaningful stage with implementation,
  validation, provenance, residuals, limitations, and remaining work.
