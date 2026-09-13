# AGENTS.md — srcSEP field-line transport application

## Applicability and precedence

These instructions apply to the entire `srcSEP` tree. Read them before editing any file. If a more deeply nested `AGENTS.md` is added later, its instructions apply to that subtree in addition to this file. Direct user instructions take precedence.

The detailed implementation and validation specification is provided outside this source directory in:

- `../docs/srcSEP_field_line_transport_improvement_plan.md`
- `../docs/srcSEP_field_line_transport_improvement_plan_with_turbulence_and_cli.docx`
- `../CURRENT_CLI_CONTRACT.md`
- `../TASK_CHECKLIST.md`

When this source tree is copied out of the handoff package, retain those documents or update the paths in this file.

## Authorization boundary

The plan describes a multi-stage program. It is not authorization to implement all stages in one turn.

- Implement only the step, priority, or bounded task explicitly requested by the user.
- Complete tasks sequentially when the user gives an ordered list.
- Stop after the requested task or requested batch, report its gate, and wait for authorization before beginning the next unrequested step.
- Diagnostic inspection, focused test execution, and documentation directly needed for the requested task are allowed.
- Do not commit, push, open a pull request, publish results, or download large observational datasets unless the user asks. Propose a commit message after a completed change when useful.

## Scientific scope

`srcSEP` is the field-line SEP transport application. Its production interface will contain exactly three particle movers:

1. `parker`: pitch-angle-averaged Parker transport along magnetic field lines.
2. `fte-dmumu`: focused transport with pitch-angle scattering specified through `D_mumu` and its required derivative.
3. `fte-mfp`: focused transport with scattering specified through parallel mean free path `lambda_parallel`.

Full three-dimensional particle transport belongs in a separate application. In particular, Cartesian/Boris motion, gradient/curvature drift motion, cross-field particle displacement, mesh-cell particle motion, and movers whose defining behavior is leaving the assigned field line are outside the final `srcSEP` production scope.

Do not confuse particle dimensionality with geometry. The following remain required in `srcSEP`:

- magnetic field lines embedded in three-dimensional space;
- three-component positions, tangents, magnetic fields, velocities, and shock normals;
- field-line curvature and focusing derivatives;
- moving field lines and evolving segment geometry;
- SWMF field-line, plasma, IMF, turbulence, and shock import;
- shock/field-line intersections and observer connectivity;
- interpolation and sampling along a spatially three-dimensional line.

The self-consistent Alfven-turbulence implementation remains a production physics subsystem. It is a turbulence provider/evolution model, not a fourth particle mover. Prescribed turbulence and SWMF-provided turbulence also remain supported choices where their ownership contracts are defined.

## Physics and data contracts

- Use SI units internally unless an interface explicitly documents a different unit.
- State units in public APIs, stored data, configuration, output metadata, and non-obvious calculations.
- Every particle step must consume one immutable background snapshot with a declared epoch, provider identity, field-line generation, validity interval, and configuration fingerprint.
- Analytic, SWCME, and SWMF providers must not overwrite each other's authoritative state.
- Imported SWMF state is immutable in read-only mode. A local-evolution handoff must copy it once and record the handoff epoch and provenance.
- Mover selection changes the transport/scattering algorithm only. It must not silently change the solar-wind, IMF, shock, turbulence source, injection normalization, units, time origin, or boundary policy.
- Coefficient selection is explicit. `kappa_parallel`, `D_mumu`, `dD_mumu/dmu`, and `lambda_parallel` must be obtained through documented providers/adapters with validity status and units.
- Accumulate particle-wave exchange in thread-local storage and apply it after deterministic reduction. Movers must not mutate shared wave state directly.
- Use one authoritative turbulence state per representation. In integrated mode, branch-integrated `E+` and `E-` are authoritative. In spectral mode, `E+(k)` and `E-(k)` bins are authoritative and integrated quantities are derived.
- Record signed energy/source ledgers for turbulence and particle-wave coupling. Limiters, floors, remaps, and rejected updates must be visible in diagnostics.
- Silent clipping, fallback, extrapolation, NaN replacement, or invalid-domain continuation is prohibited. Return an explicit status, subcycle under a documented policy, or fail with context.
- Preserve deterministic keyed random streams for fixed campaign seed, particle identity, operator, and event index. Do not use thread scheduling or MPI rank as physical randomness.

## Existing CLI is a compatibility surface

The production CLI already exists in `util/sep_cli.h` and `util/sep_cli.cpp` and is used by `main.cpp`. Preserve and test it before extending it.

- Do not introduce a second unrelated parser.
- Add test-selection fields and run modes to the existing `SEP::Util::CLI::Options` structure unless an explicitly approved refactor provides a compatibility adapter.
- Preserve current no-argument defaults until a user-approved migration changes them.
- Preserve both `--option value` and `--option=value` forms.
- Preserve case-insensitive values and the current Boolean vocabulary.
- Preserve current aliases, including `--cascase` and `--no-cascase`, during a documented deprecation interval.
- Preserve current last-occurrence-wins parsing and the current CLI-after-post-compile-input application order until a versioned precedence policy is approved.
- `-h` and `--help` must exit successfully before AMPS, field-line, particle, turbulence, MPI, or output initialization.
- Unknown/malformed options must never fall through to a long production run.
- Parser tests must verify applied production globals/function pointers, not only the temporary `Options` object.
- Preserve `--test-manager`, `--testmanager`, `--run-test-manager`, and `--no-test-manager` until their legacy behavior is mapped and tested.
- The SWMF-coupled library entry points must not parse process arguments or automatically run standalone tests.
- Keep `--turbulence-model` as the representation selector. If source/ownership selection is added, use a distinct option such as `--turbulence-source`.

Run PCLI01-PCLI12 from the plan before and after each parser, mover-selection, initialization-order, or configuration-output change. Add CLI01-CLI05 for the selectable test registry without weakening the PCLI suite.

## Required work sequence for each task

1. Inspect `git status --short`, applicable build files, configuration, code, tests, and documentation.
2. Identify existing local/user changes and do not overwrite or reformat unrelated work.
3. State the expected physics/software behavior and the focused completion gate.
4. Add or update the focused test so it can fail for the defect being addressed.
5. Implement the smallest coherent change that satisfies the requested task.
6. Add detailed comments to new code. Comments must explain physics, units, ownership, invariants, numerical limits, error paths, and compatibility decisions—not restate syntax.
7. Update relevant README sections, command examples, CLI help, configuration documentation, and output/restart schemas in the same change.
8. Run the focused test immediately.
9. Run related component tests and inspect generated evidence.
10. At the end of a requested batch or milestone, run the complete bounded regression suite.
11. Run `git diff --check`, inspect `git diff`, and review `git status --short` before reporting completion.
12. Report exact commands, PASS/FAIL/SKIP status, files changed, remaining limitations, and a proposed imperative commit subject.

A failed prerequisite gate blocks dependent work. Diagnose the failure; do not proceed by weakening assertions, tolerances, conservation limits, or statistical acceptance criteria.

## Testing interface and expectations

Preserve or implement these routine targets in the native AMPS tree:

```sh
make -j
make test-cli
make test-list
make test-case CASE=PCLI01
make test-group GROUP=turbulence
make test-turbulence
make -j test
```

Use the repository's actual target names when they differ during the transition, and document the mapping. A child test failure must propagate to `make` and the process exit status.

After relevant changes, also run the applicable strict-warning, ASan, UBSan, thread-scheduler, MPI-decomposition, restart, stress, and performance tests. Expensive stochastic/observational campaigns may be separate scheduled targets, but routine tests must remain bounded enough for development.

Every registered test must:

- have a stable ID and group;
- be individually selectable;
- report PASS, FAIL, SKIP, or ERROR;
- record metrics, tolerance, seed, configuration, duration, and artifact paths;
- return failure through the authoritative test result and process status;
- avoid hidden local failure flags and direct `exit()` calls from physics kernels;
- restore mutated global state or execute in an isolated fixture/subprocess;
- produce deterministic evidence for deterministic inputs.

Reference solutions must be independent where the plan requires independence. Do not call the production routine to generate its own expected values. Regeneration requires an explicit command, reviewed provenance, algorithm/version metadata, and a diff of old versus new metrics.

## Existing application-level baselines

The supplied `Makefile.test` identifies these relevant starting points:

- `test_SEP--Parker_spiral--FTE`: primary field-line FTE/Parker-spiral baseline.
- `test_SEP--Parker_spiral--ParkerEq`: Parker-equation baseline.
- `test_SEP--Parker_spiral--field_line`: field-line geometry/integration baseline.
- `test_SEP-SC-IH--FL--FL-attachment--Droge_diffusion_FTE`: supplemental FTE configuration after its setup is inspected.

Do not use Parker3D, Boris, Cartesian, drift, or perpendicular-diffusion particle applications as final evidence for the field-line-only production movers. Use them only to inventory or verify transfer to the separate three-dimensional application.

There is no clearly named complete self-consistent Alfven-turbulence target in the supplied test index. Implement and use TURB01-TURB20 and `make test-turbulence`; do not infer coverage from a test name alone.

## Code-change rules

- Follow the existing C++ language level and formatting in the touched subsystem unless a build-wide change is explicitly requested.
- Prefer narrow interfaces, immutable views, RAII cleanup, explicit enums/status types, and dependency injection over new mutable globals.
- Avoid function-pointer identity as a capability test; use a mover/provider descriptor.
- Keep shared field-line kernels independent of a particular background provider or mover.
- Do not duplicate a formula across movers. Move shared calculations into a tested kernel/provider.
- Validate indices, segment ownership, field-line bounds, physical domains, finite values, and configuration combinations at boundaries.
- Preserve MPI-rank-safe and thread-safe output. Only the designated writer commits shared output.
- Output files are committed transactionally: validate/preflight, write a temporary file, flush/check, and rename only after success.
- Do not make broad cosmetic changes in the same patch as a physics or numerical correction.
- Do not delete historical code until the destination/replacement is identified and its focused tests pass.
- Never use destructive Git commands to discard work. Do not modify unrelated files in a dirty worktree.

## Documentation requirements

Update the relevant documentation in every implementation change:

- root `README.md`: scope, three production movers, supported providers, configuration, examples, limitations;
- `test/README.md`: `make -j test`, CLI tests/groups, evidence, exit codes, adding tests;
- transport documentation: state, units, update order, boundaries, reproducibility;
- turbulence documentation: source/representation/ownership, equations, branch convention, boundaries, operator order, coupling ledger, restart, diagnostics, tests;
- mover documentation/comments: governing equation, required inputs, coefficient contract, frames, units, numerical scheme, limits, and validation IDs;
- migration manifest: every removed/renamed symbol, former location, new destination or replacement, and compatibility duration.

Examples must be executable and checked by tests where practical. Documentation must not claim observational validation when only synthetic or cross-model verification has been completed.

## Repository hygiene

Do not commit or package:

- object files, libraries, executables, core dumps, or profiler output;
- `.gcda`, `.gcno`, coverage HTML/data, sanitizer logs, or temporary build trees;
- generated `output/`, `test_output/`, or validation campaign results unless a small reviewed fixture is explicitly required;
- downloaded raw observational data when a downloader, provenance record, and checksum can reproduce it;
- Python caches;
- nested source archives.

Keep immutable reference inputs and expected solutions only when their license/provenance, generator version, units, schema, and checksum are documented.

## Definition of done for a requested implementation task

A task is complete only when:

- the requested behavior is implemented in the intended production path;
- new code has explanatory comments;
- focused tests fail without the fix and pass with it, or an equivalent evidence-backed test strategy is documented;
- relevant existing tests and the requested regression suite pass;
- README/CLI help/configuration documentation is current;
- no unrelated changes or generated artifacts are included;
- limitations, skipped tests, and environmental blockers are reported explicitly;
- the final response lists changed files, exact test commands/results, and a suitable commit message.
