# Upstream core integration, 2026-11-17

**Core integration, workspace rebuild and all required preservation gates are
complete.** All 843 previously passing registered CORE identities remain PASS.
The full SEP physical-development/validation plan remains open; see
[STATUS.md](STATUS.md).

## Subsequent rebase resolution onto `6bd40ea`

The user subsequently committed the integrated SEP work as
`7984c03ac5d4aaa6f9d6e7c7b4b189716effc798` and initiated its rebase onto
`6bd40eacde14b810a199151c1c884513c7fe7cc9`. This resolution supersedes the
two whitespace cleanups and appended test registrations described in the
historical delivery record below:

- `src/meshAMR/cut_cell.hpp` and `src/meshAMR/meshAMRgeneric.h` use the exact
  new upstream bytes. The conflicting versions differ only in whitespace;
  no algorithm, token, tolerance or SEP hook is lost.
- `MakefileTest/Table` uses the exact new upstream bytes. The clean automatic
  merge had appended duplicate copies of the four native core campaign
  stanzas; all four upstream originals remain, with their original profiles,
  checks and MPI settings.
- Automatically merged `src/pic/pic.h` and `srcSEP/BoundaryInjection.cpp`
  retain exactly the replayed commit's bytes, including the SEP lifecycle
  hooks and source preflight.

Focused verification passes: upstream byte comparisons, complete parsed
registry comparison (**162** entries), and the existing **4/4** table tests:

```sh
PYTHONPATH=utility/nightly PYTHONDONTWRITEBYTECODE=1 AMPS_TABLE_HOST=passion \
  python3 -B -m unittest discover -s utility/nightly/tests -p test_table.py -v
```

`passion` is this machine's declared table host. A before/after parser check
finds two registrations per affected campaign in the saved automatic merge
and one in the resolved table. Existing host-selector variants are preserved.
Evidence, original conflict stages and index backup are in
`/tmp/sep-rebase-conflicts.j5ge1br4`. These checks do not rerun the historical
full native/core campaigns or qualify the new upstream target physically.
Only the conflict resolutions, duplicate removal and these status notes are
staged; the rest of the replayed changes are preserved. The rebase remains
paused and no commit is created, as requested by the user.

## Pinned source revisions and preserved work

- Original local base: `c65494f5a6eff5988136db181abd3750b7cec8a7`.
- Verified upstream `origin/master`: `74974c2b2fc1b4cf09a73c320064a22370c0ca94`.
- The major incoming shared-core change is `873ebb5`; the selected target also
  includes subsequent upstream corrections and optional diagnostic hooks.
- Git ancestry was verified using the actual fetched commit graph.

The existing checkout contains substantial uncommitted SEP work and a separate
user linker change in `Makefile`. Before integration, 264 modified/new source
and document files, a binary Git diff and the complete working status were
backed up under `/tmp/srcsep-upstream-merge.1q31gpks`. External data, ignored
build products and unrelated user directories are not merge inputs. Source
integration is prepared in `merged-candidate` under that directory; the
original `.git` metadata is unchanged.

## Incoming scope and upstream reuse

The prepared payload has 65 upstream integration paths, one authenticated
external test fixture and one retained-local compatibility amendment
(67 delivery paths):

- **60 files are byte-identical to the pinned upstream version.** These include
  all incoming shared code under `src`, except the two whitespace-only files
  and overlapping header described below; configuration generators; the
  corrected shared species declaration; new native core fixtures/profiles;
  and the upstream macro-consistency checker. One upstream Enceladus
  `shape_manifest.json` fixture is included because the unchanged shared
  native-fixture provenance validator explicitly requires it; no Enceladus
  application implementation is changed.
- **Two files retain upstream code with four whitespace spans cleaned.** The
  upstream `cut_cell.hpp` has two trailing-space lines and
  `meshAMRgeneric.h` has two space-before-tab lines. No token, expression,
  branch, tolerance or algorithm is changed by those cleanups.
- **Two files merge local and upstream edits using the actual common base.**
  `src/pic/pic.h` retains upstream declarations and the local optional
  field-line lifecycle/source interfaces. `srcSEP/BoundaryInjection.cpp`
  retains the local source preflight and upstream whole-template zeroing.
  Both three-way merges are clean. An independent agent's `pic.h` merge has
  identical SHA256 `00c4171dddee67da8dabb287b1e592641a8b39deff43ef6c51fcdbe760ad902c`.
- **One shared test table receives four complete, exact upstream stanzas:**
  `unit_test_cutcell_tetra`, `unit_test_cut_surface_mover`,
  `unit_test_cut_surface_mover_specular`, and
  `unit_test_gas_surface_interaction`. Their native source, input profiles,
  tolerance checks, MPI ranks and negative controls accompany them. Unrelated
  application/nightly changes are outside this core-only integration.
- **One required existing external test fixture is delivered exactly.**
  `CODEX_REDUCED_SHOCK_TASK.txt` is the inert root-discovery/retained-contract
  fixture declared as untracked in `src/test/native_fixture_provenance.json`.
  The reduced-shock owner unconditionally copies it before its 97-case run.
  Its authenticated source is `/home/vtenishe/Mars2/AMPS/` with 15,578 bytes
  and SHA256 `7a65f3c0c244d01c0e44c81d9fd7194bd0298c4d892bbba5b903fab36fac3d49`.
  It was absent from the development workspace; both private full-core runs
  use these same bytes. It is test support, supplies no runtime physical
  inputs, and does not supersede the current user's scientific decisions.
- **One retained shared test-header amendment restores compatibility.** The
  local optional `ExecutionContext::testId` gains an explicit empty member
  initializer. The completed full campaign exposed the unchanged upstream
  two-field aggregate client failing `-Werror=missing-field-initializers`.
  The default preserves that aggregate API without editing the upstream
  application or changing any oracle. Strict fail-before/pass-after compilation
  and aggregate/layout/default checks pass; the layout remains 128 bytes with
  alignment 8 on this host. Header SHA256 is
  `6484ba2665ef38c9e60bf073bc3abd5742c61ea02938cb9fbb1d2988ab6c9657`.

The incoming implementation supplies exact cut-cell geometry, watertight
surface crossings, runtime gas/surface laws, mover and memory repairs, sparse
AMR particle-list sweeps, deterministic hybrid temporary-list handling and
consistent generated macro ordering. These are reused from upstream rather
than replaced by SEP implementations.

The complete selected/excluded-path list, old/local/upstream/candidate hashes
and merge method are in `candidate-integration-manifest.json` in the evidence
directory. `prepare_candidate.py` is the reviewable integration procedure.

The shared `.git` index is intentionally unchanged. Newly imported `.input`
profiles and the provenance fixture match broad existing ignore rules; the
manifest lists them explicitly. Packaging a future commit must include those
exact upstream fixtures explicitly (for example with `git add -f`), rather
than relying on an ordinary untracked-file listing alone.

## Local core contracts that still need extensions

Upstream has not changed `field_line_particle_split.cpp`, `pic_field_line.cpp`
or the eight locally modified `sep_common` files relative to the original
base. Inspection found no new upstream replacement for the following required
contracts:

- True fresh-injection and exact-split lifecycle hooks, with preflight before
  particle/RNG/list mutation, stable parent/child identities and preserved
  continuation history. Generic particle cloning remains byte-exact.
- Explicit field-line source ownership and numerical-mode preflight, complete
  fresh particle templates, and endpoint Cartesian-coordinate consistency.
- Positive-weight resampling of existing phase points with declared moments;
  the SEP relativistic six-moment policy remains application-owned.
- Persistent keyed stochastic state, validated background snapshot restore,
  coefficient contracts and the existing selectable test registry in the
  shared SEP library.

Upstream allocation/clone signatures and extension-storage alignment remain
compatible. The new upstream pre-deallocation particle observer is a passive
diagnostic; it does not replace a birth/split hook. The upstream AMR sparse
mover does not replace the existing field-line dispatch. Its fallback guards
and the upstream allocator counter are retained.

The application retains its physical source/mover/turbulence decisions. No new
scattering law, isotope, plasma input or closure-domain limit is inferred by
this software integration.

## Regression and acceptance record

Testing uses isolated checkouts because native core owners remove/reconfigure
their own generated `build` directory. They must not run against the shared
SEP development checkout. Dependencies, site configuration, executable hashes
and source identities must accompany each receipt. Separate checkouts use
separate native locks and build products.

1. Run the complete registered `CORE-*` campaign against the unmodified pinned
   upstream tree with an explicit machine profile and writable artifact paths.
2. Run the same campaign against the frozen merged candidate. Compare stable
   suite/case identities, including missing cases and PASS-to-SKIP/FAIL/ERROR
   changes. Case totals alone are insufficient.
3. Account separately for historical PASS identities and any upstream or local
   prerequisite nonpasses. Dated core documentation reports 843 PASS, one
   strict-coverage FAIL, 240 SKIP and one CME refinement-timeout ERROR. Newer
   local receipts also contain independent shared-library verification
   failures. A fresh receipt, not a documentation count, is authoritative.
4. Execute the four newly imported native core campaigns, preserving their
   original physical/numerical assertions, required one/four-rank comparisons
   and expected-failure controls. Missing or truncated reports cannot pass.
5. Rebuild the SEP application with the candidate core. Re-run its bounded
   native regression, real host/PIC timing, translated-origin coefficients,
   separate-process checkpoint continuation and both CLI campaigns.
6. Recheck source hashes, untouched local work, whitespace and the exact
   integrated payload before applying the reviewed source changes to the
   shared workspace.

The **core PASS-preservation gate passes**. Complete pristine upstream reports
843 PASS / 1 FAIL / 240 SKIP / 1 ERROR over 1,085 cases in 51 suites. The first
merged receipt captures the strict context-initializer compilation failure;
it is retained unchanged. After the narrow shared-header amendment, exact
scoped retries of the SEP-common and reduced-shock owners execute 240 cases:
153 PASS / 87 SKIP / 0 FAIL / 0 ERROR. The explicitly reconciled core receipt
then matches the pristine case/status map exactly: all 843 historical and
current baseline PASS identities remain PASS, with no missing or added cases.
The unchanged nonpasses are `T12-COVERAGE-STRICT-001` and
`CSWC-FULL-CSWC0628-MASS-TIME`, plus the 240 recorded prerequisite/qualification
SKIPs. No tolerance, test selection, coverage status or timeout is weakened.

The comparison is `core-identity-comparison.json`; the reconciled receipt is
`candidate-core-reconciled/amps-core-tests.json`, with raw attempts and exact
two-owner replacement provenance preserved. Source/dependency/profile/fixture
rehashing and independent identity/reconciliation corruption controls pass;
see `core-final-audit.md` and `core-final-independent-audit.json`. The
post-amendment SEP regression and independent audit pass, and the accepted
source payload is published.

The merged-core SEP build has succeeded with archived executable SHA256
`cc6ec6b7be9a23636eab85b94aded66a63e1c41465d4994f38e37f6568900f72`.
The combined native gate passes 121/121; HOST passes 3/3, translated-radius
cases 2/2, fresh-process checkpoint pairs 3/3 with identical final bytes, and
spectral ownership 4/4 at each of N=128 and N=129. Native CLI passes 9/9 and
production CLI software retains 18 PASS / two explicit qualification SKIPs.
These receipts are in `sep-evidence-merged` under the integration evidence
directory. The four new upstream native campaigns also pass: 32 expected
launch outcomes, 578 native assertions and 13 exact one/four-rank byte
comparisons. The two gas missing-input controls exit nonzero with their required
messages. Their strict report is `new-core-evidence-merged/checked.json`,
SHA256 `af6584c86905bb0fb108e512b6268c21ba2e4cd3ea9254719a04ba500043b464`;
the final audit confirms all 327 copied source files unchanged in that phase.
The subsequent optional-context amendment changes no source used by these
geometry/surface campaigns. The post-amendment SEP executable is separately
archived as `sep-evidence-compat/amps.merged-core-compat`, SHA256
`13dd767b86ed65174af628cba6f70a875a1b486c914bbc4802452ce89e6bf1e1`;
its five generated-configuration hashes match the preceding SEP build. All
nine native/CLI/restart/spectral campaigns and the independent raw-evidence
review pass.

## Published workspace and final executable

All 67 delivery paths match the accepted candidate hashes; see
`publication.json` and `integration-acceptance.json` in the evidence root.
The original `.git` index/refs remain unchanged and the source changes are
uncommitted. `pre-publication-core-backup` retains every replaced original.
The required external fixture is delivered as declared test support.

The existing application was refreshed and rebuilt with:

```sh
./Config.pl -application=test/sep_parker_spiral__field_line
make -j16 amps
```

The original Makefile/linker edit, Makefile.conf, `.amps.conf`, compile-time
input and effective Makefile.local settings are preserved. The preceding
workspace build is archived at
`build.before-core-74974c2-20261117T223520Z`; configuration/executable backups
are in `workspace-pre-rebuild-backup` under the evidence root. No old build is
deleted. The actual workspace `amps` and archived
`sep-evidence-workspace/amps.workspace-core` both have SHA256
`dcad7ae35263e61cf123a8b7022235c337eb34231be30e4bd4180d0ff4b3f1c6`.

All nine final workspace campaigns exit 0: 121 combined native cases, three
actual driver cases, both radius cases, three fresh-process checkpoint pairs,
nine CLI cases, 18 production-CLI software passes with two explicit
qualification SKIPs, and four spectral ownership cases at each of N=128 and
N=129. Each restart pair has identical final bytes. The reproducible commands,
child reports and hashes are in `sep-evidence-workspace/regression-summary.json`.
The independent final audit is `sep-workspace-independent-audit/summary.json`
and its raw CLI/clock supplement; it confirms all 67 published hashes and
preserved settings. Production/SWMF execution and the original broader physical
qualification requirements remain open.
