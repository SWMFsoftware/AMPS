# November 17 upstream integration contract

This work integrates upstream `4bff09213b1bfa3f59e113ba8cbb2d56f22d30ef`
with working Moon changes based on `17d49e00d1d2ccf33b06d11a3e6e5ffef34ac810`.
The user authorizes integration, prefers upstream implementations, and requires
preservation of passing core regressions. No commit, push, reference replacement,
threshold relaxation or new physics choice is authorized by this integration.

## Inventory and production paths

The upstream delta changes cut-cell measures, particle movers, runtime surface
interactions and nightly infrastructure. Existing Moon extensions cover moving
terrain in the force-enabled mover, finite competing surface release, runtime
source/reservoir policies, actual event diagnostics and sample-buffer lifecycle.
Five dirty shared C++ files additionally implement runtime PSD geometry and
passive finite-release/cohort callbacks. Moon owns the geometry, laws and
configuration; shared hooks default to null. Retain upstream samplers, mesh,
particle storage, runtime surface APIs and lifecycle fixes. Adapt Moon callers
before considering another shared change. Do not replace whole upstream files
with an older implementation. Inventory and exact deltas are retained under
`test_output/srcMoon/upstream-core-integration-20261117T185033Z/`.

## Configuration, inputs and independent oracles

Use normal Config.pl, remove the isolated generated build tree before each
configuration change, then make -j with recorded build parallelism. Preserve
the active checkout's machine settings. Carry the upstream Makefile input
prerequisite together with the existing local SHARE linker setting. Regenerate
ignored `Makefile.test` from the new upstream parser/Table so new registrations
and canonical cleanup recipes are available; retain the eleven existing local
fast-wave include-path overrides. Back up its previous bytes. Other site
configuration remains byte-identical. Native core owners retain
their existing analytical/invariant/reference oracles and thresholds. Moon U01-
U26 and I01-I10 retain literal IDs, production probes, SI units and documented
frame/time conventions. I04's reviewed reference remains SHA-256
`ed875cefc604e4ae3580fb89b8c528146929b6b9be099b41b61697faee7f0080`.
Changed output is investigated and retained, never automatically accepted.

Moon supporting runs use eight MPI ranks/one OpenMP thread, except frozen
four-rank nightly and explicit one-rank comparison fixtures. Preserve existing
core case layouts and seeds. Existing D01/D03/atomic-rate inputs and SPICE are
read-only; native labels, package QA and exact hashes control interpretation.
No observational selection, hold-out inspection or scoring is part of a merge.

## Procedure and acceptance

1. Preserve dirty tracked and new source bytes, configuration, index and source
   identities. Work outside the active tree. Pin fetched upstream independently.
2. Run the unchanged native CORE selection on pristine upstream with only
   machine dependencies and real fixtures supplied. Freeze actual case IDs,
   source/input/oracle hashes, logs and results, including nonpassing cases.
3. Integrate local changes with a three-way merge against the old revision.
   Review every overlap and upstream replacement opportunity. Verify no old
   mover, sampler or generated tree replaced upstream changes.
4. Run the same complete CORE selection on the integrated candidate. Protect
   both the earlier 843 actual PASS identities and every actual upstream PASS.
   No omitted case, PASS-to-SKIP, failed retry or missing result counts as
   preservation. Registry renames need an explicit semantic crosswalk.
   The historical 919-PASS aggregate additionally included 62 Enceladus
   application cases and fourteen Moon U aliases. Re-run the Enceladus owners
   in a separate bounded package, and bind every Moon alias to its actual
   current U result through the unchanged global adapter's literal mapping.
5. Requalify Moon local/baseline gates and the affected geometry, sampling,
   source/cohort/event-budget and atomic-rate paths. Preserve independent
   oracles, sample controls and reference values. Each development campaign
   has a 3600-second ceiling including setup, build, cleanup and reporting.
6. Deliver the tested merge to the active checkout, retaining unrelated user
   edits and machine configuration. Verify hashes and report evidence.

Expected result: regression-preservation PASS means every protected PASS stays
PASS. Existing FAIL/ERROR/SKIPPED remain individually visible; it does not mean
an all-green registry. A new evaluated scientific discrepancy is FAIL, build,
parser or execution failure is ERROR, unavailable prerequisite is SKIPPED.
No source-only or supporting result closes M1/M2/M3 or grants validation.

## Provenance and scientific contracts

The source roadmap SHA-256 is
`712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17`;
the data-root roadmap is
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
The source adds approved development-order and bounded-execution amendments;
the external document retains the baseline order. The source amendments
control implementation/support execution; shared scientific qualification
gates and authoritative external metadata remain unchanged. Formal M1 and
M2 remain open. Keep every attempt, run manifest, input snapshot, executable
identity, data hash, metric and test log under the integration evidence root.

## Final upstream advancement

Upstream advanced to `5ce92cb2feecadd11a1d60ace69153a7c7da4122` during the
package. Its two additional commits change unrelated comet Python/tooling and
Earth nightly configuration, with no core/Moon/Enceladus/native-test/input/build
changes. Exact Git-object and 656-entry native-source audits prove equivalence.
Physics receipts keep their actual 4bff092 basis; the new head's 64 nightly
unit tests and official check were actually rerun. Delivery uses the final head
without rewriting old scientific provenance or assigning new qualification.
