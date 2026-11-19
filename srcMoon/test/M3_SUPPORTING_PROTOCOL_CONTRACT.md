# M3 supporting-result and resource protocol

This supporting infrastructure extends the existing U16 orchestration controls
and U15/U17 child-result admission; it creates no new permanent U/I ID and
evaluates no lunar observations or production physics. The contract was written
before changing the supporting supervisor and global nested-result loader.

The actual entry points are `m3_supporting_services.run_services` and
`run_all_test.load_runner_result`. Their independent references are explicit
protocol/status identities and OS process lifetime, not model-generated truth.
There are no scientific units, frames, epochs, parameters or calibration cuts.
Small temporary metadata/child-process fixtures are labelled infrastructure
controls and are never substituted for qualified native tests.

All supporting groups share one supplied wall deadline, including compilation,
native children, metadata and cleanup. One existing `CampaignBudget` instance
owns every group's full descendant tree. Its finite termination allowance is
part of that same interval. Resource histories/results are retained separately
per group; the complete output tree counts toward artifact/disk limits. A timeout,
signal, leaked descendant or deadline exhausted before launch is ERROR. An
unexecuted remaining group still gets its exact runner/intended-result path and
a real diagnostic log. No later group receives a fresh wall budget.

Each child must return one of PASS/FAIL/ERROR/SKIPPED and the corresponding
0/1/2/77 exit. Normal child FAIL remains FAIL; missing prerequisites remain
SKIPPED. A status/exit contradiction or missing/corrupt child receipt is ERROR.
No directory, successful compiler or zero process exit creates a scientific PASS.

When a local result declares `supporting_checks`, each row must have a unique
nonempty check identity, valid status and nonempty runner/log strings. Exactly
four integer aggregate counters, including zeros, must reconcile with those
actual rows. Booleans are not integer counts. A parent PASS cannot contain a
non-PASS child; a parent SKIPPED cannot hide FAIL/ERROR; a parent FAIL cannot hide
ERROR. Parent ERROR can retain any child outcomes. Formal SKIPPED with local
PASS or unavailable supporting components remains SKIPPED. Contradictions become
global ERROR while preserving the actual child rows, runner/log paths and original
reported status/counts for diagnosis. Failure printing exposes nested FAIL/ERROR.

The bounded protocol tests exercise the real supervisors with a temporary child
and grandchild, plus typed metadata controls for PASS/FAIL/ERROR/SKIPPED, incorrect
counts, malformed rows, duplicate names and exit mismatches. The process fixture
must leave no owned descendants alive. No AMPS configure/build/MPI job or canonical
data access is needed. All files stay under the explicit fresh evidence directory.
Expected PASS means these infrastructure criteria were evaluated and met; it does
not close U15/U17 formal physics, I21–I24 or any prior milestone.
