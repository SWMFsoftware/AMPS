# T11 legacy-application disposition audit

This suite keeps every legacy application named by the core-test roadmap
visible while migration proceeds. It verifies three governance facts for each
of the 13 entries: the directory still exists, a bounded global suite owns its
replacement area, and the disposition is recorded in a machine-readable
receipt. It does not execute the old diagnostic program and does not turn
directory existence into a behavioral PASS.

`mapped-not-qualified` means that a production-linked global family now owns
the subject area, but each legacy application's independent oracle,
deterministic state, failure-exit propagation, and hard-coded-path behavior
still require explicit review. `diagnostic-retained` means that the old
program remains useful for manual investigation after a global residual or
invariant test replaced its visual output. Neither disposition qualifies the
legacy program itself.

The adapter requires exact ordered equality between its 13-entry table and
the central case registry before honoring an individual `--case` selection.
Missing directories are `FAIL`; unknown global owners and control-plane
problems are `ERROR`. Infrastructure failures retain all selected case IDs in
the final report through the common result protocol.

Run the suite from the repository root:

```sh
python3 src/test/run_tests.py --test CORE-LEGACY-REGISTRY-AUDIT-001 \
  --output-dir /tmp/amps-t11-legacy --jobs 1
```

The evidence file
`artifacts/CORE-LEGACY-REGISTRY-AUDIT-001/legacy-migration.json` records the
directory, global owner, and disposition of every roadmap entry. A green
receipt closes only the T11 visibility/governance check; complete legacy
behavioral migration remains an independent coverage task.
