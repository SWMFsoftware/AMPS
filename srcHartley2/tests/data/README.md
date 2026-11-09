# Source-controlled generic-comet test data

This directory is the complete data dependency of the portable test
orchestrator. `SHA256SUMS` is generated over every other file below this
directory and is verified before archived evidence is parsed.

After an intentional, reviewed fixture update, regenerate the manifest from
this directory with paths sorted in the C locale:

```sh
find . -type f ! -path ./SHA256SUMS -print0 | LC_ALL=C sort -z | \
  xargs -0 sha256sum > SHA256SUMS
sha256sum -c SHA256SUMS
```

## `archived/`

The directory names mirror the original 2026-11-04/05 validation archives so
the same strict checkers can audit either source-controlled fixtures or full
external archives.

- Stage 0 retains the final 3,010-row C2 table, Rubin comparison table, compact
  convergence records, every Stage-0 diagnostic JSON/CSV pair, and exact cell
  measures. `portable_baseline.json` authenticates the reduced final-window
  fixture and preserves its accepted numerical expectations.
- Stage 1 retains the three production policy/cadence logs. Its legacy file
  equivalence fixture is intentionally tiny and byte-identical in both arms;
  it tests the comparison contract without committing two 86-MiB plot trees.
- Stages 2--9 retain the exact validator logs and Stage-7 diagnostic sidecars
  consumed by their strict parsers.
- Stage 4/6/7/8 provenance directories retain generated inputs, definitions,
  and policy headers. Binary executable hashes belong to the complete external
  archives and are intentionally excluded from this portable code test.
- The Stage-9 checkpoint pair is the compact two-species/two-cell valid fixture
  used by the mutation suite. It exercises the complete JSON/summary/cell
  reconciliation contract without committing the 177-MiB smoke archive.

These records are regression fixtures, not evidence that the current source
has completed a new MPI calculation or achieved scientific convergence.

## `i3/Archive.zip`

This is the original supplied 3I/ATLAS archive containing the six-case PDF and
twelve H/OH MAVEN/IUVS raster panels. It has SHA-256
`344d33e69c913a5e76d57ddb7c92c1a8d1e52a92a03c06f2f5849f70e45597be`.
The portable suite regenerates the committed numerical observation table and
metadata byte for byte from this archive. No network access is required.

Do not replace or hand-edit a fixture merely to make a failing checker pass.
A fixture update must identify the originating accepted run, retain the full
external archive, update the documented reduction policy, regenerate
`SHA256SUMS`, and undergo scientific review when numerical expectations change.

