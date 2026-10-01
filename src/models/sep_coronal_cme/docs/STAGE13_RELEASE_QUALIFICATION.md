# Stage 13 release evidence and qualification

The source/package, cross-application invocation and report machinery is
implemented in `tools/release/qualification.py`. `REL3D01--03` are selected by
the shared full suite. They exercise adversarial manifests, missing products,
failed/skipped cases, independent explicit executable paths, stale reports,
archive checksums and generated-artifact rejection. These are protocol
verification, not a production release or an observation comparison.

```sh
make -C src/models/sep_coronal_cme test-stage13
python3 src/models/sep_coronal_cme/tools/qualify_release.py --profile PROFILE.json --evidence EVIDENCE.json --output release-evidence.json
```

The profile is a frozen `sccm-release-profile-v1` record, with software stage
13, mandatory test IDs, optional preregistered event/campaign-C requirements,
and `stage14_is_release_dependency=false`. The evidence record retains all
statuses and must establish clean builds of sep_common, sep_coronal_cme,
srcSEP3D and srcSEP; all D1--D10 with finite/provenance identities; integration,
MPI, convergence and cross-model evidence; domain of validity and limitations;
and qualified production coronal runtime state. Event profiles additionally
require frozen datasets, metrics/tolerances and complete candidate/N4 evidence.
Missing evidence is INCOMPLETE. It cannot update last-pass records.

For a cross-application invocation, additionally supply
`--srcsep-executable`, `--srcsep3d-executable`, `--bundle` and a fresh
`--application-output` directory. The profile owns explicit argv lists for
both applications with `{executable}`, `{bundle}`, `{report}` substitutions.
There is no shell expansion, sibling-source search or default executable.
The receiving application report identifies its application, bundle checksum,
status and evidence kind. This boundary is additive; it does not replace the
existing srcSEP test runner, which is absent from this partial source upload.

The supplied upload has the srcSEP field-line adapters, but lacks its full
application/test tree, independently built executables, complete production
coronal adapter and preregistered observational diagnostics. Consequently no
production Stage-13 qualification is asserted. Existing generic native
Parker/SWCME initialization passes remain a distinct evidence class. The
remaining application/core registrations and real campaign/MPI runs are
required before the specification's Stage-13 exit gate can close.

`source_manifest` freezes required source, build-source, example and owning
registry paths and hashes. `package_sources` copies only declared safe regular
source files and verifies archive members and hashes. `verify_sources` rejects
source drift. `update_last_pass` accepts only qualified evidence, without
automatically committing anything to git.

The delivered overlay embeds `RELEASE_SOURCE_MANIFEST.json`. It includes
all payload source hashes and file modes, owning test implementations, build
source paths and examples. `verify_archive(..., embedded_manifest_name=
"RELEASE_SOURCE_MANIFEST.json")` checks that metadata separately from its
payload to avoid a self-referential hash; extra/duplicate members fail.
The [template profile/evidence](../release/README.md) intentionally evaluates
to INCOMPLETE. A freshly requested cross-application failure also blocks the
CLI result and last-pass update, even if an older frozen evidence file passed.
