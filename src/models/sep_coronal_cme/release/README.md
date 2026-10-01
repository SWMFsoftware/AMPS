# Frozen Stage 13 release records

`profile_template.json` selects the complete 209-case baseline shared registry;
Stage 14 is deliberately excluded from production release dependencies.
`evidence_template.json` carries no fabricated build, MPI, coronal-host,
D1--D10 or observational proof. Evaluating these templates returns INCOMPLETE:

```sh
python3 src/models/sep_coronal_cme/tools/qualify_release.py --profile src/models/sep_coronal_cme/release/profile_template.json --evidence src/models/sep_coronal_cme/release/evidence_template.json --output test_output/stage13/incomplete.json
```

A production profile must additionally register all owning application/core
mandatory cases, actual machine-report checksums, supported explicit argv for
both separately built applications, and the reference checksum. Populate and
re-freeze the profile/evidence with `preprocessing.core.freeze_record` before
qualification. Null report/reference fields and empty argv are intentional
missing-evidence placeholders, not wildcard authorities. Reports and D1--D10
products are verified against their actual bytes. Observation acceptance needs
independent datasets and frozen metrics/tolerances; campaign C also needs the
complete candidate product and N4 caps.

No production `last-pass` record is supplied because no production-qualified
release evidence was available. `--last-pass` rejects INCOMPLETE evidence.
Source verification PASS and release PASS have different meanings. See
[the Stage-13 guide](../docs/STAGE13_RELEASE_QUALIFICATION.md).
