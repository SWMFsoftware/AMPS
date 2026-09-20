# V05 release evidence governance

`capabilities.json` states what is coded, component-tested, production-
integrated, scientifically validated, or blocked. `profiles.json` turns those
claims into named evidence requirements. `generate_release_evidence.py` reads
the JSON reports produced by the normal test and validation runners and emits
a JSON/Markdown matrix; it exits nonzero for every missing, skipped, failed,
or blocked required item.

The `development` profile is source-only. `native-integration` additionally
requires a linked AMPS mesh/background run. `mpi-qualification` requires the
V03 rank/thread matrix. `scientific-release` requires cross-model and reviewed
observational campaigns and is intentionally blocked by deferred R8. This is
not a documentation warning: the generator records `INCOMPLETE` and exits 1.

```bash
python3 release/generate_release_evidence.py --profile development \
  --test-report test_output/all/srcsep3d-tests.json \
  --output-dir test_output/release-development
```

For a candidate scientific release, also pass the Phase-V
`validation-summary.json`. Do not edit generated matrices into the source
tree; archive the runner reports and matrix together with the executable hash.
