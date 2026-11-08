# SEP-common global registry adapter

`run_test.py` is the T11 global adapter for all five maintained verification
scopes below `src/models/sep_common`. It copies that source subtree into the
runner-owned work directory and invokes the local `make verify` or Python
entry points there. This avoids generated archives, reference-audit products,
and parallel/perpendicular dependency builds modifying or racing in the source
checkout.

The perpendicular companion manifest declares the model Markdown document as
a bundle member, although the checkout stores the authenticated bytes one
directory above that bundle. The adapter recreates the declared layout only
inside its disposable copy. The local verifier then checks the copied bytes
against the bundle's frozen SHA-256, so this relocation neither changes the
document nor bypasses its asset-integrity gate.

The adapter imports 143 stable global cases:

- one archive membership and PIC/MPI-independence gate;
- one complete mean-free-path aggregate gate;
- 63 parallel-diffusion gates (the coefficient generator, 30 C++ records,
  and 32 advanced records);
- two perpendicular-diffusion gates, keeping the optional high-precision
  audit as a distinct PASS or SKIP; and
- 76 pitch-angle gates (build, 12 deterministic tests, 11 independent
  reference tests, and 52 Monte-Carlo acceptance identities).

The default T11 import does not spend roughly 2.3 hours running the full
pitch-angle statistical campaign. It uses the model-local
`--skip-acceptance` mode for the deterministic/reference stages, validates
the exact 52-case registry with the acceptance binary's `--list` interface,
and reports every acceptance identity as SKIP. These SKIPs are an explicit
qualification boundary, never a PASS. A future acceptance campaign adapter
may run those same identities without changing the global IDs.

All physical values, equations, tolerances, bundles, and statistical tests
remain owned by the model-local suites. The adapter validates exact native
schemas, row order, uniqueness, command exits, and terminal markers; it does
not reproduce scientific expected values.

Run through the global harness after the suite is registered in
`src/test/manifest.json`:

```sh
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-MODEL-SEP-COMMON-001 --jobs 1
```

The suite must remain exclusive because perpendicular diffusion builds the
parallel library, and every local scope shares the isolated copied tree.
