# Native test fixtures recovered during Moon/core integration

The pinned upstream revision `09fabea55e367f07a530460b488fa7a981ff6176`
registers native fixtures whose input decks were excluded by the general
`*.input` ignore rule. The `data` rule also excluded the Enceladus standalone
shape fixture. A fresh checkout therefore fails infrastructure/traceability
checks before it can reproduce the documented passing native regressions.
The reduced CME registry also requires the existing root
`CODEX_REDUCED_SHOCK_TASK.txt` task/specification document, absent from the
upstream commit. Its exact bytes are preserved with the same provenance;
it is test contract material, not an observation or a new physical closure.

The unchanged existing files were preserved from the read-only local checkout
`/home/vtenishe/Mars2/AMPS`, HEAD
`d068476edbfdecc5e2a9a23408e199898f35fb5a`, with individual SHA-256 values in
`native_fixture_provenance.json`. That source HEAD identifies the checkout,
not a claim that the ignored files were present in its commit: each record
explicitly states whether Git tracked the source file. No physical input was
invented, fitted or altered. These are numerical verification/source-parser
fixtures, never observations or validation evidence.

Both upstream and integrated core runs use the same preserved fixture bytes,
machine configuration, compiler/MPI toolchain and unchanged acceptance tests.
The initial incomplete-checkout attempt is retained separately. Input decks
exercise the registered periodic/local/global time and weight modes, chemistry
daughter/reaction cases, rotational relaxation, VHS/VSS collisions, 1/2/3-D
mesh geometry and spherical exosphere/dust sources. Their normal native test
owners still perform clean Config.pl builds and preserve generated mode,
executable and input identities. Missing capabilities/data remain typed skips;
adding a deck alone cannot make its numerical criterion pass.

The narrow `.gitignore` exceptions retain these source fixtures in future
change lists without changing the general policy for generated decks or data.
The initial Enceladus exporter checks also exposed three absolute host-path
comments in each of two generated source receipts. Their retained baseline
failures are infrastructure evidence. The separate preparation script
`srcMoon/test/prepare_core_regression_inputs.py` invokes the existing exporter
and permits changes only to `! manifest=`, `! geometry_config=` and
`! stripe_geometry=`. It backs up both receipts and verifies byte-for-byte
identity of their complete non-path payloads before replacing either file.
Hashes, source strengths, geometry and acceptance criteria are preserved. Both
upstream and candidate application campaigns receive the same preparation;
this is checkout relocation, not modification of a scientific reference.
