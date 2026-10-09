# Mean-free-path library maintenance instructions

- Treat `MEAN_FREE_PATH_MODEL.md` version 1.2 and its companion bundle as the
  scientific contract. Do not repair a printed equation or fill a missing
  parameter by analogy.
- Keep this library independent of PIC, MPI, srcSEP and srcSEP3D headers.
- All application input adapters must call the parser-neutral manager; do not
  assign `ActiveModelFunction` directly.
- Keep lambda kinds, momentum variables, operator conventions, variance
  conventions, field normalizations, published variants and runtime states
  explicit. Missing inputs are not zero and UNIT normalizations are not VALUEs.
- Preserve source-incomplete identities as typed gates. Never substitute a
  runnable model to make a test pass.
- Add independent equation fixtures and error-path tests for every new backend.
  Keep algebraic, operator, host-coupling and science-validation evidence
  separate.
- Update `README.md`, `DEVELOPMENT_LEDGER.md` and `INTEGRATION_PLAN.md` when an
  implementation or host boundary changes.
