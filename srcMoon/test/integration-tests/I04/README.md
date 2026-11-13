# I04 — canonical Moon sampler regression

The canonical target intentionally configures `-amps-test=on`. A later normal
Moon build must remove `build/` and configure `-amps-test=off` explicitly;
changing only `-application=moon` preserves the old TESTMODE. The Moon M0,
aggregate and M2 builders/restorers now make this choice explicit. Nightly
results and normal-mode supporting evidence retain their separate build
identities; a successful compile with inherited nightly mode is not a normal
production restore.

The current-core reference is now the separately user-approved November 12
replacement at SHA-256
`ed875cefc604e4ae3580fb89b8c528146929b6b9be099b41b61697faee7f0080`.
The source samplers remain those of current AMPS core. Two same-binary runs
and one delivered-checkout four-rank run have identical output; every one of
919 upstream passing regressions is preserved. `reference/provenance.json`
records the review, condition, exact evidence and old/new values; its pinned
SHA-256 is `7bc8aff99b8ad337dd2e8cd08c4af7527d638c509ea8d4edfeb35745830f651d`.
The previous corrected-Na reference/provenance/acceptance are preserved under
`reference/superseded/*pre-current-core-09fabea5*`. The strict comparison passes
with a zero-byte diff; physical coefficients and DiffNum tolerances are
unchanged. This is a reviewed whole-file regression, not independent physical
truth, observational validation or an I13 freeze. Eight protocol self-tests
also protect comparator classification, review identity and truthful NO-SPICE
kernel-dependency metadata. Their fixtures never supply science evidence.

Earlier chemistry-only review: the November 12 Na correction changed this regression because the prior
generic processor had no Na model, while the corrected production processor
uses the already-declared Moon lifetime/shadow law. Two four-rank corrected
runs agree byte-for-byte at candidate SHA-256
`4bb811b3e05d72fa708359f1bea0344efd1a9f7e89609d7e8a3ca7eade808cba`.
The user explicitly approved that earlier replacement on 2026-11-12. Its
reference hash was the candidate hash above; the archived pre-current-core
provenance records that decision, the fixture, runs and old/new values.
The earlier provenance SHA-256 is
`2facd9d35ef2017e7d0e6c8ab7da19e1f055e70e42957f9bf85e530231f2bfb0`.
The strict checker now PASSes with a zero-byte diff; DiffNum tolerances are
unchanged. The old reference remains in `reference/superseded/`.

`prepare_reference_candidate.py` can reproduce a review-only candidate from
two run directories with the same executable and identical microscopic output,
requiring a passing independent linked chemistry/budget fixture. It records
both run paths, executable/input/reference hashes and before/after values;
it never installs a reference or modifies acceptance.json. No candidate is an
independent physical truth, an I13 freeze, or observational validation.

For subsequent core-integration reviews, supply `--reason` with the actual
documented change (for example retained upstream source-sampler changes).
The default retains the original chemistry-review reason for historical
reproduction. Reports record the generator hash and never install a candidate.

To reproduce the original review after installing the replacement, supply
`--previous-reference reference/superseded/test_Moon.ref_np4.pre-na-photochemistry`
with the full I04-directory prefix. The generator version used for the actual
approval is archived under the work package's
`reference-review/prepare_reference_candidate_used.py`; its hash matches the
version recorded in reviewed provenance. Later metadata/CLI improvements do
not alter the frozen reference bytes.

Consumes the production `make -j test_Moon` regression and requires
`test_Moon.diff` to exist and contain zero bytes. The campaign runner preserves
the build/run stream and reference identity with the result.

The superseded pre-correction reference was recovered from SWMF_data Git revision
`7a1c437369a600bd248bcf8846b66b286e7109c1`, path
`PT/AMPS/data/output/test_Moon/test_Moon.ref_np4`. The checksum below was
already frozen in this test before the recovery. Reproduce its archived copy with:

```sh
python3 -B srcMoon/test/integration-tests/I04/check_reference.py \
  --import-from-repository /home/vtenishe/SWMF_data
```

`MakefileTest/Table` invokes `test_reference.py` and `check_reference.py` for
the canonical `make -j test_Moon` check. The checker requires the exact reviewed
hash, four MPI ranks, a zero-byte DiffNum result and a `test_Moon.reference.json`
identity receipt. It also verifies the pinned reviewed provenance. Corrupt/missing reference evidence is ERROR; evaluated
numerical disagreement is FAIL. Source/parser/hash tests do not constitute a
model-regression PASS. A changed reference returns a nonzero make status rather
than being ignored by the old default recipe.

The self-test executable uses unittest output buffering. Expected FAIL/ERROR
messages from deliberately invalid reference/comparator fixtures stay hidden
when those tests pass; a genuine failing self-test replays its captured output,
prints its assertion failure, and exits nonzero. The separately invoked real
production checker remains unbuffered. On success the recipe prints its real
`PASS I04`, lists the zero-byte `test_Moon.diff`, and restores the traditional
`Test test_Moon has finished: ...` completion line. This is an output/reporting
correction only; the comparison, approved reference and thresholds are unchanged.

The checker also requires the production result file before starting DiffNum.
Only DiffNum's documented completed numerical/text-difference summaries (or a
successful comparison with nonempty diagnostics) count as an evaluated FAIL.
A missing result, comparator launch failure, or comparator failure without a
completed difference summary is ERROR. Both streams and the comparator exit
code remain in the evidence; infrastructure failures cannot count as physical
disagreement.

The current-core integration also fixes the M0 aggregate's classification of
a nonzero canonical make result. A complete strict numerical comparison with
`status=FAIL`, this campaign's exact output path, and existing executable/output
is evaluated by the normal I04 scorer and remains FAIL. Missing/build/runtime/
comparator errors stay ERROR; a stale receipt from another nightly directory
cannot qualify this run. The added protocol self-test exercises these branches
without science values. This does not approve a new reference or relax DiffNum.

This is needed because `output` resolves into the shared external SWMF_data
tree. On 2026-11-09 its four-rank reference changed to SHA-256
`9a0bd0f9dca19cba820565990de846d6acc0c1efdaf3fea97c3d4a2f050cec3a`
in commit `6d2644aa2994beaedc9fc76c97c10bba3abb99ac` for later source-model
numerics. That shared file is preserved. It does not supersede this checkout's
reviewed I04 contract. The superseded archive preserves those previously
approved bytes. The current corrected reference is the separately reviewed
two-run production regression described above.

After changing the table, regenerate the generated Makefile.test with the
existing `utility/TestScripts/table_parser.cpp` utility using a copy of
MakefileTest/Table named `Table` in its working directory. Do not patch the
generated recipe as the primary change. The root Makefile remains generated
and is not a source file to commit.

The full configuration-separated M0 runner passes `TESTMoonDIR` inside its
new campaign directory when invoking the same canonical target. This keeps
prior reviewed root `run_test_Moon` executions intact. Model configuration,
four-rank MPI layout, input and comparator remain the canonical test_Moon case.

The superseded four-rank pre-Na-correction reference has SHA-256
`abc1d9c91f61a1c44e06e98bde8739b2de1650b7a6caef783bb97081f2eb80b4`.
It was accepted on 2026-10-11 only after two complete force-aware production
runs agreed byte-for-byte. It replaces the no-force reference
`9386ea2d353db638ce84bf1f202d55cb7276443f800938c6365807be8a8deef3`;
no legacy no-force configuration is retained. This reviewed regression is the
lowest available oracle in the hierarchy and does not replace the analytical,
wiring, convergence, or conservation gates that precede later validation.
