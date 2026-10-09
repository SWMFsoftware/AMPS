# UStep12Release: strict Step 12 release-contract suite

Run from the AMPS repository root:

```bash
./srcEarth/test/UStep12Release/run_test.sh
```

The suite uses Python's standard library and does not require AMPS, MPI, Geopack,
SPICE, SWMF, or network access. It creates real temporary files and SHA-256 links so
the same production release validator is exercised end to end. These manufactured
files validate decision logic only; they are not substitutes for U/I/C/F/O physics or
observation evidence.

| ID | Reference or injected defect | Required result |
|---|---|---|
| S12-U01 | Complete fixed U/I/C/F/O registry, two campaign types, seven roles each, identical finite tables | PASS with every gate and role present |
| S12-U02 | 20% field-value change and a changed trajectory terminal outcome | Both fail; discrete outcome remains exact |
| S12-U03 | Sampled-mesh `rtol` above 5% and SWMF replay above `1e-10` | Both rejected before acceptance |
| S12-U04 | Coupled build uses a different common-physics tag | Fail even when numeric files match |
| S12-U05 | Missing gate, dry-run PASS, and F8 with zero comparison rows | Each fails closed |
| S12-U06 | O4 declares post-freeze retuning | Fail; holdout cannot be normalized after inspection |
| S12-U07 | Valid restart followed by mutation of a hash-pinned artifact | Unchanged PASS is reused; mutation prevents reuse and fails |
| S12-U08 | Placeholder command and missing supported capability | Configuration fails before numeric comparison |
| S12-U09 | CLI execution of complete, validate-only, and comparator modes | Explicit RESULT line and correct exit semantics |
| S12-S01 | Production source/tag/manifest/AUXDATA/help/doc wiring | Required Step-12 interfaces present |
| S12-S02 | Validate mutable provenance separately, require portable active paths, and hash the normalized approved `test/list` structure | Every scalar or looped `last pass:` item is empty or a full 40-hex commit; Step 12 remains an independent active gate; active commands contain no home-directory paths; C9/C10 select the checkout-local Earth executable; and the normalized structural SHA-256 matches |

S12-S02 deliberately normalizes only the optional commit value following each
active or commented `last pass:` marker. Those values are written by the test
runner and must be allowed to advance after a successful validation. The marker
itself, command text, P/F expectation, active/commented state, ordering, and all
comments remain in the structural digest. This avoids the former
self-invalidating behavior—where recording a successful Step-12 pass caused the
next Step-12 run to fail—without permitting any validation entry or scientific
gate to drift unnoticed. The portable C9/C10 path requirements additionally
prevent an unrelated AMPS application in a developer-specific checkout from
being launched under an Earth validation entry.

The evaluator enforces the existing roadmap ceilings rather than fitting the synthetic
data: exact identities, `1e-10` identical-kernel relative tolerance, 2% integrated
products, 5% mesh/detector products, and exact live/replay discrete outcomes. Tests
intentionally perturb inputs beyond those gates; no production tolerance is modified
to make this suite pass.

Complete linked release requires archived PASS evidence for U-F01–U-F10 and U-F12–13;
I-F01–I-F07 and I-F10; C1–C19; F1–F17; O1–O4; both full parity campaigns; and the
frozen O4 no-retuning result. See `../../release_validation/README.md`.
