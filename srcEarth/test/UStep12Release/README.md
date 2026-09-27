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
| S12-S02 | Project current `test/list` back to Step 11 | SHA-256 equals the pre-Step-12 list, proving old commands/gates/hashes unchanged |

The evaluator enforces the existing roadmap ceilings rather than fitting the synthetic
data: exact identities, `1e-10` identical-kernel relative tolerance, 2% integrated
products, 5% mesh/detector products, and exact live/replay discrete outcomes. Tests
intentionally perturb inputs beyond those gates; no production tolerance is modified
to make this suite pass.

Complete linked release requires archived PASS evidence for U-F01–U-F10 and U-F12–13;
I-F01–I-F07 and I-F10; C1–C19; F1–F17; O1–O4; both full parity campaigns; and the
frozen O4 no-retuning result. See `../../release_validation/README.md`.
