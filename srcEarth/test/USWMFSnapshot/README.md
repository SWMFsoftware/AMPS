# USWMFSnapshot: Roadmap Step 9 validation

Run from the AMPS repository root:

```bash
./srcEarth/test/USWMFSnapshot/run_test.sh
```

The suite has no MPI, PIC, SPICE, Geopack, or SWMF dependency. It compiles the exact
portable contract used by the production live/export/replay paths with C++11 warnings
treated as errors, then audits the production wiring.

| Gate | What is compared or rejected |
| --- | --- |
| S9-U01 | Mesh revision, content fingerprint, and snapshot ID against fixed known answers |
| S9-U02 | `E=-u x B` analytic SI vector; magnetic-only default and explicit experimental-E identity |
| S9-U03 | Exact cell round trip plus byte-identical write/read/write serialization |
| S9-U04 | Reversed/interleaved records; decomposition-independent mesh and state identities |
| S9-U05 | Restart identity; B changes state only, geometry changes mesh revision |
| S9-U06 | Corrupt identity/revision, wrong units/mode/time, duplicate/missing/NaN cells, and truncation fail |
| S9-U07 | Unit-specific position/B/u/E live/replay tolerances, including required passing and failing field/mesh perturbations |
| S9-U08 | Unavailable/double-freeze/stale/corrupt states fail, including corruption queued behind a frozen batch; status is fail-closed |
| S9-S01 | Time/offset/domain/control coherence, direct-owner parity, collective I/O, batch freeze, export, strict replay, and magnetic-only default are wired into production |
| S9-S02 | Representative F4, C9, and C19 gate strings remain unchanged |

These are contract and numerical-reference tests, not substitutes for coupled runs.
The complete Step 9 validation sequence additionally requires replaying an exported
SWMF snapshot through standalone Mode3D, comparing its products with the live coupled
products, repeating the same epoch after restart, and running the applicable C-tests
with 1x1, 2x8, and the default 8x16 MPI/thread layouts. Existing acceptance thresholds
must be used unchanged.

The portable suite does not claim live/replay product agreement without a linked SWMF
run. Retain the live CSV and status JSON, replay the CSV with
`examples/standalone_step9_swmf_replay.in.template`, and compare all cutoff, access,
flux, and spectrum artifacts at the plan's existing I-F02 and I-F04–I-F07 thresholds.
A missing product or non-`PASS` status is a failure, not an excluded sample.
