# UFieldProvider — strict Step 3 reference tests

Run from the repository root:

```bash
./srcEarth/test/UFieldProvider/run_test.sh
```

It is also the second independent entry in `srcEarth/test/list`, immediately after the
Step 2 flux-numerics suite, and therefore runs under the standard test runner. Its own
`last pass:` field is updated only after a successful committed validation run.

This dependency-free C++11 suite compiles the production `FieldProvider.h` contract
and the production analytic-dipole implementation with `-Wall -Wextra -Werror`. It
does not require an AMPS executable, MPI, PIC, SPICE, Geopack, or SWMF.

The suite is intentionally not a smoke test. Its ten gates include fixed or independent
references:

| ID | Strict validation reference |
|---|---|
| U-P01 | hard-coded FNV-1a snapshot ID; identity changes for physical state but not cutoff/flux request labels |
| U-P02 | complete schema, GSM frame, interpolation, immutability, and exact SI-unit contract |
| U-P03 | exact uniform B/E component values and sample identity |
| U-P04 | explicit stale-epoch rejection |
| U-P05 | inclusive domain boundary plus distinct outside-domain and non-finite-query failures |
| U-P06 | old field values and identity remain unchanged after provider revision |
| U-P07 | exact cutoff/flux synchronization; replacement and reused-ID metadata tampering rejection |
| U-P08 | negative tests for false schema, nT units, unknown interpolation, and unavailable E |
| U-P09 | closed-form aligned equator/pole and tilted-axis dipole values |
| U-P10 | immutable dipole values after legacy-global reconfiguration and under 16 concurrent readers |

The pre-existing C/F validation cases and their thresholds are not changed by Step 3.
In particular, F4 retains its original mover, trap detector settings, numerical limits,
reference reconstruction, and pass/fail tolerances. The new unit suite supplements
those cases; it does not replace or bypass them.

Backend-linked validation still requires a configured AMPS tree:

- existing C/F cases for cutoff, flux, spectra, and external references;
- direct-versus-compact point comparisons for standalone Mode3D;
- owner-cell versus compact-array comparisons for SWMF B/u and the diagnostic
  `E=-u×B` convention (E is exposed only in explicitly experimental mode); and
- a combined cutoff+density run confirming that one snapshot ID survives both products.
