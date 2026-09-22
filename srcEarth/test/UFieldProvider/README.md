# UFieldProvider — field-provider and immutable-snapshot contract

This dependency-free unit test validates roadmap Step 3 without requiring AMPS, MPI,
Geopack, Tsyganenko libraries, or SWMF. It compiles `util/FieldProvider.h` against a
small uniform-field provider and checks the semantics required by both production
backends. It also compiles the production Boris mover, dipole evaluator, boundary
finder, and trap detector for a numerical reproduction of the F4 grid.

Run it from any directory:

```bash
./srcEarth/test/UFieldProvider/run_test.sh
```

The contract checks cover:

1. deterministic snapshot IDs that change with field-defining state;
2. complete epoch, coordinate-frame, interpolation, validity, and SI-unit metadata;
3. valid in-domain sampling with snapshot identity propagation;
4. explicit `STALE_EPOCH` failure;
5. distinct `OUTSIDE_DOMAIN` and `INVALID_REQUEST` failures;
6. immutability of an existing snapshot after the provider advances;
7. the same-snapshot gate used between cutoff and flux/spectrum products; and
8. rejection of false unit or electric-field capability claims.
9. production analytic-DIPOLE equatorial and polar point values in GSM/SI units; and
10. an owned analytic-DIPOLE parameter value remains unchanged after legacy global
    reconfiguration and under 16 concurrent readers.

The second executable then traces the exact five F4 latitudes, 32 logarithmic energies,
and 16 deterministic directions in the analytic dipole. It uses the corrected F4
policy (`BORIS`, 1-Re bounce-envelope tolerance) and fails if any point/energy bin has
zero physically resolved directions. This is a numerical regression for the original
`NaN` failure, not only a source or metadata check.

Why the F4 policy is explicit: RK4 accumulated about 0.3--0.8% momentum-magnitude
drift on long low-energy traces while the trap detector correctly required relative
energy stability of `1e-4`. Consequently, no trajectory in some bins could receive a
positive recurrent-trapping classification. Boris preserves the static-magnetic
energy invariant to roundoff. The 1-Re envelope tolerance accommodates the sparse
16-direction closure sample; time, step, and distance limits remain unresolved.

The suite tests the common API, the analytic snapshot's concurrency semantics, and the
F4 setup. Backend-linked validation is still required in a configured AMPS tree:

- point comparisons for IGRF and selected Tsyganenko models (DIPOLE is covered here);
- compact Mode3D versus source-evaluator values;
- PIC::CPLR owner-cell buffer versus compact SWMF snapshot values; and
- cutoff plus density/flux execution proving a single snapshot ID is retained.
