# AMPS SEP + corona: Stage 11 and Stage 12 source update

This archive extends the previous Stage-0--10 source package, retaining the
native `sep-corona` suite selector and the solar boundary/active corridor
validation fixes. The canonical `2026-09-28-r4` specification, application input
schema 5 and field-line bundle schema 3 are unchanged.

Stage 11 adds immutable, independently selected exterior finite-HCS and finite
moving planar sheath providers, full-vector event-split orbit stepping,
explicit relativistic shock/plasma frame transforms, restartable crossing
ledgers and global conservation audits. The qualified families and their
assumptions are documented in
`src/models/sep_coronal_cme/docs/STAGE11_DISCONTINUITY_TRANSPORT.md`.

Stage 12 adds deterministic offline parameter recovery/covariance propagation,
checksum-owned observation assets with separate construction/qualification/
withheld roles, complete preregistered formation-height candidate products,
density-conditioned radio and joint height likelihoods, reviewed ephemeris
line requests, frozen transfer protocols and matched analytic/MHD comparisons.
The detailed guide and runnable synthetic workflow are respectively
`src/models/sep_coronal_cme/docs/STAGE12_OBSERVATION_PREPROCESSING.md` and
`src/models/sep_coronal_cme/examples/stage12/README.md`.

From the AMPS source root, run every current shared-model test:

```sh
make -C src/models/sep_coronal_cme -j8 test
```

For sequential cumulative stage gates:

```sh
make -C src/models/sep_coronal_cme -j8 test-stage11
make -C src/models/sep_coronal_cme -j8 test-stage12
```

Verification on the supplied source tree:

- Clean C++17 build with `-Wall -Wextra -Wpedantic -Werror`.
- Stage 11: 202/202 canonical tests passed.
- Stage 12 aggregate: 206/206 canonical tests passed.
- Architecture/public-header and thin application adapter compilation passed.
- Deterministic model-documentation gate and synthetic preprocessing CLI
  publication/verification/selection/line-request workflow passed.
- Native boundary/corridor regression: 45 checks passed; JSON evidence passed.
- Existing application Python runner regressions: 19 tests passed.

The archive contains source, tests, examples and documentation only. Generated
objects, executables, run output and JSON/JUnit evidence are excluded. The
included `src/models/sep_coronal_cme/stages11_12.patch` describes this update
relative to the preceding source package; use either the updated files or that
patch, not both. Existing initializer/native-suite changes remain included.

The optional transport APIs and neutral preprocessing records require a host
adapter before activation in an AMPS production mover/configuration. The
manufactured planar gates do not qualify an arbitrary curved CME sheath or a
finite PFSS/SCS transition, and the offline comparison does not create an
imported-MHD runtime provider. Full linked AMPS/MPI execution and real-event
production qualification were not run because this upload is a partial source
tree. The existing corridor-outside-the-Sun contracts remain in the package.
