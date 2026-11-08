# Configured dust charging test application

`srcModelDustTest` is the native integration boundary for the legacy AMPS dust
charging model. It is selected by `input/test/model-dust.input`; it is not a
general dust injection, mover, surface, or split/merge qualification.

The registered global owner is `CORE-DUST-NATIVE-001`. Run it from the
repository root:

```sh
python3 src/test/run_tests.py --test CORE-DUST-NATIVE-001 \
  --output-dir /tmp/amps-dust-native
```

The adapter is repository-exclusive because AMPS selects one generated
application for the checkout. It verifies the root and active-process state,
removes only `build`, invokes `Config.pl -application=test/model-dust`, and
compiles through top-level `make amps` with `MAKEFLAGS=-j16`. It never edits
`.amps.conf` directly. Machine commands, rank flags and temporary/output roots
come from the layered `src/test` configuration.

## Physical and numerical contract

The application uses the production `UpdateGrainCharge` implementation at a
heliocentric distance of 2 AU. The fixed test plasma has ion temperature
`1e6 K`, number density `1e7 m^-3`, bulk velocity `(4e5,0,0) m/s`, and electron
pressure `1.5 n k_B T`. Four compact spherical grains have radii `1e-7`,
`1e-6`, `1e-5`, and `1e-4 m` and start uncharged.

For every radius, the case first solves the steady-current root and then
advances the implicit time-dependent update in `1 s` steps. It requires:

- finite positive equilibrium charge;
- finite transient iterates;
- convergence to the steady root within relative error `1e-5`;
- convergence before the explicit 200,000-iteration safety cap; and
- radius-independent `q/a` within `1e-5`.

The last criterion is an independent dimensional oracle. Each retained
collection/emission current scales with grain area `a^2` and depends on charge
through the isolated-sphere potential `q/(4*pi*epsilon0*a)`, so identical
plasma conditions imply constant equilibrium `q/a`. The tolerance matches the
production Newton solver's declared `1e-5` relative charge-increment stopping
scale; it is not a claim of roundoff-level root accuracy.

## Evidence and output

The adapter retains the clean build transcript, executable and input-deck
SHA-256 identities, MPI transcript, and full-precision numerical charging
history. It sets `AMPS_TEST_MODEL_DUST_OUTPUT` to a private artifact path, so a
nightly run cannot overwrite `PT/plots/test_model-dust.dat`. Direct legacy runs
without that variable retain the historical output location.

The application prints stable case records:

```text
[CASE] RUN  CORE-DUST-CHARGING-DYNAMICS-001 ...
[CORE-DUST-CHARGING-DYNAMICS-001] PASS
```

## Limitations

Passing this case proves configured compilation and bounded consistency of the
steady and transient production charging solvers. The case also evaluates all
four enabled signed current components at -5, 0, and +5 V against immutable
values calculated independently from equations 4 and 10--15 of Horanyi
(1996), *Annual Review of Astronomy and Astrophysics* 34, 383--418,
DOI `10.1146/annurev.astro.34.1.383`.  Its equilibrium potential is compared
with the independently bisected published-equation root
`13.64710731500525 V`; the tolerance is the production nonlinear solver's
declared `1e-5` relative stopping scale.  The exact plasma, grain, emission
parameters, signed current values, and SI units are recorded next to the
fixture in `main.cpp` and in the private numerical receipt.

The configured charging case does not exercise stochastic dust injection,
AMR ownership, MPI migration, restart, or the dust-specific size-group
split/merge operation. Those capabilities have separate T7 cases and must not
be inferred from this charging receipt.
