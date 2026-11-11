# Native ECSIM manufactured-operator tests

This exclusive suite clean-builds `test/stencil` and validates the production
ECSIM curl, grad-div, Laplacian, and curl-curl stencil builders on a periodic
unit cube.  The manufactured Fourier modes are integer multiples of `2*pi`,
so values and all derivatives are periodic.  Each operator uses five meshes
and asserts the observed final-grid order independently of its diagnostic
table: second-, fourth-, sixth-, and eighth-order schemes must achieve at
least `p-0.2`.  Curl-curl additionally asserts the discrete
`curl curl(E) = grad div(E) - laplacian(E)` residual.

Run through the portable global orchestrator:

```sh
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-ECSIM-OPERATORS-001 \
  --output-dir /tmp/amps-core-t8-operators
```

The adapter follows the native preflight, retains the executable and generated
identity, and accepts only an exact case-qualified PASS record.  These tests
cover periodic interior operators.  Physical-boundary rows, deposition,
solvers, and wave propagation have separate T8 owners.

This Cartesian fixture also requires a global no-SPICE build, independent of
the previously configured application. Its adapter selects
`Config.pl -application=test/stencil -spice-path=nospice`, verifies the saved
selection and effective make flags, and checks actual C++ build commands before
MPI. Configuration mismatches report ERROR for every requested case, with
`ecsim-operators-configuration.json`/`.txt` retained beside the build identity.
The manufactured Fourier fields and order/identity tolerances are unchanged.
