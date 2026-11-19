# Cut-surface mover regression and force-aware zero-force limits

The existing six-body scene, particle selection, independent brute-force
reference, statistics and acceptance criteria are executable with the
optional first-/second-order physical-time production tracers:

```sh
mpiexec -n 1 ./amps force-first nightly
mpiexec -n 1 ./amps force-second nightly
mpiexec -n 8 ./amps force-first nightly
mpiexec -n 8 ./amps force-second nightly
mpiexec -n 1 ./amps new nightly
mpiexec -n 8 ./amps new nightly
```

Configure/build the normal `input/test/cut-surface-mover.input` target in an
isolated checkout first; these are additional runtime selections. The modes
register continuation storage before particle-buffer allocation, install an
explicit zero analytical acceleration, and complete production MPI
continuations before comparing endpoints. They retain the original strict
reference tolerances, edge/vertex/concave adversarial cases and physical
surface callbacks. Other selections (`new`, `trajectory`, `simplified`,
`direct`) retain their existing layout and execution.

The outputs are the normal `cut-surface-mover.<selection>.<mode>.np<N>.report`
and `.digest`, plus stdout and generated geometry/mesh files. Preserve source,
compiler, input and executable identity with each report. A failed evaluated
criterion exits 1; a fatal runtime/build invariant is ERROR.

The `new` commands execute the original force-free mover and its existing
strict independent oracle. They verify the legacy five-argument geometry
service separately from the new tracers' opt-in velocity filter. Deterministic
classes' `.digest` files must agree across MPI layouts; the diffuse stochastic
class is explicitly excluded from this bitwise identity demand.

The optional `reuse-mesh /absolute/path/cut-surface-mover.mesh` argument may
reuse geometry construction only. Preserve the producing source/input/NASTRAN
and mesh SHA-256 values first, and require the regenerated NASTRAN to match.
All face-list, manifold, independent-reference and physical acceptance checks
still execute; no particle state or scoring/reference values are restored.

This tests the zero-force geometry/list limit, not Moon force sums or
observational validation. Accelerated/local-clock controls and portable
standalone math instructions are in `../srcForceCutSurfaceMoverTest/README.md`.
