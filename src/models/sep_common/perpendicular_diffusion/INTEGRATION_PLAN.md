# Parker-transport integration plan

This file records later host work; it is not evidence that either application
currently uses the perpendicular library.

1. Audit the actual srcSEP and srcSEP3D Parker/focused movers and establish for
   each whether the evolved quantity and stochastic process are forward or
   backward. Do not choose the SDE drift from the coefficient name alone.
2. Add a thin application parser that maps a documented input-file section to
   `BuildConfiguration`/`SetActiveConfiguration`. Preserve source locations in
   parser errors and broadcast one canonical validated configuration before
   worker evaluation.
3. Build one immutable local provider snapshot with the resolved mean field,
   explicit turbulence geometry, total component variances, tagged spectral
   scales/indices, provider revision/fingerprint, physical perpendicular axes
   when eigenvalues differ, and either a coherent parallel dependency or a
   paired NLGCE selection. Do not infer a 2D fraction.
4. Evaluate the selected model once per required local state. Diagnostics and
   field-line/statistical observables must be rejected by production Markov
   movers. NLGCE pair outputs must replace both coefficients together.
5. Assemble the Cartesian symmetric tensor and implement the full SDE5/SDE6
   divergence, including coefficient gradients and frame derivatives. Verify
   it independently against the SDE11 Cartesian stencil before optimizing.
6. Derive the correct forward/backward drift and covariance for the actual
   application measure. Avoid double counting resolved parallel streaming,
   pitch-angle scattering, field-line wandering, antisymmetric drift, or
   source injection width.
7. Add homogeneous covariance, zero-eigenvalue, unequal-frame rotation,
   analytic stationary-density, curved-field timestep, discontinuity, table,
   cache-invalidation, and restart tests.
8. Implement persistent frozen lines only in a trajectory-owned component,
   with explicit independent/shared line IDs, Brownian-bridge refinement,
   out-and-back reuse, checkpoint/restart, and timestep/line-resolution
   convergence.
9. Run each application's established clean configuration/build workflow and
   qualifying one-rank/four-rank cases. Record exact commands, input hashes,
   outputs, covariance/statistical uncertainty, and limitations separately
   from the standalone equation tests.

The integration must not change the generic PIC mover or silently enable
perpendicular transport in existing decks. An explicit model selection and
complete required input are mandatory.

