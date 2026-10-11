# I08 — ion Lorentz mover

Uses only the dedicated no-SPICE `Na,NA_PLUS` fixture. Four time-step levels
advance Na+ through the production Moon acceleration callback and force mover
with configured uniform electric/magnetic fields. Convergence is measured
without enabling ions in the normal neutral-Na science input.

The frozen levels are 64, 128, 256, and 512 steps per gyroperiod. An earlier
32–256 diagnostic sequence had the expected factor-of-four error decrease but
ended at 1.705% scaled error, above the unchanged 1% gate. The finer sequence
was therefore selected before an M0 pass was accepted; the tolerance was not
weakened.

The run also exercises midpoint-node wiring in
`UniformWeight_UniformTimeStep_SecondOrder`: the force callback receives the
AMR leaf containing the predicted midpoint. A missing/outside midpoint must
return the production domain-loss status before the callback is invoked; a
null node, a final-position node, or a test-only force path is not accepted.
