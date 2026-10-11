# I10 — production mover convergence

Advances the shadowed SPICE fixture at 15, 30, 60, and 120 steps. The frozen
contract requires monotonically decreasing error, every refinement ratio at
least 3 (consistent with the configured second-order mover), and at most 1%
scaled state error on the finest level against an independent RK4 reference.

