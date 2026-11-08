# Field-solver and vacuum-wave qualification

`CORE-PIC-FIELD-SOLVER-001` supplies asserted numerical receipts in place of
the old diagnostic-only linear/wave programs. It solves the manufactured
Dirichlet problem `-phi''=pi^2 sin(pi x)`, checks the actual algebraic residual,
and demonstrates second-order mesh convergence against the analytic solution.

The wave cases advance a periodic 1-D Yee pair. Their oracle is the independent
discrete relation `sin(omega*dt/2)=C sin(k*dx/2)`. Energy uses the exact
staggered Yee invariant, and a separate particle-plus-field ledger checks that
an explicitly transferred joule is counted once. This bounded suite does not
qualify the legacy multi-block GMRES assembly or particle-field feedback.

Run with:

```sh
python3 src/test/run_tests.py --test CORE-PIC-FIELD-SOLVER-001
```
