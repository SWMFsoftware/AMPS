# WND3D16

a manufactured two-zone profile reproduces the inner density for `r<=r_a`, the outer velocity for `r>=r_b`, the analytic quintic `C2` `ln(rho)` blend inside the overlap, and one invariant `eta_m` everywhere. The exact quintic-Hermite interpolant reproduces every stored node value, first derivative, and second derivative; `C2` joins, every real interval extremum, positivity, and adjacent-node no-overshoot are checked. A generic spline substitution, hidden clipping, missing derivative, nonpositive state, reversed/equal join radius, incomplete overlap support, or radius support spanning nonmonotone `r(s)` without stable segment splitting fails before a background is committed.

Run from any directory with `python3 test.py`. The launcher delegates to the global registry, so individual and cumulative gates execute the identical implementation.
