# CV08 — absorbing-boundary first passage

CV08 advances the production Parker stochastic differential equation
`ds=u dt+sqrt(2 kappa)dW` from zero until a particle reaches the absorbing
boundary `L`. For `u>0`, the arrival time has the inverse-Gaussian density
`L exp[-(L-ut)^2/(4 kappa t)]/sqrt(4 pi kappa t^3)`. The case records every
arrival or right-censored history plus the discrete-step overshoot.

Two drift/Peclet regimes, three time steps, and ten seeds test the distribution
instead of only its mean. The independent reference implements the exact CDF.
Per-seed Kolmogorov–Smirnov decisions use alpha 0.01, high-drift moments are
checked in standard-error units, overshoot must decrease under refinement, and
a deliberately wrong `1.5*kappa` reference must fit worse. Linear crossing-time
interpolation is documented model behavior; the overshoot artifact quantifies
the remaining discrete-boundary bias.

Run `python3 test/run_tests.py --amps /path/to/amps --validation-case CV08
--output-dir test_output/CV08`.
