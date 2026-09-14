# VAL02: matched focused-mover campaign

`VAL02` compares the production coefficient-driven `fte-dmumu` and event-driven
`fte-mfp` kernels under the same isotropic mean-free-path closure. It uses
`lambda_parallel=1e8 m`, particle speed `1e7 m/s`, hence `nu=v/lambda=0.1 s-1`,
4,000 particles, `dt=0.25 s`, and a 400 s duration. Seed `1502001` keys both
ensembles independently but reproducibly.

The preregistered metrics are final `mu` and `mu2`, each mover's mean-square
displacement against the persistent-random-walk diffusion limit, and detector
onset, peak time, and normalized fluence. Detector occupancy is processed with
one fixed three-cadence average before either mover is inspected. Tolerances are
stored alongside every metric in JSON/JUnit evidence.

The two algorithms have different short-time stochastic operators, so equality
of individual trajectories is neither expected nor tested. This is
cross-mover numerical closure, not proof that either mover reproduces a solar
event.
