# CV06 — Legendre-mode pitch-angle diffusion

CV06 isolates the production focused-transport `Dmumu` update in a homogeneous
segment with no spatial, focusing, or momentum evolution. For
`D_mumu=D0(1-mu^2)`, every Legendre polynomial is an eigenfunction and
`a_l(t)=a_l(0) exp[-l(l+1)D0 t]`. Modes 1–6 are initialized separately from
`f(mu)=(1+epsilon P_l(mu))/2`; all modes 0–6 are measured so cross-mode leakage
cannot hide behind a correct target coefficient.

The linked C++ stage uses the production `AdvanceFocusedTransportDmumu` mover,
its physical `mu=+-1` boundary treatment, keyed random streams, three time
steps, and three independent seeds. The case-local Python reference evaluates
the eigenvalue formula without importing production code. Acceptance covers
decay rate, leakage in standard-error units, isotropic normalization, weak
timestep behavior, and boundary diagnostics. At this practical ensemble size,
the refinement-order threshold is a regression sentinel rather than a formal
deterministic order proof; increase `particle_count` for publication studies.

Run with `python3 test/run_tests.py --amps /path/to/amps --validation-case CV06
--output-dir test_output/CV06`. The directory contains raw modal samples,
independent reference and solution CSVs, native/aggregate reports, hashes,
logs, and PNG/EPS comparison figures.
