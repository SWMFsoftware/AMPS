# IV01 — streaming plus Parker-spiral focusing

IV01 is the smallest integrated heliospheric field-line problem. An analytic
Parker spiral extends from 0.1 to 1 AU for 350 and 700 km/s solar wind. The C++
stage computes the exact spiral arclength metric, inverts that metric to sample
radius, evaluates `d ln|B|/ds`, and advances the production focused-transport
mover. Ballistic and weak `D_mumu=D0(1-mu^2)` realizations are repeated with
forward and reversed vertex storage.

The independent reference integrates the Parker characteristic at 200,000
radial midpoints. It compares physical arclength and time—not radial distance—
and uses magnetic-moment conservation for final pitch angle. Three timesteps
measure convergence. Reversed storage uses identical keyed streams and must not
change arrival or weak-scattering moments. Outputs include raw per-particle
arrival rows, exact reference/solution CSVs, JSON/JUnit, hashes, and PNG/EPS.

Run: `python3 test/run_tests.py --amps /path/to/amps --validation-case IV01
--output-dir test_output/IV01`.
