# CV12 — controlled particle-wave total-energy exchange

CV12 closes the loop between production wave-frame scattering and the
turbulence particle-exchange ledger. A narrow controlled particle population
scatters against alternating resonant branches. Each particle kinetic-energy
change is deposited with equal and opposite sign in the selected wave branch;
there are no external sources, advection, compression, damping, cascade, shock
injection, or escape.

Particle-only and wave-only controls must remain constant. The coupled run must
preserve `E_particle+E_wave` and report negligible per-step ledger residual for
three time steps and three population sizes. A suppressed-deposition negative
control must violate conservation, demonstrating that the gate can detect the
failure it targets. The independent reference uses only the closed-system
energy identity and does not reuse scatter kinematics.

Run `python3 test/run_tests.py --amps /path/to/amps --validation-case CV12
--output-dir test_output/CV12`. Review both the total-energy overlay and raw
branch-resolved exchange rows; JSON metrics remain authoritative.
