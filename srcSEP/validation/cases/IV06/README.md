# IV06 — coupled self-generated turbulence feedback

IV06 implements the plan's three-control experiment: frozen turbulence,
one-way streaming-driven growth, and two-way growth plus production pitch-angle
scattering. A narrow beam drives one registered resonant bin. Increased wave
energy raises `D_mumu`; the resulting isotropization reduces streaming and
self-limits subsequent growth.

The early one-way rate is checked against the CV11 exponential law, while the
coupled timeline compares streaming with the one-way control. Every requested
growth increment is submitted through the production turbulence core as a
pending particle-to-wave exchange. A closed reservoir removes the amount that
the production ledger reports as applied (rather than assuming that the full
request survived its positivity policy), so both the wave-only ledger residual
and the combined particle-plus-wave energy invariant are saved for audit.
Gates also require unchanged frozen waves and exact resonant-bin localization.
The PNG/EPS timeline shows the causal reduction; raw rows retain anisotropy,
streaming, `D_mumu`, wave energy, combined closure, and wave-ledger closure.

Run: `python3 test/run_tests.py --amps /path/to/amps --validation-case IV06
--output-dir test_output/IV06`.
