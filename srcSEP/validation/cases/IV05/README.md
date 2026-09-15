# IV05 — moving shock crossing and acceleration

IV05 translates a planar shock at 500 and 2000 km/s to a fixed field-line node.
Bracket endpoints are evaluated by the production shock-trajectory core; the
crossing is then interpolated continuously within the enclosing timestep.
Equivalent stationary- and moving-frame rows use identical keyed shock-cycle
histories, allowing direct momentum and spectrum comparison. Three timesteps
straddle the segment-per-step traversal scale, and exact node coincidence is
retained as an explicit edge-case flag.

The reference supplies the exact trajectory crossing and the CV09 strong-shock
`q=4` spectrum. Gates check crossing time, frame-invariant momentum, spectral
slope, and node handling. The controlled cycle model isolates crossing/frame
logic; a resolved time-dependent heliospheric shock remains a separate native
event validation rather than being inferred from IV05.

Run: `python3 test/run_tests.py --amps /path/to/amps --validation-case IV05
--output-dir test_output/IV05`.
