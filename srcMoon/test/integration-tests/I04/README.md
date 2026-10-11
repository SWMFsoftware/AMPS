# I04 — canonical Moon sampler regression

Consumes the production `make -j test_Moon` regression and requires
`test_Moon.diff` to exist and contain zero bytes. The campaign runner preserves
the build/run stream and reference identity with the result.

The current four-rank reviewed reference has SHA-256
`abc1d9c91f61a1c44e06e98bde8739b2de1650b7a6caef783bb97081f2eb80b4`.
It was accepted on 2026-10-11 only after two complete force-aware production
runs agreed byte-for-byte. It replaces the no-force reference
`9386ea2d353db638ce84bf1f202d55cb7276443f800938c6365807be8a8deef3`;
no legacy no-force configuration is retained. This reviewed regression is the
lowest available oracle in the hierarchy and does not replace the analytical,
wiring, convergence, or conservation gates that precede later validation.
