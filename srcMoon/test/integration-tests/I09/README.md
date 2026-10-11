# I09 — surface and terminal classification

Samples the installed sphere interaction callback with seed 731993 and checks
status accounting plus outward re-emission. Domain exits are classified using
the agreed Moon-relative two-body specific energy in an inertial frame,
`|v_LSO + Omega×r|²/2 - GM_Moon/|r|`: positive is physical lunar escape;
non-positive finite-box exit is numerical-domain loss.

The two terminal fixtures start immediately inside the production mesh's
sunward `+30 R_Moon` face. The receipt records the compiled
`_PARTICLE_LEFT_THE_DOMAIN_` and `_PARTICLE_MOTION_FINISHED_` values; the
runner requires both fixtures to return the former, so a normally completed
mover step cannot be mistaken for a domain exit.
