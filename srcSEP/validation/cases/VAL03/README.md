# VAL03: independent focused-transport solver

`VAL03` transports 30,000 production `fte-dmumu` particles through the same
one-zone coefficient history used by an independently implemented conservative
finite-volume pitch-angle equation. The reference solver has its own cell
indexing, face flux, no-flux boundary, and time integration and shares no
production mover, random generator, flux routine, or boundary-reflection code.

The prescribed history is `Dmumu=0.4(1-mu2) s-1` for 0.4 s with initial density
`f(mu)=0.5(1+0.8mu)`. Seed `1503001` is used only by the production particle
ensemble. The 192-cell reference is conservatively reduced to 48 comparison
bins. Acceptance requires pitch-density L1 error at most 0.08 and absolute
errors in mean `mu` and mean `mu2` at most 0.015.

This closes a controlled identical-history comparison. It does not replace a
future comparison against a separately maintained end-to-end focused-transport
code using a time-dependent heliospheric field-line history.
