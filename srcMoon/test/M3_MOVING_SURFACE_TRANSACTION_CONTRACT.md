# M3 moving-surface preview/commit and frozen-body transform contract

This extension is recorded before the Moon-side production change. It fixes
two demonstrated wiring inconsistencies without changing shared PIC or a
surface thermal/sticking law. It introduces no observation or permanent ID.

The new force tracers preview geometry repeatedly during adaptive quadrature.
`FindProductionMovingSurfaceIntersection` previously published its most
recent candidate frame and incremented the collision diagnostic. Consequently
a rejected trial/second-half candidate could overwrite the frame consumed by
the accepted impact callback. In new-tracer mode the finder must be pure:
physical collision count and committed frame remain unchanged on any preview.
The accepted interaction reconstructs the body-frame normal and both state
transforms from its actual selected LSO normal and accepted physical epoch,
then increments the accepted count once. The old legacy moving M2 path keeps
its existing finder/receipt semantics and values.

The approved new frozen profile holds body pose at SOURCE-BEGIN, with zero
material velocity relative to LSO. Body-to-LSO and inverse six-dimensional
state maps retain their rotation/velocity-rotation blocks and zero only their
lower-left derivative block. A shared production policy applies at initial
global frame setup, normal frame update, LOLA source mapping and accepted
collision mapping. Physical birth/flight clocks retain their nonzero source
offsets; an epoch is never reset merely to freeze environmental geometry.
Legacy transforms are unchanged. Thermal laws/parameters are unchanged;
body coordinates and SOURCE-BEGIN thermal caches must use the same pose.

U09 extensions use analytical closed octahedral geometry and actual furnished
SPICE transforms, never observations. Controls preview two distinct facets,
verify count/frame rollback, commit only the first selected normal, verify one
accepted collision and its independent body-frame normal/inverse transforms,
and compare nonzero-birth-offset frozen source mapping with the exact
SOURCE-BEGIN transform. A pure matrix control verifies zero derivative and
unchanged rotation for new frozen mode and exact legacy preservation. All
SPICE files actually furnished are inventoried. Missing SPICE is SKIPPED;
build/parser/runtime faults are ERROR; evaluated disagreement is FAIL. These
are local production-path checks, not I14 or M3 qualification.
