# Shared-core background-only coronal-CME event assembly

This implementation guide maps the maintained assembly to the canonical
Sections 8.4--8.6, 14, 15.4 and 17.2 of [model.md](../model.md). It does not
replace their physics, configuration or acceptance authority. The standalone
SEP integration roadmap numbers this work Stage 14; it is distinct from this
library's earlier Stage-14 research extensions and the broader `CMLU` launch
qualification roadmap.

## Public boundaries and ownership

`include/sep_coronal_cme/background_event.h` separates three operations:

1. `MaterializeBackgroundEvent` resolves a normalized, fingerprint-checked
   configuration and explicitly acquired geometry assets into typed physical
   histories. The caller supplies support, numerical controls and the actual
   ambient identity. Library code does not open files or discover sibling
   application inputs.
2. `BackgroundEvent::Create` validates the entire requested history and seals
   the physical identity from all typed physical fields. Retaining a caller's
   previous ID after changing a history cannot alias the old physical event.
3. `Prepare(time, ambientEpoch)` prepares a complete immutable epoch privately,
   then publishes front, piston, ambient ownership, support flags and local RH
   classification together. Old epoch handles remain valid. An acquisition,
   coverage, EOS, identity, topology or RH error leaves the previous parent and
   child shock generations untouched.

The ambient interface belongs to the dependency-light coronal core. The
reusable implementation around the actual canonical background belongs to
`sep_corona_swcme/coronal_event_ambient.h`; this maintains the outward
composite-to-coronal dependency direction. `BindAmbientEpoch` retains the
immutable `AmbientModel`, exact expected owner/time/generation/epoch and EOS.
It supplies actual PFSS open/closed and Parker states, signed sector, ion and
electron data, full SI tensors and numerical-stencil metadata. It does not
manufacture an ambient state for a missing provider region.

## Histories, support and attachment

All four independent values are explicit: radial center distance, radial
semiaxis and two lateral semiaxes. Each smooth law carries its own initial and
final rates, transition start/duration and reference value. The canonical law
uses `q_ref + J(t) - J(t_ref)` so the reference can occur before, within or
after acceleration. Existing pre-transition-reference arithmetic is retained.

A table supplies every value and rate. The scalar range calculation includes
interior derivative roots, checks strict time ordering and refuses temporal
extrapolation. The monotone-C1 asset admission checks derivative extrema of
each polynomial and never clips supplied slopes. Table range evaluation scans
the validated table linearly, including at the declared 100000-sample limit.
The local canonical two-knot evaluator supplies each polynomial extremum.

The `RequireSolarIntersection` policy certifies an attached history using
exact spatial radial extrema and bounded temporal interval certificates.
`AllowDetachedFront` explicitly permits the same prescribed surface to leave
the Sun. The assembly does not slow a center, change an axis or extend a time
step to force attachment. The documented short-acceleration stress case
crosses at 211.254 s and therefore rejects an attached full-run declaration.

An enabled piston has its own complete history. The conservative all-history
certificate maps the piston into normalized front coordinates and bounds its
maximum norm. This is a sufficient containment/separation certificate. Some
physically nested configurations may remain unresolved within its declared
work budget; those configurations reject instead of being guessed valid.
An independent piston is geometry only here: no contact state, ejecta field,
fluid sheath or regional material map is inferred.

## Exposed cells and complete ambient coverage

`TessellateSolarExposed` integrates actual exposed parameter-cell area using
solar cuts and azimuth topology events. Its representative may move away from
the original midpoint but retains that cell's identity. This prevents loss of
small exposed caps or off-midpoint triaxial islands. The legacy midpoint
`FixedOrientationEllipsoid::Tessellate` interface is preserved for existing
users; Stage 14 uses the new cut-cell interface for both front and piston.

The solar absorber, plasma radial support and source-qualified radius have
separate status. Strict plasma admission uses the full continuous exposed
radial extent, not its quadrature representatives. An explicitly requested
geometry-only diagnostic policy may retain unsupported surface geometry with
missing-plasma flags; it does not invent a jump there. A conservative radial
bound for each complete parameter cell uses full polar and azimuthal widths
because the exposed representative may lie anywhere within the cut cell.
Per-point support and whole-cell coverage are separate records, and an
unresolved boundary crossing remains explicitly unresolved.

The cut-cell quadrature does not certify unresolved HCS/separatrix crossings
or provide native particle transport across them. A physical PFSS/Parker
radial trace is also distinct from a branch-fixed Cartesian derivative
stencil: requesting the wrong radial side rejects. At the exact source
surface the actual canonical equality policy is Parker; no opposite trace is
synthesized by relabeling the same value.

## Local shock state and transactions

Supported points use the maintained characteristic speeds and oblique MHD RH
solver. Sub-fast geometry has no jump. Actual fast jumps must satisfy the
original entropy/compression checks and normalized residual tolerance `1e-9`.
The background-only shock inputs disable source eligibility and number/energy
rates. No source plan, particle record or random draw is needed to prepare the
background epoch.

Preparation copies the ambient descriptor, retains its owning handle and
rechecks the full descriptor before publication. Local shocks are prepared on
a private copy of the canonical provider; committing a successful private
child cannot expose a new generation before the parent is ready. Backward
publication time and regressed ambient generation are errors.

The maintained manufactured high-beta negative deliberately preserves the
canonical solver's fail-closed response when a tiny tangential magnetic stress
cannot meet its relative residual bound after Cartesian reconstruction. The
positive oblique cases use resolved fields and preserve the same tolerance;
there is no gas-dynamic fallback or gate relaxation.

## Tests and practical scope

The owning `S14BG01`--`S14BG10` records are registered in the local and shared
regression routines after all original 222 records. The ten checks cover input
and immutable ownership, independent normal/speed and area oracles, arbitrary
reference times, maximum table size, continuous support admission, delayed
shock classification, actual canonical ambient and signed sectors, full
piston history, and attached/detached histories.

The independent sphere-cap oracle starts from exact input doubles promoted to
long double and derives the exposed intersection-circle area. Normal and tiny
caps pass relative `1e-9` on 2x5, 8x17 and 16x33 grids. The actual canonical
ambient integration also advances an explicitly expanding exterior front and
requires genuine fast oblique RH patches under the unchanged residual gate
with zero particle-source area.

Run the complete maintained gate from the model directory:

```sh
make -j8 build/sep_coronal_cme_tests build/background_event_ambient_tests
python3 -B test/run_tests.py --all --expect-count 232
```

The integration executable links the actual composite ambient archive only
for this boundary test. The coronal archive itself remains independent of
composite/application/PIC/MPI code. These software gates do not claim the
broader weak-shock, integrated source-area/history, native moving-field SEP,
regional-fluid closure or observational `CMLU` qualification. The assembly
changes no AMPS time-step or sampling contract.
