# Standalone material-transport backend support

Stage 8 of [the roadmap](../INTEGRATION_ROADMAP_STANDALON.md) follows the
user's 2026-11-22 decision: **reject unsupported modes and qualify existing
supported transport**. No new plasma evolution or source prescription is
inferred from this integration request.

| Existing authority/profile | Diagnostic query/exchange | `parker` | `fte-dmumu` | `fte-mfp` | Moving restart |
| --- | --- | --- | --- | --- | --- |
| Reduced SSE front, prescribed Parker exterior, explicit finite upstream reference, net first-passage source, absorbing returns | Immutable core interval | Supported | Supported | Supported | Schema 5, same decomposition and keyed 1↔8 repartition |
| SWCME1D `SHOCK_ONLY`, maintained Leblanc density/constant radial wind | Maintained canonical adapter | Rejected | Rejected | Rejected | Unsupported |
| SWCME1D prescribed adiabatic material upstream, `swcme1d-adiabatic-outer-wind-v1` | Detached immutable core interval with full upstream gradients and local RH limits | Rejected | Rejected | Rejected | Native geometry, source and restart lifecycle pending |
| SWCME1D `FULL_ICME` | Existing diagnostics; native regional thermodynamics unavailable | Rejected | Rejected | Rejected | Unsupported |
| Coronal-CME neutral single-epoch bundle | Immutable snapshot exchange, bounded stable-ID interpolation | Rejected | Rejected | Rejected | Bundle identity alone is insufficient |
| Coronal-CME Stage 14 live ellipsoid/piston core | Immutable PFSS/Parker ambient and local shock epochs | Rejected | Rejected | Rejected | G1 and native material lifecycle pending |
| Full corona–SWCME regional volume | Completed core ambient/local shock queries only | Rejected | Rejected | Rejected | Pending qualified regional/material lifecycle |

“Supported” describes the bounded reduced profile qualified in stages 6–7
and rerun through the complete stage-9 native/application/coupling gates;
actual receipts are in the [status report](STANDALONE_INTEGRATION_STATUS.md).
This table is not an
execution receipt. It does not qualify local particle-wave evolution,
downstream propagation, PFSS material evolution, gross-emission conversion,
additional fronts or observational agreement. SWMF and 3-D Eulerian particle
transport have separate contracts and are unaffected by this standalone matrix.

## Missing physical contracts

For the maintained SWCME1D ambient law, with fixed composition,

```text
n(r) = C2/r² + C4/r⁴ + C6/r⁶,       U(r) = U0,
div(n U) = -U0 (2 C4/r⁵ + 4 C6/r⁷).
```

Thus `r² rho U` is not constant. This is a physical continuity issue for a
source-free material tube even though the empirical diagnostic profile is
intentional. Moving vertices does not remove the residual. Substituting the
reduced model's plasma, dropping the higher inverse powers, reconstructing
density silently, or adding U to relative particle displacement changes the
model. Admission needs an explicit core evolution/source contract and
independent continuity/induction qualification. `FULL_ICME` additionally needs
actual regional pressure/heating export; upstream temperature cannot be
extended throughout the sheath/ejecta.

The separate adiabatic upstream satisfies that prescribed core continuity and
induction contract. Its detached adapter exposes the controlled canonical
background while the existing standalone selection still uses the legacy
diagnostic modes. Native material admission also needs the installed moving-line
owner, source and restart lifecycle for that profile.

The neutral coronal-CME bundle has one geometry/primitive epoch. Connection
and observer history do not serialize material vertex lineage, plasma states
through time, full time-dependent Cartesian `grad(U)`, or an explicit species
number-density/composition interpretation. The executable `sep_common`
codec currently admits neutral schema **1**, while the model documentation
refers to later bundle schema **3**. This discrepancy is recorded rather than
resolved by relabeling bytes or inventing fields. Neither label establishes
material evolution.

The full corona–SWCME core leaves BG3D-4 regional closure, BG3D-5/6 material
state/handoff and BG3D-7/8 native lifecycle open. Immediate Rankine–Hugoniot
limits are not a spatial downstream volume. Adapter work cannot complete those
core qualifications.

## Enforced boundaries

[`sep_backend_admission.h`](../util/sep_backend_admission.h) returns typed,
solver-independent status. Standalone SWCME transport rejects before native
allocation. The shared `LoadParticle()` entry repeats the check for all three
movers and parser-free hosts, including diagnostic run purpose. Isolated
diagnostic operators remain available; that purpose cannot authorize actual
motion on an unsupported background. SWMF is a separate provider enum.

The original two-argument neutral `LoadAndPublish()` symbol preserves snapshot
exchange. The added `BundleUse::MaterialTransport` overload rejects before
reading/publishing a candidate or invoking an allocation callback, retaining
the previous owning view. Unknown enum values reject without fallback. No new
CLI selector advertises unavailable native support.

`SBACK01–05` cover admission, independent continuity, direct-core diagnostic
equivalence, epoch rejection, real neutral codec/identity and transactional
rejection. `SBACK06–08` require actual two-rank standalone rejection before
allocation for each mover. Passing rejection qualifies rejection only.
Reduced-profile propagation remains the separate SAPP/native/restart gates.
