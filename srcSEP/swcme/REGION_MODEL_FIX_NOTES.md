# SWCME sheath/ejecta region-model repair

## Purpose

This update repairs the optional phenomenological sheath/magnetic-ejecta model
without changing the already validated Parker field, CME/shock kinematics,
finite shock geometry, MHD jump solver, or magnetic-connectivity algorithms.
The region model is secondary to the controlled SEP study, so the principal
design goal is to make it explicit, internally consistent, and impossible for
its assumptions to leak into the SHOCK_ONLY baseline.

## Previous defects

The previous code had several independent region implementations:

1. The 1-D ejecta used `max(1.0,f_ME)` and `max(1.0,V_ME_factor)`, so documented
   values such as `f_ME=0.5` and `V_ME_factor=0.8` were silently replaced by 1.
2. The 3-D model computed sheath/ejecta widths from the apex radius and
   subtracted those absolute widths from every directional shock radius.  This
   is not self-similar at a finite SSE or ellipsoid flank.
3. LE/TE smoothing was implemented differently in 1-D and 3-D; the former LE
   blend had a reversed-side convention and the latter did not consistently use
   the configured LE/TE widths.
4. `sheath_comp_floor` still appeared in the region data model even though the
   MHD shock solver had already made an empirical compression floor physically
   inappropriate.
5. There was no explicit mode that guaranteed an untouched Parker/Leblanc
   background while retaining shock/connectivity/source bookkeeping.

## Shared region module

`swcme_regions.hpp` now owns all dimension-independent region mathematics.
`Params::region_mode` selects:

- `ShockOnly`: analytical upstream state everywhere;
- `FullICME`: RH shock boundary + phenomenological sheath/ejecta regions.

`FullICME` uses local self-similar boundaries.  If public inputs are
`d_sheath` and `d_ejecta` in AU at a 1-AU shock, their numerical values are the
fractions

```text
f_s = d_sheath / 1 AU,
f_e = d_ejecta / 1 AU.
```

For every directional shock radius `R_sh(u)`,

```text
R_LE(u) = (1-f_s) R_sh(u),
R_TE(u) = (1-f_s-f_e) R_sh(u).
```

Thus all region surfaces preserve the same local angular shape as the shock.
The public configuration is rejected when `f_s+f_e >= 1`.

## Interface smoothing

The configured LE/TE smoothing widths are also treated as local self-similar
fractions.  Each value denotes the **total** symmetric transition width centered
on the nominal boundary. CFG03 now rejects a request above 90% of the adjacent
finite layer before preparation. Every accepted effective width is the exact
requested fraction times the local shock radius; runtime construction never
caps or otherwise changes it.

The transition weight is the standard C1 polynomial

```text
s(q) = q^2 (3 - 2q),  0 <= q <= 1.
```

At the outer and inner transition endpoints both the weight and its derivative
match the corresponding pure-region branch.

## Sheath state

At `R_sh^-`, the state is exactly the shared ideal-MHD RH downstream state.  In
the sheath interior a common progress coordinate relaxes density logarithmically
and velocity/magnetic field smoothly toward the leading-edge target.  The
leading-edge radial speed is constrained between ambient `V_sw` and the radial
component of the exact RH downstream flow.  Configuration therefore requires
`V_sheath_LE_factor >= 1` for the forward-shock sheath model.

If the geometric CME surface exists but `has_shock=false`, no artificial sheath
compression is fabricated; the sheath portion remains the local ambient state.

## Magnetic ejecta

The simple target is now exactly

```text
n_ME = f_ME n_up,
V_ME = V_ME_factor V_sw.
```

Both factors may be below unity.  Negative values are invalid configuration.
The baseline ejecta field remains the Parker field; no unvalidated flux-rope
model is introduced by this remediation.

## Compatibility fields

The old `sheath_comp_floor` parameter remains in public parameter structures for
source compatibility but is ignored by physical shock and region calculations.
The corresponding `StepState::rc_floor` mirror is always 1.

**Fix 10 supersedes the original transport-side shock-smoothing statement in
this note.**  The canonical shock diagnostic is still an exact mathematical RH
discontinuity.  When `RESOLVED_COMPRESSION` is selected, however, the transport
fields now represent that jump with one finite C1 layer whose inner endpoint is
the exact RH downstream state.  In `SOURCE` mode the validated `SHOCK_ONLY`
background contains no resolved shock layer.  This distinction prevents DSA
source acceleration and `div(V)` compression acceleration from being applied to
the same particle population.

## Validation

Five tests were added:

- REG01 — SHOCK_ONLY analytical-background identity;
- REG02 — FULL_ICME RH downstream boundary;
- REG03 — exact ejecta density/velocity factors;
- REG04 — self-similar local shock/LE/TE nesting;
- REG05 — C1 artificial LE/TE transitions and 1-D/3-D consistency.

These tests are part of `test_swcme --all` and are described in detail in
`test/README.md`.
