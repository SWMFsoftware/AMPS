# Rounded affine frame inversion — pre-change regression contract

The actual approved physical Na commissioning run failed its unchanged planner
momentum gate. GDB retained the supplied SPICE rounded3x3 rotation, approximately
29.386km/s velocity offset, zero LSO parent velocity and accepted source/product
moments. AxisY residual−3.6607511516242e−37kgm/s exceeded its fixed
3.13708920657636e−37kgm/s bound. The independent prediction from
`m_parent*(I−R R^T)*offset` is−3.66078408676361e−37kgm/s. Source diagnosis:
`test_output/srcMoon/m3-p02-approved-20261119T175502Z/momentum-rounding-diagnosis.json`.
This is a representation/inversion error, not permission to alter any bound,
rate, physical state, source photon model or supplied SPICE coefficient.

Production entry: `PhysicalReactionEvents::SodiumPhotoionization::Compose`.
Replace its transpose-as-exact-inverse operation with a small long-double solve
of the ACTUAL supplied affine relation. Keep existing proper-rotation checks
and supplied matrix/offset. Explicitly verify the recovered parent's forward
roundtrip with the existing source numerical momentum bound after mass
weighting. An unrepresentable/singular solve or failed roundtrip returns a
transactional numerical error; no bath is adjusted to cancel a residual.

Extend the actual production probe with the recorded matrix/large offset and
zero event-frame parent, both source endpoints and both existing energy
conventions. Retain all prior732 checks and existing arithmetic/physics gates.
Use an independent70-digit Decimal adjugate inverse rather than the production
pivoted linear solve. Record the original transpose prediction and prove it
exceeds the original axisY commissioning bound. Include a numerical-refusal
case with unchanged zero output bounds, proving the new roundtrip diagnostic
never inserts an epsilon. New source-derived matrices are frame-arithmetic
fixtures, not observations or fake physical sodium states. Native state masses
in the component probe remain explicitly mathematical as before.

Expected: corrected production composition and independent inverse oraclePASS;
original-transpose reproducerFAIL under its original fixed bound; deliberately
unrepresentable roundtrip refuses without output mutation. Artifacts retain
all metrics, source/data SHA, CSV/PNG/EPS, exact commands and prior failed
commissioning evidence. Formal M3/observational qualification remains separate.

## Implemented scope and execution

Production now uses deterministic partial-pivoted long-double elimination;
the independently coded oracle uses a70-digit Decimal adjugate inverse.
All four quiet/active × strict/source-effective energy-convention fixtures use
the recorded rounded matrix/large offset and a zero LSO parent. In addition
to ordinary relative vector comparisons, twelve direct mass-weighted forward
roundtrips use the original1e−38+1e−12 source momentum bound; a velocity-relative
tolerance alone would miss this cancellation-sensitive error. Four independent
calculations reproduce the original axisY commissioning-bound violation of the
transpose operation. These local criteria do not establish milestone closure; the actual
final scoped campaign below records981 criteria.

The initial plan mentioned an exact-zero-bound positive-provider fixture.
No positive finite-precision photon provider is assumed to exist under exact
zero arithmetic bounds. That optional fixture was not constructed. Existing
transactional source/frame/domain refusal controls are retained; numerical
solve/roundtrip failure returns before the caller's Composition is published.
No implicit epsilon or positive zero-bound admission was introduced.

Run from the source checkout with system SciPy/matplotlib available:

```sh
/usr/bin/python3 -B srcMoon/test/test_m3_physical_reaction_events.py \
  --output-dir /absolute/NEW-rounded-affine-verification
```

The actual final local run produced981PASS,0FAIL/ERROR/SKIPPED, retaining all
732 previous criteria by exact name. Its artifact root is
`test_output/srcMoon/m3-p02-approved-20261119T175502Z/rounded-affine-inverse-final`.
The native composer regression uses the existing labelled analytical masses;
the separate corrected physical commissioning run still has to execute from
the fresh source-matched native producer. Source D04 files and supplied SPICE
coefficients are unchanged. Formal M3 status remains separate.
