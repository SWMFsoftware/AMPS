# Explicit grey-facet Diviner boundary interpretation (U10 / I16)

This supporting contract extends existing U10 and I16; it assigns no new test
ID and changes neither observation membership nor the frozen I16 15 K RMSE
gate. The user approved this interpretation on 2026-11-16. Its implementation
and local verification may precede M1 closure under the approved development
exception; formal I16 qualification still requires its declared predecessors.

## Inventory and production path

The D03 GCP processor preserves the native bolometric brightness temperature
in kelvin. `SurfaceTemperature::DivinerTable::Interpolate` interpolates that
quantity, with native quality and gap eligibility unchanged. The production
surface callback supplies kinetic temperature to sticking, release and source
physics. `MoonDivinerObservation` provides the existing grey forward operator
from kinetic temperature to equivalent bolometric brightness. `RunI16` queries
the same initialized production callback, and `integration_common.py` compares
its emitted brightness with the frozen GOOD GCP hold-out records. U10 invokes
the production operators through `common/production_kernel_probe.cpp`.

## Configuration and physical interpretation

Enabled `surface_temperature_mode = diviner-map` requires the explicit input
`diviner_map_kinetic_interpretation = isothermal-grey`. The default is empty:
old map inputs must be reviewed and amended, rather than silently interpreting
brightness as kinetic temperature. The configured
`thermal_infrared_emissivity` is dimensionless and must be in (0, 1].

This named approximation treats each resolved model facet as an isothermal
grey emitter: emitted hemispheric exitance is epsilon sigma T_kin^4. Therefore
T_kin = T_bol / epsilon^(1/4), and the matching observation operator is
T_bol = epsilon^(1/4) T_kin. Both temperatures are in kelvin. The conversion
has no frame or epoch dependence; interpolation still uses the existing lunar
body-fixed latitude/east longitude and SPICE-derived local solar time.

The approximation does not establish a physical kinetic temperature for a
heterogeneous footprint or a wavelength-dependent emissivity spectrum. Those
limitations remain explicit. This policy does not create an observation,
recover missing native uncertainty, or independently validate the boundary.

## Procedure, oracle and failure behavior

U10 tests independently known fourth-root anchors (including epsilon=1/16),
forward/inverse round trips, and rejection of nonpositive/nonfinite input,
zero/out-of-range/nonfinite emissivity, and overflow. These are deterministic
arithmetic fixtures, never observational evidence. The oracle is equality of
the grey Stefan--Boltzmann exitance, rather than current production output.
The inverse returns NaN on invalid/nonrepresentable inputs.

Native parser controls require a missing or unknown map interpretation to be
rejected and the explicit approved interpretation to be accepted. Protocol
fixtures prove parser behavior only, not dataset readiness or runtime physics.

I16 stages the named interpretation in its preserved input and records the
conversion policy. Its native CSV must retain kinetic temperature, configured
emissivity and forward-model bolometric brightness separately. The Python
scorer independently checks their algebraic consistency and compares only
brightness with native D03 brightness. The existing gap handling, immutable
processed-data hashes, hold-out membership, RMSE metric and 15 K criterion are
unchanged. Missing prerequisites remain SKIPPED; malformed receipts or
inconsistent conversion are ERROR; an evaluated excess RMSE is FAIL.

Artifacts include the preserved input, conversion-policy receipt, native model
CSV, D03 QA/hashes, comparison CSV/figures, linked stdout/stderr and run/build
manifest. The final November 16 `development-clean-7` campaign executed the
native U10 prescribed-boundary/gap checks and input controls successfully,
with all ten M0 gates and a zero-byte nightly diff. Formal I16 remains
SKIPPED for its predecessor gate; these native arithmetic/wiring results are
not observational validation. Evidence and source identities are linked from
STATUS.md section 5.14 and the completion report.

## Prescribed production source/impact consistency extension (U10)

The `prescribed-temperature` probe initializes the production state on a 2x4
logical inventory grid and a closed eight-face octahedral terrain fixture,
without AMR or MPI particle execution. It furnishes the actual compiled SPICE
kernels and fixes the start UTC to 2009-01-24T00:00:00. Both the retained
simple-cosine boundary and the explicit grey Diviner-map boundary are tested
at initialization and after a six-hour TDB/ET interval. An analytical complete
periodic brightness table (never observations) provides known interpolation
values. The chosen physical Na Arrhenius parameters are fixture inputs
E=1.85 eV and attempt frequency 1e13 s^-1, with positive inventory stock.

For every logical element, a radial production BVH query obtains the actual
terrain point. At the same explicit epoch, the probe compares cached logical
kinetic K with the ordinary point callback after an actual body-to-LSO SPICE
transform and with an independent radial-cosine or periodic-table grey
calculation. It independently checks population times nu exp(-E/k_B T)
against the normal production desorption rate. The cylindrical Earth-shadow
identity is evaluated independently for the retained simple boundary.
Acceptance is <=1e-12 relative for temperature and <=1e-11 for the exponential
rate, with nonzero changes over the two epochs required for both modes.
These are numerical production-wiring gates, not a model convergence or
observational agreement criterion. The complete table avoids real coverage
gaps; no missing observational cell is filled. Missing SPICE/Na capability is
SKIPPED and runtime/initialization failures remain ERROR.

## Read-only review seam corrections — 2026-11-16

Before the corrective patch, I16 derived its selection longitude from the
legacy geometric Sun-direction cache while the map callback used the named
`illumination_aberration_correction`. It also checked map support at the
nearest archive-bin centre and then queried the callback at actual local time.
Those geometries can differ without reference to any temperature residual.
The correction obtains the Sun's body-fixed SPICE state at the exact callback
ET with the configured correction and checks support at the actual query LST.
Nearest-bin membership logic, the 0.126-hour selection window, quality rules,
missing-data behavior and 15 K RMSE acceptance remain unchanged. This does not
inspect hold-out residuals or select a case because of agreement.

A sparse analytical U10 table demonstrates the seam: an exact native bin node
is supported, while a nearby non-node time needs a missing/flagged adjacent
corner and must be refused. Both queries use the same production interpolator;
no observational values or independent replacement interpolation are supplied.

The furnished-SPICE inventory is also hashed when a prescribed-temperature
probe later returns a metric failure, provided it emitted its actual inventory.
The science failure remains a failure; provenance success cannot promote it.
A protocol-only failed child verifies this artifact behavior separately from
native production verification. The already frozen private campaign is not
modified; the source correction requires its own subsequent native build.
