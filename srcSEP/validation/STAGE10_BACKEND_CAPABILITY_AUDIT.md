# Stage 10 backend capability audit

Baseline: commit 80466e6, with the already appended stages 10–21 roadmap and
the user's existing root Makefile edits preserved separately. Evidence:
`/data/vtenishe/TEST/sep-stages10-11-20261123T030052Z`.
The frozen source copy, baseline hashes/patch and actual process results remain
outside the source distribution.

## Executed baseline and limits

The exact five central core/geometry catalogs were executed from the frozen
baseline, independently of later interface edits. Result: **523 cases:
341 PASS, 181 SKIP, one ERROR**. The full-regional
CSWC-FULL-CSWC0628-MASS-TIME case exceeded its unchanged 180-second ceiling.
No partial refinement output was accepted. Core closure/native capability
remains unqualified; this audit records that fact rather than closing G1.

| Registered suite | PASS | SKIP | ERROR |
| --- | ---: | ---: | ---: |
| CORE-PIC-FIELD-LINE-001 | 11 | 0 | 0 |
| CORE-MODEL-SWCME-001 | 4 | 143 | 0 |
| CORE-MODEL-SEP-CORONAL-CME-001 | 219 | 3 | 0 |
| CORE-MODEL-SEP-CORONA-SWCME-FULL-001 | 45 | 0 | 1 |
| CORE-MODEL-SEP-CORONA-SWCME-REDUCED-001 | 62 | 35 | 0 |

The 46-case full-regional central catalog is the bounded executable import,
not every historical/expanded campaign listed in the core README. The 147-case
SWCME catalog reports dormant/observational cases as SKIP. Neither that status
nor the skipped native reduced-core rows is a physical PASS.

The existing three real 400-step stage-9 applications were rechecked with the
unchanged independent verifier while their pre-stage-11 compiled source/build
hashes still matched. SAPP01–03 PASS; dt=300 s, 400 complete steps, original
1e-10 bound and references. These are **revalidated prior executions**, not
new simulations. The stage-11 delivery requires new matching-build evidence.
PCLI01–12 additionally executed before edits on two MPI ranks, two repeats each:
12 PASS, actual exit 0.

## Frozen profile inventory

The identifiers below are executable profile identities, not interchangeable
labels for entire backend families. Revision 1 of each preserves its current
physical meaning. Only the first row currently qualifies material transport.

| Shared profile key | Existing software route | Material admission and missing physics |
| --- | --- | --- |
| reduced-sse-parker-exterior-v1 | Actual reduced provider/interval, native moving line, source, three movers and restart | Qualified Parker exterior, prescribed turbulence, finite-reference net source and absorbing returns; no physical downstream volume |
| swcme1d-shock-only-legacy-v1 | Canonical 1-D adapter, configuration, checked radial primitive/derivative queries | Rejected: constant radial wind with the complete Leblanc density lacks source-free material continuity |
| swcme1d-full-icme-diagnostic-v1 | Same adapter with separately resolved FULL_ICME mode | Rejected: conservative material evolution and regional thermodynamic export/qualification incomplete |
| swcme3d-sphere-shock-only-legacy-v1 | Canonical Cartesian core and mesh-provider components | Rejected: no qualified conservative material field-line lifecycle |
| swcme3d-ellipsoid-shock-only-legacy-v1 | Canonical self-similar ellipsoid geometry/local shocks | Rejected: independent 3-D material/induction/consumer qualification missing |
| swcme3d-sse-shock-only-legacy-v1 | Canonical finite-SSE cap, local shocks and connectivity | Rejected: finite cap/connection support is not material background evolution |
| swcme3d-sphere-full-icme-diagnostic-v1 | Canonical regional diagnostic evaluator | Rejected: regional material/thermal closure and field-line lifecycle unqualified |
| swcme3d-ellipsoid-full-icme-diagnostic-v1 | Canonical nonspherical regional diagnostics | Rejected: shape-specific conservative evolution/joins/consumer qualification missing |
| swcme3d-sse-full-icme-diagnostic-v1 | Canonical finite-support regional diagnostics | Rejected: regional closure/support and material consumer unqualified |
| coronal-bundle-snapshot-v1 | Neutral schema-1 codec, immutable import and bounded stable-line interpolation | Snapshot only: no temporal primitive sequence, material vertex lineage or full grad(U) |
| coronal-ellipsoid-live-planned-v1 | Shared ellipsoid/piston, ambient, local RH and source components | Planned consumer: background G1/G2 and architecture reconciliation remain prerequisites |
| full-corona-cross-ray-planned-v1 | Shared ambient/local shock and bounded piston/cross-ray/material kernels | Planned regional provider: BG3D-4/5/6 closure, handoff and native lifecycle remain unqualified |

For every row the frame/units, density basis, composition/EOS, physical
coverage, region/side and availability must be declared by its actual owner.
A table row cannot manufacture absent fields. Core diagnostics remain callable
under their established contracts; a common catalog's maximum capability mask
does not certify a particular sample or expose a missing implementation.

## Physical/numerical ownership matrix

This second matrix records the actual consumer contract separately from code
availability. “Unqualified” means no admitted material evolution/restart; it
must not be replaced by a plausible adapter law. SI/HCI is required for the
new neutral consumer; a legacy diagnostic is not relabelled automatically.

| Profile(s) | Geometry / existing evolution owner | Regions, density and thermodynamics | Derivatives, source / return / restart |
| --- | --- | --- | --- |
| Reduced Parker v1 | Finite-SSE leading surface in inertial HCI; existing core-driven radial Parker material trajectory | Qualified smooth Parker exterior; canonical nonrelativistic fluid rho, explicit e/p/alpha number densities and configured isothermal composition/EOS; local RH limits are separate | Actual full grad(U), grad(B), branch provenance; existing accelerated net-upstream-escape source, delete-on-return boundary and material schema-5 restart |
| SWCME1D SHOCK_ONLY | Canonical radial shock evaluator; constant radial wind with Leblanc density, no qualified material evolution | Upstream/local diagnostic shock state; no new proper-density claim or thermal law; canonical EOS stays its owner's responsibility | Canonical diagnostic derivative path; source/return/material restart unqualified, so no native transport |
| SWCME1D FULL_ICME | Canonical diagnostic regional evaluator, not a conservative material mapping | Diagnostic sheath/ejecta/post-event values do not qualify regional species/thermal export | Diagnostic gradients are not a production full material contract; source/return/material restart unavailable |
| SWCME3D sphere/ellipsoid/SSE SHOCK_ONLY | Canonical selected surface and local jump queries; finite SSE angular support retained | Canonical diagnostic ambient/local limits; shape identity does not supply conservative rho/U/B or thermal evolution | Existing Cartesian primitive/derivative APIs remain diagnostic; no qualified moving-line source/return/restart for any of the three shapes |
| SWCME3D sphere/ellipsoid/SSE FULL_ICME | Canonical shape-specific regional diagnostic evaluation | Regional state is unqualified for material conservation, composition/EOS and physical volume transport | No complete branch/region-aware material consumer, injection/return or native restart; all three profiles reject |
| Coronal bundle snapshot | Stored line geometry at one primitive epoch; history covers connections, not evolving plasma | Stored neutral primitives/temperature retain schema-1 semantics; no declared material species/density basis | No complete grad(U), dense physical time or trajectory; no material source/return/restart |
| Coronal live ellipsoid planned | Existing ellipsoid/piston/core ambient/local RH components; no qualified live material consumer | Intended coronal/open/exterior scope awaits G1/G2 and architecture reconciliation; no adapter-created heating/volume | Required derivative/material and source/return/restart consumer contracts remain unavailable |
| Full cross-ray regional planned | Existing core ambient, piston, cross-ray and material kernels; production evolution/joins still open | Open/closed corona, sheath/ejecta, transition and exterior require the selected BG3D closure; no guessed regional density/EOS | G1/BG3D-4/5/6 and source/boundary/handoff/restart qualification incomplete; preserve the original MASS-TIME error |

The exact twelve machine profile names/IDs split the two three-shape groups
above; none inherits a group-wide success. Reduced material propagation was
wired/qualified before this task. SWCME1D diagnostic adapter behavior is tested
by SBACK; the other SWCME shapes' presence in compiled archives is not claimed
as new shape-specific physical execution. Coronal/full/reduced core catalog
results above qualify only their actual case scopes. Neutral snapshot import
is wired, whereas the planned live coronal/full consumers are not wired.

## Source ownership and actual wiring

- SWCME1D: `src/models/swcme/swcme1d.hpp` and its resolver own model physics.
  `srcSEP/adapters/swcme1d_adapter.cpp` privately consumes them.
  `main.cpp` configures the canonical owner; `transport_common.cpp` is the
  shared three-mover entry. These sites must check the resolved profile.
- SWCME3D: `swcme3d.hpp/.cpp` supply sphere/ellipsoid/SSE and regional queries.
  The canonical archive member is linked; that is not native field-line
  qualification. The existing Cartesian provider code retains a sphere-only
  consumer check and cannot stand in for a general moving-line adapter.
- Coronal-CME: `src/models/sep_coronal_cme` owns ellipsoid/piston, upstream
  state and local RH. The current field-line consumer imports neutral bundles.
  The later launch roadmap targets a live shared provider; its architecture
  amendments and G1–G4 sequence must be implemented explicitly.
- Full composite: `src/models/sep_corona_swcme` owns the common regional
  state/contact, cross-ray/material construction and handoff. Core kernel
  verification does not close its missing production evolution.
- Shared contracts belong in `src/models/sep_common`; they must remain
  independent of PIC, MPI, either application and specific core model types.

## Decision records and frozen acceptance

**D10-1 — SWCME conservation: first controlled upstream case approved.** The
baseline retained the Leblanc density and constant radial wind as a legacy
diagnostic. The approved first material candidate is a distinct, prescribed
outer-heliospheric `SHOCK_ONLY` profile. In HCI, the wind is purely radial at a
constant speed `U0 > 0`. Its mass density is `rho_ref * (r_ref/r)^2`. The
signed field is `Br = Br_ref * (r_ref/r)^2` and
`Bphi = -Br * Omega_F * r * sin(theta) / U0`, where `theta` is colatitude
about the declared rotation axis and `Omega_F` is the declared field-line
rotation rate in rad/s. Use the declared reference field magnitude, polarity
and latitude to set `Br_ref`. The finite-source-radius offset from legacy
diagnostics does not enter the new ideal-induction expression. These laws
make `r^2 * rho * U0`, `r^2 * Br` and `r * U0 * Bphi` constant along each
smooth radial ray.

For elapsed time `tau`, a parcel starting at radius `a` on that ray reaches
`r = a + U0 * tau`. The radial deformation factor is one and both transverse
factors are `r/a`, so `J = (r/a)^2`. The density and field above must agree
with direct transport of mass and magnetic flux by that map. Check the whole
trajectory against the declared profile coverage. The full Cartesian velocity
gradient is `U0/r * (I - e_r e_r^T)`, with entries in inverse seconds. No change
to the AMPS physical timestep follows from these background relations.

The momentum equation includes the sustaining force density required by the
prescribed flow, in N/m^3. Once pressure and gravity inputs are declared,
account for the full vector balance against pressure, gravity and the Lorentz
force. Record the associated work rate `U dot f_ext` in W/m^3 separately from
thermal heating. Force and work are tested independently over the declared
support. The legacy diagnostic remains a distinct profile.

The user approved a closed provider interval of **0.30–1.20 AU** and a
controlled reference state at **1 AU**. The case takes the radial speed,
density normalization, proton temperature and field magnitude from the existing
FAST preset as explicit inputs to the new material profile. The reference
electron temperature and composition are also declared. This is a controlled
verification case. Event decks own their physical inputs and coverage.

| Quantity | Approved controlled value |
| --- | --- |
| `U0`, radial speed | `400000 m/s` |
| Number density basis and value | Electron density `n_e,ref = 6.0e6 m^-3` |
| Composition and EOS | Charge-neutral electrons and protons, zero alpha abundance, scalar multispecies ideal-gas pressure |
| Reference temperature | `T_e,ref = T_p,ref = 120000 K` |
| Thermal law | Adiabatic smooth expansion, `gamma = 5/3`, zero heat input |
| Reference magnetic magnitude and latitude | `5.0e-9 T` at the equator. Test both radial polarities. |
| Field-line rotation for the controlled case | SWCME's Carrington sidereal rate, `2.86533e-6 rad/s` |

Charge neutrality gives `n_p = n_e` and the existing SWCME mass convention
gives `rho = m_p * n_p`. The pressure is `k_B * (n_p * T_p + n_e * T_e)`.
Advect material entropy with `p/rho^gamma` constant along each smooth
trajectory. Both species temperatures scale as `(r_ref/r)^(4/3)` for this
case, so the total pressure scales as `(r_ref/r)^(10/3)`. Mechanical force
work is accounted for separately from the zero thermal source. The entropy
invariant applies up to the shock, whose one-sided states use the canonical
jump solver.

The reference values calculated from the checked-in SWCME constants are
`rho_ref = 1.003573156e-20 kg/m^3` and `p_ref = 1.988134560e-11 Pa`. For positive
radial polarity at the equator, `Br_ref = +3.411269 nT` and
`Bphi_ref = -3.655577 nT`. The opposite polarity reverses the two signed field
components. The following rounded values make useful hand checks. Use the
closed expressions and the declared constants for the actual test oracles.

| Radius (AU) | Electron density (cm^-3) | Each species temperature (K) | Equatorial field magnitude (nT) | Sound speed (km/s) |
| ---: | ---: | ---: | ---: | ---: |
| 0.30 | 66.667 | 597521 | 39.8135 | 128.221 |
| 0.50 | 24.000 | 302381 | 15.4803 | 91.2136 |
| 1.00 | 6.000 | 120000 | 5.0000 | 57.4609 |
| 1.20 | 4.1667 | 94104 | 3.8590 | 50.8845 |

Provider coverage and particle-domain boundaries are separate. The case must
keep each required vertex, insertion, shock and release query covered for the
whole fixed physical step. Check the endpoints, out-of-range requests and
trajectories that would leave the interval. Independently verify mass flux,
induction, signed field, entropy, all velocity-gradient components, sustaining
force and mechanical work, and the local fast-shock states. Gate numerical
acceptance before running the candidate. `FULL_ICME` still requires a separate
regional material and contact closure. The controlled profile has approved
inputs and awaits implementation and qualification.

**D10-2 — Full regional evolution: OPEN.** The later selected scope retains the
3-D ambient and analytical/semi-empirical cross-ray evolution without a global
MHD runtime. It supersedes the old fixed-ray completeness claim. The candidate
ambient-frame deformation is still a representability investigation; a complete
production map, inverse, force/heating/work and regional joins remain to be
selected and independently qualified in stage 15.

**D10-3 — Coronal linkage and replay: RECORDED GAP.** Current baseline architecture
requires neutral snapshot import, while the launch roadmap specifies a later
live consumer. Stage 11 can define neutral contracts without linking that
consumer. Changing linkage requires synchronized authoritative architecture and
tests at the later integration stage. Executable bundle schema 1 and documented
schema 3 must not be reconciled by relabeling bytes.

**D10-4 — New-backend physical source mapping: OPEN.** Retain the approved DSA,
upstream-frame parent isotropy, explicit finite release reference, net escape
and deletion on return. The SWCME relative-only weight is dimensionless; no
absolute rate or reference-surface conversion is inferred. A new model requires
a concrete physical mapping before native allocation in stage 20.

**D10-5 — Region/initial data and optional passage: CONDITIONAL.** Use already
specified core composition, EOS, assets and boundary laws where supplied.
Request missing ejecta inventory, topology, heat/work or boundary data before
dependent work. Default absorbing return and rejection of pass remain in force.

Original numerical/physical/runtime criteria remain controlling: fixed AMPS
dt, 1e-10 application number/energy budgets, complete physical output windows,
exact promised-deterministic continuation, original core residual/refinement
and 180-second routine ceiling. Stage-11 contract checks use exact metadata,
finite/positive supplied values, explicit absence, closed coverage and strict
catalog/schema matching. They create no new constitutive error allowance.
New closure-specific physical budgets belong to the relevant later decision,
frozen before its candidate residuals are examined.

Stage 10 completion means the baseline/catalog/matrix and these open decision
records are established. It does not assert that every audited core gate passed.


## Additional checkout-procedure findings

The broader source checks exposed PATH selecting a Python without the required
NumPy/mpmath reference dependencies. The unchanged shared SEP-common catalog
reruns with system Python: 91 PASS / 52 explicit SKIP. Two source-mode-0644
POSIX test launch recipes and the obsolete driver-only analytic ownership
marker are repaired in stage 11 without changing numerical assertions.

The separate packaging self-test target remains a baseline mismatch: its
`--root --self-test` invocation is unsupported by the current overlay helper.
Both that helper and the old recipe are unchanged from the frozen starting
revision. Preserve this ERROR until a faithful packaging/self-test contract is
restored; removing options or ignoring the exit is not qualification. These
procedure findings do not supply any missing backend physics.
