# BG3D-0 background contract and baseline audit

Status: frozen implementation contract, 2026-10-03.  This document does not
claim that BG3D-1 through BG3D-9 are implemented or qualified.  The normative
physics details remain in `model.md`, the maintained `sep_coronal_cme/model/`
documents, and revision 1.3 of the launch roadmap.

## Provenance and preservation

The writable checkout is `/home/vtenishe/Mars2/AMPS` at Git object
`fde03997f408792a544291146c5f81f057aa32a8`.  The read-only comparison checkout
`/home/vtenishe/Mars1/AMPS/BASELINE/AMPS` resolves to the same object and has a
clean tracked tree.  Mars1 records the independently supplied baseline archive
as `SEP3D_corona_coupling_compile_fix_20261003.tar.gz`, SHA-256
`8fa1e3752df9b4d790532b3587c057d417263e89e77a81ea6e563b7d950d9c2c`.
The archive bytes are not present, so that digest is provenance, not a digest
recomputed during this audit.

Mars1's later launch work is reference material only.  Its ambient provider,
ellipsoid history, matched DBM continuation, upstream sampler, and local shock
assembly may be transferred a type or hunk at a time after independent Mars2
tests.  Its aggregate launch runtime is not transferable: it makes source-plan
and wave readiness prerequisites of background publication.  It also supplies
only ambient plus a local jump, not the spatial model frozen below.

The following particle-only files in Mars2 match the baseline object byte for
byte.  These hashes are preservation sentinels for G1/G2:

| File | SHA-256 |
| --- | --- |
| `src/models/sep_coronal_cme/include/sep_coronal_cme/particle_source.h` | `5b4c07c45318132d1683353ba2025518cb480edcda8599cfdcc2f3afbb773bab` |
| `src/models/sep_coronal_cme/src/particle_source.cpp` | `fb45fe495a2d2b41c2039719c1baddd9c6a64f6cfa1f48c8bb99df4b91b1e812` |
| `srcSEP3D/adapters/source_runtime.cpp` | `215fc41631b98eb2468e78e2559fffe977d4003b242d4440ad15b72af44fb828` |
| `srcSEP3D/adapters/source_runtime.h` | `f7d965a23555729fed634a6f117f799f00037c1cbad31f9c1784bffd12b07e9b` |
| `srcSEP3D/adapters/transport_adapter.cpp` | `dcdda98541fde38cfa169fa13e0b986183bf562cd03572ab897b0967c6244923` |
| `srcSEP3D/adapters/transport_adapter.h` | `b0b898eeea5f70ffabd82ab82177685f60c7a4ad1f71c512cd6286a1db391a11` |
| `srcSEP3D/transport/focused_transport.cpp` | `6547dcf45cf90c326dc2839aed8e011b25aa86979baead608bb471f6e8dd6614` |
| `srcSEP3D/transport/focused_transport.h` | `cbea95f192082163b69b6a02dc24c77ae73d639b6904d20340b2014f6185232b` |
| `srcSEP3D/output/restart.cpp` | `c96745b8df049dc7ddaab78973f0e450a8aba46a8d9490af3661953b210ce1f7` |

Mixed files include `input/sep3d.input`, the shock provider/parser/registry,
shared and application makefiles, `srcSEP3D/main*.cpp`, `SEP3D.h`, runtime
configuration/dispatch, native mover-hook installation, test registry/stage1,
and the coronal application validator.  A change in one of these files must be
reviewed by hunk: background additions may not replace baseline particle,
legacy SWCME, CLI, or generic PIC sections.  Their complete BG3D-0 hash list is
recorded in `CODEX_CME_PLAN.md`; later stages record every intentional change.

## Supported-profile matrix

Only the first row is the selected G1 profile.  A row being parseable does not
make it qualified.

| Profile | Ambient | Front and outer motion | Spatial regions | G1 status |
| --- | --- | --- | --- | --- |
| `pfss-parker-shock-fed-map-v1` | One `sep_coronal_cme` PFSS/open/closed corona with an isothermal Parker continuation, one composition/EOS, coverage from the physical Sun through at least 1.05 AU | One fixed-orientation material ellipsoid history; full position/velocity matching into the sign-aware SWCME constant-wind DBM law; a differentiated quintic transition when direct full-surface matching is impossible | Upstream ambient; mathematical shock; shock-fed material-map sheath; distinct positive-J ejecta material map; explicit contact/transition and post-event disposition | selected; implementation begins at BG3D-1 |
| legacy SWCME `SHOCK_ONLY` | Legacy Leblanc/Parker | legacy sphere/SSE/ellipsoid and DBM | ambient plus front metadata only | compatibility control, never G1 |
| legacy SWCME `FULL_ICME_DIAGNOSTIC` | Legacy Leblanc/Parker | legacy outer geometry/DBM | phenomenological sheath/ejecta; Parker-like ejecta field | compatibility control, not the selected conservative closure |
| coronal launch ambient plus RH front | PFSS/Parker diagnostic ambient | low-coronal ellipsoid and matched DBM ingredients | no downstream volume | BG3D-2/3 intermediate only |
| spheromak or imported MHD evolution | unspecified | unspecified | unspecified | unsupported; not required and must reject if requested |

The selected profile has exactly one event identity, ambient authority,
coordinate frame, EOS/composition record, clock, and background generation.
It advertises region and observable capabilities separately.  Source planning,
particles, release surfaces, random state, and particle feedback are absent
from its background epoch.

## Frozen analytical closures

### Ambient, geometry, and shock

The ambient uses the maintained PFSS, topology, closed-plasma, ion/electron EOS,
and open-tube/Parker kernels.  The exterior is the same provider, not a second
SWCME ambient normalization.  Magnetic nulls, unsupported topology, invalid
stencils, and coverage failures are typed failures; floors do not create state.

The prescribed front is the complete material ellipsoid, including center,
three semiaxes, fixed orthonormal basis, surface identity, normal, area, point
velocity, and acceleration where used.  Solar masking does not change its
unmasked material identity.  The shock solver consumes the canonical external
upstream state and local relative normal speed.  It retains both sides,
compression, normal speed, fast Mach number, obliquity, entropy/admissibility,
and weak/subfast classification.  There is no compression floor.

The handoff crossing is found from the actual dense coronal history.  Direct
handoff is allowed only if the complete surface position and velocity match.
Otherwise the one position construction is quintic-blended and differentiated,
including the time derivative of the weight.  The outer DBM initial radius and
speed are the crossing state; shape, material labels, ambient normalization,
regional inventory, flux, and thermodynamic state are not relaunched.

### Shock-fed sheath

A sheath parcel is labelled by nonzero shock patch lineage, two material
surface coordinates, and its unique shock-crossing time `tau`.  Admission uses
the actual shock area and the relative upstream/downstream normal flux:

`dM = rho1*w1*dA(tau)*d(tau) = rho2*w2*dA(tau)*d(tau)`.

The selected motion is a prescribed regular downstream characteristic map
`x(a1,a2,tau,t)`, initialized from the RH downstream state.  Its deformation
gradient is measured relative to the crossing basis, not a current-area
surrogate.  Density and magnetic field are evaluated by the material laws

`rho = rho2/J`, `B = F*B2/J`, `U = partial_t x`, `J = det(F) > 0`,

with the reference crossing-volume factor included explicitly in `F`/`J`.
Pressure follows the declared species-composition adiabatic law
`p = p2*J^(-gamma)` unless a checksummed nonzero heating asset is selected.
The G1 profile selects zero added heating.  The map must reproduce the RH state
in the `tau -> t` limit, be one-to-one, and account for the oldest-parcel
contact plus solar/outer-domain exits.  The selected full closed surface has no
unrecorded lateral exit; any solar clipping produces an explicit absorbing
boundary flux rather than discarded mass.

### Ejecta, contact, and post-event state

The ejecta is a separate pre-existing body, not a thickened shock or quiet
ambient.  It uses a checksummed reference-domain material map `x_e(a,t)` with
`F_e = partial x_e/partial a`, `J_e > 0`,
`rho_e = rho_e0/J_e`, `B_e = F_e*B_e0/J_e`, and
`U_e = partial_t x_e`.  The reference magnetic field is supplied by a
vector-potential asset and independently differentiated; this makes a
solenoidal field without requiring a spheromak.  Composition and `gamma` are
owned by the same event record and `p_e = p_e0*J_e^(-gamma)` for the selected
zero-added-heat profile.

At the sheath/ejecta material contact the normal velocity, normal magnetic
field, and total normal traction must agree to their interface budgets.
Tangential discontinuities are allowed only when the recorded surface-current
and traction/work ledger closes.  A C1 numerical transition may sample the two
qualified limiting constructions, but direct component-wise magnetic blending
is forbidden.  The selected body covers every supported point behind the
contact down to its declared solar boundary.  Points outside that body after
its trailing passage have an explicit post-event ambient/recovery category;
they are not used to fill a missing sheath or ejecta cell.

### Prescribed force, heat, and work

This is not a global MHD solution.  The map and zero-added-heat law are frozen
before evaluation.  Momentum diagnostics independently form the inertia,
pressure, Lorentz, and gravity terms and report the remaining sustaining force
as a dimensional field and volume integral.  Energy diagnostics independently
integrate kinetic, internal, magnetic, gravitational, boundary/Poynting, and
the work of that same frozen force.  No residual-fitted pressure, field,
heating, or post hoc work term may be introduced.  Changing the map, reference
field, force convention, or heat law changes the closure fingerprint and
invalidates prior evidence.

## Frozen numerical and physical acceptance policy

Absolute normalizations belong to each fixture manifest and never floor a
physical value.  These targets are frozen before BG3D-1 implementation:

| Quantity | Requirement |
| --- | --- |
| Exact EOS/simple geometry | relative error at most `1e-12` plus declared expected-zero scale |
| Matched handoff position/velocity | normalized one-sided residual at most `1e-10` |
| Root/event time | below configured root tolerance and below one tenth of update cadence |
| Well-conditioned RH fluxes | each independently normalized residual at most `1e-9`; weak shocks additionally use strength-scaled residuals and branch/entropy checks |
| Material mass and magnetic flux | three refinements; finest integral residual at most `1e-8` for manufactured maps and observed convergence consistent with the declared order |
| Smooth divergence/induction | three refinements; finest normalized L2 at most `1e-7`, L-infinity at most `1e-6`, and decreasing at the declared order away from interfaces |
| Contact/moving-interface fluxes | independently normalized mass and normal-B residuals at most `1e-8`; traction and energy/work ledger residuals at most `1e-6` in manufactured fixtures |
| Jacobian | strictly positive at every quadrature/probe point; minimum and conditioning reported, with no floor or sign repair |
| Native identical-scalar halo transfer | bitwise equality; interpolated values use the independently qualified interpolation budget |
| Rank invariance | exact identities/categories/generations and physical values within their spatial/interpolation/reduction budgets |
| Prescribed sustaining force | report raw distribution and characteristic normalization; G1 event profile requires the volume-integrated magnitude ratio at most `2.0` and the 99th-percentile local ratio at most `5.0` |
| Force work/energy | numerical ledger residual at most `1e-6`; absolute signed work and its ratio to the event total-energy-change scale are reported and must be at most `2.0` for the G1 event profile |
| Added heating | exactly zero in the selected profile; a nonzero inferred heat residual is a physical failure, not a tunable correction |

Shock/contact position is graded as an interface observable, not by applying a
smooth-field norm across a discontinuity.  All numerical-versus-physical
residuals and dimensional values are retained.  These manufactured tolerances
are not observational accuracy claims.

## BG3D-0 evidence and remaining sequence

The baseline shared suites passed: `SEP_COMMON01`; SWCME PASS 4/4; coronal
PASS 219, SKIP 3, FAIL/ERROR 0.  The source-free legacy outer native control
passed 9/10 at one rank (one remote-ghost SKIP) and 10/10 at four ranks.  It
begins at 20 solar radii and is not launch/handoff/1-AU evidence.  The external
public-header/archive audit `ARCHSCCM01` and both application adapter compiles
pass without AMPS or MPI dependencies.  The Matplotlib compatibility defect in
`CME3D01` was repaired without changing model data/metrics; all 21 mechanics
tests now pass.

BG3D-1 next implements only the owning background input, checksummed assets,
and complete launch/history support for the selected row.  BG3D-2 through
BG3D-9 remain sequential.  G2, G3, and G4 remain closed.
