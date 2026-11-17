# Neutral material-tube kinematics and geometric conservation

`util/sep_material_tube_kinematics.*` implements two bounded value APIs. It
does not move native vertices, publish a background/wave owner, reattach
particles, solve background momentum/temperature, or select a physical MHD
domain. Every geometry candidate remains `qualifiedForProduction=false`.

## Local deformation law

Let a denote one fixed material reference and `X(a,t)` its current position.
The caller supplies `F=partial X/partial a`, `Fdot` and a declared positive
Jacobian. The module independently computes `J=det(F)`, the cofactor
`C=J F^-T`, `G=Fdot F^-1` and `Jdot=C_ij Fdot_ij`. It rejects orientation
reversal, singularity, unresolved determinant cancellation, an inconsistent
declared J, nonfinite data, and unrepresentable nonzero outputs. There is no
determinant floor or repaired deformation.

`PushForwardPoint` takes reference magnetic B [T], rho [kg/m³] and an oriented
area element [m²]. Its laws are

```text
B = F Bref/J,             rho = rhoref/J,
dA = C dAref,             B.dA = Bref.dAref,
DB/Dt = G B - trace(G) B, Drho/Dt = -trace(G) rho.
```

This is the ideal-induction/source-free-mass Cauchy solution. The primary
[MacTaggart (2012), Section2](https://www.maths.gla.ac.uk/~dmactaggart/papers/dmac12c.pdf)
states the deformation gradient, positive J and magnetic push-forward; it
separately discusses force balance rather than claiming every deformation is
an equilibrium. The area/flux identity follows directly from
`F^T cof(F)=J I`. No pressure, heating, composition, or external-force response
is inferred from this kinematic identity.

Flux preservation must also survive numerical construction and storage. A
huge shear can make `B_current` and `dA_current` individually finite while
their separately rounded dot loses a tiny reference flux completely. The
point API therefore certifies both wide-constructed and encoded binary64
states against the actual signed reference flux. It uses error-free FMA
products and TwoSum expansions only after explicit exponent checks exclude
product/residual underflow, overflow and insufficient addition headroom.
IEEE binary arithmetic and nearest rounding are required; fast-math
reassociation is an unsupported build. The allowance for nonzero flux is
`512*eps64*abs(reference_flux)` using a conservative lower magnitude enclosure.
Exact zero instead uses the absolute reference-scale criterion
`512*eps64*||Bref||*||dAref||`; it is not assigned a nonzero physical flux.
Neither allowance grows with huge cancelling transformed products. Reported
current flux derives from actual encoded B/area values. Unresolved precision
returns no point state, never a repaired field or area. Large deformations
whose states retain the invariant remain allowed; this is not a shear cutoff.

The point API also applies to explicitly supplied nonaffine local maps. It
cannot establish that a set of unrelated sampled matrices comes from a single
compatible global map. That integrability/provider contract remains with the
caller; it is not inferred from a flag or a finite determinant.

## Global affine cell interval

The cell API requires the explicit capability
`X(a,t)=F(t) a+c(t)` over the **whole** reference cell. Reference surfaces may
be curved. Each actual face supplies its outward integrated area and first
geometric moment:

```text
A_j = integral_face n_j dA              [m²]
M_ij = integral_face a_i n_j dA         [m³].
```

The cell also supplies an independent true reference volume V and integrated
material mass m. It checks closed-surface identities
`sum A=0`, `sum M_ij=V delta_ij`, and net signed magnetic flux zero. Material
side-wall magnetic flux must be explicitly zero; both cap fluxes must be
declared. Magnetic flux is a supplied integral of the actual reference field,
not `B_midpoint` times an effective planar face. The local Cauchy/Nanson law
preserves that integral even for a varying reference field.

Affine transformed surface data follow from exact integration:

```text
A_current = C A_reference
M_current = F M_reference C^T + c tensor A_current.
```

A closed curved side can have zero integrated area and nonzero first moment.
It therefore cannot be discarded or represented by an invented normal. The
outward material grid swept-volume rate on **each** surface is

```text
Q_face = Fdot_ij C_ik M_jk + cdot_i C_ik A_k       [m³/s].
```

This is exact in space for the stated affine map and supplied surface
integrals. Its time integral uses caller-supplied positive quadrature weights,
ordered times, actual F/Fdot/c/cdot samples and explicit numerical absolute
error budgets. The module checks that time weights span the interval, that
integrated Fdot/cdot agree with endpoint changes, and that

```text
V_final - V_initial = sum_caps sweptV + sum_sides sweptV.
```

Endpoint volumes are `J_endpoint*V_reference`. They are never defined by the
sweeps. No face is altered to force closure. Per-face rates also independently
sum to `Jdot*V_reference` at every quadrature node. Measured residuals and
allowed numerical bounds are returned; exceeding a declared budget rejects
without a candidate. Passing a caller-declared budget is a numerical check,
not independent proof of its accuracy or physical qualification. Arbitrary
time histories need resolved quadrature/convergence; endpoint agreement alone
does not prove the supplied samples describe the entire path.

The true material mass remains m. Mean density is exactly the declared
material average `m/(J V_reference)`. It is not substituted for a provider's
point density, nor assumed equal to `rho_midpoint*V`. All vectors, moments,
mass, volumes, time and rates use the SI units stated in the public header.

`PrepareAffineInterval` returns a private immutable candidate only after all
checks and allocations succeed. Unknown/nonaffine cell maps, missing numerical
budgets/provenance, malformed reference surfaces and unresolved quadrature
return explicit statuses. There is no mutable global, separate geometry owner,
native write, volume repair or implicit production remap.

## Independent bounded tests

The shared factory is `SEP::Testing::MaterialTubeKinematicsDescriptors()`;
cases use group `material-tube`, initialization None, deterministic repeatable
inputs. Execute `srcSEP/test/run_material_tube_kinematics_tests.sh`. Its exact
local source closure is hashed before/after compile and execution; a concurrent
edit discards the verdict. Strict C++11, AddressSanitizer memory-access checks
and UBSan are enabled. Leak detection is disabled and no leak-check claim is
made.

| Case | Independent reference |
|---|---|
| MTG01 | General triangular stretch/shear with explicit magnetic components; mapped spanning-tangent cross product for area; known full G/rates; exact90-degree rigid rotation. |
| MTG02 | Isotropic linear dilation with analytic J, six individual swept volumes, true mass and cap/side contributions. |
| MTG03 | Isochoric shear plus translation with analytic signed face antiderivatives; cap normal tilts relative to B while flux remains fixed. |
| MTG04 | Actual annular-sector extrusion, R in[2,3]m, phi in[-0.4,0.4], z in[-0.75,0.75]m. Analytic cylindrical area/first moments and elementary time antiderivatives check every cap/side sweep under affine stretch/shear; toroidal `Bphi=(5e-9 T*m)/R` checks nonuniform-field local flux. `rho_ref=(1+R/(1m))*1e-21 kg/m³` integrates to true mass `1.06e-20 kg`, separately checked against mapped mean density. |
| MTG05 | Nonaffine parabolic line `X=(a_x,a_y+q a_x²,a_z)` with explicitly expanded B/rates; whole-cell nonaffine interpolation remains unsupported. |
| MTG06 | Invalid face closure/moments/flux, wrong/singular/reversing J, missing budgets, inconsistent Fdot and an insufficient midpoint time rule reject without candidate/input mutation. Extreme-shear wide/encoded flux loss, tiny/subnormal/overflow flux ranges and unsupported rounding are checked; the same huge shear is accepted when flux remains representable. |

The reference integrals are explicitly derived from geometry or expanded
components. They do not call the production cofactor/GCL routines to generate
expected answers. For linear-in-time affine F and c, the swept integrands are
quadratic; two-node Gauss is exact up to arithmetic. A one-node midpoint rule
is an intentional negative control for isotropic dilation and must not be
accepted under a roundoff-only budget.

## Native and physical dependencies retained

This component is one part of the design in
[MATERIAL_TUBE_WKB_STAGING.md](MATERIAL_TUBE_WKB_STAGING.md). Native provider
maps/quadrature and general nonaffine cell interpolation remain unimplemented.
WKB's fixed-direction interval restriction is not lifted by this more general
local induction API. No arbitrary curved/rotating production provider gains
a WKB capability automatically.

A future outer native transaction must stage vertices/current+previous
primitives, all cell/face geometry, B/rho/U and immutable backgrounds, four wave
arrays, particle state/derived Cartesian positions, boundary reservoirs,
external accounts and ledgers atomically. Material motion requires coordinate
convection U minus grid velocity, never a second U advance. Every resident
particle's Cartesian cache must be refreshed when geometry changes even if its
field-line coordinate does not. Source sampling on a tilted material cap also
needs a guiding-center matching law or a declared aligned-cap domain; the
pointwise geometry law does not establish such a particle-source bridge.
