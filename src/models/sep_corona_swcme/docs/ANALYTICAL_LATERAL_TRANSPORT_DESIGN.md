# Analytical lateral material transport: next design and rejection gates

Date: 2026-11-10. **Research/design only. No production map is selected by
this document.** The user retained the 3-D ambient and prohibited an MHD
runtime. The implemented Cauchy, contact, geometry and conservation kernels
are necessary ingredients; none alone determines a perturbed lateral flow.
This document identifies a parameter-free reference construction to investigate
before asking the user to approve a new physical velocity ansatz.

## Preserve the actual ambient through material coordinates

Let Phi_t(a) be the material flow map of the existing canonical ambient:

```text
partial_t Phi_t(a) = Va(Phi_t(a),t), Phi_t0(a)=a
Fa = partial Phi_t / partial a
partial_t Fa = grad(Va)(Phi_t(a),t) * Fa, Fa(t0)=I.
```

This reference is fixed by already specified Va, not a guessed lateral CME
velocity. It includes corotation, field-aligned low-coronal open flow and the
Parker exterior with its time/frame conventions. It can be integrated as a
small kinematic ODE/variational system without solving MHD. Validate density
continuity and ideal induction separately before assuming that Cauchy transport
of initial rho/B reproduces the analytic ambient at later epochs. Isothermal
pressure additionally requires the ambient heating/work budget; adiabatic
transport alone is not a thermal reference solution.

Start with smooth single-topology branches. Never differentiate a parcel
through a closed/open or magnetic-sector interface as though it were smooth.
A production atlas must retain one-sided interface states, slips, traction,
flux and coverage explicitly. Inverse flow, branch identity, Jacobian,
lineage and exact time support need independent qualification.

A second variational derivative, needed by later acceleration/metric work, is
not supplied by silently setting it to zero:

```text
Ha_i,jk = partial^2 Phi_i / (partial a_j partial a_k)
partial_t Ha_i,jk = (partial_l Va_i)*Ha_l,jk
                 + (partial_l partial_m Va_i)*Fa_l,j*Fa_m,k.
```

The existing ambient gradient API supplies first derivatives only. A qualified
branch-safe second derivative or a separately validated directional-derivative
construction is required before production use of this metric term.

## Candidate reduced deformation in the ambient material frame

A possible continuation of the existing scalar piston is the composition

```text
X(a,t) = Phi_t(Psi(a,t)), Psi(a,t)=f(r0,q,t)*q, a=r0*q
F = Fa(Psi,t) * Fpsi
Fpsi = f_r*q*q^T + (f/r0)*(I-q*q^T) + q*(grad_Omega f)^T/r0
Jpsi = f_r*(f/r0)^2, J=Ja*Jpsi.
```

This is a CANDIDATE reduced ansatz, not an adopted physical law. Lateral
transport in fixed HCI rays arises from the actual reference flow/geometry;
material labels do not have to remain on a fixed observation ray. At f=r0
and f_dot=0 it reduces to the undisturbed ambient reference, avoiding the
original fixed-ray radial mass-flux inconsistency without a mass source.

With g=partial_f Phi_t(f*q), the physical velocity is

```text
v = Va(X,t) + g*f_dot.
```

The chosen corotating contact fixes only its physical normal velocity. At an
eligible intersection this gives

```text
f_dot_contact = (Vc,n - Va.n)/(g.n).
```

A vanishing/grazing g.n rejects capability; it is not divided by a floor.
No tangential surface parameterization may change the leakage condition.
The same common contact, its startup transition, signed B and manufactured
threaded-interface requirements apply.

The scalar acceleration would include the complete moving reference metric:

```text
acceleration = Da(Va) + 2*(grad Va)*g*f_dot
             + (partial_ff Phi)*f_dot^2 + g*f_ddot.
```

Only a separately specified projected momentum/pressure/heat/shock closure
could determine f_ddot. Omitting metric, transverse-field or work terms to
reuse the old radial solver would not implement this construction. Full
Cauchy B, density and material entropy belong to the same F/J. Physical
force/work residuals in directions not evolved by a scalar ansatz must remain
reported and subject to the original bounds.

## Mandatory representability check before adopting the candidate

A one-scalar perturbation spans only g. A general oblique RH downstream
velocity need not lie in `Va+span(g)`. Therefore compute the actual canonical
RH state at every required shock patch and test

```text
delta = v_RH - Va
unrepresented_velocity = delta - g*(delta.g)/abs(g)^2.
```

Use the existing applicable RH/state acceptance bounds. A nonvanishing
unrepresented component beyond those bounds rejects the candidate on that
patch. Do not replace it with ambient transverse velocity, reduce ray
coverage, mark the physical gate inapplicable, or relax the shock residual.
This necessary test can reject a proposed closure before implementing a large
production solver. Success is only necessary, not sufficient: induction,
traction, entropy, normal shock motion, contact and force/work gates remain.

If the candidate fails, a richer explicitly specified vector material map and
its physical closure will be needed. Its choice requires a concrete design
and user clarification under the instruction not to fill physical gaps with
assumptions. Conservation alone does not select that map. An MHD runtime
cannot be introduced as a fallback under the current authorization.

## Ordered next actions

1. Implement/test only the ambient reference-flow map and first variational
   gradient on smooth branches: exact rigid corotation, independent radial
   Parker mass/flux, step refinement and typed branch/support failures.
2. Establish actual spacetime ambient mass/induction consistency and inverse
   identity, then qualify second reference metrics before force use. Keep
   interface and thermal/heating work open until they have real evidence.
3. Evaluate the candidate's RH representability and startup/contact constraints
   on the actual retained fixture. Publish raw residuals/rejections.
4. If supported, present the complete reduced physical closure, sources,
   boundaries, map/inverse, shock algorithm and independent gates before
   production selection. If unsupported, document the physical conflict and
   request the needed richer closure decision.
5. Preserve the original failing CSWC0628 matrix and all acceptance gates.
   No new local PASS, reference-flow result or prescribed map closes BG3D-4
   or substitutes for native G1 propagation.

## Primary-source context and limits

[Kleimann and Roken, 2024](https://arxiv.org/abs/2207.07419v5) derive analytical
induction solutions for a prescribed flow around an axisymmetric paraboloid,
including passive fields and a mildly compressible extension. This is useful
context for independent kinematic references, not qualification of a CME
sheath with arbitrary magnetic back-reaction, upstream structure or expansion.

[Siscoe and Odstrcil, 2008](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008JA013142)
distinguish expansion sheaths from steady propagation sheaths: expansion
changes lateral deflection, accumulation and thickness. Our inference is that
a steady obstacle template cannot be adopted silently for the selected
expanding ellipsoid, and an invented flank sink cannot replace retained
material accounting. Neither reference selects the required production map.
