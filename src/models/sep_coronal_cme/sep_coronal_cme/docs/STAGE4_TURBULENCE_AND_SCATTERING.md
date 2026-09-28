# Stage 4: turbulence initialization and scattering coefficients

Stage 4 implements dependency-light wave-energy, mean-free-path, and focused
collision kernels in `turbulence_transport.{h,cpp}`.  It does not claim a
self-generated-wave solution.  Ambient wave initialization and empirical
particle scattering remain separate authorities because a scalar wave energy
does not specify the resonant spectrum, polarization, resonance broadening,
or the treatment of the 90-degree resonance.

## Directional wave state

The authoritative values are physical `w_out` and `w_in`.  After the magnetic
sector is known, the derived field-relative channels are

- positive sector: `w_+=w_out`, `w_-=w_in`;
- negative sector: `w_+=w_in`, `w_-=w_out`.

Thus an ideal HCS sign reversal exchanges field-relative labels without
changing physical energy.  `BuildDirectionalWaveState` rejects negative or
nonfinite energies and reports

`w=w_out+w_in`, `sigma_c_out=(w_out-w_in)/w`,
`sigma_c_B=sigma_B sigma_c_out`, and `delta B^2=mu_0 w`.

A direct-mean-free-path closed-field branch uses an explicitly inapplicable
wave record.  It is not represented by physical zeros.

## Prescribed and WKB initialization

The prescribed closure evaluates the configured power laws for total energy
and correlation length, then partitions the energy with the configured
physical cross helicity.  All reference scales are required to be positive.

For the nondissipative WKB branch, `EvaluateWkbOutwardWave` implements

`w_out=w_out0(A0/A)(v_A/v_A0)[(u_s0+v_A0)/(u_s+v_A)]^2`,

so `A w_out (u_s+v_A)^2/v_A` is invariant.  The baseline sets `w_in=0`; it
does not invent an inward solution by propagating the same expression through
an Alfvén critical point.  The closed-loop provider is separately named and
uses symmetric, exponentially decaying contributions from two footpoints,
normalized to reproduce the specified energy at either endpoint.

## Uncoupled-wave validity

`EvaluateWaveValidity` reports all amplitude ratios required by Section 7.5:

`epsilon_deltaB=sqrt(mu_0 w)/|B|`,
`p_w/p`, `p_w/(p+rho u_s^2)`, and
`p_w/[B^2/(2 mu_0)]=epsilon_deltaB^2`.

The normative force gate uses

`a_w=-(1/rho) dp_w/ds`

and compares its excess above a checksummed absolute derivative uncertainty
with the non-cancelling retained acceleration scale.  A large `p_w/p` with a
negligible gradient can pass; a modest pressure with a sharp gradient can
fail.  Small-amplitude enforcement is explicit and closure dependent.

## Mean free path and collision operators

The single power law and smooth broken power law implement the exact equations
in Section 7.4, including normalization at the reference/break radius and
independent rigidity dependence.  Both reject nonpositive scales and never
fall back to another coefficient model.

`IntegrateParallelMeanFreePath` evaluates

`lambda_parallel=(3v/8) integral[(1-mu^2)^2/D_mumu] dmu`

with midpoint quadrature that avoids endpoint `0/0`.  The pitch-angle SDE is
the Ito process for `D_mumu=D0(1-mu^2)`, including drift `dD/dmu=-2D0 mu` and
reflecting endpoints.  Its normal deviate is caller-supplied, keeping random
stream ownership outside the shared library.  The discrete alternative is an
isotropic Poisson reset with probability `1-exp(-v dt/lambda)`; it is a
separate collision model rather than a disguised SDE step.

Missing coefficient behavior is explicit: focused models either fail or use
the selected ballistic policy.  `none` never requires or fabricates a mean
free path.

## Validation

`make test-stage4` runs all earlier stages plus `TUR3D01--07` and
`MFP3D01--07`.  The tests cover WKB action, directional/sector identities,
nonzero center/halo-ready initialization, closed-loop symmetry, wave-force
and amplitude diagnostics, both mean-free-path asymptotes, the `D_mumu`
integral, reflected SDE endpoints, Poisson autocorrelation, and explicit
missing-coefficient dispatch.
