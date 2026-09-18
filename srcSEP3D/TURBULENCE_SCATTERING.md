# Phase T Turbulence and Scattering Inputs

Phase T makes scattering authority independent of ambient-field authority and
keeps `src/models/sep_common` as the single implementation of transport
coefficient physics.

## Provider contract

`TurbulenceProvider` publishes source/ownership metadata and evaluates a typed
`TurbulenceSample` from position plus a valid background sample. A sample
contains:

- total, along-`+B`, against-`+B`, outward, and inward magnetic variances;
- finite wave-number bounds and spectral index;
- parallel correlation length;
- generation, configuration digest, and source identity;
- separate valid and ballistic state.

The interface is separate from `BackgroundProvider`: analytic Parker fields
may use coupled turbulence, and a coupled background may use a prescribed
spectrum, provided the immutable run configuration selects that authority.

## Prescribed Kolmogorov spectrum

At radius `r`, the default amplitude is

\[
\delta B^2=(0.3|B|)^2,
\]

split equally between propagation directions. Bounds scale independently from
their 1-AU references, and correlation length has its own radial exponent.

For index `q>1`, the one-dimensional spectrum is

\[
P(k)=A k^{-q},\qquad
A=\delta B^2\frac{1-q}{k_{max}^{1-q}-k_{min}^{1-q}},
\]

so `integral_[kmin,kmax] P(k) dk = deltaB²`. `P` has units T²·m when `k`
has units m⁻¹.

## AWSoM wave convention

The coupling record is explicit:

- `wPlusJPerM3`: total Alfvén-wave energy density propagating along `+B`, with
  group velocity `U + v_A b_hat`;
- `wMinusJPerM3`: total energy propagating against `+B`, with group velocity
  `U - v_A b_hat`.

For equipartitioned Alfvén waves, magnetic energy is half the total:

\[
\frac{\delta B^2}{2\mu_0}=\frac{w}{2}
\quad\Rightarrow\quad \delta B^2=\mu_0 w.
\]

The conversion is applied independently to `w+` and `w-`. If `B·r>=0`, plus
is outward; if `B·r<0`, minus is outward. Thus a magnetic-polarity change swaps
outward/inward aliases without changing the field-aligned AWSoM record.

AWSoM supplies integrated directional energies. The finite spectral bounds,
index, and correlation length remain explicit coupling configuration; they are
not inferred from two scalar wave energies.

## Resonance range policy

`NormalizedPowerLawSpectrum::Evaluate` implements two policies:

- `Reject`: a requested `k` below `kmin` or above `kmax` returns a typed error;
- `PowerLawExtension`: evaluate the declared power law at the physical
  requested `k` and mark the result extended.

Extension never clamps to a band edge. The run configuration fingerprints the
chosen policy. Phase-P movers must propagate this state rather than silently
choosing a policy in a particle loop.

## Missing-data policy

An absent or incomplete coupled wave record fails by default. With explicit
`Ballistic` policy it produces a valid typed ballistic sample with zero
variance. `CoefficientBridge` maps that state to `ValueState::Ballistic` and a
zero Dμμ rate. A zero numerical value without this typed state is not treated
as permission for ballistic propagation.

## Shared coefficient boundary

`CoefficientBridge::ToSharedInput` converts the local 3-D background and
turbulence records into
`SEP::Transport::CoefficientPhysics::LocalInputView`. The bridge then calls:

- `EvaluateJokipiiSlab` for finite-band Dμμ;
- `KappaFromMeanFreePath` and `MeanFreePathFromKappa`;
- `IsotropicDmumuFromMeanFreePath` and its inverse.

All are linked from canonical `sep_common.a`. The bridge contains no alternate
coefficient formula, species table, relativistic conversion, or numerical
floor. Source identity, generation, and checksum accompany every local input.

The interface is split deliberately. `turbulence_models.h` declares only the
provider and normalized-spectrum API and therefore needs no `sep_common`
include path. Translation units that actually calculate transport coefficients
include `coefficient_bridge.h`, which adds the canonical shared headers. This
prevents `main_lib.cpp` from acquiring a transitive include dependency that
historic AMPS `Makefile.conf` rules omit when compiling the copied
`AMPS/build/main` application.

`EvaluateLocalScattering` is the AMPS-independent R02 entry point used by the
production resolver. At every accepted particle substep it combines the
currently pinned turbulence sample, background, species mass/charge, momentum,
and pitch cosine, then delegates to `CoefficientBridge`. It returns
`kappa_parallel`, `D_mumu`, and `dD_mumu/dmu` in one validated record. This is
why cell crossing during one AMPS call cannot reuse coefficients from the
starting cell.

## Production storage

Prescribed turbulence is evaluated for every owner-local physical cell during
`amps_init()` even though its variance can be reconstructed from the
fingerprinted configuration and background. When SWMF turbulence authority is
selected, the host must install a loaded `AwsomTurbulenceProvider`; the two
directional magnetic variances are stored at the optional Phase-M wave offset.

R03 prepares a candidate turbulence provider together with the candidate
background generation. Both are evaluated at all owner-local cells and become
active only after collective readiness. `PrescribedKolmogorovProvider` can be
re-based to an R07 saved generation so the first post-restart snapshot and
future coefficient provenance match an uninterrupted run.

Self-consistent 3-D turbulence remains reserved. Enabling it fails during
`RunConfiguration3D::Create`; it will not become a no-op or alias for the
prescribed provider.

## Evidence

| ID | Contract |
|---|---|
| `TUR3D01` | log-space quadrature closes prescribed spectral energy |
| `TUR3D02` | `deltaB²=mu0*w` and outward/inward mapping under both polarities |
| `TUR3D03` | proton/electron resonances below, inside, above band obey policy |
| `TUR3D04` | incomplete waves fail unless explicit typed ballistic mode |
| `COEF3D01` | Dμμ/mean-free-path/parallel-diffusion round trips over six decades |
| `COEF3D02` | 3-D bridge and direct shared Jokipii kernel are bitwise identical |

Run `test/run_tests.py --suite phase-t --rebuild`.
