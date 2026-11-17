# Controlled finite-amplitude, finite-time Lorentz response

This study isolates the unperturbed-orbit approximation used in weak-wave
scattering theory. It compares full relativistic test-particle trajectories
with a separately derived finite-time linear response in **prescribed**
parallel circular waves. It does not measure a Markov diffusion coefficient,
qualify broadband turbulence, evolve unstable wave modes, choose an amplitude
cutoff, or enable production. An identified discrepancy in this experiment is
conditional on its particle, phase, k, time and field ensemble.

## Field and mathematical reference

In the mean-plasma rest frame choose `B0=B0*ez` and one positive k mode per
trajectory:

```text
deltaB = a B0 (cos(psi), -h sin(psi), 0)
psi = k(z-u*t)+phi,       h=+/-1,       u=sigma*v_A
E = -u ez cross deltaB.
```

This explicitly real field has `deltaB²=a² B0²`, satisfies Faraday's law and
has a subluminal electric-field-free wave frame. The study holds it prescribed.
It does not verify its kinetic dispersion, stability, lifetime, or thermal
feedback. A mode mixture selects one fixed circular wave for each trajectory;
it is a probability mixture of separate experiments, not a coherent multimode
field. These ensembles agree at the stated linear second-moment level only;
their nonlinear dynamics are not interchanged.

For charge q and rest mass m the independent full equations are
`dp/dt=q(E+v cross B)`, `v=p/(gamma m)`. Exactly for this one-speed field,

```text
d(E_particle)/dt = u*d(p_z)/dt
E_particle-u*p_z = constant,       Delta K=u*Delta p_z.
```

Pitch is relative to the **mean** B0, not the instantaneous total field or a
fluctuating local-fluid frame. Choose initial gyrophase zero,
`p0=(p_perp,0,p_parallel)`; a uniform relative wave phase supplies the gyro
average. The unperturbed orbit has signed `Omega=q B0/(gamma0 m)`,
`p_x=p_perp cos(Omega*t)`, `p_y=-p_perp sin(Omega*t)` and `z=v*mu*t`.
Direct integration of its parallel Lorentz force gives

```text
Delta = k(v*mu-u)-h*Omega
Delta p_z^(1) = -h*q*a*B0*v_perp*T*sinc(Delta*T/2)
                *sin(phi+Delta*T/2)
Delta mu^(1) = (1-u*mu/v)*Delta p_z^(1)/p
Delta K^(1) = u*Delta p_z^(1).
```

Here `sinc(x)=sin(x)/x`, with its exact value1 at x0. Uniform phase gives
variance equal to half the squared sinusoidal amplitude. At exact resonance
the linear variance grows as T². This is a finite-time coherent response,
not a diffusion coefficient. Dividing it by2T does not establish a Markov law.
Off-resonant finite-time response can be nonzero without filling a finite-band
gap in a long-time QLT coefficient. Only the declared positive k modes enter;
there is no k floor, interpolated power, resonance broadening, or exterior
wave injection.

The distinction is consistent with the weak-wave/orbit approximations in
[Thomas & Pfrommer (2019), Section4](https://arxiv.org/html/1805.11092), while
[Isayama et al. (2023), Section2](https://arxiv.org/pdf/2204.06707) explicitly
studies coherent relativistic wave-particle dynamics. Its two-wave acceleration
model is context, not a substituted reference solution for this single-wave
experiment. The equations above are directly derived from the specified real
Lorentz force. The existing `RWN01/RWN02` normalization oracles remain
unchanged and are not called to generate these expected responses.

## One physical-input authority

`controlled_ep_maxwellian_baseline.json` and its existing `warm.Plasma` helper
supply n/B/Tp/Te, masses, charge, c, mu0 and derived scales. They are archived
with every protocol. No copy of physical constants or inferred electron
temperature is inserted into this study. The helper already distinguishes
ion-inertia and total-mass Alfvén speeds; controls must explicitly select one.
Both values and the actual wave speed are recorded. The pilot explicitly
selects `alfven_total_mass_mps`; this does not change the production model.

The thermal Maxwellian equilibrium supplies context and wave-speed/scaling
inputs. The energetic test particles are separately declared monoenergetic
probes, not samples of the thermal bath. Their species, initial kinetic energy,
mean pitch, branch and real-wave helicity are explicit. Actual binary64
momenta, u, k, amplitudes and observation coordinates are retained as hex
strings alongside raw requested decimal strings. High-precision linear
references evaluate those actual numerical coordinates, rather than silently
substituting a nominal grid. Nonzero decimal/physical-input underflow rejects
explicitly; a positive analytic response too small for binary64 remains
`UNRESOLVED` with its high-precision value retained.

## Numerical/reference identifiability

The orbit core uses dimensionless `pbar=p/(mc)`, `tau=|q|B0*t/m` and the actual
field phase. Its exact phase derivative is
`d psi/d tau=(k*c/Omega0)*(pbar_z/gamma-u/c)`. This avoids cancellation of kz
and kut without changing the Lorentz force. A vectorized RK4 integrates the
full nonlinear trajectories at three factor-two time refinements.

A separately written vector-force RHS is integrated on a declared subset by
SciPy [DOP853](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html).
This supplies an independent ODE-algorithm check, not an independent physical
model. Exact gyromotion, the analytic finite-time force integral and wave-frame
energy/momentum invariant are governing references. The full trajectory's
energy and norm changes use difference identities to avoid rest-mass
subtraction. Its moments include phase-dependent linear-response residual,
joint energy/momentum observables, and a flag for linear pitch predictions
outside[-1,1]; those predictions are never clipped.

The error accounts remain separate:

- **Orbit:** raw three-level variance differences, observed order, actual step
  sizes/counts, sampled invariant residual and the DOP853-subset differences.
- **Phase:** deterministic N/2N/4N uniform quadrature of the nonlinear response.
  These differences are quadrature estimates, not Monte Carlo confidence
  intervals.
- **Mode:** optional M/2M/4M quadrature over a finite uniform k mixture; the
  single-mode pilot records this axis as not applicable. The linear continuous
  mixture uses independent high-precision sinc-squared integration.
- **Reference/serialization:** declared high-precision refinement, exact actual
  coordinates, high-precision values and the conversion error of displayed
  binary64 values. A diagnostic error bound may be rounded upward; physical
  wave power is never floored.
- **Representation:** an explicit angular-resolution estimate from encoded
  momentum ULPs and the independent-orbit subset difference propagates to
  variance uncertainty. Two methods rounding a tiny effect away do not by
  themselves identify a physical discrepancy.
- **Finite time:** T is an experiment coordinate. It is not forced to an
  asymptotic plateau and is not assigned an invented universal bias correction.

All uncertainty values are convergence/representation **estimates**, not
certified per-trajectory rigorous bounds. `relative_response_target` applies
separately to the variance observable and DOP pitch crosscheck. The per-axis
variance criterion is `absolute_pitch_target² + relative_target*max(Vfull,Vlin)`;
the pitch crosscheck uses `absolute_pitch_target + relative_target*sqrt(max V)`.
These are declared numerical controls, never physical-model acceptance limits.
Raw estimates and their sum remain visible. Relative model discrepancy is
undefined when the linear response is zero or below numerical resolution;
absolute response, high-precision values and reasons remain in the record.
Large nonlinear discrepancy is not a failed numerical calculation. Conversely,
an unresolved phase/orbit reference is not evidence for a physical cutoff.

## Source-owned protocol and bounded pilot

The default action prepares and freezes inputs without calculating trajectories:

```sh
python3 -B srcSEP/validation/physics/run_finite_amplitude_study.py \
  --output-directory /tmp/sep-finite-amplitude-prepared
```

An explicit bounded run uses the archived runner/core/helper/input bytes:

```sh
python3 -B srcSEP/validation/physics/run_finite_amplitude_study.py \
  --output-directory /tmp/sep-finite-amplitude-run --run --plot
```

`controlled_finite_amplitude_pilot.json` declares proton1MeV, mean pitch0.5,
branch+ and h+, amplitudes[1e-4,0.01,0.1,0.3], k/resonant-k[0.5,1,2],
observations[0.5,2,8] initial gyroperiods, RK4 steps/period[64,128,256] and
phase counts[8,16,32]. The existing initializer's amplitude0.3 is a study
coordinate, not an adopted valid domain. The complete numerical controls and
resource guards are explicit in that JSON. No broad scan or tighter control
table is launched automatically after seeing a result.

`protocol.json` retains raw inputs, source/helper hashes/archives, dependency
versions and raw-to-actual coordinate mapping. Its relative calculation and
plot commands can be replayed from a relocated bundle. Source, dependency or
coordinate identity changes reject. Study output retains every requested
amplitude/k/time tuple, including unresolved rows. Exit0 means preparation or
all requested responses numerically resolved; exit1 means unresolved rows were
retained; exit2 means protocol/input/software error. None of these exit codes
certifies the physical closure.

The standalone archived plotter produces PNG/PDF and a hash/provenance JSON.
It verifies complete row coverage and actual amplitude/k/weight/time binding,
retains gaps at unresolved rows and leaves zero/undefined logarithmic values
unplotted without floors. Positive plotting underflow rejects. The1/5/10%
lines are reporting references; no accepted physical tolerance or interpolated
domain boundary is inferred.
Computed differences below their combined numerical estimate are displayed as
separate grey resolution markers, rather than identified nonlinear-error
points. The full requested amplitude range remains visible even where every
response is unresolved. Replotting with an explicitly invoked newer plotter
does not recompute or relabel the original calculation: its actual rendering
source is archived separately, and metadata records both actual and original
protocol plotter hashes. The original protocol command still reproduces its
original archived rendering. Source/protocol/study changes during rendering
invalidate the derived artifact's verdict.

### Executed development pilot: partial evidence

The approved table above was executed without modifying its numerical
controls after seeing results:36 requested responses,9 numerically resolved
and27 unresolved, using1,966,080 trajectory updates. All9 resolved responses
are at0.5 initial gyroperiods and amplitudes1e-4,0.01,0.1 across the three k
probes. Six nonzero differences exceed the combined numerical estimate;
the three weak1e-4 differences remain below numerical resolution. At amplitude
0.1 and0.5 gyroperiods the conditional relative variance differences are about
-5.99%, -10.37%, +16.97% for k/resonant-k0.5,1,2 respectively. These are
single-probe finite-time contrasts, not production cutoffs.

All27 unresolved responses exceed the predeclared sampled-invariant target
of1e-9 relative to initial kinetic energy. Some also fail the phase or DOP853
pitch criteria; the separate raw estimates remain in the artifact. A longer
observation's DOP853 subset maximum is a conservative check used for that
experiment's observations. No numerical target was relaxed and no finer or
broader scan was launched in response. This is incomplete numerical evidence
for those requested points, not a physical exclusion of those amplitudes or
times.

The scientific calculation is pinned to core
`d1465d3e682d332bcec4bed347c614e14eb3291c64344e64803af03b1c2abe68`,
runner `466eb5ee6d2fd7072f2d199ae4915b7d25033622fa56eea11fd058386e1fdda3`
and controls `e3abb9131bea2f470dec1fb3ff2698c143f100e52bc3589c9bcd245fe484e453`.
Those numerical sources remained unchanged for the reviewed replot; figures
retain the rendering source separately. Generated evidence is external to the
repository; the protocol commands and archived inputs reproduce its scope.

### Separate version2 numerical refinement

`controlled_finite_amplitude_refinement_v2.json` is a new experiment, not an
overwrite of the pilot. `finite_amplitude_refinement_v2_manifest.json` freezes
the parent protocol/study/control hashes, new control hash and exact permitted
input differences before calculation. The new time steps/gyroperiod are
[256,512,1024], phase counts[32,64,128] and total trajectory-update allowance
50,000,000. Physical inputs, observations, all acceptance criteria, DOP853
tolerances/subset, analytic-reference precision and the100,000 per-trajectory
step guard remain unchanged. Expected work is31,457,280 trajectory updates
for the same36 requested responses.

The execution uses exact parent-archived runner, orbit core, equilibrium
helper, original plotter and equilibrium bytes. Its preparation preflight
checks their hashes, resolved physical coordinates and dependency versions
against the parent; the complete old bundle is copied and readback-verified
without modifying the original. This preserves the distinction between a
numerical refinement and a new physical model or an altered analysis method.
In particular, the DOP853 subset maximum across **all** observation times is
still conservatively applied to each observation row. Replaying the archived
protocol reproduces those historical code identities; a fresh checkout must
record its actual source identities rather than relabel a different version.

No remaining unresolved point triggers another automatic refinement. Any
remaining phase failures must first be reported with their exact requested
coordinates and raw32/64/128 phase-variance values, numerical allowances and
other error accounts. Neither an unresolved calculation nor an identified
finite-time contrast selects a physical production cutoff.

The completed version2 experiment retains31 numerically resolved and5
unresolved responses. All9 previously resolved responses remain resolved;
22 of the old unresolved responses now satisfy the unchanged criteria.
Actual work is exactly31,457,280 trajectory updates. All remaining failures
are phase-quadrature failures; the time, sampled-invariant and global-DOP853
criteria pass at these points. Their phase variances are shown below
(dimensionless pitch variance; displayed values rounded, raw binary64 values
and all error accounts retained in `refinement-comparison.json`):

| a | k/resonant-k | T/initial gyroperiod | V32 | V64 | V128 | abs(V128−V64) | unchanged allowance |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 1 | 2 | 0.0322855605216072 | 0.0329758089645995 | 0.0316248399790634 | 0.00135096898553608 | 0.00118018060199995 |
| 0.1 | 2 | 2 | 0.0190867023794503 | 0.0162041450452900 | 0.0183175908146721 | 0.00211344576938216 | 0.0000366351816293443 |
| 0.1 | 2 | 8 | 0.0373033062901130 | 0.0329740128383406 | 0.0315152756473768 | 0.00145873719096380 | 0.0000630305512947536 |
| 0.3 | 2 | 2 | 0.0594247033618759 | 0.0580179644790290 | 0.0570175161026324 | 0.00100044837639651 | 0.000114035032205265 |
| 0.3 | 2 | 8 | 0.0715830100799868 | 0.0701626683533027 | 0.0674814527545214 | 0.00268121559878122 | 0.000134962905509043 |

Their fine/mid time-variance differences range from2.23e-12 to6.26e-10;
sampled-invariant errors range from1.75e-12 to2.13e-11 relative to initial
kinetic energy. These small orbit errors do not resolve the independent phase
quadrature. Some phase differences are nonmonotonic, and no convergence order
or physical cutoff is inferred from them. No further refinement was launched.

The original parent bundle was verified unchanged after computation. The
new protocol hash is
`c62c7a6025e19f82e831edd3c2a469c6c83ff2d62c5b14a4ae2d75920c6bf25d`,
new control hash
`8e136427941ba6efa16495385ee064f6c7fb08efd28b279042bd4af6df95b510`,
and new study hash
`f00601a8f4def2b87bd0a59a270a6596d660de1d19e1cc0b72d61983876d02a8`.
The comparison also verifies unchanged actual k/time/momentum coordinates,
analytic references, physical equilibrium and all execution-code identities.

Run bounded mathematical/protocol checks with
`srcSEP/test/run_finite_amplitude_study_tests.sh`. They verify exact gyro/order,
the finite-time real-wave normalization, both branch/charge/helicity signs,
nonlinear DOP853/invariant agreement, independent phase/mode quadrature,
finite-gap semantics, underflow/undefined-relative handling, raw-coordinate
binding, complete unresolved coverage, and portable pinned plotting. The
test runner hashes its full local source/input closure before and after work;
a concurrent edit invalidates its verdict. Generated evidence stays outside
the repository.

This is a partial closure-domain study. A production-domain recommendation
still requires an adopted observable error budget, broader particle/pitch/
spectral/coherence studies, the random-phase/Markov limit where invoked, kinetic
wave-dispersion evidence, and complete coupled/native convergence. The warm
study's unresolved high-k modes and the current initializer's high-k support
remain independent open issues; this prescribed nondispersive orbit experiment
cannot remove them.
