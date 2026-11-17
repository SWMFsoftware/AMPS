# Material-plane Maxwell–Jüttner kernel

This implements the user-selected relativistic thermal flux distribution as a
**full-particle, material-plane mathematical kernel**. It does not enable the
legacy native Maxwellian source. The adopted material grid, proper-density
provider mapping, guiding-center boundary condition, statistical counts, source
rate integration and native checkpoint continuation remain separate gates.
Track overall progress in [STATUS.md](../../STATUS.md). The governing derivation
and native bridge questions are in
[MAXWELL_JUTTNER_FLUX_DESIGN.md](MAXWELL_JUTTNER_FLUX_DESIGN.md).

Current gate: expanded MJ01–06 pass strict warnings/ASan/UBSan and applied
campaign/repeat/order controls. Independent review first exposed encoded cold-
boost crossing loss and ill-conditioned near-light bulk Gamma formation; the
old code failed the new regressions. The correction rejects both domains before
RNG/output commit. The expanded MJ01–06 pass **6/6 natively on two MPI ranks in
build8**, alongside their pure sanitizer/seed/order gates. Native thermal
injection remains disabled.

The native source point is build8, binary SHA256
`049a1959918dd154e5a81d75ed7ebb18a57764f29c09c3917c81b4068f9c58eb`;
its captured source manifest and structured reports are tracked in STATUS.

`util/sep_maxwell_juttner_flux.{h,cpp}` takes explicit proper species density
`n0 [m^-3]`, proper temperature `T [K]`, rest mass `m [kg]`, simulation-frame
bulk velocity `U [m/s]`, boundary velocity and inward normal direction. It
normalizes the direction; the magnitude is not face area. `v_boundary=U` must
hold componentwise, because a nonmaterial plane has a different crossing law.
Sub/superluminal invalid inputs, zero normals and nonfinite values reject before
output/RNG changes. No plasma density field is silently treated as proper
species density.

The rest distribution uses `exp(-gamma'/theta) d^3p'`, with
`theta=kB*T/(m*c^2)`. Its scalar distribution and invariant momentum measure
`d^3p/E` are described by
[Zenitani (2015), equations 22–28](https://arxiv.org/html/1504.03910v1).
[Swisdak (2013)](https://arxiv.org/html/1305.5226v1) provides the corresponding
rest kinetic-energy form. The boundary-flux simplification below is our direct
derivation, rather than adoption of either paper's volume-loading sampler.

For a material plane the Jacobian times its positive relative normal speed is
`[a·v']_+`, with `a=n-Gamma/(Gamma+1)*(beta·n)*beta`. Thus `|a|^2=1-(beta·n)^2`
and, for `x=gamma'-1`, the conditional radial law is
`x(x+2) exp(-x/theta)`. Its normalized mixture is shape-two gamma with
probability `1/(1+theta)` and shape-three gamma with probability
`theta/(1+theta)`, both with scale `theta`. Two/three sums of open-unit
exponential deviates sample these components without a momentum cutoff or a
radial rejection approximation. The incoming rest direction has
`mu=sqrt(Uopen)` about `a/|a|` and uniform azimuth.

Every event consumes six uniforms: mixture, three exponentials, cosine and
azimuth. Shape two ignores the third exponential. The fixed budget is an
explicit random-history convention. `SampleMaterialFlux` copies the complete
keyed state, including a cached normal, and commits only a successful full
sample. It preserves `x` separately from `1+x`, full proper momentum and fully
boosted simulation momentum. Any unrepresentable tail event returns
`UnresolvedCoefficient` with unchanged RNG/output. It never redraws that event,
clips a speed, truncates the distribution or replaces a small value by zero.
The native layer must handle the status before committing source counters or
allocating a particle.

Full boosted output has an additional numerical capability condition. Forming
bulk `Gamma` from the encoded input velocity is conditioned by `Gamma^2`; the
working-precision arithmetic envelope is `64 epsilon(long double) Gamma^2`.
If it exceeds `2e-12`, preparation rejects before draws. A transverse normal's
`|a|=1` cannot itself diagnose this loss. The limit is explicitly numerical,
not a newly chosen physical speed cutoff, and depends on the platform's
actual long-double precision.

After narrowing SI momenta, let `q` denote the wider boosted momentum in units
`mc`, `qd` its encoded value and `delta=qd-q`. The stable exact identities
`deltaGamma=(2q·delta+delta^2)/(Gammad+Gamma)` and
`delta(v·n)/c=n·delta/Gammad-(n·q)deltaGamma/(Gammad Gamma)` avoid subtracting
nearly equal bulk/particle speeds. The code encloses component narrowing plus
the conditioned radial/normal/boost/scaling arithmetic budget, and checks the
separately encoded gamma view too. The total allowed speed error is at most
`2e-12` times the intended **thermal crossing speed**, with no bulk-speed or
`max(1,speed)` floor. Unresolved grazing/cold events return a status without
commit; they are not repaired or resampled. This conservative numerical guard
is distinct from a formal directed-rounding interval implementation or a
guarantee for all subluminal binary64 inputs. A future explicitly rest-only
native API can retain proper momentum without requiring an unused boosted
diagnostic to be representable; the current full-output API does not silently
switch to that different contract.

The physical rate per instantaneous simulation area/time is
`n0 |a| <v'>/4`, where
`<v'>/c=2 theta(1+theta)/(exp(1/theta) K2(1/theta))`.
For `theta<=1/32`, the scaled K2 asymptotic uses
[NIST DLMF 10.40.2 and its real-argument remainder bound](https://dlmf.nist.gov/10.40).
It stops before the first omitted term after at least two terms, requiring its
relative mathematical remainder to be at most `1e-17`. That bound excludes
ordinary floating-point roundoff. Forming the series directly in `theta`
avoids inverse-temperature and exponential overflow.
An underflowing positive auxiliary remainder is rounded upward to a
representable upper bound; no physical rate/temperature/momentum is floored.

For warmer temperatures, the implementation uses the C++17 modified Bessel
function. At inverse temperature below one it forms
`z^2 exp(z) K2(z)=exp(z)[2z K1(z)+z^2 K0(z)]`; this avoids the overflowing
unscaled K2 of a very hot but finite flux. This branch has no asserted analytic
library error bound. MJ01 checks it against independently generated 90-digit
references and original positive radial integrals with relative tolerance
`2e-12`. The directly tested temperatures span `theta=1e-12` through `1e6`,
including both sides of the series transition. Special-function exceptions,
nonfinite results, a speed outside its physical bound and SI output
under/overflow return an explicit numerical failure. These checks do not claim
uniform precision throughout every representable input combination.

At zero temperature a material population has zero crossing flux; zero density
also gives zero flux. Neither has a normalized conditional positive-crossing
sample. The sampler rejects that request without draws. For positive
temperature the implemented relativistic law approaches
`n0 |a| sqrt(kB T/(2 pi m))`, with flux-weighted classical normal mean
`sqrt(pi kB T/(2m))`, normal second moment `2kB T/m` and transverse second
moment `kB T/m`. The tests compare the same sampler to this limit with explicit
finite-temperature bounds; no classical branch is implemented.

The gyroaveraged srcSEP movers do not retain gyrophase. A tilted proper cap can
receive full-particle gyro crossings that the reduced phase state cannot
represent. The present kernel therefore does not authorize projection of a
general full-flux sample onto `(p_parallel,p_perp)`. A derived guiding-center
boundary condition or an explicit properly transformed cap/B alignment domain
is required. Simulation/proper magnetic frames, density basis, material face
area, time-varying count intensity and rate-weighted birth time must also be
provided before a native source is qualified.

Run the bounded strict-warning/ASan/UBSan suite with:

```sh
SRCSEP_REPORT_DIR=/tmp/sep-mj sh srcSEP/test/run_maxwell_juttner_flux_tests.sh
```

The shared `MaxwellJuttnerFluxDescriptors()` factory supplies MJ01–06 in group
`maxwell-juttner-flux`: independent normalization, boosted moving-plane
quadrature, mixture CDF/angular moments, cold moments, atomic rejection and
keyed lifecycle. The runner additionally repeats stochastic callbacks with
different campaign/replicate seeds and reversed test order, requiring changes
in actual observables. The reference generator is explicit and read-only:

```sh
python3 -B srcSEP/test/generate_maxwell_juttner_references.py
```

No observational result, native injection test, physical source count or
production material-grid continuation follows from these kernel checks.
