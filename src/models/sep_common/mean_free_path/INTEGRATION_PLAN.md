# srcSEP3D focused-transport integration plan

This plan covers future host work only. The standalone library does not change
`srcSEP3D` configuration, movers, restart metadata, or production linkage.

## 1. Link and ownership

Build `libmean_free_path.a` once from this directory and flatten or link its
single object through the same audited mechanism used for the parallel
diffusion library. Keep all physical formulas here. The application adapter
may translate runtime state but must not copy equations. Extend archive/member
tests so duplicate strong definitions and PIC/MPI dependencies fail.

## 2. Input grammar

Add a focused-transport-only section, provisionally:

```ini
[mean_free_path]
model = PA-QFORM
lambda_parallel_m = 4.487936121e10
operator_convention = standard
q = 5/3
gap_parameter = 0.05
```

The actual section name and schema version must be selected in the srcSEP3D
configuration contract before implementation. The application parser should
retain file/line information, pass exact key/value text to
`BuildConfiguration`, and surface its typed status. It must not pre-parse
rationals, add units, choose a species, or supply defaults. Parker transport
continues to use the parallel-diffusion selector; the mean-free-path selector
is admitted only for focused transport.

## 3. Runtime adapter

At serial initialization, install the validated configuration, freeze its
fingerprint and model/runtime-state identity, and reject restart mismatch.
For each particle evaluation construct:

- `ParticleState` from the compiled species mass, signed charge, total
  momentum, and explicit nucleon count where the selected variable requires
  it;
- `LocalState` from one coherent background/turbulence snapshot, including
  field magnitude, field angle, slab quantities and provider revisions only
  when available;
- the shock side and signed normal coordinate from the reduced-front provider
  only for a selected shock-region closure; and
- a sample identity containing the background/turbulence generation.

Total turbulence variance must not be relabeled slab variance. Existing
correlation lengths must not be relabeled bend-over or slab-correlation
lengths. Missing values remain unavailable and produce the library's typed
status.

## 4. Focused mover contract

The focused mover should request `D_mumu` through
`EvaluatePitchAngleDiffusion`. For the standard operator
`d_mu(D_mumu d_mu f)`, its Ito increment must use drift
`partial_mu D_mumu` and noise `sqrt(2 D_mumu)`. For HalfD, the effective
coefficient is half the printed coefficient in both terms. The numerical
treatment must:

- maintain `-1 <= mu <= 1` and no-flux endpoints;
- explicitly impose no flux at `mu=-1` for the printed `(1-mu)^2` variant;
- document and converge the singular/discontinuous drift near `mu=0` for the
  q and epsilon forms;
- never insert an undocumented positive gap parameter; and
- treat isotropic direction-reset models as Poisson events rather than a
  pitch-angle SDE.

Focusing uses the host's signed focusing length and may apply
`MFP-FOCUS-HW13` only when explicitly selected. The library's stable scalar
helper does not choose a field model.

## 5. Geometry and perpendicular transport

Keep the particle mover in a field-aligned basis. Convert a radial SEP MFP to
parallel only with the local field angle. If a radial tensor quantity is
needed, request the existing perpendicular coefficient and apply Eq. (3);
never identify Eq. (2) with the complete tensor projection. The current
ownership of perpendicular diffusion remains unchanged.

## 6. Provider and time ordering

Freeze a coherent background/turbulence snapshot, evaluate the coefficient,
publish its identity, and only then advance particles. Any interpolation
policy belongs to the provider scheduler and must be configured. Coefficient
cache keys, if introduced, must contain particle/species state, local snapshot
and provider generations, configuration fingerprint, shock side, and formula
variant. Begin without a cache.

## 7. Failure and diagnostics

Map library statuses to srcSEP3D typed configuration/runtime failures. A source
gate, missing provider field, out-of-domain input, and numerical failure must
remain distinct. Do not replace any of them with ballistic propagation or a
different closure unless a separately documented host policy is later
approved. Route fatal initialization failures through the repository's
debugger-interceptable `exit(__LINE__, __FILE__, message)` trap.

## 8. Qualification sequence

1. Parser/component tests for every admitted model and every rejected key.
2. Restart fingerprint match/mismatch tests.
3. Fixed-background local agreement between standalone and srcSEP3D adapters.
4. Deterministic PDE/operator tests for isotropic scattering moments.
5. Stochastic convergence for mean and second moment, including endpoints.
6. One- and four-rank decomposition-independence tests with identical
   coefficient/sample identities.
7. Frozen-snapshot invariance, then time-dependent provider tests.
8. End-to-end focused-transport cases only where an independent reference
   exists; report algebraic, numerical, coupling and science evidence
   separately.

No native or MPI qualification is claimed by the standalone tests.
