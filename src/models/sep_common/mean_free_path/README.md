# Mean-free-path models for SEP and GCR transport

This directory contains the host-neutral C++17 implementation associated with
`MEAN_FREE_PATH_MODEL.md`, specification version 1.2. The library is intended
for future focused transport in `srcSEP3D`, but this delivery is deliberately
standalone: it does not modify an application parser, mover, mesh, MPI state,
turbulence provider, shock provider, or restart format. See
`INTEGRATION_PLAN.md` for that later binding.

All public dimensional values are SI. Momentum is total relativistic momentum
in kg m s^-1, rigidity is in volts, kinetic energies exposed by `Kinematics`
are in eV, length is in metres, magnetic field is in tesla, and diffusion
coefficients are in m^2 s^-1. Configuration numbers are suffix-free. The
companion-data loader retains text such as `0.3 AU` as a number plus a unit; it
never feeds that text directly to an SI evaluator.

## Delivered boundary

The library provides:

- exact-relativistic species kinematics without a default species;
- tagged parallel, SEP-radial, radial-tensor, isotropic-scattering and
  unspecified mean free paths, with explicit Eqs. (2) and (3) conversions;
- a strict parser-neutral manager, transactional configuration, a public
  active evaluator pointer, scalar/batch evaluation and provenance;
- general SEP Eq. (9) prescriptions under their source IDs;
- the confirmed `SEP-CHEN24` mapping and separately named code-derived
  `LEGACY-TENISHEV2005AIAA` law;
- q, epsilon, isotropic, printed q/epsilon, EPREM HalfD, Dröge-VA, PA-KOLMO,
  and stable He-Wan focusing calculations;
- TS2003 ion/electron and Zank turbulence-based forms;
- Bohm, Afanasiev and M-FLAMPA shock closures; and
- NWU, Corti, HelMod, Strauss, Effenberger, Wang, Tomassetti, Perugia and
  printed-Duan GCR forms.

Published source IDs install no hidden calibration. A caller supplies the
specific published value or explicitly selected member of a published range.
This is necessary where sources publish several event, sector, year, or
sensitivity values. Exact input text remains in `Configuration` for methods
and restart records.

Known but unavailable identities return typed statuses instead of another
model:

| Identifier | Reason |
| --- | --- |
| `SEP-LAITINEN16/18`, `GCR-QINSHEN17` | Need a configured pinned PARALLEL QLT/NLGCE-F dependency. |
| `SEP-MINOSHIMA26` | D-29: `Omega_n` is undefined. |
| `SHOCK-PARASOL` | U-6: `Lambda(E)` and `Delta x(E)` are absent. |
| `GCR-LUO19` | Reference data only; the formula was not legible. |
| `GCR-BOBIK12` | Normalization units are absent. |
| `GCR-JIANG23` | The printed `K0` unit is unresolved. |
| `PA-AMPS-I/II/IV` | The source/OCR record is insufficient for a source-exact normalized evaluator. |
| `PA-AMPS-III/V` | Need the pinned spectrum/QLT dependency. |

## Parser-neutral manager

Include `mean_free_path.h` and use `SEP::MeanFreePath`.
`BuildConfiguration(model, assignments, &configuration)` owns semantic
validation. The application retains file syntax and line diagnostics. Each
model has a closed key set: unknown, duplicate, missing, nonfinite,
unit-suffixed or inactive keys fail. Exact rational text such as `1/3` is
accepted only for dimensionless fields.

`SetActiveConfiguration` and `ConfigureActiveModel` install the validated
configuration together with:

```cpp
extern ModelFunction ActiveModelFunction;
```

Configure during serial initialization; concurrent reconfiguration is outside
the contract. Direct `Evaluate` calls are reentrant. `EvaluateBatch` checks
array shape before writing output.

Example reproducing Liu et al. (2025)'s central ambient prescription:

```cpp
Status status = ConfigureActiveModel(
    "SEP-MFLAMPA25",
    {{"lambda0_m", "44879361210"},
     {"lambda_kind", "parallel"},
     {"momentum_variable", "momentum_pc_ev"},
     {"x0_SI", "1e9"},
     {"momentum_exponent", "1/3"},
     {"radius0_m", "149597870700"},
     {"radial_exponent", "1"}});
```

Those values demonstrate source-explicit syntax; they are not defaults.

## Model schemas

### SEP Eq. (9)

Executable `SEP-*` power laws require `lambda0_m`, `lambda_kind`,
`momentum_variable`, `x0_SI`, `momentum_exponent`, `radius0_m`, and
`radial_exponent`. Choices are:

- `lambda_kind`: `parallel`, `radial_sep`, or `isotropic_scattering`;
  `unspecified` is rejected under U-11, while `radial_tensor` is rejected
  because Eq. (3) must be assembled explicitly with perpendicular transport;
- `momentum_variable`: `rigidity_v`, `momentum_pc_ev`, `kinetic_total_ev`, or
  `kinetic_per_nucleon_ev`.

The per-nucleon variable requires `ParticleState::nucleonCount`. It is never
inferred. `SEP-EPREM10` is refused away from 1 AU because its radial exponent
sign is unreadable.

### Confirmed legacy mappings

`SEP-CHEN24` requires `output_quantity = kappa_parallel|lambda_parallel` and
`domain_policy = error|warn_and_evaluate`. Lambda output also requires a
`species_id` matching the particle. The published output is

`kappa = 5.16e14 (r/AU)^1.17 (E/keV)^0.71 m^2 s^-1`.

Lambda is marked derived. This is the user-confirmed interpretation of the
host label `Chen2024AA`.

`LEGACY-TENISHEV2005AIAA` requires `lambda0_m`, `kinetic_energy0_J`,
`energy_exponent`, `radius0_m`, and `radial_exponent`, preserving the audited
host law `lambda=lambda0(T/T0)^alpha(r/r0)^beta`. It is labeled code-derived;
the library does not attribute the equation to an unverified publication.

### Pitch-angle models

Every pitch model requires `lambda_parallel_m` and
`operator_convention = standard|half_d`.

| Model | Additional keys |
| --- | --- |
| `PA-QFORM`, `PA-QFORM-LANG-PRINTED` | `q`, `gap_parameter` (`H`) |
| `PA-EPS`, `PA-EPS-PACHECO-PRINTED` | `gap_parameter` (`epsilon`) |
| `PA-DROGE-VA` | `q`, `alfven_to_particle_speed` |
| `PA-ISO` | none |
| `PA-EPREM` | `lambda_interpretation = printed_parameter|transport_lambda`; must use `half_d` |
| `PA-KOLMO` | explicit `q=5/3`, `gap_parameter=0`; other values are rejected |

Dröge input lambda is the nominal Eq. (16) prefactor; the returned lambda is
the Eq. (11) integral value. EPREM can reproduce either its printed parameter
or a requested transport MFP, but only after that U-8 interpretation is named.
Printed formula variants have distinct IDs and are never defaults.

### Turbulence and shock models

- `QLT-TS03-P`/`GCR-EB13`: `inertial_index`; runtime field, total slab
  variance and `k_min`.
- `QLT-TS03-E-RS`: also `dissipation_index`, `alpha_D`; runtime `k_d`, `V_A`.
- `QLT-TS03-E-DT`: also mandatory
  `dt_variant = ts_q3ms|eb13|lang24`; there is no default.
- `QLT-ZANK98`: mandatory
  `variance_convention = per_component|total`; runtime field, variance and
  slab correlation length.
- `SHOCK-BOHM`: local mean field and an explicit shock side.
- `SHOCK-AFANASIEV15`: `x0_m` and upstream signed distance, shock-frame inflow
  and Alfvén speed; `(u1-VA)(x+x0)>0` is enforced.
- `SHOCK-MFLAMPA`: `Lmax_over_radius`, `kappa_floor_m2_per_s`, and downstream
  radius, field and wave variance. The floor is diagnostic, not a Bohm bound.

The TS evaluator retains the published sum of asymptotes. No turbulence value,
slab fraction, scale conversion or dissipation onset is inferred.

### GCR models

Field-dependent GCR schemas require
`field_normalization = magnitude|radial_amplitude`, explicitly resolving U-13.
The local value must use that declared convention; the library does not build
a Parker field.

| Model | Required numerical keys beyond the field convention |
| --- | --- |
| `GCR-NWU14` | `K0_m2_per_s`, field/rigidity references, break, slopes, smoothness |
| `GCR-CORTI19` | `K0_m2_per_s`, field reference, break, slopes, smoothness |
| `GCR-HELMOD17` | `K0_AU2_per_s_numeric`, `g_low`, `normalization_convention=d21_numeric_use` |
| `GCR-HELMOD19` | HelMod17 keys plus `R_c` |
| `GCR-STRAUSS11` | `lambda0_m`, `rigidity0_V`, `radius0_m` |
| `GCR-EFFENBERGER12` | `K0_m2_per_s`, `momentum_pc0_eV`, field reference/exponent |
| `GCR-WANG19` | base scale, externally supplied `B_c_T`, mean Earth field, activity exponent |
| `GCR-TOMASSETTI17` | station-fit `a,b`, externally supplied modulation potential, field reference |
| `GCR-PERUGIA21/25` | `K0`, references, break, slopes, smoothness; time-dependent values supplied externally |
| `GCR-DUAN25` | `K0`, equatorial field, break, `a,b,c`, `formula_variant=duan25_printed_d33` |

Every `K0_m2_per_s` is a selected VALUE. A number printed only as a UNIT is not
accepted as a value. HelMod keeps the D-21 numerical-use convention visible.
Duan keeps the printed D-33 law even for `b<a`; it is not rewritten to agree
with the prose.

## Companion data and provenance

`LoadSourceRecord(bundle, relative_file, id, &record)` reads an identified JSON
object without modifying it, confines paths to the bundle and attaches the
file's recorded SHA-256 from `SHA256SUMS`. The independent Python verifier
recomputes the digest; the C++ loader does not duplicate cryptographic code.

`ParsePublishedScalar` accepts a plain number, an exact rational in a declared
dimensionless field, or a single number plus unit. It refuses ranges,
uncertainties, inequalities and unresolved expressions. This is an audit
parser; application input remains suffix-free SI.

Results record model ID, source key/location, configuration fingerprint,
runtime state and local sample identity. The FNV-1a fingerprint is a restart
identity, not a cryptographic digest.

## Build and tests

From this directory:

```sh
make clean
make
make test
make verify
```

Run only the C++ suite with:

```sh
build/test_mean_free_path --bundle MEAN_FREE_PATH_MODEL_DATA
```

The supplied independent checks require Python 3.9+ and `mpmath`:

```sh
python3 MEAN_FREE_PATH_MODEL_DATA/reference_verification.py
python3 MEAN_FREE_PATH_MODEL_DATA/independent_physics_checks.py
```

`make test` runs the C++ fixtures and both scripts. `make verify` additionally
checks that the archive has no PIC or MPI dependency. The C++ suite covers
parser rejection, transactional dispatch, kinematics, tags/geometry, SEP and
Chen fixtures, pitch normalization and operator behavior, stable focusing,
turbulence, shock, GCR, batch shape and source-record provenance.

These are coefficient-library tests. They do not validate a focused mover,
SEP event, GCR modulation, shock acceleration, MPI decomposition or evolving
provider coupling.

## Known limitations

- No srcSEP3D code consumes this library yet.
- Published rows remain in the companion archive and are selected explicitly;
  there is not yet a convenience API converting every event/year row into a
  configuration automatically.
- The C++ source loader addresses JSON records carrying `id`; CSV rows remain
  available to the independent verifier and future selection adapters.
- Wavenumber-resolved QLT/NLGCE delegation is typed but not linked here.
- Registered source gates remain unavailable; no alternate physics replaces
  them.
- Unit-test success is mathematical evidence, not observational or end-to-end
  transport validation.
