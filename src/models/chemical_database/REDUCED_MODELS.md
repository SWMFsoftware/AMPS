# Explicit-parameter reduced models: implementation and test contract

Written before implementation, 2026-11-19. This is application-independent
software implementing the user-approved reduced-model forms. It does not
authorize a database source/state, fill missing observations, select a species,
material or spectrum, or enable any application reaction. The source packet,
profile, energy zero, parameter units and physical approximation remain caller
responsibilities. Existing Runtime ABI and model defaults are untouched.

The header-only `AMPS::ChemicalDatabase::ReducedModels` functions return the
existing typed `Status`, with unchanged outputs on every refusal. No function
uses RNG, network, PIC, SPICE, files or implicit species parameters. All inputs
are finite SI scalars/vectors except dimensionless quantiles/branch fractions.
Explicit arithmetic conservation/direction bounds are caller inputs, not
physical defaults or acceptance criteria.

| Function | Declared mathematical model and domain |
|---|---|
| `IonPhysicalMass` | `m_ion=m_neutral-m_e+IE_J/c^2`; positive masses/c, nonnegative IE, neutral mass exceeds electron mass. Returns additive reference and binding contributions separately. No uncertainty/covariance is invented. |
| `RadiativeBranches` | A single supplied coherent A-set in s^-1; total `sum(A)`, lifetime `1/sum(A)`, branch `A_i/sum(A)`. Nonnegative finite A, positive total. Zero branches retained; unavailable A is not zero. |
| `OpalSecondaryEnergy` | Lower-energy-electron marginal only: `F(w)=atan(w/W)/atan(wmax/W)`, `wmax=(E_inc-I)/2`, supplied positive W and E_inc>I, quantile in [0,1]. Energy SI J. Does not define the correlated daughter angles, ion recoil or a full electron-impact event. |
| `TwoBodyProducts` | Supplied total incoming momentum P, incoming translational/available bath energy K and explicit net Q in the same frame/reference. `K_COM=K+Q-P^2/(2(m1+m2))`; products share COM velocity plus opposite relative momenta of magnitude `sqrt(2*mu*K_COM)` along supplied normalized direction. Explicit arithmetic bounds guard momentum/energy before publication. No threshold/internal-state mapping is inferred. |
| `InheritParentVelocity` | Heavy daughter takes declared parent velocity. Outgoing signed reservoir momentum/energy/charge is the accounting difference of parent plus supplied incoming bath and Q minus heavy product. Parent/heavy/bath charges are explicit C; no electron charge is assumed. A negative reservoir exchange is permitted and explicit; this is a reduced closure, not a measured joint distribution or a closed plasma bath. |
| `ThermalRelease` | First-order fixed-temperature/site Polanyi-Wigner rate `nu exp(-E_bind/(kB*T))`, exact constant-rate finite-inventory loss `N0[1-exp(-rate*dt)]`. nu, binding energy, kB, T, N0 and dt are all explicit. No replenishment, temperature evolution, site mixture or material coefficient is supplied. Log-domain arithmetic preserves exposure when a small rate times a large duration is representable; representability flags distinguish numerical underflow from a physical zero. |

Tests compile the actual header with C++17 and strict warnings. Labelled
analytical fixtures are mechanics, not observations. Independent Python
Decimal/Fraction oracles check ion-mass arithmetic, coherent lifetime/branches,
Opal CDF inversion, moving-COM two-body conservation/Q, signed reservoir budgets,
and integrated finite inventory across ordinary/adversarial scales. Tests also
exercise NaN/infinity/null output, unavailable/negative domains, explicit
conservation bounds, no mutation on refusal and genuine source-independent
relocation/compilation. Constants in the mathematical fixtures are explicit
test data, not default species/source values. All status/count/log/compile/source
identities are retained. A software PASS never becomes physical VERIFIED or
VALIDATED; current formal M3 gaps remain independent.

Run `python3 -B test/chemical_database/test_reduced_models.py --output-dir NEW`.
No application configuration or full AMPS build is needed. The shipped SDK
includes this header; any application can compile it and pass its own admitted
source/profile inputs. Model definitions and policy metadata are separately
selected by the explicit shared-database publication option.
