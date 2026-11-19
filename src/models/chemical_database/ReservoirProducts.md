# Heavy-product/reservoir moment adapter

`ReservoirProducts.h` is a C++17 header-only service in
`AMPS::ChemicalDatabase::ReservoirProducts`. It has no Moon, PIC, MPI, RNG, input
file, network or database-installation dependency. Link the existing shared
DataIO/Runtime objects when using `Status`/`CodeName`. It reuses the production
`ReducedModels::InheritParentVelocity` and `TwoBodyProducts` services; no second
rate or source formula is introduced. The preimplementation contract and
independent-oracle procedure are in
`srcMoon/test/M3_RESERVOIR_PRODUCTS_CONTRACT.md`.

## Explicit inputs and output meaning

Supply a `State` for every material participant: source/state identities,
positive isotope nuclei, integer elementary charge, kinetic/reference material
mass in kg, independently sourced physical rest mass in kg, and state energy in
J on the SAME named reference. Binding energy already counted in state Q must
not also be counted in kinetic masses. The reference material masses must close
within caller-supplied arithmetic bounds. This explicitly refuses an
incompatible compiled-species mass convention; it does not silently adjust it.

`Evidence` names pinned source, state-crosswalk, source domain, model-error and
observable documents. Callers authenticate those documents and source rights;
a name/hash string is not a scientific authorization certificate. Missing proof
or error limits returns `Unavailable`. The helper selects no species constants,
reference state, source fit, scattering support or acceptance threshold.
`EventContext` carries actual position, parent velocity, physical start/elapsed/
remaining/event time and reacted represented number. One named inertial frame
and time scale applies to every material term. No velocity rotation, source
clock, timestep ratio or new particle weight is inferred.

Both composition entry points return `Result` containing separate parent,
tracked product, incoming bath and outgoing bath `Moment` snapshots. Each moment
retains its complete source state, integer count/isotope composition, SI
momentum, kinetic/state energy and charge. Net bath changes are OUT minus IN;
represented changes multiply the SAME reacted number. Source prescription,
state Q, numerical residuals and local model discrepancy remain distinct.
There is no PIC commit, event draw or automatic physical activation.

## Prescribed-Maxwellian electron impact

`MaxwellianCollisionMeanEnergy(input,&mean)` evaluates

```
mean incoming collision energy = kBT * (3/2 + d ln k / d ln T)
```

This follows differentiation of the Maxwellian rate integral when the cross
section is T-independent and the target-rest distribution is unshifted/isotropic.
The caller supplies the actual source derivative at the requested temperature,
its supported thermal-energy interval, threshold and explicit assumptions. An
actual nonzero parent-rest bulk vector is refused. It is not legitimate to
replace a missing collision-conditioned mean by the unconditional `3*kBT/2`.
A derivative of a fitted scalar rate defines a MODEL moment, not an independent
measurement, an inverted cross section, or an outgoing electron spectrum.

`ElectronImpactInheritedVelocity(input,&result)` implements the approved reduced
`e+A -> A+ +2e` moment closure. The ion inherits the ACTUAL parent velocity.
The incoming one-electron mean is constructed from the supplied conditional
energy and the parent-to-inertial boost. The outgoing TWO-electron aggregate
uses the explicit state-Q-constrained closure, with a nonnegative realizable
second moment. Individual electron energies/directions are not sampled. The
outgoing and incoming aggregates remain separately inspectable before their net
exchange is formed. Rotationally invariant collision weighting is an additional
explicit requirement for a zero incoming first moment in the parent frame.

The local omitted ion-recoil RMS speed bound is

```
[sqrt(2*me*Ein) + sqrt(4*me*(Ein+Q))] / mIon
```

in the declared kinetic-mass model. It follows Cauchy/Minkowski and the outgoing
energy upper budget, without assuming electron correlations. The corresponding
mean absolute lab kinetic discrepancy is bounded using the actual parent speed.
This is an ENSEMBLE moment bound, not an all-event support/confidence bound.
Physical-versus-kinetic mass discrepancy is separately reported and capped.
No downstream trajectory, density, LOS or pickup-energy error follows without a
qualified transport-amplification/observation analysis.

## Forward partial charge exchange

`ForwardPartialChargeExchange(input,&result)` admits a neutral nucleus-carrying
target and singly charged nucleus-carrying projectile, with one electron
transferred. Each carrier retains its exact isotope composition. It requires an
explicit partial-channel identity, independently resolved state Q, source-domain
mapping to incoming COM-relative J, and explicit supported COM angular-error
bound. The actual incoming plasma velocity is retained. Undefined zero relative
directions and outside-domain energies are refused.

The shared conservative two-body solver sets the heavy-product COM direction
parallel to its incoming COM direction and retains Q in relative kinetic energy.
The outgoing fast neutral is an explicit bath product, and the consumed plasma
ion is an explicit incoming bath term. This is not unconditional velocity-label
exchange and is not an electron-only sink. A source-supported maximum deflection
`theta` gives a local heavy velocity discrepancy bound
`2*pRelative*sin(theta/2)/mHeavy`; the returned energy bound uses the actual
forward product velocity. An explicit full-pi mathematical bound is possible,
but never means the source measured that angular distribution.

## Failures and verification

All outputs remain unchanged on refusal. `Unavailable` means required source,
state, Q, derivative, domain or error proof is missing. `InvalidInput` means
contradictory identities, malformed hashes, integer topology/material convention
or invalid arithmetic inputs. `OutsideDomain` includes unsupported drift/reaction,
source-domain miss, impossible product energy or exceeded declared model budget.
`NumericalFailure` includes overflow, erased represented moments or arithmetic
conservation miss; `ResourceLimit` identifies allocation failure. None changes a
source value or relaxes a physical/model bound. `Evaluated` is mathematical
software success, never physical qualification.

Run the actual shared-kernel thin probe and independent Python oracles with:

```sh
python3 -B test/chemical_database/test_reservoir_products.py \
  --output-dir test_output/chemical_database/NEW-reservoir-test
```

The script preserves compiler/source/executable identities, commands and native
outputs and emits `result.json` plus exit0/1/2/77. Shared unittest discovery can
load the same methods without command-line side effects. The explicitly nominal
fixtures are not Na/He measurements, observations, a source-selection approval
or a replacement for formal I23/I24 and predecessor gates.
