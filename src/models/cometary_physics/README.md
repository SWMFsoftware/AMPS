# Generic comet collision and radiative-cooling physics

`CollisionDatabase.h` stores neutral collision data independently of comet,
geometry, and AMPS particle implementation. Application input modules resolve
species names to compile-time indices and generate immutable `cCollision`
entries.

## Lookup contract

A collision entry describes an unordered species pair. Lookup of `(a,b)` and
`(b,a)` therefore returns the same entry. Duplicate unordered pairs are
invalid.

Cross sections use m² and temperatures use K. Two laws are supported:

```text
Constant:             sigma(T) = sigma_ref
TemperaturePowerLaw:  sigma(T) = sigma_ref (T/T_ref)^alpha
```

A temperature-dependent entry also stores the closed validity interval stated
by its source. The library returns zero with `found=true` and
`temperatureInRange=false` outside that interval. The comet application treats
this as a fatal configuration error, preventing silent extrapolation.

A missing pair returns zero with `found=false`. This low-level result is for
explicit error handling; it must not be interpreted as evidence that the pair
is collisionless. The generic comet adapter stops at the first missing pair.

## Rubin et al. (2011)

`input/CometPhysicsDatabase.pm` currently provides the named `RUBIN2011`
catalog from Tables 3 and 4 of Rubin et al., Icarus 213 (2011) 655–677.
Table 3 supplies 66 unordered pairs for the 11 parent species. It does not
supply a complete daughter-species matrix, so the catalog intentionally has
no values for undocumented daughter pairs.

The H2O–H2O entry uses the Crifo relation quoted by Rubin et al.:

```text
sigma = 1.66e-19 (T/300 K)^-0.6 m²,  50 K <= T <= 300 K.
```

The remaining Table 3 parent pairs are constant hard-sphere values. The
application uses AMPS translational temperature for the H2O–H2O entry. A
future temperature-dependent mixed pair will require an explicit mixture
temperature definition before it can execute.

The application obtains that temperature with
`cDataCenterNode::GetTranslationalTemperature(species)`. The accessor averages
the completed sampling buffer and returns zero while
`PIC::LastSampleLength == 0`. Each comet input with a temperature-dependent
law must therefore supply `CollisionTemperatureFloor`. The Halley case uses
`max(T_sampled, 180 K)`, as requested; the value is case data rather than a
library default.

The generic adapter can also require a minimum completed-window mean model
particle occupancy through
`CollisionTemperatureMinimumSampleParticleNumber`. This statistic is the
same per-species value printed as `Particle Number`; a cell below the selected
cutoff uses `CollisionTemperatureFloor` instead of its statistically noisy
temperature. A zero internal default disables the filter, so a science case
must state a positive cutoff explicitly when it wants this behavior. Sampled
temperatures above a temperature-law record's cited maximum are evaluated at
that maximum rather than extrapolated; for the Crifo row this caps the
collision temperature at 300 K.

An application input can also add explicit entries. This keeps comet-specific
or newly sourced values out of C++:

```text
NeutralCollision=H2O | CO2 | CONSTANT | 3.4e-19 | source citation
NeutralCollision=H2O | H2O | TEMPERATURE_POWER | 1.66e-19 | 300 | -0.6 | 50 | 300 | Crifo (1989)
```

Explicit records may supplement a named parent database with cited daughter
pairs. Duplicate unordered pairs fail validation.

## Input-selected daughter cross-section completion

No source found in this review supplies a measured 27-by-27 elastic collision
matrix for all Rubin-network species at cold-coma temperatures. The following
traceable subsets and estimates are available:

- Combi (1996) values reproduced in Fougere's dissertation give the complete
  H2O/OH/H2/H/O/CO hard-sphere subset. Huang et al. (2016), ApJ 833:160,
  Table 2 independently prints ten of the H2O/OH/H2/H/O values.
- The standard CHEMKIN transport database supplies Lennard-Jones collision
  diameters for all 27 compiled species, including the C2 diameter of 3.621
  angstrom. These are
  transport parameters assembled for combustion modeling, and the CHEMKIN
  documentation states that some are fitted while others are estimated.
- The selected analytical hard-sphere completion uses the Lorentz diameter rule
  `d_ij=(d_i+d_j)/2` followed by `sigma_ij=pi*d_ij^2`. This is an estimator,
  not a measurement and not a choice made by Rubin et al. It also omits the
  temperature-dependent Lennard-Jones collision integral and polar effects.

The generic input keyword
`NeutralCollisionCompletion=COMBI1996_CHEMKIN_LORENTZ_HS` explicitly selects
the hierarchy. The generator retains existing entries, adds only missing
Combi rows, and fills only the pairs still absent with CHEMKIN/Lorentz
estimates. For the 27-species Halley case the resulting 378 rows comprise 66
Rubin entries, 18 additional Combi entries, and 294 estimates. Every generated
entry retains a source-tier string. Without this keyword, missing pairs remain
absent and retain the fatal lookup behavior described above.

References: Huang et al. (2016), doi:10.3847/1538-4357/833/2/160;
Fougere (2014), *The Complex Outgassing of Comets and the Resulting Coma*,
Table 2.1; Cantera `gri30.yaml`; and the Ansys Chemkin Input Manual transport
database, Table 5.1.

## Water rotational radiative cooling

`RadiativeCooling.h` contains the geometry-independent
`Gombosi1986Water` evaluator. It is the analytic approximation to the
Crovisier (1984) and Huebner (1985) H2O rotational-cooling calculations
presented by Gombosi, Nagy, and Cravens (1986), equation (48b). The positive
per-molecule loss rate is obtained from

```text
|Q| = A Trot^p nH2O exp(-tauIR)              [erg cm^-3 s^-1]
tauIR = 0.4 nH2O r sigmaIR
sigmaIR = 4.0e-15 cm^2
```

with `(A,p)=(4.4e-22,3.35)` below 52 K and
`(2.0e-20,2.47)` at and above 52 K. The function accepts K, m^-3, and m and
returns J s^-1 molecule^-1. Species selection, internal-energy storage, and
mesh traversal remain in the application adapter.

Rubin et al. (2011) state that their Halley model uses this analytic
approximation, but do not reprint it. `halley1d.input` therefore names both
the affected species and the exact model:

```text
RadiativeCooling=H2O | GOMBOSI1986_H2O
```

The adapter reads H2O density and rotational temperature from the most recent
completed sampling window and removes energy from each H2O particle's AMPS
Larsen-Borgnakke rotational-energy state after every particle iteration.
Cooling is zero until the first complete sample exists.

Run the tests with:

```sh
perl input/test_comet_physics_database.pl
g++ -std=c++17 -Wall -Wextra -pedantic \
  -Isrc/models/cometary_physics \
  src/models/cometary_physics/test_collision_database.cpp \
  -o /tmp/test_collision_database
/tmp/test_collision_database

g++ -std=c++17 -Wall -Wextra -pedantic \
  -Isrc/models/cometary_physics \
  src/models/cometary_physics/test_radiative_cooling.cpp \
  -o /tmp/test_radiative_cooling
/tmp/test_radiative_cooling
```
