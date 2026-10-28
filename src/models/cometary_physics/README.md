# Generic comet collision database

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

An application input can also add explicit entries. This keeps comet-specific
or newly sourced values out of C++:

```text
NeutralCollision=H2O | CO2 | CONSTANT | 3.4e-19 | source citation
NeutralCollision=H2O | H2O | TEMPERATURE_POWER | 1.66e-19 | 300 | -0.6 | 50 | 300 | Crifo (1989)
```

Explicit records may supplement a named parent database with cited daughter
pairs. Duplicate unordered pairs fail validation.

Run the tests with:

```sh
perl input/test_comet_physics_database.pl
g++ -std=c++17 -Wall -Wextra -pedantic \
  -Isrc/models/cometary_physics \
  src/models/cometary_physics/test_collision_database.cpp \
  -o /tmp/test_collision_database
/tmp/test_collision_database
```
