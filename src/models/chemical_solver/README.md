# Reaction-network data model

`ReactionNetwork.h` is the geometry-independent chemical-network container
used by the generic comet application. It does not own particles, choose a
time integrator, or resolve species names. The input preprocessor performs
name-to-index resolution and emits static tables; a particle or continuum
solver then decides how to execute those tables.

## Data and units

Each reaction stores slices of one flat stoichiometric-term array. A term is
an AMPS species index plus a positive coefficient. A reaction may have no
products; this represents loss to states omitted from the compiled species
list, such as photoionization in a neutral-only simulation.

The caller must use one consistent number-density system for mass-action
rates. The coefficient `A` consequently has units appropriate to reaction
molecularity. Temperature and activation temperature are in K,
heliocentric distances are in AU, and `ExcessEnergyJ` is in joules.

The implemented rate laws are:

```text
CONSTANT:      k = A
HELIOCENTRIC:  k = A (r_ref / r_h)^2
ARRHENIUS:     k = A T^beta exp(-T_a/T)
```

`EvaluateSourceLoss` calculates nonnegative source and loss arrays from the
mass-action event rate. It deliberately leaves time integration to the
solver. `IsUnitUnimolecular` identifies reactions that a one-particle decay
callback can execute without selecting additional reactant particles.

## Unit-unimolecular execution index

`cReactionNetwork` builds a compact, per-parent index of unit-unimolecular
reactions when the network object is constructed. The public
`GetUnitUnimolecularReactionCount(species)` and
`GetUnitUnimolecularReactionIndex(species,channel)` accessors expose only
reactions that `IsUnitUnimolecular` can execute. Channels for each parent keep
their original input-table order, so an application can replace a full-table
scan without changing mutually exclusive channel ordering or random-number
selection.

The index is an execution aid, not a restriction on the generic database.
Binary and higher-order reactions remain in the main table and continue to be
available through the general APIs. Index construction skips malformed table
records so that `Validate()` remains responsible for returning the detailed
diagnostic. Applications must validate the network before using the index.

The comet particle executor uses this index to build one immutable rate table
per parent species at a stage boundary. It therefore does not scan every
reaction or allocate temporary channel/rate vectors for every model particle.
If a runtime parameter used by a rate law changes, the application must rebuild
its rate cache before the next particle batch; the generic network does not
cache application state.

## Geometry-independent event diagnostics

`ChemistryDiagnostics.h` provides an optional accumulator for particle or
continuum executors. It contains no PIC, MPI, mesh, or comet dependency. A
caller configures the number of product terms for every reaction and records
local eligible visits, selected reactions/channels, destruction and sinks,
allocation pressure, product-carrier multiplicity, slot allocation/reuse, and
represented number, mass, momentum, and energy.

Product records are indexed by reaction-table term rather than distinct
species. This preserves repeated products and allows an audit to trace a
counter back to the generated input record. `CountValues()` and
`MomentValues()` expose stable flat arrays for separate integer and
floating-point MPI sum reductions; MPI itself remains an application/runtime
responsibility. Configuring the object disabled allocates no reaction table,
every record method becomes a no-op before index validation, and both output
arrays are empty.

## Input generation

`input/CometReactionNetwork.pm` parses the application record

```text
ChemistryReaction = ID | REACTANTS | PRODUCTS | LAW | A | BETA | TA_K | RREF_AU | EXCESS_EV
```

Reaction sides contain `+` separated terms. A term is `SPECIES` or
`COEFFICIENT*SPECIES`. `PRODUCTS=NONE` is an explicit sink. For example:

```text
ChemistryReaction = H2O_PHOTO | H2O | OH + H | HELIOCENTRIC | 1.0e-5 | 0 | 0 | 1 | 1.7
ChemistryReaction = H2O_ION | H2O | NONE | HELIOCENTRIC | 4.04e-7 | 0 | 0 | 1 | 0
```

The generator converts eV to joules and writes
`CometChemistry.generated.h` in the AMPS working source tree. Unknown species,
duplicate reaction identifiers, invalid numbers, and malformed sides stop
preprocessing.

## Validation and tests

`cReactionNetwork::Validate` checks array ranges, species indices,
coefficients, rate parameters, and excess energy. Application-level
validation may be stricter. The comet particle executor, for example,
requires one unit reactant, integer product multiplicities, neutral mass
balance, and individual particle-weight correction.

The standalone test also checks the parent index, preservation of channel
order, exclusion of a binary reaction, and invalid species/channel sentinels.
The comet Stage 0 test additionally checks the diagnostic accumulator with
one-, two-, repeated-, and zero-product reactions and identical reductions at
1, 4, and 12 MPI ranks.

`ChemistryDiagnostics.h` is deliberately independent of PIC geometry and MPI.
The application configures one record per network reaction and one subrecord
per product term, accumulates integer event/carrier counts plus expected and
represented number, mass, momentum, and energy, and reduces the stable
`CountValues()` and `MomentValues()` layouts collectively. Represented-number
variance is stored as a moment so the application can evaluate a statistical
production tolerance after reduction. `Reset()` clears one sampling window
without changing the configured layout. When configured disabled, no reaction
table is retained and hot-path record methods are inert.

Run the standalone tests with:

```sh
perl input/test_comet_reaction_network.pl
g++ -std=c++17 -Wall -Wextra -pedantic \
  -Isrc/models/chemical_solver \
  src/models/chemical_solver/test_reaction_network.cpp \
  -o /tmp/test_reaction_network
/tmp/test_reaction_network

srcHartley2/tests/run_stage0_tests.sh
```

