# `srcMoon` species-specific production modules

This directory contains lunar physics that belongs to one simulated species
and should not be embedded in the generic Exosphere library. Generic terrain,
thermal, inventory, mover, and chemistry services remain in their existing
modules; moving them here would incorrectly make shared physics species-local.

`SodiumCase1.{h,cpp}` is the narrowly scoped production adapter for the exact
Tenishev et al. (2013) Case-1 demonstration. It supplies input-validated MIV
rate/temperature and photoionization lifetime to the generic production
kernels, initializes and maintains the prescribed uniform surface reservoir,
and verifies the generated PSD/TD/source switches before a run. It does not
duplicate the PSD distribution, force law, impact handler, chemistry, or
samplers. Inputs without `sodium_exosphere_model=tenishev-2013-case1` do not
install any of its callbacks.

`Argon.cpp` implements the optional continuous radiogenic `40Ar` source. Its
pure parameter validator and uniform reference-sphere area partition are used
by U21. Its production callbacks use the normal AMPS user-defined-source
dispatcher, source accounting, MPI ownership, body/LSO state transform, local
surface-temperature callback, and resolved-LOLA injection adapter. It does not
duplicate a test-only injection law.

The canonical `input/moon.input` keeps `RadiogenicArgon` off because the frozen
M0/M1 executable contains Na only. An experimental Ar run requires all of:

1. `SpeciesList=Na,Ar` at configuration time;
2. `SourceProcessCode=RadiogenicArgon` set to `on` while historical `MySource`
   remains off;
3. a clean Config.pl build with SPICE;
4. `argon_source_enabled=true` in a copied runtime `.in` receipt.

For I18-I20, `srcMoon/test/build_m2_qualification.py` performs the first three
steps in a uniquely named temporary Config.pl application, hashes the inputs
and executable, verifies that generated photolytic chemistry is OFF for the
approved I19 budget, and removes every application-specific generated file.
The aggregate runner then restores the canonical Na-only build. Its linked
runner derives the fourth setting in the preserved per-test `application.in`;
neither source input file is edited in place.

The M2 build also requires individual particle weights ON to represent a final
fractional finite-release packet. The shared Exosphere finite-inventory
implementation owns depletion and departure sampling; Argon retains its
ordinary physical rate and Maxwell-flux source law. During resolved moving
source injection, early spherical MPI ownership is deferred until the final
terrain position at the sampled departure epoch. An intermediate spherical
point can have a different owner; discarding it early would lose valid atoms.
The eight-rank supporting I18/I19 checks exercise this production path, while
formal qualification and temporal/weight convergence remain open.

The runtime parser then accepts only a continuous uniform-global source and
fails if Ar is absent from the compiled species table. Potassium-weighted
geography is deliberately unavailable while D10 is regional and QA-incomplete.
See [`../SURFACE_INTERACTION_LITERATURE.md`](../SURFACE_INTERACTION_LITERATURE.md)
for the rate, residence energy, sticking interpretation, citations, and data
limitations.
