# Field-line kernels

This suite qualifies the linear interpolation now used by production segment
getters, piecewise-curved arc interpolation, unique indexed payload exchange,
exact-intrinsic particle splitting, and conservative overlap remapping when a
line is repartitioned. Remote MPI transport remains a T10 gate; T8 checks the
packing/ownership and conservation algebra on one rank.

The split and invalid-input cases also exercise the exact clone kernel used by
`cFieldLineSegment::WeightedParticleSplitting`. They check relativistic kinetic
energy and unchanged phase space, opaque-history copying, the actual
parent/child callback arguments, equal-weight selection under permuted particle
handles, duplicate/zero identity rejection, subnormal-weight skipping, and an
application preflight rejection before allocation or weight changes. The
optional `FledLine::ParticleIdentity`, `PrepareExactSplit`, and `PostExactSplit`
hooks are registered before stepping. The key supplies the stable tie break;
preflight validates application lineage/state without consuming randomness;
post-split initializes declared application metadata after the exact clone.

These dependency-light records do not qualify PIC buffer storage, native MPI
transfer, or an application's child RNG/hazard policy. `CloneParticle` remains
byte-exact for migration and restart. Native application tests must separately
verify registration and complete stochastic lifecycle behavior. The field-line
weighted merger has a different conservation contract and is not qualified by
these split checks.
