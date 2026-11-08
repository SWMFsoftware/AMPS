# Field-line kernels

This suite qualifies the linear interpolation now used by production segment
getters, piecewise-curved arc interpolation, unique indexed payload exchange,
exact-intrinsic particle splitting, and conservative overlap remapping when a
line is repartitioned. Remote MPI transport remains a T10 gate; T8 checks the
packing/ownership and conservation algebra on one rank.
