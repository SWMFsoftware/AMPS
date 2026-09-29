# Stage 8: runtime, mesh, and observer integration contracts

Stage 8 supplies dependency-free kernels used by the `srcSEP3D` AMPS adapter.
The shared library never owns an AMPS cell or particle.  Instead, the adapter
passes SI-valued native records into these functions and copies validated
results back to center/vertex storage.  This keeps physical decisions testable
without MPI while leaving allocation and halo transport under AMPS ownership.

## Physical spherical domain and refinement

`SphericalDomain` distinguishes the solar interior, the physical heliospheric
shell, and escape beyond the exact outer sphere.  The enclosing Cartesian box
must contain the complete outer sphere; its corners are classified as escaped,
not physical space.  The surface `r=R_sun` belongs to the initialized guard
shell so particles can be event-located onto the unique absorbing surface.

Solar and tube refinement are independent requests.  For target size `h_0`,
global size `h_g`, distance `d`, and decay length `L`, each request is

`h(d)=h_0+(h_g-h_0)[1-exp(-d/L)]`.

The selected size is the minimum of the solar and authoritative traced-line
requests.  No second approximate Parker spiral is constructed.  Buffered
polyline/block intersection uses the complete segment against an expanded
block, not its center; a breadth-first face-connectivity audit then proves that
all required Sun/source/observer anchors belong to one active component.

## One-sided fields and interface events

Continuous interpolation is legal only after every stencil member has the
same region, magnetic sector, mapping validity, and interface identity.
`ValidateOneSidedInterpolationStencil` rejects a stencil crossing `R_i`, a
separatrix, or the HCS before any floating-point primitive is blended.

An ideal-HCS event changes only the categorical sector coordinate.  Position,
physical velocity, momentum, and physical outward/inward wave labels are
unchanged.  Selecting drift at an ideal zero-thickness HCS returns typed
`NotImplemented`; it requires a separately qualified finite-thickness model.

## Time steps, species weights, and population control

The local step is the minimum of the configured upper bound,
`epsilon_g/Omega`, `epsilon_x*h/|v|`, and `epsilon_s*tau_sc`.  All quantities
must be positive and finite.  Tightening any accuracy factor is therefore
monotone.  `ValidateAllSpeciesWeights` requires one positive base weight for
every compiled slot and never assumes a proton in slot zero.

Split and merge operations retain the grouping key `(species, background
generation, source label)`.  A split divides represented weight while copying
phase space.  The baseline exact merge accepts identical phase-space states;
it then closes represented weight, vector momentum, and relativistic kinetic
energy exactly.  More general stochastic AMPS merge proposals must provide
the same before/after ledger before an adapter commits them.  Neither operation
can modify the physical source ledger.

## Initialization and finite output

`InitializationLedger` records ten mandatory conditions: mesh, boundary,
background, turbulence, shock, halo exchange, species weights, time steps,
observers, and output dictionary.  Output is rejected unless all ten are set
and the shock references the current background generation.  Consequently an
initialization-only run stops after physically initialized output, not after
mesh allocation.

`MakeFiniteOutput` implements one finite sentinel accompanied by typed
validity and `particlesPresent`.  A NaN or infinity is never emitted for an
empty or inapplicable quantity.

## Observer physics

Energy bins use exact linear, logarithmic, or explicit edges.  Kinetic energy
is evaluated from velocity in the detector rest frame.  Per-nucleon output
requires a validated positive integer nucleon count; electrons or ambiguous
composition cannot be silently placed in those bins.

The three geometries are discriminated:

- a volume sphere uses represented residence time divided by volume;
- a surface disk uses accepted outward crossings divided by disk area;
- a surface sphere uses accepted outward crossings divided by `4*pi*r^2`.

Empty bins retain finite zero count/intensity and an explicit false presence
flag.  Stable observer IDs are independent of vector order.

## Verification

`make test-stage8` runs 153 cumulative tests.  The 25 Stage-8 records are
`BND3D01--02`, `MESH3D01`, `COR3D01--03`, `TIM3D01--02`, `POP3D01--03`,
`MPI3D01`, `INIT3D01--04`, `NAT3D13`, `RUN3D02`, and `OBS3D01--07`.
They cover spherical clipping, independent refinement, corridor closure,
monotone time steps, all-species weights, conservative population changes,
stable reductions, initialization ordering, finite output, detector-frame
binning, per-nucleon guards, and all observer geometries.
