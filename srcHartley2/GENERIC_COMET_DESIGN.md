# Generic comet application design

## Scope

`srcHartley2` remains the application directory for build compatibility, but
new code in it is comet-generic. A comet case should be defined by input data:

1. compiled AMPS species in `SpeciesList`;
2. computational dimension and surface geometry;
3. heliocentric conditions;
4. surface and extended source descriptions for each species; and
5. an optional chemical reaction network; and
6. a named or user-supplied neutral collision database.

The AMPS core continues to own the mesh, particle mover, collisions, sampling,
and parallel distribution. Comet-specific code supplies geometry and physical
models through small interfaces. Hartley 2 values belong in input files rather
than in those interfaces.

## Implemented layers

### Geometry and physical configuration

`CometModel.h/.cpp` defines geometry and source descriptors independently of
the numerical dimension:

- spherical geometry for DIM1;
- axisymmetric geometry reserved for DIM2; and
- triangulated geometry for DIM3.

The model stores the heliocentric distance and, for each AMPS species, the
surface production rate, surface distribution, surface temperature, bulk
speed, extended production rate, and extended distribution. Input species
names are resolved to the AMPS compile-time species indices by `hartley2.pl`.

Supported surface distribution identifiers are `DISABLED`, `UNIFORM`,
`COSINE`, `TABULATED`, `LEGACYBJORN`, and `LEGACYJET`. Supported extended
source identifiers are `DISABLED`, `RADIALPROFILE`, `GRAINSUBLIMATION`, and
`USERDEFINED`. The DIM1 backend executes uniform direct surface sources. Its
extended-source executor is not implemented. The legacy DIM3 activity paths
remain available while their source maps are migrated to these descriptors.

The relevant input records are:

```text
Dimension=1
Geometry=Spherical
HeliocentricDistanceAU=1.06
SurfaceSource=H2O | UNIFORM | 2.0e27 | 180.0 | 0.0
ExtendedSource=H2O | RADIALPROFILE | 1.1e28
ChemistryMode=off
```

The `SurfaceSource` columns are species, distribution, total rate in s^-1,
temperature in K, and bulk speed in m/s. A physical quantity unavailable from
the source material can be written as `UNSPECIFIED`; the model retains `-1`
and can reject it when a runnable source backend is requested. The
`ExtendedSource` columns are species, distribution, and integrated rate in
s^-1.

For a DIM3 triangulation, input also controls the mesh instead of relying on a
Hartley 2 literal in `main.cpp`:

```text
SurfaceMeshFile=hartley2.shape.nas, KM
```

### Dimensions

DIM1 is a compiled and executed spherical radial configuration. Its domain
begins at the explicit spherical nucleus radius. The supplied surface
resolution, ideal number of radial cells, and domain extent determine a
geometric expansion profile. AMPS approximates that ideal profile with its
dyadic block refinement. `RadialSimulationMode=OFF` builds and verifies only
the mesh. `ON` additionally installs direct-source injection, species-local
cell-crossing timesteps, a thermalizing inner boundary, an absorbing outer
boundary, and the AMPS time-step loop.

DIM3 is the existing triangulated Hartley 2 DSMC application. It compiles with
the generic model and chemistry layers and reads the input-selected NASTRAN
mesh. Its legacy activity routines have not yet been replaced by the generic
per-species source executor.

DIM2 has a configuration type but no executable backend. The current core can
construct a two-dimensional planar or axial AMR mesh, but the generic
injection code contains unverified DIM2 paths. `hartley2.pl` therefore rejects
DIM2 instead of presenting it as supported. A DIM2 implementation needs an
axisymmetric surface representation, the cylindrical cell/face measures to be
selected consistently, and source sampling weighted by the ring area.

### Chemistry database

`src/models/chemical_solver/ReactionNetwork.h` is a geometry-independent,
header-only reaction database. It stores arbitrary reactant/product
stoichiometry using AMPS species indices and supports:

- constant rate coefficients;
- inverse-square heliocentric rates referenced to a supplied AU distance; and
- modified Arrhenius rates,
  `k(T)=A*T^beta*exp(-Ta/T)`.

It validates species indices, coefficients, table ranges, rate parameters,
and excess energy. It can evaluate mass-action reaction rates and separate
species source/loss terms.

`input/CometReactionNetwork.pm` parses input records and generates
`CometChemistry.generated.h` in the build tree. The syntax is:

```text
ChemistryReaction = ID | REACTANTS | PRODUCTS | LAW | A | BETA | TA_K | RREF_AU | EXCESS_EV
```

Reaction sides use `+` and optional coefficients, for example:

```text
ChemistryReaction = H2O_PHOTO | H2O | OH + H | HELIOCENTRIC | 1.03e-5 | 0 | 0 | 1 | 3.41
ChemistryReaction = OH_LOSS | 2*OH | H2O + O | ARRHENIUS | 2.0e-12 | 0.5 | 100 | 1 | 0
ChemistryReaction = H2O_ION | H2O | NONE | HELIOCENTRIC | 4.04e-7 | 0 | 0 | 1 | 0
```

The topology is arbitrary, but execution is constrained by the selected
solver. `NONE` is an explicit loss from the tracked neutral system. The comet
PIC executor accepts unit unimolecular reactions, samples exactly one channel
per event, creates all daughters together, and conserves event momentum and
the tabulated excess kinetic energy. Integer product multiplicities, neutral
mass balance, and individual particle-weight correction are validated before
the run. Bimolecular and higher-order reactions must be applied to cell
densities with a cell chemistry solver. The old
`ChemicalSolver.h`/`chemeq2.cpp` files in the core are incomplete and are not
silently treated as an operational solver.

### Neutral collision database

`src/models/cometary_physics/CollisionDatabase.h` stores unordered species
pairs, source provenance, and either a constant cross section or a bounded
temperature power law. A missing pair is observable and the comet adapter
stops instead of interpreting missing data as a collisionless pair.

`input/CometPhysicsDatabase.pm` provides named paper databases. The first is
`RUBIN2011`, containing all 66 parent-parent pairs in Rubin et al. Table 3 and
the Table 4 photoprocesses. The H2O-H2O law is evaluated only in its published
50--300 K range. The paper does not publish a complete daughter collision
matrix, so those entries remain deliberately absent.

The Rubin database also requires explicit values for its two documented
sensitivity choices:

```text
PhysicsDatabase=Rubin2011
Rubin2011OHGroundRate=<3e-6 through 6e-6 s^-1>
Rubin2011MethanolH2COBranch=<ON or OFF>
```

## Hartley 2 case from Shou et al. (2025)

The source is Shou et al., *Surface activity of H2O and CO2 on comet
103P/Hartley2 derived from EPOXI/HRI images*, Icarus 435 (2025) 116557,
doi:10.1016/j.icarus.2025.116557.

The paper supports the following entries in `hartley2.radial.input`:

| Quantity | H2O | CO2 |
| --- | ---: | ---: |
| Direct nucleus production | `2.0e27 s^-1` | `3.5e26 s^-1` |
| Night-side activity floor | `0.1` | `0.3` |
| Spherical effective temperature | `188.350171377 K` | `182.852991011 K` |
| Surface bulk drift | `0 m/s` | `0 m/s` |
| Species represented | yes | yes |

Section 4.3 reports CO2 direct production as approximately
`3-4e26 s^-1`; Fig. 5 uses `3.5e26 s^-1`, which is the value recorded in
the input. Dust and the paper's extended gas components are omitted from the
current case as requested.

The DSMC boundary uses an upward half-Maxwellian at a facet temperature from
the Tenishev et al. (2008) solar-angle profile at 1.29 AU. That profile is the
first temperature column of the existing `CG-BC.h` table; the legacy Hartley2
injector uses it with zero bulk velocity. A spherical source cannot retain
solar-angle structure, so `shou2025_spherical_parameters.py` calculates the
source-number-weighted temperature using the paper's
`G(theta)=max(cos(theta),a_night)` with `a_night=0.1` for H2O and `0.3` for
CO2. The resulting scalar preserves the half-Maxwellian thermal-energy flux
under the spherical reduction. The paper's 600 m/s is used in the analytic
inversion relation and Haser profiles, not as DSMC injection drift.

The tracked `input/hartley2.shape.nas.gz` contains 16,022 vertices and 32,040
triangles in km. `comet_shape_properties.py` finds:

```text
surface area                 5.242876570262 km^2
enclosed volume              0.809653618778 km^3
area-equivalent radius       0.645921714366 km
volume-equivalent radius     0.578189522337 km
nonmanifold edges            0
```

The requested spherical nucleus therefore uses radius `645.921714366 m`.
The paper discusses a 500-facet model, whereas the tracked mesh has 32,040
triangles. The calculation is explicitly tied to the tracked mesh and does
not claim the two triangulations have identical areas.

The supplied DIM1 science domain extends 10,000 km above that surface and uses
1,000 cells in the ideal logarithmic profile. Shou et al. do not publish a
surface mean free path or collision cross section. Using the constant
`3.4e-19 m^2` H2O/CO2 cross section already used by
`input/comet_sphere.input`, the two half-Maxwellian source fluxes imply an
outgoing surface density of `2.071008675e18 m^-3` and a conventional
hard-sphere mean free path of `1.004209133 m`. The input therefore requests a
1 m nucleus-surface cell width. This is a core-derived resolution estimate,
not a value reported by the paper.

The AMPS dyadic mesh requires maximum refinement level 21 for this combination:
with eight cells in the root block, level 20 bottoms out at `1.192 m`, while
level 21 permits `0.596 m` cells and can enforce the 1 m target.
The associated core range check now compares the refinement level with the
usable bits in `int`; its former `1 << sizeof(int)` expression incorrectly
rejected safe levels 16 through 30 on systems with four-byte integers.
The validated mesh has logarithmic ratio `1.011741915188`, ideal outer-cell
width `116057.4 m`, 189 AMR blocks, and 1,512 allocated cells.

## Halley case from Rubin et al. (2011)

The gas-only Halley configuration is kept separate from both Hartley 2 inputs:

```text
input/halley1d.input
srcHartley2/halley2011/PAPER_PARAMETERS.md
```

The input records the 6 km nucleus, 0.90 AU heliocentric distance,
`10^6 km` cometocentric outer radius, 180 K zero-drift uniform boundary, all
27 neutral species, the 11 parent rates summing to `8.54e29 s^-1`, and the
named Rubin collision/chemistry database. The supplied case uses 2,000 ideal
radial cells, a 0.5 Courant number, disabled methanol branch, a surface width
derived from the `0.3293692849882 m` H2O mean-free-path estimate, and
characteristic speeds derived from a 400 m/s parent baseline plus Table 4
excess energies. It uses species-local particle weights proportional to local
timestep so particle occupancy remains approximately uniform across the
logarithmic mesh. The selected 0.36 eV OH branch rate is `3.0e-6 s^-1` at
1 AU. Complete 500-iteration density samples are compared with a
shell-volume-weighted relative L1 norm; all species must remain within 5% for
three consecutive comparisons after the reference window, with 10,000
iterations as the hard ceiling.

## Work needed for physical DIM1 science runs

The DIM1 executor now advances particles when all numerical controls are
explicit. A case still needs:

1. a complete collision database for every pair that can coexist;
2. explicit particle-weight, timestep, sampling, and convergence controls;
3. validated internal-energy relaxation and radiative cooling where required;
   and
4. traceable reference data for plot comparisons.

Shou et al. do not state those numerical controls. Its domain extent and cell
count are supplied run settings; the surface width is derived from the
existing AMPS comet collision prescription rather than attributed to that
paper. Rubin et al. likewise omit their exact radial discretization and
numerical controls; the Halley input labels the supplied and derived values
and keeps the remaining choices explicit.
