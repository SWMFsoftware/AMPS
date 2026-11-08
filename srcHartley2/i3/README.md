# 3I/ATLAS six-case radial campaign

This directory owns the observation-derived data and deterministic input
generation for the six MAVEN/IUVS 3I/ATLAS cases supplied in `Archive.zip`.
It does not put a 3I constant into the generic comet C++ application.

Run:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 srcHartley2/i3/prepare_i3_cases.py
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest \
  srcHartley2/i3/tests/test_prepare_i3_cases.py
```

The first command reproduces `input/i3.case1.input` through
`input/i3.case6.input` from `input/i3.input` and the reviewed CSV files under
`reference/`. `input/i3.input` is a generator template, not a directly
loadable AMPS input. Each concrete case is loadable, but AMPS configuration is
compile-time: preprocess and compile one case at a time. Remove exactly the
repository `build/` directory before each recompile; never reuse a build made
for another case. The six outputs likewise require distinct run directories
or archived `PT` trees.

## Archived observation table

The supplied archive has SHA-256
`344d33e69c913a5e76d57ddb7c92c1a8d1e52a92a03c06f2f5849f70e45597be`.
Its embedded `Model_Table_6cases.pdf` has SHA-256
`a88eadce0b9dbdd3f4fda0974ff00b6a791d516ff6be497e732b57bdb7f34e0c`.
`reference/i3_cases.csv` transcribes its exposure times, heliocentric distance,
radial velocity, spatial resolution, and table production rates. Rates are
interpreted as molecules per second. H2O values inferred from the raster red
curves are retained in the separate `plot_fit_q_h2o_s-1` column for future
comparison and are not used as source rates.

The observations are in 2025. This is independently consistent with the NASA
mission report; the archive file creation date is not an observation date.

## Solar-activity correction and reaction rates

The rate calculation implements the user-selected Moore and Mapaye (2025)
per-channel law

```text
j_1AU = A0 * F10.7P^alpha,
F10.7P = (daily adjusted F10.7 + centered 81-day mean)/2.
```

For every exposure endpoint, Sun-centered J2000 ecliptic vectors for Earth
and 3I were obtained from the NASA/JPL Horizons API (3I packed designation
`K25N010`). The same solar face was found with Carrington sidereal rotation
`360/25.38 deg/day`; the resulting Earth-equivalent interval was averaged
using official NRC adjusted 20-UT daily flux. `solar_activity_audit.csv`
retains endpoint shifts, unshifted controls, final shifted values, and the
method. The source snapshots used during derivation had these hashes:

- NRC daily flux table: `5962ab4b3dda01beab622fcfd75634908aca00b716bfdb6c2f22d5d10119330a`;
- Horizons 3I vectors: `82c17649f9a06df8c6542cc7f741ead0673531d088b113a58cce18a660710aab`;
- Horizons Earth vectors: `d809d1ab766343839afbde35d2ff7a5fe8a1f589a5dd9ea742e15191e3febd93`.

The authoritative services are
`https://www.spaceweather.gc.ca/solar_flux_data/daily_flux_values/fluxtable.txt`
and `https://ssd.jpl.nasa.gov/api/horizons.api`. No runtime network access is
needed. The generated chemistry records contain the calculated 1-AU value;
the existing `HELIOCENTRIC` law then applies the case's inverse-square factor
once.

`photochemistry_channels.csv` contains the full neutral H2O descendant
network selected for this study: three H2O dissociation channels, three
state-labeled OH channels, two H2 dissociation channels, and neutral-system
photoionization/dissociative-ionization sinks for H2O, OH, H2, H, and O. O
electronic products are collapsed into the compiled O species, but reaction
IDs, rates, and excess energies remain distinct. H2O/H2 energies are from
Huebner and Mukherjee (2015) as tabulated by Carberry Mogan et al. (2022).
OH excess energies are state-labeled published values: 1.27 eV for O(3P)
from Bryans and Pesnell (2012), and 7.73/10.0 eV for O(1D)/O(1S) from
Bisikalo et al. (2015), Table 1. Bryans and Pesnell list the close 7.90/9.80
eV alternatives; the generator deliberately uses one cited table rather than
silently averaging the sources.

CO and CO2 are direct collision partners only. No CO/CO2 photochemistry is
invented.

## Source distribution, mesh, and collisions

All parents use a uniform 175-K, zero-drift surface source. This does not
sample an ordinary volume Maxwellian. The AMPS boundary injector samples

```text
p(v_vector) proportional to v_normal f_M(v_vector), v_normal > 0.
```

Consequently the speed marginal is `v` times the ordinary Maxwell speed
distribution and the angular marginal is the cosine law `p(mu)=2*mu`.
Residence weighting recovers the outgoing half-volume-Maxwellian used below.

For parent `j`, the surface number density is

```text
n_j = Q_j / [4*pi*R_n^2 * sqrt(2*k*T/(pi*m_j))].
```

For each direct parent `i`, the generator evaluates

```text
lambda_i = 1 / [sqrt(2) * sum_j(n_j*sigma_ij)]
```

with the exact collision records selected by the input at 175 K. The requested
first-cell width is the minimum of the H2O, CO, and CO2 mixture paths. Across
cases 1--6 it is respectively 1.75261509304, 1.29785287642, 0.531312387076,
0.131442336866, 0.122734766819, and 0.131442336866 m; CO2 is limiting in each
case. The domain's cometocentric outer radius is exactly 1e9 m (10^6 km), with
2,000 ideal radial cells and a 2.5-km spherical nucleus.

The collision hierarchy is explicit. Three CO2-parent pairs retain Rubin et
al. (2011) Table 3. The selected completion adds Combi (1996) H2O/OH/H2/H/O/CO
values, followed only where needed by CHEMKIN-diameter/Lorentz hard-sphere
estimates. The unsupported OH--CO2, O--CO2, H2--CO2, and H--CO2 pairs remain
labeled estimates in the generated C++ table.

The direct-source injection reference is derived for at least 800 particles
per first-cell parent population and is 114 events per shortened parent
timestep. Regular block-local adaptation keeps the validated 800--3,200
range for directly represented populations; chemically produced species use
the separate 50--200 range, and direct-source weights remain locked outside
emergency common-factor coarsening. These are statistical controls and do not
alter the physical production rates above.

## Brightness scope and limitations

Each case requests 200 logarithmic impact parameters from the nucleus
surface to the finite 10^6-km domain edge. The input selects the reusable core
models `OH_AX00_SA1988_LOWELL_UNQUENCHED` and `H_LYA_IKL2020`, and supplies
the observation's heliocentric radial velocity, explicit Sun/comet/observer
geometry, and 1-AU LASP V4 Ly-alpha photon flux. The Lowell source quantity
is an integrated energy rate; each case explicitly selects 308.5 nm to
reproduce the archived photon conversion.

Brightness is sampled from individual model particles during every AMPS
sampling iteration. At each line-of-sight quadrature point, the DIM1 local
three-component velocity is rotated into the physical radial direction and
combined with the comet velocity before evaluating the heliocentric radial
Doppler speed. Thus H Ly-alpha uses the simulated velocity distribution; the
archived table's scalar comet radial velocity only supplies the bulk vector
`V_comet=Vr*e_r`. The velocity-dependent core model is evaluated per particle
without clamping or extrapolation. A conservative line-center cross section
`9.5e-18 m2` with threshold `tau=0.1` flags where optically thin output needs
full frequency-dependent transfer. It does not attenuate the result.

The geometry vectors in `reference/i3_cases.csv` were derived at each
exposure midpoint from geometric J2000-ecliptic JPL Horizons states for the
comet and Mars center: `e_r` points from the Sun to the comet, `LOS` points
from the comet to Mars, and the positive impact direction is the projection
of `e_r` onto the observer sky plane. The vectors are explicit case data, not
compiled constants. They are unit normalized and `LOS dot IMPACT=0` is
enforced by both the generator and runtime parser. The archived scalar radial
velocities remain authoritative; current-solution Horizons velocity was not
substituted for those values.

The frozen LASP source subset used for derivation had SHA-256
`672c9577fad0019d8ee7a6a8198604afd4671a5895ef66ee3e0cc544aee9a5db`.
The authoritative dataset is the LASP/LISIRD Version-4 daily composite
(Machol et al. 2019). OH factors follow Schleicher and A'Hearn (1988), using
the Lowell Observatory service. A Sun-directed H radiation-pressure force is
not applied in exact spherical symmetry; replacing it with a cometocentric
radial force would be physically incorrect. Particle-resolved brightness is
written as ordinary Tecplot data using
`PT/plots/pic.<species>.brightness.s=<species-index>.out=<output-index>.dat`.
The integral stops at the finite domain and performs no unreported outer-tail
continuation.

### Shared velocity-dependent fluorescence data

The authoritative H and OH reference data and evaluators now live in the
application-independent core module
`src/models/resonance_fluorescence/`. Its manifest records source URLs,
citations, units, valid velocity intervals, retrieval snapshot, and SHA-256
hashes. Runtime values come from the immutable compiled representation: there
is no network or data-file access in a sampling loop.

The H model evaluates the irradiance-dependent IKL equations directly, so a
case-specific H velocity table is neither stored nor needed. The OH model
interpolates all 1,201 reviewed, non-extrapolated Lowell samples from -60 to
+60 km/s. Both models reject velocities outside their supported intervals.
Core tests reproduce the six archived scalar case values (within the stated
rounding tolerances), verify inverse-square scaling, and audit every source
hash. Run them with:

```sh
make -C src/models/resonance_fluorescence test
```

## Validation boundary

The focused generator test verifies all six shifted solar indices, reaction
count and reachability, minimum mixture paths, injection count, recoil-derived
characteristic speeds, brightness records, and byte-for-byte regeneration of
all six tracked inputs. Case 1 also has a clean C++17/MPI build, a 200,000-
sample-per-parent test of the production crossing-flux injector, and a
20,000-event-per-neutral-reaction compiled chemistry test. All eight neutral
product channels passed mass, momentum, excess-energy, and ensemble checks;
all nine neutral-system ionization channels were recognized as sinks. The
case-1 AMR mesh built 376 blocks/3,008 dyadic cells and resolved the requested
surface width.

The pure brightness tests verify the inverse DIM1 velocity rotation, exact
shell chords, nucleus occultation, Rayleigh units, tau flagging, strict input
policy, velocity-domain policy, geometry, and exact source generation.  The
repaired 12-rank case-1 run has now produced finite, 200-ray H and OH Tecplot
tables at its first completed 10,000-iteration sampling boundary.  The
companion `visualization/plot_i3_brightness_comparison.py` overlays those
profiles on the reviewed numerical raster digitization in
`visualization/references/i3_brightness_MAVEN/`.  The digitization preserves
source-pixel anchors, hashes and limitations and can be rebuilt with
`visualization/extract_i3_brightness_maven.py`; ordinary plotting does not
require `Archive.zip`.  None of this evidence establishes convergence,
optically thin H validity, or agreement with the observations. Those require
later sampling windows, original quantitative uncertainties beyond what can
be resolved in the supplied raster, and review of every H tau flag.
