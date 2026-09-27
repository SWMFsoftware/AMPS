# OV3D01: 2013 April 11 near-Earth SEP validation

This directory is a reviewed *configuration blueprint*, not a ready-to-run
input deck.  `parameters.json` deliberately has
`execution_status=blueprint-needs-event-fit`: the cited publications constrain
many physical quantities, but they do not determine every parameter of the
spherical SWCME surrogate or srcSEP3D's numerical macroparticle population.
Turning missing values into convenient defaults would make the comparison
irreproducible and physically ambiguous.

## Preparation workflow

`event_fit.template.json` enumerates every event-fit, numerical-convergence,
coordinate-transform, and provenance field that is not fixed by the
publication manifest. All distributed values are `null` intentionally. Check a
working copy with:

```bash
python3 validation/cases/2013-04-11/prepare_case.py status \
  --fit /path/to/reviewed-event-fit.json
```

The command returns nonzero and lists every incomplete field. Once the fit is
complete, render the input and validate it with the actual linked executable:

```bash
python3 validation/cases/2013-04-11/prepare_case.py render-input \
  --fit /path/to/reviewed-event-fit.json \
  --output /path/to/sep3d-2013-04-11.in \
  --verify-executable /path/to/amps
```

The renderer does not copy convenient event values into arbitrary fields. It
derives the CME, Earth, and 20-solar-radius Parker vectors through the supplied
proper rotation, requires that the current +Z Parker implementation preserve
solar north, converts the total one-AU field to the consistent radial Parker
component, installs the published radial-rigidity MFP and fixed q=5 source
shape, and calculates macroparticle weight from the declared rate contract.
It removes the example-only inner observer. Every derived deck records hashes
of both the reviewed fit and this directory's parameter manifest.

After instrument-specific calibration, quality filtering, background
subtraction, and unit conversion have produced the exact two-column grammar
`coordinate_si,value_si`, build immutable comparison evidence with:

```bash
python3 validation/cases/2013-04-11/prepare_case.py build-evidence \
  --model-csv /path/to/model.csv \
  --reference-csv /path/to/calibrated-observation.csv \
  --output-dir /path/to/evidence/OV3D01 \
  --coordinate-name time --coordinate-units s \
  --value-name differential_intensity \
  --value-units 'm^-2 s^-1 sr^-1 J^-1' \
  --model-provenance 'executable/configuration/seed/rank description' \
  --reference-provenance 'dataset/version/instrument/processing description'
```

This copies the reviewed bytes and writes their SHA-256 values into the Phase-V
manifest. It does not download, calibrate, resample, or silently repair data.

The case is nevertheless a useful first observational validation for the new
narrow Parker corridor.  The event was a moderately fast halo CME with dense
near-Earth proton coverage, and two independent event studies identify the
high-altitude western shock flank as the likely source of the particles seen
near Earth.  A one-corridor srcSEP3D calculation can therefore test transport
from that local shock/field-line intersection to Earth without pretending to
reproduce the complete three-dimensional AWSoM-R corona.

## Evidence and scope

The authoritative event study is Liu et al., *Physics-Based Simulation of the
2013 April 11 Solar Energetic Particle Event*:

- <https://arxiv.org/html/2412.07581v2>

The independent connectivity interpretation is Lario et al., *The Solar
Energetic Particle Event on 2013 April 11: An Investigation of its Solar
Origin and Longitudinal Spread*:

- <https://doi.org/10.1088/0004-637X/797/1/8>

Lario et al. report that the on-disk EUV wave did not reach the nominal
near-Earth footpoint and infer injection from the high-altitude west flank of
the CME shock.  That result is essential to this setup: the resolved Parker
line must be Earth-connected and must **not** be forced to begin at the CME
apex.  The CME direction and Parker-line start are separate physical inputs.

This case can validate:

1. local shock-source activation on the Earth-connected flank;
2. parallel Sun-to-Earth transport, onset, peak, decay, and fluence spectrum;
3. the radial-rigidity mean-free-path implementation; and
4. numerical convergence with corridor width, mesh, timestep, and particle
   population.

It cannot validate the global longitudinal distribution, the detailed
AWSoM-R shock morphology, cross-field transport, or the full published seed
boundary.  STEREO-A/B products are therefore diagnostics for a spherical
SWCME run, not release gates.  A future imported time-dependent heliospheric
background can promote those comparisons after it reproduces the published
multi-spacecraft connectivity.

## Published event constraints

Liu et al. give an M6.5 flare from AR 11719 (N09E12), beginning at 06:55 UT
and peaking at 07:16 UT.  LASCO/C2 first observed the CME at 07:24 UT.  The
study uses the DONKI-re-evaluated CME speed of 675 km/s and places the flux-rope
centre at Carrington longitude 69.5 degrees and latitude 14.5 degrees.  These
map to:

```ini
cme.launch_speed = 675 km/s
# After applying the reviewed Carrington-to-model-frame transform:
geometry.cme_direction_x = 0.33905244980934973
geometry.cme_direction_y = 0.90683696982863249
geometry.cme_direction_z = 0.25038000405444144
```

Those direction components are valid only in an event-frozen Cartesian frame
whose +X axis is Carrington longitude zero, +Y is longitude +90 degrees, and
+Z is solar north.  If `domain.coordinate_frame` uses another convention,
transform the vector and record the transform in the evidence manifest.  Do
not paste the numbers into an HCI-like frame merely because both are
Sun-centred.

For Earth, Table 2 of Liu et al. gives Carrington longitude 85.3 degrees,
latitude -5.9 degrees, radius 1 AU, and a 12-hour pre-eruption mean wind speed
of 363 km/s.  In the event-frozen convention above, the unit observer vector
is `(0.08150446536468188, 0.99135801631219989,
-0.10279253678724681)`.  Multiply it by exactly 1 AU before assigning
`observer.earth.position_[xyz]_m`.

### Earth-connected Parker start at 20 solar radii

The paper's nominal footpoint longitude, 154.0 degrees, is a near-Sun
connectivity product. srcSEP3D resolves the analytic line only from the SWCME
source radius, normally 20 solar radii. Because the initialized SWCME field
uses the source-surface correction
\(B_\phi/B_r=-\Omega(r-r_s)\sin\theta/V_{sw}\), its exact field line is

\[
  \phi(r)=\phi_s-\frac{\Omega_\odot}{V_{sw}}
  \left[(r-r_s)-r_s\ln\left(\frac r{r_s}\right)\right].
\]

For the nominal \(\Omega_\odot=2.865\times10^{-6}\ \mathrm{rad\,s^{-1}}\),
the source longitude that makes the implemented 20-solar-radius line pass
through the published Earth location is

\[
 \phi_s=85.3^\circ +
 \frac{\Omega_\odot}{363\,\mathrm{km\,s^{-1}}}
 \left[(1\,\mathrm{AU}-20R_\odot)-20R_\odot
 \ln\left(\frac{1\,\mathrm{AU}}{20R_\odot}\right)\right]
 =131.7136869974263^\circ.
\]

Use `parker_spiral.start_mode=explicit`, longitude
`2.298837508046333 rad`, and colatitude `95.9 degrees` for those nominal
values. `prepare_case.py` recomputes this longitude from the reviewed Earth
location, event-fit rotation rate, and wind speed; it does not copy the
reference number when the fit changes. Do not use
`start_mode=cme-launch-point`: the CME centre at 69.5 degrees and the
Earth-connected line at about 131.71 degrees describe the apex and west flank,
respectively.  The runtime observer/corridor connectivity check then provides
an independent guard against a sign or frame error.

## Transport setup

For the first comparison select the Parker mover and parallel transport only:

```ini
[run]
transport = parker

[transport]
spatial_diffusion_model = mean-free-path
mean_free_path_model = radial-rigidity-power-law
mean_free_path_reference_m = 4.487936121e10
mean_free_path_reference_radius_m = 1.495978707e11
mean_free_path_reference_rigidity_v = 1e9
mean_free_path_radial_exponent = 1
mean_free_path_rigidity_exponent = 0.3333333333333333
perpendicular_diffusion = none
drifts = none
```

This is the species-general representation of Liu et al. Equation 15,

\[
 \lambda_\parallel=0.3\,\mathrm{AU}
 (r/1\,\mathrm{AU})(\mathcal R/1\,\mathrm{GV})^{1/3},
 \qquad \mathcal R=pc/|q|.
\]

The implementation evaluates rigidity from the actual compiled AMPS species,
not from a proton-only kinetic-energy approximation.  Liu et al. found
0.3 AU to agree better with the event profiles than their 0.1- and 1-AU
experiments, but explicitly describe that value as tuned.  It is therefore a
registered event parameter, not a universal turbulence law.

After the Parker baseline passes, run `focused-diffusion` and
`focused-scattering` as model-form diagnostics.  They must use the same
background, source, random campaign, corridor, and observer definitions; only
the mover/coefficient contract may change.

## Shock source shape and normalization boundary

The publication uses a suprathermal boundary `f(p) proportional to p^-5`, a
10-keV injection energy, and injection coefficient `c_i=1`. It later applies
a factor 1.2 to the reported flux. The source *shape* and injection energy now
map explicitly to schema 4:

```ini
[source]
minimum_energy_j = 1.602176634e-15
spectrum_model = fixed-phase-space-power-law
phase_space_power_index = 5
```

The positive input value is q in `f(p) proportional to p^(-q)`. The sampler
therefore uses `dN/dp proportional to p^-3` after applying the isotropic
momentum-space Jacobian. `VFY3D06` independently checks that transformation and
its sampled CDF.

The normalization values must not be copied mechanically:

- `[source].injection_efficiency` is a probability-like fraction constrained
  to `[0,1]`, while the published 1.2 is an output flux scale and is never a
  legal injection efficiency; and
- `physical_particle_rate_per_s`, `samples_per_step`, and
  `macroparticle_weight` determine representation and normalization in the
  Monte Carlo implementation. The cited M-FLAMPA calculation does not
  determine those three values.

Equation 21 specifies a boundary phase-space density from local proton density
and temperature. srcSEP3D presently injects a declared physical number per
unit time. A conversion between those quantities additionally requires a
reviewed boundary-flux or source-volume operator; density and temperature alone
do not determine a birth rate. The implementation therefore closes the
published spectral-shape gap while retaining this normalization mapping as an
explicit campaign prerequisite rather than inventing a dimensional factor.

Record the rate-normalization choice in the evidence provenance. Select the three
Monte Carlo controls using a particle-count convergence study, and retain the
case registry's single global amplitude normalization.  Re-fitting a separate
amplitude for each energy channel or time interval is not allowed.

At 24 minutes after eruption the published Earth cobpoint is a weak,
quasi-parallel part of the shock: compression ratio about 1.5,
`theta_Bn` about 30 degrees, and fast-mode Mach number about 1.  Use those as
local diagnostic targets.  They are not direct inputs to the current SWCME
interface and must not be manufactured by overwriting its solved shock state.

## Why some SWCME fields remain unresolved

The Liu et al. calculation uses a Gibson-Low CME in a time-dependent AWSoM-R
background.  srcSEP3D's current released standalone path uses a spherical
SWCME/DBM surrogate.  CME speed and direction transfer; the following do not:

- DBM drag coefficient;
- sheath/ejecta thickness and smoothing widths;
- analytic ambient density, temperature, magnetic magnitude, and polarity;
- the spherical surrogate's shock timing at the Earth-connected flank; and
- numerical corridor radius and block halo.

Derive the ambient state from quality-screened OMNI plasma/MAG observations
over a declared pre-event window.  Fit drag and widths to an independently
chosen in-situ shock/ICME marker.  That background fit must be frozen before
examining SEP residuals.  Determine corridor and numerical controls through
XM3D05-style convergence.  `parameters.json` lists every unresolved item so a
campaign-preparation tool can refuse an incomplete deck.

## Active corridor and population control

The first production mesh should use `mesh.active_region.mode=parker-tube`
around the Earth-connected line.  Its radius must include the observer
collection sphere and the complete shock-flank/source intersection throughout
the injection interval. `buffer_blocks` is an integer number of AMR
touching-neighbour layers, independent of local leaf size. Start
conservatively, keep at least one such halo layer, and repeat with a wider
physical radius until onset, peak, and fluence change by less than the
registered convergence tolerance. A visually plausible narrow tube is not
convergence evidence.

Enable `population_control.mode=split-merge` only after an uncontrolled pilot
establishes count distributions.  Choose hysteretic minimum/target/maximum
counts per cell and compiled species.  Demonstrate invariance with a second,
larger target and preserve the same physical weights and source rate.  The
controller is a variance/cost device, not a source-normalization parameter.

## Observation products

The primary release comparison is near Earth:

| Product | Instrument and channels | Treatment |
|---|---|---|
| Time-intensity | SOHO/ERNE 2.0--2.5 MeV and 20.0--25.0 MeV | Background-subtracted, common cadence, compare onset/peak/decay |
| Fluence spectrum | ACE/EPAM LEMS120, SOHO/ERNE, GOES-13/EPEAD | Integrate first 259200 s; fit spectral slope only over 1--50 MeV |
| Intermediate checks | SOHO/ERNE 4--5, 8--10, 13--16, and 32--40 MeV | Diagnostic residual structure; no channel-specific amplitude fit |

Use native quality/status flags and preserve instrument energy-bin bounds.
The paper warns about SOHO/ERNE saturation at high flux and uses calibrated
GOES-13/EPEAD effective energies.  It also notes low-energy instrumental/event
structure.  Exclude or down-weight points only by a rule written before the
model comparison, and retain rejected rows plus reason codes in the evidence
bundle.

Authoritative access points are:

- SOHO/ERNE: <https://cdaweb.gsfc.nasa.gov/misc/NotesS.html#SOHO_ERNE-HED_L2-1MIN>
- ACE/EPAM: <https://cdaweb.gsfc.nasa.gov/misc/NotesA.html#AC_H2_EPM>
- GOES 1--15 instruments/data: <https://www.ncei.noaa.gov/products/goes-1-15/space-weather-instruments>
- pre-event OMNI plasma/MAG: <https://omniweb.gsfc.nasa.gov/>

The observation preparer owns downloads, calibration, background subtraction,
uncertainties, and checksums.  `validation/run_validation.py` consumes only the
prepared immutable `OV3D01/manifest.json`, `model.csv`, and `reference.csv`.
It never downloads or silently updates data during a qualification run.

## Required campaign order

1. Freeze the coordinate transform, source epoch, observation interval, and
   all data-quality rules.
2. Fit only the solar-wind/SWCME surrogate to plasma, field, and arrival data;
   do not inspect SEP scores during this step.
3. Verify the 20-solar-radius Parker curve reaches the Earth observer and the
   active corridor contains every accepted source patch.
4. Pass mesh, timestep, corridor-width, and particle-count convergence.
5. Produce time-intensity and fluence products with a single shared physical
   run configuration.
6. Build a checksum-owned OV3D01 evidence bundle and run
   `validation/run_validation.py --case OV3D01`.
7. Interpret only the pre-registered metrics.  A missing external bundle is a
   `SKIP`, malformed provenance is an `ERROR`, and an out-of-tolerance result
   is a `FAIL`.
