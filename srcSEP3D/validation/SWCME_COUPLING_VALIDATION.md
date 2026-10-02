# SWCME propagation in srcSEP3D: 20 solar radii to 1 AU

## Purpose and present implementation status

This case tests the **moving SWCME shock supplied to a running srcSEP3D/AMPS
application**, with no energetic particles. It compares unshifted shock-front
kinematics and spacecraft arrival with observations. It does not test SEP
acceleration, and it does not establish a resolved magnetohydrodynamic CME.

The comparison/figure consumer is implemented in
`run_swcme_coupling_validation.py`. `CME3D01` exercises its mechanics;
`CME3D02` is an external observational case. A reviewed bundle produces a
600-dpi PNG, vector EPS, caption, metrics and checksum report. Missing evidence
is SKIP; corrupt evidence is ERROR; a valid scientific disagreement is FAIL.
Synthetic fixtures are never an observational PASS.

The native propagation control is implemented for schema 4. The parser accepts
`run.intent=shock-propagation`, requires `source.enabled=false`, population
control off and a fresh run, and rejects incomplete launch-to-budget temporal
coverage. The canonical provider evaluates its physical front without building
injection patches. The normal AMPS step still executes; rank-zero telemetry
records the installed provider after collective clock/identity/particle checks.
The existing seven SCCM3D initialization tests remain initialization tests.

After installing these sources, **rebuild AMPS** and run from the AMPS root:

```sh
make -j16 amps
python3 srcSEP3D/validation/run_swcme_coupling_validation.py --amps ./amps --ranks 10 --input srcSEP3D/examples/sep3d_swcme_20rs_1au.in --output-dir test_output/swcme-native-control
```

The runner executes allocation-free native preflight, freezes the input and
executable identity, invokes `mpiexec -n 10`, streams diagnostics and retains
raw bytes in `native.log`. It rewrites only the three initialization output
filenames to a fresh native artifact directory. For a scheduler use
`--launcher 'srun -n {ranks}'`. Run the Python driver once, outside mpiexec.
A nonempty output directory is rejected; choose a new directory for each run.

The supplied **unfitted mechanics control** has a 1000-km/s launch, 400-km/s
ambient wind and drag `1e-8 1/km`. This smaller drag retains a physical fast-mode
shock through 1.05 AU in the canonical provider. It is not a fit to the July
2012 observations. The normal native executable can also run directly:

```sh
mpiexec -n 10 ./amps --input srcSEP3D/examples/sep3d_swcme_20rs_1au.in --output-dir test_output/swcme-direct/native
```

It writes `shock-history.csv` and, after successful close, `native-runtime.json`
in that output directory. Direct execution retains the initialization filenames
from the input. There are two stop conditions: the first completed tick at or
beyond `run.stop_shock_radius_m=1.57077764235e11` (1.05 AU), or the hard budget
of 7200 steps × 60 s = five days. The last two native samples bracket the radius
stop. `stop_shock_radius_m=0` disables that stop and retains the step budget;
positive targets must lie above launch and at/below the outer physical radius.
The 1.1-AU outer domain and 1-Rsun absorbing Sun retain the active corridor,
near-Sun sphere and x-y-corner/z-centered geometry.

The runner produces `native-propagation-control.png/.eps` on the relative clock
and checks a real 1-AU crossing, zero particles/injections, current generations,
physical shock activity and rank agreement. Successful control checks report
`native_mechanics_pass=true`; the observational case remains **SKIP** until a
reviewed absolute launch epoch and parameter fit are supplied. No observed
arrival is used to shift the clock. Native physical inactivity remains recorded
and causes a campaign FAIL, rather than being relabeled as a shock.

For an arrival comparison, independently fit the shock's 20-Rsun crossing and
launch/drag/ambient parameters using remote observations. Freeze that input and
fill `templates/swcme_event_fit.template.json`, including its original input
SHA256, launch UTC, calibration cutoff, provenance and `fit_used_holdout=false`.
Then launch with `--event-fit reviewed-event-fit.json`. The runner verifies the
binding, creates a copied evidence bundle and native manifest, and evaluates
against the installed checked Wind marker. It generates the UTC comparison PNG,
EPS and caption. This arrival-only comparison does not validate a shock-radius
track: the downloaded HELCATS traces identify brightness elongations. Native
execution can be long and requires the configured MPI host; portable regression
tests do not claim to have performed that run.

`CME3D03` tests the production parser/factory and canonical source-free DBM
provider against an independent drag-law oracle. `CME3D04` tests telemetry clock,
generation, MPI-spread, source-off, identity, close and overwrite rejection.
`CME3D01` now also exercises launch/provenance/report orchestration with explicit
OS fixtures. These tests are discovered by `test/run_tests.py --all` and
`--suite phase-v`. Native rank/time convergence remains external campaign work.

For full solar-wind/IMF comparisons, a fourth extension is required: publish
SWCME sheath/ejecta density, velocity, temperature and magnetic field through a
srcSEP3D background adapter, with a defined thermodynamic closure and derivative
contract. The current standalone background authority is analytic Parker;
SWCME supplies spherical SHOCK_ONLY/SOURCE shock geometry/compression. Plotting
that Parker field against CME sheath observations would misrepresent what the
coupling implements. This consumer therefore plots kinematics only.

## Recommended observational event

Use the **12 July 2012 Earth-directed CME**, followed at Wind near Earth on
14 July. Hu et al. (2016) provide a primary multi-spacecraft study combining
STEREO imaging, MESSENGER and Wind. Its changing deceleration makes it useful
for testing the limits of a constant-drag prescription; a miss is a meaningful
scientific result, not a reason to adjust acceptance after the run.

Do not set the 20-Rsun launch epoch equal to the flare time or first LASCO
appearance. Infer the epoch and local speed at **the tracked shock's crossing
of 20 Rsun**, with uncertainties. If the available track identifies only the
bright ejecta leading edge, either extract a separate shock track or add a
reviewed shock stand-off relation first. A flux-rope GCS fit is not automatically
a measurement of the shock radius used by this coupling.

| Observation | Use | Acquisition and cautions |
|---|---|---|
| SOHO/LASCO C3 and STEREO/SECCHI COR2 | Launch direction, width, shock radius/time near 20 Rsun | Preserve calibrated image identifiers, front selections, geometry and uncertainty. Separate shock arcs from the bright ejecta. |
| STEREO-A/B HI1/HI2 time-elongation tracks | Withheld propagation constraints beyond the launch-fit interval | Convert elongation with a declared triangulation/front geometry; retain the raw track. HELCATS HIGeoCAT is useful for event identification and geometric fits, not a set of raw measured radii or independent arrival observations. |
| Wind MFI `WI_H0_MFI` and proton SWE `WI_H1_SWE` | Locate the near-Earth shock and identify upstream/downstream jumps | MFI provides composite field data including 3-s products; SWE H1 is the 92-s definitive proton product. `WI_H0_SWE` is an older electron product, not the 2012 proton dataset. Retrieve July 12–16 to include context. |
| Wind ephemeris and shock catalog (IPShocks) | Actual spacecraft radius, shock time, normal and uncertainty | Freeze the spacecraft-specific marker. Do not merge Wind, ACE and bow-shock-shifted OMNI timestamps. |
| ACE `AC_H0_MFI`, `AC_H0_SWE` | Independent identification/cross-check | ACE products are 16-s MAG and 64-s SWEPAM. Account for location/time differences. |
| Wind/WAVES type-II emission | Additional diagnostic | Inferred heights depend on the density model and emission location. Hu et al. identify streamer/flank ambiguity; do not treat this as an independent exact shock-apex track. |

Possible second cases are the 3 April 2010 event (a high-speed stream behind the
CME challenges constant ambient wind) and the 23 July 2012 STEREO-A event
(interaction/preconditioning stress test). Neither should be advertised as a
simple isolated constant-background control. Do not use DONKI/ENLIL predictions
as observed arrival times. Energetic-particle intensities are outside this test.

## Native setup and campaign controls

| Control | Value or rule |
|---|---|
| Solar absorbing body | 1 Rsun, center at the heliocentric origin |
| CME/shock start | 20 Rsun = 1.3914e10 m, `event.launch_epoch=0`; absolute UTC epoch fitted from remote imaging |
| Outer physical extent | 1.1 AU; 1 AU is an interior diagnostic, not an escape face |
| Domain | Existing active Parker corridor plus near-Sun sphere; Sun at the x-y corner and centered in z, as in the current geometry |
| Observer | Actual Wind heliocentric position/radius at crossing; record a separate exact-1-AU crossing |
| Clock | Initial dt 60 s, native shock/diagnostic cadence 60 s; then 30 s and 15 s for convergence |
| Duration and validity | 5 days initially (432000 s); extend if front has not reached 1 AU. SWCME validity starts at launch and covers the entire run. The example's 3600–86400-s interval is insufficient. |
| Kinematics | Canonical SWCME signed DBM, with launch speed, ambient wind and drag independently frozen; retain a ballistic numerical control |
| Particles | Source disabled, no initial population, no splitting/injection; keep normal species/boundary initialization |
| MPI | Same deck on 1, 4 and 8 ranks; compare geometry, clock, fingerprints and crossings |
| AMPS fields | Analytic Parker background for the initial shock-kinematic test; label this authority in all run metadata |

The active mask does not advect with the front. Confirm it covers the full
selected trajectory and the spacecraft sampling point; expand the corridor
where necessary. A Parker field line can leave the Sun–observer radial line.
Do not infer observational coverage from the radial endpoint alone. The
near-Earth observer must reside in allocated physical blocks, outside the Sun.
Retain the existing cut-cell distributed initialization output for geometry QA.

The spherical Sun-centered provider is an approximation to a real directional
CME. Run the main comparison only for a front/observer geometry whose mapping
is explicitly reviewed. Off-axis arrivals and ejecta structure need an enhanced
front/background adapter; an active computational corridor does not create a
finite-width physical CME. SWCME prescribes front evolution analytically;
AMPS currently consumes it rather than solving the CME MHD advection equations.

## Calibration, holdout and acceptance

Freeze launch parameters from shock imaging around 20 Rsun. Freeze background
wind from a reviewed upstream context and freeze drag from an independent prior
or a declared early-imaging fit. Declare `calibration_end_utc` before processing
the later imaging and Wind arrival. Every accepted track point has a role:
`calibration` at/before the cutoff, `validation` after it. Never time-shift the
model, optimize drag against Wind arrival, or compare a reconstruction fitted
to the entire track as a withheld prediction.

Initial **proposed**, event-specific acceptance values in the template are:

* absolute spacecraft shock-arrival error at most 6 hours;
* RMS of withheld radial residuals divided by supplied one-sigma uncertainty at
  most 2, with at least five withheld points and 90% temporal coverage;
* history gaps at most 300 s; actual cadence preferably 60 s;
* MPI radius spread at most 1 m and clock spread at most 1e-9 s;
* launch radius within 1 m of 20 Rsun; continuous active shock; increasing
  provider generations; zero allocated/injected particles throughout;
* actual front crossing of 1 AU, independently of the spacecraft-radius gate.

Freeze and review these before evaluating the event. The normalized residual
is a descriptive metric, not a chi-square significance: geometrically inferred
points can have correlated errors. Arrival uncertainty is recorded separately
from the 6-hour model tolerance. The runner does not claim that these provisional
thresholds are a literature-established production acceptance standard.

Further native campaign verification (the portable prerequisites do not replace these MPI runs):

1. Ballistic control: native radius `R0 + V0*t` and exact crossing
   `(Rtarget-R0)/V0`, independently computed, including a crossing between ticks.
2. Signed-DBM control: compare the canonical native provider with an independent
   numerical integration of `dR/dt=V`, `dV/dt=-gamma*(V-w)*abs(V-w)` for faster,
   slower and equal-to-wind launches. Be explicit that gamma is in inverse
   metres (divide inverse-kilometre inputs by 1000). Check both position and speed.
3. dt convergence: 60/30/15 s; proposed crossing difference 30→15 s below 30 s.
   Exact analytic evaluation can make this a sampling/interpolation convergence
   test rather than evidence of finite-volume CME advection order.
4. Rank invariance: 1/4/8-rank histories agree within frozen tolerances; exercise
   distributed telemetry, not only repeated standalone provider evaluation.
5. Restart continuation is currently rejected for propagation mode. Supporting it
   requires retained cumulative counts/generations and an explicit history
   continuation contract; the first release requires one fresh complete run.
6. Coverage failure: short validity or a missing/outside corridor observer must
   fail explicitly. A zero-step initialization must not pass propagation.

## Evidence and native history contract

Copy `templates/swcme_coupling_evidence.template.json` into a reviewed directory
as `manifest.json`; fill every null value. Retain input bytes, raw execution log
and normalized observations beside the history. Each file descriptor contains
its relative filename and SHA256. Preserve original observation products,
quality flags, retrieval/version information, coordinate transforms, image
geometry and extraction scripts in the campaign archive; summarize them in the
provenance fields. The manifest's producer claim is auditable provenance, not
cryptographic proof of native execution. Retain the executable by checksum.

Native CSV columns are SI units unless noted:

```
time_s,tick,shock_radius_m,shock_speed_m_s,shock_active,generation,particle_count,injected_particle_count,mpi_radius_spread_m,mpi_clock_spread_s
```

The exporter must run at launch tick zero and after each successful clock advance,
using the installed canonical shock state at that same epoch. Reduce rank-min/
max radius and clock to obtain the spreads; sum owned particle populations and
cumulative injections without counting ghost copies. Verify provider identity
and configuration fingerprint on every rank before writing rank-zero rows.
Flush logs and close the file before hashing it. On restart preserve cumulative
counts/generations or expose an explicit continuation epoch; the simple v1
consumer expects one launch-to-arrival history with strictly increasing generations.

Observation CSV columns:

```
time_utc,radius_m,sigma_radius_m,quality,role
```

Use explicit UTC timestamps, positive radial uncertainties, `good`/`bad` quality
and `calibration`/`validation` roles. Bad points are counted and omitted, never
converted from fill values into physical data. Good points must be chronological.
They must identify the shock front under the declared Sun-centered spherical
mapping. The `arrival` object separately supplies the spacecraft-specific shock
timestamp, uncertainty and heliocentric radius. Solar-wind bulk speed must not
be inserted as shock-front speed. Normal shock speeds need their normal,
reference frame and geometry before any separate speed comparison is added.

## Running and publication figures

### Installed observational references

The separately distributed `SWCME_20120712_reference_data_20261001.tar.gz`
contains actual July 12–17 Wind measurements, ephemeris, catalog shock arrival,
and HELCATS time-elongation profiles. Extract it at the AMPS root, alongside
the updated source overlay. It installs `validation/reference_data/CME3D02/`
inside srcSEP3D. Both evidence runners and the dedicated consumer automatically
discover that path; an explicit evidence-root/bundle option takes precedence.
The [data preparation README](cases/2012-07-12/README.md) specifies dataset
versions, UTC intervals, unit conversions, quality flags, acknowledgements and
the reproducible downloader.

```bash
python3 srcSEP3D/validation/run_swcme_coupling_validation.py --output-dir test_output/swcme-reference-check
python3 srcSEP3D/test/run_tests.py --all --output-dir srcSEP3D/test_output/all
```

The reference manifest uses `srcsep3d-swcme-reference-v1`. The consumer verifies
every declared checksum and required CSV before reporting SKIP with
`reference_ready=true` if native history is absent. Observations no longer need
to be downloaded during a test run. The bundle includes 600-dpi PNG/vector EPS
reference previews without model curves.

Its frozen Wind shock marker is 2012-07-14 17:39:09 UTC (IPShocks Zenodo v1);
the adopted +/-20-s timing uncertainty comes from CfA's independent event
00525. Wind is at 1.00773304 AU in the NASA **Predicted Orbit** product at that
time; the position/radius source is explicit and reviewable. SWE poor fits and
gaps remain rejected: its two-hour downstream summary contains no accepted
samples. Separate 3DP ion moments provide additional context without replacing
SWE measurements or establishing a plasma-coupling acceptance test.

HELCATS brightness tracks do not yet identify an independent radial shock.
Accordingly the bundle's scope is **arrival-only**. A later `native-manifest.json`
supplies native input, log, history, executable identity, launch clock and frozen
parameters; the consumer uses only the checksum-verified observed arrival.
It retains all zero-particle, clock, MPI, generation, cadence and actual 1-AU
crossing checks, but omits unavailable radial-track gates. The resulting two-panel
comparison is explicitly an arrival diagnostic. A PASS of that narrower scope
does not establish full radial evolution, CME plasma/IMF or production readiness.
The original full-track schema and three-panel comparison remain available for
an independently reviewed shock-track campaign.

### Native comparison campaign

From the AMPS root, after the native prerequisites have produced reviewed data:

```bash
python3 srcSEP3D/validation/run_swcme_coupling_validation.py --bundle /path/to/evidence/CME3D02 --output-dir test_output/swcme-20120712-run01 --require-evidence
```

Choose a fresh output directory. Results are `result.json`,
`swcme-coupling-comparison.png`, `swcme-coupling-comparison.eps`,
`figure-caption.txt` and the copied `evidence-manifest.json`. Figures are also
written for a valid scientific FAIL, so misses remain inspectable. ERROR does
not publish comparison figures. Panels show shock radius, modeled shock speed,
and withheld radial residuals versus unshifted UTC. Observed and modeled
spacecraft crossings are marked; the exact 1-AU crossing is reported in JSON.
There are no invented measured speed, density or IMF curves.

Both external runners discover `CME3D02`:

```bash
python3 srcSEP3D/validation/run_validation.py --all --evidence-root /path/to/evidence --output-dir test_output/phase-v-run01
python3 srcSEP3D/test/run_tests.py --all --validation-data /path/to/evidence --output-dir test_output/all-run01
```

The standard test runner also runs `CME3D01` mechanics tests. The earlier
`run_coupled_sep_corona.py` aggregates shared-corona and native initialization
descriptors; its default seven SCCM3D callbacks are not this observational
campaign. Use the external evidence runners above for the figure case.

An immediately runnable **synthetic software/figure preview** is:

```bash
python3 srcSEP3D/validation/run_swcme_coupling_validation.py --demo --output-dir test_output/swcme-synthetic-preview
```

Its title/caption and JSON state SYNTHETIC/SKIP. It must not appear in a
publication as model–observation validation. A real event figure and numerical
validation result require the event fit and native export. The measured
references are supplied separately; neither a native history nor a fitted
20-Rsun shock launch is fabricated in that package.

## Primary references and data services

* Hu et al. (2016), *Sun-to-Earth Characteristics of the 2012 July 12 Coronal Mass
  Ejection and Associated Geo-effectiveness*, ApJ 829, 97:
  https://doi.org/10.3847/0004-637X/829/2/97,
  https://arxiv.org/abs/1607.06287.
* NASA CDAWeb dataset/DOI registry:
  https://cdaweb.gsfc.nasa.gov/registry/hdp/NumericalData.xql.
  Wind MFI: https://doi.org/10.48322/av38-wn55;
  Wind definitive proton SWE: https://doi.org/10.48322/nasd-j276.
* STEREO Science Center: https://stereo-ssc.nascom.nasa.gov/.
* HELCATS HIGeoCAT: https://www.helcats-fp7.eu/catalogues/wp3_cat.html;
  geometric-catalog data: https://doi.org/10.6084/m9.figshare.5803176.v1.
* Heliospheric Shock Database: https://ipshocks.helsinki.fi/.
* Mishra et al. (2020), *Probing the Thermodynamic State of a CME Up to 1 AU*,
  April 2010 event and ambient-wind limitations:
  https://doi.org/10.3389/fspas.2020.00001.
* Temmer and Nitta (2015), 23 July 2012 propagation/preconditioning:
  https://arxiv.org/abs/1411.6559.
