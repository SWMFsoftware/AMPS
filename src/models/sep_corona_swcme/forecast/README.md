# Observation-driven shock-front replay

This is the standalone archived-event prototype. It fits a finite self-similar
expansion (SSE) front to calibrated white-light sight lines. It outputs
the time evolution of the HCI front and its local shock properties for SEP
studies. The optional observer-arrival calculation uses the same history.
The prototype links the existing AMPS finite SSE, canonical PFSS/Parker
upstream, quadratic drag and
fast Rankine–Hugoniot implementations through `core_cli.cpp`. It runs without
PIC, MPI, a particle model or a production AMPS timestamp.

## Build and run

From the AMPS repository root, use C++17 and Python 3.8 or later with NumPy and
SciPy. The published-figure replay below additionally uses Pillow.

```sh
make -C src/models/sep_corona_swcme -j4 forecast-test shock-front-test PYTHON=/usr/bin/python3

/usr/bin/python3 src/models/sep_corona_swcme/forecast/replay.py /path/to/event.json \
  --issue 2012-05-17T02:45:00Z \
  --core src/models/sep_corona_swcme/build/shock_forecast_core \
  --output /path/to/replay.json
```

Repeat `--withhold-view NAME` to test the sensitivity to a specific camera.
Run the same input at earlier issue times to inspect changing coverage or
loss of identifiability. The output is written atomically after input and
core validation. A failed fit exits with a diagnostic.

## Input contract

The JSON fields below define `shock-front-replay-v1`.
The [CME event input guide](../../../../srcSEP/test/cme_shock_evolution/EVENT_INPUT.md)
describes the full observation data contract and preparing a five-hour native
srcSEP diagnostic from the same file.

| Field | Meaning |
| --- | --- |
| `schema`, `frame` | `shock-front-replay-v1`, `HCI` |
| `reference_utc` | UTC origin, ISO 8601 ending in `Z`. Replay intervals must not span a leap second. |
| `light_time_policy`, `provenance` | Describe the camera calibration, observer ephemeris, selected feature, and conversion from observation time to front time. |
| `tracks` | One record per calibrated front pick, with fields below. |
| `initial_fit` | Six starting values `[apex_radius_Rsun, apex_speed_Rsun_h, apex_accel_Rsun_h2, longitude_rad, latitude_rad, half_width_rad]`, with kinematics defined at the reference epoch. Angles stay constant over the fitted interval. |
| `fit_acceleration` | `true` for quadratic apex distance, `false` for constant speed. |
| `seed`, `ensemble_samples` | NumPy random seed and number of fit draws, 1–1024. |
| `rays_hci` | Optional unit HCI vectors to sample. The default includes the nose and 48 rays on four rings covering the fitted cap. |
| `evolution_elapsed_s` | Optional nonempty list of output epochs in elapsed SI seconds from `reference_utc`. The default is every usable front-pick epoch. Intermediate times inside the fitted interval use the fitted trajectory. Later times require `propagation` and must lie within its declared horizon. |
| `ambient_members` | Optional list of `{event_file, available_utc}` records. Paths are relative to the replay JSON or absolute. Eligible members have equal weight and were available by the replay issue time. The core verifies each member's assets and supplies its UTC validity origin and physics identity. The replay JSON supplies the trajectory, observer and drag settings. |
| `observers_hci_m` | Optional mapping of observer name to a **fixed** HCI position in metres, above the Sun and at most 1 AU. |
| `propagation` | Optional object containing `minimum_handoff_apex_m`, `wind_m_s`, `drag_m_inv`, and `horizon_s` from the last fitted epoch (maximum 14 days). Supply coefficients and minimum handoff radius supported by the event or by the population used for calibration. |

Each track record contains `view`, `observed_utc` (exposure/receive timestamp),
`available_utc` (the time the calibrated pick could have been consumed),
`elapsed_s` (front epoch in elapsed SI seconds from the origin),
`observer_hci_m` (camera position at acquisition), `los_hci` (unit vector from
the camera toward the picked feature), and `sigma_rad` (one-standard-deviation
angular uncertainty transverse to the apparent contour). Supply an
astrometrically calibrated sight line to the **outer front feature**, with
the time correction applied consistently. Use independently determined
availability timestamps for picks and upstream members in forecast
validation. A retrospective annotation or publication time is
sufficient for a reconstruction replay and has no operational lead-time claim.

## Shock evolution for SEP

Read `evolution` as a time series. Each entry gives its elapsed SI time,
apex distance and speed, and a map on the same HCI rays. Along each ray the
core evaluates the leading intersection, outward normal, local normal speed
and the upstream state at that time and position. The local MHD calculation
then reports shock classification, fast-mode Mach number, density compression
and field obliquity, with the existing validity and ensemble bounds. The ray
labels identify fixed heliocentric directions. Connecting a material field line
to the shock requires intersecting the surface with that line at the same epoch.
The top-level `map` still gives the latest fitted map for snapshot consumers.

Entries marked `fitted-history` use the observations available at the stated
issue time. Earlier points in that history describe the reconstruction at that
issue time. Run the same event with earlier issue times to reproduce how the
reconstruction itself changed as new views were received.

Request later epochs explicitly with `evolution_elapsed_s`. The
`quadratic-drag` entries use the supplied outer propagation model, joined at
the last fitted time with continuous apex position and speed. Each trajectory
realization keeps its own handoff state and must satisfy the declared minimum
handoff radius. Unsupported handoffs retain their weight as unresolved. A
request before the first fitted epoch or beyond the configured drag horizon
fails with a diagnostic. The model keeps the fitted axis and angular width
fixed and computes evolving speed and shock strength along the front.
Connecting the exported history to the SEP injection/source solver and
qualifying its predictive skill against observations are separate milestones.

## Physical and statistical bounds

The forward operator compares the angular distance from the pick to the
sphere axis with the exact finite-distance angular radius of the SSE
generating sphere. It checks that the corresponding tangent belongs to the
leading cap. This prototype models the picked outer contour as a tangent.
Brightness synthesis with Thomson scattering remains follow-on work. Poorly
conditioned fits are rejected even if the apparent contour residual is small.
The reported residual and number of views support review of a fixed-shape
fit. The ellipsoid fitting path and instrument-specific pick/WCS production
remain follow-on work.

Correlated fit draws come from the image-residual Jacobian and assume local
Gaussian, independent angular errors on the picks. Assess correlated
calibration and front-identification errors during event validation. The
draws are combined with the explicitly provided upstream
members. An invalid draw, unsupported upstream state, or unresolved local
jump remains in the denominator of the lower and upper fast-shock fraction
bounds. These fractions are exploratory and require calibration against
held-out observations. Reported Mach number, compression, and obliquity
quantiles are conditional on admissible shocks. Missing values are JSON
`null`. Geometric crossing times and times conditional on a fast-shock
solution are distinct. The optional forecast uses one supplied drag history,
matched at the last fit epoch in radius and speed. It marks crossings
unresolved when the fitted front is below the declared handoff radius or the
arrival lies outside the forecast window. Uncertainty in drag coefficients
and motion of the destination spacecraft are outside this prototype.

## 17 May 2012 observed apparent motion

[Pohjolainen, Two-Part Interplanetary Type II Solar Radio Bursts,
Figure 11](https://doi.org/10.1007/s11207-025-02449-5) separates the red COR2-A
white-light shock feature from the black CME points and the radio markers.
Section 4.1 reports approximately 1,100 km/s for that white-light feature.
The replay digitizes five centers along the authors' dashed guide. These
samples are correlated and carry no independent uncertainty on the 3-D
normal speed. The published plot supplies projected heights, without the
calibrated sight-line directions and multi-view pick times required by the
inverse model. The observational receipt therefore records a projected replay
and explicit skips for the 3-D fit, local shock strength, and arrival forecast.

```sh
curl -fL 'https://media.springernature.com/full/springer-static/image/art%3A10.1007%2Fs11207-025-02449-5/MediaObjects/11207_2025_2449_Fig11_HTML.png' \
  -o /tmp/sep-20120517-figure11.png
/usr/bin/python3 src/models/sep_corona_swcme/forecast/projected_fig11.py \
  /tmp/sep-20120517-figure11.png --output /tmp/sep-20120517-projected.json
```

The script checks SHA-256 before digitization. The source raster is retrieved
on demand. A 3-D observational gate needs calibrated picks and ephemerides
of the same outer front from [SOHO/LASCO](https://lasco-www.nrl.navy.mil/index.php?p=get_data)
and [STEREO/SECCHI](https://stereo-ssc.nascom.nasa.gov/data.shtml), with time
coverage comparable to the [Rouillard et al. reconstruction](https://arxiv.org/abs/1605.05208).
Compare it against independent observations after any fit parameters and
acceptance criteria are frozen. The observed event contains a pre-existing
magnetic cloud and a steady PFSS/Parker ambient is insufficient to assume a
qualified event-specific Mach number. The complete forecast and native SEP
qualifications remain open.

## Verification

`forecast-test` generates independent camera/sphere tangencies for a moving
surface. It checks noisy two-view recovery, an imaging gap, the exclusion of
an observation unavailable at issue time, rank rejection, the linked core's
normal speed and shock classifications, a manufactured 1-AU arrival,
atomic output on a rejected fit, time-resolved radius, velocity and local
shock properties, and forecast evolution checked against an independent
integration of the drag ODE. It also checks a nonzero-drag crossing against an
independent ODE solution, and explicit handling of late and out-of-range
upstream ensemble members and horizon-limited arrivals.
The existing `shock-front-test` independently covers the MHD jump,
conservation, geometry and shared drag laws. Synthetic evidence verifies
software and physics kernels. The Figure 11 receipt checks recovery of the
published projected motion.
