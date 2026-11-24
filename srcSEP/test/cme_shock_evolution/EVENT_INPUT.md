# Prepare a new CME shock event

Prepare calibrated front picks and a dated upstream plasma model. The
`--data` path selects the observation-event JSON. Python writes a separate
fitted event for the linked srcSEP diagnostic.

```sh
/usr/bin/python3 srcSEP/test/cme_shock_evolution/regression.py prepare /path/to/run/prepared.event \
  --data /path/to/event/observations.json --issue 2012-05-17T02:45:00Z

/path/to/amps --test CME \
  --test-input /path/to/run/prepared.event --test-output-dir /path/to/run/shock-trace
```

Use your event's issue time and a fresh output directory. The date above is
an example. `--test CME` on its own reads `prepared.event` and writes
`shock-trace/` in the working directory.

The [29 November 2020 observation package](observed_2020_11_29/README.md)
walks through image retrieval, front selection and HCI packing using actual
FITS exposures. Use the same procedure for other events and substitute your
packed JSON path. The archived event exercises rejection at the current
spherical fit boundary.

The [manufactured observations.json](observations.json) is a complete
runnable example of the native five-hour path in the legacy Table.
Run it from the AMPS root.

```sh
/usr/bin/python3 srcSEP/test/cme_shock_evolution/regression.py prepare /tmp/cme-demo/prepared.event \
  --data srcSEP/test/cme_shock_evolution/observations.json --issue 2026-10-06T03:10:00Z
```

The `fixture` subcommand generates a separate editable synthetic example.

## Observation-event JSON

The top-level object uses `shock-front-replay-v1`. The standalone replay
also accepts this file. Use JSON numbers for scalars and vector components.
The units of `initial_fit` are stated separately below.

| Key | Type, units and meaning |
| --- | --- |
| `schema` | Required string `shock-front-replay-v1`. |
| `frame` | Required string `HCI`. Use heliocentric inertial Cartesian coordinates with the Sun as origin and +Z along the declared solar rotation axis. Use one common orientation for every instrument and vector. Record the HCI realization, reference equinox and transformation version in `provenance`. |
| `reference_utc` | Required UTC origin for the elapsed-time coordinate and fit. Native form is whole seconds, `YYYY-MM-DDTHH:MM:SSZ`. Choose an origin at or before the first retained front epoch. |
| `light_time_policy` | Required nonempty string documenting the convention used to assign a physical front epoch from the observation time and geometry. Apply it before writing `elapsed_s`. |
| `provenance` | Required nonempty string identifying the event, images or catalog, ephemeris, astrometric calibration and feature-selection procedure. Include product identifiers or a source manifest. |
| `tracks` | Required array of calibrated front-pick objects. After the timestamp and camera filters, provide at least eight picks at three distinct epochs. Their angular coverage must identify the fitted geometry. |
| `holdout_tracks` | Optional array of reserved observations in the same record format. The archive packer writes exposures designated as holdout here. The fitter reads `tracks`. |
| `initial_fit` | Required array `[r0, v0, a0, lon, lat, half_width]`. The first three values are starting estimates for apex radius, speed and acceleration at the reference epoch, in solar radii, solar radii/hour and solar radii/hour². The remaining values are HCI longitude, latitude and angular half-width in radians. Longitude is measured from +X toward +Y and latitude toward +Z. Half-width is between 0 and π/2. The nominal solar radius is 695,700,000 m. |
| `fit_acceleration` | Required boolean. `true` fits constant acceleration over the observed interval. `false` uses constant speed and zero acceleration. The fit uses a single axis and width over the interval. |
| `seed` | Required nonnegative integer random seed used by the shared fitter. The native receipt carries the nominal fit. |
| `ensemble_samples` | Required integer from 1 to 1024. This controls sampling in the shared fitter. The standalone replay reports ensemble statistics. |
| `ambient_members` | Array of upstream model records below. A native trace requires an eligible member. |
| `rays_hci` | Optional nonempty array of unit HCI vectors from Sun center, identified by array order. Each ray must intersect the leading cap. The default is the fitted nose and 48 directions on four rings. |
| `propagation` | Object specifying the outer drag model. Required for a trace past the last retained fitted epoch. |
| `evolution_elapsed_s` | Optional array of snapshot times for the standalone replay. The native diagnostic schedules its own 1200 s cadence. |
| `observers_hci_m` | Optional map of observer name to fixed HCI position in metres for the standalone arrival calculation. |

### Each record in tracks

| Key | Type, units and meaning |
| --- | --- |
| `view` | String identifying the calibrated camera or viewpoint, for example `SOHO/LASCO-C2`. Repeatable `--withhold-view NAME` matches this value exactly. |
| `frame_file` | Optional local archive alias linking the pick to a specific exposure in the source manifest. |
| `observed_utc` | ISO 8601 UTC exposure or acquisition timestamp ending in `Z`. Use the timestamp associated with the calibrated camera geometry. |
| `available_utc` | Earliest UTC when the calibrated pick could have been used. It is at or after acquisition. A later annotation or processing date gives a retrospective replay cutoff. A pick is eligible when it was available by `--issue`. |
| `elapsed_s` | Elapsed SI seconds from `reference_utc` to the physical front epoch after the chosen light-time correction. It is between zero and the elapsed acquisition time. Fractional seconds are supported. |
| `observer_hci_m` | Three finite Cartesian HCI components of camera position at acquisition, in metres from Sun center. |
| `los_hci` | Three finite components of a unit HCI vector **from the camera toward the picked outer front edge**. Obtain it from calibrated image astrometry and observer geometry. |
| `sigma_rad` | Positive, finite one-standard-deviation angular uncertainty transverse to the apparent front edge, in radians. |

An image can supply several picks at the same epoch. Preserve gaps and
annotate the outer front as a distinct feature from the bright ejecta. The
fit interprets the picks as leading-cap tangencies. Image calibration,
coordinate conversion and selecting that feature precede the fitter. The
elapsed-time contract requires intervals without a leap second.

The first pick in the manufactured example shows the record format. The
camera is on the +X axis and the sight line points back toward the Sun.

```json
{
  "view": "A",
  "observed_utc": "2026-10-06T00:08:20Z",
  "available_utc": "2026-10-06T00:09:05Z",
  "elapsed_s": 0.0,
  "observer_hci_m": [149597870700.0, 0.0, 0.0],
  "los_hci": [-0.9997396826370526, -0.012733910766339725, 0.018931837664635486],
  "sigma_rad": 0.000001
}
```

### Each record in ambient_members

| Key | Type and meaning |
| --- | --- |
| `event_file` | Path to a canonical AMPS reduced-shock upstream `.event` file. A relative path is resolved against the directory of the JSON supplied to `--data`. The upstream reader resolves and checks its assets relative to the `.event` file. |
| `available_utc` | ISO 8601 UTC timestamp ending in `Z` when this realization was available. It must be no later than `--issue` to be eligible. |

If several members are eligible, choose `--ambient-member INDEX` using its
zero-based index in the original array. The native diagnostic uses one
coherent upstream realization. Its UTC origin comes from
`run.reference_epoch` in the upstream file and may differ from the
observation origin. AMPS aligns the two clocks, checks time and spatial
coverage and queries density, total thermal pressure, wind and magnetic
field at every sampled front point. The fitted trajectory supplies the front
motion. The [upstream example](../../examples/reduced-shock/positive_1au.event)
and [core contract](../../../src/models/sep_corona_swcme/shock_front/README.md)
give the upstream file fields, model assumptions, units and assets.

### Fields in propagation

| Key | Units and meaning |
| --- | --- |
| `wind_m_s` | Positive effective solar-wind speed in m/s for the drag law. |
| `drag_m_inv` | Nonnegative quadratic-drag coefficient in inverse metres. |
| `horizon_s` | Positive propagation time from the last fitted epoch in seconds, at most 14 days. |
| `minimum_handoff_apex_m` | Positive minimum heliocentric apex distance in metres for admission of the drag model. |

By default the five-hour trace starts at the first retained fitted epoch.
Use `--start-elapsed-s` for a later start. The full trace must be inside the
fitted interval or within the declared propagation coverage. The drag
history matches apex radius and speed at the last fitted epoch.
Unsupported time windows and handoffs fail before the fitted receipt is
published.

## Prepared native event

Python writes the following ordered text contract. The native diagnostic
loads it and the selected upstream file. Keep the original JSON with the
evidence.

| Line | Contents |
| --- | --- |
| 1 | Literal `cme-shock-evolution-v2`. |
| 2 | Whole-second model reference UTC. |
| 3 | Observation issue UTC. |
| 4 | Lowercase SHA-256 of the original JSON bytes. |
| 5 | Quoted absolute observation JSON path for provenance. |
| 6 | Quoted absolute path to the selected upstream `.event`. |
| 7 | First fitted epoch, last fitted epoch, and trace-start epoch in elapsed SI seconds. |
| 8 | Fitted apex radius (m), speed (m/s), acceleration (m/s²) at the reference, then the three unit HCI axis components and half-width (rad). |
| 9 | `fitted-only` if the whole trace fits within the observation span. Otherwise `quadratic-drag` followed by wind (m/s), drag (1/m), horizon (s) and minimum handoff radius (m). |
| 10 | Number of sample directions. |
| 11 onward | One whitespace-separated unit HCI vector per direction. |

C++ writes the initial map and fifteen successive 1200 s maps. Filenames
count from the **start of the trace**, so the last file is always
`shock-18000.csv`. The CSV `elapsed_s` counts from `reference_utc`. A
trace starting at 1800 s therefore ends at CSV `elapsed_s = 19800`. A
subfast or non-forward ray has blank downstream columns. Unavailable upstream
coverage, an unsupported ray or an unresolved numerical jump fails the
diagnostic with context.

The legacy Table reads the checked-in observation fixture and applies
`prepare --legacy` with a frozen observation selection and input contract.
Its five-hour comparator checks that manufactured case. A user event
produces a new trace to evaluate against its own evidence.
