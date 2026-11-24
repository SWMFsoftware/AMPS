# CME shock evolution regression

The legacy AMPS Table entry is `SEP--CME-shock-evolution` in the SEP family.
The standalone srcSEP case ID is `CME`. The descriptive CLI alias is
`--test cme-shock-evolution`.

From an installed AMPS checkout, regenerate the legacy targets after a Table
change and then run the case.

```sh
./Config.pl -install
make test_SEP--CME-shock-evolution TESTCOMPILEJOBS=16
```

The legacy run uses `/usr/bin/python3` with NumPy, SciPy, matplotlib and Pillow. Set
`AMPS_CME_PYTHON` to select another interpreter with these dependencies.
A failed admission check, native run, or numerical comparison propagates
through the Make target. A successful check leaves an empty
`test_SEP--CME-shock-evolution.diff`.

## Archived observation input

The Table runs both observation preparation scripts offline before the native
trace. Each script has its own frozen reference.

| Script | Regression output | Reference and check |
| --- | --- | --- |
| `preview_observations.py` | Three PNGs and three NPZ files in `preview/`. | [preview-reference.json](observed_2020_11_29/preview-reference.json) contains the expected displayed difference arrays, image origin, pixel extent, stretch and marks. The checker also decodes each PNG and verifies that an orange mark rendered. |
| `pack_observations.py` | `packed-observations.json` from all 59 archived marks. | [observations.json](observed_2020_11_29/observations.json) must match byte for byte before the generated file goes to the fitter. |

The packer input is from the 29 November 2020 CME. The shared fitter reaches the upper half-width
bound of the reduced spherical front and rejects the event. The expected
admission result is written to `observed-2020-11-29.json` in the run directory.
The [archival preparation guide](observed_2020_11_29/README.md) gives the
image inventory, feature selection, time and coordinate transforms, packed
input and commands for preparing another event.

## Preparation and trace

For your own event, pass the calibrated JSON path with `prepare --data`.
The [event input guide](EVENT_INPUT.md) specifies every field, its units and
time convention, the upstream input and the generated native format.

For the numerical benchmark the Table passes the checked-in
[manufactured observations.json](observations.json) directly
to Python as `--data`. It contains synthetic coronagraph picks derived from
an analytically known moving shock and three simulated cameras. The
[manufactured data generator](../../../src/models/sep_corona_swcme/forecast/manufactured_event.py)
documents that construction. The issue cutoff admits 332 exact picks and
excludes one late synthetic pick. The fitted coefficients and 49 HCI
sampling rays are written to a versioned `prepared.event`. The ambient input is the checked-in
[synthetic dipole PFSS/Parker event](../../examples/reduced-shock/positive_1au.event).
The fixture is self-contained and has exact analytical truth. It verifies
reconstruction and native shock evaluation.

The linked C++ test loads the fitted event and evaluates the finite SSE
surface, upstream state and local Rankine–Hugoniot jump. Shock evaluation
uses the fitted curve through 3 h and the configured quadratic-drag law
through 5 h. At handoff the apex position and speed are continuous. The
clock produces an initial map and fifteen completed 1200 s intervals.
Rank zero writes the CSV and completion manifest at every epoch. The cadence
is a shock-evaluation schedule for this regression.

The two run steps can also be invoked separately. Run the Python command from
the repository root. Run the executable from the directory where you want
the snapshots, with the prepared event in that directory.

```sh
/usr/bin/python3 srcSEP/test/cme_shock_evolution/regression.py prepare /path/to/run/prepared.event \
  --data srcSEP/test/cme_shock_evolution/observations.json --issue 2026-10-06T03:10:00Z --legacy
cd /path/to/run
/path/to/amps --test CME
```

The default output directory is `shock-trace`. To select explicit paths, pass
`--test-input PATH --test-output-dir DIRECTORY` together. Use a fresh output
directory for each run. Python handles preparation and comparison. The
linked executable reads the prepared event and computes the trace with the
existing C++ shock kernels in a standalone, background-only test.

## Reference and acceptance

`final-reference.csv` freezes the complete 49-ray snapshot at
`t = 18000 s`. The legacy checker verifies the manifest checksum, full
cadence, ray identity and ordering, and independent geometry and normal
speeds at **every** epoch. The geometry oracle uses the manufactured sphere
and integrates the drag ODE independently. It allows absolute errors of
2 m in position, \(2\times10^{-9}\) in normal components, and
0.005 m/s in normal speed. Solved samples must remain on the compressive,
superfast branch and satisfy the existing \(10^{-8}\) jump residual
bound. The final comparison checks exact classifications and uses relative
tolerance \(2\times10^{-8}\) and absolute tolerance \(2\times10^{-9}\)
for the frozen numeric fields. RH residuals are checked against their physical
bound at every epoch.

The initial reference came from the linked native trace after verification
against the manufactured geometry and the independent drag integration. Its
final Mach number, obliquity, compression and statuses were also cross-checked
with the standalone Python/C++ replay. There are 37 fast-shock rays and 12
subfast rays at the final epoch. The SHA-256 is
`d948bf93b10ba5f1c0059dfd801a1775a13407131e899f6d453e876619157454`.
Refresh this reference only with an explicitly reviewed physics change.
Rerun the trace in a clean directory, solve the independent checks via
`regression.py check DIRECTORY --reference CANDIDATE.csv`, and review the
numerical diff and provenance before replacing the frozen file.
