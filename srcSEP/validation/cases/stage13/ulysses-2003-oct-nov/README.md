# ulysses-2003-oct-nov

This is an observational reference package. It supplies authentic available
measurements, source provenance and explicit gaps. It does not contain a
production model comparison or an observational PASS. `case.json` is the
machine-readable contract; every artifact has a byte count and SHA256.

The package retains61 actual daily COSPIN/HET pulse-height-analysis flux
records for1October–30November2003, including39–70 and71–94MeV protons,
helium channels, statistical errors and accumulation coverage. Raw11-column
ASCII and numerical format metadata are retained. The original fractional-year
and fractional-DOY times are accumulation midpoints; values of-1 remain invalid.
No HET counting rate is converted into a physical differential flux here.

[McKibben et al.](https://doi.org/10.1029/2005JA011049) studies multiple eruptions
at approximately5.2AU and6degrees north heliographic latitude. Daily averages
cannot resolve its directional flux-tube filling or exclusion signatures.
Those require higher-cadence sector data, responses and magnetic/plasma
histories. The archive documents statistical errors separately from estimated
10–20percent systematics; this package invents no covariance or error floor.

The package leaves model/event forcing, outer/return closure, response handling
and numerical acceptance unfilled until reviewed. The approved0.30–1.20AU
upstream profile must not be extrapolated to a Voyager or Ulysses observer.

Audit the local immutable assets:

```sh
python3 reproduce.py
```

Re-fetch bounded original sources and verify/reproduce archive subsets into a
fresh directory, without changing any reference:

```sh
python3 reproduce.py --fetch --output /tmp/ulysses-2003-oct-nov-reference-reproduction
```

The fetch command requires curl, authenticates source bytes and fails if an
archive has changed. Publication constraint tables were transcribed from the
cited numerical text/labels; source document hashes are in `retrieval.json`.
They have no automatic graph-digitization path. No copyright PDF or downloaded
mission-wide data archive is vendored. Reference integrity success leaves the
comparison status `NOT_RUN`; missing data remain an explicit gate.

When network access is restricted in a child process, supply independently
retrieved original files as `source-0.bin`, `source-1.bin`, etc., in the
`retrieval.json` order:

```sh
python3 reproduce.py --source-cache /path/to/original-bytes --output /tmp/fresh-reference-copy
```

The same original SHA256 checks apply; the cache never permits altered sources.
The initial direct curl retrieval succeeded in this session. The subsequent
Python-child live fetch failed DNS and is not reported as a successful live
reproduction; source-cache extraction was tested separately.

The additional [thermodynamic input audit](THERMODYNAMIC_QUALIFICATION.md)
documents 1464 authentic hourly plasma/field/position records and 401 proton
tensor records. They remain context, with missing electron thermodynamics and
coronal mapping stated explicitly. No old observation, threshold or comparison
verdict was changed.
