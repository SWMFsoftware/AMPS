# sep-1978-09-23

This is an observational reference package. It supplies authentic available
measurements, source provenance and explicit gaps. It does not contain a
production model comparison or an observational PASS. `case.json` is the
machine-readable contract; every artifact has a byte count and SHA256.

[Reames, Ng and Tylka](https://doi.org/10.1007/s11207-012-0038-1),
Section2.4/Figure7, identifies Voyager1/2 at approximately4.2/3.9AU and
reservoir-like22–27MeV proton intensities by28September. The local table
records these published constraints, not invented intensity samples.

The actually retrieved Voyager2 daily and6-hour products contain only
22–23September in the requested22September–1October window. They omit
24–30September. `archive_coverage.json` retains the audit and bracketing
available dates. The excerpts preserve all original channel and error columns.
The wider17.9–30.0MeV archive bin is not identical to the published22–27MeV
channel. Neither gap filling nor a channel substitution is authorized.

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
python3 reproduce.py --fetch --output /tmp/sep-1978-09-23-reference-reproduction
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
