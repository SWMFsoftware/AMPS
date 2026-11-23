# sep-1979-03-01

This is an observational reference package. It supplies authentic available
measurements, source provenance and explicit gaps. It does not contain a
production model comparison or an observational PASS. `case.json` is the
machine-readable contract; every artifact has a byte count and SHA256.

[Reames, Ng and Tylka](https://doi.org/10.1007/s11207-012-0038-1),
Section2.5/Figure8, uses Voyager2 at approximately4.76AU, together with Helios1/2
and IMP8. Voyager1 was approaching Jupiter; it is not silently substituted.
Published observer distances and reservoir timing statements are retained as
constraints, not digitized intensities.

The retrieved NASA Voyager2 daily and6-hour flux archives contain ZERO rows in
the requested28February–10March window. Header-only excerpts intentionally
preserve the product schema; they are not observations or zero flux.
`archive_coverage.json` records the last earlier/first later available records.
A calibrated event product or reviewed figure digitization remains necessary.

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
python3 reproduce.py --fetch --output /tmp/sep-1979-03-01-reference-reproduction
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
