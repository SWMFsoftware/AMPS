# sep-2001-09-24

This is an observational reference package. It supplies authentic available
measurements, source provenance and explicit gaps. It does not contain a
production model comparison or an observational PASS. `case.json` is the
machine-readable contract; every artifact has a byte count and SHA256.

The package retains 360 original one-minute HED intensity/count records from
10:00–16:00 UTC and their HAPI parameter metadata. The nominal channel centers
are not the exact angular-response bands in the publication. Blank status
fields remain blank. Published path-length and pitch-angle constraints are
separate, labelled derived observations; no plot was digitized.

The return signature requires directional 17–22 and 28–54 MeV measurements,
uncertainties, acceptance and the magnetic-direction reconstruction described
by [Tan et al.](https://doi.org/10.1088/0004-637X/701/2/1753). Omnidirectional
context cannot establish a reflected beam. The quoted peak cosine is a
magnitude; this package does not infer its signed model convention. The fitted
mirror/path estimates must not be prescribed as forcing and then scored as an
independent prediction. The [mission event report](https://soho.nascom.nasa.gov/hotshots/2001_09_24/)
provides additional flare/shock context.

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
python3 reproduce.py --fetch --output /tmp/sep-2001-09-24-reference-reproduction
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
