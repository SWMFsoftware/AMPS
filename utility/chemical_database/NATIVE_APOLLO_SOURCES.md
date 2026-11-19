# Apollo water-desorption native-source access

## Contract frozen before implementation — 2026-11-19

This recipe exposes the author-deposited files supporting Jones et al. (2020),
[Investigation of Water Interactions With Apollo Lunar Regolith Grains](https://doi.org/10.1029/2019JE006147),
from fixed [Zenodo record 3822365](https://zenodo.org/records/3822365).
The native record declares CC BY 4.0. Permission admission remains separate
from parsing: the archive packet pins the actual native licence metadata.
No new physical model, synthetic observation, U/I identifier, numerical
acceptance criterion or validation claim is introduced.

`native_apollo_tpd.reproduce(source_dir, output_dir)` authenticates the complete
fourteen-file collection using maintained expected byte counts, producer MD5
and SHA-256 before creating any output. It copies every original byte into
`native/`. The thirteen data files additionally receive one reversible lexical
CSV each in `tables/`; the native `readme.txt` remains source metadata.
The independent oracle is the authenticated author file itself. This is a
source-integrity/lexical-format test, not an AMPS or adsorption-law comparison.

The explicit production CSV schema is:

```text
source_line_1based,source_field_1based,line_role,native_token,native_line_ending
```

Each tab-separated source field becomes one CSV row. Empty fields and blank
lines are retained; the source's first line is labelled `native-header` and
is included verbatim. `native_line_ending` is `LF`, `CRLF`, `CR` or `NONE`.
Indices and line-ending labels are processing metadata, not new science
columns. An independent reconstruction must reproduce every original byte,
including trailing empty fields, signed decimal strings and source endings.
There is no float parsing, invented missing header, source-series alignment,
unit conversion, clipping, smoothing, interpolation or normalization.

In particular, `figure3.txt` has two tab-separated header fields but three
numeric fields in its data lines. It remains that way. Tables with independent
measured/model temperature axes retain independent token positions. Source
error bounds, `2SD` labels, negative background-subtracted signals and blank
cells are preserved without assigning a new confidence distribution.
A source `Probability` label does not authorize treating its values as an
independently measured, unit-normalized probability mass function.

## Physical scope and units

The material cases are Apollo highland sample 14163 and mare sample 10084.
The native files include measurements, background-subtracted signals,
model-inverted activation energies/site distributions, fitted model
components and lunar-cycle predictions. Each output records this distinction;
model predictions never become observations. Native temperature is K and the
record describes desorption in monolayers per second. Energy, coverage and
probability fields whose individual headers omit units remain unresolved
at the lexical API; an application needs the actual figure/method definitions
before interpreting those fields. This tool does not guess units from
numerical magnitudes or assign a lunar epoch, coordinate frame or geometry.
Time/state conventions remain those of each native experimental/model case.

The article's inversion assumes a prefactor of 10^13 s^-1 and its model uses
unit sticking at 180 K. These assumptions are distinct from a measured
universal sticking law. The earlier Jones et al. (2018) water-cycle model
uses isotropic ground-state/unit-yield fragments, unit H2O/OH sticking,
10^13 s^-1 and 55 kJ/mol desorption on a spherical, homogeneous Moon.
Schörghofer (2023) instead supplies an idealized efficient-diffusion,
heterogeneous-site end member with illustrative Ep=0.65 eV, W=0.22 eV and
10^16 s^-1. These are bibliographic model descriptions, not executable
providers installed by this parser. Their different inversion/prefactor,
coverage, surface and state conventions must not be silently combined.

## Verification contract

Administrative lexical fixtures test only bytes/formatting, never physical
observations. The authoritative-source test class is marked
`requires_authoritative_source = True`. Its checks exercise this maintained
recipe against the real fixed deposit and independently reconstruct all
thirteen files; check source bytes/mtime, measured/inverted/model role
separation, preserved missing headers/empty/signed values, deterministic
replay, existing-output refusal, corruption/truncation/missing-member refusal
before output, and explicit missing-source CLI status. Corrupt input,
unsupported schema or runtime failure is ERROR; unavailable source is SKIPPED.
An actually evaluated mismatch is FAIL. PASS establishes only the declared
source/format criterion. No predecessor or physical validation is implied.

## Usage and artifacts

From an AMPS checkout:

```sh
python3 -B utility/chemical_database/native_apollo_tpd.py \
  --source-dir /path/to/authenticated-fourteen-native-files \
  --output-dir /path/to/NEW-reproduction
python3 -B -m unittest discover -s test/chemical_database \
  -p test_native_apollo_tpd.py -v
```

Set `AMPS_APOLLO_TPD_SOURCE_DIR` to select a retained authoritative collection
for tests. Missing source skips only that class; pure lexical fixtures still
run. Existing output directories are refused and never overwritten.
The CLI reports PASS/ERROR/SKIPPED and exits 0/2/77. It performs no network
access, AMPS configuration, native application build, external publication or
canonical-data write. Outputs include the fourteen exact native files,
thirteen long-form tables, per-file source/output hashes and roles in
`reproduction.json`. A portable installed database resolves its archived
objects; the original research paths are provenance, not runtime dependencies.
