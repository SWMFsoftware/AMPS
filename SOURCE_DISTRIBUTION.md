# SEP source-distribution contract

This archive is a source baseline for the independent `srcSEP` and `srcSEP3D`
applications and their canonical shared models.  It is not a build directory.
The only shared implementation owners are:

- `src/models/sep_common` for provider-neutral SEP kernels; and
- `src/models/swcme` for the solar-wind/CME/shock model.

The applications remain peers.  Neither application may include files from,
inspect, or invoke the other application.  Cross-application release checks
belong in `tools/`, outside both application trees.

## Current srcSEP milestone

This revision retains the B01-B05 clean-source baseline and adds improvements
D01-D03 to `srcSEP`:

- D01 makes SWCME preparation and sampling transactional and fail-closed, with
  explicit recovery policies and rank-reduced diagnostics;
- D02 centralizes the complete 1-D SWCME/shock/source key and unit schema under
  `src/models/swcme`, with deterministic normalized fingerprints; and
- D03 adds the configured AMPS build plus serial/MPI/refresh/restart campaign
  gate and records proof that each selected production mover actually returned
  from particle dispatch.

Dependency-light D01/D02 and D03-orchestrator tests are packaged. A D03 native
`PASS` still requires a configured target AMPS tree, linked executable, MPI
launcher, and reviewed site campaign; the source package never substitutes a
synthetic transport run for that evidence.

The linked srcSEP C++ registry now owns extended IDs `D01`, `D02`, and
`D03PRE`. The first two reuse the exact dependency-light callbacks; the third
checks linked mover and SWCME-refresh evidence hooks without claiming the outer
MPI/restart campaign. Python `--all` discovers all three from `--list-tests`.
Its optional `--rebuild` performs the enclosing clean/strict production build
before discovery, so a stale executable cannot hide newly registered tests.

## Authoritative manifests

Each maintained component contains `SOURCE_MANIFEST.json`.  A manifest names
the source, test, documentation, retained-input, and generated-path classes for
that component.  The manifests are release inputs: changing a source layout
requires changing the applicable manifest in the same review.

Run the complete package gate from the AMPS root:

```sh
python3 tools/sep_package_hygiene.py --root . --self-test
```

The command first checks the real tree and then runs a synthetic negative
control.  The negative control inserts a stale object into a temporary package
and must be rejected.  The real-tree check fails on:

- retired source names listed by a component manifest;
- object/dependency/archive files or native executables;
- `build`, `test_output`, or `__pycache__` trees;
- AMPS runtime binary data such as `amr.sig=*.bin`;
- source files that are not classified by exactly one manifest category; or
- a classified file whose manifest pattern no longer matches an existing file.

The complementary shared-model ownership scan is:

```sh
python3 tools/check_swcme_ownership.py --root . --self-test
```

It requires all SWCME implementation namespaces, public headers, production
source, and demonstrations to remain under `src/models/swcme`; applications
may contain only provider adapters that include bare canonical header names.
The negative control installs a disguised second implementation and proves
that the scan rejects it. Application runners retain their independence and
inspect only their own tree plus the canonical model; this cross-application
scan belongs exclusively to the AMPS-level release procedure.

Small, reviewed JSON/CSV/text reference inputs remain permitted when their
component manifest classifies them as retained inputs.  A filename extension
alone never exempts a native executable: the gate checks its magic bytes.

## Reproducible packaging

After the hygiene gate passes, create a deterministic source package with:

```sh
python3 tools/create_sep_source_package.py \
  --root . --output ../AMPS-SEP-B01-B05-D01-D03-source.tar.gz
```

The packager reruns hygiene, excludes only the generated paths declared by the
manifests plus repository metadata, sorts archive members, and normalizes
ownership and timestamps.  It never packages an existing output archive back
into itself.  Build products must be regenerated after extraction.

## Generated-path policy

Ignore rules are intentionally local to each component.  They name known
products rather than broad data extensions, because validation CSV/JSON inputs
and frozen text records are source-controlled evidence.  Build and test
commands may create ignored paths locally, but those paths must be absent from
a release archive.
