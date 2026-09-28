# `sep_coronal_cme` model specification sources

This README is a **non-normative navigation and maintenance guide**. The
generated [`model.md`](model.md) is the complete public review artifact; its
normative sources are the four modules named below plus the structured trace
registry. Physics, configuration defaults, equations, acceptance limits, API
semantics, and test definitions must not be introduced in this README.

## Why the document is generated

The model is shared by `srcSEP3D` and field-aligned `srcSEP`. Keeping one large
file as an independently edited copy would permit an equation, input enum, API
record, or release test to change in one place without the related authorities
changing with it. The generated layout makes those ownership boundaries
executable: each numbered section has exactly one maintained owner, while the
checked-in `model.md` remains convenient for review and download.

## Maintained layout

```text
sep_coronal_cme/
  README.md                         non-normative maintenance guide
  model.md                          generated canonical review artifact
  model/
    physics.md                      preamble; Sections 1-2, 4-13, 18-19
    configuration_validation.md     Section 14
    architecture_exchange.md        Sections 3 and 15
    testing_validation.md           Sections 16-17 and 20
    requirements.yaml               requirements, review dispositions, links
  tools/
    generate_model.py               deterministic generator and validator
  test/
    test_model_documentation.py     DOCSCCM01 positive and negative tests
```

The modules own noncontiguous section ranges, so simple file concatenation is
incorrect. `generate_model.py` parses actual level-two numeric headings while
ignoring headings inside fenced configuration examples, validates the declared
semantic owner, and emits Sections 1 through 20 in numeric order.

Revision `2026-09-28-r4` keeps application input schema 5 and field-line bundle
schema 3 unchanged. It corrects current diagnostic and interpretation
contracts, but the momentum-dependent release surface, return-renewal,
self-generated-wave, smooth-field drift, dynamic CME-attitude, impulsive-source,
imported-MHD-provider, and continuous wind-envelope interfaces in Section 15.8
are post-schema-5 roadmap contracts, not implemented selectors. A schema-5
deck that attempts to select them must fail as unsupported.

`requirements.yaml` deliberately uses JSON syntax. JSON is a strict subset of
YAML 1.2, which keeps the requested YAML artifact while allowing documentation
generation on a standard Python installation with no PyYAML package. It owns
the stable `SCCM-R1`--`SCCM-R10`, `SCCM-N1`--`SCCM-N7`, and
`SCCM-P1`--`SCCM-P8` records, review rows, section links, configuration-key
links, API-record links, and canonical test links. The `P` series records the
fourth technical review: its accepted corrections are therefore checked by
the same traceability gate as the earlier reviews instead of existing only as
free-form prose. The review-disposition and stable-requirement tables in
`model.md` are rendered from this registry; they are not duplicated in the
prose modules.

## Editing workflow

1. Edit only the module that owns the affected numbered section.
2. Update `model/requirements.yaml` when a review disposition, stable
   requirement, section/configuration/API link, or associated test changes.
3. Regenerate `model.md` from the model directory or any other working
   directory:

   ```sh
   python3 src/models/sep_coronal_cme/tools/generate_model.py
   ```

4. Run the byte-identity gate and its negative tests:

   ```sh
   python3 src/models/sep_coronal_cme/tools/generate_model.py --check
   python3 -B -m unittest discover \
     -s src/models/sep_coronal_cme/test \
     -p 'test_*.py'
   ```

5. Review both the maintained source diff and the regenerated `model.md` diff.
   A direct edit to `model.md` is intentionally rejected by `--check`.

For review tooling that provides Pandoc, an additional non-authoritative
Markdown parse smoke test is:

```sh
pandoc --from=gfm --to=plain \
  src/models/sep_coronal_cme/model.md \
  -o /dev/null
```

## Generator behavior

The generator is intentionally deterministic and repository-relative:

- reads UTF-8 without a byte-order mark;
- rejects carriage returns and missing final newlines instead of normalizing
  them silently;
- accepts only the four declared module paths and semantic ownership sets;
- recognizes only the two registered generated-table placeholders;
- checks stable requirement/review IDs and every declared section,
  configuration key, API record, and test link;
- requires canonical test IDs to have one definition in Section 17.2;
- checks balanced Markdown fences and sequential bibliography numbering;
- includes no timestamp, host path, environment value, or random ordering in
  generated content; and
- publishes output through an atomic same-directory replacement.

`--output PATH` writes the validated rendering to another path without
changing the canonical artifact. `--check` renders in memory and reports a
concise unified diff if committed `model.md` is stale.

## Scope boundary

This directory documents and will host the dependency-light shared coronal-CME
physics library. `srcSEP3D` links it through a thin AMPS-facing adapter.
Field-aligned `srcSEP` does not link the library; it consumes the immutable,
versioned field-line bundle through neutral `sep_common` records. The sibling
`swcme` model is an alternative backend, not a dependency of this stand-alone
model. Application, AMPS/PIC mesh, MPI, and Tecplot concerns remain outside the
shared physics layer.

## Implemented staged build

The implementation follows the roadmap as cumulative hard gates. Run
`make test-stage0` through `make test-stage6`; a later gate
always includes every earlier test. Individual canonical launchers live under
`test/individual/ID/`, and JSON/JUnit evidence is written under `build/` so it
cannot leak into a source archive.

`make test` (or `python3 test/run_tests.py --all`) is the single aggregate
Stage 0--6 gate.

Detailed implementation notes live in `docs/`. They explain the physics-to-
code mapping while leaving `model.md` as the normative equation and acceptance
authority. Stage 0 is described in
[`docs/STAGE0_CONFIGURATION_AND_ARCHITECTURE.md`](docs/STAGE0_CONFIGURATION_AND_ARCHITECTURE.md).
Stage 1 is described in [`docs/STAGE1_PFSS.md`](docs/STAGE1_PFSS.md).
Stage 2 is described in
[`docs/STAGE2_WIND_CLOSED_PLASMA_AND_INTERFACES.md`](docs/STAGE2_WIND_CLOSED_PLASMA_AND_INTERFACES.md).
Stage 3 is described in
[`docs/STAGE3_SOURCE_SURFACE_COUPLING.md`](docs/STAGE3_SOURCE_SURFACE_COUPLING.md).
Stage 4 is described in
[`docs/STAGE4_TURBULENCE_AND_SCATTERING.md`](docs/STAGE4_TURBULENCE_AND_SCATTERING.md).
Stage 5 is described in
[`docs/STAGE5_ELLIPSOID_AND_PISTON.md`](docs/STAGE5_ELLIPSOID_AND_PISTON.md).
Stage 6 is described in
[`docs/STAGE6_OBLIQUE_MHD_JUMP.md`](docs/STAGE6_OBLIQUE_MHD_JUMP.md).
