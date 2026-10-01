# Reproducible synthetic Stage-12 workflow

These are manufactured observations for software verification, with an
explicit synthetic HCI epoch `2000-01-01T00:00:00Z`. They are not a solar-event
calibration or a production background. Original source bytes, unit/frame/
support declarations, covariance, processing version and SHA-256 are specified
in `job.json`. Source generation is deterministic:

```sh
python3 tools/generate_stage12_example.py
```

Run from `src/models/sep_coronal_cme`, using a fresh output directory. The CLI
refuses to replace previously published assets or frozen records.

```sh
python3 tools/preprocess_observations.py preprocess --job examples/stage12/job.json --output build/stage12-example/assets
python3 tools/preprocess_observations.py verify --bundle build/stage12-example/assets
python3 tools/preprocess_observations.py preregister --bundle build/stage12-example/assets --input examples/stage12/preregistration.json --output build/stage12-example/frozen-preregistration.json
python3 tools/preprocess_observations.py select --bundle build/stage12-example/assets --preregistration build/stage12-example/frozen-preregistration.json --input examples/stage12/realized-candidates.json --output build/stage12-example/selection.json
python3 tools/preprocess_observations.py field-line-requests --bundle build/stage12-example/assets --ephemeris-asset ephemeris --input examples/stage12/field-line-requests.json --reviewed-by synthetic-reference-review --output build/stage12-example/field-line-requests.json
```

The job fits a known low-degree magnetogram and two density members, transforms
the synthetic observer ephemeris and publishes independent qualification and
withheld assets. The preregistration defines the complete five-factor product
with two wind/density members. Each realized row binds all construction
sources and the frozen density reconstruction recipe, then records independent
D6/topology/D1/D2 evidence and front radii. The radio-frequency likelihood is
recomputed for each density member with explicit fundamental/harmonic priors
and covariance. Inspect both rows and normalized weights in `selection.json`;
no SEP value participates in their calculation. The withheld product is only
eligible for evaluation after this selection is frozen.

Each successful command prints its frozen identity. Exit status 2 means an
invalid input, provenance/coverage failure or attempted overwrite; selection
status 1 means every candidate failed its preregistered gate, with all rows
still written for review. These neutral output records require a host adapter
before they can configure an AMPS run. Reviewed line requests contain covered
ephemeris seeds; they do not automatically select a mesh corridor.

For inference assumptions, transfer/MHD schemas and canonical test coverage,
see [the Stage-12 implementation guide](../../docs/STAGE12_OBSERVATION_PREPROCESSING.md).
