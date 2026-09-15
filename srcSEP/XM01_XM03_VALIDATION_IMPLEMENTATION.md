# XM01–XM03 cross-model validation implementation

## Scope and evidence boundary

XM01 is a genuine linked cross-solver benchmark: the selected AMPS executable
runs production focused-transport characteristics and an independent Python
PDE solver supplies the reference. XM02 and XM03 are publication-comparison
workflows. They validate and normalize a caller-supplied production export
inside the linked executable, but default to SKIP because neither a production
event export nor a signed configuration-equivalence decision is bundled.

## Publications and reference extraction

- Zhao et al., *Solar Wind with Field Lines and Energetic Particles (SOFIE)
  Model: Application to Historical Solar Energetic Particle Events*, arXiv
  2309.16903, Figure 7, supplies the XM02 MFP-sensitivity profiles.
- Liu et al., *Physics-based Simulation of the 2013 April 11 Solar Energetic
  Particle Event*, ApJ 985:82 (2025), DOI
  10.3847/1538-4357/adc4e3, Figures 14–15, supplies XM03.
- Borovikov et al., arXiv 1911.10165, documents the kinetic M-FLAMPA equation
  and numerical context; it is supporting methodology rather than a scored
  curve source.

`validation/reference/digitize_mflampa_figures.py` checks the downloaded PDF
hashes, optionally renders the exact PDF pages at 400 dpi, and applies reviewed
pixel-to-data transforms. Case-local provenance records preserve publication,
page, panel, plot box, axes, colors, sampling interval, units, uncertainty, and
limitations. Regeneration writes elsewhere by design so baseline changes can
be diffed and reviewed.

## Execution and outputs

Use `make test-xm01-xm03-unit` for the strict source/reference gate. Set
`SEP_EXECUTABLE` to add linked execution. For normal use:

```sh
python3 test/run_tests.py --amps /path/to/amps \
  --validation-case XM01 --validation-case XM02 --validation-case XM03 \
  --output-dir /evidence/XM01-XM03
```

XM01 should PASS when its statistical and moment gates pass. XM02/XM03 stage
reference-only PNG/EPS plots and SKIP until configured. A reviewed custom input
enables linked model ingestion, metrics, overlays, and—only after explicit
equivalence approval—PASS/FAIL acceptance. JSON remains authoritative.

## Remaining limitations

The exact CCMC M-FLAMPA v1 canonical-run output was not available as a numeric
table, so XM02 uses a clearly labelled published event sensitivity proxy. The
published curves are digitized rather than author-provided numeric products.
The code does not claim observational validation, and the full SWMF background,
CME/shock reproduction, multi-observer comparison, or instrument forward model
remain separate gates.
