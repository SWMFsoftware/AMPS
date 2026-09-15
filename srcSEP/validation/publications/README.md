# Publications used by XM01–XM03

The source tree stores small derived references and complete provenance, not
downloaded article PDFs. The delivery includes a separate publication packet
so repository users can inspect the exact source bytes without turning large
third-party documents into normal source-controlled inputs.

## Publication-derived model inputs

XM02 and XM03 each provide `publication_input.json`. These structured manifests
contain the model parameters that can be recovered from the papers, source
locators, figure-derived assumptions, and a concrete list of missing artifacts.
They deliberately use `reproduction_status: partial`: neither paper publishes
a complete executable SWMF run directory. The validation runner checks the
manifest schema and case identity, copies it into each evidence directory, and
hashes it in completed-comparison provenance.

`partial` remains an important provenance label even though XM02/XM03 now use
only these registered inputs. It states that the literature does not provide a
complete SWMF reproduction package; it does not create another selectable
configuration or require a command-line input override.

## Scored curve sources

1. Zhao et al., *Solar Wind with Field Lines and Energetic Particles (SOFIE)
   Model: Application to Historical Solar Energetic Particle Events*, arXiv
   2309.16903. XM02 digitizes Figure 7 on PDF page 30. Expected PDF SHA-256:
   `945361cc4c27e6481eac042eaf9a0e3f6b097ed277e1b7fab793ab9305591bbe`.
2. Liu et al., *Physics-based Simulation of the 2013 April 11 Solar Energetic
   Particle Event*, Astrophysical Journal 985:82 (2025), DOI
   `10.3847/1538-4357/adc4e3`. XM03 uses Figures 14–15 on PDF pages 23–24.
   Expected PDF SHA-256:
   `a5a613e9ad3b127c5c412366b0c4a2029339f9ac068fd9508325ab682a6dc357`.
   The article is distributed under CC BY 4.0; retain author/title/journal/DOI
   attribution with any redistributed copy.

## Method source

Borovikov et al., *Toward Quantitative Model for Simulation and Forecast of
Solar Energetic Particles Production during Gradual Events—II: Kinetic
Description of SEP*, arXiv 1911.10165, documents the M-FLAMPA kinetic method.
It is not a scored numerical reference. Expected PDF SHA-256:
`672194f51579c323061f17630f34ce21ae002a47a4a269c4668b246388825e4c`.

## Reproduction

Run the digitizer into a disposable review directory and diff its outputs with
the case references:

```sh
python3 validation/reference/digitize_mflampa_figures.py \
  --zhao-pdf /publications/zhao_2024_sofie_mflampa.pdf \
  --liu-pdf /publications/liu_2025_2013_april_11.pdf \
  --output-dir /tmp/xm-digitization-review --render-pages
```

The script checks exact PDF bytes before rendering. It never overwrites a
case reference, and reference regeneration is not part of test execution.
