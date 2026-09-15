# OV04 — 2014 January 6 connectivity-sensitive diagnostic

OV04 compares three linked first-passage spectra with the PAMELA event fluence
in Bruno et al. (2018), Figure 4 (panel labelled 2014/01/06), DOI
`10.3847/1538-4357/aacc26`.  The immutable red marker centers span 90–1100 MeV;
the plotted blue fit is not used as reference data.

Early (0 min), nominal (30 min), and late (90 min) connection-delay runs share
the same source, solar wind, mean-free-path law, source spectrum, and random
sampling budget.  One global amplitude is used across all three spectra.  The
case is diagnostic-only because the behind-limb source and anisotropic GLE
cannot be uniquely mapped onto a one-dimensional Parker line; RMSE and
correlation remain reported but do not gate the suite.  Complete reference
coverage and successful linked execution do gate it.

```sh
python3 test/run_tests.py --amps ../amps --validation-case OV04 \
  --output-dir test_output/OV04
```

The output PNG/EPS contains one panel per connection realization.  Source
times, fit values printed in the publication, source hash, and missing inputs
are preserved in the JSON provenance files.
