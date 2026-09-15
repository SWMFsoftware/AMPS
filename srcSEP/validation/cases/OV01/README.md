# OV01 — 2013 April 11 near-Earth observational benchmark

OV01 is the release-gating form of the 2013-04-11 Earth comparison.  The linked
`srcSEP/AMPS` executable advances the same event-informed Parker calculation as
XM03, while OV01 applies the tighter observational criterion from the validation
campaign.  It compares 4 h, 12 h, and 36 h spectra with 80 ACE/EPAM,
GOES-13/EPEAD, and SOHO/ERNE points extracted from Liu et al. (2025), Figure
12(a–c), DOI `10.3847/1538-4357/adc4e3`.

The reference and Figure 12(d) source history are intentionally stored once
under `../XM03/`.  `input.json` names that source case explicitly; neither the
runner nor a command-line override may replace it.  The model is normalized by
one amplitude across all times, energies, and instruments because the paper
does not publish the shock area needed for absolute one-dimensional injection.
No per-panel or per-instrument fitting is performed.

Run from the `srcSEP` directory:

```sh
python3 test/run_tests.py --amps ../amps --validation-case OV01 \
  --output-dir test_output/OV01
```

The result directory contains native arguments, linked registry JSON/JUnit,
model/reference/comparison CSV files, provenance hashes, and PNG/EPS overlays.
See `../XM03/README.md` for the vector extraction procedure and observational
cautions (background subtraction, SOHO saturation, and ACE ion contamination).
