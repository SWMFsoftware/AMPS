# OV05 — September 2017 compound-event diagnostic

OV05 runs one linked Parker transport calculation with separately tagged
September 4, 6, and 10 source pulses.  Their times come from Bruno et al.
(2019), DOI `10.1029/2018SW002085`; the comparison uses STEREO-A LET/HET
profiles digitized from Figure 2 at 4.25, 11, 31.5, and 80 MeV.  The plot title
retains both the citation and figure number.

The event is an intentionally difficult wide-longitude test: STEREO-A was
about 128 degrees east of Earth, the September 10 arrival was delayed by more
than ten hours, and the paper identifies ICMEs, high-speed streams, and a later
connection change.  A one-dimensional field line cannot reproduce those
structures.  OV05 therefore records shape RMSE/correlation as diagnostic
metrics while gating only linked execution and reference coverage.  One global
amplitude is used across every energy and time; no channel is independently
rescaled.

```sh
python3 test/run_tests.py --amps ../amps --validation-case OV05 \
  --output-dir test_output/OV05
```

The fixed input is `input/three_injection_source.csv`.  Exact publication
hash, event epoch, extraction method, uncertainty interpretation, and missing
physics are in `publication_input.json` and `reference/provenance.json`.
