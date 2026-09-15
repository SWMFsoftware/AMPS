# OV03 — 2013 May 22 interacting-CME diagnostic

OV03 compares two linked transport hypotheses with GOES-15 and STEREO-A/HET
measurements digitized from Ding et al. (2014), Figure 1(b,d), DOI
`10.1088/2041-8205/793/2/L35`.  The single-source hypothesis begins with the
13:25 UT fast CME; the twin-source hypothesis also includes the 08:48 UT slow
CME seed episode.  Both use the reported 1439 km/s CME2 speed and 444 km/s
eight-hour solar-wind average.

The 141-degree Earth/STEREO-A separation and CME interaction make this a
three-dimensional stress case.  The one-field-line reconstruction cannot
represent the required cross-field transport, so profile RMSE and correlation
are explicitly diagnostic (`gating=false` in result metrics).  Reference
coverage and successful linked execution still gate the test: a missing curve,
empty native artifact, or broken registry callback cannot pass.

```sh
python3 test/run_tests.py --amps ../amps --validation-case OV03 \
  --output-dir test_output/OV03
```

All four panels (Earth/STA × single/twin) use the identical observed baseline
for the relevant spacecraft.  The per-series unit-peak normalization preserves
onset, rise, peak time, and decay while avoiding a false absolute comparison
between integral detector channels and a reduced differential source.
