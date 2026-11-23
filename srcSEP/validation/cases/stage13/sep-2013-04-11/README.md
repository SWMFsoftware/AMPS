# 2013 April 11: preserve the existing OV01 observational comparison

This case is a continuity entry for **OV01**, not a new qualification of a
spherical or nonspherical moving shock/background. `case.json` pins the actual
existing input, publication manifest, observations, source proxy and provenance
by relative path, size and SHA256. It delegates to the maintained OV01 runner
and retains its numerical settings, normalization and release criteria.
No existing file is copied, regenerated, replaced or retuned.

The existing reference has 80 near-Earth spectral points from ACE/EPAM,
GOES-13/EPEAD and SOHO/ERNE at three times. They were extracted from Figure
12(a–c) of [Liu et al. (2025)](https://doi.org/10.3847/1538-4357/adc4e3),
with [published-paper provenance and license](https://authors.library.caltech.edu/records/qwsnd-jmx38).
The shared XM03 provenance describes the vector extraction and instrument
limitations. Original processed instrument arrays are not claimed.

The observational elapsed-hour clock starts at07:24UTC, while the shared source
history starts at06:00UTC. These origins must remain distinct. The plotted
Earth-connected thermal-energy history is a model-informed source proxy, not
independent observational evidence. The original one global amplitude across
all times, channels and instruments is retained; it does not become a license
to fit each instrument or adopt that rule for the other new cases.

Run the existing gate from `srcSEP`:

```sh
python3 -B test/run_tests.py --amps /absolute/amps --validation-case OV01 \
  --output-dir /absolute/fresh/OV01
```

The Stage13 I3-style central selector may delegate this same command. Its
result must be carried through honestly; package integrity or an alias must
not turn a failed/incomplete observational run into PASS. `comparison.status`
is `NOT_RUN` until such a run exists. New geometry qualification additionally
requires its actual event-specific shock, material field-line/background and
forward-observation contracts; the original OV01 digitization cannot supply
missing three-dimensional forcing or uncertainty.
