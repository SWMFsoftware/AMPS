# OV02 — 2020 May 29 radial PSP/STEREO-A benchmark

OV02 runs the linked application at the 0.33 AU PSP and 0.96 AU STEREO-A
locations reported by Cheng et al. (2023), DOI `10.3847/1538-4357/acac21`.
The immutable comparison points are measurements digitized from Figure 3
(PSP/EPI-Hi, 2.2 and 12.3 MeV) and Figure 6 (STEREO-A/LET, 1.8–3.6 and
4.0–6.0 MeV). Published simulation curves are explicitly excluded.

Both spacecraft use one source spectrum, one 0.0465 AU 1-GV radial mean free
path, and one 337 km/s event speed.  A single logarithmic amplitude is fitted
to every point together—never one normalization per spacecraft or channel—so
relative radial, energy, and temporal behavior remains a prediction.  The
input pulse covers the paper's 07:38–08:10 UT low-coronal shock interval.  The
reconstruction is one-dimensional and does not claim to reproduce the paper's
PFSS field, ellipsoid shock, or perpendicular diffusion.

```sh
python3 test/run_tests.py --amps ../amps --validation-case OV02 \
  --output-dir test_output/OV02
```

The PNG/EPS title names the publication and Figures 3 and 6.  Exact source
hash, included/excluded traces, time origin, and digitization uncertainty are
in `reference/provenance.json`; all numerical assumptions are in `input.json`
and `publication_input.json`.
