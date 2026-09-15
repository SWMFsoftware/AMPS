# XM03 — published M-FLAMPA 2013 April 11 event reproduction

XM03 compares linked srcSEP products with Liu et al. (2025), DOI
`10.3847/1538-4357/adc4e3`. The reference contains Earth 2.0–2.5 and 20–25 MeV
time-intensity curves for lambda0=0.3 au from Figure 15(a,b), the Earth-connected
radial mean-free-path profile from Figure 14(b), and the reported Earth fluence
spectral index -1.78 from Figure 15(c). The event interval begins at 2013-04-11
06:00 UTC. Differential intensity units are `(cm2 s sr MeV)^-1`, radius is in
solar radii, mean free path is in au, and the fitted index is dimensionless.

Export the production comparison products in long form:

```text
observable,coordinate,value
```

The required observable names are exactly those in
`reference/mflampa_2013_apr11_event.csv`. Time-series coordinates are elapsed
hours, mean-free-path coordinates are solar radii, and the spectral-index
coordinate is zero. The linked callback validates finite coordinates/values,
permits a negative value only for the signed spectral index, and writes its
normalized artifact transactionally.

Copy `input.json`, set `model_source_csv`, document the magnetic connection,
shock history, injection law/normalization, transport coefficients, energy
response, observer sampling, background evolution, and time origin, then set
`equivalence_reviewed=true` only after review. Until then the result is SKIP,
even when diagnostic metrics and figures exist.

```sh
python3 test/run_tests.py --amps /path/to/linked/srcSEP/AMPS \
  --validation-case XM03 --case-input /path/to/reviewed-xm03.json \
  --output-dir test_output/XM03
```

Provisional gates are a factor-three (0.477 dex) intensity RMSE, 3 h peak-time
error, 0.3 spectral-index error, 0.2 dex mean-free-path RMSE, and complete
series coverage. The plan's 1 h onset metric requires an agreed background and
threshold definition and is recorded in input but is not scored prematurely.
Model-to-M-FLAMPA and model-to-spacecraft results must be reported separately.
Digitization calibration and uncertainty are in `reference/provenance.json`.
