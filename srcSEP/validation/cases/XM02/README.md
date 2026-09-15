# XM02 — published M-FLAMPA Parker-spiral comparison

XM02 provides a controlled interface for comparing a production srcSEP run
with M-FLAMPA. The bundled reference is a digitization of Zhao et al. Figure 7:
M-FLAMPA >10 MeV proton intensity for far-upstream mean free paths 0.05, 0.3,
and 1 au during the 2013 April 11 event. It is useful for propagation/MFP
sensitivity, but it is not the exact quiet Parker-spiral CCMC export originally
requested by the validation plan. The case therefore cannot return PASS until
the reference and configuration-equivalence review are explicitly approved.

Prepare a production export with the exact header:

```text
elapsed_hours,series,intensity
```

`series` must use `mfp_0.05au_integral_gt10mev`,
`mfp_0.3au_integral_gt10mev`, and `mfp_1.0au_integral_gt10mev`; time is hours
since 2013-04-11 06:00 UTC and intensity is pfu. Copy `input.json`, set
`model_source_csv` to that file, and record whether injection, magnetic
connection, observer extraction, energy integration, mean-free-path law,
background evolution, normalization, and time origin are equivalent. Set
`equivalence_reviewed=true` only after that review. The linked application
validates and normalizes the table; Python never manufactures the model data.

```sh
python3 test/run_tests.py --amps /path/to/linked/srcSEP/AMPS \
  --validation-case XM02 --case-input /path/to/reviewed-xm02.json \
  --output-dir test_output/XM02
```

With no model export, the command stages the immutable reference and PNG/EPS
reference plots and returns SKIP. With an export but no equivalence approval it
computes metrics and still returns SKIP. After approval, all three series must
be present, log-intensity RMSE must be no more than 0.3 dex, and peak-time
error no more than 2 h. Publication/figure coordinates, PDF checksum, axis
calibration, and 0.12 dex extraction uncertainty are in `reference/provenance.json`.
