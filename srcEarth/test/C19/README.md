Reference generation:
The builder:

    downloads or reads the GOES-13/15 monthly five-minute directional proton files;

    extracts P4 and P5;

    applies the event-specific telemetry-head-to-physical-direction mapping;

    calculates an independent pre-event background for each detector head;

    removes invalid, fill, flagged, and low-signal measurements;

    calculates the background-subtracted physical East/West ratio;

    writes a compressed reference table and complete provenance JSON;

    optionally uses NOAA one-minute ephemeris, with nominal GEO locations as a documented fallback.

Run:

python3 srcEarth/test/C19/build_goes_reference.py --download

This creates:

srcEarth/test/C19/data/reference_C19_goes_epead_ew.csv.gz
srcEarth/test/C19/data/reference_C19_goes_epead_ew_provenance.json

C19 runner

Example GRIDDED run:

python3 srcEarth/test/C19/run_C19.py \
  --profile ROUTINE \
  --solver GRIDDED \
  --models T96,T05 \
  --reference srcEarth/test/C19/data/reference_C19_goes_epead_ew.csv.gz \
  --driver srcEarth/test/C19/data/ts05_driver_may2012.txt \
  --amps ./amps \
  -np 4 \
  -nt 16

With parallel background-field initialization:

python3 srcEarth/test/C19/run_C19.py \
  --profile ROUTINE \
  --solver GRIDDED \
  --models T96,T05 \
  --driver srcEarth/test/C19/data/ts05_driver_may2012.txt \
  --mode3d-parallel-field-init \
  --amps ./amps \
  -np 4 \
  -nt 16

The runner:


    generates one GEO trajectory location per spacecraft and epoch;

    runs T96 and/or T05;

    reads the global SM directional cutoff map;

    constructs physical eastward and westward detector apertures;

    folds cutoff energy through the P4/P5 channel range;

    assumes a configurable common incident spectrum, J(E)∝E−γJ(E)∝E−γ;

    compares modeled and observed log10(East/West);

    calculates validity, sign agreement, bias, MAE, RMSE, and correlation.

Generated plots

A completed run produces:

C19_comparison_gridded_t96.png
C19_comparison_gridded_t05.png
C19_scatter_gridded_t96.png
C19_scatter_gridded_t05.png
C19_transmission_gridded_t96.png
C19_transmission_gridded_t05.png
C19_aperture_diagnostic.png

It also writes:

C19_reference_used.csv
C19_model.csv
C19_comparison.csv
C19_metrics.csv
C19_aperture_samples.csv
C19_commands.json
C19_result.json
C19_summary.txt

The built-in synthetic runner test successfully exercises directional-map parsing, detector-aperture folding, comparison output, and plot generation:

python3 srcEarth/test/C19/run_C19.py --self-test

Remaining correction

The GRIDDED path is ready for initial use. Before using --solver GRIDLESS or --solver BOTH, make this small correction in run_C19.py: render_template() currently treats solver-specific replacement keys absent from the other template as an error.

Replace its replacement loop and validation with:

for key, value in replacements.items():
    text = text.replace(key, value)

leftovers = sorted(set(re.findall(r"__[A-Z0-9_]+__", text)))
if leftovers:
    raise ValueError(
        "template contains unresolved placeholders: %s" % leftovers
    )

A complete AMPS execution was not possible here because the environment does not contain the AMPS executable or the May 2012 TS05 driver.
