#!/usr/bin/env python3
"""Interactive comparison of the COPS reference data with LiouvilleCOPS model runs as a function of altitude.

Same window as plot_cops_model.py (NG and RG together, NG, RG, NG/RG; check boxes to switch the runs on and off),
with the altitude above the surface (log scale, increasing to the right) as the horizontal axis.
The altitude of the reference points is that of the run selected with --ref (it depends on the shape model).


DATA SOURCES
  Reference data   columns refNG, refRG of COPS_model.dat [Pa]: the COPS NG and RG pressures of the srcCG table
                   srcCG/RosinaMeasurements.h (1150 points, 2016-09-30T00:00:04 .. 10:39:26, given in mbar and
                   converted by LiouvilleCOPS with ReferenceDataUnit = mbar). The table is identical to the ROSINA COPS
                   L2 BG housekeeping pressures of the PSA/PDS3 data set RO-C-ROSINA-2-EXT3-V2.0 (N2-equivalent
                   pressures, 20 C; time = end of each housekeeping interval). 0 in the table = no valid value
                   (e.g. NG 10:24:42-10:33:14). useNG/useRG = 0 marks the points not used in the fit (no data,
                   ExcludeNG/ExcludeRG/FitInterval of the run, gauge inside the shape model).
  Model            columns pNG, pRG [Pa] (and the species columns) of each LiouvilleCOPS run, or of a fit written by
                   tools/fit_cops.py (species columns -1).
  Geometry         UTC and altitude columns: SPICE geometry of the run (altitude = distance to the nearest point of
                   the run's shape model: SHAP5, M001 or M004).

usage:
  plot_cops_model_altitude.py [run directories or COPS_model.dat files ...] [--ref RUN] [--save FILE] [--theme light|dark]
"""
import os
import sys

if "--save" in sys.argv:
    import matplotlib
    matplotlib.use("Agg")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plot_cops_model import main

if __name__ == "__main__":
    main(default_x="altitude")
