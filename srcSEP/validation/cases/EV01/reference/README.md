# EV observational-reference provenance

The CSV in this directory contains **observationally derived NASA CCMC SEP Model
Validation Challenge quantities**, not a theoretical or synthetic reference.

Source hierarchy:

1. NASA CCMC Model Input Parameters page — official numeric tables for
   `>10 MeV / 10 pfu` and `>100 MeV / 1 pfu` GOES proton measurements.
2. NASA CCMC per-event pages — unnumbered GOES Proton Measurements figures and
   adjacent threshold tables used to verify individual events.
3. NASA CCMC Data Sets page — documents that the event quantities are computed
   event-by-event with OpSEP / `operational_sep_quantities.py`.

The current figures on the CCMC event pages are **unnumbered**.  Do not cite a
fabricated figure number.  Cite the event page and the “GOES Proton
Measurements” section/threshold table.  The exact values committed here were
transcribed from the tables; the figures are supporting visualizations only.

See `provenance.json` for exact URLs, table names, processing information, and
the explicit statement that no raw NOAA/NCEI file identifier is claimed by this
pilot bundle.
