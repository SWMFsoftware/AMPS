# Utility/surface: surface elements seen by the COPS NG only or the RG only

`surface_contributions.py` writes a LiouvilleCOPS input (gauge and transport settings from a template, default
`../../LiouvilleCOPS/fit_base.in`; paths from the command line). It compiles `surface_contributions.cpp` with the
unchanged LiouvilleCOPS sources, except `main.cpp`, against CSPICE. It then runs the program for the observation
times in a UTC interval and/or an altitude interval, and prints and writes the counts and areas of the elements in
each class. Run it on a compute node: it compiles with `-march=native`.

**Definition (user decision 2026-10-27).** The contribution of an element is the gauge pressure response to a unit
uniform H2O flux at `--temperature` (Lambertian emission, the forward model's kernels and gauge conversions). With
`eps = --epsilon` (default 1e-3), as in the fine tuning (`LiouvilleCOPS/FINE_TUNING_STEP.md` §3):

* **NG only** at a sample: H_NG > 0 and |H_RG| <= eps·H_NG.
* **RG only** at a sample: |H_RG| > 0 and H_NG <= eps·|H_RG|.
* **At least once:** at one or more samples.
* **All the time, strict:** at every sample of the interval.
* **All the time, loose:** at every sample at which either gauge sees the element.

**Running counts (`samples.dat`, 2026-10-28).** After each sample k, the last 7 columns of its row count over
the samples 1..k, so one run gives the "all the time" classes for every interval end:

* `nNGonlySoFar`, `nRGonlySoFar`: elements NG only (RG only) at every sample so far, i.e. "all the time, strict"
  for an interval ending at sample k; `areaNGonlySoFar`: their area [m²].
* For the elements seen by NG at every sample so far, the largest ratio |H_RG|/H_NG so far:
  `nNGseenSoFar_ratio_le_eps`, `_le_10eps`, `_le_100eps` count those with that ratio <= eps, 10·eps, 100·eps, and
  `minMaxRatioSoFar` is the smallest of these largest ratios: how close the best element comes to "NG only all
  the time" (+inf once NG has lost every element).
* An element leaves for good at the first sample NG does not see it, so the counts never rise and
  `minMaxRatioSoFar` never falls as k grows.

Run `python3 surface_contributions.py -h` for every option. Products are described there and written to `--out`
(never overwritten). The job script used on Voyager is `$COPS_OUTPUT_ROOT/surface_contributions/surface.sbatch`.
It is checked by `check_against_base.py`, which recounts the per-sample classes from a fitting base's
`kernels/point=NNNN.bin`. For 12M (4 late samples) the counts were identical.
