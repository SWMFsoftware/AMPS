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

Run `python3 surface_contributions.py -h` for every option. Products are described there and written to `--out`
(never overwritten). The job script used on Voyager is `$COPS_OUTPUT_ROOT/surface_contributions/surface.sbatch`.
It is checked by `check_against_base.py`, which recounts the per-sample classes from a fitting base's
`kernels/point=NNNN.bin`. For 12M (4 late samples) the counts were identical.
