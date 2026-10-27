# Float-vs-double backtracing experiment

Standalone C++17 program (no AMPS dependency) that measures what single precision does
to backward trajectories as traced by the legacy tracer in a centered dipole.  It is the evidence
behind the precision/units policy of the batched engine; the analysis is in
[`../../doc/PRECISION.md`](../../doc/PRECISION.md).

## Build and run

```
make                     # g++ -O3 -march=native -std=c++17 -ffp-contract=off
make -j 3 run            # BORIS, HC4 and RK4 in parallel -> $OUT/out_<MOVER>.txt
./precision_experiment HC4 > out.txt      # one mover
```

`-ffp-contract=off` stops the compiler from fusing `a*b+c` into FMA, so the variants
differ only by their arithmetic type and the result does not depend on whether
`-march=native` provides FMA.  CPU only; single threaded.

Results are written outside srcGeospace: `OUT` defaults to
`$GEOSPACE_RUNS/results/precision_experiment/run` (`GEOSPACE_RUNS` defaults to
`<repository>/../AMPS_runs/srcGeospace`); `make run OUT=<directory>` selects another
directory.

Run time on this machine (2026-10-21, 3 movers concurrently, load average ~13 on 36
cores): BORIS 12.6 s, HC4 148.8 s, RK4 25.9 s (`time_*.txt` of the 2026-10-21 rerun below).

## What it computes

Units: length Re (6371.0 km), speed c, time Re/c = 0.021251 s, momentum m_p c, charge
e, field B_u = m_p c/(e Re) = 4.9125e-7 T (the engine itself uses amu instead of m_p,
B_u = 4.8770e-7 T; the 0.7 % difference is irrelevant here).  Dipole moment along +z
with B_eq = 3.12e-5 T (as `gridless/DipoleInterface.h`).  Proton, positive unit charge,
launched opposite to the arrival direction.  Stopping surfaces: box |x|,|y|,|z| <= 35 Re,
inner sphere r = 1 Re (segment test), 600 s, 400 Re path, 300000 steps.  Time step
min(DT_TRACE = 1 s, 0.15/Omega_gyro, time remaining) recomputed every step (the
structured legacy policy).

Movers (legacy arithmetic): `BORIS` (B at x_n, half drift, rotation, half drift),
`HC4` (Yoshida 4th-order composition of Higuera-Cary magnetic-midpoint steps, applied
as 4 substeps per step like `StepParticleChecked`), `RK4` (classical).

Variants, as template parameters `<state, accumulators, position>`:

| tag | x | p, B, mover arithmetic | time/path accumulators |
|---|---|---|---|
| DD  | double | double | double |
| DDp | double | double | double, every dt multiplied by (1+1e-7): sensitivity calibration |
| FD  | float  | float  | double |
| XD  | double | float  | double |
| FF  | float  | float  | float |

Sections of the output:

* **A** vertical cutoff by UPPER_SCAN (100-point log rigidity grid 1 MeV..20 GeV,
  downward scan to the first forbidden point, bisection to 1e-3 GV) at
  altitudes 500/9000 km and latitudes 0/30/45/60 deg, with the Stormer value.
* **E** adaptive per-trajectory dt vs. one global fixed dt (the DD workload of B).
* **B** allowed/forbidden classification of a 24 (cos theta) x 48 (phi) direction grid at
  45 deg, 9000 km, for 12 rigidities R = Rc_Stormer x (0.40..2.05); transmission T(R)
  and the number of directions whose classification differs from DD.
* **C** a trapped 10 MeV proton at 4 Re integrated for 1200 s: termination, steps,
  relative |p| spread.
* **D** the time accumulator alone: 2.4e7 additions of dt = 5e-5 s in float vs double.

## Files

* `precision_experiment.cpp` -- source (copied unchanged from the design-phase
  scratch directory).
Outputs (kept outside srcGeospace, in `$GEOSPACE_RUNS/results/precision_experiment/`,
i.e. `/data/vtenishe/GEO/AMPS_runs/srcGeospace/results/precision_experiment/` here):

* `original/out_*.txt` -- the design-phase outputs.  They were produced by an earlier
  build of the same source that did not yet print section E.
* `2026-10-21_rerun/out_*.txt`, `time_*.txt`, `make_run.log` -- rerun of 2026-10-21 with
  the Makefile settings.  Sections A-D are identical, digit for digit, to `original`;
  section E is new.
