# CV09 — planar diffusive shock acceleration

CV09 is a controlled shock-cycle benchmark executed inside the linked
application. It samples repeated crossings with constant fractional gain and
downstream escape for compression ratios 2, 3, and 4. Test-particle planar DSA
predicts the isotropic phase-space index `q=3r/(r-1)` and
`t_acc=3/(u1-u2)(kappa1/u1+kappa2/u2)`.

The scorer fits the empirical survival spectrum over a registered central
quantile interval, converts it to `q`, independently checks acceleration time,
number accounting, and resolved upstream/downstream diffusion lengths. The
50,000-particle keyed ensemble makes the slope uncertainty small without
embedding the expected answer in the model. CSV output also labels `p*v` as an
energy proxy; it is not claimed to be relativistic kinetic energy.

Scope limitation: this case verifies the analytical shock-cycle source logic,
not a resolved full PIC shock geometry. A production field-line shock run with
finite boundaries remains a separate native validation gate.

Run `python3 test/run_tests.py --amps /path/to/amps --validation-case CV09
--output-dir test_output/CV09`.
