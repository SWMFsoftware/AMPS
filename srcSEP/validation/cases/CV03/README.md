# CV03 — nonuniform diffusion and stochastic-calculus drift

## Purpose and isolated physics

CV03 verifies that spatially varying diffusion is interpreted in conservative
Ito form. It prescribes `kappa(s)=kappa0[1+0.4 sin(2 pi s/L)]` on a periodic
line. The corresponding stochastic characteristic must include `d kappa/ds`;
omitting or reversing this drift produces the wrong flux and equilibrium even
when the random diffusion amplitude looks plausible.

The linked Parker mover runs a localized transient through three timesteps
using both the analytic coefficient derivative and a centered numerical
derivative. Common keyed random numbers isolate derivative implementation
differences. It also advances a stratified uniform population, which is the
exact zero-flux equilibrium, and saves a sign-reversed-drift negative control.

## Independent reference and gate

`reference_solution.py` solves `partial_t f=partial_s(kappa partial_s f)` with
a separate 256-cell periodic finite-volume method. Harmonic face coefficients,
conservative flux differences, and a documented explicit stability factor are
used; no C++ transport routine, derivative, or random stream is imported.

The scorer measures transient relative L2 error, analytic-versus-numerical
derivative L2 difference, uniform-equilibrium bin error, probability closure,
the reconstructed equilibrium probability flux, and automatic endpoint
refinement order. The negative control must be worse
than the nominal finest solution. The native ensemble uses four independent
seeds, 6,000 particles, 64 reported bins, and 4/2/1 s timesteps. All SI values
and thresholds are reviewable in `input.json`.

## Run and evidence

```sh
python3 test/run_tests.py --amps /absolute/path/to/amps \
  --validation-case CV03 --output-dir test_output/CV03
```

Evidence includes the exact linked command, native JSON/JUnit,
`CV03_model.csv`, independent finite-volume `CV03_reference.csv`, a normalized
`CV03_solution.csv`, SHA-256 provenance, and PNG/EPS comparison/residual plots.
The linked binary is mandatory and must advertise CV03. Generic native runs
without `--test-input` and `--test-output-dir` intentionally report SKIP.

## Interpretation and limitations

Equilibrium bias that follows `kappa(s)` indicates a missing/incorrect Ito
drift. Agreement of the analytic derivative but not the numerical derivative
indicates differencing or periodic-wrap errors. Probability loss indicates a
boundary or histogram defect. Because stochastic profile error contains Monte
Carlo noise, the refinement result is a weak-order campaign diagnostic; it is
not a deterministic PDE-order proof. CV03 does not validate a heliospheric
background provider or observational data.
