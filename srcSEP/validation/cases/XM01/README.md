# XM01 — independent focused-transport PDE solver comparison

XM01 checks the linked production `fte-dmumu` characteristic against a solver
that shares neither srcSEP code nor random numbers. Five periodic uniform-line
problems isolate streaming, pitch-angle scattering, magnetic focusing,
adiabatic momentum change, and all operators together. Internal quantities are
SI: arclength in m, time in s, speed in m s^-1, momentum in kg m s^-1,
`D_mumu` in s^-1, `d ln|B|/ds` in m^-1, and `div U` in s^-1.

The production side advances samples with
`AdvanceFocusedTransportDmumu`, `D_mumu=D0(1-mu^2)`, periodic arclength, and
the reviewed initial density
`f(s,mu)=Gaussian(s;0.35L,0.06L)*(1+0.6mu)/2`. It runs `N` and `4N` samples
with deterministic keyed streams. The independent Python solver translates
the smooth periodic packet with a Fourier phase and advances zero-normal-flux
pitch transport with conservative finite-volume fluxes and Strang splitting.
Adiabatic momentum follows the independently evaluated reduced-equation
characteristic `p(t)=p0 exp[-(div U)t/3]`.

Acceptance uses only reference cells above 2% of the reference peak. The fine
resolved log-intensity RMSE must be at most 0.05 dex, the first angular moment
difference at most 0.02, the natural-log momentum error at most 0.002, and the
fine/coarse RMSE ratio at most 0.9. The threshold excludes statistically empty
tails; it does not permit dropping a populated central bin.

Run:

```sh
python3 test/run_tests.py --amps /path/to/linked/srcSEP/AMPS \
  --validation-case XM01 --output-dir test_output/XM01
```

The case writes the resolved input, native argument manifest, linked model
CSV, independent reference CSV, native and aggregate JSON/JUnit, checksums,
log, and PNG/EPS comparisons. It is an extended cross-model test: its PASS is
not observational or coupled-SWMF validation.
