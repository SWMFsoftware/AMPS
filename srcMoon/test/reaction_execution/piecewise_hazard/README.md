# Piecewise-affine optical-exposure controls

This runner compiles the maintained production
[`PiecewiseHazard`](../../../species/PiecewiseHazard.h) service with a thin
text I/O probe. It requires a C++17 compiler and Python3; Matplotlib produces
optional PNG/EPS diagnostics. It never configures AMPS or modifies `build/`.
Read [CONTRACT.md](CONTRACT.md) for all numerical and endpoint decisions.

```sh
/usr/bin/python3 -B srcMoon/test/reaction_execution/piecewise_hazard/run.py \
  --output-dir test_output/srcMoon/piecewise-hazard-NEW
```

Replace `NEW`; existing nonempty evidence directories are refused. `--cxx`
selects the compiler. Exit codes are0 PASS,1 evaluated FAIL,2 infrastructure
ERROR and77 missing compiler/SKIPPED. Public Python `run(output, compiler)`
returns the same structured result used by the parent supporting coordinator.
No permanent U/I identifier is introduced.

Fixtures declare mathematical constant/rising/falling affine functions,
competing channels, discontinuities, exact endpoints, strictly future optical
depths, gaps, zero-rate endpoints, partitions and very different rate/time
scales. They are **not observations, fitted lunar rates, a physical event
provider or permission to interpolate measured data**. Decimal closed-form
integrals/quadratics are independent of the production implementation. A
binary-exact subset additionally verifies the event and its predecessor
against the mathematical cumulative exposure, with no acceptance epsilon.

Exact segment-end events use the approved LEFT-limit rates, including an
explicit zero total when the left law ends at zero. Downstream selection must
refuse an empty CDF rather than substitute the next segment's rate. The
entire requested interval must have supplied support; an early possible
event does not excuse a later missing segment. Invalid/incomplete/numerical
queries preserve supplied output objects and all history values.

The numerical module uses an analytical affine integral and an analytical
quadratic seed. Private standard-library-only positive integer/rational
arithmetic makes threshold comparisons exact for the supplied binary inputs;
root certification stops at adjacent representable long-double timestamps.
`maximumRootRefinements` is an explicit caller work limit, never a physical
step/tolerance. Returned exposure/rates use long-double arithmetic with
exponent-scaled positive products; out-of-range final values are refused.

Artifacts include input fixtures, compiler/source/executable identity,
`compile.log`, native stdout/stderr, named `result.json` outcomes,
`comparison.csv` and `affine_hazard_comparison.png/.eps`. Figures compare
production with a mathematical reference and are not observational
comparisons. Missing Matplotlib skips visualization only. Missing/corrupt
receipts and compiler/runtime faults are ERROR; evaluated criterion misses
are FAIL. No local PASS activates the controller, debits a clock, selects a
reaction, closes I21–I24 or changes formal M1/M2 predecessors.
