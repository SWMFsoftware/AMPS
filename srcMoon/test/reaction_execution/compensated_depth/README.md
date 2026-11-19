# Full-component variable optical clocks

The [contract](CONTRACT.md) was written before this correction. The runner
links the SAME production `ReactionExposure`, `PiecewiseHazard` and
`ReactionExecution` sources, without configuring AMPS or changing physical
inputs. It extends supporting U15/U16 mechanics; it closes no M3 gate.

```sh
python3 -B srcMoon/test/reaction_execution/compensated_depth/run.py \
  --output-dir test_output/srcMoon/NEW-compensated-depth
```

The output must not exist. C++17 and matplotlib are required. Exits0/1/2/77
mean PASS/FAIL/ERROR/SKIPPED. Preserve result/manifest, commands, compiler and
source/executable hashes, raw stdout/stderr and CSV/PNG/EPS comparisons.
The global Moon supporting runner registers `compensated-depth`; this is
numerical software coverage, not a physical cohort or observation.

The exact pending value is the sum of two native binary64 components, not its
rounded long-double display. Arbitrary-size integer rational comparisons
locate events from ALL components, including1 minus double_denorm_min.
An opaque matching-origin locator receipt owns event publication and the
instantaneous CDF. No-event rational residuals are rounded high-then-low into
the existing two-double storage, with a rigorous explicit remaining-rounding
bound; positive underflow is refused. No physical acceptance epsilon or
fixed256-bit approximation replaces the decision. The old explicit scalar
projection API still refuses an unrepresentable projection.

Controls include10000ordinary7e-6s^-1 no-event flights, dyadic partitioning,
outer-interval continuation, near-boundary low components, full exponent-span
cancellation, non-dyadic residual rounding, cache identity, transactional
refusals, and the unchanged instantaneous competing-channel CDF. Python
Fraction parses C++ hex floats exactly, rather than first narrowing them to
binary64. It independently checks residuals and reported storage bounds.
Accumulated numerical error and linked MPI transport still require their
separate convergence and native tests; no qualification is inferred here.
