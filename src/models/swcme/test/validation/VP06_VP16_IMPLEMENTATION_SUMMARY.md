# VP06–VP16 implementation summary

## Delivered structure

Each validation ID has a dedicated `vpNN/` directory containing a detailed
README, `case.json`, provenance, an offline deterministic acquisition helper,
an independent `reference_solution.py`, focused oracle tests, a case runner,
and generated CSV/JSON/PNG/EPS evidence. `common_model_driver.cpp` is a thin
strict-warning probe of public SWCME shock, Parker, 1-D SEP, and 3-D SEP APIs;
it contains no acceptance thresholds or independent reference equations.

## Case outcomes

| Case | Evidence class | Focused outcome |
| --- | --- | --- |
| VP06 | Model/reference | 3/3 oblique shocks solved; worst normalized RH residual `4.01e-14`. |
| VP07 | Model/reference | Worst independent obliquity error `3.55e-14 deg`. |
| VP08 | Model/reference fixture | Hit/miss accuracy `1.0`; worst arrival error `1.5 h`. |
| VP09 | Model/reference fixture | Region classification accuracy `1.0`; boundary ordering valid. |
| VP10 | Analytical reference | Path and focusing comparisons agree at binary64 reporting precision. |
| VP11 | Literature constraint | Three observer histories; state accuracy `1.0`. |
| VP12 | Coupling benchmark | 37-field schema; source-spectrum error `1.39e-15`. |
| VP13 | Coupling benchmark | Five times; 1-D/3-D maximum relative difference `0`. |
| VP14 | Coupling benchmark | Active production source; controlled spread error `0`. |
| VP15 | Coupling benchmark | Three energy channels; profile error `1.20e-15`. |
| VP16 | Campaign benchmark | Six held-out fixtures; mean skill `0.579`; median error `5.75 h`. |

## Scope boundary

The result files deliberately distinguish model/reference, literature,
coupling, and campaign evidence. VP12–VP15 do not claim external AMPS transport
skill because the AMPS solver is not part of this repository. Likewise, the
analytical V3/V4 fixtures do not replace a future uncertainty-bearing in-situ
shock catalog. Those limitations are preserved in every case result JSON and
README so a component PASS cannot silently become an observational release
claim.

## Verification performed

- Every case-owned reference test passed in priority order.
- Every VP06–VP16 runner passed and emitted both PNG and EPS figures.
- The global runner's six orchestration tests passed.
- The complete VP01–VP16 campaign passed with existing verified external
  VP01–VP05 inputs.
- `make -j2 test` passed all 127 deterministic tests, including ASan/UBSan,
  thread reproducibility, performance guardrails, and coverage closure.
- Coverage closure reported 89.60% production lines and 57.83% branches.

