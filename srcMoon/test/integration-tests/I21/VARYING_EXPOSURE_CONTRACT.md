# I21 varying photon-exposure development cases

These cases extend the existing [I21 acceptance contract](reference/acceptance.json)
and [numerical addendum](../../M3_QUALIFICATION_ADDENDUM_20261118.md). They
exercise the quadrature and first-hit clock when the hazard genuinely varies.
The original constant-cohort FAIL stays in its original artifacts.

| Item | Contract |
| --- | --- |
| Production entry points | `TrajectoryExposure::BuildPhotonHistory`, `PiecewiseHazard::LocateEvent` with a compensated clock, and `PiecewiseHazard::Integrate`. The standalone runner links the production implementations directly. The linked scope additionally uses the actual accepted arcs, source-bound Na provider, and committed PIC event receipts. |
| Analytical inputs | Right-handed Cartesian SI axes, fixed point Sun at the origin, no occulting body, straight line starting 4 m from the source, moving outward at 1 m/s, 2 s total, reference rate 0.75 s⁻¹ at 4 m. These are dimensioned numerical verification inputs. SPICE epochs and lunar domain limits apply separately to native source-bound trajectories. |
| Numerical controls | Fixed 8, 16, 32 outer steps over the same exactly represented horizon, absolute bound 0.001 per quadrature call with zero relative contribution, and the same controls at each refinement. For a 0.25 s panel the independent curvature bound of 0.28125 s⁻³ bounds the trapezoid error by 0.000367 in optical depth. It permits a measurable step-dependent error. |
| Oracle | H(t) = k_ref*r_ref²*t/[r0*(r0 + u*t)]. The first hit for the deterministic optical-depth threshold 0.5 is 0.8 s, cumulative depth at the horizon is 1, and threshold 1.125 survives. Evaluate analytically with Decimal before comparing to reported long-double values. The standalone case is deterministic. |
| Procedure and metrics | Build exposure over each production arc, concatenate its intervals at their actual relative times, run the existing native event locator, and independently compare cumulative depth, first-hit time, and survivor classification. Report absolute depth error and maximum time error normalized by the common horizon. Check that the reported interpolation error covers the integral discrepancy. |
| Criteria | I21's 0.002 numerical ceiling for normalized first-hit errors and the frozen strictly decreasing or exact-zero rule. Exact classifications and positive accepted LEFT-limit rate. Production-reported error envelopes include explicit arithmetic representation allowance at the independently audited I/O boundary. |
| Data and identities | Standalone fixture consumes no external files. The linked source-bound case inherits exact externally pinned D04, state profile, SPICE, build and input identities from the I21 preparation. Both governing roadmap hashes are in the execution contract. |
| Artifacts and statuses | Preserve source/build/executable identities, stdout/stderr, case input, JSON and CSV/PNG/EPS numerical comparisons. Evaluated mismatch is FAIL. Tool/build/parser failure is ERROR. Missing prerequisite is SKIPPED. Scoped controls never qualify full I21. |

The prescribed source is fixed and always visible. Separate check cases cover
moving shadow boundaries and actual SPICE geometry. They remain required for
the full I21 matrix.

## Linked straight-line native scope

The approved Na P02 initiation supplies 128 original parents at fixed
`[10,000,000, 0, 0] m` in frozen LSO, with 65,536 s total, weight 2.625 and
seed 731993. Use quiet D04 and the actual SPICE Sun at the commissioning
epoch. The numerical variant declares velocity `[256, 0, 0] m/s`, Earth
occultation off, and an absolute optical-depth tolerance of `1e-10` per
production exposure call. This tolerance is fixed before launch. The
distance remains positive and its analytical axial curvature gives a
coarsest-step depth interpolation bound below `1e-10`. The same
production point-Sun qualification geometry supplies the lunar opacity.

The independent oracle requires a spherical Moon and can score the case
when every ray along the straight path lies strictly on one side of every
opaque sphere. It independently brackets the effect of the native Sun's
transverse offset and uses the original parent-keyed clock to check
classification and the exact native first-hit times and positions. A
moving-path survival probability is calculated from the integrated depth.
The constant-lifetime estimator is reserved for a homogeneous cohort.
