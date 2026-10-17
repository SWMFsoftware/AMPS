# Stage-5 DDD and neutron-equivalent validation assets

`stage5_validation_plan.v1.json` fixes the automated acceptance cases,
histories, exact reference nodes, tolerances, strict-domain case, and secondary-
neutron sensitivity before the automated campaign. One manual 1.05 MeV neutron
smoke run preceded registration and is disclosed by `registration_status`; it
is not represented as blind validation.

The two CTest gates are:

- `stage5_core`: authenticates every response table through the production
  loader and checks exact SR-NIEL nodes, log-log interpolation, a closed-form
  charged-particle interval integral, an exact ASTM group fold, the 95 MeV-mb
  conversion, unsupported species, strict domain endpoints, and case-schema
  capability boundaries.
- `stage5_transport`: drives the normal case parser, primary generator,
  multilayer Geant4 geometry, provider-defined absorber entry surface,
  event-level scorer, accuracy controller, and result/response-cell writers.
  It compares exact-node proton, helium-4, neutron DDD, and neutron `n_eq` with
  independent calculations; verifies the `1/(A*mu)` scalar-fluence weight with
  an oblique proton case; requires a 25 MeV neutron to yield
  `not_assessable`; and demonstrates sensitivity to secondary neutrons behind a
  10 cm aluminum layer.

Generated inputs, logs, results, and response cells are written below
`<build>/stage5_transport_runs`. Compact numerical evidence is written to
`<build>/validation_evidence/stage5_core.json` and
`stage5_transport.json`. Review figures are generated only after those numeric
tests pass:

```text
<build>/validation_figures/stage5_core_analytic.{png,svg}
<build>/validation_figures/stage5_transport_comparisons.{png,svg}
```

The JSON assertions are authoritative; plots are review aids. These tests
qualify implementation/reference behavior only. They do not close the formal
Stage-5 requirement for project-selected accelerator/device evidence and an
independent review record.
