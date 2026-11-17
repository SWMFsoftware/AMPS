# I20 polar comparison implementation contract — 2026-11-16

This extends the existing I20 meaning; no new permanent test ID or acceptance
threshold is introduced. The initial implementation inventory is: production
LOLA/BVH and thermal columns exist; passive cohorts record real surface/gas,
first release, censored residence and actual hazards; the numerical cluster
scorer exists, but I20 always remains SKIPPED because there is no independent
collocated polar comparison. The new module adds that missing evidence path.

| Required field | Contract |
|---|---|
| ID/name | I20, independent collocated polar environment and conditional retention comparison. |
| Purpose | Evaluate the original PSR/temperature/kinetic gates only after numerical, source and observation-operator predecessors qualify. |
| Production capability | Normal LOLA illumination, surface thermal state and Na/Ar finite-inventory cohorts. |
| Production entries | Existing I15/I17 production environment exports and `Moon::M2Cohorts`; no replacement force/thermal/source law. |
| Configuration | Existing approved Na/Ar numerical design and an explicitly reviewed/frozen external comparison contract. |
| Inputs | Hashed native-site candidate/selected registries, collocation and quantity/time-support contracts, independently processed D02/D03 rows, production model rows/manifests, actual cohort/hazard files. |
| Independent reference | Authoritative D02/D03 for environment; existing independent hazard integration for conditional first-release/residence identities. A conditional kinetic PASS does not validate the thermal history. |
| Units | Dimensionless occurrence/classification, K with explicit depth and native temperature kind, TDB/UTC according to each frozen temporal support, actual atoms and seconds for cohort evidence. |
| Frame/time | Native ME/DE421 planetocentric east-positive coordinates; exact site/quantity/depth/time-support identity required. No averaged-map date is interpreted as a simulation epoch. |
| Deterministic/stochastic controls | Approved native-pixel SHA-256 stratified selection; shared model elements retained; existing 16-seed dt/weight qualification and fixed analytical controls. |
| Procedure | Preflight without residual access; verify native site selection against candidate hashes; require all formal predecessors and independent operators; compare matched rows without interpolation or renormalization; evaluate actual cohort hazards; render only evaluated rows. |
| Metrics | Bounded visibility/PSR/temperature residuals under predeclared review criteria; geometric PSR separately from thermal retention; exact/matched row counts; first survival/occupancy/restricted residence and numerical convergence. |
| Acceptance | Original roadmap limits and explicitly supplied pre-holdout review criteria only. No new polar thermal tolerance is inferred from I16's non-PSR 15 K criterion. |
| Expected PASS | Every required real comparison and kinetic identity was evaluated and met its frozen criterion with complete formal predecessors. |
| Expected failures | Quantity/depth/time mismatch, modified native selection, missing strata hidden, shared model support called independent, unsupported observation operator, missing actual hazards, numerical/physical criterion miss. |
| Required data | D02 AVGVISIB/LPSR and D03 PRP plus exact generation/temporal/geometry metadata. D02 is currently PARTIAL. |
| Data hashes | All descriptors are path/SHA-256 pairs and are retained before science rows are opened. Altered bytes are ERROR. |
| Predecessor gates | Same-campaign formal I13-I19, exact D02/D03 READY packages, real linked runtime evidence, independent numerical/thermal/operator qualification. |
| Artifacts | Polar preflight/result JSON, immutable descriptors/identities, exact comparison CSV/PNG/EPS, conditional kinetic rows, raw-source and model manifests; no plot on a skipped run. |
| Status semantics | PASS evaluated/met; FAIL evaluated/missed; ERROR malformed/corrupt evidence; SKIPPED missing scientific capability/data/predecessor. Protocol fixtures never become observational validation. |

Source roadmap SHA-256 at inventory:
`f1a2f9265f97f36da9ef9f6a17b3c83783d9f4b51882e6d81081ef7d567757cd`.
External roadmap SHA-256:
`d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4`.
The source copy controls approved execution/development order; the common
scientific gates control qualification. Native labels and external AGENTS
control data interpretation and immutable raw files.

No real comparison will be scored in this development work. Tests use tiny
temporary metadata/arithmetic fixtures and identify them explicitly as
infrastructure verification. They do not provide a passing D02/D03 package,
new observation, guessed epoch, assumed uncertainty or threshold approval.

## Site-bound retention extension, specified before implementation

The read-only R2.7 audit found two missing reporting/binding checks. A generic
I19 cohort PASS did not prove that the cohorts started at the frozen I20
native sites, and thermal-retention labels were visible only in optional
inventory checkpoints. This extension changes neither a physical law nor an
acceptance threshold.

* **Production state:** the native I20 surface exporter must report each
  species' residence time and the production diagnostic
  `residence_time_s >= cold_trap_retention_timescale_s`, using the same
  `ResidenceTimeS` kernel and current kinetic surface state as emission. The
  explicit configured timescale is reported with those values. Only depth-zero
  rows carry these surface diagnostics. Subsurface rows keep them empty;
  subsurface temperature is not substituted into a surface release law.
* **Native CSV extension:** `surface_retention_timescale_s`,
  `na_surface_residence_time_s`, `na_surface_thermal_retention`,
  `ar_surface_residence_time_s`, `ar_surface_thermal_retention`. A positive
  infinite residence time is a valid zero-release limit. NaN, negative values,
  and unlabelled/inferred timescales are errors. Python checks the diagnostic
  relation but never reimplements a residence or binding-energy formula.
* **Site identity:** the ordered, possibly repeated logical-element list in
  every finest I19 input must match the frozen I20 site's order. Cohort index
  maps to that exact ordered site, and every actual cohort/hazard row must
  match its Na/Ar species, index and element. Shared elements retain separate
  native-site/cohort labels and are explicitly correlated support.
* **Runtime identity:** all sixteen approved finest seeds require their own
  hash-bound native manifest, receipt, input, cohort and hazard files. Their
  executable, layout, source/build identity, elapsed interval and explicit
  start epoch must match the qualified numerical design. Their interval is
  reported separately from the annual thermal/illumination supports; no
  false equality of those physical questions is imposed.
* **Time and depth:** surface kinetics are attached to exact native query
  ID/site/element/TDB time/support/depth-zero identities. The configured
  reporting timescale must match both the reviewed protocol and native input.
  Every required surface query has both species. An instantaneous retained
  label means only that its residence time exceeds that explicit timescale;
  it does not establish permanent retention or annual thermal qualification.
* **Dual classification:** geometric PSR comes from terrain-only model-ray
  occurrence. The benchmark's Earth-policy visibility occurrence is retained
  separately. Neither is used as a shortcut for the thermal-retention label.
* **Artifacts:** exact `surface_retention.csv` with per-query/per-species
  state, terrain-only PSR and benchmark visibility; exact
  `site_kinetics.csv` with per-seed/site/species actual survival, restricted
  first residence and surface occupancy. JSON preserves the identities and
  explicitly disclaims permanent trapping and independent shared-cell
  sampling. Missing new native capability is SKIPPED; malformed, mismatched
  or changed native evidence is ERROR; evaluated kinetic identities retain
  FAIL semantics. Metadata/arithmetic fixtures are never validation evidence.

The native exporter extension is coordinated with the parent agent. Until it
exists, the Python scorer must identify missing native surface-kinetics fields
as unavailable capability, rather than deriving them from observations or
assigning a default temperature/residence time.
