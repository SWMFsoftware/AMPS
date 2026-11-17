# M3 charge exchange: researched channel and implementation proposal

**Status: RESEARCH / PROPOSED; no channel activation or product-kinetics
decision is approved by this document.** Research date: 2026-11-17. The
production charge-exchange sources were not changed. Downloaded references
are research artifacts, not a READY validation-data package, observations,
or a passing U17/I24 result.

The recommended first scope is a small registry containing proton capture
from ground-state He and Na, with independently verified total-rate kernels.
Keep Ne, Ar, excited targets, alpha channels, and full daughter transport OFF
until the corresponding data and physical scope are reviewed. The available
total cross sections support a neutral-loss calculation over declared energy
ranges; they do **not** determine a unique pickup-ion/energetic-neutral
velocity distribution. Approval is needed before enabling either a reduced
neutral-loss bath or a complete two-product collision model.

## 1. Controlling contracts and the actual starting implementation

Source revision inspected: `17d49e00d1d2ccf33b06d11a3e6e5ffef34ac810`.
This is a dirty, actively developed checkout; the file hashes in the research
inventory identify the inspected bytes more precisely than HEAD alone.

| Contract | SHA-256 at inventory capture |
|---|---|
| `srcMoon/development_and_validation.md` | `712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17` |
| `/data/vtenishe/moon_validation_data/development_and_validation.md` | `d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4` |

The roadmaps differ. The source document contains approved resource,
development-order, M2 protocol and initial-state amendments, plus the approved
2026-11-17 M3 development extension; the external document retains the earlier
roadmap. Source amendments control permitted implementation timing. Their
common scientific gates, the frozen U/I meanings, and authoritative native
data metadata control qualification. Early M3 work does not close M1/M2,
authorize observational scoring, or relax formal predecessors. The exact diff
can be inspected with:

```sh
diff -u /data/vtenishe/moon_validation_data/development_and_validation.md \
  srcMoon/development_and_validation.md
```

### 1.1 U17 currently contains no documented working lunar channel

R3.4 says that an existing documented U17 channel can be enabled first. That
statement does not match the implementation. The current
[U17 contract](U17_charge_exchange/reference/acceptance.json) instead says:

```json
"implementation": "not-implemented",
"skip_reason": "No Moon charge-exchange channels or authoritative channel table are configured."
```

Its runner delegates to the common unavailable-capability path. The README
likewise lists charge exchange as SKIPPED. I24 is a planned ID in
[`test_runner.py`](test_runner.py), without a charge-exchange linked runner.
Correct the roadmap sentence when incorporating the reviewed proposal; do
not manufacture a historical channel or count a core hydrogen test as U17.

### 1.2 The generic AMPS API is hydrogen-only and returns a lifetime

[`ChargeExchange.h`](../../src/models/charge_exchange/ChargeExchange.h) declares
`MaherTinsley::LifeTime`, `LindsayStebbings::LifeTime`, and a compile-time
`ChargeExchange::LifeTime` wrapper. The header's introductory comments mention
`Init`/`GetCrossSection`, but these are not declared interfaces. The actual
Lindsay–Stebbings entry begins with:

```cpp
// src/models/charge_exchange/ChargeExchange_LindsayStebbings.cpp
if (spec != _H_SPEC_ && spec != _H_ENA_V1_SPEC_ &&
    spec != _H_ENA_V2_SPEC_ && spec != _H_ENA_V3_SPEC_)
  return 1E+100;
```

The Maher–Tinsley implementation has the same species restriction. Thus
passing a He, Ne, Ar or Na species ID does not select that atom's reaction:
it returns the established unsupported/no-reaction sentinel in seconds.

`ChargeExchange.dfn` selects Lindsay–Stebbings by default;
`ChargeExchange.pl` changes the generated selector. This is separate from
installing a Moon reaction processor or supplying projectile moments. Both
implementations accept velocities in m/s, proton temperature in K and
projectile number density in m^-3. They expose no daughter species, recoil,
angular sampling, product weights, or reaction-event ledger.

The present Lindsay–Stebbings implementation computes a Maxwellian mean
relative speed and then uses `n * mean_speed * sigma(mean_speed)`. For a
velocity-dependent cross section, that is generally different from integrating
`n * <sigma(g) g>`. Preserve this existing hydrogen mode and its regressions;
do not silently reuse this approximation for the new lunar channels.
The original Lindsay–Stebbings paper treats H/O/N2/O2 reactions, not a
He/Ne/Ar/Na parameter set. Its discussion also distinguishes total rates from
scattering and internal-state information.
[Lindsay & Stebbings 2005, DOI 10.1029/2005JA011298](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2005JA011298).

### 1.3 Existing tests and actual production connection

The optional shared-core suite `CORE-T6-REACTIONS-NATIVE-001`, alias
`CORE-CHARGE-EXCHANGE-001`, contains three hydrogen checks: Lindsay,
Maher and the generated selector. Its source
[`physics_t6_reactions.cpp`](../../test/srcCoreUnitTests/physics_t6_reactions.cpp)
explicitly states that charge exchange lacks product, weight, interpolation
and Maxwellian-quadrature APIs. These tests protect an existing hydrogen
contract; they do not verify a lunar charge-exchange channel.

The current lunar chemistry entry is
`Moon::Chemistry::ProcessParticle` in [`Moon.cpp`](../Moon.cpp). At audit time
it implements neutral-Na photolytic loss and delegates other species to the
generic photochemical processor. Its explicit refusal to create unqualified
Na+ is relevant to the required product decision:

```cpp
if (_NA_PLUS_SPEC_ >= 0)
  exit(__LINE__, __FILE__,
       "Neutral-Na chemistry cannot silently create Na+ with unqualified "
       "product weights. Use the neutral-Na production configuration; "
       "qualify Na+ product weights and kinetics before enabling that case.");
```

A new channel must be dispatched through the production chemistry path with
the same physical weights, time convention and event accounting as an
ordinary run. A lifetime library call in a standalone probe cannot establish
that connection. Other M3 work is occurring concurrently; this inventory is
the starting contract, not a claim that subsequent unreviewed changes already
provide charge exchange.

## 2. Authoritative data actually found

The reaction notation below retains the projectile identity and target
electronic state. `H*` denotes a potentially excited neutral hydrogen product;
it is not permission to assign a particular state or velocity.

| Candidate reaction | Primary numerical evidence and supported domain | Initial recommendation |
|---|---|---|
| `H+ + He -> H* + He+` | ORNL-6086/V1, A-32; total capture, degree-8 Chebyshev fit and printed table, 99–11,000,000 eV/u, cm² | First total-rate kernel candidate; original publication controls, with source QA below |
| `H+ + Na(3s) -> H* + Na+` | Morgan et al. 1985, Table 2; all-H-state capture table, Na entries at 1–200 keV/u, units 10^-16 cm² | First Na total-rate candidate over those entries only; freeze an interpolation rule before implementation |
| `H+ + Ne -> H* + Ne+` | Nakai et al. 1987 evaluation is identified by QST, but an exact numerical coefficient/range record was not acquired | OFF pending original coefficients, states, range and independent QA |
| `H+ + Ar -> H* + Ar+` | Same evaluated publication and acquisition limitation as Ne | OFF pending numerical QA; never substitute Ne or H cross sections |
| `He2+ + He -> He+ + He+` | ORNL A-98 single capture, 2.4e-4–2e6 eV/u | Later explicitly selected alpha channel; needs alpha-resolved forcing and product contract |
| `He2+ + He -> He + He2+` | ORNL A-104 double capture, 7–4e5 eV/u | Later separate channel; exchanges slow/fast neutral populations and cannot be represented as a simple He number loss |

These rows are a proposal, not an activation list. The numerical ranges are
source ranges, not asserted coverage of an arbitrarily specified plasma.
The two alpha rows are not substitutes for the proton row. No proton/alpha
abundance, helium isotope mixture, metastable fraction or solar-wind parameter
has been assumed.

### 2.1 Recommended helium reference

[Barnett 1990, ORNL-6086/V1](https://www.osti.gov/servlets/purl/6570226), A-32
(PDF page 48), provides a total capture fit and table. It gives 40% accuracy
below 5 keV/u and 20% above, and reports 10% RMS/18.7% maximum fit-to-table
deviation. These source accuracies are not stated as Gaussian 1σ errors.
A-98 documents theoretical/extrapolated low-energy content below 250 eV/u;
A-104 documents theoretical values below 100 eV/u. Preserve those distinctions
if alpha channels are later selected.

Use the original A-32 coefficients and Appendix I formula, not a transcription
of a web example. In the target frame the energy variable is per projectile
mass unit. The coefficient fit is a specific, approximate physical law; a 1%
implementation check against an independent evaluation of that fit would not
make its experimental accuracy 1%.

The official [IAEA D120144 record](https://db-amdis.org/collisiondb/datasets/120144/)
identifies A-32, Target frame, eV/u and cm², but quotes 50% uncertainty and
rounded coefficients. Its published Python example initializes the recursion
with `pcf[9]`, which its own metadata identifies as the lower energy bound;
the Fortran version starts with the ninth **coefficient**. This apparent
indexing error must not be copied. A direct degree-8 evaluation of its rounded
coefficient list agrees with its 50 retained central source samples to a
maximum relative difference of 0.001989; the original PDF remains controlling.

### 2.2 Recommended sodium reference and its limits

[Morgan, Olson, Schlachter & Gallagher 1985, DOI 10.1063/1.555752](https://srd.nist.gov/jpcrdreprint/1.555752.pdf?download=1),
Table 2 on printed page 984, distinguishes capture into all H states from the
2s-only Table 3. The Na column has eleven populated knots from 1 to 200 keV/u;
blank entries below 1 keV/u are missing data, not zero cross sections. For
example, its 1 keV/u total is `43 * 10^-16 cm²`. Experimental accuracy
assessments are reference-dependent; a uniform invented error bar is not
appropriate. The paper discusses detector acceptance and excited-state
contamination as possible systematics.

A cold proton speed of 400 km/s corresponds to approximately 0.83 keV/u,
**below this table's Na support**. This is a dimensional illustration, not a
chosen lunar input or measured wind speed. Do not extend this table to that
case merely because it is close to the boundary. Seek the exact lower-energy
evaluated fit or another reviewed primary dataset if such a case is required.

Two tempting substitutes are scientifically different quantities:

- [IAEA D107581](https://db-amdis.org/collisiondb/datasets/107581/), from the
  Igenbergs et al. database, explicitly combines single capture and direct
  ionization of Na(3s). It cannot be assigned entirely to charge exchange.
  [Igenbergs et al. 2008, DOI 10.1016/j.adt.2008.06.001](https://www2.iap.tuwien.ac.at/www/_media/atomic/publications/igenbergs_andt08.pdf)
  describes this combined proton electron-loss quantity. Its other Na-state
  records do not supply missing charge-exchange branching fractions.
- [IAEA D119990](https://db-amdis.org/collisiondb/datasets/119990/) explicitly
  labels `H(1s) + Na+(2p6)`. Its twelve native numerical samples span
  4–25 keV/u. It is not an all-H-state cross section, and the inspected page
  lacks a sufficiently explicit original reference to adopt it automatically.

### 2.3 Other references found, and reasons for deferral

[Nakai et al. 1987, DOI 10.1016/0092-640X(87)90005-2](https://www.sciencedirect.com/science/article/pii/0092640X87900052)
provides evaluated analytic total charge-transfer formulas. The
[official QST index](https://www.qst.go.jp/site/jt60-english/6485.html) explicitly
lists proton capture from He, Ne and Ar. Its linked public assets are graphs;
this audit did not recover a reviewed numerical Ne/Ar parameter table and
therefore gives no invented valid range or coefficients. The
[QST metal-vapor evaluation](https://www.qst.go.jp/site/jt60-english/6493.html)
identifies a 1988 analytic evaluation for Na; its numerical coefficients remain
an acquisition task. No curve was digitized.

[Courbin et al. 1990, DOI 10.1088/0953-4075/23/21/023](https://doi.org/10.1088/0953-4075/23/21/023)
is a candidate state-resolved/differential Na source for 0.5–5 keV. The
[authors' institutional record](https://oa.upm.es/61369/) was identified, but
the full numerical differential distribution was not acquired or qualified.
It therefore does not yet specify a production scattering sampler.
[Le, Liu & Lin 2003, DOI 10.1103/PhysRevA.68.012705](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.68.012705)
studies 1.7–40 eV and reports disagreements with another low-energy calculation.
Do not splice that region into the keV table without an explicit evaluation
and a state/energy-defect specification.

The alternative [IAEA D120235 helium fit](https://db-amdis.org/collisiondb/datasets/120235/)
also needs clarification: its comments use eV, the column uses eV/u, its
displayed lower limit and native first sample differ, and direct evaluation of
the displayed rounded log-polynomial coefficients does not reproduce its
embedded samples. This is a metadata/numerical QA problem, not permission to
alter the coefficients or choose whichever curve agrees with AMPS.
`research_fit_qa.py` reproduces these checks: the latter discrepancy reaches
240.7 in relative units over its 49 in-range native samples. Its result is
source-representation QA, not a reaction-kernel or physical validation score.

## 3. Local data readiness and retrieval evidence

There is no approved charge-exchange Dxx package in the current dependency
map. No reference was installed in the canonical data repository.

D06's `drivers/README.md`, `provenance.json`, and `qa_report.json` were read.
The QA status is **ERROR**: zero raw CDF files, zero plasma-driver rows, no
campaign acquisition intervals and no probe-selection rules. The current
CSV is schema-only. Its declared bulk velocities are GSE; AMPS's lunar mover
uses LSO and its surface uses the declared Moon body frame. A documented,
epoch-specific transformation is required before using measured vectors.
D06 deliberately supplies no alpha/proton ratio and does not infer proton
temperature from its ion thermal-velocity field. Consequently neither a
measured proton Maxwellian nor an alpha bath is presently available for a
formal event-forced calculation.

Research artifacts are retained under:

```text
test_output/srcMoon/m3-charge-exchange-research-20261117/
```

| Retained source/artifact | SHA-256 |
|---|---|
| `ornl_redbook_v1_osti.pdf`, original DOE/OSTI PDF | `3bb3fdb4be8f0f806d378083be8798cf3ce347714c3d06d58f90ee483b10c030` |
| `nist_morgan_1985_retry2.pdf`, official DOI-linked NIST PDF | `7c412cb93f827b3412b6da000e674c63237cacb0e201733a6f5a4f431359952e` |
| `collisiondb_D120144.html` | `95a103de35e543afee8c9ed23e3743fbabbbe96e2ca69ee20f3e9ff38db731e4` |
| `collisiondb_D107581.html` | `7d71f738ba904a2337abe662c38b7afaf0a35122cd2b8a211184a68c160c8d1f` |
| `collisiondb_D119990.html` | `5ed771a9dfeee2f8e42433d450d92b346c6b0ddbb5c22d45f49e54592d844313` |
| D06 `qa_report.json` | `589342d3054c080359567290591e786adfe3e464bb1105bb3b838c6a75a8a2aa` |
| D06 `processed/artemis/plasma_driver.csv` | `06179073736dc8c7207cc299d2c5ec05c21f3a81826e221acc04d82d3f08181d` |

`inventory_manifest.json` records inspected source/data hashes. The
`retrieval_manifest_*.json` files retain original URLs, HTTP status, byte
counts, content types, UTC retrieval times and hashes. Some filenames end in
`.pdf`/`.txt` but contain retained HTML failure responses; the manifest and
`%PDF` signature, not the suffix or HTTP 200 alone, determine whether a
publication was obtained. `ornl_redbook_v1.pdf` is the blocked IAEA response;
`ornl_redbook_v1_osti.pdf` is the actual publication. The unadorned NIST URL
initially returned 403; its `?download=1` form succeeded. A separate older
NIST PDF path returned a valid but different document and was not used.

IAEA API count queries found 87 H+/He and 9 H+/Na HCX records, and zero for
the exact H+/Ne and H+/Ar queries used. Counts include excited states and
composite processes; they are not counts of usable channels. Zero query
matches do not establish absence from the literature. The API archive export
and individual Download endpoints returned remote HTTP 500 with
`No space left on device`; that is the service's failure, not AMPS disk
exhaustion. Native numerical arrays were nevertheless recovered from the
retained pages' Bokeh JSON without digitizing an image, smoothing or fitting.
Their summary is explicitly research-only in
`embedded_arrays_research_summary.json`.

## 4. Physical definition required before implementation

### 4.1 Energy, rate and out-of-range handling

For a neutral velocity `v` and an explicitly declared projectile distribution
`f_p(u)` normalized to one, the total loss hazard is

```text
g = |u - v|                                      [m/s]
nu_X(v) = n_p * integral f_p(u) sigma(E(g)) g d^3u [1/s]
k_X(v) = nu_X(v) / n_p                           [m^3/s]
```

For a **nonrelativistic** projectile kinetic energy in the target rest frame,
`E_lab = 0.5*m_projectile*g*g`. A center-of-mass energy instead uses the
reduced mass. An eV/u table requires division by the projectile's declared
mass in unified atomic mass units after conversion to eV; do not use the
lunar target mass or silently treat every projectile as hydrogen. Cross
sections in cm² acquire a factor `1e-4` to obtain m²; the NIST table's
`10^-16 cm²` multiplier becomes `10^-20 m²`. An input in keV/u needs the
additional factor 1000 before an eV/u fit is queried.

The numerical cross-section kernel may accept the archive's target-rest
energy directly over its declared domain. This does not make the classical
velocity-to-energy conversion accurate throughout He A-32's 11 MeV/u upper
range. A velocity-driven application must either document a nonrelativistic
domain/error bound or implement the corresponding relativistic frame/flux
conversion before claiming the full source range. No such transport or
high-energy velocity mapping is approved by this rate-kernel proposal.

**Proposed policy:** a supported channel has explicit validity metadata and a
typed out-of-domain result. An enabled run encountering an unsupported
energy must stop with a diagnostic identifying channel, energy, source range,
and remedy. It must not silently return zero, clamp to an endpoint or extend
a fit. A deliberately disabled channel is OFF and recorded as such. Zero
physical projectile density has zero hazard; negative/nonfinite density or
invalid distribution parameters are configuration errors, not OFF modes.

A Maxwellian has unbounded tails. A finite cross-section range therefore
does not define its complete rate integral by itself. Integration over only
the available range is a partial rate, not a fully qualified Maxwellian
rate. Either obtain adequate authoritative coverage plus a documented tail
error bound, or approve an explicitly bounded physical distribution. A cold,
monoenergetic beam is suitable for a declared local reference experiment;
it must not be relabeled as a measured thermal solar wind. The electron
Maxwellian/Voronov decision does not settle this ion-distribution decision.

### 4.2 Total capture does not specify product velocities

For an actual event, the reacting projectile distribution is proportional to
`f_p(u) * sigma(g) * g`, rather than an unweighted draw from `f_p`.
State-resolved integral cross sections are needed for final electronic-state
branching, and differential cross sections for the angular distribution.
If the channel is nonradiative and its internal states are specified,
two-body kinematics must preserve center-of-mass momentum and satisfy
`K_after = K_before + Q`, with the channel's signed energy defect `Q`.
Radiative capture is a different channel with a photon energy/momentum
ledger. Isotropic scattering, forward scattering, velocity swapping, and
retaining the target velocity are different approximations; none follows
from a total-cross-section table alone.

In particular, do not assume every Na capture makes ground-state H. The
state-dependent threshold and energy defect differ between final H states.
Assigning one arbitrary `Q` to an all-state total can appear to close a
numerical energy budget while implementing the wrong physical reaction.

### 4.3 Bath versus explicit projectiles

Two scientifically distinct implementations are possible:

1. **Prescribed bath, neutral survival only.** Use a declared proton
   distribution to calculate total neutral loss. Record the removed neutral
   statistical weight and the reaction's stoichiometric transfer to named
   unresolved product reservoirs. Do not claim predicted pickup-ion/ENA
   velocities, closed kinetic-energy/momentum transfer, or a full U17/I24
   product gate. This is a reduced transport question, and requires explicit
   user approval because the existing milestone contract is broader.
2. **Full products.** For a prescribed bath, retain the consumed bath
   projectile and both products in a signed source ledger; for explicit
   projectiles, update/remove the actual projectile consistently. In both
   cases resolve state and scattering, statistical weights and collision
   sampling. There is one incoming proton and one incoming lunar atom per
   single-capture event. A virtual bath is an open-system boundary, not a
   reason to ignore its mass, charge, momentum or energy transfer.

**Recommendation:** implement/reference-test the small total-rate registry
first, keep production channels OFF, and request the scope decision below.
Do not force a same-velocity ion into the existing neutral-Na processor to
claim M3 completion. A neutral-loss development feature can be useful, but
it must retain the unclosed full-product gate explicitly.

## 5. Concrete decisions requested before channel activation

These are proposed choices, not approvals inferred from the earlier request
to research charge exchange.

| Decision | Recommendation and consequence | Relevant implementation evidence |
|---|---|---|
| Initial channels/reference law | Review `H+ + He` with the original ORNL A-32 fit, plus `H+ + Na(3s)` with the NIST all-state table. Propose log-log interpolation only between its positive populated knots; freeze this rule and supported range before coding it. Ne/Ar and alpha channels remain OFF. | The generic `LifeTime` species restriction above prevents current He/Na support; the source table distinguishes total capture from other processes. |
| First physical product scope | Prefer a rate-kernel implementation now, with production OFF until choosing either an explicitly limited neutral-loss bath or full product transport. A loss-only mode cannot close the full current U17/I24 product contract. | `Moon::Chemistry::ProcessParticle` explicitly rejects unqualified Na+ weights/kinetics. |
| Projectile distribution/forcing | For local verification, an explicit monoenergetic beam entirely inside source support is sufficient. A physical lunar run must declare its proton distribution and parameters or supply qualified measured forcing; any Maxwellian tail treatment must be reviewed. No default wind, temperature or alpha fraction is proposed. | `LifeTime` currently requires a proton temperature and density; D06 has no acquired driver rows and no alpha-resolved product. |
| Full kinetics, if selected | Obtain complete state-resolved/differential data and freeze state energies, angular law, bath/explicit treatment and product-weight strategy before implementation. Do not choose isotropy or velocity inheritance to fill missing data. | The existing `ChargeExchange.h` has no product/event API; U17's purpose explicitly includes product momentum and budgets. |

If the requested production use requires slower protons than Na's 1 keV/u
boundary, recovering a lower-energy evaluated **pure-capture** source is an
additional data prerequisite; approval cannot turn a missing table into
available physics. Likewise, lack of Ne/Ar numerical acquisition should be
reported as missing data, not filled with scaled He/Na values.

## 6. Proposed implementation and verification contract after decisions

No permanent IDs are added. Extend the existing U17/I24 meanings and keep
their currently required predecessors. The following is an implementation
plan, not a record that these tests have passed.

| Contract field | Proposed content |
|---|---|
| ID/name/purpose | U17 reaction kernels; I24 linked solar-wind-ion charge exchange; distinguish rate, wiring, products and budget subgates |
| Production capability/entries | A Moon channel registry under `srcMoon/species`, production rate evaluator and chemistry dispatch; existing generic hydrogen APIs remain unchanged |
| Configuration/inputs | Explicit OFF default; stable channel ID, reactants/products/states, source identifier/hash, law/range/energy convention, bath/explicit scope, units and distribution definition in `.in` |
| Independent reference/oracle | Original published coefficients/table knots; separately coded direct polynomial/interpolant and independent quadrature; exponential survival for constant hazard; elemental/charge and kinematic invariants where products are defined |
| Units/frame/time | SI at API boundary; table source units preserved; projectile and neutral velocities in one named frame at the same declared event epoch; no implicit GSE-to-LSO conversion |
| Deterministic/stochastic controls | Fixed recorded seeds; deterministic table/boundary and zero-density cases; at least 100000 events for a supported probabilistic branch, as the existing U17 contract requires; record represented weights and confidence intervals |
| Procedure | Verify parsing/provenance before rate evaluation; test knots and interior points, exact boundaries, unsupported channels, density scaling and relative-frame invariance; compare independent rate integration; run linked constant-bath survival; then full products only when supported |
| Metrics/acceptance | Preserve R3.4 direct-rate error <=1%, ensemble lifetime error <=3%, exact product identity and charge/elemental stoichiometry; apply predeclared channel branch and momentum/energy criteria before activation; never equate numerical tolerance with laboratory uncertainty |
| Expected PASS | Only declared evaluated subcriteria pass; full U17/I24 require all their applicable original criteria, not a source-only or rate-only result |
| Expected failure modes | Wrong initial/final state, unit/frame conversion, composite process mistaken for capture, invalid energy/tail handling, rate bias, wrong daughter weights, missing bath ledger, identity/hash corruption |
| Required data/hashes | Reviewed ORNL/NIST reference package with reproducible extraction, original metadata and frozen law; additional differential/state evidence for products; D06 only for measured-event forcing |
| Predecessors | U17 retains its U15/U16 contract; formal I24 retains the roadmap's earlier gates and R2.6 budgets; the early M3 exception permits supporting development only |
| Artifacts | stdout/stderr, exact `.in`, `result.json`, `run_manifest.json`, source/build/executable identities, source hashes, seeds/MPI/OpenMP layout, rate/survival/branch/budget CSV and PNG/EPS diagnostics |
| Status semantics | Missing channel/data/decision/production scope => SKIPPED; malformed data/hash/runtime => ERROR; evaluated scientific miss => FAIL; numerical/reference agreement alone is not observational VALIDATION |

The one-hour development campaign can contain deterministic rate tests,
bounded independent integration and a linked fixed-bath survival test, once
implemented and measured. Use eight MPI ranks for linked work. Pilot the
runtime rather than assume that all product distributions or a global
multi-species simulation fit the deadline. Larger event-forced statistics
remain an explicitly selected cluster campaign with aggregate resource
accounting. No physical horizon or confidence criterion is shortened to meet
the development budget.

Before merging a new Moon channel, rerun the protected shared hydrogen suite
and the lunar baseline. Preserve the requirement to remove `build`, configure
normally, compile with `make -j`, and obtain a zero-byte `test_Moon.diff` from
`make -j test_Moon` when a native rebuild is needed. This research package did
not configure or rebuild the application.

## 7. Reproducibility and work-package result

Commands actually used include:

```sh
git rev-parse HEAD
sha256sum srcMoon/development_and_validation.md \
  /data/vtenishe/moon_validation_data/development_and_validation.md
rg -n 'ChargeExchange|LifeTime' src/models/charge_exchange \
  test/srcCoreUnitTests/physics_t6_reactions.cpp
pdftotext -layout \
  test_output/srcMoon/m3-charge-exchange-research-20261117/nist_morgan_1985_retry2.pdf \
  test_output/srcMoon/m3-charge-exchange-research-20261117/nist_morgan_1985_text.txt
pdftotext -layout \
  test_output/srcMoon/m3-charge-exchange-research-20261117/ornl_redbook_v1_osti.pdf \
  test_output/srcMoon/m3-charge-exchange-research-20261117/ornl_redbook_v1_text.txt
/usr/bin/python3 \
  test_output/srcMoon/m3-charge-exchange-research-20261117/decode_embedded_arrays.py
/usr/bin/python3 \
  test_output/srcMoon/m3-charge-exchange-research-20261117/research_fit_qa.py
python3 -B srcMoon/test/U17_charge_exchange/test.py \
  --output-dir test_output/srcMoon/m3-charge-exchange-research-20261117/U17-current-contract
```

The last command returned **77 / SKIPPED** with its actual missing-capability
reason. U17 executed: 1; PASS: 0; FAIL: 0; ERROR: 0; SKIPPED: 1. I24 was not
executed. The PDF/text/JSON inspections are research QA, not U/I test results.
No observation residual, calibration, formal validation score, or comparison
figure was generated. No reference was approved or frozen by this work.

Research HTTP requests used Python `requests` with exact URLs and selected
IAEA query dictionaries retained in the retrieval scripts/manifests. The
official client API is documented in the
[IAEA-NDS pycollisiondb source](https://github.com/IAEA-NDS/pycollisiondb).
`decode_embedded_arrays.py` only decodes the original base64 float64 source
arrays and validates shape/finite values; it does not sample a plot image.
It is retained beside the research artifacts, not installed as a validation
processor. Any adopted processed reference must subsequently be generated by
a committed, reviewed processor with complete provenance and QA.

Changed source-tree file: this proposal only. Required follow-up documentation
after review: correct R3.4's nonexistent-U17-channel claim, add the selected
channel/data contract to the U17/I24 specifications and component inventory,
and document the approved physical scope in relevant READMEs/STATUS. Missing
production capabilities remain the lunar registry, correct rate integration,
production dispatch, and qualified product/event/weight treatment. Missing
data remain Ne/Ar numerical references, full differential/state support and
any requested measured projectile forcing.
