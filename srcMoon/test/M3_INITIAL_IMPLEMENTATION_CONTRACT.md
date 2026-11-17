# Initial M3 implementation and research contract — 2026-11-17

The user approved M3 implementation/local supporting verification before
M1/M2 closure and an explicit Maxwellian/Voronov electron rate model. The
user also requested literature proposals for charge exchange and molecular
daughter kinetics before those choices are implemented. Formal qualification
order and every existing acceptance criterion remain unchanged.

## Inventory and scope before changes

`Moon.cpp::Chemistry::ProcessParticle` implements the existing neutral-Na
photolytic loss and explicitly refuses an unqualified Na+ configuration.
`Moon.h::ExospherePhotoionizationLifeTime` retains its historical Na lifetime
and lunar/Earth shadow gates. The current normal input has no D04 loader and
does not enable the PIC electron-impact executor. The dormant
`Moon::ElectronImpactIonizationRate` uses fixed cross sections and solar-wind
speed rather than a thermal electron distribution. U16/U17 currently SKIP;
no authoritative charge-exchange channel is configured. The roadmap's phrase
"existing U17 documented channel" therefore does not describe this checkout.

The shared PIC electron-impact executor independently samples product births
and parent removal, uses product time steps in its birth weight, and shadows
the cell-index variable in its product loop. Its presence is not evidence of
qualified product transfer or conservation. This work adds/configures the
atomic rate service and its existing Moon callback; it does **not** enable
that executor or claim ensemble/product budgets. Product conversion and
kinetics require their own implementation contracts and linked checks.

Available D04 and D05 live together under the canonical `photochemistry/`
package. Native `raw/electron_impact/voronov/cfit.dat` has 406 fits with exact
atomic stages and explicit validity bounds. Its `cfit.txt` specifies eV/keV
temperature bounds, cm3/s coefficient units; `cfit.f` supplies an independent
original implementation. Processed tables contain He/Na/Ne/Ar neutral fits.
The package README/provenance/QA and all 11 processed-product hashes were
checked. D06 is not qualified electron-density/temperature forcing. Nothing
in the canonical data repository is written by this work.

## Specification

| Field | Contract |
|---|---|
| ID/name | Extend U16 electron-impact coefficient verification. No new permanent ID; I21-I24 remain planned/formally gated. |
| Purpose | Implement and verify the documented atomic Maxwellian rate without asserting product kinematics or event forcing. |
| Production capability | Shared immutable Voronov table and rate evaluator used by the existing Moon electron-impact rate callback. |
| Entry points | `species/ElectronImpact` native-table loader/coefficient/rate APIs; `MoonInput` explicit immutable rate configuration; `Moon::ElectronImpactIonizationRate`. |
| Configuration | Optional `electron_impact_rate_model=maxwellian-voronov` with native coefficients path, density m^-3 and electron temperature eV. Absent/legacy selection preserves existing behavior; no PIC executor is enabled. |
| Inputs | Native D05 table/labels/Fortran and exact selected neutral atom/stage; invalid/unavailable inputs fail explicitly. |
| Independent oracle | Independently evaluated published expression at high precision; separately compare original Fortran within its own stated floating precision. No current output is frozen as physical truth. |
| Units | Rate coefficient m3/s; density m^-3; rate s^-1; electron temperature eV (optional exact-SI K conversion); native A cm3/s converted once. |
| Frame/time | Thermal scalar coefficients have no spatial frame. The explicitly prescribed Maxwellian bath is stationary; it is not measured ARTEMIS forcing. |
| Controls | Deterministic local tests; no particles, random draw, synthetic observation or tuned parameter. |
| Procedure | Strict native parser; all four neutral coefficients across their native temperature range; independent rate/units/density/domain tests; malformed/stale data failures; actual runtime parser/install/query path; normal native clean link and retained M0 baseline. |
| Metrics | Coefficient/rate relative errors, exact schema/atomic-stage/domain checks, zero-density behavior, transactional failure behavior, native configuration identity and default regression. |
| Acceptance | Existing U16 relative error <=1e-12; no criterion relaxation. Original Fortran's rounded conversion/single precision is reported separately, never used to relax the high-precision gate. |
| Expected PASS | Only evaluated coefficient/configuration/production-rate-callback contracts. This cannot establish I23 ensemble lifetime, M3 qualification or event validation. |
| Failure modes | Wrong stage/units/fit, malformed or incomplete table, invalid/nonfinite density/temperature, silent extrapolation, stale source/executable, misleading generic-executor activation. |
| Required data/hashes | D05 native table, description, Fortran; package provenance/QA and processed/raw file hashes retained per run. Missing dependency is SKIPPED; contradiction/corruption ERROR. |
| Predecessors | Local deterministic U16 can run under approved exception; formal I23 retains original numerical/environment/budget gates. |
| Artifacts | Native probe/source/config/compiler commands; stdout/stderr; result.json, numerical comparison CSV/PNG/EPS when plotted; run manifest, input/config snapshot and hashes; research proposals and completion report. |
| Status semantics | PASS evaluated/met; FAIL computed criterion missed; ERROR parser/build/runtime/corruption; SKIPPED unavailable prerequisite. Formal I21-I24 never promoted by a local rate PASS. |

Runtime comments/documentation must distinguish a configured rate callback
from an enabled and qualified particle reaction. Legacy Na lifetimes and
reviewed nightly references are not replaced in this package. No Maxwellian
parameters are silently assigned, no intermediate solar-activity law is
invented, and no molecular/charge-exchange reaction is enabled before review.

The selected development tests remain within the user's one-hour per-campaign
limit, with eight MPI processes for linked fixtures where appropriate. Clean
builds use an isolated checkout; remove its generated `build` before normal
Config.pl/make. Root Makefile/configuration and all earlier dirty changes
remain preserved. Source/data roadmap hashes and their operational differences
are recorded in `ROADMAP_RECONCILIATION.md` and the run evidence.


## Native all-species receipt correction discovered during verification

The first all-species fixture used human-readable Na+/He+/Ne+/Ar+ labels;
AMPS requires its declared `NA_PLUS`, `HE_PLUS`, `NE_PLUS`, `AR_PLUS`
identifiers. Correct the derived fixture only; no ion is invented. Inspection
also found `input/species.input` declares `mas =4.0*_AMU_` for He_PLUS.
`ReadSpeciesBlock` initializes masses to zero and ignores this unknown key,
so compiling that species would silently give zero mass. The minimal shared
input correction changes only `mas` to `mass`, preserving its already declared
4.0-AMU value. This is an implementation/data-key fix, not a new mass model.

The native U16 callback fixture now exports actual generated neutral/ion
masses [kg] and charges [C], and refuses nonpositive/nonfinite masses or wrong
charge signs. The full core regression preservation campaign runs separately
against the corrected shared species input; all previously passing case IDs
must remain PASS. The normal Moon M0/nightly is rerun with the correction as
well. Earlier failed/misconfigured fixtures remain retained as ERROR evidence.
