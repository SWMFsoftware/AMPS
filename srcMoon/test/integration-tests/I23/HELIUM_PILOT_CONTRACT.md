# I23 stationary helium development pilot

| Item | Contract |
| --- | --- |
| Purpose | Connect the D05 rate and its analytic slope to the actual shared Maxwellian moment and conservative planner for the approved stationary-helium reference. |
| Production entry points | `HeliumElectronImpact::Prepare`, `Evaluate`, and the guarded `Compose` API, using the shared Voronov table, `ReservoirProducts` and `ReactionExecution::MakeEventPlan`. Future linked ownership uses the existing Moon controller and `ParticleTransaction`. |
| Configuration and inputs | Explicit ⁴He neutral and ⁴He⁺ ground states, free electron, 10 eV thermal energy, densities 0, 10⁶ and 10⁷ m⁻³, physical state profile and frozen D05 fit. Prescribed isotropic Maxwellian with zero target velocity and drift in one named instantaneous inertial frame. No measured D06 forcing. |
| Units and time | eV for D05 thermal energy and fit threshold, SI J internally, densities in m⁻³, coefficient m³/s, rates s⁻¹, exact charge in C and SI signed reservoir moments. Use the accepted parent event interval in seconds, with explicit frame and time-scale identities. |
| Independent reference | Decimal implementation of the original fit and a centered high-precision temperature difference for its slope. The Maxwellian identity supplies the model mean. Independent integer count and SI mass, momentum and energy checks apply to a committed packet when its model evidence is admitted. |
| Procedure and metrics | Authenticate D05 and all state profile sources. Use the production provider to evaluate rate, slope, incident conditional mean and physical versus fit threshold. Check zero density, linear density response, refusal of drift/moving-target use, and refusal to compose events when the model-error evidence is missing. Report the rate/moment discrepancies and zero activation count. |
| Controls and acceptance | Deterministic probe with no RNG or MPI. The rate uses the stricter U16 1e-12 relative limit. Slope/mean arithmetic uses the documented independent finite-difference bound. Actual I23 still requires the existing direct coefficient/lifetime, numerical, ensemble, and ledger criteria and U16/I13–I20. |
| Data and outputs | Exact D05 U16 source manifest and pinned ⁴He state profile. Preserve build/stdout/stderr logs, profile/input, source, executable, data hashes, result and numerical CSV. Missing capability is SKIPPED, infrastructure failure ERROR, evaluated criterion miss FAIL, met supporting criterion PASS. |

The physical operator needs an authenticated threshold/state and slope-error
crosswalk, and explicit observable bounds for omitted recoil and mass
convention. The guarded composer checks the whole signed bath transfer from
one incident to two outgoing electrons. A diagnostic fitted mean alone supplies
no derivative-uncertainty or model-admission evidence. Keep the separate
[formal I23 acceptance](reference/acceptance.json) and [reservoir model
contract](../../M3_RESERVOIR_PRODUCTS_CONTRACT.md) in force.

For the optional Kim source comparison, also report the fraction of the
10 eV Maxwellian helium ionization integral between the 100 and 2000 eV
incident-energy endpoints studied by Peterson, Beaty and Opal. The denominator
is the Kim fit integrated over its stated validity interval ending at 10 keV.
Evaluate below, within and above the experimental window, including a
separate quadrature of the entire interval for the closure check. Uncertainty
from those integrals is numerical. The fraction sets an upper limit on the
eligible collisions in this particular model. Individual sampled beam energies
and emitted-electron angle/energy coverage can reduce the experimental overlap.

**Update 2026-11-25.** The user approved a narrow-scope admission for this
stationary 10 eV reference: heavy-particle budgets and chronology only, with
conservative upper-end recoil tolerances. The linked coupling uses it through
`MoonHeliumImpactCohort`. Its tests follow
[ELECTRON_BATH_COUPLING_CONTRACT.md](ELECTRON_BATH_COUPLING_CONTRACT.md),
and the record is [HELIUM_ADMISSION_20261125.md](HELIUM_ADMISSION_20261125.md).
The pilot itself is unchanged; it still exercises the missing-evidence
refusal with an empty admission.
