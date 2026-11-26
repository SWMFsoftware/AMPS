# I23 electron-bath coupling — supporting test contract

This contract specifies two supporting children of the existing I23 owner. It
assigns no new test ID and does not change the formal I23 contract in
[`reference/acceptance.json`](reference/acceptance.json). Formal I23 remains
`SKIPPED / NOT QUALIFIED`.

| Item | Specification |
| --- | --- |
| ID / name | I23 supporting children: `helium-electron-bath` (MPI provider ensemble) and `helium-pic-electron-bath` (linked native PIC matrix). |
| Purpose | Verify that the prescribed electron bath is coupled to the normal reaction controller with the correct hazard `n_e k(T_e)`, exposure clock, event chronology, product publication and signed electron reservoir. |
| Production capability | `Moon::Chemistry::HeliumImpactCohort` startup owner; `HeliumElectronImpact::{Prepare,Evaluate,Compose}`; shared `ReservoirProducts::ElectronImpactInheritedVelocity`; the `MoonReactionController` trajectory-dependent channel; `ParticleTransaction`; and the normal first- and second-order force-cut-surface movers. The new He/He⁺ `.in` startup keys are used in `amps_init`. |
| Production entry points | `HeliumImpactCohort::Prepared::RegisterForNormalStartup` → `NativeRegistration::SetBeforeParticleBuffer` → `amps_init` → `Controller::Install`; the per-segment `History` callback; the per-event `Kinematics` callback → `Compose`. The ensemble probe calls `Evaluate`, `AdvanceConstantHazards` and `Compose` directly. |
| Configuration | The approved controlled bath: isotropic zero-drift Maxwellian, `k_B T_e = 10 eV`, `n_e ∈ {0, 1e6, 1e7} m⁻³`, at rest in frozen LSO. Stationary, force-free ⁴He only. A source-bound build with Na (weight reference), He and He⁺, the zero-acceleration reference problem and individual packet weights. Common-time tracing. PIC keys: `helium_ion_characteristic_speed_m_s = 3e6` (the admitted classical speed limit), He⁺/Na weight ratio 0.125, He/Na weight ratio 0.5. |
| Inputs | [`reference/helium_case.json`](reference/helium_case.json); the approval record; the frozen admission files; the state profile (`42943c73…`); D05 `cfit.dat` (`470ea307…`); the Kim comparison and CCC diagnostic pins; and the source-bound producer manifest for the PIC child. |
| Independent oracle | Decimal (70-digit) Voronov coefficient and centred-difference slope (U16 reference). First hits `−ln(U)/(n_e k)` from independently reproduced HMAC-SHA256 source draws. Endpoint-survival lifetime against `1/(n_e k)`. Reservoir moments from Decimal `E_in = k_B T (3/2 + slope)`, `E_out = E_in − I`, and the profile masses and charge. Upper-end recoil bound from the approved formula. Analytical ballistic continuation (zero displacement) for the PIC product. |
| Units | SI: m, s, m/s, kg, J, C, m⁻³, m³/s, s⁻¹. The bath temperature is given as `k_B T_e` in eV, as D05 declares. |
| Frame / time | One frozen LSO event frame for bath, admission and planner. Particle epochs are SPICE ET (TDB seconds past J2000). Ensemble intervals are relative seconds starting at 0. |
| Deterministic / stochastic controls | Seed 731993. Ensemble cohort `I23-D05-4He-10eV-electron-bath-v1`, 32,768 parents, one optical depth, parent weight 3, He⁺ weight 2. PIC cohort `I23-PIC-helium-electron-bath-v1`, a single parent conditioned on its keyed source draw plus 1/16 s, representing 2.625 atoms. MPI layouts 1 and 8, `OMP_NUM_THREADS=1`. The densities 1e6 and 1e7 reuse the same keyed draws, so they test linear density scaling, not independent statistics. |
| Procedure | Authenticate the profile, D05 and pins. Regenerate the admission and require byte equality. Compile the probe against the production sources (ensemble) or link it against the authenticated producer (PIC). Run every layout and score independently. |
| Metrics | Coefficient, rate, slope, mean incident energy, physical threshold; refusal controls; first-hit classification and time error; endpoint lifetime; per-event He⁺ identity, velocity and weights; bath electron numbers, mass, charge, energy and momentum; ledger closure; recoil bound; MPI invariance; PIC epoch, step populations, clock, lineage, weight and dt invariance; the zero-density control. |
| Acceptance criteria | Coefficient ≤ 1e-12 relative (U16) and ≤ 2% (formal criterion). Ensemble lifetime ≤ 3% (formal criterion). Explicit refusal of zero/invalid density or temperature, drift, moving target, missing frame, missing or corrupt admission, and a central (non-conservative) tolerance. Arithmetic: event-time and moment relative errors ≤ 5e-12; PIC clock ≤ 16 binary64 spacings of the epoch; exact represented numbers and charge stoichiometry. |
| Expected PASS | All checks met. This is supporting evidence only. |
| Expected failure modes | Wrong rate, units or density scaling; frozen-at-injection hazard; event on a moving or accelerated parent; product with the parent's weight correction; missing or double bath terms; bath moment taken as a residual; non-closing ledger; MPI-dependent histories; a central tolerance admitted; startup ignoring the He/He⁺ keys. |
| Required data | D05 (READY under the U16 preflight); the NIST/CODATA state profile; the pinned Kim and CCC records. No observations. D06 is not used. |
| Data hashes | Recorded in each `run_manifest.json` (`input_sha256`, `source_sha256`, SPICE kernels for PIC). |
| Predecessor gates | U16 (local coefficient). Formal I23 still requires same-campaign I13–I20, U16 and the full frozen family. |
| Artifacts | `result.json`, `run_manifest.json`, `comparison.csv/png/eps` (ensemble); `refinement.csv`, `step_populations.csv/png/eps` and per-rank native CSV tables (PIC); stdout/stderr logs; the regenerated admission; `case.json`; `independent_reference.json`. |
| Status semantics | PASS: evaluated and met. FAIL: evaluated and not met. ERROR: build, parser, corrupt-pin or launch failure. SKIPPED: missing profile, D05, compiler or MPI. The owner keeps the stronger child status and records the formal result as SKIPPED. |

Run the children through the formal owner (from the AMPS root):

```bash
/usr/bin/python3 -B srcMoon/test/run_all_test.py --only I23 --skip-build \
  --campaign-mode regression --mpi-ranks 8 --test-runtime-seconds 1800 \
  --i23-helium-profile test_output/srcMoon/m3-followup-MfThQfny/atomic-family/profile.json \
  --i23-pic-producer "$producer/normal_build_evidence.json" --i23-pic-refinement \
  --output-dir "$(mktemp -d -p test_output/srcMoon i23-formal-XXXXXXXX)"
```
