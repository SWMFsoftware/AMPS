# I23 helium electron-bath admission — user approval of 2026-11-25

This record documents the decisions approved by the user on 2026-11-25 for
coupling the prescribed electron bath to the normal Moon reaction controller.
The machine-readable authority is
[`reference/helium_admission_approval_20261125.json`](reference/helium_admission_approval_20261125.json).
The decisions resolve decision 3 of
[`M3_COMPLETION_FOLLOWUP_20261124.md`](../../M3_COMPLETION_FOLLOWUP_20261124.md)
for this coupling scope only. They grant no formal I23 PASS, observational
claim or M3 closure.

## Decisions approved by the user

Each question was presented with the derived numbers and the alternatives.
The selected answers follow.

| Question | Approved answer |
| --- | --- |
| Controlled bath | **Approve as proposed.** Isotropic zero-drift Maxwellian, `k_B T_e = 10 eV`, `n_e ∈ {0, 1e6, 1e7} m⁻³`, at rest in the frozen LSO frame, spatially uniform and constant (no wake/shadow structure, no solar-wind drift), acting only on stationary force-free ⁴He. It is a controlled verification bath, not lunar plasma. |
| Event products | **Production reservoir, narrow scope.** Use the real `HeliumElectronImpact::Compose` inherited-velocity reservoir. He⁺ inherits the parent velocity; the bath records one incoming and two outgoing electrons. Formal I23 stays SKIPPED. |
| Threshold/state and conditional mean | **Approve both.** (a) D05 fit threshold 24.6 eV against the authenticated NIST He ionization energy, all ionizations assigned to He⁺(1s); the excited-He⁺ yield is unknown, and CCC's ≤9.81% above the He⁺(n=2) onset is recorded as energy eligibility only. (b) Bound the `d ln k / d ln T` error by the published Kim (1992) minus Voronov spread, 0.120, which is a model spread and not a confidence interval. |
| Recoil, observable and He mechanical limits | **Approve both.** The observable scope is heavy-particle number, charge and mass budgets and event chronology only, with no He⁺ velocity, energy or pickup observable. He limits equal the approved Na P02 numbers: 1e-8 inertia, 1e-8 force, 1e-4 nonrelativistic, 3e6 m/s heavy-particle speed. |
| Recoil tolerance actually evaluated | **Conservative upper-end values.** The composer evaluates the omitted-recoil bound at slope + 0.120. The admitted tolerances are ≤1095 m/s RMS He⁺ speed and ≤0.0249 eV kinetic, replacing the central 1.073 km/s and 0.0239 eV first shown. |
| He/He⁺ startup | **Explicit `.in` keys,** following the tracked-H precedent, with no physical defaults. See [INPUT_FILE.md](../../../../src/INPUT_FILE.md#moon-helium-time-step-and-reference-weights). |

The approval question quoted the NIST He ionization energy as 24.587387 eV.
The authenticated profile record (`ionization energy/hc = 1.9831066637(20) ×
10⁷ m⁻¹`) gives 24.5873890 eV. The 2 × 10⁻⁶ eV quotation discrepancy does not
change the approved rule; the generator uses only the authenticated value.

## Derived evidence

[`prepare_helium_admission.py`](prepare_helium_admission.py) derives two files
from pinned sources and the approval record. Each I23 run regenerates them and
requires byte equality with the frozen copies pinned in
[`reference/helium_case.json`](reference/helium_case.json).

| File | Content | Value |
| --- | --- | --- |
| [`helium_state_crosswalk.json`](reference/helium_state_crosswalk.json) | Source threshold difference | 0.012611 eV |
| | Evaluation-rounding allowance | One binary64 spacing of the stored 24.6 eV fit threshold times e, plus 64 extended ulps; about 6 × 10⁻³⁴ J |
| [`helium_model_error.json`](reference/helium_model_error.json) | Admitted `|Δ d ln k/d ln T|` | 0.12003003713840445 (Kim − Voronov) |
| | Lower-end conditional mean | 44.22 eV, above the 24.587 eV threshold |
| | Central / upper RMS He⁺ speed bound | 1072.6 / 1094.7 m/s |
| | Upper kinetic bound, stationary parent | 0.024855 eV |

Sources: D05 `cfit.dat` SHA-256 `470ea307…`; the state profile
`42943c73…`; the retained Kim comparison `fc1cc2d9…` (produced by
`run_helium_pilot.py --kim-1992-pdf`); and the CCC diagnostic `41ab6b0e…`, an
unpublished table used for a recorded fraction only and never in a computation.
The generator refuses an approval whose tolerance falls below the composer's
upper-end bound, and it never widens a tolerance.

## Scope limits

- The reduced reservoir fixes the net bath change: +1 electron, +mₑ, −e and
  −I per event for a stationary parent, with zero momentum. The gross
  incoming/outgoing electron energies depend on the Voronov model mean and
  carry the admitted slope uncertainty.
- Omitted He⁺ recoil (up to about 1.1 km/s RMS) is comparable to thermal He
  speeds. Any velocity- or energy-sensitive He⁺ observable needs a separate
  reviewed admission.
- Moving targets, drifting or structured baths, other temperatures, measured
  D06 forcing, Ne/Ar/Na electron impact and excited-He⁺ yields are not
  admitted.
