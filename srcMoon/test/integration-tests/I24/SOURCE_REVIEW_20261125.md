# I24 first transport channel — quantitative source review (2026-11-25)

This review carries out decision 4 of
[the M3 completion follow-up](../../M3_COMPLETION_FOLLOWUP_20261124.md):
select the first charge-exchange transport case and its energy and angle
envelope from quantitative source evidence. The user asked for the review on
2026-11-25. It compares the sources for proton–sodium electron capture over
solar-wind energies. It changes no runtime selection, admits no numerical
rows and digitizes no figures. Formal I24 remains **SKIPPED**.

Energies are projectile laboratory energies per atomic mass unit with the
Na target at rest (keV/u). For this unit, `E = ½ u v²`; solar-wind speeds of
250–900 km/s correspond to 0.32–4.2 keV/u, and 1 keV/u is 438 km/s. Cross
sections are in Å² = 10⁻¹⁶ cm² = 10⁻²⁰ m².

## Status of the proposed review order

1. **IAEA D119990.** Closed by the r11 lineage review: the knots are the
   outer-shell capture sum of Knoop et al. 2004 Fig. 2, not an H(1s) partial.
   See [the r11 crosswalk](../../../../utility/chemical_database/NATIVE_M3_SODIUM_BENCHMARK_SOURCES.md).
2. **Courbin-Gaussorgues et al. 1990.** The retained PDF supplies
   state-to-state integral cross sections for capture into H(n=2) and H(n=3)
   at 0.35–5 keV/u, the H(2s) metastable fraction, and the angular extent of
   the computed differential cross sections. It has no H(1s) column.
3. **Schöffler et al. 2009 (proton–He).** Its 60–300 keV/u measurements lie
   far above solar-wind energies. Not pursued for the first Na channel.

## Finding 1 — the ground-product partial is a negligible share of capture

Shingal et al. 1986 Table 4 (A34 basis; column labels checked on the page
image) gives every H(nl) partial for H⁺ + Na(3s):

| E (keV/u) | Total (Å²) | H(1s) share | H(n=2) share | H(n=3) share |
| --- | --- | --- | --- | --- |
| 0.5 | 14.66 | 0.41% | 93.5% | 6.1% |
| 1.0 | 34.39 | 0.03% | 94.6% | 5.3% |
| 4.0 | 66.14 | 0.33% | 83.2% | 16.5% |
| 10.0 | 24.42 | 1.80% | 67.3% | 30.9% |

Knoop's MOTRIMS measurement agrees independently: at 4 keV/u, 85% of the
outer-shell recoil intensity is capture into H(n=2) (thesis Ch. 5 §5.2).
The current development case, H⁺ + Na(3s) → H(1s) + Na⁺ at 4 keV/u
(0.22 Å²), therefore carries about 1/300 of proton-capture Na⁺ production
at 4 keV/u and about 1/3400 at 1 keV/u. It remains a valid software and
kinematics fixture. It cannot represent lunar Na⁺ production.

## Finding 2 — total capture: evaluated data and theory

| Source | Kind and range | Comparison |
| --- | --- | --- |
| Morgan et al. 1985, Table 2, Na column (installed `ProtonNaMorgan1985LogLog`) | Evaluated recommended capture into all H states; knots 1–200 keV/u; no Na entry below 1 keV/u | Reference |
| Shingal et al. 1986, A34 total | Atomic-orbital theory, 0.5–20 keV/u | −20.0% at 1, −5.2% at 4, +1.8% at 10, −5.4% at 20 keV/u |
| Courbin et al. 1990, 19-state n=2 + n=3 | Molecular-orbital theory, 0.35–5 keV/u | −14% to −30% at 1–5 keV/u (no n=1 column) |

Morgan's source table (Table 1, p. 978) rates the Na total-capture inputs as
category (b), ±20–50%: Ebel 1983 (1–5 keV/u), DuBois 1984 (2–100 keV/u) and
the Fritsch 1984 atomic-orbital calculation (0.2–20 keV/u). The theory–evaluation
differences lie inside that rating. Below 0.5 keV/u Courbin reports that its
molecular-orbital H(2p) partial is twice the atomic-orbital value and that
Lyman-α experiments support the atomic-orbital result. The 0.35–0.5 keV/u
molecular values therefore carry an unresolved factor of two for H(2p).
Morgan at 4 keV/u uses the installed source-knot log-log law.

## Finding 3 — product states and energy release

With release `Q = IE(H, n) − IE(Na 3s)`, capture into n=1 releases about
8.46 eV, while capture into n=2 and n=3 absorbs about 1.74 and 3.63 eV
(Courbin R1 gives ΔE = −1.74 eV; Knoop tabulates the same magnitudes in the
opposite sign convention). Courbin Table 3 (19 states) gives the metastable
fraction H(2s)/H(n=2) as 0.245–0.44 over 0.35–5 keV/u. H(2p) decays to H(1s)
by Lyman-α emission and H(2s) is metastable. An implementation must take the
ionization energies and radiative lifetimes from a named NIST source. The
values above are only for scale.

## Finding 4 — angular evidence and Na⁺ recoil

* **Knoop, MOTRIMS, 4–25 keV/u (measured).** Outer-shell capture occurs at
  internuclear distances of 10–20 a.u. and scatters by less than 0.5 mrad
  (laboratory). Transverse-momentum peaks stay at or below 0.3 a.u. for every
  channel. The n=2 distribution peaks at 0.15 mrad at 4 keV/u. Fig. 5.10 is
  peak-normalized and plotted to 1 mrad, where the n=2 curve is a few percent
  of its peak. The figure was inspected, not digitized. An earlier
  Lyman-α coincidence experiment (Knoop's ref. [201], not yet identified)
  measured H(2p) differential cross sections at 1–5 keV/u for angles of
  0.9 mrad and above.
* **Courbin, theory, 0.5 and 2 keV/u.** The computed state-to-state
  differential cross sections extend to a center-of-mass angle of about
  0.4° at 0.5 keV/u and 0.1° at 2 keV/u (§4.3). Figs. 7–8 were not digitized.

Na⁺ recoil energy implied by each angle model, from
[`recoil.py`](../../../../test_output/srcMoon/i24-source-review-y16WTCTy/recoil.py):

| E (keV/u) | Isotropic mean / max (eV) | 0.5 mrad lab (eV) | Courbin extent (eV) | p⊥ = 1 a.u. (eV) | Longitudinal, n=2 (eV) |
| --- | --- | --- | --- | --- | --- |
| 0.35 | 28 / 57 | 3.9e-6 | — | 3.3e-4 | 1.2e-4 |
| 0.5 | 41 / 81 | 5.5e-6 | 9.9e-4 | 3.3e-4 | 8.8e-5 |
| 1.0 | 81 / 162 | 1.1e-5 | — | 3.3e-4 | 5.7e-5 |
| 2.0 | 162 / 324 | 2.2e-5 | 2.5e-4 | 3.3e-4 | 4.4e-5 |
| 4.0 | 324 / 648 | 4.4e-5 | — | 3.3e-4 | 4.2e-5 |
| 5.0 | 405 / 810 | 5.5e-5 | — | 3.3e-4 | 4.4e-5 |

Both the measurement and the theory place the Na⁺ recoil at or below about
10⁻³ eV (under 100 m/s), well below the thermal energy of the parent Na. The
approved isotropic law gives 28–810 eV and is contradicted by both sources.
No acquired source quantifies the scattering tail beyond the measured or
computed ranges.

## Leads not yet acquired

These are known only from abstracts or reference lists and were not used:
Wang, Finck, Roller-Lutz and Lutz, J. Phys. B 26 L61 (1993), 1 keV H(2p)
angle-differential alignment from Lyman-α–particle coincidences; Finck et al.,
Phys. Rev. A 38 6115 (1988), Lyman-α emission for Na(3s) and Na(3p) at
0.5–10 keV; Okasaka et al., H(2p) capture at 1–16 keV; Zapukhlyak et al.,
J. Phys. B 38 2353 (2005), the published form of Knoop's Ch. 5 data; and
Knoop's ref. [201].

## Choices that need user selection

1. **Channel:** total one-electron capture with an explicit H(n) state mix,
   H(n=2) capture only, or the existing H(1s) fixture.
2. **Energy envelope:** 1–5 keV/u (evaluated data), or an extension to
   0.5 or 0.35 keV/u that relies on theory alone.
3. **Angle law:** forward scattering with the recoil bound above, the
   isotropic sensitivity law, or acquisition of the tail data first.
4. **Product state and energy ledger:** prompt H(n) with state-specific Q,
   with or without the later radiative decay recorded as energy leaving the
   system.

## Decision — user selection, 2026-11-25

After reading this review the user selected:

1. **Channel.** Total one-electron capture H⁺ + Na(3s) → Na⁺ + H(n), with
   an explicit state mix over n = 1, 2, 3. The total rate is the installed
   Morgan et al. 1985 evaluated values. Shingal et al. 1986 A34 Table 4
   supplies the state split used for Q and the product state.
2. **Energy envelope.** 1–5 keV/u projectile energy (target at rest). Other
   energies are refused, never extrapolated.
3. **Angle law.** Zero-angle exchange in the pair center-of-mass frame with
   exact momentum and energy closure. The omitted recoil is bounded by the
   measured and computed angular extent (Finding 4). Na⁺ velocity
   observables are admitted only above that scale, and the unquantified
   scattering tail is a recorded limitation. The isotropic law remains a
   sensitivity option only.
4. **Product state and energy ledger.** Each event records the prompt H(n)
   with its state-specific Q. The later radiative decay to H(1s) is booked
   as photon energy leaving the system; the H trajectory is unchanged.
   Ionization energies and radiative data come from a named NIST source.

The existing H(1s) fixture stays as a software and kinematics test.

**Open before implementation.** The selection leaves these to be specified,
and to be approved where they are physics choices, in the implementation
contract:

* the rule for the state split between Shingal's 1.0 and 4.0 keV/u rows and
  above 4.0 keV/u, since Table 4 has no row at 5 keV/u;
* the numerical recoil bound as a function of energy inside 1–5 keV/u, given
  that the measured extent starts at 4 keV/u and the computed extent is
  given at 0.5 and 2 keV/u only;
* the H(n=3) cascade branching and the H(2s) share needed by the decay
  ledger, and whether the decay is booked at the event or timed;
* the NIST identities for the ionization energies, the H level energies and
  the radiative data.

**Resolved 2026-11-25.** The user then approved: log-log interpolation of
the shell partials between Shingal rows (no extrapolation outside
1–10 keV/u); a constant 0.4° COM bound for the omitted recoil; radiative
decay booked at the event; and NIST ASD H I levels with the statistical-weight
shell centroid, sampling only the shell n. See the
[implementation contract](TOTAL_CAPTURE_CONTRACT.md).

## Provenance and status

Evidence: `test_output/srcMoon/i24-source-review-y16WTCTy/` holds
`sources.json` (source and page-render SHA-256s), the recoil and state-fraction
scripts with their outputs, and the inspected page images. Source SHA-256s:
Shingal 1986 `58bc33f7…`, Shingal and Bransden 1987 `f63e7468…`, Courbin 1990
`030e49b5…`, Knoop Ch. 5 `7daca10c…`, Morgan 1985 `7c412cb9…`. The source
roadmap is `712e1152…` and the data roadmap `d10eca1c…`. Numerical rows admitted:
0. Figures digitized: 0. Runtime selection changed: no. Formal I24: SKIPPED.
