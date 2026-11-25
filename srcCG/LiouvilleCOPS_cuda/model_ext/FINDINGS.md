# FINDINGS.md — physics results of the 12M COPS fits, for use in publications

Kept current (user request 2026-11-12: "create FINDINGS.md and save that there. also add there new important physics
related notes in the future for using that in the publications; keep the file current"). Each entry: date, the
finding with its numbers, the evidence (job ids, files), and caveats. Newest first within each section. Paths are
relative to `COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output` unless absolute. Figures: see
`OUTPUT_VISUALIZATION.md`.

---

## A. Findings

### A25. Plausibility ranking of all new-law fits: the all-held candidate first, A second, D third; every other fit switches off most of the seen surface (2026-11-23)

* **Request** (user 2026-11-23): "update the analysys physical plasuability of each of the fit and range it". 17 fits, one
  method (job 138119, q002_plausibility/v2_138119/; the ten agree with the map-based numbers to 4e-16).
* **Rule (proposal):** tiers by the user's criterion 1b and the energy limit (<= 1 % of Q above F_E; the 1 % is mine),
  then the data cost. Flags only: convergence, old-law templates, unseen share of Q, Q over the day, uniqueness.
* **Ranking:** 1 all-held candidate (7.62e-3; 35.5 % off; max F 0.004 F_E; 95 % of Q seen; one solution from 7 starts);
  2 A (9.64e-3; 1.9 % off; 0.35 % of Q above F_E on 4 template facets; 80 % of Q unseen; not converged); 3 D (8.28e-3;
  14.7 % off, 50.2 % at 10 T; 24 % of Q above F_E); then the fits that fail 1b: within the energy limit the templates-off
  solutions 70-80, 60-70, 50-60 deg, D2, the 40-50 deg solution; over it C P_m 1e-4 (f10, f1), C f10, the 80-90 deg
  solution, C f1, softplus, B. Cliff-test fits listed apart (not candidate models).
* **Meaning.** Keeping the seen surface active costs fit: the candidate's data cost is 3 x the best of the round and its
  late NG misfit (0.209) the worst of the ranked fits. With a strict energy limit (no facet above F_E) it is the only
  fit in the first tier. No significance (no noise model); only the candidate and the templates-off model were tested
  for other solutions.
* **Evidence.** figures_cluster/12M/q002_plausibility/README.md (UPDATE 2026-11-23), v2_138119/{plausibility_v2.json,
  ranking_table.md, ranking.json}.

### A24. Cliff test: switching cliffs off costs 4.2 % in data cost, less than switching off the moderate slopes (6.7 %): UNDECIDED by the rules fixed beforehand (2026-11-22)

* **Fits** (user "submit cliff test"; plan PLAN_cliff.md, rules written 14:03 CST before any evaluation): the all-held
  option C of A23 (f 10, templates off, f_m 0 as all multipliers are held) with the cliff classes (60-180 deg) held at
  0.01, and as the control the moderate slopes (40-60 deg) held at 0.01; two starts each (the all-held baseline, D2).
  Jobs 138062-138074, none failed; both starts of each set end at the same state (objective equal to 1e-10), each a local
  minimum by the second-order test.

  | | data cost | change | objective | rms ln NG / RG | late NG / RG | seen off (1b) | Q [1/s], seen share | cliffs: share of seen production, NG / RG signal |
  |---|---|---|---|---|---|---|---|---|
  | baseline (all at 1) | 7.6201e-3 | - | 8.2409e-3 | 0.0998 / 0.0726 | 0.2087 / 0.0607 | 35.5 % | 5.66e24, 94.5 % | 13.3 %, 10.0 / 8.7 % |
  | cliffs off | 7.9373e-3 | +4.16 % | 8.5785e-3 | 0.1018 / 0.0743 | 0.2085 / 0.0599 | 36.1 % | 5.25e24, 96.4 % | 0.16 %, 0.11 / 0.09 % |
  | moderate slopes off (control) | 8.1280e-3 | +6.67 % | 8.8509e-3 | 0.1035 / 0.0745 | 0.2116 / 0.0564 | 37.1 % | 5.48e24, 95.4 % | 16.1 %, 12.1 / 10.5 % |

* **By the rules: valid (the held cliffs give 0.16 % < 1 %), UNDECIDED (1 % < 4.16 % < 5 %).** The data cost rises when
  cliffs are switched off, but by less than when about as much moderate-slope terrain is (40-60 deg gives 16.5 % of the
  baseline's seen production, the cliffs 13.3 %): nothing in this test shows that cliffs matter more than other terrain.
  The late misfit hardly changes (late NG 0.2087 to 0.2085): cliffs do not explain the late NG rise in this model.
* **Limits.** One model (no templates, no slope preference, f 10); no noise model, so the 1 % and 5 % are conventions.
  The flat terrain of the same 150-m patches takes over part of the cliffs' production (Q -7 %).
* **Evidence.** runs/shap8_012m_lategrid_dz1.5/nonneg/localnn/clifftest_q002_f10_tploff_20261122T233834Z/ (report_138074/, stationarity_138074/).

### A23. Without slope multipliers option C has one solution from every start tried, also through the chain: a stable candidate that keeps the seen surface active (follow-up c, 2026-11-22)

* **Runs** (user "submit c"; plan PLAN_ac.md): option C f 10, templates off, all multipliers held at 1, SHAP8 12M, new law:
  six starts (activity maps of the five solutions of A22 with m = 1, and a uniform map; start objectives 0.07 ... 0.49)
  and the existing fit from D2; then five cycles of fit + fine tuning. Jobs 138027-138047, none failed.
* **All seven end at the same solution:** objective 8.2408685e-3 (spread 3e-11 relative), data cost 7.6201e-3 (5e-8);
  coefficients c_j equal to 1e-7 (median) and <= 0.14 % (every c_j above 1e-3 of the largest); ftol or gtol after 30-68
  evaluations; each a local minimum by the second-order test. The six new starts shared kappa, grain and seen production;
  the D2 start did not. Seven starts are no proof that the minimum is global.
* **Chain:** unchanged through five cycles (max |h| 4.8e-4; fits 2-5: 1-3 evaluations).
* **The state:** 35.5 % of the seen area off (passes 1b; 33.1 ... 43.1 % from T/10 to 10 T), Q 5.657e24 1/s at 10:33:10
  (94.5 % seen), max F / F_E 0.004, rms ln NG / RG 0.0998 / 0.0726, late 0.2087 / 0.0607. Cliffs (>= 60 deg, 14.6 % of
  the seen area) give 13.3 % of the seen production and 10.0 / 8.7 % of the modelled NG / RG signal (median row).
* **By the plan's rules: stable, and a candidate model.** Limits: late NG is its weak point; no templates; smoothness f 10
  and the temperature law held, other settings not varied; no slope dependence by construction.
* **Evidence.** runs/shap8_012m_lategrid_dz1.5/nonneg/localnn/slopetest_c_q002_f10_pm0.01_tploff_20261122T195607Z/ (report_138047/, stationarity_138047/, chain/).
  Full figure set (2026-11-23): figures_cluster/12M/full_set_q002_C_f10_allheld/. Its facts: before the late window NG 93 % from
  the big lobe; late window small lobe NG 41 %, RG 67 %; P50 source-area union NG 6.2e5 m^2, RG 4.4e5 m^2; >= 60 deg 9 % of
  the late NG, 6 % of the late RG; overhangs 0 %.

### A22. Option C without templates has five solutions, one per slope class from 40 to 90 deg: the fits do not determine the slope class; without slope multipliers it keeps most of the seen surface active (3c + multi-start, final 2026-11-22)

* **Question and fits** (user decisions 3c with L = 0.01, 2026-11-21, and "multi-start when ready", 2026-11-22). Option C,
  f 10, P_m 0.01, templates off, new law; reference = the cycle-1 fit of the round (start D2). Eight control fits (all six
  multipliers held at 1; one class free at a time; all free >= 0.01) and four all-free (>= 1) refits from their end
  states that were better than the reference (two candidates skipped by the approved rule).
* **Five solutions of the same problem (six multipliers free, >= 1):**

  | class with the seen production | objective | data cost | m | share of seen production | rms ln NG / RG | late NG / RG | seen off (1b) | Q [1/s], seen share | max F / F_E (share of Q above) | minimum test |
  |---|---|---|---|---|---|---|---|---|---|---|
  | 40-50 deg (reference) | 6.9499e-3 | 5.6734e-3 | x102 | 92.1 % | 0.0821 / 0.0679 | 0.1879 / 0.1191 | 65.8 % fail | 7.41e24, 90.9 % | 0.172 (0 %) | local minimum |
  | 50-60 deg | 6.3746e-3 | 4.7846e-3 | x226 | 94.8 % | 0.0753 / 0.0625 | 0.1727 / 0.0952 | 69.4 % fail | 9.00e24, 88.8 % | 0.265 (0 %) | local minimum |
  | 60-70 deg | 6.1586e-3 | 4.5172e-3 | x242 | 94.2 % | 0.0722 / 0.0618 | 0.1606 / 0.0913 | 62.9 % fail | 9.05e24, 93.4 % | 0.337 (0 %) | local minimum |
  | 70-80 deg (lowest objective) | 5.9560e-3 | 3.9123e-3 | x614 | 95.9 % | 0.0674 / 0.0573 | 0.1528 / 0.0734 | 67.6 % fail | 1.138e25, 88.0 % | 0.748 (0 %) | local minimum |
  | 80-90 deg (lowest data cost) | 6.5102e-3 | 3.7569e-3 | x2166 | 97.9 % | 0.0681 / 0.0537 | 0.1630 / 0.0478 | 72.3 % fail | 1.707e25, 76.1 % | 2.78 (23.2 %) | nearly: 4e-8 of the objective left, curvature positive |

  Minimum test (job 137980): no descent on the bounded Gauss-Newton step (<= 1e-10 of the objective) and a positive
  smallest eigenvalue of the exact Hessian on the free set (checked against differences to 2e-10). All five fail 1b at
  every threshold from T/10 to 10 T (59.6 ... 82.0 % off). The seen production spans 6.7e24 ... 1.30e25 (x1.9), the total
  x2.3. The optimizer's ftol says nothing about which of the five it found. Other minima may exist.
* **What the multipliers buy, against all held at 1** (objective 8.2409e-3, data cost 7.6201e-3): data cost -25.5 %
  (40-50 deg), -37.2 % (50-60), -40.7 % (60-70), -48.7 % (70-80), -50.7 % (80-90, above the energy limit), -0.05 %
  (90-180; m 1.20; one start). Any of five classes can take the role, so the class is not measured. Which one has the
  lowest objective depends on the penalty (the lowest data cost is another one). No significance: no noise model, the
  shuffled-label test (3d) was not run.
* **All held at 1: passes 1b.** 35.5 % off (33.1 ... 43.1 % from T/10 to 10 T), a local minimum of its own problem, data
  cost below A (9.6401e-3) and D (8.2770e-3), Q 5.66e24 (94.5 % seen), max F / F_E 0.004; its late NG fit is worse
  (0.2087), its late RG fit better (0.0607). One start, no fine tuning, no templates [A23: seven starts and the chain give the same].
* **Lower bound 0.01: INCONCLUSIVE.** free_min stopped at max_nfev in the 50-60 deg solution (not a minimum; 40-50 deg at
  0.954) at objective 6.3823e-3, above three of the five solutions, which are feasible for it: it did not find its
  minimum. The bound 1 is active: in every solution with 40-50 deg at the bound, and in the two with 50-60 deg at the
  bound (60-70, 70-80), the gradient pushes it below 1.
* **How far (preview, job 138026: one bounded Gauss-Newton step with the bound 0.01 at each solution, no fit).** 40-50 deg
  solution: no change (already a minimum with the lower bound). 50-60: m(40-50) 0.95. 60-70: m(40-50) 0.85, m(50-60) 0.99.
  70-80: m(40-50) 0.79, m(50-60) 0.96. 80-90: m(40-50) 0.91. Decrease available: <= 2.8e-4 of the objective, data cost
  <= 0.1 %; the promised decrease is available at 0.75-1 of the step, so this previews the minima nearby. No class at
  60 deg or steeper goes below 1 (those not dominant sit at 1.002-1.03). So a bound below 1 changes little IN THESE FIVE
  solutions and says nothing about inactive cliffs: there almost all production comes from one class and the others
  hardly affect the data. That question needs a fit with the cliffs held inactive.
* **Fine-tuning leverage on the all-held fit (job 138026):** the 271 late facets give 2.2e-4 (mean; max 2.9e-3) of the late
  NG surface pressure, so a chain is not expected to change its late NG misfit.
* **Old law against new law, templates off: not a law effect.** The old-law fit (data cost 3.6519e-3, +55.4 % to the new
  one) sat in the 80-90 deg solution (x1300), the new one in the 40-50 deg one. In the same class the change is +2.9 %
  (3.7569e-3). The same may hold for f10 (70-80 to 50-60 deg, +12.1 %) and f1 (70-80 to 40-50 deg, +23.5 %) with
  templates and for the 'multiplier shift' of A18: not tested.
* **Not converged:** only_40_50, only_50_60, only_80_90 and free_min (max_nfev 4000; 6e-5 ... 2e-3 of the objective still
  available on one step); only_70_80 has 9e-7 left.
* **What was wrong in the interim reports (2026-11-21/22):** the reference was taken as the model's solution; "the six
  multipliers lower the data cost by 25.5 %" and "they raise the switched-off share from 35.5 to 65.8 %" described that
  solution; "nothing suggests less active cliffs" rested on the unfinished free_min fit. Not a code error: the
  reference's frozen code and direct facet sums reproduce the better states (jobs 137918, 137980; 2e-16).
* **Scope.** The other option C fits (templates on, P_m 1e-4), A, B and D were not tested for other minima (one start each).
* **Evidence.** runs/shap8_012m_lategrid_dz1.5/nonneg/localnn/slopetest_q002_f10_pm0.01_tploff_20261122T044019Z/
  (report_137981/, state_check_only_80_90_137918/) and multistart_q002_f10_pm0.01_tploff_20261122T065041Z/
  (stationarity_137980/).

### A21. Switched-off surface by decision 1b: only A and D pass; the fits did not test whether cliffs are less active (correction of A20, A19 and A18's wording; 2026-11-21)

* **Decision (user, 2026-11-21: "as recommended: 1b, 2b, 3c").** A facet is active if its activity, not counting an
  imposed floor, exceeds T = 2.114e15 = 1 % of the median srcCG AH2O of the seen facets (one value per facet; the threshold
  of localnn_report.py and nonneg_alt_summary.py since 2026-11-14). More than 50 % of the seen area switched off fails. The
  only imposed floor is B's m_b A_bg A_orig exp(h) (A_bg = 0.01 = its lower bound). D's spread penalty (no lower bound)
  and the held templates (0.001 % of the seen area) are not floors.
* **Seen area switched off:**

  | fit | data cost [1e-3] | off (1b) | floor counted | at T/10 | at 10 T | result |
  |---|---|---|---|---|---|---|
  | A | 9.640 | 1.9 % | 1.9 % | 0.1 % | 20.8 % | pass |
  | D | 8.277 | 14.7 % | 14.7 % | 3.3 % | 50.2 % | pass (on the line at 10 T) |
  | B | 5.646 | 96.2 % | 44.1 % | 96.2 % | 96.2 % | fail |
  | C f10/0.01, templates off | 5.673 | 65.8 % | 65.8 % | 63.7 % | 74.4 % | fail |
  | C f10/0.01 | 3.586 | 80.0 % | 80.0 % | 77.3 % | 87.3 % | fail |
  | C f1/0.01 | 4.173 | 81.1 % | 81.1 % | 80.2 % | 84.6 % | fail |
  | C f1/1e-4 | 2.500 | 92.3 % | 92.3 % | 89.1 % | 96.9 % | fail |
  | C f10/1e-4 | 2.499 | 93.6 % | 93.6 % | 89.0 % | 98.1 % | fail |
  | D2 | 5.671 | 93.6 % | 93.6 % | 93.6 % | 93.6 % | fail |
  | softplus | 5.216 | 96.5 % | 96.5 % | 96.5 % | 96.5 % | fail |

* **Consequence.** The two fits that pass have the two highest data costs. D has 24 % of Q above the energy limit (98
  facets); A exceeds it on the 4 template facets only (0.35 % of Q), is not converged and has 80 % of Q from unseen
  terrain. B's floor is the only activity above T on 96 % of the seen area and gives about 3 % of the seen production.
* **What was wrong.** (1) A19 counted exact zeros, although the C reports already gave the active seen area with this
  threshold (0.20 / 0.19 / 0.34 for f10, f1, templates off). (2) A20's A < k A_u was ad hoc; with it B seemed to pass at
  k = 1e-6. (3) A20's "60-80 deg = 1.00-1.07 in all five" is the imposed bound m_b >= 1 (user 2026-11-12) or within 7 %
  of it. Neither "active cliffs are not required by the data" (A19) nor "cliffs are as active as flat terrain" follows:
  no fit of this round could make a slope class less active than 0-40 deg terrain nearby. (4) A18's "NG overhang
  dominance at touchdown ... not required by the data" goes too far: the fit without templates has none, and it does
  not reproduce the NG rise of the last ~20 m.
* **Evidence.** figures_cluster/12M/q002_decisions_20261121/active_area.{json,csv} (job 137886; B rebuilt from its terms
  to 5e-16; the C and D2 shares equal the reports'); threshold_options.json (job 137884).
* **Open (decision 3c).** Control refits of C f10 / 0.01 with the templates off: all multipliers at 1; one class free at a
  time; all free with the lower bound 0.01 (user 2026-11-21). Done 2026-11-22: A22.

### A20. [THRESHOLD AND CLIFF STATEMENT SUPERSEDED, see A21] CORRECTION of A19, and activity against the angle to gravity: the C fits switch off 63-80 % of the seen surface; the slope effect in the fits equals the fitted multipliers (2026-11-21)

* **Correction of A19.** "Area with A = 0" does not measure switched-off surface: the bounded TRF leaves option C's
  switched-off coefficients slightly positive. Seen area with A < 1e-6 A_u (< 0.01 A_u; A_u = the uniform activity with
  the fit's own seen production at 10:33:10): A 0 (3) %, D 1 (29) %, B 1 (58) %, C f10/0.01 off 63 (67) %, C f10/0.01
  76 (82) %, C f1/0.01 80 (82) %, C P_m 1e-4 86-87 (96-97) %, D2 94 %, softplus 96 %. **The ranking of A19 is
  withdrawn:** its tier 1 fails the user's criterion; only A passes at every threshold tested [A21: A and D pass, B fails] (worst data cost, not
  converged). "Active cliffs are not required by the data" (A19 addendum) is withdrawn too. The "A = 0" figures of
  option C in q002_all_fits_summary and in the full sets' facts.json have the same limitation.
* **Visibility** (pass over all 4.74e9 kernel records, job 137865). Seen by >= 1 gauge 66.457 % of the area, by NG
  66.419 %, by RG 66.457 %, by both in the same observation 66.419 %; by both only at different times: none. 99.4 % of
  the seen area is seen sunlit at least once. The "5,924 NG-only / 29,210 RG-only" late facets of the visibility figures
  belong to an earlier base and do not hold for base_q002 (late window: 2 / 7 facets).
* **Angle to gravity, illumination removed** (A = F/q; each slope class against 0-40 deg in the same 150-m patch and
  diurnal-exposure quintile; 275 template / fine-tuned facets out). Matched ratio = fitted multiplier to ~2 %, A/m_b =
  0.99-1.02 in every class: C f10 x454 at 50-60, C f1 x99 and C off x102 at 40-50, A x2.34 at 40-50, B x6.7 at 50-60,
  x8.2 at 80-90, x13 overhangs. 60-80 deg = 1.00-1.07 in all five [at or near the imposed bound 1, not a measurement: A21]. F/cos SZA gives the same when matched; unmatched it
  departs by up to x1.7 (samples lit facets only, 38 % of the seen area). Same for all visibility classes.
* **Consequence.** These fits cannot give independent evidence of a slope effect: their only freedom within 150 m is
  the multiplier. Class and factor differ between fits. A test needs refits against controls (m = 1; one free class;
  shuffled slope labels): not run, the user's decision.
* **Template facets.** 275 facets (367 m2) hold 11-39 % of the seen activity of the fits with templates and make the
  overhang means look x10-270 the flat terrain (x1-8 without them).
* **Evidence.** figures_cluster/12M/q002_top5_models/slope_activity/ (jobs 137865, 137869, 137870, 137877; matched/,
  pass2/, maps_fixed/, shutdown_check/).

### A19. [RANKING AND "SURFACE SHUT DOWN" WITHDRAWN, see A20] Physical plausibility of the ten new-law fits: every fit but one exceeds the energy limit somewhere; the best-fitting fits (C, P_m 1e-4) put 83-90 % of Q on impossible fluxes, D2 and softplus switch off >90 % of the surface (2026-11-21)

* **Question** (user 2026-11-21): rank the fits by physical plausibility, with consistency with the observations and
  "not shutting down most of the surface" as criteria, plus other relevant criteria.
* **Energy limit (new criterion).** No facet can sublimate more water than the absorbed sunlight pays for:
  F_E = S0 / (r_h^2 dH) = 1.09e21 m^-2 s^-1 at r_h = 3.833 AU (S0 = 1361 W m^-2, dH = 51.0 kJ/mol, zero albedo, normal
  incidence, no thermal emission: an upper bound; free sublimation of exposed ice at the subsolar point is 3.7e20,
  0.34 F_E, with albedo 0.05, emissivity 0.9 and the Murphy & Koop vapour pressure). From the flux maps at 10:33:10:
  max F_H2O / F_E = 30 (D2), 1230-1590 (softplus, A, B, all C with templates), 5790 (D), 0.17 (C f10 / P_m 0.01 with
  the templates off, the only fit within the limit everywhere). Share of Q_H2O above F_E: 0 (templates off), 0.12 % (D2),
  0.35 % (A), 1.1-1.6 % (B, C P_m 0.01), 3.6 % (softplus; 332 m^2), 24 % (D; 511 m^2, overhang multiplier 87),
  83 % / 90 % (C f10 / f1, P_m 1e-4; 1.2e5 / 1.6e5 m^2, steep-slope multipliers 1e4-1e5).
* **Template spikes.** The largest activity of every fit with templates sits on the same two facets (11213754,
  12232930): A = 8e23-1.5e26, 10^5-10^7 x the mean; flux above F_E on 4 facets (1.1 m^2). With the templates off these
  facets have A = 2e16 and 10. The physical activity templates therefore carry fluxes ~1000 x the energy limit.
* **Surface shut down.** Area with A = 0: softplus 94 % (softplus underflow), D2 91 %, B and D 0.7 %, others 0. Area with
  A >= 10 % of the trimmed mean (mean without the most active 1 % of the area), whole surface: A 41 %, C templates off
  37 %, C f1 / 0.01 28 %, C f10 / 0.01 23 %, D 20 %, D2 8 %, B 7 %, softplus 5.5 %, C P_m 1e-4 2-4 %; B has 84 % of the
  surface below 1 % of that level.
* **Other criteria.** Seen share of Q: C P_m 0.01 88-92 %, D 41 %, D2 28 %, B 25 %, A 20 %, softplus 19 %, C P_m 1e-4
  3-4 % (Q 4.8-6.7e26 s^-1). Q max/min over the 13 map times: C P_m 0.01 1.9-2.3, A 2.8, D 3.1, C P_m 1e-4 6.5-7.1, D2
  7.9, B 11.8, softplus 15. Kappa is zero in all ten (not discriminating).
* **Proposed ranking** (figures_cluster/12M/q002_plausibility/README.md; criteria and weights are the user's decision):
  tier 1 C f10 / P_m 0.01 (data cost 3.59e-3), C f1 / P_m 0.01 (4.17e-3), C f10 / P_m 0.01 templates off (5.67e-3; the
  only energy-clean fit, but its chain did not move over five cycles); tier 2 A (widest activity, worst fit 9.64e-3, not
  converged), B (5.65e-3, activity concentrated, Q x11.8); tier 3 (fail a hard test) D2, softplus (surface off), D, and
  the two C P_m 1e-4 fits (best data cost 2.50e-3, but 83-90 % of Q above the energy limit).
* **Evidence.** figures_cluster/12M/q002_plausibility/: metrics.json / metrics.csv (job 137799, tool
  figures_cluster/scripts/plausibility_q002_metrics.py), spike_check (job 137800); fit numbers from
  model_vs_observations_nonneg_options_stage2_q002/summary.csv; seen shares and Q over the day from
  full_set_q002_<run>/checks/facts.json.
* **Caveats.** F_E is a hard bound for steady sublimation driven by the instantaneous sunlight; stored heat or other
  energy sources are not considered. The trimmed-mean thresholds are a choice (1 % and 10 %). A, B, D2 and C f1 / 1e-4
  stopped at max_nfev. No independent Q_H2O for 2016-09-30 was used.
* **Addendum 2026-11-21 (Word document on the five most plausible fits, figures_cluster/12M/q002_top5_models/).**
  Why the C f10 / P_m 0.01 chain with the templates off did not move after cycle 1 (data cost 5.6734e-3 in every cycle;
  fine tuning changed it by 1.4e-9, max |delta| 1.05e-3; fits 2-5 stopped after 3 evaluations): without the templates
  the 271 late-only facets carry only the smooth 150-m field, 2.8e-5 (NG) / 6.8e-5 (RG) of the mean late-window modelled
  pressure (templates on: 2.3-2.6 % / 3.9-4.2 %; computed for the document with the frozen localnn_problem /
  localnn_late functions, not in a run file). Even the cap (x100) on all of them would cost P_delta ~ 10 x 4.27e-3 x
  4.6^2 ~ 0.9, about 160 x the data cost, so the optimum delta is ~0. Its final NG does not follow the observed rise in
  the last ~20 m. In the C fits the multiplier penalty is the largest penalty (P_m / data cost 0.45, 0.22, 0.16 for
  f10 on, f1 on, off; c5/fit/final/fit_result.txt); each C fit puts one slope bin far above 1 (x456 at 50-60 deg,
  x99 and x102 at 40-50 deg). Free fit parameters: A 81, B 35, C 2754 (2745 nodes), D2 29.
* **Addendum 2026-11-21 (sources of the five most plausible fits, q002_top5_models section 4, numbers_sources.json from the
  full sets' facts.json).** Robust across C f10 / 0.01 on, C f1 / 0.01 on, C f10 / 0.01 off, A, B: seen production at
  10:33:10 6.22e24-1.05e25 s^-1 (factor 1.7) while the totals span 7.4e24-4.3e25 (x5.8); NG before the late window
  ~91-93 % from the big lobe; late window small lobe 26-44 % of NG, 58-69 % of RG, 30-100 m from touchdown 10-18 % /
  20-32 %; production per unit seen area enhanced on 40-60 deg slopes in every fit (x2.1 A at 40-50, x5.1 B at 50-60,
  x13 C f10 at 50-60, x10 C f1 and off at 40-50) but the bin is model-dependent; cliffs >= 60 deg (14.6 % of the seen
  area) give 1.4-3.4 % of the seen production and 0.2-2.7 % of the modelled signal in the C fits, 13 % / 4-9 % in A,
  28 % / 14-22 % in B (per unit area below average except B x1.9): active cliffs are not required by the data [WITHDRAWN, A20 and A21: not tested]. In the
  C fits 81-96 % of the pre-gap and late modelled signal comes from the one slope bin with the large multiplier. Dark
  (night + shadow) share of Q at 10:33:10: C 9.5-10.7 %, A 51 %, B 68 % (69.7 % of the area is dark). P50 source-area
  union (NG, late): C 1.0-2.9e4 m^2, B 3.9e4, A 7.1e5. Caveat: model results with held old-law templates.

### A18. CORRECTION: every result so far used H2O q_min = 0.2 q_max; with the intended 0.02 (srcCG) and the activity held fixed every final state misfits, data cost x3.3 to x53 (2026-11-18)

* **Correction notice.** The H2O night-side / cast-shadow production of EVERY fit, map, figure and finding so far
  (A1-A17 here, handoff.md, STATUS.md, all full sets, contributions, cliff and illumination products) was
  q_min = 0.2 q_max (LiouvilleCOPS Source.cpp `10.0*0.02`, the same 0.2 in fit_cops.py, cops_cuda/evaluator.py,
  build_physical_activity_templates.py and model_ext). srcCG (Comet::GetTotalProductionRateBjornNASTRAN) has 0.02. User
  decision 2026-11-18: "(b) -- use 0.02; x10 was not intended; use the same law but with decreased activity". The code
  now uses 0.02. The CO2 law (q_min = 0.1 q_max, r_h^-2) also differs from srcCG; it was not changed (not asked).
  Every number of A1-A17 holds for the 0.2 law only. Most affected: A1 (dark overhangs at the 0.2 floor), A15 / A17
  (Q_H2O and its dark / unseen shares; the x5 illumination factor becomes x50), A16 (the activities were fitted with 0.2).
* **Diagnostic, no refit** (lead 2026-11-18): the eight final states re-evaluated with A, T, kappa, grain and the
  fine-tuning factors held fixed, old law first (reproduces every recorded data cost to its printed digits; per-row model
  vs the states' tables <= 4.9e-10), then the new law:

  | state | data cost 0.2 -> 0.02 | rms ln NG | rms ln RG | late rms ln RG | Q_H2O(10:33:10) [s^-1] | seen Q | dark share of Q |
  |---|---|---|---|---|---|---|---|
  | softplus reg | 4.59e-3 -> 6.18e-2 (x13) | 0.074 -> 0.278 | 0.061 -> 0.215 | 0.042 -> 0.152 | 2.56e26 -> 3.54e25 | 1.48e25 -> 8.57e24 | 95 -> 69 % |
  | A reg | 1.07e-2 -> 3.55e-2 (x3.3) | 0.110 -> 0.243 | 0.097 -> 0.110 | 0.120 -> 0.123 | 1.84e26 -> 2.81e25 | 1.17e25 -> 5.11e24 | 93 -> 61 % |
  | B reg | 5.12e-3 -> 6.04e-2 (x12) | 0.078 -> 0.279 | 0.065 -> 0.207 | 0.052 -> 0.124 | 2.53e26 -> 3.51e25 | 1.49e25 -> 8.60e24 | 95 -> 68 % |
  | D reg | 8.12e-3 -> 4.44e-2 (x5.5) | 0.100 -> 0.263 | 0.079 -> 0.140 | 0.077 -> 0.122 | 2.70e26 -> 3.19e25 | 1.54e25 -> 5.40e24 | 98 -> 83 % |
  | C f1 / P_m 1e-4 reg | 2.44e-3 -> 1.30e-1 (x53) | 0.052 -> 0.340 | 0.047 -> 0.380 | 0.032 -> 0.535 | 3.32e27 -> 2.19e27 | 3.07e25 -> 1.65e25 | 27 -> 4 % |
  | D2 | 4.89e-3 -> 7.98e-2 (x16) | 0.078 -> 0.304 | 0.061 -> 0.259 | 0.045 -> 0.313 | 2.11e26 -> 2.75e25 | 1.21e25 -> 6.51e24 | 96 -> 73 % |
  | C f10 / P_m 0.01 reg | 3.20e-3 -> 9.61e-2 (x30) | 0.061 -> 0.311 | 0.052 -> 0.309 | 0.046 -> 0.368 | 1.51e25 -> 8.26e24 | 1.29e25 -> 7.40e24 | 44 -> 8 % |
  | softplus c9 | 4.27e-3 -> 6.48e-2 (x15) | 0.069 -> 0.262 | 0.061 -> 0.247 | 0.036 -> 0.525 | 5.38e26 -> 6.15e25 | 2.28e25 -> 7.62e24 | 98 -> 86 % |

* **What changes:** at fixed A the model is too low everywhere: mean ln(model/obs) over all used rows -0.21 ... -0.32
  (NG), -0.09 ... -0.37 (RG). Not only the dark records (7-18 % of the modelled surface pressure with 0.2, 2-21 % in the
  late window; 1-3 % with 0.02) but also the lit facets lose, since q_lit drops by 0.18 (1 - cos SZA) q_max. The seen
  production at 10:33:10 falls to 33-58 %, Q_H2O to 11-15 % (55-66 % in the two option-C states, whose Q came mostly
  from lit steep unseen terrain); the dark share of the seen production falls from 35-68 % to 6-20 %.
* **A level change alone does not restore the fits:** the best single factor on A (post-hoc bound, not a refit) is
  s = 1.17-1.43 and leaves the data cost 1.2 (A) ... 4.2 (C f1) times the old fit's. Every state needs a refit; the
  ranking of the models may change.
* **Caveats:** fixed-A numbers are not refit results (T, kappa, the multipliers and the shape of A will adapt); kappa
  is ~0 in all states, so CO2 does not compensate here; the unseen facets' activity (most of Q) is unconstrained.
* Evidence: `tests/nightfix_diagnostic_20261118/` (README.txt, table.txt / .csv / .json, results_v2/*.rows.npz; jobs
  136436-136443 v1, 136457-136464 v2; base check 136444); tool `model_ext/nightfix_diag.py`.
* **First refit with the new law: D2 (2026-11-18, user "(ii); refit D2 first, then all stage-2 chains", warm start "as
  recommended").** On `runs/shap8_012m_lategrid_dz1.5/base_q002` (the same kernels, hard links, bitwise-verified on 29
  re-run rows; COPS_model.dat of the rebuilt C++; fit_cops' base check NG 2.16e-9 / RG 4.12e-9; `--allow-legacy-base`
  because the base has no run_manifest.json), D2's settings unchanged, start = D2's final state (start data cost
  7.9817e-2 = the fixed-A diagnostic): data cost **5.671e-3** after 4000 evaluations (max_nfev, not converged;
  first-order optimality 1.2e-5), i.e. 1.16x D2's 4.886e-3 under the old law (and 0.071x D2's fixed-A value under the
  new law). rms ln NG / RG 0.0849 / 0.0643 (D2 0.0778 / 0.0610); late window (95 observations) NG 0.1763 / RG 0.0389
  (D2 0.1575 / 0.0446). c_cliff 3.727 -> 2.883; kappa and grain B at their lower bound 0 (as in D2); grain A 6.89 ->
  6.82 uPa; the harmonic coefficients moved little (22 of 25 within x0.81-1.18), active area 0.0755 -> 0.0859. So the
  intended law costs D2 about 16 % in data cost at its re-optimized state, much less than the fixed-A x16 suggested.
  Caveats: one model only; not converged; ft_c1..ft_c6 are held corrections from old-law fits; Q_H2O of the refit not yet
  computed. Evidence: `clip_refit/D2_q002_warm/` (README.txt, seg01/, provenance-136543); acceptance job 136542.
  Q_H2O of the refit (maps job 136692, 2026-11-18 13:28 CST): **3.12e25 s^-1 at 10:33:10** (3.23e25 at 10:39:26) vs
  D2's 2.11e26 with the old law: the re-optimized activity does NOT restore the total production; Q stays at ~15 % of
  D2's, close to D2's A evaluated with the new law (2.75e25). So the intended law lowers the model's Q_H2O by ~7x for D2
  while the fit to the gauges changes by only 16 %.
  Where the refit's production comes from (lead 2026-11-18 14:45 CST; source-rate figure set
  `figures_cluster/12M/source_rate_q002/D2_q002_warm/`, inspection/checks/inspect_103310.txt and _103926.txt, jobs
  136737 / 136741; new law, old-law D2 in brackets): at 10:33:10 the seen facets give 8.60e24 s^-1 (1.21e25) = 27.6 %
  of Q (5.8 %); night + cast-shadow facets give 68.8 % of Q (95.7 %) but only 6.5 % of the seen Q (39.4 %); dark area
  69.7 % of the surface, 60.2 % of the seen area; active area 8.59 % of the surface (7.55 %), 6.40 % of the seen area
  (5.31 %); 10:39:26 nearly the same (seen 28.0 %, dark 66.3 % / 6.0 %). So with the intended law the gauges' seen
  production is mostly lit (94 %), but about 2/3 of the total Q_H2O still comes from night-side facets, nearly all of
  them unseen, whose activity the data do not constrain: the total Q remains an extrapolation (A15 still applies).
  Caveats: D2 only (stage-2 chains pending); not converged (max_nfev); ft templates from old-law fits.
  First new-law stage-2 chain, softplus regularized (lead 2026-11-18 18:10 CST; chain
  nonneg/softplus/stage2r_q002_smm1_20261118T193956Z, job 136704; c5 fit ftol after 591 evaluations; figures
  source_rate_q002/softplus_q002_reg, inspect_103310.txt): data cost per cycle 5.256 / 5.229 / 5.221 / 5.218 / 5.216e-3
  (-8 % vs the new D2, +14 % vs the old-law softplus 4.593e-3). At 10:33:10: Q_H2O 6.08e25 s^-1 (old law 2.56e26,
  x0.24); seen facets 1.18e25 = 19.4 % of Q (5.8 %); night + cast shadow 76.6 % of Q (94.9 %) and 18.2 % of the seen Q
  (35.5 %); active area 5.71 % of the surface (5.99 %). Same pattern as the new D2: the seen production is mostly lit,
  the total is dominated by unseen dark facets (softplus more so than D2: 18 % vs 6.5 % dark in the seen Q).
  Option C f10 / P_m 1e-4 regularized with the new law (lead 2026-11-18 19:45 CST; chain
  nonneg/localnn/stage2r_q002_f10_pm1e-4_tplon_20261118T162407Z; c5/analysis/report.json; figures
  source_rate_q002/C_f10_pm1em4_q002_reg): data cost per cycle 2.780 / 2.551 / 2.528 / 2.511 / 2.4988e-3, -56 % vs the
  new D2 and only +1.6 % vs its old-law counterpart (2.4586e-3): the local basis absorbs the law change almost fully
  (D2 +16 %, softplus +14 %). Q_H2O(10:33:10) 4.84e26 s^-1 (old law 1.3e27), seen 2.14e25 = 4.4 % of Q; multipliers
  70-80 / 80-90 / 90-180 deg 5.46e4 / 1.79e5 / 38.2 (undetermined at f_m 1e-4, A12): the total Q is again set by
  unseen steep terrain times these factors, not by the data; kappa at its bound.
  Option C f1 / P_m 1e-4 and option D with the new law (lead 2026-11-18 21:30 CST; source_rate_q002/C_f1_pm1em4_q002_reg
  inspect_103310.txt; chains stage2r_q002_f1_pm1e-4_tplon_20261118T162407Z and logspread/stage2r_q002_f0.1_20261118T193956Z):
  C f1 final data cost 2.5005e-3 (+2.7 % vs old law 2.4358e-3); Q_H2O(10:33:10) 6.69e26 (old 3.32e27), seen 2.14e25 =
  3.2 % of Q; night + shadow only 3.7 % of Q (27 % old) and 13.4 % of the seen Q (44 %): C's total is lit unseen steep
  terrain times the undetermined multipliers, not dark facets; A > 0 on 100 % of the surface [A21: exact zeros do not measure switched-off surface]. D: every fit ftol, data cost
  8.388 / 8.359 / 8.307 / 8.283 / 8.277e-3 (+2 % vs old law 8.119e-3). D source rates at 10:33:10 (source_rate_q002/
  D_q002_reg): Q_H2O 3.25e25 (old 2.70e26), seen 41.1 % of Q (5.7 %); night + shadow 80.1 % of Q (97.6 %) and 52.9 % of
  the seen Q (67.6 %): unlike D2, softplus and C, D keeps half of its seen production on dark facets (its activity is
  A_orig exp(u), tied to the srcCG map); D is also the worst fit of the new-law runs so far (data cost 8.28e-3).
  Option C with P_m 0.01, new law (lead 2026-11-18 22:40 CST; chains stage2r_q002_*_pm0.01_*_20261118T162407Z, final
  c5/ft data cost): f10 templates on 3.5863e-3 (old law 3.1996e-3, +12 %), f1 templates on 4.1726e-3 (3.3779e-3, +24 %),
  f10 templates off 5.6734e-3 (3.6519e-3, +55 %; equal to the new D2, 5.671e-3) [A22: NOT a law effect, two different
  minima; +2.9 % within the same slope class]. So with the stronger multiplier
  penalty the law change costs 12-55 %, while with P_m 1e-4 option C absorbs it (+2-3 %) through the steep multipliers.
  Their source rates at 10:33:10 (source_rate_q002/C_*_pm0p01_*/inspection/checks/inspect_103310.txt; new law, old law in
  brackets): f10 on Q_H2O 9.71e24 (1.51e25), seen 88.4 % of Q, dark 9.6 % of Q / 6.1 % of the seen Q (43.6 / 40.8 %);
  f1 on 7.36e24 (1.52e25), seen 92.5 %, dark 10.7 / 7.0 %; f10 off 7.41e24 (1.98e25), seen 90.9 %, dark 9.9 / 5.5 %.
  So with P_m 0.01 option C's total production is ~90 % the production the gauges see, mostly lit: the one family of
  fits whose total Q (0.74-0.97e25 s^-1) is constrained by the data [A22: the templates-off model alone spans
  0.74-1.71e25 over its five solutions]. D2, softplus, D and C with P_m 1e-4 put 59-97 % of
  Q on unseen terrain (totals 3.1e25 ... 6.7e26). The seen production of all runs is 0.68-2.1e25 s^-1.
  Option B (background bg-min 1e-2) with the new law (lead 2026-11-19 01:25 CST; chain
  nonneg/background/stage2r_q002_bg1e-2_20261118T193956Z; c5 fit max_nfev, not converged; source_rate_q002/B_q002_reg):
  data cost per cycle 5.841 / 5.766 / 5.703 / 5.672 / 5.646e-3 (+10 % vs old law 5.120e-3; -0.4 % vs the new D2);
  Q_H2O(10:33:10) 4.28e25 (2.53e26), seen 1.05e25 = 24.6 % of Q (5.9 %), dark 67.7 % of Q (94.8 %) and 9.3 % of the
  seen Q (35.3 %): the D2 / softplus pattern (total dominated by unseen dark facets).
  Option A (sum of squares rank 3) with the new law, the last chain of the campaign (lead 2026-11-19 14:40 CST; chain
  nonneg/sos/stage2r_q002_r3_20261118T193956Z, every fit max_nfev; source_rate_q002/A_q002_reg): data cost per cycle
  11.57 / 10.79 / 10.36 / 9.982 / 9.640e-3, the only model that fits better than with the old law (-10 % vs 1.0694e-2),
  still the worst (+70 % vs the new D2). Q_H2O(10:33:10) 3.04e25 (1.84e26), seen 6.22e24 = 20.4 % of Q (6.4 %), dark
  50.6 % of Q (92.5 %) and 11.7 % of the seen Q (56.9 %).
  Campaign summary (all ten new-law states, 2026-11-19; comparison page source_rate_q002/cluster): data costs C 1e-4
  2.50e-3 (both f) < C f10 0.01 3.59e-3 < C f1 0.01 4.17e-3 < softplus 5.22e-3 < B 5.65e-3 ~ D2 5.67e-3 ~ C off 5.67e-3 <
  D 8.28e-3 < A 9.64e-3. The seen production at 10:33:10 is 0.62-2.1e25 s^-1 in every model; the totals span 7.4e24 ...
  6.7e26; only option C with P_m 0.01 has its total ~90 % seen (0.74-0.97e25). kappa at its bound in every fit.
  The slope multipliers move with the law (comparison set figures_cluster/12M/model_vs_observations_nonneg_options_
  stage2_q002/summary.txt, last fits, bins 40-50 / 50-60 / 60-70 / 70-80 / 80-90 / 90-180 deg, new law [old law]):
  softplus 9.68 / 10.4 / 3.03 / 1 / 6.75 / 44.2 [1 / 1 / 16.2 / 1 / 1 / 4.39]; A 2.34 / 1.05 / 1 / 1 / 1 / 1 [all 1];
  B 1 / 6.69 / 1 / 1 / 8.10 / 13.0 [1 / 1 / 11.85 / 1 / 1 / 4.37]; D 4.48 / 1 / 1 / 1 / 1 / 86.6 [14.6 / 1 / 14.9 / 1 / 1 /
  10.1]; C f1 / P_m 0.01 99.3 / 1.07 / 1.07 / 1.04 / 1.01 / 1.60 [1 / 1 / 1.02 / 290 / 1.01 / 1.66]; C f10 / P_m 0.01
  1 / 456 / 1.01 / 1.01 / 1.00 / 1.51 [1 / 1 / 1.01 / 759 / 1.00 / 1.50]; C P_m 1e-4 keep their large 70-80 / 80-90
  factors (f10 5.5e4 / 1.8e5, f1 2.5e4 / 8.1e4). So which slope bins carry extra activity (A16) is not robust to the
  night-side law: the bins with a factor > 1 change in every model with free multipliers [A22: can equally be different
  local minima; the law was not isolated]. Caveats: several last fits at
  max_nfev (A, B, C f1 / 1e-4); the multipliers of C at P_m 1e-4 are undetermined (A12); bin factors trade off against the
  harmonic map and the fine tuning.
  Full figure sets of the ten new-law states (2026-11-20, figures_cluster/12M/full_set_q002_<key>/checks/facts.json; old-law
  full sets in brackets where they exist): before the late window the big lobe (body) gives a median 87-96 % of the modelled NG
  [87-91 %] in every model; in the late window the small lobe (head) gives NG 25-64 %, RG 52-84 % [25-61 / 56-82 %]; the P50
  source-area unions of the 95 post-gap observations span 1.7e3 m^2 (C, P_m 1e-4) ... 7.1e5 m^2 (A) [1.5e3 ... 7.6e5]. The
  90-180 deg (overhang) bin gives 63-86 % of NG at the last observation (5.5 m) in softplus, A, B, D and C with templates
  [59-80 %], but only 7 % for the refitted D2 [38 %] and 0 % for C f10 / P_m 0.01 with the templates off; RG 0-1 % in all.
  So the NG overhang dominance at touchdown is a property of models with activity on overhangs (multiplier or template), not
  required by the data [TOO STRONG, A21: the fit without templates also misses the NG rise of the last ~20 m].
  Caveats: A, B, D2 and C f1 / 1e-4 stopped at max_nfev; ft templates held from old-law fits; the
  C 1e-4 multipliers are undetermined.

### A17. Night and cast shadow lower the modelled production by up to x5, and at the map times most of the production is dark (2026-11-18)

* **Question** (user 2026-11-17): the production-rate maps do not show the night-side / shadow decrease; is it there?
* **Correction 2026-11-18 (A18):** this entry describes the q_min = 0.2 q_max law of all fits so far; the intended law
  (srcCG, the code since 2026-11-18) has 0.02, i.e. a factor up to 50 (1.70 decade) between dark and overhead-Sun facets.
* **The law** (checked on every facet at 10:33:10 and 10:39:26 UTC): F_H2O = q_H2O A, q_H2O/q_max = 0.2 on night and
  cast-shadow facets, 0.2 + 0.8 cos SZA on day facets (F = q A bitwise; q_min/q_max = 0.2). At a given A illumination changes
  F by at most a factor 5 (0.70 decade), ~2.8 of the 0.25-decade colour classes of the maps.
* **Why it is not visible:** the map colour is set by A: on the seen facets with F > 0 (10:33:10, softplus reg) log10 F spans
  3.25 decades (area-weighted 2-98 %), log10 A 3.01, log10 q/q_max 0.66; A = 0 on 96.4 % of the seen area.
* **Most of the production is dark at these times:** 66.5 % of the seen area with F > 0 is night or cast shadow (35.5 % of
  the seen production); on the whole surface night / shadow is 69.7 % of the area, 85.4 % of the area with F > 0 and 94.9 % of
  Q_H2O (94.2 % of Q_H2O from facets no fitted observation sees, 98.5 % of that from their dark facets). Production-weighted
  q/q_max = production / (q_max sum A area) = 0.207 (whole), 0.347 (seen): Q_H2O is ~21 % of what the same A would give under
  an overhead Sun. The impact-area active patch is lit (72 % of its drawn area at q/q_max >= 0.8, 0.9 % dark).
* **Fit-independent geometry:** cos SZA, illumination and the seen mask are bitwise identical in the maps of all six full-set
  states; 30.3 % day / 51.6 % night / 18.1 % cast shadow by area at 10:33:10; the Sun moves 2.92 deg between the two map times
  and 4.3 % of the area changes class.
* **Caveats:** the activity of the unseen facets (which carry ~94 % of Q_H2O) is not constrained by the data; A is constant in
  time while the illumination rotates, so a facet dark at 10:33 was lit at other fitted times; the numbers are for the softplus
  reg state (the illumination itself holds for all six).
* Evidence: `figures_cluster/12M/illumination_maps_103310_103926/` (README.txt, checks/illumination_numbers.txt,
  checks/illumination_crosscase.txt, side_by_side/); jobs 136359 (numbers), 136369 / 136370 (figures).

### A16. At the same illumination the fitted source rate depends on the slope in all six final states; the size and the bins are the model's (2026-11-17)

* **Question** (user 2026-11-17): at the same illumination, does the source rate depend on the slope angle? Definitions
  0-40 (reference), 40-50 ... 80-90, >40, >60, >80 (and 90-180 separately). In the model F = A q with q = Qmin + (Qmax -
  Qmin) cos SZA on the lit side and Qmin at night / in cast shadow, so at equal SZA and shadow state the slope dependence of
  F is that of the fitted A (`figures_cluster/12M/cliff_effect_stage2/`, answer.json per case; document
  12M_stage2_cliff_effect.docx).
* **SZA-controlled ratio F(category) / F(0-40)** (area-weighted average over 10-deg SZA classes + shadowed + night of the
  per-class ratio), whole surface, 10:33:10 UTC, >40 / >60 / >80: softplus 4.66 / 8.84 / 2.48, A 1.43 / 1.39 / 1.52,
  B 4.12 / 7.55 / 2.73, D 12.6 / 10.2 / 3.83, C (f 1, P_m 1e-4) 9.86e4 / 2.16e5 / 5.44e5, D2 1.92 / 2.50 / 1.72. Single bins:
  softplus 60-70 18.6 and 90-180 4.11 (fitted multipliers 16.2, 4.389), B 15.2 / 4.38 (11.85, 4.372), D 40-50 24.9, 60-70
  21.2, 90-180 7.65 (14.64, 14.92, 10.1), C 70-80 1.9e5, 80-90 1.01e6 (7.904e4, 3.303e5); A (multipliers all at 1) 0.867-2.07
  and D2 (no multipliers) 1.32-2.93: the map's own slope trend. F and A ratios agree to <= 0.42 % (whole surface).
* **So:** in every model the >40, >60 and >80 categories emit more than the 0-40 deg facets at the same illumination
  (single bins can be below 1: A 90-180 0.867, D 70-80 0.745 and 80-90 0.586), but how much and in which bins is set by the
  fitted multipliers (where they are free) and is not robust across models: 1.4x (A) to 8-10x (softplus, B, D) for >60
  deg, and 2e5 for C, whose multipliers are not determined (A12).
* **In the signal** (P4, >60 deg, NG, 10:33:10; 15.3 % of the seen area): share / area share 0.779 A, 2.00 D, 3.73
  softplus, 4.24 D2, 4.50 B, 6.53 C. The illumination factor (0.696) and the geometry factor (0.617) are the same in every
  model (vs the 0-40 deg facets the >60 deg facets are farther from the gauge, 1/r^2 x 0.41, but face it better, cos theta
  x 1.48); the activity factor (1.13 A ... 6.53 C) and the covariance (1.6-2.4) make the difference.
* **Caveats:** a model property, not a measurement: the slope dependence comes from multipliers and maps the data constrain
  only on the seen surface; A, B, C and D2 stopped at max_nfev (not converged); C's multipliers are undetermined at f_m 1e-4;
  kappa at its lower bound everywhere. On the gauge-seen surface the 0-40 reference is sparse at night / high SZA and at 20 m
  and 5.5 m in several classes; there the specified ratio is unstable (it can exceed 1e100), so the whole-surface numbers
  above are the robust ones (variant ratios over well-sampled classes in answer.csv).

### A15. Full figure sets of the six final states: what is common to all models and what is not (2026-11-17; corrected 19:00 CST)

* **States:** the regularized stage-2 chains softplus, A (sum of squares r3), B (background 1e-2), D (log spread f 0.1),
  the best regularized C (f 1, P_m 1e-4) and D2 (`figures_cluster/12M/full_set_<key>/`, checks/facts.json, computed by
  full_set_facts.py from the maps, contributions and contrib_phase products of each state).
* **The SEEN production is similar across models; the total is set by the unseen extrapolation.** H2O production at
  10:33:10 UTC of the facets that contribute to a fitted observation (seen) and of all facets (total), s^-1 (full sets:
  checks/facts.json from the maps; option C: c5/analysis/report.json production_1033, incl. its c9 / D2 rows):
  softplus reg. 1.49e25 of 2.56e26 (5.8 %), B 1.49e25 / 2.53e26 (5.9 %), A 1.17e25 / 1.84e26 (6.3 %), D 1.54e25 / 2.70e26
  (5.7 %), D2 1.21e25 / 2.11e26 (5.8 %), c9 2.28e25 / 5.38e26 (4.2 %), C reg. f 1 / P_m 1e-4 3.07e25 / 3.32e27 (0.9 %),
  C reg. f 10 / 1e-4 2.98e25 / 1.30e27 (2.3 %); but the option C chains with f_m 0.01: C reg. f 10 / 0.01 1.30e25 / 1.51e25
  (86 %), f 1 / 0.01 1.30e25 / 1.52e25 (86 %), f 10 / 0.01 templates off 1.47e25 / 1.98e25 (74 %), unregularized f 10 / 0.01
  3.21e25 / 3.37e25 (95 %), f 1 / 0.01 2.98e25 / 3.13e25 (95 %). So the seen production is about 1.2-3.2e25 s^-1 in every
  model, while Q_H2O varies by more than 100x (1.5e25 ... 3.3e27) with the activity each model extrapolates onto terrain no
  fitted observation sees: ~94 % unseen in the six full-set states except C f 1 / 1e-4 (99 %), most of it from unseen facets
  >= 60 deg there (the large undetermined 70-80 / 80-90 deg multipliers, 7.9e4 / 3.3e5, times the smoothness extension of
  the c-field), whereas the multiplier penalty f_m 0.01 keeps the unseen part small. Q_H2O is therefore not constrained by
  COPS; the seen production is. (The seen share does not depend on whether the seen surface is mostly switched off - A = 0
  on 96.4 % / 94.7 % of the seen area in softplus / D2 - or active everywhere - A, B, D: A = 0 on <= 0.9 %.) [A21: exact zeros do not measure switched-off surface.]
* **Robust across the six:** before the late window the big lobe (body) gives a median 87-91 % of the modelled NG; at the
  last observation (5.5 m) the 90-180 deg (overhang) bin gives 38 % (D2) to 80 % (D) of NG but only 1-2 % of RG (NG sees
  through the shield onto facets behind the gauge; A1, A11).
* **Model-dependent:** the late-window small-lobe (head) share (median NG 25 % D2 ... 61 % C, RG 56 % ... 82 %) and the size
  of the P50 source areas of the post-gap observations (NG union 1.5e3 m^2 C, 1.9e4 softplus, 2.8e4 B, 1.4e5 D2, 1.7e5 D,
  7.6e5 m^2 A): the data constrain the integrated signal, not where on the near surface it comes from.
* **kappa** is at its lower bound in all six (A14).
* **Caveats:** C's multipliers are not determined (A12), so its Q_H2O is not a measurement; B, A and C's last fits and D2
  stopped at max_nfev; Q_H2O and the unseen share depend on the model's extension of the activity to unseen terrain (and on
  the H2O law); the D2 state is the export of the clipped fit (A re-derived bitwise from its restart; forward model = D2's
  per-row model to 3.1e-11 / 2.4e-15, data cost to 1e-15; clip_refit/D2/tools_state_verify.json).

### A14. Stage 2 (5 fit / fine-tuning cycles): option C stays best by ~47 %; the regularized fine tuning costs 5–13 % data cost, mostly in the last ~100 m; κ is at 0 in every fit (2026-11-17)

* **Setup** (user / lead 2026-11-15/16): each chain = 5 cycles of a fit (max_nfev 4000; slope multipliers 40-50 ... 90-180
  deg ≥ 1; D2's ft_c1..ft_c6 held at 1 unless stated) and a late-facet fine tuning (271 facets seen only in the last 5 min).
  Regularized chains (`stage2r_*`): P_delta = f_delta C_ref mean_w (h + delta)^2, f_delta 10, cap |h + delta| ≤ ln 100,
  templates on. Unregularized chains (`stage2_*`; softplus: chain_5cycles_smm1) as context. Final state = c5 fine tuning.
* **Final data cost, regularized (unregularized):** C f 1 / P_m 1e-4 2.4358e-3 (2.2973e-3); C f 10 / P_m 1e-4 2.4586e-3
  (2.2543e-3); C f 10 / P_m 0.01 3.1996e-3 (3.0375e-3); C f 1 / P_m 0.01 3.3779e-3 (2.9946e-3); C f 10 / P_m 0.01 templates
  off 3.6519e-3; softplus 4.5925e-3 (4.2861e-3 after 5 cycles, 4.2688e-3 after 9); B bg 1e-2 5.1201e-3 (4.7834e-3); D f 0.1
  8.1192e-3 (7.2138e-3); A rank 3 1.0694e-2 (9.4386e-3); D2 4.8864e-3. The best regularized C chain is 47 % below the
  regularized softplus and 50 % below D2; the regularization raises the data cost by +5.3 % (C f 10 / 0.01), +6.0 %
  (C f 1 / 1e-4), +7.0 % (B), +7.1 % (softplus), +9.1 % (C f 10 / 1e-4), +12.6 % (D), +12.8 % (C f 1 / 0.01), +13.3 % (A).
  rms ln NG / RG (late window): best C reg. 0.0516 / 0.0470 (0.116 / 0.032), softplus reg. 0.0739 / 0.0610 (0.141 /
  0.042), B 0.0779 / 0.0646, D 0.0999 / 0.0791, A 0.1096 / 0.0968.
* **Where the regularization costs:** in the last ~100 m. Lowest NG residual ln(model/obs) between 100 and 10 m,
  regularized / unregularized: softplus -0.31 / -0.15, C f 10 / 1e-4 -0.24 / -0.18, C f 10 / 0.01 -0.30, B -0.36, D -0.40,
  A -0.48 (unreg. B -0.22, D -0.22, A -0.33); RG below 20 m regularized up to +0.11 (softplus, C P_m 1e-4), +0.19 (C P_m
  0.01, B), +0.34 (D), +0.57 (A); unregularized within about +0.1 (A +0.18). The capped, penalized late facets can no
  longer supply the steep NG rise and the RG level near touchdown at the same time.
* **Multipliers:** P_m 1e-4 lets 70-80 / 80-90 deg run to 7.9e4-8.7e4 / 3.3e5-3.4e5 (regularized; unregularized up to
  5.3e5) - not determined (A12 degeneracy); P_m 0.01 keeps one large bin (70-80 deg 290-759, templates off 80-90 deg 1305).
  softplus reg. 60-70 deg 16.2, 90-180 deg 4.39 (the others at 1); B reg. 11.9 / 4.37; D reg. 40-50 14.6, 60-70 14.9,
  90-180 10.1; A all six at 1.
* **Convergence:** every last fine tuning stopped on ftol (regularized 2-16, unregularized 4-46 evaluations). Last fits: ftol except B reg. (max_nfev),
  A reg. (max_nfev, second segment after a time limit), C f 1 / 1e-4 reg. (max_nfev); unregularized softplus c5 / c9 and
  A at max_nfev, C f 10 / 1e-4 by xtol. Reaching max_nfev is not convergence.
* **κ (CO2 scale on the srcCG CO2 map; 1 = the map, lower bound 0) is at its lower bound in all 59 fits checked**
  (`kappa_check/kappa_table.txt`): the last fits of the stage-2 chains: softplus reg. 1.9e-40, B reg. 4.2e-36, D reg.
  1.2e-45, A reg. 4.6e-12, C reg. 4.7e-22 (f 10 / 0.01), 4.2e-24 (f 1 / 0.01), 1.5e-20 (templates off), 1.7e-37 (f 10 /
  1e-4), 1.8e-51 (f 1 / 1e-4); unregularized: softplus c5 2.7e-37, c9 1.7e-35, B 1.7e-45, D 7.1e-37, A 4.0e-13, C 5.6e-23,
  5.6e-35, 5.9e-25, 6.5e-34. **It was already so before stage 2:** D2 2.4e-36 (fit_result.txt kappa row = fit_restart
  value), the unbounded softplus chain c5 3.7e-37, the bounded softplus c1 / c5 / c9 4.1e-37 / 2.7e-37 / 1.7e-35, every
  stage-1 fit and continuation of B, A, D (3.6e-47 ... 1.3e-12) and of C (1.2e-23 ... 7.2e-9; the C report's criterion
  κ ≤ 1e-8); each first fit of a chain starts at κ ≈ 0 already (A5). Caveat (as A5 / A12): the gauges measure total
  pressure, CO2 and H2O differ only through the production law, so κ = 0 says the fits prefer the H2O law for the whole
  signal, not that CO2 is measured absent; several fits ended at max_nfev.
* **Evidence:** `figures_cluster/12M/model_vs_observations_nonneg_options_stage2_final/` (summary.txt / .csv,
  12M_stage2_*_residuals, kappa_check/; stage2_figures.py, CPU job 136130); the chain roots and their README.md /
  HOW-TO-RERUN.md (`runs/shap8_012m_lategrid_dz1.5/nonneg/<model>/stage2r_*`, `stage2_*`).
* **Caveats:** one chain per setting; C on CPU (the P_m 1e-4 chains: cycles 1-3 and part of 4 on CPU, the rest on the
  GPU); the C multipliers are not determined; the fine tuning acts on 271 facets only.

### A13. Forcing more of the seen surface to be active costs fit quality; the best fit (local basis) concentrates production on steep facets (2026-11-15, stage 1)

* Four non-negative activity parameterizations, all with the 7 slope-bin multipliers (0-40 held 1, six >= 1), D2's
  ft templates (sum >= 0 on 265 facets) and D2's setup, stage-1 single fits from D2 (B/A/D continued once from their
  end points with the optimizer scale rebuilt; D with the fixed-level spread penalty). Data cost vs the active seen-area
  fraction (A > 1 % of the srcCG median of the seen facets):

  | model | data cost | vs softplus c9 | active seen area | converged |
  |---|---|---|---|---|
  | C local hat basis 150 m, smoothness f = 1 | **2.689e-3** | −37 % | 12.5 % | max_nfev |
  | C, f = 10 / f = 0.1 | 2.831e-3 / 2.833e-3 | −34 % | 20.1 / 13.6 % | max_nfev |
  | C, f = 10 + multiplier penalty 0.01 | 3.266e-3 | −23 % | 24.5 % | ftol |
  | softplus c9 (reference) | 4.269e-3 | — | 2.6 % | — |
  | B background >= 1e-2 A_orig | 5.366e-3 | +26 % | 55 % | ftol |
  | D log form, f = 0.1 (fixed level) | 8.390e-3 (+1.7e-3 penalty) | +97 % | 88 % | ftol |
  | A sum of squares rank 3 | 1.271e-2 | ×3.0 | 96 % | max_nfev |
  | B background >= 1e-1 | 6.617e-2 | ×15 | 88 % | max_nfev |

* **The data allow a fit 37 % better than the softplus**, but only with a local basis, and that fit still makes 50 % of
  the seen production on 0.06 % of the seen area (steep 60-90 deg facets, A12). Models that force activity over most of
  the seen surface (B, D, A with 55-96 % active) fit 1.3-15x worse. So "more seen surface contributing" is not favoured
  by the data under the present physics (production law with q_min = 0.2 q_max, original T law, gauge response without
  spacecraft occlusion).
* In the last ~100 m every new model is below the NG data (C −0.25, A −0.47 in ln) and above the last RG points; the
  softplus c9 is better there because its 9 cycles include late-facet fine tuning (stage-1 fits have none).
* Caveats: stage 1 only (no fine-tuning cycles); most fits not converged (max_nfev); multipliers not determined (A12).
* Evidence: `figures_cluster/12M/model_vs_observations_nonneg_options_stage1_final/` (summary.txt/csv, figures,
  activity maps), run roots under `runs/shap8_012m_lategrid_dz1.5/nonneg/{background,sos,logspread,localnn}/`,
  `checks/stage1_diag_20261114/`.

### A12. A local non-negative basis (option C) fits better than D2 and c9 and puts the activity on 60–90° facets; the multipliers' absolute values are not determined (2026-11-14, preliminary)

* **Setup** (lead / user 2026-11-14): A = m_b (sum_j c_j phi_j + A_s sum ft_ck), c_j >= 0 on 2745 hat functions of a 150 m
  coarse triangulation, the 7 slope-bin multipliers (0-40 deg held at 1, six >= 1), D2's templates held at 1, D2's data,
  W10x5, grain and crossing; smoothness lambda_C x_c^T Lg x_c with lambda_C = f C_ref / R_pilot (f = 0.1, 1, 10; C_ref
  = 4.2688e-3, R_pilot = 517.46 of a 200-evaluation pilot); start = the NNLS projection of D2; max_nfev 4000 (not
  converged). f = 1 is still running (resubmitted after a node failure).
* **Fit quality:** data cost f 0.1: 2.8331e-3, f 10: 2.8306e-3, i.e. -34 % vs the bounded softplus c9 (4.2688e-3) and
  -42 % vs D2 (4.886e-3); rms ln NG / RG 0.0553 / 0.0510 and 0.0558 / 0.0505 (c9: 0.0690 / 0.0614); late window NG 0.123
  / 0.121, RG 0.064 / 0.050. Already the 200-evaluation pilot reached 3.80e-3.
* **More of the seen surface is active, but the production is more concentrated:** A > 1 % of the srcCG median on 13.6 %
  (f 0.1) / 20.1 % (f 10) of the seen area (D2 5.3 %, c9 2.6 %); yet 50 % of the seen H2O production at 10:33:10 comes
  from 0.05 % / 0.09 % of the seen area (D2 0.71 %, c9 < 0.01 %) and 90 % from 0.33 % / 0.66 %.
* **Where:** the 60-90 deg bins give 92-94 % of the late NG and 86-90 % of the late RG model signal (whole descent:
  NG 96-98 %, RG 95-97 %); mean activity on the seen facets relative to 0-40 deg: 60-70 deg 17 / 134, 70-80 deg 345 /
  979, 80-90 deg 1067 / 2156 (f 0.1 / f 10); flat terrain is ~10-40x less active than in D2.
* **Not determined:** the c-field is nearly the same on every bin at 150 m (it cannot follow the slope pattern), so the
  multipliers alone make the steep / flat contrast (m 70-80 = 675 / 1608, m 80-90 = 2269 / 3617). No 150 m node is
  dominated by a steep bin, but the global level of the c-field is free under the Laplacian penalty: lowering all c_j
  and raising the six multipliers by the same factor changes the residual vector by only 2.3-4.3 % of the multiplier
  step alone, and each multiplier is compensated by the c_j to 0.2-0.3 %. So the multipliers and the flat-terrain
  activity level are weakly constrained (only by the 0-40 deg share of the signal, 4-6 % late NG, 9-13 % late RG); the
  steep-facet activity is what the data fix.
* **With the multiplier penalty** (user 2026-11-14 ~17:30 CST: P_m = 0.01 C_ref sum (ln m_b)^2; one continuation of
  the f = 10 fit, job 125360, converged by ftol after 219 evaluations): data cost 3.2655e-3 (+15 % vs the f = 10 end
  point, still -24 % vs c9), P_m 1.96e-3; the 60-70 and 80-90 deg multipliers drop to 1.01 / 1.003, 40-60 deg stay at
  the bound 1, and only 70-80 deg keeps a large one (827; 90-180 deg 2.4): one steep bin then carries 86 % of the seen
  production and 88 % / 80 % of the late NG / RG signal. So the data accept either several steep bins (no P_m) or
  essentially one (with P_m) at a 15 % higher data cost; which steep bin carries the signal is not fixed by the data.
* **Stage 2, chain f = 10, f_m = 0.01 complete** (2026-11-15, `runs/.../nonneg/localnn/stage2_f10_pm0.01_20261115T190019Z`,
  5 fit / fine-tuning cycles, every fit by ftol: 219, 78, 33, 18, 12 evaluations; fine tunings 243 (ftol), 500
  (max_nfev), 5, 5, 5 (ftol)). The full objective (data + smoothness + P_m; crossing 0) falls monotonically:
  9.3738e-3 (fit start) -> 5.3947e-3 (c1 fit) -> 5.0178e-3 (c1 ft) -> 4.9440e-3 -> 4.9255e-3 -> 4.9191e-3 -> 4.9161e-3
  -> 4.9141e-3 -> 4.9125e-3 -> 4.9112e-3 (c5 fit) -> 4.9099e-3 (c5 ft); each fit starts exactly at the previous fine
  tuning's objective. The data cost alone rises after c1: 2.8886e-3 (c1 ft) -> 3.0122 / 2.9937 -> 3.0333 / 3.0303 ->
  3.0388 / 3.0372 -> 3.0387 / 3.0375e-3 (fit / ft): the fits give back data cost for smoothness and P_m (P_m 1.96e-3 ->
  1.65e-3, smoothness 1.69e-4 -> 2.22e-4). Final multipliers 40-50, 50-60 at the bound 1, 60-70 1.013, 70-80 473,
  80-90 1.008, 90-180 2.35; kappa 0. The c1 fine tuning set delta in [-10.7, +11.3] on the 271 late facets (5 with
  h > 5, 17 with h < -5; later cycles change h by <= 0.9); final max h 11.8 (x 1.4e5). c5 fit state: Q_H2O(10:33:10)
  3.37e25 s^-1, active seen area 30.6 %, 50 % / 90 % of the seen production from 3.0 m^2 / 1.16e5 m^2 (a few
  fine-tuned late facets); final state (c5 ft): seen production at 10:33:10 by bin 63 % 90-180 deg, 35 % 70-80 deg,
  1.2 % 0-40 deg; late-window signal NG 81 % / RG 72 % from 70-80 deg (0-40 deg 11 % / 22 %, 90-180 deg 7 % / 4 %);
  NG P50 source area 6.9e3 m^2 at 10:33:10, 2.0 m^2 (2 facets) at 10:39:26. Caveat: the fine tuning has no lower
  bound and lambda 0 (user decision), so the extreme per-facet factors are allowed by construction; one chain of
  four; GPU cross-check of the final state pending.
* **Regularizing the fine tuning** (test 2026-11-15, the c5 fine tuning re-run with P_delta = f_delta C_ref
  area-mean (h+delta)^2 and the cap |h+delta| <= ln 100; `tests/localnn_20261114/ftreg/summary_134507.txt`): without
  it the accumulated per-facet log factors reach |h| 11.8-12.3 (C, B, D) and 50 % of the seen production at 10:33:10
  sits on 1-3 m^2. The cap alone clips 20-37 of the 271 tuned facets at the start (start data cost +35-50 %). Data
  cost after the regularized fine tuning vs unregularized: C f1 / 0.01 2.995e-3 -> 3.355 / 3.423 / 3.825e-3 (f_delta
  0.1 / 1 / 10), C f10 / 0.01 3.037e-3 -> 3.239 / 3.283 / 3.568e-3, B 4.783e-3 -> 4.928 / 4.954 / 5.239e-3, D 7.214e-3 ->
  7.737 / 7.833 / 8.305e-3; late rms ln NG 0.113 -> 0.131-0.149 (C f10), RG up to 0.05-0.11; max |h+delta| 4.6 (cap)
  / 2.6-3.2 / 0.6-1.1. Only f_delta = 10 spreads the production: seen area holding 50 % 1.03 m^2 (one facet) for
  f_delta 0.1 and 1, 3.3e4-1.8e5 m^2 for 10 (C Q_H2O(10:33:10) 3.4e25 -> 1.4e25; B, D change < 10 %). The tuned facets
  carry 3-13 % of the mean late-window surface signal but 80-90 % (NG) / 30-70 % (RG) at the last observation
  (10:39:26), regularized or not. All runs ftol; no chain was rerun.
  The 1.03 m^2 holder is ONE tuned facet, 12232930 (90-180 deg bin), in C as in B and D (agent 1): in C its A is
  m_90-180 (2.3-2.7) x TS (D2's fine-tuning template value 7.894e23, u ~ 0) x exp(h), with h only 2.2 (C f1) / 2.9
  (C f10) after the chains; the concentration comes from D2's templates held at 1, not from large h. With the
  area-weighted penalty this small facet is cheap: at f_delta 0.1 / 1 its h + delta RISES to 3.3 / 3.1-3.2 (A x 2-3);
  only f_delta 10 brings it to 0.1-0.16 and the largest A moves to facet 11213754.
* **Stage-2 reruns with the regularized fine tuning (f_delta 10, cap ln 100), D2's templates on vs off** (user
  decisions 2026-11-15/16; f = 10, f_m = 0.01, 5 cycles, all ftol; `tests/localnn_20261114/stage2r/chain_summary_135107.txt`,
  `chain_summary_135114.txt`, `share_c5_*.txt`). Templates ON (`stage2r_f10_pm0.01_tplon_20261116T044825Z`): data cost
  3.2655e-3 (c1 fit) -> 3.1700e-3 (c1 ft) -> ... -> 3.1996e-3 (c5 ft), objective (data + smoothness + P_m) monotone
  5.3947e-3 -> 5.2623e-3, P_delta 2.6-3.7e-5, max |h+delta| 0.93; late rms ln NG / RG 0.130 / 0.046; multipliers 1, 1,
  1.01, 759 (70-80 deg), 1.004, 1.50 (90-180 deg, 2.42 at c1); Q_H2O(10:33:10) 1.51e25 s^-1; active seen area 25.4 %;
  50 % / 90 % of the seen production on 4.0e4 / 1.8e5 m^2; facet 12232930 holds 4.1 % (43.6 % unregularized). The
  tuned late facets carry 2.4 % of the mean late NG signal and 75 % at the last observation. Templates OFF
  (`stage2r_f10_pm0.01_tploff_20261116T044825Z`, c1 fit from the f = 10 end point minus the template part: start data
  5.458e-3): c1 fit 238 evaluations -> data 3.6519e-3, objective 6.0628e-3 (P_m 2.20e-3, smoothness 2.14e-4), late rms
  ln NG / RG 0.161 / 0.049; then nothing moves (fine tunings 2-3 evaluations, |delta| <= 5e-4, fits 3 evaluations):
  without the templates the 271 late facets carry 0.01 % of the mean late signal and 0.13 % (NG) / 0.27 % (RG) at the
  last observation, so the fine tuning has no leverage. Multipliers 1, 1, 1.013, 1.014, 1305 (80-90 deg), 1.003: the
  steep activity moves from 70-80 to 80-90 deg; Q_H2O 1.98e25; active seen area 28.5 %; 50 % / 90 % of the seen
  production on 2.1e4 / 1.1e5 m^2 (largest facet 0.6 %). Data cost: templates OFF 3.65e-3 > ON (regularized) 3.20e-3 >
  ON unregularized 3.04e-3; all below the softplus c9 4.27e-3. So the templates (D2's late fine-tuning facets) are what
  the late signal uses; without them the late window fits worse (rms ln NG 0.161 vs 0.130) and the fine tuning cannot
  help.
  Same templates-ON regularized setup for f = 1, f_m = 0.01 (`stage2r_f1_pm0.01_tplon_20261116T055923Z`, 22 min, all
  ftol; `tests/localnn_20261114/stage2r/chain_summary_135382.txt`): data 3.4438e-3 (c1 fit) -> 3.3779e-3 (c5 ft),
  objective 5.0205e-3 -> 4.9098e-3 monotone, late rms ln NG / RG 0.140 / 0.044, multipliers 1, 1, 1.018, 290, 1.008,
  1.66; Q_H2O 1.52e25; active 22.4 %; 50 % of the seen production on 3.8e4 m^2 (facet 12232930 4.6 %); unregularized
  f1 / 0.01: data 2.9946e-3, objective 4.4228e-3, 50 % on 3.0 m^2 (two facets, 27 % / 30 %).
* **Regularized f_m = 1e-4 chains** (templates on, f_delta 10 + cap; 2026-11-17; c1 = the original stage-2 c1 fits,
  c4 partly and c5 on the GPU; `tests/localnn_20261114/stage2r/chain_summary_136129.txt`): f = 10: data 2.6315e-3 (c1
  fit) -> 2.4586e-3 (c5 ft), objective 2.7299e-3 -> 2.5930e-3; f = 1: 2.5521e-3 -> 2.4358e-3, 2.6525e-3 -> 2.5703e-3
  (c1-c3 fits and f1 c4 / c5 at max_nfev 4000: not converged); late rms ln NG / RG 0.115 / 0.032 (both); multipliers
  40-70 deg ~1, 70-80 deg 8.7e4 / 7.9e4, 80-90 deg 3.4e5 / 3.3e5, 90-180 deg 3.5; Q_H2O(10:33:10) 1.3e27 / 3.3e27 s^-1;
  active seen area 4.4 % / 4.8 %; 50 % / 90 % of the seen production on 9.0e3 / 5.3e4 m^2 (f 10) and 7.7e3 / 5.2e4 m^2
  (f 1); facet 12232930 holds 4.2 % / 4.1 %. Unregularized f_m = 1e-4 (same c1 fits): data 2.2543e-3 / 2.2973e-3,
  50 % on 443 / 895 m^2, facet 12232930 27 % / 23 %. So with the weak multiplier penalty the data cost is the lowest of
  all option C runs (regularized 2.44-2.46e-3 vs 3.20-3.38e-3 at f_m = 0.01), bought with steep-bin multipliers of
  1e5 and Q_H2O 1e27+: Q is then set by unseen steep terrain multiplied by those factors, not by the data (caveat), and
  the multipliers are only weakly determined (the global-level degeneracy).
* **Night law 0.02 (2026-11-18, A18): option C design** (job 136507): with H2O q_min = 0.02 q_max the design columns of
  the hat functions fall to 0.1 - 1 of the 0.2-law values (median 0.79 NG, 0.80 RG; the template columns 0.71; CO2
  unchanged) and the late-facet kernels to a median 0.1 (most records of the 271 late facets are night / shadow), so the
  late fine tuning acts on activity that now produces about 10 x less. The 0.2-law D2 / c9 states evaluated with the
  new design: data cost 0.0798 / 0.0648 (were 4.89e-3 / 4.27e-3); all option C results above are 0.2-law results.
* **kappa = 0 in every option C stage-2 state - and in every other fit** (c5 fits, 2026-11-17; widened 2026-11-17 15:00
  CST: kappa is at its lower bound in all 59 fits checked, all options, the softplus chains, D2 and the stage-1 fits, and
  was so before stage 2; values in A14 and `figures_cluster/12M/model_vs_observations_nonneg_options_stage2_final/
  kappa_check/kappa_table.txt`; the ninth C state, stage2r f1 / 1e-4: 1.8e-51). In C, kappa (the CO2 scale, bound >= 0) is
  5.6e-23 (stage2 f10 / 0.01), 5.6e-35 (f10 / 1e-4), 5.9e-25 (f1 / 0.01), 6.5e-34 (f1 / 1e-4), 1.5e-20 (stage2r f10 /
  0.01 templates off), 4.7e-22 (templates on), 1.7e-37 (stage2r f10 / 1e-4), 4.2e-24 (stage2r f1 / 0.01): at the bound
  in all eight (Q_CO2 ~ 0), as in the earlier fits (A5). With f_m = 1e-4 the steep multipliers run to 1e4 - 5e5
  (stage2r f10 / 1e-4: 70-80 deg 8.67e4, 80-90 deg 3.40e5; stage2 f10 / 1e-4: 1 / 5.3e5); with f_m = 0.01 one bin
  carries 2e2 - 8e2 (templates off: 80-90 deg 1305). Caveats: the gauges measure total pressure, so CO2 and H2O differ
  only through the production law (q_max ~ r_h^-2 vs r_h^-4.21, the night floor 0.1 vs 0.2 q_max) and kappa = 0 says
  that the data prefer the H2O law for the whole signal, not that CO2 is measured to be absent; several fits ended at
  max_nfev (not converged); the f_m = 1e-4 multipliers are weakly determined (the global-level degeneracy, above).
* **Evidence:** `runs/shap8_012m_lategrid_dz1.5/nonneg/localnn/stage1_anch_f{0.1,10}_*/analysis/report.txt`,
  `analysis_degeneracy/degeneracy.txt`, `gpu_eval/gpu_check.json` (f 10: GPU standard evaluation = CPU to 1.2e-15 in
  the cost); jobs 125107/125108, 125119/125120, pilot 125071; README.md "Option C".
* **Caveats:** not converged (max_nfev); f = 1 pending; Q_H2O differences between the models are dominated by unseen
  terrain (the smoothness extension), not by data.

### A11. Spacecraft occlusion of the gauges: none for RG, up to ~45 % of the modelled NG signal at the very end (2026-11-13)

* The fit has no spacecraft occlusion (LiouvilleCOPS README §4.5: the NG shield table "is an instrument-response
  approximation, not a spacecraft CAD occlusion model"). Test with the SPICE spacecraft model (bus box from
  ROS_SPACECRAFT_C1..C8, ten solar-array panels, 2.2 m HGA disc; SA and HGA C-kernels, held after 10:39:19.851) seen
  from the only documented ROSINA location, ROS_ROSINA (-226300) = (0.8928, −0.7447, 2.6058) m in ROS_SPACECRAFT, which
  lies 4.9 cm inside the SPICE bus box and was lifted to 1 mm above the +Z face (the COPS boom geometry is undocumented).
* Fraction of the gauge's front half-space blocked by the spacecraft: RG (boresight +Z, toward the comet) 0.000 at
  every row (the wings and HGA lie behind the deck plane); NG (boresight −Y) 0.498, essentially constant (the bus
  behind the deck plane).
* Share of the modelled signal (c9 bounded softplus, the fit's weights; late window with W10x5) that falls in those
  blocked directions: RG ≤ 1.7 % (mean 0.03 %, maximum at 5.5 m); NG 0 before the gap; in the late window up to
  **45.3 % at 5.5 m** (mean 2.7 %), starting below ~37 m. Separately, NG's back-side shield part (in the fit) is up to
  38.5 % of the late NG signal (mean 3.0 %, below ~50 m) and 18.5 % mean / 1.24 % median before the gap.
* Consequence: in the last ~40 m the modelled NG pressure may be overestimated if the bus shields NG as this geometry
  suggests; the size depends on the real NG position (boom), which is unknown here. RG is not affected.
* Evidence: jobs 124942 (scrows), 124944 (`fov_share.py --hemi-sc`),
  `figures_cluster/12M/movies_softplus_smm1_9cycles_detach300_pan_hemi/validation/hemi_share/`, `rosina_point.c/.log`.

### A10. Most of the modelled NG signal comes from the big lobe; near the end the local east rise takes over (2026-11-13)

* Bounded softplus fit, c9 final state; per-facet NG contributions as the fit computes them (frozen nonneg_contrib
  path; the 95 late sums match the contribution data to 2.6e-14). Lobes split by a plane through the neck (narrowest
  cross-section): head 18.5 km², body 34.9 km². **The touchdown point (lon 12.16°, lat 31.69°) is on the head (small
  lobe), 1018 m from the neck plane.**
* Mean per-observation share of the modelled NG surface signal: all 934 NG observations: **body 76.2 %**, head 23.8 %
  (72 % from 1–3 km away); 95 late observations (W10x5, 297 → 5.5 m): body 46.4 %, head 53.6 % (44 % from 30–100 m).
  The body's share per late observation falls with altitude: 89 % at 297 m, 50 % at 141 m, 12 % at 101 m, 0 below
  ~88 m.
* Dominant regions: one 100 m cell 79 m east of the touchdown point (the rise east of the impact, 25 m higher;
  lon 14.7°, lat 32.1°) gives 44.9 % of the late NG signal; on the body, 1.67 km away at azimuth 105° (lon 66.6°,
  lat 12.6°) 18.7 % of the late and 11.6 % of the all-observation signal, and 2.67 km away at azimuth 46° (lon 149°,
  lat 55°) 12.6 % of the all-observation signal.
* Topography: the impact site lies in a depression; seen from 2 m above the touchdown point the highest skyline is the
  east rise (24.7° elevation at azimuth 91°, 37–39 m high within 78–82 m); this rise hides the body from NG below
  ~88 m. Higher but farther walls: 79–80 m at 280 m to the north-west, 106–127 m at 450–860 m to the south and
  south-west.
* RG near the end (2026-11-14, exact kernels at the movie frame times, job 124972): below ~8 m the modelled RG surface
  flux comes from a south-facing slope ~17 m north of the touchdown point. Extrapolated past the last observation along
  the SPICE trajectory, the modelled RG flux drops from 5.2e15 to ~1e-18 m⁻² s⁻¹ between 2.9 m and 2.1 m, when RG sees
  only facets with (numerically) zero activity — model-only, beyond the data.
* The fitted active regions have sharp edges (softplus with s = 1e-5 A_orig switches from active to ~0 within a
  narrow band of L), so contribution maps show crisp-edged patches, e.g. a near-circular active patch on the rise east
  of the touchdown point (visible from ~300 m to 40 m); these edges are activity edges, not instrument fields of view.
* Caveat: model-dependent (bounded softplus c9); facet slopes on the 12M mesh are not used as cliff evidence (55–83 %
  of facets are steeper than 60°: surface roughness).
* Evidence: job 124913 (`figures_cluster/scripts/ng_source_map.py`, `ng_source_regions.py`,
  `ng_lobe_share_vs_altitude.py`); `figures_cluster/12M/movies_softplus_smm1_9cycles_detach300_pan/analysis/`
  (`ng_source_regions_c9_map.png`, `ng_lobe_share_vs_altitude_c9.png/.tsv`).

### A9. Late in the descent the modelled RG signal comes mostly from outside the RG IK field of view (2026-11-13)

* The fit integrates RG over the whole front half-space of its boresight (RGAcceptance = hemisphere) and NG with the
  srccg_shield acceptance (B). Share of the modelled surface signal that comes from directions OUTSIDE the SPICE IK
  fields of view (ROS_ROSINA_V12.TI: RG -226332 rectangle 45° × 45° half-angles, NG -226331 77.5° × 77.5°), bounded
  softplus fit (A8), final state:

  | gauge | before the NG gap (to 10:24:38) | in the NG gap (700–297 m) | 95 late observations |
  |---|---|---|---|
  | RG | median 0 % (max 58 %) | median 78 % (58–85 %) | median 60 % (16–100 %) |
  | NG | median 17 % (100 % above ~5 km) | — | median 0.1 % (0–59 %) |

  RG by altitude: 0 % above 2 km, 21 % at 1–2 km, 53 % at 100–297 m, 80 % at 30–100 m, 91 % at 10–30 m, 99 % below
  10 m. Above 7.3 km the nucleus lies entirely outside the NG IK FOV.
* **Update 2026-11-13 (c9 state, circular fields of view):** the RG entrance is a round 6 mm opening into a 60 mm
  equilibrium sphere (Pestoni et al. 2021, A&A, arXiv:2012.01495), so its acceptance is axisymmetric; the IK rectangle
  is only SPICE's encoding (no angular width is documented besides the IK's 90° / 155° full widths). With circular
  cones (RG 45°, NG 77.5° half-angle) the shares outside are: RG before the gap median 0 % (up to 87 %, mean 20 %), in
  the NG gap median 88 % (75–90 %), late window median 62 % (18–100 %); NG before the gap median 17 % (0.6–100 %),
  late median 0.03 %. Evidence: job 124885,
  `figures_cluster/12M/movies_softplus_smm1_9cycles_detach200/validation/fov_share/`. A thin round opening by itself
  accepts the whole front hemisphere with a cosine weight, which is what the fit assumes; the 45° cone is the IK's
  nominal width only.
* **Consequence:** below ~1 km the modelled RG pressure depends on the assumed RG angular acceptance; with an IK-FOV
  acceptance most of the late RG signal of this model would be lost. Which acceptance is right is open (C). The IK
  itself is inconsistent: the same file lists FOV_ANGULAR_SIZE = π rad (180°) for RG and 5.41 rad (310°) for NG.
* Evidence: `figures_cluster/scripts/fov_share.py` (jobs 124731 etc.),
  `figures_cluster/12M/movies_softplus_smm1_5cycles/validation/fov_share/fov_share_smm1_c5_summary.txt`; mask check vs
  getfov 0 mismatches (8 frames × 336,769 directions per gauge).

### A8. With all slope multipliers ≥ 1 the fit is better: switching off 40–60° and 70–90° was not required by the data (2026-11-13)

* **Update: continued to 9 cycles (2026-11-13, job 124772; user OK).** c6–c9 fit / fine tuning: 4.2774e-3 /
  4.2750e-3, 4.2736e-3 / 4.2734e-3, 4.2694e-3 / 4.2694e-3, 4.2688e-3 / **4.2688e-3**; every fit stage max_nfev 4000
  (not converged), every fine tuning ftol; gain per cycle 0.26 → 0.04 → 0.10 → 0.01 %. Final c9: rms ln NG 0.0690, RG
  0.0614; −0.40 % vs c5, −2.2 % vs the unbounded fit, −12.6 % vs D2. Multipliers c9 (c5): 40–50° 33.5 (27.2), 60–70°
  14.1 (15.6), 80–90° 17.0 (10.5), 90–180° 5.94 (5.62); 50–60° and 70–80° still at the bound 1; κ and grain B at 0. The
  40–50° and 80–90° multipliers are still drifting, so their values are not determined; the pattern (40–50° and 80–90°
  up, 60–70° and overhangs moderate, 50–60° and 70–80° at the bound) is stable from c5 to c9. A = 0 on 97.4 % of the
  seen area. P50 source areas c9 vs c5: area Jaccard 0.83 (NG) / 0.82 (RG), ~98.5 % of the c9 union inside the c5
  union (15–16 % smaller): stable; vs the unbounded fit 0.11 / 0.09. Slope ≥ 60° share of the modelled pressure:
  before the late window NG 22 %, RG 21 %; late NG 27 %, RG 23 %; 40–50° gives 58 % before the late window. Overhangs
  below 20 m NG 68 %, RG 70 %; at 5.5 m NG 75 %, RG 47 %. Q_H2O at 10:33:10 5.38e26 s^-1. Model NG = RG crossing
  472.6 m. Evidence: `figures_cluster/12M/model_vs_observations_softplus_smm1_9cycles/chain_stages.txt`,
  `full_set_softplus_smm1_9cycles/comparison_vs_c5/`, `comparison_vs_unbounded/`; maps 124775, contributions 124776
  (`<chain>/maps_c9`, `contributions_c9`). The c5 numbers below are kept as the 5-cycle state.
* **Setup** (user 2026-11-12: the multipliers may not decrease below 1, start 1, may grow): softplus × slope
  multipliers exactly as the 5-cycle run of A3/A6 (start D2, fine tuning late_only 5 min, 5 cycles), but the six free
  multipliers have lower bound 1 (was 0). Option `nonneg_fit.py --slope-mult-min 1`.
* **Cost:** final 4.2861e-3 (rms ln NG 0.0688, RG 0.0620): 1.8 % below the unbounded fit (4.3668e-3) and 12.3 % below
  D2 (4.886e-3). Per cycle (fit / fine tuning): c1 4.6610e-3 (max_nfev) / 4.4415e-3; c2 4.3872e-3 (max_nfev) /
  4.3820e-3; c3 4.3648e-3 (ftol) / 4.3608e-3; c4 4.3353e-3 (ftol) / 4.3211e-3; c5 4.3098e-3 (**max_nfev**) /
  4.2861e-3 (ftol). **Not converged** (the last fit stopped at max_nfev 4000).
* **So the unbounded fit's ftol stop with four multipliers at 1e-22 … 1e-31 (A3) was not the best solution of this
  model.** Keeping those bins ≥ 1 led the optimizer to a lower cost with a different pattern.
* **Final multipliers** (0–40° held at 1): 40–50° **27.2**, 50–60° 1 (at the bound), 60–70° **15.6**, 70–80° 1 (at the
  bound), 80–90° **10.5**, 90–180° **5.62**; κ and the grain parameter B at their bound 0 (κ 2.7e-37). The values are
  not settled: across cycles 40–50° went 1 → 1.29 → 4.25 → 17.6 → 27.2 and 80–90° 1 → 1 → 1 → 6.15 → 10.5.
* **Activity:** A = 0 on 97.2 % of the seen area (unbounded 96.5 %), 95.7 % of the facets, 95.1 % of the area.
* **Production:** Q_H2O at 10:33:10 = 4.75e26 s^-1 (unbounded 2.41e26). The 40–50° bin gives 28 % of the seen
  production and 55 % of the whole-surface production; most of the total is on unseen surface and unconstrained (A4).
* **Which slopes carry the modelled pressure** (mean shares of the surface signal):

  | | 0–40° | 40–50° | 60–70° | ≥ 60° total | ≥ 60° unbounded fit |
  |---|---|---|---|---|---|
  | NG before the late window | 21 % | 54 % | 19 % | 24 % (15–35 %) | 54 % |
  | NG late window (95 obs) | 15 % | 54 % | 18 % | 29 % | 69 % |
  | RG late window (95 obs) | 19 % | 56 % | 14 % | 24 % | 61 % |

  The conclusion "cliffs ≥ 60° carry most of the signal" (A3) therefore depends on the bound: the data are fitted at
  least as well with the enhancement on 40–50° slopes.
* **Overhangs near the surface (A1 holds):** 90–180° share below 20 m NG 66 %, RG 64 % (unbounded 58 % / 54 %); at
  5.5 m NG 73 %, RG 45 % (67 % / 40 %); dark overhangs produce 0.2 × 5.62 = 1.12 × a sunlit flat facet.
* **The P50 source areas moved:** union area overlap (Jaccard) with the unbounded fit 0.16 (NG) / 0.14 (RG), P90
  0.43 / 0.42; NG P50 union 4,036 facets, 1.13e4 m² (unbounded 3,114), RG 3,770 facets, 1.03e4 m² (unbounded 4,413).
  The dense core is again 1–6° east of the impact (lon 13.3–17.0°, lat 31.5–32.9°). Specific facets are model
  dependent; the general region is not.
* Model NG = RG crossing at 475.5 m, inside the NG data gap (unbounded 496.8 m).
* Evidence: chain `runs/shap8_012m_lategrid_dz1.5/nonneg/softplus/chain_5cycles_smm1_20261113T041215Z` (job 124664,
  code `model_ext_frozen/nonneg_20261113T035423Z`, chain code `model_ext_frozen/chain_onejob_20261113T041215Z`), maps
  124667, contributions 124668; `figures_cluster/12M/model_vs_observations_softplus_smm1/summary.txt`,
  `chain_stages.txt`; `figures_cluster/12M/full_set_softplus_smm1_5cycles/comparison_vs_unbounded/` (job 124727).

### A1. Overhangs contribute strongly only in the last ~10 m, through the night-side production floor (2026-11-12)

Question (user): an effective slope > 90° means an overhang, mostly in shadow, so it should contribute little; why is
it a significant part close to the surface?
* **Overhangs are mostly, not always, dark.** The effective gravitational slope is the angle between the inward facet
  normal and the local effective gravity (self-gravity at 533 kg/m³ + centrifugal term); > 90° = the facet faces
  partly downward relative to gravity, not away from the Sun. Seen overhang facets at 10:33:10 UTC: 27 % sunlit, 54 %
  night side (cos SZA ≤ 0), 19 % in cast shadow (flat 0–40° facets: 61 % / 26 % / 14 %).
* **Dark facets still produce in the model.** The production law kept from srcCG (as in D2) sets the night-side and
  cast-shadow rate to q_min = 0.2 q_max for H2O. With the fitted overhang multiplier 4.42, a dark overhang produces
  0.2 × 4.42 = 0.88 of a fully sunlit (cos SZA = 1) flat facet with the same harmonic activity.
* **Only the last seconds are dominated by overhangs.** Overhang share of the model's facet signal after the gap:
  1–2 % for most of the 95 post-gap observations (NG 0.024 at 10:33:10); NG 0.672 at 10:39:26 (all night side), RG
  maximum 0.608 at 10:39:18 and 0.397 at 10:39:26 (night side); summed over the 95 observations NG 15 %, RG 6 %.
  At 10:39:26 the remaining flat-facet RG signal is mostly from cast-shadowed facets (0.544).
* **Mechanism:** with the gauge a few metres above the surface next to the cliff at the impact site, the nearest
  facets dominate (contribution ∝ facet area / distance², times the gauge's angular response), and those are dark
  overhangs producing at the night-side floor times the fitted multiplier.
* **Correction 2026-11-18 (A18):** the 0.2 floor was not intended; srcCG and the fixed code use 0.02, which makes dark
  overhangs produce 10x less at the same A. To be re-evaluated after a refit with the 0.02 law.
* **Caveat for interpretation:** this is a consequence of the 20 % night/shadow floor plus a free overhang multiplier;
  the fit uses the overhang multiplier mainly to reproduce the last seconds of NG and RG. Not yet tested: a fit with
  the overhang multiplier held at 1 or 0, or without a night/shadow floor for overhangs.
* Evidence: softplus 5-cycle final state (chain `runs/shap8_012m_lategrid_dz1.5/nonneg/softplus/chain_5cycles_20261112T012553Z`,
  job 124310); facet contributions (job 124549, `<chain>/contributions/`); illumination split
  `analysis/overhang_illumination_12m/out_124578/summary.txt` (CPU job 124578); illumination from the maps geometry
  (job 124324).

### A2. Where the post-gap signal comes from: P% source areas (2026-11-12)

* **Definition:** the P% source area S_P(k, g) of observation k and gauge g is the smallest set of facets, taken in
  decreasing order of their contribution to the modelled surface signal, that supplies P % of it (after the footprint
  "source area" of Schmid 1994); its source-area occurrence N_P(e) is the number of the 95 post-gap observations
  (10:33:10–10:39:26 UTC) whose S_P contains facet e. All facet contributions are ≥ 0; they sum to the model values to
  ≤ 4.7e-10 (print precision).
* **The post-gap NG and RG signal comes from a compact area east of the impact point,** about lon 13.5–17°, lat
  31–33°, around a cliff near lon 15.3°; RG's area is wider than NG's.
* Union over the 95 observations: P50 NG 3114 facets (9.1e3 m², 0.017 % of the nucleus), RG 4413 facets (1.19e4 m²);
  P90 NG 37,228 facets (8.1e4 m², 0.15 %), RG 37,256 (7.7e4 m²). Median facets per observation: P50 NG 1515, RG 1580;
  P90 NG 9211, RG 9399. S_50 shrinks from ~3000 (NG) / ~4000 (RG) facets at 10:33:10 to 1 (NG) / 3 (RG) facets at
  10:39:26 (gauge ~5 m above the surface).
* **NG dust-grain term** (not attributed to facets): 29–75 % of the modelled NG pressure after the gap (median 71 %).
* **Correction 2026-11-13 (the "compact area" above is incomplete):** by AREA the P50 unions are not compact. Measured
  from the last gauge position, only ~10 % of the union area lies within 500 m (unbounded NG: 20–100 m 7 %, 100–500 m
  3 %, 0.5–2 km 49 %, > 2 km 40 %); the 5–95 % longitude range is 15–156°. The region east of the impact point is the
  dense core seen in the close-up occurrence maps; the rest is spread over distant facets. The bounded fit (A8) gives
  the same picture (> 2 km: 50 % NG, 46 % RG). Source: `full_set_softplus_smm1_5cycles/comparison_vs_unbounded/
  source_areas_and_slope_shares.txt`.
* Evidence: job 124549 (`<chain>/contributions/README.txt`, `source_areas_{NG,RG}.npz`); figures
  `figures_cluster/12M/full_set_softplus_5cycles_v2/source_areas/`.

### A3. Steep slopes and overhangs carry the fitted production; 40–60° and 70–90° are switched off (2026-11-10 … 12)

* **Correction 2026-11-13:** the switch-off of 40–60° and 70–90° is NOT required by the data. With all multipliers ≥ 1
  the fit is 1.8 % better, with 40–50° ×27 and 80–90° ×10 and a ≥ 60° share of the signal of 24–29 % instead of
  54–69 % (A8). The numbers below describe the unbounded fit only.
* Softplus × slope multipliers, 5 cycles, final multipliers: 0–40° 1 (held), 40–50° 1.3e-22, 50–60° 1.1e-25,
  **60–70° 27.8**, 70–80° 4.3e-29, 80–90° 1.8e-31, **90–180° 4.42**. Area shares of the nucleus: 68.3, 9.6, 7.7, 6.1,
  4.1, 2.3, 2.0 %.
* Instantaneous H2O production of the seen surface at 10:33:10 / 10:39:26: 60–70° 0.440 / 0.448, overhangs 0.400 /
  0.395, 0–40° 0.160 / 0.158, the other bins ~0.
* Share of the modelled surface pressure: 60–70° facets 34–71 % of the RG signal before the gap; overhangs ≤ 2 %
  before the gap and up to ~60–70 % in the last ~10 m (A1).
* Earlier, clipped refit D2 (2026-11-10): cliffs ≥ 60° are 14.6 % of the seen area but carry 44.7 % of D2's seen
  production; the 45–60° band ×3 activity in the clipped fits, ×~1 in the positive fits; the cliff coefficient is
  significant in every fit.
* The 3-cycle and 5-cycle squared-harmonics fits show the same pattern (50–60° ×16.7, overhangs ×3.5, others ≤ 8e-4).
* Evidence: `<chain>/c5/fit/final/nonneg_bins.txt`, contribution job 124549 (`slope_bins.txt`); figures
  `full_set_softplus_5cycles_v2/cliffs/`.

### A4. Most of the seen surface has exactly zero activity in the best fits (2026-11-10 … 12)

* Clipped refit D2: activity exactly 0 on 95.2 % of the seen facets (the linear field before the clip is negative on
  94.7 % of the seen area; its median −1.1e20).
* Softplus (s = 1e-5 A_orig, start from D2): A = 0 on 96.4 % of the seen area (softplus underflow; the start is
  numerically D2's clip). The non-zero activity forms a few smooth patches (two near the south pole, an elongated one
  around lon −45°, small ones near 60° and 150° E, one just east of the impact point).
* The zero-activity surface of D2 is NOT mainly dark: never sunlit while seen = 0.6 % of its area; 95.8 % (NG) and
  96.4 % (RG) of its kernel weight is sunlit. The fit switches off sunlit surface, not night or shadow.
* The unseen surface carries ~92 % of the softplus fit's total production (8.4 × the original srcCG map's total);
  it is not constrained by the data.
* Evidence: `analysis/shadow_zero_activity_12m/` (job 123614), `nonneg_bins.txt`, maps job 124324.

### A5. No CO2 is needed by the best fits (2026-11-09 … 12)

* The CO2 scale κ (on the fixed srcCG CO2 map ACO2) goes to its lower bound 0 in D2 (2e-36), softplus (3.7e-37),
  squared (9.7e-35), bumps (7e-17) and is 9e-6 in the best positive fit; only the slope multipliers × harmonics fit
  (floor 1e-5) kept κ = 0.466.
* Softplus final state: Q_CO2 ≈ 1.4e-12 s^-1, Q_CO2/Q_H2O ≈ 5.9e-39. Bounded softplus (A8, 2026-11-13): κ 2.7e-37,
  CO2/H2O ≈ 2.2e-39.
* Caveat: COPS measures total density / ram pressure and does not separate species; "no CO2" means the CO2 term is
  not needed to fit the pressures with this CO2 map and law.

### A6. Model comparison on the same data and base (2026-11-08 … 12)

Data cost (lower is better; all on the lategrid base with the gauge +1.5 m and the W10x5 time response):

| model | data cost | rms ln NG / RG | stop |
|---|---|---|---|
| softplus × slope multipliers ≥ 1, 9 cycles (A8, 2026-11-13) | **4.269e-3** | 0.069 / 0.061 | every fit max_nfev 4000 (not converged); gain 0.01 %/cycle at c9 |
| softplus × slope multipliers ≥ 1, 5 cycles (A8, 2026-11-13) | 4.286e-3 | 0.069 / 0.062 | last fit max_nfev 4000 (not converged) |
| softplus × slope multipliers ≥ 0, 5 cycles | 4.367e-3 | 0.071 / 0.061 | ftol (a lower cost exists: A8) |
| softplus, 3 cycles | 4.370e-3 | 0.071 / 0.061 | ftol |
| clipped refit D2 (cliff + 6 ft templates) | 4.886e-3 | 0.078 / 0.061 | max_nfev 4000 |
| squared harmonics × slope multipliers, 5 cycles | 6.649e-3 | 0.088 / 0.074 | max_nfev 4000 (not converged) |
| positive (exp) activity, ratio 1e-5, λ 1e-4, 3 cycles | 8.537e-3 | 0.107 / 0.075 | ftol |
| slope multipliers × harmonics, floor 1e-5 | 8.723e-3 | 0.105 / 0.080 | max_nfev 4000 |
| Gaussian bumps × slope multipliers, 3 cycles | 1.413e-2 | 0.119 / 0.118 | ftol |

* Models that allow (near-)zero activity over most of the seen surface fit best; strictly positive (exp) activity
  fits 1.7–2.2 × worse.
* Gauge treatment: the clipped c6 model refitted without the late-window time response and gauge offset (D0, flat
  base) reached 4.759e-3; with them (D2) 4.886e-3: the response/offset choice changes the cost by 3–4 %.

### A7. Observational facts used in the figures (2026-11-12)

* **Gap in the NG data:** 700 → 297 m altitude, 10:24:38 → 10:33:10 UTC; RG kept measuring.
* The model's NG = RG crossing (softplus 5 cycles) is at 496.8 m (bounded softplus, A8: 475.5 m), inside the NG gap;
  the data crossing must lie in the gap as well.
* **Rosetta's own shadow** (2026-11-13, SPICE geometry, used in the movies only, not in any fit): in the late window
  the Sun is 30° from the anti-nadir (Sun in ROS_SPACECRAFT ≈ (0.495, 0, −0.869)), so the spacecraft shadow falls on the
  surface ~30° off the sub-spacecraft point, ~29.7° from the RG axis (inside the RG IK FOV) and ~90° from the NG axis.
  Bus 2.10 × 2.00 × 2.58 m and two solar-array wings 14.36 × 2.25 m (32.6 m span) from ROS_STRUCT_V8.BSP;
  at r_h ≈ 3.8 AU the umbra of a 2-m body is ~820 m long. In the movie geometry the shadow first reaches the surface at
  1330 m altitude and has an umbra from 652 m on. Unknown: the offset between the structure origin (-226000) and the
  trajectory point (-226), ~1 m; no loaded kernel gives it (0 used).
* Late window: 10:33:10–10:39:26 UTC (95 observations), gauge altitude 297–5.5 m (spacecraft reference point
  298–7.0 m). Rosetta impact (extrapolated from the last two observations): 10:39:33.8 UTC, lon 12.164°, lat 31.685°
  (facet 11,194,754 of the 12M mesh); check 2026-11-13: the SPICE RORB trajectory of ROSETTA crosses the 12M surface on
  the same facet, 7.4 mm from this point, at 10:39:33.874 UTC (the COPS track, +1.5 m, 8.7 mm away, 1.70 s earlier).
  With the SPICE spacecraft model (bus box, solar-array wings, 2.2 m HGA; SA and HGA orientation held after
  10:39:19.851) the first contact with the 12M surface is the −Y solar-array wing at 10:39:30.624 UTC (±1 cm), i.e.
  ~3.2 s before the centre of mass reaches the surface; the bus +Z face ~10:39:30.8.
* None of the models reproduces the NG dip around 10:35–10:37.

## B. Method facts for the publication's methods section
* Gauge acceptance as used in the fit (LiouvilleCOPS src/Gauges.cpp, Liouville.h, Geometry.cpp; base settings
  NGAcceptance = srccg_shield, NGHeadOnWeight = on, RGAcceptance = hemisphere): RG accepts the full front half-space of
  its boresight (+Z of ROS_ROSINA_COPS = spacecraft +Z, toward the comet), 2 pi sr, flux weight -l.b_RG; NG (boresight
  -Y) accepts 2.95 pi sr (74 % of the sphere): front half-space weight 1 - l.b_NG, the back half-space only through the
  shield (fraction 1.00 at 90 deg from the axis, 0.65 at 60, 0.22 at 37, 0 within 17.5 deg of +Y). Both gauges therefore
  receive gas from essentially the whole visible nucleus; the IK fields of view (RG 45 deg, NG 77.5 deg) are not used.
  Over the descent the RG boresight is 2-27 deg from the nucleus centre; the NG boresight 66-94 deg (the nucleus
  straddles NG's front/back boundary: 55-81 % in front). Evidence: analysis/movie_geometry_check/ (2026-11-12).

* Collisionless (Liouville) transport from the 12M SHAP8 shape-model surface (12,364,290 facets) to the gauges;
  per-facet transport kernels Kn (NG), KΦ and KV (RG) (LiouvilleCOPS src/Liouville.h), cos SZA and cast shadow per
  record; 1551 base rows (1150 observations + 401 model-only 1-s grid rows in the late window).
* Production law (srcCG): q = q_min (night or cast shadow) else q_min + (q_max − q_min) cos SZA; H2O q_max =
  600 / r_h^p (r_h in AU, p = 4.2143229), q_min = 0.02 q_max as srcCG (since 2026-11-18; EVERY fit up to 2026-11-18
  used 0.2 q_max, not intended: A18); CO2 q_max = 600 / r_h², q_min = 0.1 q_max; flux
  F = q × A [m^-2 s^-1].
* Gauge time response in the late window (W10x5): mean of five 10-s running averages ending at T, T−2, …, T−8 s,
  T = the L2 STOP time; gauge 1.5 m along +Z from the trajectory point; dust-grain term in NG; NG = RG crossing
  constraint (weight 10).
* Effective gravitational slope: polyhedral self-gravity at 533 kg/m³ plus the centrifugal term (SPICE rotation).
* Fits: TRF least squares, max_nfev 4000; per-facet fine tuning on the facets seen only in the last 5 min;
  3–5 fit / fine-tuning cycles.
* Option C fine tuning on the CPU (2026-11-15): the GPU fine-tuning code run unchanged with a CPU engine on precomputed
  late kernels; on the softplus chain it selects the same 271 facets and reproduces the start costs to 8e-16 and a
  4-evaluation case bitwise; longer runs differ by up to ~1e-8 in the final cost, the same spread as the CPU code
  itself under 1e-15 start perturbations (truncated-LSMR trust-region steps amplify rounding on the flat end).

## C. Open questions (physics)

* (2026-11-18, A18) How do the fits change with the intended q_min = 0.02 q_max: costs, ranking of the models, the late
  NG/RG, where the activity goes (seen vs unseen, steep vs flat)? Needs refits (and, first, bases whose COPS_model.dat
  is consistent with the new law: tests/nightfix_diagnostic_20261118/README.txt).
* Should dark overhangs produce at the night-side floor times a free multiplier (A1)? Tests: overhang multiplier
  held at 1 or 0; no floor for overhangs; merge the 90–180° bin into 80–90°.
* Is the switch-off of sunlit surface (A4) physical (dust mantle) or a basis limitation?
* The NG dip at 10:35–10:37 is not reproduced by any model.
* The slope multipliers are not determined by the data (A8: still changing between cycles; the bound changes which
  slopes carry the signal). Which slope pattern, if any, is physical?
* RG angular acceptance: hemisphere (as fitted) or the IK field of view? Late in the descent 53–99 % of the modelled RG
  signal comes from outside the IK FOV (A9).
* Does Rosetta's own shadow (A7) change the late NG/RG? Not estimated (user 2026-11-13: movies only).
* Spacecraft occlusion of NG in the last ~40 m (A11): where exactly is the NG (boom tip) relative to the bus, and
  should the fit block directions through the spacecraft?
* Option C (A12): is the steep / flat activity contrast physical, or the only way a 150 m basis can make local sources?
  The multipliers trade off with the global level of the c-field; a penalty on the level (or m held) would decide it.

## D. Change log of this file

* 2026-11-21 (21:10 CST): A20: correction of A19 (switched-off surface, ranking, cliffs); visibility by gauge; activity vs slope with illumination removed.
* 2026-11-21 (21:30 CST): A21: switched-off surface by the user's decision 1b; cliff wording of A20, A19 and the overhang wording of A18 corrected.
* 2026-11-21 (17:50 CST): A19 addendum: sources of the five most plausible fits (robust vs model-dependent).
* 2026-11-21 (15:15 CST): A19 addendum: why the templates-off C chain did not move; C penalty shares; parameter counts. A19: C f1 share above F_E 1.6 % (was 1.7 %, rounding).
* 2026-11-21 (14:40 CST): A19: plausibility of the ten new-law fits (energy limit, active area, seen share, ranking).
* 2026-11-20 (16:00 CST): A18: the full figure sets of the ten new-law states (lobe, shell, P50 and overhang shares vs the old law).
* 2026-11-19 (17:30 CST): A18: the slope multipliers of the new-law fits vs the old law (comparison set stage2_q002).
* 2026-11-19 (14:40 CST): A18 addendum: option A with the new law and the campaign summary of all ten new-law states.
* 2026-11-19 (01:25 CST): A18 addendum: option B with the new law (5.646e-3; Q_H2O 4.28e25, 25 % seen).
* 2026-11-18 (23:10 CST): A18 addendum: source rates of the three option C chains with P_m 0.01 (total Q ~90 % seen).
* 2026-11-18 (22:40 CST): A18 addendum: the three option C chains with P_m 0.01 under the new law (+12 / +24 / +55 %).
* 2026-11-18 (21:30 CST): A18 addendum: option C f1 / P_m 1e-4 and option D with the new law.
* 2026-11-18 (19:45 CST): A18 addendum: option C f10 / P_m 1e-4 with the new law (2.4988e-3, +1.6 % vs old law;
  Q_H2O 4.84e26 set by the undetermined steep multipliers).
* 2026-11-18 (18:10 CST): A18 addendum: the first new-law stage-2 chain (softplus regularized): costs per cycle and its
  source rates (Q_H2O 6.08e25 at 10:33:10; dark 76.6 % of Q, 18.2 % of the seen Q).
* 2026-11-18 (14:45 CST): A18 addendum: where the D2 refit's production comes from (seen 27.6 % of Q; dark 68.8 % of Q
  but 6.5 % of the seen Q; source_rate_q002/D2_q002_warm).
* 2026-11-18 (13:35 CST): A18 addendum: Q_H2O of the D2 refit (3.12e25 at 10:33:10 vs D2's 2.11e26).
* 2026-11-18 (04:40 CST): A18 addendum: the D2 refit with the new law on base_q002 (data cost 5.671e-3 vs 4.886e-3).
* 2026-11-18 (01:45 CST): A18 correction notice (all results used H2O q_min = 0.2 q_max; intended 0.02) with the
  fixed-A diagnostic of the eight final states; correction notes in A1 and A17; B production law; C open question.
* 2026-11-18 (00:40 CST): A17 (night / cast-shadow factor of the production: present, up to x5, hidden by the A range; most of
  the production is dark at the map times; illumination identical in the six states).
* 2026-11-18 (00:10 CST): A16 (slope at equal illumination: SZA-controlled source-rate ratios of the six final states).
* 2026-11-18 (02:35 CST): A12: option C design matrices of the 0.02 night law.
* 2026-11-17 (14:50 CST): A12: the regularized f_m = 1e-4 chains (f = 10, f = 1).
* 2026-11-17 (18:30 CST): A15 (full figure sets of the six final states: the seen surface gives ~6 % of Q_H2O in every
  model; robust and model-dependent quantities).
* 2026-11-17 (19:00 CST): A15 corrected / extended (lead): the ~6 % seen share holds for the six full-set states, not for
  the option C chains with f_m 0.01 (74-95 % seen); the seen production at 10:33:10 is about 1.2-3.2e25 s^-1 in every model
  while Q_H2O varies by > 100x.
* 2026-11-17 (15:00 CST): A14 (final stage-2 comparison: costs reg. vs unreg., late residuals, multipliers,
  convergence, kappa of all 59 fits); A12 kappa bullet widened to every fit (with the ninth C state).
* 2026-11-17 (14:05 CST): A12: kappa = 0 in all eight C stage-2 states.
* 2026-11-16 (00:05 CST): A12 update: stage-2 reruns, regularized fine tuning, templates on vs off.
* 2026-11-15 (21:40 CST): A12 update: fine-tuning regularization test (C, B, D c5).
* 2026-11-15 (14:15 CST): A12 update: the complete stage-2 chain f = 10, f_m = 0.01.
* 2026-11-15 (13:20 CST): A12 update (stage 2: the first option C fine tuning); B: CPU fine tuning of option C.
* 2026-11-15: A13 (stage-1 comparison of the four non-negative options: active area vs fit quality).
* 2026-11-14 (20:50 CST): A12 update: the multiplier-penalty continuation (125360).
* 2026-11-14 (16:40 CST): A12 (option C local basis: better fits, steep-facet activity, undetermined multipliers;
  preliminary, f = 1 pending); C new question.

* 2026-11-13 (late night): A11 (spacecraft occlusion: RG none, NG up to 45 % of the late signal at 5.5 m; shield
  part); C new question.
* 2026-11-13 (night): A10 (where the NG signal comes from: big lobe 76 % over the descent, east rise near the end).
* 2026-11-13 (evening): A9 update: RG acceptance axisymmetric (round 6 mm opening); fov_share with circular cones;
  A7 contact times of the spacecraft model with the 12M surface on the RORB trajectory (first contact −Y wing
  10:39:30.624 UTC).
* 2026-11-13 (later): A8 update with the 9-cycle continuation (c9 state; source areas stable vs c5); A6 row; A7 the
  touchdown point check: the RORB trajectory meets the 12M surface 7.4 mm from the ground-track impact point
  (10:39:33.874 UTC, facet 11,194,754).
* 2026-11-13: A8 (bounded softplus fit, multipliers ≥ 1: better cost, different slope pattern, source areas moved),
  A9 (share of the modelled signal outside the IK FOVs); corrections in A2 (P50 unions are not compact by area) and A3
  (the switch-off is not required); A1, A5, A6, A7 updated with the bounded fit; A7 Rosetta's shadow geometry; C three
  new questions.

* 2026-11-12: created; A1 (overhangs) and the findings of 2026-11-08 … 12 (A2–A7), methods (B), open questions (C).
* 2026-11-12: B: gauge acceptance as used in the fit (RG 2 pi sr, NG 2.95 pi sr; IK FOVs unused).
