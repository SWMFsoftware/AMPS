# FINDINGS.md — physics results of the 12M COPS fits, for use in publications

Kept current (user request 2026-11-12: "create FINDINGS.md and save that there. also add there new important physics
related notes in the future for using that in the publications; keep the file current"). Each entry: date, the
finding with its numbers, the evidence (job ids, files), and caveats. Newest first within each section. Paths are
relative to `COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output` unless absolute. Figures: see
`OUTPUT_VISUALIZATION.md`.

---

## A. Findings

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
  on 96.4 % / 94.7 % of the seen area in softplus / D2 - or active everywhere - A, B, D: A = 0 on <= 0.9 %.)
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
  600 / r_h^p (r_h in AU, p = 4.2143229), q_min = 0.2 q_max; CO2 q_max = 600 / r_h², q_min = 0.1 q_max; flux
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
