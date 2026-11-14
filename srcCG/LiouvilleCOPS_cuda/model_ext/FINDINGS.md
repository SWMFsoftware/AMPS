# FINDINGS.md — physics results of the 12M COPS fits, for use in publications

Kept current (user request 2026-11-12: "create FINDINGS.md and save that there. also add there new important physics
related notes in the future for using that in the publications; keep the file current"). Each entry: date, the
finding with its numbers, the evidence (job ids, files), and caveats. Newest first within each section. Paths are
relative to `COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output` unless absolute. Figures: see
`OUTPUT_VISUALIZATION.md`.

---

## A. Findings

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

## D. Change log of this file

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
