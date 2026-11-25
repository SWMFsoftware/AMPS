# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/q002_top5_content.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/q002_top5_content.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: f5e745befe880ef32e6d64469a4501978328768ec939eff1e07bc379765c850e
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: text of 12M_q002_top5_models.docx: summary, tables, captions, section 4 (ground track, source rates, areas,
#   slopes and cliffs, robust vs model-dependent), the methods sections with their sources (2026-11-21)
# RUN: imported by q002_top5_docx.py
# CHANGED: 2026-11-21 user decisions 1b / 2b (FINDINGS A21)
# ==== end of header ====
"""Text of 12M_q002_top5_models.docx (2026-11-21): used by q002_top5_docx.py.

Every statement below is taken from the code that ran or from the run files. Main sources (cited in the text):
  fit_cops.py (FC) = /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS/tools/fit_cops.py (sha256 = the runs' fit_manifest)
  LiouvilleCOPS src (Source.cpp, Gauges.cpp, Liouville.cpp) and README.md; base_q002/input.snapshot.in
  frozen option A/B code: cops_output/model_ext_frozen/nonneg_20261118T072729Z/ (nonneg_alt.py, nonneg_mult.py,
    nonneg_problem.py, nonneg_finetune.py, positive_problem.py, positive_finetune.py, grain_crossing.py, slope_mult.py)
  frozen option C code: cops_output/model_ext_frozen/localnn_q002_20261118T162327Z/ (localnn_problem.py, localnn_fit.py,
    localnn_mesh.py, localnn_design.py, localnn_late.py, localnn_ft.py, localnn_chain.sbatch.copy)
  run folders under cops_output/runs/shap8_012m_lategrid_dz1.5/nonneg/ (README.md, commands.sh, fit_result.txt,
    finetune_result.txt, c5/analysis/report.txt), the comparison set summary.csv, q002_plausibility/metrics.json,
    FINDINGS.md (A12, A18-A21, B); decisions/active_area.json
Numbers of the five fits come from numbers.json (q002_top5_figures.py) where available; the others are quoted from the
run files named next to them.
"""

LAW = "q_min = 0.02 q_max"
TITLE = "Five selected fits of the new-law round (SHAP8 12M)"
SUBTITLE = ("Summary, summary plots, how the five activity models differ numerically, their assumptions, and what the "
            "fitting procedure minimizes. Rosetta/COPS end of mission, comet 67P, 2016-09-30; H2O night/shadow law "
            "q_min = 0.02 q_max. Prepared 2026-11-21; corrected 2026-11-22.")
FRONT = [
    "**Request (user, 2026-11-21):** \"prepare a summary word doc for the 5 models in the table; in addition to the "
    "summary plots, include a detailed description of how those model numerically are difference from each other, what "
    "are ther assumptions, what the fitting procedure is actuly minimizes\". The five = ranks 1-5 of the plausibility "
    "ranking as first made on 2026-11-21 (withdrawn, see the correction), the user's choice.",
    "**Corrections (2026-11-21).** (1) The first ranking counted only facets with activity exactly zero as switched "
    "off; it is withdrawn. (2) The measure used in the 21:00 CST version (activity below a fraction of the fit's own "
    "uniform-equivalent activity) was ad hoc. It is replaced by the project's existing threshold with the imposed "
    "floor not counted (user decision 1b): option B now fails, and A is the only one of the five that passes. (3) "
    "That version said that 60-80° cliffs are as active as flat terrain. This is withdrawn: the slope multipliers "
    "have an imposed lower bound of 1 and these two sit at or near it, so the fits could not represent less active "
    "cliffs. The opposite earlier statement, that active cliffs are not required by the data, is withdrawn too. (4) "
    "The gauge visibility of the late window is corrected (Section 4.1). Also corrected: the number of observations "
    "against which option C's parameters are counted (2084, not 1551) and the range of the dust-grain term.",
    "**Correction (2026-11-22): the C fit with templates off is one of five solutions of its own problem.** It stopped "
    "on ftol in a local minimum, not at the best solution. From other starts the same model puts 94-98 % of the seen "
    "water on 50-60°, 60-70°, 70-80° or 80-90°, with data costs of 4.78, 4.52, 3.91 and 3.76×10^{-3} (this fit: "
    "5.67×10^{-3}). Its slope class, multiplier and switched-off share describe that one solution, and the slope class "
    "is not measured. Earlier versions presented it as converged and 'a 40-60° class raised' as common to all five "
    "fits; as conclusions both are withdrawn (Section 10). The other fits were not tested for other solutions.",
    "**Folder:** /data001/highresolu_vt/cops_output/figures_cluster/12M/q002_top5_models/ (this document, figures/, "
    "numbers.json, build_manifest.json). Paths below are relative to /data001/highresolu_vt/cops_output unless absolute.",
    "**Order used throughout** (that of the withdrawn ranking, kept as labels only): (1) C f10, P_m 0.01, templates on; (2) C f1, P_m 0.01, templates on; "
    "(3) C f10, P_m 0.01, templates off; (4) A sum of squares, rank 3; (5) B background, bg 1e-2.",
]
ORDER_LABEL = {"C_f10_pm0p01_q002": "C f10 / P_m 0.01 / on", "C_f1_pm0p01_q002": "C f1 / P_m 0.01 / on",
               "C_f10_pm0p01_tploff_q002": "C f10 / P_m 0.01 / off", "A_q002": "A (rank 3)", "B_q002": "B (bg 1e-2)"}
# from the run files (sources in the methods text)
FREE = {"C_f10_pm0p01_q002": 2754, "C_f1_pm0p01_q002": 2754, "C_f10_pm0p01_tploff_q002": 2754, "A_q002": 81,
        "B_q002": 35}
CZERO = {"C_f10_pm0p01_q002": "1814 (66.1 %)", "C_f1_pm0p01_q002": "1867 (68.0 %)",
         "C_f10_pm0p01_tploff_q002": "1531 (55.8 %)"}
C_PEN = {"C_f10_pm0p01_q002": ("2.290e-4", "1.607e-3", "5.4228e-3"), "C_f1_pm0p01_q002": ("1.592e-4", "9.125e-4", "5.2448e-3"),
         "C_f10_pm0p01_tploff_q002": ("3.621e-4", "9.144e-4", "6.9499e-3")}
FIT_NFEV = {"C_f10_pm0p01_q002": "3859 (3719 + 140 after a wall-time restart), 40, 14, 8, 6",
            "C_f1_pm0p01_q002": "3861 (3798 + 63), 33, 9, 9, 7",
            "C_f10_pm0p01_tploff_q002": "3765 (3729 + 36), 3, 3, 3, 3",
            "A_q002": "4000 in every cycle (max_nfev)",
            "B_q002": "2723 (ftol), 4000, 4000, 2751 (ftol), 4000"}
MAXHD = {"C_f10_pm0p01_q002": "3.36", "C_f1_pm0p01_q002": "3.42", "C_f10_pm0p01_tploff_q002": "1.05e-3",
         "A_q002": "3.97", "B_q002": "1.33"}
GRAIN_A = {"A_q002": "7.635 µPa", "B_q002": "6.990 µPa"}


def e(v, d=2):
    m, x = ("%%.%de" % d % v).split("e")
    return "%s×10^{%d}" % (m, int(x))


def pc(v, d=1):
    return ("%%.%df %%%%" % d) % (100 * v)


def summary(R):
    c10, c1, off, A, B = (R[k] for k in ("C_f10_pm0p01_q002", "C_f1_pm0p01_q002", "C_f10_pm0p01_tploff_q002",
                                         "A_q002", "B_q002"))
    return [
        ("p", "The five fits use the same data, the same forward model and the same H2O law (%s). They differ only in how "
              "the H2O surface activity A of the 12,364,290 facets is parameterized, kept non-negative and regularized, "
              "and in where their chains started." % LAW),
        ("b", "**Two model families.** Options A and B describe the activity with global spherical harmonics of degree "
              "≤ 4 (25 functions) inside a non-negative wrapper: a sum of three squared harmonic fields (A, 81 free "
              "parameters) or a fixed srcCG background map plus a softplus of one harmonic field (B, 35). Option C "
              "describes it with 2745 local, non-negative hat functions on a 150-m mesh (2754 free parameters), with a "
              "Laplacian smoothness penalty of strength f and a penalty P_m on the slope multipliers."),
        ("b", "**What is minimized.** Every fit minimizes the same data cost, ½ Σ_gauges (1/N_g) Σ_i ln²(p_model/p_obs) "
              "over 934 NG and 1150 RG observations (equal total weight per gauge), plus an NG/RG crossing penalty "
              "(zero at every final state). C adds its smoothness and multiplier penalties. Between fits, a fine "
              "tuning of 271 late-seen facets minimizes the data cost plus a penalty P_δ on their log factors. Five "
              "fit + fine-tuning cycles; the final state is the fifth fine tuning (Section 6)."),
        ("b", "**Fit to the data.** Final data cost: C f10 %.3g, C f1 %.3g, B %.3g, C templates off %.3g, A %.3g "
              "(×10^{-3}); late-window rms ln(model/data) NG / RG: %.3f / %.3f (C f10) ... %.3f / %.3f (A). C f10 fits "
              "best; A worst (2.7 × C f10)." % (c10["data_cost"] * 1e3, c1["data_cost"] * 1e3, B["data_cost"] * 1e3,
                                                  off["data_cost"] * 1e3, A["data_cost"] * 1e3, c10["rms_ln_NG_late"],
                                                  c10["rms_ln_RG_late"], A["rms_ln_NG_late"], A["rms_ln_RG_late"])),
        ("b", "**Switched-off surface (user decision 1b, 2026-11-21).** A facet is active if its activity, not "
              "counting an imposed floor, exceeds T = 2.11×10^{15}, which is 1 %% of the median srcCG activity of the "
              "seen facets. A fit fails if more than half of the seen area is switched off. Switched off: A %s, "
              "C templates off %s, C f10 %s, C f1 %s, B %s. Only A passes. With its floor counted B would have "
              "%s off: most of what B keeps above T is that background floor, which the fit pushed to its lower limit. "
              "Of the ten fits of the round only A and D pass, and they have the two highest data costs (Section 9)."
              % (tuple(pc(x["off_seen"], 0) for x in (A, off, c10, c1, B)) + (pc(B["off_seen_floor_counted"], 0),))),
        ("b", "**Total water production.** The C fits give Q_{H2O} = %s ... %s s^{-1} at 10:33:10 UTC, about 90 %% of it "
              "from terrain the gauges saw; A and B give %s and %s s^{-1}, 80 %% and 75 %% of it from unseen terrain "
              "(extrapolation of the global harmonics)." % (e(min(c1["Q_H2O_1033"], off["Q_H2O_1033"]), 1),
                                                           e(c10["Q_H2O_1033"], 1), e(A["Q_H2O_1033"], 1),
                                                           e(B["Q_H2O_1033"], 1))),
        ("b", "**Energy limit.** Only C with the templates off stays below the sunlight-powered maximum flux everywhere "
              "(max F/F_E = %.2f). In the other four the flux exceeds it on 4 facets (1.1 m^{2}) carrying D2's held "
              "fine-tuning templates, which hold %s ... %s of Q." % (off["F_max_over_F_E"],
                                                                      pc(A["Q_share_over_F_E"], 2),
                                                                      pc(c1["Q_share_over_F_E"], 1))),
        ("b", "**Convergence.** Every fit and fine tuning of the three C chains ended on ftol. All five fits of A and "
              "three of the five fits of B (including the last) stopped at max_nfev = 4000: A and B are not converged. "
              "The templates-off C chain did not change after cycle 1: its fine-tuned facets carry only 3×10^{-5} of the "
              "late NG signal, so the fine tuning had nothing to adjust (Section 8.4). ftol "
              "means that a fit stopped improving near its start's minimum, not that it found the best one "
              "(Section 10)."),
        ("b", "**Not unique (Section 10).** Option C with the templates off has five solutions, one for each slope "
              "class from 40° to 90°; the fit shown here is the worst of them in both objective and data cost. All "
              "five switch off most of the seen surface. The same model with every slope multiplier held at 1 does "
              "not (35.5 % off), at a data cost of 7.62×10^{-3}, below A and D, and reaches the same solution from seven "
              "starts. Switching its cliffs off costs 4.2 % in data cost, less than switching off the 40-60° slopes "
              "(6.7 %): undecided by the test's rules."),
        ("b", "**Main caveats.** The six held templates come from old-law (q_min = 0.2 q_max) fits; the gas temperature "
              "law is held fixed; κ (CO2) is 0 in all five; the slope multipliers trade off against the activity map "
              "and are not well determined. The switched-off criterion is the user's decision 1b; the other plausibility "
              "criteria and their thresholds are proposals."),
    ]


def glance_table(R, order):
    hdr = ["#", "Fit", "Data cost (10^{-3})", "Late rms ln NG / RG", "Seen area off (floor not counted / counted)",
           "Q_{H2O} 10:33:10 (s^{-1})", "Seen share of Q", "Q above F_E", "Q max/min", "Last fit"]
    rows = []
    for i, k in enumerate(order):
        r = R[k]
        rows.append([str(i + 1), ORDER_LABEL[k], "%.3f" % (r["data_cost"] * 1e3),
                     "%.3f / %.3f" % (r["rms_ln_NG_late"], r["rms_ln_RG_late"]),
                     "%s / %s" % (pc(r["off_seen"], 0), pc(r["off_seen_floor_counted"], 0)),
                     e(r["Q_H2O_1033"]), pc(r["seen_share"], 0), pc(r["Q_share_over_F_E"], 2),
                     "%.2f" % r["Q_day_ratio"],
                     "ftol" if r["fit_termination"].startswith("`ftol`") else "max_nfev"])
    widths = [0.6, 2.7, 1.4, 1.9, 2.0, 1.9, 1.3, 1.3, 1.2, 1.4]
    cap = ("The five fits (%s): final data cost (0.5 r·r of the data rows at the last fine tuning), rms of ln(model/data) "
           "in the late window 10:33:10-10:39:26 UTC, share of the seen area switched off (activity below T = 2.11×10^{15}; decision 1b "
           "does not count B's imposed floor, the second value does), Q_H2O and its seen share at 10:33:10 UTC, share of Q from facets above the "
           "energy limit F_E = 1.09×10^{21} m^{-2} s^{-1}, ratio of the highest to the lowest Q over the 13 map times, "
           "termination of the last fit. Sources: summary.csv of the comparison set, q002_plausibility/metrics.json, "
           "full_set_q002_<run>/checks/facts.json (numbers.json)." % LAW)
    return hdr, rows, widths, cap


def settings_table(R, order):
    hdr = ["Fit", "Activity model", "Non-negativity", "Free fit parameters", "Penalties in the fit",
           "Largest multiplier", "Start of cycle 1", "Fit evaluations per cycle"]
    act = {"C": "local hat basis, 150 m (2745 nodes)", "A": "3 squared harmonic fields (degree ≤ 4)",
           "B": "srcCG map × A_bg + softplus(harmonics)"}
    nn = {"C": "bounds c_j ≥ 0", "A": "squares", "B": "floor + softplus"}
    pen = {"C_f10_pm0p01_q002": "crossing; S (f = 10); P_m (f_m 0.01)", "C_f1_pm0p01_q002": "crossing; S (f = 1); P_m (f_m 0.01)",
           "C_f10_pm0p01_tploff_q002": "crossing; S (f = 10); P_m (f_m 0.01)", "A_q002": "crossing", "B_q002": "crossing"}
    start = {"C": "D2's activity projected onto the hats (NNLS)", "A": "rank-3 Gram fit of D2's activity",
             "B": "D2's harmonics, A_bg = 0.01"}
    rows = []
    for k in order:
        r = R[k]
        o = r["option"]
        mb = max(r["multipliers"].items(), key=lambda kv: kv[1])
        rows.append([ORDER_LABEL[k], act[o] + ("" if o != "C" else (", templates " + ("off" if "tploff" in k else "on"))),
                     nn[o], str(FREE[k]), pen[k], "%.3g (%s deg)" % (mb[1], mb[0].replace("_", "-")), start[o],
                     FIT_NFEV[k]])
    widths = [2.2, 3.0, 1.6, 1.3, 2.3, 1.7, 2.3, 2.6]
    cap = ("Settings and convergence of the five fits. Free fit parameters include the six free slope multipliers, κ "
           "and the grain A, B; each fine tuning adds 271 facet factors + grain A, B. Fit evaluations: the five cycles "
           "of each chain (max_nfev 4000; 'ftol' or 'max_nfev' as recorded in fit_result.txt). Sources: the run folders' "
           "README.md, commands.sh, fit_result.txt and summary.csv.")
    return hdr, rows, widths, cap


PLOTS_INTRO = ("Figures 1-5 compare the five fits directly; Figures 6-10 show each fit's activity and production-rate "
               "maps. Colours and line styles are the same in every comparison figure: blue solid C f10 / P_m 0.01 / "
               "templates on, vermillion dashed C f1 / P_m 0.01 / templates on, green dotted C f10 / P_m 0.01 / "
               "templates off, pink dash-dotted A, orange solid B. Every figure is also in figures/ as .png and .eps "
               "(300 dpi).")


def summary_figures(R):
    data = ("Black (NG) and grey (RG) dots: the COPS observations used in the fits (934 NG, 1150 RG; no NG between "
            "10:24:38 and 10:33:10 UTC). The 95 observations from 10:33:10 UTC are compared with the gauges' response "
            "average W10x5 of the 1-s model grid, the NG model including the fitted grain term; the shaded band is that "
            "late window (297-5.5 m). Lines: the final state of each fit (label: data cost).")
    return [
        ("12M_q002_top5_altitude_NG_RG", 15.5,
         "Model vs COPS against altitude, NG (solid) and RG (dashed) of the five fits, H2O law %s. %s Hatched: the NG "
         "gap. Circles: where the model NG first equals RG in descent order (the crossing that the crossing penalty "
         "constrains). Made with model_obs_figures.py from the comparison set's model tables." % (LAW, data)),
        ("12M_q002_top5_time", 16.5,
         "Model vs COPS against time on 2016-09-30 (UTC), the same fits and data as Figure 1: (a) NG, (b) RG, (c) the "
         "last hour. %s" % data),
        ("12M_q002_top5_residuals", 16.5,
         "Log residuals ln(model/data) of the five fits against altitude (NG and RG; the late window shaded). Made with "
         "options_residual_figure.py from the comparison set's model tables and summary.csv."),
        ("12M_q002_top5_criteria", 17.0,
         "The plausibility criteria of the five fits (q002_plausibility, FINDINGS A21): (a) final data cost, (b) late "
         "rms ln(model/data) NG filled, RG open, (c) seen area switched off at T = 2.11×10^{15}: filled without B's imposed floor (decision 1b), open with it, limit 50 %, "
         "(d) share of Q_H2O above the energy limit F_E = S_0/(r_h^2 ΔH) = 1.09×10^21 m^-2 s^-1 at 3.833 AU, "
         "(e) share of Q_H2O from seen terrain, (f) highest / lowest Q_H2O over the 13 map times (activity fixed, "
         "illumination changing)."),
        ("12M_q002_top5_convergence", 12.5,
         "Data cost per chain cycle: circles after the fit of the cycle, squares after its fine tuning (the state the "
         "next cycle starts from). The templates-off C chain (green) is flat from cycle 1; A (pink) still decreases at "
         "cycle 5 (every A fit stopped at max_nfev)."),
    ]


MAPS_INTRO = ("Left: the H2O activity factor A of the final state; right: the H2O production rate F_H2O = q A "
              "(m^-2 s^-1) at 10:33:10 UTC, the first observation of the late window. Each image: (a) global "
              "longitude/latitude map with Rosetta's ground track, (b) the impact area. Colour only on the facets seen "
              "during the fit (66.5 % of the area); unseen facets grey. From the full sets full_set_q002_<run>/ "
              "(activity/ and production_rate/, pub_maps_v2), unchanged; colour classes differ between fits (see each "
              "full set's class tables).")


def map_caption(r, fs):
    return ("%s: H2O activity A (left) and H2O production rate at 10:33:10 UTC (right), final state, %s. Q_H2O at "
            "10:33:10 = %s s^-1, %s from seen terrain; seen area switched off %s (decision 1b; %s with an imposed floor counted). "
            "Source: %s/{activity,production_rate}/." % (r["label"], LAW, e(r["Q_H2O_1033"]), pc(r["seen_share"], 0),
                                                         pc(r["off_seen"], 0),
                                                         pc(r["off_seen_floor_counted"], 0), fs.name))


def methods(R, order, TEN):
    c10, c1, off, A, B = (R[k] for k in ("C_f10_pm0p01_q002", "C_f1_pm0p01_q002", "C_f10_pm0p01_tploff_q002",
                                         "A_q002", "B_q002"))
    S = []
    # ------------------------------------------------------------------ 4
    S.append((("5. The common forward model and its assumptions", 1), [
        ("p", "All five fits evaluate the gauge pressures with the same LiouvilleCOPS forward model, on the same fitting "
              "base (runs/shap8_012m_lategrid_dz1.5/base_q002) and with the same data. Only the activity map A_e "
              "differs (Section 7)."),
        ("h3", "5.1 Gas transport and geometry"),
        ("b", "Collisionless (free-molecular) flow from the surface to the gauges, Liouville theorem; only H2O and CO2; "
              "no collisions, no gravity on the gas, no photochemistry, no distributed coma (LiouvilleCOPS README "
              "§4.1)."),
        ("b", "Static kernels: the nucleus orientation at the observation time, no rotation during the flight "
              "(FlightTimeRotation = none)."),
        ("b", "Shape: SHAP8 12M (cg-spg_mspcd-shap8-v0.4cmb-12Mfacets), 12,364,290 facets; each facet subdivided "
              "4^{level} times until its edge ≤ 0.05 × the distance to the gauge (≤ 8 levels); a facet contributes only "
              "where the ray from the gauge hits it first (BVH; Liouville.cpp)."),
        ("b", "Emission: cosine (\"flat\") angular distribution with a Maxwellian flux at the surface gas temperature "
              "T_e: h(v) = 2β^{4} v^{3} exp(−β^{2} v^{2}), β = (m/2kT)^{1/2}."),
        ("b", "Both gauges at one point: the spacecraft position + 1.5 m along the spacecraft +Z axis, towards the "
              "comet (GaugeOffsetSC = 0 0 1.5; user decision 2026-11-07)."),
        ("h3", "5.2 Gauges"),
        ("p", "For each visible sub-element with solid-angle weight w = ΔA cos θ/(π r^{2}) the kernels are "
              "K_n = Σ w W_NG(l), K_Φ = Σ w W_RG(l)(−l·b_RG), K_V = Σ w W_RG(l) (Liouville.cpp), and"),
        ("eq", "n_NG = Σ_e F_e (√π/2) β_e K_n,e,     Φ_RG = Σ_e F_e [K_Φ,e + (V·b_RG)(√π/2) β_e K_V,e]", "1"),
        ("eq", "p_NG = Σ_s n_NG,s k T_NG / S_s + g(t),     p_RG = Σ_s Φ_RG,s (2π m_s k T_RG)^{1/2} / S_s", "2"),
        ("b", "NG (nude gauge, density-type): T_NG = 293 K; acceptance 'srccg_shield' with head-on weight: front "
              "(μ = l·b_NG < 0) weight 1 − μ, back the tabulated shield fraction, only for molecules arriving from the "
              "RG front; 2.95π sr in all (Gauges.cpp; FINDINGS B)."),
        ("b", "RG (ram gauge, flux-type): full front hemisphere (2π sr), flux weight −l·b_RG, equilibrium-sphere "
              "pressure model, T_RG from the COPS sensor-temperature housekeeping, spacecraft velocity V included."),
        ("b", "Sensitivities S_H2O = 0.893, S_CO2 = 0.704 (N2-equivalent pressures). The SPICE IK fields of view and "
              "spacecraft occlusion are not used."),
        ("h3", "5.3 Illumination and the source laws"),
        ("p", "cos SZA = n·ŝ with the apparent Sun direction; a facet is in night if cos SZA ≤ 0 and in cast shadow if "
              "a ray from its centroid towards the Sun hits the surface (Source.cpp; geometry fixed per base record). The "
              "H2O flux of facet e at time t is (Source.cpp, fit_cops.py; q_min = 0.02 q_max since 2026-11-18, FINDINGS "
              "A18):"),
        ("eq", "F_H2O,e(t) = A_e q_e(t),   q = q_min (night or cast shadow),   q = q_min + (q_max − q_min) cos SZA (lit)",
         "3"),
        ("eq", "q_max = 600 r_h^{−p},   p = 4.2143229,   q_min = 0.02 q_max,   r_h ≈ 3.83 AU (per record, SPICE)", "4"),
        ("p", "CO2 uses the fixed srcCG CO2 map A_CO2 with one scale κ ≥ 0 and the law q_CO2 = 600 r_h^{−2} × (0.1 in "
              "night/shadow; 0.1 + 0.9 cos SZA lit) (Source.cpp). κ ended at its lower bound 0 in all five fits."),
        ("h3", "5.4 Gas temperature (held)"),
        ("eq", "T = max(T_min, a / cos SZA + b) for lit facets with SZA ≤ 75°, T = T_min otherwise; "
               "a = (T_max − T_75)/(1 − 1/cos 75°), b = T_max − a", "5"),
        ("p", "Held in all five fits at T_max = 182.09922 K, T_75 = 136.11149 K, T_min = 133 K (--fix in commands.sh; "
              "the C design matrices are built with the same law). The temperature enters through β (Eq. 1)."),
        ("h3", "5.5 Dust grain in NG"),
        ("eq", "g(t) = max(0, A_g − B_g (t − t_end)) for t ≥ t_end = 10:33:10 UTC, 0 before;   A_g, B_g ≥ 0", "6"),
        ("p", "An additive NG term for the late NG excess (grain_crossing.py); RG has none. Final A_g: A 7.635, B 6.990, "
              "C f10 6.652, C f1 6.627, C templates off 7.206 µPa; B_g at its bound 0 in all five."),
        ("h3", "5.6 Data and the gauge time response"),
        ("b", "COPS NG and RG pressures of 2016-09-30 00:00:04-10:39:26 UTC (srcCG/RosinaMeasurements.h, = the L2 "
              "housekeeping pressures): 1150 observation rows (60-s cadence to 10:03:02, then 4 s); selection profile "
              "'all_available': 934 NG (216 missing or non-positive) and 1150 RG; no per-point weights."),
        ("b", "Late window: the 95 observations from 10:33:10 UTC (297-5.5 m) are compared with W10x5 of the model on "
              "401 model-only 1-s grid rows (10:32:46-10:39:26 UTC); earlier observations with the instantaneous model:"),
        ("eq", "p_mod(T) = (1/5) Σ_{k=0..4} W10(T − 2k s),   W10 = 10-s running mean of the surface model", "7"),
        ("p", "(the electrometer 10-s average and the DPU mean of the last five readouts; positive_problem.py). The grain "
              "term is added at the time stamp, unaveraged."),
        ("h3", "5.7 Seen facets, slope bins, templates"),
        ("b", "Seen facets: those with K_n, K_Φ or K_V > 0 for at least one fitted row: 9,211,996 facets, 66.5 % of the "
              "area. Only these are constrained by the data; the activity elsewhere is the model's extrapolation."),
        ("b", "Slope bins: θ = angle between the inward normal and the local effective gravity (homogeneous polyhedron, "
              "533 kg m^{-3}, plus rotation). Bins [0,40) (held at 1; 68.3 % of the area), [40,50) 9.6 %, [50,60) 7.7 %, "
              "[60,70) 6.1 %, [70,80) 4.1 %, [80,90) 2.3 %, [90,180] 2.0 %; one multiplier m_b ≥ 1 per bin (no upper "
              "bound) scales the whole activity of its facets."),
        ("b", "Templates: D2's six fine-tuning templates ft_c1..ft_c6, per-facet activity corrections from earlier "
              "old-law (q_min = 0.2 q_max) chains, held as fixed added terms with amplitude 1. They are non-zero on the "
              "same 265 seen facets (351 m^{2}); their sum reaches A_s ΣT = 7.894×10^{23}. Used in A, B and C "
              "'templates on'; not in C 'templates off'. D2's binary cliff template (θ ≥ 60°) is in none of the five."),
        ("h3", "5.8 Numerics"),
        ("p", "About 4.7×10^{9} kernel records; FP32 records with FP64 arithmetic ('mixed' precision) in the A/B GPU fits "
              "(one node, 8 × A100) and in the C design matrices; C fits on CPU from the design matrices (exact to ≤ "
              "9.4×10^{-15} against direct facet sums, re-checked after every fit: data cost relative difference ≤ "
              "8.9×10^{-16})."),
    ]))
    # ------------------------------------------------------------------ 5
    S.append((("6. What the fitting procedure minimizes", 1), [
        ("h3", "6.1 Data cost (all five)"),
        ("eq", "C_data = ½ Σ_{g = NG, RG} (1/N_g) Σ_{i ∈ g} [ln p_mod,g,i − ln p_obs,g,i]^{2},   N_NG = 934, N_RG = 1150",
         "8"),
        ("p", "Residuals r_i = N_g^{−1/2} (ln p_mod − ln p_obs) (fit_cops.Problem.data_residuals; model pressures floored "
              "at 10^{−300}). Each gauge has the same total weight; every point the same weight within its gauge. This is "
              "the 'data cost' of every table and figure."),
        ("h3", "6.2 Crossing penalty (all five)"),
        ("eq", "C_cross = ½ Σ_k [μ (w_NG w_RG)^{1/2} max(0, σ_k (ln p_NG,k − ln p_RG,k))]^{2},   μ = 10", "9"),
        ("p", "σ = +1 before the NG gap (model NG ≤ RG) and −1 from 10:33:10 UTC (NG ≥ RG), as observed; 839 + 95 pairs "
              "(grain_crossing.py). Zero at the final state of all five fits; not included in the data cost."),
        ("h3", "6.3 Penalties of option C only"),
        ("eq", "S = λ_C x_c^{T} L_g x_c,   x_c = c / A_s,   λ_C = f C_ref / R_pilot", "10"),
        ("p", "L_g = the P1 stiffness matrix (cotangent Laplacian) of the 150-m mesh, so S ∝ ∫ |∇(Σ_j c_j φ_j)|^{2}: a "
              "penalty on the gradient of the activity field. C_ref = 4.2688×10^{−3} (data cost of the old-law bounded "
              "softplus c9 fit), R_pilot = 517.46 (x_c^{T} L_g x_c of a 200-evaluation pilot fit): at the pilot state S = "
              "f C_ref. f = 1 → λ_C = 8.2495×10^{−6}; f = 10 → λ_C = 8.2495×10^{−5} (lambda_rule_anchored_20261114.json; "
              "user decision 2026-11-14)."),
        ("eq", "P_m = f_m C_ref Σ_{b = 6 free bins} (ln m_b)^{2},   f_m = 0.01 in the three C fits", "11"),
        ("p", "P_m keeps the slope multipliers as small as the data allow; it picks one point of the degeneracy between "
              "the multipliers and the overall level of the activity field (FINDINGS A12). A and B have no smoothness "
              "and no multiplier penalty (λ = 0: no Tikhonov term in any fit)."),
        ("h3", "6.4 Objective of a fit"),
        ("eq", "Φ_fit = C_data + C_cross  (A, B);     Φ_fit = C_data + C_cross + S + P_m  (C)", "12"),
        ("p", "Free in the fit: the activity parameters of Section 7, the six slope multipliers, κ and the grain A_g, B_g "
              "(T law and templates held; the fine-tuning factors h_e held from the previous cycle)."),
        ("h3", "6.5 Fine tuning and its penalty (all five)"),
        ("p", "After every fit, the 271 'late_only' facets (seen only in the last 5 minutes, 10:34:26-10:39:26 UTC, and "
              "never before) get per-facet factors: A_e → A_e exp(δ_e), with all other parameters fixed except the "
              "grain A_g, B_g (273 parameters). Pressures are linear in A at fixed T and κ, so dp/dA of these facets is "
              "computed once. The fine tuning minimizes"),
        ("eq", "Φ_ft = C_data + C_cross + P_δ,   P_δ = f_δ C_ref Σ_e w_e (h_e + δ_e)^{2} / Σ_e w_e,   |h_e + δ_e| ≤ ln 100",
         "13"),
        ("p", "with f_δ = 10, C_ref = 4.2688×10^{−3}, w_e = facet area, h_e = the accumulated log factor of the earlier "
              "fine tunings (nonneg_finetune.py). The cap ln 100 is a bound (no facet reached it). P_δ at the final "
              "state: %s (C f10), %s (C f1), %s (C off), %s (A), %s (B)." % (
                  e(c10["P_delta"]), e(c1["P_delta"]), e(off["P_delta"], 1), e(A["P_delta"]), e(B["P_delta"]))),
        ("h3", "6.6 The stage-2 chain"),
        ("p", "Five cycles of (fit, fine tuning). The fit of cycle k > 1 starts from the previous fit's parameters, takes "
              "the grain from the previous fine tuning and holds h = h_previous + δ. The final state is the fifth fine "
              "tuning. Because each stage omits only terms that are constant in its own variables, the chain is a "
              "block-coordinate descent on C_data + C_cross (+ S + P_m for C) + P_δ: each fit starts exactly at the "
              "data cost of the preceding fine tuning (checked in the run files)."),
        ("h3", "6.7 Optimizer and bounds"),
        ("table", ["Item", "Fits", "Fine tunings"],
         [["Method", "scipy least_squares TRF; A/B: x_scale 'jac', analytic Jacobian (GPU forward model); C: x_scale "
                     "'jac', exact trust-region solver, analytic Jacobian (checked against finite differences)",
           "TRF, tr_solver lsmr (≤ 50 iterations), x_scale 1"],
          ["Tolerances", "ftol = xtol = gtol = 10^{−10}", "10^{−10}"],
          ["max_nfev", "4000 per fit (a wall-time stop resumes as a new TRF run from the best point)", "500"],
          ["Bounds", "harmonic coefficients and A's L free; c_j ≥ 0 (C); m_b ≥ 1; κ ≥ 0; A_g, B_g ≥ 0; A_bg ≥ 0.01 (B)",
           "|h_e + δ_e| ≤ ln 100; A_g, B_g ≥ 0"],
          ["Hardware", "A, B: 8 × A100 (gpu002); C: CPU (16 threads, ≈ 11 s per evaluation)", "GPU for dp/dA"]],
         [2.4, 9.6, 5.0], 8, "Optimizer settings (fit_cops.py, fit_manifest.json, localnn_fit.py, positive_finetune.py)."),
    ]))
    # ------------------------------------------------------------------ 6
    S.append((("7. The activity models: how they differ numerically", 1), [
        ("p", "In every formula b(e) is the slope bin of facet e, m_b its multiplier (m_[0,40) = 1), T_k the held "
              "templates (amplitude 1), h_e the held fine-tuning log factor (non-zero only on the 271 late facets), y_i(e) "
              "the 25 real spherical harmonics of degree ≤ 4 at the facet centroid (srcCG functions, Source.cpp) and A_s "
              "= 2.2228×10^{18} the area-weighted RMS of the srcCG linear field (the scale of the fitter's coordinates)."),
        ("h3", "7.1 Option A: sum of squares, rank 3"),
        ("eq", "A_e = m_b(e) [ ‖L^{T} y(e)‖^{2} / A_s + A_s Σ_k T_k(e) ] exp(h_e),   L: 25 × 3 lower-trapezoidal", "14"),
        ("b", "The activity is the sum of the squares of three harmonic fields (the three columns of L), so it is "
              "non-negative by construction and zero only where all three fields vanish at once; it is a polynomial of "
              "degree ≤ 8 in the harmonics. Column 0 of L = the fitter's a00..a24, columns 1 and 2 = 24 + 23 extra "
              "coefficients (lower-trapezoidal: removes the rotation freedom L → LQ)."),
        ("b", "Free parameters: 72 (L) + 6 multipliers + κ + A_g, B_g = 81. Start: the rank-3 Gram (PSD) fit of the new "
              "D2's activity, refined by Levenberg-Marquardt (relative rms 0.127); start data cost 9.404."),
        ("b", "Behaviour: smooth and global; keeps the surface active (2 % of the seen area switched off, decision "
              "1b) but cannot produce sharp local sources, hence the largest data cost; Q mostly from unseen terrain "
              "(global extrapolation)."),
        ("h3", "7.2 Option B: background, bg 1e-2"),
        ("eq", "A_e = m_b(e) [ A_bg A_orig,e + s_e ln(1 + exp(L_e / s_e)) ] exp(h_e),   L_e = Σ_i a_i y_i(e) + A_s Σ_k T_k(e),"
               "   s_e = 10^{−5} A_orig,e", "15"),
        ("b", "A_orig = the srcCG H2O map (elements.bin). The softplus is a smooth version of max(0, L_e) (s_e is tiny), "
              "i.e. D2's clipped field; the background A_bg A_orig keeps every facet with A_orig > 0 active at a floor. "
              "Where A_orig = 0 (0.7 % of the area) the term is max(0, L) and the floor 0."),
        ("b", "Free parameters: 25 harmonics + 6 multipliers + κ + A_bg + A_g, B_g = 35. A_bg ≥ 0.01 ended at its bound "
              "0.01: the data push the background down as far as allowed. Start: D2's harmonics, A_bg = 0.01; start data "
              "cost 6.606×10^{-2}."),
        ("b", "Behaviour: the softplus part behaves like D2 (most of the field clipped). On 96 % of the seen area "
              "nothing above T remains once the imposed background is taken out; that area gives about 3 % of the seen "
              "production. With the floor counted 44 % would be switched off. Q swings ×11.8 over the rotation."),
        ("h3", "7.3 Option C: local non-negative basis"),
        ("eq", "A_e = m_b(e) max(0, Σ_j c_j φ_j(e) + TS_e) exp(h_e),   c_j ≥ 0,   TS_e = A_s Σ_k T_k(e) (on) or 0 (off)",
         "16"),
        ("b", "Mesh: 2745 nodes, n = round(2S/(√3 h^{2})) for S = 5.35×10^{7} m^{2} and h = 150 m, placed by area-weighted "
              "farthest-point sampling and 6 Lloyd iterations of a geodesic Voronoi partition on the 12M facets; 5486 "
              "coarse triangles (median edge 141 m). 2207 nodes touch seen facets; the other 538 are set only by the "
              "smoothness term (localnn_mesh.py; grid_150m/README.txt)."),
        ("b", "Basis: piecewise-linear hat functions φ_j (barycentric in the coarse triangle that contains the facet "
              "centroid's projection), three non-zero per facet, φ ≥ 0, Σ_j φ_j = 1. With c ≥ 0, TS ≥ 0 and m ≥ 1 the "
              "activity is non-negative without clipping (0 clipped facets)."),
        ("b", "Each c_j changes the activity only within ≈ 150 m of node j: the model can switch regions on and off "
              "locally. Non-negativity is a box bound in the optimizer; at the final states %s, %s and %s of the 2745 "
              "coefficients are at 0 (f10 on, f1 on, off; c5/analysis/report.txt). These are the regions the fit "
              "switches off: 66-81 %% of the seen area has A below T (decision 1b), although A is formally positive." % (CZERO["C_f10_pm0p01_q002"], CZERO["C_f1_pm0p01_q002"], CZERO["C_f10_pm0p01_tploff_q002"])),
        ("b", "Free parameters: 2745 c_j + 6 multipliers + κ + A_g, B_g = 2754, against 2084 observations: the problem is "
              "under-determined without the smoothness S (Eq. 10), which sets the unseen and weakly seen nodes."),
        ("b", "Fast evaluation: with T fixed the pressures are linear in the activity, p = Σ_b m_b (G^{(b)} c + Tcol^{(b)}) "
              "+ κ Ccol, with design matrices G (1551 rows × 7 bins × 2745 nodes) computed once from the 4.7×10^{9} "
              "kernel records (localnn_design.py, new law)."),
        ("b", "Start (cycle 1): the area-weighted non-negative least-squares projection of D2's activity onto the hats "
              "(minus the templates when they are on); 799 (on) / 811 (off) nodes start at 0."),
        ("b", "Settings of the three: f = 10 or 1 (smoothness strength, Eq. 10), f_m = 0.01 (Eq. 11), templates on or "
              "off. Not to be confused with f_δ = 10 of the fine tuning, the same in all five."),
        ("h3", "7.4 Reference: D2 (clipped harmonics; start of the chains)"),
        ("eq", "A_e = max(0, Σ_i a_i y_i(e) + A_s [c_cliff 1(θ_e ≥ 60°) + Σ_k T_k(e)])", "17"),
        ("p", "29 free parameters; the hard clip sets 91.4 % of the area to zero (FINDINGS A19), which is what the "
              "positive models A, B and C were built to avoid (user 2026-11-14)."),
        ("h3", "7.5 Summary of the numerical differences"),
        ("table", ["", "A (rank 3)", "B (bg 1e-2)", "C (three settings)"],
         [["Activity A_e", "m_b [‖L^{T}y‖^{2}/A_s + A_s ΣT] e^{h}", "m_b [A_bg A_orig + softplus(Ya + A_s ΣT)] e^{h}",
           "m_b max(0, Σ c_j φ_j + TS) e^{h}"],
          ["Basis", "25 global harmonics, 3 fields", "25 global harmonics + fixed srcCG map", "2745 local hats, 150 m"],
          ["Non-negativity", "squares (exact)", "floor + softplus", "bounds c_j ≥ 0 (TRF)"],
          ["Locality", "global: each coefficient changes A everywhere", "global", "local: ≈ 150 m per coefficient"],
          ["Free fit parameters", "81", "35", "2754"],
          ["Extra fit penalties", "none", "none", "smoothness S (f = 1 or 10); P_m (f_m = 0.01)"],
          ["Templates", "held (on)", "held (on)", "on or off"],
          ["Can switch regions off", "only where 3 fields vanish together", "yes, down to the floor (96 % of the seen area)", "yes, node by node "
                                                                                                  "(56-68 % of nodes at 0)"],
          ["Start of cycle 1", "rank-3 Gram fit of D2", "D2 harmonics, A_bg 0.01", "NNLS projection of D2"],
          ["Solver hardware", "GPU (8 × A100)", "GPU", "CPU (design matrices)"]],
         [3.0, 4.4, 4.6, 5.0], 8, "The activity models side by side (sources in Sections 7.1-7.4)."),
    ]))
    # ------------------------------------------------------------------ 7
    mrows = []
    for k in order:
        r = R[k]
        m = r["multipliers"]
        mrows.append([ORDER_LABEL[k]] + ["%.3g" % m[b] for b in ("40_50", "50_60", "60_70", "70_80", "80_90", "90_180")]
                     + ["%.1e" % r["kappa"], MAXHD[k]])
    prow = []
    for k in ("C_f10_pm0p01_q002", "C_f1_pm0p01_q002", "C_f10_pm0p01_tploff_q002"):
        s_, pm, obj = C_PEN[k]
        r = R[k]
        prow.append([ORDER_LABEL[k], "%.4e" % r["data_cost"], s_, pm, obj, "%.2f" % (float(pm) / r["data_cost"])])
    S.append((("8. How the five fits differ in their results", 1), [
        ("h3", "8.1 Fitted multipliers, κ and fine-tuning factors"),
        ("table", ["Fit", "40-50°", "50-60°", "60-70°", "70-80°", "80-90°", "90-180°", "κ", "max |h+δ|"], mrows,
         [3.0, 1.4, 1.4, 1.4, 1.4, 1.4, 1.5, 1.7, 1.6], 8,
         "Slope multipliers of the final states (bin [0,40) held at 1; lower bound 1), κ (lower bound 0) and the largest "
         "total fine-tuning log factor (cap ln 100 = 4.61). Sources: summary.csv, finetune_result.txt."),
        ("p", "Each C fit puts one bin far above 1 (×456 at 50-60° for f10, ×99 and ×102 at 40-50° for f1 and "
              "templates off) and keeps the others at 1.0-1.6; for f10 on, the 50-60° bin gives 83 % of the seen "
              "production and 91 % / 84 % of the late NG / RG signal (degeneracy.txt). B uses the 50-60°, 80-90° and "
              "overhang bins (6.7, 8.1, 13); A barely uses multipliers (2.3 at 40-50°). The bins that carry extra "
              "activity change between models and between the starts of one model (Section 10): they are not "
              "measured. The change between the old-law and new-law fits can equally be a change of local minimum "
              "(FINDINGS A22)."),
        ("h3", "8.2 Penalties of the C fits at the final state"),
        ("table", ["C fit", "Data cost", "S", "P_m", "Φ_fit (c5 fit)", "P_m / data"], prow,
         [3.4, 2.4, 2.4, 2.4, 2.8, 2.0], 8,
         "Option C, final (cycle-5) fit: data cost, smoothness S, multiplier penalty P_m, fit objective and the ratio "
         "P_m/data (c5/fit/final/fit_result.txt)."),
        ("p", "P_m is the largest penalty: the big multiplier costs 45 % of the data cost in C f10, so f_m directly limits "
              "how much activity can be moved onto one slope bin. Going from f = 10 to f = 1 (less smoothing) lowers S "
              "but raises the data cost from 3.59 to 4.17×10^{-3}. Interpretation (not tested): the two chains ended in "
              "different local optima (f = 10 puts ×456 on 50-60°, f = 1 ×99 on 40-50°), so the weaker smoothing did not "
              "buy a better fit; the total fit objective is lower for f = 1 (5.245 vs 5.423×10^{-3})."),
        ("h3", "8.3 Convergence"),
        ("p", "C: every fit and fine tuning ended on ftol; the cycle-1 fits needed a wall-time restart (≈ 3,800 "
              "evaluations), later fits 3-40. A: all five fits stopped at max_nfev 4000 and the data cost still fell by "
              "3.4 % in the last cycle (9.98 → 9.64×10^{-3}); in its last fine tuning the data cost rose slightly "
              "(9.634 → 9.640×10^{-3}) while data + P_δ fell. B: three of five fits at max_nfev (c5 included, optimality "
              "9.1×10^{-7}); its data cost changed by 0.45 % in the last cycle (5.6716 → 5.6459×10^{-3})."),
        ("h3", "8.4 Why the templates-off C chain did not move"),
        ("p", "Its cycle-1 fit stopped on ftol at 5.6734×10^{-3}, a local minimum but not the best (Section 10); the following fine tuning changed the data cost by 1.4×10^{-9} "
              "(max |δ| = 1.05×10^{-3}), and fits 2-5 each stopped after 3 evaluations with the same objective to 9 "
              "digits. Without the templates the 271 late facets hold only the smooth 150-m field: they carry 2.8×10^{-5} "
              "(NG) and 6.8×10^{-5} (RG) of the mean late-window pressure (templates on: 2.3-2.6 % and 3.9-4.2 %; "
              "computed for this document with the frozen localnn_problem / localnn_late functions, not in a run file). "
              "Even the maximum factor 100 on all of them would cost P_δ ≈ 10 × 4.27×10^{-3} × 4.6^{2} ≈ 0.9, about 160 × "
              "the data cost, so the optimum δ is ≈ 0 and each later fit starts at its own optimum. Its data cost "
              "equals the new D2's (5.6707×10^{-3}, +0.05 %). Its late NG (Figure 1) does not follow the observed rise "
              "in the last ≈ 20 m."),
    ]))
    # ------------------------------------------------------------------ 8
    S.append((("9. Plausibility of the five (criteria and caveats)", 1), [
        ("p", "The switched-off criterion is the user's decision 1b of 2026-11-21. The other criteria are the "
              "proposals of q002_plausibility/README.md; their thresholds are not decided."),
        ("b", "**Fit to the data:** C f10 < C f1 < B ≈ C off < A (data cost)."),
        ("b", "**Surface not shut down (decision 1b):** active = activity without an imposed floor above T = 2.11×10^{15} "
              "(1 % of the median srcCG activity of the seen facets; a facet at T emits 4.4×10^{15} m^{-2} s^{-1} under "
              "the subsolar Sun, about 10^{-5} of free ice sublimation). More than 50 % of the seen area switched off "
              "fails. The only imposed floor is B's background A_bg A_orig. D's spread penalty has no lower bound "
              "and is not counted as a floor, and neither are the held templates (0.001 % of the seen area)."),
        ("table", ["Fit", "Data cost (10^{-3})", "Off (decision 1b)", "Off, floor counted", "Off at T/10",
                   "Off at 10 T", "Result", "Q above F_E"],
         [[t["label"], "%.3f" % (1e3 * t["data_cost"]), pc(t["off_T"]), pc(t["off_T_floor_counted"]),
           pc(t["off_T0p1"]), pc(t["off_T10"]), "pass" if t["passes"] else "fail", pc(t["Q_share_over_F_E"], 2)]
          for t in TEN],
         [3.2, 1.8, 2.0, 2.1, 1.8, 1.8, 1.5, 2.8], 8,
         "All ten fits of the round against decision 1b: share of the seen area switched off, its sensitivity to the "
         "threshold T, and the share of Q above the energy limit. At 10 T option D is on the 50 % line. Source: "
         "q002_decisions_20261121/active_area.json, summary.csv."),
        ("b", "**Energy limit** (F_E = S_0/(r_h^{2} ΔH) = 1.09×10^{21} m^{-2} s^{-1}, all sunlight into sublimation; free "
              "sublimation of exposed ice ≈ 0.34 F_E): only C off is below everywhere; the others exceed it on the 4 "
              "template facets (0.35 % of Q in A to 1.6 % in C f1)."),
        ("b", "**Q set by the data:** C 88-92 % of Q from seen terrain; A 20 %, B 25 %."),
        ("b", "**Q stable over the rotation:** C 1.9-2.3, A 2.8, B 11.8 (max/min over the 13 map times)."),
        ("b", "**Overhang dependence:** NG at the last observation (5.5 m) comes 70-84 % from overhang facets (90-180°) in "
              "A, B and C with templates and 0 % in C with the templates off (FINDINGS A18). This is the work of the held "
              "templates; without them the NG rise of the last ≈ 20 m is not reproduced, so the data need some extra "
              "late source, but not necessarily overhangs."),
        ("p", "Assessment. Under decision 1b two of the ten fits pass, A and D, and they have the two highest data "
              "costs. D (not among the five) has 24 % of Q above the energy limit, on 98 facets. That leaves A, which "
              "exceeds the limit only on the 4 template facets (0.35 % of Q), has the worst data cost of the ten, is "
              "not converged and takes 80 % of its Q from terrain no gauge saw. The C fits, B, D2 and softplus switch "
              "off most of the seen surface. No fit of this round both keeps the surface active and fits the data as "
              "well as those do. Outside the round, the control fit of option C with all "
              "multipliers held at 1 and no templates passes (35.5 % off) with a lower data cost than A and D and no "
              "flux near the energy limit (Section 10). What would sharpen this: an independent Q_H2O for 2016-09-30, template amplitudes "
              "capped at the energy limit, and a converged fit that penalizes switched-off seen area."),
        ("h3", "Caveats common to all five"),
        ("b", "The six templates are corrections from old-law fits, held at 1 under the new law (FINDINGS A18)."),
        ("b", "The gas temperature law is held (not fitted); κ = 0 means the data are fitted by H2O alone, not that "
              "CO2 is absent (FINDINGS A14)."),
        ("b", "W10x5, the +1.5 m gauge offset and T_RG = sensor temperature are modelling decisions, not instrument "
              "documentation; the real COPS position and the RG angular acceptance are open (FINDINGS C)."),
        ("b", "The 12M shape has genus 4; its slope bins were admitted with allow_genus 4 as an experiment."),
    ]))
    return S


SOURCES_INTRO = ("Inputs read by the document generator (figures, maps, numbers), with their SHA-256 (first 16 hex "
                 "digits). The method text quotes the code and run files listed below.")
EXTRA_SOURCES = [
    "Code that ran: model_ext_frozen/nonneg_20261118T072729Z/ (A, B, fine tuning); model_ext_frozen/"
    "localnn_q002_20261118T162327Z/ (C); /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS/tools/fit_cops.py and src/ "
    "(forward model); model_ext_frozen/chain_onejob_20261118T191837Z/ (A, B chains).",
    "Runs: runs/shap8_012m_lategrid_dz1.5/nonneg/sos/stage2r_q002_r3_20261118T193956Z (A), "
    "nonneg/background/stage2r_q002_bg1e-2_20261118T193956Z (B), nonneg/localnn/stage2r_q002_{f10_pm0.01_tplon, "
    "f1_pm0.01_tplon, f10_pm0.01_tploff}_20261118T162407Z (C); clip_refit/D2_q002_warm (D2).",
    "Generators: figures_cluster/scripts/q002_top5_figures.py (figures, numbers.json), q002_top5_docx.py and "
    "q002_top5_content.py (this document); copies in AMPS/srcCG/LiouvilleCOPS_cuda/model_ext/.",
]


# ====================================================================== section 4 (added 2026-11-21, user: "enhance the
# summary: that also should include ground tracks, source rates, quntification of the contribution of different area to
# the onservations, importance of the cliffs, etc -- all useful information about the comet, sources, etc")
BINS = ["[0,40)", "[40,50)", "[50,60)", "[60,70)", "[70,80)", "[80,90)", "[90,180]"]
BIN_LAB = ["0-40°", "40-50°", "50-60°", "60-70°", "70-80°", "80-90°", "90-180°"]
FSK = {"C_f10_pm0p01_q002": "C_f10_pm0p01", "C_f1_pm0p01_q002": "C_f1_pm0p01",
       "C_f10_pm0p01_tploff_q002": "C_f10_pm0p01_tploff", "A_q002": "A", "B_q002": "B"}


def summary_sources(NS):
    """Extra summary bullets on the comet and its sources (numbers from numbers_sources.json)."""
    S = NS["runs"]
    seenQ = {k: r["Q_time"]["10:33:10"] * r["production_seen_unseen"]["10:33:10"]["seen_share"] for k, r in S.items()}
    tot = [r["Q_time"]["10:33:10"] for r in S.values()]
    big = [r["regions"]["big_lobe_NG_median_pre"] for r in S.values()]
    sng = [r["regions"]["small_lobe_median_late"][0] for r in S.values()]
    srg = [r["regions"]["small_lobe_median_late"][1] for r in S.values()]
    return [
        ("b", "**Sources (Section 4).** The production of the terrain the gauges saw is %s ... %s s^{-1} at 10:33:10 "
              "UTC in all five fits (within a factor 1.7), while the totals differ by ×%.1f: the data constrain the seen "
              "production much more than the total (one model alone, C with templates off, spans ×1.9 in seen production "
              "over its solutions, Section 10). Before the late window NG sees mostly the big lobe (median %s ... %s); in "
              "the last 6 minutes the small lobe, where Rosetta touched down, gives %s ... %s of NG and %s ... %s of "
              "RG." % (e(min(seenQ.values()), 2), e(max(seenQ.values()), 2), max(tot) / min(tot), pc(min(big), 0),
                       pc(max(big), 0), pc(min(sng), 0), pc(max(sng), 0), pc(min(srg), 0), pc(max(srg), 0))),
        ("b", "**Angle to the local gravity (Section 4.7).** With the illumination divided out (A = F/q) and facets "
              "compared within the same 150-m patch and exposure class, each slope class differs from 0-40° terrain by "
              "exactly its fitted multiplier, to within about 2 %. The slope dependence is therefore what the multipliers "
              "put in, not independent evidence. These five end states raise 40-50° or 50-60° (×2.3 to ×454), but that "
              "is not a result: from other starts C with templates off raises 60-70°, 70-80° or 80-90° and fits better "
              "(Section 10). The 60-80° multipliers are at or within 7 % of their imposed lower bound of 1, so whether such cliffs "
              "are less active than flat terrain was not tested: the fits cannot represent it. The result is the same for "
              "facets seen by at least one gauge and by both gauges in the same observation: these are 99.9 % the same "
              "facets."),
    ]


def sources(R, NS, F12):
    S = NS["runs"]
    order = list(FSK)
    blocks = []
    lab = {k: R[k]["label"] for k in order}

    def fs(k, rel):
        return F12 / ("full_set_q002_%s" % FSK[k]) / rel

    # -------------------------------------------------- 4.1 ground track
    blocks += [
        ("h2", "4.1 Ground track, visibility and the observing geometry"),
        ("p", "Rosetta descended onto the small lobe (head) of 67P on 2016-09-30; the impact point on the SHAP8 12M shape "
              "is at longitude 12.164°, latitude 31.685°. COPS observed every 60 s until 10:03:02 UTC and every 4 s "
              "after; NG has no data from 10:24:38 to 10:33:10 UTC (700-297 m). The late window 10:33:10-10:39:26 UTC "
              "(297-5.5 m above the shape) holds the 95 observations that see the touchdown region at close range; "
              "the last observation is 5.5 m above the surface. Of the surface, 66.5 % of the area (9,211,996 facets) "
              "is seen by at least one fitted observation. Over the whole descent NG and RG see the same terrain: 99.9 %% of "
              "the seen area is seen by both, always also in a common observation. Over the late window 2 facets are "
              "seen only by NG and 7 only by RG. At single instants near the ground the two differ more (10:39:26: 22,088 "
              "facets seen by both, 4,680 by NG only, 1,893 by RG only). Earlier versions of this document quoted 5,924 "
              "NG-only and 29,210 RG-only late facets: those masks belong to an earlier base (other gauge position) and "
              "do not apply to these fits.".replace("%%", "%")),
        ("fig", fs("A_q002", "production_rate/12M_groundtrack_h2o_production_rate_103926UTC_latlon.png"), 13.0,
         "Ground track of Rosetta on 2016-09-30 and its impact point, over the H2O source rate at 10:39:26 UTC of option "
         "A (the only one of the five that passes the switched-off criterion): (a) global longitude/latitude map, (b) the impact area. The track "
         "is the same for all fits. From full_set_q002_A/production_rate/."),
        ("fig", fs("A_q002", "activity/12M_groundtrack_activity_3d.png"), 17.0,
         "The ground track in 3-D over option A's H2O activity: (a) the nucleus seen down the final approach, (b) the "
         "impact region, (c) the late descent; scale bars in each panel. From full_set_q002_A/activity/."),
    ]
    # -------------------------------------------------- 4.2 source rates
    qrows = []
    for k in order:
        r = S[k]
        ps = r["production_seen_unseen"]["10:33:10"]
        il = r["illumination"]["103310"]
        q = r["Q_time"]
        qrows.append([ORDER_LABEL[k], e(q["10:33:10"]), e(q["10:39:26"]),
                      e(q["10:33:10"] * ps["seen_share"]), pc(ps["seen_share"], 0), pc(ps["unseen_ge60deg_share"], 1),
                      "%s / %s / %s" % (pc(1 - il["dark_share_of_Q"], 0), pc(il["night_share_of_Q"], 0),
                                        pc(il["cast_shadow_share_of_Q"], 0)),
                      "%s – %s" % (e(min(q.values()), 1), e(max(q.values()), 1))])
    blocks += [
        ("h2", "4.2 Source rates: how much water, from where and when"),
        ("p", "The H2O production rate of a facet is F = q A (Eq. 3): the fitted activity A times the illumination "
              "law. Q_H2O = Σ F × area. The activity maps and the 10:33:10 production maps of each fit are Figures 6-10 "
              "(Section 3.6). At 10:33:10 UTC 69.7 % of the surface area is dark (night or cast shadow), the same in all "
              "fits; dark facets produce q_min = 0.02 q_max per unit activity."),
        ("table", ["Fit", "Q 10:33:10 (s^{-1})", "Q 10:39:26", "Seen production 10:33:10", "Seen share",
                   "Unseen ≥ 60° share", "Lit / night / shadow share of Q", "Q range over the day"], qrows,
         [2.6, 1.9, 1.9, 2.1, 1.3, 1.6, 2.8, 2.6], 7.8,
         "Water production of the five fits (%s): total Q_H2O at the two map times, the production of the seen terrain "
         "and its share, the share from unseen cliffs (≥ 60°), the shares from lit, night and cast-shadow facets at "
         "10:33:10 UTC, and the range of Q over the 13 map times of the day. Sources: full_set_q002_<run>/checks/"
         "facts.json, source_rate_q002/<run>/inspection/checks/inspect_103310.json (numbers_sources.json)." % LAW),
        ("fig", "figures_sources/12M_q002_top5_Q_time.png", 13.0, None),
        ("fig", "figures_sources/12M_q002_top5_illumination.png", 17.0, None),
        ("b", "**Seen production is constrained by the data much more than the total.** The seen terrain produces %s ... %s s^{-1} "
              "in all five fits; the total ranges from %s (C f1) to %s (B). In A and B 75-80 %% of Q comes from "
              "terrain no observation saw, much of it on the night side (dark share 51 %% and 68 %%), extrapolated by the global "
              "harmonics; in the C fits ~90 %% of Q is seen and ~90 %% is lit." % (
                  e(min(S[k]["Q_time"]["10:33:10"] * S[k]["production_seen_unseen"]["10:33:10"]["seen_share"] for k in order), 2),
                  e(max(S[k]["Q_time"]["10:33:10"] * S[k]["production_seen_unseen"]["10:33:10"]["seen_share"] for k in order), 2),
                  e(S["C_f1_pm0p01_q002"]["Q_time"]["10:33:10"], 1), e(S["B_q002"]["Q_time"]["10:33:10"], 1))),
        ("b", "**Q over the day.** With the activity fixed, Q changes with the illumination as the nucleus rotates: ×1.9-2.3 "
              "in the C fits (maximum 03:00-04:00, minimum 08:00 UTC), ×2.8 in A, ×11.8 in B (up to 4.2×10^{26} s^{-1} "
              "at 07:00): a likely reason (not tested) is that B's activity sits in small areas that rotate into sunlight."),
        ("b", "**Not tested here:** an independent Q_H2O for 2016-09-30 (e.g. ROSINA/MIRO) would discriminate the totals; "
              "none was used."),
    ]
    # -------------------------------------------------- 4.3 contribution of areas
    rrows = []
    for k in order:
        r = S[k]
        rg, su = r["regions"], r["source_area_unions"]
        rrows.append([ORDER_LABEL[k], pc(rg["big_lobe_NG_median_pre"], 0),
                      "%s / %s" % (pc(rg["small_lobe_median_late"][0], 0), pc(rg["small_lobe_median_late"][1], 0)),
                      "%s / %s" % (pc(rg["shell_30_100m_median_late"][0], 0), pc(rg["shell_30_100m_median_late"][1], 0)),
                      "%s / %s" % (e(su["NG_P50"]["area_m2"], 1), e(su["RG_P50"]["area_m2"], 1)),
                      "%s / %s" % (e(su["NG_P90"]["area_m2"], 1), e(su["RG_P90"]["area_m2"], 1))])
    blocks += [
        ("h2", "4.3 Which areas the gauges' signal comes from"),
        ("p", "Each facet's contribution to each observation is the forward model term F_e × kernel_e (Eqs. 1-2); "
              "summing contributions over facet groups gives the share of the modelled signal from each group. A P % "
              "source area of an observation is the smallest set of facets, in decreasing order of contribution, that "
              "gives P % of its modelled signal; the union over the 95 late observations measures how large an area the "
              "late signal comes from."),
        ("table", ["Fit", "Big lobe, NG before late (median)", "Small lobe, late NG / RG (median)",
                   "30-100 m from touchdown, late NG / RG", "P50 union NG / RG (m^{2})", "P90 union NG / RG (m^{2})"],
         rrows, [2.8, 2.3, 2.5, 2.6, 3.2, 3.2], 7.8,
         "Contribution of areas to the modelled signal (%s): lobes are separated by the neck plane; the touchdown "
         "point is on the small lobe; source-area unions over the 95 late observations. Source: full_set_q002_<run>/"
         "checks/facts.json (regions, source_area_unions)." % LAW),
        ("fig", "figures_sources/12M_q002_top5_regions_sources.png", 17.0, None),
        ("b", "**Same in the five:** before the late window ~92 % of NG comes from the big lobe (body) in every fit; in the late "
              "window the small lobe gives 26-44 % of NG and 58-69 % of RG, and facets 30-100 m from the touchdown "
              "point 10-18 % of NG and 20-32 % of RG."),
        ("b", "**Model-dependent:** the late signal comes from 1.0-2.9×10^{4} m^{2} (P50 union, NG) in the C fits, "
              "3.9×10^{4} m^{2} in B and 7.1×10^{5} m^{2} in A: the local C basis concentrates the sources, the global "
              "harmonics of A spread them."),
        ("grid", [(fs(k, "contribution_density/12M_contrib_density_late_NG_RG.png"), lab[k]) for k in order], 8.2,
         "Contribution per unit surface area to the gas density at NG (left half of each image) and to the flux into "
         "RG (right half), mean over the 95 late observations (10:33:10-10:39:26 UTC, each the W10x5 response average), "
         "as each fitted model computes it; longitude/latitude around the impact. From full_set_q002_<run>/"
         "contribution_density/ (pre-gap and single-time maps there as well)."),
        ("grid", [(fs(k, "source_areas/12M_source_area_union_NG_RG.png"), lab[k]) for k in order], 8.2,
         "Union of the P % source areas (P = 50 darkest ... 90 lightest) of the 95 late observations, NG (left) and RG "
         "(right): global map (top) and close-up around the impact (bottom). From full_set_q002_<run>/source_areas/."),
    ]
    for k in order:
        blocks.append(("fig", fs(k, "regions/12M_region_share_lobe_NG_RG.png"), 15.5,
                       "%s: share of the modelled NG (left) and RG (right) surface pressure from the small lobe "
                       "(head, with the touchdown point) and the big lobe (body) along the descent (altitude, log); box: "
                       "the late window. From full_set_q002_%s/regions/ (distance shells in the same folder)." % (
                           lab[k], FSK[k])))
    # -------------------------------------------------- 4.4 cliffs and slopes
    area = S[order[0]]["seen_area_share"]
    prow = [["Seen area share"] + [pc(area[b], 1) for b in BINS]]
    for k in order:
        sh = S[k]["production_share_seen_103310"]
        prow.append([ORDER_LABEL[k]] + ["%.2f" % (sh[b] / area[b]) for b in BINS])
    blocks += [
        ("h2", "4.4 Importance of the slopes and cliffs"),
        ("p", "Slope = effective gravitational slope (Section 5.7). Cliffs (≥ 60°) are 14.6 % of the seen area, "
              "overhangs (90-180°) 1.9 %. Each fit has one multiplier per slope bin (Section 8.1); the table "
              "below gives the resulting production per unit area of each bin relative to the mean of the seen surface "
              "(1 = proportional to area)."),
        ("table", ["Fit"] + BIN_LAB, prow, [3.2] + [1.97] * 7, 7.8,
         "Seen H2O production at 10:33:10 UTC per unit seen area, by slope bin, relative to the seen mean (production "
         "share / area share); first row: the area shares. Source: full_set_q002_<run>/checks/facts.json "
         "(production_share_seen_103310, seen_area_share)."),
        ("fig", "figures_sources/12M_q002_top5_slope_production.png", 15.5, None),
        ("fig", "figures_sources/12M_q002_top5_slope_pressure.png", 17.0, None),
        ("b", "**These are production numbers: illumination and templates are still in them.** Steeper terrain receives "
              "less sunlight (diurnal mean q/q_max 0.29 at 0-40°, 0.19 at 80-90°, 0.15 for overhangs), so its production "
              "per unit area is lower at equal activity. Section 4.7 removes both effects."),
        ("b", "**The class with the large multiplier dominates.** Per unit area the production of one 40-60° class is "
              "raised in every fit (×2.1 at 40-50° in A, ×5.1 at 50-60° in B, ×13 at 50-60° in C f10, ×10 at 40-50° in "
              "C f1 and C templates off). This describes these end states only: which class dominates depends on the "
              "start of the fit (Section 10)."),
        ("b", "**Cliffs (≥ 60°)** give 1.4-3.4 % of the seen production and 0.2-2.7 % of the modelled signal in the C "
              "fits, 13 % and 4-9 % in A, 28 % and 14-22 % in B. These fits do not measure how active cliffs are. "
              "Their slope multipliers cannot go below 1 (imposed 2026-11-12), so no slope class can be less active "
              "than 0-40° terrain in the same neighbourhood, and at 60-80° they sit at that bound or within 7 % of "
              "it. Neither 'active cliffs are not required by the data' nor 'cliffs are as active as flat terrain' "
              "follows from these fits; both earlier statements are withdrawn."),
        ("b", "**Overhangs at touchdown.** At the last observation (5.5 m) overhang facets give 70-84 % of the modelled "
              "NG in the fits with templates and 0 % without, and ≤ 1 % of RG in all. This comes from 275 template and "
              "fine-tuned facets (367 m^{2}) that hold 11-39 % of the seen activity of those fits; without them the late NG rise is not reproduced (Section 8.4)."),
        ("b", "**Earlier cliff analysis.** At equal illumination the fitted source rate depended on slope in all six "
              "old-law final states (FINDINGS A16; figures_cluster/12M/cliff_effect_stage2/); those numbers hold for the "
              "old law q_min = 0.2 q_max only and were not redone for the new law."),
    ]
    for k in order:
        blocks.append(("fig", fs(k, "cliffs/12M_cliff_pressure_fraction_NG_RG.png"), 15.5,
                       "%s: share of the modelled NG (left) and RG (right) surface pressure from each slope bin along "
                       "the descent (altitude, log; every observation row); box: the late window. From "
                       "full_set_q002_%s/cliffs/ (with 12M_cliff_slope_bins: area, multiplier and production per bin)." % (
                           lab[k], FSK[k])))
    # -------------------------------------------------- 4.5 other properties
    blocks += [
        ("h2", "4.5 Other properties of the comet in these fits"),
        ("b", "**CO2:** the CO2 scale κ is 0 in all five (κ ≤ 6×10^{−18}): the total pressures are fitted by H2O alone. "
              "This does not measure CO2 as absent: the gauges measure total pressure only (FINDINGS A14)."),
        ("b", "**Late NG excess (dust grain term):** the five fits add 6.6-7.6 µPa to NG from 10:33:10 UTC, without decay "
              "(B_g = 0); attributed to a dust grain in NG in the model."),
        ("b", "**Gas temperature:** held at T_max = 182.1 K (subsolar), T_75 = 136.1 K, T_min = 133 K; not fitted."),
        ("b", "**Activity map:** exact zeros are rare (≤ 0.7 % of the area), but that is formal: by decision 1b "
              "66-81 % of the seen area is switched off in the C fits, 96 % in B (44 % if its imposed floor is "
              "counted) and 2 % in A."),
        ("h2", "4.6 What is robust and what is model-dependent"),
        ("table", ["Quantity", "Same in all five", "Model-dependent"],
         [["Production", "seen production 0.62-1.05×10^{25} s^{-1} at 10:33:10", "total 7.4×10^{24} ... 4.3×10^{25} s^{-1}; "
                                                                                  "share from unseen terrain 8-80 %"],
          ["Areas seen by the gauges", "NG before late: ~92 % big lobe; late: small lobe 26-44 % NG, 58-69 % RG",
           "size of the late source areas (×70 between C f10 and A)"],
          ["Slopes", "nothing: the raised class depends on the start of the fit (Section 10)", "which class and how much (×2.3-454); 80-90° in B, overhangs in B and C with templates; all of it set by the multipliers"],
          ["Switched-off surface", "-", "2 % (A) to 96 % (B) of the seen area (decision 1b)"],
          ["Touchdown", "RG at 5.5 m mostly from 0-40° and 40-60° slopes", "NG at 5.5 m from overhangs (templates on) "
                                                                            "or from 40-50° (templates off, in this one of its five solutions)"],
          ["Daily variation of Q", "-", "×1.9 (C) ... ×11.8 (B)"],
          ["CO2", "κ = 0 (pressures fitted by H2O)", "-"]],
         [3.0, 6.6, 7.4], 8, "What the five fits agree on and where they differ (Sections 4.1-4.5)."),
    ]
    return blocks

SOURCES_SECTION_INTRO = (
    "What the five fits say about the comet: where Rosetta flew and what the gauges saw (4.1), how much water the "
    "surface produced and from which illumination (4.2), which areas the observed gas came from (4.3), the role of "
    "slopes, cliffs and overhangs (4.4), other properties (4.5), what is robust across the five (4.6), and the activity against the angle to gravity with the illumination removed, with maps around the touchdown (4.7). All numbers "
    "are the fits' modelled values (final states, H2O law q_min = 0.02 q_max) from the files named in each caption; "
    "they are model results, not measurements.")


def slope_section(P2, M, out, AA):
    """Section 4.7 (user 2026-11-21: touchdown-vicinity source rate; activity against the angle to gravity with the
    illumination accounted for; facets seen by one gauge / by both). P2 = slope_activity/pass2, M = slope_activity/matched."""
    d = "slope_activity/"
    keys = list(M["runs"])
    short = [ORDER_LABEL[k.replace("_reg", "")] for k in keys]
    V = M["visibility"]
    rows = []
    for k, lab in zip(keys, short):
        B = P2["runs"][k]["any"]["bins"]
        m = M["runs"][k]["any"]
        f = lambda v: "%.3g" % v
        rows += [[lab, "mean A (10^{17})"] + [f(b["A_mean"] / 1e17) for b in B],
                 ["", "A / A(0-40°)"] + [f(b["A_mean"] / B[0]["A_mean"]) for b in B],
                 ["", "F/cos / same at 0-40°"] + [f(b["M1_F_over_cos_mean"] / B[0]["M1_F_over_cos_mean"]) for b in B],
                 ["", "matched A (proposed)"] + [f(x["ratio"]) for x in m["A"]],
                 ["", "matched F/cos"] + [f(x["ratio"]) for x in m["F_over_cos"]],
                 ["", "fitted multiplier"] + [f(x) for x in P2["runs"][k]["multipliers"]],
                 ["", "switched off (%)"] + ["%.0f" % (100 * x) for x in AA["runs"][k]["off_by_bin"]]]
    crow = []
    for k, lab in zip(keys, short):
        crow.append([lab, "any"] + ["%.3g" % x["ratio"] for x in M["runs"][k]["any"]["A"]])
        crow.append(["", "both, any time"] + ["%.3g" % x["ratio"] for x in M["runs"][k]["both"]["A"]])
        crow.append(["", "both, same observation"] + ["%.3g" % x["ratio"] for x in M["runs"][k]["both_concurrent"]["A"]])
    vrow = [["seen by at least one gauge", "%.3f %%" % (100 * V["any"])],
            ["seen by NG", "%.3f %%" % (100 * V["ng_obs"])], ["seen by RG", "%.3f %%" % (100 * V["rg_obs"])],
            ["seen by both, at the same or different times", "%.3f %%" % (100 * V["both_any_time"])],
            ["seen by both in the same observation", "%.3f %%" % (100 * V["both_obs"])],
            ["seen by both, only at different times", "%d facets (no measurable area)" % V["both_only_different_times_facets"]],
            ["seen while sunlit at least once (share of the seen area)",
             "%.1f %%" % (100 * V["seen_lit_at_least_once_of_seen"])]]
    w9 = [2.0, 2.7] + [1.76] * 7
    return [
        ("h2", "4.7 Source rate around the touchdown, and activity against the angle to the local gravity"),
        ("h3", "Source rate in the vicinity of the touchdown"),
        ("fig", d + "maps_fixed/12M_q002_top5_touchdown_source_rate_103926.png", 17.0, None),
        ("fig", d + "maps_fixed/12M_q002_top5_touchdown_source_rate_103310.png", 17.0, None),
        ("b", "**What the maps show.** Black is ≤ 10^{15} m^{-2} s^{-1}, very low production (switched off or poorly lit); pale grey, facets no observation saw. In "
              "C f10 and C f1 the touchdown point is a few tens of metres inside the edge of an active patch, and most of "
              "the terrain to its west and south has very low production. The straight edges are the footprints of the 150-m "
              "hat functions, a property of the model and not of the comet. The brightest facets are those of the slope "
              "class with the large multiplier, on the scarp 30-100 m east of the touchdown (and, with the templates "
              "off, also on the scarp to the west). A is smooth, with the 40-50° facets slightly brighter. B has an "
              "elliptical spot about 120 m long on the eastern scarp, starting about 20 m from the touchdown."),
        ("h3", "How the illumination is accounted for: the metrics"),
        ("eq", "F_e(t) = A_e q_e(t),   q = q_min (night, shadow),   q = q_min + (q_max − q_min) cos SZA (lit)", "18"),
        ("b", "**F / cos SZA** (as requested) is undefined at night and in shadow, which leaves %.0f %% of the seen area "
              "out at 10:33:10 UTC. It also retains a term A q_min / cos SZA that grows at grazing incidence (+18 %% at "
              "cos = 0.1). It is given for comparison, for lit facets with cos SZA ≥ 0.1."
              % (100 * (1 - P2["runs"][keys[0]]["any"]["M1_included_area"]))),
        ("b", "**A = F / q** removes the illumination exactly under the model's law, at every time and on every facet."),
        ("b", "**Proposed metric: matched A.** Each slope class is compared with 0-40° facets in the same 150-m patch "
              "and the same quintile of exposure over the rotation (diurnal mean q/q_max):"),
        ("eq", "R_b = median over cells c of  ⟨A⟩_{b,c} / ⟨A⟩_{0-40°,c},   cell = patch × exposure class", "19"),
        ("p", "Location (regional activity) and illumination history are then the same on both sides of the ratio, so "
              "what remains is the effect of the angle to gravity. The 275 template and fine-tuned facets (367 m^{2}) "
              "are left out: with them, overhangs would look ×10-270 as active as 0-40° terrain in the fits with templates, instead of ×1-8 (both unmatched)."),
        ("fig", d + "pass2/12M_q002_top5_slope_illumination.png", 15.5, None),
        ("h3", "Facets seen by at least one gauge"),
        ("table", ["Fit", "Statistic"] + BIN_LAB, rows, w9, 7.2,
         "Activity against the angle to the local effective gravity, facets seen by at least one gauge (area weighted). "
         "Mean A; the ratio to 0-40° without matching, for A and for F/cos SZA; the matched ratios; the fitted "
         "multiplier; the share of each class switched off by decision 1b. All other statistics for all visibility "
         "classes: slope_activity/matched/slope_statistics.csv."),
        ("fig", d + "matched/12M_q002_top5_slope_matched.png", 17.0, None),
        ("fig", d + "pass2/12M_q002_top5_slope_strata.png", 17.0, None),
        ("h3", "Facets seen by both gauges, concurrently or at different times"),
        ("table", ["Visibility class", "Share of the surface area"], vrow, [11.0, 6.0], 8,
         "Gauge visibility over the fitted observations, from a pass over all 4.74×10^{9} kernel records "
         "(slope_activity/seen_gauges.npz). Concurrent = the same facet in the same row, with data from both gauges."),
        ("table", ["Fit", "Class"] + BIN_LAB, crow, w9, 7.2,
         "The proposed (matched) metric for the three visibility classes. A class 'only at different times' cannot be "
         "analysed: it is empty."),
        ("h3", "What follows"),
        ("b", "**The choice of visibility class does not matter.** Over the descent both gauges saw the same terrain in "
              "common observations. They differ only at single instants near the ground."),
        ("b", "**The slope effect in these fits is the fitted multiplier.** The matched ratio reproduces m_b to within "
              "about 2 %, and A/m_b gives 0.99-1.02 in every class: inside a 150-m neighbourhood neither the hats nor "
              "the degree-4 harmonics can distinguish slopes."),
        ("b", "**Illumination does not create it.** The pattern persists within each illumination group, and matched "
              "F/cos and matched A agree within a few per cent. Without matching, F/cos differs from A by up to about "
              "×1.7 on steep terrain, because it samples only the lit part of each class."),
        ("b", "**What can and cannot be claimed.** These five end states put more activity on one class between 40° and "
              "60°; other starts of the same model put it anywhere from 40° to 90° (Section 10). Also "
              "the class and the factor (×2.3 to ×454) differ between the fits, so neither is measured. The multipliers "
              "of 60-80° are at or near their lower bound of 1: 'no excess' is all that can be said, and whether "
              "these cliffs are less active than flat terrain was not tested. 80-90° is raised only in B (×8.2); "
              "overhangs ×13 in B, ×1.5-1.6 in C with templates, not in A or in C without templates."),
        ("b", "**Test chosen by the user (decision 3c), done 2026-11-22:** Section 10. Several classes can take the "
              "role, so the slope class is not determined; whether the gain of any class exceeds chance needs the "
              "shuffled-label null test (3d, not run)."),
    ]


def multistart_section(T, ST):
    """Section 10 (2026-11-22, user decisions 3c and 'multi-start when ready'; FINDINGS A22)."""
    R = {r["variant"]: r for r in T["rows"]}
    lab = (("reference", "reference", "40-50° (the fit of this report)"), ("ms_only_50_60", "ms_50_60", "50-60°"),
           ("ms_only_60_70", "ms_60_70", "60-70°"), ("ms_only_70_80", "ms_70_80", "70-80° (lowest objective)"),
           ("ms_only_80_90", "ms_80_90", "80-90° (lowest data cost)"), ("none", "none", "all multipliers held at 1"))
    rows = []
    for k, s, name in lab:
        r = R[k]
        m = max(r["m_" + b] for b in ("40_50", "50_60", "60_70", "70_80", "80_90", "90_180"))
        test = "minimum" if ST[s]["local_minimum"] else "nearly"
        rows.append([name, "%.3f" % (1e3 * r["objective"]), "%.3f" % (1e3 * r["data_cost"]),
                     "-" if k == "none" else "×%.0f: %s" % (m, pc(r["top_class_share_of_seen_production"], 1)),
                     "%.3f / %.3f" % (r["late_rms_NG"], r["late_rms_RG"]),
                     "%s %s" % (pc(r["off_seen_1b"], 1), "pass" if r["passes_1b"] else "fail"),
                     "%s; %s" % (e(r["Q_103310"], 2), pc(r["Q_seen_share"], 0)),
                     "%.2g (%s)" % (r["F_max_over_F_E"], pc(r["Q_share_over_F_E"], 1)), test])
    f = R["free_min"]
    return (("10. One model, five solutions: the slope-effect test", 1), [
        ("p", "The user's test of 2026-11-21/22 took option C with the templates off and refitted it with all slope "
              "multipliers held at 1, with one class free at a time, and with all free down to 0.01. Four of the "
              "restricted fits beat the fit of this report, which a correct all-free fit cannot allow. The all-free "
              "fit (multipliers ≥ 1) was then restarted from each of their end states."),
        ("table", ["Class with the seen water", "Objective (10^{-3})", "Data cost (10^{-3})", "Multiplier: share of seen water",
                   "Late rms NG / RG", "Seen area off (1b)", "Q (s^{-1}); seen share", "Max F/F_E (Q above)",
                   "Minimum test"], rows,
         [2.4, 1.5, 1.5, 2.2, 1.9, 1.9, 2.3, 1.8, 1.5], 7.5,
         "Solutions of option C, f 10, P_m 0.01, templates off, SHAP8 12M, new H2O law: the five found with all "
         "multipliers free (≥ 1), and the same model with all held at 1. Minimum test: no descent on the bounded "
         "Gauss-Newton step and positive curvature of the exact Hessian on the free parameters ('nearly': 4×10^{-8} "
         "of the objective is left, curvature positive). Sources: report_137981, stationarity_137980 (FINDINGS A22)."),
        ("b", "**The slope class is not measured.** Each of five classes can supply almost all the seen water. The "
              "optimizer ends in the valley nearest its start, and ftol does not tell which valley that is. Other "
              "solutions may exist."),
        ("b", "**Better fit, no more active surface.** All five solutions fail the switched-off criterion at every "
              "threshold from T/10 to 10 T. The best-fitting one (80-90°) also puts 23 % of Q above the energy limit."),
        ("b", "**Without slope multipliers the surface stays active, with the same solution from every start tried.** 35.5 % of the seen area "
              "is off (33-43 % from T/10 to 10 T), with a data cost below A and D and no flux near the energy limit. Seven "
              "different starts end at the same solution, and five cycles with fine tuning leave it unchanged (FINDINGS "
              "A23); seven starts do not prove it is the best. Its late NG fit is the weak point."),
        ("b", "**Cliff test (FINDINGS A24).** In that model, cliffs (≥ 60°) held at 0.01 raise the data cost by 4.2 %; "
              "the moderate slopes (40-60°) held at 0.01 instead, a slightly larger part of the production, raise it by "
              "6.7 %. By the rules fixed before the test this is undecided (between 1 and 5 %). Nothing in it shows "
              "that cliffs matter more than other terrain, and the late NG misfit does not change."),
        ("b", "**Are cliffs less active than flat terrain? Not answered.** The fit allowed down to 0.01 stopped at "
              "max_nfev with an objective of %.3f×10^{-3}, above three solutions that were open to it. In every "
              "solution with 40-50° at the lower bound of 1 the gradient pushes it below 1: the bound is active, and "
              "how far that class would go is unknown." % (1e3 * f["objective"])),
        ("b", "**Old against new H2O law.** The old-law fit of this model was in the 80-90° solution and the new one "
              "in the 40-50° one. Their 55 % difference in data cost is mostly that change of solution: within "
              "80-90° it is 2.9 %."),
        ("b", "**Not tested:** other starts of the fits with templates, of A, B and D; the significance of any slope "
              "gain (no noise model; the shuffled-label test 3d was not run)."),
    ])
