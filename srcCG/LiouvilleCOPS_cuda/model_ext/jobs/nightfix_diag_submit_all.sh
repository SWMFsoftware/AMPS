#!/bin/bash
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/nightfix_diagnostic_20261118/jobs/submit_all.sh (on Voyager: /data001/highresolu_vt/cops_output/tests/nightfix_diagnostic_20261118/jobs/submit_all.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 38a3b4a84466df092aff5acc23cebd4ceb70afce0b6e4b068e1684e508c5cfcd
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: launcher (2026-11-18). Submits the eight final states (v1) and the base consistency check.
# RUN: as run from tests/nightfix_diagnostic_20261118/jobs/ (jobs 136436-136444).
# ==== end of header ====
# H2O night-side fix diagnostic (2026-11-18): submit the eight final states and the base consistency check (CPU jobs,
# normal128cores, 16 cores, 96 GB each). NPYDIS=1 everywhere: the CPU nodes (AMD EPYC 9554) have AVX-512, the fits' hosts
# (EPYC 7513) do not; without it the A re-derivation of the state checks is not bitwise (first attempt 136407 stopped
# at "state checks failed": A_fit_reproduced_bitwise false; with it: true, srun check 2026-11-18).
set -eu
cd "$(dirname "$0")"
R=/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5
F=/data001/highresolu_vt/cops_output/model_ext_frozen
CODE=$F/nightfix_diag_20261118T073237Z
G1=$R/maps_geometry/job_136152/geometry G2=$R/maps_geometry/job_136166/geometry G3=$R/maps_geometry/job_136219/geometry
G4=$R/maps_geometry/job_124775/geometry
sub() { local key=$1; shift; echo -n "$key: "; sbatch --parsable --job-name=nightfix_diag_$key --export=ALL,CODE=$CODE,NPYDIS=1,KEY=$key,"$@" nightfix_diag.sbatch; }
N=$R/nonneg
sub stage2_softplus_reg FT=$N/softplus/stage2r_smm1_20261116T160802Z/c5/ft/final,FIT=$N/softplus/stage2r_smm1_20261116T160802Z/c5/fit/final,GEOM=$G1,RECORDED=4.5925244867e-03
sub stage2_A_reg FT=$N/sos/stage2r_r3_20261116T160802Z/c5/ft/final,FIT=$N/sos/stage2r_r3_20261116T160802Z/c5/fit/final,GEOM=$G2,RECORDED=1.0693676519e-02
sub stage2_B_reg FT=$N/background/stage2r_bg1e-2_20261116T160802Z/c5/ft/final,FIT=$N/background/stage2r_bg1e-2_20261116T160802Z/c5/fit/final,GEOM=$G2,RECORDED=5.1201315445e-03
sub stage2_D_reg FT=$N/logspread/stage2r_f0.1_20261116T160802Z/c5/ft/final,FIT=$N/logspread/stage2r_f0.1_20261116T160802Z/c5/fit/final,GEOM=$G3,RECORDED=8.1191532728e-03
L1=$N/localnn/stage2r_f1_pm1e-4_tplon_20261116T055923Z/c5/tools_state_ft
sub stage2_C_f1_pm1em4_reg FT=$L1,FIT=$L1/fit_state,GEOM=$G4,ADAPTER=$F/localnn_s2jc2_20261117T080730Z/localnn_tools.py,RECORDED=2.4358002153e-03
sub D2 FT=$R/clip_refit/D2/tools_state,FIT=$R/clip_refit/D2/tools_state,GEOM=$G4,ADAPTER=$F/d2_tools_20261117T225752Z/d2_tools.py,RECORDED=4.886437e-03
L2=$N/localnn/stage2r_f10_pm0.01_tplon_20261116T044825Z/c5/tools_state_ft
sub stage2_C_f10_pm0.01_reg FT=$L2,FIT=$L2/fit_state,GEOM=$G4,ADAPTER=$F/localnn_s2r_20261116T044801Z/localnn_tools.py,RECORDED=3.1996268075e-03
sub softplus_c9 FT=$N/softplus/chain_5cycles_smm1_20261113T041215Z/c9/ft/final,FIT=$N/softplus/chain_5cycles_smm1_20261113T041215Z/c9/fit/final,GEOM=$G4,RECORDED=4.2688066575e-03
echo -n "base_run: "; sbatch --parsable --job-name=nightfix_diag_base_run --export=ALL,CODE=$CODE,NPYDIS=1,BASE_RUN=1 nightfix_diag.sbatch
