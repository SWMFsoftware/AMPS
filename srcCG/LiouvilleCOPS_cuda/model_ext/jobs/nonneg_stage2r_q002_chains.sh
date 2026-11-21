#!/bin/bash
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/nonneg_stage2r_q002_chains.sh (on Voyager: /data001/highresolu_vt/cops_output/jobs/nonneg_stage2r_q002_chains.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: a2fde1b7f64906dae8c0e810684004ee242fc6808213dbe2d0158e98b0b858cc
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: launcher (PREPARED 2026-11-18, not yet submitted). The four regularized nonneg chains (softplus, sos r3,
#   background 1e-2, logspread f 0.1 + u0) on base_q002 (H2O q_min 0.02), every cycle-1 fit from the new D2
#   (clip_refit/D2_q002_warm); read-only pre-checks; dry run by default. Plan model_ext/logs/plan_q002_stage2r_nonneg.md.
# 2026-11-18 04:40 CST: --node-map (lead: gpu002 softplus -> B -> A, gpu001 D); pre-checks run, all OK.
# 2026-11-18 10:30 CST: --serial / --after-any-map / --post; SUBMITTED 136617-136643 (plan section 11).
# 2026-11-18 13:40 CST: posts with nonneg_maps/nonneg_contrib_20261118T192126Z, SHAPE_MANIFEST/SHAPE_REFERENCE, GEOM=job_136691;
#   submitted 136702-136727 (plan section 12).
# ==== end of header ====
# Stage-2 REGULARIZED chains of the four nonneg models on the H2O q_min 0.02 base (user decisions 2026-11-18 ~02:40 CST
# "as recommended", via the lead; plan model_ext/logs/plan_q002_stage2r_nonneg.md). PREPARED, NOT SUBMITTED until the
# lead / user OK. One one-job chain (jobs/chain_onejob.sh -> jobs/chain_onejob.sbatch -> chain_onejob.py v7.4, kind
# nonneg) per model, 5 fit / fine-tuning cycles, every cycle-1 fit DIRECTLY from the NEW D2 (no stage 1, no adoption):
#   softplus   --start d2 (default)                               D2_q002 restart
#   sos        --sos-rank 3, --start gram with the LM start file  nonneg/bundle_q002/sos_start_r3_lm.npz (from the NEW D2)
#   background --bg-min 0.01, --start d2 (default)
#   logspread  --spread-f 0.1 --spread-cref 4.2688066575e-03 --spread-u0 -1.864126788612443, --start uniform (default)
# all: --base runs/shap8_012m_lategrid_dz1.5/base_q002; --source / --ft-templates-from = clip_refit/D2_q002_warm/final;
# --slope-mult-min 1; fit NFEV 4000, CHECK 1; fine tuning late_only, WINDOW_MIN 5, MAXNFEV 500, --ft-delta-penalty-f 10
# --ft-delta-cref 4.2688066575e-03 --ft-delta-cap 4.605170185988092 (ln 100); L_THREADS 32; MODEL_EXT
# model_ext_frozen/nonneg_20261118T072729Z; 08:00:00, 3 attempts, any GPU node but gpu003, MIN_FREE_GB 20.
# Chain roots: runs/shap8_012m_lategrid_dz1.5/nonneg/<model>/stage2r_q002_<setting>_<UTC>/ (README.md, HOW-TO-RERUN.md).
#   jobs/nonneg_stage2r_q002_chains.sh --chain-code <frozen v7.4> [--utc U] [--after JOBID]
#       [--models "softplus background sos logspread"] [--node-map "softplus:gpu002,background:gpu002,sos:gpu002,logspread:gpu001"]
#       [--submit]
# --node-map (lead 2026-11-18 04:37 CST: gpu002 runs softplus -> B -> A, gpu001 runs D after agent C's chains): pins
# each chain's attempts to that node (--nodelist); models without an entry go to any node but gpu003.
# --serial "m1 m2 ..." (lead 2026-11-18 05:15 CST: gpu002 softplus -> B -> A strictly): each listed chain's attempts
#   also depend on afterany of ALL attempts of the previous listed chain (a cancelled spare attempt counts as ended).
# --after-any-map "model:ID[:ID...],..." extra afterany jobs per chain (D after agent C's two GPU attempts on gpu001).
# --post (user 2026-11-18 "include the maps in post-processing"): per chain a guard (jobs/positive_chain_guard.sbatch
#   KIND=ft TARGET=<root>/c5/ft, afterok of ANY attempt) -> jobs/nonneg_maps.sbatch (CHAINS=<root>, MAPS_NAME=maps_c5,
#   EXPECT_CYCLES=5, BASE=base_q002, SHAPE_MANIFEST=base_q002/cpp_run_manifest.json, SHAPE_REFERENCE=base/run_manifest.json,
#   MODEL_EXT model_ext_frozen/nonneg_maps_20261118T192126Z) -> jobs/nonneg_contrib.sbatch (CYCLE=5,
#   CONTRIB_NAME=contributions_c5, GEOM=maps_geometry/job_136691/geometry (validated on base_q002), BASE=base_q002,
#   MODEL_EXT model_ext_frozen/nonneg_contrib_20261118T192126Z); CPU jobs (normal128cores), afterok, --kill-on-invalid-dep=yes.
# DRY RUN by default: the pre-checks and the sbatch lines; nothing is created. --submit refuses unless every pre-check
# passes: D2_q002_warm/final complete; base_q002 complete (COPS_model.dat, elements.bin, kernels); the new D2's ft_c*
# template bundles = the old D2's (names, paths, sha256: "D2's ft templates kept"); the LM start file computed from
# the new D2's restart (sha256); at least 20 GiB free per chain still to submit (each ~17-21 GB).
set -uo pipefail
ROOT_OUT=/data001/highresolu_vt/cops_output
J=$ROOT_OUT/jobs; RUN=$ROOT_OUT/runs/shap8_012m_lategrid_dz1.5; NN=$RUN/nonneg
SNAP=$ROOT_OUT/model_ext_frozen/nonneg_20261118T072729Z
BASEQ=$RUN/base_q002; D2N=$RUN/clip_refit/D2_q002_warm/final; D2O=$RUN/clip_refit/D2/final
LM=$NN/bundle_q002/sos_start_r3_lm.npz
DELTA=(--ft-delta-penalty-f 10 --ft-delta-cref 4.2688066575e-03 --ft-delta-cap 4.605170185988092)
PY="/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python -W ignore"
UTC=$(date -u +%Y%m%dT%H%M%SZ); AFTER=""; MODELS="softplus background sos logspread"; SUBMIT=0; CC=""; NODEMAP=""
SERIAL=""; AAMAP=""; POST=0; PREV=""; LASTSERIAL=""
# posts (lead 2026-11-18 13:35 CST): the base_q002-capable snapshots (validate-geometry --shape-manifest); contributions on
# the validated base_q002 geometry of job 136691
MX=$ROOT_OUT/model_ext_frozen/nonneg_maps_20261118T192126Z; CX=$ROOT_OUT/model_ext_frozen/nonneg_contrib_20261118T192126Z
GEOMQ=$RUN/maps_geometry/job_136691/geometry; SHAPEM=$RUN/base_q002/cpp_run_manifest.json; SHAPER=$RUN/base/run_manifest.json
LOG=$ROOT_OUT/model_ext/logs/submissions.txt
die() { echo "nonneg_stage2r_q002_chains.sh: $*" >&2; exit 2; }
while [ $# -gt 0 ]; do
    case $1 in
        --utc) UTC=$2; shift 2;; --after) AFTER=$2; shift 2;; --models) MODELS=$2; shift 2;;
        --chain-code) CC=$2; shift 2;; --submit) SUBMIT=1; shift;; --node-map) NODEMAP=$2; shift 2;;
        --serial) SERIAL=$2; shift 2;; --after-any-map) AAMAP=$2; shift 2;; --post) POST=1; shift;;
        -h|--help) sed -n '2,31p' "$0"; exit 0;;
        *) die "unknown argument $1";;
    esac
done
[ -n "$CC" ] && [ -f "$CC/chain_onejob.py" ] || die "--chain-code <frozen chain_onejob v7.4 snapshot> is required"
case $UTC in [0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]T[0-9][0-9][0-9][0-9][0-9][0-9]Z) ;; *) die "--utc YYYYmmddTHHMMSSZ";; esac
case $AFTER in ""|*[!0-9]*) [ -z "$AFTER" ] || die "--after takes ONE job id";; esac
# ---- pre-checks (all read-only)
BAD=0
note() { echo "# CHECK $*"; }
for f in COPS_model.dat elements.bin kernels; do [ -e "$BASEQ/$f" ] || { note "FAIL base_q002 has no $f"; BAD=1; }; done
$PY - "$D2N" "$D2O" "$LM" <<'PYEOF' || BAD=1
import hashlib, json, sys
from pathlib import Path
import numpy as np
d2n, d2o, lm = (Path(x) for x in sys.argv[1:4])
bad = 0
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""):
            h.update(b)
    return h.hexdigest()
def say(ok, msg):
    global bad
    print("# CHECK %s %s" % ("ok  " if ok else "FAIL", msg)); bad |= (not ok)
r = d2n / "fit_restart.npz"
if not r.is_file():
    say(False, "new D2 %s: no fit_restart.npz (refit not done)" % d2n); sys.exit(1)
m = json.loads(str(np.load(r)["meta_json"]))
o = m.get("optimizer") or {}
say(m.get("status") == "complete", "new D2 %s: status %s, cost %s, %s" % (d2n.resolve(), m.get("status"), o.get("cost"),
                                                                        o.get("message")))
def fts(man):
    t = (json.loads((man / "fit_manifest.json").read_text()).get("activity_enhancement") or {}).get("activity_templates") or []
    return [(x["name"], x["bundle"]["path"], x["bundle"]["sha256"]) for x in t if str(x.get("name", "")).startswith("ft_c")]
a, b = fts(d2n), fts(d2o)
say(a == b and len(a) == 6, "ft templates of the new D2 = the old D2's (%d: %s)" % (len(a), [x[0] for x in a]))
if not lm.is_file():
    say(False, "sos LM start %s missing (jobs/nonneg_start.sbatch MODEL=sos RANK=3 REFINE=lm D=<new D2> OUT=...)" % lm)
else:
    sm = json.loads(str(np.load(lm)["meta_json"]))
    say((sm.get("d2_restart") or {}).get("sha256") == sha(r),
        "sos LM start computed from the new D2 restart (%s)" % (sm.get("d2_restart") or {}).get("path"))
sys.exit(1 if bad else 0)
PYEOF
NCH=$(echo $MODELS | wc -w); FREE=$(df -P -B1G "$ROOT_OUT" | awk 'NR == 2 {print $4}')
if [ "$FREE" -lt $((20 * NCH + 20)) ]; then note "FAIL $FREE GiB free < $((20 * NCH + 20)) (20 per chain + the 20 GiB guard)"; BAD=1
else note "ok   $FREE GiB free (>= 20 per chain + 20)"; fi
[ "$SUBMIT" = 1 ] && [ "$BAD" = 1 ] && die "pre-checks failed: nothing submitted"
AMPSC=$(git -C /data001/highresolu_vt/AMPS rev-parse HEAD 2>/dev/null)
for m in $MODELS; do
    case $m in
        softplus)   SET=smm1; TOK=();        LABEL="softplus (A = m_b s ln(1 + exp(L/s)), s = 1e-5 A_orig), the reference";;
        sos)        SET=r3; TOK=(--sos-rank 3 --start-file "$LM"); LABEL="option A (sum of squares, rank 3), LM start from the new D2";;
        background) SET=bg1e-2; TOK=(--bg-min 0.01); LABEL="option B (background floor A_bg >= 0.01)";;
        logspread)  SET=f0.1; TOK=(--spread-f 0.1 --spread-cref 4.2688066575e-03 --spread-u0 -1.864126788612443)
                    LABEL="option D (log form, spread penalty f 0.1, C_ref 4.2688066575e-03, fixed level u0 -1.864126788612443)";;
        *) die "--models: softplus, sos, background, logspread";;
    esac
    NODE=""
    for kv in ${NODEMAP//,/ }; do [ "${kv%%:*}" = "$m" ] && NODE=${kv#*:}; done
    case $NODE in ""|gpu001|gpu002) ;; *) die "--node-map $m:$NODE: gpu001 or gpu002";; esac
    AA=""
    for kv in ${AAMAP//,/ }; do [ "${kv%%:*}" = "$m" ] && AA=${kv#*:}; done
    INSERIAL=0; for x in $SERIAL; do [ "$x" = "$m" ] && INSERIAL=1; done
    [ "$INSERIAL" = 1 ] && [ -n "$PREV" ] && AA=${AA:+$AA:}$PREV
    [ "$SUBMIT" != 1 ] && [ "$INSERIAL" = 1 ] && [ -z "$PREV" ] && [ -n "$LASTSERIAL" ] && AA=${AA:+$AA:}"<attempts of $LASTSERIAL>"
    L=$NN/$m/stage2r_q002_${SET}_$UTC
    [ -e "$L" ] && die "exists: $L"
    ARGS=(--kind nonneg --model "$m" "${TOK[@]}" --slope-mult-min 1 --base "$BASEQ" --source "$D2N/fit_restart.npz"
          --ft-templates-from "$D2N/fit_manifest.json" --chain-root "$L" --model-ext "$SNAP" --chain-code "$CC"
          --cycles 5 --nfev 4000 --check 1 --ft-maxnfev 500 --ft-facets late_only --ft-window-min 5 --l-threads 32
          "${DELTA[@]}" --min-free-gb 20 --time 08:00:00 --attempts 3 ${AFTER:+--after "$AFTER"} ${NODE:+--node "$NODE"})
    if [ -n "$AA" ]; then case $AA in *"<"*) echo "# (dry run) --after-any $AA";; *) ARGS+=(--after-any "$AA");; esac; fi
    echo "# ===== $m: $L"
    if [ "$SUBMIT" != 1 ]; then
        $J/chain_onejob.sh "${ARGS[@]}" | grep -E '^  sbatch|^# (one-job|cycle-1|attempt|NOTE)'
        [ "$INSERIAL" = 1 ] && LASTSERIAL=$m
        [ "$POST" = 1 ] && echo "# post: guard (afterok of any attempt) -> nonneg_maps.sbatch -> nonneg_contrib.sbatch (BASE=$BASEQ)"
        continue
    fi
    mkdir -p "$L"
    cat > "$L/README.md" <<EOF
# Stage 2 (regularized), H2O q_min 0.02 law: $m $SET, 5 fit / fine-tuning cycles (one job per chain)

Created $UTC by \`jobs/nonneg_stage2r_q002_chains.sh\` (agent 3). Plan: \`model_ext/logs/plan_q002_stage2r_nonneg.md\`. Job ids: HOW-TO-RERUN.md.

- **Model:** $LABEL; slope-bin multipliers >= 1, D2's ft_c1..ft_c6 templates held at 1; D2's scheme otherwise.
- **Law:** H2O night/shadow q_min = 0.02 q_max (user 2026-11-18, the x10 bug fixed; FINDINGS A18); base \`$BASEQ\`.
- **Decisions** (user 2026-11-18 ~02:40 CST "as recommended"): constants of the 0.2-law fits kept as numbers (C_ref 4.2688066575e-03,
  u0 -1.864126788612443); cycle-1 fit directly from the NEW D2 (\`$D2N\`) with the earlier start rules; D2's ft templates kept;
  regularized fine tuning (late_only 5, MAXNFEV 500, f_delta 10, C_ref 4.2688066575e-03, cap ln 100), 5 cycles, --slope-mult-min 1,
  NFEV 4000, CHECK 1.
- **Code:** model \`$SNAP\`; chain driver \`$CC\` (v7.4); AMPS $AMPSC (+ the night-fix edits in the tracked tree).
- Reaching max_nfev is not convergence; compare fits by final cost, termination and active bounds.
EOF
    { echo "# How to rerun: $m $SET on base_q002 (5 cycles, one job)"; echo; echo '```bash'
      echo "export COPS_OUTPUT_ROOT=$ROOT_OUT"; printf '%q ' "$J/chain_onejob.sh" "${ARGS[@]}" --submit; echo; echo '```'; } > "$L/HOW-TO-RERUN.md"
    OUT=$($J/chain_onejob.sh "${ARGS[@]}" --submit) || { echo "$OUT"; die "chain_onejob.sh failed for $m"; }
    echo "$OUT" | grep -v '^#   '
    echo "$OUT" | grep -E '^attempt [0-9]+: job [0-9]+$' | sed 's/^/- /' >> "$L/HOW-TO-RERUN.md"
    IDS=$(echo "$OUT" | sed -n 's/^attempt [0-9]*: job \([0-9]*\)$/\1/p' | paste -sd: -)
    [ -n "$IDS" ] || die "no attempt ids for $m"
    [ "$INSERIAL" = 1 ] && PREV=$IDS
    if [ "$POST" = 1 ]; then
        OR=$(echo "$IDS" | sed 's/:/?afterok:/g')
        G=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency=afterok:$OR --job-name=guard_q002_${m} \
            --export=ALL,KIND=ft,TARGET=$L/c5/ft $J/positive_chain_guard.sbatch) || die "guard sbatch failed for $m"
        echo "$(date -u +%FT%TZ) $G jobs/positive_chain_guard.sbatch normal128cores q002 $m post guard (KIND=ft TARGET=$L/c5/ft): afterok:$OR" >> "$LOG"
        MJ=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency=afterok:$G --job-name=maps_q002_${m} \
            --export=ALL,CHAINS=$L,MAPS_NAME=maps_c5,EXPECT_CYCLES=5,BASE=$BASEQ,MODEL_EXT=$MX,SHAPE_MANIFEST=$SHAPEM,SHAPE_REFERENCE=$SHAPER \
            $J/nonneg_maps.sbatch) || die "maps sbatch failed for $m"
        echo "$(date -u +%FT%TZ) $MJ jobs/nonneg_maps.sbatch normal128cores q002 $m post maps: CHAINS=$L MAPS_NAME=maps_c5 EXPECT_CYCLES=5 BASE=$BASEQ MODEL_EXT=$MX SHAPE_MANIFEST=$SHAPEM SHAPE_REFERENCE=$SHAPER afterok:$G" >> "$LOG"
        CJ=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency=afterok:$MJ --job-name=contrib_q002_${m} \
            --export=ALL,CHAIN=$L,CYCLE=5,CONTRIB_NAME=contributions_c5,GEOM=$GEOMQ,BASE=$BASEQ,MODEL_EXT=$CX \
            $J/nonneg_contrib.sbatch) || die "contrib sbatch failed for $m"
        echo "$(date -u +%FT%TZ) $CJ jobs/nonneg_contrib.sbatch normal128cores q002 $m post contributions: CHAIN=$L CYCLE=5 CONTRIB_NAME=contributions_c5 GEOM=$GEOMQ BASE=$BASEQ MODEL_EXT=$CX afterok:$MJ" >> "$LOG"
        echo "post $m: guard $G, maps $MJ, contributions $CJ"
        { echo; echo "Post-processing: guard $G (afterok of any attempt), maps $MJ (\`maps_c5\`), contributions $CJ (\`contributions_c5\`)."; } >> "$L/HOW-TO-RERUN.md"
    fi
done
[ "$SUBMIT" = 1 ] || echo "# dry run: nothing created or submitted (add --submit after the OK)"
