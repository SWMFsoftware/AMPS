#!/bin/bash
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/test_chain_deps.sh (on Voyager: /data001/highresolu_vt/cops_output/jobs/test_chain_deps.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 595b7dc0f7e0085f6aa1bc185479dcdf1ca9578cf4827e7e0cf2838f2b8feb85
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job (test, CPU). Dummy-job test of chain dependency layouts on this Slurm 22.05 (2026-11-10, after the flaw
#   124014 -> 124015): the old layout reproduces the flaw; the guard layout cancels every attempt behind a failed gate
#   or stage and lets retries run; tiny CPU jobs only.
# RUN: test_chain_deps.sh --submit [--dir DIR]; test_chain_deps.sh --check DIR   (jobs 124035-124058,
#   checks/chain_deps_20261110T214653Z: 24 of 24 PASS; they were submitted with the version before the SUMMARY fix,
#   .before_verdictfix_*)
# ==== end of header ====
# Dummy-job test of chain dependency layouts on this Slurm (22.05), short CPU jobs only (normal128cores, 1 CPU,
# 200 MB, 2 min limit; each job sleeps 5 s and exits with a given code; nothing touches a GPU or a product).
# Written 2026-11-10 after the 3a chain flaw: c1ft.1 124014 (afterok:124013) was cancelled by --kill-on-invalid-dep
# when the GPU test 124013 failed, and c1ft.2 124015 (afternotok:124014) then STARTED.
#
# Scenarios (prefix depT_<tag>_):
#   old_gate_fail   the CURRENT layout of slope_cycles.sh / positive_cycles.sh behind a failed upstream: T exits 1;
#                   B1 afterok:T, B2 afternotok:B1. Expected (the flaw): B1 cancelled, never started; B2 STARTS.
#   new_gate_fail   the PROPOSED layout behind a failed single upstream: T exits 1; B1 afterok:T,
#                   B2 afterok:T,afternotok:B1. Expected: B1 and B2 cancelled, neither started.
#   new_stage_fail  PROPOSED, upstream stage of 2 attempts that both fail: U1 exit 1, U2 afternotok:U1 exit 1;
#                   guard G afterok:U1?afterok:U2 (checks U's done marker); B1 afterok:G, B2 afterok:G,afternotok:B1,
#                   B3 afterok:G,afternotok:B1:B2. Expected: U1, U2 failed; G, B1, B2, B3 cancelled, never started.
#   new_retry       PROPOSED, retries used: U1 exit 75, U2 exit 0 (writes the marker); G runs (marker present) and
#                   completes; B1 exit 75, B2 exit 0, B3 cancelled. Expected states accordingly.
#   new_first_ok    PROPOSED, first attempts succeed: U1 exit 0, U2 cancelled; G completes; B1 exit 0; B2, B3 cancelled.
# Usage:
#   jobs/test_chain_deps.sh --submit [--dir DIR]   submit all scenarios (logged in model_ext/logs/submissions.txt);
#                                                  prints DIR
#   jobs/test_chain_deps.sh --check DIR            wait until every job left the queue, compare with the expectations
#                                                  (state; started = sacct Start is a time), print PASS/FAIL per job
set -uo pipefail
ROOT=/data001/highresolu_vt/cops_output; LOG=$ROOT/model_ext/logs/submissions.txt
MODE=""; DIR=""
while [ $# -gt 0 ]; do
    case $1 in --submit) MODE=submit; shift;; --check) MODE=check; DIR=$2; shift 2;; --dir) DIR=$2; shift 2;;
               *) echo "usage: $0 --submit [--dir DIR] | --check DIR" >&2; exit 2;; esac
done
[ -n "$MODE" ] || { echo "usage: $0 --submit [--dir DIR] | --check DIR" >&2; exit 2; }
if [ "$MODE" = submit ]; then
    DIR=${DIR:-$ROOT/checks/chain_deps_$(date -u +%Y%m%dT%H%M%SZ)}
    [ -e "$DIR" ] && { echo "exists: $DIR" >&2; exit 2; }
    mkdir -p "$DIR"; cp -p "$0" "$DIR/"
    : > "$DIR/jobs.tsv"
    sub() {   # name expected_state expected_started exit_code [marker_to_write] -- sbatch dependency options...
        local name=$1 est=$2 estart=$3 code=$4 mark=$5; shift 5
        local cmd="echo \"\$(date -u +%FT%TZ) $name start\"; sleep 5"
        [ -n "$mark" ] && [ "$code" = 0 ] && cmd+="; touch $DIR/$mark"
        cmd+="; exit $code"
        local id
        id=$(sbatch --parsable --partition=normal128cores --time=00:02:00 --mem=200M --cpus-per-task=1 \
             --kill-on-invalid-dep=yes --job-name="depT_$name" --output="$DIR/%x-%j.out" "$@" --wrap="$cmd") || {
            echo "sbatch failed for $name" >&2; exit 1; }
        printf '%s\t%s\t%s\t%s\n' "$name" "$id" "$est" "$estart" >> "$DIR/jobs.tsv"
        echo "$(date -u +%FT%TZ) $id test_chain_deps.sh $name (dummy CPU job, dependency test; $DIR): $*" >> "$LOG"
        ID=$id
    }
    guard() {   # name expected_state expected_started marker -- dependency options
        local name=$1 est=$2 estart=$3 mark=$4; shift 4
        local id
        id=$(sbatch --parsable --partition=normal128cores --time=00:02:00 --mem=200M --cpus-per-task=1 \
             --kill-on-invalid-dep=yes --job-name="depT_$name" --output="$DIR/%x-%j.out" "$@" \
             --wrap="echo \"\$(date -u +%FT%TZ) $name guard\"; test -e $DIR/$mark || { echo 'upstream marker $mark missing'; exit 1; }") || {
            echo "sbatch failed for $name" >&2; exit 1; }
        printf '%s\t%s\t%s\t%s\n' "$name" "$id" "$est" "$estart" >> "$DIR/jobs.tsv"
        echo "$(date -u +%FT%TZ) $id test_chain_deps.sh $name (dummy CPU guard job, dependency test; $DIR): $*" >> "$LOG"
        ID=$id
    }
    # old layout behind a failed gate
    sub old_gate_fail_T FAILED yes 1 ""; T=$ID
    sub old_gate_fail_B1 CANCELLED no 0 "" --dependency=afterok:$T; B1=$ID
    sub old_gate_fail_B2 COMPLETED yes 0 "" --dependency=afternotok:$B1          # the flaw: it starts
    # proposed layout behind a failed gate
    sub new_gate_fail_T FAILED yes 1 ""; T=$ID
    sub new_gate_fail_B1 CANCELLED no 0 "" --dependency=afterok:$T; B1=$ID
    sub new_gate_fail_B2 CANCELLED no 0 "" --dependency=afterok:$T,afternotok:$B1
    # proposed layout, upstream stage fails
    sub new_stage_fail_U1 FAILED yes 1 sf_U.done; U1=$ID
    sub new_stage_fail_U2 FAILED yes 1 sf_U.done --dependency=afternotok:$U1; U2=$ID
    guard new_stage_fail_G CANCELLED no sf_U.done "--dependency=afterok:$U1?afterok:$U2"; G=$ID
    sub new_stage_fail_B1 CANCELLED no 0 "" --dependency=afterok:$G; B1=$ID
    sub new_stage_fail_B2 CANCELLED no 0 "" --dependency=afterok:$G,afternotok:$B1; B2=$ID
    sub new_stage_fail_B3 CANCELLED no 0 "" --dependency=afterok:$G,afternotok:$B1:$B2
    # proposed layout, retries used
    sub new_retry_U1 FAILED yes 75 rt_U.done; U1=$ID
    sub new_retry_U2 COMPLETED yes 0 rt_U.done --dependency=afternotok:$U1; U2=$ID
    guard new_retry_G COMPLETED yes rt_U.done "--dependency=afterok:$U1?afterok:$U2"; G=$ID
    sub new_retry_B1 FAILED yes 75 "" --dependency=afterok:$G; B1=$ID
    sub new_retry_B2 COMPLETED yes 0 "" --dependency=afterok:$G,afternotok:$B1; B2=$ID
    sub new_retry_B3 CANCELLED no 0 "" --dependency=afterok:$G,afternotok:$B1:$B2
    # proposed layout, first attempts succeed
    sub new_first_ok_U1 COMPLETED yes 0 fo_U.done; U1=$ID
    sub new_first_ok_U2 CANCELLED no 0 fo_U.done --dependency=afternotok:$U1; U2=$ID
    guard new_first_ok_G COMPLETED yes fo_U.done "--dependency=afterok:$U1?afterok:$U2"; G=$ID
    sub new_first_ok_B1 COMPLETED yes 0 "" --dependency=afterok:$G; B1=$ID
    sub new_first_ok_B2 CANCELLED no 0 "" --dependency=afterok:$G,afternotok:$B1; B2=$ID
    sub new_first_ok_B3 CANCELLED no 0 "" --dependency=afterok:$G,afternotok:$B1:$B2
    echo "$DIR"
    exit 0
fi
# check
[ -f "$DIR/jobs.tsv" ] || { echo "no $DIR/jobs.tsv" >&2; exit 2; }
IDS=$(cut -f2 "$DIR/jobs.tsv" | paste -sd, -)
while [ -n "$(squeue -h -j "$IDS" 2>/dev/null)" ]; do sleep 10; done
sleep 5
fail=0
while IFS=$'\t' read -r name id est estart; do
    line=$(sacct -X -n -P -j "$id" -o State,Start,ExitCode | head -1)
    st=${line%%|*}; rest=${line#*|}; start=${rest%%|*}; ec=${rest#*|}
    st=${st%% *}
    case $start in Unknown|None|"") started=no;; *) started=yes;; esac
    ok=PASS; { [ "$st" = "$est" ] && [ "$started" = "$estart" ]; } || { ok=FAIL; fail=1; }
    printf '%-22s %-8s %-10s started %-3s exit %-6s expected %-10s started %-3s %s\n' "$name" "$id" "$st" "$started" \
        "$ec" "$est" "$estart" "$ok"
done < "$DIR/jobs.tsv" | tee "$DIR/check.txt"
grep -q ' FAIL$' "$DIR/check.txt" && { echo "SUMMARY: FAIL"; exit 1; }   # the verdict column (not the state FAILED)
echo "SUMMARY: PASS"
