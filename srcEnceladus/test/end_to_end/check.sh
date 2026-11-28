#!/bin/sh
# Gating check of the Enceladus end-to-end Table tests (MakefileTest/Table,
# srcEnceladus/test/end_to_end/README.md).
#
# usage: sh srcEnceladus/test/end_to_end/check.sh NAME RUN_DIR REQUEST ZONES RANKS DIFFNUM
#   NAME     Table test name, e.g. Enceladus--E3-INMS (prefix of its references)
#   RUN_DIR  the test's run directory (run.log, PT/plots, PT/thread*)
#   REQUEST  observation block of the deck: inms_e3_h2o or uvis_2010_e10_solar
#   ZONES    expected sampling windows: total_iterations / sample_length_steps
#   RANKS    MPI ranks the references belong to; empty = the ranks of the run
#   DIFFNUM  share/Scripts/DiffNum.pl
# Run from the AMPS root. Every gating difference is written to stdout, which
# the Table redirects into test_NAME.diff: the test passes when it is empty.
#
# Gates, in order:
#   1. the run finished (run.log has the final receipt and "End of the run:");
#   2. the run used RANKS MPI ranks (one PT/thread* entry per rank);
#   3. the observation file has exactly ZONES zones (one per sampling window);
#   4. the particle ledger closes: injected == expected and ledger_residual == 0;
#   5. regression: products/<label>.dat equals
#      output/test_Enceladus/NAME--<label>.ref_np<RANKS>.gz under DiffNum.pl's
#      default tolerances (relative 1e-12, absolute 1e-30). The decks set
#      reproducible_parallel_run=true (without it, load rebalancing by
#      wall-clock time and the MPI message-arrival order of exchanged particles
#      make the same seed and rank count diverge); repeated 8-rank runs are then
#      bit-identical (2026-11-28), so any difference is a code, input, build or
#      rank-count change, not Monte Carlo noise.
# Products (RUN_DIR/products):
#   observation.dat  PT/plots/observation.REQUEST.dat without its '#' comment
#                    lines (they hold absolute input paths of this checkout);
#   ledger.dat       "key value" lines of the final receipt and species ledger
#                    from run.log, plus the global time step.
# The agreement with INMS/UVIS is reported by compare.sh and never gates.

NAME=$1; D=$2; REQUEST=$3; ZONES=$4; RANKS=$5; DIFFNUM=$6
if [ -z "$NAME" ] || [ -z "$D" ] || [ -z "$REQUEST" ] || [ -z "$ZONES" ] || [ -z "$DIFFNUM" ]; then
    echo "usage: check.sh NAME RUN_DIR REQUEST ZONES RANKS DIFFNUM"; exit 0
fi
ROOT=$(pwd)
rm -rf "$D/check" "$D/products"; mkdir -p "$D/check" "$D/products"

# 1. completion
LOG=$D/run.log
if [ ! -f "$LOG" ]; then echo "run.log is missing in $D"; exit 0; fi
grep -q "final_receipt completed_iterations=" "$LOG" || echo "run.log has no final_receipt: the run did not finish"
grep -q "End of the run:" "$LOG" || echo "run.log has no 'End of the run:': the run did not finalize"

# 2. rank count
NP=$(ls -d "$D"/PT/thread* 2>/dev/null | wc -l | tr -d ' ')
[ -n "$RANKS" ] || RANKS=$NP
[ "$NP" = "$RANKS" ] || echo "the run used $NP MPI ranks; the references are for $RANKS"

# 3. observation product and its zone count
OBS=$D/PT/plots/observation.$REQUEST.dat
if [ -f "$OBS" ]; then
    grep -v '^#' "$OBS" > "$D/products/observation.dat"
    NZ=$(grep -c '^ZONE' "$D/products/observation.dat")
    [ "$NZ" = "$ZONES" ] || echo "observation.$REQUEST.dat has $NZ zones; expected $ZONES"
else
    echo "observation product $OBS is missing"
fi

# 4. ledger product and its closure
awk '
  /\[srcEnceladus\] forward_model / {
    for (i = 1; i <= NF; i++) if ($i ~ /^global_dt_s=/) { split($i, kv, "="); print "global_dt_s", kv[2] } }
  /\[srcEnceladus\] final_receipt / {
    for (i = 1; i <= NF; i++) if ($i ~ /^[a-z_]+=-?[0-9]/) { split($i, kv, "="); print kv[1], kv[2] } }
  /\[srcEnceladus\] final_species=/ {
    species = ""
    for (i = 1; i <= NF; i++) if ($i ~ /^final_species=/) { split($i, kv, "="); species = kv[2] }
    for (i = 1; i <= NF; i++) if ($i ~ /^[a-z_]+=-?[0-9]/) { split($i, kv, "="); print species "." kv[1], kv[2] } }
' "$LOG" > "$D/products/ledger.dat"
if [ -s "$D/products/ledger.dat" ]; then
    awk '
      $1 ~ /\.ledger_residual$/ && $2 != 0 { print "particle ledger does not close: " $1 "=" $2 }
      $1 ~ /\.injected_macroparticles$/ { injected[substr($1, 1, index($1, ".") - 1)] = $2 }
      $1 ~ /\.expected_injected_macroparticles$/ { expected[substr($1, 1, index($1, ".") - 1)] = $2 }
      END { for (s in expected) if (injected[s] != expected[s])
              print s ": injected " injected[s] " macroparticles, expected " expected[s] }
    ' "$D/products/ledger.dat"
else
    echo "no final receipt or species ledger found in run.log"
fi

# 5. regression against the references
for LABEL in observation ledger; do
    REF=$ROOT/output/test_Enceladus/$NAME--$LABEL.ref_np$RANKS.gz
    if [ ! -f "$D/products/$LABEL.dat" ]; then
        echo "product $LABEL.dat is missing"
    elif [ ! -f "$REF" ]; then
        echo "reference $REF is missing"
    else
        gunzip -c "$REF" > "$D/check/$LABEL.ref"
        # DiffNum.pl keeps scratch files in its working directory: run it in check/
        (cd "$D/check" && perl "$DIFFNUM" "../products/$LABEL.dat" "$LABEL.ref") > "$D/check/$LABEL.diff" 2>&1
        cat "$D/check/$LABEL.diff"
    fi
done
exit 0
