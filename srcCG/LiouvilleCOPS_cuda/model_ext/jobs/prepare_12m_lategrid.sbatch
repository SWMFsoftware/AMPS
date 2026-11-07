#!/bin/bash
# ==== version-controlled copy: header added 2026-11-07, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/prepare_12m_lategrid.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/prepare_12m_lategrid.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: e89c87fb06a3e6e2de50302e4ef9d995e497942cdca6394818cd38c873a638ed
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job. Fitting base of the positive fits (user decisions 2026-11-07): observation file with the 1-s late grid
#   (late_time_grid.py), base with GaugeOffsetSC = 0 0 1.5 (native/liouville_cops_cuda), bundle links, verification
#   against the 12M base and the late_response base_dz1.5; STOP_AFTER=base (no record store: disk); free-space guard.
# RUN: sbatch prepare_12m_lategrid.sbatch   (job 122291)
# ==== end of header ====
# New 12M fitting base + FP32 record store: all 1150 original observations plus a uniform 1-s grid in the late window,
# gauge offset GaugeOffsetSC = 0 0 1.5 m.
# Why (user decisions 2026-11-07): the COPS time response (10-s running average and 10 s + DPU 5-readout average, ROSINA
# COPS instrument modes RO-ROS-MAN-1019 sect. 4.2) can be fitted only if the fitting base has the model on a 1-s grid in
# the late window; the gauge position is taken 1.5 m along +Z of the spacecraft (toward the comet) from the trajectory
# point (COPS position relative to the s/c reference point not known; late_response.sbatch job 122289 sensitivity check).
# 1. model_ext/late_time_grid.py --keep-from 2016-09-30T00:00:00 --allow-missing-l2: the 1150 observations of
#    srcCG/RosinaMeasurements.h (times and reference values unchanged) + the grid 10:32:46.282 .. 10:39:26.282 step 1 s
#    (401 points, reference 0 = gap, never fitted), merged in time order (1551 points) -> <run>/late_grid.h, and the
#    point table with the L2 START/STOP of every observation -> <run>/late_grid_points.txt (observations 553, 554 =
#    9:13:04, 9:14:04 have no L2 bg file in the HK table: L2 "-").
# 2. base: the existing 12M base's boundary.in (runs/shap8_012m_flat_earlier/base/boundary.in, built by the runner with
#    rebuild_cops_boundary.py) with only ObservationFile, OutputDirectory, LastPoint, GaugeOffsetSC and Threads changed
#    (checked), run with the same executable native/liouville_cops_cuda (sha256 checked) from ../LiouvilleCOPS, one GPU.
# 3. store: tools/c15_record_store.py --run <base> --out <store> --threads 8, the command of the runner's store stage
#    (runs/shap8_012m_flat_earlier/stamps/store.json).
# 4. bundle: the 12M restart template bundle is mesh-based; linked (symlinks) only if this base's elements.bin is
#    byte-identical to the one the bundle records (then the fitter's mesh identity = sha256(centroid, area) matches).
# 5. verification: 1551 points, originals keep reference values and use flags of the fitting base, grid rows unused,
#    altitude shift of the dz offset, agreement with the late-window base dz1.5 of job 122289, store record count.
# Products: runs/shap8_012m_lategrid_dz1.5/{late_grid.h,late_grid_points.txt,base,store,bundle,provenance/prepare-<job>}
# Disk (estimated from the 12M base and job 122289): base about 172 GB, store about 83 GB. The job refuses to start a
# stage when the file system has less than NEED_BASE_GB / NEED_STORE_GB free.
# Resumable: a complete base (run_manifest process complete) or store (store_manifest.json) is kept; an incomplete one is
# moved aside as *.incomplete_<UTC>. STOP_AFTER=base stops after the base and its verification.
#   sbatch prepare_12m_lategrid.sbatch                       (base + store + bundle links + verification)
#   sbatch --export=ALL,STOP_AFTER=base prepare_12m_lategrid.sbatch
#SBATCH --job-name=cops12m_lgrid
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_lgrid-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
DZ=1.5; THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-48}}; STORE_THREADS=8; STOP_AFTER=${STOP_AFTER:-}
NEED_BASE_GB=${NEED_BASE_GB:-185}; NEED_STORE_GB=${NEED_STORE_GB:-95}
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier; RUN=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz$DZ
X=$COPS_OUTPUT_ROOT/model_ext; PY="$APP/.venv/bin/python -W ignore"
BIN=$APP/native/liouville_cops_cuda; BIN_SHA=34a7551fc5d1381d3c2d5e59f7ab56a2cbbdb8c0222811f126bc72039c0dc0bd
B=$RUN/base; S=$RUN/store; LATE=$R/late_response/122289/base_dz1.5
OBS=(--observations "$AMPS/srcCG/RosinaMeasurements.h" --hk "$COPS_OUTPUT_ROOT/data/COPS_HK_20160930.txt"
     --start 2016-09-30T10:32:46.282 --stop 2016-09-30T10:39:26.282 --step 1 --keep-from 2016-09-30T00:00:00 --allow-missing-l2)
now() { date -u +%FT%TZ; }
stamp() { date -u +%Y%m%dT%H%M%SZ; }
freegb() { df -P -B1G "$COPS_OUTPUT_ROOT" | awk 'NR==2{print $4}'; }
base_complete() {
    [ -f "$B/run_manifest.json" ] && [ -f "$B/COPS_model.dat" ] && grep -q '^Done: ' "$B/build.log" 2>/dev/null &&
    $PY -c 'import json,sys; s=json.load(open(sys.argv[1]))["states"]; sys.exit(0 if s.get("process")=="complete" and s.get("cache")=="verified" else 1)' "$B/run_manifest.json"
}
store_complete() { [ -f "$S/store_manifest.json" ]; }
[ "$(sha256sum < "$BIN" | cut -d' ' -f1)" = "$BIN_SHA" ] || { echo "liouville_cops_cuda differs from the executable of the 12M base" >&2; exit 2; }
# a new run directory, or one this script made (resume); an earlier job that stopped before the observation file left
# only provenance/
RESUME=0
if [ -e "$RUN" ]; then
    if [ -f "$RUN/late_grid.h" ] && [ -f "$RUN/late_grid_points.txt" ]; then RESUME=1
    elif [ "$(ls -A "$RUN")" != provenance ]; then echo "exists and not made by this script: $RUN" >&2; exit 2; fi
fi
need=0; base_complete || need=$((need + NEED_BASE_GB)); [ "$STOP_AFTER" = base ] || store_complete || need=$((need + NEED_STORE_GB))
[ "$(freegb)" -ge "$need" ] || { echo "$(now) refused: $(freegb) GB free on $COPS_OUTPUT_ROOT, this job needs about $need GB" >&2; exit 6; }
P=$RUN/provenance/prepare-$SLURM_JOB_ID; mkdir -p "$P"; cd "$AMPS/srcCG/LiouvilleCOPS" || exit 1   # rebuild.json cwd
cp -p "$0" "$P/prepare_12m_lategrid.sbatch"; cp -p "$X/late_time_grid.py" "$P/"
sha256sum "$BIN" "$APP/tools/c15_record_store.py" "$X/late_time_grid.py" > "$P/executables.sha256"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(CUDA_[A-Z_]*|CUPY_[A-Z_]*|COPS_[A-Z_]*|THREADS|STOP_AFTER|NEED_[A-Z_]*|OMP_[A-Z_]*|OPENBLAS_[A-Z_]*|MKL_[A-Z_]*|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH|SLURM_JOB_GPUS|CUDA_VISIBLE_DEVICES)=' | sort > "$P/environment.txt"
# commands.sh: every command of this job with all arguments resolved, runnable as is (HOW-TO-RERUN.md)
{ echo '#!/bin/bash'; echo "# job $SLURM_JOB_ID on $(hostname), $(now); AMPS $(git -C "$AMPS" rev-parse HEAD)"; grep '^export ' "$0"; echo "cd $AMPS/srcCG/LiouvilleCOPS"; } > "$P/commands.sh"
rec() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; }
recto() { local o=$1; shift; printf '%q ' "$@" >> "$P/commands.sh"; printf '> %q\n' "$o" >> "$P/commands.sh"; }
free -g > "$P/free_at_start.txt"; df -h "$COPS_OUTPUT_ROOT" > "$P/df_at_start.txt"
if [ -n "${CUDA_VISIBLE_DEVICES:-}" ]; then NGPU=$(echo "$CUDA_VISIBLE_DEVICES" | tr , "\n" | wc -l); else NGPU=$(nvidia-smi -L | wc -l); fi
echo "$(now) job $SLURM_JOB_ID on $(hostname): $NGPU GPU(s) visible (CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}), threads $THREADS, resume $RESUME, $(freegb) GB free (needs about $need GB)"

# 1. observation file and point table
recto "$P/late_time_grid.log" $PY "$X/late_time_grid.py" "${OBS[@]}" --out "$RUN/late_grid.h" --table "$RUN/late_grid_points.txt"
if [ $RESUME = 0 ]; then
    $PY "$X/late_time_grid.py" "${OBS[@]}" --out "$RUN/late_grid.h" --table "$RUN/late_grid_points.txt" > "$P/late_time_grid.log" 2>&1 || { cat "$P/late_time_grid.log"; exit 3; }
else   # the files of the first job must be what the helper writes now
    $PY "$X/late_time_grid.py" "${OBS[@]}" --out "$P/late_grid.h" --table "$P/late_grid_points.txt" > "$P/late_time_grid.log" 2>&1 || { cat "$P/late_time_grid.log"; exit 3; }
    cmp "$P/late_grid.h" "$RUN/late_grid.h" && cmp "$P/late_grid_points.txt" "$RUN/late_grid_points.txt" || { echo "observation file changed since the first job" >&2; exit 3; }
fi
cat "$P/late_time_grid.log"
NP=$(grep -vc '^#' "$RUN/late_grid_points.txt"); [ "$NP" -eq 1551 ] || { echo "expected 1551 points, got $NP" >&2; exit 3; }

# 2. base
if base_complete; then
    echo "$(now) base: complete, kept ($B)"
else
    [ -e "$B" ] && { mv "$B" "$B.incomplete_$(stamp)"; echo "$(now) moved an incomplete base aside: $B.incomplete_*"; }
    mkdir -p "$B"
    SED=(sed -e "s|^ObservationFile *=.*|ObservationFile = $RUN/late_grid.h|"
         -e "s|^OutputDirectory *=.*|OutputDirectory = $B|"
         -e "s|^FirstPoint *=.*|FirstPoint = 0|" -e "s|^LastPoint *=.*|LastPoint  = $((NP - 1))|"
         -e "s|^GaugeOffsetSC *=.*|GaugeOffsetSC  = 0 0 $DZ            # prepare_12m_lategrid.sbatch: COPS 1.5 m along +Z of the s/c (toward the comet) [m], user decision 2026-11-07|"
         -e "s|^Threads *=.*|Threads = $THREADS|" "$R/base/boundary.in")
    recto "$B/boundary.in" "${SED[@]}"
    "${SED[@]}" > "$B/boundary.in"
    diff "$R/base/boundary.in" "$B/boundary.in" > "$B/boundary.in.diff"
    KEYS='^(ObservationFile|OutputDirectory|FirstPoint|LastPoint|GaugeOffsetSC|Threads) *='
    nexp=8; [ "$THREADS" = 64 ] || nexp=10       # 4 or 5 changed lines (Threads unchanged only with 64)
    [ "$(grep -c '^[<>]' "$B/boundary.in.diff")" -eq $nexp ] || { echo "unexpected boundary.in changes:"; cat "$B/boundary.in.diff"; exit 4; }
    diff <(grep -vE "$KEYS" "$R/base/boundary.in") <(grep -vE "$KEYS" "$B/boundary.in") > /dev/null || { echo "boundary.in differs outside the intended keys" >&2; exit 4; }
    for want in "ObservationFile = $RUN/late_grid.h" "OutputDirectory = $B" "FirstPoint = 0" "LastPoint  = $((NP - 1))" "Threads = $THREADS"; do
        [ "$(grep -cxF "$want" "$B/boundary.in")" -eq 1 ] || { echo "boundary.in: missing '$want'" >&2; exit 4; }
    done
    [ "$(grep -cE "^GaugeOffsetSC  = 0 0 $DZ " "$B/boundary.in")" -eq 1 ] || { echo "boundary.in: GaugeOffsetSC not set" >&2; exit 4; }
    echo "$(now) base: boundary.in checked (changes:"; sed 's/^/    /' "$B/boundary.in.diff"; echo "  )"
    echo "env COPS_CUDA_DEVICE=0 bash -c 'echo 1000 > /proc/self/oom_score_adj; exec \"\$@\"' _ $BIN $B/boundary.in > $B/build.log 2>&1" >> "$P/commands.sh"
    echo "$(now) base: start ($NP points, device 0, $(freegb) GB free)"
    COPS_CUDA_DEVICE=0 bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "$BIN" "$B/boundary.in" > "$B/build.log" 2>&1; rc=$?
    echo "$(now) base: exit $rc"; tail -4 "$B/build.log" | cut -c1-220
    [ $rc -eq 0 ] && base_complete || { echo "$(now) base build failed" >&2; exit 5; }
    du -sh "$B" | sed 's/^/  /'
fi

# 3. bundle links (mesh identity: byte-identical elements.bin)
EL_REC=$($PY -c 'import json,sys; print(json.load(open(sys.argv[1]))["elements"]["sha256"])' "$R/bundle/restart_templates.json")
EL_NEW=$(sha256sum < "$B/elements.bin" | cut -d' ' -f1)
echo "elements.bin sha256: new $EL_NEW, bundle records $EL_REC" | tee "$P/elements_sha256.txt"
if [ "$EL_NEW" = "$EL_REC" ]; then
    if [ ! -e "$RUN/bundle" ]; then
        mkdir -p "$RUN/bundle"
        for f in restart_templates.npz restart_templates.json; do rec ln -s "$R/bundle/$f" "$RUN/bundle/$f"; ln -s "$R/bundle/$f" "$RUN/bundle/$f"; done
        { echo "$(now) job $SLURM_JOB_ID: symlinks to $R/bundle (made by the runner's bundle stage, c05_write_bundle.py, job 120818)."
          echo "Valid for this base: elements.bin is byte-identical (sha256 $EL_NEW), so the mesh identity (sha256 of the facet"
          echo "centroids and areas, activity_templates.load_bundle) matches. Columns cliff_ge_45deg and cliff_ge_60deg depend only"
          echo "on the mesh and the gravity cache. late_ng_not_rg and late_rg_not_ng are visibility masks of the dz=0 base"
          echo "(1150 points), not recomputed for this base (GaugeOffsetSC 0 0 1.5, 1551 points)."; } > "$RUN/bundle/LINKED_FROM.txt"
    fi
    echo "$(now) bundle: linked ($RUN/bundle)"
else
    echo "$(now) bundle: NOT linked (elements.bin differs from the bundle's base)"
fi

# 4. verification of the base
cat > "$P/verify_base.py" <<'EOF'
"""prepare_12m_lategrid.sbatch: check the new base against the point table, the 12M fitting base and job 122289's dz1.5 base."""
import json, sys
from datetime import datetime, timedelta
run, old, late, table = sys.argv[1:5]
def utc(s):
    d, t = s.split("T"); h, m, sec = t.split(":")
    return datetime.fromisoformat(d) + timedelta(hours=int(h), minutes=int(m), seconds=float(sec))
def model(path):
    cols, rows = None, []
    for line in open(path):
        if line.startswith("#iPoint"):
            cols = {c: i for i, c in enumerate(line[1:].split())}
        elif line.strip() and not line.startswith("#"):
            rows.append(line.split())
    return cols, rows
C, new = model(run + "/base/COPS_model.dat")
C0, ref = model(old + "/base/COPS_model.dat")
CL, lat = model(late + "/COPS_model.dat")
assert C == C0 == CL, "column layout differs"
pts = [l.split() for l in open(table) if not l.startswith("#")]
bad = []
def check(ok, msg):
    if not ok:
        bad.append(msg)
check(len(new) == len(pts) == 1551, "points: model %d, table %d" % (len(new), len(pts)))
g = lambda r, c: r[C[c]]
obs = [(j, int(p[3])) for j, p in enumerate(pts) if p[2] == "obs"]
grid = [j for j, p in enumerate(pts) if p[2] == "grid"]
check(len(obs) == 1150 and len(grid) == 401 and [i for _, i in obs] == list(range(1150)), "obs/grid split")
for j, p in enumerate(pts):
    r = new[j]
    check(int(g(r, "iPoint")) == j, "iPoint %s at row %d" % (g(r, "iPoint"), j))
    check(abs((utc(g(r, "UTC")) - utc(p[1])).total_seconds()) < 1e-3, "UTC row %d: %s vs %s" % (j, g(r, "UTC"), p[1]))
for j, i in obs:
    for c in ("refNG", "refRG", "useNG", "useRG"):
        check(float(g(new[j], c)) == float(g(ref[i], c)), "%s of original %d: %s vs %s" % (c, i, g(new[j], c), g(ref[i], c)))
    check(g(new[j], "UTC") == g(ref[i], "UTC"), "UTC of original %d" % i)
for j in grid:
    for c in ("refNG", "refRG", "useNG", "useRG"):
        check(float(g(new[j], c)) == 0.0, "%s of grid row %d = %s" % (c, j, g(new[j], c)))
check(all(g(r, "inside") == "0" for r in new), "a point inside the nucleus")
s = lambda rows, c, idx: sum(int(rows[k][C[c]]) for k in idx)
print("points %d (obs %d, grid %d); useNG %d useRG %d (fitting base %d %d); grid useNG+useRG %d" % (
    len(new), len(obs), len(grid), s(new, "useNG", range(len(new))), s(new, "useRG", range(len(new))),
    s(ref, "useNG", range(len(ref))), s(ref, "useRG", range(len(ref))),
    s(new, "useNG", grid) + s(new, "useRG", grid)))
check(s(new, "useNG", range(len(new))) == s(ref, "useNG", range(len(ref))) and
      s(new, "useRG", range(len(new))) == s(ref, "useRG", range(len(ref))), "use counts differ")
d = sorted(float(g(new[j], "altitude")) - float(g(ref[i], "altitude")) for j, i in obs)
print("altitude new - fitting base (dz 0) at the %d originals: min %.3f median %.3f max %.3f m" % (len(d), d[0], d[len(d) // 2], d[-1]))
jl, il = obs[-1]
print("last original %s: altitude %s m (fitting base %s m); nVisible %s (%s); pNG %s pRG %s Pa (fitting base %s %s)" % (
    g(new[jl], "UTC"), g(new[jl], "altitude"), g(ref[il], "altitude"), g(new[jl], "nVisible"), g(ref[il], "nVisible"),
    g(new[jl], "pNG"), g(new[jl], "pRG"), g(ref[il], "pNG"), g(ref[il], "pRG")))
check(abs(float(g(new[jl], "altitude")) - 5.47) < 0.05, "last original altitude %s, expected about 5.47 m" % g(new[jl], "altitude"))
by = {g(r, "UTC"): r for r in new}
worst, same, n = {}, 0, 0
for r in lat:
    q = by.get(r[C["UTC"]])
    if q is None:
        check(False, "122289 dz1.5 time %s missing" % r[C["UTC"]]); continue
    n += 1; same += r[2:] == q[2:]
    for c in ("altitude", "nVisible", "pNG", "pRG", "pNG_H2O", "pRG_H2O", "pNG_CO2", "pRG_CO2"):
        a, b = float(r[C[c]]), float(q[C[c]])
        worst[c] = max(worst.get(c, 0.0), abs(a - b) / max(abs(a), abs(b), 1e-300))
print("vs job 122289 base_dz1.5 at its %d points: %d rows identical in all columns but iPoint/UTC; max relative difference %s" % (
    n, same, ", ".join("%s %.2e" % kv for kv in worst.items())))
check(n == 496, "122289 common points %d" % n)
m = json.load(open(run + "/base/run_manifest.json"))
k = m["outputs"]["kernels"]
print("run_manifest: states %s; kernels count %d records %d size %.1f GB set %s" % (
    m["states"], k["count"], k["records"], k["size_bytes"] / 1e9, k["set_sha256"][:16]))
check(m["states"].get("process") == "complete" and m["states"].get("cache") == "verified" and k["count"] == 1551, "run_manifest")
if bad:
    print("FAILED (%d):" % len(bad)); print("\n".join(bad[:40])); sys.exit(1)
print("base verification: OK")
EOF
rec $PY "$P/verify_base.py" "$RUN" "$R" "$LATE" "$RUN/late_grid_points.txt"
$PY "$P/verify_base.py" "$RUN" "$R" "$LATE" "$RUN/late_grid_points.txt" > "$P/verify_base.log" 2>&1; rc=$?
cat "$P/verify_base.log"; [ $rc -eq 0 ] || { echo "$(now) base verification failed" >&2; exit 7; }
[ "$STOP_AFTER" = base ] && { echo "$(now) STOP_AFTER=base: done"; exit 0; }

# 5. store (the runner's store stage command; cwd as in stamps/store.json)
if store_complete; then
    echo "$(now) store: complete, kept ($S)"
else
    [ -e "$S" ] && { mv "$S" "$S.incomplete_$(stamp)"; echo "$(now) moved an incomplete store aside: $S.incomplete_*"; }
    [ "$(freegb)" -ge "$NEED_STORE_GB" ] || { echo "$(now) refused: $(freegb) GB free, the store needs about $NEED_STORE_GB GB" >&2; exit 6; }
    echo "cd $APP" >> "$P/commands.sh"; cd "$APP" || exit 1
    rec bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ $PY "$APP/tools/c15_record_store.py" --run "$B" --out "$S" --threads $STORE_THREADS
    echo "$(now) store: start ($(freegb) GB free)"
    bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ $PY "$APP/tools/c15_record_store.py" --run "$B" --out "$S" --threads $STORE_THREADS > "$P/store.log" 2>&1; rc=$?
    echo "$(now) store: exit $rc"; tail -3 "$P/store.log" | cut -c1-300
    [ $rc -eq 0 ] && store_complete || { echo "$(now) store build failed" >&2; exit 8; }
fi
$PY - "$S/store_manifest.json" "$B/run_manifest.json" "$R/store/store_manifest.json" <<'EOF' 2>&1 | tee "$P/verify_store.log"
import json, sys
s, m, o = (json.load(open(p)) for p in sys.argv[1:4])
k = m["outputs"]["kernels"]
print("store: %d points, %d records (base %d), %.1f GB FP32, kernel set %s verified against run_manifest %s, mesh %s (12M store %s), "
      "nfaces %d, max relative rounding %.3g" % (len(s["points"]), s["records"], k["records"], s["records"] * 20 / 1e9,
      s["kernel_set_sha256"][:16], s["kernel_set_verified_against_run_manifest"], s["mesh_identity"][:12], o["mesh_identity"][:12],
      s["nfaces"], s["conversion"]["max_relative_rounding"]))
ok = (len(s["points"]) == 1551 and s["records"] == k["records"] and s["kernel_set_sha256"] == k["set_sha256"]
      and s["kernel_set_verified_against_run_manifest"] and s["mesh_identity"] == o["mesh_identity"] and s["nfaces"] == o["nfaces"])
print("store verification: %s" % ("OK" if ok else "FAILED")); sys.exit(0 if ok else 1)
EOF
rc=${PIPESTATUS[0]}
du -sh "$B" "$S" | sed 's/^/  /'; df -h "$COPS_OUTPUT_ROOT" > "$P/df_at_end.txt"; free -g > "$P/free_at_end.txt"
echo "$(now) done: exit $rc"
exit $rc
