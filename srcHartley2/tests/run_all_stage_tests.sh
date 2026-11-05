#!/usr/bin/env bash
# Run the complete focused and archived Hartley2 chemistry-development test
# stack for Stages 0 through 8.
#
# The lower-stage integrated tests intentionally consume frozen logs: their
# own runners define those logs as the authoritative, reproducible runtime
# protocol.  Stages 4 through 8 additionally support fresh MPI matrices when
# the caller supplies the appropriate executable and a new output root.  This
# orchestrator never configures or compiles AMPS, never edits an input, and
# never removes a build tree.  Keeping build mutation outside the test process
# guarantees that a failed fresh run cannot leave the canonical workspace
# correction-on.
#
# All optional inputs of the individual runners are made mandatory here.  A
# complete all-stage pass must not hide missing production/runtime evidence
# behind one of the per-stage runners' deliberate standalone SKIP paths.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source_root=$(cd -- "$script_dir/../.." && pwd)

evidence_root="$source_root/validation_runs"
stage4_mode=archived
stage4_executable=
stage4_output_root=
stage5_mode=archived
stage5_executable=
stage5_output_root=
stage6_mode=archived
stage6_executable=
stage6_output_root=
stage7_mode=archived
stage7_executable=
stage7_output_root=
stage8_mode=archived
stage8_executable=
stage8_output_root=
mpi_launcher=mpiexec
validation_particles=200000
chemistry_validation_events=200000
plan_only=false
requested_log_root=

usage() {
  cat <<'EOF'
usage: run_all_stage_tests.sh [options]

Run every focused/pure test and every archived runtime-evidence checker used
to close Hartley2 chemistry-development Stages 0 through 8.

Options:
  --evidence-root PATH
      Archive bundle root. Default: REPOSITORY/validation_runs. The expected
      dated subdirectory layout is documented below and preflighted strictly.
  --stage4-mode archived|fresh
      archived (default) revalidates the frozen correction-on 1/4/12 matrix.
      fresh invokes run_stage4_integrated_tests.sh and then compares its
      all-c=1 arm with the frozen correction-off reference.
  --stage4-executable PATH
      Required in fresh mode. Must be an already-built correction-on general-
      NTC executable. This runner deliberately does not compile it.
  --stage4-output-root PATH
      Required in fresh mode and must not exist. The Stage-4 runtime runner
      creates this evidence archive.
  --stage5-mode archived|fresh
      archived (default) revalidates the frozen Stage-5 correction-window and
      exact-chemistry 1/4/12 matrix. fresh invokes the Stage-5 MPI runner.
  --stage5-executable PATH
      Required in fresh Stage-5 mode. Must be an already-built DIM1,
      correction-window-enabled executable.
  --stage5-output-root PATH
      Required in fresh Stage-5 mode and must not exist.
  --stage6-mode archived|fresh
      archived (default) revalidates the frozen Stage-6 packet, Stage-5
      correction-window, and exact-chemistry 1/4/12 matrix. fresh invokes the
      Stage-6 MPI runner with a caller-supplied correction/packet-on binary.
  --stage6-executable PATH
      Required in fresh Stage-6 mode. Must be an already-built DIM1 binary
      with correction windows, Bernoulli residuals, and packets enabled.
  --stage6-output-root PATH
      Required in fresh Stage-6 mode and must not exist.
  --stage7-mode archived|fresh
      archived (default) revalidates the frozen Stage-7 channel-importance,
      Stage-6 packet, Stage-5 correction-window, and exact-chemistry 1/4/12
      matrix. fresh invokes the Stage-7 MPI runner.
  --stage7-executable PATH
      Required in fresh Stage-7 mode. Must be an already-built DIM1 binary
      with correction windows, Bernoulli residuals, packets, diagnostics, and
      channel importance enabled.
  --stage7-output-root PATH
      Required in fresh Stage-7 mode and must not exist.
  --stage8-mode archived|fresh
      archived (default) revalidates the frozen Stage-8 deterministic forced-
      segment matrix and its Stage-4/5/6/7 dependencies. fresh invokes the
      Stage-8 MPI runner.
  --stage8-executable PATH
      Required in fresh Stage-8 mode. Must be an already-built DIM1 binary
      with corrected chemistry, windows, packets, channel importance, and
      deterministic two-branch total-reaction forcing enabled.
  --stage8-output-root PATH
      Required in fresh Stage-8 mode and must not exist.
  --mpiexec COMMAND
      MPI launcher for Stage-0 focused and fresh Stage-4/5/6/7/8 tests.
  --validation-particles N
      Fresh Stage-4 and Stage-8 dependency-replay mover particles per case.
      Default: 200000.
  --chemistry-validation-events N
      Fresh Stage-4/5/6/7/8 chemistry events per case. Default: 200000.
  --log-root PATH
      Preserve every per-test log and summary under a new directory. Without
      this option, successful-run logs are temporary; failed/skipped-run logs
      are retained automatically and their absolute paths are reported.
  --plan
      Preflight all tools/artifacts and print the resolved test plan without
      executing a compiler, Python/Perl test, or MPI process.
  -h, --help
      Show this help.

Required evidence layout beneath --evidence-root:
  halley1d_16x_particles_40000iter_sampling10000_20261104_18mpi
  stage0_closure_20261105
  stage1_population_policy_20261105
  stage2_component_carrier_20261105
  stage3_systematic_20261105
  stage4-correction-on-full-20261105
  stage4-correction-off-reference-20261105
  stage4-correction-off-regression-20261105
  stage4-unsupported-mf-20261105b
  stage5-correction-window-full-20261105
  stage6-daughter-packet-full-20261105-r4
  stage7-channel-importance-full-20261105-r2
  stage8-total-reaction-forcing-full-20261105-r4

The current repository build is used only as a read-only Stage-3 parser
fixture and must be configured with individual corrections off. Rebuild it
before invoking this runner if that precondition is not satisfied; as always,
remove exactly REPOSITORY/build before recompiling.
EOF
}

while (( $# > 0 )); do
  case "$1" in
    --evidence-root)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      evidence_root=$2
      shift 2
      ;;
    --stage4-mode)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      stage4_mode=$2
      shift 2
      ;;
    --stage4-executable)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      stage4_executable=$2
      shift 2
      ;;
    --stage4-output-root)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      stage4_output_root=$2
      shift 2
      ;;
    --stage5-mode)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      stage5_mode=$2
      shift 2
      ;;
    --stage5-executable)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      stage5_executable=$2
      shift 2
      ;;
    --stage5-output-root)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      stage5_output_root=$2
      shift 2
      ;;
    --stage6-mode)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      stage6_mode=$2
      shift 2
      ;;
    --stage6-executable)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      stage6_executable=$2
      shift 2
      ;;
    --stage6-output-root)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      stage6_output_root=$2
      shift 2
      ;;
    --stage7-mode)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      stage7_mode=$2
      shift 2
      ;;
    --stage7-executable)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      stage7_executable=$2
      shift 2
      ;;
    --stage7-output-root)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      stage7_output_root=$2
      shift 2
      ;;
    --stage8-mode)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      stage8_mode=$2
      shift 2
      ;;
    --stage8-executable)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      stage8_executable=$2
      shift 2
      ;;
    --stage8-output-root)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      stage8_output_root=$2
      shift 2
      ;;
    --mpiexec)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      mpi_launcher=$2
      shift 2
      ;;
    --validation-particles)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      validation_particles=$2
      shift 2
      ;;
    --chemistry-validation-events)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      chemistry_validation_events=$2
      shift 2
      ;;
    --log-root)
      [[ $# -ge 2 ]] || { usage >&2; exit 2; }
      requested_log_root=$2
      shift 2
      ;;
    --plan)
      plan_only=true
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      printf 'unknown all-stage option: %s\n' "$1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ -n "$requested_log_root" ]]; then
  if [[ "$requested_log_root" != /* ]]; then
    requested_log_root="$source_root/$requested_log_root"
  fi
  if [[ -e "$requested_log_root" ]]; then
    printf 'all-stage log root already exists: %s\n' \
      "$requested_log_root" >&2
    exit 2
  fi
fi

case "$stage4_mode" in
  archived|fresh) ;;
  *)
    printf 'invalid --stage4-mode value: %s\n' "$stage4_mode" >&2
    exit 2
    ;;
esac
case "$stage5_mode" in
  archived|fresh) ;;
  *)
    printf 'invalid --stage5-mode value: %s\n' "$stage5_mode" >&2
    exit 2
    ;;
esac
case "$stage6_mode" in
  archived|fresh) ;;
  *)
    printf 'invalid --stage6-mode value: %s\n' "$stage6_mode" >&2
    exit 2
    ;;
esac
case "$stage7_mode" in
  archived|fresh) ;;
  *)
    printf 'invalid --stage7-mode value: %s\n' "$stage7_mode" >&2
    exit 2
    ;;
esac
case "$stage8_mode" in
  archived|fresh) ;;
  *)
    printf 'invalid --stage8-mode value: %s\n' "$stage8_mode" >&2
    exit 2
    ;;
esac
for numeric_option in validation_particles chemistry_validation_events; do
  numeric_value=${!numeric_option}
  if [[ ! "$numeric_value" =~ ^[1-9][0-9]*$ ]]; then
    printf 'invalid positive integer for %s: %s\n' \
      "$numeric_option" "$numeric_value" >&2
    exit 2
  fi
done

if [[ ! -d "$evidence_root" ]]; then
  printf 'evidence root is not a directory: %s\n' "$evidence_root" >&2
  exit 2
fi
evidence_root=$(cd -- "$evidence_root" && pwd -P)

if [[ "$stage4_mode" == fresh ]]; then
  if [[ -z "$stage4_executable" || -z "$stage4_output_root" ]]; then
    printf '%s\n' \
      'fresh Stage-4 mode requires --stage4-executable and --stage4-output-root' >&2
    exit 2
  fi
  if [[ ! -f "$stage4_executable" || ! -x "$stage4_executable" ]]; then
    printf 'Stage-4 executable is not executable: %s\n' \
      "$stage4_executable" >&2
    exit 2
  fi
  stage4_executable=$(readlink -f -- "$stage4_executable")
  if [[ "$stage4_output_root" != /* ]]; then
    stage4_output_root="$source_root/$stage4_output_root"
  fi
  if [[ -e "$stage4_output_root" ]]; then
    printf 'fresh Stage-4 output already exists: %s\n' \
      "$stage4_output_root" >&2
    exit 2
  fi
elif [[ -n "$stage4_executable" || -n "$stage4_output_root" ]]; then
  printf '%s\n' \
    '--stage4-executable/--stage4-output-root are valid only in fresh mode' >&2
  exit 2
fi

if [[ "$stage5_mode" == fresh ]]; then
  if [[ -z "$stage5_executable" || -z "$stage5_output_root" ]]; then
    printf '%s\n' \
      'fresh Stage-5 mode requires --stage5-executable and --stage5-output-root' >&2
    exit 2
  fi
  if [[ ! -f "$stage5_executable" || ! -x "$stage5_executable" ]]; then
    printf 'Stage-5 executable is not executable: %s\n' \
      "$stage5_executable" >&2
    exit 2
  fi
  stage5_executable=$(readlink -f -- "$stage5_executable")
  if [[ "$stage5_output_root" != /* ]]; then
    stage5_output_root="$source_root/$stage5_output_root"
  fi
  if [[ -e "$stage5_output_root" ]]; then
    printf 'fresh Stage-5 output already exists: %s\n' \
      "$stage5_output_root" >&2
    exit 2
  fi
elif [[ -n "$stage5_executable" || -n "$stage5_output_root" ]]; then
  printf '%s\n' \
    '--stage5-executable/--stage5-output-root are valid only in fresh mode' >&2
  exit 2
fi

if [[ "$stage6_mode" == fresh ]]; then
  if [[ -z "$stage6_executable" || -z "$stage6_output_root" ]]; then
    printf '%s\n' \
      'fresh Stage-6 mode requires --stage6-executable and --stage6-output-root' >&2
    exit 2
  fi
  if [[ ! -f "$stage6_executable" || ! -x "$stage6_executable" ]]; then
    printf 'Stage-6 executable is not executable: %s\n' \
      "$stage6_executable" >&2
    exit 2
  fi
  stage6_executable=$(readlink -f -- "$stage6_executable")
  if [[ "$stage6_output_root" != /* ]]; then
    stage6_output_root="$source_root/$stage6_output_root"
  fi
  if [[ -e "$stage6_output_root" ]]; then
    printf 'fresh Stage-6 output already exists: %s\n' \
      "$stage6_output_root" >&2
    exit 2
  fi
elif [[ -n "$stage6_executable" || -n "$stage6_output_root" ]]; then
  printf '%s\n' \
    '--stage6-executable/--stage6-output-root are valid only in fresh mode' >&2
  exit 2
fi

if [[ "$stage7_mode" == fresh ]]; then
  if [[ -z "$stage7_executable" || -z "$stage7_output_root" ]]; then
    printf '%s\n' \
      'fresh Stage-7 mode requires --stage7-executable and --stage7-output-root' >&2
    exit 2
  fi
  if [[ ! -f "$stage7_executable" || ! -x "$stage7_executable" ]]; then
    printf 'Stage-7 executable is not executable: %s\n' \
      "$stage7_executable" >&2
    exit 2
  fi
  stage7_executable=$(readlink -f -- "$stage7_executable")
  if [[ "$stage7_output_root" != /* ]]; then
    stage7_output_root="$source_root/$stage7_output_root"
  fi
  if [[ -e "$stage7_output_root" ]]; then
    printf 'fresh Stage-7 output already exists: %s\n' \
      "$stage7_output_root" >&2
    exit 2
  fi
elif [[ -n "$stage7_executable" || -n "$stage7_output_root" ]]; then
  printf '%s\n' \
    '--stage7-executable/--stage7-output-root are valid only in fresh mode' >&2
  exit 2
fi

if [[ "$stage8_mode" == fresh ]]; then
  if [[ -z "$stage8_executable" || -z "$stage8_output_root" ]]; then
    printf '%s\n' \
      'fresh Stage-8 mode requires --stage8-executable and --stage8-output-root' >&2
    exit 2
  fi
  if [[ ! -f "$stage8_executable" || ! -x "$stage8_executable" ]]; then
    printf 'Stage-8 executable is not executable: %s\n' \
      "$stage8_executable" >&2
    exit 2
  fi
  stage8_executable=$(readlink -f -- "$stage8_executable")
  if [[ "$stage8_output_root" != /* ]]; then
    stage8_output_root="$source_root/$stage8_output_root"
  fi
  if [[ -e "$stage8_output_root" ]]; then
    printf 'fresh Stage-8 output already exists: %s\n' \
      "$stage8_output_root" >&2
    exit 2
  fi
elif [[ -n "$stage8_executable" || -n "$stage8_output_root" ]]; then
  printf '%s\n' \
    '--stage8-executable/--stage8-output-root are valid only in fresh mode' >&2
  exit 2
fi

resolved_mpi_launcher=$(command -v -- "$mpi_launcher" || true)
if [[ -z "$resolved_mpi_launcher" ]]; then
  printf 'MPI launcher is not on PATH: %s\n' "$mpi_launcher" >&2
  exit 2
fi
resolved_mpi_launcher=$(readlink -f -- "$resolved_mpi_launcher")

# Every tool below is used by at least one child runner.  Preflight the full
# toolchain before the first test so a missing Perl/compiler cannot turn a
# long suite into a late partial result.
for required_tool in mpicxx c++ g++ perl python3 grep sha256sum tee; do
  if ! command -v -- "$required_tool" >/dev/null 2>&1; then
    printf 'required all-stage test tool is not on PATH: %s\n' \
      "$required_tool" >&2
    exit 2
  fi
done

stage0_baseline="$evidence_root/halley1d_16x_particles_40000iter_sampling10000_20261104_18mpi"
stage0_production="$evidence_root/stage0_closure_20261105/diagnostics_on_10000iter_14mpi/production"
stage1_root="$evidence_root/stage1_population_policy_20261105"
stage1_policy_log="$stage1_root/canonical_policy/run.log"
stage1_alternate_log="$stage1_root/alternate_policy/run.log"
stage1_cadence_log="$stage1_root/cadence/run.log"
stage1_legacy_baseline="$evidence_root/stage0_closure_20261105/diagnostics_off/production/PT/plots"
stage1_legacy_candidate="$stage1_root/legacy_candidate/PT/plots"
stage2_root="$evidence_root/stage2_component_carrier_20261105"
stage3_root="$evidence_root/stage3_systematic_20261105"
stage4_on_archive="$evidence_root/stage4-correction-on-full-20261105"
stage4_off_archive="$evidence_root/stage4-correction-off-reference-20261105"
stage4_off_regression="$evidence_root/stage4-correction-off-regression-20261105"
stage4_mf_archive="$evidence_root/stage4-unsupported-mf-20261105b"
stage5_archive="$evidence_root/stage5-correction-window-full-20261105"
stage6_archive="$evidence_root/stage6-daughter-packet-full-20261105-r4"
stage7_archive="$evidence_root/stage7-channel-importance-full-20261105-r2"
stage8_archive="$evidence_root/stage8-total-reaction-forcing-full-20261105-r4"

require_file() {
  if [[ ! -f "$1" ]]; then
    printf 'required all-stage evidence file is missing: %s\n' "$1" >&2
    exit 2
  fi
}

require_directory() {
  if [[ ! -d "$1" ]]; then
    printf 'required all-stage evidence directory is missing: %s\n' "$1" >&2
    exit 2
  fi
}

require_directory "$stage0_baseline/production"
require_directory "$stage0_production/PT/diagnostics"
require_file "$stage1_policy_log"
require_file "$stage1_alternate_log"
require_file "$stage1_cadence_log"
require_directory "$stage1_legacy_baseline"
require_directory "$stage1_legacy_candidate"
for stage_rank in 1 4 12; do
  require_file "$stage2_root/mpi${stage_rank}/run.log"
  require_file "$stage3_root/mpi${stage_rank}/stage3.mpi${stage_rank}.log"
  require_file \
    "$stage4_off_archive/mpi${stage_rank}/run.log"
  require_file \
    "$stage4_off_regression/mpi${stage_rank}/weight-adaptation/run.log"
  require_file \
    "$stage4_off_regression/mpi${stage_rank}/systematic-residual/run.log"
  if [[ "$stage4_mode" == archived ]]; then
    require_file \
      "$stage4_on_archive/mpi${stage_rank}/local-timestep/run.log"
    require_file \
      "$stage4_on_archive/mpi${stage_rank}/local-timestep-unit/run.log"
    require_file \
      "$stage4_on_archive/mpi${stage_rank}/weight-adaptation/run.log"
    require_file \
      "$stage4_on_archive/mpi${stage_rank}/chemistry/run.log"
  fi
  if [[ "$stage5_mode" == archived ]]; then
    require_file \
      "$stage5_archive/mpi${stage_rank}/correction-windows/run.log"
    require_file \
      "$stage5_archive/mpi${stage_rank}/chemistry/run.log"
  fi
  if [[ "$stage6_mode" == archived ]]; then
    require_file \
      "$stage6_archive/mpi${stage_rank}/correction-windows/run.log"
    require_file \
      "$stage6_archive/mpi${stage_rank}/chemistry/run.log"
    require_file \
      "$stage6_archive/mpi${stage_rank}/daughter-packets/run.log"
  fi
  if [[ "$stage7_mode" == archived ]]; then
    require_file \
      "$stage7_archive/mpi${stage_rank}/correction-windows/run.log"
    require_file \
      "$stage7_archive/mpi${stage_rank}/chemistry/run.log"
    require_file \
      "$stage7_archive/mpi${stage_rank}/daughter-packets/run.log"
    require_file \
      "$stage7_archive/mpi${stage_rank}/channel-importance/run.log"
    require_file \
      "$stage7_archive/mpi${stage_rank}/channel-importance/PT/diagnostics/comet_stage7_importance.out=$((7000000+stage_rank)).csv"
    require_file \
      "$stage7_archive/mpi${stage_rank}/channel-importance/PT/diagnostics/comet_stage7_importance.out=$((7000000+stage_rank)).json"
  fi
  if [[ "$stage8_mode" == archived ]]; then
    require_file \
      "$stage8_archive/mpi${stage_rank}/local-timestep/run.log"
    require_file \
      "$stage8_archive/mpi${stage_rank}/local-timestep-unit/run.log"
    require_file \
      "$stage8_archive/mpi${stage_rank}/weight-adaptation/run.log"
    require_file \
      "$stage8_archive/mpi${stage_rank}/correction-windows/run.log"
    require_file \
      "$stage8_archive/mpi${stage_rank}/chemistry/run.log"
    require_file \
      "$stage8_archive/mpi${stage_rank}/daughter-packets/run.log"
    require_file \
      "$stage8_archive/mpi${stage_rank}/channel-importance/run.log"
    require_file \
      "$stage8_archive/mpi${stage_rank}/channel-importance/PT/diagnostics/comet_stage7_importance.out=$((7000000+stage_rank)).csv"
    require_file \
      "$stage8_archive/mpi${stage_rank}/channel-importance/PT/diagnostics/comet_stage7_importance.out=$((7000000+stage_rank)).json"
    require_file \
      "$stage8_archive/mpi${stage_rank}/total-reaction-forcing/run.log"
  fi
done
require_file "$stage3_root/mpi1/PT/diagnostics/stage3_cells.C2.out=0.csv"
require_file "$stage3_root/mpi1/PT/diagnostics/cell_measures.csv"
require_file "$stage4_on_archive/provenance/SHA256SUMS"
require_file "$stage4_mf_archive/run.log"
require_file "$stage4_mf_archive/error.log"
if [[ "$stage4_mode" == archived ]]; then
  require_file \
    "$stage4_on_archive/mpi1/unsupported-simplified-ntc/run.log"
  require_file \
    "$stage4_on_archive/mpi1/unsupported-simplified-ntc/error.log"
fi
if [[ "$stage6_mode" == archived ]]; then
  require_file "$stage6_archive/provenance/SHA256SUMS"
fi
if [[ "$stage7_mode" == archived ]]; then
  require_file "$stage7_archive/provenance/SHA256SUMS"
fi
if [[ "$stage8_mode" == archived ]]; then
  require_file "$stage8_archive/provenance/SHA256SUMS"
fi

# The Stage-3 parser test copies generated sources from the current build but
# never modifies them.  Require the complete fixture and correction-off mode;
# using a correction-on Stage-4 build here would make an all-stage result
# dependent on the caller's previous local configuration.
require_directory "$source_root/build/main"
require_file "$source_root/build/pic/picGlobal.dfn"
require_file "$source_root/halley1d.input.Assembled.Block"
require_file "$source_root/halley1d.input.Assembled.Species"
if ! grep -Eq -- \
    '^#define _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ _INDIVIDUAL_PARTICLE_WEIGHT_OFF_$' \
    "$source_root/build/pic/picGlobal.dfn"; then
  printf '%s\n' \
    'current Stage-3 parser fixture is not a correction-off AMPS build' >&2
  exit 2
fi

print_plan() {
  printf 'Hartley2 all-stage test plan\n'
  printf '  repository: %s\n' "$source_root"
  printf '  evidence:   %s\n' "$evidence_root"
  printf '  MPI:        %s\n' "$resolved_mpi_launcher"
  printf '  Stage 4:    %s\n' "$stage4_mode"
  printf '  Stage 5:    %s\n' "$stage5_mode"
  printf '  Stage 6:    %s\n' "$stage6_mode"
  printf '  Stage 7:    %s\n' "$stage7_mode"
  printf '  Stage 8:    %s\n' "$stage8_mode"
  if [[ -n "$requested_log_root" ]]; then
    printf '  logs:       %s\n' "$requested_log_root"
  else
    printf '%s\n' \
      '  logs:       temporary on success; automatically retained on failure/skip'
  fi
  if [[ "$stage4_mode" == fresh ]]; then
    printf '  Stage-4 executable: %s\n' "$stage4_executable"
    printf '  Stage-4 output:     %s\n' "$stage4_output_root"
  fi
  if [[ "$stage5_mode" == fresh ]]; then
    printf '  Stage-5 executable: %s\n' "$stage5_executable"
    printf '  Stage-5 output:     %s\n' "$stage5_output_root"
  fi
  if [[ "$stage6_mode" == fresh ]]; then
    printf '  Stage-6 executable: %s\n' "$stage6_executable"
    printf '  Stage-6 output:     %s\n' "$stage6_output_root"
  fi
  if [[ "$stage7_mode" == fresh ]]; then
    printf '  Stage-7 executable: %s\n' "$stage7_executable"
    printf '  Stage-7 output:     %s\n' "$stage7_output_root"
  fi
  if [[ "$stage8_mode" == fresh ]]; then
    printf '  Stage-8 executable: %s\n' "$stage8_executable"
    printf '  Stage-8 output:     %s\n' "$stage8_output_root"
  fi
  printf '%s\n' \
    '  tests:' \
    '    Stage 0 focused MPI, analyzer, frozen baseline, production report' \
    '    Stage 1 parser/network, runtime policy, cadence, legacy equivalence' \
    '    Stage 2 pure planner and integrated 1/4/12 transaction evidence' \
    '    Stage 3 accumulator, parser, adaptive bins, diagnostic semantics,' \
    '            synthetic science gate, CLI smoke, integrated 1/4/12 evidence' \
    '    Stage 4 correction-off Stage-2/3 regressions, provenance, MF and' \
    '            simplified-NTC rejection, correction-on matrix, all-c=1 identity' \
    '    Stage 5 input/planner tests and correction-window/exact-chemistry' \
    '            integrated 1/4/12 evidence' \
    '    Stage 6 packet-policy/planner tests, Stage-5 dependency replay,' \
    '            production transaction, and integrated 1/4/12 evidence' \
    '    Stage 7 channel-policy/planner tests, Stage-4/5/6 dependency replay,' \
    '            weighted channel moments, sinks, ESS, durable diagnostics,' \
    '            and decomposition-identical integrated 1/4/12 evidence' \
    '    Stage 8 input/planner and analytic forced-segment references,' \
    '            strict parser mutation tests, Stage-4/5/6/7 dependency replay,' \
    '            and mover-integrated deterministic 1/4/12 evidence'
}

print_plan
if [[ "$plan_only" == true ]]; then
  printf '%s\n' 'PASS all-stage preflight and plan'
  exit 0
fi

all_stage_tmp=$(mktemp -d "${TMPDIR:-/tmp}/hartley2-all-stages.XXXXXX")
scratch_root="$all_stage_tmp/scratch"
mkdir -p -- "$scratch_root"
explicit_log_root=false
if [[ -n "$requested_log_root" ]]; then
  mkdir -p -- "$requested_log_root"
  log_root=$(cd -- "$requested_log_root" && pwd -P)
  explicit_log_root=true
else
  log_root="$all_stage_tmp/logs"
  mkdir -p -- "$log_root"
fi

# Keep default logs only when they are needed to diagnose an incomplete run.
# Explicit --log-root output is always retained.  On a signal or unexpected
# shell error suite_success remains false, so the default logs survive too.
suite_success=false
cleanup_all_stage() {
  rm -rf -- "$scratch_root"
  if [[ "$suite_success" == true || "$explicit_log_root" == true ]]; then
    rm -rf -- "$all_stage_tmp"
  else
    printf 'All-stage diagnostic logs preserved at: %s\n' "$log_root" >&2
  fi
}
trap cleanup_all_stage EXIT

test_names=()
test_stages=()
test_results=()
test_logs=()
test_exit_codes=()
passed_tests=0
skipped_tests=0
failed_tests=0
executed_tests=0

run_test() {
  local stage=$1
  local label=$2
  shift 2
  local started=$SECONDS
  local test_log="$log_root/${label//[^A-Za-z0-9_.-]/_}.log"
  local command_status=0
  local result=
  executed_tests=$((executed_tests+1))
  printf 'START all-stage test=%s stage=%s log=%s\n' \
    "$label" "$stage" "$test_log"
  # Preserve live child output while retaining a temporary receipt.  A child
  # runner's standalone SKIP is not a pass for this complete orchestrator,
  # even if that child intentionally returns zero for partial checkouts.  Do
  # not fail fast: independent later groups should still run and contribute to
  # the final diagnostic summary.
  if "$@" 2>&1 | tee "$test_log"; then
    command_status=0
  else
    command_status=$?
  fi

  if (( command_status != 0 )); then
    result=FAIL
    failed_tests=$((failed_tests+1))
  elif grep -Eq -- '^SKIP ' "$test_log"; then
    result=SKIP
    skipped_tests=$((skipped_tests+1))
  else
    result=PASS
    passed_tests=$((passed_tests+1))
  fi

  test_names+=("$label")
  test_stages+=("$stage")
  test_results+=("$result")
  test_logs+=("$test_log")
  test_exit_codes+=("$command_status")

  printf '%s all-stage test=%s stage=%s exit=%d elapsed-seconds=%d log=%s\n' \
    "$result" "$label" "$stage" "$command_status" \
    "$((SECONDS-started))" "$test_log" | tee -a "$test_log"
}

run_stage3_adaptive_smoke() {
  local output_prefix="$scratch_root/stage3-adaptive-smoke"
  local input_table=
  input_table="$stage3_root/mpi1/PT/diagnostics/stage3_cells.C2.out=0.csv"
  local sampling_length=
  # The focused Stage-3 archive and a production science archive need not use
  # the same sampling window.  Read the value embedded in the authoritative
  # table instead of importing the production example's value into this test.
  sampling_length=$(PYTHONDONTWRITEBYTECODE=1 python3 -c \
    'import csv,sys; print(next(csv.DictReader(open(sys.argv[1], newline="")))["sampling_length"])' \
    "$input_table")
  if [[ ! "$sampling_length" =~ ^[1-9][0-9]*$ ]]; then
    printf 'invalid Stage-3 archived sampling length: %s\n' \
      "$sampling_length" >&2
    return 1
  fi
  PYTHONDONTWRITEBYTECODE=1 python3 \
    "$source_root/srcHartley2/visualization/adaptive_radial_bins.py" \
    "$input_table" \
    --cell-measures "$stage3_root/mpi1/PT/diagnostics/cell_measures.csv" \
    --sampling-length "$sampling_length" --minimum-exposure 100000 \
    --output "$output_prefix"
  [[ -s "$output_prefix.csv" && -s "$output_prefix.json" ]]
}

check_stage4_provenance() {
  (
    cd -- "$stage4_on_archive/provenance"
    sha256sum -c SHA256SUMS
  )
}

check_rejection_message() {
  local label=$1
  local expected=$2
  shift 2
  if ! grep -Fq -- "$expected" "$@"; then
    printf 'missing %s rejection message: %s\n' "$label" "$expected" >&2
    return 1
  fi
}

check_stage4_mf_rejection() {
  check_rejection_message MF \
    'individual corrections with collisions require the validated general NTC model' \
    "$stage4_mf_archive/run.log" "$stage4_mf_archive/error.log"
}

check_stage4_simplified_rejection() {
  check_rejection_message simplified-NTC \
    'simplified NTC does not support individual particle-weight corrections' \
    "$stage4_on_root/mpi1/unsupported-simplified-ntc/run.log" \
    "$stage4_on_root/mpi1/unsupported-simplified-ntc/error.log"
}

check_stage4_matrix_and_identity() {
  local stage4_arguments=()
  for stage_rank in 1 4 12; do
    stage4_arguments+=(
      --local-log "$stage_rank" \
        "$stage4_on_root/mpi${stage_rank}/local-timestep/run.log"
      --unit-local-log "$stage_rank" \
        "$stage4_on_root/mpi${stage_rank}/local-timestep-unit/run.log"
      --off-unit-local-log "$stage_rank" \
        "$stage4_off_archive/mpi${stage_rank}/run.log"
      --adaptation-log "$stage_rank" \
        "$stage4_on_root/mpi${stage_rank}/weight-adaptation/run.log"
      --chemistry-log "$stage_rank" \
        "$stage4_on_root/mpi${stage_rank}/chemistry/run.log"
    )
  done
  PYTHONDONTWRITEBYTECODE=1 python3 \
    "$script_dir/test_stage4_integrated_logs.py" \
    "${stage4_arguments[@]}"
}

check_stage5_matrix() {
  local stage5_arguments=()
  for stage_rank in 1 4 12; do
    stage5_arguments+=(
      --window-log "$stage_rank" \
        "$stage5_on_root/mpi${stage_rank}/correction-windows/run.log"
      --chemistry-log "$stage_rank" \
        "$stage5_on_root/mpi${stage_rank}/chemistry/run.log"
    )
  done
  PYTHONDONTWRITEBYTECODE=1 python3 \
    "$script_dir/test_stage5_integrated_logs.py" \
    "${stage5_arguments[@]}"
}

check_stage6_matrix() {
  local stage5_arguments=()
  local packet_arguments=(--expected-events "$chemistry_validation_events")
  for stage_rank in 1 4 12; do
    # A Stage-6 archive is accepted only if its prerequisite Stage-5 windows
    # and exact-chemistry validator still pass under the same frozen binary.
    stage5_arguments+=(
      --window-log "$stage_rank" \
        "$stage6_on_root/mpi${stage_rank}/correction-windows/run.log"
      --chemistry-log "$stage_rank" \
        "$stage6_on_root/mpi${stage_rank}/chemistry/run.log"
    )
    packet_arguments+=(
      --packet-log "$stage_rank" \
        "$stage6_on_root/mpi${stage_rank}/daughter-packets/run.log"
    )
  done
  PYTHONDONTWRITEBYTECODE=1 python3 \
    "$script_dir/test_stage5_integrated_logs.py" \
    "${stage5_arguments[@]}" &&
  PYTHONDONTWRITEBYTECODE=1 python3 \
    "$script_dir/test_stage6_daughter_packet_logs.py" \
    "${packet_arguments[@]}"
}

check_stage6_provenance() {
  (
    cd -- "$stage6_archive/provenance"
    sha256sum -c SHA256SUMS
  )
}

check_stage7_matrix() {
  local stage5_arguments=()
  local packet_arguments=(--expected-events "$chemistry_validation_events")
  local importance_arguments=(--expected-events "$chemistry_validation_events")
  for stage_rank in 1 4 12; do
    # Revalidate the complete correction/packet dependency stack from the
    # same frozen Stage-7 binary before accepting its importance estimator.
    stage5_arguments+=(
      --window-log "$stage_rank" \
        "$stage7_on_root/mpi${stage_rank}/correction-windows/run.log"
      --chemistry-log "$stage_rank" \
        "$stage7_on_root/mpi${stage_rank}/chemistry/run.log"
    )
    packet_arguments+=(
      --packet-log "$stage_rank" \
        "$stage7_on_root/mpi${stage_rank}/daughter-packets/run.log"
    )
    importance_arguments+=(
      --importance-log "$stage_rank" \
        "$stage7_on_root/mpi${stage_rank}/channel-importance/run.log"
      --diagnostics-root "$stage_rank" \
        "$stage7_on_root/mpi${stage_rank}/channel-importance"
    )
  done
  PYTHONDONTWRITEBYTECODE=1 python3 \
    "$script_dir/test_stage5_integrated_logs.py" \
    "${stage5_arguments[@]}" &&
  PYTHONDONTWRITEBYTECODE=1 python3 \
    "$script_dir/test_stage6_daughter_packet_logs.py" \
    "${packet_arguments[@]}" &&
  PYTHONDONTWRITEBYTECODE=1 python3 \
    "$script_dir/test_stage7_channel_importance_logs.py" \
    "${importance_arguments[@]}"
}

check_stage7_provenance() {
  (
    cd -- "$stage7_archive/provenance"
    sha256sum -c SHA256SUMS
  )
}

check_stage8_matrix() {
  local stage4_arguments=()
  local stage5_arguments=()
  local packet_arguments=(--expected-events "$chemistry_validation_events")
  local importance_arguments=(--expected-events "$chemistry_validation_events")
  local forcing_arguments=(--expected-events "$chemistry_validation_events")
  for stage_rank in 1 4 12; do
    # Stage 8 changes the mover transaction, so accept its evidence only when
    # the same frozen binary still satisfies every corrected-chemistry layer.
    # The Stage-4 replay exercises the callback-enabled first-order mover's
    # ordinary local/MPI face, continuation, correction, collision, and
    # adaptation paths.  Stage-7 sidecars then prove that the reacted branch
    # continues to transport L, while Stage 6 proves the subsequent /N.
    stage4_arguments+=(
      --local-log "$stage_rank" \
        "$stage8_on_root/mpi${stage_rank}/local-timestep/run.log"
      --unit-local-log "$stage_rank" \
        "$stage8_on_root/mpi${stage_rank}/local-timestep-unit/run.log"
      --adaptation-log "$stage_rank" \
        "$stage8_on_root/mpi${stage_rank}/weight-adaptation/run.log"
      --chemistry-log "$stage_rank" \
        "$stage8_on_root/mpi${stage_rank}/chemistry/run.log"
    )
    stage5_arguments+=(
      --window-log "$stage_rank" \
        "$stage8_on_root/mpi${stage_rank}/correction-windows/run.log"
      --chemistry-log "$stage_rank" \
        "$stage8_on_root/mpi${stage_rank}/chemistry/run.log"
    )
    packet_arguments+=(
      --packet-log "$stage_rank" \
        "$stage8_on_root/mpi${stage_rank}/daughter-packets/run.log"
    )
    importance_arguments+=(
      --importance-log "$stage_rank" \
        "$stage8_on_root/mpi${stage_rank}/channel-importance/run.log"
      --diagnostics-root "$stage_rank" \
        "$stage8_on_root/mpi${stage_rank}/channel-importance"
    )
    forcing_arguments+=(
      --forcing-log "$stage_rank" \
        "$stage8_on_root/mpi${stage_rank}/total-reaction-forcing/run.log"
    )
  done
  PYTHONDONTWRITEBYTECODE=1 python3 \
    "$script_dir/test_stage4_integrated_logs.py" \
    "${stage4_arguments[@]}" &&
  PYTHONDONTWRITEBYTECODE=1 python3 \
    "$script_dir/test_stage5_integrated_logs.py" \
    "${stage5_arguments[@]}" &&
  PYTHONDONTWRITEBYTECODE=1 python3 \
    "$script_dir/test_stage6_daughter_packet_logs.py" \
    "${packet_arguments[@]}" &&
  PYTHONDONTWRITEBYTECODE=1 python3 \
    "$script_dir/test_stage7_channel_importance_logs.py" \
    "${importance_arguments[@]}" &&
  PYTHONDONTWRITEBYTECODE=1 python3 \
    "$script_dir/test_stage8_total_reaction_forcing_logs.py" \
    "${forcing_arguments[@]}"
}

check_stage8_provenance() {
  (
    cd -- "$stage8_archive/provenance"
    sha256sum -c SHA256SUMS
  )
}

run_test stage0 stage0.complete env \
  HARTLEY2_MPIEXEC="$resolved_mpi_launcher" \
  HARTLEY2_STAGE0_BASELINE_ROOT="$stage0_baseline" \
  HARTLEY2_STAGE0_PRODUCTION_ROOT="$stage0_production" \
  "$script_dir/run_stage0_tests.sh"

run_test stage1 stage1.complete env \
  HARTLEY2_STAGE1_POLICY_LOG="$stage1_policy_log" \
  HARTLEY2_STAGE1_ALTERNATE_LOG="$stage1_alternate_log" \
  HARTLEY2_STAGE1_CADENCE_LOG="$stage1_cadence_log" \
  HARTLEY2_STAGE1_LEGACY_BASELINE="$stage1_legacy_baseline" \
  HARTLEY2_STAGE1_LEGACY_CANDIDATE="$stage1_legacy_candidate" \
  "$script_dir/run_stage1_tests.sh"

run_test stage2 stage2.planner "$script_dir/run_stage2_planner_tests.sh"
run_test stage2 stage2.integrated env \
  HARTLEY2_STAGE2_LOG_1="$stage2_root/mpi1/run.log" \
  HARTLEY2_STAGE2_LOG_4="$stage2_root/mpi4/run.log" \
  HARTLEY2_STAGE2_LOG_12="$stage2_root/mpi12/run.log" \
  "$script_dir/run_stage2_integrated_tests.sh"

run_test stage3 stage3.accumulator "$script_dir/run_stage3_accumulator_tests.sh"
run_test stage3 stage3.adaptive-unit env PYTHONDONTWRITEBYTECODE=1 \
  python3 "$script_dir/test_stage3_adaptive_bins.py"
run_test stage3 stage3.diagnostic-semantics env PYTHONDONTWRITEBYTECODE=1 \
  python3 "$script_dir/test_stage3_diagnostic_semantics.py"
run_test stage3 stage3.parser "$script_dir/run_stage3_parser_tests.sh"
run_test stage3 stage3.science-gate "$script_dir/run_stage3_science_gate_tests.sh"
run_test stage3 stage3.adaptive-cli run_stage3_adaptive_smoke
run_test stage3 stage3.integrated env \
  HARTLEY2_STAGE3_LOG_1="$stage3_root/mpi1/stage3.mpi1.log" \
  HARTLEY2_STAGE3_LOG_4="$stage3_root/mpi4/stage3.mpi4.log" \
  HARTLEY2_STAGE3_LOG_12="$stage3_root/mpi12/stage3.mpi12.log" \
  "$script_dir/run_stage3_integrated_tests.sh"

# These two replays prove that the final Stage-4 source retained correction-off
# Stage-2/3 behavior; the earlier original Stage-2/3 archives above remain
# useful independent implementation evidence and are deliberately also run.
run_test stage4 stage4.correction-off-stage2-regression env \
  HARTLEY2_STAGE2_LOG_1="$stage4_off_regression/mpi1/weight-adaptation/run.log" \
  HARTLEY2_STAGE2_LOG_4="$stage4_off_regression/mpi4/weight-adaptation/run.log" \
  HARTLEY2_STAGE2_LOG_12="$stage4_off_regression/mpi12/weight-adaptation/run.log" \
  "$script_dir/run_stage2_integrated_tests.sh"
run_test stage4 stage4.correction-off-stage3-regression env \
  HARTLEY2_STAGE3_LOG_1="$stage4_off_regression/mpi1/systematic-residual/run.log" \
  HARTLEY2_STAGE3_LOG_4="$stage4_off_regression/mpi4/systematic-residual/run.log" \
  HARTLEY2_STAGE3_LOG_12="$stage4_off_regression/mpi12/systematic-residual/run.log" \
  "$script_dir/run_stage3_integrated_tests.sh"
run_test stage4 stage4.provenance check_stage4_provenance
run_test stage4 stage4.mf-rejection check_stage4_mf_rejection

if [[ "$stage4_mode" == fresh ]]; then
  run_test stage4 stage4.fresh-runtime \
    "$script_dir/run_stage4_integrated_tests.sh" \
    --executable "$stage4_executable" \
    --output-root "$stage4_output_root" \
    --mpiexec "$resolved_mpi_launcher" \
    --validation-particles "$validation_particles" \
    --chemistry-validation-events "$chemistry_validation_events"
  stage4_on_root=$stage4_output_root
else
  stage4_on_root=$stage4_on_archive
fi

run_test stage4 stage4.simplified-ntc-rejection check_stage4_simplified_rejection
run_test stage4 stage4.matrix-and-identity check_stage4_matrix_and_identity

run_test stage5 stage5.input-and-planner \
  "$script_dir/run_stage5_correction_window_tests.sh"
if [[ "$stage5_mode" == fresh ]]; then
  run_test stage5 stage5.fresh-runtime \
    "$script_dir/run_stage5_integrated_tests.sh" \
    --executable "$stage5_executable" \
    --output-root "$stage5_output_root" \
    --mpiexec "$resolved_mpi_launcher" \
    --chemistry-validation-events "$chemistry_validation_events"
  stage5_on_root=$stage5_output_root
else
  stage5_on_root=$stage5_archive
fi
run_test stage5 stage5.matrix check_stage5_matrix

run_test stage6 stage6.input-and-planner \
  "$script_dir/run_stage6_daughter_packet_tests.sh"
if [[ "$stage6_mode" == archived ]]; then
  run_test stage6 stage6.provenance check_stage6_provenance
fi
if [[ "$stage6_mode" == fresh ]]; then
  run_test stage6 stage6.fresh-runtime \
    "$script_dir/run_stage6_integrated_tests.sh" \
    --executable "$stage6_executable" \
    --output-root "$stage6_output_root" \
    --mpiexec "$resolved_mpi_launcher" \
    --chemistry-validation-events "$chemistry_validation_events"
  stage6_on_root=$stage6_output_root
else
  stage6_on_root=$stage6_archive
fi
run_test stage6 stage6.matrix check_stage6_matrix

run_test stage7 stage7.input-and-planner \
  "$script_dir/run_stage7_channel_importance_tests.sh"
if [[ "$stage7_mode" == archived ]]; then
  run_test stage7 stage7.provenance check_stage7_provenance
fi
if [[ "$stage7_mode" == fresh ]]; then
  run_test stage7 stage7.fresh-runtime \
    "$script_dir/run_stage7_integrated_tests.sh" \
    --executable "$stage7_executable" \
    --output-root "$stage7_output_root" \
    --mpiexec "$resolved_mpi_launcher" \
    --chemistry-validation-events "$chemistry_validation_events"
  stage7_on_root=$stage7_output_root
else
  stage7_on_root=$stage7_archive
fi
run_test stage7 stage7.matrix check_stage7_matrix

run_test stage8 stage8.input-and-planner \
  "$script_dir/run_stage8_total_reaction_forcing_tests.sh"
run_test stage8 stage8.analytic-reference env PYTHONDONTWRITEBYTECODE=1 \
  python3 "$script_dir/test_stage8_forced_segment_reference.py"
run_test stage8 stage8.parser-unit env PYTHONDONTWRITEBYTECODE=1 \
  python3 "$script_dir/test_stage8_total_reaction_forcing_parser.py"
if [[ "$stage8_mode" == archived ]]; then
  run_test stage8 stage8.provenance check_stage8_provenance
fi
if [[ "$stage8_mode" == fresh ]]; then
  run_test stage8 stage8.fresh-runtime \
    "$script_dir/run_stage8_integrated_tests.sh" \
    --executable "$stage8_executable" \
    --output-root "$stage8_output_root" \
    --mpiexec "$resolved_mpi_launcher" \
    --validation-particles "$validation_particles" \
    --chemistry-validation-events "$chemistry_validation_events"
  stage8_on_root=$stage8_output_root
else
  stage8_on_root=$stage8_archive
fi
run_test stage8 stage8.matrix check_stage8_matrix

stage_order=(stage0 stage1 stage2 stage3 stage4 stage5 stage6 stage7 stage8)
group_results=()
passed_groups=0
skipped_groups=0
failed_groups=0
executed_groups=0

# A test group is one roadmap stage.  A group fails when any contained test
# fails, is skipped only when it has no failures and at least one skip, and
# passes only when every contained test passes.  This gives the requested
# stage-level result without hiding the finer per-test diagnostics.
for summary_stage in "${stage_order[@]}"; do
  stage_passed=0
  stage_skipped=0
  stage_failed=0
  stage_total=0
  for ((summary_index=0; summary_index<${#test_names[@]}; summary_index++)); do
    if [[ "${test_stages[summary_index]}" != "$summary_stage" ]]; then
      continue
    fi
    stage_total=$((stage_total+1))
    case "${test_results[summary_index]}" in
      PASS) stage_passed=$((stage_passed+1)) ;;
      SKIP) stage_skipped=$((stage_skipped+1)) ;;
      FAIL) stage_failed=$((stage_failed+1)) ;;
    esac
  done

  if (( stage_total == 0 )); then
    continue
  fi
  executed_groups=$((executed_groups+1))
  if (( stage_failed != 0 )); then
    stage_result=FAIL
    failed_groups=$((failed_groups+1))
  elif (( stage_skipped != 0 )); then
    stage_result=SKIP
    skipped_groups=$((skipped_groups+1))
  else
    stage_result=PASS
    passed_groups=$((passed_groups+1))
  fi
  group_results+=(
    "$summary_stage result=$stage_result tests-run=$stage_total passed=$stage_passed skipped=$stage_skipped failed=$stage_failed"
  )
done

if (( failed_tests == 0 && skipped_tests == 0 )); then
  overall_result=PASS
else
  overall_result=FAIL
fi

summary_file="$log_root/summary.txt"
{
  printf '\nHartley2 all-stage test summary status=%s stage4-mode=%s stage5-mode=%s stage6-mode=%s stage7-mode=%s stage8-mode=%s\n' \
    "$overall_result" "$stage4_mode" "$stage5_mode" "$stage6_mode" \
    "$stage7_mode" "$stage8_mode"
  printf '  tests: run=%d passed=%d skipped=%d failed=%d\n' \
    "$executed_tests" "$passed_tests" "$skipped_tests" "$failed_tests"
  printf '  test-groups: run=%d passed=%d skipped=%d failed=%d\n' \
    "$executed_groups" "$passed_groups" "$skipped_groups" "$failed_groups"
  printf '%s\n' '  group-results:'
  for group_record in "${group_results[@]}"; do
    printf '    %s\n' "$group_record"
  done

  printf '%s\n' '  failed-tests:'
  if (( failed_tests == 0 )); then
    printf '%s\n' '    none'
  else
    for ((summary_index=0; summary_index<${#test_names[@]}; summary_index++)); do
      if [[ "${test_results[summary_index]}" == FAIL ]]; then
        printf '    test=%s stage=%s exit=%s log=%s\n' \
          "${test_names[summary_index]}" "${test_stages[summary_index]}" \
          "${test_exit_codes[summary_index]}" "${test_logs[summary_index]}"
      fi
    done
  fi

  printf '%s\n' '  skipped-tests:'
  if (( skipped_tests == 0 )); then
    printf '%s\n' '    none'
  else
    for ((summary_index=0; summary_index<${#test_names[@]}; summary_index++)); do
      if [[ "${test_results[summary_index]}" == SKIP ]]; then
        printf '    test=%s stage=%s log=%s\n' \
          "${test_names[summary_index]}" "${test_stages[summary_index]}" \
          "${test_logs[summary_index]}"
      fi
    done
  fi

  if [[ "$explicit_log_root" == true || "$overall_result" != PASS ]]; then
    printf '  logs-root: %s\n' "$log_root"
    printf '  summary-log: %s\n' "$summary_file"
  else
    printf '%s\n' '  logs: temporary; removed after successful completion'
  fi
} | tee "$summary_file"

if [[ "$overall_result" == PASS ]]; then
  printf 'PASS Hartley2 all-stage suite tests=%d test-groups=%d stage4-mode=%s stage5-mode=%s stage6-mode=%s stage7-mode=%s stage8-mode=%s\n' \
    "$executed_tests" "$executed_groups" "$stage4_mode" "$stage5_mode" \
    "$stage6_mode" "$stage7_mode" "$stage8_mode"
  suite_success=true
  exit 0
fi

printf 'FAIL Hartley2 all-stage suite tests=%d test-groups=%d stage4-mode=%s stage5-mode=%s stage6-mode=%s stage7-mode=%s stage8-mode=%s\n' \
  "$executed_tests" "$executed_groups" "$stage4_mode" "$stage5_mode" \
  "$stage6_mode" "$stage7_mode" "$stage8_mode" >&2
exit 1
