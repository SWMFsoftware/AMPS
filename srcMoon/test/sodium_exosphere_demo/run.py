#!/usr/bin/env python3
"""Run the bounded Tenishev-2013 Case-1 sodium demonstration.

This runner executes the normal srcMoon production executable and preserves
its density/column/brightness products.  It is intentionally outside the
frozen U01--U26/I01--I46 registry: a successful process is an EXPERIMENTAL
meeting demonstration, not PASS for observational validation or convergence.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys


TEST_ROOT = Path(__file__).resolve().parents[1]
REPO = TEST_ROOT.parents[1]
sys.path.insert(0, str(TEST_ROOT / "common"))
from moon_testlib import run_one  # noqa: E402


NOTICE = (
    "MEETING DEMONSTRATION - EXPERIMENTAL; NO OBSERVATIONAL VALIDATION OR "
    "NUMERICAL/STATISTICAL CONVERGENCE CLAIM"
)
PAPER = Path("/home/vtenishe/lunar-exosphere.pdf")
PAPER_SHA256 = "75d1e76039268c6f013ec80c53dfd7f0fb70f4779daff92461bd15185ff12b49"
SOURCE_ROADMAP = REPO / "srcMoon/development_and_validation.md"
DATA_ROADMAP = Path("/data/vtenishe/moon_validation_data/development_and_validation.md")


def sha256(path: Path) -> str:
    """Return the byte-level identity used in the immutable run manifest."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def replace_setting(text: str, key: str, value: str) -> str:
    """Replace one unambiguous flat Moon setting in a derived receipt."""
    pattern = rf"(?m)^({re.escape(key)}\s*=\s*)([^!\n]+?)(\s*)$"
    matches = list(re.finditer(pattern, text))
    if len(matches) != 1:
        raise ValueError(f"expected exactly one {key} in the input template")
    return re.sub(pattern, lambda match: match.group(1) + value + match.group(3),
                  text)


def parse_observation_utc(value: str) -> str:
    """Require an explicit UTC on the paper's stated date.

    The paper provides no time of day.  Requiring this argument prevents the
    template's midnight placeholder from silently becoming observation
    geometry.  Python validates civil syntax here; CSPICE remains authoritative
    for UTC-to-ET conversion inside production.
    """
    normalized = value[:-1] + "+00:00" if value.endswith("Z") else value
    try:
        parsed = dt.datetime.fromisoformat(normalized)
    except ValueError as exception:
        raise ValueError("--simulation-start-utc must be ISO-8601 UTC") from exception
    if parsed.tzinfo is None or parsed.utcoffset() != dt.timedelta(0):
        raise ValueError("--simulation-start-utc must include UTC designator Z or +00:00")
    if parsed.date() != dt.date(1991, 9, 30):
        raise ValueError(
            "the accepted near-Moon product is restricted to 1991-09-30; "
            "another date is a different campaign")
    return parsed.astimezone(dt.timezone.utc).isoformat(
        timespec="microseconds").replace("+00:00", "Z")


def git_text(*arguments: str) -> str:
    """Read repository identity without changing the checkout."""
    completed = subprocess.run(["git", *arguments], cwd=REPO, text=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, check=False)
    return completed.stdout.strip() if completed.returncode == 0 else "UNAVAILABLE"


def clean_build() -> None:
    """Apply the project-approved clean SPICE configuration/build sequence."""
    build = REPO / "build"
    if build.exists():
        shutil.rmtree(build)
    commands = [
        [str(REPO / "Config.pl"), "-application=moon",
         "-spice-path=/home/vtenishe/SPICE/cspice",
         "-spice-kernels=/home/vtenishe/SPICE/Kernels"],
        ["make", "-j"],
    ]
    for command in commands:
        completed = subprocess.run(command, cwd=REPO, check=False)
        if completed.returncode != 0:
            raise RuntimeError(
                f"clean build command exited {completed.returncode}: {command}")


def preflight_build(executable: Path, build_dir: Path) -> None:
    """Reject missing/stale or physically incompatible build inputs."""
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise RuntimeError(f"production executable is unavailable: {executable}")
    if not (build_dir / "libAMPS.a").is_file():
        raise RuntimeError(f"generated production library is unavailable: {build_dir}")
    moon_build = (REPO / "input/moon.input").read_text(encoding="utf-8")
    required_fragments = (
        "Source:ImpactVaporization=on",
        "Source:ThermalDesorption=on",
        "Source:PhotonStimulatedDesorption=on",
        "Source:Sputtering=off",
        "SaveParticleRestartFile=off",
    )
    missing = [item for item in required_fragments if item not in moon_build]
    if missing:
        raise RuntimeError(
            "input/moon.input does not satisfy the Case-1 build contract: " +
            ", ".join(missing))


def discover_product_groups(runtime: Path) -> dict[str, list[Path]]:
    """Return the production Na products required by this demonstration.

    The bounded Case-1 run advances only about 19 hours at the current global
    time step.  It therefore cannot trigger the separate
    ``Moon::Sampling::SubsolarLimbColumnIntegrals::PrintDataFile`` tables,
    whose production trigger is completion of a lunar orbital-direction
    cycle.  The generic exosphere observation operator runs at every 500-step
    output and writes the same required altitude coordinate, Na column
    density, and 5891.58/5897.56-A brightness variables into
    ``pic.LimbColumnDensity.*``.  Requiring those per-stage files exercises the
    normal bounded-run production path without pretending that a month-scale
    orbit-accumulated product was generated.
    """
    model_files = sorted(path for path in runtime.rglob("*") if path.is_file())
    return {
        "volume_number_density": [
            path for path in model_files
            if path.name.startswith("pic.NA.s=") and path.suffix == ".dat"
        ],
        "limb_altitude_column_and_d_line_brightness": [
            path for path in model_files
            if path.name.startswith("pic.LimbColumnDensity.") and
            path.suffix == ".dat"
        ],
        "spherical_column_and_d_line_brightness_map": [
            path for path in model_files
            if path.name.startswith("pic.ColumnDensityMap.spherical.") and
            path.suffix == ".dat"
        ],
        "antisunward_column_and_d_line_brightness_map": [
            path for path in model_files
            if path.name.startswith("pic.Moon.Anti-sunwardColumnIntegrals.out=")
            and path.suffix == ".dat"
        ],
    }


def write_success_receipts(arguments: argparse.Namespace, output: Path,
                           observation_utc: str, command: list[str],
                           verification: dict[str, str],
                           recovery_mode: bool = False) -> dict[str, object]:
    """Validate production artifacts and write immutable success receipts.

    ``recovery_mode`` exists only to re-run deterministic post-processing when
    an already completed expensive MPI run was rejected by an older runner
    contract.  It never launches AMPS and never changes a model artifact.
    """
    runtime = output / "runtime"
    application = output / "application.in"
    stdout = output / "stdout.log"
    stderr = output / "stderr.log"
    template = arguments.input.resolve()
    executable = arguments.executable.resolve()
    build_dir = arguments.build_dir.resolve()

    product_groups = discover_product_groups(runtime)
    missing_products = [
        name for name, paths in product_groups.items()
        if not paths or not any(path.stat().st_size > 0 for path in paths)
    ]
    if missing_products:
        raise RuntimeError(
            "production completed without nonempty required Na products: " +
            ", ".join(missing_products))
    scientific_files = sorted({
        path for paths in product_groups.values() for path in paths
        if path.stat().st_size > 0
    })
    model_files = sorted(path for path in runtime.rglob("*") if path.is_file())

    formal_vv = {
        "I11_numerical_convergence": "NOT RUN",
        "I13_baseline_acceptance": "NOT RUN",
        "I32_observational_validation": "SKIPPED - Flynn/Mendillo 1991 "
            "machine-readable observations, uncertainty, and exact UTC are unavailable",
    }
    manifest = {
        "schema": "srcMoon-tenishev-2013-case1-demo-manifest-v1",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "notice": NOTICE,
        "postprocessing_mode": (
            "existing-completed-run-recovery" if recovery_mode else "normal"),
        "source_revision": git_text("rev-parse", "HEAD"),
        "source_dirty_paths": git_text("status", "--short").splitlines(),
        "command": command, "mpi_ranks": arguments.mpi_ranks,
        "openmp_threads": 1, "simulation_start_utc": observation_utc,
        "paper": {"path": str(PAPER), "sha256": sha256(PAPER)},
        "input": {"path": str(application), "sha256": sha256(application)},
        "template": {"path": str(template), "sha256": sha256(template)},
        "executable": {"path": str(executable), "sha256": sha256(executable)},
        "build_library": {"path": str(build_dir / "libAMPS.a"),
                          "sha256": sha256(build_dir / "libAMPS.a")},
        "roadmaps": {
            "source": {"path": str(SOURCE_ROADMAP),
                       "sha256": sha256(SOURCE_ROADMAP)},
            "validation_data": {"path": str(DATA_ROADMAP),
                                "sha256": sha256(DATA_ROADMAP)},
        },
        "local_kernel_gates": verification,
        "formal_vv": formal_vv,
        "required_product_groups": {
            name: [str(path) for path in paths]
            for name, paths in product_groups.items()
        },
        "model_artifacts": [
            {"path": str(path), "sha256": sha256(path),
             "bytes": path.stat().st_size} for path in model_files],
    }
    manifest_path = output / "run_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    result = {
        "schema": "srcMoon-tenishev-2013-case1-demo-result-v1",
        "execution_status": "PASS", "capability_label": "EXPERIMENTAL",
        "notice": NOTICE, "simulation_start_utc": observation_utc,
        "iterations_requested": 5000, "sample_length_steps": 500,
        "postprocessing_mode": manifest["postprocessing_mode"],
        "scientific_artifact_count": len(scientific_files),
        "scientific_artifacts": [str(path) for path in scientific_files],
        "required_product_groups": {
            name: [str(path) for path in paths]
            for name, paths in product_groups.items()
        },
        "local_kernel_gates": verification, "formal_vv": formal_vv,
        "artifacts": {
            "application": {"path": str(application),
                            "sha256": sha256(application)},
            "stdout": {"path": str(stdout), "sha256": sha256(stdout)},
            "stderr": {"path": str(stderr), "sha256": sha256(stderr)},
            "run_manifest": {"path": str(manifest_path),
                             "sha256": sha256(manifest_path)},
        },
    }
    (output / "result.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def finalize_existing(arguments: argparse.Namespace) -> dict[str, object]:
    """Finalize a preserved run after a post-processing-contract failure.

    Recovery is fail-closed: it requires the normal AMPS completion trailer,
    the exact derived UTC receipt, all six local PASS artifacts, and all four
    nonempty production product families.  It does not infer observational or
    convergence status from process completion.
    """
    output = arguments.output_dir.resolve()
    application = output / "application.in"
    stdout = output / "stdout.log"
    runtime = output / "runtime"
    for required in (application, stdout, output / "stderr.log", runtime):
        if not required.exists():
            raise RuntimeError(f"cannot finalize; required artifact is absent: {required}")
    observation_utc = parse_observation_utc(arguments.simulation_start_utc)
    application_text = application.read_text(encoding="utf-8")
    if f"simulation_start_utc = {observation_utc}" not in application_text:
        raise RuntimeError(
            "cannot finalize; application.in does not contain the requested UTC")
    stdout_text = stdout.read_text(encoding="utf-8", errors="replace")
    if "Total AMPS runtime is" not in stdout_text:
        raise RuntimeError(
            "cannot finalize; stdout lacks the normal AMPS completion trailer")

    verification: dict[str, str] = {}
    for test_id in ("U03", "U04", "U05", "U12", "U14", "U22"):
        result_path = output / "verification" / test_id / "result.json"
        if not result_path.is_file():
            raise RuntimeError(f"cannot finalize; missing {result_path}")
        result = json.loads(result_path.read_text(encoding="utf-8"))
        verification[test_id] = str(result.get("status", "ERROR"))
        if verification[test_id] != "PASS":
            raise RuntimeError(
                f"cannot finalize; {test_id} is {verification[test_id]}")

    command = ["mpiexec", "-n", str(arguments.mpi_ranks),
               str(arguments.executable.resolve()), "-input", str(application)]
    return write_success_receipts(arguments, output, observation_utc, command,
                                  verification, recovery_mode=True)


def run(arguments: argparse.Namespace) -> dict[str, object]:
    """Execute local gates and one bounded eight-rank production run."""
    if arguments.rebuild:
        clean_build()
    output = arguments.output_dir.resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"refusing to overwrite non-empty output: {output}")
    output.mkdir(parents=True, exist_ok=True)
    if arguments.mpi_ranks <= 0:
        raise ValueError("--mpi-ranks must be positive")
    observation_utc = parse_observation_utc(arguments.simulation_start_utc)
    executable = arguments.executable.resolve()
    build_dir = arguments.build_dir.resolve()
    preflight_build(executable, build_dir)
    if not PAPER.is_file() or sha256(PAPER) != PAPER_SHA256:
        raise RuntimeError(
            "the reviewed /home/vtenishe/lunar-exosphere.pdf is missing or "
            "its SHA-256 differs")

    # These gates cover forces, radiation pressure/shadowing, impacts,
    # photoionization and Na source kernels.  They do not replace linked
    # convergence or observational tests.
    verification_root = output / "verification"
    verification: dict[str, str] = {}
    for test_id in ("U03", "U04", "U05", "U12", "U14", "U22"):
        result = run_one(test_id, build_dir, verification_root)
        verification[test_id] = str(result.get("status", "ERROR"))
        if verification[test_id] != "PASS":
            raise RuntimeError(
                f"Case-1 runner requires {test_id} PASS; observed "
                f"{verification[test_id]}: {result.get('reason', '')}")

    template = arguments.input.resolve()
    text = template.read_text(encoding="utf-8")
    text = replace_setting(text, "simulation_start_utc", observation_utc)
    text = replace_setting(
        text, "surface_cea_file", str((output / "sphere-unused.cea").resolve()))
    text = replace_setting(
        text, "surface_tecplot_file", str((output / "sphere-unused.dat").resolve()))
    application = output / "application.in"
    application.write_text(text, encoding="utf-8")

    runtime = output / "runtime"
    runtime.mkdir()
    stdout = output / "stdout.log"
    stderr = output / "stderr.log"
    command = ["mpiexec", "-n", str(arguments.mpi_ranks), str(executable),
               "-input", str(application)]
    environment = dict(os.environ)
    environment["OMP_NUM_THREADS"] = "1"
    with stdout.open("w", encoding="utf-8") as out, \
            stderr.open("w", encoding="utf-8") as err:
        completed = subprocess.run(command, cwd=runtime, env=environment,
                                   stdout=out, stderr=err, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            f"production executable exited {completed.returncode}; see "
            f"{stdout} and {stderr}")

    return write_success_receipts(arguments, output, observation_utc, command,
                                  verification)


def main() -> int:
    timestamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    parser = argparse.ArgumentParser(
        description=("Run the bounded production Tenishev et al. (2013) "
                     "Case-1 Na exosphere demonstration."),
        epilog=("The publication gives 1991-09-30 but no observation time. "
                "Supply an independently justified exact UTC; the runner "
                "does not interpret successful execution as validation."))
    parser.add_argument("--simulation-start-utc", required=True,
                        help="exact ISO-8601 UTC on 1991-09-30; no default")
    parser.add_argument("--mpi-ranks", type=int, default=8)
    parser.add_argument("--input", type=Path,
                        default=REPO / "srcMoon/examples/tenishev2013_case1.in")
    parser.add_argument("--executable", type=Path, default=REPO / "amps")
    parser.add_argument("--build-dir", type=Path, default=REPO / "build")
    parser.add_argument("--rebuild", action="store_true",
                        help="remove build, configure SPICE Moon, and make -j")
    parser.add_argument(
        "--finalize-existing", action="store_true",
        help=("do not launch AMPS; re-check and write receipts for an existing "
              "completed output after a runner post-processing failure"))
    parser.add_argument("--output-dir", type=Path,
                        default=REPO / "test_output/srcMoon/sodium-exosphere-demo" /
                        timestamp)
    arguments = parser.parse_args()
    try:
        if arguments.rebuild and arguments.finalize_existing:
            raise ValueError("--rebuild and --finalize-existing are incompatible")
        result = (finalize_existing(arguments) if arguments.finalize_existing
                  else run(arguments))
    except (OSError, ValueError, RuntimeError) as exception:
        print(f"ERROR: {exception}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
