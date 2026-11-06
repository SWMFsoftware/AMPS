#!/usr/bin/env python3
"""Deterministic mutation suite for Stage-10 merge-readiness collection."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
TOOL = HERE / "stage10_merge_readiness.py"
POLICY = HERE / "stage10_release_policy.json"


def load_tool():
    spec = importlib.util.spec_from_file_location("stage10_readiness", TOOL)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot import readiness tool")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


READINESS = load_tool()


class TestFailure(RuntimeError):
    pass


def require(condition, message):
    if not condition:
        raise TestFailure(message)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(command, cwd, environment=None):
    return subprocess.run(
        list(command), cwd=cwd, env=environment, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
    )


def git(real_git, repository, *arguments):
    completed = run([real_git, *arguments], repository)
    if completed.returncode:
        raise TestFailure(f"fixture git {' '.join(arguments)} failed: {completed.stderr}")
    return completed.stdout.strip()


def make_git_shim(directory, real_git, transcript):
    directory.mkdir()
    shim = directory / "git"
    shim.write_text(
        "#!/usr/bin/env python3\nimport json,subprocess,sys\n"
        f"real={real_git!r}\nlog={str(transcript)!r}\n"
        "with open(log,'a',encoding='utf-8') as stream:\n"
        " stream.write(json.dumps(sys.argv[1:])+'\\n')\n"
        "raise SystemExit(subprocess.run([real,*sys.argv[1:]]).returncode)\n",
        encoding="utf-8",
    )
    shim.chmod(0o755)


def counts(number):
    return {"run": number, "passed": number, "skipped": 0, "failed": 0}


def write_receipt(root, repository, dev, head_tree):
    """Write a valid 39-test/11-group receipt and referenced files."""
    root.mkdir()
    per_stage = [1, 1, 2, 7, 6, 2, 3, 3, 5, 8, 1]
    tests, groups = [], []
    for stage_number, test_count in enumerate(per_stage):
        stage = f"stage{stage_number}"
        groups.append({"id": stage, "result": "PASS", "counts": counts(test_count)})
        for index in range(test_count):
            test_id = f"{stage}.fixture-{index}"
            log = root / f"{stage}.{index}.log"
            log.write_text(f"PASS {test_id}\n", encoding="utf-8")
            tests.append(
                {"id": test_id, "stage": stage, "result": "PASS", "exit_code": 0,
                 "log": log.name, "log_sha256": sha256(log)}
            )
    summary = root / "summary.txt"
    summary.write_text("PASS deterministic all-stage fixture\n", encoding="utf-8")
    receipt = {
        "schema": "amps.comet.all-stage-test-result.v1",
        "generated_at_utc": "2026-11-05T12:00:00Z",
        "repository": {"root": str(repository.resolve()), "head_commit": dev,
                       "head_tree": head_tree, "branch": "1d_comet_dev",
                       "detached_head": False},
        "modes": {f"stage{stage}": "archived" for stage in range(4, 10)},
        "tests": tests, "groups": groups,
        "counts": {"tests": counts(39), "groups": counts(11)},
        "summary": {"path": summary.name, "sha256": sha256(summary)},
        "interpretation": {
            "focused_implementation": {"result": "PASS", "accepted": True,
                                       "counts": counts(39), "meaning": "not science"},
            "stage9_science": {"result": "NOT_EVALUATED", "accepted": False,
                               "manifest": None, "manifest_sha256": None,
                               "meaning": "science remains distinct"},
        },
        "overall_result": "PASS", "overall_scope": "focused_implementation_only",
    }
    result = root / "all_stage_result.json"
    result.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def collect_command(repository, target, base, dev, receipt, output, *,
                    target_ref="refs/heads/master", makefile_digest=None, science=None):
    command = [
        sys.executable, str(TOOL), "collect", "--repo", str(repository),
        "--target-ref", target_ref, "--target-commit", target,
        "--base-commit", base, "--dev-commit", dev,
        "--all-stage-result", str(receipt),
        "--all-stage-result-sha256", sha256(receipt),
        "--output-dir", str(output), "--policy", str(POLICY),
    ]
    if makefile_digest:
        command.extend(["--dirty-exclusion", f"Makefile={makefile_digest}"])
    if science:
        command.extend(
            ["--stage9-science-manifest", str(science["manifest"]),
             "--stage9-science-manifest-sha256", sha256(science["manifest"]),
             "--stage9-science-result", str(science["result"]),
             "--stage9-science-result-sha256", sha256(science["result"]),
             "--rubin-comparison", str(science["rubin"]),
             "--rubin-comparison-sha256", sha256(science["rubin"])]
        )
    return command


def expect_failure(command, cwd, environment, label):
    completed = run(command, cwd, environment)
    require(completed.returncode != 0, f"{label} unexpectedly succeeded")


def forged_science(root):
    """Create PASS-shaped evidence that the strict evaluator must reject."""
    root.mkdir()
    seeds = [104729, 1299709, 15485863, 32452843]
    runs = [{"method": method, "seed": seed, "summary": "absent.json"}
            for seed in seeds for method in ("ANALOG", "FORCED")]
    manifest = root / "RUN_MANIFEST.json"
    manifest.write_text(json.dumps(
        {"schema": "amps.comet.stage9-science-manifest.v1", "policy": {}, "runs": runs},
        sort_keys=True) + "\n", encoding="utf-8")
    result = root / "SCIENCE_RESULT.json"
    result.write_text(json.dumps(
        {"schema": "amps.comet.stage9-science-result.v1", "outcome": "PASS",
         "accepted": True, "seeds": seeds, "run_count": 8, "failures": []},
        sort_keys=True) + "\n", encoding="utf-8")
    rubin = root / "comparison.json"
    rubin.write_text(json.dumps(
        {"schema": "amps.comet.stage9-rubin-comparison.v1", "stage9_outcome": "PASS",
         "seeds": seeds, "run_count": 8,
         "convergence": {"all_runs_converged_before_hard_limit": True},
         "provenance": {"manifest_sha256": sha256(manifest)}, "profiles": {"C2": {}}},
        sort_keys=True) + "\n", encoding="utf-8")
    for suffix in (".csv", ".png", ".eps", ".pdf"):
        rubin.with_suffix(suffix).write_bytes(b"forged\n")
    return {"manifest": manifest, "result": result, "rubin": rubin}


def main():
    real_git = shutil.which("git")
    if real_git is None:
        raise TestFailure("git is not available")
    with tempfile.TemporaryDirectory(prefix="stage10-readiness-selftest-") as temporary:
        root = Path(temporary)
        repository = root / "repository"
        repository.mkdir()
        # Support the older Git installed on some AMPS validation hosts, where
        # ``git init -b`` is unavailable and the default unborn branch is
        # already named ``master``.
        git(real_git, repository, "init")
        git(real_git, repository, "config", "user.name", "Stage10 Selftest")
        git(real_git, repository, "config", "user.email", "stage10@example.invalid")
        (repository / "Makefile").write_text("all:\n\t@true\n", encoding="utf-8")
        (repository / "srcHartley2").mkdir()
        source = repository / "srcHartley2" / "Comet.cpp"
        source.write_text("int comet = 0;\n", encoding="utf-8")
        git(real_git, repository, "add", "Makefile", "srcHartley2/Comet.cpp")
        git(real_git, repository, "commit", "-m", "baseline master")
        target = git(real_git, repository, "rev-parse", "HEAD")
        git(real_git, repository, "checkout", "-b", "1d_comet_dev")
        source.write_text("int comet = 1;\n", encoding="utf-8")
        git(real_git, repository, "add", "srcHartley2/Comet.cpp")
        git(real_git, repository, "commit", "-m", "develop comet")
        dev = git(real_git, repository, "rev-parse", "HEAD")
        tree = git(real_git, repository, "rev-parse", "HEAD^{tree}")
        base = git(real_git, repository, "merge-base", target, dev)
        receipt = write_receipt(root / "receipt", repository, dev, tree)

        transcript = root / "git-commands.jsonl"
        shim_directory = root / "shim"
        make_git_shim(shim_directory, real_git, transcript)
        environment = dict(os.environ)
        environment["PATH"] = f"{shim_directory}{os.pathsep}{environment.get('PATH', '')}"

        expect_failure(
            collect_command(repository, dev, base, dev, receipt, root / "wrong-target"),
            repository, environment, "wrong target")
        expect_failure(
            collect_command(repository, target, base, dev, receipt, root / "wrong-target-ref",
                            target_ref="refs/remotes/origin/master"),
            repository, environment, "wrong target ref")
        git(real_git, repository, "checkout", "-b", "wrong_branch")
        expect_failure(
            collect_command(repository, target, base, dev, receipt, root / "wrong-branch"),
            repository, environment, "wrong branch")
        git(real_git, repository, "checkout", "1d_comet_dev")

        missing_command = collect_command(
            repository, target, base, dev, receipt, root / "missing-tests")
        receipt_index = missing_command.index("--all-stage-result") + 1
        missing_command[receipt_index] = str(root / "absent.json")
        expect_failure(missing_command, repository, environment, "missing test receipt")

        failed_receipt = receipt.parent / "failed_result.json"
        failed = json.loads(receipt.read_text(encoding="utf-8"))
        failed["overall_result"] = "FAIL"
        failed_receipt.write_text(json.dumps(failed, sort_keys=True) + "\n", encoding="utf-8")
        expect_failure(
            collect_command(repository, target, base, dev, failed_receipt, root / "failed-tests"),
            repository, environment, "failed tests")
        wrong_tree_receipt = receipt.parent / "wrong_tree_result.json"
        wrong_tree = json.loads(receipt.read_text(encoding="utf-8"))
        wrong_tree["repository"]["head_tree"] = "0" * 40
        wrong_tree_receipt.write_text(
            json.dumps(wrong_tree, sort_keys=True) + "\n", encoding="utf-8")
        expect_failure(
            collect_command(repository, target, base, dev, wrong_tree_receipt,
                            root / "wrong-tree"),
            repository, environment, "wrong receipt tree")

        untracked_source = repository / "srcHartley2" / "UntrackedModel.cpp"
        untracked_source.write_text("int omitted_from_commit = 1;\n", encoding="utf-8")
        expect_failure(
            collect_command(repository, target, base, dev, receipt,
                            root / "untracked-source"),
            repository, environment, "untracked reviewable source")
        untracked_source.unlink()

        # The PT exclusions are delimiter-qualified (`PT/` and `PT.`), not a
        # broad `PT` prefix that could hide an unrelated root source file.
        prefix_lookalike = repository / "PTmalicious.py"
        prefix_lookalike.write_text("raise RuntimeError('not a PT artifact')\n", encoding="utf-8")
        expect_failure(
            collect_command(repository, target, base, dev, receipt,
                            root / "artifact-prefix-lookalike"),
            repository, environment, "generated-artifact prefix lookalike")
        prefix_lookalike.unlink()

        # A policy-listed raw artifact is excluded and summarized, not hashed
        # file-by-file.  It remains present through the successful collection
        # so the positive report path is exercised as well.
        generated = repository / "build" / "runtime.cache"
        generated.parent.mkdir()
        generated.write_bytes(b"generated runtime artifact\n")

        (repository / "Makefile").write_text(
            "all:\n\t@true\n# machine-local linker change\n", encoding="utf-8")
        makefile_digest = sha256(repository / "Makefile")
        expect_failure(
            collect_command(repository, target, base, dev, receipt, root / "unhashed-dirt"),
            repository, environment, "unhashed Makefile")
        source.write_text("int comet = 2;\n", encoding="utf-8")
        expect_failure(
            collect_command(repository, target, base, dev, receipt, root / "extra-dirt",
                            makefile_digest=makefile_digest),
            repository, environment, "extra tracked dirt")
        git(real_git, repository, "checkout", "--", "srcHartley2/Comet.cpp")

        bundle = root / "focused-bundle"
        command = collect_command(repository, target, base, dev, receipt, bundle,
                                  makefile_digest=makefile_digest)
        completed = run(command, repository, environment)
        require(completed.returncode == 0, f"valid collection failed: {completed.stderr}")
        report = json.loads((bundle / "readiness.json").read_text(encoding="utf-8"))
        require(report["status"] == "FOCUSED_ONLY_NOT_SCIENCE_READY",
                "focused-only status was not preserved")
        require(report["authorization"]["state"] == "NOT_GRANTED",
                "collector granted merge authorization")
        exclusions = report["worktree"]["tracked_dirty_exclusions"]
        require(len(exclusions) == 1 and exclusions[0]["worktree_sha256"] == makefile_digest,
                "Makefile exclusion is not hash-bound")
        artifact_summary = report["worktree"]["untracked_artifact_exclusions"]
        require(
            artifact_summary["path_count"] == 1
            and artifact_summary["groups"][0]["rule"] == "prefix:build/",
            "allowed untracked artifact was not summarized",
        )
        verify = [sys.executable, str(TOOL), "verify", "--bundle", str(bundle)]
        completed = run(verify, repository, environment)
        require(completed.returncode == 0, f"bundle verify failed: {completed.stderr}")
        expect_failure(command, repository, environment, "output overwrite")
        markdown = bundle / "readiness.md"
        markdown.write_text(markdown.read_text(encoding="utf-8") + "tamper\n", encoding="utf-8")
        expect_failure(verify, repository, environment, "bundle tamper")

        forged = forged_science(root / "forged-science")
        expect_failure(
            collect_command(repository, target, base, dev, receipt, root / "forged-ready",
                            makefile_digest=makefile_digest, science=forged),
            repository, environment, "forged Stage-9 PASS")

        # Commit an unclassified path last.  All provenance is otherwise
        # regenerated and valid, so failure specifically exercises semantic
        # path classification rather than stale-HEAD rejection.
        unknown = repository / "mystery.bin"
        unknown.write_bytes(b"unknown merge semantics\n")
        git(real_git, repository, "add", "mystery.bin")
        git(real_git, repository, "commit", "-m", "unknown path mutation")
        unknown_dev = git(real_git, repository, "rev-parse", "HEAD")
        unknown_tree = git(real_git, repository, "rev-parse", "HEAD^{tree}")
        unknown_receipt = write_receipt(
            root / "unknown-receipt", repository, unknown_dev, unknown_tree)
        expect_failure(
            collect_command(repository, target, base, unknown_dev, unknown_receipt,
                            root / "unknown-path", makefile_digest=makefile_digest),
            repository, environment, "unclassified diff path")

        records = [json.loads(line) for line in transcript.read_text(encoding="utf-8").splitlines()]
        require(records, "git command transcript is empty")
        for arguments in records:
            require(arguments[0] not in READINESS.ReadOnlyGit.FORBIDDEN_SUBCOMMANDS,
                    f"forbidden production git verb: {arguments}")
            READINESS.ReadOnlyGit.validate_arguments(arguments)

    print(
        "PASS Stage-10 merge-readiness selftest cases=target-ref,branch,receipt,"
        "untracked-source,artifact-exclusion,dirty-exclusion,extra-dirt,overwrite,"
        "tamper,forged-science,unknown-path,read-only-git"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except TestFailure as exc:
        print(f"stage10_merge_readiness_selftest.py: {exc}", file=sys.stderr)
        raise SystemExit(1)
