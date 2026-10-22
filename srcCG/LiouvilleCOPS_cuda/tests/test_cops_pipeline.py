"""The campaign runner tools/cluster/cops_pipeline.py (Voyager; C15/P0), CPU only.

Every stage command is a fake tool (``FAKE_TOOL`` below, installed under the real tool names in a
temporary application tree) that logs its arguments and writes the products the runner checks. Covered:

* stage order, the exact command lines per fit configuration (templates, record store, devices), and a
  second ``run`` that does nothing;
* ``plan`` writes nothing; ``status --check`` exit codes (0 complete, 75 resumable, 1 failed);
* a failed stage is moved aside (never deleted) and runs again; a restart directory refuses another mesh
  or emission law; external base/store/bundles are identity-checked and skipped;
* wall time: a fit stopped at the deadline continues as segment k+1 from the interrupted restart with
  the remaining evaluations; a failed segment runs again from the last clean one; segments are recorded
  in the stamp and pipeline.json; a stage whose recorded duration exceeds the time left is not started;
* the memory model behind ``--precision auto`` (FP64 resident, else FP32 with streaming, host check);
* ``--gravity-precision`` is independent of ``--precision``; both policies are recorded; the gravity
  time estimate comes from recorded runs or is reported unknown;
* ``submit``: GPUs per node from ``sinfo``, task lines from --nodes, afterok chain or self-resubmission;
  ``sbatch`` is never run; fit_mpi.slurm hardwires no rank count;
* ``--restart-from`` template checks.

The GPU and MPI end-to-end runs of the runner are pending (after the 20M chain; STATUS.md C15).
"""

import json
import os
from pathlib import Path
import re
import stat
import sys

import numpy as np
import pytest

APP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP / "tools" / "cluster"))
import cops_pipeline as cp  # noqa: E402

FACES = 12
TOOLS = ("rebuild_cops_boundary.py", "c15_record_store.py", "c05_write_bundle.py", "c04_gravity_equivalence.py",
         "c04_write_gpu_cache.py", "gpu_fit.py", "fine_tune_gpu.py", "restart_tool.py")

FAKE_TOOL = r'''
import hashlib, json, os, struct, sys
from pathlib import Path
name = Path(__file__).name
argv = sys.argv[1:]
keep = ("COPS_CUDA_RECORD_STORE", "COPS_OUTPUT_ROOT", "COPS_CUDA_DEVICE", "LIOUVILLE_COPS_REFERENCE")
with open(os.environ["FAKE_CALLS"], "a") as log:
    log.write(json.dumps({"tool": name, "argv": argv, "env": {k: os.environ.get(k) for k in keep}}) + "\n")

def opt(flag, default=None):
    return argv[argv.index(flag) + 1] if flag in argv else default

def fail(key):
    if os.environ.get("FAKE_FAIL") == key:
        print("fake failure of " + key, flush=True)
        sys.exit(3)

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

if name == "rebuild_cops_boundary.py":
    out = Path(opt("--out")); out.mkdir(parents=True); (out / "kernels").mkdir()
    (out / "elements.bin").write_bytes(struct.pack("<iii", 1, int(os.environ.get("FAKE_FACES", "12")), 25))
    for i in range(3):
        (out / "kernels" / ("point=%04d.bin" % i)).write_bytes(b"\0" * (92 + 40 * 100))
    (out / "run_manifest.json").write_text(json.dumps({"states": {"process": "complete", "cache": "verified"}}))
    fail("base")
    (out / "rebuild.json").write_text(json.dumps({"status": "complete", "surface_bc": opt("--surface-bc"),
        "shape_model": opt("--shape-model"), "run_manifest": {"sha256": sha(out / "run_manifest.json")}}))
elif name == "c15_record_store.py":
    run, out = Path(opt("--run")), Path(opt("--out")); out.mkdir()
    fail("store")
    (out / "store_manifest.json").write_text(json.dumps({"base_run": str(run), "records": 300,
        "run_manifest": {"sha256": sha(run / "run_manifest.json")}, "mesh_identity": "m" * 64,
        "record_precision": "fp32", "nfaces": int(os.environ.get("FAKE_FACES", "12"))}))
elif name == "c05_write_bundle.py":
    out = Path(opt("--out"))
    names = ["late_ng_not_rg", "late_rg_not_ng"] + ([] if "--visibility-only" in argv else
                                                    ["cliff_ge_45deg", "cliff_ge_60deg"])
    fail("visibility" if "--visibility-only" in argv else "bundle")
    out.write_bytes(b"npz")
    out.with_suffix(".json").write_text(json.dumps({"base_run": opt("--run"), "mesh_identity": "m" * 64,
                                                    "names": names}))
elif name == "c04_gravity_equivalence.py":
    eq = Path(opt("--output-root")) / "cuda_c04" / "equivalence_20990101T000000Z"; eq.mkdir(parents=True)
    fail("gravity")
    (eq / "gpu_field.npz").write_bytes(b"field")
    (eq / "gpu_manifest.json").write_text("{}")
elif name == "c04_write_gpu_cache.py":
    fail("cache")
    Path(opt("--out")).write_bytes(b"cache")
elif name == "gpu_fit.py":
    split = argv.index("--")
    out = Path(argv[argv.index("--out", split) + 1]); out.mkdir()
    fit = os.environ.get("FAKE_FIT_NAME", out.parent.name.split("_")[0])
    fail("fit:" + out.parent.name.split("_")[0])
    budget, segment = int(opt("--max-nfev")), int(opt("--segment", "1"))
    step = int(os.environ.get("FAKE_FIT_STEP", "0"))
    if "--deadline" in argv and os.environ.get("FAKE_FIT_EXHAUSTED"):
        sys.exit(75)
    stopped = "--deadline" in argv and step and budget > step
    nfev = step if stopped else budget
    record = {"index": segment, "status": "interrupted" if stopped else "complete", "nfev": nfev,
              "max_nfev": budget, "cost": 1.0 / (segment + 1)}
    (out / "fit_restart.npz").write_bytes(b"restart %d" % segment)
    (out / "fit_result.txt").write_text("final_cost %.6e\n" % record["cost"])
    (out / "fit_manifest.json").write_text(json.dumps({"states": {"process": "complete"},
        "compute_backend": {"wall_time": {"segments": [record]}}}))
    sys.exit(75 if stopped else 0)
elif name == "fine_tune_gpu.py":
    out = Path(opt("--out")); out.mkdir()
    fail("finetune")
    (out / "fine_tuning.json").write_text(json.dumps({"restart": opt("--restart-from")}))
elif name == "restart_tool.py":
    Path(opt("--save")).write_bytes(b"translated")
'''


@pytest.fixture
def world(tmp_path, monkeypatch):
    """A fake application and reference tree, a mesh, and the fake tools' call log."""
    app, ref = tmp_path / "app", tmp_path / "LiouvilleCOPS"
    (app / "tools").mkdir(parents=True)
    (app / "cops_cuda").mkdir()
    (app / "cops_cuda" / "x.py").write_text("x = 1\n")
    (app / "native").mkdir()
    (app / "native" / "liouville_cops_cuda").write_text("#!/bin/sh\n")
    (ref / "tools").mkdir(parents=True)
    for name in TOOLS:
        target = (ref if name == "rebuild_cops_boundary.py" else app) / "tools" / name
        target.write_text(FAKE_TOOL)
    (ref / "tools" / "fit_cops.py").write_text("# reference fitter\n")
    (ref / "fit_base.in").write_text("ObservationFile = /nonexistent/obs.h\nKernelDirectory = /nonexistent\n"
                                     "Kernel = lsk/NAIF0011.TLS\nSurfaceBC = flat\n")
    shape = tmp_path / "mesh.bdf"
    shape.write_text("$ test mesh\n" + "".join("GRID %d\n" % k for k in range(8)) +
                     "".join("CTRIA3 %d 1 1 2 3\n" % k for k in range(FACES)))
    calls = tmp_path / "calls.jsonl"
    calls.write_text("")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    write_exe(bin_dir / "sbatch", "#!/bin/sh\necho called >> %s\nexit 1\n" % (tmp_path / "sbatch_called"))
    monkeypatch.setenv("PATH", "%s:%s" % (bin_dir, os.environ["PATH"]))
    monkeypatch.setenv("FAKE_CALLS", str(calls))
    for name in ("SLURM_JOB_ID", "SLURM_JOB_NUM_NODES", "SLURM_NNODES", "SLURM_GPUS_ON_NODE", "SLURM_NTASKS",
                 "SLURM_MEM_PER_NODE", "SLURM_CPUS_ON_NODE", "SLURM_JOB_END_TIME", "FAKE_FAIL", "FAKE_FIT_STEP",
                 "FAKE_FIT_EXHAUSTED"):
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    monkeypatch.delenv("COPS_OUTPUT_ROOT", raising=False)
    return {"tmp": tmp_path, "app": app, "ref": ref, "shape": shape, "calls": calls, "bin": bin_dir,
            "restart": tmp_path / "campaign"}


def write_exe(path, text):
    path.write_text(text)
    path.chmod(path.stat().st_mode | stat.S_IEXEC)


def base_args(world, *extra, restart=None, shape=None, bc="flat", gpus="2"):
    return (["--restart-dir", str(restart or world["restart"]), "--shape", str(shape or world["shape"]),
             "--surface-bc", bc, "--app-root", str(world["app"]), "--reference", str(world["ref"]),
             "--python", sys.executable] + (["--gpus-per-node", gpus] if gpus else []) +
            ["--gpu-memory-gb", "48", "--host-memory-gb", "256", "--nodes", "1", "--threads", "4"] + list(extra))


def pipeline(command, world, *extra, **kw):
    return cp.main([command] + base_args(world, *extra, **kw))


def calls(world):
    return [json.loads(line) for line in world["calls"].read_text().splitlines()]


def tool_order(world):
    out = []
    for c in calls(world):
        if c["tool"] == "gpu_fit.py":
            out.append("fit:" + Path(c["argv"][c["argv"].index("--out", c["argv"].index("--")) + 1]).parent.name.split("_")[0])
        elif c["tool"] == "c05_write_bundle.py":
            out.append("visibility" if "--visibility-only" in c["argv"] else "bundle")
        else:
            out.append(c["tool"].replace(".py", ""))
    return out


def fit_call(world, fit, index=-1):
    found = [c for c in calls(world) if c["tool"] == "gpu_fit.py" and
             Path(c["argv"][c["argv"].index("--out", c["argv"].index("--")) + 1]).parent.name.startswith(fit + "_")]
    return found[index]


def test_plan_writes_nothing_and_lists_the_stages_in_order(world, capsys):
    assert pipeline("plan", world, "--precision", "mixed", "--fine-tune") == 0
    out = capsys.readouterr().out
    assert not world["restart"].exists() and calls(world) == []
    order = [line.split()[0] for line in out.splitlines() if re.match(r"^\S+\s+(prepare|fit|finetune)\s", line)]
    assert order == ["base", "store", "visibility", "gravity", "cache", "bundle", "fit:B29_n2000", "fit:N30_n2000",
                     "fit:C30_n2000", "fit:NC31_n2000", "finetune:B29_n2000", "finetune:N30_n2000",
                     "finetune:C30_n2000", "finetune:NC31_n2000"]
    assert "MISSING       base input file: ObservationFile = /nonexistent/obs.h" in out
    assert "records mixed, fp32_records_fp64_arithmetic_v1" in out
    assert "gravity       mixed (mixed_near_fp64_exact_far_fp32_quad3_v1)" in out
    assert "%d facets" % FACES in out


def test_run_in_order_with_the_production_commands_then_nothing(world, capsys):
    assert pipeline("run", world, "--precision", "mixed", "--fine-tune", "--max-nfev", "50") == 0
    assert tool_order(world) == ["rebuild_cops_boundary", "c15_record_store", "visibility", "c04_gravity_equivalence",
                                 "c04_write_gpu_cache", "bundle", "fit:B29", "fit:N30", "fit:C30", "fit:NC31",
                                 "fine_tune_gpu", "fine_tune_gpu", "fine_tune_gpu", "fine_tune_gpu"]
    d = world["restart"]
    vis, bundle = str(d / "visibility" / "visibility_templates.npz"), str(d / "bundle" / "restart_templates.npz")
    expect = {"B29": [], "N30": ["--activity-template", vis, "late_ng_not_rg"],
              "C30": ["--activity-template", bundle, "cliff_ge_60deg"],
              "NC31": ["--activity-template", bundle, "late_ng_not_rg", "--activity-template", bundle, "cliff_ge_60deg"]}
    for fit, templates in expect.items():
        c = fit_call(world, fit)
        argv = c["argv"]
        split = argv.index("--")
        assert argv[:split] == ["--devices", "0,1", "--precision", "mixed", "--record-store", str(d / "store"),
                                "--max-nfev", "50", "--segment", "1"]
        out = argv.index("--out", split)
        assert argv[split + 1] == str(d / "base")
        assert argv[split + 2:out] == cp.DEFAULT_FIT_ARGS.split()
        assert argv[out + 1] == str(d / "fits" / ("%s_n50" % fit) / "seg01") and argv[out + 2:] == templates
        assert c["env"]["COPS_CUDA_RECORD_STORE"] == str(d / "store")
    base_call = calls(world)[0]
    assert base_call["env"]["COPS_CUDA_DEVICE"] == "0"
    assert base_call["argv"][base_call["argv"].index("--surface-bc") + 1] == "flat"
    gravity = next(c for c in calls(world) if c["tool"] == "c04_gravity_equivalence.py")["argv"]
    assert gravity[gravity.index("--devices") + 1] == "0,1"
    assert gravity[gravity.index("--precision"):gravity.index("--precision") + 4] == ["--precision", "mixed",
                                                                                      "--near-factor", "6"]
    ft = [c["argv"] for c in calls(world) if c["tool"] == "fine_tune_gpu.py"]
    assert ft[0][ft[0].index("--restart-from") + 1] == str(d / "fits" / "B29_n50" / "seg01" / "fit_restart.npz")
    assert ft[0][ft[0].index("--minutes") + 1] == "5" and ft[0][ft[0].index("--method") + 1] == "cpu"
    manifest = json.loads((d / "pipeline.json").read_text())
    assert manifest["policies"]["record_precision"]["policy"] == "fp32_records_fp64_arithmetic_v1"
    assert manifest["policies"]["gravity_precision"]["policy"] == "mixed_near_fp64_exact_far_fp32_quad3_v1"
    assert manifest["fits"]["NC31_n50"]["segments"][0]["nfev"] == 50
    n = len(calls(world))
    capsys.readouterr()
    assert pipeline("run", world) == 0 and len(calls(world)) == n        # stored request; nothing left to do
    assert cp.main(["status", "--restart-dir", str(d), "--check"]) == 0


def test_auto_precision_and_needed_stages(world):
    assert pipeline("run", world, "--fits", "B29") == 0          # 300 records, 12 facets: FP64 fits resident
    assert tool_order(world) == ["rebuild_cops_boundary", "fit:B29"]
    argv = fit_call(world, "B29")["argv"]
    assert argv[argv.index("--precision") + 1] == "fp64" and "--record-store" not in argv
    resolved = json.loads((world["restart"] / "pipeline.json").read_text())["resolved"]["precision"]
    assert resolved["value"] == "fp64" and resolved["estimate"]["records"] == 300
    assert resolved["hardware"]["gpu_memory_bytes"] == 48e9


def test_gravity_precision_is_its_own_setting(world, capsys):
    assert pipeline("run", world, "--fits", "C30", "--precision", "fp64", "--gravity-precision", "fp64") == 0
    gravity = next(c for c in calls(world) if c["tool"] == "c04_gravity_equivalence.py")["argv"]
    assert gravity[gravity.index("--precision") + 1] == "fp64" and "--near-factor" not in gravity
    policies = json.loads((world["restart"] / "pipeline.json").read_text())["policies"]
    assert policies["gravity_precision"] == {"value": "fp64", "policy": "fp64_reference_v1", "reason": "requested"}
    assert policies["record_precision"]["policy"] == "fp64_reference_v1"
    assert pipeline("plan", world, "--gravity-precision", "mixed") == cp.EXIT_REFUSED
    assert "gravity_precision" in capsys.readouterr().err


def test_policy_names_match_the_modules():
    precision = (APP / "cops_cuda" / "precision.py").read_text()
    gravity = (APP / "cops_cuda" / "gravity.py").read_text()
    assert 'POLICY_NAME = "%s"' % cp.RECORD_POLICIES["fp64"] in precision
    assert 'MIXED_POLICY_NAME = "%s"' % cp.RECORD_POLICIES["mixed"] in precision
    assert 'MIXED_POLICY = "%s"' % cp.GRAVITY_POLICIES["mixed"] in gravity
    assert cp.GRAVITY_POLICIES["fp64"] == cp.RECORD_POLICIES["fp64"]     # gravity fp64 records precision.POLICY_NAME


def test_gravity_time_from_recorded_runs_or_unknown(world, capsys, tmp_path):
    pipeline("plan", world, "--fits", "C30")
    out = capsys.readouterr().out
    assert "gravity time *mixed unknown (no recorded mixed gravity run)" in out
    root = tmp_path / "output_root"
    eq = root / "cuda_c04" / "equivalence_20261022T005551Z"
    eq.mkdir(parents=True)
    (eq / "gpu_manifest.json").write_text(json.dumps({
        "targets": 20, "seconds": 100.0, "pairs_per_second": 8.0, "devices": [0, 1],
        "evaluator": {"precision": {"policy": "mixed_near_fp64_exact_far_fp32_quad3_v1"}},
        "manifest": {"device": {"name": "NVIDIA RTX A6000"}}}))
    pipeline("plan", world, "--fits", "C30", "--output-root", str(root))
    out = capsys.readouterr().out
    # 12 x 12 pairs at 4 pairs/s per GPU on 2 GPUs = 18 s
    assert "gravity time *mixed ~18 s on 2 GPU(s), scaled from equivalence_20261022T005551Z" in out
    assert "NVIDIA RTX A6000" in out and "gravity time  fp64  unknown" in out


def test_failed_stage_is_moved_aside_and_runs_again(world, monkeypatch):
    monkeypatch.setenv("FAKE_FAIL", "gravity")
    assert pipeline("run", world, "--fits", "C30", "--precision", "mixed") == cp.EXIT_FAILED
    d = world["restart"]
    assert json.loads((d / "stamps" / "gravity.json").read_text())["status"] == "failed"
    assert cp.main(["status", "--restart-dir", str(d), "--check"]) == cp.EXIT_FAILED
    monkeypatch.delenv("FAKE_FAIL")
    assert pipeline("run", world) == 0
    moved = sorted(p.name for p in d.iterdir() if p.name.startswith("gravity.incomplete_"))
    assert len(moved) == 1 and (d / moved[0] / "cuda_c04").is_dir()
    assert any(p.name.startswith("gravity.json.incomplete_") for p in (d / "stamps").iterdir())
    assert tool_order(world).count("c04_gravity_equivalence") == 2 and tool_order(world).count("rebuild_cops_boundary") == 1


def test_other_mesh_or_emission_law_is_refused(world, capsys, tmp_path):
    assert pipeline("run", world, "--fits", "B29") == 0
    capsys.readouterr()
    assert pipeline("run", world, "--fits", "B29", bc="rough") == cp.EXIT_REFUSED
    assert "surface_bc differ" in capsys.readouterr().err
    other = tmp_path / "other.bdf"
    other.write_text(world["shape"].read_text() + "$ changed\n")
    assert pipeline("plan", world, shape=other) == cp.EXIT_REFUSED
    assert "shape differ" in capsys.readouterr().err


def test_external_products_are_checked_and_skipped(world, capsys, tmp_path):
    first = tmp_path / "first"
    assert pipeline("run", world, "--fits", "NC31,N30", "--precision", "mixed", restart=first) == 0
    n = len(calls(world))
    ext = ["--base", str(first / "base"), "--store", str(first / "store"),
           "--visibility-bundle", str(first / "visibility" / "visibility_templates.npz"),
           "--bundle", str(first / "bundle" / "restart_templates.npz")]
    capsys.readouterr()
    assert pipeline("plan", world, "--fits", "NC31,N30", *ext) == 0
    out = capsys.readouterr().out
    for sid in ("base", "store", "visibility", "bundle"):
        assert re.search(r"^%s\s+prepare\s+external" % sid, out, re.M), sid
    assert re.search(r"^gravity\s+prepare\s+not needed", out, re.M)
    assert pipeline("run", world, "--fits", "NC31,N30", *ext) == 0
    assert [c["tool"] for c in calls(world)[n:]] == ["gpu_fit.py", "gpu_fit.py"]
    argv = fit_call(world, "NC31")["argv"]
    assert argv[argv.index("--record-store") + 1] == str(first / "store")
    capsys.readouterr()
    assert pipeline("plan", world, *ext, restart=tmp_path / "rough", bc="rough") == cp.EXIT_REFUSED
    assert "SurfaceBC=flat, not rough" in capsys.readouterr().err


def test_wall_time_segments_resume_with_the_remaining_evaluations(world, monkeypatch, capsys):
    monkeypatch.setenv("FAKE_FIT_STEP", "7")
    args = ("--fits", "B29", "--max-nfev", "20", "--wall-time", "01:00:00", "--fine-tune")
    assert pipeline("run", world, *args) == cp.EXIT_RESUMABLE
    d = world["restart"]
    assert cp.main(["status", "--restart-dir", str(d), "--check"]) == cp.EXIT_RESUMABLE
    assert pipeline("run", world) == cp.EXIT_RESUMABLE
    assert pipeline("run", world) == 0
    fits = [c["argv"] for c in calls(world) if c["tool"] == "gpu_fit.py"]
    assert [a[a.index("--max-nfev") + 1] for a in fits] == ["20", "13", "6"]
    assert [a[a.index("--segment") + 1] for a in fits] == ["1", "2", "3"]
    assert "--restart-from" not in fits[0]
    for k in (1, 2):
        assert fits[k][fits[k].index("--restart-from") + 1] == str(d / "fits" / "B29_n20" / ("seg%02d" % k) /
                                                                     "fit_restart.npz")
    assert all("--deadline" in a and a[a.index("--save-seconds") + 1] == "600" for a in fits)
    stamp = json.loads((d / "stamps" / "fit_B29_n20.json").read_text())
    assert stamp["status"] == "complete"
    assert [(s["index"], s["status"], s["nfev"]) for s in stamp["segments"]] == [
        (1, "interrupted", 7), (2, "interrupted", 7), (3, "complete", 6)]
    recorded = json.loads((d / "pipeline.json").read_text())["fits"]["B29_n20"]
    assert [s["nfev"] for s in recorded["segments"]] == [7, 7, 6] and recorded["nfev_used"] == 20
    ft = [c["argv"] for c in calls(world) if c["tool"] == "fine_tune_gpu.py"]
    assert len(ft) == 1 and ft[0][ft[0].index("--restart-from") + 1] == str(d / "fits" / "B29_n20" / "seg03" /
                                                                             "fit_restart.npz")
    capsys.readouterr()
    cp.main(["status", "--restart-dir", str(d)])
    assert "segment 3 complete" in capsys.readouterr().out


def test_failed_segment_runs_again_from_the_last_clean_one(world, monkeypatch):
    monkeypatch.setenv("FAKE_FIT_STEP", "7")
    args = ("--fits", "B29", "--max-nfev", "20", "--wall-time", "01:00:00")
    assert pipeline("run", world, *args) == cp.EXIT_RESUMABLE
    monkeypatch.setenv("FAKE_FAIL", "fit:B29")
    assert pipeline("run", world) == cp.EXIT_FAILED
    monkeypatch.delenv("FAKE_FAIL")
    assert pipeline("run", world) == cp.EXIT_RESUMABLE
    d = world["restart"] / "fits" / "B29_n20"
    assert any(p.name.startswith("seg02.incomplete_") for p in d.iterdir())
    last = [c["argv"] for c in calls(world) if c["tool"] == "gpu_fit.py"][-1]
    assert last[last.index("--segment") + 1] == "2" and last[last.index("--max-nfev") + 1] == "13"
    assert last[last.index("--restart-from") + 1] == str(d / "seg01" / "fit_restart.npz")


def test_no_progress_segment_and_deadline_deferral(world, monkeypatch, capsys):
    monkeypatch.setenv("FAKE_FIT_EXHAUSTED", "1")
    assert pipeline("run", world, "--fits", "B29", "--wall-time", "01:00:00") == cp.EXIT_RESUMABLE
    stamp = json.loads((world["restart"] / "stamps" / "fit_B29_n2000.json").read_text())
    assert stamp["status"] == "deferred" and stamp["segments"] == []
    # a stage whose earlier attempt took longer than the time left is not started
    d = world["restart"]
    (d / "stamps" / "base.json").rename(d / "stamps" / "base.json.incomplete_x")
    data = json.loads((d / "stamps" / "base.json.incomplete_x").read_text())
    data.update(status="interrupted", seconds=1e6)
    (d / "stamps" / "base.json.incomplete_x").write_text(json.dumps(data))
    capsys.readouterr()
    n = len(calls(world))
    assert pipeline("run", world) == cp.EXIT_RESUMABLE
    assert len(calls(world)) == n and "an earlier attempt took 1000000 s" in capsys.readouterr().out


def test_memory_model_of_the_precision_decision():
    gib = 1 << 30
    # 20M on the two local A6000 (headroom 0.12, host activity): FP64 does not fit; FP32 streams a part
    est = cp.memory_estimate(5057053157, 19999504, 2, 1, 48 * gib, 251 * gib, 0.12, "host")
    assert not est["fp64_fits"] and est["gpu_fits"] and est["host_fits"]
    assert 2.9e9 < est["fp32_resident_records"] < 3.6e9            # measured: 2 x 1.768e9 resident
    assert 1.4e9 < est["fp32_streamed_records"] < 2.2e9            # measured: 2 x 0.761e9 streamed
    assert 45e9 < est["host_bytes_rank0_node"] < 80e9              # measured: about 54 GB
    # 8 GPUs of 40 GB: FP32 records all resident, no pinned buffers
    est = cp.memory_estimate(5057053157, 19999504, 8, 1, 40e9, 512e9, 0.20, "device")
    assert not est["fp64_fits"] and est["fp32_streamed_records"] == 0 and est["pinned_bytes_per_node"] == 0
    # small mesh: FP64 resident
    assert cp.memory_estimate(3e8, 1194250, 2, 1, 48 * gib, 251 * gib, 0.2, "host")["fp64_fits"]
    # host memory too small for the pinned buffers: refused
    assert not cp.memory_estimate(5057053157, 19999504, 2, 1, 48 * gib, 32 * gib, 0.12, "host")["host_fits"]


def test_memory_refusal_is_reported_and_stops_the_fit(world, capsys):
    assert pipeline("run", world, "--fits", "B29", "--precision", "mixed") == 0
    capsys.readouterr()
    assert pipeline("plan", world, "--host-memory-gb", "0.001") == 0
    assert "REFUSED       the host memory of a node" in capsys.readouterr().out
    n = len(calls(world))
    assert pipeline("run", world, "--host-memory-gb", "0.001", "--tag", "small") == cp.EXIT_FAILED
    assert len(calls(world)) == n                                  # no fit was started
    stamp = json.loads((world["restart"] / "stamps" / "fit_B29_n2000_small.json").read_text())
    assert stamp["status"] == "failed" and stamp["result"].startswith("refused: the host memory of a node")


def test_status_check_and_lock(world, capsys):
    assert pipeline("run", world, "--fits", "B29") == 0
    import fcntl
    with open(world["restart"] / ".pipeline.lock", "a+") as held:
        fcntl.flock(held.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        assert pipeline("run", world) == cp.EXIT_REFUSED
    assert "another cops_pipeline.py" in capsys.readouterr().err


def test_restart_from_template_checks(world, capsys, tmp_path):
    restart = tmp_path / "nc31_restart.npz"
    meta = {"mesh": {"identity": "z" * 64, "faces": FACES}, "templates": {"names": ["late_ng_not_rg", "cliff_ge_60deg"]}}
    np.savez(restart, meta_json=np.array(json.dumps(meta)))
    assert pipeline("plan", world, "--restart-from", "B29=%s" % restart) == cp.EXIT_REFUSED
    assert "templates ['late_ng_not_rg', 'cliff_ge_60deg'] that B29 does not fit" in capsys.readouterr().err
    assert pipeline("plan", world, "--restart-from", str(restart)) == cp.EXIT_REFUSED
    assert "give FIT=FILE" in capsys.readouterr().err
    assert pipeline("plan", world, "--fits", "NC31", "--restart-from", str(restart)) == 0
    out = capsys.readouterr().out
    assert "NC31: restart mesh zzzzzzzzzzzz; checked against the base when the fit starts" in out
    assert "--restart-from %s" % restart in out


def test_submit_scripts_chain_and_never_run_sbatch(world, capsys):
    write_exe(world["bin"] / "sinfo", "#!/bin/sh\necho 'gpu:a100:4(S:0-1)|64'\n")
    assert cp.main(["submit"] + base_args(world, "--fits", "B29", "--nodes", "3", "--gpus-per-node", "auto",
                                          gpus=None)) == 0
    out = capsys.readouterr().out
    d = world["restart"]
    prepare, fit = (d / "slurm" / "prepare.sbatch").read_text(), (d / "slurm" / "fit.sbatch").read_text()
    assert "GPUs per node from sinfo: 4, 64 CPUs per node" in out
    assert "#SBATCH --nodes=1" in prepare and "#SBATCH --gpus-per-node=4" in prepare
    assert "#SBATCH --nodes=3" in fit and "#SBATCH --ntasks-per-node=4" in fit and "#SBATCH --gpus-per-task=1" in fit
    assert "#SBATCH --cpus-per-task=16" in fit and "##SBATCH --partition=SITE" in fit
    assert "run --restart-dir %s --group fit" % d in fit
    assert "--dependency=afterok:$J_PREPARE" in out and not (world["tmp"] / "sbatch_called").exists()
    assert "afterany" not in fit
    # resubmission needs a wall time; then the scripts continue themselves
    capsys.readouterr()
    assert cp.main(["submit", "--restart-dir", str(d), "--resubmit", "3", "--gpus-per-node", "4"]) == cp.EXIT_REFUSED
    assert "--resubmit needs --wall-time" in capsys.readouterr().err
    assert cp.main(["submit", "--restart-dir", str(d), "--resubmit", "3", "--wall-time", "slurm",
                    "--gpus-per-node", "4"]) == 0
    out = capsys.readouterr().out
    prepare = (d / "slurm" / "prepare.sbatch").read_text()
    assert "sbatch --dependency=afterany:$SLURM_JOB_ID --export=ALL,COPS_PIPELINE_RESUBMIT_LEFT=$((left - 1)) %s" \
        % (d / "slurm" / "prepare.sbatch") in prepare
    assert "status --restart-dir %s --group prepare --check; state=$?" % d in prepare
    assert "sbatch --export=ALL,COPS_PIPELINE_RESUBMIT_LEFT=3 %s" % (d / "slurm" / "fit.sbatch") in prepare
    assert "elif [ $state -eq 75 ]" in prepare and "left=${COPS_PIPELINE_RESUBMIT_LEFT:-3}" in prepare
    assert any(p.name.startswith("prepare.sbatch.incomplete_") for p in (d / "slurm").iterdir())
    assert "afterok" not in out and not (world["tmp"] / "sbatch_called").exists()


def test_submit_without_sinfo_needs_gpus_per_node(world, capsys, monkeypatch):
    monkeypatch.setenv("PATH", str(world["bin"]))                 # no sinfo
    assert cp.main(["submit"] + base_args(world, "--fits", "B29", gpus=None)) == cp.EXIT_REFUSED
    assert "give --gpus-per-node" in capsys.readouterr().err


def test_fit_mpi_template_hardwires_no_ranks():
    text = (APP / "tools" / "cluster" / "fit_mpi.slurm").read_text()
    active = [line for line in text.splitlines() if line.startswith("#SBATCH")]
    assert not any("--ntasks-per-node" in line or "--nodes" in line for line in active)
    assert "SLURM_NTASKS_PER_NODE" in text


def test_hardware_detection(world, monkeypatch):
    args = cp.parse_args(["plan", "--restart-dir", "x"])
    request = {"python": sys.executable}
    monkeypatch.setenv("SLURM_JOB_NUM_NODES", "3")
    monkeypatch.setenv("SLURM_GPUS_ON_NODE", "8")
    monkeypatch.setenv("SLURM_NTASKS", "24")
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "512000")
    monkeypatch.setenv("SLURM_JOB_ID", "5")
    hw = cp.detect_hardware(args, request, live=False)
    assert (hw["nodes"], hw["gpus_per_node"], hw["ranks"], hw["launcher"]) == (3, 8, 24, "srun")
    assert hw["host_memory_bytes"] == 512000 * (1 << 20) and hw["gpu_memory_bytes"] is None
    args = cp.parse_args(["plan", "--restart-dir", "x", "--nodes", "2", "--gpus-per-node", "4",
                          "--gpu-memory-gb", "80", "--host-memory-gb", "1000"])
    hw = cp.detect_hardware(args, request, live=False)
    assert (hw["nodes"], hw["gpus_per_node"], hw["ranks"]) == (2, 4, 8)
    assert hw["gpu_memory_bytes"] == 80e9 and hw["host_memory_bytes"] == 1000e9
    for name in ("SLURM_JOB_NUM_NODES", "SLURM_GPUS_ON_NODE", "SLURM_NTASKS", "SLURM_MEM_PER_NODE", "SLURM_JOB_ID"):
        monkeypatch.delenv(name)
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "2,3")
    hw = cp.detect_hardware(cp.parse_args(["plan", "--restart-dir", "x"]), request, live=False)
    assert hw["gpus_per_node"] == 2 and hw["gpus_source"] == "CUDA_VISIBLE_DEVICES"
    assert hw["host_memory_source"] == "/proc/meminfo" and hw["host_memory_bytes"] > 0


@pytest.mark.parametrize("chunk", [7, 8, 13, 64, 1 << 20])
def test_facet_count_is_exact_across_read_blocks(tmp_path, monkeypatch, chunk):
    """A card ending exactly at a block boundary is counted once (the 20M BDF read 19,999,506 before)."""
    mesh = tmp_path / "mesh.bdf"
    cards = "".join("CTRIA3 %d 1 1 2 3\n" % k for k in range(50))
    mesh.write_text("$ header\nGRID 1 0 0. 0. 0.\n" + cards + "ENDDATA\n")
    monkeypatch.setattr(cp, "CHUNK", chunk)
    assert cp.file_identity(mesh, count_facets=True)["facets"] == 50
