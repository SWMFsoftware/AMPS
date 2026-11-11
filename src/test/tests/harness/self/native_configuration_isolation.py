"""Control-plane regression of native adapters using isolated external tools.

The child adapters exercise their real configuration, verification and result
paths. Their compiler/MPI stand-ins never compile or execute AMPS physics.
Synthetic native markers are private fixture inputs, not globally reported
physics-case results. Real native qualification is run separately.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from typing import Dict


TOOL = r'''import json, os, shutil, subprocess, sys
from pathlib import Path
root=Path(os.environ["AMPS_NATIVE_ISOLATION_FIXTURE"])
spec=json.loads((root/"fixture.json").read_text())
name=Path(sys.argv[0]).name
if name=="Config.pl": role,args="configure",sys.argv[1:]
elif name.startswith("amps"): role,args="native",sys.argv[1:]
else: role,args=sys.argv[1],sys.argv[2:]
with (root/"events.jsonl").open("a") as stream:
 stream.write(json.dumps({"role":role,"argv":args})+"\n")
def selection():
 return dict(line.split("=",1) for line in (root/".amps.conf").read_text().splitlines() if "=" in line)
variant=spec["variant"]
if role=="configure":
 if variant=="configure-failure":
  print("SYNTHETIC_CONFIG_FAILURE: deliberately refused selection");sys.exit(7)
 state=selection()
 for arg in args:
  if arg.startswith("-application="):
   app=arg.split("=",1)[1].rsplit("/",1)[-1]
   state["APPLICATION"]=app if variant!="wrong-application" else "wrong-application"
   state["InputFileAMPS"]=app+".input"
   shutil.copyfile(root/spec["deck"],root/(app+".input"))
  if arg.startswith("-spice-path=") and variant!="ignored-mode":state["SPICE"]=arg.split("=",1)[1]
 (root/".amps.conf").write_text("".join(k+"="+v+"\n" for k,v in state.items()))
 print("SYNTHETIC configuration selection");sys.exit(0)
if role=="make":
 state=selection()
 if any(arg.startswith("amps_test_native_configuration_") for arg in args):
  if variant=="probe-failure":print("SYNTHETIC_PROBE_FAILURE");sys.exit(9)
  mode=state["SPICE"] if variant!="contradictory-make" else "/synthetic/other-spice"
  flags="-D _NO_SPICE_CALLS_" if mode=="nospice" else "-I/synthetic/spice"
  if variant=="compiler-undefines":flags+=" -U_NO_SPICE_CALLS_"
  print("AMPS_NATIVE_CONFIG_SPICE="+mode)
  print("AMPS_NATIVE_CONFIG_SEARCH_C="+flags)
  print("AMPS_NATIVE_CONFIG_FLAGCC="+flags+" -c -g -std=c++17")
  sys.exit(0)
 if args!=["amps"]:raise RuntimeError("unexpected synthetic make target")
 if state["SPICE"]!="nospice":print("SYNTHETIC_BUILD_WRONG_MODE");sys.exit(8)
 with (root/"events.jsonl").open("a") as stream:stream.write(json.dumps({"role":"build"})+"\n")
 flag="" if variant=="wrong-compile-mode" else "-D _NO_SPICE_CALLS_"
 print("fixture-cxx "+flag+" -c synthetic.cpp -o synthetic.o")
 if variant=="unexpected-library":print("fixture-link synthetic.o /synthetic/cspice/lib/cspice.a")
 shutil.copyfile(root/"tool.py",root/"amps")
 (root/"amps").chmod(0o700)
 sys.exit(0)
if role=="mpiexec":
 if args[:2]!=["-n","1"]:raise RuntimeError("fixture must remain one rank")
 sys.exit(subprocess.call(args[2:]))
if role=="native":
 print("SYNTHETIC CONTROL-PLANE RECEIPT; NO AMPS PHYSICS EXECUTED")
 if spec["kind"]=="exosphere":
  print("exosphere_source_counts impact=1 psd=1 thermal=1 sputtering=1 vertical=1")
 if spec["kind"]=="dust":
  out=Path(os.environ["AMPS_TEST_MODEL_DUST_OUTPUT"])
  out.write_text("SYNTHETIC HARNESS DATA\n"+"".join("Charging of a grain with radius "+r+"\nEquilibrium charge: 1\n" for r in ["1e-7","1e-6","1e-5","1e-4"]))
 print("["+spec["case"]+"] PASS")
 sys.exit(0)
raise RuntimeError("unknown synthetic tool role")
'''


FIXTURES = [
    # These trailing strings are synthetic artifact labels, not physical-case
    # descriptions. Keep a comment between case ID and label so the catalog's
    # prose-tuple recognizer cannot attribute mock metadata to a real oracle.
    ("exosphere", "models/exosphere_native", "test/core-exosphere-sphere",
     "CORE-EXOSPHERE-SPHERE-NATIVE-001", "CORE-EXOSPHERE-SPHERE-INTEGRATION-001",  # label below
     "exosphere-native"),
    ("dust", "models/dust_native", "test/model-dust", "CORE-DUST-NATIVE-001",
     "CORE-DUST-CHARGING-DYNAMICS-001",  # label below
     "dust-native"),
    ("dust-injection", "models/dust_injection_native", "test/model-dust-injection",
     "CORE-DUST-INJECTION-NATIVE-001", "CORE-DUST-INJECTION-INTEGRATION-001",  # label below
     "dust-injection-native"),
    ("stencil", "pic/ecsim_operators", "test/stencil", "CORE-ECSIM-OPERATORS-001",
     "ECSIM-STENCIL-CURL-001",  # label below
     "ecsim-operators"),
]


def test_native_configuration_isolation(context: Dict) -> None:
    """Exercise state leakage and rejected mismatches through real adapters.

    Every row has its own fake repository and external-tool state. A fresh
    SPICE-enabled start before dust is essential: using the exosphere result
    as its initial state would let an uncorrected dust adapter pass unnoticed.
    Success here proves control flow/evidence only, not physical correctness.
    """
    repo = Path(context["repo_root"])
    tests = repo / "src/test/tests"
    artifacts = Path(context["artifact_dir"]) / "native-configuration-isolation"
    artifacts.mkdir(parents=True, exist_ok=True)
    rows = []
    real_selection = {name: (repo / name).read_bytes()
                      for name in (".amps.conf", "Makefile.conf", "Makefile.local")
                      if (repo / name).is_file()}
    scenarios = [("enabled", "/synthetic/previous-spice", "ok"),
                 ("disabled", "nospice", "ok")]
    negative = ["ignored-mode", "contradictory-make", "compiler-undefines",
                "configure-failure", "probe-failure", "wrong-application",
                "wrong-compile-mode", "unexpected-library"]
    for kind, runner, application, suite, identifier, label in FIXTURES:
        cases = list(scenarios)
        if kind in ("exosphere", "dust"):
            cases += [(variant, "/synthetic/previous-spice", variant)
                      for variant in negative]
        for name, initial, variant in cases:
            evidence = artifacts / (kind + "-" + name)
            evidence.mkdir()
            # TemporaryDirectory owns only this row, never the real checkout.
            # Copy its resulting logs/receipts before automatic fixture cleanup.
            with TemporaryDirectory(prefix="amps-native-mode-", dir=context["work_dir"]) as tmp:
                private = Path(tmp)
                (private / "src").mkdir()
                (private / "input/test").mkdir(parents=True)
                (private / "work").mkdir()
                (private / "artifacts").mkdir()
                deck = "input/" + application + ".input"
                (private / deck).write_bytes((repo / deck).read_bytes())
                for filename in ("Makefile", "Makefile.conf", "Makefile.def", "Makefile.local"):
                    (private / filename).write_text("# isolated synthetic build fixture\n")
                (private / ".amps.conf").write_text(
                    "APPLICATION=previous\nInputFileAMPS=previous.input\nSPICE=" + initial + "\n")
                spec = {"kind": kind, "case": identifier, "deck": deck, "variant": variant}
                (private / "fixture.json").write_text(json.dumps(spec))
                code = "#!" + sys.executable + "\n" + TOOL
                for filename in ("Config.pl", "tool.py"):
                    (private / filename).write_text(code)
                    (private / filename).chmod(0o700)
                child = {
                    "schema_version": 1, "repo_root": str(private),
                    "test": {"id": suite, "expected_cases": [identifier],
                             "selected_cases": [identifier], "timeout_seconds": 10},
                    "work_dir": str(private / "work"),
                    "artifact_dir": str(private / "artifacts"),
                    "result_file": str(private / "result.json"), "paths": {},
                    "tools": {"make": [sys.executable, str(private / "tool.py"), "make"],
                              "mpiexec": [sys.executable, str(private / "tool.py"), "mpiexec"]},
                    "execution": {"mpi_rank_flag": "-n"},
                }
                config = private / "context.json"
                config.write_text(json.dumps(child))
                environment = dict(os.environ, AMPS_TEST_CONTEXT=str(config),
                                   AMPS_NATIVE_ISOLATION_FIXTURE=str(private))
                argv = [sys.executable, "-B", str(tests / runner / "run_test.py"),
                        "--context", str(config)]
                completed = subprocess.run(argv, cwd=private, env=environment,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    universal_newlines=True, timeout=15)
                (evidence / "adapter.log").write_text(completed.stdout)
                for path in (private / "artifacts").iterdir():
                    if path.is_file():
                        (evidence / path.name).write_bytes(path.read_bytes())
                assert (private / "result.json").is_file(), (
                    kind, variant, completed.returncode, completed.stdout)
                result = json.loads((private / "result.json").read_text())
                (evidence / "synthetic-result.json").write_text(json.dumps(result, indent=2))
                mode = json.loads((private / "artifacts" / (label + "-configuration.json")).read_text())
                events = [json.loads(line) for line in (private / "events.jsonl").read_text().splitlines()]
                (evidence / "events.json").write_text(json.dumps(events, indent=2))
                expected = "PASS" if variant == "ok" else "ERROR"
                assert result["status"] == expected, (kind, variant, result, completed.stdout)
                assert [c["id"] for c in result["cases"]] == [identifier]
                assert all(c["status"] == expected for c in result["cases"])
                assert (completed.returncode == 0) == (expected == "PASS")
                assert mode["configure_command"][-1] == "-spice-path=nospice"
                assert mode["before_selection"]["SPICE"] == initial
                assert all(Path(path).is_file() for path in result["evidence"])
                executed = [e for e in events if e["role"] == "native"]
                built = [e for e in events if e["role"] == "build"]
                if expected == "PASS":
                    assert len(executed) == len(built) == 1
                    assert all(mode["checks"].values()) and all(mode["build_checks"].values())
                    identity = json.loads((private / "artifacts" /
                                          (label + "-build-identity.json")).read_text())
                    assert identity["configure_command"] == mode["configure_command"]
                    assert identity["build_environment"] == {"MAKEFLAGS": "-j16"}
                else:
                    assert not executed, (kind, variant, events)
                    assert len(built) == (1 if variant in
                                         ("wrong-compile-mode", "unexpected-library") else 0)
                    if variant == "configure-failure":
                        assert "SYNTHETIC_CONFIG_FAILURE" in (
                            private / "artifacts" / (label + "-configuration.txt")).read_text()
                rows.append({"fixture": kind, "initial_spice": initial,
                             "scenario": variant, "expected_status": expected,
                             "observed_status": result["status"],
                             "native_physics_executed": False})
    assert all((repo / name).read_bytes() == value for name, value in real_selection.items())
    (artifacts / "summary.json").write_text(json.dumps({
        "scope": "synthetic configuration-control checks only",
        "native_physics_executed": False, "rows": rows, "status": "PASS",
    }, indent=2) + "\n")
