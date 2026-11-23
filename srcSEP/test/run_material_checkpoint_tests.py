#!/usr/bin/env python3
"""Real-MPI material restart and repartition gate, using a frozen native build.

First prepare the ordinary native integration runner with --prepare-only.
This runner consumes that compiler/object/header/archive provenance rather than
inventing a second build configuration. --prepare-only emits dependent jobs;
--collect requires receipts containing their actual process return codes.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import re
import os
import signal
import struct
import subprocess
import time
import xml.etree.ElementTree as ET

MOVERS = {"MCP01": "parker", "MCP02": "fte-dmumu", "MCP03": "fte-mfp"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def physical(condition, message):
    """Separate a computed contract failure from missing execution evidence."""
    if not condition:
        raise AssertionError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fields(text, prefix):
    return [dict(re.findall(r"([a-z_]+)=([^\s]+)", line))
            for line in text.splitlines() if line.startswith(prefix)]


def particles(directory, ranks):
    """Compare every persisted physical/history byte, excluding only links."""
    result = {}
    for rank in range(ranks):
        raw = (directory / f"material-particles-rank{rank}.bin").read_bytes()
        require(len(raw) >= 8, "missing particle-record header")
        width, = struct.unpack_from("<Q", raw)
        require(width > 0 and (len(raw)-8) % (width+8) == 0, "truncated particle history")
        for offset in range(8, len(raw), width+8):
            identity, = struct.unpack_from("<Q", raw, offset)
            require(identity != 0 and identity not in result, "duplicate physical particle identity")
            result[identity] = raw[offset+8:offset+8+width]
    require(result, "empty population cannot qualify stochastic-history continuation")
    return result


def ledger(text, steps):
    """Recompute signed SI balances independently of printed residual fields."""
    rows = fields(text, "REDUCED_SOURCE_LEDGER ")
    require(len(rows) == len(steps), "missing physical-step ledger")
    maximum = 0.0
    for row, step in zip(rows, steps):
        values = {key: float(value) for key, value in row.items()}
        require(all(math.isfinite(value) for value in values.values()), "nonfinite ledger")
        physical(values["epoch_s"] == 12000 + 300*step, "wrong physical clock or skipped step")
        born, energy = values["born"], values["born_kinetic_j"]
        require(born > 0 and energy > 0, "missing physical accelerated source")
        terms = [values[name] for name in ("retained", "returned", "domain_exit", "geometry_retired")]
        birth = [values[name] for name in ("retained_birth_j", "returned_birth_j", "domain_birth_j", "retired_birth_j")]
        events = [values[name] for name in ("retained_kinetic_j", "returned_event_j", "domain_event_j", "retired_event_j")]
        work = values["background_work_j"] + values["wave_work_j"]
        residuals = (abs(born-math.fsum(terms))/born,
                     abs(energy-math.fsum(birth))/energy,
                     abs(math.fsum([energy, work, -math.fsum(events)])) /
                     (energy+abs(values["background_work_j"])+abs(values["wave_work_j"])))
        maximum = max(maximum, *residuals)
        physical(max(residuals) <= 1e-10, "represented number/birth/signed energy budget failed")
    return rows[-1], maximum


def compare_output(reference, actual, repartition):
    physical(reference.with_suffix(reference.suffix+".sampling").read_bytes() ==
            actual.with_suffix(actual.suffix+".sampling").read_bytes(),
            "material labels/exposure/cadence changed after restart")
    a, b = reference.read_text(), actual.read_text()
    if not repartition:
        physical(a == b, "same-decomposition completed output differs")
        return 0.0
    # Native Tecplot prints six decimal places in scientific notation. A
    # changed MPI summation tree can cross the last printed rounding boundary;
    # the fixed 2e-6 relative bound covers two independently rounded values.
    # This serialization bound does not replace the 1e-10 SI budget gate.
    left = [line.split() for line in a.splitlines() if line and line[0] in "-+.0123456789"]
    right = [line.split() for line in b.splitlines() if line and line[0] in "-+.0123456789"]
    require(len(left) == len(right) and left, "missing/repartitioned output rows")
    maximum = 0.0
    for ra, rb in zip(left, right):
        require(len(ra) == len(rb), "output schema changed with MPI decomposition")
        for sa, sb in zip(ra, rb):
            x, y = float(sa), float(sb)
            if math.isnan(x) or math.isnan(y):
                require(math.isnan(x) and math.isnan(y), "vacuum/unmeasured distinction changed")
            else:
                require(math.isfinite(x) and math.isfinite(y), "nonfinite measured output")
                scale = max(abs(x), abs(y))
                maximum = max(maximum, abs(x-y)/scale if scale else 0.0)
    physical(maximum <= 2e-6, "repartitioned output exceeds native text precision")
    return maximum


def prepare(native, output):
    # Reuse the actual native executable and its frozen build hashes. The
    # continuous and restarted jobs must differ only in the checkpoint
    # handoff and, for the explicit repartition control, the MPI layout.
    frozen = json.loads(native.read_text())
    hashes = dict(frozen["hashes"])
    hashes[str(Path(__file__).resolve())] = digest(__file__)
    first = frozen["jobs"][0]
    executable = first["argv"][3]
    deck, event = first["argv"][4:6]
    ranks = frozen["ranks"]
    require(set((1, 2, 4, 8)).issubset(ranks), "site gate requires MPI 1, 2, 4 and 8")
    jobs = []
    for identity, mover in MOVERS.items():
        for count in ranks:
            for role in ("continuous", "write", "resume"):
                directory = output / f"{identity}-ranks{count}-{role}"
                directory.mkdir()
                checkpoint = output / f"{identity}-ranks{count}-write/material.checkpoint"
                command = ["mpiexec", "-n", str(count), "--bind-to", "none", executable,
                           deck, event, f"checkpoint-{mover}-{role}", str(checkpoint)]
                jobs.append(dict(id=identity, mover=mover, ranks=count, role=role,
                                 cwd=str(directory), argv=command, timeout_s=300))
        for writer, reader in ((1, 8), (8, 1)):
            directory = output / f"{identity}-repartition{writer}to{reader}"
            directory.mkdir()
            checkpoint = output / f"{identity}-ranks{writer}-write/material.checkpoint"
            command = ["mpiexec", "-n", str(reader), "--bind-to", "none", executable,
                       deck, event, f"checkpoint-{mover}-resume", str(checkpoint)]
            jobs.append(dict(id="MCP0"+str(int(identity[-1])+3), mover=mover, ranks=reader,
                             writer_ranks=writer, role="repartition", cwd=str(directory),
                             argv=command, timeout_s=300))
        # A separate resume from the continuous process's step-eight save
        # restores a different native segment count and then advances through
        # the third output. This supplements the odd-buffer-phase step-four
        # saves; equality of positions alone cannot qualify changing topology.
        for count in ranks:
            directory = output/f"{identity}-topology-ranks{count}"
            directory.mkdir()
            checkpoint = output/f"{identity}-ranks{count}-continuous/material-topology.checkpoint"
            jobs.append(dict(id=f"MCP{int(identity[-1])+8:02d}", mover=mover, ranks=count,
                             role="topology", cwd=str(directory), timeout_s=300,
                             argv=["mpiexec", "-n", str(count), "--bind-to", "none", executable,
                                   deck, event, f"checkpoint-{mover}-resume8", str(checkpoint)]))
    for count in (1, 8):
        # A CRC-valid populated save reaches allocation/table preparation.
        # This complements early corrupt-file rejection; each MPI job injects
        # exactly one late local failure and requires atomic all-rank rollback.
        directory = output / f"MCP12-ranks{count}"
        directory.mkdir()
        jobs.append(dict(id="MCP12", ranks=count, role="late-reject", cwd=str(directory),
                         argv=["mpiexec", "-n", str(count), "--bind-to", "none", executable,
                               deck, event, "checkpoint-parker-late-reject",
                               str(output/"MCP01-ranks1-write/material.checkpoint")],timeout_s=300))
        directory = output / f"MCP07-ranks{count}"
        directory.mkdir()
        jobs.append(dict(id="MCP07", ranks=count, role="reject", cwd=str(directory),
                         argv=["mpiexec", "-n", str(count), "--bind-to", "none", executable,
                               deck, event, "checkpoint-parker-reject", str(output/"corrupt.checkpoint")],
                         timeout_s=300))
    plan = dict(schema="srcsep-material-checkpoint-plan-v1", hashes=hashes, jobs=jobs)
    (output/"prepared.json").write_text(json.dumps(plan, indent=2)+"\n")
    return plan


def collect(plan, output):
    for path, expected in plan["hashes"].items():
        require(digest(path) == expected, "source/configuration/object/archive identity changed: "+path)
    results = []
    for job in plan["jobs"]:
        start = time.monotonic()
        row = dict(job, status="ERROR", metrics={})
        directory = Path(job["cwd"])
        try:
            receipt = json.loads((directory/"process-result.json").read_text())
            require(receipt["argv"] == job["argv"], "receipt command differs from prepared command")
            require(type(receipt["exit_code"]) is int and receipt["exit_code"] == 0,
                    "MPI process failed or has no actual return code")
            require(not receipt.get("launch_error"), "MPI launch/timeout error: "+str(receipt.get("launch_error")))
            text = (directory/"native.rank0.log").read_text()
            marker = job["id"] if job["role"] in ("reject","late-reject") else next(key for key, value in MOVERS.items() if value == job["mover"])
            require(f"[{marker}] PASS" in text and " ERROR " not in text, "missing native PASS or native error")
            if job["role"] == "late-reject":
                late=fields(text,"MATERIAL_LATE_REJECTION ")
                physical(len(late)==1 and int(late[0]["prepared_particles"])>0 and
                         int(late[0]["failed_rank"])==job["ranks"]-1 and int(late[0]["retry_step"])==4,
                         "late rejection did not exercise staged particles, one failing owner and successful retry")
                row["metrics"]["prepared_particles"]=int(late[0]["prepared_particles"])
            elif job["role"] != "reject":
                steps = range(1, 5) if job["role"] == "write" else [9] if job["role"] == "topology" else range(5, 10) if job["role"] in ("resume", "repartition") else range(1, 10)
                _, maximum = ledger(text, steps)
                row["metrics"]["maximum_relative_budget_error"] = maximum
                actual = particles(directory, job["ranks"])
                row["metrics"]["retained_macroparticles"] = len(actual)
                if job["role"] in ("resume", "repartition", "topology"):
                    identity = next(key for key, value in MOVERS.items() if value == job["mover"])
                    reference = output/f"{identity}-ranks{job['ranks']}-continuous"
                    physical(actual == particles(reference, job["ranks"]), "particle state/history differs from continuous execution")
                    repartition = job["role"] == "repartition"
                    if not repartition:
                        for rank in range(job["ranks"]):
                            physical((directory/f"material-final-rank{rank}.bin").read_bytes() ==
                                    (reference/f"material-final-rank{rank}.bin").read_bytes(),
                                    "same-decomposition complete native continuation image differs")
                    number = 2 if job["role"] == "topology" else 1
                    name = f"material-output-{number}.dat"
                    error = compare_output(reference/name, directory/name, repartition)
                    if job["role"] == "topology":
                        restored = fields(text, "MATERIAL_TOPOLOGY_RESTORE ")
                        require(len(restored)==1 and int(restored[0]["original_segments"])!=int(restored[0]["restored_segments"]),
                                "missing changed-topology execution proof")
                        # Only output two is generated after a step-eight
                        # restart; previous outputs belong to the prior job.
                        row["metrics"]["original_segment_count"] = int(restored[0]["original_segments"])
                        row["metrics"]["restored_segment_count"] = int(restored[0]["restored_segments"])
                    row["metrics"]["maximum_output_relative_difference"] = error
            row["status"] = "PASS"
        except (OSError, ValueError, KeyError, struct.error) as error:
            row["message"] = str(error)
        except AssertionError as error:
            row["status"],row["message"] = "FAIL",str(error)
        row["duration_s"] = time.monotonic()-start
        results.append(row)
    report = dict(schema="srcsep-material-checkpoint-v1", results=results, hashes=plan["hashes"],
                  acceptance=dict(physical_dt_s=300, same_decomposition="exact complete image",
                                  repartition="exact particle histories and exposure; SI balances <=1e-10",
                                  output_text_relative_tolerance=2e-6,
                                  reference="separate uninterrupted native process; no regenerated reference solution"))
    (output/"summary.json").write_text(json.dumps(report, indent=2)+"\n")
    suite = ET.Element("testsuite", name="material-checkpoint", tests=str(len(results)))
    for row in results:
        case = ET.SubElement(suite, "testcase", classname=row["id"], name=Path(row["cwd"]).name,
                             time=str(row["duration_s"]))
        if row["status"] != "PASS":
            ET.SubElement(case, "failure" if row["status"] == "FAIL" else "error",
                          message=row.get("message", "missing evidence"))
    ET.ElementTree(suite).write(output/"junit.xml", encoding="utf-8", xml_declaration=True)
    print("Material checkpoint:", {status: sum(row["status"] == status for row in results)
                                   for status in ("PASS", "FAIL", "SKIP", "ERROR")})
    return int(any(row["status"] != "PASS" for row in results))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native-prepared", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--collect", action="store_true")
    parser.add_argument("--prepare-corruption", action="store_true")
    args = parser.parse_args()
    output = args.output.resolve()
    if args.prepare_corruption:
        raw = bytearray((output/"MCP01-ranks1-write/material.checkpoint").read_bytes())
        require(len(raw) > 128, "missing validated writer checkpoint")
        raw[-1] ^= 1  # Deliberate CRC failure, never a modified production file.
        require(not (output/"corrupt.checkpoint").exists(), "negative-control evidence already exists")
        (output/"corrupt.checkpoint").write_bytes(raw)
        return 0
    if args.collect:
        return collect(json.loads((output/"prepared.json").read_text()), output)
    require(args.native_prepared, "supply --native-prepared from the ordinary native integration runner")
    require(not output.exists() or not any(output.iterdir()), "choose a fresh output directory")
    output.mkdir(parents=True, exist_ok=True)
    plan = prepare(args.native_prepared.resolve(), output)
    if args.prepare_only:
        print("Prepared dependent MPI jobs; no execution or PASS claimed:", output/"prepared.json")
        return 0
    for job in plan["jobs"]:
        if job["role"] == "reject" and not (output/"corrupt.checkpoint").exists():
            raw = bytearray((output/"MCP01-ranks1-write/material.checkpoint").read_bytes())
            raw[-1] ^= 1
            (output/"corrupt.checkpoint").write_bytes(raw)
        directory = Path(job["cwd"])
        launch_error = ""
        with (directory/"native.rank0.log").open("w") as log:
            try:
                process = subprocess.Popen(job["argv"], cwd=directory, stdout=log,
                                           stderr=subprocess.STDOUT, start_new_session=True)
                try:
                    code = process.wait(timeout=job["timeout_s"])
                except subprocess.TimeoutExpired:
                    # Terminate the MPI launcher AND its ranks before another
                    # case starts. Retain the actual reaped process code; a
                    # timeout is an execution error, never a physics pass.
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
                    code,launch_error = process.returncode,"MPI job exceeded the declared 300-s deadline"
            except OSError as error:
                code,launch_error = None,str(error)
        (directory/"process-result.json").write_text(json.dumps(
            dict(argv=job["argv"], exit_code=code, launch_error=launch_error))+"\n")
    return collect(plan, output)


if __name__ == "__main__":
    raise SystemExit(main())
