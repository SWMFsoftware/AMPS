#!/usr/bin/env python3
"""Unified srcSEP3D test runner.

The command vocabulary intentionally follows ``srcSEP/test/run_tests.py``:
``--list``, repeatable ``--test``/``--group``/``--suite``, ``--routine``,
``--all``, and ``--output-dir`` have the same meaning.  This runner joins
the foundation plus Phase-M/B/T/P/A/O/V evidence classes without pretending they are
interchangeable:

* standalone C++ tests compile core/mesh/background/turbulence/runtime without
  AMPS or MPI;
* source/ABI checks inspect the actual production manifest and AMPS pic.h;
* R1 gates audit the canonical sep_common/SWCME archives and relocated runner;
* the strict production build runs only with a real configured AMPS checkout.

Missing production configuration is reported as SKIP, never as PASS.  That
distinction keeps a source-only development run useful while ensuring it
cannot satisfy the AMPS production-build acceptance gate.
"""

from __future__ import annotations

import argparse
import ast
import datetime as _datetime
import hashlib
from dataclasses import asdict, dataclass
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
from xml.sax.saxutils import escape as _xml_escape


ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "test" / "stage1"


@dataclass(frozen=True)
class TestDefinition:
    test_id: str
    group: str
    name: str
    kind: str
    routine: bool = True


@dataclass
class Result:
    test_id: str
    group: str
    status: str
    message: str
    elapsed_seconds: float
    command: List[str]


# The C++ registry remains authoritative for callbacks and result semantics.
# This compact manifest exists only so --list and unknown-ID validation work
# before a binary is built, matching srcSEP's discovery-first CLI behavior.
TESTS: Tuple[TestDefinition, ...] = (
    TestDefinition("BLD01", "BLD", "Standalone binary excludes AMPS/MPI symbols", "cpp"),
    TestDefinition("HARN01", "HARN", "Empty-registry contract", "cpp"),
    TestDefinition("HARN02", "HARN", "Failure/error exit-code contract", "cpp"),
    TestDefinition("HARN03", "HARN", "Skip exit-code contract", "cpp"),
    TestDefinition("HARN04", "HARN", "JSON/JUnit writer contract", "cpp"),
    TestDefinition("LAY01", "LAY", "L0/L1 AMPS dependency exclusion", "cpp"),
    TestDefinition("LAY02", "LAY", "Layering negative control", "cpp"),
    TestDefinition("LIFE3D01", "LIFE3D", "Legal lifecycle and counters", "cpp"),
    TestDefinition("LIFE3D02", "LIFE3D", "Illegal transition matrix", "cpp"),
    TestDefinition("LIFE3D03", "LIFE3D", "Frozen layout and fingerprint", "cpp"),
    TestDefinition("LIFE3D04", "LIFE3D", "Standalone/SWMF adapter parity", "cpp"),
    TestDefinition("R3D01", "R3D", "AMPS mover selection hook", "cpp"),
    TestDefinition("R3D02", "R3D", "Complete requested-time advance", "cpp"),
    TestDefinition("R3D03", "R3D", "Transactional snapshot update", "cpp"),
    TestDefinition("R3D04", "R3D", "Authoritative integer clock", "cpp"),
    TestDefinition("R3D05", "R3D", "Shock source lifecycle", "cpp"),
    TestDefinition("R3D06", "R3D", "Observer publication transaction", "cpp"),
    TestDefinition("R3D07", "R3D", "Complete restart contract", "cpp"),
    TestDefinition("R3D08", "R3D", "Canonical initialization source", "cpp"),
    TestDefinition("R3D09", "R3D", "Finite empty-cell Tecplot output", "cpp"),
    TestDefinition("CFG3D01", "CFG3D", "Input schema and CLI", "cpp"),
    TestDefinition("CFG3D02", "CFG3D", "Complete typed contracts", "cpp"),
    TestDefinition("CFG3D03", "CFG3D", "Domain and boundaries", "cpp"),
    TestDefinition("CFG3D04", "CFG3D", "Shared Parker geometry", "cpp"),
    TestDefinition("CFG3D05", "CFG3D", "Mesh and memory preflight", "cpp"),
    TestDefinition("CFG3D06", "CFG3D", "Finite Parker initialization schema", "cpp"),
    TestDefinition("CFG3D07", "CFG3D", "Complete compiled AMPS species binding", "cpp"),
    TestDefinition("CFG3D08", "CFG3D", "Complete schema-3 initialization", "cpp"),
    TestDefinition("CFG3D09", "CFG3D", "Background and turbulence selection", "cpp"),
    TestDefinition("CFG3D10", "CFG3D", "CME and Parker start linkage", "cpp"),
    TestDefinition("CFG3D11", "CFG3D", "Transport/control schema", "cpp"),
    TestDefinition("CFG3D12", "CFG3D", "Corner/sphere input and endpoint normalization", "cpp"),
    TestDefinition("DOM3D01", "DOM3D", "Whole-corridor corner bounds", "source"),
    TestDefinition("DOM3D02", "DOM3D", "Solar sphere and photospheric coarsening", "source"),
    TestDefinition("DOM3D03", "DOM3D", "Connected sphere/corridor AMR allocation", "source"),
    TestDefinition("DOM3D04", "DOM3D", "X-y corner with Sun-centered z and complete corridor", "source"),
    TestDefinition("OUT3D01", "OUT3D", "Empty output stencil and byte-slice safety", "source"),
    TestDefinition("OUT3D02", "OUT3D", "Excluded-volume output rows and MPI callback branches", "source"),
    TestDefinition("MSH3D01", "MSH3D", "Resolution bounds", "cpp"),
    TestDefinition("MSH3D02", "MSH3D", "Radial closed forms", "cpp"),
    TestDefinition("MSH3D03", "MSH3D", "Parker tube centreline", "cpp"),
    TestDefinition("MSH3D04", "MSH3D", "Tube-distance convergence", "cpp"),
    TestDefinition("MSH3D05", "MSH3D", "Shoulder and 2:1 balance", "cpp"),
    TestDefinition("MSH3D06", "MSH3D", "Rotation invariance", "cpp"),
    TestDefinition("MSH3D07", "MSH3D", "Octree budget and ownership", "cpp"),
    TestDefinition("MSH3D08", "MSH3D", "Earth and Mars presets", "cpp"),
    TestDefinition("MSH3D09", "MSH3D", "Refinement-boundary gradients", "cpp"),
    TestDefinition("MSH3D10", "MSH3D", "Finite Parker line and origin-relative AMR", "cpp"),
    TestDefinition("MSH3D11", "MSH3D", "Initialization Tecplot output", "cpp"),
    TestDefinition("MSH3D12", "MSH3D", "Active Parker corridor", "cpp"),
    TestDefinition("MSH3D13", "MSH3D", "Parker geometry authority", "cpp"),
    TestDefinition("MSH3D14", "MSH3D", "Hole-free active mask", "cpp"),
    TestDefinition("MSH3D15", "MSH3D", "Solar internal-boundary geometry", "cpp"),
    TestDefinition("BGP3D01", "BGP3D", "Divergence-free Parker field", "cpp"),
    TestDefinition("BGP3D02", "BGP3D", "Parker component laws", "cpp"),
    TestDefinition("BGP3D03", "BGP3D", "Field-line tangency", "cpp"),
    TestDefinition("BGP3D04", "BGP3D", "Focusing length", "cpp"),
    TestDefinition("BGP3D05", "BGP3D", "Velocity derivatives", "cpp"),
    TestDefinition("BGP3D06", "BGP3D", "Polar limits", "cpp"),
    TestDefinition("BGP3D07", "BGP3D", "SWCME ambient closure", "cpp"),
    TestDefinition("SNAP3D01", "SNAP3D", "Required-field completeness", "cpp"),
    TestDefinition("SNAP3D02", "SNAP3D", "Finite-value policy", "cpp"),
    TestDefinition("SNAP3D03", "SNAP3D", "SWMF unit conversion", "cpp"),
    TestDefinition("SNAP3D04", "SNAP3D", "SWMF epoch consistency", "cpp"),
    TestDefinition("SNAP3D05", "SNAP3D", "Atomic snapshot publication", "cpp"),
    TestDefinition("SNAP3D06", "SNAP3D", "Snapshot time interpolation", "cpp"),
    TestDefinition("SNAP3D07", "SNAP3D", "Batch per-sample status", "cpp"),
    TestDefinition("SNAP3D08", "SNAP3D", "Coordinate-frame rejection", "cpp"),
    TestDefinition("TUR3D01", "TUR3D", "Prescribed spectrum normalization", "cpp"),
    TestDefinition("TUR3D02", "TUR3D", "AWSoM wave mapping", "cpp"),
    TestDefinition("TUR3D03", "TUR3D", "Resonance bounds", "cpp"),
    TestDefinition("TUR3D04", "TUR3D", "Missing turbulence policy", "cpp"),
    TestDefinition("TUR3D05", "TUR3D", "Selectable prescribed models", "cpp"),
    TestDefinition("TUR3D06", "TUR3D", "Amplitude models and Tecplot energy", "cpp"),
    TestDefinition("COEF3D01", "COEF3D", "Coefficient conversions", "cpp"),
    TestDefinition("COEF3D02", "COEF3D", "Shared coefficient kernel", "cpp"),
    TestDefinition("COEF3D03", "COEF3D", "Parallel tensor assembly", "cpp"),
    TestDefinition("COEF3D04", "COEF3D", "Complete Ito drift", "cpp"),
    TestDefinition("COEF3D05", "COEF3D", "Invalid coefficient status", "cpp"),
    TestDefinition("COEF3D06", "COEF3D", "Parallel coefficient gradient", "cpp"),
    TestDefinition("COEF3D07", "COEF3D", "Selectable coefficient models", "cpp"),
    TestDefinition("PRK3D01", "PRK3D", "Parallel diffusion moments", "cpp"),
    TestDefinition("PRK3D02", "PRK3D", "Advection", "cpp"),
    TestDefinition("PRK3D03", "PRK3D", "Orientation invariance", "cpp"),
    TestDefinition("PRK3D04", "PRK3D", "Nonuniform equilibrium", "cpp"),
    TestDefinition("PRK3D05", "PRK3D", "Adiabatic cooling", "cpp"),
    TestDefinition("PRK3D06", "PRK3D", "First passage", "cpp"),
    TestDefinition("PRK3D07", "PRK3D", "Radial PDE comparison", "cpp"),
    TestDefinition("PRK3D08", "PRK3D", "Named substep limits", "cpp"),
    TestDefinition("FTE3D01", "FTE3D", "Ballistic streaming", "cpp"),
    TestDefinition("FTE3D02", "FTE3D", "Magnetic focusing", "cpp"),
    TestDefinition("FTE3D03", "FTE3D", "Pitch-angle eigenmodes", "cpp"),
    TestDefinition("FTE3D04", "FTE3D", "Pitch boundaries", "cpp"),
    TestDefinition("FTE3D05", "FTE3D", "Momentum characteristic", "cpp"),
    TestDefinition("FTE3D06", "FTE3D", "Strong-scattering reduction", "cpp"),
    TestDefinition("FTE3D07", "FTE3D", "Zero perpendicular identity", "cpp"),
    TestDefinition("FTE3D08", "FTE3D", "Event-driven mean-free-path scattering", "cpp"),
    TestDefinition("FTE3D09", "FTE3D", "Scattering-frame energy", "cpp"),
    TestDefinition("RNG3D01", "RNG3D", "Thread reproducibility", "cpp"),
    TestDefinition("RNG3D02", "RNG3D", "Order independence", "cpp"),
    TestDefinition("RNG3D03", "RNG3D", "Purpose isolation", "cpp"),
    TestDefinition("V1D01", "V1D", "Gyrotropic tensor and Ito drift", "cpp"),
    TestDefinition("V1D02", "V1D", "Perpendicular diffusion moments", "cpp"),
    TestDefinition("V1D03", "V1D", "Guiding-centre drift direction", "cpp"),
    TestDefinition("V1D04", "V1D", "Focused perpendicular transport", "cpp"),
    TestDefinition("V1D05", "V1D", "Tensor diffusion timestep", "cpp"),
    TestDefinition("POP3D01", "POP3D", "Relativistic split/merge conservation", "cpp"),
    TestDefinition("V2D01", "V2D", "Distinct 1-D/3-D production-core parity", "source"),
    TestDefinition("V5D01", "V5D", "Validation/release governance contracts", "source"),
    TestDefinition("ADP3D01", "ADP3D", "Production mover dispatch", "cpp"),
    TestDefinition("NAT3D04", "NAT3D", "Boundary dispositions", "cpp"),
    TestDefinition("NAT3D05", "NAT3D", "Particle ledger closure", "cpp"),
    TestDefinition("NAT3D06", "NAT3D", "Sampling isolation", "cpp"),
    TestDefinition("NAT3D07", "NAT3D", "Output schema and publication", "cpp"),
    TestDefinition("NAT3D08", "NAT3D", "Shock crossing dispatch", "cpp"),
    TestDefinition("SHK3D01", "SHK3D", "Dimensional source identity", "cpp"),
    TestDefinition("SHK3D02", "SHK3D", "Expanding shock geometry", "cpp"),
    TestDefinition("SHK3D03", "SHK3D", "Source ownership guards", "cpp"),
    TestDefinition("SHK3D04", "SHK3D", "Source weight normalization", "cpp"),
    TestDefinition("RST3D01", "RST3D", "Complete restart round trip", "cpp"),
    TestDefinition("RST3D02", "RST3D", "Transactional restart rejection", "cpp"),
    TestDefinition("RST3D03", "RST3D", "Snapshot restart policy", "cpp"),
    TestDefinition("INT3D01", "INT3D", "Deterministic rank gather", "cpp"),
    TestDefinition("INT3D02", "INT3D", "Global conservation audit", "cpp"),
    TestDefinition("INT3D03", "INT3D", "Load and resource budgets", "cpp"),
    TestDefinition("VFY3D01", "VFY3D", "Profile metric normalization", "cpp"),
    TestDefinition("VFY3D02", "VFY3D", "Coverage and malformed evidence", "cpp"),
    TestDefinition("VFY3D03", "VFY3D", "Parker Green-function validation", "cpp"),
    TestDefinition("VFY3D04", "VFY3D", "Focused convergence", "cpp"),
    TestDefinition("VFY3D05", "VFY3D", "SWCME DSA distribution", "cpp"),
    TestDefinition("VFY3D06", "VFY3D", "Fixed phase-space distribution", "cpp"),
    TestDefinition("UTIL02", "UTIL", "Shared-kernel frozen record", "cpp"),
    TestDefinition("HARN02-EXITCODE", "HARN_SHELL", "Outer failure exit code", "shell", False),
    TestDefinition("HARN03-EXITCODE", "HARN_SHELL", "Outer skip exit code", "shell", False),
    TestDefinition("RUN3D01", "RUNNER", "Python runner CLI contract", "source"),
    TestDefinition("VALRUN3D01", "RUNNER", "Phase-V evidence runner contract", "source"),
    TestDefinition("BLDL3D01", "BLDL3D", "Configured enclosing AMPS build", "source"),
    TestDefinition("BLDL3D02", "BLDL3D", "Retired production-symbol exclusion", "source"),
    TestDefinition("BLDL3D03", "BLDL3D", "AMPS mover return-code mapping", "source"),
    TestDefinition("BLDL3D04", "BLDL3D", "AMPS Pi-macro namespace hygiene", "source"),
    TestDefinition("BLDL3D05", "BLDL3D", "Source/build makefile path resolution", "source"),
    TestDefinition("BLDL3D06", "BLDL3D", "Production transitive-header boundary", "source"),
    TestDefinition("BLDL3D07", "BLDL3D", "Application-object ABI freshness", "source"),
    TestDefinition("BLDL3D08", "BLDL3D", "Initialized native background output ordering", "source"),
    TestDefinition("BLDL3D09", "BLDL3D", "Active-region and population-control wiring", "source"),
    TestDefinition("BLDL3D10", "BLDL3D", "Solar internal-boundary wiring", "source"),
    TestDefinition("BLDL3D11", "BLDL3D", "Coronal-CME native-test wiring", "source"),
    TestDefinition("BLDL3D12", "BLDL3D", "Permanent PIC source/build header guards", "source"),
    TestDefinition("ARCH3D02", "ARCH3D", "Canonical shared-archive ownership", "source"),
    TestDefinition("SWCME3D01", "SWCME3D", "Relocated SWCME common runner", "source"),
    # Linked and external-evidence cases are intentionally non-routine.  They
    # remain visible through --all/--suite phase-v and report SKIP when their
    # configured executable or reviewed bundle is unavailable.
    TestDefinition("NAT3D01", "NAT3D", "AMPS mesh integration", "validation", False),
    TestDefinition("NAT3D02", "NAT3D", "Cell background integration", "validation", False),
    TestDefinition("NAT3D03", "NAT3D", "AMR gradient integration", "validation", False),
    TestDefinition("NAT3D09", "NAT3D", "Particle load balance", "validation", False),
    TestDefinition("NAT3D10", "NAT3D", "Production resource budgets", "validation", False),
    TestDefinition("NAT3D11", "NAT3D", "Coupled snapshot/shock schedule", "validation", False),
    TestDefinition("NAT3D12", "NAT3D", "Production product grammar", "validation", False),
    TestDefinition("MPI3D01", "MPI3D", "Multi-rank sampling reproducibility", "validation", False),
    TestDefinition("MPI3D02", "MPI3D", "Multi-rank restart continuation", "validation", False),
    TestDefinition("SCCM3D01", "SCCM3D", "Coronal-CME initialization ledger", "validation", False),
    TestDefinition("SCCM3D02", "SCCM3D", "Coronal-CME species numerics", "validation", False),
    TestDefinition("SCCM3D03", "SCCM3D", "Coronal-CME source species binding", "validation", False),
    TestDefinition("SCCM3D04", "SCCM3D", "Coronal-CME mesh boundary", "validation", False),
    TestDefinition("SCCM3D05", "SCCM3D", "Coronal-CME provider generations", "validation", False),
    TestDefinition("SCCM3D06", "SCCM3D", "Coronal-CME initialization output", "validation", False),
    TestDefinition("SCCM3D07", "SCCM3D", "Coronal-CME MPI identity", "validation", False),
    TestDefinition("XM3D01", "XM3D", "Parker cross-model profiles", "validation", False),
    TestDefinition("XM3D02", "XM3D", "Focused cross-model profiles", "validation", False),
    TestDefinition("XM3D03", "XM3D", "Longitudinal displacement diagnostic", "validation", False),
    TestDefinition("XM3D04", "XM3D", "SWCME shock-source reduction", "validation", False),
    TestDefinition("XM3D05", "XM3D", "Resolution convergence", "validation", False),
    TestDefinition("XM3D06", "XM3D", "Independent 3-D PDE reference", "validation", False),
    TestDefinition("OV3D01", "OV3D", "2013 April 11 near-Earth event", "validation", False),
    TestDefinition("OV3D02", "OV3D", "2020 May 29 PSP/STEREO-A event", "validation", False),
    TestDefinition("OV3D03", "OV3D", "2014 January 6 PAMELA diagnostic", "validation", False),
    TestDefinition("OV3D04", "OV3D", "Electron multi-spacecraft diagnostic", "validation", False),
    TestDefinition("SWMF3D01", "SWMF3D", "Live SWMF coupling replay (blocked by R8)", "validation", False),
)

BY_ID: Dict[str, TestDefinition] = {item.test_id: item for item in TESTS}
GROUPS: Dict[str, List[TestDefinition]] = {}
for _item in TESTS:
    GROUPS.setdefault(_item.group, []).append(_item)

SUITES: Dict[str, Tuple[str, ...]] = {
    "standalone": tuple(item.test_id for item in TESTS
                        if item.kind in ("cpp", "shell") or
                        item.test_id in ("RUN3D01", "VALRUN3D01")),
    "r0": ("BLDL3D01", "BLDL3D02", "BLDL3D03", "BLDL3D04", "BLDL3D05",
           "BLDL3D06", "BLDL3D07", "BLDL3D08", "BLDL3D09", "BLDL3D10",
           "BLDL3D11",
           "RUN3D01", "LAY01", "BLD01"),
    "r1": ("ARCH3D02", "SWCME3D01", "UTIL02"),
    "r2": ("LIFE3D01", "LIFE3D02", "LIFE3D03", "LIFE3D04"),
    "improvements-r": tuple(item.test_id for item in TESTS
                             if item.group == "R3D"),
    "improvements-c": tuple(item.test_id for item in TESTS
                            if item.group == "CFG3D"),
    "improvements-v": tuple(item.test_id for item in TESTS
                            if item.group in ("V1D", "V2D", "V5D")),
    "phase-m": tuple(item.test_id for item in TESTS if item.group in ("MSH3D", "DOM3D")),
    "phase-b": tuple(item.test_id for item in TESTS
                     if item.group in ("BGP3D", "SNAP3D")),
    "phase-t": tuple(item.test_id for item in TESTS
                     if item.group == "TUR3D" or
                     item.test_id in ("COEF3D01", "COEF3D02", "COEF3D06",
                                      "COEF3D07")),
    "phase-p": tuple(item.test_id for item in TESTS
                     if item.group in ("PRK3D", "FTE3D", "RNG3D") or
                     item.test_id in ("COEF3D03", "COEF3D04", "COEF3D05",
                                      "POP3D01")),
    "phase-a": tuple(item.test_id for item in TESTS
                     if item.group in ("ADP3D", "SHK3D") or
                     item.test_id in ("NAT3D04", "NAT3D05", "NAT3D08")),
    "phase-o": tuple(item.test_id for item in TESTS
                     if item.group in ("RST3D", "OUT3D") or
                     item.test_id in ("NAT3D06", "NAT3D07")),
    "phase-v": tuple(item.test_id for item in TESTS
                     if item.group in ("V1D", "V2D", "V5D", "INT3D", "VFY3D", "MPI3D", "SCCM3D", "XM3D", "OV3D", "SWMF3D") or
                     (item.group == "NAT3D" and item.kind == "validation") or
                     item.test_id == "VALRUN3D01"),
    "production": ("BLDL3D01", "BLDL3D02", "BLDL3D03", "BLDL3D04",
                   "BLDL3D05", "BLDL3D06", "BLDL3D07", "BLDL3D08",
                   "BLDL3D09", "BLDL3D10", "BLDL3D11", "BLDL3D12"),
}


class RunnerError(RuntimeError):
    """Usage/setup error that must exit 2 rather than look like physics FAIL."""


class _HelpFormatter(argparse.ArgumentDefaultsHelpFormatter,
                     argparse.RawDescriptionHelpFormatter):
    pass


EPILOG = """
Examples:
  python3 test/run_tests.py --list
  python3 test/run_tests.py --routine --amps-source /path/to/AMPS
  python3 test/run_tests.py --test BLDL3D02 --output-dir test_output/r0-source
  python3 test/run_tests.py --group HARN --group BLDL3D
  python3 test/run_tests.py --suite standalone --suite production
  python3 test/run_tests.py --all --amps-source /path/to/AMPS

For a configured AMPS checkout, either place srcSEP3D in its normal
application location or provide --make-config /path/to/Makefile.conf.  A
source-only archive can run the standalone, R1, R2, M/B/T/P/A/O/V, and source-only
build tests. BLDL3D01 is recorded as SKIP until a real production configuration is
present; BLDL3D03 is SKIP if the actual AMPS pic.h is unavailable.
"""


def _utc_stamp() -> str:
    return _datetime.datetime.now(_datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run srcSEP3D standalone, R0-R2/C01-C05/M/B/T/P/A/O/V, and production-build tests.",
        formatter_class=_HelpFormatter,
        epilog=EPILOG)
    parser.add_argument("--amps", default=os.environ.get("SEP3D_EXECUTABLE"),
                        help="linked srcSEP3D/AMPS executable for Phase-V native cases")
    parser.add_argument("--validation-data", type=Path,
                        default=os.environ.get("SEP3D_VALIDATION_DATA"),
                        help="root containing CASE_ID/manifest.json evidence bundles")
    parser.add_argument(
        "--validation-input", type=Path,
        default=os.environ.get("SEP3D_VALIDATION_INPUT"),
        help="complete immutable input deck for Phase-V linked cases")
    parser.add_argument("--validation-launch-prefix", default="",
                        help="optional argv prefix for linked cases, e.g. 'mpiexec -n 8'")
    parser.add_argument("--amps-source", type=Path,
                        default=os.environ.get("AMPS_SOURCE_ROOT"),
                        help="AMPS root containing src/pic/pic.h for ABI checks")
    parser.add_argument("--sep1d-source", type=Path,
                        help="explicit srcSEP root for V2D01; never inferred")
    parser.add_argument("--make-config", type=Path,
                        default=os.environ.get("AMPS_MAKE_CONFIG"),
                        help="configured AMPS Makefile.conf for BLDL3D01")
    parser.add_argument("--list", action="store_true",
                        help="list registered tests and suites, then exit")
    parser.add_argument("--test", dest="tests", action="append", default=[],
                        metavar="ID", help="run one stable test ID; repeatable")
    parser.add_argument("--group", dest="groups", action="append", default=[],
                        metavar="GROUP", help="run one test group; repeatable")
    parser.add_argument("--routine", action="store_true",
                        help="run the bounded routine set")
    parser.add_argument("--all", action="store_true",
                        help="run every public test separately and continue after failures")
    parser.add_argument("--suite", dest="suites", action="append", default=[],
                        choices=sorted(SUITES), help="run a named suite; repeatable")
    parser.add_argument("--output-dir", type=Path,
                        default=Path("test_output") / f"srcsep3d-tests-{_utc_stamp()}",
                        help="directory for logs and JSON/JUnit summaries")
    parser.add_argument("--sep-common-dir", type=Path,
                        help="canonical sep_common source/header directory")
    parser.add_argument("--sep-common-archive", type=Path,
                        help="path to sep_common.a")
    parser.add_argument("--cxx", default="g++", help="standalone C++ compiler")
    parser.add_argument("--rebuild", action="store_true",
                        help="force rebuilding test/stage1")
    parser.add_argument("--no-build", action="store_true",
                        help="use an existing test/stage1 binary")
    parser.add_argument("--timeout", type=float, default=300.0,
                        help="per-test command timeout in seconds")
    parser.add_argument("--verbose", action="store_true",
                        help="print complete commands and subprocess output")
    parser.add_argument("--no-color", action="store_true",
                        help="accepted for CLI parity; output is already plain text")
    return parser


def _resolve_existing(explicit: Optional[Path], candidates: Iterable[Path],
                      description: str) -> Path:
    if explicit is not None:
        path = explicit.expanduser().resolve()
        if not path.exists():
            raise RunnerError(f"{description} does not exist: {path}")
        return path
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    rendered = "\n  ".join(str(item) for item in candidates)
    raise RunnerError(f"cannot locate {description}; checked:\n  {rendered}")


def _sep_common_paths(args: argparse.Namespace) -> Tuple[Path, Path]:
    common_dir = _resolve_existing(
        args.sep_common_dir,
        (ROOT.parent / "src" / "models" / "sep_common",
         ROOT.parent / "sep_common",
         ROOT.parent.parent / "sep_common"),
        "sep_common directory")

    archive = (args.sep_common_archive.expanduser().resolve()
               if args.sep_common_archive is not None
               else common_dir / "sep_common.a")
    if not archive.is_file():
        make = shutil.which("make")
        if make is None:
            raise RunnerError(
                f"sep_common archive is absent and make is unavailable: {archive}")
        completed = subprocess.run(
            [make, "-C", str(common_dir), "all"], text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
        if completed.returncode != 0:
            raise RunnerError("sep_common build failed:\n" + completed.stdout)
    if not archive.is_file():
        raise RunnerError(f"sep_common build produced no archive: {archive}")
    return common_dir, archive


def _build_standalone(args: argparse.Namespace) -> None:
    if args.no_build:
        if not BINARY.is_file():
            raise RunnerError("--no-build was requested but test/stage1 is absent")
        return

    source, archive = _sep_common_paths(args)
    compiler = shutil.which(args.cxx)
    if compiler is None:
        raise RunnerError(f"C++ compiler not found: {args.cxx}")
    make = shutil.which("make")
    if make is None:
        raise RunnerError("make is required to build the standalone registry")

    # Use the application makefile instead of one monolithic compiler command.
    # This keeps the runner and documented build manifest identical and lets
    # GNU Make honor inherited MAKEFLAGS=-jN for independent object files.
    # Resolve SWCME from the explicitly supplied AMPS root as well.  Relying on
    # the process environment here makes an otherwise detached srcSEP3D tree
    # accidentally search its incomplete sibling layout during --all.
    _, swcme = _canonical_model_dirs(args)
    variables = [f"CXX={compiler}", f"SEP_COMMON_DIR={source}",
                 f"SEP_COMMON_ARCHIVE={archive}", f"SWCME_DIR={swcme}"]
    commands: List[List[str]] = []
    if args.rebuild:
        commands.append([make, "-C", str(ROOT), "clean-standalone", *variables])
    commands.append([make, "-C", str(ROOT), "test/stage1", *variables])
    for command in commands:
        if args.verbose:
            print("RUN:", shlex.join(command))
        completed = subprocess.run(
            command, cwd=ROOT, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, timeout=args.timeout, check=False)
        if completed.returncode != 0:
            raise RunnerError("standalone build failed:\n" + completed.stdout)
        if args.verbose and completed.stdout:
            print(completed.stdout, end="")


def _run_command(command: List[str], cwd: Path, timeout: float,
                 verbose: bool) -> Tuple[int, str, float]:
    if verbose:
        print("RUN:", shlex.join(command))
    started = time.monotonic()
    try:
        completed = subprocess.run(command, cwd=cwd, text=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT,
                                   timeout=timeout, check=False)
        elapsed = time.monotonic() - started
        if verbose and completed.stdout:
            print(completed.stdout, end="")
        return completed.returncode, completed.stdout, elapsed
    except subprocess.TimeoutExpired as exc:
        output = (exc.stdout or "") + (exc.stderr or "")
        return 2, output + f"\nTimed out after {timeout:g} s", time.monotonic() - started


def _parse_cpp_result(report: Path, definition: TestDefinition,
                      command: List[str], elapsed: float, output: str) -> Result:
    if not report.is_file():
        return Result(definition.test_id, definition.group, "ERROR",
                      "C++ test wrote no JSON report. Output: " + output[-2000:],
                      elapsed, command)
    try:
        payload = json.loads(report.read_text(encoding="utf-8"))
        records = payload.get("results", [])
        record = next(item for item in records
                      if str(item.get("id", "")).upper() == definition.test_id)
    except (OSError, ValueError, StopIteration) as exc:
        return Result(definition.test_id, definition.group, "ERROR",
                      f"cannot parse C++ report: {exc}", elapsed, command)
    return Result(definition.test_id, definition.group,
                  str(record.get("status", "ERROR")).upper(),
                  str(record.get("message", "")), elapsed, command)


def _run_cpp(definition: TestDefinition, args: argparse.Namespace,
             output_dir: Path) -> Result:
    report = output_dir / f"{definition.test_id}.json"
    command = [str(BINARY), "--test", definition.test_id,
               "--test-json", str(report)]
    code, output, elapsed = _run_command(command, ROOT, args.timeout, args.verbose)
    result = _parse_cpp_result(report, definition, command, elapsed, output)
    if code == 2 and result.status not in ("ERROR",):
        result.status = "ERROR"
        result.message += " (binary exit code 2)"
    if definition.test_id == "LIFE3D04" and result.status == "PASS":
        # The C++ half proves both adapters traverse the same Runtime calls.
        # Audit the actual coupled boundary as well: a fixture alone cannot
        # prove main_lib.cpp has not gained an independent parser.
        coupled = _strip_cpp_comments(
            (ROOT / "main_lib.cpp").read_text(encoding="utf-8"))
        forbidden = ("argc", "argv", "AMPS_PARAM.in", "std::ifstream",
                     "getenv(", "ParseCommandLine")
        observed = [token for token in forbidden if token in coupled]
        if observed:
            result.status = "FAIL"
            result.message = (
                "coupled entry points contain internal configuration parsing: " +
                ", ".join(observed))
        else:
            result.message += (
                "; production coupled entry points contain no process-argument, "
                "AMPS parameter-file, stream, or environment parser")
    return result


def _run_shell(definition: TestDefinition, args: argparse.Namespace) -> Result:
    script_name = {
        "HARN02-EXITCODE": "run_harn02.sh",
        "HARN03-EXITCODE": "run_harn03.sh",
    }[definition.test_id]
    command = ["sh", str(ROOT / "test" / "individual-test" / script_name)]
    code, output, elapsed = _run_command(command, ROOT, args.timeout, args.verbose)
    return Result(definition.test_id, definition.group,
                  "PASS" if code == 0 else "FAIL",
                  output.strip() or f"script exited {code}", elapsed, command)


def _production_files() -> List[Path]:
    files = [ROOT / "SEP3D.h", ROOT / "main_lib.cpp", ROOT / "main.cpp"]
    for directory in (ROOT / "core", ROOT / "background", ROOT / "runtime",
                      ROOT / "mesh", ROOT / "turbulence", ROOT / "transport",
                      ROOT / "adapters", ROOT / "output", ROOT / "validation",
                      ROOT / "amps"):
        files.extend(sorted(directory.glob("*.h")))
        files.extend(sorted(directory.glob("*.cpp")))
    return files


def _strip_cpp_comments(text: str) -> str:
    """Remove comments before checking for retired *identifiers*.

    Migration comments are allowed to explain a removed name.  The production
    gate is concerned with compiled tokens, so both // and /*...*/ comments
    are removed while quoted strings are retained (a runtime lookup by a
    retired name must still fail the check).
    """
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
    return re.sub(r"//.*?$", "", text, flags=re.MULTILINE)


def _check_retired_sources(definition: TestDefinition) -> Result:
    started = time.monotonic()
    errors: List[str] = []
    retired_file = ROOT / "SEP3D.cpp"
    if retired_file.exists():
        errors.append(
            "retired SEP3D.cpp is still present; remove this exact stale file "
            "from the installed srcSEP3D tree before rebuilding")

    makefile = (ROOT / "makefile").read_text(encoding="utf-8")
    # Strip comments before checking the active production manifest.
    active_make = "\n".join(line.split("#", 1)[0] for line in makefile.splitlines())
    if "SEP3D.cpp" in active_make or "SEP3D.o" in active_make:
        errors.append("makefile production manifest still references SEP3D.cpp")

    forbidden = (
        "Mover_Axisymmetric_SecondOrder", "TotalParticleAcceleration",
        "GlobalEnergyDistribution", "inject_particle_onto_field_line",
        "8.760e+08", "9.445e+08",
    )
    for path in _production_files():
        text = _strip_cpp_comments(path.read_text(encoding="utf-8"))
        for token in forbidden:
            if token in text:
                errors.append(f"{path.relative_to(ROOT)} contains retired token {token}")

    main_lib = (ROOT / "main_lib.cpp").read_text(encoding="utf-8")
    # OutputDistributedDataTECPLOT is the cut-cell-aware native initialization
    # evidence path. Prepopulation and binary mesh saves remain retired
    # prototype operations and must not return to the production driver.
    for call in ("PrepopulateDomain", "saveMeshFile"):
        if call in main_lib:
            errors.append(f"main_lib.cpp retains prototype operation {call}")

    status = "FAIL" if errors else "PASS"
    message = "; ".join(errors) if errors else (
        "retired mover/sampler source is absent; application sources are explicit "
        "and M/B/T/P/A/O/V implementations remain in layered modules; wedge and "
        "prepopulation operations are absent")
    return Result(definition.test_id, definition.group, status, message,
                  time.monotonic() - started, [])


def _check_initialized_native_background(definition: TestDefinition) -> Result:
    """Guard the data-bearing initialization product's ordering contract.

    A linked AMPS qualification run remains the numerical authority, but this
    source gate catches both output-bridge regressions seen in production:
    filling only srcSEP3D's private cache left native Bx/By/Bz at zero, while
    omitting an InterpolateCenterNode callback left application-owned
    turbulence at zero in AMPS' temporary vertex nodes.
    """
    started = time.monotonic()
    main_lib = (ROOT / "main_lib.cpp").read_text(encoding="utf-8")
    output_sampling = (ROOT / "output" / "sampling.cpp").read_text(
        encoding="utf-8")
    run_configuration = (ROOT / "runtime" / "run_configuration.cpp").read_text(
        encoding="utf-8")
    required = (
        "PIC::CPLR::DATAFILE::CenterNodeAssociatedDataOffsetBegin",
        "Offset::PlasmaNumberDensity",
        "Offset::PlasmaBulkVelocity",
        "Offset::PlasmaTemperature",
        "Offset::PlasmaIonPressure",
        "Offset::MagneticField",
        "Offset::ElectricField",
        "Offset::MagneticFieldGradient",
        "ZeroApplicationStateOnOwnedCells();",
        "ZeroNativeAmpsBackgroundOnOwnedCells();",
        "StoreNativeAmpsBackground(cells[i].cell, snapshot->samples()[i]);",
        "StoreTurbulenceAtCellCenter(cells[i].cell, waves)",
        "CompleteNativeAmpsBackgroundInstallation();",
        "PIC::Mesh::InterpolateCenterNode.push_back(",
        "InterpolateInitializationCellData",
        "WriteInitializationDataTecplotAfterBackground();",
        "OutputDistributedDataTECPLOT(",
        "radiusM >= innerRadiusM && radiusM <= outerRadiusM",
        "TurbulenceTecplotVariableList()",
        "PrepareTurbulenceTecplotPresentation(",
        "mandatory directional turbulence storage is absent",
    )
    missing = [token for token in required if token not in main_lib]
    required_output = (
        "turbulence_wave_energy_density_J_per_m3",
        "turbulence_wave_energy_plus_J_per_m3",
        "turbulence_wave_energy_minus_J_per_m3",
        "candidate.waveEnergyJPerM3 = candidate.deltaB2T2 /",
        "Core::Status InterpolateStaticCenterState(",
    )
    missing.extend(token for token in required_output
                   if token not in output_sampling)
    if "AppendField(2, &cursor, &layout.waveEnergyOffset);" not in \
            run_configuration:
        missing.append("unconditional wave-energy storage allocation")
    if missing:
        return Result(
            definition.test_id, definition.group, "FAIL",
            "native initialization bridge is incomplete: " +
            "; ".join(missing), time.monotonic() - started, [])

    try:
        fill = main_lib[main_lib.index("void FillAndPublishBackground()"):
                        main_lib.index("void WriteInitializationDataTecplotAfterBackground()")]
        writer = main_lib[
            main_lib.index("void WriteInitializationDataTecplotAfterBackground()"):
            main_lib.index("void RefreshBackgroundAtBoundary()")]
        initialize = main_lib[main_lib.index("void amps_init()"):
                              main_lib.index("int amps_time_step()")]
        fill_order = [
            fill.index("ZeroApplicationStateOnOwnedCells();"),
            fill.index("ZeroNativeAmpsBackgroundOnOwnedCells();"),
            fill.index("StoreNativeAmpsBackground("),
            fill.index("StoreTurbulenceAtCellCenter("),
            fill.rindex("CompleteNativeAmpsBackgroundInstallation();"),
        ]
        initialize_order = [
            initialize.index("FillAndPublishBackground();"),
            initialize.index("WriteInitializationDataTecplotAfterBackground();"),
        ]
    except ValueError as error:
        return Result(
            definition.test_id, definition.group, "FAIL",
            f"cannot identify native initialization ordering boundary: {error}",
            time.monotonic() - started, [])

    if fill_order != sorted(fill_order) or initialize_order != sorted(initialize_order):
        return Result(
            definition.test_id, definition.group, "FAIL",
            "application/native zero/fill/turbulence/halo sequence or final writer call is out of order",
            time.monotonic() - started, [])
    if ("gNativeAmpsBackgroundReady" not in writer or
            not re.search(r"PIC::Mesh::mesh->OutputDistributedDataTECPLOT\(\s*"
                          r"path\.c_str\(\),\s*true,\s*species\.ampsIndex\s*\)", writer) or
            "SetAssembleDistributedOutputFileFlag(true);" not in writer):
        return Result(
            definition.test_id, definition.group, "FAIL",
            "initialization writer must assemble cut-cell data for every species "
            "after the completed native-background boundary",
            time.monotonic() - started, [])

    return Result(
        definition.test_id, definition.group, "PASS",
        "validated background/turbulence are stored at physical center nodes, "
        "the application slice participates in AMPS vertex interpolation, wave "
        "energy is emitted as mandatory total/directional SI columns, halos are "
        "exchanged, and completion is checked before assembled cut-cell initialization output",
        time.monotonic() - started, [])


def _check_active_population_wiring(definition: TestDefinition) -> Result:
    """Guard AMPS lifecycle placement and split/merge safety invariants.

    This source test complements the portable geometry and conservation unit
    tests.  It deliberately checks the production driver because neither an
    isolated classifier test nor an isolated relativistic merge test can prove
    that AMPS consumes those implementations at a safe, globally joined point.
    A configured BLDL3D01 run remains the authority for the concrete AMPS ABI.
    """
    started = time.monotonic()
    main_lib = (ROOT / "main_lib.cpp").read_text(encoding="utf-8")
    mesh_model = (ROOT / "mesh" / "mesh_model.cpp").read_text(
        encoding="utf-8")
    adapter = (ROOT / "amps" / "amps_particle_adapter.cpp").read_text(
        encoding="utf-8")
    population = (ROOT / "transport" / "population_control.cpp").read_text(
        encoding="utf-8")
    core_split_path = (ROOT.parent / "src" / "pic" /
                       "pic_particle_spliting.cpp")
    if not core_split_path.is_file():
        return Result(
            definition.test_id, definition.group, "FAIL",
            "required AMPS core source is absent from the tested tree: "
            "src/pic/pic_particle_spliting.cpp",
            time.monotonic() - started, [])
    core_split = core_split_path.read_text(encoding="utf-8")

    required_by_file = {
        "main_lib.cpp": (
            "BuildActiveRegionPlan(",
            "SetTreeNodeActiveUseFlag(",
            "GetNeibFace(",
            "VerifyActiveRegionAllocation();",
            "PIC::ParticleSplitting::SetMode(PIC::ParticleSplitting::_disactivated);",
            "node != nullptr && node->IsUsedInCalculationFlag",
            "ApplyPopulationControlAtBoundary();",
        ),
        "mesh/mesh_model.cpp": (
            "double ActiveTubeRadiusM(",
            "bool BlockIntersectsActiveRegion(",
        ),
        "amps/amps_particle_adapter.cpp": (
            "MergeRelativisticThreeToTwo(",
            "NewPopulationStableId(",
            "ResolvePopulationMagneticDirection(",
            "PIC::ParticleBuffer::CloneParticle(",
        ),
        "transport/population_control.cpp": (
            "RelativisticTotalEnergyJ(",
            "MergeRelativisticThreeToTwo(",
        ),
    }
    texts = {
        "main_lib.cpp": main_lib,
        "mesh/mesh_model.cpp": mesh_model,
        "amps/amps_particle_adapter.cpp": adapter,
        "transport/population_control.cpp": population,
    }
    missing: List[str] = []
    for name, tokens in required_by_file.items():
        missing.extend(f"{name}: {token}" for token in tokens
                       if token not in texts[name])

    # These checks are intentionally scoped to the two public functions and
    # validate ordering as well as spelling.  A matching token in a comment,
    # or a guard placed after reserve(), particles[0], or std::next(end()),
    # would not repair the undefined-behaviour path that this gate protects.
    try:
        merge_start = core_split.index(
            "void PIC::ParticleSplitting::MergeParticleList(")
        split_start = core_split.index(
            "void PIC::ParticleSplitting::SplitParticleList(", merge_start)
        merge_body = core_split[merge_start:split_start]
        split_body = core_split[split_start:]

        merge_target = merge_body.index(
            "if (nRequestedParticleNumber < 1) return;")
        merge_traversal = merge_body.index("while (particleIndex != -1)")
        merge_empty = merge_body.index(
            "if (particles.empty() || particles.size() <=")
        merge_first_record = merge_body.index(
            "double vMin[3] = {particles[0]")
        merge_tuple = merge_body.index("if (particles.size() < 3) return;")
        merge_empty_bin = merge_body.index("if (binList.empty())")
        merge_successor = merge_body.index("if (std::next(currentBinIt)")

        split_target = split_body.index(
            "if (nRequestedParticleNumber < 1) return;")
        split_reserve = split_body.index("particles.reserve(")
        split_tuple = split_body.index("if (particles.size() < 2) return;")
        split_first_record = split_body.index(
            "double vMin[3] = {particles[0]")
        split_empty_bin = split_body.index("if (binList.empty())")
        split_successor = split_body.index("if (std::next(currentBinIt)")
    except ValueError as error:
        missing.append(
            "pic_particle_spliting.cpp: incomplete guarded entry points "
            f"({error})")
    else:
        if merge_target > merge_traversal or split_target > split_reserve:
            missing.append(
                "pic_particle_spliting.cpp: non-positive target guards must "
                "precede list traversal/reserve")
        if merge_empty > merge_first_record:
            missing.append(
                "pic_particle_spliting.cpp: empty/no-op merge guard must "
                "precede particles[0]")
        if merge_tuple > merge_successor or merge_empty_bin > merge_successor:
            missing.append(
                "pic_particle_spliting.cpp: merge tuple/bin guards must "
                "precede iterator successor access")
        if split_tuple > split_first_record or split_empty_bin > split_successor:
            missing.append(
                "pic_particle_spliting.cpp: split tuple/bin guards must "
                "precede record/iterator access")

    # The core fix must travel with the SEP overlay.  Previously the README
    # and this test described the guarded implementation while the release
    # allowlist silently omitted the actual .cpp, leaving an old destination
    # copy active after installation.
    package_policy = ROOT.parent / "tools" / "sep_package_hygiene.py"
    required_release_members: Tuple[str, ...] = ()
    if package_policy.is_file():
        # Read the literal policy declaration instead of accepting the path
        # anywhere in the file.  A path mentioned only in a comment or in the
        # permissive allowlist does not prove that a missing member fails the
        # package audit (the defect this gate is intended to prevent).
        try:
            policy_tree = ast.parse(
                package_policy.read_text(encoding="utf-8"),
                filename=str(package_policy))
            for statement in policy_tree.body:
                target = None
                value = None
                if isinstance(statement, ast.AnnAssign):
                    target = statement.target
                    value = statement.value
                elif isinstance(statement, ast.Assign) and \
                     len(statement.targets) == 1:
                    target = statement.targets[0]
                    value = statement.value
                if isinstance(target, ast.Name) and \
                   target.id == "TOP_LEVEL_REQUIRED" and value is not None:
                    literal = ast.literal_eval(value)
                    if isinstance(literal, tuple) and \
                       all(isinstance(item, str) for item in literal):
                        required_release_members = literal
                    break
        except (OSError, SyntaxError, ValueError):
            required_release_members = ()
    if "src/pic/pic_particle_spliting.cpp" not in required_release_members:
        missing.append(
            "sep_package_hygiene.py: guarded AMPS core source is not a "
            "required release member")
    if missing:
        return Result(
            definition.test_id, definition.group, "FAIL",
            "active-region/population integration is incomplete: " +
            "; ".join(missing), time.monotonic() - started, [])

    try:
        init_mesh = main_lib[main_lib.index("void amps_init_mesh()"):
                             main_lib.index("void amps_init()")]
        init_order = [
            init_mesh.index("PIC::Mesh::mesh->buildMesh();"),
            init_mesh.index("ApplyActiveRegionMask(resolution);"),
            init_mesh.index("PIC::Mesh::mesh->SetParallelLoadMeasure("),
            init_mesh.index("PIC::Mesh::mesh->CreateNewParallelDistributionLists();"),
            init_mesh.index("PIC::Mesh::mesh->AllocateTreeBlocks();"),
            init_mesh.index("VerifyActiveRegionAllocation();"),
        ]
        time_step = main_lib[main_lib.index("int amps_time_step()") :]
        boundary_order = [
            time_step.index("SEP3D::AMPS::Movers::InjectParticles(plan);"),
            time_step.index("ApplyPopulationControlAtBoundary();"),
            time_step.index("PublishObserversAtBoundary();"),
            time_step.index("WriteCheckpointAtBoundary();"),
        ]
        population_function = adapter[
            adapter.index("PopulationControlReport ApplyPopulationControl("):
            adapter.index("int MoveParticle(")]
        species_loop = population_function.index(
            "for (int species = 0; species < PIC::nTotalSpecies; ++species)")
        semantic_ordinal = population_function.index(
            "std::uint64_t operationIndex = 0;")
        first_population_action = population_function.index(
            "const bool mergeTriggered")
    except ValueError as error:
        return Result(
            definition.test_id, definition.group, "FAIL",
            f"cannot identify production lifecycle boundary: {error}",
            time.monotonic() - started, [])

    if init_order != sorted(init_order):
        return Result(
            definition.test_id, definition.group, "FAIL",
            "active-use flags must be applied after AMR construction and before "
            "load measurement and distribution, with allocation verified "
            "immediately after block creation",
            time.monotonic() - started, [])
    if boundary_order != sorted(boundary_order):
        return Result(
            definition.test_id, definition.group, "FAIL",
            "population control must follow shock injection and precede observer "
            "publication and checkpointing",
            time.monotonic() - started, [])
    if not (species_loop < semantic_ordinal < first_population_action) or \
       population_function.count("std::uint64_t operationIndex = 0;") != 1:
        return Result(
            definition.test_id, definition.group, "FAIL",
            "resampling operation ordinals must restart inside each physical "
            "cell/species population so MPI block ownership cannot change "
            "post-resampling random keys",
            time.monotonic() - started, [])

    return Result(
        definition.test_id, definition.group, "PASS",
        "the finite Parker-corridor planner uses native AMR neighbours and "
        "AMPS active-use flags before distribution, then audits allocation; "
        "inactive source patches are rejected; the legacy automatic splitter is "
        "disabled; SEP-aware relativistic resampling runs after injection and "
        "before observers/checkpoints; resampling keys are independent of MPI "
        "block traversal; and generic AMPS empty/singleton guards are present",
        time.monotonic() - started, [])


def _check_solar_boundary_wiring(definition: TestDefinition) -> Result:
    """Guard the AMPS solar-sphere registration and cut-cell lifecycle.

    Portable MSH3D15 owns the geometry arithmetic.  This source gate checks
    the AMPS-only API and ordering that cannot be linked by the standalone
    test binary; BLDL3D01 remains the final concrete-ABI authority.
    """
    started = time.monotonic()
    main_lib = (ROOT / "main_lib.cpp").read_text(encoding="utf-8")
    mesh_header = (ROOT / "mesh" / "mesh_model.h").read_text(
        encoding="utf-8")
    mesh_model = (ROOT / "mesh" / "mesh_model.cpp").read_text(
        encoding="utf-8")
    run_configuration = (ROOT / "runtime" / "run_configuration.cpp").read_text(
        encoding="utf-8")

    required_main = (
        "#if _INTERNAL_BOUNDARY_MODE_ != _INTERNAL_BOUNDARY_MODE_ON_",
        "_USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ !=",
        "_USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ON_",
        "void RegisterSolarSurfaceBoundary()",
        "PIC::BC::InternalBoundary::Sphere::Init();",
        "PIC::BC::InternalBoundary::Sphere::RegisterInternalSphere();",
        "SetSphereGeometricalParameters(",
        "gSolarSurfaceBoundary->localResolution = localResolution;",
        "gSolarSurfaceBoundary->InjectionRate = nullptr;",
        "gSolarSurfaceBoundary->InjectionBoundaryCondition = nullptr;",
        "gSolarSurfaceBoundary->ParticleSphereInteraction =",
        "return _PARTICLE_DELETED_ON_THE_FACE_;",
        "MPI_Allreduce(localGeometry, minimumGeometry",
        "MPI_Allreduce(localGeometry, maximumGeometry",
        "AxisAlignedBoxEntirelyInsideSolarBoundary(",
        "cell->Measure = 0.0;",
        "-_GHOST_CELLS_X_",
        "-_GHOST_CELLS_Y_",
        "-_GHOST_CELLS_Z_",
    )
    required_mesh = (
        "struct SolarBoundaryGeometry",
        "SolarBoundaryGeometry MakeSolarBoundary(",
        "result.radiusM = Core::Const::R_sun;",
        "bool AxisAlignedBoxEntirelyInsideSolarBoundary(",
    )
    missing = [f"main_lib.cpp: {token}" for token in required_main
               if token not in main_lib]
    combined_mesh = mesh_header + "\n" + mesh_model
    missing.extend(f"mesh_model: {token}" for token in required_mesh
                   if token not in combined_mesh)
    if "normalized.innerRadiusM < Core::Const::R_sun" not in run_configuration:
        missing.append("run_configuration.cpp: source shell below R_sun guard")
    if "sep3d-physics-v8" not in run_configuration:
        missing.append("run_configuration.cpp: solar-boundary physics identity")
    for token in (";solar_boundary=amps-absorbing-sphere-v1",
                  ";solar_radius_m=", "Core::Const::R_sun"):
        if token not in run_configuration:
            missing.append(
                f"run_configuration.cpp: solar-boundary provenance {token}")
    singleton_calls = (
        "PIC::BC::InternalBoundary::Sphere::Init();",
        "PIC::BC::InternalBoundary::Sphere::RegisterInternalSphere();",
        "RegisterSolarSurfaceBoundary();",
    )
    for call in singleton_calls:
        if main_lib.count(call) != 1:
            missing.append(f"main_lib.cpp: exactly one call to {call}")
    if "PIC::Mesh::mesh->RegisterInternalBoundary" in main_lib:
        missing.append("main_lib.cpp: duplicate direct internal-boundary registration")

    try:
        callback = main_lib[
            main_lib.index("int AbsorbParticleAtSolarSurface("):
            main_lib.index("void RegisterSolarSurfaceBoundary()")]
        registration = main_lib[
            main_lib.index("void RegisterSolarSurfaceBoundary()"):
            main_lib.index("double InitLoadMeasure(")]
        init_mesh = main_lib[main_lib.index("void amps_init_mesh()"):
                             main_lib.index("void amps_init()")]
        correction = main_lib[
            main_lib.index("void CorrectSolarInteriorCellMeasures()"):
            main_lib.index("bool TrajectoryTrackingCondition(")]
        active_mask = main_lib[
            main_lib.index("void ApplyActiveRegionMask("):
            main_lib.index("void VerifyActiveRegionAllocation()")]
        geometry = mesh_model[
            mesh_model.index("SolarBoundaryGeometry MakeSolarBoundary("):
            mesh_model.index("double TubeRadiusM(")]
        call_site = init_mesh[
            init_mesh.index("SEP3D::Init_BeforeParser();"):
            init_mesh.index("PIC::Mesh::initCellSamplingDataBuffer();")]
        init_order = [
            init_mesh.index("PIC::Init_BeforeParser();"),
            init_mesh.index("SEP3D::Init_BeforeParser();"),
            init_mesh.index("RegisterSolarSurfaceBoundary();"),
            init_mesh.index("PIC::Mesh::initCellSamplingDataBuffer();"),
            init_mesh.index("PIC::Mesh::mesh->init("),
            init_mesh.index("PIC::Mesh::mesh->buildMesh();"),
            init_mesh.index("ApplyActiveRegionMask(resolution);"),
            init_mesh.index("PIC::Mesh::mesh->AllocateTreeBlocks();"),
            init_mesh.index("VerifyActiveRegionAllocation();"),
            init_mesh.index(
                "PIC::DomainBlockDecomposition::UpdateBlockTable();"),
            init_mesh.index("PIC::Mesh::mesh->InitCellMeasure();"),
            init_mesh.index("CorrectSolarInteriorCellMeasures();"),
        ]
    except ValueError as error:
        missing.append(f"production boundary/lifecycle slice: {error}")
    else:
        if "DeleteParticle" in callback:
            missing.append("absorbing callback must not delete the particle twice")
        if "PIC::ThisThread == 0" in registration.split(
                "PIC::BC::InternalBoundary::Sphere::Init();")[0]:
            missing.append("sphere registration is incorrectly rank-zero-only")
        if "PIC::ThisThread" in call_site or "if (" in call_site:
            missing.append(
                "sphere registration call site is conditionally rank-local")
        if init_order != sorted(init_order):
            missing.append("sphere/cut-cell initialization calls are out of order")
        if "innerRadiusM" in geometry:
            missing.append("physical solar radius is coupled to innerRadiusM")
        registration_authority = (
            "SEP3D::Mesh::MakeSolarBoundary(Configuration().options())")
        if registration_authority not in registration or \
                "centerM, geometry.radiusM" not in registration or \
                "geometry.centerM.x" not in registration:
            missing.append(
                "sphere registration bypasses authoritative center/radius")
        active_contract = (
            "SEP3D::Mesh::MakeSolarBoundary(Configuration().options())",
            "AxisAlignedBoxEntirelyInsideSolarBoundary(",
            "solarInterior[ordinal]",
            "plan.leafClass[ordinal] == "
            "SEP3D::Mesh::ActiveLeafClass::Inactive ||",
            "plan.leafClass[i] != SEP3D::Mesh::ActiveLeafClass::Inactive &&",
            "solarInterior[i] == 0",
        )
        for token in active_contract:
            if token not in active_mask:
                missing.append(
                    f"active mask omits solar/tube composition: {token}")
        if "AxisAlignedBoxEntirelyInsideSolarBoundary(" not in correction or \
                "cell->Measure = 0.0;" not in correction:
            missing.append("fully interior AMPS cell measures are not corrected")
        else:
            predicate = correction.index(
                "AxisAlignedBoxEntirelyInsideSolarBoundary(")
            preserve = correction.index("continue;", predicate)
            zero = correction.index("cell->Measure = 0.0;")
            if not predicate < preserve < zero:
                missing.append(
                    "fractional cut-cell measures are not preserved before "
                    "fully-interior zeroing")

    if missing:
        return Result(
            definition.test_id, definition.group, "FAIL",
            "solar internal-boundary integration is incomplete: " +
            "; ".join(missing), time.monotonic() - started, [])
    return Result(
        definition.test_id, definition.group, "PASS",
        "the fixed R_sun sphere is registered once on every rank before mesh "
        "construction, uses the application surface resolution and absorbing "
        "callback without double deletion, composes with leaf pruning, and "
        "repairs only fully solid AMPS cell measures after cut-cell setup",
        time.monotonic() - started, [])


def _check_coronal_cme_native_wiring(definition: TestDefinition) -> Result:
    """Audit the linked test boundary without claiming an AMPS execution.

    BLDL3D01 and SCCM3D01..07 remain the concrete build/runtime authorities.
    This fast source gate prevents the executable CLI, production archive, and
    read-only lifecycle probe from drifting apart in a source-only checkout.
    """
    started = time.monotonic()
    files = {
        "main.cpp": ROOT / "main.cpp",
        "main_lib.cpp": ROOT / "main_lib.cpp",
        "configuration_io.cpp": ROOT / "runtime" / "configuration_io.cpp",
        "standalone CLI": ROOT / "runtime" / "standalone_command_line.cpp",
        "native evaluator": ROOT / "validation" /
            "coronal_cme_application_test.cpp",
        "makefile": ROOT / "makefile",
    }
    model_root = ROOT.parent / "src" / "models" / "sep_coronal_cme"
    files["coronal-CME particle source"] = model_root / "src" / "particle_source.cpp"
    missing_files = [f"{name}: {path}" for name, path in files.items()
                     if not path.is_file()]
    if missing_files:
        return Result(
            definition.test_id, definition.group, "FAIL",
            "coronal-CME native integration file is absent: " +
            "; ".join(missing_files), time.monotonic() - started, [])

    text = {name: path.read_text(encoding="utf-8")
            for name, path in files.items()}
    required = {
        "main.cpp": (
            "CoronalCmeNativeTests()",
            "SelectCoronalCmeNativeTests(",
            "CaptureNativeApplicationState(",
            "EvaluateCoronalCmeNativeTests(",
            "WriteNativeTestJson(",
            "amps_time_step()",
            "MPI_Barrier(MPI_GLOBAL_COMMUNICATOR)",
        ),
        "main_lib.cpp": (
            "CaptureNativeApplicationState(",
            "ValidateCompleteSample(sample, capabilities)",
            "GetLocalTimeStep(species.ampsIndex)",
            "GetLocalParticleWeight(species.ampsIndex)",
            "FileHasOnlyFiniteNumericTokens(path)",
            "MPI_Allreduce(&localFingerprint",
        ),
        "standalone CLI": (
            'argument == "--test-suite"',
            'argument == "--test-input"',
            'argument == "--test-json"',
            'argument == "--artifact-directory"',
            'argument == "--test-steps"',
            'argument == "--expect-mpi-ranks"',
        ),
        "configuration_io.cpp": (
            "nativeTestRun ? candidate.commandLine.testInputPath",
        ),
        "native evaluator": (
            "srcsep-component-tests-v1",
            '"SCCM3D01"', '"SCCM3D02"', '"SCCM3D03"',
            '"SCCM3D04"', '"SCCM3D05"', '"SCCM3D06"',
            '"SCCM3D07"',
            "ValidateInitializationForOutput(",
            "ValidateAllSpeciesWeights(",
            "ValidateSourceSpecies(",
        ),
        "makefile": (
            "SEP_CORONAL_CME_DIR",
            "SEP_CORONAL_CME_OBJECTS",
            "CORONAL_CME_VALIDATION_OBJ",
            "SEP::CoronalCME::ValidateInitializationForOutput(",
        ),
        "coronal-CME particle source": (
            "item.nucleonCount < 0",
        ),
    }
    missing = [f"{name}: {token}" for name, tokens in required.items()
               for token in tokens if token not in text[name]]
    if missing:
        return Result(
            definition.test_id, definition.group, "FAIL",
            "coronal-CME native integration is incomplete: " +
            "; ".join(missing), time.monotonic() - started, [])
    if "no simulation was started" in text["main.cpp"]:
        return Result(
            definition.test_id, definition.group, "FAIL",
            "main.cpp still rejects native linked tests before AMPS starts",
            time.monotonic() - started, [])
    return Result(
        definition.test_id, definition.group, "PASS",
        "srcSEP3D discovers and executes named linked tests through the normal "
        "AMPS lifecycle, captures finite collective state without mutation, "
        "evaluates public coronal-CME contracts, and emits runner-compatible "
        "JSON evidence from the production archive",
        time.monotonic() - started, [])


def _find_pic_header(args: argparse.Namespace) -> Optional[Path]:
    roots: List[Path] = []
    if args.amps_source is not None:
        roots.append(args.amps_source.expanduser().resolve())
    roots.extend((ROOT.parent, ROOT.parent.parent))
    for root in roots:
        for candidate in (root / "src" / "pic" / "pic.h",
                          root / "pic" / "pic.h", root / "pic.h"):
            if candidate.is_file():
                return candidate
    return None


def _macro_int(text: str, name: str) -> Optional[int]:
    match = re.search(rf"^\s*#\s*define\s+{re.escape(name)}\s+(-?\d+)\b",
                      text, re.MULTILINE)
    return int(match.group(1)) if match else None


def _check_return_codes(definition: TestDefinition,
                        args: argparse.Namespace) -> Result:
    started = time.monotonic()
    adapter = ROOT / "amps" / "amps_mover_status.h"
    adapter_text = adapter.read_text(encoding="utf-8")
    required_adapter_tokens = (
        "static_assert(_PARTICLE_DELETED_ON_THE_FACE_ == 0",
        "static_assert(_PARTICLE_LEFT_THE_DOMAIN_ == 2",
        "static_assert(_PARTICLE_MOTION_FINISHED_ == 3",
        "Core::ParticleMotionOutcome::Advanced",
        "return _PARTICLE_MOTION_FINISHED_",
        "Core::ParticleMotionOutcome::LeftDomain",
        "return _PARTICLE_LEFT_THE_DOMAIN_",
    )
    missing = [token for token in required_adapter_tokens if token not in adapter_text]
    if missing:
        return Result(definition.test_id, definition.group, "FAIL",
                      "adapter contract is incomplete: " + ", ".join(missing),
                      time.monotonic() - started, [])

    pic_header = _find_pic_header(args)
    if pic_header is None:
        return Result(definition.test_id, definition.group, "SKIP",
                      "AMPS pic.h not found; provide --amps-source to verify ABI values",
                      time.monotonic() - started, [])
    pic_text = pic_header.read_text(encoding="utf-8", errors="replace")
    expected = {
        "_PARTICLE_DELETED_ON_THE_FACE_": 0,
        "_PARTICLE_LEFT_THE_DOMAIN_": 2,
        "_PARTICLE_MOTION_FINISHED_": 3,
    }
    observed = {name: _macro_int(pic_text, name) for name in expected}
    mismatches = [f"{name}: expected {value}, found {observed[name]}"
                  for name, value in expected.items() if observed[name] != value]
    if mismatches:
        return Result(definition.test_id, definition.group, "FAIL",
                      "; ".join(mismatches), time.monotonic() - started, [])
    return Result(definition.test_id, definition.group, "PASS",
                  f"adapter mappings match {pic_header}",
                  time.monotonic() - started, [])


def _check_macro_hygiene(definition: TestDefinition,
                         args: argparse.Namespace) -> Result:
    """Compile the application constants after a synthetic AMPS Pi macro."""
    compiler = shutil.which(args.cxx)
    if compiler is None:
        return Result(definition.test_id, definition.group, "ERROR",
                      f"C++ compiler not found: {args.cxx}", 0.0, [])

    # SEP3D.h is included from generated pic.h, so a concrete source/restart
    # include here leaks sep_common into unrelated AMPS interface and mesh
    # translation units.  The particle-adapter declaration is reachable from
    # the same public boundary and must keep InjectionPlan incomplete until its
    # implementation file.  Audit these precise dependency edges first.
    umbrella = (ROOT / "SEP3D.h").read_text(encoding="utf-8")
    particle_adapter = (ROOT / "amps" / "amps_particle_adapter.h").read_text(
        encoding="utf-8")
    forbidden_umbrella = (
        '"adapters/source_runtime.h"',
        '"output/restart.h"',
        '"amps/amps_particle_adapter.h"',
    )
    leaked = [item for item in forbidden_umbrella if item in umbrella]
    if '"../adapters/source_runtime.h"' in particle_adapter:
        leaked.append('amps_particle_adapter.h -> source_runtime.h')
    if leaked:
        return Result(
            definition.test_id, definition.group, "FAIL",
            "AMPS-facing headers expose canonical-model-only dependencies: " +
            ", ".join(leaked), 0.0, [])
    source = ROOT / "test" / "compile_macro_hygiene.cpp"
    command = [compiler, "-std=c++17", "-Wall", "-Wextra", "-Wpedantic",
               "-Werror", "-fsyntax-only", str(source)]
    code, output, elapsed = _run_command(command, ROOT, args.timeout, args.verbose)
    if code != 0:
        return Result(definition.test_id, definition.group, "FAIL",
                      "AMPS Pi macro still corrupts srcSEP3D headers:\n" + output[-3000:],
                      elapsed, command)
    return Result(definition.test_id, definition.group, "PASS",
                  "sep3d_types.h compiles after the AMPS Pi macro is defined",
                  elapsed, command)


def _check_turbulence_header_boundary(definition: TestDefinition,
                                      args: argparse.Namespace,
                                      output_dir: Path) -> Result:
    """Check both transitive headers and the installed AMPS object recipe.

    AMPS compiles the copied build/main/main_lib.cpp with a generic rule from
    Makefile.conf. Several deployed revisions of that rule do not consume the
    CPPFLAGS/CXXFLAGS/INCLUDE variables appended by the application makefile.
    Consequently transitive public headers must stay dependency-light, while
    application-owned objects that legitimately consume sep_common/SWCME must
    receive the canonical directories through the target-scoped compiler
    environment.  The final probe deliberately uses a generic recipe with no
    include variables; it reproduces the command shape used by those AMPS
    installations instead of merely inspecting makefile text.
    """
    compiler = shutil.which(args.cxx)
    if compiler is None:
        return Result(definition.test_id, definition.group, "ERROR",
                      f"C++ compiler not found: {args.cxx}", 0.0, [])

    try:
        common_dir, swcme_dir = _canonical_model_dirs(args)
    except RunnerError as error:
        return Result(definition.test_id, definition.group, "ERROR",
                      str(error), 0.0, [])

    fixture = output_dir / "BLDL3D06-header-boundary"
    fixture.mkdir(parents=True, exist_ok=True)
    provider_source = fixture / "provider_only.cpp"
    bridge_source = fixture / "coefficient_bridge.cpp"
    provider_source.write_text(
        '#include "turbulence/turbulence_models.h"\n'
        'int main() {\n'
        '  SEP3D::Turbulence::PrescribedKolmogorovConfiguration value;\n'
        '  return value.spectralIndex > 1.0 ? 0 : 1;\n'
        '}\n', encoding="utf-8")
    bridge_source.write_text(
        '#include "turbulence/coefficient_bridge.h"\n'
        'int main() { return 0; }\n', encoding="utf-8")

    common_flags = [compiler, "-std=c++17", "-Wall", "-Wextra",
                    "-Wpedantic", "-Werror", "-fsyntax-only",
                    f"-I{ROOT}"]
    provider_command = [*common_flags, str(provider_source)]
    code, output, elapsed = _run_command(
        provider_command, ROOT, args.timeout, args.verbose)
    if code != 0:
        return Result(
            definition.test_id, definition.group, "FAIL",
            "production-facing turbulence_models.h requires an undeclared "
            "sep_common include path:\n" + output[-3000:],
            elapsed, provider_command)

    bridge_command = [*common_flags, f"-I{common_dir}", str(bridge_source)]
    code, output, bridge_elapsed = _run_command(
        bridge_command, ROOT, args.timeout, args.verbose)
    elapsed += bridge_elapsed
    if code != 0:
        return Result(
            definition.test_id, definition.group, "FAIL",
            "coefficient_bridge.h does not compile with canonical "
            "SEP_COMMON_DIR:\n" + output[-3000:], elapsed, bridge_command)

    # Reproduce an installed Makefile.conf whose generic C++ recipe ignores
    # CPPFLAGS, CXXFLAGS, and INCLUDE.  Name the target like a real srcSEP3D
    # production object so the makefile's target-scoped CPLUS_INCLUDE_PATH is
    # exercised.  The probe must find both canonical model headers even though
    # neither directory appears on the compiler command line.
    make = shutil.which("make")
    if make is None:
        return Result(definition.test_id, definition.group, "ERROR",
                      "make is unavailable for the generic-recipe probe",
                      elapsed, provider_command + bridge_command)
    generic_root = fixture / "generic-recipe"
    generic_root.mkdir(parents=True, exist_ok=True)
    # main_lib.o intentionally has no application-local explicit recipe; in
    # production it is exactly the root-level object compiled by the generic
    # Makefile.conf rule shown in AMPS build logs.
    generic_source = generic_root / "main_lib.cpp"
    generic_source.write_text(
        '#include "sep_injection_spectrum.h"\n'
        '#include "swcme_sep_source.hpp"\n'
        'int sep3d_model_header_probe() {\n'
        '  SEP::Injection::Configuration injection;\n'
        '  swcme::sep::SEPSourceState source;\n'
        '  return injection.macroparticlesPerEvent == 0 && !source.active;\n'
        '}\n', encoding="utf-8")
    generic_config = generic_root / "Makefile.conf"
    generic_config.write_text(
        '%.o: %.cpp\n'
        '\t$(CXX) -std=c++17 -Wall -Wextra -Werror -c $< -o $@\n',
        encoding="utf-8")
    generic_command = [
        make, "-f", str(ROOT / "makefile"),
        "main_lib.o",
        f"AMPS_ROOT={ROOT.parent}",
        f"AMPS_CONFIG={generic_config}",
        f"SEP_COMMON_DIR={common_dir}",
        f"SWCME_DIR={swcme_dir}",
        f"CXX={compiler}",
    ]
    code, output, generic_elapsed = _run_command(
        generic_command, generic_root, args.timeout, args.verbose)
    elapsed += generic_elapsed
    if code != 0:
        return Result(
            definition.test_id, definition.group, "FAIL",
            "an AMPS generic recipe that ignores make include variables cannot "
            "compile a srcSEP3D-owned canonical-model consumer:\n" +
            output[-3000:], elapsed, generic_command)

    return Result(
        definition.test_id, definition.group, "PASS",
        "SEP3D.h and the AMPS adapter hide source/restart model dependencies; "
        "the turbulence provider remains sep_common-independent; and an AMPS "
        "generic object recipe that ignores CPPFLAGS/CXXFLAGS/INCLUDE still "
        "receives canonical sep_common/SWCME headers for srcSEP3D objects",
        elapsed, provider_command + bridge_command + generic_command)


def _check_makefile_relocation(definition: TestDefinition,
                               args: argparse.Namespace,
                               output_dir: Path) -> Result:
    """Verify path discovery from both source and AMPS/build/main locations."""
    make = shutil.which("make")
    if make is None:
        return Result(definition.test_id, definition.group, "ERROR",
                      "make is not available", 0.0, [])

    compiler = shutil.which(args.cxx)
    if compiler is None:
        return Result(definition.test_id, definition.group, "ERROR",
                      f"C++ compiler is not available: {args.cxx}", 0.0, [])

    # The temporary tree reproduces the paths involved in discovery and the
    # archive contract checked after an enclosing AMPS build.  It deliberately
    # compiles tiny fixture translation units instead of any production code;
    # BLDL3D01 remains authoritative for a real configured AMPS compilation.
    fixture = output_dir / "BLDL3D05-layout"
    if fixture.exists():
        shutil.rmtree(fixture)
    source_dir = fixture / "srcSEP3D"
    build_dir = fixture / "build" / "main"
    common_dir = fixture / "src" / "models" / "sep_common"
    swcme_dir = fixture / "src" / "models" / "swcme"
    coronal_cme_dir = fixture / "src" / "models" / "sep_coronal_cme"
    source_dir.mkdir(parents=True)
    build_dir.mkdir(parents=True)
    common_dir.mkdir(parents=True)
    swcme_dir.mkdir(parents=True)
    (coronal_cme_dir / "build").mkdir(parents=True)
    (source_dir / "amps").mkdir(parents=True)
    (fixture / "build" / "pic").mkdir(parents=True)
    (fixture / "Makefile.conf").write_text(
        "# Empty BLDL3D05 configuration fixture.\n", encoding="utf-8")

    # The production audit invokes nm on every archive member.  Consequently,
    # zero-byte placeholders are not adequate: each fixture member must be a
    # valid object file.  The mesh fixture also exports the exact C01/C05 ABI
    # names required by audit-production-symbols.  Return types are immaterial
    # to these ordinary C++ mangled names, so minimal local structures keep the
    # routing test independent of all AMPS and SEP3D headers.
    empty_source = fixture / "fixture_empty.cpp"
    mesh_source = fixture / "fixture_mesh_abi.cpp"
    empty_object = fixture / "fixture_empty.o"
    mesh_object = fixture / "fixture_mesh_abi.o"
    empty_source.write_text(
        "// Valid object with no production symbols; used for archive shape.\n"
        "namespace { constexpr int kArchiveFixture = 0; }\n",
        encoding="utf-8")
    mesh_source.write_text(
        "namespace SEP3D {\n"
        "namespace RuntimeModel {\n"
        "struct RunConfiguration3DOptions {};\n"
        "struct StorageLayout {};\n"
        "}\n"
        "namespace Mesh {\n"
        "struct DomainBounds {};\n"
        "struct ResolutionConfiguration {};\n"
        "struct RefinementPreflight {};\n"
        "struct SolarBoundaryGeometry {};\n"
        "DomainBounds MakeDomain(\n"
        "    const RuntimeModel::RunConfiguration3DOptions&) { return {}; }\n"
        "SolarBoundaryGeometry MakeSolarBoundary(\n"
        "    const RuntimeModel::RunConfiguration3DOptions&) { return {}; }\n"
        "int BuildRefinementPreflight(\n"
        "    const DomainBounds&, const ResolutionConfiguration&,\n"
        "    const RuntimeModel::StorageLayout&, RefinementPreflight*) {\n"
        "  return 0;\n"
        "}\n"
        "}\n"
        "}\n"
        "namespace SEP {\n"
        "namespace CoronalCME {\n"
        "struct InitializationLedger {};\n"
        "int ValidateInitializationForOutput(const InitializationLedger&) {\n"
        "  return 0;\n"
        "}\n"
        "}\n"
        "}\n",
        encoding="utf-8")

    commands: List[str] = []
    elapsed = 0.0
    for source, target in ((empty_source, empty_object),
                           (mesh_source, mesh_object)):
        command = [compiler, "-std=c++17", "-c", str(source),
                   "-o", str(target)]
        code, output, duration = _run_command(
            command, fixture, args.timeout, args.verbose)
        commands.extend(command)
        elapsed += duration
        if code != 0:
            return Result(
                definition.test_id, definition.group, "ERROR",
                "could not compile a valid archive fixture object:\n" +
                output[-3000:], elapsed, commands)

    # Keep this list explicit.  A production-manifest change must update this
    # relocation fixture as well, otherwise the archive audit should fail and
    # expose the disagreement rather than silently accepting a partial build.
    application_members = (
        "parker_geometry.o domain_geometry.o mesh_model.o "
        "bg_provider.o bg_parker.o bg_swmf.o background_snapshot.o "
        "turbulence_models.o keyed_random.o time_step.o perpendicular_transport.o "
        "parker_transport.o focused_transport.o population_control.o "
        "run_configuration.o configuration_io.o standalone_command_line.o runtime.o runtime_adapters.o "
        "transport_adapter.o particle_ledger.o swcme_source_adapter.o source_runtime.o "
        "sampling.o observer_runtime.o publication.o restart.o output_coordinator.o "
        "validation_metrics.o coronal_cme_application_test.o main_lib.o "
        "amps_particle_adapter.o")
    # Stage 14 adds the coherent-transport kernel to the coronal archive.
    # Include its flattened object explicitly here as in the production
    # makefile; the layout probe must exercise the current archive contract.
    shared_members = (
        "sep_transport_common.o sep_coefficient_physics.o "
        "sep_coefficient_registry.o sep_background_snapshot.o "
        "sep_test_registry.o sep_injection_spectrum.o sep_species_source.o "
        "swcme3d.o common_sep_field_line_exchange.o "
        "common_sep_field_line_bundle_io.o common_sep_coherent_transport.o")
    (fixture / "Makefile").write_text(
        f"APPLICATION_MEMBERS := {application_members}\n"
        f"SHARED_MEMBERS := {shared_members}\n"
        "amps:\n"
        "\t@mkdir -p build/main\n"
        "\t@for member in $(APPLICATION_MEMBERS) $(SHARED_MEMBERS); do "
        "cp fixture_empty.o build/main/$$member; done\n"
        "\t@cp fixture_mesh_abi.o build/main/mesh_model.o\n"
        "\t@cd build/main && rm -f mainlib.a && "
        "ar -rcs mainlib.a $(APPLICATION_MEMBERS) $(SHARED_MEMBERS)\n"
        "\t@cp fixture_empty.o build/main/main.o\n"
        "\t@cd build/main && rm -f main.a && ar -rcs main.a main.o\n",
        encoding="utf-8")
    shutil.copy2(ROOT / "makefile", source_dir / "makefile")
    shutil.copy2(ROOT / "makefile", build_dir / "makefile")
    shutil.copy2(ROOT / "amps" / "install_mover_hook.py",
                 source_dir / "amps" / "install_mover_hook.py")
    (fixture / "build" / "pic" / "picGlobal.dfn").write_text(
        "#ifndef _PIC_GLOBAL_DEFINITIONS_H_\n"
        "#define _PIC_GLOBAL_DEFINITIONS_H_\n"
        "#define _PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_(ptr,LocalTimeStep,node) "
        "PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory_SecondOrder"
        "(ptr,LocalTimeStep,node);\n"
        "#endif\n", encoding="utf-8")

    expected = {
        f"AMPS_ROOT={fixture.resolve()}",
        f"AMPS_CONFIG={(fixture / 'Makefile.conf').resolve()}",
        f"SEP_COMMON_DIR={common_dir.resolve()}",
        f"SWCME_DIR={swcme_dir.resolve()}",
        f"SEP_CORONAL_CME_DIR={coronal_cme_dir.resolve()}",
    }
    for location in (source_dir, build_dir):
        # Deliberately run from the fixture root, not from either makefile
        # directory.  A relative-to-CWD implementation would fail this probe.
        command = [make, "--no-print-directory", "-f",
                   str(location / "makefile"), "print-layout-paths",
                   f"SWCME_DIR={swcme_dir.resolve()}",
                   f"SEP_CORONAL_CME_DIR={coronal_cme_dir.resolve()}"]
        code, output, duration = _run_command(
            command, fixture, args.timeout, args.verbose)
        commands.extend([f"cwd={fixture}", *command])
        elapsed += duration
        observed = set(output.splitlines())
        expected_for_location = expected | {
            f"SEP3D_MAKEFILE_DIR={location.resolve()}"
        }
        if code != 0:
            return Result(
                definition.test_id, definition.group, "FAIL",
                f"makefile path probe failed from {location}:\n" + output[-3000:],
                elapsed, commands)
        missing = sorted(expected_for_location - observed)
        if missing:
            return Result(
                definition.test_id, definition.group, "FAIL",
                f"makefile resolved incorrect paths from {location}; missing: " +
                ", ".join(missing), elapsed, commands)

    production_command = [
        make, "--no-print-directory", "-f", str(source_dir / "makefile"),
        "strict-production", f"AMPS_ROOT={fixture.resolve()}",
        f"AMPS_CONFIG={(fixture / 'Makefile.conf').resolve()}",
        f"SWCME_DIR={swcme_dir.resolve()}",
        f"SEP_CORONAL_CME_DIR={coronal_cme_dir.resolve()}",
    ]
    code, output, duration = _run_command(
        production_command, output_dir, args.timeout, args.verbose)
    commands.extend([f"cwd={output_dir}", *production_command])
    elapsed += duration
    if code != 0:
        return Result(
            definition.test_id, definition.group, "FAIL",
            "strict-production did not delegate to the enclosing AMPS target "
            "and audit build/main archives:\n" + output[-3000:],
            elapsed, commands)

    return Result(
        definition.test_id, definition.group, "PASS",
        "source srcSEP3D and copied build/main makefiles resolve the same "
        "AMPS root, Makefile.conf, sep_common, SWCME, and sep_coronal_cme "
        "directories; production orchestration delegates to enclosing make amps",
        elapsed, commands)


def _canonical_model_dirs(args: argparse.Namespace) -> Tuple[Path, Path]:
    """Resolve R1 model roots without depending on the process CWD."""
    roots: List[Path] = []
    if args.amps_source is not None:
        roots.append(args.amps_source.expanduser().resolve())
    roots.append(ROOT.parent.resolve())
    for root in roots:
        common = root / "src" / "models" / "sep_common"
        swcme = root / "src" / "models" / "swcme"
        if common.is_dir() and swcme.is_dir():
            return common, swcme
    raise RunnerError(
        "cannot locate canonical src/models/sep_common and src/models/swcme")


def _check_shared_archives(definition: TestDefinition,
                           args: argparse.Namespace) -> Result:
    """Audit canonical R1 archives and this application's ownership manifest.

    srcSEP and srcSEP3D are separately selectable AMPS applications.  A
    srcSEP3D test must therefore remain runnable in a checkout that does not
    install srcSEP at all.  Cross-application ownership is checked by each
    application's own runner, never by reaching sideways into a sibling tree.
    """
    started = time.monotonic()
    try:
        common, swcme = _canonical_model_dirs(args)
    except RunnerError as error:
        return Result(definition.test_id, definition.group, "ERROR", str(error),
                      time.monotonic() - started, [])
    make = shutil.which("make")
    ar = shutil.which("ar")
    if make is None or ar is None:
        return Result(definition.test_id, definition.group, "ERROR",
                      "make and ar are required for the R1 archive audit",
                      time.monotonic() - started, [])

    commands: List[str] = []
    for directory in (common, swcme):
        command = [make, "-C", str(directory), "verify"]
        code, output, _ = _run_command(
            command, ROOT, args.timeout, args.verbose)
        commands.extend(command)
        if code != 0:
            return Result(
                definition.test_id, definition.group, "FAIL",
                f"canonical archive verification failed in {directory}:\n" +
                output[-4000:], time.monotonic() - started, commands)

    expected = {
        common / "sep_common.a": (
            "sep_transport_common.o", "sep_coefficient_physics.o",
            "sep_coefficient_registry.o", "sep_background_snapshot.o",
            "sep_test_registry.o", "sep_injection_spectrum.o",
            "sep_species_source.o"),
        swcme / "swcme.a": ("swcme3d.o",),
    }
    for archive, members in expected.items():
        completed = subprocess.run(
            [ar, "t", str(archive)], text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, check=False)
        observed = tuple(line.strip() for line in completed.stdout.splitlines()
                         if line.strip())
        if completed.returncode != 0 or observed != members:
            return Result(
                definition.test_id, definition.group, "FAIL",
                f"{archive} members {observed!r} do not equal {members!r}",
                time.monotonic() - started, commands)

    makefile = ROOT / "makefile"
    text = makefile.read_text(encoding="utf-8")
    required = ("SEP_COMMON_ARCHIVE", "SWCME_ARCHIVE",
                "SEP_COMMON_OBJECTS", "SWCME_OBJECTS")
    missing = [token for token in required if token not in text]
    if missing:
        return Result(
            definition.test_id, definition.group, "FAIL",
            f"{makefile} does not consume canonical R1 objects: " +
            ", ".join(missing), time.monotonic() - started, commands)

    return Result(
        definition.test_id, definition.group, "PASS",
        "sep_common.a and swcme.a have exact canonical membership and "
        "srcSEP3D consumes their canonical objects without inspecting srcSEP",
        time.monotonic() - started, commands)


def _run_swcme_common(definition: TestDefinition, args: argparse.Namespace,
                      output_dir: Path) -> Result:
    """Run the relocated common 1-D/3-D SWCME acceptance set."""
    try:
        _, swcme = _canonical_model_dirs(args)
    except RunnerError as error:
        return Result(definition.test_id, definition.group, "ERROR", str(error),
                      0.0, [])
    runner = swcme / "test" / "run_tests.py"
    if not runner.is_file():
        return Result(definition.test_id, definition.group, "FAIL",
                      f"relocated SWCME runner is absent: {runner}", 0.0, [])
    command = [sys.executable, str(runner), "--routine", "--output-dir",
               str(output_dir / "swcme-r1")]
    if args.rebuild:
        command.append("--rebuild")
    code, output, elapsed = _run_command(
        command, swcme, args.timeout, args.verbose)
    return Result(
        definition.test_id, definition.group,
        "PASS" if code == 0 else ("FAIL" if code == 1 else "ERROR"),
        ("relocated SWCME 1-D/3-D common suite passed"
         if code == 0 else output[-5000:]), elapsed, command)


def _find_make_config(args: argparse.Namespace) -> Optional[Path]:
    if args.make_config is not None:
        path = args.make_config.expanduser().resolve()
        return path if path.is_file() else None
    for candidate in (ROOT / ".." / ".." / "Makefile.conf",
                      ROOT.parent / "Makefile.conf",
                      ROOT.parent.parent / "Makefile.conf"):
        resolved = candidate.resolve()
        if resolved.is_file():
            return resolved
    return None


def _find_amps_root(args: argparse.Namespace,
                    config: Path) -> Optional[Path]:
    """Find the enclosing checkout that owns the production `make amps` target."""
    candidates: List[Path] = []
    if args.amps_source is not None:
        candidates.append(args.amps_source.expanduser().resolve())
    candidates.extend((config.parent.resolve(), ROOT.parent.resolve()))
    for candidate in candidates:
        if (candidate / "Makefile").is_file():
            return candidate
    return None


def _check_production_build(definition: TestDefinition,
                            args: argparse.Namespace) -> Result:
    config = _find_make_config(args)
    if config is None:
        return Result(definition.test_id, definition.group, "SKIP",
                      "configured AMPS Makefile.conf not found; provide --make-config",
                      0.0, [])
    amps_root = _find_amps_root(args, config)
    if amps_root is None:
        return Result(
            definition.test_id, definition.group, "ERROR",
            "enclosing AMPS Makefile not found; provide --amps-source pointing "
            "to the AMPS root", 0.0, [])
    command = ["make", "-f", str(ROOT / "makefile"), "strict-production",
               f"AMPS_ROOT={amps_root}", f"AMPS_CONFIG={config}"]
    code, output, elapsed = _run_command(command, ROOT, args.timeout, args.verbose)
    if code == 0:
        return Result(definition.test_id, definition.group, "PASS",
                      "enclosing `make amps` completed and build/main archives "
                      "passed retired-symbol and exact shared-member audits",
                      elapsed, command)
    return Result(definition.test_id, definition.group, "FAIL",
                  "configured enclosing AMPS build failed:\n" + output[-4000:],
                  elapsed, command)


def _check_domain_geometry(definition: TestDefinition, args: argparse.Namespace,
                           output_dir: Path) -> Result:
    """Compile real geometry/AMR kernels without missing host/model headers.

    All geometry IDs share one content-addressed probe. Its source/header hash
    makes archive overlays safe even when source timestamps are normalized.
    Results describe portable geometry verification, not an MPI host run.
    """
    sources = [ROOT / name for name in (
        "core/parker_geometry.cpp", "core/domain_geometry.cpp",
        "mesh/mesh_model.cpp", "test/individual-test/domain_geometry_probe.cpp")]
    headers = [ROOT / name for name in (
        "core/parker_geometry.h", "core/domain_geometry.h", "core/sep3d_types.h",
        "mesh/mesh_model.h", "runtime/run_configuration.h")]
    work = output_dir / "domain-geometry-probe"
    work.mkdir(parents=True, exist_ok=True)
    executable = work / "probe"
    signature = work / "source.sha256"
    digest = hashlib.sha256()
    digest.update(str(args.cxx).encode("utf-8"))
    for path in sources + headers:
        digest.update(path.relative_to(ROOT).as_posix().encode("utf-8"))
        digest.update(path.read_bytes())
    identity = digest.hexdigest()
    elapsed = 0.0
    commands = []
    if not executable.is_file() or not signature.is_file() or signature.read_text().strip() != identity:
        command = [args.cxx, "-std=c++17", "-O2", "-Wall", "-Wextra",
                   "-Wpedantic", "-Werror", *map(str, sources), "-o", str(executable)]
        code, output, duration = _run_command(command, ROOT, args.timeout, args.verbose)
        elapsed += duration
        commands.extend(command)
        if code != 0:
            return Result(definition.test_id, definition.group, "ERROR",
                          "portable production-kernel compilation failed:\n" + output[-4000:],
                          elapsed, commands)
        signature.write_text(identity + "\n", encoding="utf-8")
    command = [str(executable), definition.test_id]
    code, output, duration = _run_command(command, ROOT, args.timeout, args.verbose)
    elapsed += duration
    commands.extend(command)
    return Result(definition.test_id, definition.group, "PASS" if code == 0 else "FAIL",
                  output.strip(), elapsed, commands)


def _check_output_boundary(definition: TestDefinition, args: argparse.Namespace,
                           output_dir: Path) -> Result:
    """Compile the actual native callbacks with portable host-service doubles.

    Extraction keeps the crash regression tied to main_lib.cpp, including its
    input guards and owner/root control flow. This is not a native MPI run.
    """
    work = output_dir / "output-boundary-probe"
    work.mkdir(parents=True, exist_ok=True)
    native = (ROOT / "main_lib.cpp").read_text(encoding="utf-8")
    callbacks = []
    for name in ("InterpolateInitializationCellData", "PrintInitializationCellData"):
        match = re.search(r"^void " + name + r"\s*\(", native, re.M)
        if match is None:
            raise RunnerError(f"missing native output callback {name}")
        start = native.index("{", match.start())
        # These functions have balanced braces and no brace-containing strings.
        # Stop at their closing brace, never at an unrelated later declaration.
        depth, end = 1, start + 1
        while depth and end < len(native):
            depth += (native[end] == "{") - (native[end] == "}")
            end += 1
        if depth:
            raise RunnerError(f"unclosed native output callback {name}")
        callbacks.append(native[match.start():end])
    extracted = "\n\n".join(callbacks) + "\n"
    (work / "native_output_callbacks.inc").write_text(extracted, encoding="utf-8")
    sources = [ROOT / "output/sampling.cpp",
               ROOT / "test/individual-test/output_boundary_probe.cpp"]
    digest = hashlib.sha256(str(args.cxx).encode() + extracted.encode())
    # Include transitive headers so a portable probe cannot survive a stale ABI.
    for path in sorted(set(sources + list(ROOT.rglob("*.h")))):
        digest.update(path.relative_to(ROOT).as_posix().encode())
        digest.update(path.read_bytes())
    identity = digest.hexdigest()
    executable, signature = work / "probe", work / "source.sha256"
    elapsed, commands = 0.0, []
    if (not executable.is_file() or not signature.is_file() or
            signature.read_text().strip() != identity):
        command = [args.cxx, "-std=c++17", "-O2", "-Wall", "-Wextra",
                   "-Wpedantic", "-Werror", "-I", str(work),
                   *map(str, sources), "-o", str(executable)]
        code, output, duration = _run_command(command, ROOT, args.timeout, args.verbose)
        elapsed += duration
        commands.extend(command)
        if code != 0:
            return Result(definition.test_id, definition.group, "ERROR",
                          "portable output-callback compilation failed:\n" + output[-4000:],
                          elapsed, commands)
        signature.write_text(identity + "\n", encoding="utf-8")
    command = [str(executable), definition.test_id]
    code, output, duration = _run_command(command, ROOT, args.timeout, args.verbose)
    return Result(definition.test_id, definition.group, "PASS" if code == 0 else "FAIL",
                  output.strip(), elapsed + duration, commands + command)


def _check_application_object_freshness(definition: TestDefinition) -> Result:
    """Protect the copied build/main tree from timestamp-stale C++ ABIs.

    Deterministic source packages normalize mtimes.  Without a force edge,
    make may retain an object compiled against an earlier mesh_model.h when a
    package is overlaid on an existing configured tree.  The final linker then
    sees new callers but an old provider object.  This source gate verifies
    both the forced-rebuild edge and the archive-time ABI/member checks.
    """
    started = time.monotonic()
    makefile = (ROOT / "makefile").read_text(encoding="utf-8")
    required = (
        ".PHONY: FORCE_SEP3D_APPLICATION_OBJECTS",
        "$(MAINLIBOBJ) $(MAINOBJ): FORCE_SEP3D_APPLICATION_OBJECTS",
        "ar -rcs mainlib.a",
        "$(SEP_CORONAL_CME_OBJECTS)",
        "MakeDomain(SEP3D::RuntimeModel::RunConfiguration3DOptions const&)",
        "SEP3D::Mesh::BuildRefinementPreflight(",
        "MakeSolarBoundary(SEP3D::RuntimeModel::RunConfiguration3DOptions const&)",
    )
    missing = [token for token in required if token not in makefile]
    if missing:
        return Result(
            definition.test_id, definition.group, "FAIL",
            "production stale-object guard is incomplete: " +
            "; ".join(missing), time.monotonic() - started, [])
    return Result(
        definition.test_id, definition.group, "PASS",
        "production rebuilds application-owned objects and audits every "
        "mainlib member plus the domain/preflight/solar-boundary mesh ABIs "
        "before link",
        time.monotonic() - started, [])


def _header_guard_range(text: str, macro: str) -> Tuple[int, int, int]:
    """Locate a macro guard's matching endif despite nested conditionals."""
    offset, depth = 0, 0
    start = None
    for line in text.splitlines(keepends=True):
        if start is None:
            if re.match(r"\s*#\s*ifndef\s+" + re.escape(macro) + r"\b", line):
                start, depth = offset, 1
        else:
            directive = re.match(r"\s*#\s*(if|ifdef|ifndef|endif)\b", line)
            if directive:
                depth += -1 if directive.group(1) == "endif" else 1
                if depth == 0:
                    return start, offset, offset + len(line)
        offset += len(line)
    raise ValueError(f"missing or unclosed permanent guard {macro}")


def _check_pic_header_guards(definition: TestDefinition,
                             args: argparse.Namespace,
                             output_dir: Path) -> Result:
    """Verify permanent source guards without any separately shipped script.

    Use the same AMPS-source discovery as the ABI gates. Compile the actual
    domain declaration and PIC model-header tail through distinct source/build
    copies: pragma once alone cannot protect those copies. Dependency stubs
    reproduce recursive halo/gyro includes without requiring MPI or a full PIC
    build. Each negative control removes one correction and must be rejected.
    Generated probes remain in the report directory for failure inspection.
    """
    started = time.monotonic()
    pic_header = _find_pic_header(args)
    if pic_header is None:
        return Result(definition.test_id, definition.group, "SKIP",
                      "AMPS pic.h not found; provide --amps-source to verify guards",
                      time.monotonic() - started, [])
    domain_header = pic_header.parent / "ecsim/domain_bc.h"
    domain_guard = "AMPS_PIC_ECSIM_DOMAIN_BC_H_INCLUDED"
    command: List[str] = []
    try:
        domain = domain_header.read_text(encoding="utf-8", errors="replace")
        pic = pic_header.read_text(encoding="utf-8", errors="replace")
        start, end, _ = _header_guard_range(domain, domain_guard)
        if (not start < domain.index("class cDomainBC") < end or
                not re.search(r"^\s*#\s*define\s+" + domain_guard + r"\b",
                              domain[start:end], re.M)):
            raise ValueError("cDomainBC is not enclosed by its permanent macro guard")
        start, end, _ = _header_guard_range(pic, "_PIC_")
        marker = pic.index("//include headers for individual physical models")
        if not start < marker < end:
            raise ValueError("PIC model-header tail lies outside _PIC_")
        for include in ('#include "ecsim/domain_bc.h"',
                        '#include "gyro/gyro_mover.h"'):
            if include not in pic[marker:end]:
                raise ValueError(f"guarded PIC model-header tail is missing {include}")
        tail = pic[marker:]
        work = (output_dir / "pic-header-guards-probe").resolve()
        work.mkdir(parents=True, exist_ok=True)
        translation = work / "probe.cpp"
        command = [args.cxx, "-std=c++17", "-Wall", "-Wextra", "-Werror",
                   "-fsyntax-only", str(translation)]

        def compile_probe(label: str, success: bool,
                          diagnostic: Optional[str] = None) -> None:
            code, output, _ = _run_command(command, ROOT, args.timeout, args.verbose)
            (work / f"{label}.log").write_text(output, encoding="utf-8")
            if (code == 0) != success or (diagnostic and diagnostic not in output):
                raise ValueError(f"{label} compile contract failed (exit={code}); "
                                 f"log={work / (label + '.log')}:\n" + output[-4000:])

        prefix = "#ifndef _PIC_\n#define _PIC_\nstruct PicBodyMarker {};\n"
        trees = ("src/pic", "build/pic")
        for tree in trees:
            directory = work / tree
            (directory / "ecsim").mkdir(parents=True, exist_ok=True)
            (directory / "pic.h").write_text(prefix + tail, encoding="utf-8")
            (directory / "ecsim/domain_bc.h").write_text(
                domain + f"\n// physical copy: {tree}\n", encoding="utf-8")
            # Distinct file contents prevent compiler pragma-once shortcuts.
            # Only halo/gyro stubs recurse; the declarations expose an early
            # PIC guard endpoint even when the domain's own guard is correct.
            for relative in ("mesh/domain_boundary_bc_collectors.h",
                             "ecsim/halo_sync.h", "parallel/halo_sync.h",
                             "gyro/drift_velocity.h", "gyro/gyro_mover.h"):
                path = directory / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                content = f"#pragma once\n// physical copy: {tree}\n"
                if relative == "ecsim/halo_sync.h":
                    content += (f'#include "{work / "build/pic/pic.h"}"\n'
                                "void SyncJ(int tag=43000);\n")
                elif relative == "gyro/gyro_mover.h":
                    content += (f'#include "{work / "build/pic/pic.h"}"\n'
                                "void GyroMover(int order=1);\n")
                elif relative == "parallel/halo_sync.h":
                    content += "class HaloManager {};\n"
                path.write_text(content, encoding="utf-8")

        fixture = (
            "#include <vector>\n"
            "namespace PIC { namespace Mesh { struct cBoundaryCellInfo {}; "
            "struct cBoundaryCornerNodeInfo {}; } }\n"
            "#define _PIC_MODEL__DUST__MODE_ 0\n"
            "#define _PIC_MODEL__DUST__MODE__ON_ 1\n"
        )
        translation.write_text(fixture + '#include "src/pic/pic.h"\n'
                               '#include "build/pic/pic.h"\n', encoding="utf-8")
        compile_probe("guarded-model-tail", True)

        # Reproduce the historical premature _PIC_ endif; the actual guarded
        # domain remains unchanged so repeated defaults independently expose it.
        _, tail_end, tail_after = _header_guard_range(prefix + tail, "_PIC_")
        broken = (prefix + "#endif\n" +
                  (prefix + tail)[len(prefix):tail_end] + (prefix + tail)[tail_after:])
        for tree in trees:
            (work / tree / "pic.h").write_text(broken, encoding="utf-8")
        compile_probe("unguarded-model-tail", False, "default argument")
        for tree in trees:
            (work / tree / "pic.h").write_text(prefix + tail, encoding="utf-8")
        compile_probe("restored-model-tail", True)

        translation.write_text(fixture + '#include "src/pic/ecsim/domain_bc.h"\n'
                               '#include "build/pic/ecsim/domain_bc.h"\n',
                               encoding="utf-8")
        compile_probe("guarded-domain", True)
        domain_start, domain_end, domain_after = _header_guard_range(domain, domain_guard)
        # Remove just the outer macro pair and define, preserving pragma once
        # and every nested conditional in the production declaration.
        unguarded = domain[:domain_start] + domain[
            domain.index("\n", domain_start) + 1:domain_end] + domain[domain_after:]
        unguarded = re.sub(r"^\s*#\s*define\s+" + domain_guard + r"[^\n]*\n",
                           "", unguarded, flags=re.M)
        for tree in trees:
            (work / tree / "ecsim/domain_bc.h").write_text(
                unguarded + f"\n// physical copy: {tree}\n", encoding="utf-8")
        compile_probe("unguarded-domain", False, "redefinition")
    except (OSError, ValueError) as error:
        return Result(definition.test_id, definition.group, "FAIL", str(error),
                      time.monotonic() - started, command)
    return Result(definition.test_id, definition.group, "PASS",
                  "permanent domain/PIC guards compile through source/build copies; "
                  f"both negative controls fail; probe logs={work}",
                  time.monotonic() - started, command)


def _run_source(definition: TestDefinition, args: argparse.Namespace,
                output_dir: Path) -> Result:
    if definition.test_id == "BLDL3D01":
        return _check_production_build(definition, args)
    if definition.test_id == "BLDL3D02":
        return _check_retired_sources(definition)
    if definition.group == "DOM3D":
        return _check_domain_geometry(definition, args, output_dir)
    if definition.group == "OUT3D":
        return _check_output_boundary(definition, args, output_dir)
    if definition.test_id == "BLDL3D03":
        return _check_return_codes(definition, args)
    if definition.test_id == "BLDL3D04":
        return _check_macro_hygiene(definition, args)
    if definition.test_id == "BLDL3D05":
        return _check_makefile_relocation(definition, args, output_dir)
    if definition.test_id == "BLDL3D06":
        return _check_turbulence_header_boundary(definition, args, output_dir)
    if definition.test_id == "BLDL3D07":
        return _check_application_object_freshness(definition)
    if definition.test_id == "BLDL3D08":
        return _check_initialized_native_background(definition)
    if definition.test_id == "BLDL3D09":
        return _check_active_population_wiring(definition)
    if definition.test_id == "BLDL3D10":
        return _check_solar_boundary_wiring(definition)
    if definition.test_id == "BLDL3D11":
        return _check_coronal_cme_native_wiring(definition)
    if definition.test_id == "BLDL3D12":
        return _check_pic_header_guards(definition, args, output_dir)
    if definition.test_id == "ARCH3D02":
        return _check_shared_archives(definition, args)
    if definition.test_id == "SWCME3D01":
        return _run_swcme_common(definition, args, output_dir)
    if definition.test_id == "RUN3D01":
        command = [sys.executable, str(ROOT / "test" / "test_python_runner.py")]
        code, output, elapsed = _run_command(
            command, ROOT, args.timeout, args.verbose)
        return Result(definition.test_id, definition.group,
                      "PASS" if code == 0 else "FAIL",
                      output.strip() or f"runner unit test exited {code}",
                      elapsed, command)
    if definition.test_id == "VALRUN3D01":
        command = [sys.executable,
                   str(ROOT / "test" / "test_validation_runner.py")]
        code, output, elapsed = _run_command(
            command, ROOT, args.timeout, args.verbose)
        return Result(definition.test_id, definition.group,
                      "PASS" if code == 0 else "FAIL",
                      output.strip() or
                      f"Phase-V runner unit test exited {code}",
                      elapsed, command)
    if definition.test_id == "V2D01":
        if args.sep1d_source is None:
            return Result(definition.test_id, definition.group, "SKIP",
                          "provide --sep1d-source explicitly; srcSEP3D never searches a sibling application",
                          0.0, [])
        command = [sys.executable, str(ROOT / "test" / "run_cross_model_parity.py"),
                   "--sep1d-root", str(args.sep1d_source.expanduser().resolve()),
                   "--output-dir", str(output_dir / "V2D01")]
        code, output, elapsed = _run_command(command, ROOT, args.timeout,
                                             args.verbose)
        return Result(definition.test_id, definition.group,
                      "PASS" if code == 0 else "FAIL",
                      output.strip() or f"cross-model parity exited {code}",
                      elapsed, command)
    if definition.test_id == "V5D01":
        command = [sys.executable, str(ROOT / "test" / "test_v_program.py")]
        code, output, elapsed = _run_command(command, ROOT, args.timeout,
                                             args.verbose)
        return Result(definition.test_id, definition.group,
                      "PASS" if code == 0 else "FAIL",
                      output.strip() or f"V03-V05 governance tests exited {code}",
                      elapsed, command)
    raise RunnerError(f"no source-test implementation for {definition.test_id}")


def _run_validation(definition: TestDefinition, args: argparse.Namespace,
                    output_dir: Path) -> Result:
    """Run one external Phase-V descriptor and preserve its evidence class.

    The nested campaign runner performs checksum, schema, and linked-registry
    checks.  This front end only translates its one-case result into the common
    srcSEP3D summary so ``--all`` can continue after an unavailable dataset or
    a failed event without losing the other development evidence.
    """
    campaign_output = output_dir / "phase-v"
    command = [sys.executable, str(ROOT / "validation" / "run_validation.py"),
               "--case", definition.test_id,
               "--output-dir", str(campaign_output),
               "--timeout", str(args.timeout)]
    if args.amps:
        command.extend(("--amps", str(Path(args.amps).expanduser().resolve())))
    if args.validation_data is not None:
        command.extend(("--evidence-root",
                        str(args.validation_data.expanduser().resolve())))
    if args.validation_input is not None:
        command.extend(("--test-input",
                        str(args.validation_input.expanduser().resolve())))
    if args.validation_launch_prefix:
        command.extend(("--launch-prefix", args.validation_launch_prefix))
    code, output, elapsed = _run_command(
        command, ROOT, args.timeout + 5.0, args.verbose)
    report = campaign_output / definition.test_id / "result.json"
    if not report.is_file():
        return Result(definition.test_id, definition.group, "ERROR",
                      "Phase-V runner wrote no case report: " + output[-2000:],
                      elapsed, command)
    try:
        record = json.loads(report.read_text(encoding="utf-8"))
        status = str(record.get("status", "ERROR")).upper()
        message = str(record.get("message", ""))
    except (OSError, ValueError) as error:
        return Result(definition.test_id, definition.group, "ERROR",
                      f"cannot parse Phase-V case report: {error}",
                      elapsed, command)
    if status not in ("PASS", "FAIL", "SKIP", "ERROR"):
        status = "ERROR"
        message = "Phase-V case report has an invalid status"
    expected_code = 2 if status == "ERROR" else (1 if status == "FAIL" else 0)
    if code != expected_code:
        status = "ERROR"
        message += f" (runner exit {code}, expected {expected_code})"
    return Result(definition.test_id, definition.group, status, message,
                  elapsed, command)


def _select(args: argparse.Namespace) -> List[TestDefinition]:
    modes = sum((bool(args.list), bool(args.tests or args.groups),
                 bool(args.routine), bool(args.all), bool(args.suites)))
    if modes != 1:
        raise RunnerError("choose exactly one mode: --list, --test/--group, "
                          "--routine, --all, or --suite")
    if args.all and (args.tests or args.groups or args.routine or args.suites):
        raise RunnerError("--all cannot be combined with another selector")

    selected: List[TestDefinition] = []
    if args.routine:
        selected = [item for item in TESTS if item.routine]
    elif args.all:
        selected = list(TESTS)
    elif args.suites:
        for suite in args.suites:
            selected.extend(BY_ID[test_id] for test_id in SUITES[suite])
    else:
        for test_id in args.tests:
            key = test_id.upper()
            if key not in BY_ID:
                raise RunnerError(f"unknown test ID '{test_id}'; use --list")
            selected.append(BY_ID[key])
        for group in args.groups:
            key = group.upper()
            if key not in GROUPS:
                raise RunnerError(f"unknown group '{group}'; use --list")
            selected.extend(GROUPS[key])

    # Stable-ID de-duplication mirrors the srcSEP C++ registry.
    unique: Dict[str, TestDefinition] = {}
    for item in selected:
        unique.setdefault(item.test_id, item)
    return [unique[key] for key in sorted(unique)]


def _print_list() -> None:
    print(f"{'ID':<20} {'GROUP':<14} {'KIND':<10} NAME")
    print("-" * 88)
    for item in TESTS:
        print(f"{item.test_id:<20} {item.group:<14} {item.kind:<10} {item.name}")
    print("\nSuites:")
    for name, ids in sorted(SUITES.items()):
        print(f"  {name:<12} {' '.join(ids)}")


def _write_reports(results: Sequence[Result], output_dir: Path) -> None:
    counts = {status: sum(item.status == status for item in results)
              for status in ("PASS", "FAIL", "SKIP", "ERROR")}
    exit_code = 2 if counts["ERROR"] else (1 if counts["FAIL"] else 0)
    payload = {
        "schema": "srcsep3d-tests-v1",
        "generated_utc": _datetime.datetime.now(_datetime.timezone.utc).isoformat(),
        "root": str(ROOT),
        "exit_code": exit_code,
        "totals": {key.lower(): value for key, value in counts.items()},
        "results": [asdict(item) for item in results],
    }
    (output_dir / "srcsep3d-tests.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    cases = []
    for item in results:
        body = ""
        if item.status == "FAIL":
            body = f'<failure message="{_xml_escape(item.message)}"/>'
        elif item.status == "ERROR":
            body = f'<error message="{_xml_escape(item.message)}"/>'
        elif item.status == "SKIP":
            body = f'<skipped message="{_xml_escape(item.message)}"/>'
        cases.append(
            f'  <testcase classname="srcSEP3D.{item.group}" '
            f'name="{item.test_id}" time="{item.elapsed_seconds:.9f}">{body}</testcase>')
    xml = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (f'<testsuite name="srcSEP3D" tests="{len(results)}" '
         f'failures="{counts["FAIL"]}" errors="{counts["ERROR"]}" '
         f'skipped="{counts["SKIP"]}">'),
        *cases, "</testsuite>", "",
    ]
    (output_dir / "srcsep3d-tests.xml").write_text("\n".join(xml), encoding="utf-8")


def _print_summary(results: Sequence[Result]) -> int:
    for item in results:
        print(f"[{item.test_id}] {item.status} ({item.elapsed_seconds:.3f}s) {item.message}")
    counts = {status: sum(item.status == status for item in results)
              for status in ("PASS", "FAIL", "SKIP", "ERROR")}
    print("\nSummary: " + " ".join(f"{key}={value}" for key, value in counts.items()))
    for status in ("FAIL", "ERROR"):
        ids = [item.test_id for item in results if item.status == status]
        print(f"{status}: " + (", ".join(ids) if ids else "none"))
    return 2 if counts["ERROR"] else (1 if counts["FAIL"] else 0)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    if args.timeout <= 0.0:
        raise RunnerError("--timeout must be positive")
    if args.list:
        # Selection validation is still performed so contradictory modes fail
        # consistently instead of silently preferring --list.
        _select(args)
        _print_list()
        return 0

    selected = _select(args)
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if any(item.kind in ("cpp", "shell") for item in selected):
        _build_standalone(args)

    results: List[Result] = []
    for definition in selected:
        if definition.kind == "cpp":
            result = _run_cpp(definition, args, output_dir)
        elif definition.kind == "shell":
            result = _run_shell(definition, args)
        elif definition.kind == "source":
            result = _run_source(definition, args, output_dir)
        elif definition.kind == "validation":
            result = _run_validation(definition, args, output_dir)
        else:
            raise RunnerError(f"unknown test kind for {definition.test_id}: "
                              f"{definition.kind}")
        results.append(result)
        print(f"[{result.test_id}] {result.status}", flush=True)

    _write_reports(results, output_dir)
    return _print_summary(results)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RunnerError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
