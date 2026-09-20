#!/usr/bin/env python3
"""Audit the Stage 1–2 integration boundaries without requiring an AMPS build.

These checks intentionally inspect the production entry points and native test
fixtures.  They complement numerical unit tests: the numerical suites cannot
detect a second driver-level call, an unconditional provider initializer, or a
fixture that bypasses the immutable snapshot phase when the linked executable
is unavailable in a source-only package.
"""

from __future__ import annotations

import pathlib
import re
import sys


ROOT = pathlib.Path(__file__).resolve().parents[1]


def fail(message: str) -> None:
    """Terminate with one stable diagnostic suitable for Make and CI logs."""

    raise AssertionError(message)


def code_without_comments(source: str) -> str:
    """Remove C/C++ comments before counting executable call expressions.

    The audit is about live source, so an API name retained in an explanatory
    comment must not be mistaken for a second invocation.  String literals are
    preserved because none of the audited call spellings are user-facing text.
    """

    source = re.sub(r"/\*.*?\*/", "", source, flags=re.DOTALL)
    return re.sub(r"//[^\n]*", "", source)


def function_body(source: str, signature: str) -> str:
    """Return one brace-balanced C++ function body selected by signature."""

    start = source.find(signature)
    if start < 0:
        fail(f"missing function signature: {signature}")
    opening = source.find("{", start)
    if opening < 0:
        fail(f"missing opening brace after: {signature}")

    depth = 0
    for index in range(opening, len(source)):
        if source[index] == "{":
            depth += 1
        elif source[index] == "}":
            depth -= 1
            if depth == 0:
                return source[opening : index + 1]
    fail(f"unterminated function body: {signature}")
    return ""  # Unreachable, but documents the function's return type.


def require_order(body: str, earlier: str, later: str, context: str) -> None:
    """Require two contract markers to occur in their safety-critical order."""

    earlier_index = body.find(earlier)
    later_index = body.find(later)
    if earlier_index < 0 or later_index < 0 or earlier_index >= later_index:
        fail(f"{context}: expected {earlier!r} before {later!r}")


def main() -> int:
    main_cpp = (ROOT / "main.cpp").read_text(encoding="utf-8")
    main_lib_cpp = (ROOT / "main_lib.cpp").read_text(encoding="utf-8")
    tests_cpp = (ROOT / "tests.cpp").read_text(encoding="utf-8")
    source_compile_scripts = {
        "CV01": ROOT / "test" / "run_cv01_tests.sh",
        "CV02-CV05": ROOT / "test" / "run_cv02_cv05_tests.sh",
        "CV06-CV12": ROOT / "test" / "run_cv06_cv12_tests.sh",
        "IV01-IV06": ROOT / "test" / "run_iv01_iv06_tests.sh",
    }

    # Recommendation 1: the common application timestep owns exactly one live
    # turbulence transaction and the standalone driver owns none.
    advance_pattern = r"SEP::Turbulence::PICAdapter::Advance\s*\("
    if len(re.findall(advance_pattern, code_without_comments(main_lib_cpp))) != 1:
        fail("main_lib.cpp must contain exactly one live PICAdapter::Advance call")
    if re.search(advance_pattern, code_without_comments(main_cpp)):
        fail("main.cpp must not advance turbulence outside amps_time_step()")

    # sep.h already imports the canonical coefficient registry from
    # src/models/sep_common.  Including the retired application-local header in
    # the same translation unit creates a second set of Provider, Ownership,
    # snapshot, and coefficient declarations because the two physical headers
    # intentionally have different include guards.  Keep the production driver
    # on the canonical dependency surface only.
    if '#include "util/sep_coefficient_registry.h"' in main_cpp:
        fail("main.cpp includes the retired local coefficient registry")

    # Recommendation 2: startup initialization is explicitly source-aware.
    executable_main = code_without_comments(main_cpp)
    if "ModelInit::Init()" in executable_main:
        fail("legacy unconditional ModelInit::Init remains executable")
    for marker in (
        "SelfConsistentIntegrated",
        "SelfConsistentSpectral",
        "ConfiguredProvider()",
        "Provider::Analytic",
        "Provider::Swcme",
    ):
        if marker not in executable_main:
            fail(f"provider/source ownership marker is absent: {marker}")

    # Recommendation 4: each native mover fixture enters the immutable read
    # phase after its temporary setup and before its first production move.
    fixture_contracts = (
        ("bool ParkerModelMoverTest_const_plasma_field()", "ParticleMover_Parker"),
        ("bool ParkerModelMoverTest_convection()", "ParticleMover_Parker"),
        ("bool FTE_Convectoin()", "ParticleMover_FocusedTransport_Dmumu"),
    )
    for signature, mover in fixture_contracts:
        body = function_body(tests_cpp, signature)
        require_order(body, "BeginParticleRead", mover, signature)

    # Recommendation 5: every detached controlled/integrated native compile
    # must resolve the canonical repository-root include path.  The scripts use
    # either src_root or root for their absolute srcSEP directory variable.
    for label, path in source_compile_scripts.items():
        script = path.read_text(encoding="utf-8")
        if '-I"$src_root/.."' not in script and '-I"$root/.."' not in script:
            fail(f"{label} compile gate omits the canonical repository root")

    print("Stage 1–2 source ownership and fixture contracts: PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AssertionError as error:
        print(f"Stage 1–2 contract audit: FAIL: {error}", file=sys.stderr)
        raise SystemExit(1)
