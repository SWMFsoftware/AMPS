#!/usr/bin/env python3
"""Generate thin individual launchers and human-readable test metadata."""

# Keep this maintenance utility usable with the same Python 3.7+ installation
# as the global runner.  Native ``list[str]`` evaluation requires Python 3.9,
# so annotations are postponed on older HEC login nodes.
from __future__ import annotations

from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
ROADMAP = ROOT / "model" / "testing_validation.md"

def registry():
    # The executable/Python release registry owns stage and backend metadata.
    # Do not maintain another ID-prefix/range table that can omit future tests.
    import importlib.util
    import sys
    spec = importlib.util.spec_from_file_location("sccm_layout_registry", ROOT / "test" / "run_tests.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    ok, message = module.validate_registry()
    if not ok:
        raise ValueError(message)
    return module.TESTS


def description(identifier: str) -> str:
    pattern = re.compile(rf"^- `{re.escape(identifier)}`: (.+?)(?=\n- `[A-Z]|\n\n)",
                         re.MULTILINE | re.DOTALL)
    match = pattern.search(ROADMAP.read_text(encoding="utf-8"))
    return " ".join(match.group(1).split()) if match else "Canonical roadmap gate."

def implementation_notes(test) -> str:
    """Keep generated launchers auditable without duplicating acceptance limits.

    Normative descriptions above come from the specification; these notes name
    the executable evidence and the supplied-family assumptions. They do not
    promote future event-campaign IDs into the software release registry.
    """
    if test.stage == 11:
        return ("\nImplementation: `test/tests_stage11.cpp` calls the public "
                "`discontinuity_transport.h` owning providers. It uses independent "
                "manufactured field/orbit/finite-volume references. Supplied-family "
                "domains, frame invariants, passive-wave policy and integration "
                "limits are documented in the [Stage-11 guide]"
                "(../../../docs/STAGE11_DISCONTINUITY_TRANSPORT.md). "
                "Passing these planar-family checks does not qualify a curved "
                "CME sheath or finite composite transition.\n")
    if test.stage == 12:
        return ("\nImplementation: the identically named unittest class in "
                "`test/test_stage12.py` exercises `tools/preprocessing/`. It uses "
                "synthetic maps, geometry, covariance and candidate products; no "
                "event observation is hardcoded in C++. The [Stage-12 guide]"
                "(../../../docs/STAGE12_OBSERVATION_PREPROCESSING.md) explains "
                "inference assumptions, immutable metadata and role ownership. "
                "The [synthetic CLI example](../../../examples/stage12/README.md) "
                "provides a complete runnable asset/selection workflow.\n")
    if test.stage == 13:
        return ("\nImplementation: `test/test_stage13.py` verifies the release "
                "machinery in [the Stage-13 guide](../../../docs/STAGE13_RELEASE_QUALIFICATION.md). "
                "Actual production qualification additionally requires clean owning-app "
                "builds, MPI/convergence/campaign evidence and D1--D10.\n")
    if test.stage == 14:
        disposition=(" The canonical campaign result is SKIP until actual "
                     "independently registered evidence is available; the synthetic "
                     "contract verification is retained separately." if test.identifier in {"EVT3D01","XMD3D01","SLM3D01"} else "")
        return ("\nImplementation: `test/test_stage14.py` calls public research kernels "
                "and immutable offline producers, with synthetic independent references "
                "and negative fixtures. [The Stage-14 guide]"
                "(../../../docs/STAGE14_RESEARCH_EXTENSIONS.md) records the implemented "
                "domains and outstanding host/campaign gates. A software PASS is not "
                "an observational-campaign or production-adapter qualification."+disposition+"\n")
    return ""

def main() -> int:
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first-stage", type=int, default=0)
    args = parser.parse_args()
    selected = [test for test in registry() if test.stage >= args.first_stage]
    root = ROOT / "test" / "individual"
    for test in selected:
        identifier, stage = test.identifier, test.stage
        folder = root / identifier
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "test.py").write_text("#!/usr/bin/env python3\n"
            "from pathlib import Path\nimport subprocess\nimport sys\n\n"
            "root = Path(__file__).resolve().parents[3]\n"
            f"raise SystemExit(subprocess.run([sys.executable, str(root / 'test' / 'run_tests.py'), '--test', '{identifier}'], cwd=root).returncode)\n",
            encoding="utf-8")
        (folder / "README.md").write_text(f"# {identifier}\n\n{description(identifier)}\n\n"
            "Run from any directory with `python3 test.py`. The launcher delegates "
            "to the global registry, so individual and cumulative gates execute "
            "the identical implementation.\n" + implementation_notes(test), encoding="utf-8")
        expected="skip-without-actual-campaign-evidence" if identifier in {"EVT3D01","XMD3D01","SLM3D01"} else "pass"
        (folder / "reference.json").write_text(f'{{\n  "id": "{identifier}",\n  "first_stage": {stage},\n'
            f'  "expected": "{expected}",\n  "tolerance_authority": "model/testing_validation.md"\n}}\n', encoding="utf-8")
    print(f"generated {len(selected)} individual test launchers")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
