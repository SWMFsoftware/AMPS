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

def identifiers() -> list[str]:
    return (["DOCSCCM01", "ARCHSCCM01"] +
        [f"CFG3D{x}" for x in range(12, 21)] + ["RST3D04", "RST3D05"] +
        [f"PFSS3D{x:02d}" for x in range(1, 10)] +
        [f"WND3D{x:02d}" for x in range(1, 21)] +
        [f"CLS3D{x:02d}" for x in range(1, 11)] + ["STR3D01"] +
        [f"SCS3D{x:02d}" for x in range(1, 10)] +
        [f"HCS3D{x:02d}" for x in range(1, 4)] +
        [f"CPL3D{x:02d}" for x in (*range(1, 10), 11, 12)] +
        [f"PLS3D{x:02d}" for x in range(1, 5)] +
        ["OFX3D01", "LOS3D01"] +
        [f"TUR3D{x:02d}" for x in range(1, 8)] +
        [f"MFP3D{x:02d}" for x in range(1, 8)] +
        [f"ELL3D{x:02d}" for x in range(1, 11)] +
        [f"RH3D{x:02d}" for x in range(1, 12)])

def description(identifier: str) -> str:
    pattern = re.compile(rf"^- `{re.escape(identifier)}`: (.+?)(?=\n- `[A-Z]|\n\n)",
                         re.MULTILINE | re.DOTALL)
    match = pattern.search(ROADMAP.read_text(encoding="utf-8"))
    return " ".join(match.group(1).split()) if match else "Canonical roadmap gate."

def main() -> int:
    root = ROOT / "test" / "individual"
    for identifier in identifiers():
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
            "the identical implementation.\n", encoding="utf-8")
        if identifier.startswith(("DOC", "ARCH", "CFG", "RST")):
            stage = 0
        elif identifier.startswith("PFSS"):
            stage = 1
        elif identifier.startswith(("WND", "CLS", "STR")):
            stage = 2
        elif identifier.startswith(("SCS", "HCS", "CPL", "PLS", "OFX", "LOS")):
            stage = 3
        elif identifier.startswith(("TUR", "MFP")):
            stage = 4
        elif identifier.startswith("ELL"):
            stage = 5
        else:
            stage = 6
        (folder / "reference.json").write_text(f'{{\n  "id": "{identifier}",\n  "first_stage": {stage},\n'
            '  "expected": "pass",\n  "tolerance_authority": "model/testing_validation.md"\n}\n', encoding="utf-8")
    print(f"generated {len(identifiers())} individual test launchers")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
