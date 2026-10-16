#!/usr/bin/env python3
"""Run linked M2 test I17; policy and artifact handling are shared."""
import sys
import json
import subprocess
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from integration_common import main

# Run launcher/CSV contract checks even when scientific prerequisites are
# missing. Their success is only infrastructure evidence; main() still owns
# I17's scientific PASS/FAIL/SKIPPED result. Keep the checks in the individual
# runner so the aggregate run_all_test.py registry includes them automatically.
checked = subprocess.run([
    sys.executable, "-B", str(Path(__file__).with_name("test_infrastructure.py"))],
    check=False)
if checked.returncode:
    # Preserve the aggregate runner's expected ERROR result instead of leaving
    # an absent result.json to be misread as a scientific failure or skip.
    if "--output-dir" in sys.argv:
        output = Path(sys.argv[sys.argv.index("--output-dir") + 1]) / "I17"
        output.mkdir(parents=True, exist_ok=True)
        (output / "result.json").write_text(json.dumps({
            "id": "I17", "status": "ERROR", "metrics": {},
            "reason": "I17 infrastructure contract failed; see runner log",
        }, indent=2) + "\n")
    raise SystemExit(2)

raise SystemExit(main())
