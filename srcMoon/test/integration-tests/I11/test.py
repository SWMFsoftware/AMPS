#!/usr/bin/env python3
"""Run linked M1 AMR-convergence test I11."""
import sys
import subprocess
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from integration_common import main

# Configuration checks are automatically included by run_all_test's existing
# I11 registration. A failure is infrastructure ERROR, never science FAIL/PASS.
checked = subprocess.run([sys.executable, "-B",
    str(Path(__file__).with_name("test_infrastructure.py"))], check=False)
if checked.returncode:
    raise SystemExit(2)

raise SystemExit(main())
