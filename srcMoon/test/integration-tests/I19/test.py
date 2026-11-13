#!/usr/bin/env python3
"""Run linked M2 test I19; policy and artifact handling are shared."""
import sys
import subprocess
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from integration_common import main

code=subprocess.run([sys.executable,"-B",str(Path(__file__).with_name("test_budget.py"))],check=False).returncode
if code:raise SystemExit(code)
raise SystemExit(main())
