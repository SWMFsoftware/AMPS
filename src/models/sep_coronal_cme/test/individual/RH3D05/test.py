#!/usr/bin/env python3
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parents[3]
raise SystemExit(subprocess.run([sys.executable, str(root / 'test' / 'run_tests.py'), '--test', 'RH3D05'], cwd=root).returncode)
