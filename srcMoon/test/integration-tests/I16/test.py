#!/usr/bin/env python3
"""Run linked M2 test I16; policy and artifact handling are shared."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from integration_common import main
from integration_protocol_checks import run_protocol_checks

# This local arithmetic/staging contract runs even when I13 or D03 prevents a
# scientific run. Its PASS is infrastructure evidence only; main() still owns
# every frozen I16 science/data/predecessor gate and the final scientific SKIP.
if not run_protocol_checks("I16", (
        Path(__file__).resolve().parents[2] / "test_i16_grey_boundary.py",)):
    raise SystemExit(2)

raise SystemExit(main())
