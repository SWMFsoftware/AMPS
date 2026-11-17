#!/usr/bin/env python3
"""Run I20 infrastructure contracts, then the linked scientific dispatcher."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from integration_common import main
from integration_protocol_checks import run_protocol_checks

# Temporary metadata/arithmetic/export-control fixtures run from the ordinary
# I20 entry, so run_all_test.py includes them even when I13/D02 keeps science
# SKIPPED. Their success never becomes native or observational qualification.
# The shared helper preserves commands, source hashes, stdout/stderr and an
# actual per-ID infrastructure ERROR result if a subordinate check fails.
test_root=Path(__file__).resolve().parents[2]
if not run_protocol_checks("I20", (
        test_root / "test_m2_polar.py",
        test_root / "test_m2_operator_checks.py")):
    raise SystemExit(2)

raise SystemExit(main())
