#!/usr/bin/env python3
"""Run linked M2 test I15; policy and artifact handling are shared."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from integration_common import main
from integration_protocol_checks import run_protocol_checks

# These deterministic control-plane fixtures never read lunar observations.
# Include them in the individual runner so the literal global registry runs
# them even while missing D02 timing metadata keeps scientific I15 SKIPPED.
test_root=Path(__file__).resolve().parents[2]
if not run_protocol_checks("I15", (
        test_root / "test_m2_illumination.py",
        test_root / "data_processing/test_i15_spatial_support.py",
        test_root / "test_m2_contracts.py")):
    raise SystemExit(2)

raise SystemExit(main())
