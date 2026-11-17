#!/usr/bin/env python3
"""Run I19 infrastructure contracts, then the linked scientific dispatcher."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from integration_common import main
from integration_protocol_checks import run_protocol_checks

# Preserve the signed-budget checks and register the formal dispatch/import
# fixtures under the SAME I19 entry. The literal global runner executes these
# checks even when formal I19 later SKIPs. Their PASS is infrastructure only;
# it cannot replace the numerical design, native receipts or predecessor gates.
test_root=Path(__file__).resolve().parents[2]
if not run_protocol_checks("I19", (
        Path(__file__).with_name("test_budget.py"),
        test_root / "test_m2_protocol.py",
        test_root / "test_m2_formal_campaign.py")):
    raise SystemExit(2)

raise SystemExit(main())
