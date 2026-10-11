#!/usr/bin/env python3
"""Run frozen linked test I07 through the shared M0 harness."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from integration_common import main
raise SystemExit(main())
