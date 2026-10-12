#!/usr/bin/env python3
"""Run linked M1 angular surface-refinement test I12."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from integration_common import main

raise SystemExit(main())
