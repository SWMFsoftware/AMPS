#!/usr/bin/env python3
"""Source-tree entry for explicit, offline chemical-release publication."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent / "python"))
from amps_chemical_database.build import main

if __name__ == "__main__":
    raise SystemExit(main())
