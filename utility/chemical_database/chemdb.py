#!/usr/bin/env python3
"""Run the shared source-tree chemical catalogue CLI."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent / "python"))
from amps_chemical_database.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
