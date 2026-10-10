#!/usr/bin/env python3
"""Thin U08 CLI; policy/oracles live in acceptance.json and common runner."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from moon_testlib import run_one_cli

raise SystemExit(run_one_cli("U08"))
