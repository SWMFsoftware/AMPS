#!/usr/bin/env python3
"""Thin U25 CLI; policy/oracles live in acceptance.json, packages.json and the verifier.

Without options U25 evaluates every dataset declared by the tests that list U25
as a predecessor; --dataset Dxx and --for-test UXX select a specific target.
"""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from package_verification import main

raise SystemExit(main())
