#!/usr/bin/env python3
"""Offline global-runner entry for maintained P05 NIST source controls.

The default points only to the retained research acquisition. It does not
download data, write canonical data, activate a cascade or normalize missing
branches. A checkout without that packet reports SKIPPED, never PASS.
"""
from pathlib import Path
from state_cascades.test_prepare_sources import main

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PACKET = ROOT / 'test_output/srcMoon/m3-oxygen-cascade-sources-20261118T170845Z/acquired'

if __name__ == '__main__':
    raise SystemExit(main(default_packet=DEFAULT_PACKET, entrypoint=Path(__file__).resolve()))
