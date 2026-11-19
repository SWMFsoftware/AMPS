#!/usr/bin/env python3
"""Run independent P02 arithmetic/source/owner controls without fake events.

This bounded supporting group is registered with the existing U15 owner. It
tests the independent verifier, not a duplicate production model, and supplies
no physical defaults. Actual event assessment remains the separate explicit
verify_p02.py --binding workflow with authentic completed native artifacts.
"""
import argparse
import json
from pathlib import Path

from test_m3_source_build_protocol import run_control_suite

HERE=Path(__file__).resolve().parent


def run(output):
    return run_control_suite(output,HERE/'sodium_cohort/test_p02.py',
        [HERE/'sodium_cohort/verify_p02.py',HERE/'sodium_cohort/verify_ledger.py',
         HERE/'sodium_cohort/audit.py',HERE/'sodium_cohort/prepare_physical_case.py',
         HERE/'M3_P02_NUMERICAL_BUDGET.json'],runner=__file__,
        schema='srcMoon-P02-independent-verifier-controls-v1',
        scope='independent verifier arithmetic/source/owner controls only; actual events require a completed native owner')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',required=True,type=Path)
    value=run(parser.parse_args().output_dir)
    print(json.dumps(value,indent=2))
    raise SystemExit({'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[value['status']])
