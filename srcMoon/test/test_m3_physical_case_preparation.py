#!/usr/bin/env python3
"""Register explicit Na case-preparation controls without running particles.

The maintained unittest module tests parser/provenance/transactional refusal,
including the actual native scalar grammar. Its fixtures are administrative
records, never validation observations. An incomplete scientific configuration
continues to return SKIPPED from the preparer; this runner supplies no defaults.
"""
import argparse
import json
from pathlib import Path

from test_m3_source_build_protocol import run_control_suite

HERE=Path(__file__).resolve().parent


def run(output):
    """Preserve method results/logs and hashes using the normal runner schema."""
    return run_control_suite(output,HERE/'sodium_cohort/test_prepare_physical_case.py',
        [HERE/'sodium_cohort/prepare_physical_case.py',HERE/'sodium_cohort/profile_projection.py',
         HERE/'m2_initial_state.py'],
        runner=__file__,schema='srcMoon-physical-case-preparation-controls-v1',
        scope='explicit case parser/provenance controls only; NO particle execution or physical qualification')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',required=True,type=Path)
    value=run(parser.parse_args().output_dir)
    print(json.dumps(value,indent=2))
    raise SystemExit({'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[value['status']])
