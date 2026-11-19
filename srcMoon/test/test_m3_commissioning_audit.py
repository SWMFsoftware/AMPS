#!/usr/bin/env python3
"""Register analytical/provenance controls for the physical commissioning audit.

No successful native event fixture is constructed. Actual commissioning audits
require separately retained completed owners; these bounded controls verify the
independent arithmetic and fail-closed parser/identity paths under existing U15.
"""
import argparse
import json
from pathlib import Path
from test_m3_source_build_protocol import run_control_suite

HERE=Path(__file__).resolve().parent


def run(output):
    return run_control_suite(output,HERE/'sodium_cohort/test_audit_commissioning.py',
        [HERE/'sodium_cohort/audit_commissioning.py',HERE/'sodium_cohort/audit.py',
         HERE/'sodium_cohort/verify_ledger.py',HERE/'sodium_cohort/verify_numerical.py',
         HERE/'m3_qualification.py',HERE/'sodium_cohort/verify_p02.py'],
        runner=__file__,schema='srcMoon-commissioning-audit-controls-v1',
        scope='independent analytical/provenance controls; no native physical events or formal qualification')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',required=True,type=Path)
    value=run(parser.parse_args().output_dir)
    print(json.dumps(value,indent=2))
    raise SystemExit({'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[value['status']])
