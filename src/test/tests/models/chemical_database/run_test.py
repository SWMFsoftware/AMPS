#!/usr/bin/env python3
"""Import the maintained shared data-access suite without site/science data.

The one CORE identity owns the complete protocol software scope, including
actual native compilation and Python API/CLI controls. Its child records stay
in result.json and logs; real release/observational qualification is separate.
"""
from __future__ import annotations
import json
from pathlib import Path
import sys

TEST_ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(TEST_ROOT))
from lib.runner_support import (RunnerError,case,command,error_result,load_context,
                                main_guard,run_command,write_result)

EXPECTED_CASES=['CHEMICAL-DATABASE-PROTOCOL-001']


def main():
    context=load_context('Shared immutable chemical database software controls')
    try:
        configured=context['test'].get('expected_cases',[])
        selected=context['test'].get('selected_cases',configured)
        if configured!=EXPECTED_CASES or selected!=EXPECTED_CASES:
            raise RunnerError('chemical database complete case registry/selection drift')
        repo=Path(context['repo_root']);output=Path(context['artifact_dir'])/'shared-protocol'
        argv=command(context,'python')+['-B',str(repo/'test/chemical_database/run_all_tests.py'),
            '--output-dir',str(output),'--protocol-only','--cxx-json',json.dumps(command(context,'cxx'))]
        completed=run_command(argv,repo,timeout=int(context['test'].get('timeout_seconds',300)))
        path=output/'result.json'
        if not path.is_file():raise RunnerError('shared protocol runner did not publish result.json')
        result=json.loads(path.read_text())
        if result.get('schema')!='amps-chemical-database-software-tests-v1' or result.get('protocol_only') is not True:
            raise RunnerError('wrong shared protocol receipt schema/scope')
        records=result.get('tests');counts=result.get('counts');status=result.get('status')
        if not isinstance(records,list) or not records or not isinstance(counts,dict):
            raise RunnerError('shared protocol receipt lacks complete case outcomes/counts')
        actual={name:sum(row.get('status')==name for row in records) for name in ('PASS','FAIL','ERROR','SKIPPED')}
        if (any(row.get('status') not in actual for row in records) or actual!=counts or
            result.get('executed')!=len(records)):
            raise RunnerError('shared protocol child record/count contradiction')
        expected='ERROR' if counts['ERROR'] else 'FAIL' if counts['FAIL'] else 'PASS' if counts['PASS'] else 'SKIPPED'
        if status!=expected or completed.returncode!={'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[status]:
            raise RunnerError('shared protocol status/exit contradicts evaluated records')
        # In this explicit source-only selection, a dependency SKIP remains a
        # CORE SKIP even if another available software method passed.
        core_status='SKIP' if status=='SKIPPED' or (status=='PASS' and counts['SKIPPED']) else status
        return write_result(context,core_status,'shared software/API checks only; no physical validation',
            cases=[case(EXPECTED_CASES[0],core_status,'actual maintained protocol: '+str(counts))],
            evidence=[str(path),str(output/'unittest.log'),str(output/'build.json')],
            metrics={'child_records':len(records),'child_counts':counts,'observational_validation':False})
    except Exception as error:return error_result(context,error)


if __name__=='__main__':main_guard(main)
