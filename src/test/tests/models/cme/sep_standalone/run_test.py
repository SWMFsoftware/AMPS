#!/usr/bin/env python3
"""Expose the production srcSEP/core integration cases to the shared runner."""
from pathlib import Path
import json
import os
import sys

TEST_ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(TEST_ROOT))
from lib.runner_support import case, load_context, run_command, write_result, error_result, main_guard


def main():
    context=load_context('Standalone SEP/core integration')
    if context['test']['id']=='CORE-SEP-MATERIAL-CHECKPOINT-001':
        return material_checkpoint(context)
    if context['test']['id']=='CORE-SEP-STAGE7-CHECKPOINT-001':
        return stage7_checkpoint(context)
    native=context['test']['id']=='CORE-SEP-STANDALONE-NATIVE-001'
    try:
        local=Path(context['artifact_dir'])/'local'
        script='run_native_standalone_integration.py' if native else 'run_standalone_integration_tests.py'
        argv=[sys.executable,str(Path(context['repo_root'])/'srcSEP/test'/script),'--output-dir',str(local)]
        if native:
            argv+=['--build-root',os.environ.get('AMPS_FIELD_LINE_NATIVE_BUILD',context['repo_root'])]
        result=run_command(argv,Path(context['repo_root']),timeout=1800 if native else 600)
        payload=json.loads((local/'summary.json').read_text())
        # Freeze the complete executed catalog, including source/boundary
        # cases. Accepting only the earlier stages 0–4 subset would discard
        # new failures at the shared-regression boundary. Keep this list and
        # manifest.expected_cases synchronized; partial reports remain errors.
        expected=(['SINT08','SINT09','SINT13N','SINT14N','SINT15N','SINT16N',
                   'SINT25N','SINT26N','SINT27N','SINT28N','SINT29N','SINT30N','SINT31N','SINT32N','MCP08']
                  if native else
                  ['SINT01','SINT02','SINT03','SINT04','SINT05','SINT06','SINT07','SINT07R',
                   'SINT10','SINT11','SINT12','SINT13','SINT14','SINT15','SINT16',
                   'SINT17','SINT18','SINT19','SINT20','SINT21','SINT22','SINT23','SINT24',
                   'SINT26','SINT27','SINT28','SINT29','SINT30','SINT31','SINT32'])
        rows=payload['results']
        schema='srcsep-native-standalone-integration-v1' if native else 'srcsep-standalone-integration-v1'
        if payload.get('schema')!=schema or [r['test_id'] for r in rows]!=expected:
            raise RuntimeError('integration case registry/schema mismatch')
        if any(r['status'] not in ['PASS','FAIL','ERROR','SKIP'] for r in rows):
            raise RuntimeError('integration case has an unknown result status')
        failed=any(r['status'] in ['FAIL','ERROR'] for r in rows)
        if result.returncode!=int(failed):
            raise RuntimeError('integration process exit disagrees with its results')
        selected=context['test'].get('selected_cases',expected)
        imported=[case(r['test_id'],r['status'],r['message']) for r in rows if r['test_id'] in selected]
        status='ERROR' if any(r['status']=='ERROR' for r in imported) else 'FAIL' if any(r['status']=='FAIL' for r in imported) else 'SKIP' if any(r['status']=='SKIP' for r in imported) else 'PASS'
        return write_result(context,status,payload.get('scope','production adapter/kernel verification; native scope is separate'),cases=imported,
                            evidence=[str(local/'summary.json'),str(local/'commands.json')])
    except Exception as error:
        return error_result(context,error)


def material_checkpoint(context):
    """Expose each material lifecycle ID without discarding any MPI role.

    Each ID aggregates its full requested decomposition/role catalog. A missing
    role is an infrastructure error, even when every remaining role passed.
    The maintained material reader enforces that catalog and preserves ERROR
    separately from a measured numerical FAIL. No old reference is generated.
    """
    try:
        repo=Path(context['repo_root'])
        local=Path(context['artifact_dir'])/'local'
        native=local/'native'
        campaign=local/'campaign'
        first=run_command([sys.executable,str(repo/'srcSEP/test/run_native_standalone_integration.py'),
            '--build-root',os.environ.get('AMPS_FIELD_LINE_NATIVE_BUILD',str(repo)),
            '--output-dir',str(native),'--ranks','1','2','4','8','--prepare-only'],repo,timeout=600)
        if first.returncode:
            raise RuntimeError('native material fixture preparation failed: '+first.stdout)
        result=run_command([sys.executable,str(repo/'srcSEP/test/run_material_checkpoint_tests.py'),
            '--native-prepared',str(native/'prepared.json'),'--output',str(campaign)],repo,timeout=1800)
        # Import the same catalog validator as the shared srcSEP orchestrator;
        # this path must not accept a different subset of the 58 real MPI jobs.
        sys.path.insert(0,str(repo/'tools'))
        from sep_test_orchestrator import material_checkpoint_rows
        rows=material_checkpoint_rows(campaign)
        failed=any(r['status'] in ['FAIL','ERROR'] for r in rows)
        if result.returncode!=int(failed):
            raise RuntimeError('material process exit disagrees with results')
        identities=['MCP01','MCP02','MCP03','MCP04','MCP05','MCP06','MCP07','MCP09','MCP10','MCP11','MCP12']
        selected=context['test'].get('selected_cases',identities)
        imported=[]
        for identity in identities:
            related=[row for row in rows if row['id']==identity]
            status='ERROR' if any(r['status']=='ERROR' for r in related) else 'FAIL' if any(r['status']=='FAIL' for r in related) else 'SKIP' if any(r['status']=='SKIP' for r in related) else 'PASS'
            if identity in selected:
                imported.append(case(identity,status,'; '.join(r.get('message','') for r in related if r.get('message')) or 'all required real-MPI roles verified'))
        status='ERROR' if any(r['status']=='ERROR' for r in imported) else 'FAIL' if any(r['status']=='FAIL' for r in imported) else 'SKIP' if any(r['status']=='SKIP' for r in imported) else 'PASS'
        return write_result(context,status,'moving material checkpoint/repartition; fixed physical dt',cases=imported,
            evidence=[str(campaign/'summary.json'),str(campaign/'prepared.json'),str(campaign/'junit.xml')])
    except Exception as error:
        return error_result(context,error)


def stage7_checkpoint(context):
    """Keep every role visible through the existing shared validation entry.

    Configured worker builds are explicit prerequisites. Missing builds or MPI
    startup failure are infrastructure errors, never an optional physics pass.
    The same catalog reader is used by the application-wide orchestrator.
    """
    try:
        repo=Path(context['repo_root']);local=Path(context['artifact_dir'])/'local'
        command=[sys.executable,str(repo/'srcSEP/test/run_stage7_checkpoint_tests.py'),
                 '--build-root',os.environ.get('AMPS_FIELD_LINE_NATIVE_BUILD',str(repo)),
                 '--output',str(local)]
        result=run_command(command,repo,timeout=2400)
        sys.path.insert(0,str(repo/'tools'))
        from sep_test_orchestrator import stage7_checkpoint_rows
        rows=stage7_checkpoint_rows(local)
        failed=any(r['status'] in ('FAIL','ERROR') for r in rows)
        if result.returncode!=int(failed):raise RuntimeError('stage7 process status disagrees with native role results')
        selected=context['test'].get('selected_cases',sorted({r['id'] for r in rows}))
        imported=[]
        for identity in selected:
            related=[r for r in rows if r['id']==identity]
            if not related:raise RuntimeError('missing selected stage7 case '+identity)
            status='ERROR' if any(r['status']=='ERROR' for r in related) else 'FAIL' if any(r['status']=='FAIL' for r in related) else 'SKIP' if any(r['status']=='SKIP' for r in related) else 'PASS'
            imported.append(case(identity,status,'; '.join(r.get('message','') for r in related if r.get('message')) or 'all required native roles verified'))
        status='ERROR' if any(r['status']=='ERROR' for r in imported) else 'FAIL' if any(r['status']=='FAIL' for r in imported) else 'SKIP' if any(r['status']=='SKIP' for r in imported) else 'PASS'
        return write_result(context,status,'stage7 material/native/thread continuation; isolated wave codec is separate from physical transport',
                            cases=imported,evidence=[str(local/'summary.json'),str(local/'prepared.json'),str(local/'junit.xml')])
    except Exception as error:return error_result(context,error)


if __name__=='__main__':
    main_guard(main)
