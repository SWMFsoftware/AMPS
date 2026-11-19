#!/usr/bin/env python3
"""Run actual registered M3 supporting services without promoting qualification.

Each child owns its source/data identity, independent oracle and results. This
module only checks the CLI/result protocol and retains runner/log/artifact
locations. It never synthesizes science metrics or replaces a failed child.
"""
from __future__ import annotations
import json
from pathlib import Path
import subprocess
import sys
import time

from campaign_budget import CampaignBudget

ROOT=Path(__file__).resolve().parents[2]
TEST=Path(__file__).resolve().parent
RUNNERS={
    'reaction-primitives':TEST/'test_m3_reaction_execution.py',
    'runtime-input':TEST/'test_m3_runtime_configuration.py',
    'photon-selection':TEST/'common/m3_rates/test_photon_selection.py',
    'rate-registry':TEST/'common/m3_rates/test.py',
    'piecewise-exposure':TEST/'reaction_execution/piecewise_hazard/run.py',
    'exposure-clock-bridge':TEST/'reaction_execution/exposure_bridge/run.py',
    'compensated-depth':TEST/'reaction_execution/compensated_depth/run.py',
    'common-time-clock':TEST/'run_m3_transport_clock.py',
    'qualification-protocol-controls':TEST/'test_m3_qualification_protocol.py',
    # User-approved M3 completion additions retain their local/experimental
    # scope. Neither an actual trajectory-controller check nor an approved
    # statistical design supplies a missing physical source/state provider.
    'trajectory-exposure':TEST/'test_m3_trajectory_exposure.py',
    'physical-reaction-events':TEST/'test_m3_physical_reaction_events.py',
    'profile-state':ROOT/'utility/chemical_database/state_profiles/run.py',
    # Opt-in real atomic-state/crosswalk preparation and new shared product/
    # cascade services extend local verification only. No source availability
    # or numerical PASS changes the formal I21--I24 predecessor contract.
    'atomic-state-profile':ROOT/'utility/chemical_database/state_profiles/run_atomic.py',
    'source-build-protocol':TEST/'test_m3_source_build_protocol.py',
    'physical-case-preparation':TEST/'test_m3_physical_case_preparation.py',
    'p02-independent-controls':TEST/'test_m3_p02.py',
    'commissioning-audit-controls':TEST/'test_m3_commissioning_audit.py',
    'reservoir-products':ROOT/'test/chemical_database/test_reservoir_products.py',
    'oxygen-cascade':TEST/'test_m3_oxygen_cascade.py',
    'qualification-infrastructure':TEST/'test_m3_qualification.py',
    # Pure administrative/independent source-admission controls. Heavy native
    # configuration/runs remain explicit isolated packages, never automatic
    # active-checkout reconfiguration from this global supporting registry.
    'normal-source-registration-controls':TEST/'test_m3_sodium_cohort_support.py',
    'coefficient-stage-controls':TEST/'test_m3_coefficient_stage.py',
    'sodium-ledger':TEST/'test_m3_sodium_ledger.py',
    'sodium-numerical':TEST/'test_m3_sodium_numerical.py',
    'source-compile-recipe':ROOT/'utility/chemical_database/state_profiles/compile_recipe/run.py',
    # Approved P01--P09 component checks exercise maintained production
    # kernels. Their explicit experimental/source-only scopes never promote
    # a full U/I owner whose physical inputs or predecessors remain absent.
    'reduced-daughter-kinematics':TEST/'test_m3_reduced_kinematics.py',
    'ideal-ion-collector':TEST/'test_m3_ion_collector.py',
    'helium-silica-conditional':TEST/'test_m3_helium_silica.py',
    'sodium-partial-capture':TEST/'test_m3_sodium_partial_capture.py',
    'oxygen-state-source-controls':TEST/'test_m3_oxygen_sources.py',
}
EXIT={'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}


def run_services(output,groups,wall_seconds=300):
    """Execute selected groups under ONE supplied remaining wall interval."""
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=False)
    started=time.monotonic();records=[]
    # Reuse ONE monotonic deadline across every group. A subprocess.run
    # timeout only kills its direct Python child; a compiler/native grandchild
    # could otherwise survive and consume resources after this function returns.
    # The existing supervisor tracks verified process birth identities and
    # reaps nested sessions without signalling another test's descendants.
    budget=CampaignBudget(wall_seconds,min(5.0,wall_seconds*.1))
    for name in groups:
        runner=RUNNERS[name];child=output/name;log=output/(name+'.log')
        command=[sys.executable,'-B',str(runner),'--output-dir',str(child)]
        if name=='common-time-clock':command+=['--local-only']
        row={'check':name,'runner':str(runner),'log':str(log),
             'result':str(child/'result.json'),'scope':'supporting-production-service-only',
             'command':command,'status':'ERROR','execution':'not-run'}
        # The intended runner/log remain discoverable even if the ONE budget
        # is exhausted before this group launches or its setup metadata fail.
        log.write_text('SUPPORTING GROUP '+name+'\n')
        try:
            if name=='rate-registry':
                # Preserve the configured primary-publication location across
                # isolated checkouts. Never search for a more convenient source.
                from m2_initial_state import assignments
                record=assignments(TEST.parent/'examples/m3_reaction_rates.in')['charge_exchange_research_directory']
                selected=Path(record['value'])
                if not selected.is_absolute():selected=Path(record['file']).parent/selected
                command+=['--cex-reference-dir',str(selected.resolve())]
            if budget.remaining()<=0:raise RuntimeError('single supporting-service deadline exhausted before launch')
            row['execution']='launched'
            # Observe the WHOLE owned output tree, not only a compiler's tiny
            # receipt directory. Move observer-owned records after each group
            # so later groups cannot overwrite earlier resource evidence. The
            # child's own output directory remains absent until it creates it.
            execution=budget.run(command,ROOT,log,output=output)
            resources=output/'resource-budgets'/name;resources.mkdir(parents=True)
            for filename in ('budget_result.json','budget_resources.csv'):
                (output/filename).rename(resources/filename)
            row.update(returncode=execution['returncode'],
                resource_result=str(resources/'budget_result.json'),
                resource_history=str(resources/'budget_resources.csv'))
            if execution.get('reason') or execution['status']=='ERROR':
                raise RuntimeError(execution.get('reason') or 'abnormal supporting child exit '+str(execution['returncode']))
            value=json.loads((child/'result.json').read_text())
            if not isinstance(value,dict):raise ValueError('supporting child result must be an object')
            status=value.get('status')
            if status not in EXIT or execution['returncode']!=EXIT[status]:
                raise ValueError('child status/exit contradiction')
            row.update(status=status,
                reason=value.get('reason'),child_counts=value.get('counts',{}))
        except (OSError,ValueError,RuntimeError,KeyError,TypeError,ImportError,subprocess.TimeoutExpired) as error:
            row.update(status='ERROR',reason=str(error))
            with log.open('a') as stream:stream.write('ERROR: '+str(error)+'\n')
        records.append(row)
        # A timed-out process cannot borrow another per-child wall budget.
        # Remaining groups are recorded explicitly, never silently omitted.
    status=('ERROR' if any(r['status']=='ERROR' for r in records) else
            'FAIL' if any(r['status']=='FAIL' for r in records) else
            'SKIPPED' if any(r['status']=='SKIPPED' for r in records) else 'PASS')
    result={'schema':'srcMoon-M3-supporting-services-v1','status':status,
        'supporting_checks':records,'wall_seconds':time.monotonic()-started,
        'formal_M3_status':'SKIPPED / NOT QUALIFIED',
        'counts':{s:sum(r['status']==s for r in records) for s in EXIT}}
    (output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
