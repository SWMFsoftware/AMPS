"""Approved I19/I20 numerical design and independent diagnostic evaluation.

No physics is advanced here: the executable owns forces, sources, stochastic
surface interactions and weighted packets. This module constructs explicit
runtime receipts and evaluates actual ledgers. Short development differences
are diagnostics, never sixteen-seed confidence or observational validation.
"""
from __future__ import annotations
import csv
from datetime import datetime, timezone
import importlib.util
import itertools
import json
import math
import re
from pathlib import Path
import sys

from integration_common import (prepare_na_ar_input, replace_moon_setting,
                                 moon_setting, parse_receipt, sha256)
from m2_initial_state import (physics_identity,
                             physics_identity_sha256, stage_cluster_template)


def load_protocol(path=None):
    """Read the approved contract; source presence alone is not qualification."""
    path=Path(path) if path else Path(__file__).with_name('M2_I19_I20_PROTOCOL.json')
    value=json.loads(path.read_text())
    if value.get('status')!='APPROVED':raise ValueError('numerical protocol is not approved')
    return value


def grid_pairs():
    """The frozen union of adjacent axes and non-finest-to-finest contrasts."""
    grid = list(itertools.product((1,2,4),repeat=2));pairs=set()
    for i,(dt,w) in enumerate(grid):
        for j,(other_dt,other_w) in enumerate(grid):
            if (dt==other_dt and other_w==2*w) or (w==other_w and other_dt==2*dt):
                pairs.add(tuple(sorted((i,j))))
        if (dt,w)!=(4,4):pairs.add(tuple(sorted((i,8))))
    return grid, sorted(pairs)


def cases(tier='development'):
    """Return every selected case before launching: no result-based omission."""
    if tier not in ('development','cluster'):
        raise ValueError('unknown M2 execution tier')
    grid,_=grid_pairs();count=2 if tier=='development' else 16
    protocol=load_protocol()
    base=protocol['development_design']['seed_base'] if tier=='development' else protocol['seed_base']
    return [dict(key=f'{name}-d{dt}-w{w}-s{j:02d}',case=name,dt_divisor=dt,
                 weight_divisor=w,seed_index=j,seed=base+1009*j,grid_index=index)
            for name in ('source-fed','cohort') for index,(dt,w) in enumerate(grid)
            for j in range(count)]


def without_verification_controls(text):
    """Keep the original positive-coverage adapters separate from the sweep.

    An explicit numerical template may select I19-only optional controls.
    Passing them to I14/I18 would correctly be rejected by production, and
    passing them to the old I19 fixture would silently select the new case.
    Remove only this complete named group in a separately retained receipt;
    all physical surface/thermal/source settings remain in that receipt.
    """
    return re.sub(r"(?m)^\s*m2_verification_[a-z_]+\s*=.*\n?","",text)


def write_case(template, destination, case, tier='development', reference=None,
               elements=None, horizon_s=2550240.0):
    """Preserve .in grammar/physics and add the complete optional control group.

    ``auto`` computes native baseline dt/weights before overrides. The launcher
    then freezes the first full-precision receipt across the remaining grid.
    Physical cohort stock always uses the unrefined native reference weight.
    """
    if tier == 'cluster':
        # Formal coupling must use the reviewed AMR/source/thermal problem.
        # The bounded I18 helper is intentionally development-only: enabling
        # Ar or silently replacing volume cells would invalidate I13 and the
        # reviewed initial-state evidence even if the template hash matched.
        stage_cluster_template(template,destination)
    elif tier == 'development':
        prepare_na_ar_input(template,destination)
    else:
        raise ValueError('unknown M2 execution tier')
    text=destination.read_text()
    # These are named numerical controls: preserve explicit template choices
    # rather than forcing the approved example defaults on every .in receipt.
    prior={m.group(1):m.group(2).strip() for m in re.finditer(
        r"(?m)^\s*(m2_verification_[a-z_]+)\s*=\s*([^!\n]+)",text)}
    text=without_verification_controls(text)
    changes=[('random_seed',str(case['seed']))]
    if tier=='development':
        changes += [('simulation_environment_mode','evolving'),
                    ('surface_inventory_initial_state','zero'),
                    ('surface_inventory_checkpoint_file','none')]
    for key,value in changes:
        text=replace_moon_setting(text,key,value,destination)
    settings={'m2_verification_case':case['case'], 'm2_verification_profile':tier,
        'm2_verification_coarse_steps':str(int(reference['coarse_steps'])) if tier=='cluster' and reference is not None else prior.get('m2_verification_coarse_steps','24'),'m2_verification_dt_divisor':str(case['dt_divisor']),
        'm2_verification_weight_divisor':str(case['weight_divisor']),
        'm2_verification_horizon_s':'short' if tier=='development' else format(horizon_s,'.17g'),
        # The first case must honor an explicit input reference. Later cases
        # use its actual receipt, independently pinning W when dt is refined.
        # Replacing a user-supplied positive reference with auto would silently
        # discard the very numerical control this input group provides.
        'm2_verification_coarse_timestep_s':prior.get('m2_verification_coarse_timestep_s','auto') if reference is None else format(reference['coarse_timestep_s'],'.17g'),
        'm2_verification_na_weight_atoms':prior.get('m2_verification_na_weight_atoms','auto') if reference is None else format(reference['reference_na_weight_atoms'],'.17g'),
        'm2_verification_ar_weight_atoms':prior.get('m2_verification_ar_weight_atoms','auto') if reference is None else format(reference['reference_ar_weight_atoms'],'.17g'),
        'm2_verification_cohort_packets':prior.get('m2_verification_cohort_packets','32'),
        'm2_verification_source_integral_tolerance':prior.get('m2_verification_source_integral_tolerance','1e-8'),
        # A frozen cluster registry is authoritative when explicitly supplied.
        # Development may instead retain its input fixture's site selection.
        # Reusing the template's site 3000 for a different registry would make
        # the normalization and actual initial stocks describe different cases.
        'm2_verification_surface_elements':(','.join(str(e) for e in elements)
            if elements is not None else prior.get('m2_verification_surface_elements','3000'))}
    if tier=='cluster':
        # Actual full-column native state is part of every formal case's
        # evidence. A separate reviewed spin-up proof cannot stand in for the
        # state the process really installed before injecting its first atom.
        settings['m2_verification_record_initial_thermal_state']='true'
    elif 'm2_verification_record_initial_thermal_state' in prior:
        # A supporting/pilot caller can explicitly request the native state
        # without selecting formal execution. Do not strip that request with
        # the rest of the old verification group. Absent retains the native
        # false default; malformed supplied values remain a parser ERROR.
        settings['m2_verification_record_initial_thermal_state']=prior[
            'm2_verification_record_initial_thermal_state']
    lines=text.splitlines();end=[i for i,line in enumerate(lines) if line.strip().startswith('#end application')]
    if len(end)!=1:raise ValueError('expected one application closing boundary')
    lines[end[0]:end[0]]=['! Explicit M2 numerical verification controls; not observational data.']+[
        f'{key} = {value}' for key,value in settings.items()]
    destination.write_text('\n'.join(lines)+'\n')
    if tier=='cluster' and physics_identity(destination)!=physics_identity(template):
        raise ValueError('formal case numerical staging changed frozen physics')


def accounting_module():
    path=Path(__file__).parent/'integration-tests/I19/analyze_budget.py'
    spec=importlib.util.spec_from_file_location('m2_signed_budget',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def evaluate_case(directory, case, tier='development', reference=None):
    """Check actual controls, clocks, exact accounting and passive tag identities."""
    receipt=parse_receipt(directory/'receipt.txt')
    result=accounting_module().evaluate(directory/'m2_budget.csv',require_coverage=False,
                                       allow_numerical_loss=tier=='cluster')
    failures=result['failures'];result.update(check=case['key'],id='I19',scope=f'{tier}-M2-numerical-case',
        formal_i19_status='SKIPPED',runner=str(Path(__file__)),log=str(directory/'stdout.log'))
    expected_steps=int(moon_setting(directory/'application.in','m2_verification_coarse_steps'))*case['dt_divisor']
    if receipt.get('seed')!=case['seed']:failures.append({'control':'seed not applied'})
    if tier=='development' and receipt.get('steps')!=expected_steps:failures.append({'control':'wrong step count'})
    dt=float(receipt['synchronized_timestep_s']);coarse=float(receipt['coarse_timestep_s'])
    for spec,key in (('Na','na'),('Ar','ar')):
        if not math.isclose(float(receipt[key+'_weight_atoms'])*case['weight_divisor'],float(receipt['reference_'+key+'_weight_atoms']),rel_tol=1e-14):
            failures.append({'control':key+' weight changed with dt'})
    if not math.isclose(dt*case['dt_divisor'],coarse,rel_tol=1e-14):failures.append({'control':'dt divisor'})
    if reference:
        for key in ('coarse_timestep_s','reference_na_weight_atoms','reference_ar_weight_atoms'):
            if not math.isclose(float(receipt[key]),float(reference[key]),rel_tol=1e-14):failures.append({'control':'baseline drift '+key})
    rows=list(csv.DictReader((directory/'m2_budget.csv').open()))
    for row in rows:
        if float(row['species_elapsed_s'])!=float(row['elapsed_s']):failures.append({'clock':'species not synchronized'})
        error=abs(float(row['actual_epoch_tdb_s'])-float(row['epoch_tdb_s']))
        bound=16*sys.float_info.epsilon*max(1,abs(float(row['epoch_tdb_s'])))
        if error>bound:failures.append({'clock':'epoch differs from indexed physical time'})
        if case['case']=='cohort' and abs(float(row['emission_atoms'])-float(row['surface_release_atoms']))>float(row['roundoff_bound_atoms']):
            failures.append({'source':'external emission in source-free cohort'})
    # Independently expose cumulative EXTERNAL supply separately from gross
    # adsorption/release. Internal recycling cannot enlarge this denominator.
    external_relative=0.0
    for row in rows:
        external=float(row['emission_atoms'])-float(row['surface_release_atoms'])
        residual=abs(float(row['combined_residual_atoms']))
        if case['case']=='source-fed' and external>0:
            external_relative=max(external_relative,residual/external)
            if residual/external>=.005:failures.append({'external_supply_budget':'imbalance >=0.5%'})
    result['maximum_external_source_relative_imbalance']=external_relative
    write_windows(rows,directory/'window_budget.csv')
    cohorts=directory/'m2_cohorts.csv'
    if case['case']=='cohort':
        cohort_result=evaluate_cohorts(cohorts,directory/'m2_hazards.csv',directory/'application.in')
        result['cohort_metrics']=cohort_result
        failures.extend(cohort_result['failures'])
    if receipt.get('budget_valid')!=1:failures.append({'ledger':'production invariant failed'})
    result.update(status='FAIL' if failures else 'PASS',receipt=receipt,
                  formal_reason='short supporting scope or formal predecessor/data/qualification gate not supplied')
    return result


def write_windows(rows,path):
    """Reconstruct completed-window signed transfers from cumulative receipts.

    Decimal arithmetic preserves serialized physical atoms; differences never
    turn cumulative recycling into an external source. Original rows remain
    intact. This is a reproducible diagnostic, not an independent observation.
    """
    from decimal import Decimal,localcontext
    output=[]
    with localcontext() as context:
        context.prec=60
        for species in ('Na','Ar'):
            group=[r for r in rows if r['species']==species]
            for before,after in zip(group,group[1:]):
                delta=lambda key:Decimal(after[key])-Decimal(before[key])
                residual=delta('gas_atoms')+delta('surface_atoms')-delta('emission_atoms')+delta('adsorption_atoms')-delta('surface_deposit_atoms')+delta('surface_release_atoms')-delta('reservoir_supply_atoms')+delta('escape_atoms')+delta('numerical_removal_atoms')
                output.append({'species':species,'begin_s':before['elapsed_s'],'end_s':after['elapsed_s'],
                    'external_supply_atoms':str(delta('emission_atoms')-delta('surface_release_atoms')),
                    'adsorbed_atoms':str(delta('adsorption_atoms')),'deposited_atoms':str(delta('surface_deposit_atoms')),
                    'released_atoms':str(delta('surface_release_atoms')),'escaped_atoms':str(delta('escape_atoms')),
                    'numerical_loss_atoms':str(delta('numerical_removal_atoms')),'combined_residual_atoms':str(residual)})
    if output:
        with path.open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(output[0]));writer.writeheader();writer.writerows(output)


def evaluate_cohorts(path, hazard_path, application=None):
    """Independent mass/survival/first-residence integration of actual hazards.

    This is conditional on the physical release hazards recorded at runtime,
    not independent validation of the temperature/source law. Existing U12/
    I18 analytical kernels separately check those physical laws. No average
    temperature is substituted into an Arrhenius exponential.
    """
    rows=list(csv.DictReader(path.open()));hazards=list(csv.DictReader(hazard_path.open()))
    failures=[];maximum=0.0
    if application is not None:
        # Independently evaluate the already approved species residence laws
        # using explicit SI constants and the actual diagnostic temperature.
        # This checks rate wiring; it does not validate the temperature itself.
        energy=float(moon_setting(application,'sodium_activation_energy_ev'))*1.602176634e-19
        frequency=float(moon_setting(application,'sodium_attempt_frequency_hz'))
        argon_energy=float(moon_setting(application,'argon_desorption_energy_j_mol'))
        for h in hazards:
            T=float(h['temperature_K']);thermal=float(h['thermal_release_hazard_s-1'])
            expected=(frequency*math.exp(-energy/(1.380649e-23*T)) if h['species']=='Na' else
                      (1.380649e-23*T/6.62607015e-34)*math.exp(-argon_energy/(8.31446261815324*T)))
            if not math.isclose(thermal,expected,rel_tol=1e-12,abs_tol=0):
                failures.append({'cohort':'thermal hazard differs from independent species law'})
            total=float(h['total_release_hazard_s-1']);photon=float(h['photon_release_hazard_s-1'])
            if not math.isclose(total,thermal+photon,rel_tol=1e-14,abs_tol=0):
                failures.append({'cohort':'competing release hazards do not sum'})
    for row in rows:
        n=float(row['initial_atoms']);t=float(row['elapsed_s'])
        # Independent oracle integrates each complete constant-coefficient
        # interval and retains the censored tail, including zero-hazard limits.
        survival=1.;mean=0.
        for h in hazards:
            if h['species']!=row['species'] or h['cohort']!=row['cohort']:continue
            begin=float(h['begin_elapsed_s']);dt=float(h['dt_s']);rate=float(h['total_release_hazard_s-1'])
            if begin>=t:continue
            duration=min(dt,t-begin)
            mean+=survival*(duration if rate==0 else -math.expm1(-rate*duration)/rate)
            survival*=math.exp(-rate*duration)
        actual=float(row['first_survival']);actual_mean=float(row['restricted_first_residence_s'])
        mass_error=abs(float(row['tag_residual_atoms']))/n;maximum=max(maximum,mass_error)
        # Diagnostics accumulate real double packet weights and replicated
        # stock updates; the forward arithmetic bound uses counted operations,
        # including MPI reductions, rather than observed agreement. This is
        # distinct from the approved 1% resolution error allocations.
        if 'arithmetic_operation_bound' in row:
            operations=int(row['arithmetic_operation_bound']);nu=operations*sys.float_info.epsilon
            if not 0<=nu<1:raise ValueError('invalid passive-tracer forward arithmetic bound')
            magnitude=1+sum(abs(float(row[key])) for key in ('surface_fraction','gas_fraction','escape_fraction','numerical_loss_fraction'))
            arithmetic_bound=nu/(1-nu)*magnitude
            if not math.isclose(float(row['tag_roundoff_bound_atoms'])/n,arithmetic_bound,rel_tol=1e-12,abs_tol=0):
                raise ValueError('reported passive-tag roundoff bound differs from independent evaluation')
        else:
            # Tiny exact protocol fixtures contain no real particle sums.
            arithmetic_bound=64*sys.float_info.epsilon*(len(hazards)+1)
        if mass_error>arithmetic_bound:failures.append({'cohort':'tag mass not conserved','relative':mass_error})
        if abs(actual-survival)>arithmetic_bound or abs(actual_mean-mean)>arithmetic_bound*max(1,t):
            failures.append({'cohort':'first residence does not match recorded hazard integral'})
        for key in ('first_survival','surface_fraction','gas_fraction','escape_fraction','numerical_loss_fraction'):
            v=float(row[key])
            if not -arithmetic_bound<=v<=1+arithmetic_bound:failures.append({'cohort':'unbounded '+key})
        if not -arithmetic_bound<=actual_mean<=t+arithmetic_bound*max(1,t):failures.append({'cohort':'invalid censored residence'})
    return {'status':'FAIL' if failures else 'PASS','rows':len(rows),
            'maximum_relative_tag_residual':maximum,'failures':failures,
            'scope':'conditional production hazard integration and actual passive-tag accounting'}


def render_case(directory):
    """Render only existing evaluated production CSVs; never invent skipped data."""
    accounting_module().render(directory/'m2_budget.csv',directory)
    # Keep the common per-stage artifact names as exact copies of the actual
    # evaluated ledger/plot. This is conservation evidence, not observations.
    import shutil
    for source,destination in (('m2_budget.csv','comparison.csv'),
                               ('budget.png','comparison.png'),
                               ('budget.eps','comparison.eps')):
        shutil.copy2(directory/source,directory/destination)
    if not (directory/'m2_cohorts.csv').is_file():return
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows=list(csv.DictReader((directory/'m2_cohorts.csv').open()))
    fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
    for species in ('Na','Ar'):
        group=[r for r in rows if r['species']==species and r['cohort']=='0']
        for key in ('first_survival','surface_fraction','gas_fraction'):
            axes[0].plot([float(r['elapsed_s']) for r in group],[float(r[key]) for r in group],label=species+' '+key)
        axes[1].plot([float(r['elapsed_s']) for r in group],[float(r['restricted_first_residence_s']) for r in group],label=species)
    axes[0].set(xlabel='physical seconds',ylabel='fraction');axes[1].set(xlabel='physical seconds',ylabel='restricted first residence [s]')
    for ax in axes:ax.legend(fontsize=7)
    fig.suptitle('SUPPORTING COHORT VERIFICATION — NOT FORMAL I20 / NOT OBSERVATIONAL VALIDATION')
    for ext in ('png','eps'):fig.savefig(directory/('cohorts.'+ext),dpi=150)
    plt.close(fig)


def execute_design(executable,template,output,*,tier,mpi_ranks,workers,
                   remaining_seconds,on_result=None,elements=None,horizon_s=2550240.0,
                   memory_limit_bytes=200000000000):
    """Run the fixed design through production MPI, with common phase budget.

    Call inside the enclosing process-tree supervisor. Each grid cell completes
    its declared replicates before its diagnostics are published; selected
    cases are never omitted based on their residuals. Every failure keeps its
    actual delegated runner/log. This function does not qualify predecessors.
    """
    from integration_common import manifest
    from m1_execution import execute_levels
    reference=None;results=[];selected=cases(tier)
    for grid_index in range(9):
        group=[case for case in selected if case['grid_index']==grid_index]
        batches=([group[:1],group[1:]] if reference is None else [group])
        for batch in batches:
            jobs={}
            for case in batch:
                directory=output/case['key'];directory.mkdir(parents=True);(directory/'runtime').mkdir()
                application=directory/'application.in'
                write_case(template,application,case,tier,reference,elements,horizon_s)
                command=['mpiexec','-n',str(mpi_ranks),str(executable),
                    '--moon-integration-test','I19','--m2-budget-check',
                    '--integration-output',str(directory/'receipt.txt'),'-input',str(application)]
                identity=manifest('I19',executable,command,case['seed'],application,mpi_ranks)
                identity.update(execution_scope=tier+'-numerical-design',case=case,
                    protocol={'path':str(Path(__file__).with_name('M2_I19_I20_PROTOCOL.json')),
                              'sha256':sha256(Path(__file__).with_name('M2_I19_I20_PROTOCOL.json'))},
                    formal_test_status='SKIPPED',formal_reason='full formal gate assembly is separate')
                if tier=='cluster':
                    # Freeze chronology and actual STAGED physics before the
                    # executor may launch. Native state proof verifies both;
                    # the original template's full hash alone is insufficient.
                    identity['started_utc']=datetime.now(timezone.utc).isoformat().replace('+00:00','Z')
                    identity['physics_identity_sha256']=physics_identity_sha256(application)
                (directory/'run_manifest.json').write_text(json.dumps(identity,indent=2,sort_keys=True)+'\n')
                jobs[case['key']]={'directory':directory,'command':command}
            monitor=output/'monitors'/f'grid-{grid_index}-batch-{len(batch)}';monitor.mkdir(parents=True)
            print(f'RUNNING {tier} M2 grid {grid_index+1}/9 ({len(batch)} cases)',flush=True)
            codes=execute_levels(jobs,monitor,workers,remaining_seconds(),
                memory_gib=memory_limit_bytes/1024**3,render_previews=False)
            for case in batch:
                directory=jobs[case['key']]['directory']
                result={'check':case['key'],'id':'I19','status':'ERROR','scope':tier+'-numerical-design',
                    'runner':str(Path(__file__)),'log':str(directory/'stdout.log')}
                try:
                    if codes.get(case['key'])!=0:raise RuntimeError(f"MPI exit {codes.get(case['key'])}; inspect retained monitor")
                    result.update(evaluate_case(directory,case,tier,reference))
                    if reference is None:reference=result['receipt']
                    render_case(directory)
                except (OSError,ValueError,KeyError,RuntimeError) as error:
                    result.update(status='ERROR',reason=str(error))
                (directory/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
                results.append(result)
                if on_result is not None:on_result(result)
            if any(codes.get(case['key'])!=0 for case in batch):
                raise RuntimeError('selected numerical cases did not complete; retained ERROR evidence')
    return results
