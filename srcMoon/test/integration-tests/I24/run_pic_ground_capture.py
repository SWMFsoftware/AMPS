#!/usr/bin/env python3
"""Run one actual Na charge-exchange transaction with each mover at 1 and 8 ranks."""
import argparse
import csv
from decimal import Decimal as D, localcontext
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import time
import traceback

from run_ground_capture import (CASE, DEFAULT_INVENTORY, EXIT, HERE, MissingPrerequisite,
                                ROOT, TEST, dump, prepare_inputs, sha, source_uniform)

sys.path.insert(0, str(TEST))
from build_m2_qualification import production_source_identity, sha256
from campaign_budget import CampaignBudget
from run_m3_moon_force_tracers import native_sources
from common.moon_testlib import _probe_compile_command

COHORT = 'I24-PIC-ground-capture-v1'
TEMPLATE = HERE / 'reference/pic_ground_capture.in'
# These exactly representable weights deliberately differ by species. Copying
# the Na packet correction to either product would fail the population checks.
WEIGHT = 2.625
POSITION = (1.e7,1001.,1003.)
MAPPING = {'Na':'23Na:ground','NA_PLUS':'23Na+:ground','H':'1H:ground'}


def binary_spacing(value):
    """Spacing toward increasing magnitude for a finite normal binary64 SI quantity."""
    return math.ldexp(1., math.frexp(abs(value))[1]-53)


def checked_build(path):
    """Require the actual clean-build receipts and the effective three-species ABI."""
    if not path.is_file():
        raise MissingPrerequisite('source-bound producer manifest missing: '+str(path))
    producer = json.loads(path.read_text())
    recipe = producer['source_qualification_recipe']
    if (producer['schema'] != 'srcMoon-actual-normal-restore-producer-v1' or producer['status'] != 'PASS'
            or recipe['packet_layout'] != 'conservative-products'
            or recipe['reference_problem'] != 'homogeneous-zero-acceleration'
            or recipe['mapping']['states'] != MAPPING):
        raise ValueError('expected a clean force-free Na, Na+, H build with individual particle weights')
    if native_sources(producer['production_source_hashes']) != native_sources(production_source_identity()):
        raise ValueError('current authoritative production sources differ from the producer')
    # Hash the retained evidence first. Expiry of the disposable configured
    # checkout must not hide a changed build log, application, or ABI snapshot.
    for record in (producer['build_log'], producer['executable'], recipe, recipe['mapping'], recipe['profile']):
        if sha256(Path(record['path'])) != record['sha256']:
            raise ValueError('changed build evidence: '+record['path'])
    commands = Path(producer['build_log']['path']).read_text()
    if re.findall(r'^EXIT: (-?\d+)$', commands, re.M) != ['0','0']:
        raise ValueError('clean configure/build pair did not succeed')
    library = Path(producer['library']['path'])
    build = library.parent
    if not producer['retained_generated_snapshots']:
        raise ValueError('missing generated ABI evidence')
    members = []
    for name, digest in producer['retained_generated_snapshots'].items():
        generated = build / Path(name).relative_to(path.parent / 'generated')
        if sha(name) != digest:
            raise ValueError('changed retained generated ABI: '+name)
        members.append((generated,digest))
    # A reboot can discard /tmp while the producer's published evidence still
    # verifies. Such a library cannot be relinked. A present artifact whose
    # bytes differ is a corrupt input and remains an infrastructure ERROR.
    if not re.fullmatch(r'[0-9a-f]{64}',producer['library']['sha256']):
        raise ValueError('malformed configured library SHA-256 in producer receipt')
    missing = []
    for member, digest in [(library,producer['library']['sha256']), *members]:
        if not member.exists():
            missing.append(str(member))
        elif sha256(member) != digest:
            raise ValueError('changed configured build member: '+str(member))
    if missing:
        raise MissingPrerequisite('disposable configured build member unavailable: '+', '.join(missing)+
                                  '; recreate the source-bound clean producer and supply its new manifest')
    return build


def rows(folder, table, ranks):
    records = []
    for rank in range(ranks):
        with (folder / f'{table}-rank-{rank}.csv').open() as stream:
            for record in csv.DictReader(stream):
                if None in record or any(value is None for value in record.values()):
                    raise ValueError('malformed '+table+' row')
                records.append(record)
    return records


def audit(folder, ranks, hit, dt, reference, steps=1, weight_scale=1):
    """Independently reconstruct the event and the continued physical packets.

    The test conditions its one trial on a specified source random key to
    ensure both products can be inspected before they leave the mesh. Its
    lifetime estimator belongs to the separately sampled 32,768-parent test.
    """
    live, initial, events = [rows(folder, table, ranks) for table in ('particles','products','events')]
    checks = {}
    checks['one_event_two_committed_and_two_live_particles'] = [len(events),len(initial),len(live)] == [1,2,2]
    if not checks['one_event_two_committed_and_two_live_particles']:
        return checks, {}, {}
    event, = events
    by_species = lambda records: {int(r['species']):r for r in records}
    births, final = map(by_species, (initial, live))
    checks['only_sodium_ion_and_neutral_hydrogen_remain'] = set(births) == set(final) == {1,2}
    if not checks['only_sodium_ion_and_neutral_hydrogen_remain']:
        return checks, {}, {}
    # Native clocks are absolute ET binary64 seconds. Subtraction at a large
    # epoch can lose low bits. Bound only this rounding and endpoint arithmetic.
    origin, epoch = map(float,(event['origin_et'], event['epoch_et']))
    time_bound = 16*max(binary_spacing(value) for value in (origin, epoch, origin+dt, dt))
    # elapsed_s is elapsed time WITHIN the native step that contained the
    # event. Only epoch_et is absolute; the independent exponential draw
    # measures elapsed time from injection. A fresh optical depth drawn on a
    # new native step would move or suppress this event.
    increment = dt/steps
    hit_step = math.ceil(hit/increment)
    hit_within_step = hit-(hit_step-1)*increment
    absolute_epoch_error = abs((epoch-origin)-hit)
    checks['source_key_and_event_epoch'] = (int(event['parent']) == 0
        and abs(float(event['elapsed_s'])-hit_within_step) <= time_bound and absolute_epoch_error <= time_bound)
    by_step = {}
    for sample in rows(folder,'steps',ranks):
        by_step.setdefault(int(sample['step']),[]).append(sample)
    checks['committed_step_history'] = set(by_step)==set(range(1,steps+1))
    # Count actual packets after exchange across all ranks. The interval
    # endpoints in this frozen design lie strictly on either side of t_hit.
    # Each weight population must be unchanged by numerical reference scaling.
    for number, samples in by_step.items():
        before = number < hit_step
        counts = tuple(sum(int(row[key]) for row in samples) for key in
                       ('committed_events','neutral_parents','sodium_ions','hydrogen'))
        population = tuple(sum(float(row[key]) for row in samples) for key in
                           ('represented_Na','represented_Na_ion','represented_H'))
        checks['committed_step_history'] &= (len(samples)==ranks and
            all(float(row['relative_end_s'])==number*increment for row in samples) and
            counts==((0,1,0,0) if before else (1,0,1,1)) and
            population==((WEIGHT,0,0) if before else (0,WEIGHT,WEIGHT)))
    rate = float(reference['rate'])
    checks['single_source_rate_and_exothermic_Q'] = (abs(float(event['rate_s'])/rate-1) < 5e-12
        and abs(float(event['Q_J'])/float(reference['q'])-1) < 5e-12)
    velocity = {}
    max_position_error = 0.
    for species, wref in ((1,56.*weight_scale),(2,.875*weight_scale)):
        begin, end = births[species], final[species]
        v = [float(end['v'+axis]) for axis in 'xyz']
        velocity[species] = v
        prefix = f'species_{species}_'
        checks[prefix+'lineage'] = (int(begin['parent']) == int(end['parent']) == 0
            and int(begin['ordinal']) == int(end['child']))
        checks[prefix+'weights'] = all(float(r['represented']) == WEIGHT
            and float(r['reference_weight']) == wref and float(r['correction'])*wref == WEIGHT for r in (begin,end))
        checks[prefix+'velocity_after_continuation'] = all(float(begin['v'+a]) == v[i] for i,a in enumerate('xyz'))
        checks[prefix+'clock'] = (abs(float(begin['epoch_et'])-epoch) <= time_bound
            and abs(float(begin['remaining_s'])-(hit_step*increment-hit)) <= time_bound
            and abs(float(end['epoch_et'])-(origin+dt)) <= time_bound)
        # Bound the effect of rounding ET on x = x0 + v*(dt-hit) in metres.
        checks[prefix+'continued_position'] = all(float(begin[a]) == POSITION[i]
            and abs(float(end[a])-POSITION[i]-v[i]*(dt-hit)) <= abs(v[i])*time_bound+8*binary_spacing(POSITION[i])
            for i,a in enumerate('xyz'))
        max_position_error=max(max_position_error,*(abs(float(end[a])-POSITION[i]-v[i]*(dt-hit))
                               for i,a in enumerate('xyz')))
    mp, mi, mh, speed = map(float,(reference['proton'],reference['ion'],reference['h'],reference['speed']))
    incident = [mp*speed*f for f in (.6,.8,0.)]
    kinetic = .5*mp*speed**2
    q = float(reference['q'])
    residual = math.dist([mi*a+mh*b for a,b in zip(velocity[1],velocity[2])],incident)/(mp*speed)
    energy_ratio = abs((mi*sum(a*a for a in velocity[1])+mh*sum(a*a for a in velocity[2]))*.5-kinetic-q)/kinetic
    checks['independent_product_momentum_energy'] = max(residual,energy_ratio) < 5e-12
    # The packet count is a magnitude. The ledger stores signed contributions
    # to system changes: a consumed prescribed proton has negative charge,
    # mass, momentum and kinetic-plus-internal energy in the bath. The
    # native charge fields store integer stoichiometry for *one* physical
    # event. Multiply by the separately measured represented populations to
    # check charge conservation for this noninteger packet.
    checks['signed_proton_moments'] = (float(event['bath_number']) == WEIGHT
        and int(event['bath_charge_e']) == -1 and int(event['tracked_charge_e']) == 1
        and float(event['bath_number'])*int(event['bath_charge_e'])
            + float(final[1]['represented'])*int(event['tracked_charge_e']) == 0
        and abs(float(event['bath_mass_kg'])/(-WEIGHT*mp)-1) < 5e-12
        and math.dist([float(event['bath_p'+a]) for a in 'xyz'],[-WEIGHT*v for v in incident])/(WEIGHT*mp*speed) < 5e-12)
    checks['closed_committed_ledger'] = (event['integer_ok'] == event['kinematic_ok'] == '1'
        and abs(float(event['bath_energy_J'])+float(event['tracked_energy_J'])) < WEIGHT*kinetic*5e-12)
    metrics = dict(event_elapsed_s=float(event['elapsed_s']),continued_s=dt-hit,clock_bound_s=time_bound,
                   step_count=steps,step_s=increment,reference_weight_scale=weight_scale,horizon_s=dt,
                   physical_represented_number=WEIGHT,event_epoch_et_s=epoch,source_relative_hit_s=hit,
                   absolute_event_epoch_error_s=absolute_epoch_error,maximum_position_error_m=max_position_error,
                   relative_momentum_error=residual,relative_energy_error=energy_ratio)
    # Ignore MPI ownership and local packet handles, which have no physical
    # meaning across decompositions. Compare semantic IDs and phase-space data.
    physical = [tuple(row[k] for k in ('parent','child','species','represented','x','y','z','vx','vy','vz','epoch_et'))
                for _,row in sorted(final.items())]
    return checks, metrics, physical


def same_endpoint(baseline, current, time_bound):
    """Require the same physical products and recoil in every weight/dt case."""
    if len(baseline)!=len(current):
        return False
    for expected,actual in zip(baseline,current):
        if actual[:4] != expected[:4]:
            return False
        for axis in range(3):
            v=float(actual[7+axis]);v0=float(expected[7+axis])
            x=float(actual[4+axis]);x0=float(expected[4+axis])
            if abs(v-v0) > 64*binary_spacing(max(1.,abs(v),abs(v0))):
                return False
            if abs(x-x0) > 2*max(abs(v),abs(v0))*time_bound + 16*binary_spacing(max(1.,abs(x),abs(x0))):
                return False
        if abs(float(actual[10])-float(expected[10]))>2*time_bound:
            return False
    return True


def write_comparison(folder, records, duration):
    """Plot Cartesian displacements from actual mesh packets and ballistic references.

    Continuous positions in metres are preserved alongside the figure. A
    co-located marker is expected for a successful force-free continuation.
    The coordinates describe this controlled event in LSO.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    with (folder/'continuation.csv').open('w',newline='') as stream:
        out=csv.writer(stream)
        out.writerow(('species','axis','actual_displacement_m','analytical_displacement_m'))
        figure,axes=plt.subplots(1,2,figsize=(8,3.6))
        for axis,row in zip(axes,records):
            # phase-space output field order is defined by audit(). The last
            # value is an absolute ET epoch; displacement uses the local flight.
            species=int(row[2])
            actual=[float(row[4+i])-POSITION[i] for i in range(3)]
            target=[float(row[7+i])*duration for i in range(3)]
            for coord,measured,expected in zip('xyz',actual,target):
                out.writerow((species,coord,measured,expected))
            axis.plot(range(3),target,'o',mfc='none',color='#264b79',ms=10,label='Analytical')
            axis.plot(range(3),actual,'x',color='#ac4615',ms=6,label='PIC')
            axis.set_xticks(range(3))
            axis.set_xticklabels(list('xyz'))
            axis.set_xlabel('LSO component')
            axis.set_title('Na⁺' if species==1 else 'H')
        axes[0].set_ylabel('Displacement after capture (m)')
        axes[0].legend(frameon=False)
        figure.tight_layout()
        for ext in ('png','eps'):
            figure.savefig(folder/('continuation.'+ext),dpi=150)
        plt.close(figure)


def main(args):
    output=args.output_dir.resolve()
    output.mkdir(parents=True,exist_ok=False)
    start=time.monotonic()
    result={'schema':'srcMoon-I24-native-PIC-ground-capture-v1','status':'ERROR','checks':{},'metrics':{},
            'formal_I24_status':'SKIPPED / NOT QUALIFIED'}
    manifest={'source_revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
              'runner':str(Path(__file__).resolve()),'source_cohort':COHORT,'seed':731993,'mpi_layouts':[1,8],
              'openmp_threads':1,'commands':[]}
    try:
        if not 0 < args.wall_seconds <= 3600:
            raise ValueError('development limit is one hour')
        budget=CampaignBudget(args.wall_seconds,min(60,args.wall_seconds/10))
        def execute(command,label,cwd,cores,keep_independent_cases=False):
            folder=output/(label+'-resources')
            folder.mkdir()
            receipt=budget.run(command,cwd,output/(label+'.log'),output=folder,allocated_cores=cores,
                environment=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
            manifest['commands'].append(receipt)
            if receipt.get('reason') or (receipt['returncode'] and not keep_independent_cases):
                raise RuntimeError(label+' did not complete; see '+str(output/(label+'.log')))
            return receipt
        producer=args.producer_manifest.resolve()
        build=checked_build(producer)
        manifest['producer']={'path':str(producer),'sha256':sha(producer)}
        with localcontext() as context:
            context.prec=70
            argv,ref,hashes=prepare_inputs(args.lock,args.mass_inventory,output)
            hit = -D.from_float(source_uniform(731993,COHORT,0)).ln()/ref['rate']
            dt = float(hit+D(1)/16)
            # In the stationary-target LSO fixture energy conservation bounds
            # the fastest possible outgoing H even without an angular law.
            h_speed = ((ref['proton']*ref['speed']**2+2*ref['q'])/ref['h']).sqrt()
        manifest['sources']=hashes
        manifest['test_files']={str(p):sha(p) for p in (Path(__file__).resolve(),HERE/'pic_capture_probe.cpp',
            HERE/'ground_capture_setup.h',HERE/'run_ground_capture.py',TEMPLATE,CASE,TEST/'sodium_cohort/verify_numerical.py')}
        manifest['analytical_first_event_s']=str(hit)
        manifest['interval_s']=dt
        executable=output/'native-pic-ground-capture'
        command=_probe_compile_command(build,executable)
        command.insert(1,'-DMPI_ON')
        command[command.index(str(TEST/'common/production_kernel_probe.cpp'))]=str(HERE/'pic_capture_probe.cpp')
        execute(command,'compile',ROOT,1)
        manifest['probe_sha256']=sha(executable)
        baseline=None
        layouts=[(order,ranks,1,1.,f'{order}-{ranks}mpi') for order in ('first','second') for ranks in (1,8)]
        if args.refinement:
            # Existing second-order runs provide both unrefined anchors. One
            # numerical control varies at a time and the exact binary divisions
            # of the horizon/weights preserve the source and represented count.
            for steps,scale,level in ((2,1.,'steps2'),(4,1.,'steps4'),(1,.5,'weight-half'),(1,.25,'weight-quarter')):
                layouts.extend(('second',ranks,steps,scale,f'second-{level}-{ranks}mpi') for ranks in (1,8))
        manifest['refinement_requested']=args.refinement
        manifest['layouts']=[dict(label=label,order=order,mpi_ranks=ranks,steps=steps,weight_scale=scale)
                             for order,ranks,steps,scale,label in layouts]
        for order,ranks,steps,weight_scale,label in layouts:
            folder=output/label
            folder.mkdir()
            (folder/'PT').mkdir()
            input_file=folder/'moon.in'
            text=TEMPLATE.read_text().replace('force-cut-surface-second-order','force-cut-surface-'+order+'-order')
            text=text.replace('#end application',
                f'hydrogen_characteristic_speed_m_s = {h_speed}\n'
                'hydrogen_to_sodium_reference_weight_ratio = 0.125\n#end application')
            input_file.write_text(text.replace('tracing_maximum_substep_s = 8192','tracing_maximum_substep_s = '+repr(dt)))
            receipt=execute(['mpiexec','-n',str(ranks),str(executable),*argv,str(input_file),str(folder),
                             repr(dt),str(steps),str(weight_scale)],label,folder,ranks,keep_independent_cases=args.refinement)
            if receipt['returncode']:
                # Per-case native aborts are ERRORs. Independent later cases
                # still run within the original resource allowance, preserving
                # the weight evidence when a time partition fails (and vice versa).
                result.setdefault('layout_errors',{})[label]=dict(returncode=receipt['returncode'],log=receipt['log'])
                result['checks'][label]={'completed_native_launch':False}
                result['metrics'][label]={}
                continue
            # Hash actual kernel bytes reported by CSPICE, including the
            # frame kernels needed for this particular LSO initialization.
            paths=(folder/'spice-kernels.txt').read_text().splitlines()
            if not paths:
                raise ValueError('empty list of loaded SPICE kernels')
            manifest.setdefault('spice_sha256',{})[label]={p:sha256(Path(p)) for p in paths}
            checks,metrics,physical=audit(folder,ranks,float(hit),dt,ref,steps,weight_scale)
            if physical:
                write_comparison(folder,physical,dt-float(hit))
            if physical:
                if label=='first-1mpi':
                    baseline=physical
                elif baseline is None:
                    # Preserve the identity of the declared reference even
                    # when its native trial is absent or malformed.
                    checks['first_order_one_rank_anchor_available']=False
                elif steps==1 and weight_scale==1:
                    checks['phase_space_matches_first_order_one_rank']=physical == baseline
                else:
                    # Changing the incoming numerical reference weight must not
                    # select new recoil directions. Compare the physical packet
                    # at a common endpoint, allowing ET and binary64 rounding
                    # in differing step partitions; native handles are local.
                    checks['common_endpoint_matches_unrefined_packets'] = same_endpoint(
                        baseline,physical,metrics['clock_bound_s'])
            result['checks'][label]=checks
            result['metrics'][label]=metrics
        if args.refinement:
            # Retain a row even when an evaluated check fails with missing
            # event metrics. The status column preserves that evaluated FAIL.
            path=output/'refinement.csv'
            columns=('case','order','ranks','steps','weight_scale','status','source_relative_hit_s','event_epoch_et_s',
                     'absolute_event_epoch_error_s','clock_bound_s','maximum_position_error_m','step_s','horizon_s',
                     'physical_represented_number','relative_energy_error','relative_momentum_error')
            with path.open('w',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=columns,extrasaction='ignore')
                writer.writeheader()
                for order,ranks,steps,scale,label in layouts:
                    status=('ERROR' if label in result.get('layout_errors',{}) else
                            'PASS' if all(result['checks'][label].values()) else 'FAIL')
                    writer.writerow(dict(result['metrics'][label],case=label,order=order,ranks=ranks,steps=steps,
                                         weight_scale=scale,status=status))
        result['status']='PASS' if all(all(c.values()) for c in result['checks'].values()) else 'FAIL'
        if result.get('layout_errors'):
            result.update(status='ERROR',reason='native launch error in '+', '.join(result['layout_errors']))
    except MissingPrerequisite as exc:
        result.update(status='SKIPPED',reason=str(exc))
    except Exception as exc:
        (output/'exception.log').write_text(traceback.format_exc())
        result.update(status='ERROR',reason=str(exc))
    result['wall_seconds']=time.monotonic()-start
    result['counts']={s:int(result['status']==s) for s in EXIT}
    manifest['result_files']={str(p.relative_to(output)):sha(p) for p in output.rglob('*')
        if p.is_file() and (p.suffix in ('.csv','.in','.log','.png','.eps','.txt') or p.name=='exception.log')}
    dump(output/'run_manifest.json',manifest)
    dump(output/'result.json',result)
    print(json.dumps(result,indent=2))
    return EXIT[result['status']]


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--producer-manifest',required=True,type=Path)
    parser.add_argument('--lock',required=True,type=Path)
    parser.add_argument('--mass-inventory',type=Path,default=DEFAULT_INVENTORY)
    parser.add_argument('--output-dir',required=True,type=Path)
    parser.add_argument('--refinement',action='store_true',help='add second-order fixed-horizon timestep and fixed-population weight controls at 1 and 8 ranks')
    parser.add_argument('--wall-seconds',type=float,default=900)
    sys.exit(main(parser.parse_args()))
