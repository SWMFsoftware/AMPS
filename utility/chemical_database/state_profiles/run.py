#!/usr/bin/env python3
"""Bounded real-source controls for the production M3 state-profile loader.

This supporting component can PASS while actual physical reaction families
remain SKIPPED. It never configures AMPS, changes compiled masses or the
installed archive, fetches data, or substitutes its arithmetic for physics.
"""
from __future__ import annotations

import argparse
import copy
import csv
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

import prepare


ROOT = Path(__file__).resolve().parents[3]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def arguments(executable, profile):
    command = [str(executable)]
    for name in ('isotopes','constants','states'):
        command.extend([profile[name]['path'],profile[name]['sha256']])
    for source in profile['sources']:
        command.extend([source[k] for k in ('id','path','sha256','role')])
    return command


def run(args):
    started = time.monotonic()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError('test output must not exist')
    output.mkdir(parents=True)
    tests = []
    manifest = dict(schema='srcMoon-M3-state-profile-run-v1',commands=[],
        source_revision=subprocess.run(['git','rev-parse','HEAD'],cwd=ROOT,text=True,
                                       capture_output=True).stdout.strip(),
        application_configuration_modified=False, physical_activation=False)

    def check(name, passed):
        tests.append(dict(name=name,status='PASS' if passed else 'FAIL'))

    def execute(command, name):
        logfile = output / (name+'.log')
        manifest['commands'].append(dict(argv=command,cwd=str(ROOT),log=str(logfile)))
        result = subprocess.run(command,cwd=ROOT,text=True,capture_output=True,timeout=60)
        logfile.write_text(result.stdout+result.stderr)
        return result

    try:
        profile = prepare.prepare(args.raw_dir,args.decision_register,output/'profile',args.canonical_root)
        native = [ROOT/'srcMoon/species/PhysicalStateProfile.cpp',
                  ROOT/'srcMoon/species/ReducedKinematics.cpp',
                  ROOT/'src/models/chemical_database/DataIO.cpp',
                  Path(__file__).with_name('probe.cpp')]
        identities = native + [ROOT/'srcMoon/species/PhysicalStateProfile.h',
            ROOT/'srcMoon/species/ReducedKinematics.h',ROOT/'srcMoon/species/ReactionExecution.h',ROOT/'srcMoon/species/Photochemistry.h',
            ROOT/'src/models/chemical_database/DataIO.h']
        manifest['native_source_sha256'] = {str(p):digest(p) for p in identities}
        data_files = [Path(s['path']) for s in profile['sources']]
        before = {str(p):(digest(p),p.stat().st_mtime_ns) for p in data_files}
        compiler = shutil.which('g++')
        if not compiler:
            raise FileNotFoundError('g++ is unavailable')
        executable = output/'production-profile-probe'
        command = [compiler,'-std=c++17','-O2','-Wall','-Wextra','-Werror','-pedantic',
                   '-I',str(ROOT)] + [str(p) for p in native] + ['-o',str(executable)]
        compiled = execute(command,'compile')
        if compiled.returncode:
            raise RuntimeError('production profile compiler failed; see compile.log')
        manifest['compiler'] = subprocess.run([compiler,'--version'],text=True,capture_output=True).stdout.splitlines()[0]
        manifest['executable_sha256'] = digest(executable)
        actual = execute(arguments(executable,profile),'native-real-profile')
        if actual.returncode not in (0,1):
            raise RuntimeError('production profile rejected real source; see native-real-profile.log')
        lines = [line.split('\t') for line in actual.stdout.splitlines()]
        native_checks = [row for row in lines if row[0]=='CHECK']
        for _,name,status in native_checks:
            check('native:'+name,status=='PASS')
        check('native control exit agrees with declared checks',
              actual.returncode==0 and bool(native_checks) and
              all(r[2]=='PASS' for r in native_checks))
        epsilon = Decimal(next(row[1] for row in lines if row[0]=='EPSILON'))
        # The independent prior extractor is a research reference, not a
        # managed copy. Select ONLY its permitted isotope/constant quantities;
        # no ASD threshold is copied to this profile or used as a state energy.
        legacy = prepare.strict_json(args.reference_inventory)
        manifest['independent_reference'] = dict(path=str(args.reference_inventory.resolve()),
                                                  sha256=digest(args.reference_inventory))
        expected = {r['identity']+':ground':Decimal(r['reference_neutral_mass_kg'])
                    for r in legacy['isotopes']}
        expected['electron:free'] = Decimal(legacy['constants']['electron mass']['value'])
        comparisons = []
        for row in lines:
            if row[0] in ('STATE','PHYSICAL'):
                value, reference = Decimal(row[2]), expected[row[1]]
                with localcontext() as context:
                    context.prec=80
                    relative=abs(value-reference)/abs(reference)
                passed=relative<=2*epsilon
                check('independent-Decimal:'+row[0]+':'+row[1],passed)
                comparisons.append([row[0],row[1],str(reference),row[2],str(relative),str(2*epsilon),'PASS' if passed else 'FAIL'])
            elif row[0]=='CONSTANT':
                reference=Decimal(legacy['constants'][row[1]]['value'])
                with localcontext() as context:
                    context.prec=80;relative=abs(Decimal(row[2])-reference)/abs(reference)
                check('independent-constant:'+row[1],relative<=2*epsilon and
                      row[3]==legacy['constants'][row[1]]['unit'])
        with (output/'comparison.csv').open('w',newline='') as stream:
            writer=csv.writer(stream);writer.writerow(['quantity','state_id','independent_kg','native_kg','relative_error','two_machine_epsilon_limit','status']);writer.writerows(comparisons)
        check('actual profile retains all fourteen unresolved states',
              len([row for row in lines if row[0]=='UNAVAILABLE'])==14)
        with Path(profile['states']['path']).open(newline='') as stream:
            state_rows=list(csv.reader(stream))
        with Path(profile['isotopes']['path']).open(newline='') as stream:
            isotope_rows=list(csv.reader(stream))
        variants=[]
        def variant(name,table,rows,edit):
            changed=copy.deepcopy(rows);edit(changed)
            path=output/(name+'.csv')
            with path.open('w',newline='') as stream:
                csv.writer(stream).writerows(changed)
            staged=copy.deepcopy(profile);staged[table]=dict(path=str(path),sha256=digest(path))
            variants.append((name,staged))
        variant('duplicate-isotope','isotopes',isotope_rows,lambda r:r.append(r[1]))
        variant('missing-state-energy','states',state_rows,lambda r:r[1].__setitem__(5,''))
        variant('unknown-energy-source','states',state_rows,lambda r:r[1].__setitem__(6,'UNKNOWN'))
        variant('excited-state-cannot-use-reference-zero','states',state_rows,lambda r:r[1].__setitem__(2,'excited'))
        variant('molecule-cannot-use-atomic-formation-zero','states',state_rows,lambda r:r[1].__setitem__(3,'1H:2'))
        variant('unavailable-state-cannot-contain-invented-zero','states',state_rows,lambda r:r[8].__setitem__(5,'0'))
        variant('fractional-charge','states',state_rows,lambda r:r[1].__setitem__(4,'0.5'))
        variant('unknown-isotope','states',state_rows,lambda r:r[1].__setitem__(3,'999X:1'))
        variant('negative-physical-mass','states',state_rows,lambda r:r[1].__setitem__(11,'-1'))
        variant('duplicate-state','states',state_rows,lambda r:r.append(r[1]))
        for name,staged in variants:
            rejected=execute(arguments(executable,staged),'refusal-'+name)
            check('refusal:'+name,rejected.returncode==2 and
                  rejected.stdout.startswith('LOAD-REFUSED\t'))
        no_physical=copy.deepcopy(state_rows)
        no_physical[1][11:14]=['','','']
        path=output/'explicitly-missing-physical-mass.csv'
        with path.open('w',newline='') as stream:csv.writer(stream).writerows(no_physical)
        staged=copy.deepcopy(profile);staged['states']=dict(path=str(path),sha256=digest(path))
        missing_mass=execute(arguments(executable,staged),'missing-physical-mass')
        check('missing physical mass never falls back to reference mass',
              missing_mass.returncode==0 and 'PHYSICAL\t1H:ground\t' not in missing_mass.stdout)
        check('native sources unchanged during compilation and tests',
              manifest['native_source_sha256']=={str(p):digest(p) for p in identities})
        check('original data bytes and timestamps unchanged',
              before=={str(p):(digest(p),p.stat().st_mtime_ns) for p in data_files})
        failures=sum(t['status']=='FAIL' for t in tests)
        result=dict(status='FAIL' if failures else 'PASS',passed=len(tests)-failures,
                    failed=failures,errors=0,skipped=0)
    except FileNotFoundError as error:
        result=dict(status='SKIPPED',reason=str(error),passed=0,failed=0,errors=0,skipped=1)
    except (OSError,ValueError,KeyError,StopIteration,RuntimeError,subprocess.TimeoutExpired) as error:
        result=dict(status='ERROR',reason=str(error),passed=sum(t['status']=='PASS' for t in tests),
                    failed=sum(t['status']=='FAIL' for t in tests),errors=1,skipped=0)
    # Retain the original scalar totals for direct users while supplying the
    # common four-state summary consumed by the global Moon support runner.
    # These are the SAME evaluated controls, never extra tests or promotions.
    result['counts']={status:result[key] for status,key in
        (('PASS','passed'),('FAIL','failed'),('ERROR','errors'),('SKIPPED','skipped'))}
    result.update(schema='srcMoon-M3-state-profile-support-v1',
        scope='source-bound-profile-and-parser-only',formal_M3_status='SKIPPED / NOT QUALIFIED',
        physical_activation=False,checks=tests,runner=str(Path(__file__).resolve()),
        elapsed_s=time.monotonic()-started,output_dir=str(output))
    (output/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    packet=ROOT/'test_output/srcMoon/m3-approved-products-20261118T163142Z/source-acquisition-v2'
    parser.add_argument('--raw-dir',type=Path,default=packet/'raw')
    parser.add_argument('--reference-inventory',type=Path,default=packet/'reference_mass_inventory.json')
    parser.add_argument('--decision-register',type=Path,default=ROOT/'srcMoon/test/M3_DAUGHTER_AND_ION_SURFACE_DECISIONS.json')
    parser.add_argument('--canonical-root',type=Path,default=Path('/data/vtenishe/moon_validation_data/photochemistry'))
    args=parser.parse_args(argv)
    try:value=run(args)
    except (OSError,ValueError) as error:
        print(json.dumps(dict(status='ERROR',reason=str(error))));return 2
    print(json.dumps(value,indent=2));return {'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[value['status']]


if __name__=='__main__':
    raise SystemExit(main())
