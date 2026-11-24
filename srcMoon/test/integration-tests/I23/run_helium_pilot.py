#!/usr/bin/env python3
"""Check the source-bound stationary 4He provider before event admission."""
import argparse
import csv
from decimal import Decimal, localcontext
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
TEST = HERE.parents[1]
ROOT = TEST.parents[1]
sys.path.insert(0, str(ROOT/'utility/chemical_database/state_profiles'))
from run import arguments
sys.path.insert(0, str(TEST/'U16_electron_impact'))
from voronov_verification import (MissingPrerequisite, independent_rows,
                                 preflight, reference, reference_log_slope)
from chemical_database_build import shared_dependencies, shared_sources
from m3_qualification import sha, strict_object

EXIT = {'PASS': 0, 'FAIL': 1, 'ERROR': 2, 'SKIPPED': 77}
KIM_SOURCE = 'https://nvlpubs.nist.gov/nistpubs/jres/097/jresv97n6p689_A1b.pdf'
KIM_SHA256 = 'b777277069851715e5810c17f36ad8fcbaf005fd6c0071c2b9baba28082be961'

def kim_model_comparison(pdf, profile, reference_coefficient, reference_mean_j, manifest):
    """Integrate the published total-cross-section fit as a second model.

    This measures fit-to-fit discrepancy within the paper's finite support.
    Its authors used laboratory data in the fit; the quadrature is not a
    state-resolved observation or an uncertainty bound for event admission.
    """
    from scipy.integrate import quad
    if not pdf.is_file() or shutil.which('pdftotext') is None:
        raise MissingPrerequisite('optional NIST source PDF or pdftotext unavailable')
    if sha(pdf)!=KIM_SHA256:
        raise ValueError('optional NIST Kim 1992 source SHA-256 mismatch')
    manifest['input_sha256'][str(pdf.resolve())]=KIM_SHA256
    command=['pdftotext','-layout','-f','2','-l','2',str(pdf),'-']
    native=subprocess.run(command,text=True,capture_output=True,check=False,timeout=30)
    manifest['commands'].append({'argv':command,'returncode':native.returncode})
    if native.returncode:
        raise RuntimeError('could not read the original NIST fit table')
    pattern=r'^\s+He\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(-\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*$'
    rows=re.findall(pattern,native.stdout,re.MULTILINE)
    if len(rows)!=1:
        raise ValueError('expected one neutral-He row in the pinned native Table 1')
    fit_threshold,a,b,c,d=map(float,rows[0])
    energy_max_ev,temp_ev=10_000,10
    with Path(profile['constants']['path']).open() as stream:
        constants={r['name']:r['value'] for r in csv.DictReader(stream)}
    electron_mass=float(constants['electron mass'])
    elementary_charge=float(constants['elementary charge'])

    def moment_integrand(z,order):
        # Dimensionless z=(E-I)/T. Eq. (2) uses dimensionless x and the
        # tabulated D in x+D; convert its cross section from 1e-16 cm² to m².
        energy=fit_threshold+temp_ev*z
        x=energy/fit_threshold
        log_x=math.log1p(temp_ev*z/fit_threshold)
        cross_section=1e-20*(a*log_x+(b*log_x+c*(x-1))/(x+d))/x
        if cross_section<0 or not math.isfinite(cross_section):
            raise ArithmeticError('published fit evaluation is outside its stated positive domain')
        return temp_ev*cross_section*energy**(order+1)*math.exp(-energy/temp_ev)

    domain_end=(energy_max_ev-fit_threshold)/temp_ev
    # Split the dimensionless support so a long cold tail cannot hide the
    # near-threshold integral from the adaptive quadrature.
    cuts=(0,5,25,100,domain_end)
    moments=[]
    for order in (0,1):
        pieces=[quad(moment_integrand,lo,hi,args=(order,),epsabs=1e-32,epsrel=1e-11)
                for lo,hi in zip(cuts,cuts[1:])]
        moments.append([math.fsum(piece[i] for piece in pieces) for i in (0,1)])
    factor=math.sqrt(8*elementary_charge/(math.pi*electron_mass))/temp_ev**1.5
    result_coefficient=factor*moments[0][0]
    result_mean_ev=moments[1][0]/moments[0][0]
    result_mean_j=result_mean_ev*elementary_charge
    # Peterson et al. report discrete incident electron energies within
    # [100, 2000] eV. Only calculate the fraction that could be compared by
    # energy. Neither continuous angular coverage nor a ground-ion fraction
    # follows from these limits. Ionization weighting includes v*sigma(E),
    # so integrating the Maxwellian number distribution alone would be wrong.
    window_low,window_high=100,2000
    boundaries=(fit_threshold,window_low,window_high,energy_max_ev)
    bands=[]
    for lower,upper in zip(boundaries,boundaries[1:]):
        z_low,z_high=(lower-fit_threshold)/temp_ev,(upper-fit_threshold)/temp_ev
        # Subdivide the long exponential tail in dimensionless units, using
        # relative error for even the very small above-window contribution.
        # The whole-domain integration above uses different quadrature panels.
        panels=[z_low,*(z_low+s for s in (5,25,100) if z_low+s<z_high),z_high]
        pieces=[quad(moment_integrand,lo,hi,args=(0,),epsabs=0,epsrel=1e-11)
                for lo,hi in zip(panels,panels[1:])]
        integral,error=(math.fsum(piece[i] for piece in pieces) for i in (0,1))
        bands.append({'incident_interval_eV':[lower,upper],'fraction':integral/moments[0][0],
                      'raw_integral':integral,'estimated_absolute_quadrature_error':error})
    residual=abs(math.fsum(band['raw_integral'] for band in bands)-moments[0][0])
    error_budget=moments[0][1]+math.fsum(band['estimated_absolute_quadrature_error'] for band in bands)
    return {
        'scope':'finite-support fit-to-fit diagnostic only; no activation or source-error bound',
        'source_url':KIM_SOURCE,'source_sha256':KIM_SHA256,
        'locator':'printed p. 690, Table 1 neutral-He row and equation (2)',
        'source_coefficients':dict(zip(('threshold_eV','A','B','C','D'),rows[0])),
        'thermal_energy_eV':temp_ev,'integration_range_eV':[fit_threshold,energy_max_ev],
        'laboratory_incident_energy_coverage':{
            'study':'https://doi.org/10.1103/PhysRevA.5.712','locator':'Abstract, incident beam energy range',
            'source_interval_eV':[window_low,window_high],'bands':bands,
            'fraction_in_studied_energy_interval':bands[1]['fraction'],
            'absolute_integral_closure_error':residual,'estimated_numerical_error_budget':error_budget,
            'interpretation':'Conditional on ionization within the Kim fit domain. The collision fraction counts all incident energies between the experimental endpoints. It is an upper bound on coverage by energy in this model, not a measured fraction or a model error bound.'},
        'raw_moments_and_estimated_absolute_quadrature_errors':moments,
        'kim_rate_coefficient_m3_s':result_coefficient,
        'voronov_rate_coefficient_m3_s':float(reference_coefficient),
        'signed_relative_coefficient_difference':result_coefficient/float(reference_coefficient)-1,
        'kim_conditional_mean_eV':result_mean_ev,
        'voronov_conditional_mean_eV':float(reference_mean_j)/elementary_charge,
        'signed_mean_difference_J':result_mean_j-float(reference_mean_j),
        'analytic_slope_equivalent_difference':(result_mean_j-float(reference_mean_j))/(elementary_charge*temp_ev)}


def run(profile_path, expected_sha, output, data_root, kim_pdf=None):
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    manifest = {'source_revision': subprocess.check_output(
        ['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'MPI': 'none', 'OpenMP': 'none', 'seed': None, 'commands': []}
    checks = []
    try:
        if not profile_path.is_file():
            raise MissingPrerequisite('prepared ground-ion state profile unavailable')
        if sha(profile_path) != expected_sha:
            raise ValueError('state profile digest mismatch')
        profile = strict_object(profile_path)
        shutil.copyfile(profile_path, output/'profile.snapshot.json')
        package, hashes = preflight(data_root)
        table = package/'raw/electron_impact/voronov/cfit.dat'
        manifest['input_sha256'] = {str(profile_path.resolve()): expected_sha, **hashes}
        for descriptor in [profile[name] for name in ('isotopes','constants','states')]+profile['sources']:
            native = Path(descriptor['path'])
            if not native.is_file():
                raise MissingPrerequisite('profile source missing: '+str(native))
            if sha(native) != descriptor['sha256']:
                raise ValueError('profile source changed: '+str(native))
            manifest['input_sha256'][str(native)] = descriptor['sha256']
        compiler = shutil.which('g++')
        if compiler is None:
            raise MissingPrerequisite('C++17 compiler unavailable')

        def invoke(argv, label, timeout):
            with (output/(label+'.stdout.log')).open('x') as stdout, (output/(label+'.stderr.log')).open('x') as stderr:
                process = subprocess.run(argv, cwd=ROOT, timeout=timeout, stdout=stdout, stderr=stderr)
            manifest['commands'].append({'argv': argv, 'returncode': process.returncode})
            return process.returncode

        executable = output/'helium-pilot'
        species = ROOT/'srcMoon/species'
        components = ('electron-impact','photochemistry')
        sources = [HERE/'helium_pilot_probe.cpp',
            *[species/(name+'.cpp') for name in ('HeliumElectronImpact','PhysicalReactionEvents',
               'PhysicalStateProfile','ReducedKinematics','IonCollector','ReactionExecution','PiecewiseHazard')],
            ROOT/'srcMoon/MoonTransportClock.cpp',
            *shared_sources(ROOT,components)]
        manifest['source_sha256'] = {str(p):sha(p) for p in (
            Path(__file__).resolve(), HERE/'HELIUM_PILOT_CONTRACT.md',
            *sources,*[p for p in species.glob('*.h')],
            *shared_dependencies(ROOT,components),
            ROOT/'src/models/chemical_database/ReservoirProducts.h',
            ROOT/'src/models/chemical_database/ReducedModels.h')}
        argv = [compiler,'-std=c++17','-O2','-Wall','-Wextra','-Werror','-pedantic',
                '-I',str(ROOT),*[str(source) for source in sources],'-o',str(executable)]
        if invoke(argv,'compile',120):
            raise RuntimeError('helium provider build failed; see compile.stderr.log')
        manifest['executable_sha256']=sha(executable)
        argv = arguments(executable,profile)
        argv[1:1] = [str(table),hashes[str(table)]]
        status = invoke(argv,'native',30)
        if status == 77:
            raise MissingPrerequisite((output/'native.stderr.log').read_text().strip())
        if status not in (0,1):
            raise RuntimeError('helium provider probe failed; see native.stderr.log')
        rows=[line.split('\t') for line in (output/'native.stdout.log').read_text().splitlines()]
        checks=[{'check':row[1],'status':row[2]} for row in rows if row[0]=='CHECK']
        guard_names={'drift-refusal','moving-target-refusal','invalid-density-preserves-output',
            'unreviewed-temperature-refusal','missing-evidence-preserves-composition',
            'one-incoming-two-outgoing-bath-electrons'}
        if (len(checks)!=len(guard_names) or {row['check'] for row in checks}!=guard_names
                or any(row['status'] not in ('PASS','FAIL') for row in checks)
                or status!=int(any(row['status']=='FAIL' for row in checks))):
            raise ValueError('native guard records contradict the required checks or exit')

        with localcontext() as context:
            context.prec=70
            native_fit = independent_rows(package)['He']
            energy_ev = Decimal(10)
            coefficient = reference(native_fit,energy_ev)
            slope = reference_log_slope(native_fit,energy_ev)
            mean = energy_ev*Decimal('1.602176634e-19')*(Decimal('1.5')+slope)
            with Path(profile['states']['path']).open() as stream:
                ion = [r for r in csv.DictReader(stream) if r['state_id']=='4He+:ground']
            if len(ion)!=1:
                raise ValueError('unique physical ground-ion state required')
            threshold = Decimal(ion[0]['formation_J'])
            data = [list(map(Decimal,row[1:])) for row in rows if row[0]=='DATA']
            if [r[0] for r in data] != [0,10**6,10**7]:
                raise ValueError('native pilot omitted or reordered prescribed densities')
            for density,*values in data:
                expected=[coefficient,density*coefficient,slope,mean,native_fit['energy'],threshold]
                names=('coefficient_m3_s','rate_s','log_slope','mean_incoming_energy_j',
                       'fit_threshold_ev','physical_threshold_j')
                for name,actual,wanted in zip(names,values,expected):
                    relative=abs(actual/wanted-1) if wanted else abs(actual)
                    limit=Decimal('1e-12') if wanted else Decimal(0)
                    checks.append({'check':f'{name}@{density}/m3','actual':str(actual),
                        'reference':str(wanted),'relative_or_zero_error':str(relative),
                        'limit':str(limit),'status':'PASS' if relative<=limit else 'FAIL'})
            with (output/'comparison.csv').open('x',newline='') as stream:
                fields=['check','actual','reference','relative_or_zero_error','limit','status']
                writer=csv.DictWriter(stream,fieldnames=fields)
                writer.writeheader()
                writer.writerows(row for row in checks if 'actual' in row)
        if kim_pdf is not None:
            diagnostic=kim_model_comparison(kim_pdf,profile,coefficient,mean,manifest)
            (output/'kim_model_comparison.json').write_text(json.dumps(diagnostic,indent=2,allow_nan=False)+'\n')
            coverage=diagnostic['laboratory_incident_energy_coverage']
            closed=coverage['absolute_integral_closure_error']<=coverage['estimated_numerical_error_budget']
            valid=all(0<=band['fraction']<=1 for band in coverage['bands'])
            checks.append({'check':'Kim-incident-energy-window-partition','status':'PASS' if closed and valid else 'FAIL',
                           'absolute_integral_closure_error':coverage['absolute_integral_closure_error'],
                           'estimated_numerical_error_budget':coverage['estimated_numerical_error_budget']})
        result_status='FAIL' if any(row['status']=='FAIL' for row in checks) else 'PASS'
    except (MissingPrerequisite,ImportError) as error:
        result_status='SKIPPED'
        checks.append({'check':'prerequisite','status':result_status,'reason':str(error)})
    except (OSError,ValueError,KeyError,RuntimeError,ArithmeticError,subprocess.TimeoutExpired) as error:
        result_status='ERROR'
        checks.append({'check':'infrastructure','status':result_status,'reason':str(error)})
    result={'status':result_status,'scope':'I23 stationary helium rate, model moment and admission guards',
        'linked_status':'SKIPPED','physical_activation':False,'checks':checks,
        'counts':{s:sum(row['status']==s for row in checks) for s in EXIT}}
    for name,payload in [('run_manifest.json',manifest),('result.json',result)]:
        (output/name).write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile',type=Path,required=True)
    parser.add_argument('--profile-sha256',required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--data-root',type=Path,default=Path('/data/vtenishe/moon_validation_data'))
    parser.add_argument('--kim-1992-pdf',type=Path,
                        help='optional exact NIST paper for the independent finite-support fit comparison')
    args=parser.parse_args()
    result=run(args.profile,args.profile_sha256,args.output_dir,args.data_root,args.kim_1992_pdf)
    print(json.dumps(result,indent=2))
    raise SystemExit(EXIT[result['status']])
