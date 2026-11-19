#!/usr/bin/env python3
"""Exact-Fraction controls of production compensated variable clock semantics."""
from __future__ import annotations
import argparse
import csv
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
NAMES=set('''origin ten_thousand_realistic_no_event_outer_intervals
dyadic_ordinary_debits_exact_in_native_pair equivalent_dyadic_partition_same_full_clock
near_endpoint_pair_constructed legacy_exact_scalar_query_still_refuses_loss
full_components_no_early_event_at_predecessor retained_positive_tail_has_its_real_event
tail_event_exact_time full_pair_first_representable_event_not_projected_early
legacy_scalar_integrated_api_retains_refusal full_double_exponent_span_pair_constructed
cache_detects_low_component_hidden_in_scalar cache_identical_full_pair
full_exponent_span_no_event full_exponent_span_tail_event non_dyadic_integral_supported
non_dyadic_storage_bound_explicit invalid_component_pair_transactional
nonfinite_component_pair_transactional positive_native_storage_underflow_refused
unproven_default_receipt_refused real_receipt_available receipt_matches_all_origin_bits
instantaneous_CDF_retained query_and_commit_do_not_mutate_origin
consumed_components_query_transactional'''.split())
FILES=[ROOT/'srcMoon/species'/name for name in ('ReactionExecution.h','ReactionExecution.cpp',
 'PiecewiseHazard.h','PiecewiseHazard.cpp','ReactionExposure.h','ReactionExposure.cpp')]+[
 HERE/'probe.cpp',HERE/'run.py',HERE/'CONTRACT.md',HERE/'README.md']


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def binary(text):
    """Exact C++ hex-float conversion, including LONG DOUBLE and subnormals.

    float.fromhex would first discard long-double low bits. Integer numerator
    and power-of-two denominator preserve the entire printed binary input.
    """
    match=re.fullmatch(r'(-?)0x([0-9a-fA-F]+)(?:\.([0-9a-fA-F]*))?p([+-]?\d+)',text)
    if not match:raise ValueError('malformed binary diagnostic: '+text)
    sign,integer,fraction,exponent=match.groups();fraction=fraction or ''
    numerator=int(integer+fraction,16)*(-1 if sign else 1)
    power=int(exponent)-4*len(fraction)
    return Fraction(numerator*(2**power),1)if power>=0 else Fraction(numerator,2**(-power))


def run(output,cxx='g++'):
    started=time.monotonic();output=Path(output).resolve();output.mkdir(parents=True,exist_ok=False)
    result=dict(schema='srcMoon-compensated-depth-support-v1',status='ERROR',
        scope='production numerical clocks only; no physical chemistry qualification',
        formal_I21_I24_status='SKIPPED / NOT QUALIFIED',counts={})
    manifest=dict(commands=[],external_data=[],MPI='none',random_seed=None,
                  stochastic_controls='explicit uniforms; no RNG')
    try:
        manifest['source_hashes']={str(p):sha(p)for p in FILES}
        manifest['source_revision']=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
        compiler=shutil.which(cxx)
        if not compiler:raise FileNotFoundError('C++17 compiler unavailable')
        os.environ.setdefault('MPLCONFIGDIR',str(output/'matplotlib-cache'))
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
        except ImportError as error:raise FileNotFoundError('matplotlib required for declared diagnostics')from error
        manifest['compiler']=subprocess.check_output([compiler,'--version'],text=True,timeout=15)
        executable=output/'compensated-depth-probe'
        command=[compiler,'-std=c++17','-O2','-Wall','-Wextra','-Werror','-pedantic','-I',str(ROOT)]+[
            str(ROOT/'srcMoon/species'/name)for name in('ReactionExecution.cpp','PiecewiseHazard.cpp','ReactionExposure.cpp')]+[
            str(HERE/'probe.cpp'),'-o',str(executable)]
        manifest['commands'].append(dict(argv=command,cwd=str(ROOT),log=str(output/'compile.log')))
        with(output/'compile.log').open('w')as stream:
            compiled=subprocess.run(command,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT,timeout=120)
        if compiled.returncode:raise RuntimeError('production compilation failed; see compile.log')
        manifest['executable']=dict(path=str(executable),sha256=sha(executable))
        manifest['commands'].append(dict(argv=[str(executable)],cwd=str(ROOT),log=str(output/'stdout.log')))
        native=subprocess.run([str(executable)],cwd=ROOT,text=True,capture_output=True,timeout=180)
        (output/'stdout.log').write_text(native.stdout);(output/'stderr.log').write_text(native.stderr)
        rows=[dict(status=status,name=name)for status,name in re.findall(r'^(PASS|FAIL) (.+)$',native.stdout,re.M)]
        summaries=re.findall(r'^SUMMARY checks=(\d+) failed=(\d+)$',native.stdout,re.M)
        if len(summaries)!=1 or {r['name']for r in rows}!=NAMES or len(rows)!=len(NAMES):
            raise ValueError('incomplete/duplicate native compensated-depth receipt')
        count,failed=map(int,summaries[0])
        if count!=len(rows)or failed!=sum(r['status']=='FAIL'for r in rows)or native.returncode!=int(failed>0):
            raise ValueError('native count/status/exit contradiction')
        source=re.findall(r'^INPUT rate (\S+) predecessor (\S+) gap (\S+)$',native.stdout,re.M)
        if len(source)!=1:raise ValueError('missing binary source diagnostics')
        rate,predecessor,gap=map(binary,source[0])
        metrics=re.findall(r'^CLOCK (\S+) (\S+) (\S+) (\S+) (\S+) (\S+)$',native.stdout,re.M)
        expected={'ordinary10000','near_endpoint','full_exponent_span','third'}
        if len(metrics)!=len(expected)or {r[0]for r in metrics}!=expected:
            raise ValueError('missing/duplicate Fraction-comparison metric')
        comparison=[]
        for name,*numbers in metrics:
            ah,al,bh,bl,bound=map(binary,numbers);before,after=ah+al,bh+bl
            debit={'ordinary10000':10000*rate,'near_endpoint':predecessor,
                   'full_exponent_span':predecessor,'third':Fraction(1,3)}[name]
            exact=before-debit;difference=abs(after-exact)
            passed=exact>0 and bound>=0 and difference<=bound
            if name!='third':passed=passed and difference==0 and bound==0
            rows.append(dict(name='independent_Fraction_'+name,status='PASS'if passed else'FAIL'))
            comparison.append([name,str(exact),str(after),str(difference),str(bound),'PASS'if passed else'FAIL'])
        unchanged={str(p):sha(p)for p in FILES}==manifest['source_hashes']
        rows.append(dict(name='production_and_test_source_unchanged',status='PASS'if unchanged else'FAIL'))
        with(output/'comparison.csv').open('w',newline='')as stream:
            writer=csv.writer(stream);writer.writerow(['case','exact_rational_residual','native_pair_rational_residual','absolute_storage_error','rigorous_storage_error_bound','status']);writer.writerows(comparison)
        fig,ax=plt.subplots();indices=list(range(len(comparison)))
        ax.plot(indices,[float(Fraction(r[1]))for r in comparison],'o',label='Independent exact Fraction')
        ax.plot(indices,[float(Fraction(r[2]))for r in comparison],'x',label='Compiled native pair')
        ax.set_xticks(indices);ax.set_xticklabels([r[0]for r in comparison],rotation=20)
        ax.set(ylabel='Pending dimensionless optical depth',title='Mathematical clock controls — no physical validation');ax.legend()
        for suffix in('png','eps'):fig.savefig(output/('comparison.'+suffix),bbox_inches='tight')
        plt.close(fig)
        counts={name:sum(r['status']==name for r in rows)for name in('PASS','FAIL','ERROR','SKIPPED')}
        result.update(status='FAIL'if counts['FAIL']else'PASS',counts=counts,results=rows)
    except FileNotFoundError as error:
        result.update(status='SKIPPED',reason=str(error),counts={'PASS':0,'FAIL':0,'ERROR':0,'SKIPPED':1})
    except(OSError,ValueError,RuntimeError,subprocess.SubprocessError)as error:
        result.update(status='ERROR',reason=str(error),counts={'PASS':0,'FAIL':0,'ERROR':1,'SKIPPED':0})
    result.update(elapsed_s=time.monotonic()-started,runner=str(HERE/'run.py'),output_dir=str(output))
    (output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    (output/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--cxx',default='g++');args=parser.parse_args()
    value=run(args.output_dir,args.cxx);print(json.dumps(value,indent=2))
    return {'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[value['status']]


if __name__=='__main__':raise SystemExit(main())
