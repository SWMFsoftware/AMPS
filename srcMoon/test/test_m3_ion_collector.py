#!/usr/bin/env python3
"""Compile production P08 collector/parser and compare analytical book budgets.

No observations, source rates, lunar ion retention, or daughter recoil model is
created here. Normal linked PIC callbacks need a separately recorded native run.
Each requested output directory is new; failing artifacts are never replaced.
"""
from __future__ import annotations
import argparse,csv,datetime,hashlib,json,math,re,shutil,subprocess,time,traceback
from pathlib import Path
from chemical_database_build import shared_dependencies
ROOT=Path(__file__).resolve().parents[2]
NAMES=(
 'collect-explicit','one-publication','identity-state','integer-nuclei','actual-two-epochs',
 'frame-normal-retained','energy-provenance','moving-facet-galilean',
 'reject-missing-formation','reject-missing-excitation','reject-negative-excitation',
 'reject-missing-energy-zero','reject-missing-evidence','reject-missing-state','reject-neutral',
 'reject-uncompiled','reject-zero-mass','reject-invalid-nuclei','reject-electron-without-nuclei','reject-zero-weight',
 'reject-negative-weight','reject-nonfinite-weight','reject-missing-charge-unit',
 'reject-missing-impact-epoch','reject-missing-geometry-epoch','reject-nonfinite-position',
 'reject-nonfinite-velocity','reject-nonfinite-surface-velocity','reject-zero-normal',
 'reject-outgoing','reject-missing-geometry','reject-missing-transforms',
 'reject-nonfinite-transform','reject-noninvertible-transform','explicit-zero-excitation',
 'distinct-explicit-states','reject-duplicate-state','reject-mixed-energy-zero','reject-ledger-energy-zero-change',
 'reject-empty-inventory','reject-null-ledger','negative-ion-signed-charge','explicit-zero-formation')
EXPECTED={'nucleus_A':5,'nucleus_B':2.5,'charge_elementary':5,'charge_C':2.5,'mass_kg':10,
 'momentum_x':20,'momentum_y':-20,'momentum_z':-60,'kinetic_facet_J':220,
 'kinetic_lso_J':250,'formation_J':-17.5,'excitation_J':7.5,'total_facet_J':210}

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def change(text,key,value):
 pattern=r'(?m)^\s*'+re.escape(key)+r'\s*=[^\n]*$'
 if re.search(pattern,text):return re.sub(pattern,lambda _:key+' = '+str(value),text)
 return text.replace('#end application',key+' = '+str(value)+'\n#end application')

def run_checks(output):
 output=Path(output).resolve();output.mkdir(parents=True,exist_ok=False)
 started=time.monotonic();result={'schema':'srcMoon-M3-ion-collector-controls-v1','status':'ERROR',
  'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
  'scope':'analytical production-kernel and actual parser support only; NOT lunar ion-surface validation',
  'formal_status':'SKIPPED','checks':[],'linked_callback_execution':False,'observations_used':False}
 sources=[ROOT/'srcMoon'/n for n in ['test/ion_collector/probe.cpp','species/IonCollector.cpp',
   'species/IonCollector.h','species/ReactionExecution.cpp','species/ReactionExecution.h','MoonInput.cpp',
   'MoonInput.h','species/ElectronImpact.cpp','species/Photochemistry.cpp','species/ChargeExchange.cpp']]
 sources += shared_dependencies(ROOT,'runtime')
 sources += [ROOT/'srcMoon/species'/(name+'.h') for name in
             ('ElectronImpact','Photochemistry','ChargeExchange')]
 try:
  result['source_revision']=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
  result['source_hashes']={str(p.relative_to(ROOT)):sha(p) for p in sources}
  compiler=shutil.which('g++')
  if not compiler:
   result.update(status='SKIPPED',reason='native C++ compiler unavailable')
  else:
   exe=output/'ion-collector-probe'
   cpp=[p for p in sources if p.suffix=='.cpp']
   cmd=[compiler,'-std=c++17','-O2','-g','-Wall','-Wextra','-pedantic','-I',str(ROOT),'-I',str(ROOT/'srcMoon'),*(str(p) for p in cpp),'-o',str(exe)]
   p=subprocess.run(cmd,cwd=ROOT,text=True,capture_output=True,timeout=180)
   (output/'compile.stdout.log').write_text(p.stdout);(output/'compile.stderr.log').write_text(p.stderr)
   result['compile']={'command':cmd,'returncode':p.returncode}
   if p.returncode:raise RuntimeError('production collector/parser probe did not compile')
   result['executable']={'path':str(exe),'sha256':sha(exe)}
   p=subprocess.run([str(exe),str(output/'production_budget.csv')],cwd=ROOT,text=True,capture_output=True,timeout=30)
   (output/'stdout.log').write_text(p.stdout);(output/'stderr.log').write_text(p.stderr)
   rows=[line.split() for line in p.stdout.splitlines() if line.startswith('CHECK ')]
   names=[r[1] for r in rows if len(r)==3]
   if len(rows)!=len(NAMES) or len(names)!=len(NAMES) or len(set(names))!=len(names) or set(names)!=set(NAMES):
    raise RuntimeError('native collector controls missing/duplicate/malformed; zero checks never imply PASS')
   if p.returncode not in (0,1):raise RuntimeError('native collector crashed or exited outside its check contract')
   failures=sum(r[2]!='1' for r in rows)
   if p.returncode!=int(bool(failures)) or p.stdout.count('SUMMARY failed='+str(failures)+'\n')!=1:
    raise RuntimeError('native collector summary/exit does not match its fixed named controls')
   for r in rows:result['checks'].append({'name':r[1],'status':'PASS' if r[2]=='1' else 'FAIL'})
   raw=list(csv.DictReader((output/'production_budget.csv').open()))
   if len(raw)!=len(EXPECTED) or {r['quantity'] for r in raw}!=set(EXPECTED):raise RuntimeError('native budget artifact has missing/duplicate quantities')
   comparisons=[]
   for row in raw:
    actual=float(row['actual']);expected=EXPECTED[row['quantity']]
    # The fixture values are exactly representable binary fractions/integers.
    # Exact equality therefore tests the number and energy bookkeeping without
    # selecting a new physical tolerance or hiding the signed formation term.
    ok=math.isfinite(actual) and actual==expected
    result['checks'].append({'name':'budget-'+row['quantity'],'status':'PASS' if ok else 'FAIL'})
    comparisons.append({'quantity':row['quantity'],'actual':actual,'independent_expected':expected,'difference':actual-expected})
   with (output/'comparison.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(comparisons[0]));w.writeheader();w.writerows(comparisons)
   try:
    import matplotlib
   except ImportError:
    result['visualization']={'status':'SKIPPED','reason':'optional matplotlib unavailable; analytical CSV retained'}
   else:
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(10,4))
    ax.bar([r['quantity'] for r in comparisons],[abs(r['difference']/r['independent_expected']) for r in comparisons])
    ax.set_ylabel('Absolute relative difference (dimensionless)')
    ax.set_title('Analytical ideal-collector bookkeeping — not lunar observations')
    ax.tick_params(axis='x',rotation=70);fig.tight_layout()
    products=[]
    for suffix in ('png','eps'):
     path=output/('comparison.'+suffix);fig.savefig(path);products.append({'path':str(path),'sha256':sha(path)})
    plt.close(fig);result['visualization']={'status':'PASS','products':products}
   templates={'legacy':ROOT/'srcMoon/examples/m1_regression.in','first':ROOT/'srcMoon/examples/m3_common_time_first_order.in','second':ROOT/'srcMoon/examples/m3_common_time_second_order.in'}
   variants=[('default-unsupported','legacy',None,0,'unsupported'),('explicit-unsupported','legacy','unsupported',0,'unsupported'),
    ('reject-unreviewed-boundary','legacy','neutral-sticking',3,None),('reject-collector-legacy','legacy','ideal-collector',3,None),
    ('collector-first-explicit','first','ideal-collector',0,'ideal-collector'),('collector-second-explicit','second','ideal-collector',0,'ideal-collector')]
   for name,which,selector,code,expected in variants:
    directory=output/name;directory.mkdir()
    text=templates[which].read_text()
    for key in ('surface_cea_file','surface_tecplot_file'):text=change(text,key,str(directory/(key+'.unused')))
    if selector is not None:text=change(text,'ion_surface_boundary',selector)
    application=directory/'application.in';application.write_text(text)
    process=subprocess.run([str(exe),'parse',str(application)],cwd=ROOT,text=True,capture_output=True,timeout=30)
    (directory/'stdout.log').write_text(process.stdout);(directory/'stderr.log').write_text(process.stderr)
    values=dict(line.split('=',1) for line in process.stdout.splitlines() if '=' in line)
    # A negative parser case must fail for THIS selector, not incidentally
    # because this machine lacks the examples' external SPICE/LOLA paths.
    missing_paths='spice_path must contain' in process.stderr or 'lola_image_file and lola_label_file' in process.stderr
    ok=process.returncode==code and (values.get('ion_surface_boundary')==expected if code==0 else values.get('preserved')=='1' and 'ion_surface_boundary' in process.stderr)
    result['checks'].append({'name':name,'status':'SKIPPED' if missing_paths else 'PASS' if ok else 'ERROR' if process.returncode<0 else 'FAIL','returncode':process.returncode,'expected_returncode':code,'reason':'example path prerequisites absent; no parser result inferred' if missing_paths else ''})
   after={str(p.relative_to(ROOT)):sha(p) for p in sources}
   if after!=result['source_hashes']:raise RuntimeError('production/test source changed during controls; retain this attempt and rerun a matched snapshot')
   result['status']='ERROR' if any(c['status']=='ERROR' for c in result['checks']) else 'FAIL' if any(c['status']=='FAIL' for c in result['checks']) else 'SKIPPED' if any(c['status']=='SKIPPED' for c in result['checks']) else 'PASS'
 except Exception as exc:
  result.update(status='ERROR',reason=str(exc));(output/'infrastructure-error.log').write_text(traceback.format_exc())
 result['elapsed_seconds']=time.monotonic()-started
 result['finished_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
 result['counts']={s:sum(c['status']==s for c in result['checks']) for s in ['PASS','FAIL','SKIPPED','ERROR']}
 (output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({'status':result['status'],'counts':result['counts'],'result':str(output/'result.json')}))
 return result
run=run_checks
if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output-dir',type=Path,required=True);args=parser.parse_args()
 result=run_checks(args.output_dir);raise SystemExit(0 if result['status'] in ['PASS','SKIPPED'] else 1)
