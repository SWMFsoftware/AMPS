#!/usr/bin/env python3
"""Read-only real-release QA against acquired originals; run explicitly.

This is not unittest discovery and not a physical-validation campaign. It opens
an explicitly pinned real database lock through the maintained API, then checks
lexical preservation against the original source packets/native bytes. Missing
original evidence is SKIPPED; malformed/corrupt infrastructure is ERROR; an
available, evaluated preservation mismatch is FAIL. No chemistry law is chosen.
"""
from __future__ import annotations
import argparse,bz2,csv,hashlib,importlib.util,io,json,subprocess,sys,traceback,zipfile
from datetime import datetime,timezone
from decimal import Decimal
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'utility/chemical_database/python'))
from amps_chemical_database import Database,DatabaseError

class MissingOriginal(Exception):pass

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def plain(value):
 if hasattr(value,'items'):return {key:plain(item)for key,item in value.items()}
 if isinstance(value,(tuple,list)):return [plain(item)for item in value]
 return value

def same(actual,expected,message):
 if actual!=expected:raise AssertionError(message)

def csv_data(raw):
 """Decode source CSV without numeric conversion, reordering or gap filling."""
 values=list(csv.reader(io.StringIO(raw.decode('utf-8'))))
 if not values:raise ValueError('Native projection is empty')
 return tuple(values[0]),tuple(tuple(row)for row in values[1:])

class Audit:
 """Retain real input identities, typed outcomes and native command receipts."""
 def __init__(self,database,output,native_probe=None):
  self.db=database;self.output=output;self.rows=[];self.inputs={};self.log=[];self.native_probe=Path(native_probe).resolve()if native_probe else None;self.native_checked=set()
  if self.native_probe:
   if not self.native_probe.is_file():raise MissingOriginal('Native catalogue probe unavailable')
   self.inputs[str(self.native_probe)]=sha(self.native_probe);(output/'native-queries').mkdir()
 def original(self,path,expected):
  """Read an authoritative oracle only after its declared SHA agrees.

  Missing originals make the requested comparison unavailable; changed bytes
  are corruption. Neither condition authorizes a replacement expected value.
  """
  path=Path(path)
  if not path.is_file():raise MissingOriginal('Original source unavailable: '+str(path))
  actual=sha(path)
  if actual!=expected:raise ValueError('Original source SHA mismatch: '+str(path))
  self.inputs[str(path)]=actual
  return path.read_bytes()
 def check(self,name,function):
  """Distinguish evaluated mismatch, missing evidence and infrastructure.

  AssertionError means an available preservation criterion missed (FAIL).
  MissingOriginal is SKIPPED. Parser/hash/runtime faults are ERROR. All actual
  successful subchecks remain visible if a later required check cannot run.
  """
  try:metrics=function()or{};row={'name':name,'status':'PASS','metrics':metrics}
  except MissingOriginal as e:row={'name':name,'status':'SKIPPED','reason':str(e)}
  except AssertionError as e:row={'name':name,'status':'FAIL','reason':str(e)}
  except Exception as e:
   row={'name':name,'status':'ERROR','reason':str(e)};self.log.append(name+'\n'+traceback.format_exc())
  self.rows.append(row);self.log.append(row['status']+' '+name+' '+json.dumps(plain(row.get('metrics',row.get('reason'))),sort_keys=True));print(row['status']+' '+name,flush=True);return row
 def packet(self,dataset):
  """Authenticate original packet JSON; never discover a more convenient version."""
  pin=dataset.get('source_packet')
  if not pin:raise MissingOriginal('Explicit original packet unavailable for '+dataset['id'])
  path=Path(pin['path']);raw=self.original(path,pin['sha256']);value=json.loads(raw)
  if value.get('schema')!='amps-chemical-source-packet-v1' or value.get('physical_activation')is not False:raise ValueError('Malformed/activating original packet')
  return path.parent,value,{row['id']:row for row in value['sources']}
 def native(self,base,row):return self.original(base/row['path'],row['sha256'])
 def stored_equals(self,source_id,expected):
  path=self.db.source_path(source_id);same(path.read_bytes(),expected,'Stored source bytes differ: '+source_id)
 def table_equals(self,table_id,columns,rows):
  """Compare every original token/column/order through the actual public API.

  Numbers remain strings in native units. Exact identity is appropriate for
  lossless source copying, not an invented physical-agreement tolerance. The
  optional C++ adapter reads the SAME pinned runtime export. Its TableRows
  returns data rows without a header; successful native loading also checks
  the authenticated native CSV header against the runtime column declaration.
  """
  actual=self.db.table(table_id);same(actual.columns,tuple(columns),'Column identity/order changed: '+table_id)
  same(actual.rows,tuple(tuple(row)for row in rows),'Original lexical values/order changed: '+table_id)
  if self.native_probe and table_id not in self.native_checked:
   # Invoke the existing real Catalogue adapter; no new native implementation
   # or formula is compiled. Every command/output is retained independently.
   identity=self.db.identity;command=[str(self.native_probe),str(self.db.root),identity['export_path'],identity['export_sha256'],identity['release_id'],identity['release_manifest_sha256'],identity['profile_id'],identity['profile_sha256'],'table-access','table',table_id]
   stem=self.output/'native-queries'/hashlib.sha256(table_id.encode()).hexdigest()[:16]
   process=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=30)
   stem.with_suffix('.stdout.json').write_bytes(process.stdout);stem.with_suffix('.stderr.log').write_bytes(process.stderr)
   receipt={'command':command,'table_id':table_id,'executable_sha256':self.inputs[str(self.native_probe)],'returncode':process.returncode,'stdout_sha256':hashlib.sha256(process.stdout).hexdigest(),'stderr_sha256':hashlib.sha256(process.stderr).hexdigest()}
   stem.with_suffix('.command.json').write_text(json.dumps(receipt,indent=2)+'\n')
   if process.returncode!=0:raise RuntimeError('Native catalogue probe process failed: '+table_id)
   native=json.loads(process.stdout)
   same(native.get('open_ok'),True,'Native Catalogue refused authenticated release: '+table_id)
   same(native.get('query_ok'),True,'Native Catalogue refused allowed actual table: '+table_id)
   same(native.get('rows'),[list(row)for row in actual.rows],'Native/Python/original lexical rows differ: '+table_id)
   self.native_checked.add(table_id)
  return actual

def packet_sources(a,dataset):
 """Preserve source bytes, declared domains/units and policy independently.

 Held originals are authenticated as research evidence but are never copied
 or read through a managed numerical-payload API by this check. Their denial
 is exercised separately. Metadata preservation grants no physical meaning.
 """
 base,packet,source=a.packet(dataset);count=0;held=0
 canonical={row['id']:row for row in dataset['files']}
 for row in packet['sources']:
  raw=a.native(base,row);current=canonical[row['id']]
  same(current.get('physical_activation'),False,'Source unexpectedly activates physics')
  same(plain(current.get('quantity_metadata')),row.get('quantity_metadata'),'Native quantity/domain metadata changed: '+row['id'])
  if current['storage_permission']=='allowed':a.stored_equals(row['id'],raw);count+=1
  else:held+=1
 return {'packet_id':packet['id'],'verified_available_native_or_derived_sources':count,'held_sources':held}

def collisions(a,dataset):
 """Use the native HTML coefficient table as oracle for all21 actual fits.

 The maintained lexical reader is bound to the packet's exact recipe SHA.
 Fit domains/reaction labels remain the source strings, including absent
 limits. Parameter dimensions are expression-specific; no blanket eV, cm2
 or dimensionless interpretation, frame choice or rate evaluation is made.
 """
 base,packet,sources=a.packet(dataset)
 spec=importlib.util.spec_from_file_location('real_collision_native',ROOT/'utility/chemical_database/native_collision_sources.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 tables=[row for row in packet['sources']if row.get('table')]
 same(len(tables),21,'IAEA acquired table membership changed');total=0
 canonical={row['id']:row for row in dataset['files']}
 for row in tables:
  def one(row=row):
   parent=sources[row['derived_from'][0]];raw=a.native(base,parent);ident=parent['quantity_metadata']['reaction_as_rendered']
   parsed=module.parse_iaea_fit_page(raw.decode('utf-8'));same(parsed['native_metadata']['reaction_as_rendered'],ident,'IAEA reaction source mismatch')
   a.original(base/row['recipe']['path'],row['recipe']['sha256'])
   same(sha(ROOT/'utility/chemical_database/native_collision_sources.py'),row['recipe']['sha256'],'Maintained collision recipe differs from frozen acquired recipe')
   a.table_equals(row['table']['id'],('parameter_label','value_native'),parsed['parameters'])
   original=parent['quantity_metadata'];stored=plain(canonical[parent['id']]['quantity_metadata'])
   same(stored,original,'IAEA source state/domain metadata changed')
   same(stored.get('x_range'),parsed['native_metadata']['x_range'],'IAEA native domain changed')
   return {'rows':len(parsed['parameters']),'domain_as_published':stored.get('x_range'),'reaction_as_published':ident,'qualified_evaluator':False}
  result=a.check('IAEA-native-fit:'+row['table']['id'],one)
  if result['status']=='PASS':total+=result['metrics']['rows']
 same(total,172,'IAEA 172 native coefficients were not all evaluated/preserved')
 return {'tables':len(tables),'exact_coefficient_rows':total,'full_joint_pdf_qualified':False}

def nasa(a,dataset):
 """Check NASA19/56 selections against actual native fixed-column bytes.

 Original one-based line locators address the pinned thermo.inp. The source
 A15 name, metadata columns53-65/66-80, interval temperature fields, and
 original16-character coefficient slots are compared independently of the
 catalogue builder. Fortran-D exponents and signed zero markers remain exact
 strings. Native enthalpy stays J/mol, temperature K and source mass u; no
 thermodynamic formula, isotope rest mass, zero-of-energy or reactionQ is
 inferred. The empty/ignored coefficient slot is not turned into a value.
 """
 base,packet,sources=a.packet(dataset);native=sources['nasa-cea:nasa-cea-data-thermo.inp'];lines=a.native(base,native).decode('ascii').splitlines();out={}
 for ident,expected in(('nasa-selected-species',19),('nasa-selected-coefficients',56)):
  row=sources[ident];columns,rows=csv_data(a.native(base,row));actual=a.table_equals(ident,columns,rows);same(len(rows),expected,'NASA native selection count changed');out[ident]=len(rows)
  for values in actual.dictionaries():
   if ident=='nasa-selected-species':
    start=int(values['native_start_line'])-1;meta=lines[start+1].ljust(80)
    same(lines[start][:15].strip(),values['name'],'Native NASA species line mismatch')
    same(meta[52:65].strip(),values['molecular_weight_native_u'],'NASA source mass lexeme changed')
    same(meta[65:80].strip(),values['hf298_native_J_mol'],'NASA source enthalpy lexeme changed')
   else:
    start=int(values['native_interval_start_line'])-1;domain,left,right=[v.ljust(80)for v in lines[start:start+3]]
    same((domain[:11].strip(),domain[11:22].strip()),(values['temperature_min_K'],values['temperature_max_K']),'NASA temperature domain changed')
    tokens=[left[k:k+16].strip()for k in range(0,80,16)]+[right[k:k+16].strip()for k in(0,16,48,64)]
    same(tokens,[values[key]for key in('a1_native','a2_native','a3_native','a4_native','a5_native','a6_native','a7_native','b1_native','b2_native')],'NASA original Fortran coefficient spellings changed')
 return dict(out,oracle='Original native fixed-column source lines; no thermo model evaluated')

def exomol(a,dataset):
 """Compare ExoMol's original whitespace/bz2 tokens, including signed data.

 Bzip2 decompression is lossless. The actual OH definition controls column
 names/count. nm, cm2/molecule, cm-1, s and K remain source conventions;
 spectroscopic energies are not formation energies. Negative H2O samples
 and their exact indices survive; they are not clipped into a probability
 distribution. No molecular reference frame or channel product is chosen.
 """
 base,packet,sources=a.packet(dataset);metrics={}
 for row in packet['sources']:
  if not row.get('table'):continue
  ident=row['table']['id'];parent=sources[row['derived_from'][0]];raw=a.native(base,parent)
  if parent['path'].endswith('.bz2'):raw=bz2.decompress(raw)
  native=tuple(tuple(line.split())for line in raw.decode('utf-8').splitlines())
  if any(not values for values in native):raise ValueError('Native ExoMol blank row requires explicit interpretation')
  actual=a.table_equals(ident,row['table']['columns'],native);item={'native_rows':len(native),'original_tokens_preserved':True}
  if 'h2o-dtu-'in ident:
   same(len(native),7933,'ExoMol H2O native row count changed');negative=[i for i,values in enumerate(native)if Decimal(values[1])<0];zero=sum(Decimal(values[1])==0 for values in native)
   same([i for i,values in enumerate(actual.rows)if Decimal(values[1])<0],negative,'ExoMol negative positions changed')
   same(len(negative),row['quantity_metadata']['source_checks']['negative_sample_count'],'ExoMol negative count metadata changed')
   same(zero,row['quantity_metadata']['source_checks']['zero_sample_count'],'ExoMol zero count metadata changed');item.update(negative_samples=len(negative),zero_samples=zero)
  elif 'states' in ident:
   definition=json.loads(a.native(base,sources['exomol-20261119:oh-mythos-def-json']))['dataset']['states'];same(len(native),10606,'OH native states count changed');same(len(native),definition['number_of_states'],'OH source definition count changed');same(actual.columns,tuple(field['name']for field in definition['states_file_fields']),'OH source definition columns changed');item['spectroscopic_energy_not_formation_Q']=True
  elif ident.endswith(':oh-mythos-pf:table'):same(len(native),10000,'OH partition source rows changed')
  metrics[ident]=item
 same(len(metrics),4,'Required four ExoMol source tables not all checked');return metrics

def surface(a,dataset):
 """Use original ZIP member bytes as measured/model table oracle.

 Four named experimental files have17 rows;14 explicitly model-output files
 have2033. Source degrees, amu/ion and atoms/ion columns are not recalculated.
 Error columns are not assigned a confidence level. These are sputtered
 material yields, not a measured incident-ion retention or charge-state law.
 """
 base,packet,sources=a.packet(dataset);archive_row=sources['surface.tuwien.hjn8d-hyh84.archive'];archive=zipfile.ZipFile(io.BytesIO(a.native(base,archive_row)));measured=simulated=mr=sr=0
 for row in packet['sources']:
  if not row.get('table'):continue
  member=row['quantity_metadata']['native_archive_member'];raw=archive.read(member);same(a.native(base,row),raw,'Surface source projection differs from ZIP member');a.stored_equals(row['id'],raw);columns,rows=csv_data(raw);a.table_equals(row['table']['id'],columns,rows)
  role=row['quantity_metadata']['role'];expected='measured-sputter-yield'if'Experiment'in member else'model-output-sputter-yield'
  same(role,expected,'Measured versus simulated surface role changed')
  if 'Experiment'in member:measured+=1;mr+=len(rows)
  else:simulated+=1;sr+=len(rows)
 same((measured,simulated,mr),(4,14,17),'Native surface experimental/simulation membership changed')
 return {'measured_tables':measured,'measured_rows':mr,'simulation_tables':simulated,'simulation_rows':sr,'incident_retention_qualified':False}

def archived_packets(a):
 # Exact explicit declarations are part of the immutable release DAG, rather
 # than historical filesystem locators alone. Old staging builds lack this
 # capability and are SKIPPED; a present but wrong pin is corruption/ERROR.
 prefix=Path(a.db.identity['release_manifest_path']).parent
 release=json.loads((a.db.root/a.db.identity['release_manifest_path']).read_text())
 provenance=a.db.root/prefix/'provenance/import.json'
 details=json.loads(provenance.read_text());entries=details.get('explicit_source_packets')
 if not isinstance(entries,list)or not entries:raise MissingOriginal('Release does not archive explicit source packet declarations')
 same(len(entries),6,'Final acquired release must archive all six explicit source packets')
 expected={(d['source_packet']['path'],d['source_packet']['sha256'])for d in a.db.datasets()if d.get('source_packet')}
 same({(r['path'],r['sha256'])for r in entries},expected,'Archived packet set differs from admitted catalogue declarations')
 for row in entries:
  relative=row.get('archived_path')
  if not isinstance(relative,str)or Path(relative).is_absolute()or'..'in Path(relative).parts:raise ValueError('Invalid archived packet path')
  if not relative.startswith(str(prefix/'provenance/source_packets')+'/'):raise ValueError('Archived packet not in declared release provenance')
  if release['members'].get(relative)!=row['sha256']:raise ValueError('Archived packet lacks exact release member pin')
  archived=a.original(a.db.root/relative,row['sha256']);original=a.original(Path(row['path']),row['sha256'])
  same(archived,original,'Archived packet bytes differ from original declaration')
 return {'archived_original_packets':len(entries),'release_member_sha_pins_verified':len(entries)}

def denial(a):
 """Verify every held source/table refuses actual numerical payload access.

 Public metadata may be inspected. Managed source paths must remain empty
 and both real APIs must return PermissionRequired. Unavailable data is not
 silently recast as a permissible zero-row table or a successful source.
 """
 count=0;table_count=0
 for source in a.db.sources():
  if source['storage_permission']!='allowed' or source['evaluation_permission']!='allowed':
   same(source['path'],'','Held numerical source has an imported payload path')
   try:a.db.source_path(source['source_id'])
   except DatabaseError as e:same(e.code,'PermissionRequired','Held source returned wrong denial');count+=1
   else:raise AssertionError('Held source unexpectedly accessible')
 for table in a.db.tables():
  if table['evaluation_permission']!='allowed':
   same(table['path'],'','Held table has an imported payload path')
   try:a.db.table(table['table_id'])
   except DatabaseError as e:same(e.code,'PermissionRequired','Held table returned wrong denial');table_count+=1
   else:raise AssertionError('Held table unexpectedly accessible')
 return {'permission_denied_sources':count,'permission_denied_tables':table_count}

def unchanged_inputs(a):
 """Close the read-only audit by rehashing every original/lock/executable."""
 for name,digest in a.inputs.items():
  if sha(Path(name))!=digest:raise ValueError('Original input changed during QA: '+name)
 return {'unchanged_original_inputs':len(a.inputs)}

def main(argv=None):
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--database-lock',type=Path,required=True);parser.add_argument('--output-dir',type=Path,required=True);parser.add_argument('--native-probe',type=Path,help='Existing compiled test/chemical_database/probe.cpp; no compilation performed');args=parser.parse_args(argv);args.output_dir.mkdir(parents=True,exist_ok=False);a=None;started=datetime.now(timezone.utc);result={'status':'ERROR'}
 try:
  lock=args.database_lock.resolve();before=sha(lock);db=Database.from_lock(lock);a=Audit(db,args.output_dir,args.native_probe);a.inputs[str(lock)]=before;a.check('actual-deep-release-integrity',lambda:db.verify(deep=True));datasets={d['id']:d for d in db.datasets()}
  required=('collision-input-acquisition-20261119','nasa-cea-thermochemistry-2f79a647','exomol-physical-source-inputs-20261119','surface.tuwien.regolith.hjn8d-hyh84')
  for ident in required:
   if ident not in datasets:raise MissingOriginal('Required acquired dataset absent from release: '+ident)
  for dataset in db.datasets():
   if dataset.get('source_packet'):a.check('original-packet-bytes:'+dataset['id'],lambda dataset=dataset:packet_sources(a,dataset))
  a.check('IAEA21-native-coefficient-tables',lambda:collisions(a,datasets[required[0]]));a.check('NASA19-and56-native-fixed-column-strings',lambda:nasa(a,datasets[required[1]]));a.check('ExoMol-original-token-and-negative-preservation',lambda:exomol(a,datasets[required[2]]));a.check('surface-original-ZIP-and-measured-model-separation',lambda:surface(a,datasets[required[3]]));a.check('all-six-original-packets-archived-and-pinned',lambda:archived_packets(a));a.check('all-held-payload-access-denied',lambda:denial(a));
  if a.native_probe:a.check('all45-new-CSV-native-Python-original-comparisons',lambda:(same(len(a.native_checked),45,'Not all45native table queries were evaluated')or {'native_tables':len(a.native_checked),'native_executable_sha256':a.inputs[str(a.native_probe)]}))
  a.check('lock-and-original-evidence-unchanged',lambda:unchanged_inputs(a))
  counts={s:sum(r['status']==s for r in a.rows)for s in('PASS','FAIL','ERROR','SKIPPED')};status=next((s for s in('ERROR','FAIL','SKIPPED')if counts[s]),'PASS');result={'status':status,'checks':a.rows,'counts':dict(counts,TOTAL=len(a.rows)),'identity':plain(db.identity),'inputs':a.inputs,'scope':'Real acquired source/lexical/permission QA only; physical model qualification and observational validation NOT EVALUATED','physical_activation':False,'native_comparison':{'requested':bool(a.native_probe),'query_count':len(a.native_checked),'executable':str(a.native_probe)if a.native_probe else None},'optional_extra_source_packets':'Accepted and original/hash/quantity checked when present; no ExoAtom physical model is selected'}
 except MissingOriginal as e:result={'status':'SKIPPED','reason':str(e),'checks':a.rows if a else[]}
 except Exception as e:result={'status':'ERROR','reason':str(e),'checks':a.rows if a else[]};(args.output_dir/'exception.log').write_text(traceback.format_exc())
 if 'counts'not in result:
  result['checks'].append({'name':'required-real-release-prerequisite','status':result['status'],'reason':result.get('reason','')});result['counts']={status:sum(row['status']==status for row in result['checks'])for status in('PASS','FAIL','ERROR','SKIPPED')};result['counts']['TOTAL']=len(result['checks'])
 result.update(schema='amps-real-acquired-release-qa-v1',started_utc=started.isoformat(),finished_utc=datetime.now(timezone.utc).isoformat(),source_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),script_sha256=sha(Path(__file__)),command=[sys.executable,'-B',str(Path(__file__).resolve()),'--database-lock',str(args.database_lock),'--output-dir',str(args.output_dir)]+(['--native-probe',str(args.native_probe)]if args.native_probe else[]))
 (args.output_dir/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');(args.output_dir/'audit.log').write_text('\n'.join(a.log if a else[])+'\n');print(json.dumps({'status':result['status'],'counts':result.get('counts'),'result':str(args.output_dir/'result.json')},indent=2));return {'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[result['status']]
if __name__=='__main__':raise SystemExit(main())
