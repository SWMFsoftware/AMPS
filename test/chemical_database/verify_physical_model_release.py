#!/usr/bin/env python3
"""Explicit, read-only QA for a real physical-source/model database release.

Every available CSV cell is compared with its pinned native/original table bytes,
through production Python and the actual C++ Catalogue probe. Selected explicit
packets additionally bind original members, parent links, recipes, provenance,
rights and native quantity metadata. Model definitions are retrieved as inert
metadata, never evaluated. This is source/software QA, not physical validation.
"""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time
import traceback
from verify_acquired_release import (Audit,Database,DatabaseError,MissingOriginal,
                                     ROOT,csv_data,plain,same,sha,unchanged_inputs)


def strict_json(raw):
    """Refuse duplicate/nonfinite metadata; do not repair an authority chain."""
    def unique(rows):
        result={}
        for key,value in rows:
            if key in result:raise ValueError('duplicate JSON key '+key)
            result[key]=value
        return result
    def bad(value):raise ValueError('nonfinite JSON '+value)
    return json.loads(raw,object_pairs_hook=unique,parse_constant=bad)


class ReleaseAudit(Audit):
    """Use the maintained query adapter, retaining per-command stdout/stderr."""
    def __init__(self,database,output,native_probe,deadline):
        super().__init__(database,output,native_probe)
        self.deadline=deadline
        self.packet_table_ids=set()
        self.denied_native_ids=set()
    def remaining(self):
        left=self.deadline-time.monotonic()
        if left<=0:raise RuntimeError('original release-QA wall budget exhausted; results are incomplete')
        return left
    def query(self,operation,identifier,capability):
        """Run one real Catalogue query; signal/timeout/non-JSON is ERROR.

        A valid PermissionRequired result is assessed by the caller. The probe
        itself returns0 after reporting typed query refusal; process failure is
        never counted as successful refusal, even if output includes that text.
        """
        identity=self.db.identity
        command=[str(self.native_probe),str(self.db.root),identity['export_path'],identity['export_sha256'],
                 identity['release_id'],identity['release_manifest_sha256'],identity['profile_id'],identity['profile_sha256'],
                 capability,operation,identifier]
        stem=self.output/'native-queries'/hashlib.sha256((operation+'\0'+identifier).encode()).hexdigest()[:20]
        process=subprocess.run(command,capture_output=True,timeout=min(30,self.remaining()))
        stem.with_suffix('.stdout.json').write_bytes(process.stdout);stem.with_suffix('.stderr.log').write_bytes(process.stderr)
        receipt={'command':command,'returncode':process.returncode,'executable_sha256':self.inputs[str(self.native_probe)],
                 'stdout_sha256':hashlib.sha256(process.stdout).hexdigest(),'stderr_sha256':hashlib.sha256(process.stderr).hexdigest()}
        stem.with_suffix('.command.json').write_text(json.dumps(receipt,indent=2)+'\n')
        if process.returncode!=0:raise RuntimeError('C++ probe process failed '+operation+':'+identifier)
        result=strict_json(process.stdout)
        if result.get('open_ok')is not True:raise RuntimeError('C++ Catalogue failed to open pinned release')
        return result

    def table_equals(self,identifier,columns,rows):
        """Compare lexical views with the same original absolute wall deadline.

        The shared Audit helper uses a fixed per-process timeout; this override
        uses query() so even a late table cannot start a fresh30-second budget.
        No interpretation or numerical conversion is added by either reader.
        """
        self.remaining()
        actual=self.db.table(identifier)
        same(actual.columns,tuple(columns),'CSV columns changed '+identifier)
        same(actual.rows,tuple(tuple(row)for row in rows),'CSV source tokens/order changed '+identifier)
        if identifier not in self.native_checked:
            native=self.query('table',identifier,'table-access')
            same(native.get('query_ok'),True,'native CSV query refused '+identifier)
            same(native.get('rows'),[list(row)for row in actual.rows],'native/Python/source CSV cells differ '+identifier)
            self.native_checked.add(identifier)
        return actual


def selected_packet(a,path,datasets,release,archived):
    """Cross-bind selected original packet, admitted metadata and archived bytes.

    All declared source bytes are authenticated; a held source is only examined
    in its original research packet and must remain inaccessible in the managed
    database. Derived parent/recipe identities are checked without executing a
    source's scientific model or selecting a more convenient packet version.
    """
    path=Path(path).resolve()
    if not path.is_file():raise MissingOriginal('selected original source packet unavailable '+str(path))
    pins=[d['source_packet']for d in datasets.values()if d.get('source_packet')and
          Path(d['source_packet']['path']).resolve()==path]
    if len(pins)!=1:raise ValueError('selected original packet has no unique admitted digest')
    digest=pins[0]['sha256'];value=strict_json(a.original(path,digest))
    if value.get('schema')!='amps-chemical-source-packet-v1':raise ValueError('unsupported source packet schema')
    if value.get('physical_activation',False)is not False:raise ValueError('packet attempts physical activation')
    identifier=value['id']
    if identifier not in datasets:raise MissingOriginal('selected packet dataset is absent from release '+identifier)
    dataset=datasets[identifier];pin=dataset.get('source_packet')
    if not isinstance(pin,dict)and not hasattr(pin,'get'):raise ValueError('admitted packet identity missing')
    same((pin.get('path'),pin.get('sha256')),(str(path),digest),'selected packet identity differs from admitted catalogue '+identifier)
    matches=[row for row in archived if row['path']==str(path)and row['sha256']==digest]
    same(len(matches),1,'selected original packet missing/duplicated in release provenance '+identifier)
    retained=matches[0]['archived_path']
    relative=Path(retained)
    if relative.is_absolute()or '..'in relative.parts or 'provenance/source_packets/'not in retained:
        raise ValueError('archived packet path escapes declared provenance')
    if release['members'].get(retained)!=digest:raise ValueError('archived packet missing exact release member hash')
    same(a.original(a.db.root/retained,digest),path.read_bytes(),'archived packet differs from original')
    expected={row['id']:row for row in value['sources']}
    if len(expected)!=len(value['sources']):raise ValueError('duplicate original source ID')
    actual={row['id']:row for row in dataset['files']}
    same(set(actual),set(expected),'admitted source membership differs from selected packet')
    tables=0;parents=0;allowed=0
    for key,row in expected.items():
        a.remaining();native=a.original(path.parent/row['path'],row['sha256']);stored=actual[key]
        same(stored.get('observed_sha256'),row['sha256'],'admitted source SHA changed '+key)
        same(plain(stored.get('quantity_metadata')),row.get('quantity_metadata'),'native quantity/role/domain metadata changed '+key)
        same(stored.get('policy'),row['policy'],'source policy changed '+key)
        for item in row.get('rights_evidence',[]):a.original(path.parent/item['path'],item['sha256'])
        if stored['storage_permission']=='allowed'and stored['evaluation_permission']=='allowed':
            same(a.db.source_path(key).read_bytes(),native,'stored source bytes changed '+key);allowed+=1
        if row['role']=='derived-table':
            recipe=row['recipe'];a.original(path.parent/recipe['path'],recipe['sha256'])
            for dependency in row.get('quantity_metadata',{}).get('recipe_dependencies',[]):
                a.original(path.parent/dependency['path'],dependency['sha256'])
            same(plain(stored.get('derived_from')),row['derived_from'],'derived parent IDs changed '+key)
            for parent in row['derived_from']:
                if parent not in expected or expected[parent]['role']=='derived-table':raise ValueError('derived source has absent/non-native parent')
                a.original(path.parent/expected[parent]['path'],expected[parent]['sha256']);parents+=1
        if row.get('table'):
            table=row['table'];a.packet_table_ids.add(table['id']);tables+=1
            if stored['evaluation_permission']=='allowed'and table.get('format','csv')=='csv':
                columns,rows=csv_data(native);a.table_equals(table['id'],columns,rows)
    return {'packet_id':identifier,'native_or_derived_members':len(expected),'allowed_source_bytes_checked':allowed,
            'declared_tables':tables,'direct_parent_bindings':parents,'archived_packet_exact':True}


def archived_declarations(a,release,provenance):
    """Verify exact portable original-packet AND declared-recipe DAG members.

    Archived recipes are compared byte-for-byte without executing their models.
    Common helper dependencies remain independently pinned release members.
    Missing archive capability is unavailable; present malformed/mismatched pins
    are corruption. Current source files never replace an old expected digest.
    """
    packets=provenance.get('explicit_source_packets');recipes=provenance.get('packet_recipes')
    if not isinstance(packets,list)or not isinstance(recipes,list):raise MissingOriginal('portable packet/recipe archives unavailable')
    expected_packets={(d['source_packet']['path'],d['source_packet']['sha256'])for d in a.db.datasets()if d.get('source_packet')}
    same({(p['path'],p['sha256'])for p in packets},expected_packets,'archived packet set differs from admitted selection')
    expected_recipes={};count=0
    for row in packets:
        relative=Path(row['archived_path'])
        if relative.is_absolute()or '..'in relative.parts or '/provenance/source_packets/'not in str(relative):raise ValueError('invalid archived packet path')
        same(release['members'].get(str(relative)),row['sha256'],'archived packet not release-member pinned')
        original=a.original(Path(row['path']),row['sha256'])
        same(a.original(a.db.root/relative,row['sha256']),original,'archived packet bytes changed')
        document=strict_json(original)
        for source in document['sources']:
            declarations=([source['recipe']]if source.get('recipe')else[])
            declarations+=source.get('quantity_metadata',{}).get('recipe_dependencies',[])
            for recipe in declarations:
                expected_recipes[(row['sha256'],recipe['path'])]=(Path(row['path']).parent/recipe['path'],recipe['sha256'])
        count+=1
    observed={}
    for row in recipes:
        key=(row['packet_sha256'],row['original_path'])
        if key in observed:raise ValueError('duplicate archived packet recipe')
        observed[key]=row
    same(set(observed),set(expected_recipes),'declared recipe archive membership differs')
    for key,(original,digest)in expected_recipes.items():
        row=observed[key];relative=Path(row['archived_path'])
        if relative.is_absolute()or '..'in relative.parts or '/provenance/packet_recipes/'not in str(relative):raise ValueError('invalid archived recipe path')
        same(row['sha256'],digest,'archived recipe pin differs from packet')
        same(release['members'].get(str(relative)),digest,'archived recipe lacks exact release-member pin')
        same(a.original(a.db.root/relative,digest),a.original(original,digest),'archived recipe bytes changed')
    return {'exact_original_packets_archived':count,'exact_declared_recipes_archived':len(expected_recipes)}


def available_tables(a):
    """Read every available CSV without float conversion or units changes.

    Original source records are used whenever present; otherwise the immutable
    release CSV is the lexical oracle. Both cases are identified in metrics.
    Exact cell identity is appropriate for copying/access software, never an
    invented physical-agreement tolerance. Opaque tables are not recast as CSV.
    """
    originals={}
    for dataset in a.db.datasets():
        for row in dataset.get('files',()):
            if row.get('table'):originals[row['table']['id']]=row
    count=0;cells=0;original_count=0
    for table in a.db.tables():
        a.remaining()
        if table['format']!='csv'or table['availability']!='available'or table['evaluation_permission']!='allowed':continue
        identifier=table['table_id'];path=a.db.root/table['path']
        raw=a.original(path,table['sha256'])
        row=originals.get(identifier)
        if row and row.get('path')and Path(row['path']).is_file():
            original=a.original(Path(row['path']),row.get('observed_sha256',row.get('sha256')))
            same(raw,original,'release CSV differs from available original '+identifier);original_count+=1
        columns,rows=csv_data(raw);a.table_equals(identifier,columns,rows);count+=1;cells+=sum(len(row)for row in rows)
    if count==0:raise MissingOriginal('release has no allowed available CSV tables')
    same(len(a.native_checked),count,'not every actual available CSV was compared through C++')
    return {'available_CSV_tables':count,'exact_data_cells':cells,'available_original_table_byte_oracles':original_count,
            'native_Python_original_equality':True,'scientific_formulas_evaluated':0}


def models(a):
    """Cross-check inert model API records against the published native CSV.

    This verifies accessibility/immutability and serialization only. Equations,
    source citations, limitations and input requirements are documentary strings;
    no rate, cascade, state assignment or approximation is activated here.
    """
    if not hasattr(a.db,'models')or not hasattr(a.db,'model'):raise MissingOriginal('model descriptor API unavailable')
    values=list(a.db.models())
    if not values:raise MissingOriginal('explicit model definitions not selected in this release')
    table=a.db.table('amps-reduced-model-definitions-v1');rows={row['model_id']:row for row in table.dictionaries()}
    same(set(rows),{row['id']for row in values},'model API and published CSV membership differ')
    for value in values:
        identifier=value['id'];actual=plain(a.db.model(identifier));same(actual,plain(value),'model query changed list record')
        same(plain(a.db.describe(identifier)),actual,'model describe differs from query')
        same(actual.get('physical_activation'),False,'model descriptor activates physics')
        row=rows[identifier]
        for key in ('implementation','evidence_kind'):same(actual[key],row[key],'model scalar differs from native table')
        for key in ('equations','inputs','outputs','assumptions','limitations','citations'):
            same(actual[key],strict_json(row[key+'_json']),'model documentary JSON differs from native table')
        try:a.db.model(identifier)['physical_activation']=True
        except TypeError:pass
        else:raise AssertionError('model descriptor is mutable')
    return {'model_descriptors':len(values),'inert_immutable_views':True,'model_table_compared_in_native_sweep':True}


def denied(a):
    """Test every held numerical source/table through both real permission gates."""
    requests=[];source_count=table_count=0
    for row in a.db.sources():
        if row['storage_permission']=='allowed'and row['evaluation_permission']=='allowed':continue
        same(row['path'],'','held source carries managed payload path')
        try:a.db.source_path(row['source_id'])
        except DatabaseError as error:same(error.code,'PermissionRequired','held source wrong denial code')
        else:raise AssertionError('held source was numerically accessible')
        requests.append(('source-path',row['source_id'],'source-access'));source_count+=1
    for row in a.db.tables():
        if row['evaluation_permission']=='allowed':continue
        same(row['path'],'','held table carries managed payload path')
        try:a.db.table(row['table_id'])
        except DatabaseError as error:same(error.code,'PermissionRequired','held table wrong denial code')
        else:raise AssertionError('held table was numerically accessible')
        requests.append(('table',row['table_id'],'table-access'));table_count+=1
    if not requests:raise MissingOriginal('no historical held inputs available for preservation check')
    # One thin native operation calls the SAME production guard for every
    # held metadata ID. This avoids hundreds of redundant immutable opens;
    # no ID sampling or alternative permission implementation is involved.
    result=a.query('denial-sweep','', 'source-access')
    same(result.get('denial_sweep_ok'),True,'native held-query sweep failed')
    native_rows=result.get('denials')
    if not isinstance(native_rows,list):raise ValueError('native denial sweep result missing')
    expected={(operation,identifier)for operation,identifier,_ in requests}
    observed=[(row.get('kind'),row.get('id'))for row in native_rows]
    same(len(observed),len(set(observed)),'native held-query sweep duplicates IDs')
    same(set(observed),expected,'native held-query sweep omitted/added metadata IDs')
    for row in native_rows:
        same(row.get('query_ok'),False,'C++ held query unexpectedly succeeded')
        same(row.get('query_code'),'PermissionRequired','C++ held query wrong denial code')
        same(row.get('output_preserved'),True,'C++ held query changed refused output')
        a.denied_native_ids.add(row['kind']+':'+row['id'])
    return {'denied_sources_both_APIs':source_count,'denied_tables_both_APIs':table_count}


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--database-lock',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--source-packet',type=Path,action='append',default=[])
    parser.add_argument('--native-probe',type=Path)
    parser.add_argument('--wall-seconds',type=float,default=600)
    args=parser.parse_args(argv)
    if not math.isfinite(args.wall_seconds)or not 0<args.wall_seconds<=3600:parser.error('wall seconds must be positive and at most3600')
    args.output_dir.mkdir(parents=True,exist_ok=False)
    started=datetime.now(timezone.utc);deadline=time.monotonic()+args.wall_seconds;a=None;result={}
    try:
        lock=args.database_lock.resolve()
        if not lock.is_file():raise MissingOriginal('explicit database lock unavailable')
        db=Database.from_lock(lock)
        if args.native_probe is None:raise MissingOriginal('actual C++ Catalogue --native-probe is required for this QA scope')
        a=ReleaseAudit(db,args.output_dir,args.native_probe,deadline);a.inputs[str(lock)]=sha(lock)
        a.check('actual-deep-release-integrity',lambda:db.verify(deep=True))
        release=strict_json((db.root/db.identity['release_manifest_path']).read_bytes())
        prefix=Path(db.identity['release_manifest_path']).parent
        provenance=strict_json((db.root/prefix/'provenance/import.json').read_bytes())
        archived=provenance.get('explicit_source_packets')
        if not isinstance(archived,list):raise MissingOriginal('explicit source packet provenance unavailable')
        datasets={row['id']:row for row in db.datasets()};seen=set()
        a.check('all-admitted-original-packets-and-declared-recipes-archived',lambda:archived_declarations(a,release,provenance))
        if not args.source_packet:raise MissingOriginal('explicit source packets are required as original new-source oracles')
        provided={path.resolve()for path in args.source_packet}
        admitted={Path(row['path']).resolve()for row in archived}
        if provided-admitted:raise ValueError('requested packet is not an admitted release selection')
        if admitted-provided:raise MissingOriginal('all admitted original --source-packet arguments are required; missing '+str(len(admitted-provided)))
        for path in args.source_packet:
            resolved=path.resolve()
            if resolved in seen:raise ValueError('duplicate requested source packet')
            seen.add(resolved)
            a.check('selected-packet:'+str(path),lambda path=path:selected_packet(a,path,datasets,release,archived))
        a.check('all-available-CSV-Python-native-original-cells',lambda:available_tables(a))
        a.check('inert-model-descriptor-retrieval',lambda:models(a))
        a.check('all-historical-held-permissions-both-APIs',lambda:denied(a))
        a.check('lock-source-and-native-executable-unchanged',lambda:unchanged_inputs(a))
        counts={key:sum(row['status']==key for row in a.rows)for key in('PASS','FAIL','ERROR','SKIPPED')}
        status=next((key for key in('ERROR','FAIL','SKIPPED')if counts[key]),'PASS')
        result={'status':status,'checks':a.rows,'counts':dict(counts,TOTAL=len(a.rows)),
                'identity':plain(db.identity),'inputs':a.inputs,'native_table_queries':len(a.native_checked),
                'native_permission_queries':len(a.denied_native_ids),'selected_packet_tables':sorted(a.packet_table_ids)}
    except MissingOriginal as error:result={'status':'SKIPPED','reason':str(error),'checks':a.rows if a else[]}
    except Exception as error:
        result={'status':'ERROR','reason':str(error),'checks':a.rows if a else[]}
        (args.output_dir/'exception.log').write_text(traceback.format_exc())
    if 'counts'not in result:
        result['checks'].append({'name':'required-release-prerequisite','status':result['status'],'reason':result.get('reason','')})
        result['counts']={key:sum(row['status']==key for row in result['checks'])for key in('PASS','FAIL','ERROR','SKIPPED')}
        result['counts']['TOTAL']=len(result['checks'])
    if time.monotonic()>deadline:
        result['checks'].append({'name':'original-wall-budget','status':'ERROR','reason':'release QA exceeded original wall budget'})
        result['counts']['ERROR']+=1;result['counts']['TOTAL']+=1;result['status']='ERROR'
    result.update(schema='amps-physical-model-release-qa-v1',started_utc=started.isoformat(),finished_utc=datetime.now(timezone.utc).isoformat(),
                  scope='Real source/lexical/API/permission QA only; physical qualification/observational scoring NOT EVALUATED',physical_activation=False,
                  script_sha256=sha(Path(__file__)),audit_helper_sha256=sha(Path(__file__).with_name('verify_acquired_release.py')),
                  source_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                  command=[sys.executable,'-B',str(Path(__file__).resolve())]+sys.argv[1:])
    (args.output_dir/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    (args.output_dir/'audit.log').write_text('\n'.join(a.log if a else[])+'\n')
    print(json.dumps({'status':result['status'],'counts':result['counts'],'result':str(args.output_dir/'result.json')},indent=2))
    return {'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[result['status']]

if __name__=='__main__':raise SystemExit(main())
