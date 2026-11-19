"""Production native/Python data-access controls, never physical observations.

Fixtures contain administrative strings only. Their digest chains are rebuilt
only to exercise a specific parser/permission guard; real release bytes are
never rewritten. Numerical service tests remain the independent maintained
Moon/chemistry controls rather than duplicated formulas here.
"""
import csv
import hashlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import unittest
from types import SimpleNamespace

SOURCE_HEADER='source_id,path,sha256,media_type,availability,storage_permission,evaluation_permission,capabilities,source_uri,notes'.split(',')
TABLE_HEADER='table_id,path,sha256,format,schema_id,source_ids,availability,evaluation_permission,capabilities,columns,row_count,notes'.split(',')
CONTEXT=None


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def csv_bytes(header,rows):
    stream=io.StringIO(newline='')
    writer=csv.writer(stream,lineterminator='\n');writer.writerow(header);writer.writerows(rows)
    return stream.getvalue()


class Fixture:
    """A complete explicitly labelled administrative export, no science data."""
    def __init__(self,root):
        self.root=Path(root);self.root.mkdir(parents=True,exist_ok=False)
        (self.root/'source-z.txt').write_text('ADMINISTRATIVE PARSER FIXTURE Z; NOT OBSERVATIONS\n')
        (self.root/'source-a.txt').write_text('ADMINISTRATIVE PARSER FIXTURE A; NOT OBSERVATIONS\n')
        (self.root/'table.csv').write_text(csv_bytes(['marker','label'],[['z','quoted, marker'],['a','line1\nline2']]))
        (self.root/'profile.json').write_text(json.dumps({'schema':'amps-chemical-profile-v1',
            'id':'fixture-profile','capabilities':['parser-control'],
            'scope':'administrative access controls only; NOT observations'})+'\n')
        (self.root/'catalogue.json').write_text(json.dumps({'schema':'amps-chemical-catalogue-v1',
            'datasets':[{'id':'fixture-z','scope':'administrative parser fixture'},
                        {'id':'fixture-a','scope':'administrative parser fixture'}],
            'species':[],'reactions':[],'table_metadata':{'table-z':{'scope':'NOT science data'}}})+'\n')
        self.sources=[['source-z','source-z.txt',digest(self.root/'source-z.txt'),'text/plain','available','allowed','allowed','parser-control','urn:fixture:z','not science'],
                      ['source-a','source-a.txt',digest(self.root/'source-a.txt'),'text/plain','available','allowed','allowed','parser-control','urn:fixture:a','not science']]
        self.tables=[['table-z','table.csv',digest(self.root/'table.csv'),'csv','administrative-fixture-v1','source-z;source-a','available','allowed','parser-control','marker;label','2','not science']]
        self.fields={'schema':'amps-chemical-runtime-v1','codec':'csv-block-v1','release_id':'fixture-release',
            'release_manifest_path':'release.json','release_manifest_sha256':'',
            'profile_id':'fixture-profile','profile_path':'profile.json','profile_sha256':digest(self.root/'profile.json'),
            'source_index_path':'sources.csv','table_index_path':'tables.csv',
            'profile_capabilities':'parser-control','topology_fingerprint':'none'}
        self.publish()

    def publish(self,extra=''):
        (self.root/'sources.csv').write_text(csv_bytes(SOURCE_HEADER,self.sources))
        (self.root/'tables.csv').write_text(csv_bytes(TABLE_HEADER,self.tables))
        self.fields['source_index_sha256']=digest(self.root/'sources.csv')
        self.fields['table_index_sha256']=digest(self.root/'tables.csv')
        # These are administrative authority-chain fixtures, not repaired
        # source releases. Rebuilding their graph permits a specific semantic
        # guard to be reached without first failing an unrelated digest.
        release={'schema':'amps-chemical-release-v1','release_id':'fixture-release',
            'catalogue':{'path':'catalogue.json','sha256':digest(self.root/'catalogue.json')},
            'members':{name:digest(self.root/name) for name in
                       ('profile.json','catalogue.json','sources.csv','tables.csv')}}
        (self.root/'release.json').write_text(json.dumps(release,sort_keys=True)+'\n')
        self.fields['release_manifest_sha256']=digest(self.root/'release.json')
        (self.root/'export.in').write_text('#begin chemical_database: runtime\n'+
            ''.join(key+' = '+value+'\n' for key,value in self.fields.items())+extra+'#end chemical_database\n')


class NativeProtocol(unittest.TestCase):
    def setUp(self):
        self.fixture=Fixture(CONTEXT.output/self.id().split('.')[-1])

    def call(self,operation='catalogue',identifier='',capability='parser-control',fixture=None,**overrides):
        return CONTEXT.native(fixture or self.fixture,operation,identifier,capability,**overrides)

    def rejected(self,result):
        self.assertFalse(result['open_ok'],result)
        self.assertFalse(result['loaded'],result)

    def test_exact_source_order_and_quoted_rows(self):
        result=self.call('table','table-z')
        self.assertTrue(result['open_ok'],result);self.assertTrue(result['query_ok'],result)
        self.assertEqual(result['source_ids'],['source-z','source-a'])
        self.assertEqual(result['table_ids'],['table-z'])
        self.assertEqual(result['rows'],[['z','quoted, marker'],['a','line1\nline2']])

    def test_rows_outlive_catalogue(self):
        result=self.call('lifetime','table-z')
        self.assertTrue(result['query_ok'],result)
        self.assertEqual(result['rows'][0],['z','quoted, marker'])

    def test_null_outputs_and_missing_metadata_preserved(self):
        result=self.call('guards','table-z')
        self.assertTrue(result['open_ok'],result);self.assertTrue(result['guards_ok'],result)

    def test_concurrent_immutable_owned_reads(self):
        result=self.call('concurrent','table-z')
        self.assertTrue(result['open_ok'],result);self.assertTrue(result['concurrent_reads_ok'],result)

    def test_published_sha256_empty_and_abc_vectors(self):
        result=self.call('sha-vectors')
        self.assertTrue(result['open_ok'],result);self.assertTrue(result['known_sha_vectors_ok'],result)

    def test_rejected_reload_preserves_immutable_storage(self):
        result=self.call('rollback')
        self.assertTrue(result['open_ok'],result);self.assertFalse(result['reopen_ok'])
        self.assertTrue(result['old_identity_preserved']);self.assertTrue(result['copy_identity_shared'])

    def test_wrong_export_identity(self):
        self.rejected(self.call(expected_export='0'*64))

    def test_wrong_release_identity(self):
        self.rejected(self.call(expected_release='0'*64))

    def test_wrong_profile_identity(self):
        self.rejected(self.call(expected_profile='0'*64))

    def test_changed_authenticated_index(self):
        with (self.fixture.root/'sources.csv').open('a') as stream:stream.write('\n')
        self.rejected(self.call())

    def test_unknown_schema(self):
        self.fixture.fields['schema']='unsupported-administrative-schema';self.fixture.publish()
        self.rejected(self.call())

    def test_unknown_codec(self):
        self.fixture.fields['codec']='unknown';self.fixture.publish();self.rejected(self.call())

    def test_unknown_key(self):
        self.fixture.publish('unowned_key = never_ignored\n');self.rejected(self.call())

    def test_duplicate_key(self):
        self.fixture.publish('schema = amps-chemical-runtime-v1\n');self.rejected(self.call())

    def test_include_forbidden_even_if_target_exists(self):
        (self.fixture.root/'empty.in').write_text('! explicitly empty grammar fixture\n')
        self.fixture.publish('#include empty.in\n');self.rejected(self.call())

    def test_duplicate_source_identifier(self):
        self.fixture.sources.append(self.fixture.sources[0][:]);self.fixture.publish();self.rejected(self.call())

    def test_duplicate_table_identifier(self):
        self.fixture.tables.append(self.fixture.tables[0][:]);self.fixture.publish();self.rejected(self.call())

    def test_unknown_table_source(self):
        self.fixture.tables[0][5]='not-declared';self.fixture.publish();self.rejected(self.call())

    def test_parent_path_refused(self):
        self.fixture.sources[0][1]='../escape.txt';self.fixture.publish();self.rejected(self.call())

    def test_absolute_path_refused(self):
        self.fixture.sources[0][1]=str(self.fixture.root/'source-z.txt');self.fixture.publish();self.rejected(self.call())

    def test_backslash_path_refused(self):
        self.fixture.sources[0][1]='sub\\source-z.txt';self.fixture.publish();self.rejected(self.call())

    def test_symlink_escape_refused(self):
        outside=self.fixture.root.parent/(self.fixture.root.name+'-outside.txt')
        outside.write_text('administrative outside-root path fixture\n')
        (self.fixture.root/'escaped.txt').symlink_to(outside)
        self.fixture.sources[0][1]='escaped.txt';self.fixture.sources[0][2]=digest(outside)
        self.fixture.publish();self.rejected(self.call())

    def test_explicit_resource_limit(self):
        result=self.call(file_limit=1);self.rejected(result)
        self.assertEqual(result['open_code'],'ResourceLimit')

    def test_explicit_total_resource_limit(self):
        result=self.call(total_limit=1);self.rejected(result)
        self.assertEqual(result['open_code'],'ResourceLimit')

    def test_profile_capability_not_inferred(self):
        self.rejected(self.call(capability='unapproved-capability'))

    def test_unknown_table_preserves_output(self):
        result=self.call('table','not-declared')
        self.assertTrue(result['open_ok'],result);self.assertFalse(result['query_ok'],result)
        self.assertTrue(result['output_preserved'])

    def blocked_table(self,fixture=None):
        result=self.call('table','table-z',fixture=fixture)
        self.assertTrue(result['open_ok'],result);self.assertFalse(result['query_ok'],result)
        self.assertTrue(result['output_preserved'],result)

    def test_pending_source_evaluation_permission(self):
        self.fixture.sources[0][6]='pending';self.fixture.publish();self.blocked_table()

    def test_unknown_source_storage_permission(self):
        self.fixture.sources[0][5]='unknown';self.fixture.publish();self.blocked_table()

    def test_prohibited_table_evaluation_permission(self):
        self.fixture.tables[0][7]='prohibited';self.fixture.publish();self.blocked_table()

    def test_metadata_only_source_not_zero_table(self):
        self.fixture.sources[0][4]='metadata-only';self.fixture.publish();self.blocked_table()

    def test_missing_payload_can_remain_declared_unavailable(self):
        self.fixture.sources[0][1:3]=['',''];self.fixture.sources[0][4]='unavailable'
        self.fixture.publish();self.blocked_table()

    def test_source_capability_access_guard(self):
        self.fixture.sources[0][7]='another-declared-capability';self.fixture.publish()
        result=self.call('source-path','source-z')
        self.assertTrue(result['open_ok'],result);self.assertFalse(result['query_ok'],result)
        self.assertTrue(result['output_preserved'])

    def test_table_capability_access_guard(self):
        self.fixture.tables[0][8]='another-declared-capability';self.fixture.publish();self.blocked_table()

    def test_changed_table_bytes_refused_before_use(self):
        with (self.fixture.root/'table.csv').open('a') as stream:stream.write('extra,fixture\n')
        result=self.call('table','table-z')
        # Implementations may eagerly check all payloads while opening or
        # lazily authenticate only requested data. Both must refuse its use.
        self.assertFalse(result['open_ok'] and result.get('query_ok',False),result)
        if result['open_ok']:self.assertTrue(result['output_preserved'])

    def test_declared_row_count_mismatch_refused(self):
        self.fixture.tables[0][10]='3';self.fixture.publish()
        result=self.call('table','table-z')
        self.assertFalse(result['open_ok'] and result.get('query_ok',False),result)

    def test_relocation_uses_explicit_new_root(self):
        moved=self.fixture.root.parent/(self.fixture.root.name+'-moved')
        shutil.move(str(self.fixture.root),moved);self.fixture.root=moved
        result=self.call('table','table-z')
        self.assertTrue(result['query_ok'],result)
        self.assertEqual(result['rows'][1],['a','line1\nline2'])


class PythonProtocol(unittest.TestCase):
    """Actual Python service against the same administrative authority graph."""
    def setUp(self):
        self.api=CONTEXT.python_api()
        self.fixture=Fixture(CONTEXT.output/self.id().split('.')[-1])

    def open(self):
        f=self.fixture.fields
        return self.api.Database.open(self.fixture.root,f['release_id'],f['release_manifest_sha256'],
            f['profile_id'],f['profile_sha256'],'export.in',digest(self.fixture.root/'export.in'))

    def test_native_python_exact_order_and_rows(self):
        database=self.open();table=database.table('table-z','parser-control')
        native=CONTEXT.native(self.fixture,'table','table-z','parser-control')
        self.assertEqual([r['source_id'] for r in database.sources()],native['source_ids'])
        self.assertEqual([r['table_id'] for r in database.tables()],native['table_ids'])
        self.assertEqual(table.rows,tuple(tuple(row) for row in native['rows']))
        self.assertEqual(table.columns,('marker','label'));self.assertEqual(table.row_count,2)

    def test_frozen_recursive_metadata_and_table(self):
        database=self.open();before=repr(database.datasets());dataset=database.datasets()[0]
        with self.assertRaises(TypeError):dataset['id']='changed'
        table=database.table('table-z','parser-control')
        with self.assertRaises(TypeError):table.metadata['scope']='changed'
        with self.assertRaises((AttributeError,TypeError)):table.rows=(('edited',),)
        self.assertEqual(repr(database.datasets()),before)

    def test_owned_rows_survive_database_lifetime(self):
        database=self.open();table=database.table('table-z','parser-control');del database
        self.assertEqual(table.rows[0],('z','quoted, marker'))

    def test_failed_open_leaves_existing_object_unchanged(self):
        database=self.open();before=database.identity
        with self.assertRaises(self.api.DatabaseError):
            self.api.Database.open(self.fixture.root,'fixture-release','0'*64,'fixture-profile',
                self.fixture.fields['profile_sha256'],'export.in',digest(self.fixture.root/'export.in'))
        self.assertIs(database.identity,before)
        self.assertEqual(database.table('table-z','parser-control').rows[1],('a','line1\nline2'))

    def test_wrong_schema_and_duplicate_key_refused(self):
        self.fixture.fields['schema']='unknown';self.fixture.publish()
        with self.assertRaises(self.api.DatabaseError):self.open()
        self.fixture.fields['schema']='amps-chemical-runtime-v1';self.fixture.publish('codec = csv-block-v1\n')
        with self.assertRaises(self.api.DatabaseError):self.open()

    def test_hash_corruption_and_unknown_rights_fail_closed(self):
        (self.fixture.root/'sources.csv').write_text('changed administrative bytes\n')
        with self.assertRaises(self.api.DatabaseError):self.open()
        self.fixture.sources[0][5]='unknown';self.fixture.publish()
        with self.assertRaises(self.api.DatabaseError):self.open()

    def test_path_traversal_and_symlink_escape_refused(self):
        self.fixture.sources[0][1]='../escape';self.fixture.publish()
        with self.assertRaises(self.api.DatabaseError):self.open()
        outside=self.fixture.root.parent/(self.fixture.root.name+'-outside');outside.write_text('ADMINISTRATIVE ONLY\n')
        (self.fixture.root/'escape').symlink_to(outside)
        self.fixture.sources[0][1]='escape';self.fixture.sources[0][2]=digest(outside);self.fixture.publish()
        with self.assertRaises(self.api.DatabaseError):self.open()

    def test_missing_table_cell_retained_and_number_refused(self):
        (self.fixture.root/'table.csv').write_text(csv_bytes(['marker','label'],[['','missing administrative cell'],['1','present']]))
        self.fixture.tables[0][2]=digest(self.fixture.root/'table.csv');self.fixture.publish()
        table=self.open().table('table-z','parser-control')
        self.assertEqual(table.rows[0][0],'')
        with self.assertRaises(self.api.DatabaseError):table.number(0,'marker','undeclared-unit')

    def test_exact_padded_source_column_labels_preserved(self):
        # Real archives can pad column labels. This administrative analogue
        # checks exact lexical preservation, never normalization/renaming.
        columns=['marker           ',' label ']
        (self.fixture.root/'table.csv').write_text(csv_bytes(columns,[['z','quoted, marker'],['a','line1\nline2']]))
        self.fixture.tables[0][2]=digest(self.fixture.root/'table.csv')
        self.fixture.tables[0][9]=';'.join(columns);self.fixture.publish()
        native=CONTEXT.native(self.fixture,'table','table-z','parser-control')
        self.assertTrue(native['query_ok'],native)
        table=self.open().table('table-z','parser-control')
        self.assertEqual(table.columns,tuple(columns))
        self.assertEqual(table.rows,tuple(tuple(row) for row in native['rows']))

    def test_capability_missing_not_inferred(self):
        database=self.open()
        with self.assertRaises(self.api.DatabaseError):database.table('table-z','unapproved-capability')

    def test_permission_missing_not_inferred(self):
        self.fixture.sources[0][6]='pending';self.fixture.publish()
        with self.assertRaises(self.api.DatabaseError):self.open().table('table-z','parser-control')

    def test_relocation_without_original_path_dependency(self):
        moved=self.fixture.root.parent/(self.fixture.root.name+'-moved')
        shutil.move(str(self.fixture.root),moved);self.fixture.root=moved
        self.assertEqual(self.open().table('table-z','parser-control').rows[0],('z','quoted, marker'))

    def test_lock_explicit_identity_and_cli_order(self):
        lock=self.fixture.root/'database-lock.json';f=self.fixture.fields
        lock.write_text(json.dumps({'schema':'amps-chemical-database-lock-v1','database_root':'.',
            'release_id':f['release_id'],'release_manifest_sha256':f['release_manifest_sha256'],
            'profile_id':f['profile_id'],'profile_sha256':f['profile_sha256'],
            'export_path':'export.in','export_manifest_sha256':digest(self.fixture.root/'export.in')}))
        database=self.api.Database.from_lock(lock)
        self.assertEqual([r['id'] for r in database.datasets()],['fixture-z','fixture-a'])
        command=[CONTEXT.python_executable,'-B',str(CONTEXT.repo/'utility/chemical_database/chemdb.py'),
                 '--lock',str(lock),'list']
        process=CONTEXT.command(command,'python-cli-list')
        self.assertEqual(process.returncode,0,process.stderr)
        self.assertEqual([r['id'] for r in json.loads(process.stdout)],['fixture-z','fixture-a'])

    def test_duplicate_json_authority_fields_refused(self):
        path=self.fixture.root/'profile.json'
        path.write_text('{"schema":"amps-chemical-profile-v1","id":"fixture-profile","id":"duplicate"}\n')
        self.fixture.fields['profile_sha256']=digest(path);self.fixture.publish()
        with self.assertRaises(self.api.DatabaseError):self.open()

    def test_malformed_lock_shapes_are_typed_refusals(self):
        # Shape-corrupt metadata must not escape as AttributeError/KeyError;
        # callers rely on the declared data-access diagnostic contract.
        path=self.fixture.root/'malformed-lock.json'
        for value in (None,[],{}, {'schema':'amps-chemical-database-lock-v1'}):
            path.write_text(json.dumps(value)+'\n')
            with self.subTest(value=value):
                with self.assertRaises(self.api.DatabaseError):self.api.Database.from_lock(path)

    def test_runner_preserves_failed_subcases_and_statuses(self):
        import run_all_tests as runner
        runner.OUTPUT=CONTEXT.output
        class AdministrativeResultFixtures(unittest.TestCase):
            def test_pass(self):pass
            def test_fail(self):self.fail('DECLARED result-protocol fixture; not a scientific failure')
            def test_error(self):raise RuntimeError('DECLARED result-protocol fixture')
            def test_skip(self):self.skipTest('DECLARED missing prerequisite fixture')
            def test_subcase(self):
                with self.subTest(administrative_case='deliberate failure'):
                    self.fail('DECLARED subcase fixture')
        stream=io.StringIO()
        suite=unittest.defaultTestLoader.loadTestsFromTestCase(AdministrativeResultFixtures)
        result=unittest.TextTestRunner(stream=stream,resultclass=runner.RecordedResult).run(suite)
        (CONTEXT.output/'runner-protocol-control.log').write_text('ADMINISTRATIVE EXPECTED FAILURE FIXTURES ONLY\n'+stream.getvalue())
        self.assertEqual({status:sum(row['status']==status for row in result.records)
                          for status in ('PASS','FAIL','ERROR','SKIPPED')},
                         {'PASS':1,'FAIL':2,'ERROR':1,'SKIPPED':1})


class PythonProduction(unittest.TestCase):
    """Real staged builder/API checks; unavailable implementation is SKIPPED."""
    @classmethod
    def setUpClass(cls):
        cls.database=CONTEXT.python_database()

    def test_actual_api_lists_stable_metadata(self):
        for method in ('datasets','species','reactions'):
            first=getattr(self.database,method)();second=getattr(self.database,method)()
            self.assertEqual(first,second)
        self.assertTrue(self.database.verify(deep=True))

    def test_metadata_views_cannot_mutate_database(self):
        before=repr(self.database.datasets());view=self.database.datasets()
        if isinstance(view,list):view.clear()
        elif isinstance(view,dict):view.clear()
        self.assertEqual(repr(self.database.datasets()),before)

    def test_real_source_table_exact_native_python_agreement(self):
        identifier='iaea-collisiondb-research:e33e73d1967727c64e6096db3df8a57f20af4197639e07fb77ee5062d0ab4af1:table'
        try:table=self.database.table(identifier)
        except CONTEXT.python_api().DatabaseError as error:
            if error.code in ('Unavailable','PermissionRequired'):self.skipTest(str(error))
            raise
        self.assertEqual(table.columns,('energy_eV_per_u','cross_section_cm2'))
        self.assertEqual(table.row_count,12)
        self.assertEqual(table.metadata['sha256'],'2a887ad43497a7f9e2b3d7b6173ae1b2597cc893d897926762c181946c92db93')
        fixture=SimpleNamespace(root=self.database.root,fields=dict(self.database.identity),
                                export_path=self.database.identity['export_path'])
        native=CONTEXT.native(fixture,'table',identifier,'table-access')
        self.assertTrue(native['open_ok'],native);self.assertTrue(native['query_ok'],native)
        self.assertEqual(table.rows,tuple(tuple(row) for row in native['rows']))
        self.assertEqual([row['source_id'] for row in self.database.sources()],native['source_ids'])
        self.assertEqual([row['table_id'] for row in self.database.tables()],native['table_ids'])

    def test_installed_sdk_and_cli_after_relocation(self):
        from verify_deployment import run_checks
        receipt=run_checks(CONTEXT,self.database)
        self.assertEqual(receipt['status'],'PASS')
        self.assertFalse(receipt['physical_activation'])
