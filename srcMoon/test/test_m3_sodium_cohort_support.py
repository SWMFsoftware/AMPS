#!/usr/bin/env python3
"""Administrative projection/receipt controls, never physical Na qualification."""
import argparse,hashlib,importlib.util,io,json,tempfile,time,traceback,unittest
from pathlib import Path
import test_m3_normal_cohort as normal
TEST=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('profile_projection',TEST/'sodium_cohort/profile_projection.py');projection=importlib.util.module_from_spec(spec);spec.loader.exec_module(projection)
spec=importlib.util.spec_from_file_location('source_preflight',TEST/'sodium_cohort/prepare_source_preflight.py');preparation=importlib.util.module_from_spec(spec)
import sys
sys.modules['profile_projection']=projection;spec.loader.exec_module(preparation)
class ProjectionControls(unittest.TestCase):
 """Tiny administrative bytes below are not observations or physical states."""
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name)
  self.profile=self.root/'profile.json';self.value={'schema':'srcMoon-M3-state-profile-v1','physical_activation':False,'sources':[]}
  for name in('isotopes','constants','states'):
   p=self.root/(name+'.csv');p.write_text('ADMINISTRATIVE PARSER FIXTURE ONLY\n');self.value[name]={'path':str(p),'sha256':projection.sha(p)}
  for i in range(2):
   p=self.root/('source'+str(i)+'.txt');p.write_text('NOT OBSERVATIONS OR PHYSICAL SOURCE DATA '+str(i));self.value['sources'].append({'id':'ADMIN-'+str(i),'path':str(p),'sha256':projection.sha(p),'role':'administrative-fixture'})
  self.save()
 def save(self):self.profile.write_text(json.dumps(self.value))
 def test_exact_ordered_files_and_hashes(self):
  out=self.root/'native.in';result=projection.project(self.profile,out);text=out.read_text()
  self.assertEqual(result['profile_json_sha256'],projection.sha(self.profile));self.assertEqual(result['sha256'],projection.sha(out))
  self.assertLess(text.index('source: ADMIN-0'),text.index('source: ADMIN-1'))
  for name in('isotopes','constants','states'):self.assertIn(self.value[name]['sha256'],text)
 def test_corruption_is_error_and_no_projection(self):
  Path(self.value['states']['path']).write_text('changed');out=self.root/'bad.in'
  with self.assertRaises(ValueError):projection.project(self.profile,out)
  self.assertFalse(out.exists())
 def test_source_only_profile_cannot_self_activate(self):
  self.value['physical_activation']=True;self.save()
  with self.assertRaises(ValueError):projection.project(self.profile,self.root/'bad.in')
 def test_duplicate_metadata_keys_are_error(self):
  self.profile.write_text('{"schema":"x","schema":"y"}')
  with self.assertRaises(ValueError):projection.project(self.profile,self.root/'bad.in')
 def test_duplicate_source_identity_is_error(self):
  self.value['sources'][1]['id']=self.value['sources'][0]['id'];self.save()
  with self.assertRaises(ValueError):projection.project(self.profile,self.root/'bad.in')
 def test_existing_projection_is_preserved(self):
  out=self.root/'native.in';projection.project(self.profile,out);old=out.read_bytes()
  with self.assertRaises(ValueError):projection.project(self.profile,out)
  self.assertEqual(out.read_bytes(),old)
 def test_hashes_cannot_relocate_another_file(self):
  self.value['sources'][0]['path']=self.value['sources'][1]['path'];self.save()
  with self.assertRaises(ValueError):projection.project(self.profile,self.root/'bad.in')
 def test_preflight_stages_no_missing_physics_choices(self):
  app=self.root/'application.in';app.write_text('#begin application: moon\n#end application\n')
  result=preparation.prepare(self.profile,app,'EXPLICIT-PARENT','EXPLICIT-MISSING-ION','EXPLICIT-ELECTRON',self.root/'request')
  text=Path(result['input']).read_text();self.assertIn('ion_state_id = EXPLICIT-MISSING-ION',text)
  for absent in('maximum_relative_inertial_mass_error','physical_step_s','source_position_lso_m','transport_policy_identity'):self.assertNotIn(absent,text)

def run(output):
 output=Path(output).resolve();output.mkdir(parents=True,exist_ok=False);start=time.monotonic();stream=io.StringIO()
 suite=unittest.TestSuite([unittest.defaultTestLoader.loadTestsFromTestCase(normal.ReceiptControls),unittest.defaultTestLoader.loadTestsFromTestCase(ProjectionControls)])
 try:
  tested=unittest.TextTestRunner(stream=stream,verbosity=2).run(suite);errors=len(tested.errors);failed=len(tested.failures);skipped=len(tested.skipped)
  status='ERROR'if errors else'FAIL'if failed else'PASS';result={'schema':'srcMoon-M3-source-admission-controls-v1','status':status,'passed':tested.testsRun-errors-failed-skipped,'failed':failed,'errors':errors,'skipped':skipped,'checks':tested.testsRun,'scope':'administrative parser/projection/receipt controls only','formal_M3_status':'SKIPPED / NOT QUALIFIED','elapsed_wall_seconds':time.monotonic()-start}
 except Exception:
  stream.write(traceback.format_exc());result={'status':'ERROR','passed':0,'failed':0,'errors':1,'skipped':0,'checks':1}
 result['counts']={'PASS':result['passed'],'FAIL':result['failed'],'ERROR':result['errors'],'SKIPPED':result['skipped']}
 (output/'controls.log').write_text(stream.getvalue());paths=[Path(__file__).resolve(),Path(normal.__file__).resolve(),TEST/'run_m3_normal_cohort.py',TEST/'sodium_cohort/profile_projection.py',TEST/'sodium_cohort/prepare_source_preflight.py']
 (output/'run_manifest.json').write_text(json.dumps({'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},'scope':result.get('scope','infrastructure')},indent=2)+'\n');(output/'result.json').write_text(json.dumps(result,indent=2)+'\n');return result
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',required=True);a=p.parse_args();r=run(a.output_dir);print(json.dumps(r,indent=2));raise SystemExit({'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[r['status']])
