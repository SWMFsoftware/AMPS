#!/usr/bin/env python3
"""Source-only tests for V03-V05 campaign and governance contracts."""
import json, pathlib, subprocess, sys, tempfile, unittest
ROOT=pathlib.Path(__file__).resolve().parents[1]
class VProgramTests(unittest.TestCase):
 def test_v04_keeps_deferred_r8_blocked(self):
  campaign=json.loads((ROOT/'validation/v04_campaign.json').read_text()); live=campaign['rungs'][-1]
  self.assertEqual(live['blocked_by'],'R8'); self.assertEqual(live['status'],'BLOCKED'); self.assertIn('not accepted',live['reason'])
 def test_native_profiles_are_monotone(self):
  p={x['name']:x for x in json.loads((ROOT/'validation/native_profiles.json').read_text())['profiles']}
  self.assertTrue(set(p['small']['ranks'])<=set(p['medium']['ranks'])<=set(p['production']['ranks']))
 def test_release_generator_rejects_missing_and_accepts_complete_development(self):
  with tempfile.TemporaryDirectory() as td:
   d=pathlib.Path(td); report=d/'tests.json'; required=['V1D01','V1D02','V1D03','V1D04','V1D05','V2D01','V5D01']
   report.write_text(json.dumps({'results':[{'test_id':x,'status':'PASS'} for x in required]}))
   cmd=[sys.executable,str(ROOT/'release/generate_release_evidence.py'),'--profile','development','--test-report',str(report),'--output-dir',str(d/'ok')]
   self.assertEqual(subprocess.run(cmd).returncode,0)
   report.write_text(json.dumps({'results':[{'test_id':'V1D01','status':'PASS'}]}))
   self.assertEqual(subprocess.run(cmd[:-1]+[str(d/'bad')]).returncode,1)
if __name__=='__main__': unittest.main(verbosity=2)
