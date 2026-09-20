#!/usr/bin/env python3
"""Generate the V05 capability matrix and enforce a named release profile.

Only machine reports count. Missing, skipped, stale, or failed required IDs
remain visible and make the selected profile fail. The tool cannot promote a
capability merely because its source exists, and R8 remains an explicit block
for the scientific-release profile.
"""
import argparse, datetime, json, pathlib, sys
ROOT=pathlib.Path(__file__).resolve().parents[1]
def load(path): return json.loads(path.read_text())
def main():
 p=argparse.ArgumentParser(); p.add_argument('--profile',default='development'); p.add_argument('--test-report',type=pathlib.Path,required=True); p.add_argument('--validation-report',type=pathlib.Path); p.add_argument('--output-dir',type=pathlib.Path,required=True); a=p.parse_args()
 profiles=load(ROOT/'release/profiles.json')['profiles']; caps=load(ROOT/'release/capabilities.json')['capabilities']
 if a.profile not in profiles: raise RuntimeError('unknown release profile')
 required=[]; cursor=a.profile; seen=set()
 while cursor:
  if cursor in seen: raise RuntimeError('release profile inheritance cycle')
  seen.add(cursor); spec=profiles[cursor]; required.extend(spec.get('required',[])); cursor=spec.get('inherits')
 reports=[load(a.test_report)]
 if a.validation_report: reports.append(load(a.validation_report))
 status={}
 for report in reports:
  for row in report.get('results',[]): status[str(row.get('test_id',row.get('id',''))).upper()]=str(row.get('status','ERROR')).upper()
 rows=[{'id':item,'status':status.get(item,'MISSING'),'required':True} for item in dict.fromkeys(required)]
 blocked=profiles[a.profile].get('blocked_by'); passed=not blocked and all(r['status']=='PASS' for r in rows)
 payload={'schema':'srcsep3d-release-evidence-v1','generated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'profile':a.profile,'status':'PASS' if passed else 'INCOMPLETE','blocked_by':blocked,'requirements':rows,'capabilities':caps}
 out=a.output_dir.resolve(); out.mkdir(parents=True,exist_ok=True); (out/'release-evidence.json').write_text(json.dumps(payload,indent=2)+'\n')
 lines=['# srcSEP3D release evidence','',f'- Profile: `{a.profile}`',f'- Status: **{payload["status"]}**',f'- Blocked by: `{blocked}`' if blocked else '- Blocked by: none','','| Evidence | Required | Status |','|---|---:|---|']+[f'| `{r["id"]}` | yes | {r["status"]} |' for r in rows]
 (out/'release-evidence.md').write_text('\n'.join(lines)+'\n')
 return 0 if passed else 1
if __name__=='__main__':
 try: sys.exit(main())
 except Exception as e: print('ERROR:',e,file=sys.stderr); sys.exit(2)
