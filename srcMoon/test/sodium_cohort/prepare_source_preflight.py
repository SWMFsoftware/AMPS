#!/usr/bin/env python3
"""Stage data-only source inspection without selecting missing physics values."""
import argparse,json
from pathlib import Path
from profile_projection import project,pairs,sha

def prepare(profile,application,parent,ion,electron,output):
 profile=Path(profile).resolve();application=Path(application).resolve();output=Path(output).resolve()
 if output.exists():raise ValueError('output exists; preserve earlier source inspection')
 value=json.loads(profile.read_text(),object_pairs_hook=pairs)
 output.mkdir(parents=True);projection=project(profile,output/'profile.native.in')
 def quote(path):
  text=str(Path(path).resolve())
  if any(c in text for c in '\r\n!"'):raise ValueError('unrepresentable native path')
  return '"'+text+'"'
 def record_path(record):
  p=Path(record['path']);return p if p.is_absolute()else profile.parent/p
 lines=['! Source availability only: numerical/physical transport choices are NOT supplied.','#begin sodium_qualification: source-availability','  schema = srcMoon-native-sodium-cohort-v1','  application_input = '+quote(application),'  profile_metadata_file = '+quote(profile),'  profile_metadata_sha256 = '+sha(profile),'  profile_native_manifest_file = '+quote(projection['path']),'  profile_native_manifest_sha256 = '+projection['sha256']]
 for name in('isotopes','constants','states'):
  record=value[name];lines+=['  '+name+'_file = '+quote(record_path(record)),'  '+name+'_sha256 = '+record['sha256']]
 for name,identifier in(('parent',parent),('ion',ion),('electron',electron)):
  if not identifier or any(c in identifier for c in '\n\r!'):raise ValueError('explicit representable state ID required')
  lines+=['  '+name+'_state_id = '+identifier]
 for source in value['sources']:
  lines+=['  #begin source: '+source['id'],'    path = '+quote(record_path(source)),'    sha256 = '+source['sha256'],'    role = '+source['role'],'  #end source']
 lines+=['#end sodium_qualification','']
 path=output/'source-preflight.in';path.write_text('\n'.join(lines));return {'input':str(path),'sha256':sha(path),'scope':'data-only source availability; no invented error budget, geometry, cohort or physical state'}

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--profile-json',required=True);p.add_argument('--application-input',required=True);p.add_argument('--parent-state-id',required=True);p.add_argument('--ion-state-id',required=True);p.add_argument('--electron-state-id',required=True);p.add_argument('--output-dir',required=True);a=p.parse_args();print(json.dumps(prepare(a.profile_json,a.application_input,a.parent_state_id,a.ion_state_id,a.electron_state_id,a.output_dir),indent=2))
