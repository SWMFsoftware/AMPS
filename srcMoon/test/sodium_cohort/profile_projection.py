#!/usr/bin/env python3
"""Deterministic administrative projection, not state/energy invention.

Native production code authenticates this complete projection and compares its
ordered files/sources against the explicit Request. The original JSON remains
hash-bound too. No READY/activation field is generated or interpreted.
"""
import hashlib,json
from pathlib import Path

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def pairs(items):
 value={}
 for key,item in items:
  if key in value:raise ValueError('duplicate profile JSON key: '+key)
  value[key]=item
 return value

def project(profile_path,output):
 profile_path=Path(profile_path).resolve();output=Path(output)
 value=json.loads(profile_path.read_text(),object_pairs_hook=pairs,parse_constant=lambda x:(_ for _ in ()).throw(ValueError('nonfinite JSON token '+x)))
 if not isinstance(value,dict)or value.get('schema')!='srcMoon-M3-state-profile-v1' or value.get('physical_activation')is not False:raise ValueError('requires actual source-only state profile schema')
 def pin(record):
  if not isinstance(record,dict)or not isinstance(record.get('path'),str)or not isinstance(record.get('sha256'),str):raise ValueError('profile file requires path/SHA')
  path=Path(record['path']);path=path if path.is_absolute()else profile_path.parent/path;path=path.resolve()
  if sha(path)!=record['sha256']:raise ValueError('profile byte identity mismatch: '+str(path))
  if any(c in str(path)for c in '\r\n!"'):raise ValueError('path cannot be represented without changing AMPS grammar')
  return path,record['sha256']
 lines=['#begin physical_state_profile: native','  schema = srcMoon-M3-state-profile-native-v1','  profile_json_file = "'+str(profile_path)+'"','  profile_json_sha256 = '+sha(profile_path)]
 for key in('isotopes','constants','states'):
  path,digest=pin(value[key]);lines+=['  '+key+'_file = "'+str(path)+'"','  '+key+'_sha256 = '+digest]
 sources=value.get('sources');seen=set()
 if not isinstance(sources,list)or not sources:raise ValueError('ordered sources required')
 for item in sources:
  if not isinstance(item,dict)or set(item)!={'id','path','sha256','role'}:raise ValueError('exact source shape required')
  identifier=item['id'];role=item['role']
  if not isinstance(identifier,str)or not identifier or any(not(c.isascii()and(c.isalnum()or c in'_-'))for c in identifier)or identifier in seen:raise ValueError('unique AMPS source name required')
  if not isinstance(role,str)or not role or any(c in role for c in '\r\n!'):raise ValueError('unrepresentable source role')
  seen.add(identifier);path,digest=pin(item)
  lines+=['  #begin source: '+identifier,'    path = "'+str(path)+'"','    sha256 = '+digest,'    role = '+role,'  #end source']
 lines+=['#end physical_state_profile','']
 if output.exists():raise ValueError('projection output exists; preserve earlier artifact')
 output.write_text('\n'.join(lines))
 return {'path':str(output.resolve()),'sha256':sha(output),'profile_json_path':str(profile_path),'profile_json_sha256':sha(profile_path),'scope':'administrative exact projection only; not physical activation'}

if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('profile_json');p.add_argument('output_in');a=p.parse_args();print(json.dumps(project(a.profile_json,a.output_in),indent=2))
