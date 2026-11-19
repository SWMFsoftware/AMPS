#!/usr/bin/env python3
"""Bounded independent event-bookkeeping controls; no physical activation."""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import hashlib
import io
import json
from pathlib import Path
import sys
import time
import traceback
import unittest

ROOT=Path(__file__).resolve().parents[2]
COMPONENT=Path(__file__).with_name('sodium_cohort')
sys.path.insert(0,str(COMPONENT))
import test_ledger
from verify_ledger import verify,write_artifacts,digest

class Recorded(unittest.TextTestResult):
    """Keep one result per executed unittest method, including subtest failures."""
    def startTest(self,test):
        super().startTest(test);self.rows=getattr(self,'rows',[]);self.current={'name':test.id(),'status':'PASS'}
    def addFailure(self,test,error):super().addFailure(test,error);self.current['status']='FAIL'
    def addError(self,test,error):super().addError(test,error);self.current['status']='ERROR'
    def addSkip(self,test,reason):super().addSkip(test,reason);self.current.update(status='SKIPPED',reason=reason)
    def addSubTest(self,test,subtest,error):
        super().addSubTest(test,subtest,error)
        if error is not None:self.current['status']='FAIL' if issubclass(error[0],test.failureException) else 'ERROR'
    def stopTest(self,test):self.rows.append(self.current.copy());super().stopTest(test)


def run(output):
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    rows=[];manifest={'created_utc':datetime.now(timezone.utc).isoformat(),'scope':'independent verifier software + explicitly analytical toy ledger; NO physical source activation','command':[sys.executable,'-B',str(Path(__file__).resolve()),'--output-dir',str(output)],'cwd':str(Path.cwd()),'source_sha256':{}}
    files=[Path(__file__),COMPONENT/'verify_ledger.py',COMPONENT/'test_ledger.py',COMPONENT/'LEDGER_CONTRACT.md',COMPONENT/'LEDGER_README.md']
    manifest['source_sha256']={str(p.resolve()):digest(p) for p in files}
    try:
        stream=io.StringIO();suite=unittest.defaultTestLoader.loadTestsFromModule(test_ledger)
        actual=unittest.TextTestRunner(stream=stream,verbosity=2,resultclass=Recorded).run(suite)
        (output/'unittest.log').write_text(stream.getvalue());rows=list(actual.rows)
        # Retain one complete, clearly labelled arithmetic fixture and its
        # independently recomputed CSV/figures. It is never called sodium data.
        fixture=output/'analytical-software-fixture';fixture.mkdir()
        profile,mapping,bounds,values=test_ledger.fixture(fixture)
        events=test_ledger.write_events(fixture,values)
        (fixture/'mapping.json').write_text(json.dumps(mapping,indent=2)+'\n')
        (fixture/'bounds.json').write_text(json.dumps(bounds,indent=2)+'\n')
        diagnostic=verify([events],profile,mapping,bounds);write_artifacts(diagnostic,fixture/'independent-reconstruction')
        rows.append({'name':'retained-analytical-ledger-artifacts','status':diagnostic['status']})
        if any(digest(path)!=expected for path,expected in manifest['source_sha256'].items()):raise RuntimeError('verifier source changed during run')
    except ImportError as error:rows.append({'name':'artifact-dependency','status':'SKIPPED','reason':str(error)})
    except Exception as error:
        rows.append({'name':'verifier-infrastructure','status':'ERROR','reason':str(error)});(output/'exception.log').write_text(traceback.format_exc())
    counts={s:sum(r['status']==s for r in rows) for s in ('PASS','FAIL','ERROR','SKIPPED')}
    status=next(s for s in ('ERROR','FAIL','SKIPPED','PASS') if counts[s])
    result={'schema':'srcMoon-M3-independent-ledger-tests-v1','status':status,'counts':counts,'results':rows,'scope':manifest['scope'],'formal_I21_I24_status':'SKIPPED / NOT QUALIFIED','elapsed_wall_seconds':time.monotonic()-started}
    manifest['artifacts']={str(p.relative_to(output)):digest(p) for p in output.rglob('*') if p.is_file()}
    (output/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');(output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output-dir',required=True,type=Path)
    args=parser.parse_args()
    if args.output_dir.exists():parser.error('output must not exist; preserve earlier evidence')
    r=run(args.output_dir);print(json.dumps(r,indent=2));raise SystemExit({'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[r['status']])
