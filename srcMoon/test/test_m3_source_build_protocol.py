#!/usr/bin/env python3
"""Administrative guards for the isolated source build; never run Config.pl.

The real build remains an explicit opt-in workflow. This bounded global-runner
entry tests path/source/mapping/refusal policy using the maintained protocol
controls, and preserves ordinary per-method status/log/source identities.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import time
import traceback
import unittest

HERE=Path(__file__).resolve().parent


class Recorded(unittest.TextTestResult):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs);self.records=[]
    def startTest(self,test):
        super().startTest(test);self.current={'name':test.id(),'status':'PASS'}
    def addFailure(self,test,error):
        super().addFailure(test,error);self.current['status']='FAIL'
    def addError(self,test,error):
        super().addError(test,error);self.current['status']='ERROR'
    def addSkip(self,test,reason):
        super().addSkip(test,reason);self.current.update(status='SKIPPED',reason=reason)
    def addSubTest(self,test,subtest,error):
        super().addSubTest(test,subtest,error)
        if error is not None:self.current['status']='FAIL' if issubclass(error[0],test.failureException)else'ERROR'
    def stopTest(self,test):
        self.records.append(dict(self.current));super().stopTest(test)


def run_control_suite(output,source,dependencies,*,runner,scope,schema):
    """Record an administrative unittest suite with the global-runner protocol.

    ``source`` is an actual maintained test module, not generated test results.
    Every dependency is hashed before and after execution, and unittest owns
    PASS/FAIL/ERROR/SKIPPED classification. This helper performs no Config.pl,
    physical-model selection, particle execution, or qualification promotion.
    Separate callers supply their own schema, runner and explicit limited scope.
    """
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=False)
    started=time.monotonic();rows=[];log=io.StringIO()
    source=Path(source).resolve();runner=Path(runner).resolve()
    paths=[Path(__file__).resolve(),runner,source,*map(Path,dependencies)]
    digests={str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths}
    try:
        spec=importlib.util.spec_from_file_location('recorded_m3_administrative_controls',source)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        result=unittest.TextTestRunner(stream=log,verbosity=2,resultclass=Recorded).run(
            unittest.defaultTestLoader.loadTestsFromModule(module))
        rows=result.records
        if any(hashlib.sha256(Path(p).read_bytes()).hexdigest()!=digest for p,digest in digests.items()):
            raise RuntimeError('protocol/build source changed during verification')
    except Exception as error:
        log.write(traceback.format_exc());rows.append({'name':'protocol-infrastructure','status':'ERROR','reason':str(error)})
    (output/'unittest.log').write_text(log.getvalue())
    counts={s:sum(r['status']==s for r in rows)for s in ('PASS','FAIL','ERROR','SKIPPED')}
    status='ERROR'if counts['ERROR']else'FAIL'if counts['FAIL']else'SKIPPED'if counts['SKIPPED']else'PASS'
    value={'schema':schema,'status':status,'counts':counts,'checks':rows,
           'scope':scope,
           'elapsed_seconds':time.monotonic()-started,'formal_M3_status':'SKIPPED / NOT QUALIFIED'}
    (output/'result.json').write_text(json.dumps(value,indent=2)+'\n')
    (output/'run_manifest.json').write_text(json.dumps({'created_utc':datetime.now(timezone.utc).isoformat(),
        'source_sha256':digests,'command':[sys.executable,'-B',str(runner),'--output-dir',str(output)]},indent=2)+'\n')
    return value


def run(output):
    return run_control_suite(output,HERE/'source_qualification_build/test_protocol.py',
        [HERE/'build_m3_source_qualification.py'],runner=__file__,
        schema='srcMoon-source-build-protocol-controls-v1',
        scope='administrative isolation/provenance guards only; NO build or physical qualification')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output-dir',required=True,type=Path)
    args=parser.parse_args();value=run(args.output_dir);print(json.dumps(value,indent=2))
    raise SystemExit({'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[value['status']])
