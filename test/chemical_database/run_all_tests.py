#!/usr/bin/env python3
"""Single shared chemical-database software runner; no application build.

Every method outcome and compile/query log is preserved. A missing production
API/source is SKIPPED, evaluated mismatch FAIL and build/runtime failure ERROR.
Administrative fixtures are explicitly labelled, never validation evidence.
"""
import argparse
from collections import Counter
from datetime import datetime,timezone
import hashlib
import importlib
import inspect
import json
import os
import shlex
from pathlib import Path
import shutil
import subprocess
import sys
import time
import traceback
import unittest

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
import test_protocol as tests


class Context:
    def __init__(self,output,lock,cxx=None):
        self.output=output;self.lock=lock;self.commands=[];self.executable=output/'catalogue-probe'
        self.repo=ROOT;self.python_executable=sys.executable
        self.cxx=cxx or shlex.split(os.environ.get('CXX','g++'))
        self.compile()

    def compile(self):
        compiler=shutil.which(self.cxx[0]) if self.cxx else None
        if not compiler:raise unittest.SkipTest('native C++ compiler unavailable')
        base=ROOT/'src/models/chemical_database'
        self.sources=[HERE/'probe.cpp',HERE/'test_protocol.py',HERE/'verify_deployment.py',Path(__file__).resolve(),base/'Runtime.cpp',base/'Runtime.h',
            base/'DataIO.cpp',base/'DataIO.h',ROOT/'src/general/input/amps_input.h']
        self.sources += [ROOT/'srcMoon/species'/(name+'.h') for name in
                         ('Photochemistry','ElectronImpact','ChargeExchange')]
        self.sources += [base/(name+'.h') for name in ('Photochemistry','ElectronImpact','ChargeExchange')]
        # New mathematical helpers and their thin executable adapters are
        # shared production code too, even though no AMPS particle build is
        # needed. Bind their actual bytes in the same campaign identity.
        self.sources += [base/'ReducedModels.h']
        self.sources += sorted((base/'examples').glob('*.cpp'))
        self.sources += sorted(HERE.glob('reduced_models_probe.cpp'))
        # Python queries are production paths too. Bind their actual code,
        # including the CLI imported by a subprocess, rather than identifying
        # only the native half of a cross-language result.
        self.sources += sorted((ROOT/'utility/chemical_database/python/amps_chemical_database').glob('*.py'))
        self.sources += [ROOT/'utility/chemical_database/chemdb.py']
        # Acquisition/import and native projection code is part of the source
        # identity even though it is Python rather than the compiled reader.
        # New source-parser tests are discovered below, with honest SKIPs for
        # explicitly unavailable authoritative packets in another checkout.
        self.sources += sorted((ROOT/'utility/chemical_database').glob('*.py'))
        self.sources += sorted(HERE.glob('test_*.py'))
        self.sources = list(dict.fromkeys(self.sources))
        missing=[str(p) for p in self.sources if not p.is_file()]
        if missing:raise unittest.SkipTest('shared runtime implementation unavailable: '+', '.join(missing))
        compiler_command=[compiler]+self.cxx[1:]
        command=compiler_command+['-std=c++17','-O2','-g','-Wall','-Wextra','-Werror','-pedantic','-pthread',
            str(HERE/'probe.cpp'),str(base/'Runtime.cpp'),str(base/'DataIO.cpp'),'-o',str(self.executable)]
        (self.output/'compile_command.json').write_text(json.dumps(command,indent=2)+'\n')
        process=subprocess.run(command,cwd=ROOT,text=True,capture_output=True,timeout=180)
        (self.output/'compile.stdout.log').write_text(process.stdout);(self.output/'compile.stderr.log').write_text(process.stderr)
        self.commands.append({'command':command,'returncode':process.returncode,'log':str(self.output/'compile.stderr.log')})
        if process.returncode:raise RuntimeError('actual shared Catalogue compilation failed; see compile.stderr.log')
        self.build={'source_revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
            'compiler':subprocess.check_output(compiler_command+['--version'],text=True).splitlines()[0],
            'source_sha256':{str(p.relative_to(ROOT)):tests.digest(p) for p in self.sources},
            'executable_sha256':tests.digest(self.executable),'command':command,
            'MPI':'not used; PIC-independent C++17 API','OMP_NUM_THREADS':os.environ.get('OMP_NUM_THREADS','unspecified')}
        (self.output/'build.json').write_text(json.dumps(self.build,indent=2,sort_keys=True)+'\n')

    def native(self,fixture,operation,identifier,capability,**changes):
        directory=self.output/'queries'/('query-%04d'%len(self.commands));directory.mkdir(parents=True)
        fields=fixture.fields
        export_relative=getattr(fixture,'export_path','export.in')
        command=[str(self.executable),str(fixture.root),export_relative,changes.get('expected_export',tests.digest(fixture.root/export_relative)),
            fields['release_id'],changes.get('expected_release',fields['release_manifest_sha256']),
            fields['profile_id'],changes.get('expected_profile',fields['profile_sha256']),capability or 'none',operation,identifier,
            str(changes.get('file_limit',0)),str(changes.get('total_limit',0))]
        process=subprocess.run(command,cwd=directory,text=True,capture_output=True,timeout=20)
        (directory/'stdout.log').write_text(process.stdout);(directory/'stderr.log').write_text(process.stderr)
        self.commands.append({'command':command,'returncode':process.returncode,'log':str(directory/'stdout.log')})
        if process.returncode:raise RuntimeError('native API adapter crashed or returned infrastructure error: '+str(directory))
        result=json.loads(process.stdout)
        (directory/'result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
        return result

    def python_database(self):
        api=self.python_api()
        if self.lock is None:raise unittest.SkipTest('real staged production lock not supplied; no science release fabricated')
        return api.Database.from_lock(self.lock)

    def python_api(self):
        package=ROOT/'utility/chemical_database/python'
        if not (package/'amps_chemical_database').exists():raise unittest.SkipTest('production Python package not available yet')
        sys.path.insert(0,str(package))
        api=importlib.import_module('amps_chemical_database')
        return api

    def command(self,command,label,timeout=30):
        directory=self.output/'queries'/('query-%04d-%s'%(len(self.commands),label));directory.mkdir(parents=True)
        process=subprocess.run(command,cwd=directory,text=True,capture_output=True,timeout=timeout)
        (directory/'stdout.log').write_text(process.stdout);(directory/'stderr.log').write_text(process.stderr)
        self.commands.append({'command':command,'returncode':process.returncode,'log':str(directory/'stdout.log')})
        return process


class RecordedResult(unittest.TextTestResult):
    def __init__(self,*args,**kwargs):super().__init__(*args,**kwargs);self.records=[]
    def addSuccess(self,test):super().addSuccess(test);self.record(test,'PASS')
    def addFailure(self,test,error):super().addFailure(test,error);self.record(test,'FAIL')
    def addError(self,test,error):super().addError(test,error);self.record(test,'ERROR')
    def addSkip(self,test,reason):super().addSkip(test,reason);self.record(test,'SKIPPED',reason)
    def addSubTest(self,test,subtest,error):
        super().addSubTest(test,subtest,error)
        # unittest does not call addFailure/addError for a failed subTest.
        # Preserve each failed subcase explicitly; otherwise its method could
        # disappear from our aggregate while the text log correctly failed.
        if error is not None:
            self.record(subtest,'FAIL' if issubclass(error[0],test.failureException) else 'ERROR')
    def record(self,test,status,reason=None):
        owner = getattr(test, 'test_case', test)
        runner = inspect.getsourcefile(owner.__class__) or str(HERE/'test_protocol.py')
        self.records.append({'test':test.id(),'status':status,'reason':reason,'runner':str(runner),
            'log':str(OUTPUT/'unittest.log')})


OUTPUT=None


def select_additional_tests(loaded, protocol_only, excluded):
    """Select the declared software scope before execution, never by outcome.

    CORE's protocol case promises no external-data prerequisite. Real-source
    classes remain registered for the full runner, where missing originals
    report SKIPPED. Omitting them from protocol scope is explicit in its receipt
    and does not convert an attempted/skipped physical test into PASS.
    """
    selected = unittest.TestSuite()
    for item in loaded:
        if isinstance(item, unittest.TestSuite):
            selected.addTests(select_additional_tests(item, protocol_only, excluded))
        elif protocol_only and getattr(type(item), 'requires_authoritative_source', False):
            excluded.append(item.id())
        else:
            selected.addTest(item)
    return selected


def select_owner_suites(owners, protocol_only, excluded):
    """Load the protocol owner classes for the declared scope, never by outcome.

    ``owners`` is the ordered (NativeProtocol, PythonProtocol, PythonProduction)
    triple. PythonProduction needs a real staged/installed release
    (--database-lock), so --protocol-only does not run it. Its method IDs are
    recorded in ``excluded`` (the receipt's outside_selected_scope), like the
    authoritative-source classes, so the receipt names every registered method
    that this scope leaves out.
    """
    native, python_protocol, production = owners
    loader = unittest.defaultTestLoader
    suite = unittest.TestSuite(loader.loadTestsFromTestCase(owner) for owner in (native, python_protocol))
    production_tests = loader.loadTestsFromTestCase(production)
    if protocol_only:
        excluded.extend(test.id() for test in production_tests)
    else:
        suite.addTest(production_tests)
    return suite


def main():
    global OUTPUT
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--database-lock',type=Path,help='real staged/installed pinned lock for additional Python API controls')
    parser.add_argument('--protocol-only',action='store_true',help='select only site-independent administrative/native/Python API controls')
    parser.add_argument('--cxx-json',help='explicit compiler argv JSON for the global core scheduler, preserving wrappers/flags')
    args=parser.parse_args();OUTPUT=args.output_dir.resolve()
    cxx=json.loads(args.cxx_json) if args.cxx_json else None
    if cxx is not None and (not isinstance(cxx,list) or not cxx or any(not isinstance(item,str) or not item for item in cxx)):
        parser.error('--cxx-json must be a nonempty JSON array of nonempty strings')
    if OUTPUT.exists():parser.error('output directory must be new to preserve prior evidence')
    OUTPUT.mkdir(parents=True)
    started=time.monotonic();records=[];context=None;excluded=[]
    try:
        context=Context(OUTPUT,args.database_lock,cxx)
        tests.CONTEXT=context
        with (OUTPUT/'unittest.log').open('w') as log:
            suite=select_owner_suites((tests.NativeProtocol,tests.PythonProtocol,tests.PythonProduction),
                                      args.protocol_only,excluded)
            # These registered parser/source-reproduction controls extend the
            # existing shared database case; they do not allocate a lunar U/I ID
            # or make missing source files/physical providers pass.
            for path in sorted(HERE.glob('test_*.py')):
                if path.name == 'test_protocol.py':
                    continue
                module = importlib.import_module(path.stem)
                suite.addTests(select_additional_tests(
                    unittest.defaultTestLoader.loadTestsFromModule(module), args.protocol_only, excluded))
            result=unittest.TextTestRunner(stream=log,verbosity=2,resultclass=RecordedResult).run(suite)
        records=result.records
        changed=[str(p) for p in context.sources if tests.digest(p)!=context.build['source_sha256'][str(p.relative_to(ROOT))]]
        if changed:
            records.append({'test':'immutable-source-basis','status':'ERROR','reason':'source changed during run: '+', '.join(changed),
                'runner':str(Path(__file__).resolve()),'log':str(OUTPUT/'build.json')})
    except unittest.SkipTest as error:
        records=[{'test':'native-build-prerequisite','status':'SKIPPED','reason':str(error),'runner':str(Path(__file__)), 'log':str(OUTPUT/'compile.stderr.log')}]
    except Exception as error:
        (OUTPUT/'runner-error.log').write_text(traceback.format_exc())
        records=[{'test':'native-build-or-runner','status':'ERROR','reason':str(error),'runner':str(Path(__file__)), 'log':str(OUTPUT/'runner-error.log')}]
    counts={status:sum(row['status']==status for row in records) for status in ('PASS','FAIL','SKIPPED','ERROR')}
    status='ERROR' if counts['ERROR'] else 'FAIL' if counts['FAIL'] else 'PASS' if counts['PASS'] else 'SKIPPED'
    summary={'schema':'amps-chemical-database-software-tests-v1','status':status,
        'scope':'production data-access/integrity software verification; NOT physical validation',
        'protocol_only':args.protocol_only,
        'outside_selected_scope':excluded,
        'created_utc':datetime.now(timezone.utc).isoformat(),'executed':len(records),'counts':counts,'tests':records,
        'failed_or_error_tests':[row for row in records if row['status'] in ('FAIL','ERROR')],
        'wall_seconds':time.monotonic()-started,'commands':context.commands if context else []}
    (OUTPUT/'result.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    print('Shared chemical database test summary\n  executed: %d'%len(records))
    for name,value in counts.items():print('  %-8s %d'%(name.lower()+':',value))
    for row in summary['failed_or_error_tests']:
        print('%s %s\n  runner: %s\n  log: %s'%(row['status'],row['test'],row['runner'],row['log']))
    print('Aggregate result: '+str(OUTPUT/'result.json'))
    return {'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[status]


if __name__=='__main__':raise SystemExit(main())
