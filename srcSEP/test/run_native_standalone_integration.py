#!/usr/bin/env python3
"""Qualify stages 2–4 against a matching, configured real-MPI AMPS archive.

A configured build is an explicit prerequisite; this test never reconfigures
another checkout. Every requested rank count must execute every diagnostic
mode. Missing native output/MPI startup is ERROR, never PASS or optional SKIP.
Source, generated objects, archive and executable identities accompany results.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import shlex
import subprocess
import time
import tempfile
import xml.etree.ElementTree as ET

REPO = Path(__file__).resolve().parents[2]
MODES = {'geometry': ['SINT08', 'SINT09'], 'connectivity': ['SINT13N'],
         'transport': ['SINT14N', 'SINT15N', 'SINT16N']}
MODES.update({'source-parker':['SINT25N'],'source-fte-dmumu':['SINT26N'],
              'source-fte-mfp':['SINT27N'],'source-mono':['SINT28N'],
              'source-zero':['SINT29N'],'source-invalid':['SINT30N'],
              'source-disconnected':['SINT31N'],'source-invalid-energy':['SINT32N']})
MODES['source-merge']=['MCP08']
IDS = [identity for identities in MODES.values() for identity in identities]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pic_dependency_closure(build):
    """Reject stale PIC template instantiations outside the selected movers.

    An old object can supply an inline/weak pool or particle-buffer method at
    link time even when every selected mover was rebuilt. Inspect every MMD
    closure that consumes PIC headers, including application/model objects.
    Dependency paths are resolved in the compiler's working directory; nested
    util/foo.d files commonly contain targets relative to build/main.
    """
    generated = build / "build"
    paths = []
    for depfile in sorted(generated.rglob("*.d")):
        rule = depfile.read_text().split("\n\n", 1)[0].replace("\\\n", " ")
        if ":" not in rule:
            continue
        target, text = rule.split(":", 1)
        tokens = shlex.split(text)
        if not any(Path(token).name in ("pic.h", "picGlobal.dfn") for token in tokens):
            continue
        first = depfile.relative_to(generated).parts[0]
        candidates = [generated / first, depfile.parent]
        bases = [base for base in candidates if tokens and (base / tokens[0]).is_file()]
        if not bases:
            raise ValueError("cannot resolve native compiler working directory: " + str(depfile))
        base = bases[0]
        targets = shlex.split(target)
        if len(targets) != 1:
            raise ValueError("ambiguous native compiler target: " + str(depfile))
        obj = (base / targets[0]).resolve()
        dependencies = [(base / token).resolve() for token in tokens]
        if not obj.is_file() or any(not path.is_file() or path.stat().st_mtime_ns > obj.stat().st_mtime_ns
                                    for path in dependencies):
            raise ValueError("stale native PIC template/header closure: " + str(obj))
        # A copied build must not silently compile consumers against another
        # configured PIC tree. Model headers outside this build can be valid;
        # the actual pic.h/picGlobal.dfn authority must be this generated tree.
        for path in dependencies:
            if path.name in ("pic.h", "picGlobal.dfn") and path.parent != generated / "pic":
                raise ValueError("native PIC dependency belongs to another configured build: " + str(path))
        paths.extend([depfile, obj] + dependencies)
    if not paths:
        raise ValueError("no native PIC compiler dependency evidence")
    return paths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-root', type=Path, default=os.environ.get('AMPS_FIELD_LINE_NATIVE_BUILD'))
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--ranks', type=int, nargs='+', default=[1, 4])
    parser.add_argument('--prepare-only', action='store_true', help='compile/hash fixtures and write explicit launch jobs without claiming test results')
    parser.add_argument('--collect', action='store_true', help='validate prepared hashes and collect externally launched jobs with captured process return codes')
    args = parser.parse_args()
    if args.prepare_only and args.collect:parser.error('--prepare-only and --collect are mutually exclusive')
    reports=Path(os.environ['SRCSEP_REPORT_DIR']).resolve() if os.environ.get('SRCSEP_REPORT_DIR') else None
    if reports:reports.mkdir(parents=True,exist_ok=True)
    output=args.output_dir.resolve() if args.output_dir else Path(tempfile.mkdtemp(
        prefix='standalone-native-',dir=str(reports) if reports else None))
    output.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    commands, runs, hashes = [], [], {}
    rows = []
    try:
        if not args.build_root or any(n < 1 for n in args.ranks) or len(set(args.ranks)) != len(args.ranks):
            raise ValueError('supply a configured --build-root and distinct positive MPI rank counts')
        build = Path(args.build_root).resolve()
        archive = build/'build/libAMPS.a'
        paths = [archive, build/'Makefile',build/'Makefile.conf',build/'Makefile.def',build/'Makefile.local',
                 build/'build/main/makefile',build/'build/pic/makefile',
                 REPO/'srcSEP/test/native_standalone_integration.cpp',
                 REPO/'srcSEP/test/native_standalone_integration.mk', Path(__file__),
                 REPO/'input/test/sep_reduced_front__field_line.input',
                 REPO/'input/test/sep_reduced_front_protons.input',
                 REPO/'srcSEP/examples/reduced-shock/field_line.in', REPO/'srcSEP/examples/reduced-shock/positive_1au.event']
        paths += pic_dependency_closure(build)
        # Compare the actual generated consumer with this checkout. Ignore only
        # the repository's explicit generated macro-check comments, never code.
        def normalized(path):
            value=re.sub(r'/\* AMPS_CHECKMACRO_BEGIN \*/.*?/\* AMPS_CHECKMACRO_END \*/', '', path.read_text(), flags=re.S)
            # The maintained test profile expands this console-only prefix.
            # No numerical/configuration token is normalized; both original
            # files are hashed below, including their actual emitted prefix.
            return re.sub(r'\s+', '', value.replace('$PREFIX:', 'AMPS:'))
        for name in ['sep', 'physics_checkpoint','util/sep_population_production_adapter','util/sep_conservative_resampling',
                     'reduced_shock_geometry', 'transport_common', 'parker_mover',
                     'reduced_shock_particles','reduced_shock_boundary','util/sep_reduced_source','util/sep_absorbing_bridge',
                     'focused_transport_dmumu', 'focused_transport_mfp', 'main_lib',
                     'adapters/reduced_shock_background_adapter', 'util/sep_reduced_geometry', 'util/sep_parker_core',
                     'util/sep_focused_transport_core', 'util/sep_focused_transport_mfp_core']:
            source, generated = REPO/'srcSEP'/(name+'.cpp'), build/'build/main'/(name+'.cpp')
            obj = generated.with_suffix('.o')
            if normalized(source) != normalized(generated):
                raise ValueError('native generated source differs from checkout: '+name)
            if not obj.is_file() or obj.stat().st_mtime_ns < generated.stat().st_mtime_ns or archive.stat().st_mtime_ns < obj.stat().st_mtime_ns:
                raise ValueError('native object/archive is stale: '+name)
            paths += [source, generated, obj]
            # Follow the compiler's actual MMD closure rather than assuming
            # cpp timestamps prove an ABI/header update was rebuilt. Only the
            # first dependency rule is physical input; later -MP rules are
            # empty make convenience targets, not additional source files.
            dependency_file=obj.with_suffix('.d')
            dependency_text=dependency_file.read_text().split('\n\n',1)[0].replace('\\\n',' ')
            dependencies=[(build/'build/main'/name).resolve() for name in shlex.split(dependency_text.split(':',1)[1])]
            if any(not path.is_file() or path.stat().st_mtime_ns>obj.stat().st_mtime_ns for path in dependencies):
                raise ValueError('native header/dependency closure is stale: '+name)
            paths += [dependency_file]+dependencies
        # sep.h legitimately contains post-configuration Cartesian force flags.
        # Compare the new public output contract verbatim rather than removing
        # or canonicalizing those physical flags. The full generated header is
        # already frozen in the compiler dependency closure, and is hashed with
        # its source here; a change during execution remains an error.
        source,generated=REPO/'srcSEP/sep.h',build/'build/main/sep.h'
        contract=r'struct ShockDiagnostic\s*\{.*?QueryShockDiagnostic\([^;]*;'
        a,b=re.search(contract,source.read_text(),re.S),re.search(contract,generated.read_text(),re.S)
        if not a or not b or re.sub(r'\s+','',a.group())!=re.sub(r'\s+','',b.group()):
            raise ValueError('native generated shock-output interface differs from checkout')
        paths += [source,generated]
        for name in ['physics_checkpoint.h','transport_common.h','reduced_shock_geometry.h','adapters/reduced_shock_background_adapter.h',
                     'reduced_shock_particles.h','reduced_shock_boundary.h','util/sep_reduced_source.h',
                     'util/sep_reduced_particle_options.h','util/sep_absorbing_bridge.h',
                     'util/sep_reduced_geometry.h','util/sep_parker_core.h','util/sep_focused_transport_core.h',
                     'util/sep_focused_transport_mfp_core.h']:
            source,generated=REPO/'srcSEP'/name,build/'build/main'/name
            if normalized(source)!=normalized(generated):
                raise ValueError('native generated interface differs from checkout: '+name)
            paths += [source,generated]
        for relative in ['src/models/sep_common/sep_background_snapshot.cpp',
                         'src/models/sep_corona_swcme/src/ambient_state.cpp',
                         'src/models/sep_corona_swcme/include/sep_corona_swcme/ambient_state.h',
                         'src/models/sep_coronal_cme/src/particle_source.cpp',
                         'src/models/sep_coronal_cme/src/mhd_jump_solver.cpp',
                         'src/models/sep_coronal_cme/include/sep_coronal_cme/particle_source.h',
                         'src/models/sep_corona_swcme/shock_front/provider.cpp']:
            if (REPO/relative).read_bytes() != (build/relative).read_bytes():
                raise ValueError('native core physics source differs: '+relative)
            paths += [REPO/relative, build/relative]
        hashes = {str(path): digest(path) for path in paths}
        executable = output/'native-integration'
        environment = dict(os.environ, AMPS_INTEGRATION_REPO=str(REPO), AMPS_INTEGRATION_OUTPUT=str(executable))
        argv = ['make', '-f', str(REPO/'srcSEP/test/native_standalone_integration.mk'), 'native-standalone-integration']
        commands.append({'argv': argv, 'cwd': str(build), 'log': str(output/'build.log')})
        if not args.collect:
            with (output/'build.log').open('w') as stream:
                compiled = subprocess.run(argv, cwd=build, env=environment, stdout=stream, stderr=subprocess.STDOUT, timeout=300)
            if compiled.returncode:
                raise RuntimeError('native fixture compile/link failed; see build.log')
        hashes[str(executable)] = digest(executable)
        if args.collect:
            prepared=json.loads((output/'prepared.json').read_text())
            if prepared['hashes']!=hashes or prepared['ranks']!=args.ranks:
                raise RuntimeError('prepared source/build identity or requested MPI ranks changed before collection')
        if args.prepare_only:
            # Some execution environments allow an explicit MPI launcher but
            # prohibit socket creation in Python children. Preparation is not
            # a test PASS. The external launcher must record the REAL return
            # code of every job; collection verifies all source/binary hashes
            # again and requires the ordinary complete structured case output.
            jobs=[]
            for ranks in args.ranks:
                for mode in MODES:
                    directory=output/(mode+'-ranks'+str(ranks));directory.mkdir()
                    jobs.append({'ranks':ranks,'mode':mode,'cwd':str(directory),'log':str(directory/'native.log'),
                        'argv':['mpiexec','-n',str(ranks),str(executable),
                            str(REPO/'srcSEP/examples/reduced-shock/field_line.in'),
                            str(REPO/'srcSEP/examples/reduced-shock/positive_1au.event'),mode]})
            (output/'prepared.json').write_text(json.dumps({'hashes':hashes,'ranks':args.ranks,'jobs':jobs},indent=2)+'\n')
            print('Prepared native MPI jobs (not executed): '+str(output/'prepared.json'),flush=True)
            return 0
        for ranks in args.ranks:
            for mode, expected in MODES.items():
                directory = output/(mode+'-ranks'+str(ranks))
                directory.mkdir(exist_ok=args.collect)
                argv = ['mpiexec', '-n', str(ranks), str(executable),
                        str(REPO/'srcSEP/examples/reduced-shock/field_line.in'),
                        str(REPO/'srcSEP/examples/reduced-shock/positive_1au.event'), mode]
                log = directory/'native.log'
                commands.append({'argv': argv, 'cwd': str(directory), 'log': str(log)})
                run_started=time.monotonic()
                if args.collect:
                    captured=json.loads((directory/'process-result.json').read_text())
                    code=captured.get('exit_code')
                    if type(code) is not int:raise RuntimeError('external MPI job has no captured completed process return code: '+mode)
                    result=subprocess.CompletedProcess(argv,code)
                else:
                    with log.open('w') as stream:
                        result = subprocess.run(argv, cwd=directory, stdout=stream, stderr=subprocess.STDOUT, timeout=180)
                # Stage-7 merger cases use MCP IDs; retain the exact expected
                # per-mode list check below so an unrelated PASS cannot qualify
                # the requested native fixture or conceal a missing result.
                records = re.findall(r'^\[((?:SINT\d+N?)|(?:MCP\d+))\] (PASS|FAIL) (.*)$', log.read_text(), re.M)
                status = 'ERROR' if [r[0] for r in records] != expected else 'FAIL' if result.returncode or any(r[1] == 'FAIL' for r in records) else 'PASS'
                runs.append({'mode': mode, 'ranks': ranks, 'returncode': result.returncode,
                             'status': status, 'records': records, 'log': str(log),
                             'seconds':captured.get('wall_time_seconds',0) if args.collect else time.monotonic()-run_started})
        for path, before in hashes.items():
            if digest(Path(path)) != before:
                raise RuntimeError('native source/build changed during qualification: '+path)
        for mode, identities in MODES.items():
            related = [run for run in runs if run['mode'] == mode]
            status = 'ERROR' if any(r['status'] == 'ERROR' for r in related) else 'FAIL' if any(r['status'] == 'FAIL' for r in related) else 'PASS'
            rows += [{'test_id': identity, 'status': status, 'message': 'actual native '+mode+'; MPI ranks '+str(args.ranks),
                      'seconds':sum(r['seconds'] for r in related),'seed':0,
                      'configuration':{'ranks':args.ranks,'mode':mode,'physical_dt':'unchanged runtime input','units':'SI/HCI'},
                      'acceptance':'Frozen assertions in hashed native_standalone_integration.cpp; every requested MPI mode/rank must execute',
                      'metrics':[line for r in related for line in Path(r['log']).read_text().splitlines() if '_METRIC ' in line or line.startswith('REDUCED_')],
                      'artifacts':[r['log'] for r in related]} for identity in identities]
    except Exception as error:
        rows = [{'test_id': identity, 'status': 'ERROR', 'message': str(error)} for identity in IDS]
    # Structured evidence exists even for missing build or startup failure. The
    # shared registries import these individual rows and preserve ERROR/FAIL.
    for row in rows:
        print(f"[{row['test_id']}] {row['status']} {row['message']}", flush=True)
    seconds = time.monotonic()-started
    (output/'commands.json').write_text(json.dumps(commands, indent=2)+'\n')
    (output/'source-build-sha256.json').write_text(json.dumps(hashes, indent=2)+'\n')
    (output/'summary.json').write_text(json.dumps({'schema': 'srcsep-native-standalone-integration-v1',
        'scope': 'real MPI/native; actual reduced net source plus labelled manufactured cases; no downstream volume or moving restart qualification',
        'seconds': seconds, 'runs': runs, 'results': rows}, indent=2)+'\n')
    suite = ET.Element('testsuite', name='native-standalone-integration', tests=str(len(rows)), time=str(seconds))
    for row in rows:
        node = ET.SubElement(suite, 'testcase', name=row['test_id'])
        if row['status'] != 'PASS':
            ET.SubElement(node, 'error' if row['status'] == 'ERROR' else 'failure', message=row['message'])
    ET.ElementTree(suite).write(output/'junit.xml', encoding='unicode')
    if reports:
        # The existing srcSEP source-suite merger expects id/results records,
        # and scans only top-level JSON. Keep command/hash/run details in the
        # owned child directory and publish this sole compatibility report.
        imported=[dict(row,id=row['test_id'],group='standalone-core-native',artifacts=[str(output/'summary.json')]) for row in rows]
        (reports/'standalone-native-results.json').write_text(json.dumps({'schema':'srcsep-component-tests-v1','results':imported},indent=2)+'\n')
    return int(any(row['status'] != 'PASS' for row in rows))


if __name__ == '__main__':
    raise SystemExit(main())
