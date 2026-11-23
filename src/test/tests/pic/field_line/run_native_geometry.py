#!/usr/bin/env python3
"""Registered native-container gate; serial and real MPI have distinct IDs."""
from pathlib import Path
import hashlib
import json
import math
import os
import re
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from lib.runner_support import case, command, error_result, load_context, main_guard, run_command, write_result

def main():
    context = load_context('Native material geometry and unequal-lifetime sampling')
    repo = Path(context['repo_root'])
    work = Path(context['work_dir'])
    artifact = Path(context['artifact_dir'])
    serial = context['test']['id'].endswith('SERIAL-001')
    identity = 'PIC-FIELDLINE-NATIVE-SERIAL-001' if serial else 'PIC-FIELDLINE-NATIVE-MPI1-001'
    scope = 'single-process MPI dependency adapter; NOT MPI qualification' if serial else 'native one-rank MPI; NOT a decomposition/coupling qualification'
    try:
        build = Path(os.environ.get('AMPS_FIELD_LINE_NATIVE_BUILD', str(repo))).resolve()
        native = build/'build/pic'
        required = [build/'build/libAMPS.a', native/'pic.h', native/'pic_field_line_geometry.h']
        if not all(p.is_file() for p in required):
            return write_result(context, 'ERROR', 'Declared native-field-line capability lacks a configured, clean native build', cases=[case(identity, 'ERROR', scope)])
        # Configuration adds macro guards and species tables to generated
        # sources. Compare the actual field-line class layout separately from
        # those unrelated species substitutions; hashing an arbitrary old
        # library alone would not establish that this checkout was tested.
        def normalized(text):
            text = re.sub(r'/\* AMPS_CHECKMACRO_BEGIN \*/.*?/\* AMPS_CHECKMACRO_END \*/', '', text, flags=re.S)
            return re.sub(r'\s+', '', text)
        def field_line_layout(text):
            return text[text.index('namespace FieldLine{'):text.index('namespace ParticleSplitting {')]
        if normalized(field_line_layout((repo/'src/pic/pic.h').read_text())) != normalized(field_line_layout((native/'pic.h').read_text())):
            raise ValueError('Native field-line class layout differs from the checkout; regenerate the build')
        for name in ['pic_field_line.cpp', 'field_line_parallel.cpp', 'pic_field_line_geometry.h']:
            generated = native/name
            if not generated.is_file() or normalized(generated.read_text()) != normalized((repo/'src/pic'/name).read_text()):
                raise ValueError('Native generated source differs from the checkout: '+name)
            required.append(generated)
        # These ordinary objects and their archive must postdate the relevant
        # generated layout/kernel. This detects stale incremental builds; the
        # documented clean-build prerequisite still covers other dependencies.
        for name in ['pic_field_line', 'field_line_parallel']:
            obj = native/(name+'.o')
            dependencies = [native/(name+'.cpp'), native/'pic.h', native/'pic_field_line_geometry.h']
            if not obj.is_file() or obj.stat().st_mtime_ns < max(p.stat().st_mtime_ns for p in dependencies) or required[0].stat().st_mtime_ns < obj.stat().st_mtime_ns:
                raise ValueError('Native field-line object/archive is stale: '+name)
            required.append(obj)
        sources = [repo/'src/pic/pic.h', repo/'src/pic/pic_field_line.cpp', repo/'src/pic/field_line_parallel.cpp', repo/'src/pic/pic_field_line_geometry.h', Path(__file__), Path(__file__).with_name('native_geometry.cpp'), Path(__file__).with_name('serial_mpi_dependency.cpp')]
        before = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources+required}
        executable = work/'native-field-line-fixture'
        environment = {'AMPS_NATIVE_TEST_REPO': str(repo), 'AMPS_NATIVE_TEST_OUTPUT': str(executable)}
        target = 'native-field-line-fixture-serial' if serial else 'native-field-line-fixture'
        result = run_command(command(context, 'make')+['-f', str(Path(__file__).with_name('native_geometry.mk')), target], build, timeout=300, env=environment)
        if result.returncode:
            return write_result(context, 'FAIL', 'Native fixture failed to compile/link; clean matching field-line build is required', cases=[case(identity, 'FAIL', scope)])
        argv = [str(executable)] if serial else command(context, 'mpiexec')+[context['execution']['mpi_rank_flag'], '1', str(executable)]
        result = run_command(argv, work, timeout=120)
        artifact.mkdir(parents=True, exist_ok=True)
        transcript = artifact/'native-output.txt'
        transcript.write_text(result.stdout or '')
        if result.returncode or 'PASS native material insertion/removal/particle-coordinate transaction' not in (result.stdout or ''):
            # MPI startup is infrastructure ERROR, not a physics/reference FAIL.
            status = 'ERROR' if 'NATIVE_CONTAINER_STARTED' not in (result.stdout or '') else 'FAIL'
            return write_result(context, status, 'Native execution did not complete: '+scope, cases=[case(identity, status, scope)], evidence=[str(transcript)])
        data = work/'native-material-sampling.dat'
        lines = data.read_text().splitlines()
        columns = re.findall(r'"([^"]+)"', next(line for line in lines if line.startswith('VARIABLES=')))
        rows = []
        for line in lines:
            try:
                values = [float(value) for value in line.split()]
            except ValueError:
                continue
            if values:
                rows.append(values)
        column = columns.index('Number Density[1/m^3]')
        if len(rows)!=5 or any(len(row)!=len(columns) for row in rows):
            raise ValueError('Native field-line output has malformed dimensions')
        if any(not math.isfinite(row[column]) for row in rows) or rows[0][column]!=2.25 or rows[2][column]!=2.25:
            raise ValueError('Fresh/old segment densities fail independent equal-density oracle')
        exposure = Path(str(data)+'.sampling')
        history = [line.split() for line in exposure.read_text().splitlines() if line and not line.startswith('#')]
        if len(history)!=4 or [int(row[4]) for row in history]!=[1,1,2,2] or [float(row[6]) for row in history]!=[0.75]*4:
            raise ValueError('Native exposure metadata does not preserve actual segment lifetimes/constant dt')
        after = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources+required}
        if before!=after:
            raise ValueError('Source or native build changed during execution; evidence discarded')
        provenance = artifact/'provenance.json'
        provenance.write_text(json.dumps({'scope': scope, 'sha256': before, 'fresh_density': rows[0][column], 'old_density': rows[2][column]}, indent=2)+'\n')
        return write_result(context, 'PASS', scope, cases=[case(identity, 'PASS', 'Actual container, particle loss, shock intersection, fluence and unequal-lifetime native output')], evidence=[str(transcript), str(provenance), str(data), str(exposure)])
    except Exception as error:
        return error_result(context, error)

if __name__ == '__main__':
    main_guard(main)
