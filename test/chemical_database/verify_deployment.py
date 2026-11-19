"""Read-only real-release deployment/relocation checks, no fabricated data.

The caller already opened an authenticated prepared database. Copy that small
published archive into owned test scratch, relocate the copy, run its OWN
installed CLI and compile its OWN C++ SDK. All original archive bytes remain
read-only, and no AMPS source-tree include/library is supplied to compilation.
"""
import json
from pathlib import Path
import shutil


def run_checks(context,database):
    output=context.output/'installed-deployment';output.mkdir(exist_ok=False)
    source=database.root
    # Copies, never hardlinks: an accidental test write cannot modify a
    # canonical or installed source object through shared inode ownership.
    copied=output/'initial-copy';shutil.copytree(source,copied)
    moved=output/'relocated-copy';copied.rename(moved)
    if copied.exists():raise RuntimeError('relocation retained stale original copy path')
    lock=moved/'database-lock.json'
    relocated=context.python_api().Database.from_lock(lock)
    if relocated.table('iaea-collisiondb-research:e33e73d1967727c64e6096db3df8a57f20af4197639e07fb77ee5062d0ab4af1:table').rows != database.table('iaea-collisiondb-research:e33e73d1967727c64e6096db3df8a57f20af4197639e07fb77ee5062d0ab4af1:table').rows:
        raise AssertionError('relocated real native lexical table changed')
    invocation=[context.python_executable,'-B',str(moved/'bin/chemdb'),'--lock',str(lock),'verify']
    cli=context.command(invocation,'installed-cli')
    if cli.returncode:raise RuntimeError('actual relocated installed CLI failed: '+cli.stderr)
    verified=json.loads(cli.stdout)
    if verified.get('status')!='PASS' or verified.get('physical_activation') is not False:
        raise AssertionError('installed integrity CLI scope/status changed')
    module=moved/'sdk/src/models/chemical_database'
    names=('Runtime','DataIO','Photochemistry','ElectronImpact','ChargeExchange','RateServices')
    sources=[module/(name+'.cpp') for name in names]
    executable=output/'sdk-inspect'
    command=context.cxx+['-std=c++17','-O2','-g','-Wall','-Wextra','-Werror','-pedantic',
        *(str(path) for path in sources),str(module/'example.cpp'),'-o',str(executable)]
    compiled=context.command(command,'installed-sdk-compile',timeout=180)
    if compiled.returncode:raise RuntimeError('actual relocated SDK did not compile: '+compiled.stderr)
    pin=relocated.identity
    native=context.command([str(executable),str(moved),pin['export_path'],pin['export_sha256'],
        pin['release_id'],pin['release_manifest_sha256'],pin['profile_id'],pin['profile_sha256'],'catalogue'],
        'installed-sdk-query')
    if native.returncode:raise RuntimeError('actual relocated native SDK rejected published release: '+native.stdout+native.stderr)
    import hashlib
    def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
    receipt={'scope':'real lexical source/integrity and deployed software relocation; NOT physical validation',
        'original_database_root':str(source),'relocated_root':str(moved),'lock_sha256':digest(lock),
        'export_sha256':pin['export_sha256'],'native_executable_sha256':digest(executable),
        'installed_sdk_sources':{str(path.relative_to(moved)):digest(path) for path in
             [*sources,*(module/(name+'.h') for name in names),module/'example.cpp',moved/'sdk/src/general/input/amps_input.h']},
        'installed_python_sources':{str(path.relative_to(moved)):digest(path) for path in
             sorted((moved/'python/amps_chemical_database').glob('*.py'))},
        'cli_integrity':verified,'source_tree_compiler_includes_or_libraries':[],
        'status':'PASS','physical_activation':False}
    (output/'result.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    return receipt
