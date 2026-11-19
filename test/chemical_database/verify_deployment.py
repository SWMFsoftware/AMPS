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
    model_deployment = {"selected": False}
    if (module/'ReducedModels.h').is_file():
        # Compile ONLY the relocated header/example, with no checkout include.
        # The explicit unit-coefficient fixture is an analytical software
        # control, never a species/material parameter set or observation.
        example=moved/'sdk/examples/reduced_models.cpp'
        helper=output/'reduced-models-example'
        command=context.cxx+['-std=c++17','-O2','-g','-Wall','-Wextra','-Werror','-pedantic',
            '-I',str(module),str(example),'-o',str(helper)]
        compiled=context.command(command,'installed-reduced-models-compile',timeout=180)
        if compiled.returncode:raise RuntimeError('relocated reduced-model SDK failed to compile: '+compiled.stderr)
        evaluated=context.command([str(helper),'1','0','1','1','1','1'],'installed-reduced-models-query')
        if evaluated.returncode:raise RuntimeError('relocated reduced-model example rejected analytical inputs')
        from decimal import Decimal, localcontext
        values=json.loads(evaluated.stdout,parse_float=Decimal,parse_int=Decimal)
        with localcontext() as arithmetic:
            arithmetic.prec=60
            # A platform-independent roundoff bound for this order-one kernel
            # composition; not a physical-model or observational tolerance.
            roundoff=Decimal('1e-15')
            if values['rate_s_inverse']!=1 or abs(values['remaining_population']-Decimal(-1).exp())>roundoff or abs(values['released_population']-(1-Decimal(-1).exp()))>roundoff:
                raise AssertionError('relocated thermal helper violates independent analytic inventory solution')
        refused=context.command([str(helper),'1','0','1','0','1','1'],'installed-reduced-models-refusal')
        if refused.returncode!=2:raise AssertionError('relocated helper did not refuse zero kinetic temperature')
        listing=context.command([context.python_executable,'-B',str(moved/'bin/chemdb'),'--lock',str(lock),
            'list','--kind','models'],'installed-model-metadata')
        if listing.returncode:raise RuntimeError('relocated installed model CLI failed: '+listing.stderr)
        if len(json.loads(listing.stdout))!=len(relocated.models()):
            raise AssertionError('relocated installed CLI model selection differs from pinned metadata')
        model_deployment={"selected":True,"models":len(relocated.models()),
            "header_and_example_used_only_relocated_SDK":True,"analytical_inventory_and_refusal":"PASS"}
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
        'reduced_model_deployment':model_deployment,
        'status':'PASS','physical_activation':False}
    (output/'result.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    return receipt
