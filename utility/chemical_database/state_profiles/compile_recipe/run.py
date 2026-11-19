#!/usr/bin/env python3
"""Real-source compile-table checks, without Config or AMPS application build."""
from __future__ import annotations
import argparse,importlib.util,json,subprocess,time,traceback,re,shutil
from pathlib import Path
import prepare
ROOT=Path(__file__).resolve().parents[4]

def parser_species_symbol_check(recipe,output,manifest):
    """Run the exact maintained ReadSpeciesBlock body, then compile CEX names.

    Only the normal parser's startup is replaced by explicit disposable file
    handles. The production subroutine and ampsConfigLib perform every species
    interpretation/macro write; there is no reimplementation of that logic.
    """
    work=output/'actual-species-parser';pic=work/'generated/pic';pic.mkdir(parents=True)
    source=(ROOT/'ampsConfig.pl').read_text()
    body=source[source.index('sub ReadSpeciesBlock {'):source.index('#=============================== Process AMPS\' settings')]
    for name in ('pic.h','pic_mol.cpp','picSpeciesMacro.dfn'):
        shutil.copyfile(ROOT/'src/pic'/name,pic/name)
    assembled=work/'species.Assembled'
    lines=(recipe/'source_species.input').read_text().splitlines()
    assembled.write_text(''.join(str(i+1)+' source_species.input\n'+line+'\n' for i,line in enumerate(lines)))
    harness=work/'species_parser.pl'
    harness.write_text('use strict;use warnings;use lib "'+str(ROOT)+'";use ampsConfigLib;\n'
        'my ($line,$InputFileLineNumber,$FileName,$InputLine,$InputComment,$s0,$s1);\n'
        'my $InputFileName="'+str(work/'normal')+'";my @SpeciesList=("NA","ELECTRON");my $TotalSpeciesNumber=2;\n'
        '$ampsConfigLib::WorkingSourceDirectory="'+str(work/'generated')+'";\n'
        'open(InputFile,"<'+str(assembled)+'") or die $!;ReadSpeciesBlock();close(InputFile);\n'+body)
    command=['perl',str(harness)];child=subprocess.run(command,cwd=work,capture_output=True,text=True,timeout=15)
    (work/'parser.stdout.log').write_text(child.stdout);(work/'parser.stderr.log').write_text(child.stderr)
    manifest['commands'].append({'argv':command,'returncode':child.returncode})
    if child.returncode:raise RuntimeError('normal ReadSpeciesBlock failed; see actual-species-parser logs')
    cex_sources=[ROOT/'src/models/charge_exchange'/name for name in
                 ('ChargeExchange_MaherTinsley.cpp','ChargeExchange_LindsayStebbings.cpp')]
    cex='\n'.join(p.read_text() for p in cex_sources)
    symbols=sorted(set(re.findall(r'_H_ENA_V[123]_SPEC_',cex)))
    if len(symbols)!=3:raise ValueError('expected explicit three canonical CEX ENA symbol dependencies')
    probe=work/'symbols.cpp';probe.write_text('#include "'+str(pic/'picSpeciesMacro.dfn')+'"\n'+
        '\n'.join('static_assert('+s+'<0,"unselected CEX dependency remains inactive");' for s in symbols)+'\nint main(){return 0;}\n')
    command=['g++','-std=c++17','-Wall','-Wextra','-Werror','-pedantic',str(probe),'-o',str(work/'symbols')]
    child=subprocess.run(command,capture_output=True,text=True,timeout=15)
    (work/'compiler.stdout.log').write_text(child.stdout);(work/'compiler.stderr.log').write_text(child.stderr)
    manifest['commands'].append({'argv':command,'returncode':child.returncode})
    manifest['species_parser_source_sha256']={str(p):prepare.sha(p) for p in
        (ROOT/'ampsConfig.pl',ROOT/'ampsConfigLib.pm',ROOT/'input/species.input',*cex_sources)}
    if child.returncode:raise RuntimeError('generated unselected CEX symbols do not compile')
    return True

def parser_compiler_reference_check(recipe,output,manifest):
    """Exercise the maintained normal Exosphere DEFINE parser in isolation.

    Copied model files are disposable parser outputs only, never a primary
    production fix. No Config.pl, AMPS application build or force law executes.
    The generated macro is independently stringified by the C++ preprocessor.
    """
    work=output/'actual-exosphere-parser';model=work/'generated/models/exosphere';model.mkdir(parents=True)
    files=[ROOT/'exosphere.pl',ROOT/'ampsConfigLib.pm']+[ROOT/'src/models/exosphere'/('Exosphere.'+suffix) for suffix in ('cpp','h','dfn')]+[ROOT/'src/models/exosphere/ExosphereDefinitions.dfn',ROOT/'srcMoon/test/sodium_cohort/zero_force_recipe.h',ROOT/'src/meshAMR/meshAMR_UserDefinitions.h']
    manifest['parser_source_sha256']={str(p):prepare.sha(p) for p in files}
    shutil.copyfile(ROOT/'ampsConfigLib.pm',work/'ampsConfigLib.pm')
    for name in ('Exosphere.cpp','Exosphere.h','Exosphere.dfn','ExosphereDefinitions.dfn'):shutil.copyfile(ROOT/'src/models/exosphere'/name,model/name)
    mesh=work/'generated/meshAMR';mesh.mkdir()
    shutil.copyfile(ROOT/'src/meshAMR/meshAMR_UserDefinitions.h',mesh/'meshAMR_UserDefinitions.h')
    text=(recipe/'moon_source_profile.input').read_text()
    definitions=[line for line in text.splitlines() if re.match(r'\s*define\s+(_PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_|_MOON_COHORT_ZERO_FORCE_REFERENCE_)',line)]
    if len(definitions)!=2:raise ValueError('exact two reference DEFINE directives required')
    block=work/'reference.Block';block.write_text(''.join(str(i+1)+' '+str(recipe/'moon_source_profile.input')+'\n'+line+'\n' for i,line in enumerate(definitions))+'3 native-reference\n#endblock\n')
    (work/'reference.Species').write_text('2\nNa\nELECTRON\n')
    command=['perl',str(ROOT/'exosphere.pl'),str(block),str(work/'generated')]
    child=subprocess.run(command,cwd=work,capture_output=True,text=True,timeout=15)
    (work/'parser.stdout.log').write_text(child.stdout);(work/'parser.stderr.log').write_text(child.stderr)
    manifest['commands'].append({'argv':command,'cwd':str(work),'returncode':child.returncode})
    if child.returncode:raise RuntimeError('normal Exosphere parser failed; see parser logs')
    probe=work/'macro_check.cpp';probe.write_text('#include <string>\n#include "'+str(model/'Exosphere.dfn')+'"\n#include "'+str(ROOT/'srcMoon/test/sodium_cohort/zero_force_recipe.h')+'"\n#define SHOW_RAW(x) #x\n#define SHOW(x) SHOW_RAW(x)\nstatic_assert(_MOON_COHORT_ZERO_FORCE_REFERENCE_==1,"explicit reference marker");\nint main(){const std::string actual=SHOW(_PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_(a,s,p,x,v,n));return actual=="PIC::Mover::TotalParticleAcceleration_default(a,s,p,x,v,n)" && Moon::SodiumCohortInput::ExactZeroAccelerationMacro(actual)?0:1;}\n')
    executable=work/'macro_check';command=['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror','-pedantic',str(probe),'-o',str(executable)]
    child=subprocess.run(command,cwd=work,capture_output=True,text=True,timeout=15)
    (work/'compiler.stdout.log').write_text(child.stdout);(work/'compiler.stderr.log').write_text(child.stderr)
    manifest['commands'].append({'argv':command,'cwd':str(work),'returncode':child.returncode})
    if child.returncode:raise RuntimeError('parsed reference macro compilation failed')
    child=subprocess.run([str(executable)],cwd=work,capture_output=True,text=True,timeout=5)
    manifest['commands'].append({'argv':[str(executable)],'cwd':str(work),'returncode':child.returncode})
    if child.returncode:raise RuntimeError('parsed reference binding differs from declared production zero-function call')
    return True

def run(args):
    output=args.output_dir.resolve();output.mkdir(parents=True,exist_ok=False);started=time.monotonic();checks=[]
    manifest={'scope':'normal compile-input/table recipe ONLY; no physical provider activated','root_configuration_modified':False,'commands':[]}
    def check(name,ok):checks.append({'name':name,'status':'PASS' if ok else 'FAIL'})
    try:
        spec=importlib.util.spec_from_file_location('m3_source_profile_prepare',Path(__file__).resolve().parents[1]/'prepare.py');source_prepare=importlib.util.module_from_spec(spec);spec.loader.exec_module(source_prepare)
        profile=source_prepare.prepare(args.raw_dir,args.decision_register,output/'source-profile',args.canonical_root)
        profile_path=output/'source-profile/profile.json'
        protected=[ROOT/p for p in ('src/general/constants.h','input/species.input','input/moon.input','Makefile','Makefile.local') if (ROOT/p).is_file()]
        before={str(p):(prepare.sha(p),p.stat().st_mtime_ns) for p in protected}
        mapping=output/'mapping.json';mapping.write_text(json.dumps({'Na':'23Na:ground','ELECTRON':'electron:free'}))
        receipt=prepare.prepare(ROOT,profile_path,mapping,ROOT/'input/moon.input',output/'available')
        species=(output/'available/source_species.input').read_text()
        normal=(output/'available/moon_source_profile.input').read_text();original=(ROOT/'input/moon.input').read_text()
        force=lambda text:[line for line in text.splitlines() if re.match(r'\s*define\s+_PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_',line)]
        check('default normal-template preserves exact existing force definition',force(normal)==force(original) and '_MOON_COHORT_ZERO_FORCE_REFERENCE_' not in normal and receipt['reference_problem']=='normal-template')
        check('default packet layout preserves exact template bytes',prepare.packet_template(original)==original and receipt['packet_layout']=='template')
        conservative=prepare.packet_template(original,'conservative-products')
        weight=lambda text:re.findall(r'(?mi)^\s*ParticleWeightCorrectionMode\s*=\s*(on|off)',text)
        check('explicit conservative product packet enables normal individual weights',weight(conservative)==['on'] and prepare.packet_template(conservative,'conservative-products')==conservative)
        timestep=lambda text:[line for line in text.splitlines() if re.match(r'\s*TimeStep',line,re.I)]
        check('conservative packet leaves timestep policy unchanged',timestep(conservative)==timestep(original))
        for name,text in [('missing',re.sub(r'(?mi)^.*ParticleWeightCorrectionMode.*\n','',original)),('duplicate',original+'\nParticleWeightCorrectionMode=on\n')]:
            try:prepare.packet_template(text,'conservative-products');ok=False
            except ValueError:ok=True
            check('ambiguous '+name+' packet layout refuses',ok)
        zero=prepare.prepare(ROOT,profile_path,mapping,ROOT/'input/moon.input',output/'homogeneous',reference_problem='homogeneous-zero-acceleration')
        zero_text=(output/'homogeneous/moon_source_profile.input').read_text()
        check('named homogeneous experiment explicitly selects existing core zero function and marker',zero['reference_problem']=='homogeneous-zero-acceleration' and 'PIC::Mover::TotalParticleAcceleration_default(a,s,p,x,v,n)' in zero_text and 'define _MOON_COHORT_ZERO_FORCE_REFERENCE_ 1' in zero_text)
        check('normal Exosphere parser and compiler preserve actual zero reference binding',parser_compiler_reference_check(output/'homogeneous',output,manifest))
        for name,changed in [('missing',re.sub(r'(?m)^.*define[ \t]+_PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_.*\n','',original)),('duplicate',original+'\n'+force(original)[0]+'\n'),('marker',original+'\n define _MOON_COHORT_ZERO_FORCE_REFERENCE_ 1\n')]:
            try:prepare.reference_template(changed,'homogeneous-zero-acceleration');ok=False
            except ValueError:ok=True
            check('ambiguous '+name+' reference template refuses',ok)
        check('actual available source states stage explicit compile entries',len(receipt['native_receipts'])==3 and '#component=Na' in species and '#component=ELECTRON' in species)
        original_species=(ROOT/'input/species.input').read_text()
        markers=lambda text:re.findall(r'(?mi)^\s*#component\s*=\s*(\S+)',text)
        check('complete canonical component inventory and inactive ENA names preserved',markers(species)==markers(original_species) and all('#component='+name in species for name in ('H_ENA_V1','H_ENA_V2','H_ENA_V3')))
        unchanged=lambda text:re.sub(r'(?mi)^\s*(mass|charge)\s*=.*$','',text)
        check('canonical non-mass-charge fields and component order unchanged',unchanged(species)==unchanged(original_species))
        check('actual ReadSpeciesBlock emits compilable inactive CEX dependency symbols',parser_species_symbol_check(output/'available',output,manifest))
        check('source SI unit exactly matches authenticated CODATA value',receipt['compiled_charge_unit_C']=='1.6021766339999999e-19')
        check('native emitted nonzero charge expression independently compiled',all(r['returncode']==0 for r in receipt['commands']) and any(row[0]=='STATE' and row[1]=='electron:free' and row[4]=='-1' and row[5]!='-1' for row in receipt['native_receipts']))
        check('unrelated background species selection preserved','SpeciesList=CO2,O' in (output/'available/moon_source_profile.input').read_text())
        check('application parser and kernel references use selected basename','#block ../moon_source_profile.input.pl' in (output/'available/moon_source_profile.input').read_text() and '#include moon_source_profile.input.spicekernels' in (output/'available/moon_source_profile.input').read_text())
        check('normal kernel companion preserved',prepare.sha(output/'available/moon_source_profile.input.spicekernels')==prepare.sha(ROOT/'input/moon.input.spicekernels'))
        check('normal application compiler companion preserved',prepare.sha(output/'available/moon_source_profile.input.pl')==prepare.sha(ROOT/'input/moon.input.pl'))
        missing=output/'missing.json';missing.write_text(json.dumps({'Na':'23Na:ground','NA_PLUS':'23Na+:unresolved'}))
        try:prepare.prepare(ROOT,profile_path,missing,ROOT/'input/moon.input',output/'missing-state');ok=False
        except LookupError as error:ok=not (output/'missing-state').exists();(output/'missing-state-refusal.log').write_text(str(error)+'\n')
        check('actual missing ion state refuses before output',ok)
        try:prepare.prepare(ROOT,profile_path,missing,ROOT/'input/moon.input',output/'missing-homogeneous',reference_problem='homogeneous-zero-acceleration');ok=False
        except LookupError as error:ok=not (output/'missing-homogeneous').exists();(output/'missing-homogeneous-refusal.log').write_text(str(error)+'\n')
        check('named reference does not bypass missing ion source refusal',ok)
        corrupt=json.loads(profile_path.read_text());corrupt['sources'][0]['sha256']='0'*64;changed=output/'corrupt.json';changed.write_text(json.dumps(corrupt))
        try:prepare.prepare(ROOT,changed,mapping,ROOT/'input/moon.input',output/'corrupt');ok=False
        except RuntimeError as error:ok=not (output/'corrupt').exists();(output/'corrupt-refusal.log').write_text(str(error)+'\n')
        check('changed source digest refuses before output',ok)
        try:prepare.prepare(ROOT,profile_path,mapping,ROOT/'input/moon.input',output/'available');ok=False
        except FileExistsError:ok=True
        check('existing output evidence preserved',ok)
        try:prepare.prepare(ROOT,profile_path,mapping,ROOT/'input/moon.input',output/'space unsafe');ok=False
        except ValueError:ok=not (output/'space unsafe').exists()
        check('unsupported native include path refuses before output',ok)
        check('root source/configuration bytes and mtimes unchanged',before=={str(p):(prepare.sha(p),p.stat().st_mtime_ns) for p in protected})
        manifest['recipe']=receipt
    except FileNotFoundError as error:checks.append({'name':'required-source-or-compiler','status':'SKIPPED','reason':str(error)})
    except Exception as error:
        checks.append({'name':'recipe-infrastructure','status':'ERROR','reason':str(error)});(output/'exception.log').write_text(traceback.format_exc())
    files=[Path(__file__),Path(__file__).with_name('prepare.py'),Path(__file__).with_name('export.cpp'),Path(__file__).with_name('CONTRACT.md'),Path(__file__).with_name('README.md')]
    manifest['source_sha256']={str(p):prepare.sha(p) for p in files}
    counts={s:sum(c['status']==s for c in checks) for s in ('PASS','FAIL','ERROR','SKIPPED')}
    status=next(s for s in ('ERROR','FAIL','SKIPPED','PASS') if counts[s])
    result={'schema':'srcMoon-M3-compile-species-recipe-controls-v1','status':status,'counts':counts,'checks':checks,'scope':manifest['scope'],'formal_M3_status':'SKIPPED / NOT QUALIFIED','elapsed_wall_seconds':time.monotonic()-started}
    (output/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');(output/'result.json').write_text(json.dumps(result,indent=2)+'\n');return result
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',required=True,type=Path)
    packet=ROOT/'test_output/srcMoon/m3-approved-products-20261118T163142Z/source-acquisition-v2'
    p.add_argument('--raw-dir',type=Path,default=packet/'raw');p.add_argument('--decision-register',type=Path,default=ROOT/'srcMoon/test/M3_DAUGHTER_AND_ION_SURFACE_DECISIONS.json');p.add_argument('--canonical-root',type=Path,default=Path('/data/vtenishe/moon_validation_data/photochemistry'));a=p.parse_args()
    if a.output_dir.exists():p.error('output must not exist')
    r=run(a);print(json.dumps(r,indent=2));raise SystemExit({'PASS':0,'FAIL':1,'ERROR':2,'SKIPPED':77}[r['status']])
