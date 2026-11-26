#!/usr/bin/env python3
"""Stage an application-owned species include through the real native resolver.

No AMPS configuration, global constant edits or profile qualification occurs.
The caller's output does not exist until all native identity/mapping checks and
an independent compilation of the emitted charge expressions have succeeded.
"""
from pathlib import Path
import argparse,hashlib,json,re,subprocess,tempfile,shutil

def strict(path):
    def unique(items):
        d={}
        for k,v in items:
            if k in d:raise ValueError('duplicate JSON key '+k)
            d[k]=v
        return d
    def bad(x):raise ValueError('nonfinite JSON '+x)
    return json.loads(Path(path).read_text(),object_pairs_hook=unique,parse_constant=bad)

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def packet_template(text,packet_layout='template'):
    """Opt in to normal individual-weight storage for conservative products.

    This is a particle-record ABI choice, not an altered physical particle
    population or timestep policy. The caller must clean-build this input.
    The default returns the caller's exact original bytes.
    """
    if packet_layout=='template':return text
    if packet_layout!='conservative-products':raise ValueError('unknown explicit packet layout')
    pattern=re.compile(r'(?mi)^([ \t]*ParticleWeightCorrectionMode[ \t]*=[ \t]*)(on|off)([ \t]*(?:![^\r\n]*)?)$')
    matches=list(pattern.finditer(text))
    if len(matches)!=1:raise ValueError('conservative product layout requires one normal ParticleWeightCorrectionMode assignment')
    return pattern.sub(lambda m:m.group(1)+'on'+m.group(3),text)

def source_species_include(canonical,mapping,native):
    """Replace ONLY selected mass/charge fields in the complete species owner.

    ReadSpeciesBlock emits negative-index macros for every unselected component.
    Those declarations are required even when a legacy shared source module is
    not active in this particular simulation. Dropping unselected components
    therefore changes the normal compile API; it must not be used to stage an
    application-owned physical mass table. All other canonical bytes survive.
    """
    selected={symbol.upper():(symbol,state) for symbol,state in mapping.items()}
    markers=list(re.finditer(r'(?mi)^[ \t]*#component[ \t]*=[ \t]*([A-Za-z0-9_]+)[^\r\n]*$',canonical))
    found=set();pieces=[];end=0
    for i,match in enumerate(markers):
        name=match.group(1).upper();stop=markers[i+1].start() if i+1<len(markers) else len(canonical)
        block=canonical[match.start():stop]
        if name in selected:
            if name in found:raise ValueError('duplicate canonical selected component '+name)
            found.add(name);symbol,state=selected[name];row=native[state]
            for key,value in (('mass',row[2]),('charge',row[5])):
                pattern=re.compile(r'(?mi)^([ \t]*'+key+r'[ \t]*=[ \t]*)[^!\r\n]*')
                block,count=pattern.subn(lambda m:m.group(1)+value,block)
                if count!=1:raise ValueError('canonical selected component requires exactly one '+key+': '+symbol)
        pieces.extend((canonical[end:match.start()],block));end=stop
    pieces.append(canonical[end:])
    if found!=set(selected):raise ValueError('mapped component absent from canonical species table: '+','.join(sorted(set(selected)-found)))
    return ''.join(pieces)

# Named reference experiments: each substitutes ONE existing production
# function for the acceleration macro and adds its marker. No force law is
# written here. lunar-point-mass-gravity (I21 decision 1, approved 2026-11-26)
# applies only the production lunar point-mass term.
REFERENCE_PROBLEMS={
    'homogeneous-zero-acceleration':('PIC::Mover::TotalParticleAcceleration_default(a,s,p,x,v,n)',
                                     '_MOON_COHORT_ZERO_FORCE_REFERENCE_'),
    'lunar-point-mass-gravity':('Moon::LunarPointMassReferenceAcceleration(a,s,p,x,v,n)',
                                '_MOON_COHORT_LUNAR_GRAVITY_REFERENCE_'),
}


def reference_template(text,reference_problem='normal-template'):
    """Select an explicit reference experiment in NORMAL AMPS input syntax.

    No source-state availability or physical approval is inferred here. The
    normal default is exactly the caller's template. An opt-in experiment
    names an existing production function; it never substitutes a test force law.
    """
    if reference_problem!='normal-template' and reference_problem not in REFERENCE_PROBLEMS:
        raise ValueError('unknown explicit reference problem')
    if reference_problem=='normal-template':return text
    function,marker=REFERENCE_PROBLEMS[reference_problem]
    for other in REFERENCE_PROBLEMS.values():
        if re.search(r'(?mi)^\s*define\s+'+other[1]+r'\b',text):
            raise ValueError('reference template already contains a cohort marker')
    pattern=re.compile(r'(?m)^([ \t]*define[ \t]+)_PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_\(([^()\r\n]*)\)[ \t]+([^\r\n]+)$')
    matches=list(pattern.finditer(text))
    if len(matches)!=1:raise ValueError('reference experiment requires one simple acceleration definition')
    match=matches[0];arguments=match.group(2).split(',')
    if len(arguments)!=6 or any(not re.fullmatch(r'[A-Za-z_][A-Za-z0-9_]*',x.strip()) for x in arguments):
        raise ValueError('reference acceleration definition requires six simple arguments')
    blocks=list(re.finditer(r'(?m)^\s*#block[ \t]+\.\./exosphere\.pl(?:[ \t]*![^\r\n]*)?[ \t]*$',text))
    if len(blocks)!=1:raise ValueError('normal Exosphere compile-input block must be explicit and unique')
    end=re.search(r'(?mi)^\s*#endblock\b',text[blocks[0].end():])
    if end is None or not blocks[0].end()<match.start()<blocks[0].end()+end.start():
        raise ValueError('acceleration definition is outside normal Exosphere input block')
    replacement=(match.group(1)+'_PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_(a,s,p,x,v,n) '
        +function+'\n'+match.group(1)+marker+' 1')
    return text[:match.start()]+replacement+text[match.end():]


def prepare(amps_root,profile_path,mapping_path,template_path,output,*,reference_problem='normal-template',packet_layout='template'):
    amps_root,profile_path,mapping_path,template_path,output=map(lambda p:Path(p).resolve(),(amps_root,profile_path,mapping_path,template_path,output))
    if output.exists():raise FileExistsError('output already exists')
    if any(c.isspace() or c in '!#' for c in str(output)):raise ValueError('compile-time include path cannot contain whitespace or comment/directive characters')
    profile,mapping=strict(profile_path),strict(mapping_path)
    if profile.get('schema')!='srcMoon-M3-state-profile-v1' or not isinstance(mapping,dict) or not mapping:raise ValueError('explicit profile schema/nonempty mapping required')
    if any(not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*',k) or not isinstance(v,str) or not v for k,v in mapping.items()):raise ValueError('unsafe compiled symbol or empty explicit state')
    if len({k.upper() for k in mapping})!=len(mapping) or len(set(mapping.values()))!=len(mapping):raise ValueError('duplicate compiled symbol or mapped state')
    symbols=list(mapping)
    original=packet_template(reference_template(template_path.read_text(),reference_problem),packet_layout)
    marker=re.search(r'(?mi)^\s*#endmain\s*$',original)
    if not marker:raise ValueError('template main section terminator missing')
    main_section=original[:marker.end()]
    main=re.compile(r'(?mi)^(\s*SpeciesList\s*=)[^\n]*$')
    include=re.compile(r'(?mi)^\s*#include\s+species\.input\s*$')
    if len(main.findall(main_section))!=1 or len(include.findall(original))!=1:raise ValueError('template requires exactly one SpeciesList and species.input include')
    compiler=shutil.which('g++')
    if not compiler:raise FileNotFoundError('g++ is required')
    sources=[Path(__file__).with_name('export.cpp'),amps_root/'srcMoon/species/PhysicalStateProfile.cpp',amps_root/'srcMoon/species/ReducedKinematics.cpp',amps_root/'src/models/chemical_database/DataIO.cpp']
    headers=[amps_root/p for p in ('srcMoon/species/PhysicalStateProfile.h','srcMoon/species/ReducedKinematics.h','srcMoon/species/ReactionExecution.h','srcMoon/species/Photochemistry.h','src/models/chemical_database/DataIO.h','src/general/constants.h','src/general/constants.AtomicData.h','src/general/constants.PlanetaryData.h')]
    commands=[]
    with tempfile.TemporaryDirectory(prefix='amps-profile-export-') as work:
        work=Path(work);exe=work/'export'
        argv=[compiler,'-std=c++17','-O2','-Wall','-Wextra','-Werror','-pedantic','-I',str(amps_root),*[str(s) for s in sources],'-o',str(exe)]
        r=subprocess.run(argv,capture_output=True,text=True,timeout=60);commands.append({'argv':argv,'returncode':r.returncode,'stdout':r.stdout,'stderr':r.stderr})
        if r.returncode:raise RuntimeError('native profile exporter compile failed: '+r.stderr)
        argv=[str(exe)]
        for name in ('isotopes','constants','states'):argv.extend(profile[name][k] for k in ('path','sha256'))
        for record in profile['sources']:argv.extend(record[k] for k in ('id','path','sha256','role'))
        argv.extend(['--states',*mapping.values()]);r=subprocess.run(argv,capture_output=True,text=True,timeout=60);commands.append({'argv':argv,'returncode':r.returncode,'stdout':r.stdout,'stderr':r.stderr})
        if r.returncode==77:raise LookupError('required source state/rest mass unavailable: '+r.stderr.strip())
        if r.returncode:raise RuntimeError('native profile refused: '+r.stderr.strip())
        rows=[x.split('\t') for x in r.stdout.splitlines()]
        constant=rows[0]
        native={row[1]:row for row in rows[1:]}
        if constant[0]!='CONSTANT' or len(native)!=len(mapping) or set(native)!=set(mapping.values()):raise ValueError('malformed native exporter receipt')
        # Compile the emitted charge expressions separately. This exercises the
        # exact ordinary AMPS expression shape, including native double literal
        # rounding; a Python division is never accepted as an oracle for C++.
        checks=['#include "src/general/constants.h"','int main(){']
        for state in mapping.values():
            row=native[state]
            checks.append('if(('+row[5]+'*ElectronCharge)!=('+row[4]+'.0*'+constant[1]+'))return 1;')
        checks.append('return 0;}');source=work/'verify.cpp';source.write_text('\n'.join(checks)+'\n')
        argv=[compiler,'-std=c++17','-O2','-Wall','-Wextra','-Werror','-pedantic','-I',str(amps_root),str(source),'-o',str(work/'verify')]
        r=subprocess.run(argv,capture_output=True,text=True,timeout=60);commands.append({'argv':argv,'returncode':r.returncode,'stdout':r.stdout,'stderr':r.stderr})
        if r.returncode:raise RuntimeError('charge expression compilation failed: '+r.stderr)
        r=subprocess.run([str(work/'verify')],capture_output=True,text=True,timeout=10)
        if r.returncode:raise ValueError('emitted native charge expression differs from source SI charge')
        canonical_species=amps_root/'input/species.input'
        species=source_species_include(canonical_species.read_text(),mapping,native)
        # An absolute include is intentional: the native assembler resolves
        # includes from process cwd, not from the including file's directory.
        staged=main.sub(lambda m:m.group(1)+','.join(symbols),main_section)+original[marker.end():]
        staged=include.sub('#include '+str(output/'source_species.input'),staged)
        # Match Config.pl's application-basename staging. Depending on old
        # moon.input.pl/kernel files left by an earlier configure is unsafe.
        staged,n=re.subn(r'(?m)^#block[ \t]+\.\./moon\.input\.pl[^\r\n]*$', '#block ../moon_source_profile.input.pl ! source-bound application companion',staged)
        if n!=1:raise ValueError('normal Moon parser block missing or duplicate')
        staged,n=re.subn(r'(?mi)^\s*#include\s+moon\.input\.spicekernels\s*$', '#include moon_source_profile.input.spicekernels',staged)
        if n!=1:raise ValueError('normal Moon kernel include missing or duplicate')
        manifest={'schema':'srcMoon-M3-source-species-recipe-v1','reference_problem':reference_problem,'physical_activation':False,'formal_M3_status':'SKIPPED / NOT QUALIFIED','profile':{'path':str(profile_path),'sha256':sha(profile_path)},'mapping':{'path':str(mapping_path),'sha256':sha(mapping_path),'states':mapping},'template':{'path':str(template_path),'sha256':sha(template_path)},'sources':{str(s):sha(s) for s in sources+headers},'constants_header':{'path':str(amps_root/'src/general/constants.h'),'sha256':sha(amps_root/'src/general/constants.h')},'native_receipts':rows,'commands':commands,'compiled_charge_unit_C':constant[1],'legacy_charge_unit_C':constant[2],'charge_source_identity':constant[3],'charge_source_record':constant[4],'availability':'software recipe only; no complete physical family admitted'}
        manifest['canonical_species']={'path':str(canonical_species),'sha256':sha(canonical_species),'operation':'preserve all components and fields; replace only explicitly mapped mass/charge values'}
        manifest['packet_layout']=packet_layout
        companion=Path(str(template_path)+'.pl')
        if not companion.is_file():raise ValueError('normal application compiler companion is missing')
        manifest['application_compiler']={'path':str(companion),'sha256':sha(companion)}
        kernels=Path(str(template_path)+'.spicekernels')
        if not kernels.is_file():raise ValueError('normal application kernel companion is missing')
        manifest['application_kernels']={'path':str(kernels),'sha256':sha(kernels)}
        output.mkdir(parents=True,exist_ok=False)
        (output/'source_species.input').write_text(species)
        (output/'moon_source_profile.input').write_text(staged)
        shutil.copyfile(companion,output/'moon_source_profile.input.pl')
        (output/'moon_source_profile.input.pl').chmod(companion.stat().st_mode)
        shutil.copyfile(kernels,output/'moon_source_profile.input.spicekernels')
        (output/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
        return manifest
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('amps-root','profile','mapping','template','output-dir'):p.add_argument('--'+key,required=True,type=Path)
    p.add_argument('--reference-problem',choices=('normal-template',*REFERENCE_PROBLEMS),default='normal-template',help='explicit reference experiment; normal-template preserves the supplied force binding')
    p.add_argument('--packet-layout',choices=('template','conservative-products'),default='template',help='explicit particle-record ABI; conservative-products enables normal individual weight correction storage')
    a=p.parse_args()
    try:prepare(a.amps_root,a.profile,a.mapping,a.template,a.output_dir,reference_problem=a.reference_problem,packet_layout=a.packet_layout)
    except LookupError as error:
        print(json.dumps({'status':'SKIPPED','reason':str(error),'formal_M3_status':'SKIPPED / NOT QUALIFIED'}));raise SystemExit(77)
    except (OSError,ValueError,KeyError,RuntimeError,subprocess.TimeoutExpired) as error:
        print(json.dumps({'status':'ERROR','reason':str(error)}));raise SystemExit(2)
