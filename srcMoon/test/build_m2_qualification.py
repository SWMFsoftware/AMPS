#!/usr/bin/env python3
"""Build the dedicated Na+Ar srcMoon M2 qualification executable.

The canonical ``input/moon.input`` is the frozen Na-only M0/M1 application and
must not be edited merely to run an M2 species test.  This script derives a
temporary Config.pl application with reviewed compile-time changes: add
neutral Ar to ``SpeciesList``, register the dedicated ``RadiogenicArgon``
source, disable photolytic chemistry for the I19 accounting contract, and
enable individual particle weights for fractional finite-release packets.
It verifies that each applicable edit occurs once,
performs a clean SPICE build, copies the resulting executable into the
campaign directory, and removes only the three temporary input files it made.
It also removes the uniquely prefixed root copies and `.Assembled*`
intermediates that Config.pl derives from those inputs.

Runtime activation remains separate: I18--I20 derive a preserved ``.in``
receipt with ``argon_source_enabled=true``.  This two-key arrangement makes it
impossible for an accidental Na-only build or an accidental runtime switch to
silently masquerade as a two-species qualification run.

The optional purpose ``sampling`` applies only the first two species/source
changes, retaining canonical chemistry and individual-weight controls,
for I11's deterministic storage regression. Its executable is named separately
and the build receipt explicitly omits the I19 chemistry contract. The default
M2 purpose still requires the original strict I19 generated-configuration gate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time


REPO = Path(__file__).resolve().parents[2]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def replace_once(text: str, pattern: str, replacement: str,
                 description: str) -> str:
    updated, count = re.subn(pattern, replacement, text, count=0,
                             flags=re.MULTILINE)
    if count != 1:
        raise ValueError(
            f"expected exactly one {description}; found {count}")
    return updated


def derive_compile_input(source: Path, application_name: str,
                         purpose: str = "m2") -> str:
    """Return the reviewed two-species Config.pl input without writing it."""
    text = source.read_text(encoding="utf-8")
    text = replace_once(
        text, r"^(\s*SpeciesList\s*=\s*)Na(\s*!.*)$",
        r"\1Na,Ar\2", "active canonical SpeciesList=Na declaration")

    # Anchor the source edit to its unique process code.  A broad replacement
    # of ``Source:UserDefined=off`` could enable the historical MySource block,
    # which the M0 decision explicitly keeps disabled.
    block_pattern = (
        r"(^\s*Source:UserDefined\s*=\s*)off"
        r"(?=(?:[^\n]*\n){1,3}\s*SourceProcessCode\s*=\s*"
        r"RadiogenicArgon\b)")
    text = replace_once(
        text, block_pattern, r"\1on", "RadiogenicArgon source block")
    if purpose == "m2":
        # Explicitly choose the existing long-period generator for M2's
        # multi-seed/monthly design. This is a derived build setting only;
        # canonical M0/M1 and the sampling-only fixture retain their RNG.
        # RND belongs to ampsConfig.pl's ReadGeneralBlock, not ReadMainBlock.
        # Keeping the option in its actual native block is required even
        # though both sections use the same compile-time key=value spelling.
        text=replace_once(text,r"^(#General[ \t]*)$",
            r"\1\nRND=MERSENNETWISTER ! explicit M2 independent-replicate policy",
            "General block for explicit M2 random generator")
        # I19 explicitly excludes chemistry. The generated macro check below
        # protects this physical contract; make the corresponding source
        # transformation rather than failing after an expensive wrong build.
        text = replace_once(text, r"^(\s*PhotolyticReactions\s*=\s*)on\b",
                            r"\1off", "canonical photolytic reaction switch")
        # The Exosphere application block also defines this macro explicitly.
        # Upstream now inserts such definitions inside their include guard,
        # after the global definition. Select the SAME chemistry contract in
        # both authoritative controls; otherwise Exosphere.h can turn the
        # apparently disabled global mode back ON in a production TU.
        text = replace_once(text,
            r"^(\s*define\s+_PIC_PHOTOLYTIC_REACTIONS_MODE_\s+)"
            r"_PIC_PHOTOLYTIC_REACTIONS_MODE_ON_\b",
            r"\1_PIC_PHOTOLYTIC_REACTIONS_MODE_OFF_",
            "explicit application photolytic reaction macro")
        text = replace_once(text, r"^(\s*ParticleWeightCorrectionMode\s*=\s*)off\b",
                            r"\1on", "finite-release individual particle weights")

    # Config.pl copies files matching the selected application basename into
    # the repository root.  Point the derived input at its equally temporary
    # parser/kernel companions so it never depends on a stale generated
    # ``moon.input.pl`` from an earlier canonical configuration.
    text = replace_once(
        text, r"^#block\s+\.\./moon\.input\.pl\s*.*$",
        f"#block ../{application_name}.input.pl ! generated M2 qualification parser",
        "moon.input.pl block reference")
    text = replace_once(
        text, r"^\s*#include\s+moon\.input\.spicekernels\s*$",
        f"  #include {application_name}.input.spicekernels",
        "moon.input.spicekernels include")
    return text


def production_source_identity():
    """Hash authoritative native inputs, never generated build copies.

    New untracked Moon sources must be included before the user's eventual
    commit. Generated application names and private build paths are excluded
    so a retained binary can be checked after temporary inputs are removed.
    """
    paths=set()
    for base in (REPO/"src",REPO/"srcMoon"):
        if base.exists():
            paths.update(p for p in base.rglob("*") if p.is_file() and
                (p.suffix.lower() in (".cpp",".h",".hpp",".c",".f",".f90",".dfn") or p.name=="makefile"))
    # The shared species include supplies physical masses/charges even though
    # it is not a C++ source file. A changed declaration must invalidate a
    # retained executable just as a changed kernel does; generated header
    # snapshots alone cannot authenticate the original physical input.
    # ampsConfig.pl delegates preprocessing to ampsConfigLib.pm. In particular
    # its macro-guard handling changes generated compilation behavior even
    # when the .pl launcher and physical input bytes are unchanged. Bind this
    # authoritative dependency too; generated snapshots alone cannot identify
    # which preprocessor produced a retained executable.
    paths.update(REPO/name for name in ("Config.pl","ampsConfig.pl","ampsConfigLib.pm",
        "srcMoon/test/build_m2_qualification.py","input/moon.input",
        "input/moon.input.pl","input/moon.input.spicekernels","input/species.input")
        if (REPO/name).is_file())
    return {str(p.relative_to(REPO)):sha256(p) for p in sorted(paths)}


def verify_retained_executable(executable,*,require_current_generated_tree=True):
    """Require bytes, native source and generated chemistry/layout identity."""
    receipt=executable.parent/"build-m2-na-ar.json"
    value=json.loads(receipt.read_text())
    if not isinstance(value,dict) or not isinstance(value.get("executable"),dict):
        raise ValueError("malformed retained native build receipt")
    if value.get("status")!="PASS" or value.get("executable",{}).get("sha256")!=sha256(executable):
        raise ValueError("retained executable does not match its successful build receipt")
    if value.get("production_source_hashes")!=production_source_identity():
        raise ValueError("retained executable native source identity is stale; clean rebuild required")
    contract=value.get("i19_chemistry_contract") or {}
    if not isinstance(contract,dict) or not isinstance(value.get("rng_contract",{}),dict):
        raise ValueError("malformed retained chemistry/RNG configuration contract")
    # A multi-configuration campaign restores ordinary Na-only production
    # after building Na/Ar. Verify immutable generated-header snapshots for
    # a retained binary; do not compare it with an unrelated later build.
    # Every effective-mode proof now requires all four preserved headers.
    # In particular, a live Exosphere.dfn override must never escape the
    # strict probe-ABI check through an empty legacy snapshot map.
    preserved=value.get("retained_generated_files",{})
    if not isinstance(preserved,dict):raise ValueError("malformed retained generated-header map")
    for relative,record in preserved.items():
        if not isinstance(relative,str) or not isinstance(record,dict) or not isinstance(record.get("path"),str):
            raise ValueError("malformed retained generated-header descriptor")
    mandatory={"pic/picGlobal.dfn","general/rnd.h","pic/pic.h","models/exosphere/Exosphere.dfn"}
    if not mandatory<=set(preserved):
        raise ValueError("retained executable lacks the complete generated chemistry/layout/RNG header set")
    def generated(relative,expected):
        record=preserved.get(relative)
        if record:
            path=Path(record["path"])
            if not path.is_absolute():path=receipt.parent/path
            if not path.is_file() or sha256(path)!=record.get("sha256") or record["sha256"]!=expected:
                raise ValueError("retained generated configuration altered: "+relative)
        else:
            raise ValueError("retained executable lacks immutable generated configuration: "+relative)
        if require_current_generated_tree:
            current=REPO/"build"/relative
            if not current.is_file() or sha256(current)!=expected:
                raise ValueError("current generated tree differs from retained configuration: "+relative)
    generated("pic/picGlobal.dfn",contract.get("generated_definitions_sha256"))
    generated("general/rnd.h",value.get("rng_contract",{}).get("generated_header_sha256"))
    for relative,record in preserved.items():
        generated(relative,record.get("sha256"))
    # Textual macro matches in an older receipt did not establish effective
    # production modes after Exosphere.h. Reject that stale ABI in BOTH modes.
    proof=value.get("effective_configuration")
    if (not isinstance(proof,dict) or proof.get("schema")!="srcMoon-M2-effective-configuration-v1" or
            proof.get("status")!="PASS" or proof.get("photolytic_reactions")!="OFF" or
            proof.get("individual_particle_weights")!="ON" or proof.get("rng")!="MERSENNE_TWISTER" or
            proof.get("executable_sha256")!=value["executable"]["sha256"]):
        raise ValueError("retained executable lacks effective compiled M2 configuration proof")
    for name in ("probe_source","macro_dump"):
        record=proof.get(name)
        if not isinstance(record,dict) or not isinstance(record.get("path"),str):
            raise ValueError("malformed effective configuration proof: "+name)
        path=Path(record["path"])
        if not path.is_absolute():path=receipt.parent/path
        if not path.is_file() or sha256(path)!=record.get("sha256"):
            raise ValueError("effective compiled configuration proof altered: "+name)
    return value


def native_preprocessor_prefix(output: Path) -> list:
    """Expand the same compiler and FLAGCC recipe as production make.

    The extra target has no prerequisites and builds nothing. Structured
    subprocess arguments preserve literal make syntax without shell quoting.
    Remove only object-compilation -c; retain optimization, debug, OpenMP,
    wrapper, include and user-supplied preprocessor options unchanged.
    """
    output.mkdir(parents=True,exist_ok=True)
    makefile=output/"configuration-flags.mk"
    makefile.write_text('.PHONY: __moon_m2_configuration_flags\n__moon_m2_configuration_flags:\n\t$(info AMPS_CXX=$(COMPILE.mpicxx))\n\t$(info AMPS_FLAGCC=$(FLAGCC))\n\t@:\n')
    command=["make","--no-print-directory","-s","-f","Makefile","-f",str(makefile.resolve()),
             "OPT3=-O3 -g","__moon_m2_configuration_flags"]
    completed=subprocess.run(command,cwd=REPO,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,check=False)
    (output/"make-expanded-flags.log").write_text(completed.stdout+completed.stderr)
    (output/"make-command.json").write_text(json.dumps(command,indent=2)+"\n")
    rows={}
    for line in completed.stdout.splitlines():
        if line.startswith(("AMPS_CXX=","AMPS_FLAGCC=")):
            key,value=line.split("=",1)
            if key in rows:raise RuntimeError("duplicate expanded native compiler setting: "+key)
            rows[key]=value
    if completed.returncode or not rows.get("AMPS_CXX") or not rows.get("AMPS_FLAGCC"):
        raise RuntimeError("cannot expand native compiler configuration; inspect "+str(output))
    prefix=shlex.split(rows["AMPS_CXX"])+[arg for arg in shlex.split(rows["AMPS_FLAGCC"]) if arg!="-c"]
    (output/"preprocessor-prefix.json").write_text(json.dumps(prefix,indent=2)+"\n")
    return prefix


def effective_configuration(build: Path, spice_path: Path, output: Path,
                            *, compiler: str = "mpicxx", compiler_prefix=None) -> dict:
    """Ask the real preprocessor for modes AFTER the production include chain.

    No executable, model trajectory or scientific oracle is produced here.
    #error evaluates the same aliases/guard ordering as Moon's production TUs;
    an early OFF followed by a later ON fails even though a regex finds OFF.
    """
    output.mkdir(parents=True,exist_ok=True)
    probe=output/"effective-configuration.cpp"
    required=("_PIC_PHOTOLYTIC_REACTIONS_MODE_","_PIC_PHOTOLYTIC_REACTIONS_MODE_OFF_",
        "_INDIVIDUAL_PARTICLE_WEIGHT_MODE_","_INDIVIDUAL_PARTICLE_WEIGHT_ON_",
        "_RND_MODE_","_RND_MODE_MERSENNE_TWISTER_")
    guards="".join('#ifndef '+name+'\n#error "missing effective configuration macro '+name+'"\n#endif\n' for name in required)
    probe.write_text('''// Compiler configuration probe; not a physical test.\n#include "pic.h"\n#include "Moon.h"\n'''+guards+'''#if _PIC_PHOTOLYTIC_REACTIONS_MODE_ != _PIC_PHOTOLYTIC_REACTIONS_MODE_OFF_\n#error "M2 effective photolytic chemistry must be OFF after Moon.h"\n#endif\n#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ != _INDIVIDUAL_PARTICLE_WEIGHT_ON_\n#error "M2 effective individual particle weights must be ON"\n#endif\n#if _RND_MODE_ != _RND_MODE_MERSENNE_TWISTER_\n#error "M2 effective random generator must be Mersenne Twister"\n#endif\n''',encoding="utf-8")
    directories=[build/name for name in ("pic","main","meshAMR","interface","general",
        "models/exosphere","models/electron_impact","models/sputtering","models/dust",
        "models/charge_exchange","models/photolytic_reactions","species")]
    directories += [REPO/"srcInterface", REPO/"share/Library/src", REPO,
                    spice_path/"include"]
    command=list(compiler_prefix) if compiler_prefix is not None else [compiler,"-DMPI_ON","-std=c++17"]
    # The fallback exists for isolated minimal-header protocol fixtures only;
    # the real builder always supplies make's actual expanded recipe above.
    if compiler_prefix is None:
        for directory in directories:command.extend(["-I",str(directory)])
    command.extend(["-E","-dM"])
    command.append(str(probe))
    completed=subprocess.run(command,cwd=REPO,stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,text=True,check=False)
    dump=output/"effective-macros.txt";dump.write_text(completed.stdout)
    (output/"stderr.log").write_text(completed.stderr)
    (output/"command.json").write_text(json.dumps(command,indent=2)+"\n")
    if completed.returncode:
        raise RuntimeError("effective M2 compile-configuration rejected; inspect "+str(output/"stderr.log"))
    return {"schema":"srcMoon-M2-effective-configuration-v1","status":"PASS",
        "photolytic_reactions":"OFF","individual_particle_weights":"ON",
        "rng":"MERSENNE_TWISTER","command":command,
        "probe_source":{"path":str(probe.relative_to(output.parent)),"sha256":sha256(probe)},
        "macro_dump":{"path":str(dump.relative_to(output.parent)),"sha256":sha256(dump)}}


def run_logged(command: list[str], log: Path) -> tuple[int, float]:
    started = time.monotonic()
    with log.open("a", encoding="utf-8") as stream:
        stream.write("command: " + " ".join(command) + "\n")
        stream.flush()
        completed = subprocess.run(
            command, cwd=REPO, stdout=stream, stderr=subprocess.STDOUT,
            check=False)
    return completed.returncode, time.monotonic() - started


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--spice-path", type=Path, required=True)
    parser.add_argument("--spice-kernels", type=Path, required=True)
    parser.add_argument("--purpose", choices=("m2", "sampling"), default="m2",
                        help="sampling builds a two-species storage fixture without claiming the I19 chemistry contract")
    args = parser.parse_args()

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    executable = output / ("amps-moon-m2-na-ar" if args.purpose == "m2" else
                           "amps-moon-sampling-na-ar")
    log = output / "build-m2-na-ar.log"
    result_path = output / "build-m2-na-ar.json"
    # Sandbox/container launches can reuse the same small PID, so a PID alone
    # is not a safe application identity.  Bind the generated name to the
    # absolute campaign path; the runner already requires that path to be new
    # or empty, making the first 12 SHA-256 digits deterministic and collision
    # resistant without introducing an unrecorded random value.
    application_token = hashlib.sha256(
        str(output).encode("utf-8")).hexdigest()[:12]
    application_name = f"moon_m2_qualification_{application_token}"
    source_input = REPO / "input/moon.input"
    staged = {
        REPO / f"input/{application_name}.input": None,
        REPO / f"input/{application_name}.input.pl":
            REPO / "input/moon.input.pl",
        REPO / f"input/{application_name}.input.spicekernels":
            REPO / "input/moon.input.spicekernels",
    }
    # Config.pl copies the selected application's three inputs from input/
    # into the repository root.  They are generated build intermediates, not
    # source files.  Track and remove these exact uniquely named copies as
    # well as the staged inputs so M2 qualification leaves no Python/build
    # debris for a later commit.
    generated_copies = tuple(
        REPO / f"{application_name}{suffix}"
        for suffix in (".input", ".input.pl", ".input.spicekernels")
    )
    assembled_pattern = f"{application_name}.input.Assembled*"
    result: dict[str, object] = {
        "schema": "srcMoon-M2-Na-Ar-build-v1",
        "purpose": args.purpose,
        "status": "ERROR",
        "source_input": {"path": str(source_input),
                         "sha256": sha256(source_input)},
        "compile_time_changes": [
            "SpeciesList=Na -> SpeciesList=Na,Ar",
            "RadiogenicArgon Source:UserDefined=off -> on",
        ],
        "commands": [],
        "log": str(log),
    }
    # Cleanup ownership starts only after every temporary name is verified
    # absent. If a previous/user-created file causes the preflight to refuse,
    # finally must preserve it rather than remove the very collision reported.
    cleanup_owned = False
    if args.purpose == "m2":
        result["compile_time_changes"].extend([
            "PhotolyticReactions=on -> off (I19 chemistry contract)",
            "ParticleWeightCorrectionMode=off -> on (exact partial release packets)"])
    try:
        for path in (*staged, *generated_copies,
                     *tuple(REPO.glob(assembled_pattern))):
            if path.exists():
                raise FileExistsError(
                    f"refusing to overwrite temporary application file {path}")
        cleanup_owned = True
        derived = derive_compile_input(source_input, application_name,args.purpose)
        # The staged files are build inputs, not validation evidence.  Their
        # exact hashes are retained before cleanup in the machine-readable
        # build receipt.
        (REPO / f"input/{application_name}.input").write_text(
            derived, encoding="utf-8")
        for destination, source in staged.items():
            if source is not None:
                shutil.copy2(source, destination)
        result["staged_input_hashes"] = {
            path.name: sha256(path) for path in staged
        }

        build = REPO / "build"
        if build.exists():
            shutil.rmtree(build)
        commands = [
            # Select production mode explicitly: Config.pl otherwise retains
            # a previous nightly build's sorting/debug/iteration policies.
            [str(REPO / "Config.pl"), f"-application={application_name}", "-amps-test=off",
             f"-spice-path={args.spice_path.resolve()}",
             f"-spice-kernels={args.spice_kernels.resolve()}"],
            # This is the same source-generation command used by make.
            # Generate first so the compiler can reject contradictory modes
            # before any expensive object compilation occurs.
            [str(REPO/"ampsConfig.pl"),"-input",f"{application_name}.input","-no-compile"],
            # The current core native owners use sixteen build jobs. Keep
            # debug symbols/optimization explicit while bounding compilers.
            ["make", "-j16", "OPT3=-O3 -g"],
        ]
        log.write_text("", encoding="utf-8")
        for command in commands:
            returncode, duration = run_logged(command, log)
            result["commands"].append({
                "command": command, "returncode": returncode,
                "duration_seconds": duration})
            if returncode != 0:
                raise RuntimeError(
                    f"M2 qualification build command failed ({returncode}): "
                    + " ".join(command))
            if args.purpose=="m2" and Path(command[0]).name=="ampsConfig.pl":
                result["prebuild_effective_configuration"]=effective_configuration(
                    build,args.spice_path.resolve(),output/"prebuild-effective-configuration",
                    compiler_prefix=native_preprocessor_prefix(output/"prebuild-native-flags"))

        # I19's approved M2 accounting case excludes chemistry so its
        # interspecies conversion terms have a physical expectation of
        # exactly zero.  Verify the *generated* contract instead of assuming
        # that a human-readable Config.pl input spelling maps to a particular
        # macro.  A future configuration change that enables photochemistry
        # must fail this dedicated build rather than silently alter I19's
        # closure equation.
        generated_definitions = REPO / "build/pic/picGlobal.dfn"
        definitions_text = generated_definitions.read_text(encoding="utf-8")
        chemistry_pattern = (
            r"^#define\s+_PIC_PHOTOLYTIC_REACTIONS_MODE_\s+"
            r"_PIC_PHOTOLYTIC_REACTIONS_MODE_OFF_\s*$"
        )
        if (args.purpose == "m2" and
                re.search(chemistry_pattern, definitions_text,
                          flags=re.MULTILINE) is None):
            raise RuntimeError(
                "dedicated M2 build violates the approved I19 contract: "
                "photolytic chemistry is not compiled OFF")
        if args.purpose == "m2" and re.search(
                r"^#define\s+_INDIVIDUAL_PARTICLE_WEIGHT_MODE_\s+_INDIVIDUAL_PARTICLE_WEIGHT_ON_\s*$",
                definitions_text,flags=re.MULTILINE) is None:
            raise RuntimeError("finite M2 release requires generated individual particle weights ON")
        # A sampling-only fixture needs the real two-species storage layout,
        # but performs no chemistry or transport advances. It must not claim
        # I19's chemistry-off configuration when the normal build enables it.
        result["i19_chemistry_contract"] = None if args.purpose == "sampling" else {
            "photolytic_reactions": "OFF",
            "converted_in_atoms": 0.0,
            "converted_out_atoms": 0.0,
            "generated_definitions": str(generated_definitions),
            "generated_definitions_sha256": sha256(generated_definitions),
        }
        if args.purpose=="m2":
            result["effective_configuration"]=effective_configuration(
                build,args.spice_path.resolve(),output/"effective-configuration",
                compiler_prefix=native_preprocessor_prefix(output/"native-flags"))
            rng=REPO/"build/general/rnd.h"
            if not re.search(r"^#define\s+_RND_MODE_\s+_RND_MODE_MERSENNE_TWISTER_\s*$",rng.read_text(),re.MULTILINE):
                raise RuntimeError("M2 generated RNG is not the explicitly selected Mersenne Twister")
            result["rng_contract"]={"mode":"MERSENNE_TWISTER","generated_header_sha256":sha256(rng),
                "rank_seed_policy":"configured positive base + rank; distinct replicate/rank labels; reviewed statistical approximation, not proof of independence"}
            # Preserve the actual generated chemistry/layout and RNG headers
            # beside this executable before another profile removes build/.
            # Later formal imports verify these immutable bytes as well as
            # current authoritative source hashes and executable identity.
            result["retained_generated_files"]={}
            for relative in ("pic/picGlobal.dfn","general/rnd.h","pic/pic.h","models/exosphere/Exosphere.dfn"):
                source=REPO/"build"/relative
                if not source.is_file():raise RuntimeError("missing generated header: "+relative)
                target=output/"generated-configuration"/relative
                target.parent.mkdir(parents=True,exist_ok=True)
                shutil.copy2(source,target)
                result["retained_generated_files"][relative]={"path":str(target.relative_to(output)),"sha256":sha256(target)}
        result["production_source_hashes"]=production_source_identity()
        shutil.copy2(REPO / "amps", executable)
        if args.purpose=="m2":
            result["effective_configuration"]["executable_sha256"]=sha256(executable)
        result.update(
            status="PASS",
            executable={"path": str(executable),
                        "sha256": sha256(executable)},
        )
    except (OSError, ValueError, RuntimeError) as error:
        result["error"] = str(error)
    finally:
        cleanup_errors: list[str] = []
        # Config.pl determines the exact set of `.Assembled*` derivatives, so
        # enumerate only this uniquely named application's prefix after the
        # build. No broad repository cleanup or unrelated generated file is
        # touched.
        cleanup_paths = ((*staged, *generated_copies,
                          *tuple(REPO.glob(assembled_pattern)))
                         if cleanup_owned else ())
        for path in cleanup_paths:
            try:
                path.unlink(missing_ok=True)
            except OSError as error:
                cleanup_errors.append(f"{path}: {error}")
        if cleanup_errors:
            result["status"] = "ERROR"
            result["cleanup_errors"] = cleanup_errors
        result_path.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n",
            encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == "PASS" else 2

if __name__ == "__main__":
    raise SystemExit(main())
