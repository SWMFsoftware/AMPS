#!/usr/bin/env python3
"""Exercise M2 control boundaries through the compiled PRODUCTION .in parser.

These are configuration fixtures, not observations, force laws or scientific
validation. The thin probe calls ParseApplicationInput and exits before model
initialization. In particular an accepted path is not proof its data are READY.
The aggregate development driver runs this after the clean Na/Ar build.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import subprocess
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/"common"))
from moon_testlib import compile_probe
from m2_illumination import setting,digest,resolve


def change(text,key,value):
    """Edit a fixture; this is not a second production syntax parser."""
    pattern=r"(?m)^\s*"+re.escape(key)+r"\s*=[^\n]*$"
    if re.search(pattern,text):return re.sub(pattern,key+" = "+str(value),text)
    endings=list(re.finditer(r"(?m)^\s*#end application[^\n]*",text))
    if len(endings)!=1:raise ValueError("fixture requires one Moon block")
    point=endings[0].start()
    return text[:point]+key+" = "+str(value)+"\n"+text[point:]


def classify_invalid_psd(returncode, output):
    """Separate an evaluated guard miss from a broken fatal-path adapter.

    A signal, an MPI-before-init rejection or a missing declared diagnostic is
    infrastructure ERROR, even when some expected words appeared before the
    failure. Only the actual PSD refusal, with initialized MPI and a nonzero
    process status, is a successful negative production-path check.
    """
    if returncode<0:
        return "ERROR","negative PSD child terminated by signal "+str(-returncode)
    if "before MPI_INIT" in output or "before MPI_Init" in output:
        return "ERROR","negative PSD adapter reached MPI_Abort before initializing MPI"
    if returncode==77 and "requires compiled positive Na PSD coefficients" in output:
        return "SKIPPED","positive compiled Na PSD coefficients unavailable"
    if returncode==0 or "psd_nonfinite_factor_was_accepted" in output:
        return "FAIL","invalid PSD photon factor was accepted or the process reported success"
    if ("direct photon-flux factor must be finite and nonnegative" in output and
            "native_negative_probe_mpi_initialized=1" in output):
        return "PASS","production PSD guard rejected the factor through the initialized fatal path"
    return "ERROR","negative PSD child lacks the declared initialized-MPI fatal diagnostic"


def run_checks(build,input_file,output):
    output=Path(output).resolve();output.mkdir(parents=True,exist_ok=True)
    input_file=Path(input_file).resolve()
    probe,compilation=compile_probe(Path(build).resolve(),output)
    result={"schema":"srcMoon-M2-native-control-tests-v1","scope":"supporting-production-parser-not-validation",
            "status":"ERROR","checks":[],"compile":compilation,
            "input_template":{"path":str(input_file),"sha256":digest(input_file)}}
    if probe is None:
        result.update(status=compilation["status"],reason="compiled native probe unavailable")
    else:
        text=input_file.read_text()
        for key in ("spice_path","lola_image_file","lola_label_file"):
            text=change(text,key,str(resolve(input_file,setting(input_file,key))))
        finite=change(text,"illumination_terrain_solar_model","uniform-disc")
        finite=change(finite,"illumination_disc_radial_points",4)
        finite=change(finite,"illumination_disc_azimuth_points",16)
        # This file is only an existence/parser fixture: the native probe
        # never initializes a model or reads a temperature value from it.
        # It must never be promoted to a Diviner observation or READY dataset.
        table=output/"protocol-only-not-observations.csv"
        table.write_text("PROTOCOL FIXTURE ONLY; NO OBSERVATIONS\n")
        diviner=change(text,"surface_temperature_mode","diviner-map")
        diviner=change(diviner,"diviner_temperature_table_file",str(table))
        diviner_missing=re.sub(r"(?m)^\s*diviner_map_kinetic_interpretation\s*=[^\n]*\n?","",diviner)
        variants=[("legacy-default-preserved",text,True),
                  ("explicit-disc-valid",finite,True),
                  ("disc-grid-missing",change(text,"illumination_terrain_solar_model","uniform-disc"),False),
                  ("disc-negative-ring",change(finite,"illumination_disc_radial_points",-1),False),
                  ("disc-zero-ring",change(finite,"illumination_disc_radial_points",0),False),
                  ("disc-asymmetric-grid",change(finite,"illumination_disc_azimuth_points",7),False),
                  ("disc-noninteger-grid",change(finite,"illumination_disc_radial_points",1.5),False),
                  ("ignored-disc-controls",change(finite,"illumination_terrain_solar_model","point-sun"),False),
                  ("unknown-disc-physics",change(finite,"illumination_terrain_solar_model","limb-darkening"),False),
                  ("disabled-environment",change(finite,"m2_environment","disabled"),False),
                  ("unknown-aberration",change(text,"illumination_aberration_correction","guess"),False),
                  ("unreviewed-holdout",change(text,"diviner_eclipse_selection","holdout"),False),
                  ("unknown-selection",change(text,"diviner_eclipse_selection","december-fit"),False),
                  ("diviner-map-missing-interpretation",diviner_missing,False),
                  ("diviner-map-unknown-interpretation",change(diviner,"diviner_map_kinetic_interpretation","brightness-is-kinetic"),False),
                  ("diviner-map-explicit-grey-valid",change(diviner,"diviner_map_kinetic_interpretation","isothermal-grey"),True),
                  ("case-preserving-relative-contract",change(text,"m2_formal_contract_file","Evidence/M2.json"),True)]
        for name,body,accepted in variants:
            directory=output/name;directory.mkdir()
            source=directory/"application.in";source.write_text(body)
            command=[str(probe),"parse-m2-input",str(source)]
            process=subprocess.run(command,text=True,capture_output=True,check=False)
            (directory/"stdout.log").write_text(process.stdout)
            (directory/"stderr.log").write_text(process.stderr)
            expected=0 if accepted else 4
            ok=process.returncode==expected
            if name=="explicit-disc-valid":
                ok &= "radial_points=4\n" in process.stdout and "azimuth_points=16\n" in process.stdout
            if name=="case-preserving-relative-contract":
                ok &= "formal_contract="+str(directory/"Evidence/M2.json")+"\n" in process.stdout
            result["checks"].append({"name":name,"status":"PASS" if ok else "FAIL",
                "expected_returncode":expected,"returncode":process.returncode,
                "command":command,"input_sha256":digest(source),"log":str(directory/"stdout.log")})
        # Exercise the real core PSD refusal, not a copied Python validator.
        # A crash without the declared fatal diagnostic is not a successful
        # negative test; no observation or physical reference is fabricated.
        invalid=output/"invalid-psd-factor";invalid.mkdir()
        command=[str(probe),"psd-illumination-invalid"]
        negative={"name":"invalid-psd-factor","status":"ERROR","command":command,
                  "working_directory":str(invalid),"log":str(invalid/"stderr.log")}
        try:
            # Config.pl substitutes the diagnostic filename in this generated
            # source. Read the ACTUAL compiled-build convention, rather than
            # guessing that a bare kernel probe needs no runtime directories.
            # The retained failed build used PT/amps_error.log. Keep relative
            # log paths inside this child-owned directory; an absolute/external
            # log would violate isolation and needs an explicit adapter change.
            log_source=Path(build).resolve()/"general/specfunc.cpp"
            log_names=set(re.findall(r'fopen\(\s*"([^"\\\n]+)"\s*,\s*"a\+"\s*\)',log_source.read_text()))
            if len(log_names)!=1:
                raise ValueError("cannot identify one compiled diagnostic log path in generated specfunc.cpp")
            log_name=Path(next(iter(log_names)))
            if log_name.is_absolute() or ".." in log_name.parts:
                raise ValueError("negative probe requires a relative, isolated configured diagnostic log path")
            log_path=invalid/log_name
            log_path.parent.mkdir(parents=True,exist_ok=True)
            # Prove the diagnostic can be opened before asking a fatal handler
            # to write it. The ordinary core logger is intentionally unchanged.
            with log_path.open("a"):
                pass
            negative["diagnostics_configuration"]={"source":str(log_source),
                "source_sha256":digest(log_source),"relative_log_path":str(log_name)}
            process=subprocess.run(command,cwd=invalid,text=True,capture_output=True,check=False)
            (invalid/"stdout.log").write_text(process.stdout)
            (invalid/"stderr.log").write_text(process.stderr)
            status,reason=classify_invalid_psd(process.returncode,process.stdout+process.stderr)
            negative.update(status=status,reason=reason,returncode=process.returncode)
        except (OSError,ValueError) as error:
            negative.update(reason=str(error))
            (invalid/"stderr.log").write_text(str(error)+"\n")
        result["checks"].append(negative)
        states={x["status"] for x in result["checks"]}
        result["status"]=("ERROR" if "ERROR" in states else "FAIL" if "FAIL" in states else
                          "SKIPPED" if "SKIPPED" in states else "PASS")
        result["counts"]={state:sum(x["status"]==state for x in result["checks"]) for state in ("PASS","FAIL","ERROR","SKIPPED")}
        result["probe"]={"path":str(probe),"sha256":digest(probe)}
    (output/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir",required=True,type=Path)
    parser.add_argument("--input",required=True,type=Path)
    parser.add_argument("--output-dir",required=True,type=Path)
    args=parser.parse_args()
    if args.output_dir.exists() and any(args.output_dir.iterdir()):parser.error("output must be new/empty")
    result=run_checks(args.build_dir,args.input,args.output_dir)
    print(result["status"],"M2 native input controls",result.get("counts",{}))
    return {"PASS":0,"FAIL":1,"ERROR":2,"SKIPPED":77}[result["status"]]


if __name__=="__main__":raise SystemExit(main())
