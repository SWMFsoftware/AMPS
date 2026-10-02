#!/usr/bin/env python3
"""Compile the actual AMPS scheduler/getter bodies without the MPI mesh build.

The fixture replaces mesh allocation, MPI and file I/O, not the code being
tested. It extracts the production functions from the supplied checkout and
compiles them twice (file interpolation OFF and ON). Empty-schedule accesses
are caught by vector assertions and optionally address/undefined sanitizers.
This is a portable regression check, not evidence of a native MPI run.
"""

import argparse
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile


def cpp_block(text, pattern):
    """Return one unique C++ definition, ignoring braces in strings/comments."""
    matches = list(re.finditer(pattern, text))
    if len(matches) != 1:
        raise ValueError("expected one production anchor for " + pattern)
    start = matches[0].start()
    brace = text.index("{", matches[0].end())
    # Tokenization prevents a diagnostic string containing a brace from
    # changing the boundary of the extracted production definition.
    tokens = re.finditer(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"|'
                         r"'(?:\\.|[^'\\])*'|[{}]", text[brace:])
    depth = 0
    for token in tokens:
        value = token.group()
        if value == "{":
            depth += 1
        elif value == "}":
            depth -= 1
            if depth == 0:
                return text[start:brace + token.end()]
    raise ValueError("unclosed production definition")


FIXTURE = r'''
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#define _TARGET_HOST_
#define _TARGET_DEVICE_
#define _CUDA_MANAGED_
#define _PIC_MODE_ON_ 1
#define _PIC_MODE_OFF_ 0
#define _PIC_COUPLER_MODE_ 2
#define _PIC_COUPLER_MODE__DATAFILE_ 2
#define _PIC_FIELD_LINE_MODE_ 1
#define _PIC_TIMESTEP_RETURN_CODE__END_SIMULATION_ -7
#define _PIC_DEBUGGER_MODE_ 1
#define _PIC_DEBUGGER_MODE_ON_ 1
#define _PIC_DEBUGGER_MODE__CHECK_FINITE_NUMBER_ 1
#include "pic_background_update_mode.h"
using std::vector;
using std::isfinite;
int assertions=0, imports=0, schedules=0, copies=0, fieldLines=0, finiteChecks=0;
void Check(bool value,const char* message) {
  ++assertions;
  if (!value) throw std::runtime_error(message);
}
// Production AMPS uses its own three-argument exit overload. Throw in the
// fixture so a rejected direct loader call can be checked without ending it.
void exit(int,const char*,const char* message) { throw std::runtime_error(message); }
namespace PIC {
int ThisThread=1;
namespace SimulationTime {
double now=0;
double Get() { return now; }
void SetInitialValue(double value) { now=value; }
void Update() { now+=60; }
}
namespace FieldLine { void Update() { ++fieldLines; } }
namespace Debugger {
void CatchOutLimitValue(double* value,int n,int,const char*) {
  ++finiteChecks;
  for (int i=0;i<n;++i) if (!isfinite(value[i])) throw std::runtime_error("nonfinite field");
}
}
namespace Mesh {
struct cDataCenterNode {
  // Two aligned 4-double slots suffice for three-component B/U/E reads.
  double storage[8]={2,3,4,0,10,11,12,0};
  char* GetAssociatedDataBufferPointer() { return reinterpret_cast<char*>(storage); }
};
}
namespace CPLR { namespace DATAFILE {
int CenterNodeAssociatedDataOffsetBegin=0, nTotalBackgroundVariables=4;
void ImportData(const char*) { ++imports; }
namespace MULTIFILE {
struct cScheduleItem { double Time; const char* FileName; };
vector<cScheduleItem> Schedule;
int nFile=0,iFileLoadNext=-1,CurrDataFileOffset=0,NextDataFileOffset=-1;
bool ReachedLastFile=false,BreakAtLastFile=true;
void GetSchedule() {
  ++schedules; Schedule={{0,"first"},{120,"second"},{240,"third"}}; nFile=3;
}
void CopyCurrDataFile2NextDataFile() { ++copies; }
void Init(bool,int);
void UpdateDataFile();
@IS_TIME@
}
@GETTER@
} }
}
@POLICY@
@INIT@
@UPDATE@
int Step() {
  // The extracted coupler block is used intact, including EOF termination and
  // the field-line callback. Model ownership must bypass only file lifecycle.
  PIC::SimulationTime::Update();
  @STEP@
  return 0;
}
int main() {
  try {
    using namespace PIC::CPLR::DATAFILE;
    using namespace PIC::CPLR::DATAFILE::MULTIFILE;
    Check(UsesFileSchedule(),"legacy default changed");
    Init(true,0);
    const int initialImports=(_PIC_DATAFILE__TIME_INTERPOLATION_MODE_==1)?2:1;
    Check(schedules==1 && imports==initialImports && copies==1,"file initialization changed");
    Check(PIC::SimulationTime::Get()==0,"file epoch initialization changed");
    // A real file sequence still advances and updates field lines.
    PIC::SimulationTime::now=70;
    Check(Step()==0,"file step failed");
    Check(imports==initialImports+1 && fieldLines==1,"file loading/field lines changed");
    ReachedLastFile=true; BreakAtLastFile=true;
    Check(Step()==-7,"legacy EOF termination changed");
    Check(fieldLines==1,"EOF ordering changed");
    BreakAtLastFile=false;
    Check(Step()==0 && fieldLines==2,"legacy continue-after-EOF changed");

    PIC::Mesh::cDataCenterNode cell;
    double values[3]={};
    CurrDataFileOffset=0; NextDataFileOffset=4*sizeof(double);
    Schedule={{0,"first"},{120,"second"}}; nFile=2; iFileLoadNext=2;
    ReachedLastFile=false;
    GetBackgroundValue(values,3,0,&cell,60);
    Check(values[0]==(_PIC_DATAFILE__TIME_INTERPOLATION_MODE_==1?6:2),"file interpolation changed");
    ReachedLastFile=true;
    GetBackgroundValue(values,3,0,&cell,60);
    Check(values[0]==(_PIC_DATAFILE__TIME_INTERPOLATION_MODE_==1?6:2),"file EOF interpolation changed");

    BackgroundUpdatePolicy=BackgroundUpdateMode::RuntimeProvider;
    Schedule.clear(); nFile=0; iFileLoadNext=-1;
    CurrDataFileOffset=0; NextDataFileOffset=-1;
    // Poison the next slot. Runtime reads must use the complete current slot,
    // even if file interpolation was compiled ON or a caller supplies NaN time.
    for (int i=4;i<8;++i) cell.storage[i]=NAN;
    for (int reached=0;reached<2;++reached) {
      for (int shouldBreak=0;shouldBreak<2;++shouldBreak) {
        ReachedLastFile=reached; BreakAtLastFile=shouldBreak;
        PIC::SimulationTime::now=0;
        const int oldImports=imports,oldLines=fieldLines;
        cell.storage[0]=2; cell.storage[1]=3; cell.storage[2]=4;
        Check(!IsTimeToUpdate(),"runtime consulted a file schedule");
        Check(Step()==0,"runtime step 1 ended at file EOF");
        GetBackgroundValue(values,3,0,&cell,NAN);
        Check(values[0]==2 && values[1]==3 && values[2]==4,"runtime current field read failed");
        // Emulate the application's joined publication at the step boundary;
        // the next read must observe the new epoch, not interpolate file times.
        cell.storage[0]=7; cell.storage[1]=8; cell.storage[2]=9;
        Check(Step()==0,"runtime step 2 ended at file EOF");
        GetBackgroundValue(values,3,0,&cell,120);
        Check(values[0]==7 && values[1]==8 && values[2]==9,"runtime publication invisible");
        Check(PIC::SimulationTime::Get()==120,"runtime clock progression changed");
        Check(imports==oldImports && fieldLines==oldLines+2,"runtime loader/field-line ownership failed");
      }
    }
    const int oldSchedules=schedules,oldImports=imports;
    bool rejectedInit=false,rejectedUpdate=false;
    try { Init(true,0); } catch (const std::runtime_error& e) {
      rejectedInit=std::string(e.what()).find("runtime-owned")!=std::string::npos;
    }
    try { UpdateDataFile(); } catch (const std::runtime_error& e) {
      rejectedUpdate=std::string(e.what()).find("runtime-owned")!=std::string::npos;
    }
    Check(rejectedInit && rejectedUpdate,"direct runtime file calls were not rejected");
    Check(schedules==oldSchedules && imports==oldImports,"rejected loader touched file I/O");
    Check(PIC::SimulationTime::Get()==120 && CurrDataFileOffset==0 && NextDataFileOffset==-1,
          "rejected loader changed clock/slots");
    Check(finiteChecks==10,"field finite checks were bypassed");
    std::cout << "PASS interpolation=" << _PIC_DATAFILE__TIME_INTERPOLATION_MODE_
              << " assertions=" << assertions << "\n";
  } catch (const std::exception& e) {
    std::cerr << "FAIL: " << e.what() << "\n"; return 1;
  }
}
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amps-root", type=Path,
                        default=Path(__file__).resolve().parents[2])
    parser.add_argument("--cxx", default="g++")
    parser.add_argument("--sanitize", action="store_true",
                        help="enable address and undefined-behavior sanitizers")
    parser.add_argument("--output", type=Path, help="write machine-readable results")
    args = parser.parse_args()
    root = args.amps_root.resolve()
    h = (root / "src/pic/pic.h").read_text()
    data = (root / "src/pic/pic_datafile.cpp").read_text()
    step = (root / "src/pic/pic.cpp").read_text()
    app = (root / "srcSEP3D/main_lib.cpp").read_text()
    # Ordering is a production integration contract: native getters used during
    # initial publication/output must already know there is no file schedule.
    init = cpp_block(app, r"void amps_init\(\)")
    if not (init.index("PIC::Init_AfterParser();") <
            init.index("BackgroundUpdateMode::RuntimeProvider") <
            init.index("\n  FillAndPublishBackground();")):
        raise ValueError("SEP3D did not claim ownership before publication")
    replacements = {
        "@IS_TIME@": cpp_block(h, r"inline bool IsTimeToUpdate\(\)"),
        "@GETTER@": cpp_block(h, r"inline void GetBackgroundValue\([^\n]+\)"),
        "@INIT@": cpp_block(data, r"void PIC::CPLR::DATAFILE::MULTIFILE::Init\([^\n]+\)"),
        "@UPDATE@": cpp_block(data, r"void PIC::CPLR::DATAFILE::MULTIFILE::UpdateDataFile\(\)"),
        "@STEP@": cpp_block(step, r"if \(_PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__DATAFILE_\)"),
    }
    policy = re.findall(r"PIC::CPLR::DATAFILE::BackgroundUpdateMode\s+_TARGET_DEVICE_\s+"
                        r"_CUDA_MANAGED_\s+PIC::CPLR::DATAFILE::BackgroundUpdatePolicy\s*=\s*"
                        r"PIC::CPLR::DATAFILE::BackgroundUpdateMode::FileSchedule;", data)
    if len(policy) != 1:
        raise ValueError("missing production file-schedule default definition")
    replacements["@POLICY@"] = policy[0]
    fixture = FIXTURE
    for anchor, replacement in replacements.items():
        fixture = fixture.replace(anchor, replacement)
    results = []
    # Build in temporary storage; this check neither configures AMPS nor leaves
    # generated compiler products in the production source tree.
    with tempfile.TemporaryDirectory(prefix="sep3d-runtime-coupling-") as tmp:
        source = Path(tmp) / "fixture.cpp"
        source.write_text(fixture)
        for interpolation in (0, 1):
            executable = Path(tmp) / ("test-" + str(interpolation))
            command = shlex.split(args.cxx) + [
                "-std=c++11", "-O1", "-g", "-Wall", "-Wextra", "-Werror",
                "-Wno-unused-parameter",
                "-D_GLIBCXX_ASSERTIONS",
                "-D_PIC_DATAFILE__TIME_INTERPOLATION_MODE_=" + str(interpolation),
                "-I" + str(root / "src/pic"), str(source), "-o", str(executable)]
            if args.sanitize:
                command += ["-fsanitize=address,undefined", "-fno-omit-frame-pointer"]
            subprocess.run(command, check=True)
            environment = os.environ.copy()
            if args.sanitize:
                # The fixture checks invalid access/undefined behavior. Leak
                # scanning is unrelated and cannot inspect /proc task lists in
                # some managed runners; keep ASan's memory-access checks active.
                environment["ASAN_OPTIONS"] = environment.get("ASAN_OPTIONS", "") + ":detect_leaks=0"
            run = subprocess.run([str(executable)], capture_output=True, text=True,
                                 env=environment)
            if run.returncode:
                # Preserve the actual assertion/sanitizer diagnostic rather
                # than reporting only the temporary executable's exit status.
                raise RuntimeError(run.stdout + run.stderr)
            print(run.stdout, end="")
            results.append({"interpolation": interpolation, "status": "PASS",
                            "output": run.stdout.strip(), "sanitizers": args.sanitize})
    report = {"scope": "extracted production core functions; MPI/mesh/file-I/O stubs",
              "native_mpi_executed": False, "results": results}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
