#!/usr/bin/env python3
"""Reproduce unsupported-configuration REFUSAL controls, never those capabilities.

Only a new compile-input file is written. Configuration/build remain the
normal AMPS commands and must be run in an isolated checkout when another
application is configured. The alternative check action interprets the real
native exit/log, rejecting unrelated crashes and incomplete refusals.
"""
import argparse,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
CONTROLS={
    'periodic':('_PIC_BC__PERIODIC_MODE_ _PIC_BC__PERIODIC_MODE_ON_','periodic wrapping is not implemented'),
    'backward':('_PIC_PARTICLE_MOVER__BACKWARD_TIME_INTEGRATION_MODE_ _PIC_PARTICLE_MOVER__BACKWARD_TIME_INTEGRATION_MODE__ENABLED_','backward integration is not implemented'),
    'field-solver':('_PIC_FIELD_SOLVER_MODE_ _PIC_FIELD_SOLVER_MODE__ELECTROMAGNETIC__ECSIM_','self-consistent PIC field solving/current-charge deposition is not implemented'),
}
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',required=True,choices=sorted(CONTROLS))
    action=parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--output',type=Path,help='NEW compile input; existing files are never overwritten')
    action.add_argument('--check-log',type=Path,help='actual native stdout/stderr from the matching producer')
    parser.add_argument('--exit-code',type=int,help='actual mpiexec return code, mandatory with --check-log')
    args=parser.parse_args();definition,diagnostic=CONTROLS[args.mode]
    if args.output:
        template=ROOT/'input/test/force-cut-surface-mover.input';text=template.read_text()
        if text.count('#endGeneral')!=1:raise ValueError('expected exactly one General block terminator')
        text=text.replace('#endGeneral','  ! Unsupported-configuration refusal control; no capability qualification.\n  define '+definition+'\n#endGeneral')
        args.output.parent.mkdir(parents=True,exist_ok=True)
        with args.output.open('x') as stream:stream.write(text)
        print(json.dumps(dict(scope='configuration refusal only',mode=args.mode,input=str(args.output.resolve()),template_sha256=hashlib.sha256(template.read_bytes()).hexdigest(),input_sha256=hashlib.sha256(args.output.read_bytes()).hexdigest(),expected_native_exit_code=1,expected_diagnostic=diagnostic),indent=2));return 0
    if args.exit_code is None:parser.error('--exit-code is mandatory with --check-log')
    text=args.check_log.read_text(errors='replace')
    accepted=args.exit_code==1 and diagnostic in text and 'Segmentation fault' not in text
    print(json.dumps(dict(status='PASS' if accepted else 'ERROR',scope='configuration refusal only; capability remains unsupported',mode=args.mode,native_exit_code=args.exit_code,log=str(args.check_log.resolve()),expected_diagnostic=diagnostic),indent=2))
    return 0 if accepted else 2
if __name__=='__main__':raise SystemExit(main())
