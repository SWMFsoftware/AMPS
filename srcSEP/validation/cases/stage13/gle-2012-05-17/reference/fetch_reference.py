#!/usr/bin/env python3
"""Reacquire this case's small official products without changing references.

Only artifacts with an explicit HTTPS download_url in the reviewed case.json
are requested. The caller must provide a fresh evidence directory; this tool
never overwrites committed observations, substitutes a response, interpolates
missing samples or changes validation readiness. An HTTP200 response is not
sufficient: expected size and SHA256 must match before atomic publication.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--source-cache',type=Path,help='verified prior direct-download files; records cache reproduction, never a new network acquisition')
    args=parser.parse_args()
    case_dir=Path(__file__).resolve().parents[1]
    manifest_path=case_dir/'case.json'
    # Pin the exact declaration consumed by this acquisition. The receipt
    # distinguishes source availability from a successful model comparison.
    manifest_bytes=manifest_path.read_bytes()
    manifest=json.loads(manifest_bytes)
    output=args.output_dir.resolve()
    if output.exists():
        parser.error('a fresh output directory is required; existing evidence is preserved')
    output.mkdir(parents=True)
    records=[]
    for artifact in manifest['artifacts']:
        url=artifact.get('download_url')
        if url is None:
            continue
        if not url.startswith('https://'):
            raise ValueError('reference downloads require an explicit HTTPS source')
        destination=output/Path(artifact['path']).name
        temporary=destination.with_name(destination.name+'.partial')
        # The2MB ceiling enforces the authorized small-reference scope. It is
        # not a scientific cutoff: larger products require separate acquisition.
        command=['curl','-fsSL','--proto','=https','--proto-redir','=https',
                 '--max-time','60','--max-filesize','2000000',url,'-o',str(temporary)]
        if args.source_cache is None:
            completed=subprocess.run(command,text=True,stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT,check=False)
            actual_code=completed.returncode
            diagnostic=completed.stdout
            payload=temporary.read_bytes() if temporary.exists() else b''
            method='direct-HTTPS-download'
        else:
            # Restricted runners may allow a top-level curl acquisition while
            # blocking DNS in Python children. Replay the separately acquired
            # bytes explicitly; this is not reported as a successful network
            # request. The same size/hash requirements apply before publication.
            cached=args.source_cache.resolve()/destination.name
            payload=cached.read_bytes() if cached.is_file() else b''
            actual_code=None
            diagnostic='explicit source-cache replay: '+str(cached)
            method='checksum-verified-source-cache-replay'
        digest=hashlib.sha256(payload).hexdigest()
        verified=((actual_code==0 or args.source_cache is not None)
                  and len(payload)==artifact['bytes'] and digest==artifact['sha256'])
        if verified and args.source_cache is not None:
            temporary.write_bytes(payload)
        record={'url':url,'command':command if args.source_cache is None else None,
                'method':method,'actual_exit_code':actual_code,
                'bytes':len(payload),'sha256':digest,'expected_bytes':artifact['bytes'],
                'expected_sha256':artifact['sha256'],'status':'PASS' if verified else 'ERROR',
                'diagnostic':diagnostic,'published_path':str(destination) if verified else None}
        if verified:
            temporary.replace(destination)
        # A failed/changed download stays explicitly unqualified. Do not turn
        # missing data or a new archive version into a silently updated baseline.
        records.append(record)
    unchanged=manifest_path.read_bytes()==manifest_bytes
    good=bool(records) and unchanged and all(r['status']=='PASS' for r in records)
    receipt={'schema':'srcsep-stage13-reference-acquisition-v1','case_id':manifest['case_id'],
             'manifest_sha256':hashlib.sha256(manifest_bytes).hexdigest(),
             'manifest_unchanged':unchanged,'results':records,
             'status':'PASS' if good else 'ERROR','comparison_status':'NOT_RUN'}
    (output/'acquisition.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'case_id':manifest['case_id'],'status':receipt['status'],
                      'comparison_status':'NOT_RUN','receipt':str(output/'acquisition.json')}))
    return 0 if good else 2


if __name__=='__main__':
    raise SystemExit(main())
