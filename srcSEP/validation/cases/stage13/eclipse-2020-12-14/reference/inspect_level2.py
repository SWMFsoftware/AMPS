#!/usr/bin/env python3
"""Verify and inspect the exact small K-Cor exposure named by this case.

This is a strict header inspector, not a density inversion or a general FITS
science reader. It consumes an independently acquired compressed exposure and
pins its original bytes before extracting the primary header's80-byte cards.
No image pixel is transformed, clipped, interpolated or treated as a density.
"""
from __future__ import annotations
import argparse
import gzip
import hashlib
import json
from pathlib import Path


def main() -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    case=Path(__file__).resolve().parents[1]
    metadata_path=case/'reference/kcor_level2_exposure.json'
    metadata_bytes=metadata_path.read_bytes()
    metadata=json.loads(metadata_bytes)
    output=args.output_dir.resolve()
    if output.exists():
        parser.error('a fresh output directory is required; evidence is never overwritten')
    compressed=args.input.read_bytes()
    digest=hashlib.sha256(compressed).hexdigest()
    if digest!=metadata['compressed_sha256'] or len(compressed)!=metadata['compressed_bytes']:
        parser.error('compressed exposure bytes differ from the reviewed product; no new version is accepted silently')
    # gzip.decompress verifies the transport CRC and trailer. The full expected
    # length additionally rejects a truncated or accidentally different product.
    raw=gzip.decompress(compressed)
    if len(raw)!=metadata['decompressed_bytes'] or len(raw)%2880:
        parser.error('invalid or unexpected FITS block length')
    cards=[]
    for offset in range(0,len(raw),80):
        card=raw[offset:offset+80].decode('ascii')
        cards.append(card)
        if card[:8].strip()=='END':
            break
    else:
        parser.error('FITS primary header has no END card')
    reproduced=('\n'.join(cards)+'\n').encode('ascii')
    expected=json.loads((case/metadata['raw_header_artifact']).read_text())
    if expected['card_width_bytes']!=80 or cards!=expected['cards'] or metadata_path.read_bytes()!=metadata_bytes:
        parser.error('original FITS header or its interpretation changed; no metadata is published')
    # The source image is intentionally not repackaged or copied to this
    # output. The exact cards preserve unit, WCS, timing and processing version
    # without manufacturing a calibrated error model or a quiet-corona label.
    output.mkdir(parents=True)
    (output/'header.txt').write_bytes(reproduced)
    receipt={'schema':'srcsep-stage13-level2-header-check-v1',
             'source_sha256':digest,'source_bytes':len(compressed),
             'header_sha256':hashlib.sha256(reproduced).hexdigest(),
             'header_cards':len(cards),'status':'PASS_METADATA_INTEGRITY_ONLY',
             'comparison_status':'NOT_RUN','quiet_baseline_verified':False,
             'electron_density_derived':False}
    temporary=output/'inspection.json.partial'
    temporary.write_text(json.dumps(receipt,indent=2)+'\n')
    temporary.replace(output/'inspection.json')
    print(json.dumps(receipt))
    return 0


if __name__=='__main__':
    raise SystemExit(main())
