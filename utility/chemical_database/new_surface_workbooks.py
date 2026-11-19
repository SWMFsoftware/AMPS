#!/usr/bin/env python3
"""Source-pinned single-XLSX entry point for the existing lexical decoder.

The native Li2025 workbooks include sparse cells, Chinese/English labels,
pressure logs, measured water contents and source-calculated columns. This
wrapper retains all worksheets/cells and does not decide their physical role,
linearize logarithms, evaluate formulas, infer units or adopt a surface law.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import importlib.util
import json
import xml.etree.ElementTree as ET
import zipfile


def project(source, output, *, expected_sha256, decoder_sha256):
    """Bind both native bytes and maintained decoder before lexical export.

    The decoder must be the authenticated sibling ``native_workbooks.py``.
    A packaged recipe therefore keeps both maintained files, rather than
    introducing a new XLSX implementation or relying on a workstation package.
    """
    source, output=Path(source),Path(output)
    if output.exists():raise FileExistsError("CSV output already exists")
    raw=source.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=expected_sha256:raise ValueError("native source SHA256 mismatch")
    decoder=Path(__file__).with_name("native_workbooks.py")
    if hashlib.sha256(decoder.read_bytes()).hexdigest()!=decoder_sha256:
        raise ValueError("maintained decoder SHA256 mismatch")
    spec=importlib.util.spec_from_file_location("authenticated_native_workbooks",decoder)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    rows,metadata=module.decode_xlsx(raw,source.name)
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open("x",encoding="utf-8",newline="")as stream:
        writer=csv.writer(stream,lineterminator="\n");writer.writerow(module.COLUMNS);writer.writerows(rows)
    with output.open(encoding="utf-8",newline="")as stream:
        if list(csv.reader(stream))!=[module.COLUMNS]+rows:raise ValueError("CSV changed source cells")
    return dict(metadata,status="PASS",physical_activation=False,scope="lexical source access only",
                columns=module.COLUMNS,recipe_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                decoder_sha256=decoder_sha256,output_sha256=hashlib.sha256(output.read_bytes()).hexdigest())


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-file",type=Path,required=True)
    parser.add_argument("--source-sha256",required=True)
    parser.add_argument("--decoder-sha256",required=True)
    parser.add_argument("--output-file",type=Path,required=True)
    args=parser.parse_args(argv)
    try:
        value=project(args.source_file,args.output_file,expected_sha256=args.source_sha256,
                      decoder_sha256=args.decoder_sha256)
    except FileNotFoundError as error:
        print(json.dumps(dict(status="SKIPPED",reason=str(error))));return 77
    except (OSError,ValueError,TypeError,zipfile.BadZipFile,ET.ParseError)as error:
        print(json.dumps(dict(status="ERROR",reason=str(error))));return 2
    print(json.dumps(value,indent=2));return 0


if __name__=="__main__":raise SystemExit(main())
