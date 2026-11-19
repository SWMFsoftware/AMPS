#!/usr/bin/env python3
"""Lexical views of Zhao etal2021 OH(A)+H2 source-data archives, DOI730-1.

This is a source reader, not a quenching evaluator. Figure2 mixes reproduced
experimental data and zero-total-angular-momentum quantum calculations;
Figure3 contains calculated H2 marginal distributions; Figure4 contains
wave-function densities. None is a Maxwellian rate or complete joint product
law. The three elastic/inelastic/reactive outcomes must not be silently
renormalized into a two-channel destruction probability.

Duplicate native column headings in Fig2ab and literal ``--`` missing values
are scientific source content. Long-cell CSVs preserve row/column position and
every token, including headers. No units, Unicode labels or missingness are
inferred from the values. Figure4's six columns are defined only by its native
README: two angles in radians and four adiabatic wave-function densities.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import stat
import zipfile

COLUMNS = ["native_member", "row", "column", "value_token"]
MEMBERS = {
    5: {"Fig2ab.txt": ("tab", 12), "Fig2cd.txt": ("tab", 9),
        "Fig2ef.txt": ("tab", 9), "Fig2gh.txt": ("tab", 9)},
    6: {"Fig3"+letter+".txt": ("tab", 4) for letter in "abcdef"},
    7: {"Fig4"+letters+".txt": ("whitespace", 6) for letters in ("ab", "cd", "ef")},
}


def decode_table(raw, member, mode, width):
    """Keep source ordering/headers/gaps; reject rows that cannot be located.

    TSV blank fields are preserved. The fixed-width Figure4 format has six
    whitespace-delimited tokens and no header. Its native blank lines are
    retained as administrative records with column="0", value_token=""; data
    columns remain1--6. This reserves the source row identity without inventing
    six missing values, joining separated blocks, or dropping a real gap. No
    numerical conversion is used, so zero/signed/``--``/exponent tokens retain
    their exact spelling. A wholly blank source is still an ERROR.
    """
    text = raw.decode("utf-8")
    if mode == "tab":
        data = list(csv.reader(io.StringIO(text, newline=""), delimiter="\t", strict=True))
    elif mode == "whitespace":
        data = [line.split() for line in text.splitlines()]
    else:
        raise ValueError("unknown explicit native row delimiter")
    if (not data or not any(data) or
        any(len(row) != width and not (mode == "whitespace" and not row) for row in data)):
        raise ValueError("empty/ragged native quenching table; no gap filling")
    result=[]
    for i,row in enumerate(data,1):
        if not row:
            result.append([member,str(i),"0",""])
        else:
            result.extend([member,str(i),str(j),token]for j,token in enumerate(row,1))
    return result


def project(source, output_dir, *, archive_number, expected_sha256):
    """Authenticate and decode ALL fixed source members before publishing.

    Native archives are never extracted to filesystem paths. Exact member-set,
    duplicate/symlink/encryption and expanded-size checks prevent silently
    incomplete or unsafe projections. The 32MiB cap is an explicit parser
    resource guard, not a scientific sampling/cropping limit; refusal leaves
    all source bytes intact. A NEW output directory is mandatory.
    """
    source, output = Path(source), Path(output_dir)
    if output.exists():
        raise FileExistsError("output directory already exists")
    if archive_number not in MEMBERS:
        raise ValueError("only source archives5,6,7 are supported")
    raw = source.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != expected_sha256:
        raise ValueError("native source SHA256 mismatch")
    expected = MEMBERS[archive_number]
    prepared = []
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        names = archive.namelist()
        if len(names) != len(set(names)) or set(names) != set(expected) | {"README.txt"}:
            raise ValueError("native archive member set changed")
        if sum(item.file_size for item in archive.infolist()) > 32*1024*1024:
            raise ValueError("native archive exceeds32MiB parser budget; review explicitly")
        for item in archive.infolist():
            if item.flag_bits & 1 or stat.S_ISLNK(item.external_attr >> 16):
                raise ValueError("encrypted/symlink native member unsupported")
        readme = archive.read("README.txt")
        if b"definitions of the data columns" not in readme or b"radians" not in readme:
            raise ValueError("native units/column README changed")
        for member, (mode, width) in expected.items():
            content = archive.read(member)
            rows = decode_table(content, member, mode, width)
            prepared.append((member, rows, hashlib.sha256(content).hexdigest()))
    output.mkdir(parents=True)
    results=[]
    for member, rows, member_sha in prepared:
        path = output / (member+".cells.csv")
        with path.open("x", encoding="utf-8", newline="") as stream:
            writer=csv.writer(stream, lineterminator="\n")
            writer.writerow(COLUMNS); writer.writerows(rows)
        with path.open(encoding="utf-8", newline="") as stream:
            if list(csv.reader(stream, strict=True)) != [COLUMNS]+rows:
                raise ValueError("long-cell CSV changed native tokens")
        results.append(dict(member=member, member_sha256=member_sha,
                            cell_count=len(rows), output=path.name,
                            output_sha256=hashlib.sha256(path.read_bytes()).hexdigest(), columns=COLUMNS))
    receipt=dict(schema="amps-oh-h2-quenching-source-projection-v1", status="PASS",
                 scope="lexical source reproduction only", physical_activation=False,
                 source_sha256=digest, recipe_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                 archive_number=archive_number, tables=results)
    (output/"projection.json").write_text(json.dumps(receipt,indent=2)+"\n")
    return receipt


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-file",type=Path,required=True)
    parser.add_argument("--source-sha256",required=True)
    parser.add_argument("--archive-number",type=int,choices=sorted(MEMBERS),required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args(argv)
    try:
        value=project(args.source_file,args.output_dir,archive_number=args.archive_number,
                      expected_sha256=args.source_sha256)
    except FileNotFoundError as error:
        print(json.dumps(dict(status="SKIPPED",reason=str(error)))); return 77
    except (OSError,ValueError,TypeError,zipfile.BadZipFile,csv.Error) as error:
        print(json.dumps(dict(status="ERROR",reason=str(error)))); return 2
    print(json.dumps(value,indent=2)); return 0


if __name__=="__main__":
    raise SystemExit(main())
