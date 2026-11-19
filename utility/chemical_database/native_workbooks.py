#!/usr/bin/env python3
"""Lexical XLSX access using only ZIP/XML; never evaluate source formulae.

See NATIVE_WORKBOOK_SOURCES.md for the preimplementation contract. All values
are strings: the module does not infer dates, units, probabilities, uncertainty,
measured/model roles, missing observations or numerical normalization. Sparse
coordinates stay sparse and explicit blank cells stay explicit. The original
ZIP/XLSX remains the full formatting/metadata oracle; CSV is a detached view.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import posixpath
import re
import sys
import xml.etree.ElementTree as ET
import zipfile

MAIN='http://schemas.openxmlformats.org/spreadsheetml/2006/main'
REL='http://schemas.openxmlformats.org/officeDocument/2006/relationships'
PACKAGE='http://schemas.openxmlformats.org/package/2006/relationships'
COLUMNS=['workbook','sheet','sheet_state','cell','cell_type','value_token',
         'resolved_text','formula_token','formula_attributes_json','style_id']
# These are parser resource guards, not physical acceptance thresholds. A caller
# receiving a refusal must preserve the source and increase/review its parser
# budget explicitly; the parser never samples/crops a large scientific source.
MAX_EXPANDED_BYTES=128*1024*1024
MAX_MEMBER_BYTES=64*1024*1024
CELL=re.compile(r'[A-Z]+[1-9][0-9]*\Z')


def sha256(data):
    """Return the exact byte identity; no text/newline normalization."""
    return hashlib.sha256(data).hexdigest()


def _archive(raw):
    """Read a bounded archive without extraction, refusing ambiguous member IDs."""
    try:
        archive=zipfile.ZipFile(io.BytesIO(raw))
    except (zipfile.BadZipFile, OSError) as error:
        raise ValueError('invalid source ZIP/XLSX') from error
    seen=set();total=0
    for info in archive.infolist():
        name=info.filename
        if (name in seen or '\\' in name or name.startswith('/') or
                '..' in PurePosixPath(name).parts or '\x00' in name):
            archive.close();raise ValueError('duplicate or escaping archive member: '+name)
        seen.add(name);total+=info.file_size
        if info.flag_bits&1 or info.file_size>MAX_MEMBER_BYTES or total>MAX_EXPANDED_BYTES:
            archive.close();raise ValueError('encrypted/oversized source archive; review parser resource budget')
    try:
        if archive.testzip() is not None:raise ValueError('source ZIP member CRC failure')
    except zipfile.BadZipFile as error:
        archive.close();raise ValueError('source ZIP member CRC failure') from error
    return archive


def _xml(archive,name):
    """Parse one native XML part; no DTD/entity expansion or external retrieval."""
    try:
        raw=archive.read(name)
    except KeyError as error:
        raise ValueError('required XLSX part missing: '+name) from error
    if b'<!DOCTYPE' in raw.upper() or b'<!ENTITY' in raw.upper():
        raise ValueError('DTD/entities forbidden in workbook XML')
    try:return ET.fromstring(raw)
    except ET.ParseError as error:raise ValueError('malformed workbook XML: '+name) from error


def _target(base,target):
    """Resolve an OPC relationship, without following URLs or escaping package."""
    if not target or '\\' in target or '://' in target or '?' in target or '#' in target:
        raise ValueError('invalid/external workbook relationship target')
    resolved=posixpath.normpath(target.lstrip('/') if target.startswith('/') else
                               posixpath.join(posixpath.dirname(base),target))
    if resolved.startswith('../') or resolved in ('..','.'):
        raise ValueError('workbook relationship escapes package')
    return resolved


def _relations(archive,part):
    """Return exact relationship IDs; external relationships are never fetched."""
    name=posixpath.join(posixpath.dirname(part),'_rels',posixpath.basename(part)+'.rels')
    root=_xml(archive,name);result={}
    if root.tag!='{'+PACKAGE+'}Relationships':raise ValueError('unsupported relationship namespace')
    for row in root:
        if row.tag!='{'+PACKAGE+'}Relationship':raise ValueError('unexpected relationship element')
        identifier=row.get('Id')
        if not identifier or identifier in result:raise ValueError('duplicate/missing relationship ID')
        result[identifier]=dict(row.attrib)
    return result


def _display_text(item):
    """Resolve source display runs, excluding phonetic annotation (rPh) text.

    The original XML retains annotations and formatting. Joining every nested
    <t> would incorrectly append phonetic guides to the scientific cell label.
    Only direct <t> and rich-text <r>/<t> form the displayed source string.
    """
    if item is None:return ''
    pieces=[]
    for child in item:
        if child.tag=='{'+MAIN+'}t':pieces.append(child.text or '')
        elif child.tag=='{'+MAIN+'}r':pieces.extend(t.text or ''for t in child if t.tag=='{'+MAIN+'}t')
    return ''.join(pieces)


def decode_xlsx(raw,workbook='workbook.xlsx'):
    """Return ordered long-cell rows plus every worksheet's source metadata.

    ``value_token`` is the literal <v> text, including cached formula values or
    shared-string indices. ``resolved_text`` resolves only source text strings;
    it never converts number/date/boolean/error tokens. Formula text/attributes
    are retained separately and never evaluated. Shared-formula followers with
    empty formula text remain empty: no hidden reconstruction is performed.
    """
    rows=[];sheets=[]
    with _archive(raw) as archive:
        root=_xml(archive,'xl/workbook.xml')
        if root.tag!='{'+MAIN+'}workbook':raise ValueError('unsupported workbook namespace')
        relationships=_relations(archive,'xl/workbook.xml')
        shared=[]
        for item in relationships.values():
            if item.get('Type')==REL+'/sharedStrings':
                if item.get('TargetMode')=='External':raise ValueError('external shared strings')
                strings=_xml(archive,_target('xl/workbook.xml',item.get('Target')))
                if strings.tag!='{'+MAIN+'}sst':raise ValueError('unsupported shared string namespace')
                shared=[_display_text(item)for item in strings]
        sheet_root=root.find('{'+MAIN+'}sheets')
        if sheet_root is None:raise ValueError('workbook has no sheets declaration')
        names=set()
        for declaration in sheet_root:
            name=declaration.get('name');identifier=declaration.get('{'+REL+'}id')
            if not name or name in names or identifier not in relationships:
                raise ValueError('duplicate/missing sheet name or relationship')
            names.add(name);relationship=relationships[identifier]
            if relationship.get('TargetMode')=='External':raise ValueError('external worksheet relationship')
            kind=relationship.get('Type','')
            if kind!=REL+'/worksheet':raise ValueError('unsupported non-worksheet sheet; no silent omission')
            part=_target('xl/workbook.xml',relationship.get('Target'))
            sheet=_xml(archive,part)
            if sheet.tag!='{'+MAIN+'}worksheet':raise ValueError('unsupported worksheet namespace')
            state=declaration.get('state','visible');seen=set();before=len(rows)
            data=sheet.find('{'+MAIN+'}sheetData')
            if data is None:raise ValueError('worksheet lacks sheetData')
            for source_row in data:
                if source_row.tag!='{'+MAIN+'}row':raise ValueError('unexpected sheetData element')
                for cell in source_row:
                    if cell.tag!='{'+MAIN+'}c':raise ValueError('unexpected worksheet cell element')
                    coordinate=cell.get('r')
                    if not coordinate or not CELL.fullmatch(coordinate) or coordinate in seen:
                        raise ValueError('duplicate/invalid source cell coordinate')
                    seen.add(coordinate)
                    cell_type=cell.get('t','n');style=cell.get('s','')
                    if cell_type not in ('n','s','inlineStr','str','b','e','d'):
                        raise ValueError('unsupported cell type: '+cell_type)
                    values=cell.findall('{'+MAIN+'}v');formulas=cell.findall('{'+MAIN+'}f')
                    if len(values)>1 or len(formulas)>1:raise ValueError('duplicate cell value/formula')
                    token=(values[0].text or '')if values else ''
                    resolved=''
                    if cell_type=='s':
                        if not re.fullmatch(r'[0-9]+',token) or int(token)>=len(shared):
                            raise ValueError('invalid shared string index')
                        resolved=shared[int(token)]
                    elif cell_type=='inlineStr':
                        inline=cell.find('{'+MAIN+'}is')
                        resolved=_display_text(inline)
                    elif cell_type=='str':resolved=token
                    formula=(formulas[0].text or '')if formulas else ''
                    attributes=json.dumps(formulas[0].attrib,sort_keys=True,separators=(',',':'))if formulas else ''
                    rows.append([workbook,name,state,coordinate,cell_type,token,resolved,formula,attributes,style])
            sheets.append({'name':name,'state':state,'part':part,'cell_count':len(rows)-before,
                           'dimension':sheet.find('{'+MAIN+'}dimension').get('ref')if sheet.find('{'+MAIN+'}dimension')is not None else None})
    return rows,{'workbook':workbook,'sha256':sha256(raw),'worksheets':sheets,'cell_count':len(rows)}


def project_zip(source_zip,output_dir):
    """Export ALL XLSX members as lexical CSVs, without overwriting any output.

    Files are emitted only after all members decode successfully. This does not
    promise transactionality against disk exhaustion; partial write evidence
    remains on IO failure and a new output directory is required for any retry.
    """
    source=Path(source_zip);output=Path(output_dir)
    if output.exists():raise FileExistsError('output directory already exists')
    raw=source.read_bytes();prepared=[]
    with _archive(raw)as archive:
        for name in archive.namelist():
            if name.lower().endswith('.xlsx'):
                rows,metadata=decode_xlsx(archive.read(name),name)
                prepared.append((rows,metadata))
    if not prepared:raise ValueError('source ZIP contains no XLSX workbooks')
    output.mkdir(parents=True)
    result=[]
    for ordinal,(rows,metadata)in enumerate(prepared):
        filename='workbook-%03d.cells.csv'%ordinal
        with (output/filename).open('x',encoding='utf-8',newline='')as stream:
            writer=csv.writer(stream,lineterminator='\n');writer.writerow(COLUMNS);writer.writerows(rows)
        with (output/filename).open(encoding='utf-8',newline='')as stream:
            if list(csv.reader(stream))!=[COLUMNS]+rows:raise ValueError('CSV lexical round-trip changed source cells')
        result.append(dict(metadata,output=filename,output_sha256=sha256((output/filename).read_bytes()),columns=COLUMNS))
    receipt={'schema':'amps-native-workbook-projection-v1','status':'PASS','scope':'source lexical access only; not physical qualification',
             'source_sha256':sha256(raw),'recipe_sha256':sha256(Path(__file__).read_bytes()),'workbooks':result,
             'physical_activation':False,'formula_evaluation':False,'sparse_cells_filled':False}
    (output/'projection.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt


def project_text(source, output, *, skip_rows, source_width, selected_columns, columns):
    """Project an explicitly declared whitespace-token table, without conversion.

    The caller must bind the native header/units/column-role contract and source
    identity. This helper only enforces row width and exact token preservation;
    it does not infer a header, skip unexpected gaps, interpret source zeros as
    detections, manufacture uncertainties or evaluate a model column. Separate
    calls can expose observed and model columns without mixing their roles.
    """
    source,output=Path(source),Path(output)
    if output.exists():raise FileExistsError('output already exists')
    if (type(skip_rows)is not int or skip_rows<0 or type(source_width)is not int or source_width<1 or
        not selected_columns or len(columns)!=len(selected_columns) or len(set(columns))!=len(columns) or
        any(type(i)is not int or i<0 or i>=source_width for i in selected_columns)):
        raise ValueError('explicit source token schema is invalid')
    lines=source.read_text(encoding='utf-8').splitlines()
    if skip_rows>=len(lines):raise ValueError('source has no declared data rows')
    rows=[]
    for line in lines[skip_rows:]:
        tokens=line.split()
        if len(tokens)!=source_width:raise ValueError('native text row width/gap differs from declared schema')
        rows.append([tokens[i]for i in selected_columns])
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x',encoding='utf-8',newline='')as stream:
        writer=csv.writer(stream,lineterminator='\n');writer.writerow(columns);writer.writerows(rows)
    with output.open(encoding='utf-8',newline='')as stream:
        if list(csv.reader(stream))!=[columns]+rows:raise ValueError('CSV text token round-trip changed')
    return {'scope':'source lexical access only','source_sha256':sha256(source.read_bytes()),
            'output_sha256':sha256(output.read_bytes()),'row_count':len(rows),'columns':columns,
            'skip_rows':skip_rows,'source_width':source_width,'selected_columns':list(selected_columns),
            'native_header':lines[:skip_rows],'physical_activation':False}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    source=parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--source-zip',type=Path)
    source.add_argument('--source-text',type=Path)
    parser.add_argument('--output-dir',type=Path)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--skip-rows',type=int)
    parser.add_argument('--source-width',type=int)
    parser.add_argument('--select',help='explicit zero-based native column indices, comma separated')
    parser.add_argument('--columns',help='explicit CSV labels, comma separated')
    args=parser.parse_args()
    try:
        if args.source_zip:
            if args.output_dir is None:parser.error('--source-zip requires --output-dir')
            result=project_zip(args.source_zip,args.output_dir)
        else:
            if args.output is None or args.select is None or args.columns is None:
                parser.error('--source-text requires --output/--select/--columns')
            result=project_text(args.source_text,args.output,skip_rows=args.skip_rows,
                                source_width=args.source_width,
                                selected_columns=[int(x)for x in args.select.split(',')],
                                columns=args.columns.split(','))
    except FileNotFoundError as error:print('SKIPPED: '+str(error),file=sys.stderr);return 77
    except (ValueError,OSError,zipfile.BadZipFile)as error:print('ERROR: '+str(error),file=sys.stderr);return 2
    print(json.dumps(result,indent=2));return 0


if __name__=='__main__':raise SystemExit(main())
