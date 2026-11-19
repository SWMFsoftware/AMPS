"""Workbook lexical controls: administrative fixtures are never observations.

The real-source group independently reads native member XML and authenticates
retained author ZIPs. It is explicitly data-dependent so CORE --protocol-only
can omit it without mistaking administrative strings for validation evidence.
"""
import csv
import hashlib
import io
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
import zipfile
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'utility/chemical_database'))
import native_workbooks as native


def zip_bytes(parts):
    stream=io.BytesIO()
    with zipfile.ZipFile(stream,'w')as archive:
        for name,payload in parts:archive.writestr(name,payload)
    return stream.getvalue()


def workbook(cells=None,second=True,shared=True,external=False):
    """Administrative text/numbers exercise storage, not physical quantities."""
    if cells is None:
        cells=('<row r="1"><c r="A1" t="inlineStr"><is><t>NOT OBSERVATIONS</t></is></c>'
               '<c r="C1" t="s"><v>0</v></c><c r="D1" s="7"/></row>'
               '<row r="3"><c r="B3"><v>-0.0000</v></c>'
               '<c r="D3"><f t="shared" si="2">1/0</f><v>NA-CACHE</v></c>'
               '<c r="F3" t="e"><v>#DIV/0!</v></c></row>')
    sheets='<sheet name="Visible" sheetId="1" r:id="r1"/>'
    relationships='<Relationship Id="r1" Type="'+native.REL+'/worksheet" Target="worksheets/a.xml"'+(' TargetMode="External"'if external else '')+'/>'
    if second:
        sheets+='<sheet name="Hidden" state="hidden" sheetId="2" r:id="r2"/>'
        relationships+='<Relationship Id="r2" Type="'+native.REL+'/worksheet" Target="worksheets/b.xml"/>'
    if shared:relationships+='<Relationship Id="s" Type="'+native.REL+'/sharedStrings" Target="sharedStrings.xml"/>'
    parts=[('xl/workbook.xml','<workbook xmlns="'+native.MAIN+'" xmlns:r="'+native.REL+'"><sheets>'+sheets+'</sheets></workbook>'),
           ('xl/_rels/workbook.xml.rels','<Relationships xmlns="'+native.PACKAGE+'">'+relationships+'</Relationships>'),
           ('xl/worksheets/a.xml','<worksheet xmlns="'+native.MAIN+'"><sheetData>'+cells+'</sheetData></worksheet>')]
    if second:parts.append(('xl/worksheets/b.xml','<worksheet xmlns="'+native.MAIN+'"><sheetData><row r="9"><c r="Z9" t="inlineStr"><is><r><t>--</t></r><r><t xml:space="preserve"> keep </t></r></is></c></row></sheetData></worksheet>'))
    if shared:parts.append(('xl/sharedStrings.xml','<sst xmlns="'+native.MAIN+'"><si><r><t>shared</t></r><r><t xml:space="preserve"> text</t></r></si></sst>'))
    return zip_bytes(parts)


class WorkbookProtocol(unittest.TestCase):
    def setUp(self):self.temp=tempfile.TemporaryDirectory();self.work=Path(self.temp.name)
    def tearDown(self):self.temp.cleanup()
    def test_all_sheets_including_hidden_are_ordered(self):
        rows,meta=native.decode_xlsx(workbook())
        self.assertEqual([(x['name'],x['state'])for x in meta['worksheets']],[('Visible','visible'),('Hidden','hidden')])
        self.assertEqual(rows[-1][1:4],['Hidden','hidden','Z9'])
    def test_missing_sparse_cells_are_not_filled(self):
        rows,_=native.decode_xlsx(workbook(second=False));self.assertEqual([r[3]for r in rows],['A1','C1','D1','B3','D3','F3'])
    def test_explicit_blank_and_literal_missing_marker_survive(self):
        rows,_=native.decode_xlsx(workbook());self.assertEqual(rows[2][5:9],['','','','']);self.assertEqual(rows[-1][6],'-- keep ')
    def test_source_number_lexemes_are_not_floats(self):
        rows,_=native.decode_xlsx(workbook());self.assertEqual(rows[3][5],'-0.0000')
    def test_shared_strings_retain_index_and_rich_text(self):
        rows,_=native.decode_xlsx(workbook());self.assertEqual(rows[1][5:7],['0','shared text'])
    def test_formula_and_cache_are_never_evaluated(self):
        rows,_=native.decode_xlsx(workbook());self.assertEqual(rows[4][5],'NA-CACHE');self.assertEqual(rows[4][7],'1/0');self.assertEqual(rows[4][8],'{"si":"2","t":"shared"}')
    def test_error_token_and_style_remain_source_tokens(self):
        rows,_=native.decode_xlsx(workbook());self.assertEqual(rows[5][4:6],['e','#DIV/0!']);self.assertEqual(rows[2][9],'7')
    def test_phonetic_annotation_is_not_appended_to_source_label(self):
        raw=workbook('<row><c r="A1" t="inlineStr"><is><r><t>LABEL</t></r><rPh sb="0" eb="5"><t>guide</t></rPh></is></c></row>')
        rows,_=native.decode_xlsx(raw);self.assertEqual(rows[0][6],'LABEL')
    def test_corrupt_unused_member_crc_refused(self):
        raw=bytearray(zip_bytes([('unused.txt','UNIQUE-DAMAGE-TARGET')]))
        offset=raw.index(b'UNIQUE-DAMAGE-TARGET');raw[offset]^=1
        with self.assertRaisesRegex(ValueError,'CRC'):native.decode_xlsx(bytes(raw))
    def test_bad_shared_reference_refused(self):
        with self.assertRaisesRegex(ValueError,'shared string'):native.decode_xlsx(workbook(shared=False))
    def test_duplicate_coordinate_refused(self):
        with self.assertRaisesRegex(ValueError,'coordinate'):native.decode_xlsx(workbook('<row><c r="A1"/><c r="A1"/></row>'))
    def test_external_sheet_refused_without_network(self):
        with self.assertRaisesRegex(ValueError,'external worksheet'):native.decode_xlsx(workbook(external=True))
    def test_archive_traversal_refused(self):
        with self.assertRaisesRegex(ValueError,'escaping'):native.decode_xlsx(zip_bytes([('../escape','x')]))
    def test_bad_xml_and_doctype_refused(self):
        for xml in ['<bad>','<!DOCTYPE x [<!ENTITY a "x">]><x/>']:
            with self.assertRaises(ValueError):native.decode_xlsx(zip_bytes([('xl/workbook.xml',xml)]))
    def test_duplicate_archive_member_refused(self):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore');raw=zip_bytes([('xl/workbook.xml','x'),('xl/workbook.xml','y')])
        with self.assertRaisesRegex(ValueError,'duplicate'):native.decode_xlsx(raw)
    def test_unsupported_cell_type_not_silently_dropped(self):
        with self.assertRaisesRegex(ValueError,'cell type'):native.decode_xlsx(workbook('<row><c r="A1" t="unknown"/></row>'))
    def test_archive_expansion_guard_refuses_not_truncates(self):
        old=native.MAX_EXPANDED_BYTES
        try:
            native.MAX_EXPANDED_BYTES=10
            with self.assertRaisesRegex(ValueError,'resource budget'):native.decode_xlsx(workbook())
        finally:native.MAX_EXPANDED_BYTES=old
    def test_projection_roundtrip_and_no_overwrite(self):
        source=self.work/'source.zip';source.write_bytes(zip_bytes([('a.xlsx',workbook()),('b.xlsx',workbook(second=False))]));before=source.read_bytes()
        result=native.project_zip(source,self.work/'out');self.assertEqual(len(result['workbooks']),2)
        with (self.work/'out/workbook-000.cells.csv').open(newline='')as f:rows=list(csv.reader(f))
        self.assertEqual(rows[0],native.COLUMNS);self.assertEqual(rows[1][6],'NOT OBSERVATIONS')
        with self.assertRaises(FileExistsError):native.project_zip(source,self.work/'out')
        self.assertEqual(source.read_bytes(),before)
    def test_missing_source_cli_is_skipped(self):
        r=subprocess.run([sys.executable,'-B',str(ROOT/'utility/chemical_database/native_workbooks.py'),'--source-zip',str(self.work/'missing'),'--output-dir',str(self.work/'out')],capture_output=True,text=True)
        self.assertEqual(r.returncode,77,r.stderr);self.assertFalse((self.work/'out').exists())

    def test_explicit_text_columns_keep_observed_model_signs_and_markers(self):
        source=self.work/'source.txt';source.write_text('DECLARED HEADER; NOT OBSERVATIONS\n1.00 -2.0 -- 9e-3\n')
        result=native.project_text(source,self.work/'observed.csv',skip_rows=1,source_width=4,selected_columns=[0,1,2],columns=['coordinate','observed','missing'])
        self.assertEqual(result['row_count'],1)
        with(self.work/'observed.csv').open(newline='')as stream:self.assertEqual(list(csv.reader(stream))[1],['1.00','-2.0','--'])
    def test_unexpected_text_gap_refused_without_filling(self):
        source=self.work/'source.txt';source.write_text('HEADER\n1 2\n\n3 4\n')
        with self.assertRaisesRegex(ValueError,'width/gap'):
            native.project_text(source,self.work/'out',skip_rows=1,source_width=2,selected_columns=[0,1],columns=['a','b'])
        self.assertFalse((self.work/'out').exists())


class RealAuthorWorkbooks(unittest.TestCase):
    requires_authoritative_source=True
    @classmethod
    def setUpClass(cls):
        base=ROOT/'test_output/srcMoon/m3-physics-literature-20261119T072531Z'
        cls.sources=[base/'molecular-laboratory/41467_2021_22824_MOESM3_ESM.zip',base/'molecular-source-data/chang2021-excited-H2-source-data.zip']
        if any(not p.is_file()for p in cls.sources):raise unittest.SkipTest('author-deposited Chang ZIPs unavailable')
    def test_every_native_cell_coordinate_and_value_matches_xml(self):
        for source in self.sources:
            with zipfile.ZipFile(source)as outer:
                for name in outer.namelist():
                    if not name.endswith('.xlsx'):continue
                    raw=outer.read(name);rows,meta=native.decode_xlsx(raw,name)
                    with zipfile.ZipFile(io.BytesIO(raw))as archive:
                        expected=[]
                        for sheet in meta['worksheets']:
                            root=ET.fromstring(archive.read(sheet['part']))
                            for cell in root.findall('.//{'+native.MAIN+'}sheetData/{'+native.MAIN+'}row/{'+native.MAIN+'}c'):
                                value=cell.find('{'+native.MAIN+'}v')
                                expected.append((sheet['name'],cell.get('r'),value.text or ''if value is not None else ''))
                    self.assertEqual([(r[1],r[3],r[5])for r in rows],expected)
    def test_all_five_author_workbooks_are_included(self):
        observed=[]
        for source in self.sources:
            with zipfile.ZipFile(source)as archive:observed += [n for n in archive.namelist()if n.endswith('.xlsx')]
        self.assertEqual(observed,['Source Data/f1.xlsx','Source Data/f2.xlsx','Source Data/f4.xlsx','Fig4.xlsx','Fig3.xlsx'])
    def test_real_source_projection_replay_preserves_every_csv_byte(self):
        with tempfile.TemporaryDirectory()as t:
            for i,source in enumerate(self.sources):
                before=(hashlib.sha256(source.read_bytes()).hexdigest(),source.stat().st_mtime_ns)
                a=native.project_zip(source,Path(t)/('a'+str(i)));b=native.project_zip(source,Path(t)/('b'+str(i)))
                self.assertEqual(a,b)
                for record in a['workbooks']:self.assertEqual((Path(t)/('a'+str(i))/record['output']).read_bytes(),(Path(t)/('b'+str(i))/record['output']).read_bytes())
                self.assertEqual(before,(hashlib.sha256(source.read_bytes()).hexdigest(),source.stat().st_mtime_ns))

if __name__=='__main__':unittest.main()
