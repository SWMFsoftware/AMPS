#!/usr/bin/env python3
"""Reproduce Tellez2012 TableI from the exact acquired publisher PDF.

The extractor reads native PDF character text; no OCR, plot digitization,
manual science-value repair or fitted parameter is generated. Two distinct
pdftotext layouts must agree. The known font/CMap glyph mappings are scoped
strictly to this hash-verified source's table and documented in its manifest.
"""
from __future__ import annotations
import argparse
import csv
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess

PDF_SHA256='bf554c667bb01d03952d151f86e4e348b706f20630c8fe896eedf15f6fa30ed2'
SOURCE_URL='https://pure.tue.nl/ws/portalfiles/portal/3503988/583864115491353.pdf'
COLUMNS=['fit_id','target','mechanism','energy_min_eV','energy_min_inclusive',
 'energy_max_eV','energy_max_inclusive','published_energy_inequality',
 'intercept','intercept_reported_uncertainty','characteristic_velocity_m_s',
 'vc_reported_uncertainty_m_s','r_squared','O_Si_ratio','O_Si_ratio_reported_uncertainty',
 'source_vc_in_1e5_m_s','source_vc_uncertainty_in_1e5_m_s']

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def parse_native_table(text):
    """Decode only the frozen native TableI, retaining missing cells as blanks.

    In this PDF's character map U+0006 denotes <= in the Si energy cell,
    while the parenthesized uncertainty glyph is extracted as ASCII '6'.
    These mappings were checked against the rendered native table. They are
    not generic substitutions or guesses about arbitrary scientific PDFs.
    """
    if 'TABLE I.' not in text or 'Ei range (keV)' not in text:
        raise ValueError('expected source TableI/energy units missing')
    table=text[text.index('Mechanism',text.index('TABLE I.')):]
    if not re.search(r'Vc\s*\(\x02105\s*m/s\)',table):
        raise ValueError('native TableI velocity scale is not the verified1e5m/s')
    table=table.split('Downloaded')[0]
    table=re.sub(r'\b(16|18)\s+AN\s+O\b',lambda m:'AN O'+m[1],table)
    table=re.sub(r'\b(16|18)\s+O\b',lambda m:'O'+m[1],table)
    table=table.replace('\x06','<=').replace('\x08\x08\x08','MISSING')
    table=re.sub(r'O\s*\(mean\)','O_mean',table)
    number=r'([0-9]+\.[0-9]+)'
    pair=number+r'\s*\((?:6|±)\s*'+number+r'\)'
    pattern=re.compile(r'(?:(CIN|AN)\s+)?(Si|O16|O18|O_mean)\s+Ei\s*(<=|>|<)\s*'+number+r'\s+'+pair+r'\s+'+pair+r'\s+('+number+r'|MISSING)(?:\s+('+number+r'|MISSING)(?:\s*\((?:6|±)\s*'+number+r'\))?)?')
    # Keep row-boundary parsing unambiguous: process one physical row after
    # joining the native superscript/isotope line, not a whitespace-flattened
    # table that could borrow the following row's numbers.
    result=[];mechanism=None
    for line in table.splitlines():
        if not re.search(r'\bEi\s*(?:<=|>|<)',line):continue
        match=pattern.fullmatch(line.strip())
        if not match:raise ValueError('unexpected native TableI row: '+repr(line))
        g=match.groups()
        # Repeated captures in the optional columns are named by position;
        # explicit validation below prevents shifting a missing cell into a
        # neighbouring scientific quantity.
        supplied,target,relation,limit,a,da,vc,dvc,r2,_r2dup,ratio,_ratiodup,dratio=g
        mechanism=supplied or mechanism
        if mechanism not in ('CIN','AN'):raise ValueError('mechanism not declared before row')
        if target=='Si':
            if mechanism!='CIN' or relation!='<=' or limit!='7.0':raise ValueError('Si source range changed')
            lower,lo,upper,hi='600','true','7000','true'
        elif mechanism=='CIN':
            if relation!='>' or limit!='2.0':raise ValueError('oxygen CIN source range changed')
            lower,lo,upper,hi='2000','false','7000','true'
        else:
            if relation!='<' or limit!='0.8':raise ValueError('oxygen AN source range changed')
            lower,lo,upper,hi='600','true','800','false'
        clean=lambda value:'' if value in (None,'MISSING') else value
        result.append(dict(zip(COLUMNS,[target+'_'+mechanism,target,mechanism,
            lower,lo,upper,hi,'Ei'+relation+limit+'keV',a,da,
            format(Decimal(vc)*100000,'.0f'),format(Decimal(dvc)*100000,'.0f'),
            clean(r2),clean(ratio),clean(dratio),vc,dvc])))
    expected=['Si_CIN','O16_CIN','O18_CIN','O_mean_CIN','O16_AN','O18_AN','O_mean_AN']
    if [row['fit_id'] for row in result]!=expected:raise ValueError('complete ordered seven-row family missing')
    return result

def prepare(pdf,output):
    pdf=Path(pdf).resolve();output=Path(output).resolve()
    if sha(pdf)!=PDF_SHA256:raise ValueError('publisher PDF source SHA256 mismatch')
    executable=shutil.which('pdftotext')
    if not executable:raise FileNotFoundError('pdftotext dependency unavailable')
    output.mkdir(parents=True,exist_ok=False);parsed=[];commands=[]
    for mode in ('raw','layout'):
        target=output/('table_page4_'+mode+'.txt')
        command=[executable,'-f','4','-l','4','-'+mode,str(pdf),str(target)]
        subprocess.run(command,check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=30)
        commands.append(command);parsed.append(parse_native_table(target.read_text()))
    if parsed[0]!=parsed[1]:raise ValueError('independent native raw/layout TableI extractions differ')
    table=output/'table_i.csv'
    with table.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=COLUMNS,lineterminator='\n');writer.writeheader();writer.writerows(parsed[0])
    version=subprocess.run([executable,'-v'],capture_output=True,text=True,timeout=15)
    manifest={'schema':'srcMoon-Tellez2012-source-table-v1','scope':'published conditional charge-survival fits only; not observations or lunar retention',
      'source_pdf':{'path':str(pdf),'sha256':PDF_SHA256,'url':SOURCE_URL,'doi':'10.1063/1.4758699'},
      'source_locator':'publisher article151602-3/TableI, PDFpage4; Eq3 on151602-2/PDFpage3',
      'commands':commands,'pdftotext_version_stdout_stderr':version.stdout+version.stderr,
      'processor':{'path':str(Path(__file__).resolve()),'sha256':sha(__file__)},
      'output_table':{'path':str(table),'sha256':sha(table),'rows':len(parsed[0])},
      'operations':['native PDF text extraction in raw and layout order','source-version-specific verified character-map normalization','TableI Vc and reported uncertainty multiply exactly1e5','source row inequalities intersect actual600..7000eVstudy support'],
      'uncertainty':'published +/- only; confidence/covariance not supplied, no confidence label assigned',
      'missing_cell_policy':'blank retained; no value synthesized','interpolation':'none','canonical_data_modified':False,
      'physical_velocity_mapping':'UNAVAILABLE: exact TableI fitted velocity-coordinate convention not established',
      'formal_or_observational_validation':'NOT EVALUATED'}
    (output/'source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (output/'qa.json').write_text(json.dumps({'status':'PASS','scope':'source identity, units, complete rows and two native text layouts only','row_count':7,'raw_layout_identical':True,'source_sha256':PDF_SHA256,'table_sha256':sha(table),'synthetic_observations':False},indent=2)+'\n')
    return manifest

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source-pdf',type=Path,required=True);parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args();print(json.dumps(prepare(args.source_pdf,args.output_dir),indent=2));return 0
if __name__=='__main__':raise SystemExit(main())
