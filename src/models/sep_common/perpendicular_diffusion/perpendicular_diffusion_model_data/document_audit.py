#!/usr/bin/env python3
"""Check equation/reference identifiers, delimiters, contents, and tables."""
import argparse
import hashlib
import json
import re
import unicodedata
from pathlib import Path

from document_fixtures import tables

ROOT = Path(__file__).resolve().parent


def slug(title):
    cleaned = "".join(ch for ch in title.lower()
                      if ch in "-_ " or unicodedata.category(ch)[0] in "LN")
    return cleaned.replace(" ","-")


def verify(path):
    text = path.read_text()
    assert "\u200b" not in text,"Unexpected zero-width character"
    assert text.count("```")%2==0,"Unclosed code fence"
    assert text.count("$$")%2==0,"Unclosed display-math delimiter"
    tags = re.findall(r"\\tag\{([^}]+)\}",text)
    assert len(tags)==len(set(tags)),"Duplicate equation tag"
    defined = set(tags)
    # PoS(ICRC2015) is a conference identifier, not an equation reference.
    mentions = set(re.findall(r"(?<!PoS)\(([A-Z]+\d+)\)",text))
    assert mentions<=defined,f"Unknown equation identifiers {mentions-defined}"
    refs = re.findall(r"^\| (R\d\d) \|",text,re.M)
    assert len(refs)==len(set(refs)),"Duplicate reference identifier"
    used_refs = set(re.findall(r"\bR\d\d\b",text))
    assert used_refs<=set(refs),"Undefined reference"
    displays = re.findall(r"\$\$(.*?)\$\$",text,re.S)
    for number,math_text in enumerate(displays,1):
        cleaned = re.sub(r"\\[{}]","",math_text)
        level = 0
        for character in cleaned:
            if character=="{":level+=1
            if character=="}":level-=1
            assert level>=0,f"Brace mismatch in display {number}"
        assert level==0,f"Brace mismatch in display {number}"
        assert re.findall(r"\\begin\{([^}]+)\}",math_text)==re.findall(r"\\end\{([^}]+)\}",math_text),"Environment mismatch"
    headings = re.findall(r"^## (.+)$",text,re.M)
    anchors = {slug(h) for h in headings}
    links = re.findall(r"\]\(#([^)]*)\)",text)
    assert set(links)<=anchors,f"Broken contents anchors {set(links)-anchors}"
    assert len(links)==20,"Expected twenty contents entries"
    assert headings[-1]=="20. Implementation roadmap for Codex in AMPS","Roadmap is not final section"
    fixture = json.loads((ROOT/"benchmark_points.json").read_text())
    for name,content in tables(fixture).items():
        block = f"<!-- BEGIN GENERATED {name} -->\n{content}\n<!-- END GENERATED {name} -->"
        assert text.count(block)==1,f"Table mismatch {name}"
    assert not re.search(r"<!-- TABLE_[A-Z_]+ -->",text),"Unfilled table marker"
    # Enforce equal columns for ordinary Markdown table blocks.
    blocks = re.findall(r"(?:^\|.*\|\s*\n){2,}",text,re.M)
    for block in blocks:
        counts = {line.count("|") for line in block.splitlines() if line.strip()}
        assert len(counts)==1,"Inconsistent Markdown table columns"
    return {"state":"passed","sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
            "lines":len(text.splitlines()),"unique_equation_tags":len(tags),
            "display_math_blocks":len(displays),"defined_references":len(refs),
            "verified_contents_anchors":len(links),"verified_generated_tables":len(tables(fixture)),
            "checked_markdown_table_blocks":len(blocks),
            "scope":"Structural and fixture-table checks; source and numerical checks are separate."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--document",type=Path,default=ROOT.parent/"PERPENDICULAR_DIFFUSION_COEFFICIENT_MODEL.md")
    parser.add_argument("--write-report",action="store_true")
    args = parser.parse_args()
    result = verify(args.document)
    if args.write_report:
        (ROOT/"document_audit.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
