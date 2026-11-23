"""Lossless cell views of Burgess (2016), DOI 10.6028/jres.121.005.

The input is the separately authenticated PMC7339710 article HTML.  This is a
lexical adapter, not an evaluation of thermochemistry: Table 3 mixes methods and
dates, and Table 13 mixes energy quantities and native units.  In particular we
never turn an unlabelled thermochemical species into an isotopic ground state.
The original cell HTML retains superscripts, links and footnotes that plain text
cannot encode uniquely.  Spanning cells are retained once; no gap is filled.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
from pathlib import Path
import re

TABLES = ("T3", "T13")
COLUMNS = ("table_id", "row_index", "cell_index", "tag", "rowspan", "colspan",
           "native_text", "native_html")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def cells(document: str, table_id: str):
    """Return literal native cells; indices describe markup, not physical axes.

    This deliberately narrow parser accepts the article's explicit closed HTML
    table markup.  It does not repair an incomplete HTML response, infer omitted
    cells from row spans, normalize whitespace or remove a source's footnotes.
    """
    if table_id not in TABLES:
        raise ValueError("unsupported article table ID")
    starts = list(re.finditer(r'<section\b[^>]*\bid=[\"\']' + table_id +
                             r'[\"\'][^>]*>', document, re.I))
    if len(starts) != 1:
        raise ValueError("missing/duplicate native table section: " + table_id)
    tail = document[starts[0].end():]
    opening = re.search(r'<table\b[^>]*>', tail, re.I)
    closing_section = re.search(r'</section\s*>', tail, re.I)
    if opening is None or (closing_section and opening.start() > closing_section.start()):
        raise ValueError("native section has no table: " + table_id)
    end = re.search(r'</table\s*>', tail[opening.end():], re.I)
    if end is None or (closing_section and opening.end() + end.end() > closing_section.start()):
        raise ValueError("unclosed native table: " + table_id)
    body = tail[opening.end():opening.end() + end.start()]
    if re.search(r'<table\b', body, re.I):
        raise ValueError("nested tables are outside the declared adapter")
    rows = re.findall(r'<tr\b[^>]*>(.*?)</tr\s*>', body, re.I | re.S)
    if not rows or len(re.findall(r'<tr\b', body, re.I)) != len(rows):
        raise ValueError("missing/malformed native rows")
    result = []
    for row_index, row in enumerate(rows):
        found = list(re.finditer(r'<(td|th)\b([^>]*)>(.*?)</\1\s*>', row, re.I | re.S))
        if not found or len(re.findall(r'<t[dh]\b', row, re.I)) != len(found):
            raise ValueError("missing/malformed native cells")
        for index, match in enumerate(found):
            tag, attrs, inner = match.groups()
            spans = []
            for name in ("rowspan", "colspan"):
                attr = re.findall(r'\b' + name + r'\s*=\s*[\"\']([^\"\']*)[\"\']', attrs, re.I)
                present = len(re.findall(r'\b' + name + r'\s*=', attrs, re.I))
                if present != len(attr) or len(attr) > 1 or (attr and (not attr[0].isdigit() or int(attr[0]) < 1)):
                    raise ValueError("invalid/duplicate native cell span")
                spans.append(attr[0] if attr else "")
            text = html.unescape(re.sub(r'<[^>]*>', '', inner))
            result.append((table_id, str(row_index), str(index), tag.lower(),
                           *spans, text, inner))
    return tuple(result)


def prepare(source, expected_sha256, output_dir):
    """Hash-check and parse both tables before creating a new output directory.

    Units, statistical confidence and state meanings remain in native cells and
    the parent article.  A successful receipt certifies byte/view integrity only.
    """
    source, output = Path(source).resolve(), Path(output_dir)
    if not re.fullmatch(r'[0-9a-f]{64}', expected_sha256):
        raise ValueError("explicit lowercase SHA-256 is required")
    if sha256(source) != expected_sha256:
        raise ValueError("native article hash mismatch")
    if output.exists():
        raise ValueError("output already exists; preserve prior evidence")
    document = source.read_text(encoding="utf-8")
    if "PMC7339710" not in document or "10.6028/jres.121.005" not in document:
        raise ValueError("wrong article or challenged/error response")
    parsed = {name: cells(document, name) for name in TABLES}
    output.mkdir(parents=True, exist_ok=False)
    receipt = {"schema": "amps-nist-hxoy-lexical-view-v1", "status": "PASS",
               "scope": "native-source-and-cell-integrity-only",
               "source": {"path": str(source), "sha256": expected_sha256},
               "physical_activation": False, "tables": []}
    for name, values in parsed.items():
        path = output / (name + "_cells.csv")
        with path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(COLUMNS)
            writer.writerows(values)
        receipt["tables"].append({"table_id": name, "path": path.name,
                                  "row_count": len(values), "columns": list(COLUMNS),
                                  "sha256": sha256(path)})
    (output / "qa_report.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--expected-sha256", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    try:
        print(json.dumps(prepare(args.source, args.expected_sha256, args.output_dir), indent=2))
        return 0
    except (OSError, UnicodeError, ValueError) as error:
        print(json.dumps({"status": "ERROR", "message": str(error), "physical_activation": False}))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
