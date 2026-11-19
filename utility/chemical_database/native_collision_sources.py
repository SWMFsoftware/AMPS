"""Strict lexical access to native IAEA CollisionDB fit-parameter pages.

A successful parse means that published coefficient *strings* were retained;
no fit is evaluated and no event distribution, source state, validity bound or
permission is inferred. In particular, the plotted sample grid must never be
substituted for author knots or a missing fit-domain limit. Source-specific
licence admission belongs to source_packets.py, separately from this reader.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
from decimal import Decimal, InvalidOperation
from html.parser import HTMLParser
from pathlib import Path


class NativeCollisionSourceError(ValueError):
    """Missing, inconsistent or corrupt native fit metadata; no fallback value."""


class _Text(HTMLParser):
    """Decode native HTML text while excluding executable/visualization scripts."""

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.depth = 0
        self.parts = []

    def handle_starttag(self, tag, attrs):
        if tag in ("script", "style", "nav"):
            self.depth += 1

    def handle_endtag(self, tag):
        if tag in ("script", "style", "nav"):
            self.depth -= 1

    def handle_data(self, data):
        if not self.depth:
            self.parts.append(data)


def _text(fragment):
    parser = _Text()
    parser.feed(fragment)
    return " ".join(" ".join(parser.parts).split())


def _field(html, name):
    match = re.search(r"<td[^>]*>\s*" + re.escape(name) +
                      r"\s*</td>\s*<td[^>]*>(.*?)</td>", html, re.S | re.I)
    return _text(match[1]) if match else None


def parse_iaea_fit_page(html, expected_dataset_id=None):
    """Return exact ordered native coefficient lexemes and explicit metadata.

    The source template places coefficients in one ``fit-coeff-table``. Its
    decimal strings are not rounded to Python float. Parameter dimensions are
    formula-specific and remain unspecified where the source omits them. Blank
    formula/domain metadata is retained as missing, rather than filled from a
    plot or from a different dataset. This routine does not execute source code.
    """
    if not isinstance(html, str):
        raise NativeCollisionSourceError("Native HTML must be decoded text")
    ids = re.findall(r"<h2[^>]*>\s*DataSet\s+(D[0-9]+)\s*</h2>", html, re.I)
    if len(ids) != 1 or (expected_dataset_id is not None and ids[0] != expected_dataset_id):
        raise NativeCollisionSourceError("Native dataset identity missing, ambiguous or mismatched")
    headings = re.findall(r"<h[1-6][^>]*>(.*?)</h[1-6]>", html, re.S | re.I)
    reactions = [_text(item) for item in headings if "\\mathrm" in item and "\\rightarrow" in item]
    if len(reactions) != 1:
        raise NativeCollisionSourceError("Exactly one explicit native reaction heading is required")
    tables = re.findall(r'<table\s+class="fit-coeff-table"[^>]*>(.*?)</table>', html, re.S | re.I)
    if len(tables) != 1:
        raise NativeCollisionSourceError("Exactly one native coefficient table is required")
    rows = []
    for row in re.findall(r"<tr[^>]*>(.*?)</tr>", tables[0], re.S | re.I):
        cells = re.findall(r"<td[^>]*>(.*?)</td>", row, re.S | re.I)
        if len(cells) != 2:
            raise NativeCollisionSourceError("Coefficient row must have exactly two native cells")
        label, value = [_text(cell) for cell in cells]
        if not label or not value:
            raise NativeCollisionSourceError("Coefficient label/value cannot be blank")
        try:
            number = Decimal(value)
        except InvalidOperation as error:
            raise NativeCollisionSourceError("Non-numerical native coefficient: " + value) from error
        if not number.is_finite():
            raise NativeCollisionSourceError("Nonfinite native coefficient: " + label)
        rows.append((label, value))
    if not rows or len({label for label, _ in rows}) != len(rows):
        raise NativeCollisionSourceError("Empty coefficient table or duplicate native labels")
    formula = re.search(r"<td[^>]*>\s*Fit Function.*?</td>\s*<td[^>]*>(.*?)</td>", html, re.S | re.I)
    formula_text = _text(formula[1]) if formula else None
    details = re.findall(r'href="(/collisiondb/fit-functions/[^"#]+)"', html)
    metadata = {key: _field(html, label) for key, label in (
        ("process", "Process"), ("columns", "Columns"), ("threshold", "Threshold"),
        ("uncertainty", "Uncertainty"), ("method", "Method"), ("comment", "Comment"),
        ("x_range", "x-range"))}
    metadata.update(reaction_as_rendered=reactions[0], formula_as_rendered=formula_text,
                    function_details_urls=["https://db-amdis.org" + value for value in details])
    missing = []
    if not formula_text or formula_text in (r"\[\]", r"\[ \]"):
        missing.append("NATIVE_FIT_FORMULA_NOT_RENDERED")
    if not metadata["x_range"]:
        missing.append("NATIVE_FIT_DOMAIN_NOT_SUPPLIED")
    return {"dataset_id": ids[0], "parameters": rows, "native_metadata": metadata,
            "scientific_status": "SOURCE_INPUT_NOT_PHYSICAL_QUALIFICATION",
            "physical_activation": False, "missing_metadata": missing}


def parameters_csv(parsed):
    """Project native strings only; the header labels are parser metadata.

    Units remain parameter-specific in the source expression. No universal
    dimensionless, eV, cm2, SI or uncertainty convention is assigned to them.
    CSV quoting is standard and row order/decimal spelling remain unchanged.
    """
    output = io.StringIO(newline="")
    writer = csv.writer(output)
    writer.writerow(("parameter_label", "value_native"))
    writer.writerows(parsed["parameters"])
    return output.getvalue().encode("utf-8")


def project(input_path, expected_sha256, output_path, dataset_id=None):
    """Authenticate original bytes before producing a new deterministic table.

    Existing outputs are never overwritten. An HTTP error document cannot pass
    the required native identity/table/reaction grammar. HTTP status/rights must
    additionally be checked by the packet loader; this parser grants neither.
    """
    path, target = Path(input_path), Path(output_path)
    raw = path.read_bytes()
    actual = hashlib.sha256(raw).hexdigest()
    if (not isinstance(expected_sha256, str) or
            not re.fullmatch(r"[0-9a-f]{64}", expected_sha256) or actual != expected_sha256):
        raise NativeCollisionSourceError("Native source SHA-256 mismatch")
    parsed = parse_iaea_fit_page(raw.decode("utf-8"), dataset_id)
    payload = parameters_csv(parsed)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("xb") as stream:
        stream.write(payload)
    return dict(parsed, source_sha256=actual, output_sha256=hashlib.sha256(payload).hexdigest(),
                row_count=len(parsed["parameters"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--expected-sha256", required=True)
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = project(args.input, args.expected_sha256, args.output, args.dataset_id)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
