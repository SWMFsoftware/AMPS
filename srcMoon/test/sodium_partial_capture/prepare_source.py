#!/usr/bin/env python3
"""Reproduce the actual D119990 partial-channel numbers from native HTML arrays.

This parser decodes the retained typed numeric arrays; it never digitizes a
curve, follows an arbitrary renderer, fits a law, interpolates or renormalizes.
The native response/hash, reaction, renderer and units fix source meaning.
"""
import argparse
import ast
import base64
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import struct

SOURCE_SHA256 = "ee669388646f4dfc54720e9d53ded92c3a853a041e5b37d2ea9c39bfed61f590"
CSV_SHA256 = "2a887ad43497a7f9e2b3d7b6173ae1b2597cc893d897926762c181946c92db93"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def walk(value):
    if isinstance(value, dict):
        yield value
        for child in value.values():
            yield from walk(child)
    elif isinstance(value, list):
        for child in value:
            yield from walk(child)


def extract_native_rows(source):
    """Decode the unique source Scatter x/y arrays with native float64 order."""
    source = Path(source)
    if sha(source) != SOURCE_SHA256:
        raise ValueError("native D119990 response hash changed; review source identity before importing another response")
    text = source.read_text()
    matches = re.findall(r"const docs_json = ('(?:\\.|[^'])*');", text)
    if len(matches) != 1 or "DataSet D119990" not in text:
        raise ValueError("native D119990 data/JSON identity missing or ambiguous")
    document = json.loads(ast.literal_eval(matches[0]))
    objects = list(walk(document))
    renderers = [row for row in objects if row.get("name") == "GlyphRenderer"]
    sources = [row for row in objects if row.get("name") == "ColumnDataSource"]
    if len(renderers) != 1 or len(sources) != 1 or renderers[0]["attributes"]["glyph"]["name"] != "Scatter":
        raise ValueError("expected one measured source Scatter; no arbitrary line/polygon selection")
    attributes = sources[0]["attributes"]["data"]
    if attributes.get("type") != "map":
        raise ValueError("unsupported native source map")
    axes = {row["attributes"]["axis_label"] for row in objects if row.get("name") in ("LinearAxis", "LogAxis")}
    if axes != {"E /eV.u-1", "sigma /cm2"}:
        raise ValueError("native energy/cross-section units changed")
    arrays = {}
    for key, value in attributes["entries"]:
        if key not in ("x", "y") or value.get("type") != "ndarray" or value.get("dtype") != "float64" or value.get("shape") != [12] or value.get("order") not in ("little", "big"):
            raise ValueError("unsupported/extra native numeric array")
        if key in arrays:
            raise ValueError("duplicate source numeric axis")
        payload = base64.b64decode(value["array"]["data"], validate=True)
        if len(payload) != 12 * 8:
            raise ValueError("native float64 byte count mismatch")
        arrays[key] = struct.unpack(("<" if value["order"] == "little" else ">") + "12d", payload)
    if set(arrays) != {"x", "y"}:
        raise ValueError("missing native numeric axis")
    rows = list(zip(arrays["x"], arrays["y"]))
    if any(not math.isfinite(x) or not math.isfinite(y) or x <= 0 or y <= 0 for x, y in rows) or any(rows[i][0] >= rows[i+1][0] for i in range(len(rows)-1)):
        raise ValueError("invalid native values/order; no repair is permitted")
    return rows


def prepare(source, output):
    source, output = Path(source).resolve(), Path(output).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError("output must be new or empty")
    rows = extract_native_rows(source)
    output.mkdir(parents=True, exist_ok=True)
    raw = output / "raw" / "collisiondb_Na_H1s_119990.html"
    raw.parent.mkdir();shutil.copy2(source, raw)
    target = output / "collisiondb_119990_embedded_arrays.csv"
    with target.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["energy_eV_per_u", "cross_section_cm2"])
        writer.writerows(rows)
    if sha(target) != CSV_SHA256:
        raise ValueError("reproduced native numeric CSV hash mismatch")
    record = {"dataset_id": 119990, "source_url": "https://db-amdis.org/collisiondb/datasets/119990/",
        "native_file": str(raw.relative_to(output)), "native_sha256": SOURCE_SHA256,
        "output": target.name, "output_sha256": CSV_SHA256, "rows": len(rows),
        "reaction": "H+ + Na -> H(1s) + Na+(2p6); measured ground-product partial capture",
        "native_initial_Na_state": "not_supplied_by_native_label",
        "approved_reference_initial_Na_state": "3s; explicit reference model, not recovered experimental metadata",
        "physical_initial_state_crosswalk_qualified": False,
        "extraction": "native typed float64 arrays; no interpolation/digitization/fitting",
        "status": "SOURCE_EXTRACTION_ONLY_NOT_FULL_KINETIC_READY", "script_sha256": sha(__file__)}
    (output / "native_inventory.json").write_text(json.dumps([record], indent=2) + "\n")
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-html", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        record = prepare(args.source_html, args.output_dir)
    except FileNotFoundError as error:
        print(json.dumps({"status": "SKIPPED", "reason": str(error)}));return 77
    except (OSError, ValueError, KeyError, TypeError) as error:
        print(json.dumps({"status": "ERROR", "reason": str(error)}));return 2
    print(json.dumps({"status": "PASS", "scope": "source extraction only", "source": record}, indent=2));return 0


if __name__ == "__main__":
    raise SystemExit(main())
