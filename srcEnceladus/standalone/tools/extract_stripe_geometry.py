#!/usr/bin/env python3
"""Extract the seven 2014 stripe polylines from the existing AMPS table."""
import argparse
import json
import re
from pathlib import Path


NAME_MAP = {
    "ALEXANDRIA": "alexandria",
    "BAGDADHOT": "baghdad_extra",
    "BAGDAD": "baghdad",
    "CAIROHOT": "cairo_extra",
    "CAIRO": "cairo",
    "DAMASCUSHOT": "damascus_extra",
    "DAMASCUS": "damascus",
}


def extract(text):
    pattern = re.compile(
        r"TigerStripeGeometry__(\w+)\s*=\s*\{(\d+),\{\{(.*?)\}\}\};",
        re.DOTALL)
    result = {}
    for match in pattern.finditer(text):
        raw_name, count_text, body = match.groups()
        if raw_name not in NAME_MAP:
            continue
        count = int(count_text)
        number = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
        values = [float(value) for value in re.findall(number, body)]
        if len(values) < 3*count:
            raise ValueError("Not enough coordinate values for "+raw_name)
        points = [values[index:index+3] for index in range(0, 3*count, 3)]
        if len(points) != count:
            raise ValueError("Point count mismatch for "+raw_name)
        result[NAME_MAP[raw_name]] = points
    if set(result) != set(NAME_MAP.values()):
        raise ValueError("Did not recover exactly the seven published fissures")
    if sum(map(len, result.values())) != 102:
        raise ValueError("Expected 102 published source nodes")
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source")
    parser.add_argument("output")
    args = parser.parse_args()
    source = Path(args.source)
    output = Path(args.output)
    document = {
        "schema_version": 1,
        "frame": "IAU_ENCELADUS",
        "units": "m",
        "provenance": {
            "paper": "Tenishev et al. (2014), doi:10.1002/2014JE004700",
            "source_file": str(source),
            "source_symbol_prefix": "TigerStripeGeometry__",
            "note": "102 nodes in four normal and three extra fissures",
        },
        "stripes": extract(source.read_text()),
    }
    with output.open("x") as stream:
        json.dump(document, stream, indent=2, allow_nan=False)
        stream.write("\n")


if __name__ == "__main__":
    main()
