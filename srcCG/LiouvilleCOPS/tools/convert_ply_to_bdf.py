#!/usr/bin/env python3
"""Convert a binary triangular PLY to the solver's metre-unit NASTRAN surface.

No welding, smoothing, triangulation, winding repair, or frame transformation
is performed. Source units and their evidence must be explicitly supplied:
PLY itself does not standardize coordinate units. Scalar vertex attributes
other than xyz are ignored; unsupported element layouts are refused.
"""
import argparse
import json
from pathlib import Path
import tempfile

import numpy as np

from convert_dsk_to_bdf import sha256, nastran_real, continuation, write_triangle
from gravity_topology import require_outward_closed_body

TYPES = {"char": "i1", "uchar": "u1", "short": "i2", "ushort": "u2",
         "int": "i4", "uint": "u4", "float": "f4", "double": "f8",
         "int8": "i1", "uint8": "u1", "int16": "i2", "uint16": "u2",
         "int32": "i4", "uint32": "u4", "float32": "f4", "float64": "f8"}


def read_ply(path):
    """Memory-map fixed triangular records; reject truncation and extra data.

    The list count is still checked for EVERY face: assuming three entries
    based only on the first face could silently misalign a polygonal file.
    """
    with Path(path).open("rb") as stream:
        lines = []
        while True:
            line = stream.readline(4097)
            if not line or len(line) > 4096 or len(lines) > 4096:
                raise ValueError("missing/oversized PLY header")
            text = line.decode("ascii").strip()
            lines.append(text)
            if text == "end_header":
                break
        offset = stream.tell()
    if lines[0] != "ply":
        raise ValueError("not a PLY file")
    formats = {"format binary_big_endian 1.0": ">",
               "format binary_little_endian 1.0": "<"}
    if lines[1] not in formats:
        raise ValueError("only binary PLY 1.0 is supported")
    endian = formats[lines[1]]
    elements = []
    for line in lines[2:-1]:
        words = line.split()
        if not words or words[0] in ("comment", "obj_info"):
            continue
        if words[0] == "element" and len(words) == 3:
            elements.append([words[1], int(words[2]), []])
        elif words[0] == "property" and elements:
            elements[-1][2].append(words[1:])
        else:
            raise ValueError("unsupported PLY header: " + line)
    if [e[0] for e in elements] != ["vertex", "face"]:
        raise ValueError("expected vertex then face elements only")
    (_, nv, vp), (_, nf, fp) = elements
    if nv < 4 or nf < 4:
        raise ValueError("surface too small")
    if any(len(p) != 2 or p[0] not in TYPES for p in vp):
        raise ValueError("vertex properties must be supported scalars")
    vd = np.dtype([(p[1], endian + TYPES[p[0]]) for p in vp])
    if not all(name in vd.names for name in ("x", "y", "z")):
        raise ValueError("missing xyz")
    if (len(fp) != 1 or len(fp[0]) != 4 or fp[0][0] != "list"
            or fp[0][3] not in ("vertex_index", "vertex_indices")):
        raise ValueError("faces must contain only a vertex-index list")
    _, count_type, index_type, _ = fp[0]
    count_dtype = np.dtype(endian + TYPES[count_type])
    index_dtype = np.dtype(endian + TYPES[index_type])
    if count_dtype.kind not in "iu" or index_dtype.kind not in "iu":
        raise ValueError("face counts and indices must be integers")
    fd = np.dtype([("count", count_dtype), ("indices", index_dtype, (3,))])
    expected = offset + nv * vd.itemsize + nf * fd.itemsize
    if Path(path).stat().st_size != expected:
        raise ValueError("PLY payload size differs from fixed triangle layout")
    vertices = np.memmap(path, mode="r", dtype=vd, offset=offset, shape=(nv,))
    faces = np.memmap(path, mode="r", dtype=fd, offset=offset + nv * vd.itemsize,
                      shape=(nf,))
    if np.any(faces["count"] != 3):
        raise ValueError("non-triangular face")
    xyz = np.column_stack([vertices[k].astype(np.float64) for k in ("x", "y", "z")])
    indices = faces["indices"]
    if not np.isfinite(xyz).all():
        raise ValueError("nonfinite coordinate")
    if np.any(indices < 0) or np.any(indices >= nv):
        raise ValueError("face index outside vertex table")
    if np.any((indices[:, 0] == indices[:, 1]) |
              (indices[:, 0] == indices[:, 2]) |
              (indices[:, 1] == indices[:, 2])):
        raise ValueError("repeated vertex in face")
    return xyz, indices, lines


def convert(source, output, units, frame, evidence, allow_failed_gravity_audit=False):
    """Publish only complete files; explicitly distinguish conversion and audit.

    A surface may be readable by the forward solver but fail the stricter
    genus-zero gravity admission. Such a product requires an explicit override
    and is marked NOT gravity-admitted; no geometry is repaired automatically.
    """
    source, output = Path(source), Path(output)
    manifest = output.with_suffix(output.suffix + ".conversion.json")
    if output.exists() or manifest.exists():
        raise ValueError("refusing existing output or manifest")
    xyz, faces, header = read_ply(source)
    xyz *= {"m": 1., "km": 1000.}[units]
    print("Read %d vertices and %d triangles" % (len(xyz), len(faces)), flush=True)
    try:
        audit = {"passed": True, "result": require_outward_closed_body(xyz, faces)}
    except ValueError as exc:
        if not allow_failed_gravity_audit:
            raise
        audit = {"passed": False, "reason": str(exc),
                 "warning": "NOT admitted for the current gravity builder"}
    print("Gravity topology audit: %s" % audit, flush=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="ascii",
                                         dir=str(output.parent), delete=False,
                                         prefix=output.name + ".", suffix=".tmp") as stream:
            temporary = Path(stream.name)
            stream.write("$ PLY conversion: coordinates in metres; original vertex/face order.\n")
            stream.write("$ No frame transform, winding repair, welding or decimation.\n")
            if not audit["passed"]:
                stream.write("$ WARNING: gravity topology admission FAILED; see conversion JSON.\n")
            for node, point in enumerate(xyz, 1):
                x, y, z = [nastran_real(v) for v in point]
                cont = continuation(node)
                stream.write("%-8s%16d%16d%16s%16s%-8s\n" %
                             ("GRID*", node, 0, x, y, cont))
                stream.write("%-8s%16s\n" % (cont, z))
            for element, nodes in enumerate(faces, 1):
                write_triangle(stream, element, tuple(int(i) + 1 for i in nodes))
            stream.write("ENDDATA\n")
        report = {"source": str(source.resolve()), "source_sha256": sha256(source),
                  "source_header": header, "source_units": units, "output_units": "m",
                  "frame": frame, "units_and_frame_evidence": evidence,
                  "frame_transform": "none", "vertices": len(xyz), "faces": len(faces),
                  "bounds_m": [xyz.min(axis=0).tolist(), xyz.max(axis=0).tolist()],
                  "geometry_changes": "none except unit scaling and decimal formatting",
                  "gravity_topology_audit": audit, "output_sha256": sha256(temporary),
                  "coordinate_format": "GRID* %.9E, 10 significant decimal digits",
                  "converter_sha256": sha256(__file__)}
        temporary.rename(output)
        temporary = None
        manifest.write_text(json.dumps(report, indent=2) + "\n")
        print("Wrote %s and %s" % (output, manifest), flush=True)
        return report
    finally:
        if temporary is not None:
            temporary.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--input-units", required=True, choices=("m", "km"))
    parser.add_argument("--frame", required=True)
    parser.add_argument("--units-and-frame-evidence", required=True)
    parser.add_argument("--allow-failed-gravity-audit", action="store_true",
                        help="preserve geometry but mark output NOT admitted for gravity")
    args = parser.parse_args()
    convert(args.source, args.output, args.input_units, args.frame,
            args.units_and_frame_evidence, args.allow_failed_gravity_audit)


if __name__ == "__main__":
    main()
