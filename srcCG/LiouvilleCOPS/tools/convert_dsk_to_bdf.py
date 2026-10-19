#!/usr/bin/env python3
"""Convert single-segment SPICE type-2 DSK models to NASTRAN BDF.

LiouvilleCOPS accepts DSK directly, but the historical SHAP5 input is a
NASTRAN bulk-data file containing ``GRID*`` and ``CTRIA3`` cards.  This tool
uses the NAIF N0067 utilities as the authoritative DSK reader, preserves the
relative order of retained DSK plates, changes vertex units from kilometres
to metres, and writes the surface-only subset of the historical BDF convention.

The conversion is deliberately strict.  ``dskbrief`` must report exactly one
type-2 segment in the requested frame, the counts exported by ``dskexp`` must
equal the DSK metadata, and every face must contain three distinct in-range
node indices. Reversed exact duplicate fin pairs are removed; every other
topology defect is rejected. Existing BDF/manifest paths are refused. A temporary BDF is
atomically renamed only after conversion and hashing succeed.
"""

import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np


SCHEMA = "liouville_cops.dsk_to_nastran"
SCHEMA_VERSION = 1
KM_TO_M = 1000.0


def sha256(path):
    """Return a streaming SHA-256 identity without loading a large mesh."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_naif_tool(name, explicit, cspice_root):
    """Resolve a NAIF executable without embedding a machine-specific path."""
    if explicit:
        candidate = Path(explicit).expanduser().resolve()
    elif cspice_root:
        candidate = Path(cspice_root).expanduser().resolve() / "exe" / name
    else:
        found = shutil.which(name)
        if not found:
            raise ValueError(
                "%s is not on PATH; pass --%s or --cspice-root" %
                (name, name))
        candidate = Path(found).resolve()
    if not candidate.is_file() or not os.access(str(candidate), os.X_OK):
        raise ValueError("NAIF executable is missing or not executable: %s" % candidate)
    return candidate


def inspect_dsk(dskbrief, frame_kernel, dsk, expected_frame):
    """Return strict segment metadata parsed from NAIF dskbrief output."""
    command = [str(dskbrief), "-full", str(frame_kernel), str(dsk)]
    completed = subprocess.run(
        command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, check=False)
    if completed.returncode:
        raise ValueError("dskbrief failed for %s:\n%s" % (dsk, completed.stdout))
    text = completed.stdout
    segments = re.findall(r"^Segment number\s+\d+\s+of file", text, re.M)
    if len(segments) != 1:
        raise ValueError("%s must contain exactly one DSK segment (found %d)" %
                         (dsk, len(segments)))
    data_type = re.search(r"^\s*Data type:\s+(\d+)\b", text, re.M)
    frame = re.search(r"^\s*Reference frame:\s+(.+?)\s*$", text, re.M)
    vertices = re.search(r"^\s*Number of vertices:\s+(\d+)\b", text, re.M)
    plates = re.search(r"^\s*Number of plates:\s+(\d+)\b", text, re.M)
    if not all((data_type, frame, vertices, plates)):
        raise ValueError("cannot parse type/frame/counts from dskbrief output for %s" % dsk)
    if int(data_type.group(1)) != 2:
        raise ValueError("%s has DSK data type %s, not type 2" %
                         (dsk, data_type.group(1)))
    actual_frame = frame.group(1).strip()
    if actual_frame != expected_frame:
        raise ValueError("%s uses frame %s, expected %s" %
                         (dsk, actual_frame, expected_frame))
    return {
        "segments": 1,
        "data_type": 2,
        "frame": actual_frame,
        "vertices": int(vertices.group(1)),
        "plates": int(plates.group(1)),
        "dskbrief_command": command,
    }


def nastran_real(value):
    """Format a metre coordinate in a 16-character large-field real."""
    result = "%.9E" % value
    if len(result) > 16:
        raise ValueError("coordinate cannot fit a NASTRAN 16-character field: %s" % result)
    return result


def continuation(node_id):
    """Make an eight-character large-field continuation identifier."""
    if node_id < 1 or node_id > 0xFFFFFFF:
        raise ValueError("node ID is outside the supported continuation range")
    return "*%07X" % node_id


def write_grid(stream, node_id, xyz_km):
    """Write one GRID* card, converting the DSK coordinates to metres."""
    x, y, z = (nastran_real(value * KM_TO_M) for value in xyz_km)
    cont = continuation(node_id)
    stream.write("%-8s%16d%16d%16s%16s%-8s\n" %
                 ("GRID*", node_id, 0, x, y, cont))
    stream.write("%-8s%16s\n" % (cont, z))


def write_triangle(stream, element_id, nodes):
    """Write one fixed-small-field CTRIA3 card with property ID 1."""
    if element_id > 99999999 or any(value > 99999999 for value in nodes):
        raise ValueError("CTRIA3 or node ID does not fit an 8-character field")
    stream.write("%-8s%8d%8d%8d%8d%8d\n" %
                 (("CTRIA3", element_id, 1) + tuple(nodes)))


def permutation_sign(face):
    """Return the orientation parity of a triangle relative to sorted nodes."""
    inversions = sum(face[i] > face[j] for i in range(3) for j in range(i + 1, 3))
    return 1 if inversions % 2 == 0 else -1


def sanitize_faces(faces):
    """Remove reversed duplicate fins and prove closed, consistent topology.

    Four archived decimations contain a small number of zero-thickness fins:
    the same triangle occurs twice with opposite winding and is attached to
    the exterior along one edge. Such a pair contributes no closed volume but
    makes that attachment edge non-manifold. Removing both exact copies is
    deterministic; any same-winding duplicate or other topology defect is
    rejected rather than guessed away.
    """
    canonical = np.sort(faces, axis=1)
    key_type = np.dtype([("a", "<i4"), ("b", "<i4"), ("c", "<i4")])
    keys = np.ascontiguousarray(canonical).view(key_type).reshape(-1)
    order = np.argsort(keys, kind="mergesort")
    ordered = keys[order]
    same = np.zeros(len(ordered), dtype=bool)
    if len(ordered) > 1:
        same[1:] = ordered[1:] == ordered[:-1]
    duplicate_positions = np.flatnonzero(same)
    remove = np.zeros(len(faces), dtype=bool)
    removed_pairs = []
    for position in duplicate_positions:
        if position + 1 < len(ordered) and ordered[position + 1] == ordered[position]:
            raise ValueError("a triangle occurs more than twice in the DSK export")
        first, second = int(order[position - 1]), int(order[position])
        if remove[first] or remove[second]:
            raise ValueError("overlapping duplicate-triangle groups in DSK export")
        if permutation_sign(faces[first]) == permutation_sign(faces[second]):
            raise ValueError("same-winding duplicate triangle in DSK export")
        remove[first] = remove[second] = True
        removed_pairs.append((first + 1, second + 1))
    del canonical, keys, order, ordered, same

    retained = faces[~remove]
    # Prove that every undirected edge occurs exactly twice and that its two
    # retained faces traverse it in opposite directions. This is the topology
    # contract enforced by LiouvilleCOPS before normals enter the physics.
    edge_count = 3 * len(retained)
    edges = np.empty((edge_count, 2), dtype=np.int32)
    edges[0::3] = retained[:, (0, 1)]
    edges[1::3] = retained[:, (1, 2)]
    edges[2::3] = retained[:, (2, 0)]
    direction = edges[:, 0] < edges[:, 1]
    edges.sort(axis=1)
    edge_type = np.dtype([("a", "<i4"), ("b", "<i4")])
    edge_keys = np.ascontiguousarray(edges).view(edge_type).reshape(-1)
    edge_order = np.argsort(edge_keys, kind="mergesort")
    sorted_keys = edge_keys[edge_order]
    if (edge_count % 2 or
            np.any(sorted_keys[0::2] != sorted_keys[1::2]) or
            (edge_count > 2 and np.any(sorted_keys[1:-1:2] == sorted_keys[2::2]))):
        raise ValueError("surface remains open or non-manifold after exact-fin removal")
    sorted_direction = direction[edge_order]
    if np.any(sorted_direction[0::2] == sorted_direction[1::2]):
        raise ValueError("surface winding remains inconsistent after exact-fin removal")
    return retained, removed_pairs


def convert_export(export_path, output_path, expected_vertices, expected_plates,
                   source_name, source_digest):
    """Convert a NAIF vertex-facet export into an atomic BDF candidate."""
    vertices = 0
    plates = 0
    faces_started = False
    face_table = np.empty((expected_plates, 3), dtype=np.int32)
    output_path = Path(output_path)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
                mode="w", encoding="ascii", newline="\n", delete=False,
                dir=str(output_path.parent), prefix=output_path.name + ".",
                suffix=".tmp") as target:
            temporary = Path(target.name)
            target.write("$ LiouvilleCOPS NASTRAN surface mesh\n")
            target.write("$ Converted from SPICE DSK: %s\n" % source_name)
            target.write("$ Source SHA256: %s\n" % source_digest)
            target.write("$ DSK vertices were converted from km to m; plate order is preserved.\n")
            with Path(export_path).open("r", encoding="ascii") as source:
                for line_number, line in enumerate(source, 1):
                    fields = line.split()
                    if not fields:
                        continue
                    kind = fields[0]
                    if kind == "v":
                        if faces_started or len(fields) != 4:
                            raise ValueError("invalid/out-of-order vertex at export line %d" %
                                             line_number)
                        try:
                            xyz = tuple(float(value) for value in fields[1:])
                        except ValueError:
                            raise ValueError("invalid vertex number at export line %d" %
                                             line_number) from None
                        if not all(value == value and abs(value) != float("inf")
                                   for value in xyz):
                            raise ValueError("non-finite vertex at export line %d" % line_number)
                        vertices += 1
                        write_grid(target, vertices, xyz)
                    elif kind == "f":
                        faces_started = True
                        if len(fields) != 4:
                            raise ValueError("non-triangular face at export line %d" % line_number)
                        try:
                            nodes = tuple(int(value) for value in fields[1:])
                        except ValueError:
                            raise ValueError("invalid face index at export line %d" %
                                             line_number) from None
                        if len(set(nodes)) != 3 or any(value < 1 or value > vertices
                                                      for value in nodes):
                            raise ValueError("degenerate/out-of-range face at export line %d" %
                                             line_number)
                        if plates >= expected_plates:
                            raise ValueError("export contains more plates than dskbrief reported")
                        face_table[plates] = nodes
                        plates += 1
                    else:
                        raise ValueError("unexpected dskexp record %r at line %d" %
                                         (kind, line_number))
        if vertices != expected_vertices or plates != expected_plates:
            raise ValueError(
                "export count mismatch: got %d vertices/%d plates, expected %d/%d" %
                (vertices, plates, expected_vertices, expected_plates))
        retained, removed_pairs = sanitize_faces(face_table)
        # Append faces only after the complete surface passes the topology
        # audit. Element IDs are compacted, while relative order of every
        # retained DSK plate remains unchanged.
        with temporary.open("a", encoding="ascii", newline="\n") as target:
            if removed_pairs:
                target.write("$ Removed %d reversed duplicate fin pairs (%d plates).\n" %
                             (len(removed_pairs), 2 * len(removed_pairs)))
            for element_id, nodes in enumerate(retained, 1):
                write_triangle(target, element_id, tuple(int(value) for value in nodes))
            target.write("ENDDATA\n")
        return temporary, vertices, plates, len(retained), removed_pairs
    except Exception:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
        raise


def atomic_json(path, payload):
    """Write a conversion manifest without exposing a partial JSON file."""
    path = Path(path)
    with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", newline="\n", delete=False,
            dir=str(path.parent), prefix=path.name + ".", suffix=".tmp") as stream:
        temporary = Path(stream.name)
        json.dump(payload, stream, indent=2, sort_keys=True)
        stream.write("\n")
    os.replace(str(temporary), str(path))


def convert_one(dsk, output, dskbrief, dskexp, frame_kernel, expected_frame):
    """Inspect, export, convert, hash, manifest, and publish one DSK."""
    dsk = Path(dsk).resolve()
    output = Path(output).resolve()
    manifest = output.with_suffix(output.suffix + ".conversion.json")
    label = dsk.with_suffix(".lbl")
    if not dsk.is_file():
        raise ValueError("DSK does not exist: %s" % dsk)
    if output.exists() or manifest.exists():
        raise ValueError("refusing existing conversion output: %s or %s" %
                         (output, manifest))
    output.parent.mkdir(parents=True, exist_ok=True)
    metadata = inspect_dsk(dskbrief, frame_kernel, dsk, expected_frame)
    dsk_digest = sha256(dsk)
    with tempfile.TemporaryDirectory(prefix="liouville_cops_dskexp_") as directory:
        exported = Path(directory) / (dsk.stem + ".obj")
        command = [str(dskexp), "-dsk", str(dsk), "-text", str(exported),
                   "-format", "vertex-facet", "-prec", "17"]
        completed = subprocess.run(
            command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, check=False)
        if completed.returncode or not exported.is_file():
            raise ValueError("dskexp failed for %s:\n%s" % (dsk, completed.stdout))
        temporary, vertices, input_plates, output_plates, removed_pairs = convert_export(
            exported, output, metadata["vertices"], metadata["plates"],
            dsk.name, dsk_digest)
    try:
        output_digest = sha256(temporary)
        payload = {
            "schema": SCHEMA,
            "schema_version": SCHEMA_VERSION,
            "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "source": {
                "path": str(dsk), "bytes": dsk.stat().st_size,
                "sha256": dsk_digest,
                "pds_label": ({"path": str(label), "bytes": label.stat().st_size,
                               "sha256": sha256(label)} if label.is_file() else None),
            },
            "dsk": metadata,
            "conversion": {
                "dskexp_command": command,
                "coordinate_scale_km_to_m": KM_TO_M,
                "vertices": vertices,
                "input_plates": input_plates,
                "output_plates": output_plates,
                "removed_reversed_duplicate_fin_pairs": removed_pairs,
                "removed_plates": input_plates - output_plates,
                "bdf_cards": ["GRID*", "CTRIA3"],
                "retained_plate_order_preserved": True,
                "volume_elements_written": False,
            },
            "output": {
                "path": str(output), "bytes": temporary.stat().st_size,
                "sha256": output_digest,
            },
        }
        os.replace(str(temporary), str(output))
        atomic_json(manifest, payload)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    return payload


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dsk", nargs="+", help="single-segment type-2 DSK file(s)")
    parser.add_argument("--out-dir", required=True,
                        help="existing/new directory for surface-only .bdf files")
    parser.add_argument("--frame-kernel", required=True,
                        help="SPICE frame kernel resolving the DSK frame code")
    parser.add_argument("--expected-frame", default="67P/C-G_CK")
    parser.add_argument("--cspice-root",
                        help="CSPICE toolkit root containing exe/dskbrief and exe/dskexp")
    parser.add_argument("--dskbrief", help="explicit NAIF dskbrief executable")
    parser.add_argument("--dskexp", help="explicit NAIF dskexp executable")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        dskbrief = resolve_naif_tool("dskbrief", args.dskbrief, args.cspice_root)
        dskexp = resolve_naif_tool("dskexp", args.dskexp, args.cspice_root)
        frame_kernel = Path(args.frame_kernel).expanduser().resolve()
        if not frame_kernel.is_file():
            raise ValueError("frame kernel does not exist: %s" % frame_kernel)
        out_dir = Path(args.out_dir).expanduser().resolve()
        outputs = [out_dir / (Path(name).stem + ".bdf") for name in args.dsk]
        if len(set(outputs)) != len(outputs):
            raise ValueError("input basenames produce duplicate output names")
        # Preflight every destination before performing the first conversion.
        for output in outputs:
            manifest = output.with_suffix(output.suffix + ".conversion.json")
            if output.exists() or manifest.exists():
                raise ValueError("refusing existing conversion output: %s or %s" %
                                 (output, manifest))
        for source, output in zip(args.dsk, outputs):
            print("Converting %s -> %s" % (source, output), flush=True)
            report = convert_one(
                source, output, dskbrief, dskexp, frame_kernel,
                args.expected_frame)
            print("  %d vertices, %d input / %d retained plates, %d bytes, SHA256 %s" %
                  (report["conversion"]["vertices"],
                   report["conversion"]["input_plates"],
                   report["conversion"]["output_plates"],
                   report["output"]["bytes"],
                   report["output"]["sha256"]), flush=True)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Interrupted", file=sys.stderr)
        sys.exit(130)
