#!/usr/bin/env python3
"""Repair exactly collinear PLY triangles by surface-preserving edge flips.

The input BDF must be the sequential-node/face conversion of the input PLY.
All BDF vertices and all unaffected records are retained byte-for-byte. This
does not smooth, weld, fill holes, change frames, or relax gravity admission.
The output still requires validation by the production C++ mesh reader.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from convert_ply_to_bdf import read_ply
from convert_dsk_to_bdf import sha256, write_triangle


def plan_repairs(vertices, faces, chunk=500000):
    """Only accept isolated flips with an absent new diagonal.

    A zero-area (a,b,c), with c strictly inside its longest edge ab, and
    adjacent (b,a,d) can be replaced by (c,a,d),(c,d,b). These partition
    exactly the same triangle abd; no new physical surface is introduced.
    Exact collinearity is required, not a tolerance-based near-zero test.
    """
    v = np.asarray(vertices)
    zeros = []
    for start in range(0, len(faces), chunk):
        f = faces[start:start+chunk]
        cross = np.cross(v[f[:, 1]]-v[f[:, 0]], v[f[:, 2]]-v[f[:, 0]])
        zeros.extend((start+np.flatnonzero(np.all(cross == 0, axis=1))).tolist())
    changes, records = {}, []
    for index in zeros:
        old = np.asarray(faces[index], dtype=np.int64)
        lengths = [np.dot(v[old[(k+1)%3]]-v[old[k]],
                          v[old[(k+1)%3]]-v[old[k]]) for k in range(3)]
        k = int(np.argmax(lengths))
        a, b, c = map(int, np.roll(old, -k))
        ab, ac = v[b]-v[a], v[c]-v[a]
        fraction = float(np.dot(ac, ab)/np.dot(ab, ab))
        if not 0 < fraction < 1:
            raise ValueError("collinear vertex is not strictly interior")
        neighbors = []
        for start in range(0, len(faces), chunk):
            f = faces[start:start+chunk]
            mask = np.any(f == a, axis=1) & np.any(f == b, axis=1)
            neighbors.extend((start+np.flatnonzero(mask)).tolist())
        if len(neighbors) != 2 or index not in neighbors:
            raise ValueError("longest edge is not two-sided")
        neighbor = next(i for i in neighbors if i != index)
        other = list(map(int, faces[neighbor]))
        j = other.index(b)
        if other[(j+1)%3] != a:
            raise ValueError("inconsistent edge orientation")
        d = other[(j+2)%3]
        if d in (a, b, c) or index in changes or neighbor in changes:
            raise ValueError("overlapping or collapsed repair")
        for start in range(0, len(faces), chunk):
            f = faces[start:start+chunk]
            if np.any(np.any(f == c, axis=1) & np.any(f == d, axis=1)):
                raise ValueError("new diagonal already exists")
        new = [(c, a, d), (c, d, b)]
        normal = np.cross(v[a]-v[b], v[d]-v[b])
        normals = [np.cross(v[y]-v[x], v[z]-v[x]) for x, y, z in new]
        if any(np.dot(n, normal) <= 0 for n in normals):
            raise ValueError("replacement is degenerate or reverses orientation")
        if not np.allclose(sum(normals), normal, rtol=1e-12, atol=1e-18):
            raise ValueError("oriented area is not conserved")
        changes[index], changes[neighbor] = new
        records.append(dict(face_index=index, neighbor_index=neighbor,
                            old_faces=[old.tolist(), other],
                            new_faces=[list(t) for t in new],
                            interior_fraction=fraction))
    return changes, records


def repair(ply, bdf, output):
    """Check source identity before streaming a new BDF; never overwrite."""
    ply, bdf, output = map(Path, (ply, bdf, output))
    manifest = output.with_suffix(output.suffix+".repair.json")
    if output.exists() or manifest.exists():
        raise ValueError("output already exists")
    conversion = json.loads(bdf.with_suffix(bdf.suffix+".conversion.json").read_text())
    # Compare hashes using the converter's actual schema, not file basenames.
    ply_hash, bdf_hash = sha256(ply), sha256(bdf)
    if (conversion.get("source_sha256") != ply_hash
            or conversion.get("output_sha256") != bdf_hash):
        raise ValueError("conversion manifest does not match input hashes")
    xyz, faces, _ = read_ply(ply)
    changes, records = plan_repairs(xyz, faces)
    print("Planned %d surface-preserving flips" % len(records), flush=True)
    if not records:
        raise ValueError("no exactly collinear triangles found")
    # Validate every original triangle against the source, including records
    # not being changed. GRID records are copied verbatim from hashed input.
    count = 0
    partial = output.with_suffix(output.suffix+".partial")
    with bdf.open() as source, partial.open("x") as target:
        for line in source:
            if line[:8].strip() != "CTRIA3":
                target.write(line)
                continue
            fields = [int(line[8+8*k:16+8*k]) for k in range(5)]
            if fields[0] != count+1 or fields[2:] != (faces[count]+1).tolist():
                raise ValueError("BDF triangle order differs from original PLY")
            if count in changes:
                write_triangle(target, count+1, tuple(i+1 for i in changes[count]))
            else:
                target.write(line)
            count += 1
    if count != len(faces):
        raise ValueError("incomplete BDF triangle table")
    partial.rename(output)
    result = dict(source_ply=str(ply), source_ply_sha256=ply_hash,
                  source_bdf=str(bdf), source_bdf_sha256=bdf_hash,
                  output=str(output), output_sha256=sha256(output),
                  vertex_count=len(xyz), face_count=count,
                  vertices_moved=0, changed_faces=len(changes),
                  repairs=records, indices="zero-based",
                  validation="pending production C++ closed-geometry probe",
                  gravity_admitted=False,
                  gravity_note="Existing genus-zero gravity policy is unchanged.")
    manifest.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ply", required=True)
    parser.add_argument("--bdf", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    repair(args.ply, args.bdf, args.output)
