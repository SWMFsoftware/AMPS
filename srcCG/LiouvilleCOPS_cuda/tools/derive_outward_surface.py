#!/usr/bin/env python3
"""Write an outward-wound, surface-only copy of a consistently inward-wound BDF.

Decision of 2026-10-19 (ROADMAP §8). ``SHAP5_Volume.bdf`` is a Gmsh volume mesh. Its
100,000 CTRIA3 boundary triangles form a closed, consistently oriented manifold
(Euler 2) whose winding is uniformly *inward*: the signed volume is negative.
The repaired CPU builder correctly refuses such a body, because it never
re-orients a mesh. Reversing *every* triangle is an exact, global change of
orientation, unlike HEAL's per-face flips, which depend on rounding. It yields
the outward body that the builder and the GPU evaluator accept.

The derived file is written so that the reference reader
(``build_physical_activity_templates.read_nastran_surface``) returns:

* **vertices:** byte-identical to the source's. Every ``GRID``/``GRID*`` card
  and its continuation is copied verbatim, so no coordinate is re-printed or
  re-rounded.
* **faces:** exactly ``source_faces[:, [0, 2, 1]]``, in the same face order. Each
  ``CTRIA3`` keeps its element and property ids; nodes n2 and n3 are swapped.

All other cards (CTETRA volume elements, comments) are omitted; the reader
ignores them anyway. Before anything is written:

* the source must pass the audit except for its sign of volume (a closed,
  consistently oriented manifold with negative volume);
* the result must pass both the CPU reference audit
  (``gravity_topology.require_outward_closed_body``) and ``cops_cuda.gravity_body``.

The output and its provenance JSON are created exclusively; an existing file
is never replaced.

Usage (CPU or CUDA interpreter)::

    python3 tools/derive_outward_surface.py SOURCE.bdf OUTPUT.bdf
"""

import argparse
import hashlib
from pathlib import Path
import sys

import numpy as np

APP_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_ROOT))
from cops_cuda import gravity_body, provenance, reference  # noqa: E402


def reversed_ctria3(line):
    """Swap nodes 2 and 3 of one free-field CTRIA3 card (the reader splits on whitespace)."""
    fields = line.split()
    if len(fields) < 6 or fields[0] != "CTRIA3":
        raise ValueError("unsupported CTRIA3 card: %r" % line)
    if len(fields) > 6:
        # Extra fields (theta/material angle, offsets) would be dropped by a
        # naive rewrite; the source has none, so refuse rather than guess.
        raise ValueError("CTRIA3 card with optional fields is not supported: %r" % line)
    eid, pid, n1, n2, n3 = fields[1:6]
    return "CTRIA3  %s %s %s %s %s\n" % (eid, pid, n1, n3, n2)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("source", help="consistently inward-wound closed NASTRAN surface")
    parser.add_argument("output", help="new outward-wound surface BDF (must not exist)")
    args = parser.parse_args()
    source, output = Path(args.source).resolve(), Path(args.output).resolve()
    sidecar = output.with_suffix(output.suffix + ".provenance.json")
    if output.exists() or sidecar.exists():
        raise SystemExit("refusing to replace %s or %s" % (output, sidecar))

    builder = reference.load("build_physical_activity_templates")
    topology = reference.load("gravity_topology")
    V, F = builder.read_nastran_surface(source)
    before = gravity_body.topology(V, F)
    if not (before["consistently_oriented_closed_manifold"] and
            before["euler_characteristic"] == 2 and before["signed_volume_m3"] < 0.0 and
            before["zero_area_faces"] == 0 and before["unreferenced_vertices"] == 0):
        raise SystemExit("source is not a closed, consistently INWARD-wound manifold: %s" % before)

    # Copy GRID cards verbatim; rewrite CTRIA3 cards; drop everything else.
    counts = {"grid_cards": 0, "ctria3_cards": 0, "dropped_lines": 0}
    output.parent.mkdir(parents=True, exist_ok=True)
    with source.open("r", errors="strict") as stream, output.open("x") as out:
        out.write("$ Outward-wound surface derived from %s (sha256 %s)\n" %
                  (source.name, provenance.sha256_file(source)))
        out.write("$ by LiouvilleCOPS_cuda/tools/derive_outward_surface.py: CTRIA3 nodes 2,3 swapped\n")
        lines = iter(stream)
        for line in lines:
            if line.startswith("GRID*"):
                out.write(line)
                out.write(next(lines))          # the GRID* continuation card
                counts["grid_cards"] += 1
            elif line.startswith("GRID"):
                out.write(line)
                counts["grid_cards"] += 1
            elif line.startswith("CTRIA3"):
                out.write(reversed_ctria3(line))
                counts["ctria3_cards"] += 1
            else:
                counts["dropped_lines"] += 1

    # Verify through the reference reader and both audits.
    V2, F2 = builder.read_nastran_surface(output)
    checks = {"vertices_bitwise_identical": bool(V2.tobytes() == V.tobytes()),
              "faces_equal_source_with_nodes_2_3_swapped":
                  bool(np.array_equal(F2, F[:, [0, 2, 1]]))}
    after = gravity_body.GravityBody(V2, F2, 1.0).audit        # raises if inadmissible
    cpu_audit = topology.require_outward_closed_body(V2, F2)
    checks["volume_is_negated_source_volume"] = bool(
        np.isclose(after["signed_volume_m3"], -before["signed_volume_m3"], rtol=1e-12, atol=0.0))
    if not all(checks.values()):
        output.rename(output.with_suffix(output.suffix + ".rejected"))
        raise SystemExit("verification failed: %s" % checks)

    record = {"schema": "liouville_cops_cuda.derived_surface", "schema_version": 1,
              "created_utc": provenance.utc_now(),
              "decision": "user 2026-10-19: test the gravity model on SHAP5_Volume via an "
                          "exact global winding reversal (ROADMAP §8)",
              "source": provenance.file_identity(source, "source_volume_mesh"),
              "output": provenance.file_identity(output, "derived_outward_surface"),
              "transformation": "CTRIA3 (n1,n2,n3) -> (n1,n3,n2) for every face; GRID cards "
                                "copied verbatim; all other cards dropped",
              "cards": counts, "verification": checks,
              "topology_source": before, "topology_output": after,
              "cpu_reference_audit": cpu_audit if isinstance(cpu_audit, dict) else str(cpu_audit),
              "tool": provenance.file_identity(Path(__file__), "tool"),
              "reference_root": str(reference.bind())}
    provenance.write_json_exclusive(sidecar, record)
    print("wrote %s (%d vertices, %d faces, volume %.6g m^3)" %
          (output, after["vertices"], after["faces"], after["signed_volume_m3"]))
    print("provenance %s" % sidecar)
    return 0


if __name__ == "__main__":
    sys.exit(main())
