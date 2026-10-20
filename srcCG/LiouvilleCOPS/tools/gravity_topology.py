"""Independent combinatorial audit before disabling the package's HEAL check.

HEAL's absolute ray-start epsilon is smaller than an ulp on kilometre-scale
meshes expressed in metres. Never let that heuristic rewrite production faces.
This audit checks a single outward closed genus-zero triangulated manifold.
It does not prove geometric non-self-intersection; shape provenance remains
required. No winding correction or vertex welding is performed here.
"""
import hashlib

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

ORIENTATION_POLICY = "DISABLE after independent outward closed manifold audit v1"


def require_outward_closed_body(vertices, faces, label="gravity mesh"):
    """Refuse bad topology, including disconnected vertex links (pinches).

    Sparse connectivity keeps the audit linear after the edge sort, without
    any ray casting or quadratic face-pair search. Signed volume is computed
    relative to a nearby origin to reduce translation-induced cancellation.
    """
    v = np.asarray(vertices, dtype=np.float64)
    raw = np.asarray(faces)
    def refuse(reason):
        raise ValueError("%s refused: %s" % (label, reason))
    if v.ndim != 2 or v.shape[1] != 3 or len(v) < 4 or not np.isfinite(v).all():
        refuse("vertices must be a finite (n,3) table")
    if raw.ndim != 2 or raw.shape[1] != 3 or len(raw) < 4 or raw.dtype.kind not in "iu":
        refuse("faces must be an integer (n,3) table")
    f = np.asarray(raw, dtype=np.int64)
    if np.any(f < 0) or np.any(f >= len(v)):
        refuse("face vertex index outside vertex table")
    if len(np.unique(f)) != len(v):
        refuse("unreferenced vertices")
    a, b, c = (v[f[:, i]] for i in range(3))
    area2 = np.linalg.norm(np.cross(b - a, c - a), axis=1)
    if not np.isfinite(area2).all() or np.any(area2 <= 0):
        refuse("nonfinite or zero-area faces")
    if len(np.unique(np.sort(f, axis=1), axis=0)) != len(f):
        refuse("duplicate faces")
    directed = np.stack((f, np.roll(f, -1, axis=1)), axis=-1).reshape(-1, 2)
    undirected = np.sort(directed, axis=1)
    order = np.lexsort((undirected[:, 1], undirected[:, 0]))
    ordered = undirected[order]
    starts = np.r_[0, np.flatnonzero(np.any(ordered[1:] != ordered[:-1], axis=1)) + 1]
    counts = np.diff(np.r_[starts, len(ordered)])
    if np.any(counts != 2):
        refuse("not closed: every edge must have exactly two incident faces")
    left, right = order[0::2], order[1::2]
    if np.any(directed[left] != directed[right, ::-1]):
        refuse("inconsistent face winding across shared edges")
    nfaces = len(f)
    face_graph = coo_matrix((np.ones(len(left)), (left // 3, right // 3)),
                            shape=(nfaces, nfaces)).tocsr()
    bodies = connected_components(face_graph, directed=False, return_labels=False)
    if bodies != 1:
        refuse("expected one connected body, found %d" % bodies)
    # Each face corner is a node in the link of its mesh vertex. Across a
    # shared edge the starting corner meets the opposite ending corner.
    # With two faces per edge every link has degree two; one component per
    # vertex therefore proves each link is a single circle, not a pinch.
    next_corner = (np.arange(3 * nfaces) // 3) * 3 + (np.arange(3 * nfaces) + 1) % 3
    rows = np.r_[left, next_corner[left]]
    cols = np.r_[next_corner[right], right]
    corner_graph = coo_matrix((np.ones(len(rows)), (rows, cols)),
                              shape=(3 * nfaces, 3 * nfaces)).tocsr()
    links = connected_components(corner_graph, directed=False, return_labels=False)
    if links != len(v):
        refuse("nonmanifold vertex links (pinched surface)")
    euler = len(v) - len(left) + nfaces
    if euler != 2:
        refuse("expected genus-zero surface, Euler characteristic %d" % euler)
    origin = np.mean(v, axis=0)
    volume = float(np.sum(np.einsum("ij,ij->i", a-origin,
                                    np.cross(b-origin, c-origin))) / 6.)
    if not np.isfinite(volume) or volume <= 0:
        refuse("winding is inward or volume is nonpositive/nonfinite")
    return {"orientation_policy": ORIENTATION_POLICY, "closed_consistent": True,
            "connected_bodies": int(bodies), "vertex_links_manifold": True,
            "euler_characteristic": int(euler), "signed_volume_m3": volume,
            "zero_area_faces": 0, "unreferenced_vertices": 0,
            "vertex_count": len(v), "face_count": nfaces,
            "face_order_sha256": hashlib.sha256(np.asarray(f, dtype="<i8").tobytes()).hexdigest(),
            "vertices_sha256": hashlib.sha256(np.asarray(v, dtype="<f8").tobytes()).hexdigest(),
            "self_intersection_test": "not performed; requires trusted shape provenance"}
