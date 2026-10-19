"""The homogeneous gravity polyhedron: reading, topology audit, identity.

Python 3.8 compatible; NumPy only. This module is original code of this
project. The GPLv3-derived parts of the gravity port are confined to
``gravity.py`` and ``kernels/gravity.cu`` (see ``GRAVITY_REFERENCE.md`` §7).

Why a topology audit instead of the package's HEAL
--------------------------------------------------
``polyhedral_gravity``'s HEAL/VERIFY ray test misjudges about 46% of the
faces of a kilometre-scale mesh given in metres (``GRAVITY_REFERENCE.md``
§5.1). The port therefore never re-orients faces. It evaluates the winding
as given, and admits a body only if an exact combinatorial and geometric
audit proves that winding outward and consistent:

* every undirected edge is shared by exactly two faces (closed surface);
* every directed edge occurs exactly once (consistent orientation: the two
  faces of an edge traverse it in opposite directions);
* Euler characteristic V - E + F = 2 (a sphere-like surface);
* no face has zero area;
* signed volume > 0 (the consistent winding is outward, not inward).

This is the C00 ``mesh_topology`` audit, which admitted the 788k body, moved
here unchanged. It is exact integer bookkeeping plus one volume sum, so,
unlike the ray test, it does not depend on coordinate magnitude.

In addition, ``GravityBody`` runs the CPU reference's own audit,
``gravity_topology.require_outward_closed_body``. That audit is stricter: it also
requires unique triangles, one connected body and manifold vertex links (no
pinches). A body is admitted only if both pass. The reference result is kept
verbatim (``reference_audit``), because the CPU builder reuses a gravity cache
only when the cache's recorded ``topology`` equals that result and its
``orientation_policy`` equals the reference constant. A GPU-produced cache
(C04) therefore records exactly these values, and the unchanged CPU workflow
accepts it.
"""

import numpy as np

from . import provenance, reference
from .errors import IdentityError


def topology(vertices, faces):
    """Closedness/orientation facts, independent of any ray casting."""
    v = np.asarray(vertices, dtype=np.float64)
    f = np.asarray(faces, dtype=np.int64)
    directed = np.concatenate((f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]))
    undirected = np.sort(directed, axis=1)
    _, und_counts = np.unique(undirected, axis=0, return_counts=True)
    _, dir_counts = np.unique(directed, axis=0, return_counts=True)
    a, b, c = v[f[:, 0]], v[f[:, 1]], v[f[:, 2]]
    signed_volume = float(np.sum(np.einsum("ij,ij->i", a, np.cross(b, c))) / 6.0)
    area = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
    used = np.unique(f)
    return {"vertices": int(len(v)), "faces": int(len(f)),
            "unreferenced_vertices": int(len(v) - len(used)),
            "edges": int(len(und_counts)),
            "edges_not_shared_by_exactly_two_faces": int(np.count_nonzero(und_counts != 2)),
            "directed_edges_repeated": int(np.count_nonzero(dir_counts != 1)),
            "consistently_oriented_closed_manifold": bool(
                np.all(und_counts == 2) and np.all(dir_counts == 1)),
            "euler_characteristic": int(len(v) - len(und_counts) + len(f)),
            "signed_volume_m3": signed_volume, "area_min_m2": float(area.min()),
            "zero_area_faces": int(np.count_nonzero(area == 0.0)),
            "max_abs_coordinate_m": float(np.max(np.abs(v)))}


def admissible(audit):
    """Reasons a body fails the audit (empty list: admitted)."""
    reasons = []
    if not audit["consistently_oriented_closed_manifold"]:
        reasons.append("not a closed, consistently oriented manifold "
                       "(%d edges not shared by two faces, %d repeated directed edges)" %
                       (audit["edges_not_shared_by_exactly_two_faces"],
                        audit["directed_edges_repeated"]))
    if audit["euler_characteristic"] != 2:
        reasons.append("Euler characteristic %d != 2" % audit["euler_characteristic"])
    if audit["zero_area_faces"]:
        reasons.append("%d zero-area faces" % audit["zero_area_faces"])
    if not audit["signed_volume_m3"] > 0.0:
        reasons.append("signed volume %.6g <= 0: winding is inward" % audit["signed_volume_m3"])
    if audit["unreferenced_vertices"]:
        reasons.append("%d unreferenced vertices" % audit["unreferenced_vertices"])
    return reasons


class GravityBody:
    """A homogeneous polyhedron admitted for exact gravity evaluation.

    ``vertices`` are in metres, and ``faces`` index them in the evaluated
    winding. The constructor runs the audit and refuses an inadmissible body
    with ``IdentityError``. There is no repair mode.
    """

    def __init__(self, vertices, faces, density, source=None):
        self.vertices = np.ascontiguousarray(vertices, dtype=np.float64)
        self.faces = np.ascontiguousarray(faces, dtype=np.int64)
        if self.vertices.ndim != 2 or self.vertices.shape[1] != 3:
            raise ValueError("vertices must have shape (n, 3)")
        if self.faces.ndim != 2 or self.faces.shape[1] != 3:
            raise ValueError("faces must have shape (m, 3)")
        if len(self.faces) == 0 or self.faces.min() < 0 or self.faces.max() >= len(self.vertices):
            raise IdentityError("face indices outside the vertex table")
        if not np.isfinite(self.vertices).all():
            raise IdentityError("non-finite vertex coordinates")
        self.density = float(density)
        if not self.density > 0.0:
            raise ValueError("density must be positive")
        self.source = source
        self.audit = topology(self.vertices, self.faces)
        reasons = admissible(self.audit)
        if reasons:
            raise IdentityError("gravity body refused by the topology audit: " + "; ".join(reasons))
        gravity_topology = reference.load("gravity_topology")
        try:
            self.reference_audit = gravity_topology.require_outward_closed_body(
                self.vertices, self.faces)
        except ValueError as exc:
            raise IdentityError("gravity body refused by the reference audit: %s" % exc) from None
        self.orientation_policy = gravity_topology.ORIENTATION_POLICY

    @classmethod
    def from_bdf(cls, path, density):
        """Read a NASTRAN surface with the reference builder's own reader.

        The reader is ``build_physical_activity_templates.read_nastran_surface``.
        It numbers vertices in sorted GRID-id order and keeps the CTRIA3 face
        order, which is the order in which the campaign evaluated faces.
        """
        builder = reference.load("build_physical_activity_templates")
        vertices, faces = builder.read_nastran_surface(path)
        return cls(vertices, faces, density, source=provenance.file_identity(path, "gravity_shape"))

    @property
    def nfaces(self):
        return len(self.faces)

    def identity(self):
        """Content identity: the evaluated arrays, density and audit."""
        return {"source": self.source, "density_kg_m3": self.density,
                "vertices_sha256": provenance.array_sha256(self.vertices),
                "faces_sha256": provenance.array_sha256(self.faces),
                "orientation": "OUTWARDS (input winding, admitted by topology audit; no HEAL)",
                "orientation_policy": self.orientation_policy,
                "topology": self.audit, "reference_audit": self.reference_audit}
