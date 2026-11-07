# Core-mesh cut-cell fixtures

`sparse_out_of_order_tetrahedron.bdf` is a deterministic free-field subset of
NASTRAN used by `MESH-CUTCELL-NASTRAN-ROUNDTRIP-001`. Its four `GRID*` records
use the sparse IDs `42`, `900`, `7`, and `305` in deliberately nonnumeric
order. The four outward-wound `CTRIA3` records form the closed unit right
tetrahedron and carry property/AMPS face attributes 11 through 14.

AMPS accepts this bounded surface-file contract:

- the first record is a comment/header, followed by two-line `GRID*` records
  and then `CTRIA3` records;
- a `GRID*` record is `GRID* id coordinate-system x y`, followed by
  `* z`; IDs are keys and need not be dense or ordered;
- a `CTRIA3` record is `CTRIA3 element-id property-id node-0 node-1 node-2`;
  the node order defines the face normal by
  `(x1-x0) cross (x2-x0)`, and every referenced node must exist in the same
  input file;
- coordinates and the reader's conversion factor must use consistent length
  units. The fixture is dimensionless; the test reads it with factor 1;
- duplicate node IDs and undefined face references are typed fatal input
  errors. They are not repaired or silently renumbered.

`CutCell::SaveNastranSurfaceMeshLongFormat` writes the same subset. It assigns
deterministic consecutive one-based serialization IDs, so read/write/read
preserves double-precision coordinates, topology, winding, and face
attributes. Sparse input IDs, input record ordering, element IDs, comments,
and per-file `MeshFileID` provenance are intentionally not round-trip state.
The test compares the promised state and does not treat textual identity as a
format requirement.
