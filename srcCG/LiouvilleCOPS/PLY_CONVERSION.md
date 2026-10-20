# Binary PLY surface conversion

`tools/convert_ply_to_bdf.py` reads binary big- or little-endian PLY 1.0 with
scalar vertex properties and one triangular vertex-index list per face. It
writes the existing LiouvilleCOPS NASTRAN representation: large-field `GRID*`
coordinates in metres, one-based `CTRIA3` indices, and `ENDDATA`.

The converter preserves vertex order, face order, winding, and coordinates
except for the explicitly requested unit scaling and decimal formatting.
There is no frame transform, decimation, smoothing, welding, hole filling,
duplicate deletion, or inferred correction of topology. Non-coordinate scalar
vertex attributes are not transferred to BDF. A PLY is not a SPICE frame
definition: the frame argument records the provenance, not a verified
alignment with the mission trajectory.

Example (replace placeholders with actual paths and documented metadata):

```bash
python3 tools/convert_ply_to_bdf.py /path/source.ply /path/output.bdf \
  --input-units km --frame CHEOPS \
  --units-and-frame-evidence 'Reference to source metadata and unit checks'
```

`--input-units` and evidence are mandatory because PLY does not define physical
units. Do not use `km` simply because a file is named after a comet.

The default conversion requires the existing audited, outward, closed,
single-body, genus-zero gravity topology. It rejects malformed/truncated
payloads, non-triangular records, bad vertex indices, and nonfinite coordinates.
Existing outputs are refused. The accompanying `.bdf.conversion.json` records
source/output/converter SHA256, original header, units, frame evidence, counts,
bounds, and topology results.

`--allow-failed-gravity-audit` is for faithful format conversion of a surface
that fails the current gravity admission policy. It prominently records the
failure in both BDF comments and JSON; it does NOT authorize running gravity or
prove that the C++ forward model accepts the surface. Investigate any failure
before scientific use. In particular, a genus-zero restriction and actual open
or non-manifold geometry are different issues and must not be conflated.

The converter uses the same 16-character coordinate formatting as the DSK
converter. Binary float32 inputs are read into float64 before scaling;
conversion does not create extra physical precision.

Tests:

```bash
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 \
  python3 -m unittest discover -s tests -p test_convert_ply_to_bdf.py -v
```

New shape files require new transport kernels, activity-template masks and
gravity products. A successful format conversion is not a validation of the
shape's epoch, local surface accuracy, or frame alignment for the final descent.

## SHAP8 combined 12M product

The MiARD product is `cg-spg_mspcd-shap8-v0.4cmb-12Mfacets.ply`, with
6,182,139 vertices and 12,364,290 triangles. D1.1 v7 section 2.2 documents the
Cheops frame. Its binary PLY header does not declare units. Raw bounds span
5.076575 coordinate units, signed volume is 18.559092 coordinate-units cubed,
and a sampled median edge length is 0.002559 coordinate units. Interpreting
these as kilometres matches the known nucleus dimensions and the report's
approximately 2 m sampling. This unit inference must remain explicit in the
conversion provenance, not be presented as a literal PLY header declaration.

The original mesh passes the audit's edge closure, consistent winding,
connectivity, and manifold-vertex-link stages, but has Euler characteristic -6
(genus 4), not the genus zero currently required by `gravity_topology.py`.
Consequently, a faithful conversion requires the explicit
`--allow-failed-gravity-audit` option. Do not automatically relax the gravity
policy, fill its handles, or apply HEAL. Acceptance by the C++ surface reader
and acceptance by the gravity builder are separate checks.

The converted 12M file passed the actual `src/Mesh.cpp` `ReadNastran` and
`FinalizeGeometry(true, true)` checks on 2026-10-20: all 6,182,139 vertices
and 12,364,290 triangles were read, the surface was closed and consistently
oriented, and the original winding was kept. Computed area was
5.349535e7 m² and volume 1.855909e10 m³. This does not remove the separate
genus-zero gravity-policy restriction above.

The retained data directory is
`/data/vtenishe/rosetta_data_end_of_mission_cops/shape_models/SHAP8_MIARD`.
Its README and `download_provenance.json` identify the archived download,
all retained hashes, and validation details. The original live MiARD server
was unavailable, so the original ZIP was retrieved from its 2020-09-07
Internet Archive capture. All ZIP members passed CRC validation.

## 20M and 44M extensions and eight-digit IDs

The 20M product has 9,999,746 vertices and 19,999,504 triangles. Like the 12M
product, its audited topology has Euler characteristic -6 rather than the
current gravity policy's required 2.

The 44M product has 21,983,227 vertices and 43,966,474 triangles. Its original
PLY contains nine exactly zero-area triangles. The source-geometry report
beside the conversion lists their zero-based facet indices and coordinates.
The faithful BDF is not a geometry repair and must not be treated as a
simulation-ready shape until those degeneracies are addressed and revalidated.

For node IDs at or above 10,000,000, adjacent standard eight-character CTRIA3
fields no longer have intervening whitespace. `src/Mesh.cpp` and the Python
`read_nastran_surface` now fall back to fixed-column parsing when whitespace
parsing fails. The writer and mesh geometry are unchanged. Rebuild
`liouville_cops` before using the 44M file with older binaries.

The sparse-eight-digit fixture in `test_convert_ply_to_bdf.py` validates the
Python path; set `COPS_MESH_PROBE` to an executable calling the production
C++ `ReadNastran` and `FinalizeGeometry(true,true)` to exercise that path too.

The actual 20M file passed that C++ check. The actual 44M file parsed all
21,983,227 vertices and 43,966,474 faces, then correctly failed geometry
finalization with `Error: degenerate surface triangle`. No degeneracy was
silently repaired and no existing scientific configuration was switched to
either model.

## Separate, surface-preserving 44M repair

`tools/repair_collinear_ply_bdf.py` addresses exactly collinear triangles,
not general mesh defects. It identifies the longest edge of each zero-area
source triangle and requires its third vertex to be strictly inside that
edge. The edge must have precisely one other incident triangle with opposite
winding. Replacing the two triangles with the other diagonal splits that
neighbor at the existing collinear vertex: the geometric surface is unchanged
in the source coordinates. No vertex is moved or deleted.

The script refuses an already-existing new diagonal, overlapping repairs,
wrong winding, open/nonmanifold repair edges, and degenerate replacements.
It checks conservation of oriented area. Near-collinear but nonzero triangles
are not modified. All original BDF GRID records are copied verbatim; the
existing ten-significant-digit decimal formatting is retained. Consequently,
exact geometric equivalence refers to the binary PLY source; rounding can
leave an extremely small difference between original and repaired BDF
triangulations. The actual output must pass C++ validation separately.

The original BDF and PLY remain unchanged. The repair checks their SHA256
against the conversion manifest, verifies every BDF triangle against the
original PLY ordering, and writes a separate `.bdf.repair.json` containing
every changed triangle (zero-based indices), hashes, and admission status.
Existing output is refused. An interrupted write leaves a `.partial` file,
not a completed BDF.

```bash
SHAP8_DIR=/data/vtenishe/rosetta_data_end_of_mission_cops/shape_models/SHAP8_MIARD
python3 tools/repair_collinear_ply_bdf.py \
  --ply "$SHAP8_DIR/cg-spg_mspcd-shap8-v0.4cmb-44Mfacets.ply" \
  --bdf "$SHAP8_DIR/cg-spg_mspcd-shap8-v0.4cmb-44Mfacets.bdf" \
  --output "$SHAP8_DIR/cg-spg_mspcd-shap8-v0.4cmb-44Mfacets_repaired.bdf"
```

Reproducible production-reader check (set `CSPICE_DIR` to your installation):

```bash
CSPICE_DIR=/home/vtenishe/SPICE/cspice
g++ -O2 -std=c++11 -I src -I "$CSPICE_DIR/include" \
  tests/mesh_geometry_probe.cpp src/Mesh.cpp \
  "$CSPICE_DIR/lib/cspice.a" -lm -o /tmp/cops_mesh_geometry_probe
/tmp/cops_mesh_geometry_probe \
  "$SHAP8_DIR/cg-spg_mspcd-shap8-v0.4cmb-44Mfacets_repaired.bdf"
```

After successful validation, select the repaired forward-simulation mesh in
your input file:

```ini
ShapeModel = /data/vtenishe/rosetta_data_end_of_mission_cops/shape_models/SHAP8_MIARD/cg-spg_mspcd-shap8-v0.4cmb-44Mfacets_repaired.bdf
ShapeModelFormat = nastran
```

Rebuild old binaries to include the eight-digit CTRIA3 fix. Do not reuse
transport matrices, per-facet activity files, masks, or gravity caches from
SHAP5, another SHAP8 resolution, or the unrepaired connectivity. This is a
44-million-face model: mesh validation alone is not a full transport/fit
performance or memory-capacity test. The existing strict gravity policy is
unchanged; C++ geometry acceptance does not imply gravity/cliff readiness.

Tests:

```bash
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 \
  python3 -m unittest discover -s tests -p test_repair_collinear_ply_bdf.py -v
```

The actual repaired 44M BDF passed the production C++ check on 2026-10-20:
21,983,227 nodes, 43,966,474 faces, closed and consistently oriented, original
winding retained, area `5.379692e7 m²`, volume `1.855906e10 m³`. Exactly nine
flips changed 18 face records; zero vertices moved. Its SHA256 is
`d80d5915345c209819d29d4369b640a97c36f7f0ecedc23c36d1237ec8f097e0`.
The adjacent `_repaired.cpp_mesh_validation.txt` and `.bdf.repair.json`
retain the validation and detailed repair evidence. Original BDF comments
remain unchanged and describe the original conversion's audit failure;
consult the repair manifest for the repaired product's current status.

The subsequent strict source-coordinate topology audit passed zero-area,
duplicate-face, edge-closure, winding, connected-body, and manifold-vertex-link
checks, then rejected Euler characteristic **−10**, corresponding to genus 6,
rather than the required genus zero. Thus the repair resolves forward-model
geometry admission, **not gravity/cliff admission**. A separate reviewed
extension of gravity validation is needed before using this shape in the
cliff-enhanced pipeline. Do not fill topological handles merely to satisfy
the current policy or bypass it without establishing physical validity.
Neither audit performs a geometric self-intersection test.
