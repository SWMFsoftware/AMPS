# Stage-10 shape examples

These six JSON files are complete, directly runnable `shieldSim.case.v2`
verification inputs. They demonstrate both non-slab geometry providers without
claiming that either provider is qualified for mission or anomaly analysis.
Every file defines its absorber material, ordered shielding, and requested
calculation functions directly in JSON. No one of those selections is compiled
into the shape provider.

The separate Boolean CSG enclosure in
`../stage10_analytic_csg_case.json` uses the same selection fields and defaults.

| Case | Provider | Sensitive geometry | Exact volume (cm³) | Entry definition |
|---|---|---:|---:|---|
| `analytic_box.json` | `analytic_csg.v1` | 20×12×2 mm box | 0.48 | local −z face |
| `analytic_cylinder.json` | `analytic_csg.v1` | r=5 mm, L=4 mm cylinder | π/10 | local −z face |
| `analytic_sphere.json` | `analytic_csg.v1` | r=5 mm sphere | π/6 | entire exterior |
| `mesh_box.json` | `tessellated_mesh.v1` | same 20×12×2 mm box, 12 triangles | 0.48 | two tagged −z triangles |
| `mesh_tetrahedron.json` | `tessellated_mesh.v1` | 10 mm right tetrahedron, 4 triangles | 1/6 | one tagged z=0 triangle |
| `mesh_octahedron.json` | `tessellated_mesh.v1` | axis radius 5 mm, 8 triangles | 1/6 | entire exterior |

Every committed case defaults to natural silicon, a 1 mm Al plus 2 mm HDPE
shield, and all five calculations shown in the supplied specification:

```json
"quantities": ["TID", "DDD", "n_eq", "LET", "H100/10"]
```

### Runtime selections in each input file

Edit the following paths in any of the six case files.

1. **Shielding:** `geometry.shield_stack.layers` is the ordered incident-to-
   absorber stack. Each of one through five layers defines `layer_id`,
   `material`, and either geometric or areal-density thickness:

   ```json
   "layers": [
     {
       "layer_id": "outer_al",
       "material": "Al",
       "thickness": {"type": "geometric", "value": 1.0, "unit": "mm"}
     },
     {
       "layer_id": "inner_hdpe",
       "material": "HDPE",
       "thickness": {"type": "areal_density", "value": 0.19, "unit": "g/cm2"}
     }
   ]
   ```

   Usable aliases from the screenshot are `Al`, `Cu`, `W`, `Ta`, `HDPE`,
   `BPE`, `Kevlar`, `Kapton`, `Water`, and `Si`. `CFRP`, `SiCComposite`, and
   lunar/Mars regolith are deliberately disabled templates until their actual
   project properties are supplied through `material_overrides`. Layer order
   is physical; each layer is limited to 100 mm derived geometric thickness.

2. **Absorber/detector material:** change the sensitive region's `material` at
   `geometry.absorbers[0].geometry.regions[*].material`. The region with
   `"role":"sensitive"` is scored; passive-region materials, if present, are
   selected independently. For example:

   ```json
   {
     "region_id": "detector",
     "role": "sensitive",
     "material": "Si"
   }
   ```

   Input-library aliases corresponding to the screenshot are `Skin`,
   `EyeLens`, `BFO`, `CNS`, `SoftTissue`, `Si`, `SiO2`, `SiC`, `GaAs`,
   `InGaAs`, `Ge`, and `Water`. Density/composition changes belong in the
   top-level `material_overrides` array; changing only a display label never
   changes physics.

3. **Calculation functions:** edit the top-level `quantities` array. The five
   names are exactly `TID`, `DDD`, `n_eq`, `LET`, and `H100/10`. `DDD` or
   `n_eq` requires the top-level `damage_response` object. `LET` requires
   `let_response`. Remove an unused response object when its corresponding
   calculation is removed; the strict parser rejects stale hidden settings.

The parser enforces the currently validated compatibility boundary:

| Calculation | Accepted absorber material in the current model |
|---|---|
| `TID` | `Al`, `Si`, or `Water`; proton/He-4 source domain |
| `DDD` | natural `Si` only; registered p/He-4/neutron response |
| `n_eq` | natural `Si` only; registered Stage-5 convention |
| `LET` | `Al`, `Si`, or `Water`; registered p/He-4/C-12 domain |
| `H100/10` | primary-proton shield-exit quantity; absorber-independent |

Thus silicon is the default because it is the only absorber common to all five
requested functions. The other screenshot materials are valid explicit
transport materials, but their TID/DDD/LET proxy domains have not been
qualified. The model rejects those unsupported combinations rather than
silently extending silicon or water response data.

Candidate Ge, SiO2, GaAs, Ge-DDD, and InGaAs-DDD reference packages are now
checksum-pinned for future qualification. Their presence does not change this
table. The exact data-ready, partial, undefined, and not-applicable states are
recorded in `tests/reference_material_functions/material_quantity_data_status.v1.json`.

One case file represents one shielding stack and one source definition. To
compare alternative shielding configurations, copy the case, give it a unique
`case_id` and `output.prefix`, and edit its layer array. This keeps each result
traceable to one unambiguous physical configuration.

The analytic and mesh boxes are coextensive and are compared through the full
same-seed Geant4/TID path by `stage10_tessellated_mesh`. The tetrahedron and
octahedron add non-box convex topology with closed-form volume, centroid, and
surface area. They are geometry/transport fixtures, not spacecraft models.

From the repository root, validate all examples without transport:

```bash
source /data/vtenishe/G4/install/bin/geant4.sh
for case_file in examples/stage10_shapes/*.json; do
  test_output/build/shieldSim --validate-case="$case_file" || exit 1
done
```

Run an individual example with outputs isolated under `test_output/`:

```bash
mkdir -p test_output/stage10_shapes/mesh_octahedron
cd test_output/stage10_shapes/mesh_octahedron
../../build/shieldSim \
  --case=../../../examples/stage10_shapes/mesh_octahedron.json
```

Every file uses only 20 histories for quick execution. Do not interpret its
Monte Carlo value as a converged physical result. Run the authoritative bounded
verification gates from `test_output/build`:

```bash
ctest --output-on-failure -R '^stage10_(analytic_csg|tessellated_mesh)$'
```

The mesh files use explicit millimetre coordinates and outward triangle
winding. Provider v1 accepts closed, connected, convex triangular shells only;
it does not import or repair STL/CAD. See the Stage-10 section of `RUN.md` for
the full contract and limitations.

The JSON gallery is generated by
`scripts/generate_stage10_shape_examples.py`. After an intentional generator
change, rerun it from the repository root and review every JSON diff:

```bash
python3 scripts/generate_stage10_shape_examples.py
```
