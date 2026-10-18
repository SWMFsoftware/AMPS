# Runnable material, shield, species, and shape configurations

These seven complete `shieldSim.case.v2` files show how the four independent
input selections fit together:

1. the incident particle species and energy are in
   `source.spectrum.nodes`;
2. the ordered shielding materials, thicknesses, and common transverse size
   are in `geometry.shield_stack`;
3. the scored absorber material and its dimensions/mesh are in the sensitive
   region under `geometry.absorbers[0].geometry`; and
4. the calculation functions are selected by the top-level `quantities`
   array, with mandatory `damage_response` and `let_response` blocks when
   applicable.

Every case is intentionally limited to 20 histories so it is quick to run.
That is enough to prove that the complete input reaches Geant4 and all selected
writers, but it is **not** enough for a converged radiation assessment. The
files use `purpose: "verification_test"` because the Stage-10 non-slab shape
providers have not passed an independent nontrivial-shape qualification gate.

## Configuration matrix

| File | Incident species | Sensitive absorber and shape | Shield, upstream to downstream | Selected calculations |
|---|---|---|---|---|
| `01_si_box_proton_all_quantities.json` | 50/150 MeV protons | Si analytic 20×12×2 mm box | 60×40 mm: 1.5 mm Al; 0.19 g/cm² HDPE; 0.1 mm Ta | TID, DDD, n_eq, LET, H100/10 |
| `02_si_cylinder_neutron_damage.json` | 1.05 MeV neutron | Si analytic r=5 mm, L=4 mm cylinder | 45×45 mm: 0.05 mm Al | DDD, n_eq |
| `03_si_sphere_alpha_tid_ddd_let.json` | 100 MeV total-energy He-4 | Si analytic r=5 mm sphere | 50×50 mm: 0.05 mm Kapton; 0.10 mm Al; 0.10 mm HDPE; 0.01 mm Ta; 0.05 mm water | TID, DDD, LET |
| `04_water_cylinder_proton_tid_let_h100.json` | 50/150 MeV protons | water analytic r=5 mm, L=4 mm cylinder | 40×30 mm: 0.5 mm Cu; 4.0 mm HDPE | TID, LET, H100/10 |
| `05_al_mesh_box_proton_tid_let_h100.json` | 50/150 MeV protons | Al 20×12×2 mm closed triangle-mesh box | 70×45 mm: 0.10 mm Kapton; 0.75 mm Al; 0.50 mm Kevlar; 0.095 g/cm² HDPE | TID, LET, H100/10 |
| `06_si_mesh_tetrahedron_neutron_damage.json` | 1.05 MeV neutron | Si 10 mm right tetrahedron, four triangles | 55×35 mm: 0.05 mm Al; 0.05 mm HDPE | DDD, n_eq |
| `07_si_mesh_octahedron_case_material_override.json` | 50/150 MeV protons | Si r=5 mm closed octahedral mesh | 65×65 mm: 0.27 g/cm² case-local Al; 2.0 mm HDPE | TID, DDD, n_eq, LET, H100/10 |

The first and seventh files demonstrate that all five names can be requested
together for natural silicon. Their proton-only `n_eq` is intentionally zero:
the model defines `n_eq` as the ASTM E722 neutron contribution and never folds
proton DDD into that number. Files 02 and 06 provide the complementary neutron
case in which `n_eq` is physically active. A single source cannot make all five
quantities active because TID's qualified source domain is proton/He-4 while
`n_eq` is neutron-only. Do not add a neutron node to a TID case; the parser
rejects that unsupported combination.

## Edit the input selections

### Incident species

Each discrete node declares its Geant4 particle identity, energy coordinate,
direction, and probability explicitly. For example, the neutron definition is:

```json
{
  "incident_id": "neutron_1p05MeV",
  "particle": {"name": "neutron", "Z": 0, "A": 1},
  "energy": {
    "coordinate": "total_kinetic_energy",
    "value": 1.05,
    "unit": "MeV"
  },
  "direction": {"mu": 1.0, "phi": {"value": 0.0, "unit": "rad"}},
  "weight": 1.0
}
```

Allowed exact response-node identities are proton `(1,1)`, alpha/He-4
`(2,4)`, neutron `(0,1)`, gamma `(0,0)`, and `ion` with positive explicit
`Z,A`. Gamma and general ions are raw-transport verification inputs, not TID
inputs. Node weights must be positive and sum to one within `1e-12`.

### Shield material, number of layers, thickness, and footprint

Layer array order is physical: element zero is the upstream layer first seen by
the incident radiation. One through five layers are accepted. Each layer must
use either a geometric thickness or an areal density:

```json
"shield_stack": {
  "size_x": {"value": 60.0, "unit": "mm"},
  "size_y": {"value": 40.0, "unit": "mm"},
  "transform": {
    "translation": {
      "x": {"value": 0.0, "unit": "mm"},
      "y": {"value": 0.0, "unit": "mm"},
      "z": {"value": -20.0, "unit": "mm"}
    },
    "rotation_quaternion_local_to_world": {"w": 1.0, "x": 0.0, "y": 0.0, "z": 0.0}
  },
  "layers": [
    {
      "layer_id": "outer_al",
      "material": "Al",
      "thickness": {"type": "geometric", "value": 1.5, "unit": "mm"}
    },
    {
      "layer_id": "middle_hdpe",
      "material": "HDPE",
      "thickness": {"type": "areal_density", "value": 0.19, "unit": "g/cm2"}
    }
  ]
}
```

Every layer's density comes from the selected material record and converts
areal density to geometric thickness as
`thickness_mm = 10 * areal_density_g_cm2 / density_g_cm3`. The resolved value
and material hash are archived in `*_resolved_case.json`. Each derived layer
thickness must be greater than zero and no more than 100 mm (10 cm).

### Absorber material and shape

For a non-slab provider, change only the region whose `role` is `sensitive`:

```json
{
  "region_id": "detector",
  "role": "sensitive",
  "material": "Water",
  "solid": {
    "node_id": "detector_cylinder",
    "type": "cylinder",
    "inner_radius": {"value": 0.0, "unit": "mm"},
    "outer_radius": {"value": 5.0, "unit": "mm"},
    "length": {"value": 4.0, "unit": "mm"}
  }
}
```

Analytic primitives are box, sphere, and z-axis cylinder. Mesh files define
every vertex in millimetres and every outward-wound triangular facet inline.
The current mesh provider accepts one connected, closed, convex shell; it does
not import or repair STL/CAD. `entry_surface` is part of the response: analytic
boxes/cylinders can select `local_negative_z_face`; a sphere uses
`entire_exterior`; a mesh uses tagged coplanar facets or `entire_exterior`.

The currently qualified quantity/material boundary remains:

| Quantity | Sensitive absorber and incident boundary |
|---|---|
| TID | natural Al, natural Si, or water; proton or He-4 |
| DDD | natural Si; registered proton, He-4, or neutron damage response |
| n_eq | natural Si; ASTM E722 neutron contribution only |
| LET | natural Al, natural Si, or water; registered proton, He-4, or C-12 domain |
| H100/10 | primary-proton shield-exit ratio; independent of absorber material |

Downloaded data for additional materials are qualification candidates, not a
license to apply silicon response functions to them. The parser fails an
unsupported pairing explicitly.

### Case-local material properties

File 07 demonstrates a case-local shield density record. It extends and pins
the exact parent as `Al_natural_G4_11_4@1`, assigns a new stable material ID,
and supplies provenance. The value is deliberately labeled synthetic and must
be replaced by project evidence before science use. Case-local records are
selected by their new `material_id`; parent aliases are never silently
redirected. A density or composition change changes material, geometry, case,
and response hashes.

### Calculation functions

Select the exact names at top level:

```json
"quantities": ["TID", "DDD", "n_eq", "LET", "H100/10"]
```

Keep `damage_response` when either DDD or `n_eq` is selected, and remove it
when neither is selected. Keep `let_response` exactly when LET is selected.
The strict parser rejects stale response blocks so an old hidden setting cannot
leak into a revised case.

## Validate and run

From the repository root, after building as described in `../../RUN.md`:

```bash
source /data/vtenishe/G4/install/bin/geant4.sh

for case_file in examples/configuration_cases/*.json; do
  test_output/build/shieldSim --validate-case="$case_file" || exit 1
done
```

Run one case while keeping generated products under `test_output/`:

```bash
mkdir -p test_output/configuration_cases/02_si_cylinder_neutron_damage
cd test_output/configuration_cases/02_si_cylinder_neutron_damage
../../build/shieldSim \
  --case=../../../examples/configuration_cases/02_si_cylinder_neutron_damage.json
cd ../../..
```

Run every example:

```bash
for case_file in examples/configuration_cases/*.json; do
  case_name=$(basename "$case_file" .json)
  output_dir="test_output/configuration_cases/$case_name"
  mkdir -p "$output_dir"
  (
    cd "$output_dir"
    ../../build/shieldSim --case="../../../$case_file"
  ) || exit 1
done
```

The registered regression gate runs those same seven complete transports and
checks the resolved selectors, output set, case-local density conversion, zero
proton-only `n_eq`, and positive neutron-case `n_eq`:

```bash
cd test_output/build
ctest --output-on-failure -R '^configuration_examples$'
cd ../..
```

The generated inputs are maintained by
`scripts/generate_stage10_shape_examples.py`; edit that generator and review
the regenerated JSON diff when changing the gallery.
