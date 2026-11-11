# CME model-local registry adapters

These adapters execute model-owned registries and translate only their typed
result records into the global schema. They do not duplicate physics expected
values. Every local registry is validated in full before a global `--case`
selection is applied, so selecting one case cannot conceal local registry
drift.

The staged coronal-CME adapter imports 222 rows as **software verification**;
its three campaign-evidence cases remain typed `SKIP`, and the local
`production_qualified=false` and `observational_campaign_qualified=false`
claims are retained. The full-volume Corona--SWCME adapter imports 46 bounded
code markers. BG3D-4/G1 implementation resumed on 2026-11-10, but remains
unqualified. The reduced
adapter imports 97 scope-qualified rows. Routine global execution uses the
model runner's `--skip-native` mode, producing 62 executed portable rows and
35 explicit native `SKIP` rows rather than silently dropping native scope.

All Makefile-driven sources are copied below the suite work directory before
building. This excludes stale checkout objects and supports read-only source
mounts. The authoritative local JSON/JUnit/text products are retained below
the global artifact directory.

Run the registered suites with:

```sh
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-MODEL-SEP-CORONAL-CME-001
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-MODEL-SEP-CORONA-SWCME-FULL-001
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-MODEL-SEP-CORONA-SWCME-REDUCED-001
```

A complete native reduced-front campaign remains a separate model-local run:

```sh
python3 srcSEP3D/test/run_reduced_shock_front.py --output-dir /absolute/evidence/path
```

That campaign is zero-particle and reduced-front-plus-ambient only. Neither a
geometric observer arrival nor its PASS summary asserts an accepted fast shock
or a spatial sheath/contact/ejecta volume.

## Current analytical cross-ray development

The user retained the 3-D ambient and selected conservative analytical/
semi-empirical lateral material transport without an MHD runtime, a threaded
ideal material contact, and a 1200-s startup-corotation taper. The full-volume
owner now has 46 rows: its original 15 plus 31 new model checks (29 cross-ray/
material/contact cases and two numerical controls). It remains unqualified as
a production spatial CME provider. The complete core campaign preserves all
792 old PASS identities and passes all 32 new global cases; its sole remaining
FAIL/ERROR are strict incomplete coverage and the original refinement timeout.

Use [the owner README](sep_corona_swcme_full/README.md) for exact execution,
selection/concurrency, timeout and failure/report rules; [the model test README](../../../../models/sep_corona_swcme/test/README.md)
for quick binaries; and [the API README](../../../../models/sep_corona_swcme/README_CROSS_RAY.md)
for equations, SI/frame conventions, immutable state and physical limitations.
The individual descriptions below supply what/how/why/where without becoming
a second physics oracle.

## Numerical controls, fluxes and static geometry

| Global case | What / how / why / where |
| --- | --- |
| `CSWC-FULL-CMBGU05-DERIVATIVE-CONTROLS` | Uses independent expanding-sphere and translated-plane trajectories to verify the actual CMBGU05 fourth-order position derivative and dimensional normal-flux calculation. Tests fourth-order convergence, both signs of imposed slip, noncancelling absolute leakage, wrong normals, nonfinite data and unsupported endpoints. Prevents a numerical derivative repair from bypassing physical leakage rejection. Source: `src/models/sep_corona_swcme/test/test_bg3d4_contact_derivative.cpp`. |
| `CSWC-FULL-CMBGU05-SHELL-GEOMETRY-CONTROLS` | Independently integrates exact spherical volume and first moments with two-node Gaussian quadrature, including broad, thin and one-ulp shells from small SI scales to 1 AU; rejects invalid/unrepresentable geometry. Validates the factored production moments without changing the finite-volume equations or physical acceptance limits. Source: `src/models/sep_corona_swcme/test/test_spherical_shell.cpp`. |
| `CSWC-FULL-CSWC0628-XRAY-FLUX` | Hand-evaluated SI ideal-MHD fluxes, zero-field acoustic limits, parallel/perpendicular and degenerate fast characteristics, primitive/conservative energy and both upwind branches. Validates the new shared cross-ray kernel before coupling. Source: `src/models/sep_corona_swcme/test/test_cross_ray_transport.cpp`. |
| `CSWC-FULL-CSWC0628-XRAY-ALE` | Independent moving-face Reynolds flux values, normal-orientation reversal and Cartesian rotation. Prevents missing pressure/Maxwell energy work or confusing moving induction flux with a stationary face. Same source. |
| `CSWC-FULL-CSWC0628-XRAY-CONSERVATION` | Equal-and-opposite shared-face transfer on unequal cell volumes conserves mass, all Cartesian momentum components, total energy and volume-integrated B; checks actual extensive SI ledgers. Same source. |
| `CSWC-FULL-CSWC0628-XRAY-REJECTION` | Rejects invalid density/pressure/normal/gamma, inconsistent normal B, unresolved thermal energy, unqualified timestep and unimplemented moving-volume update; immutable inputs remain unchanged. Same source. |
| `CSWC-FULL-CSWC0628-XRAY-ADVECTION` | Advances an analytic lateral periodic density contact at 32/64/128 cells through actual shared flux/transfer routines, checks first-order convergence, positivity and global conserved inventories. Same source. |
| `CSWC-FULL-CSWC0628-XRAY-GEOMETRY` | Analytic octant volume/centroid/signed face vectors, HCI rotation and uniform MHD flux closure on curved static faces. Source: `src/models/sep_corona_swcme/test/test_spherical_sector.cpp`. |
| `CSWC-FULL-CSWC0628-XRAY-MESH-CLOSURE` | Exact whole-sphere shell volume/area and all 792 interfaces of a 336-cell explicit sector reference mesh. Same source. |
| `CSWC-FULL-CSWC0628-XRAY-GEOMETRY-REJECTION` | One-ulp heliocentric shells remain positive; invalid radii, angles and left-handed/nonorthonormal frames reject without clipping. Same source. |
| `CSWC-FULL-CSWC0628-XRAY-TOPOLOGY` | Unique shared nodes/faces/edges, exact integer boundary-of-boundary cancellation, polar/seam handling and the analytic shell Euler characteristic. Source: `src/models/sep_corona_swcme/test/test_cross_ray_mesh.cpp`. |
| `CSWC-FULL-CSWC0628-XRAY-EDGE-GEOMETRY` | Actual curved edge samples meet canonical vertices with correct radial, meridional and latitude arc-length line elements; unsupported queries reject. Same source. |
| `CSWC-FULL-CSWC0628-XRAY-MESH-IDENTITY` | Immutable event/geometry/epoch/generation identity, changed-geometry fingerprints, cell-budget and periodicity rejection, and transactional failed preparation. Same source. |

The reference helper forms antisymmetric position differences before large
weighted sums. Its formula, steps, physical bound, and the original CMBGU05
convergence condition are unchanged. The API's reported zero flux is never the
independent reference. The refinement owner retains its 180-second execution
limit and all five production advances. These software checks do not qualify
the resumed, still-unqualified full-volume physical closure.


## Static constrained-transport reference

| Global case | What / how / why / where |
| --- | --- |
| `CSWC-FULL-CSWC0628-XRAY-CT-INITIALIZATION` | Independent Gaussian integration of a known vector potential on actual curved edges; compares its circulation against analytic constant-field face flux and retains raw discrete divergence. Source: `src/models/sep_corona_swcme/test/test_cross_ray_ct.cpp`. |
| `CSWC-FULL-CSWC0628-XRAY-CT-FARADAY` | Manufactured Faraday evolution followed by twenty common random edge circulations on the full shared mesh; checks every face flux and divergence with explicit arithmetic bounds. Same source. |
| `CSWC-FULL-CSWC0628-XRAY-CT-AUTHORITY` | Rejects initial monopoles, nonfinite or incomplete edge data, foreign identity and stale epochs/generations; verifies immutable committed states and zero-EMF publication. These static CT tests do not qualify fluid-energy or moving-region coupling. Same source. |

## Field reconstruction and atomic fluid-energy references

| Global case | What / how / why / where |
| --- | --- |
| `CSWC-FULL-CSWC0628-XRAY-RECONSTRUCTION-CONSTANT` | Recovers known Cartesian fields on curved faces under HCI rotation and magnetic SI scaling without changing CT face flux. Source: `src/models/sep_corona_swcme/test/test_cross_ray_reconstruction.cpp`. |
| `CSWC-FULL-CSWC0628-XRAY-RECONSTRUCTION-CONVERGENCE` | Uses independent edge-potential integration of a divergence-free affine field; checks volume-average field error over 4/8/16 cells per axis and retains nonzero fit residuals. Same source. |
| `CSWC-FULL-CSWC0628-XRAY-RECONSTRUCTION-REJECTION` | Rejects a geometrically rank-deficient half-circle cell without regularization; accepts zero magnetic field on resolvable cells. Same source. |
| `CSWC-FULL-CSWC0628-XRAY-FLUID-UNIFORM` | Twelve static curved-mesh epochs with hand-integrated uniform ideal-MHD flux and electric circulation preserve all fluid/field values and common generation. Source: `src/models/sep_corona_swcme/test/test_cross_ray_fluid.cpp`. |
| `CSWC-FULL-CSWC0628-XRAY-FLUID-CONSERVATION` | Checks analytic opposite cross-ray transfers, all Cartesian momentum components and energy, untouched cells, and the outward sign of inner-boundary mass/work. Same source. |
| `CSWC-FULL-CSWC0628-XRAY-FLUID-ENERGY` | Prescribes a manufactured external curl-E, keeps total energy conservative with the new CT field, and separately checks explicitly supplied electromagnetic source work. This is kernel validation, not ideal-MHD edge-solver qualification. Same source. |
| `CSWC-FULL-CSWC0628-XRAY-FLUID-REJECTION` | Incomplete face/source data, stale CT authority, nonfinite transverse flux, invalid EOS and nonpositive thermal energy reject the whole epoch; no energy reset or input mutation. Same source. |

## Corotating geometry and threaded-contact diagnostics

| Global case | What / how / why / where |
| --- | --- |
| `CSWC-FULL-CSWC0628-XRAY-CONTACT-STARTUP` | Independently checks the selected 1200-s quintic angular transition and actual four-ray ambient/contact startup leakage while retaining threaded magnetic flux. Source: `src/models/sep_corona_swcme/test/test_cross_ray_contact.cpp`. |
| `CSWC-FULL-CSWC0628-XRAY-CONTACT-DERIVATIVES` | Differentiates actual material-point positions and fixed-HCI ray intersections independently through startup, handoff and late propagation; checks normal speed, radial acceleration and implicit geometry. Same source. |
| `CSWC-FULL-CSWC0628-XRAY-CONTACT-INTERFACE` | Manufactured density-jump material contact, signed leakage, oblique-frame tangential electric mismatch, Maxwell traction and moving-interface pressure work; continuous-field zero-Bn limit. Same source. |
| `CSWC-FULL-CSWC0628-XRAY-CONTACT-REJECTION` | Explicit duration and owning contact fingerprint, unsupported epochs, invalid geometry/normal, and immutable failed queries. These checks do not qualify moving-mesh or actual two-sided ejecta coupling. Same source. |

## Analytical Cauchy material transport

| Global case | What / how / why / where |
| --- | --- |
| `CSWC-FULL-CSWC0628-XRAY-MATERIAL-INVARIANTS` | Independent affine triangular map, deformed face cross products and signed flux/mass/entropy references verify full Cauchy transport rather than radial field scaling. Source: `src/models/sep_corona_swcme/test/test_material_transport.cpp`. |
| `CSWC-FULL-CSWC0628-XRAY-MATERIAL-RATES` | Independently differentiates a closed-form transported field and Eulerian velocity, checking material density/pressure/induction rates and full velocity-gradient tensor. Same source. |
| `CSWC-FULL-CSWC0628-XRAY-MATERIAL-NONLINEAR` | Nonlinear spatially varying lateral deformation checks global mass and divergence with explicit inverse-map sampling at 24 points. No production algorithm generates expected fields. Same source. |
| `CSWC-FULL-CSWC0628-XRAY-MATERIAL-REJECTION` | Rejects folds, unresolved rank, unrepresentable Jacobian and nonfinite map rates; distinguishes resolved small-J kernel behavior from the separate frozen physical minimum-J gate. Same source. |
