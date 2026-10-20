# HEAL corruption of the gravity, slope and cliff products
      
      Status: defect found 2026-10-19 during the CUDA C00 reference capture
      (`../LiouvilleCOPS_cuda/GRAVITY_REFERENCE.md` §5.1). This document describes the
      problem and a fix for `tools/build_physical_activity_templates.py`.

      Update 2026-10-19: the repair is now implemented in both gravity paths,
      with independent topology validation and a new cache algorithm/orientation
      identity. The five known corrupted production files were preserved with
      `_heal_defect` appended to their filenames. The template selector also
      refuses cliff columns without repaired gravity provenance. See
      `GRAVITY_REPAIR.md` for tests and remaining production-validation gates.
      Corrected full-field production products still need regeneration; the
      historical evidence and proposed acceptance checks below are retained.
      
      ## 1. Summary
     
      `build_physical_activity_templates.py` builds the gravity polyhedron with
      `integrity_check=pg.PolyhedronIntegrity.HEAL` and `metric_unit=pg.MetricUnit.METER`.
      For a kilometre-scale body given in metres, HEAL's orientation test misjudges
      about 46% of correctly oriented faces as inward and flips them. The flipped
      polyhedron is what gets evaluated, so the self-gravity, effective gravity, slopes
      and cliff masks are wrong.
      
      The mesh itself is correct. The SHAP5 788k BDF is a closed, consistently
      outward-oriented 2-manifold: every edge is shared by exactly two faces and every
      directed edge occurs once. Its Euler characteristic is 2, its signed volume is
      +1.876e10 m³, and no face has zero area. Nothing needed healing; HEAL damaged it.

      ## 2. Mechanism (polyhedral-gravity 3.3)
      
      `Polyhedron::checkPlaneUnitNormalOrientation` / `countRayPolyhedronIntersections`
      (`Polyhedron.cpp`):
  
      1. For each face, a ray starts at `centroid + 1e-14 * N` and points along the
         plane normal `N`.
      2. It is intersected with every face (Möller–Trumbore, eps = 1e-14, hit if
         `t > 1e-14`). The distinct intersection points are counted, and an **odd**
         count marks the face as inward.
      3. Majority rule: if at most n/2 faces are flagged, the orientation is OUTWARDS.
         HEAL then swaps vertices 0 and 1 of every flagged face.
         
      The defect: `EPSILON_ZERO_OFFSET = 1e-14` is **absolute, in mesh units**. At
      coordinates of about 2,600 m, one rounding unit of a double is about 4.5e-13 m,
      far larger than the 1e-14 m offset. `centroid + 1e-14*N` therefore rounds back
      onto the face's own plane. Möller–Trumbore then reports the starting face itself
      at `t ≈ 1e-13 > 1e-14` for about half of the faces, depending on rounding. That
      makes the count odd, and a correct face is flagged and flipped. The flags are
      effectively random, so they cannot be undone after the fact.

      ## 3. Evidence
      
      | Mesh | Independent topology | Package result | Effect |
      |---|---|---|---|
      | icosphere L3, r = 1, 10, 100 m | outward | OUTWARDS, 0 flagged | correct |
      | icosphere L3, r = 1000 m | outward | INWARDS, 612/1280 flagged | gravity error up to 121% |
      | icosphere L3, r = 2000 m | outward | INWARDS, 588/1280 flagged | gravity error up to 127% |
      | SHAP5 788k BDF | closed, outward, Euler 2 | exact replica: 138/300 sampled faces flagged (46%), all own-face hits | about 362k faces flipped |
      
      On 48 evenly spaced 788k facets (0.1 m offset), flipping 46% of faces at random gives:
      - |g| error: median 92%;
      - slope error: median 81°, p95 144°;
      - facets ≥ 60°: 38/48 with flips, against 2/48 with correct gravity.
      
      The finished 788k template build of 2026-10-19
      (`physical_templates_shap5_788k.json`) shows the same signature:
      - slope median 71.96°, maximum 179.77° (a slope above 90° means gravity points
        out of the surface);
      - `cliff_ge_45deg` area fraction 0.710 and `cliff_ge_60deg` area fraction 0.595;
      - |g| from 1.4e-5 to 6.1e-5 m/s². Correct gravity is about 1.4e-4 to 2.3e-4 m/s²
        with a median slope of about 18°.

      Note that the build's own offset-halving audit **passed** (0.108° < 0.25°). The
      audit only compares two evaluations of the same corrupted polyhedron, so it
      cannot detect a wrong orientation.

      References:
      - C00 capture: `cops_output/cuda_reference/c00_20261019T064536Z/gravity/`
        (`gravity_fixtures.json`, `gravity_788k.npz`: correct DISABLE gravity at 832
        targets plus 32 on-surface points, and the replica flags).
      - C03 GPU evaluation: `cops_output/cuda_c03/gravity_20261019T171225Z/` (`report.json`;
        `full_field.npz` holds the full 787,754-facet DISABLE field, sample-validated
        against the package at ≤ 2.3e-13 relative).

 
      ## 4. Affected products
     
      Invalid, to be quarantined or regenerated:
      - `cops_output/physical_templates_shap5_788k.gravity.npz` (self-gravity cache);
      - `cops_output/physical_templates_shap5_788k.npz`: the `cliff_ge_45deg` and
        `cliff_ge_60deg` columns;
      - `cops_output/physical_templates_shap5_788k.png` and the gravity/slope fields of
        `physical_templates_shap5_788k.json`;
      - `physical_templates_shap5_788k.gravity_100k_rejected.npz` and the "rejected
        100k" decision. Its 5° resolution-gate failure was probably this artifact,
        because `validate_gravity_resolution` also uses HEAL;
      - any fit that used the cliff templates.

      Unaffected: the `late_ng_not_rg` template (geometry only, no gravity), and all
      fits without cliff templates.


      ### 5.1 Code changes in `tools/build_physical_activity_templates.py`
     
      1. Add a topology audit, i.e. exact combinatorics plus one volume sum, which
         does not depend on coordinate scale:
     
      ```python
      def mesh_topology(vertices, faces):
          """Closedness and orientation facts without ray casting."""
          v = np.asarray(vertices, dtype=float)
          f = np.asarray(faces, dtype=np.int64)
          directed = np.concatenate((f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]))
          _, und_counts = np.unique(np.sort(directed, axis=1), axis=0, return_counts=True)
          _, dir_counts = np.unique(directed, axis=0, return_counts=True)
          a, b, c = v[f[:, 0]], v[f[:, 1]], v[f[:, 2]]
          volume = float(np.sum(np.einsum("ij,ij->i", a, np.cross(b, c))) / 6.0)
          area = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
          return {
              "closed_consistent": bool(np.all(und_counts == 2) and np.all(dir_counts == 1)),
              "euler_characteristic": int(len(v) - len(und_counts) + len(f)),
              "signed_volume_m3": volume,
              "zero_area_faces": int(np.count_nonzero(area == 0.0)),
              "unreferenced_vertices": int(len(v) - len(np.unique(f))),
          }

      def require_outward_closed_body(vertices, faces, label):
          t = mesh_topology(vertices, faces)
          problems = []
          if not t["closed_consistent"]:
              problems.append("not a closed, consistently oriented manifold")
          if t["euler_characteristic"] != 2:
              problems.append("Euler characteristic %d" % t["euler_characteristic"])
          if t["zero_area_faces"]:
              problems.append("%d zero-area faces" % t["zero_area_faces"])
          if not t["signed_volume_m3"] > 0.0:
              problems.append("winding is inward (volume %.6g)" % t["signed_volume_m3"])
          if t["unreferenced_vertices"]:
              problems.append("%d unreferenced vertices" % t["unreferenced_vertices"])
          if problems:
              raise ValueError("%s refused: %s" % (label, "; ".join(problems)))
          return t
      ```


      2. In `effective_slope` **and** `validate_gravity_resolution`, call the audit
         before building the polyhedron, and replace
              
         `integrity_check=pg.PolyhedronIntegrity.HEAL`
         
         with 
         
         `integrity_check=pg.PolyhedronIntegrity.DISABLE`.
         
         Do not "repair" a mesh that fails the audit; refuse it. A consistently
         *inward* mesh could be reversed explicitly (`faces[:, [1, 0, 2]]`), but only
         as a documented, opt-in step.
         
      3. Record the audit result and `"integrity_check": "DISABLE after topology audit"`
         in the returned metadata and in the bundle JSON.


      ### 5.2 Invalidate the corrupted cache (essential) 
         
      The gravity-cache definition (`gravity_cache.make_definition`) contains the
      algorithm, geometry identity, shape SHA-256, density and offset, but **not** the
      integrity mode. Without a definition change, the fixed builder would find the
      existing HEAL cache, judge it valid, and silently reuse the corrupted
      self-gravity. Do both of the following:
      - change `gravity_cache.ALGORITHM` to a new value, for example
        `"homogeneous_polyhedron_self_gravity_at_outward_offset_v2_disable_after_topology_audit"`,
        or add an `"orientation_policy"` field to the definition, so old caches no
        longer match;
      - rename the existing products with a `_heal_defect` suffix rather than
        deleting them, so they stay available as evidence.


      ### 5.3 Add tests (`tests/test_activity_templates.py` or a new file)
        
      - **Regression:** an icosphere L3 of radius 1000 m in metres. HEAL flags about
        600 of 1280 faces; DISABLE after the audit gives g within 1e-12 of the same
        mesh scaled to 1 m (rescaled by r).
      - **Refusal:** a cube with one face flipped is refused by the audit, as is an
        open mesh (one face removed).
      - **Reference agreement:** on the C00 targets (`gravity_788k.npz`,
        `points_offset`), the fixed builder's self-gravity equals
        `acceleration_offset`. It should be bitwise identical, since it is the same
        binary, DISABLE and the same face order.
      
      ## 6. Alternatives considered (not recommended)
        
      - **HEAL on coordinates in km** (`metric_unit=KILOMETER`). At about 2.6 km one
        rounding unit is about 4.4e-16, below 1e-14, so the ray test would work. But it
        keeps an O(n²) ray test (6.2e11 tests for 788k faces), adds a unit conversion
        (outputs in km/s², which must be multiplied by 1000), and still depends on a
        scale-sensitive heuristic. The topology audit is exact and fast.
      - **Signed volume only.** It detects a globally inward mesh but not a mixed
        winding. The edge-count audit is needed for consistency.
      - **Fixing the package.** That is upstream work, out of scope here.

      ## 7. Rebuild and acceptance
     
      1. Apply 5.1–5.3. Run the build with 36 threads; expect several hours on CPU,
         since the previous 788k evaluation took 13,826 s.
      2. Acceptance checks:
         - the offset-halving audit max is about 0.0050° (C00 value 0.00502388°);
         - plausibility: |g| about 1.4e-4–2.3e-4 m/s², median slope about 18°, no
           slope > 90° except at genuine overhangs, which should be reviewed;
         - the self-gravity agrees with the C00 fixture at its 832 targets, and with
           the C03 GPU full field (`full_field.npz`) to ≤ 1e-12 relative at every facet;
         - cliff masks: report the area fractions for 45° and 60°, and list the facets
           within 1e-6° of a threshold. 
      3. Re-run the 100k resolution gate with the fixed code before relying on the
         "rejected" decision.
      4. Re-run any fit that used cliff templates.
      
      ## 8. Physical assumptions unchanged by the fix
      
      These are unchanged by the fix:
      - homogeneous density 533 kg/m³;
      - targets at centroid + 0.1 m·outward normal, centrifugal term at the centroid;
      - SPICE omega in the body frame;
      - G = 6.67430e-11;
      - cliff thresholds 45° and 60°.

      Only the polyhedron orientation handling changes.
