// demo3d_1.cpp — supported-domain swcme3d SSE output example
//
// Build:
//   g++ -std=c++17 -O3 -march=native demo3d_1.cpp swcme3d.cpp -o demo3d_1
//
// Run:
//   ./demo3d_1
//
// Outputs:
//   • Time series at 1 AU (CSV):        ts_cone.csv
//   • Shock kinematics (CSV):           shock_cone.csv
//   • Tecplot dataset (4 zones):        sse_apex_bundle_tecplot.dat
//
// Notes:
//   - Densities are in m^-3; velocities in m/s.
//   - The Tecplot file has the production order: surface_cells
//     (FETRIANGLE/BLOCK), surface_nodal (FEPOINT), volume_box
//     (structured POINT), and box_face_minX (structured POINT).
//   - The volume is constructed by default_apex_box() around the 6-hour shock
//     apex.  Unlike the historical origin-crossing box, every requested point
//     is inside the analytical model domain and the checked writer can finish.

#include "swcme3d.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <string>

using namespace swcme3d;

static inline double hours(double h){ return h * 3600.0; }

int main(){
  try {
    // -----------------------------
    // Configure a finite SSE shock
    // -----------------------------
    Params P;
    P.shape = ShockShape::SSE;

    // Apex direction along +Z (the symmetry axis of the finite SSE cap)
    P.cme_dir[0] = 0.0;
    P.cme_dir[1] = 0.0;
    P.cme_dir[2] = 1.0;

    // The production Parker field uses the actual angle to the solar-rotation
    // axis.  Choose +Y so the existing +Z CME/observer direction is explicitly
    // equatorial rather than accidentally lying on the rotation pole.
    P.solar_rotation_axis[0] = 0.0;
    P.solar_rotation_axis[1] = 1.0;
    P.solar_rotation_axis[2] = 0.0;
    P.sin_theta = 1.0;  // B1AU_nT reference normalization at the equator

    // Finite SSE opening: the ray at this half width is tangent to the spherical cap.
    P.half_width_rad   = 40.0 * (3.14159265358979323846/180.0);

    // (Other Params defaults in swcme3d.hpp control DBM, smoothing, etc.)

    // Build model
    Model model(P);

    // -------------------------------------------------------
    // Time series at 1 AU on +Z axis (observer at x=y=0, z=1 AU)
    // -------------------------------------------------------
    // Five-minute output remains well resolved for this usage example while
    // keeping OUT07's full demonstration run compact and deterministic.
    const double t0 = 0.0, t1 = hours(72.0), dt = 300.0;
    const std::size_t Nt = static_cast<std::size_t>((t1 - t0)/dt) + 1;

    double x_obs[1] = { 0.0 };
    double y_obs[1] = { 0.0 };
    double z_obs[1] = { AU   };

    std::ofstream ts("ts_cone.csv");
    ts << "t_s,n_m3,Vx_ms,Vy_ms,Vz_ms,V_mag_ms\n";

    std::ofstream sh("shock_cone.csv");
    sh << "t_s,R_sh_AU,V_sh_km_s,rc,R_LE_AU,R_TE_AU\n";

    // OUT07 runs examples in an empty directory and requires every declared
    // product.  Fail immediately if either CSV cannot be opened rather than
    // allowing the later Tecplot file to make the demo appear successful.
    if (!ts.is_open() || !sh.is_open())
      throw std::runtime_error("cannot open demo3d_1 CSV outputs");

    for (std::size_t k=0; k<Nt; ++k){
      const double t = t0 + k*dt;
      StepState S = model.prepare_step(t);

      double n,Vx,Vy,Vz;
      model.evaluate_cartesian_fast(S, x_obs, y_obs, z_obs, &n,&Vx,&Vy,&Vz, 1);
      const double Vmag = std::sqrt(Vx*Vx + Vy*Vy + Vz*Vz);

      ts << std::setprecision(10) << t << ','
         << std::setprecision(10) << n << ','
         << std::setprecision(10) << Vx << ','
         << std::setprecision(10) << Vy << ','
         << std::setprecision(10) << Vz << ','
         << std::setprecision(10) << Vmag << '\n';

      sh << std::setprecision(10) << t << ','
         << std::setprecision(10) << (S.r_sh_m/AU) << ','
         << std::setprecision(10) << (S.V_sh_ms/1e3) << ','
         << std::setprecision(10) << S.rc << ','
         << std::setprecision(10) << (S.r_le_m/AU) << ','
         << std::setprecision(10) << (S.r_te_m/AU) << '\n';
    }

    // Explicit close/status checks turn delayed filesystem failures into a
    // nonzero demo exit, matching the production Tecplot writer's lifecycle
    // contract and making OUT07 meaningful on full or interrupted storage.
    ts.close();
    sh.close();
    if (!ts || !sh)
      throw std::runtime_error("failed to finish demo3d_1 CSV outputs");

    // -------------------------------------------------------
    // Visualization snapshot at 6 h:
    // - build shock surface
    // - compute cell metrics
    // - define a supported apex-centered volume and its min-X face
    // -------------------------------------------------------
    const double t_mesh = hours(6.0);
    StepState Smesh = model.prepare_step(t_mesh);

    // Demonstration-scale resolution keeps the executable and OUT07 fast while
    // retaining the unique apex, periodic rings, and finite SSE boundary.
    ShockMesh surf = model.build_shock_mesh(Smesh, /*nTheta=*/24, /*nPhi=*/48);

    TriMetrics tri;
    model.compute_triangle_metrics(surf, tri);

    // default_apex_box() applies the same OUT04 structural rules as the writer
    // and places this compact diagnostic grid away from the unsupported solar
    // interior.  Twelve points per dimension keep the file useful but small.
    const BoxSpec B=model.default_apex_box(Smesh,/*half_AU=*/0.02,/*N=*/12);

    // Use the checked API so model-domain, mesh, and I/O failures are visible
    // in stderr and force a nonzero exit instead of the old false-success path.
    const swcme::ModelStatus output_status=
        model.write_tecplot_dataset_bundle_checked(
            surf,tri,Smesh,B,"sse_apex_bundle_tecplot.dat");
    if (!output_status.ok()) {
      std::cerr << "Dataset output failed: " << output_status.summary() << '\n';
      return 1;
    }

    std::cout << "Done. Wrote ts_cone.csv, shock_cone.csv, "
                 "sse_apex_bundle_tecplot.dat\n";
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }
}
