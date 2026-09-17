// demo3d_2.cpp — End-to-end demonstration of the solar-wind + CME shock model
// ============================================================================
// What this example shows
// -----------------------
// 1) **How CME “strength” is defined and measured** in this model:
//      • Apex shock compression ratio rc(t)
//      • Fast-mode Mach number from the production ideal-MHD shock state
//      • Shock-frame upstream normal speed and shock normal speed
//      • Actual downstream/upstream |B| ratio and magnetic-field rotation
//      • Normalized Rankine-Hugoniot conservation residuals
// 2) **Time series at 1 AU** (n, V) + apex shock kinematics.
// 3) **Snapshot at t = 36 h**:
//      • Triangulated shock surface (SSE) + per-triangle metrics
//      • A supported apex-centered structured volume box and min-X face
//      • Point sampling of plasma + B (with sheath tangential amplification)
//      • 10 random points from the **first triangle** are printed to stdout
//      • 10 random points **per triangle** written to Tecplot (POINT)
//
// Files written
// -------------
//   strength_summary.csv                 : apex “strength” metrics vs time
//   ts_cone.csv                          : time series at 1 AU (n, V components, |V|)
//   shock_cone.csv                       : apex shock radius/speed + sheath/ejecta edges
//   sse_apex_bundle_tecplot.dat          : four-zone surface/volume/face bundle
//   predefined_points_tecplot.dat        : Tecplot POINT zone (predefined samples)
//   surface_random_samples_tecplot.dat   : Tecplot POINT zone (10 samples per triangle)
//
// Build & run
// -----------
//   g++ -std=c++17 -O3 -march=native demo3d_2.cpp swcme3d.cpp -o demo3d_2
//   ./demo3d_2
//
// Units & conventions
// -------------------
//   Distance: m;  Time: s;  Velocity: m/s;  Density: m^-3;  Magnetic field: Tesla.
//   Inputs ending with *_kms are km/s; n1AU is in cm^-3; B1AU in nT; angles in radians.
//   Sun center is (0,0,0). Parker spiral uses Ω_sun ≈ 2.86533e-6 rad/s.
//
// Physics used (short recap; see swcme3d.hpp/.cpp head comments for equations)
// ----------------------------------------------------------------------------
// • Apex kinematics via Drag-Based Model (DBM):
//     u(t) = (V0 − Vsw) / (1 + Γ (V0 − Vsw) t)
//     r(t) = r0 + Vsw t + (1/Γ) ln(1 + Γ (V0 − Vsw) t)
//     Vsh  = Vsw + u(t)
// • Ambient density: Leblanc et al. (1998) scaled to n(1 AU) = n1AU_cm3.
// • B field: Parker spiral normalized to |B|(1 AU) = B1AU_nT.
// • Sheath/ejecta: shock jump → compressed sheath → depleted ejecta → ambient
//   blended with C^1 smoothsteps (independent edge smoothness).
// • The full ideal-MHD Rankine-Hugoniot solver returns compression, fast Mach
//   number, and conservative downstream rho/p/V/B together with normalized
//   mass, magnetic, electric, momentum, and energy residuals.  This example
//   reports those values directly rather than inverting a gas-dynamic proxy.
//
// Magnetic field evaluator used here
// ----------------------------------
// `evaluate_cartesian_with_B(...)` returns upstream Parker B everywhere but in
// the **sheath** it amplifies the **tangential** component toward rc at the shock,
// relaxing to upstream at the sheath leading edge. This captures first-order
// rotation/strength right behind a forward shock.
//
// ============================================================================

#include "swcme3d.hpp"   // model API (Params, Model, AU, etc.)
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cstdio>
#include <random>
#include <stdexcept>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

using namespace swcme3d;

static inline double hours(double h){ return h * 3600.0; }

// -----------------------------------------------------------------------------
// Utility: sample a point uniformly inside a triangle using barycentric coords
// -----------------------------------------------------------------------------
struct Vec3 { double x,y,z; };

static inline Vec3 bary_sample(const Vec3& A, const Vec3& B, const Vec3& C,
                               std::mt19937_64& rng)
{
  std::uniform_real_distribution<double> U(0.0, 1.0);
  double r1 = U(rng), r2 = U(rng);
  double s  = std::sqrt(r1);
  double wA = 1.0 - s;
  double wB = s * (1.0 - r2);
  double wC = s * r2;
  return Vec3{ wA*A.x + wB*B.x + wC*C.x,
               wA*A.y + wB*B.y + wC*C.y,
               wA*A.z + wB*B.z + wC*C.z };
}

// -----------------------------------------------------------------------------
// Utility: write a Tecplot POINT zone with plasma + B + shock diagnostics
// -----------------------------------------------------------------------------
static swcme::ModelStatus write_points_tecplot_checked(
                                 const char* path,
                                 const std::vector<double>& X,
                                 const std::vector<double>& Y,
                                 const std::vector<double>& Z,
                                 const std::vector<double>& n,
                                 const std::vector<double>& Vx,
                                 const std::vector<double>& Vy,
                                 const std::vector<double>& Vz,
                                 const std::vector<double>& Bx,
                                 const std::vector<double>& By,
                                 const std::vector<double>& Bz,
                                 const std::vector<double>& rc,
                                 const std::vector<double>& Vsh_n)
{
  const std::size_t N = X.size();
  if (!path)
    return swcme::ModelStatus::make(
        swcme::StatusCode::NullPointer,"demo3d_2 point-cloud path");
  if (N==0 || Y.size()!=N || Z.size()!=N || n.size()!=N || Vx.size()!=N ||
      Vy.size()!=N || Vz.size()!=N || Bx.size()!=N || By.size()!=N ||
      Bz.size()!=N || rc.size()!=N || Vsh_n.size()!=N)
    return swcme::ModelStatus::make(
        swcme::StatusCode::InvalidConfiguration,
        "demo3d_2 point-cloud parallel arrays");

  // Validate every precomputed field before creating a transaction.  The
  // point-cloud helper belongs only to this example, but it follows the same
  // no-partial-output rule as the public writers exercised by OUT02-OUT05.
  for (std::size_t i=0; i<N; ++i) {
    const double values[]={X[i],Y[i],Z[i],n[i],Vx[i],Vy[i],Vz[i],Bx[i],By[i],
                           Bz[i],rc[i],Vsh_n[i]};
    for (double value : values)
      if (!std::isfinite(value))
        return swcme::ModelStatus::make_value(
            swcme::StatusCode::NonFiniteResult,
            "demo3d_2 point-cloud value",value,i);
  }

  // CheckedTextFile gives auxiliary demo files the same exact-write,
  // flush/error/close, staging cleanup, and atomic-commit behavior as the
  // model products.  OUT07 then parses the committed bytes independently.
  swcme::output::CheckedTextFile output(
      swcme::output::stdio_file_operations());
  if (!output.open_transactional(path))
    return swcme::ModelStatus::make(
        swcme::StatusCode::FileOpenFailure,"demo3d_2 point-cloud open");
  output.print("demo3d_2 point-cloud title",swcme::ModelStatus::npos,
      "TITLE = \"Point cloud (plasma + B + shock diagnostics)\"\n");
  output.print("demo3d_2 point-cloud variables",swcme::ModelStatus::npos,
      "VARIABLES = \"X[m]\",\"Y[m]\",\"Z[m]\",\"n[m^-3]\","
      "\"Vx[m/s]\",\"Vy[m/s]\",\"Vz[m/s]\",\"Bx[T]\",\"By[T]\","
      "\"Bz[T]\",\"rc[-]\",\"Vsh_n[m/s]\"\n");
  output.print("demo3d_2 point-cloud zone",swcme::ModelStatus::npos,
      "ZONE T=\"points\", N=%zu, F=POINT\n",N);
  for (std::size_t i=0; i<N && output.good(); ++i)
    output.print("demo3d_2 point-cloud row",i,
        "%.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e\n",
        X[i],Y[i],Z[i],n[i],Vx[i],Vy[i],Vz[i],Bx[i],By[i],Bz[i],rc[i],
        Vsh_n[i]);
  return output.finish(
      "demo3d_2 point-cloud flush","demo3d_2 point-cloud stream error",
      "demo3d_2 point-cloud close","demo3d_2 point-cloud commit");
}

// -----------------------------------------------------------------------------
// Helper: write the production ideal-MHD apex state over time.  The residuals
// expose solver quality directly and avoid the historical, invalid inference
// of fast Mach number from a gas-dynamic compression formula.
// -----------------------------------------------------------------------------
static bool write_strength_summary_csv(Model& model,double t0,double t1,
                                       double dt,const char* path)
{
  std::ofstream out(path);
  if (!out.is_open()) return false;
  out << "t_s,has_shock,rc_apex,M_fast,U1n_ms,Vsh_n_ms,theta_Bn_rad,"
         "B2overB1,field_rotation_rad,mass_residual,normal_B_residual,"
         "electric_residual,momentum_residual,energy_residual\n";

  for (double t=t0; t<=t1+1e-9; t+=dt){
    StepState S = model.prepare_step(t);
    double u_apex[3] = { S.e1[0], S.e1[1], S.e1[2] };
    LocalShockState shock;
    const swcme::ModelStatus status=
        model.shock_state_direction_checked(S,u_apex,shock);
    if (!status.ok()) return false;

    const double upstream_normal_velocity=
        shock.upstream.velocity_m_s[0]*shock.normal[0]+
        shock.upstream.velocity_m_s[1]*shock.normal[1]+
        shock.upstream.velocity_m_s[2]*shock.normal[2];
    const double U1n=shock.Vsh_n_m_s-upstream_normal_velocity;
    const auto magnitude=[](const swcme::shock::Vec3& vector) {
      return std::sqrt(vector[0]*vector[0]+vector[1]*vector[1]+
                       vector[2]*vector[2]);
    };
    const double B1=magnitude(shock.upstream.magnetic_T);
    const double B2=magnitude(shock.downstream.magnetic_T);
    const double dot=shock.upstream.magnetic_T[0]*shock.downstream.magnetic_T[0]+
        shock.upstream.magnetic_T[1]*shock.downstream.magnetic_T[1]+
        shock.upstream.magnetic_T[2]*shock.downstream.magnetic_T[2];
    const double cosine=(B1>0.0 && B2>0.0)
        ? std::max(-1.0,std::min(1.0,dot/(B1*B2))) : 1.0;

    out << std::setprecision(10)
        << t << ',' << (shock.has_shock ? 1 : 0) << ',' << shock.compression
        << ',' << shock.fast_mach << ',' << U1n << ',' << shock.Vsh_n_m_s
        << ',' << shock.theta_Bn_rad << ',' << (B1>0.0 ? B2/B1 : 1.0)
        << ',' << std::acos(cosine) << ',' << shock.mass_residual << ','
        << shock.normal_B_residual << ',' << shock.electric_residual << ','
        << shock.momentum_residual << ',' << shock.energy_residual << '\n';
  }
  out.close();
  return static_cast<bool>(out);
}

// ============================================================================
// MAIN
// ============================================================================
int main(){
  try{
    // =========================================================================
    // 1) Define CME/shock parameters (SSE) — with explanatory comments
    // =========================================================================
    Params P;

    // --- Shock front geometry ------------------------------------------------
    P.shape = ShockShape::SSE;           // choices: Sphere, Ellipsoid, SSE
    P.cme_dir[0] = 0.0;                      // apex direction (unit vector)
    P.cme_dir[1] = 0.0;                      // here: +Z (0,0,1)
    P.cme_dir[2] = 1.0;

    // SSE angular half-width (radians). Only directions with theta <= half_width intersect the finite cap.
    P.half_width_rad   = 40.0 * (3.14159265358979323846/180.0); // 40°

    // --- Apex kinematics (DBM) ----------------------------------------------
    // r0_Rs    : DBM reference distance in solar radii; ~15-20 Rs is the recommended science range
    // V0_sh_kms: initial shock apex speed (km/s)
    // V_sw_kms : ambient solar-wind speed (km/s)
    // Gamma_kmInv: drag parameter Γ (1/km). Larger Γ → stronger deceleration/acceleration to V_sw.
    P.r0_Rs       = 20.0;
    P.V0_sh_kms   = 1500.0;
    P.V_sw_kms    = 400.0;
    P.Gamma_kmInv = 0.2e-7;

    // --- Ambient plasma & thermodynamics ------------------------------------
    // n1AU_cm3 rescales Leblanc so that n(1 AU) matches this value (cm^-3).
    P.n1AU_cm3 = 5.0;
    // Adiabatic sound speed uses T_K and gamma_ad.
    P.T_K      = 1.5e5;
    P.gamma_ad = 5.0/3.0;

    // --- Magnetic field (Parker) --------------------------------------------
    // B1AU_nT fixes |B|(1 AU) at the reference colatitude given by
    // sin_theta.  The actual 3-D winding is computed from the local position
    // and the explicit solar-rotation axis.  Choose +Y here so this example's
    // +Z CME/observer direction lies in the solar equatorial plane while the
    // existing visualization geometry remains unchanged.
    P.B1AU_nT   = 5.0;
    P.sin_theta = 1.0;  // reference normalization is equatorial
    P.solar_rotation_axis[0] = 0.0;
    P.solar_rotation_axis[1] = 1.0;
    P.solar_rotation_axis[2] = 0.0;

    // --- Sheath / magnetic ejecta parameterization --------------------------
    // Thickness values are defined at 1 AU and scale self-similarly with the current apex radius.
    P.sheath_thick_AU_at1AU = 0.10;   // thicker sheath → longer relaxation
    P.ejecta_thick_AU_at1AU = 0.20;   // ME (flux-rope) radial thickness

    // C1 transition widths (AU at a 1-AU shock).  The shock width is active
    // only in RESOLVED_COMPRESSION; SOURCE mode uses SHOCK_ONLY and no velocity
    // compression layer.  Keep the resolved shock layer relatively narrow.
    P.edge_smooth_shock_AU_at1AU = 0.01;
    P.edge_smooth_le_AU_at1AU    = 0.03;
    P.edge_smooth_te_AU_at1AU    = 0.03;

    // sheath_comp_floor is deprecated and ignored; the MHD RH solution alone
    // sets the physical compression.  sheath_ramp_power shapes post-shock relaxation.
    P.sheath_comp_floor = 1.5;  // compatibility-only; has no physical effect
    P.sheath_ramp_power = 1.5;  // >1 → faster drop from rc at the shock

    // Target speeds inside sheath/ejecta relative to ambient.
    P.V_sheath_LE_factor = 1.05; // sheath relaxes toward, but not below, ambient
    P.V_ME_factor        = 0.8; // ejecta typically slower

    // Ejecta density factor (depleted compared to upstream).
    P.f_ME = 0.6;

    // Build the model with these parameters
    Model model(P);

    // =========================================================================
    // 2) “Strength” summary over time (apex metrics): 0–72 h, Δt=5 min
    // =========================================================================
    {
      const double t0=0.0, t1=hours(72.0), dt=300.0;
      if (!write_strength_summary_csv(
              model,t0,t1,dt,"strength_summary.csv"))
        throw std::runtime_error("failed to write strength_summary.csv");

      // Also print a quick apex strength readout at t=36 h (for convenience)
      StepState S = model.prepare_step(hours(36.0));
      double u_apex[3] = { S.e1[0], S.e1[1], S.e1[2] };
      LocalShockState shock;
      const swcme::ModelStatus shock_status=
          model.shock_state_direction_checked(S,u_apex,shock);
      if (!shock_status.ok())
        throw std::runtime_error("36-hour apex shock state failed: "+
                                 shock_status.summary());

      std::cout << "[Strength @ apex, t=36h] "
                << "rc=" << shock.compression
                << ", M_fast=" << shock.fast_mach
                << ", Vsh_n=" << shock.Vsh_n_m_s << " m/s"
                << ", theta_Bn=" << shock.theta_Bn_rad << " rad"
                << ", max_RH_residual="
                << std::max({std::abs(shock.mass_residual),
                             std::abs(shock.normal_B_residual),
                             std::abs(shock.electric_residual),
                             std::abs(shock.momentum_residual),
                             std::abs(shock.energy_residual)}) << '\n';
    }

    // =========================================================================
    // 3) Time series at 1 AU (+Z). Observer at (0,0,1 AU)
    // =========================================================================
    {
      const double t0=0.0, t1=hours(72.0), dt=300.0;
      const std::size_t Nt = static_cast<std::size_t>((t1-t0)/dt)+1;

      double x_obs[1]={0.0}, y_obs[1]={0.0}, z_obs[1]={AU};

      std::ofstream ts("ts_cone.csv");    ts << "t_s,n_m3,Vx_ms,Vy_ms,Vz_ms,V_mag_ms\n";
      std::ofstream sh("shock_cone.csv"); sh << "t_s,R_sh_AU,V_sh_km_s,rc,R_LE_AU,R_TE_AU\n";
      if (!ts.is_open() || !sh.is_open())
        throw std::runtime_error("cannot open demo3d_2 time-series outputs");

      for (std::size_t k=0;k<Nt;++k){
        // The explicit cast documents the intentional discrete-to-continuous
        // conversion and keeps the demonstration inside OUT08's warning gate.
        const double t=t0+static_cast<double>(k)*dt;
        StepState S = model.prepare_step(t);

        double n,Vx,Vy,Vz;
        model.evaluate_cartesian_fast(S, x_obs,y_obs,z_obs, &n,&Vx,&Vy,&Vz, 1);
        const double Vmag = std::sqrt(Vx*Vx+Vy*Vy+Vz*Vz);

        ts << std::setprecision(10) << t << ','
           << n << ',' << Vx << ',' << Vy << ',' << Vz << ',' << Vmag << '\n';

        sh << std::setprecision(10) << t << ','
           << (S.r_sh_m/AU) << ',' << (S.V_sh_ms/1e3) << ',' << S.rc << ','
           << (S.r_le_m/AU) << ',' << (S.r_te_m/AU) << '\n';
      }
      // CSV streams can report delayed errors only at flush/close.  Checking
      // both here prevents OUT07 from accepting a nominally successful demo
      // whose diagnostics were truncated.
      ts.close();
      sh.close();
      if (!ts || !sh)
        throw std::runtime_error("failed to finish demo3d_2 time-series outputs");
    }

    // =========================================================================
    // 4) Snapshot at t = 36 h: surface triangulation + volume box + bundle
    // =========================================================================
    const double t_mesh = hours(36.0);
    StepState Smesh = model.prepare_step(t_mesh);

    // Build the surface (unique-apex, periodic-ring triangulation of the finite SSE cap) and per-triangle metrics
    // Demonstration-grade resolution exercises finite-cap topology and random
    // sampling without the historical multi-minute, hundreds-of-megabytes run.
    ShockMesh  surf = model.build_shock_mesh(Smesh, /*nTheta=*/24, /*nPhi=*/48);
    TriMetrics tri;  model.compute_triangle_metrics(surf, tri);

    // Construct a small apex-centered grid through the validated factory.  The
    // old box touched the solar origin and was outside the model domain, so its
    // writer failed after the demo had already generated expensive samples.
    const BoxSpec B=model.default_apex_box(Smesh,/*half_AU=*/0.02,/*N=*/12);

    const swcme::ModelStatus bundle_status=
        model.write_tecplot_dataset_bundle_checked(
            surf,tri,Smesh,B,"sse_apex_bundle_tecplot.dat");
    if (!bundle_status.ok())
      throw std::runtime_error("bundle output failed: "+bundle_status.summary());

    // =========================================================================
    // 5) Evaluate plasma + B at a set of predefined points (for quick sanity)
    // =========================================================================
    {
      std::vector<Vec3> Pts = {
        {0.0,        0.0,        0.80*AU},   // inside 1 AU, sunward of shock
        {0.0,        0.0,        0.95*AU},
        {0.0,        0.0,        1.05*AU},   // just beyond 1 AU
        {0.02*AU,    0.00,       1.00*AU},   // slight off-axis
        {-0.02*AU,   0.00,       1.00*AU},
        {0.00,       0.02*AU,    1.00*AU}
      };

      const std::size_t Np = Pts.size();
      std::vector<double> X(Np),Y(Np),Z(Np), n(Np),Vx(Np),Vy(Np),Vz(Np), Bx(Np),By(Np),Bz(Np), rc(Np),Vsh_n(Np);

      for (std::size_t i=0;i<Np;++i){
        X[i]=Pts[i].x; Y[i]=Pts[i].y; Z[i]=Pts[i].z;

        // Plasma + B (Tesla) with sheath tangential amplification
        model.evaluate_cartesian_with_B(Smesh, &X[i],&Y[i],&Z[i],
                                        &n[i],&Vx[i],&Vy[i],&Vz[i],
                                        &Bx[i],&By[i],&Bz[i], 1);

        // Local shock diagnostics (along the LOS direction of this point)
        const double r = std::sqrt(X[i]*X[i]+Y[i]*Y[i]+Z[i]*Z[i]);
        const double invr = (r>0)? 1.0/r : 0.0;
        double udir[3]={X[i]*invr, Y[i]*invr, Z[i]*invr};
        double Rdir, n_hat[3], rc_loc, Vsh_n_loc;
        model.diagnose_direction(Smesh, udir, Rdir, n_hat, rc_loc, Vsh_n_loc);
        rc[i]=rc_loc; Vsh_n[i]=Vsh_n_loc;
      }

      const swcme::ModelStatus point_status=write_points_tecplot_checked(
          "predefined_points_tecplot.dat",X,Y,Z,n,Vx,Vy,Vz,Bx,By,Bz,rc,Vsh_n);
      if (!point_status.ok())
        throw std::runtime_error("predefined point output failed: "+
                                 point_status.summary());
    }

    // =========================================================================
    // 6) Surface sampling: (a) print 10 random points in the **first triangle**
    // =========================================================================
    if (!surf.tri_i.empty()){
      std::mt19937_64 rng(42);
      // Generated mesh connectivity is validated, positive, and one-based.
      // Convert it once to the vector index type instead of allowing implicit
      // signed-to-unsigned conversions at every coordinate access.
      const std::size_t ia=static_cast<std::size_t>(surf.tri_i[0]-1);
      const std::size_t ib=static_cast<std::size_t>(surf.tri_j[0]-1);
      const std::size_t ic=static_cast<std::size_t>(surf.tri_k[0]-1);
      Vec3 A{surf.x[ia], surf.y[ia], surf.z[ia]};
      Vec3 Bv{surf.x[ib], surf.y[ib], surf.z[ib]};
      Vec3 Cv{surf.x[ic], surf.y[ic], surf.z[ic]};

      std::cout << "Ten random points in first triangle (A="<<ia<<",B="<<ib<<",C="<<ic<<"):\n";
      std::cout << "idx, x, y, z [m], n[m^-3], Vx, Vy, Vz[m/s], Bx, By, Bz[T], rc, Vsh_n[m/s]\n";

      for (int m=0;m<10;++m){
        Vec3 Pp = bary_sample(A,Bv,Cv, rng);

        double n_,Vx_,Vy_,Vz_,Bx_,By_,Bz_;
        model.evaluate_cartesian_with_B(Smesh, &Pp.x,&Pp.y,&Pp.z,
                                        &n_,&Vx_,&Vy_,&Vz_,
                                        &Bx_,&By_,&Bz_, 1);

        const double rP = std::sqrt(Pp.x*Pp.x+Pp.y*Pp.y+Pp.z*Pp.z);
        const double invrP = (rP>0)?1.0/rP:0.0;
        double udir[3]={Pp.x*invrP, Pp.y*invrP, Pp.z*invrP};
        double Rdir, n_hat[3], rc_loc, Vsh_n_loc;
        model.diagnose_direction(Smesh, udir, Rdir, n_hat, rc_loc, Vsh_n_loc);

        std::cout << m << ", "
                  << std::setprecision(10) << Pp.x << ", " << Pp.y << ", " << Pp.z << ", "
                  << n_ << ", " << Vx_ << ", " << Vy_ << ", " << Vz_ << ", "
                  << Bx_ << ", " << By_ << ", " << Bz_ << ", "
                  << rc_loc << ", " << Vsh_n_loc << "\n";
      }
    }

    // =========================================================================
    // 7) Surface sampling: (b) 10 random points **per triangle** → Tecplot POINT zone
    // =========================================================================
    {
      std::mt19937_64 rng(2025);
      const std::size_t Ne = surf.tri_i.size();
      const int SAMPLES_PER_TRI = 10;
      const std::size_t Ns = Ne * SAMPLES_PER_TRI;

      std::vector<double> SX, SY, SZ, Sn, SVx, SVy, SVz, SBx, SBy, SBz, Src, SVshn;
      SX.reserve(Ns); SY.reserve(Ns); SZ.reserve(Ns);
      Sn.reserve(Ns); SVx.reserve(Ns); SVy.reserve(Ns); SVz.reserve(Ns);
      SBx.reserve(Ns); SBy.reserve(Ns); SBz.reserve(Ns);
      Src.reserve(Ns); SVshn.reserve(Ns);

      for (std::size_t e=0;e<Ne;++e){
        const std::size_t i=static_cast<std::size_t>(surf.tri_i[e]-1);
        const std::size_t j=static_cast<std::size_t>(surf.tri_j[e]-1);
        const std::size_t k=static_cast<std::size_t>(surf.tri_k[e]-1);
        Vec3 A{surf.x[i], surf.y[i], surf.z[i]};
        Vec3 Bv{surf.x[j], surf.y[j], surf.z[j]};
        Vec3 Cv{surf.x[k], surf.y[k], surf.z[k]};

        for (int m=0;m<SAMPLES_PER_TRI;++m){
          Vec3 Pp = bary_sample(A,Bv,Cv, rng);

          double n_,Vx_,Vy_,Vz_,Bx_,By_,Bz_;
          model.evaluate_cartesian_with_B(Smesh, &Pp.x,&Pp.y,&Pp.z,
                                          &n_,&Vx_,&Vy_,&Vz_,
                                          &Bx_,&By_,&Bz_, 1);

          const double rP = std::sqrt(Pp.x*Pp.x+Pp.y*Pp.y+Pp.z*Pp.z);
          const double invrP = (rP>0)?1.0/rP:0.0;
          double udir[3]={Pp.x*invrP, Pp.y*invrP, Pp.z*invrP};
          double Rdir, n_hat[3], rc_loc, Vsh_n_loc;
          model.diagnose_direction(Smesh, udir, Rdir, n_hat, rc_loc, Vsh_n_loc);

          SX.push_back(Pp.x);  SY.push_back(Pp.y);  SZ.push_back(Pp.z);
          Sn.push_back(n_);    SVx.push_back(Vx_);  SVy.push_back(Vy_);  SVz.push_back(Vz_);
          SBx.push_back(Bx_);  SBy.push_back(By_);  SBz.push_back(Bz_);
          Src.push_back(rc_loc); SVshn.push_back(Vsh_n_loc);
        }
      }

      const swcme::ModelStatus sample_status=write_points_tecplot_checked(
          "surface_random_samples_tecplot.dat",SX,SY,SZ,Sn,SVx,SVy,SVz,
          SBx,SBy,SBz,Src,SVshn);
      if (!sample_status.ok())
        throw std::runtime_error("surface sample output failed: "+
                                 sample_status.summary());
    }

    std::cout << "Done. Wrote strength_summary.csv, CSV time series, and Tecplot datasets.\n";
    return 0;

  } catch (const std::exception& e){
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }
}
