// swcme3d.hpp
//
// ─────────────────────────────────────────────────────────────────────────────
// PUBLIC API: Lightweight solar-wind + CME forward-shock model
// ─────────────────────────────────────────────────────────────────────────────
//
// What this model does (high level):
//   • Prescribes a time-dependent CME-driven shock front that expands through a
//     steady Parker-spiral solar wind. The shock kinematics follow the drag-based
//     model (DBM). Inside the shock we blend a compressed sheath and a CME ejecta
//     (“magnetic cloud”) toward ambient downstream values using C^1 smoothstep
//     transitions with user-controlled widths.
//   • Computes upstream density from the Leblanc et al. (1998) profile, normalized
//     to a chosen n(1 AU). The background magnetic field follows the Parker spiral,
//     normalized to |B|(1 AU) = B1AU_nT.
//   • Approximates the local oblique shock compression ratio (rc) using a fast-mode
//     Mach-number proxy (Edmiston & Kennel 1984; see also Priest 2014). We cap
//     rc ≤ 4 for γ = 5/3 (hydrodynamic limit).
//   • Returns plasma n, V, B fields and a radial-centered estimate of ∇·Vsw,
//     and can export these together with a triangulated shock surface to Tecplot.
//
// Main assumptions / scope:
//   • Kinematics: a single DBM apex; anisotropy via geometry (sphere/ellipsoid/cone).
//   • Flow: purely radial (bulk V aligned with r̂); azimuthal/latitudinal flows
//     are neglected. This keeps evaluation very fast (sub-ms per grid slice).
//   • MHD: rc is an oblique fast-shock proxy based on upstream Parker B, cs, vA.
//     We do not integrate Rankine–Hugoniot relations; we use a minimal closure
//     suited for quick forward-modeling and visualization.
//
// References (suggested starting points):
//   • Parker, E.N. (1958), ApJ 128:664 — Parker spiral.
//   • Leblanc, Y., Dulk, G.A., Bougeret, J.-L. (1998), Solar Phys. 183:165 — density model.
//   • Vršnak, B., et al. (2013), Solar Phys. 285:295 — drag-based model (DBM).
//   • Edmiston, J.P., Kennel, C.F. (1984), J. Plasma Phys. 32:429 — oblique shock speeds.
//   • Priest, E. (2014), Magnetohydrodynamics of the Sun — textbook treatment.
//   • Russell, C.T., Mulligan, T. (2002), PSS 50:527 — CME shocks & sheath properties.
//
// Build example:
//   g++ -std=c++17 -O3 -march=native demo3d_2.cpp swcme3d.cpp -o demo
//
#ifndef SWCME3D_HPP
#define SWCME3D_HPP

#include <cstddef>
#include <vector>

namespace swcme3d {

// ---------- Physical constants (defined in swcme3d.cpp) ----------
extern const double AU;  // Astronomical Unit [m]
extern const double Rs;  // Solar radius [m]
extern const double PI;  // Pi

// ---------- Geometry / shape ----------
enum class ShockShape { Sphere, Ellipsoid, ConeSSE };

// ---------- User parameters ----------
struct Params {
  // Shock geometry
  ShockShape shape = ShockShape::ConeSSE;
  double cme_dir[3] = {1.0, 0.0, 0.0};   // apex direction (unit vector)
  double half_width_rad = 35.0 * PI/180.0; // Cone half-width (for ConeSSE)
  double flank_slowdown_m = 0.0;         // cosine exponent for flank slowdown
  double axis_ratio_y = 1.0;             // Ellipsoid: b/a
  double axis_ratio_z = 1.0;             // Ellipsoid: c/a

  // Kinematics (DBM) — Vršnak et al. (2013)
  // dV/dt = -Γ (V - Vsw) |V - Vsw| → closed-form r(t), V(t)
  double r0_Rs       = 1.05; // initial apex radius [R_sun]
  double V0_sh_kms   = 1500; // initial apex speed [km/s]
  double V_sw_kms    = 400;  // background solar wind speed [km/s]
  double Gamma_kmInv = 0.2e-7; // DBM drag parameter [1/km]

  // Ambient plasma
  double n1AU_cm3 = 5.0;      // Leblanc density normalized to this at 1 AU [cm^-3]
  double T_K      = 1.5e5;    // proton temperature [K] (for cs)
  double gamma_ad = 5.0/3.0;  // adiabatic index

  // Parker spiral parameters (near-ecliptic: sinθ ~ 1)
  double B1AU_nT  = 5.0;      // |B|(1 AU) in nT
  double sin_theta = 1.0;     // sin(colatitude) factor in B_phi term

  // CME radial structure (self-similar scaling with r_sh)
  double sheath_thick_AU_at1AU     = 0.08; // sheath thickness at 1 AU [AU]
  double ejecta_thick_AU_at1AU     = 0.25; // ejecta thickness at 1 AU [AU]
  double edge_smooth_shock_AU_at1AU= 0.01; // shock smoothing half-width [AU]
  double edge_smooth_le_AU_at1AU   = 0.02; // LE smoothing half-width [AU]
  double edge_smooth_te_AU_at1AU   = 0.03; // TE smoothing half-width [AU]

  // Sheath ramp and floors
  double sheath_ramp_power = 1.0;  // ramp exponent (≥1 for sharper shock)
  double sheath_comp_floor = 1.2;  // min compression in sheath

  // Target speeds inside regions (fractions of V_sw)
  double V_sheath_LE_factor = 0.9; // speed at inner edge of sheath (LE)
  double V_ME_factor        = 0.8; // speed in ejecta

  // Ejecta density factor (relative to upstream)
  double f_ME = 0.5;
};

// ---------- Time-dependent state (derived from params for each time step) ----------
struct StepState {
  // Apex-aligned orthonormal basis (e1 along CME direction)
  double e1[3], e2[3], e3[3];

  // Apex shock radius and speed
  double r_sh_m = 0.0;   // [m]
  double V_sh_ms = 0.0;  // [m/s]
  double a_m = 0.0;      // characteristic size (using r_sh)

  // Region thicknesses and smoothing widths
  double dr_sheath_m = 0.0;
  double dr_me_m     = 0.0;
  double w_shock_m   = 0.0;
  double w_le_m      = 0.0;
  double w_te_m      = 0.0;

  // Radii of internal transitions
  double r_le_m = 0.0;   // leading edge (sheath→ejecta)
  double r_te_m = 0.0;   // trailing edge (ejecta→ambient)

  // Target speeds (absolute) for blending
  double V_sheath_LE_ms = 0.0;
  double V_ME_ms        = 0.0;
  double V_dn_ms        = 0.0; // downstream/ambient target

  // Diagnostics
  double rc = 1.0;          // apex compression ratio
  double inv_dr_sheath = 0.0;
  double rc_floor = 1.0;
};

// ---------- Shock surface mesh ----------
struct ShockMesh {
  // Nodal fields
  std::vector<double> x, y, z;          // node positions
  std::vector<double> n_hat_x, n_hat_y, n_hat_z; // nodal outward normals
  std::vector<double> rc;               // nodal compression ratio
  std::vector<double> Vsh_n;            // nodal normal shock speed

  // Connectivity (Tecplot is 1-based)
  std::vector<int> tri_i, tri_j, tri_k; // triangle indices
};

struct TriMetrics {
  // Per-triangle (cell-centered) metrics
  std::vector<double> area;             // [m^2]
  std::vector<double> nx, ny, nz;       // unit normal (from triangle geometry)
  std::vector<double> cx, cy, cz;       // centroid position [m]
  std::vector<double> rc_mean;          // mean of nodal rc at triangle vertices
  std::vector<double> Vsh_n_mean;       // mean of nodal Vsh_n at triangle vertices
};

// ---------- Volume sampling box ----------
struct BoxSpec {
  // Box center and half-sizes (axis-aligned in global coordinates)
  double cx=0, cy=0, cz=0; // [m]
  double hx=0, hy=0, hz=0; // [m]

  // Grid resolution
  int Ni=0, Nj=0, Nk=0;
};

// ---------- Main model ----------
class Model {
public:
  explicit Model(const Params& P);

  // Build state for a given time (seconds since launch epoch)
  StepState prepare_step(double t_s) const;

  // Directional shock radius and outward normal for unit direction u=(ux,uy,uz)
  void shape_radius_normal(const StepState& S,
                           double ux, double uy, double uz,
                           double& Rdir_m, double n_hat[3]) const;

  // Local oblique-shock proxy:
  //  rc_out  : compression ratio (≤4 for γ=5/3)
  //  Vsh_n   : shock propagation speed projected on local normal
  //  thetaBn : obliquity angle between upstream B and shock normal
  void local_oblique_rc(const StepState& S, const double u[3], const double n_hat[3],
                        double Rdir_m, double r_eval_m,
                        double& rc_out, double& Vsh_n_out, double& thetaBn_out) const;

  // Evaluate density and velocity (no magnetic field)
  void evaluate_cartesian_fast(const StepState& S,
                               const double* x_m, const double* y_m, const double* z_m,
                               double* n_m3, double* Vx_ms, double* Vy_ms, double* Vz_ms,
                               std::size_t N) const;

  // Evaluate density, velocity, and magnetic field (Parker upstream with sheath Bt amplification)
  void evaluate_cartesian_with_B(const StepState& S,
                                 const double* x_m, const double* y_m, const double* z_m,
                                 double* n_m3, double* Vx_ms, double* Vy_ms, double* Vz_ms,
                                 double* Bx_T, double* By_T, double* Bz_T,
                                 std::size_t N) const;

  // Convenience: same as above but also returns div(Vsw) at each point
  void evaluate_cartesian_with_B_div(const StepState& S,
                                     const double* x_m, const double* y_m, const double* z_m,
                                     double* n_m3, double* Vx_ms, double* Vy_ms, double* Vz_ms,
                                     double* Bx_T, double* By_T, double* Bz_T,
                                     double* divVsw,
                                     std::size_t N,
                                     double dr_frac = 1e-3) const;

  // Radial centered approximation for ∇·Vsw (works with the radially-aligned flow)
  void compute_divV_radial(const StepState& S,
                           const double* x_m, const double* y_m, const double* z_m,
                           double* divV, std::size_t N, double dr_frac=1e-3) const;

  // Quick diagnostic in a given unit direction u
  void diagnose_direction(const StepState& S, const double u[3],
                          double& Rdir_m, double n_hat[3],
                          double& rc_loc, double& Vsh_n) const;

  // Build a lat–lon triangulated shock surface (with connectivity)
  ShockMesh build_shock_mesh(const StepState& S, std::size_t nTheta, std::size_t nPhi) const;

  // Compute per-triangle metrics (area, centroid, geometric normals, mean rc, mean Vsh_n)
  void compute_triangle_metrics(const ShockMesh& M, TriMetrics& T) const;

  // An outward-shifted cubic box centered near the apex (for volume visualization)
  BoxSpec default_apex_box(const StepState& S, double half_AU, int N) const;

  // FE surface writer (cell-centered metrics + nodal rc, Vsh_n); includes connectivity
  bool write_shock_surface_center_metrics_tecplot(const ShockMesh& M,
                                                  const TriMetrics& T,
                                                  const char* path) const;

  // Full dataset writer (surface FE nodal + FE cell metrics + structured volume).
  // Optionally may include a structured 2-D face zone (controlled in .cpp).
  bool write_tecplot_dataset_bundle(const ShockMesh& M, const TriMetrics& T,
                                    const StepState& S, const BoxSpec& B,
                                    const char* path) const;


  bool write_box_face_minX_tecplot_FE(const StepState& S,
                                    const BoxSpec& B,
                                    int I, int J,
                                    const char* path) const;

  // Structured 2-D min-X box face writer (no connectivity required).
  // Use this to create files like "cone_face_origin_tecplot.dat" safely.
  bool write_box_face_minX_tecplot_structured(const StepState& S,
                                              const BoxSpec& B,
                                              const char* path) const;

private:
  Params P_;
};

} // namespace swcme3d

#endif // SWCME3D_HPP

