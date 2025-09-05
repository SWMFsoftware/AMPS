// swcme3d.hpp
//
// ─────────────────────────────────────────────────────────────────────────────
// SOLAR WIND + CME FORWARD-SHOCK SEMI-ANALYTICAL MODEL — PUBLIC API
// ─────────────────────────────────────────────────────────────────────────────
//
// What this header declares
// -------------------------
// • Public constants (AU, Rs, PI) with external linkage.
// • Parameter bundle (Params) with physically-named fields and units.
// • Time-dependent state (StepState) computed from Params at time t.
// • Simple geometric containers: ShockMesh, TriMetrics, BoxSpec.
// • Model class with:
//   - Directional shock geometry: radius and normal.
//   - Oblique MHD proxy for compression ratio rc and normal shock speed Vsh_n.
//   - Field evaluators returning density (n), bulk velocity (V), magnetic
//     field (B) and divergence ∇·V at arbitrary Cartesian points.
//   - Shock surface meshing + per-triangle metrics (area, centroid, normal,
//     rc_mean, Vsh_n_mean).
//   - Tecplot writers (surface, volume box, an optional box face).
//
// Units & conventions
// -------------------
// • Positions in SI meters (m). Velocities in m/s. Density in m^-3.
// • B field in Tesla (1 nT = 1e-9 T).
// • Angles in radians. Time in seconds.
// • The Sun is at the origin (0,0,0).
// • The CME apex direction is the unit vector Params::cme_dir (global frame).
//
// Physics model (brief)
// ---------------------
// • Ambient density n_up(r): Leblanc et al. (1998), scaled to n(1 AU).
// • Parker spiral upstream B (Parker 1958), with |B|(1 AU) set by B1AU_nT.
// • CME apex kinematics via Drag-Based Model (DBM): Vršnak et al. (2013).
// • Shock shape: sphere, ellipsoid, or cone-like SSE (with flank slowdown).
// • Local rc proxy from oblique fast Mach number (Edmiston & Kennel 1984;
//   Priest 2014), limited to ≤ 4.
// • Sheath / ejecta blends using C^1 smoothsteps and independent edge widths.
// • Tangential B amplified smoothly in the sheath (Bn continuous).
// • ∇·V by robust radial finite-difference: (1/r^2) ∂(r^2 V_r)/∂r.
//
// Output (Tecplot)
// ----------------
// VARIABLES order (used consistently across all zones):
//   1: X [m]   2: Y [m]   3: Z [m]
//   4: n [m^-3]      5..7: Vx,Vy,Vz [m/s]
//   8..10: Bx,By,Bz [T]   11: divVsw [1/s]
//   12: rc [-]        13: Vsh_n [m/s]
//   14..16: nx,ny,nz (geom normal; cell-centered in surface_cells)
//   17: area [m^2]    18: rc_mean [-]  19: Vsh_n_mean [m/s]
//   20..22: tnx,tny,tnz (reserved)     23..25: cx,cy,cz [m] (centroid)
//
// References (short list)
// -----------------------
// • Parker, E.N. (1958) ApJ 128:664 — Parker spiral.
// • Leblanc, Y.; Dulk, G.A.; Bougeret, J.-L. (1998) Sol. Phys. 183:165 — n(r).
// • Vršnak, B. et al. (2013) Sol. Phys. 285:295 — DBM CME kinematics.
// • Edmiston, J.P.; Kennel, C.F. (1984) J. Plasma Phys. 32:429 — oblique shocks.
// • Priest, E. (2014) Magnetohydrodynamics of the Sun, CUP.
// • Russell, C.T.; Mulligan, T. (2002) PSS 50:527 — CME sheath/ICME structure.
// • Manchester, W.B. IV et al. (2005) ApJ 622:1225 — CME-driven shocks.
//
#pragma once
#include <cstddef>
#include <vector>

namespace swcme3d {

// ─────────────────────────────────────────────────────────────────────────────
// Public constants (external linkage; defined in swcme3d.cpp)
// ─────────────────────────────────────────────────────────────────────────────
extern const double AU;   // Astronomical Unit [m]
extern const double Rs;   // Solar radius [m]
extern const double PI;   // π

// ─────────────────────────────────────────────────────────────────────────────
// Shock shape type
// ─────────────────────────────────────────────────────────────────────────────
enum class ShockShape { Sphere, Ellipsoid, ConeSSE };

// ─────────────────────────────────────────────────────────────────────────────
// Parameter pack (user inputs)
// ─────────────────────────────────────────────────────────────────────────────
struct Params {
  // --- CME geometry & kinematics ---
  ShockShape shape = ShockShape::ConeSSE; // default: cone-like
  double cme_dir[3] = {1,0,0};            // unit vector (global) apex direction
  double axis_ratio_y = 1.0;              // Ellipsoid only: b/a
  double axis_ratio_z = 1.0;              // Ellipsoid only: c/a
  double half_width_rad = 45.0 * (PI/180.0); // Cone: half-angle
  double flank_slowdown_m = 1.0;          // Cone: radius∝cos^m(theta)

  // DBM apex initial conditions & drag
  double r0_Rs = 1.05;     // initial shock apex radius in solar radii
  double V0_sh_kms = 1500; // initial shock speed [km/s]
  double V_sw_kms  = 400;  // ambient solar wind speed [km/s]
  double Gamma_kmInv = 0.5e-7; // drag parameter Γ [km^-1] (typical ~1e-7)

  // Ambient plasma & magnetic field at 1 AU (for scaling)
  double n1AU_cm3 = 5.0;   // number density at 1 AU [cm^-3]
  double B1AU_nT  = 5.0;   // |B| at 1 AU [nT]
  double sin_theta = 1.0;  // sin(colatitude) for Parker pitch (≈ sin heliographic latitude)

  // Thermodynamics (for cs, fast speed, Mach)
  double gamma_ad = 5.0/3.0; // adiabatic index
  double T_K      = 1.5e5;   // effective temperature [K] for sound speed

  // Sheath & ejecta thicknesses (self-similar, relative to r_sh)
  double sheath_thick_AU_at1AU  = 0.08; // sheath thickness at 1 AU
  double ejecta_thick_AU_at1AU  = 0.15; // ejecta (ICME) thickness at 1 AU

  // Edge smoothing widths (C^1), scaled ~ r_sh
  double edge_smooth_shock_AU_at1AU = 0.01; // shock ramp width
  double edge_smooth_le_AU_at1AU    = 0.02; // leading edge (sheath→ejecta)
  double edge_smooth_te_AU_at1AU    = 0.03; // trailing edge (ejecta→ambient)

  // Target velocities within regions (factors of V_sw)
  double V_sheath_LE_factor = 1.15; // sheath LE magnitude ~ 1.15 V_sw
  double V_ME_factor        = 0.9;  // ejecta magnitude ~ 0.9 V_sw

  // Ejecta density fraction (simple cavity or enhancement)
  double f_ME = 0.5; // ejecta density = f_ME * n_up

  // Sheath compression shaping
  double sheath_ramp_power = 3.0; // steeper front as power increases
  double sheath_comp_floor = 1.2; // minimum compression at shock onset
};

// ─────────────────────────────────────────────────────────────────────────────
// Time-dependent state (built once per time step by prepare_step)
// ─────────────────────────────────────────────────────────────────────────────
struct StepState {
  // Apex-aligned orthonormal frame (e1 = apex direction)
  double e1[3], e2[3], e3[3];

  // Shock radii and speeds
  double r_sh_m;      // apex shock radius [m]
  double V_sh_ms;     // apex shock speed [m/s]
  double a_m;         // semimajor axis (for sphere/ellipsoid; here = r_sh_m)

  // Region geometry (self-similar thicknesses measured radially)
  double dr_sheath_m; // sheath thickness [m]
  double dr_me_m;     // ejecta thickness [m]
  double w_shock_m;   // smoothing width at shock
  double w_le_m;      // smoothing width at leading edge
  double w_te_m;      // smoothing width at trailing edge
  double r_le_m;      // radius of sheath→ejecta transition
  double r_te_m;      // radius of ejecta→ambient transition

  // Region target speeds
  double V_sheath_LE_ms; // sheath speed at leading edge [m/s]
  double V_ME_ms;        // ejecta speed [m/s]
  double V_dn_ms;        // downstream ambient speed [m/s]

  // Diagnostics
  double rc;             // apex rc (diagnostic only)
  double inv_dr_sheath;  // 1 / dr_sheath (if positive)
  double rc_floor;       // from Params::sheath_comp_floor
};

// ─────────────────────────────────────────────────────────────────────────────
// Shock surface mesh (nodal attributes + triangle list)
// ─────────────────────────────────────────────────────────────────────────────
struct ShockMesh {
  // Nodal positions
  std::vector<double> x, y, z;         // size Nv

  // Nodal outward normals and nodal shock diagnostics
  std::vector<double> n_hat_x, n_hat_y, n_hat_z; // nodal unit normals
  std::vector<double> rc;              // nodal compression ratio
  std::vector<double> Vsh_n;           // nodal normal shock speed [m/s]

  // Connectivity (1-based indices for Tecplot FE)
  std::vector<int> tri_i, tri_j, tri_k; // size 3*Ne
};

// ─────────────────────────────────────────────────────────────────────────────
// Per-triangle (cell-centered) metrics
// ─────────────────────────────────────────────────────────────────────────────
struct TriMetrics {
  std::vector<double> area;      // [m^2]
  std::vector<double> nx, ny, nz;// geometric normal (unit)
  std::vector<double> cx, cy, cz;// centroid [m]
  std::vector<double> rc_mean;   // mean of nodal rc
  std::vector<double> Vsh_n_mean;// mean of nodal Vsh_n [m/s]
};

// ─────────────────────────────────────────────────────────────────────────────
// Structured volume box spec
// ─────────────────────────────────────────────────────────────────────────────
struct BoxSpec {
  // Center and half-sizes (axis-aligned in global frame)
  double cx=0, cy=0, cz=0;
  double hx=0, hy=0, hz=0;
  // Grid resolution (I×J×K)
  int Ni=16, Nj=16, Nk=16;
};

// ─────────────────────────────────────────────────────────────────────────────
// Model (stateless w.r.t. time; per-step state is StepState)
// ─────────────────────────────────────────────────────────────────────────────
class Model {
public:
  explicit Model(const Params& P);

  // Build all time-dependent quantities at time t [s]
  StepState prepare_step(double t_s) const;

  // Directional shock radius and outward normal for unit direction u=(ux,uy,uz)
  void shape_radius_normal(const StepState& S,
                           double ux,double uy,double uz,
                           double& Rdir_m,double n_hat[3]) const;

  // Oblique-shock proxy: compression ratio rc and normal shock speed Vsh_n
  // Evaluated using upstream quantities at r_eval_m.
  void local_oblique_rc(const StepState& S,
                        const double u[3], const double n_hat[3],
                        double Rdir_m, double r_eval_m,
                        double& rc_out, double& Vsh_n_out, double& thetaBn_out) const;

  // Evaluate n and radial V (returned as Cartesian V) at N points
  void evaluate_cartesian_fast(const StepState& S,
                               const double* x_m,const double* y_m,const double* z_m,
                               double* n_m3,double* Vx_ms,double* Vy_ms,double* Vz_ms,
                               std::size_t N) const;

  // Evaluate n, V, Parker-based B with sheath Bt amplification
  void evaluate_cartesian_with_B(const StepState& S,
                                 const double* x_m,const double* y_m,const double* z_m,
                                 double* n_m3,double* Vx_ms,double* Vy_ms,double* Vz_ms,
                                 double* Bx_T,double* By_T,double* Bz_T,
                                 std::size_t N) const;

  // Same as above + ∇·V (radial finite-difference)
  void evaluate_cartesian_with_B_div(const StepState& S,
                                     const double* x_m,const double* y_m,const double* z_m,
                                     double* n_m3,double* Vx_ms,double* Vy_ms,double* Vz_ms,
                                     double* Bx_T,double* By_T,double* Bz_T,double* divVsw,
                                     std::size_t N, double dr_frac=1e-3) const;

  // Radial ∇·V helper (public so callers can recompute with different dr_frac)
  void compute_divV_radial(const StepState& S,
                           const double* x_m,const double* y_m,const double* z_m,
                           double* divV,std::size_t N,double dr_frac) const;

  // Quick diagnostic along a given unit direction u
  void diagnose_direction(const StepState& S,
                          const double u[3],
                          double& Rdir_m,double n_hat[3],
                          double& rc_loc,double& Vsh_n) const;

  // Build triangulation of shock surface (θ×φ sampling)
  ShockMesh build_shock_mesh(const StepState& S,
                             std::size_t nTheta, std::size_t nPhi) const;

  // Per-triangle metrics
  void compute_triangle_metrics(const ShockMesh& M, TriMetrics& T) const;

  // A convenience box near the apex (shifted outward along e1)
  BoxSpec default_apex_box(const StepState& S,double half_AU,int N) const;

  // Tecplot writers
  bool write_shock_surface_center_metrics_tecplot(const ShockMesh& M,
                                                  const TriMetrics& T,
                                                  const char* path) const;

  bool write_tecplot_dataset_bundle(const ShockMesh& M,const TriMetrics& T,
                                    const StepState& S,const BoxSpec& B,
                                    const char* path) const;

  bool write_box_face_minX_tecplot_structured(const StepState& S,
                                              const BoxSpec& B,
                                              const char* path) const;

private:
  Params P_;
};

} // namespace swcme3d

