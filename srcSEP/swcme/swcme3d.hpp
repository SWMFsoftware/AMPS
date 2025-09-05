#pragma once
/*
  swcme3d.hpp — Solar-wind + CME kinematic/shock model (1.05 Rsun → interplanetary)

  # Purpose
  Lightweight forward model for solar wind plasma and a CME-driven shock, with
  simple geometry (sphere / ellipsoid / cone-cap SSE), drag-based kinematics (DBM),
  ambient Parker-spiral magnetic field, and a parameterized sheath/ejecta profile
  that blends smoothly from upstream → sheath → ejecta → ambient.

  # Outputs
  - Plasma number density n [m^-3], bulk velocity vector V [m/s]
  - Magnetic field vector B [Tesla] (Parker upstream; tangentially amplified in sheath)
  - Shock diagnostics (local radius, normal, compression ratio rc, normal shock speed Vsh_n)
  - Surface triangulation + per-triangle metrics (area, normal, rc_mean, Vsh_n_mean)
  - Structured-box (Tecplot POINT) volume with n, V, B, and divV (= ∇·V)

  # Core ingredients and equations
  (1) Ambient density — Leblanc et al. (1998) empirical profile:
      n_e(r) [cm^-3] = 3.3e5 (Rs/r)^2 + 4.1e6 (Rs/r)^4 + 8.0e7 (Rs/r)^6
      Here we scale the coefficients by a single factor so that n(1 AU) = n1AU_cm3 (user param).
      Output density n = n_e * 1e6 [m^-3].

  (2) Ambient magnetic field — Parker spiral (Parker 1958):
      Br ∝ r^-2 ; Bphi = - Br (Ω r sinθ / Vsw)
      We set Br(1 AU) so that |B|(1 AU) = B1AU_nT (user param), then construct
      B = Br r̂ + Bphi ϕ̂ with the rotation axis along +z.

  (3) Kinematics — Drag-Based Model (DBM; Vršnak et al. 2013):
      dV/dt = - Γ (V - Vsw) |V - Vsw|  →  V(t) = Vsw + (V0 - Vsw) / (1 + Γ (V0 - Vsw) t)
      r(t) = r0 + Vsw t + ln(1 + Γ (V0 - Vsw) t) / Γ
      Γ in [km^-1] is input as Gamma_kmInv; we convert to [m^-1].

  (4) Shock geometry (choose one):
      • Sphere:               R(û) = r_sh
      • Ellipsoid:            (x/a)^2 + (y/b)^2 + (z/c)^2 = 1, apex along ê1 = CME direction
      • Cone-cap (SSE-like):  R(θ) = r_sh cos^m θ for θ ≤ half_width; flattened flanks via m.

  (5) Sheath/ejecta parameterization (phenomenological):
      - Behind the shock we ramp the compression ratio from a floor (≥1) up to local rc at
        the shock, with a tunable power p and width w_sh (smoothstep C^1).
      - Through the sheath into the ejecta (ME), we blend to target values (n, V) at the
        leading edge (width w_LE) and relax at the trailing edge (width w_TE).
      - The ejecta number-density factor f_ME and speed factor V_ME_factor set targets.

      Blending uses the C^1 smoothstep s(x)=x^2(3−2x) on normalized offsets. Each interface
      uses its own width so you can choose a sharper shock and a softer ejecta tail.

  (6) MHD shock proxy (oblique; Priest 2014; Edmiston & Kennel 1984):
      We estimate the (fast-mode) normal Mach number using upstream sound & Alfven speeds,
      and obliquity θ_Bn from the Parker field vs. the local surface normal. Then:
        M_fn = U1n / c_f  ,  c_f^2 = ½[(v_A^2 + c_s^2) + sqrt((v_A^2 + c_s^2)^2
                                                      − 4 c_s^2 v_A^2 cos^2 θ_Bn)]
        rc = ((γ+1)M_fn^2) / ((γ−1)M_fn^2 + 2)   (capped at 4 for γ=5/3)
      This rc is used as the **target** compression at the shock. In the sheath interior,
      the tangential B is amplified toward rc and then relaxes toward ambient across the
      sheath thickness.

  (7) Divergence of V — efficient radial estimate at a point:
      For a locally radial flow along r̂ with scalar speed V(r) (and slow angular variation),
      ∇·(V r̂) = (1/r^2) d/dr (r^2 V).
      We evaluate V at r±δr along the same direction r̂ and use a centered difference on
      r^2 V, i.e.
        divV ≈ {[(r+δr)^2 V(r+δr) − (r−δr)^2 V(r−δr)] / (2 δr)} / r^2
      with δr = max(dr_frac·r, δr_min). This captures the jump/smoothing across the shock.

  # Units
    • distance [m], time [s], velocity [m/s], density [m^-3], magnetic field [Tesla]

  # References (by topic)
    Parker spiral:           Parker, E.N. (1958), ApJ 128:664.
    Solar rotation rate:     Howard et al. (1984), ARA&A 22:131 (for Ω; we use 2.865e-6 rad/s).
    Density (quiet wind):    Leblanc, Dulk & Bougeret (1998), Sol. Phys. 183:165.
    DBM (CME kinematics):    Vršnak et al. (2013), Sol. Phys. 285:295.
    MHD shocks (oblique):    Edmiston & Kennel (1984), J. Plasma Phys. 32:429;
                             Priest (2014), "Magnetohydrodynamics of the Sun".
    CME sheath phenomenology: Russell & Mulligan (2002), Planet. Space Sci. 50:527;
                             Manchester et al. (2005), ApJ 622:1225.

  -------------------------------------------------------------------------------
  This header defines the public API. See swcme3d.cpp for implementation details.
*/

#include <cstddef>
#include <vector>

namespace swcme3d {

// ------------------------- physical constants (ODR-defined in .cpp) ----------
extern const double AU;   // 1.495978707e11 m
extern const double Rs;   // 6.957e8 m
extern const double PI;   // 3.14159...

// ------------------------- configuration & state structures -------------------
enum class ShockShape { Sphere, Ellipsoid, ConeSSE };

struct Params {
  // Geometry
  ShockShape shape = ShockShape::Sphere;
  double axis_ratio_y = 1.0;   // ellipsoid b/a
  double axis_ratio_z = 1.0;   // ellipsoid c/a
  double half_width_rad = 50.0 * (PI/180.0);  // cone half-angle
  double flank_slowdown_m = 2.0;              // cone slowdown exponent m

  // Kinematics (DBM)
  double r0_Rs = 1.05;       // starting apex radius in Rs
  double V0_sh_kms = 1500.0; // initial apex shock speed [km/s]
  double V_sw_kms  = 400.0;  // ambient solar-wind speed [km/s]
  double Gamma_kmInv = 0.2e-7; // DBM drag [km^-1]

  // Orientation (ê1 = CME apex direction)
  double cme_dir[3] = {1.0, 0.0, 0.0};

  // Ambient density scaling (Leblanc → n(1 AU)=n1AU_cm3)
  double n1AU_cm3 = 5.0;

  // Thermodynamics for sound speed and MHD shock proxy
  double T_K = 1.5e5;        // proton temperature [K]
  double gamma_ad = 5.0/3.0; // adiabatic index

  // Parker spiral parameters
  double sin_theta = 1.0;    // effective sin(latitude) of footpoint (≈1 near equator)
  double B1AU_nT   = 5.0;    // magnitude |B|(1 AU) in nT

  // Sheath / ejecta thickness and smooth widths (self-similar ∝ r_sh)
  double sheath_thick_AU_at1AU = 0.06;
  double ejecta_thick_AU_at1AU = 0.20;
  double edge_smooth_shock_AU_at1AU = 0.01; // sharp near the shock
  double edge_smooth_le_AU_at1AU    = 0.03; // softer at leading edge into ME
  double edge_smooth_te_AU_at1AU    = 0.05; // softest at trailing edge

  // Sheath ramp and floor for compression
  double sheath_ramp_power = 2.0;
  double sheath_comp_floor = 1.2;

  // Target speeds inside sheath and ME (fractions of V_sw)
  double V_sheath_LE_factor = 1.2; // at sheath leading edge
  double V_ME_factor        = 0.9; // inside ejecta

  // Ejecta density factor relative to ambient upstream at same r
  double f_ME = 0.6;
};

struct StepState {
  // Local orthonormal basis aligned with apex direction
  double e1[3]{1,0,0}, e2[3]{0,1,0}, e3[3]{0,0,1};

  // Apex kinematics
  double r_sh_m = 0.0;   // apex shock radius [m]
  double V_sh_ms = 0.0;  // apex shock speed [m/s]
  double a_m = 0.0;      // alias for semi-axis a

  // Self-similar thicknesses and edge widths at current scale
  double dr_sheath_m = 0.0, dr_me_m = 0.0;
  double w_shock_m = 0.0, w_le_m = 0.0, w_te_m = 0.0;

  // Convenience radii along apex
  double r_le_m = 0.0, r_te_m = 0.0;

  // Target speeds
  double V_sheath_LE_ms = 0.0; // speed at sheath leading edge
  double V_ME_ms        = 0.0; // speed in ejecta (ME)
  double V_dn_ms        = 0.0; // downstream far

  // Apex compression ratio (approximate diagnostic)
  double rc = 1.0;

  // Other cached helpers
  double inv_dr_sheath = 0.0;
  double rc_floor = 1.0;
};

struct ShockMesh {
  // Nodal positions
  std::vector<double> x, y, z;
  // Nodal surface normals (unit, outward)
  std::vector<double> n_hat_x, n_hat_y, n_hat_z;
  // Nodal shock properties
  std::vector<double> rc;     // compression ratio
  std::vector<double> Vsh_n;  // normal shock speed [m/s]
  // Connectivity (1-based Tecplot indexing)
  std::vector<int> tri_i, tri_j, tri_k;
};

struct TriMetrics {
  // Cell-centered (per triangle)
  std::vector<double> area;        // [m^2]
  std::vector<double> nx, ny, nz;  // triangle unit normals
  std::vector<double> cx, cy, cz;  // centroids
  std::vector<double> rc_mean;     // mean nodal rc per triangle
  std::vector<double> Vsh_n_mean;  // mean nodal Vsh_n per triangle
};

struct BoxSpec {
  // Box center and half-sizes (axis-aligned in global coordinates)
  double cx=0, cy=0, cz=0;
  double hx=0, hy=0, hz=0;
  // Structured sampling counts
  int Ni=41, Nj=41, Nk=41;
};

// ------------------------- model class ----------------------------------------
class Model {
public:
  explicit Model(const Params&);

  // Time-step preparation (build basis, DBM apex position/speed, thicknesses, etc.)
  StepState prepare_step(double t_s) const;

  // Directional shock radius R(û) and outward normal n̂ for the selected shape
  void shape_radius_normal(const StepState& S,
                           double ux, double uy, double uz,
                           double& Rdir_m, double n_hat[3]) const;

  // Local oblique-shock proxy at direction û (returns rc, Vsh_n, θ_Bn)
  void local_oblique_rc(const StepState& S, const double u[3], const double n_hat[3],
                        double Rdir_m, double r_eval_m,
                        double& rc_out, double& Vsh_n_out, double& thetaBn_out) const;

  // Evaluate n and V at Cartesian points (no magnetic field)
  void evaluate_cartesian_fast(const StepState& S,
                               const double* x_m, const double* y_m, const double* z_m,
                               double* n_m3, double* Vx_ms, double* Vy_ms, double* Vz_ms,
                               std::size_t N) const;

  // Evaluate n, V, B **and** div(Vsw) at Cartesian points.
  // This is a convenience wrapper around evaluate_cartesian_with_B + compute_divV_radial.
  void evaluate_cartesian_with_B_div(const StepState& S,
                                     const double* x_m, const double* y_m, const double* z_m,
                                     double* n_m3, double* Vx_ms, double* Vy_ms, double* Vz_ms,
                                     double* Bx_T, double* By_T, double* Bz_T,
                                     double* divVsw,            // <- new output
                                     std::size_t N,
                                     double dr_frac = 1e-3) const;

  // Evaluate n, V, and B (with shock-tangential amplification in the sheath)
  void evaluate_cartesian_with_B(const StepState& S,
                                 const double* x_m, const double* y_m, const double* z_m,
                                 double* n_m3, double* Vx_ms, double* Vy_ms, double* Vz_ms,
                                 double* Bx_T, double* By_T, double* Bz_T,
                                 std::size_t N) const;

  // Radial divergence estimate at points (centered 1D derivative along local r̂)
  // dr_frac ∈ [~1e-4 .. 1e-2] controls δr = max(dr_frac·r, δr_min)
  void compute_divV_radial(const StepState& S,
                           const double* x_m, const double* y_m, const double* z_m,
                           double* divV, std::size_t N,
                           double dr_frac = 1e-3) const;

  // Build a lat–lon triangulation of the visible cap; returns nodal arrays + 1-based tris
  ShockMesh build_shock_mesh(const StepState& S,
                             std::size_t nTheta, std::size_t nPhi) const;

  // Per-triangle metrics from a ShockMesh
  void compute_triangle_metrics(const ShockMesh& M, TriMetrics& T) const;

  // Default apex-centered box (slightly shifted outward along ê1)
  BoxSpec default_apex_box(const StepState& S, double half_AU, int N) const;

  // Tecplot writers
  // 1) Surface with nodal rc/Vsh_n and cell-centered metrics
  bool write_shock_surface_center_metrics_tecplot(const ShockMesh& M,
                                                  const TriMetrics& T,
                                                  const char* path) const;

  // 2) Bundle: surface nodal, surface cell metrics, and **volume POINT** zone
  //    Volume zone now includes n, V, **B**, and **divV** at each grid point.
  bool write_tecplot_dataset_bundle(const ShockMesh& M, const TriMetrics& T,
                                    const StepState& S, const BoxSpec& B,
                                    const char* path) const;

  // Directional diagnostic helper (Rdir, n̂, rc, Vsh_n)
  void diagnose_direction(const StepState& S, const double u[3],
                          double& Rdir_m, double n_hat[3],
                          double& rc_loc, double& Vsh_n) const;

private:
  Params P_;
};

} // namespace swcme3d

