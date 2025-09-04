#ifndef SWCME3D_HPP
#define SWCME3D_HPP
/*
swcme3d — Solar wind + CME forward-shock model (1D & 3D) with Tecplot writers
=============================================================================

WHAT
----
Lightweight C++17 model that returns solar-wind **density [m^-3]** and **bulk
velocity [m/s]** at arbitrary points and times, including a 3D CME forward
shock with oblique MHD compression. Includes triangulated shock surface and
Tecplot ASCII writers (surface + volume).

PHYSICS (key formulas)
----------------------
(E1) Leblanc–Dulk–Bougeret (1998) density (cm^-3, r in R_sun), scaled to match n(1 AU):
     n_cm3(r) = 3.3e5 r^-2 + 4.1e6 r^-4 + 8.0e7 r^-6;  n_m3 = 1e6 * scale * n_cm3.
(E2) Parker spiral (Parker 1958): |B|(r) ∝ (AU/r)^2 √[1+(Ω r sinθ / V_sw)^2], with −φ̂ component.
(E3) c_s = √(γ k_B T / m_p),  ρ = m_p n,  v_A = |B| / √(μ0 ρ).
(E5) Fast speed (oblique): c_f(θ)^2 = ½[(c_s^2+v_A^2)+√((c_s^2+v_A^2)^2−4 c_s^2 v_A^2 cos^2θ)].
(E7) DBM apex (Vršnak & Žic 2007; Vršnak+ 2013):
     du/dt = −Γ u|u|, u(t)=u0/(1+Γ u0 t), u0=V0−V_sw.
(E8) Apex radius: R_sh(t)=R0+V_sw t + (1/Γ) ln(1+Γ u0 t).
(E9) U1n = V_sh^n − V_sw^n;  M_fn = U1n / c_f(θ_Bn).
(E10) r_c = ((γ+1) M_fn^2)/((γ−1) M_fn^2 + 2) ≤ 4 (γ=5/3).
(S) Regions: shock → sheath (graded compression) → ejecta (lower density) → ambient, with
    independent C¹ smoothstep widths for the shock jump, sheath LE, and ejecta TE.

COORDINATES & UNITS
-------------------
• Heliocentric Cartesian; **Sun center is (0,0,0)**. r=|x| [m], t [s].
• Inputs/outputs: density [m^-3], velocity [m/s].

REFERENCES (titles)
-------------------
• Parker (1958) “Dynamics of the Interplanetary Gas and Magnetic Fields,” ApJ.
• Leblanc, Dulk & Bougeret (1998) “Tracing the Electron Density from the Corona to 1 AU,” Sol. Phys.
• Vršnak & Žic (2007) A&A; Vršnak et al. (2013) Sol. Phys. — Drag-Based Model.
• Priest (2014) *Magnetohydrodynamics of the Sun*, Cambridge UP.

TIPS
----
• For near real-time stepping, call prepare_step(t) every ~0.8 s; per-point eval is O(1).
• Sharper shock vs softer ejecta: decrease w_shock, increase w_te.
• Start at 1.05 R_sun via Params::r0_Rs=1.05.
*/

#include <cstddef>
#include <cstdint>
#include <vector>
#include <cstdio>

namespace swcme3d {

// ---------- constants (SI) ----------
extern const double PI, kB, mp, mu0, Rs, AU, Omega_sun;

// ---------- runtime I/O toggle ----------
void EnableIO();
void DisableIO();
bool IsIOEnabled();

// ---------- enums, params, state ----------
enum class ShockShape { Sphere, Ellipsoid, ConeSSE };

struct Params {
  // Ambient
  double V_sw_kms    = 400.0;
  double n1AU_cm3    = 7.0;
  double T_K         = 1.2e5;
  double B1AU_nT     = 5.0;

  // Shock initial (DBM) at r0_Rs * Rs
  double r0_Rs       = 1.05;
  double V0_sh_kms   = 1800.0;
  double Gamma_kmInv = 8.0e-8;

  // Thermo/geometry
  double gamma_ad    = 5.0/3.0;
  double sin_theta   = 1.0;

  // Sheath/ejecta shaping
  double sheath_thick_AU_at1AU = 0.12;
  double ejecta_thick_AU_at1AU = 0.25;
  double sheath_comp_floor     = 1.2;
  double sheath_ramp_power     = 1.0;
  double V_sheath_LE_factor    = 1.05;
  double f_ME                  = 0.5;
  double V_ME_factor           = 1.0;

  // Independent smoothstep half-widths (scale ∝ R_sh)
  double edge_smooth_shock_AU_at1AU = 0.005;
  double edge_smooth_le_AU_at1AU    = 0.010;
  double edge_smooth_te_AU_at1AU    = 0.020;

  // 3D shock shape
  ShockShape shape = ShockShape::Sphere;
  double cme_dir[3] = {1.0, 0.0, 0.0};
  double axis_ratio_y = 0.85;  // Ellipsoid b/a
  double axis_ratio_z = 0.70;  // Ellipsoid c/a
  double half_width_rad = 40.0 * (3.14159265358979323846/180.0); // Cone SSE half-width
  double flank_slowdown_m = 2.0; // Cone SSE f(θ)=cos^m θ
};

struct StepState {
  // Apex shock
  double t_s=0, r_sh_m=0, V_sh_ms=0, rc=1;

  // Up/down (apex)
  double V_up_ms=0, V_dn_ms=0;

  // Region geometry (apex-based)
  double dr_sheath_m=0, dr_me_m=0, r_le_m=0, r_te_m=0, inv_dr_sheath=0;

  // Smooth widths
  double w_shock_m=0, w_le_m=0, w_te_m=0;

  // Targets
  double rc_floor=1.2, V_sheath_LE_ms=0, V_ME_ms=0;

  // 3D frame & axes
  double e1[3]{}, e2[3]{}, e3[3]{};
  double a_m=0, b_m=0, c_m=0; // semi-axes (a=a_m=r_sh_m)
};

// Surface mesh + nodal fields
struct ShockMesh {
  std::vector<double> x,y,z;                    // node positions
  std::vector<double> rc, Vsh_n;                // nodal rc, normal shock speed
  std::vector<double> n_hat_x,n_hat_y,n_hat_z;  // nodal outward normals
  std::vector<int> tri_i,tri_j,tri_k;           // 1-based triangle connectivity
};

// Per-triangle metrics
struct TriMetrics {
  std::vector<double> area, rc_mean, Vsh_n_mean;
  std::vector<double> nx, ny, nz; // outward per-cell normal
  std::vector<double> cx, cy, cz; // centroid
};

// Volume sampling box
struct BoxSpec { double cx=0,cy=0,cz=0, hx=0,hy=0,hz=0; int Ni=0,Nj=0,Nk=0; };

// ---------- model ----------
class Model {
public:
  explicit Model(const Params&);

  // time-step state (~0.8s recommended for tight sync)
  StepState prepare_step(double t_s) const;

  // evaluators
  void evaluate_radii_fast(const StepState& S, const double* r_m,
                           double* n_m3, double* V_ms, std::size_t N) const;

  void evaluate_cartesian_fast(const StepState& S,
                               const double* x_m, const double* y_m, const double* z_m,
                               double* n_m3, double* Vx_ms, double* Vy_ms, double* Vz_ms,
                               std::size_t N) const;

  // shock surface & metrics
  ShockMesh build_shock_mesh(const StepState& S, std::size_t nTheta, std::size_t nPhi) const;
  void compute_triangle_metrics(const ShockMesh& M, TriMetrics& T) const;

  // Tecplot writers (single zones)
  bool write_shock_surface_tecplot(const ShockMesh& M, const char* path) const;
  bool write_shock_surface_center_metrics_tecplot(const ShockMesh& M, const TriMetrics& T, const char* path) const;
  bool write_box_tecplot(const StepState& S, const BoxSpec& B, const char* path) const;

  // Tecplot bundle (3 zones in one dataset: nodal surface, cell-centered surface, volume)
  bool write_tecplot_dataset_bundle(const ShockMesh& M, const TriMetrics& T,
                                    const StepState& S, const BoxSpec& B,
                                    const char* path) const;

  // convenience
  BoxSpec default_apex_box(const StepState& S, double half_AU, int N) const;

  // optional I/O helpers
  void write_step_csv_header(std::FILE* fp=stdout) const;
  void write_step_csv_line(const StepState& S, std::FILE* fp=stdout) const;

  // accessor
  double V_sw() const;

private:
  // internal helpers
  double leblanc_cm3(double r_m) const;
  double parker_norm(double r_m) const;
  void   make_frame(const double d[3], double e1[3], double e2[3], double e3[3]) const;
  void   shock_position_speed(double t, double& r_sh_m, double& V_sh_ms) const;

  void   B_parker_vec(const double u[3], double r_m, double V_sw_ms,
                      double B0_T, double Bout[3]) const;
  double fast_speed_oblique(double a, double vA, double cosTheta) const;
  double B_parker_mag(double r_m) const;
  double v_A(double r_m) const;
  double c_fast(double r_m) const;

  void shape_radius_normal(const StepState& S, double ux, double uy, double uz,
                           double& R, double n_hat[3]) const;

  void local_oblique_rc(const StepState& S,
                        const double u[3], const double n_hat[3],
                        double Rdir_m, double r_eval_m,
                        double& rc_out, double& Vsh_n_out, double& thetaBn_out) const;

private:
  Params P_;
  double leb_scale_{1.0};
  double B0_T_{0.0};
  double c_s_{0.0};
};

// clamp helper
double clamp01(double x);

} // namespace swcme3d

#endif // SWCME3D_HPP

