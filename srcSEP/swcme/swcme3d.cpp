/*
swcme3d.cpp
===========

Implementation for the solar wind + CME forward-shock model.

Updates in this revision
------------------------
1) Surface mesh no longer collapses to (0,0,0) for Cone/SSE:
   - The triangulation for Cone/SSE now covers only the physical cap
     0 <= theta <= half_width, avoiding any nodes outside the cone.
   - Previously, directions outside the cone were assigned R=0, which
     placed nodes at the origin and led to triangles effectively
     "connected" to (0,0,0). This is now eliminated.

2) Volume box enlarged & biased outward:
   - default_apex_box(...) now offsets the box center by +0.6 * half_AU
     along the apex direction (e1), so the box samples more of the region
     ahead (outside) of the forward shock.

3) Tecplot BLOCK zones keep the **wrapped output** (8 values per line)
   to avoid parsing errors due to extremely long lines.

Units & coordinates
-------------------
• Distances [m], time [s], density [m^-3], velocities [m/s].
• Sun center is (0,0,0).
• See swcme3d.hpp for full equations and references.
*/

#include "swcme3d.hpp"
#include <cmath>
#include <algorithm>
#include <array>
#include <cstdio>

namespace swcme3d {

//=============================================================================
// 1) Physical constants (SI)
//=============================================================================
const double PI        = 3.14159265358979323846;
const double kB        = 1.380649e-23;
const double mp        = 1.67262192369e-27;
const double mu0       = 4.0e-7 * PI;
const double Rs        = 6.957e8;
const double AU        = 1.495978707e11;
const double Omega_sun = 2.86533e-6;

//=============================================================================
// 2) Runtime I/O flag (printf-like progress dumps)
//=============================================================================
static bool& io_flag(){ static bool f=false; return f; }
void EnableIO()  { io_flag() = true;  }
void DisableIO() { io_flag() = false; }
bool IsIOEnabled(){ return io_flag(); }

//=============================================================================
// 3) Small vector helpers (local to this TU)
//=============================================================================
static inline double dot3(const double a[3], const double b[3]) {
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
static inline void cross3(const double a[3], const double b[3], double c[3]){
  c[0]=a[1]*b[2]-a[2]*b[1]; c[1]=a[2]*b[0]-a[0]*b[2]; c[2]=a[0]*b[1]-a[1]*b[0];
}
static inline double norm3(const double a[3]) { return std::sqrt(dot3(a,a)); }
static inline void unit3(double a[3]) { double n=norm3(a); if(n>0){a[0]/=n;a[1]/=n;a[2]/=n;} }
double clamp01(double x){ return x<0.0?0.0:(x>1.0?1.0:x); }

//=============================================================================
// 4) Model ctor: precompute scales for density and Parker-B normalization
//=============================================================================
Model::Model(const Params& p) : P_(p) {
  // Scale Leblanc so that n(1 AU) matches user-selected n1AU_cm3 (in cm^-3)
  leb_scale_ = P_.n1AU_cm3 / leblanc_cm3(AU);

  // Normalize Parker spiral magnitude so that |B|(1 AU) = B1AU_nT (in nT)
  B0_T_      = (P_.B1AU_nT * 1e-9) / parker_norm(AU);

  // Sound speed from selected (scalar) temperature
  c_s_       = std::sqrt(P_.gamma_ad * kB * P_.T_K / mp);
}

//=============================================================================
// 5) Time step: compute apex-shock kinematics + derived region widths
//=============================================================================
StepState Model::prepare_step(double t_s) const {
  StepState S{};
  S.t_s     = t_s;

  // Ambient solar wind speed (constant, radial)
  S.V_up_ms = V_sw();

  // DBM apex position/speed
  shock_position_speed(t_s, S.r_sh_m, S.V_sh_ms);

  // Apex compression (parallel proxy). Local oblique compression is computed
  // later, per direction (for 3D fields & nodal values).
  const double cf   = c_fast(S.r_sh_m);
  const double Vrel = std::max(S.V_sh_ms - S.V_up_ms, 0.0);
  if (cf <= 1.0 || Vrel <= 0.0) {
    S.rc = 1.0; S.V_dn_ms = S.V_up_ms;
  } else {
    const double Mn = Vrel / cf;
    const double g  = P_.gamma_ad;
    double rc = ((g+1.0)*Mn*Mn)/((g-1.0)*Mn*Mn + 2.0);
    if (!std::isfinite(rc)) rc = 1.0;
    S.rc     = std::min(4.0, std::max(1.0, rc)); // γ=5/3 ⇒ r_c ≤ 4
    S.V_dn_ms = S.V_sh_ms + (S.V_up_ms - S.V_sh_ms)/S.rc;
  }

  // Sheath + ejecta thickness scale self-similarly with R_sh
  S.dr_sheath_m   = std::max(0.0, P_.sheath_thick_AU_at1AU) * (S.r_sh_m / AU) * AU;
  S.dr_me_m       = std::max(0.0, P_.ejecta_thick_AU_at1AU) * (S.r_sh_m / AU) * AU;
  S.inv_dr_sheath = (S.dr_sheath_m > 0.0) ? 1.0 / S.dr_sheath_m : 0.0;

  // Leading/trailing edges (apex-radial)
  S.r_le_m = S.r_sh_m - S.dr_sheath_m;  // sheath LE
  S.r_te_m = S.r_le_m - S.dr_me_m;      // ejecta TE

  // Independent smoothstep half-widths (shock / sheath-LE / ejecta-TE)
  S.w_shock_m = std::max(0.0, P_.edge_smooth_shock_AU_at1AU) * (S.r_sh_m / AU) * AU;
  S.w_le_m    = std::max(0.0, P_.edge_smooth_le_AU_at1AU)    * (S.r_sh_m / AU) * AU;
  S.w_te_m    = std::max(0.0, P_.edge_smooth_te_AU_at1AU)    * (S.r_sh_m / AU) * AU;

  // Floor compression inside the sheath + target speeds at inner interfaces
  S.rc_floor       = std::max(1.0, std::min(S.rc, P_.sheath_comp_floor));
  S.V_sheath_LE_ms = P_.V_sheath_LE_factor * S.V_up_ms;  // toward ambient
  S.V_ME_ms        = P_.V_ME_factor        * S.V_up_ms;  // ejecta mean

  // Build an orthonormal frame (e1 along CME apex direction) and ellipsoid axes
  make_frame(P_.cme_dir, S.e1, S.e2, S.e3);
  S.a_m = S.r_sh_m;
  S.b_m = S.a_m * ((P_.axis_ratio_y>0)? P_.axis_ratio_y : 1.0);
  S.c_m = S.a_m * ((P_.axis_ratio_z>0)? P_.axis_ratio_z : 1.0);

  return S;
}

//=============================================================================
// 6) Evaluators: 1D (radial) and 3D (cartesian)
//=============================================================================
void Model::evaluate_radii_fast(const StepState& S, const double* r_m,
                                double* n_m3, double* V_ms, std::size_t N) const
{
  // Precompute Leblanc coefficients with scale so we don't redo it per point
  const double lebA = leb_scale_ * 3.3e5;
  const double lebB = leb_scale_ * 4.1e6;
  const double lebC = leb_scale_ * 8.0e7;

  const double Vup = S.V_up_ms, Vdn = S.V_dn_ms;
  const double rsh=S.r_sh_m, rle=S.r_le_m, rte=S.r_te_m;
  const double rc=S.rc, rc_floor=S.rc_floor, inv_drs=S.inv_dr_sheath;
  const double Vshe_LE=S.V_sheath_LE_ms, Vme=S.V_ME_ms;
  const double p = P_.sheath_ramp_power;
  const double w_sh=S.w_shock_m, w_le=S.w_le_m, w_te=S.w_te_m;

  // C¹ smoothstep mapped by radius around a center with half-width w
  auto sstep = [](double r, double r0, double w)->double {
    if (w <= 0.0) return (r <= r0) ? 1.0 : 0.0;
    double t = 0.5 + (r0 - r) / (2.0 * w); // t=1 at r=r0−w, t=0 at r=r0+w
    if (t <= 0.0) return 0.0;
    if (t >= 1.0) return 1.0;
    return t*t*(3.0 - 2.0*t);
  };

  for (std::size_t i=0;i<N;++i){
    const double r = r_m[i];

    // Ambient upstream density (Leblanc) in m^-3
    const double s = Rs / r;
    const double inv2=s*s, inv4=inv2*inv2, inv6=inv4*inv2;
    const double n_up = (lebA*inv2 + lebB*inv4 + lebC*inv6) * 1e6;

    // Sheath graded compression between rc at shock and rc_floor at sheath LE
    const double xi_raw = (S.dr_sheath_m > 0.0) ? (rsh - r) * inv_drs : 0.0;
    const double xi = (xi_raw<=0.0) ? 0.0 : (xi_raw>=1.0 ? 1.0 : xi_raw);
    const double ramp = (p==1.0) ? (1.0 - xi) : std::pow(1.0 - xi, p);
    const double Csheath = rc_floor + (rc - rc_floor)*ramp;
    const double n_sheath = Csheath * n_up;
    const double V_sheath = Vdn + (Vshe_LE - Vdn)*xi;

    // Ejecta target
    const double n_ejecta = std::max(0.0, P_.f_ME) * n_up;
    const double V_ejecta = S.V_ME_ms;

    // Blend shock → sheath region (around rsh)
    double a = (rc>1.0) ? sstep(r, rsh, w_sh) : 0.0;
    double n_mix = (1.0 - a)*n_up + a*n_sheath;
    double V_mix = (1.0 - a)*Vup  + a*V_sheath;

    // Blend sheath LE → ejecta (around rle)
    a = (rc>1.0) ? sstep(r, rle, w_le) : 0.0;
    n_mix = (1.0 - a)*n_mix + a*n_ejecta;
    V_mix = (1.0 - a)*V_mix + a*V_ejecta;

    // Blend ejecta TE → ambient (around rte)
    a = (rc>1.0) ? sstep(r, rte, w_te) : 0.0;
    n_mix = (1.0 - a)*n_mix + a*n_up;
    V_mix = (1.0 - a)*V_mix + a*Vup;

    n_m3[i] = n_mix;
    V_ms[i] = V_mix;
  }
}

void Model::evaluate_cartesian_fast(const StepState& S,
                                    const double* x_m, const double* y_m, const double* z_m,
                                    double* n_m3, double* Vx_ms, double* Vy_ms, double* Vz_ms,
                                    std::size_t N) const
{
  // Leblanc precompute with scale
  const double lebA = leb_scale_ * 3.3e5;
  const double lebB = leb_scale_ * 4.1e6;
  const double lebC = leb_scale_ * 8.0e7;

  const double Vup = S.V_up_ms;
  const double rc_floor=S.rc_floor, inv_drs=S.inv_dr_sheath;
  const double Vshe_LE=S.V_sheath_LE_ms, Vme=S.V_ME_ms;
  const double p = P_.sheath_ramp_power;
  const double w_sh=S.w_shock_m, w_le=S.w_le_m, w_te=S.w_te_m;

  auto sstep = [](double r, double r0, double w)->double {
    if (w <= 0.0) return (r <= r0) ? 1.0 : 0.0;
    double t = 0.5 + (r0 - r) / (2.0 * w);
    if (t <= 0.0) return 0.0;
    if (t >= 1.0) return 1.0;
    return t*t*(3.0 - 2.0*t);
  };

  for (std::size_t i=0; i<N; ++i) {
    // Sample position
    const double x=x_m[i], y=y_m[i], z=z_m[i];
    const double r2 = x*x + y*y + z*z;
    const double r  = std::sqrt(r2);
    const double r_safe = (r>0.0)? r : 1e-9;
    const double invr = 1.0/r_safe;
    const double u[3] = {x*invr, y*invr, z*invr}; // local radial unit

    // Local shock radius and outward normal for this direction
    double Rdir, n_hat[3];
    shape_radius_normal(S, u[0], u[1], u[2], Rdir, n_hat);

    // Upstream density at the sample radius
    const double s = Rs / r_safe;
    const double inv2=s*s, inv4=inv2*inv2, inv6=inv4*inv2;
    const double n_up = (lebA*inv2 + lebB*inv4 + lebC*inv6) * 1e6;

    // Local oblique compression ratio at the shock in this direction
    double rc_loc=1.0, Vsh_n=0.0, thBn=0.0;
    local_oblique_rc(S, u, n_hat, Rdir, std::max(r_safe, Rdir),
                     rc_loc, Vsh_n, thBn);

    // Sheath graded state with local rc
    const double xi_raw = (S.dr_sheath_m > 0.0) ? (Rdir - r_safe) * inv_drs : 0.0;
    const double xi = (xi_raw<=0.0) ? 0.0 : (xi_raw>=1.0 ? 1.0 : xi_raw);
    const double ramp = (p==1.0) ? (1.0 - xi) : std::pow(1.0 - xi, p);
    const double Csheath = rc_floor + (rc_loc - rc_floor)*ramp;
    const double n_sheath = Csheath * n_up;
    const double V_sheath = (S.V_dn_ms) + (Vshe_LE - S.V_dn_ms)*xi;

    // Ejecta
    const double n_ejecta = std::max(0.0, P_.f_ME) * n_up;
    const double V_ejecta = Vme;

    // Blends based on the *directional* shock radius Rdir
    double a = (rc_loc>1.0) ? sstep(r_safe, Rdir, w_sh) : 0.0;
    double n_mix = (1.0 - a)*n_up + a*n_sheath;
    double Vmag  = (1.0 - a)*Vup  + a*V_sheath;

    a = (rc_loc>1.0) ? sstep(r_safe, Rdir - S.dr_sheath_m, w_le) : 0.0;
    n_mix = (1.0 - a)*n_mix + a*n_ejecta;
    Vmag  = (1.0 - a)*Vmag + a*V_ejecta;

    a = (rc_loc>1.0) ? sstep(r_safe, Rdir - S.dr_me_m - S.dr_sheath_m, w_te) : 0.0;
    n_mix = (1.0 - a)*n_mix + a*n_up;
    Vmag  = (1.0 - a)*Vmag + a*Vup;

    // Return a radial velocity vector with magnitude Vmag
    Vx_ms[i] = Vmag * u[0];
    Vy_ms[i] = Vmag * u[1];
    Vz_ms[i] = Vmag * u[2];
    n_m3[i]  = n_mix;
  }
}

//=============================================================================
// 7) Shock surface triangulation (lat-lon grid) + per-triangle metrics
//=============================================================================
ShockMesh Model::build_shock_mesh(const StepState& S,
                                  std::size_t nTheta, std::size_t nPhi) const
{
  // Ensure minimal grid
  if (nTheta < 2) nTheta=2; if (nPhi < 3) nPhi=3;

  ShockMesh M;

  // --- KEY CHANGE ---
  // For Cone/SSE, only mesh the physical cap: 0 <= theta <= half_width.
  // This prevents any node from being placed at the origin (R=0).
  const double theta_max =
    (/*Cone cap*/ (P_.shape == ShockShape::ConeSSE))
      ? std::max(1e-6, std::min(P_.half_width_rad, PI))  // safe clamp
      : PI;                                              // full sphere/ellipsoid

  // With theta in [0, theta_max], we create (nTheta+1) rings and nPhi longitudes.
  const std::size_t Nv = (nTheta+1)*nPhi;
  M.x.resize(Nv); M.y.resize(Nv); M.z.resize(Nv);
  M.rc.resize(Nv); M.Vsh_n.resize(Nv);
  M.n_hat_x.resize(Nv); M.n_hat_y.resize(Nv); M.n_hat_z.resize(Nv);

  auto idx=[nPhi](std::size_t it,std::size_t ip){return it*nPhi + (ip % nPhi);};

  // Nodes (longitude wraps; the last theta ring is exactly theta_max)
  for (std::size_t it=0; it<=nTheta; ++it){
    const double th = (double(it)/double(nTheta))*theta_max;
    const double sth=std::sin(th), cth=std::cos(th);
    for (std::size_t ip=0; ip<nPhi; ++ip){
      const double ph=(double(ip)/double(nPhi))*2.0*PI;
      const std::size_t k=idx(it,ip);
      double u[3]={sth*std::cos(ph), sth*std::sin(ph), cth};

      double Rdir, n_hat[3];
      shape_radius_normal(S, u[0], u[1], u[2], Rdir, n_hat);

      // Rdir is guaranteed > 0 in the meshed angular range (no origin nodes)
      M.x[k]=Rdir*u[0]; M.y[k]=Rdir*u[1]; M.z[k]=Rdir*u[2];
      M.n_hat_x[k]=n_hat[0]; M.n_hat_y[k]=n_hat[1]; M.n_hat_z[k]=n_hat[2];

      double rc_loc=1.0, Vsh_n=0.0, thBn=0.0;
      local_oblique_rc(S, u, n_hat, Rdir, Rdir, rc_loc, Vsh_n, thBn);
      M.rc[k]=rc_loc; M.Vsh_n[k]=Vsh_n;
    }
  }

  // Two triangles per cell; for Cone cap the "rim" is at theta_max (open edge).
  for (std::size_t it=0; it<nTheta; ++it){
    for (std::size_t ip=0; ip<nPhi; ++ip){
      const std::size_t k00=idx(it  ,ip  );
      const std::size_t k01=idx(it  ,ip+1);
      const std::size_t k10=idx(it+1,ip  );
      const std::size_t k11=idx(it+1,ip+1);
      M.tri_i.push_back((int)k00+1); M.tri_j.push_back((int)k10+1); M.tri_k.push_back((int)k11+1);
      M.tri_i.push_back((int)k00+1); M.tri_j.push_back((int)k11+1); M.tri_k.push_back((int)k01+1);
    }
  }
  return M;
}

void Model::compute_triangle_metrics(const ShockMesh& M, TriMetrics& T) const {
  const std::size_t Ne = M.tri_i.size();
  T.area.resize(Ne);
  T.rc_mean.resize(Ne);
  T.Vsh_n_mean.resize(Ne);
  T.nx.resize(Ne); T.ny.resize(Ne); T.nz.resize(Ne);
  T.cx.resize(Ne); T.cy.resize(Ne); T.cz.resize(Ne);

  for (std::size_t e=0; e<Ne; ++e){
    const int i = M.tri_i[e]-1, j = M.tri_j[e]-1, k = M.tri_k[e]-1;
    const double p0[3]={M.x[i], M.y[i], M.z[i]};
    const double p1[3]={M.x[j], M.y[j], M.z[j]};
    const double p2[3]={M.x[k], M.y[k], M.z[k]};

    // Centroid
    double c[3]={(p0[0]+p1[0]+p2[0])/3.0, (p0[1]+p1[1]+p2[1])/3.0, (p0[2]+p1[2]+p2[2])/3.0};
    T.cx[e]=c[0]; T.cy[e]=c[1]; T.cz[e]=c[2];

    // Area and raw normal
    double a01[3]={p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]};
    double a02[3]={p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]};
    double nvec[3]; cross3(a01, a02, nvec);
    double area = 0.5 * norm3(nvec);
    T.area[e] = area;

    // Orient normal outward (away from the Sun)
    double nn = norm3(nvec);
    if (nn>0){ nvec[0]/=nn; nvec[1]/=nn; nvec[2]/=nn; }
    double rhat[3]={c[0],c[1],c[2]}; double rn=norm3(rhat); if(rn>0){rhat[0]/=rn;rhat[1]/=rn;rhat[2]/=rn;}
    if (dot3(nvec, rhat) < 0.0) { nvec[0]*=-1; nvec[1]*=-1; nvec[2]*=-1; }
    T.nx[e]=nvec[0]; T.ny[e]=nvec[1]; T.nz[e]=nvec[2];

    // Cell-centered means
    T.rc_mean[e]    = (M.rc[i]    + M.rc[j]    + M.rc[k])    / 3.0;
    T.Vsh_n_mean[e] = (M.Vsh_n[i] + M.Vsh_n[j] + M.Vsh_n[k]) / 3.0;
  }
}

//=============================================================================
// 8) Tecplot writers
//=============================================================================

// a) Nodal surface writer (FEPOINT). One node per line: naturally wrapped.
bool Model::write_shock_surface_tecplot(const ShockMesh& M, const char* path) const {
  const std::size_t Nv=M.x.size(), Ne=M.tri_i.size();
  std::FILE* fp=std::fopen(path,"w"); if(!fp) return false;

  std::fprintf(fp, "TITLE = \"Shock surface (nodal)\"\n");
  std::fprintf(fp, "VARIABLES = \"X\",\"Y\",\"Z\",\"rc\",\"Vsh_n\",\"nx\",\"ny\",\"nz\"\n");
  std::fprintf(fp, "ZONE T=\"surface\", N=%zu, E=%zu, F=FEPOINT, ET=TRIANGLE\n", Nv, Ne);

  for (std::size_t i=0;i<Nv;++i){
    std::fprintf(fp, "%.9e %.9e %.9e %.6e %.6e %.6e %.6e %.6e\n",
      M.x[i],M.y[i],M.z[i], M.rc[i],M.Vsh_n[i], M.n_hat_x[i],M.n_hat_y[i],M.n_hat_z[i]);
  }
  for (std::size_t e=0;e<Ne;++e)
    std::fprintf(fp, "%d %d %d\n", M.tri_i[e], M.tri_j[e], M.tri_k[e]);

  std::fclose(fp); return true;
}

// b) Cell-centered surface writer (BLOCK). **WRAPPED output** to avoid mega-lines.
bool Model::write_shock_surface_center_metrics_tecplot(
  const ShockMesh& M, const TriMetrics& T, const char* path) const
{
const std::size_t Nv = M.x.size();
const std::size_t Ne = M.tri_i.size();

std::FILE* fp = std::fopen(path,"w");
if(!fp) return false;

// Now includes rc and Vsh_n as NODAL variables, plus cell-centered metrics.
std::fprintf(fp, "TITLE = \"Shock surface (cell-centered metrics + nodal rc)\"\n");
std::fprintf(fp, "VARIABLES = "
                 "\"X\",\"Y\",\"Z\",\"rc\",\"Vsh_n\","
                 "\"area\",\"rc_mean\",\"Vsh_n_mean\",\"tnx\",\"tny\",\"tnz\",\"cx\",\"cy\",\"cz\"\n");
// [1-5] nodal (XYZ + rc + Vsh_n), [6-14] cell-centered
std::fprintf(fp, "ZONE T=\"surface_cells\", N=%zu, E=%zu, ZONETYPE=FETRIANGLE, DATAPACKING=BLOCK,\n", Nv, Ne);
std::fprintf(fp, "VARLOCATION=([1-5]=NODAL, [6-14]=CELLCENTERED)\n");

// Wrapped writers to avoid very long lines in BLOCK format.
const int PER_LINE = 8;
auto dumpN = [&](const std::vector<double>& a){
  int cnt=0; for (double v : a){ std::fprintf(fp, "%.9e ", v); if(++cnt==PER_LINE){ std::fprintf(fp, "\n"); cnt=0; } }
  if (cnt) std::fprintf(fp, "\n");
};
auto dumpE = dumpN;

// BLOCK order matches VARIABLES order:
// NODAL: X, Y, Z, rc, Vsh_n
dumpN(M.x); dumpN(M.y); dumpN(M.z); dumpN(M.rc); dumpN(M.Vsh_n);
// CELLCENTERED: area, rc_mean, Vsh_n_mean, tnx, tny, tnz, cx, cy, cz
dumpE(T.area); dumpE(T.rc_mean); dumpE(T.Vsh_n_mean);
dumpE(T.nx);   dumpE(T.ny);      dumpE(T.nz);
dumpE(T.cx);   dumpE(T.cy);      dumpE(T.cz);

// Connectivity (1-based)
for (std::size_t e=0; e<Ne; ++e)
  std::fprintf(fp, "%d %d %d\n", M.tri_i[e], M.tri_j[e], M.tri_k[e]);

std::fclose(fp);
return true;
}

// c) Volume box writer (POINT structured). One sample per line.
bool Model::write_box_tecplot(const StepState& S, const BoxSpec& B, const char* path) const {
  std::FILE* fp=std::fopen(path,"w"); if(!fp) return false;

  std::fprintf(fp, "TITLE = \"SW+CME box\"\n");
  std::fprintf(fp, "VARIABLES = \"X\",\"Y\",\"Z\",\"n\",\"Vx\",\"Vy\",\"Vz\"\n");
  std::fprintf(fp, "ZONE T=\"t=%.1fs\", I=%d, J=%d, K=%d, DATAPACKING=POINT\n", S.t_s, B.Ni, B.Nj, B.Nk);

  for (int k=0;k<B.Nk;++k){
    const double zk = B.cz + (-B.hz + (2.0*B.hz)*(k/(double)(B.Nk-1)));
    for (int j=0;j<B.Nj;++j){
      const double yj = B.cy + (-B.hy + (2.0*B.hy)*(j/(double)(B.Nj-1)));
      for (int i=0;i<B.Ni;++i){
        const double xi = B.cx + (-B.hx + (2.0*B.hx)*(i/(double)(B.Ni-1)));
        double n,Vx,Vy,Vz;
        evaluate_cartesian_fast(S, &xi,&yj,&zk, &n,&Vx,&Vy,&Vz, 1);
        std::fprintf(fp, "%.9e %.9e %.9e %.9e %.9e %.9e %.9e\n", xi,yj,zk, n,Vx,Vy,Vz);
      }
    }
  }

  std::fclose(fp);
  return true;
}

// d) Combined dataset: 3 zones in a single file (with wrapped BLOCK zone).
bool Model::write_tecplot_dataset_bundle(const ShockMesh& M, const TriMetrics& T,
  const StepState& S, const BoxSpec& B,
  const char* path) const
{
const std::size_t Nv = M.x.size();
const std::size_t Ne = M.tri_i.size();

std::FILE* fp = std::fopen(path,"w");
if(!fp) return false;
auto p = [&](const char* fmt, auto... args){ std::fprintf(fp, fmt, args...); };

// Shared VARIABLE list: superset across all zones
p("TITLE = \"SW+CME dataset\"\n");
p("VARIABLES = "
"\"X\",\"Y\",\"Z\",\"n\",\"Vx\",\"Vy\",\"Vz\",\"rc\",\"Vsh_n\","
"\"nx\",\"ny\",\"nz\",\"area\",\"rc_mean\",\"Vsh_n_mean\","
"\"tnx\",\"tny\",\"tnz\",\"cx\",\"cy\",\"cz\"\n");

// -------- Zone 1: surface (nodal, FEPOINT) --------
p("ZONE T=\"surface_nodal\", N=%zu, E=%zu, F=FEPOINT, ET=TRIANGLE\n", Nv, Ne);
for (std::size_t i=0;i<Nv;++i){
p("%.9e %.9e %.9e %.9e %.9e %.9e %.9e "
"%.6e %.6e %.6e %.6e %.6e "
"%.9e %.9e %.9e %.6e %.6e %.6e %.9e %.9e %.9e\n",
M.x[i],M.y[i],M.z[i],
0.0,0.0,0.0,0.0,          // n,Vx,Vy,Vz not defined here
M.rc[i],M.Vsh_n[i],       // nodal shock quantities
M.n_hat_x[i],M.n_hat_y[i],M.n_hat_z[i],
0.0,0.0,0.0,              // area, rc_mean, Vsh_n_mean
0.0,0.0,0.0,              // tnx,tny,tnz
0.0,0.0,0.0               // cx,cy,cz
);
}
for (std::size_t e=0;e<Ne;++e)
p("%d %d %d\n", M.tri_i[e], M.tri_j[e], M.tri_k[e]);

// -------- Zone 2: surface (BLOCK) with nodal rc/Vsh_n + cell-centered metrics --------
p("ZONE T=\"surface_cells\", N=%zu, E=%zu, ZONETYPE=FETRIANGLE, DATAPACKING=BLOCK,\n", Nv, Ne);
// Mark X,Y,Z **and** rc,Vsh_n as NODAL; others remain CELLCENTERED
p("VARLOCATION=([1-3,8-9]=NODAL, [4-7,10-21]=CELLCENTERED)\n");

const int PER_LINE = 8;
auto dumpN_wrap = [&](const std::vector<double>& a){
int cnt=0; for (double v : a){ std::fprintf(fp, "%.9e ", v); if(++cnt==PER_LINE){ std::fprintf(fp, "\n"); cnt=0; } }
if (cnt) std::fprintf(fp, "\n");
};
auto dumpE_wrap = dumpN_wrap;
auto dumpEzeros_wrap = [&](){
int cnt=0; for (std::size_t e=0;e<Ne;++e){ std::fprintf(fp, "0 "); if(++cnt==PER_LINE){ std::fprintf(fp, "\n"); cnt=0; } }
if (cnt) std::fprintf(fp, "\n");
};

// VARIABLES order:
// 1:X(N) 2:Y(N) 3:Z(N) 4:n(E) 5:Vx(E) 6:Vy(E) 7:Vz(E) 8:rc(N) 9:Vsh_n(N)
// 10:nx(E) 11:ny(E) 12:nz(E) 13:area(E) 14:rc_mean(E) 15:Vsh_n_mean(E)
// 16:tnx(E) 17:tny(E) 18:tnz(E) 19:cx(E) 20:cy(E) 21:cz(E)

// NODAL first: X,Y,Z
dumpN_wrap(M.x);
dumpN_wrap(M.y);
dumpN_wrap(M.z);

// CELLCENTERED placeholders for n, Vx, Vy, Vz (unused in this surface zone)
dumpEzeros_wrap(); // n
dumpEzeros_wrap(); // Vx
dumpEzeros_wrap(); // Vy
dumpEzeros_wrap(); // Vz

// NODAL shock properties: rc and Vsh_n
dumpN_wrap(M.rc);
dumpN_wrap(M.Vsh_n);

// CELLCENTERED placeholders for (nodal) normals nx,ny,nz in this zone
dumpEzeros_wrap(); // nx
dumpEzeros_wrap(); // ny
dumpEzeros_wrap(); // nz

// CELLCENTERED per-triangle metrics
dumpE_wrap(T.area);
dumpE_wrap(T.rc_mean);
dumpE_wrap(T.Vsh_n_mean);
dumpE_wrap(T.nx);  // triangle normal x
dumpE_wrap(T.ny);  // triangle normal y
dumpE_wrap(T.nz);  // triangle normal z
dumpE_wrap(T.cx);  // triangle centroid x
dumpE_wrap(T.cy);  // triangle centroid y
dumpE_wrap(T.cz);  // triangle centroid z

// Connectivity
for (std::size_t e=0;e<Ne;++e)
p("%d %d %d\n", M.tri_i[e], M.tri_j[e], M.tri_k[e]);

// -------- Zone 3: volume (POINT structured) --------
p("ZONE T=\"volume_box\", I=%d, J=%d, K=%d, DATAPACKING=POINT\n", B.Ni, B.Nj, B.Nk);
for (int kk=0; kk<B.Nk; ++kk){
const double zk = B.cz + (-B.hz + (2.0*B.hz)*(kk/(double)(B.Nk-1)));
for (int jj=0; jj<B.Nj; ++jj){
const double yj = B.cy + (-B.hy + (2.0*B.hy)*(jj/(double)(B.Nj-1)));
for (int ii=0; ii<B.Ni; ++ii){
const double xi = B.cx + (-B.hx + (2.0*B.hx)*(ii/(double)(B.Ni-1)));
double n,Vx,Vy,Vz;
evaluate_cartesian_fast(S, &xi,&yj,&zk, &n,&Vx,&Vy,&Vz, 1);
p("%.9e %.9e %.9e %.9e %.9e %.9e %.9e "
"%.9e %.9e %.9e %.9e %.9e "
"%.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e\n",
xi,yj,zk, n,Vx,Vy,Vz,
0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0);
}
}
}

std::fclose(fp);
return true;
}


//=============================================================================
// 9) Convenience + internals
//=============================================================================
BoxSpec Model::default_apex_box(const StepState& S, double half_AU, int N) const {
  // --- KEY CHANGE ---
  // Bias the box center OUTWARD along the apex direction (e1) so that
  // more of the box volume lies AHEAD of the forward shock.
  // Adjust OUTSIDE_BIAS if you want a different skew (0..1).
  constexpr double OUTSIDE_BIAS = 0.6; // 60% of half-size outward shift

  const double h = half_AU * AU;
  const double shift = OUTSIDE_BIAS * h;

  BoxSpec B{};
  B.cx = S.a_m*S.e1[0] + shift*S.e1[0];
  B.cy = S.a_m*S.e1[1] + shift*S.e1[1];
  B.cz = S.a_m*S.e1[2] + shift*S.e1[2];

  // Axis-aligned half sizes
  B.hx = h; B.hy = h; B.hz = h;

  // Resolution
  B.Ni = N; B.Nj = N; B.Nk = N;
  return B;
}

void Model::write_step_csv_header(std::FILE* fp) const {
  if (!IsIOEnabled()) return; if(!fp) fp=stdout;
  std::fprintf(fp,"t_s,R_sh_AU,V_sh_km_s,rc,R_LE_AU,R_TE_AU\n");
}
void Model::write_step_csv_line(const StepState& S, std::FILE* fp) const {
  if (!IsIOEnabled()) return; if(!fp) fp=stdout;
  std::fprintf(fp,"%.3f,%.9f,%.6f,%.6f,%.9f,%.9f\n",
    S.t_s, S.r_sh_m/AU, S.V_sh_ms/1e3, S.rc, S.r_le_m/AU, S.r_te_m/AU);
}

double Model::V_sw() const { return P_.V_sw_kms * 1e3; }

// Leblanc density (cm^-3) without scale; r in meters ⇒ r/Rs dimensionless inside
double Model::leblanc_cm3(double r_m) const {
  const double s = Rs / r_m;
  const double inv2=s*s, inv4=inv2*inv2, inv6=inv4*inv2;
  return 3.3e5*inv2 + 4.1e6*inv4 + 8.0e7*inv6;
}

// Parker spiral |B| normalization, relative to B0 at 1 AU.
double Model::parker_norm(double r_m) const {
  const double Br_over_B0 = (AU / r_m) * (AU / r_m);
  const double x = Omega_sun * r_m * P_.sin_theta / V_sw();
  return Br_over_B0 * std::sqrt(1.0 + x*x);
}

// Orthonormal frame from an arbitrary direction (robust even if nearly axis-aligned)
void Model::make_frame(const double d[3], double e1[3], double e2[3], double e3[3]) const {
  double nrm = std::sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]); if(nrm==0) nrm=1.0;
  e1[0]=d[0]/nrm; e1[1]=d[1]/nrm; e1[2]=d[2]/nrm;

  // Pick a helper axis roughly orthogonal to e1, then Gram–Schmidt
  double h0 = (std::fabs(e1[0])<0.9)?1.0:0.0, h1 = (h0==1.0)?0.0:1.0, h2=0.0;
  double dot = h0*e1[0]+h1*e1[1]+h2*e1[2];
  double ex=h0-dot*e1[0], ey=h1-dot*e1[1], ez=h2-dot*e1[2];
  double rn=std::sqrt(ex*ex+ey*ey+ez*ez); if(rn==0) rn=1.0;
  e2[0]=ex/rn; e2[1]=ey/rn; e2[2]=ez/rn;

  // e3 = e1 × e2
  e3[0]= e1[1]*e2[2]-e1[2]*e2[1];
  e3[1]= e1[2]*e2[0]-e1[0]*e2[2];
  e3[2]= e1[0]*e2[1]-e1[1]*e2[0];
}

// DBM apex position/speed (closed form using ln(1 + Γ u0 t))
void Model::shock_position_speed(double t, double& r_sh_m, double& V_sh_ms) const {
  const double r0 = P_.r0_Rs * Rs;
  const double Vsw = V_sw();
  const double V0  = P_.V0_sh_kms * 1e3;
  const double G   = P_.Gamma_kmInv * 1e-3;  // [1/m]
  const double u0  = V0 - Vsw;

  if (std::abs(u0) < 1e-12) { r_sh_m=r0+Vsw*t; V_sh_ms=Vsw; return; }

  if (u0 > 0.0) {
    // decelerating (drag opposes motion)
    const double denom=1.0 + G*u0*t;
    const double u=u0/denom;
    const double integral=(1.0/G)*std::log(denom);
    r_sh_m=r0+Vsw*t+integral; V_sh_ms=Vsw+u;
  } else {
    // accelerating toward Vsw (u0 negative)
    double denom=1.0 - G*u0*t; if(denom<=0.0) denom=1e-30;
    const double u=u0/denom;
    const double integral=-(1.0/G)*std::log(denom);
    r_sh_m=r0+Vsw*t+integral; V_sh_ms=Vsw+u;
  }
}

// Parker spiral vector at (direction u, radius r_m)
void Model::B_parker_vec(const double u[3], double r_m, double V_sw_ms, double B0_T, double Bout[3]) const {
  // Radial magnitude scales as (AU/r)^2
  const double scale = (AU / r_m)*(AU / r_m);
  const double Br = B0_T * scale;

  // φ̂ direction ≈ Ω̂ × r̂; normalize to get unit azimuthal
  const double OMEGA_HAT[3]={0,0,1};
  double eph[3]; cross3(OMEGA_HAT, u, eph);
  double ephn = norm3(eph);
  if (ephn < 1e-16) { eph[0]=0.0; eph[1]=1.0; eph[2]=0.0; ephn=1.0; }
  eph[0]/=ephn; eph[1]/=ephn; eph[2]/=ephn;

  const double sinTheta = ephn; // ≡ |Ω̂×r̂|
  const double Bphi = - Br * (Omega_sun * r_m * sinTheta / V_sw_ms);

  // Total B = Br r̂ + Bphi φ̂
  Bout[0] = Br*u[0] + Bphi*eph[0];
  Bout[1] = Br*u[1] + Bphi*eph[1];
  Bout[2] = Br*u[2] + Bphi*eph[2];
}

// Fast-mode speed for oblique propagation: c_f(θ)
double Model::fast_speed_oblique(double a, double vA, double cosTheta) const {
  const double c2 = a*a + vA*vA;
  const double D  = (c2*c2) - 4.0*a*a*vA*vA*(cosTheta*cosTheta);
  const double root = (D>0.0)? std::sqrt(D) : 0.0;
  const double cf2 = 0.5*( c2 + root );
  return std::sqrt(std::max(0.0, cf2));
}

double Model::B_parker_mag(double r_m) const { return B0_T_ * parker_norm(r_m); }
double Model::v_A(double r_m) const {
  const double B = B_parker_mag(r_m);
  const double n = leb_scale_ * leblanc_cm3(r_m) * 1e6; // to m^-3
  const double rho = mp * n;
  return B / std::sqrt(mu0 * rho);
}
double Model::c_fast(double r_m) const { const double vA=v_A(r_m); return std::sqrt(c_s_*c_s_ + vA*vA); }

// Shock radius & normal for selected shape and direction u=(ux,uy,uz)
void Model::shape_radius_normal(const StepState& S, double ux, double uy, double uz,
                                double& R, double n_hat[3]) const
{
  if (P_.shape == ShockShape::Sphere) {
    // Sphere: radius is constant; normal aligns with radial direction
    R = S.r_sh_m; n_hat[0]=ux; n_hat[1]=uy; n_hat[2]=uz; return;
  }

  // Rotate direction into local CME frame (e1,e2,e3)
  const double u1 = ux*S.e1[0] + uy*S.e1[1] + uz*S.e1[2];
  const double u2 = ux*S.e2[0] + uy*S.e2[1] + uz*S.e2[2];
  const double u3 = ux*S.e3[0] + uy*S.e3[1] + uz*S.e3[2];

  if (P_.shape == ShockShape::Ellipsoid) {
    // Ellipsoid with semi-axes (a,b,c) in (e1,e2,e3) directions.
    // Directional radius is solution to (x/a)^2 + (y/b)^2 + (z/c)^2 = 1.
    const double invR2 = (u1*u1)/(S.a_m*S.a_m) + (u2*u2)/(S.b_m*S.b_m) + (u3*u3)/(S.c_m*S.c_m);
    R = 1.0 / std::sqrt(invR2);

    // Surface normal ∝ (x/a^2, y/b^2, z/c^2) evaluated at the point
    const double px = R*ux, py=R*uy, pz=R*uz;
    double nx =  px/(S.a_m*S.a_m);
    double ny =  py/(S.b_m*S.b_m);
    double nz =  pz/(S.c_m*S.c_m);
    double nn = std::sqrt(nx*nx+ny*ny+nz*nz); if(nn>0){nx/=nn;ny/=nn;nz/=nn;}
    n_hat[0]=nx; n_hat[1]=ny; n_hat[2]=nz;
    return;
  }

  // Cone/SSE (self-similar expansion; simple angular mask with cos^m scaling)
  const double cosT = std::max(-1.0, std::min(1.0, ux*S.e1[0]+uy*S.e1[1]+uz*S.e1[2]));
  const double theta = std::acos(cosT);
  const bool inside = (theta <= P_.half_width_rad);

  // For mesh generation we now ONLY sample inside the cone (handled by caller).
  // Here we provide the local radius when inside; caller avoids outside ranges.
  const double f = inside ? std::pow(std::max(0.0, cosT), P_.flank_slowdown_m) : 0.0;
  R = S.r_sh_m * f;

  // Normal ~ radial direction for simple cone front
  double nx=ux, ny=uy, nz=uz; double nn=std::sqrt(nx*nx+ny*ny+nz*nz); if(nn>0){nx/=nn;ny/=nn;nz/=nn;}
  n_hat[0]=nx; n_hat[1]=ny; n_hat[2]=nz;
}

// Local oblique compression and normal shock speed for direction u at radius r_eval_m
void Model::local_oblique_rc(const StepState& S,
                             const double u[3], const double n_hat[3],
                             double Rdir_m, double r_eval_m,
                             double& rc_out, double& Vsh_n_out, double& thetaBn_out) const
{
  const double Vsw = V_sw();

  // Parker magnetic field (vector) at r_eval in direction u
  double Bv[3]; B_parker_vec(u, r_eval_m, Vsw, B0_T_, Bv);
  const double Bmag = norm3(Bv);
  double b_hat[3] = { (Bmag>0)? Bv[0]/Bmag : 0.0,
                      (Bmag>0)? Bv[1]/Bmag : 0.0,
                      (Bmag>0)? Bv[2]/Bmag : 0.0 };
  const double cosBn = clamp01(std::fabs(dot3(b_hat, n_hat)));
  thetaBn_out = std::acos(cosBn);

  // Upstream density & Alfvén speed at r_eval
  const double n_up_m3 = leb_scale_ * leblanc_cm3(r_eval_m) * 1e6;
  const double rho1    = mp * n_up_m3;
  const double vA      = (Bmag>0 && rho1>0) ? (Bmag/std::sqrt(mu0*rho1)) : 0.0;
  const double a       = c_s_;
  const double cf_n    = fast_speed_oblique(a, vA, cosBn);

  // Self-similar directional shock speed magnitude
  const double a_axis  = (S.a_m>0.0)? S.a_m : S.r_sh_m;
  const double Vsh_mag = (a_axis>0.0) ? (S.V_sh_ms * (Rdir_m / a_axis)) : S.V_sh_ms;

  // Normal inflow speed into the shock (projected radial SW speed on the normal)
  const double Vsw_n   = Vsw * (u[0]*n_hat[0] + u[1]*n_hat[1] + u[2]*n_hat[2]);
  const double U1n     = std::max(Vsh_mag - Vsw_n, 0.0);

  if (cf_n <= 1e-6 || U1n <= 0.0) { rc_out=1.0; Vsh_n_out=Vsh_mag; return; }

  // Fast-mode normal Mach number and RH compression (capped at 4)
  const double g   = P_.gamma_ad;
  const double Mfn = U1n / cf_n;
  double rc = ((g+1.0)*Mfn*Mfn)/((g-1.0)*Mfn*Mfn + 2.0);
  if (!std::isfinite(rc)) rc = 1.0;
  rc_out    = std::min(4.0, std::max(1.0, rc));
  Vsh_n_out = Vsh_mag;
}

} // namespace swcme3d

