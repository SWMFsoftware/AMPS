/***************************************************************************************
 * curl_curl_e.cpp  (C++11/C++14 compatible)
 * -------------------------------------------------------------------------------------
 * PURPOSE
 *   Validate and cross-verify the discrete operator
 *
 *       curl_curl(E) = ∇×(∇×E) = ∇(∇·E) − ∇²E
 *
 *   at CORNER nodes (collocated stencils) for orders:
 *   • 2nd (compact & wide)  • 4th  • 6th  • 8th
 *
 * WHAT THIS TEST COMPUTES
 *   Numerical operators, all applied to the same periodic vector field E(x,y,z):
 *     1) CC_num    := curl_curl(E)        (built from cCurlCurlEStencil rows)
 *     2) GD_num    := grad_div(E)         (built from cGradDivEStencil rows)
 *     3) Lap_num   := ∇²E                 (built from cLaplacianStencil, per component)
 *     4) ID_num    := GD_num − Lap_num    (pure numerical identity reconstruction)
 *
 *   Analytic references at each grid point (x,y,z):
 *     • CC_ana   = ∇(∇·E) − ∇²E
 *     • GD_ana   = ∇( (a+b+c) cos(ax)cos(by)cos(cz) )
 *     • Lap_ana  = −(a²+b²+c²) * E    (component-wise)
 *
 * FIELD (periodic on [0,L]^3)
 *   E(x,y,z) = ( sin(ax) cos(by) cos(cz),
 *                cos(ax) sin(by) cos(cz),
 *                cos(ax) cos(by) sin(cz) )
 *   with a=2π, b=3π, c=5π (unequal to avoid trivial cancellations).
 *
 * STENCIL APPLICATION
 *   Mirrors curl_b.cpp: each sub-stencil is exported via
 *     ExportStencil(cStencil::cStencilData*)
 *   then applied by looping integer offsets and multiplying by stored weights.
 *   All metric factors (1/dx, 1/dx^2, etc.) are already baked into the
 *   exported coefficients — do NOT rescale during application.
 *
 * GLOBAL METRICS (printed per variant/order)
 *   • CC vs analytic   : Linf,  RelL2   of (CC_num − CC_ana)
 *   • CC vs (GD−Lap)   : Linf,  RelL2   of (CC_num − ID_num)
 *   • GD vs analytic   : Linf,  RelL2   of (GD_num − GD_ana)
 *   • Lap vs analytic  : Linf,  RelL2   of (Lap_num − Lap_ana)
 *   (RelL2 := ||error||_2 / ||reference||_2 over the whole grid.)
 *
 * PER-POINT TABLES (component-wise, at an interior node)
 *   1) CurlCurl comparison (direct & identity, both vs analytic)
 *
 *     Component          Analytic            Num CurlCurl         Num (GD-Lap)       |Err(CC)|    |Err(ID)|
 *                         CC_ana               CC_num               ID_num          |CC-CC_ana|  |ID-CC_ana|
 *
 *     Meaning:
 *       “Analytic”     = CC_ana component from closed-form ∇(∇·E) − ∇²E
 *       “Num CurlCurl” = CC_num component (direct curl_curl stencil)
 *       “Num (GD-Lap)” = ID_num component (grad_div − laplacian, both numerical)
 *       “|Err(CC)|”    = |CC_num − CC_ana|,  “|Err(ID)|” = |ID_num − CC_ana|
 *
 *   2) GradDiv vs Analytic
 *
 *     Component          Analytic GD         Numerical GD         |Err(GD)|
 *                         GD_ana               GD_num            |GD-GD_ana|
 *
 *   3) Laplacian vs Analytic
 *
 *     Component          Analytic Lap        Numerical Lap        |Err(Lap)|
 *                         Lap_ana              Lap_num           |Lap-Lap_ana|
 *
 * EXPECTED BEHAVIOR
 *   • |Err(CC)| and |Err(ID)| decrease with order and with grid refinement.
 *   • GD and Lap errors independently converge at their formal orders.
 *   • CC_num and ID_num should be very close (consistency of the discrete identity).
 *
 * CONVERGENCE SUMMARY (added at end of run)
 *   After running across Ns = {16,24,32,48,64}, we print a compact table with
 *   L∞ and L2 (RMS) errors of CC vs analytic for 2nd-compact (“Second”), 4th, 6th, 8th,
 *   plus observed orders computed as:
 *       p = log(e_prev/e_curr) / log(N_curr/N_prev),   with h ~ 1/N.
 ***************************************************************************************/

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <utility>
#include <functional>

#include "test_harness.h"
#include "test_register.h"
#include "test_force_link_all.h"
#include "pic.h"  // brings stencil types and ExportStencil

using namespace PIC::FieldSolver::Electromagnetic::ECSIM::Stencil;

namespace {

// ---------- small helpers ----------
inline int wrap(int i, int N){ int r = i % N; return (r<0)? r+N : r; }
struct Vec3 { double x,y,z; };

// Analytic field
inline Vec3 analyticE(double x, double y, double z, double a, double b, double c){
  Vec3 v;
  v.x = std::sin(a*x)*std::cos(b*y)*std::cos(c*z);
  v.y = std::cos(a*x)*std::sin(b*y)*std::cos(c*z);
  v.z = std::cos(a*x)*std::cos(b*y)*std::sin(c*z);
  return v;
}

// Analytic pieces
inline Vec3 analyticCurlCurlE(double x,double y,double z,double a,double b,double c){
  const double sx = std::sin(a*x), cx = std::cos(a*x);
  const double sy = std::sin(b*y), cy = std::cos(b*y);
  const double sz = std::sin(c*z), cz = std::cos(c*z);

  const double Ex = sx*cy*cz;
  const double Ey = cx*sy*cz;
  const double Ez = cx*cy*sz;

  const double K   = (a + b + c);
  const double dDiv_dx = K * (-a*sx)*cy*cz;
  const double dDiv_dy = K * cx*(-b*sy)*cz;
  const double dDiv_dz = K * cx*cy*(-c*sz);

  const double lam = (a*a + b*b + c*c);
  const double lapx = -lam * Ex;
  const double lapy = -lam * Ey;
  const double lapz = -lam * Ez;

  return { dDiv_dx - lapx, dDiv_dy - lapy, dDiv_dz - lapz };
}

inline Vec3 analyticGradDivE(double x,double y,double z,double a,double b,double c){
  const double sx = std::sin(a*x), cx = std::cos(a*x);
  const double sy = std::sin(b*y), cy = std::cos(b*y);
  const double sz = std::sin(c*z), cz = std::cos(c*z);
  const double K = (a + b + c);
  return { K * (-a*sx)*cy*cz,
           K * cx*(-b*sy)*cz,
           K * cx*cy*(-c*sz) };
}

inline Vec3 analyticLaplacianE(double x,double y,double z,double a,double b,double c){
  const double sx = std::sin(a*x), cx = std::cos(a*x);
  const double sy = std::sin(b*y), cy = std::cos(b*y);
  const double sz = std::sin(c*z), cz = std::cos(c*z);
  const double Ex = sx*cy*cz;
  const double Ey = cx*sy*cz;
  const double Ez = cx*cy*sz;
  const double lam = (a*a + b*b + c*c);
  return { -lam*Ex, -lam*Ey, -lam*Ez };
}

// Exported stencil application (curl_b.cpp style)
static inline double apply_exported(const cStencil::cStencilData& S,
                                    const std::vector<double>& F,
                                    int i, int j, int k,
                                    int Nx, int Ny, int Nz)
{
  double acc = 0.0;
  for (int n=0; n<S.Length; ++n) {
    const int ii = wrap(i + S.Data[n].i, Nx);
    const int jj = wrap(j + S.Data[n].j, Ny);
    const int kk = wrap(k + S.Data[n].k, Nz);
    const size_t idx = (size_t)kk*Ny*Nx + (size_t)jj*Nx + (size_t)ii;
    acc += S.Data[n].a * F[idx];
  }
  return acc;
}

// ----- pretty printers for interior-point comparisons -----

static void print_point_curlcurl(int N, double L,
                                 const std::vector<double>& CCx,
                                 const std::vector<double>& CCy,
                                 const std::vector<double>& CCz,
                                 const std::vector<double>& CCxA,
                                 const std::vector<double>& CCyA,
                                 const std::vector<double>& CCzA,
                                 const std::vector<double>& IDx,
                                 const std::vector<double>& IDy,
                                 const std::vector<double>& IDz,
                                 const char* flavor_label)
{
  const double dx = L/N;
  const int ii = N/4, jj = N/3, kk = (2*N)/5;   // interior, nontrivial
  const int Nx=N, Ny=N, Nz=N;
  const size_t idx = (size_t)kk*Ny*Nx + (size_t)jj*Nx + (size_t)ii;
  const double x = ii*dx, y = jj*dx, z = kk*dx;

  auto line = [&](const char* name, double a_, double n_, double n2){
    const double err1 = std::abs(n_  - a_);
    const double err2 = std::abs(n2 - a_);
    std::cout << "    " << std::left << std::setw(16) << name
              << std::right << std::scientific << std::setprecision(8)
              << std::setw(16) << a_  << "   "
              << std::setw(16) << n_  << "   "
              << std::setw(16) << n2  << "   "
              << std::setprecision(3) << std::setw(10) << err1
              << std::setw(12) << err2 << "\n";
  };

  std::cout << "\n[" << flavor_label << "] curl_curl(E) @ interior point\n"
            << "  N=" << N << ", (i,j,k)=(" << ii << "," << jj << "," << kk << ")"
            << ", (x,y,z)=(" << std::fixed << std::setprecision(6)
            << x << ", " << y << ", " << z << ")\n"
            << "  --------------------------------------------------------------------------------------------------------------\n"
            << "    Component          Analytic            Num CurlCurl         Num (GD-Lap)       |Err(CC)|    |Err(ID)|\n"
            << "                        CC_ana               CC_num               ID_num          |CC-CC_ana|  |ID-CC_ana|\n"
            << "  --------------------------------------------------------------------------------------------------------------\n";

  line("(CurlCurlE)_x", CCxA[idx], CCx[idx], IDx[idx]);
  line("(CurlCurlE)_y", CCyA[idx], CCy[idx], IDy[idx]);
  line("(CurlCurlE)_z", CCzA[idx], CCz[idx], IDz[idx]);
}

static void print_point_gd(int N, double L,
                           const std::vector<double>& GDx,
                           const std::vector<double>& GDy,
                           const std::vector<double>& GDz,
                           const std::vector<double>& GDAx,
                           const std::vector<double>& GDAy,
                           const std::vector<double>& GDAz,
                           const char* flavor_label)
{
  const double dx = L/N;
  const int ii = N/4, jj = N/3, kk = (2*N)/5;
  const int Nx=N, Ny=N, Nz=N;
  const size_t idx = (size_t)kk*Ny*Nx + (size_t)jj*Nx + (size_t)ii;
  const double x = ii*dx, y = jj*dx, z = kk*dx;

  auto line = [&](const char* name, double a_, double n_){
    const double err = std::abs(n_  - a_);
    std::cout << "    " << std::left << std::setw(16) << name
              << std::right << std::scientific << std::setprecision(8)
              << std::setw(16) << a_  << "   "
              << std::setw(16) << n_  << "   "
              << std::setprecision(3) << std::setw(10) << err << "\n";
  };

  std::cout << "\n[" << flavor_label << "] grad_div(E) @ interior point\n"
            << "  --------------------------------------------------------------------------------------------------------------\n"
            << "    Component          Analytic GD         Numerical GD         |Err(GD)|\n"
            << "                         GD_ana              GD_num            |GD-GD_ana|\n"
            << "  --------------------------------------------------------------------------------------------------------------\n";

  line("(GradDivE)_x", GDAx[idx], GDx[idx]);
  line("(GradDivE)_y", GDAy[idx], GDy[idx]);
  line("(GradDivE)_z", GDAz[idx], GDz[idx]);
}

static void print_point_lap(int N, double L,
                            const std::vector<double>& LPx,
                            const std::vector<double>& LPy,
                            const std::vector<double>& LPz,
                            const std::vector<double>& LPAx,
                            const std::vector<double>& LPAy,
                            const std::vector<double>& LPAz,
                            const char* flavor_label)
{
  const double dx = L/N;
  const int ii = N/4, jj = N/3, kk = (2*N)/5;
  const int Nx=N, Ny=N, Nz=N;
  const size_t idx = (size_t)kk*Ny*Nx + (size_t)jj*Nx + (size_t)ii;
  const double x = ii*dx, y = jj*dx, z = kk*dx;

  auto line = [&](const char* name, double a_, double n_){
    const double err = std::abs(n_  - a_);
    std::cout << "    " << std::left << std::setw(16) << name
              << std::right << std::scientific << std::setprecision(8)
              << std::setw(16) << a_  << "   "
              << std::setw(16) << n_  << "   "
              << std::setprecision(3) << std::setw(10) << err << "\n";
  };

  std::cout << "\n[" << flavor_label << "] laplacian(E) @ interior point\n"
            << "  --------------------------------------------------------------------------------------------------------------\n"
            << "    Component          Analytic Lap        Numerical Lap        |Err(Lap)|\n"
            << "                         Lap_ana             Lap_num           |Lap-Lap_ana|\n"
            << "  --------------------------------------------------------------------------------------------------------------\n";

  line("(LaplacianE)_x", LPAx[idx], LPx[idx]);
  line("(LaplacianE)_y", LPAy[idx], LPy[idx]);
  line("(LaplacianE)_z", LPAz[idx], LPz[idx]);
}

// ---------- error structure ----------
// NOTE: l2rel is still used for on-screen diagnostics; l2abs (RMS) is added for convergence table.
struct ErrStats { double linf=0.0, l2abs=0.0, l2rel=0.0; };

// ---------- core runner ----------
template<typename BuildGradDiv, typename BuildLap, typename BuildCurlCurl>
ErrStats run_one_order(const char* label,
                       BuildGradDiv build_grad_div,
                       BuildLap     build_lap,
                       BuildCurlCurl build_curl_curl,
                       int N, double L,
                       double a, double b, double c)
{
  const double dx=L/N, dy=L/N, dz=L/N;
  const int Nx=N, Ny=N, Nz=N;
  const size_t NT = (size_t)Nx*Ny*Nz;

  // initialize E at corners
  std::vector<double> Ex(NT), Ey(NT), Ez(NT);
  for (int k=0; k<Nz; ++k){
    const double z=k*dz;
    for (int j=0; j<Ny; ++j){
      const double y=j*dy;
      for (int i=0; i<Nx; ++i){
        const double x=i*dx;
        const Vec3 E = analyticE(x,y,z,a,b,c);
        const size_t id = (size_t)k*Ny*Nx + (size_t)j*Nx + (size_t)i;
        Ex[id]=E.x; Ey[id]=E.y; Ez[id]=E.z;
      }
    }
  }

  // build stencils
  cGradDivEStencil  G[3];
  cLaplacianStencil Ls;
  cCurlCurlEStencil CC[3];

  build_grad_div(G, dx,dy,dz);
  build_lap(&Ls, dx,dy,dz);
  build_curl_curl(CC, dx,dy,dz);

  // export taps
  cStencil::cStencilData
    GxEx,GxEy,GxEz, GyEx,GyEy,GyEz, GzEx,GzEy,GzEz,
    LEx,LEy,LEz,
    CCxEx,CCxEy,CCxEz, CCyEx,CCyEy,CCyEz, CCzEx,CCzEy,CCzEz;

  G[0].Ex.ExportStencil(&GxEx); G[0].Ey.ExportStencil(&GxEy); G[0].Ez.ExportStencil(&GxEz);
  G[1].Ex.ExportStencil(&GyEx); G[1].Ey.ExportStencil(&GyEy); G[1].Ez.ExportStencil(&GyEz);
  G[2].Ex.ExportStencil(&GzEx); G[2].Ey.ExportStencil(&GzEy); G[2].Ez.ExportStencil(&GzEz);

  Ls.Ex.ExportStencil(&LEx); Ls.Ey.ExportStencil(&LEy); Ls.Ez.ExportStencil(&LEz);

  CC[0].Ex.ExportStencil(&CCxEx); CC[0].Ey.ExportStencil(&CCxEy); CC[0].Ez.ExportStencil(&CCxEz);
  CC[1].Ex.ExportStencil(&CCyEx); CC[1].Ey.ExportStencil(&CCyEy); CC[1].Ez.ExportStencil(&CCyEz);
  CC[2].Ex.ExportStencil(&CCzEx); CC[2].Ey.ExportStencil(&CCzEy); CC[2].Ez.ExportStencil(&CCzEz);

  // arrays for diagnostics
  std::vector<double> CCx(NT), CCy(NT), CCz(NT);          // numerical curl_curl
  std::vector<double> IDx(NT), IDy(NT), IDz(NT);          // GD - Lap
  std::vector<double> GDx(NT), GDy(NT), GDz(NT);          // numerical grad_div
  std::vector<double> LPx(NT), LPy(NT), LPz(NT);          // numerical laplacian
  std::vector<double> CCxA(NT), CCyA(NT), CCzA(NT);       // analytic curl_curl
  std::vector<double> GDAx(NT), GDAy(NT), GDAz(NT);       // analytic grad_div
  std::vector<double> LPAx(NT), LPAy(NT), LPAz(NT);       // analytic laplacian

  // error accumulators
  double linfA=0.0, l2A=0.0, l2refA=0.0;   // CC vs analytic
  double linfI=0.0, l2I=0.0, l2refI=0.0;   // CC vs ID
  double linfGD=0.0, l2GD=0.0, l2refGD=0.0;// GD vs analytic
  double linfLP=0.0, l2LP=0.0, l2refLP=0.0;// Lap vs analytic

  for (int k=0; k<Nz; ++k){
    const double z=k*dz;
    for (int j=0; j<Ny; ++j){
      const double y=j*dy;
      for (int i=0; i<Nx; ++i){
        const double x=i*dx;
        const size_t id = (size_t)k*Ny*Nx + (size_t)j*Nx + (size_t)i;

        // --- numerical curl_curl rows ---
        const double CCx_num =
            apply_exported(CCxEx, Ex,i,j,k,Nx,Ny,Nz)
          + apply_exported(CCxEy, Ey,i,j,k,Nx,Ny,Nz)
          + apply_exported(CCxEz, Ez,i,j,k,Nx,Ny,Nz);
        const double CCy_num =
            apply_exported(CCyEx, Ex,i,j,k,Nx,Ny,Nz)
          + apply_exported(CCyEy, Ey,i,j,k,Nx,Ny,Nz)
          + apply_exported(CCyEz, Ez,i,j,k,Nx,Ny,Nz);
        const double CCz_num =
            apply_exported(CCzEx, Ex,i,j,k,Nx,Ny,Nz)
          + apply_exported(CCzEy, Ey,i,j,k,Nx,Ny,Nz)
          + apply_exported(CCzEz, Ez,i,j,k,Nx,Ny,Nz);
        CCx[id]=CCx_num; CCy[id]=CCy_num; CCz[id]=CCz_num;

        // --- numerical grad_div rows ---
        const double GDx_num =
            apply_exported(GxEx, Ex,i,j,k,Nx,Ny,Nz)
          + apply_exported(GxEy, Ey,i,j,k,Nx,Ny,Nz)
          + apply_exported(GxEz, Ez,i,j,k,Nx,Ny,Nz);
        const double GDy_num =
            apply_exported(GyEx, Ex,i,j,k,Nx,Ny,Nz)
          + apply_exported(GyEy, Ey,i,j,k,Nx,Ny,Nz)
          + apply_exported(GyEz, Ez,i,j,k,Nx,Ny,Nz);
        const double GDz_num =
            apply_exported(GzEx, Ex,i,j,k,Nx,Ny,Nz)
          + apply_exported(GzEy, Ey,i,j,k,Nx,Ny,Nz)
          + apply_exported(GzEz, Ez,i,j,k,Nx,Ny,Nz);
        GDx[id]=GDx_num; GDy[id]=GDy_num; GDz[id]=GDz_num;

        // --- numerical laplacian component-wise ---
        const double LPx_num = apply_exported(LEx, Ex,i,j,k,Nx,Ny,Nz);
        const double LPy_num = apply_exported(LEy, Ey,i,j,k,Nx,Ny,Nz);
        const double LPz_num = apply_exported(LEz, Ez,i,j,k,Nx,Ny,Nz);
        LPx[id]=LPx_num; LPy[id]=LPy_num; LPz[id]=LPz_num;

        // --- numerical identity ---
        IDx[id] = GDx_num - LPx_num;
        IDy[id] = GDy_num - LPy_num;
        IDz[id] = GDz_num - LPz_num;

        // --- analytic references at (x,y,z) ---
        const Vec3 CCa  = analyticCurlCurlE(x,y,z,a,b,c);
        const Vec3 GDa  = analyticGradDivE (x,y,z,a,b,c);
        const Vec3 LPa  = analyticLaplacianE(x,y,z,a,b,c);
        CCxA[id]=CCa.x; CCyA[id]=CCa.y; CCzA[id]=CCa.z;
        GDAx[id]=GDa.x; GDAy[id]=GDa.y; GDAz[id]=GDa.z;
        LPAx[id]=LPa.x; LPAy[id]=LPa.y; LPAz[id]=LPa.z;

        // --- accumulate errors ---
        // CC vs analytic
        {
          const double ex = CCx_num-CCa.x, ey = CCy_num-CCa.y, ez = CCz_num-CCa.z;
          linfA = std::max(linfA, std::max(std::abs(ex), std::max(std::abs(ey), std::abs(ez))));
          l2A  += ex*ex + ey*ey + ez*ez;
          l2refA += CCa.x*CCa.x + CCa.y*CCa.y + CCa.z*CCa.z;
        }
        // CC vs ID (numerical identity)
        {
          const double ex = CCx_num-IDx[id], ey = CCy_num-IDy[id], ez = CCz_num-IDz[id];
          linfI = std::max(linfI, std::max(std::abs(ex), std::max(std::abs(ey), std::abs(ez))));
          l2I  += ex*ex + ey*ey + ez*ez;
          const double rx = IDx[id], ry = IDy[id], rz = IDz[id];
          l2refI += rx*rx + ry*ry + rz*rz;
        }
        // GD vs analytic
        {
          const double ex = GDx_num-GDa.x, ey = GDy_num-GDa.y, ez = GDz_num-GDa.z;
          linfGD = std::max(linfGD, std::max(std::abs(ex), std::max(std::abs(ey), std::abs(ez))));
          l2GD  += ex*ex + ey*ey + ez*ez;
          l2refGD += GDa.x*GDa.x + GDa.y*GDa.y + GDa.z*GDa.z;
        }
        // Lap vs analytic
        {
          const double ex = LPx_num-LPa.x, ey = LPy_num-LPa.y, ez = LPz_num-LPa.z;
          linfLP = std::max(linfLP, std::max(std::abs(ex), std::max(std::abs(ey), std::abs(ez))));
          l2LP  += ex*ex + ey*ey + ez*ez;
          l2refLP += LPa.x*LPa.x + LPa.y*LPa.y + LPa.z*LPa.z;
        }
      }
    }
  }

  // norms (relative L2 kept for per-variant diagnostics)
  const double relL2A  = (l2refA  > 0.0) ? std::sqrt(l2A  / l2refA ) : std::sqrt(l2A);
  const double relL2I  = (l2refI  > 0.0) ? std::sqrt(l2I  / l2refI ) : std::sqrt(l2I);
  const double relL2GD = (l2refGD > 0.0) ? std::sqrt(l2GD / l2refGD) : std::sqrt(l2GD);
  const double relL2LP = (l2refLP > 0.0) ? std::sqrt(l2LP / l2refLP) : std::sqrt(l2LP);

  // absolute L2 (RMS) of CC vs analytic over the whole grid (used for convergence table)
  const double l2absA = std::sqrt(l2A / double(NT));

  // print per-variant diagnostics (unchanged)
  std::cout << "  " << std::left << std::setw(12) << label
            << " CC vs analytic:   Linf=" << std::scientific << std::setprecision(3) << linfA
            << "   RelL2=" << relL2A << "\n";
  std::cout << "  " << std::left << std::setw(12) << ""
            << " CC vs (GD-Lap):   Linf=" << linfI
            << "   RelL2=" << relL2I << "\n";
  std::cout << "  " << std::left << std::setw(12) << ""
            << " GD vs analytic:   Linf=" << linfGD
            << "   RelL2=" << relL2GD << "\n";
  std::cout << "  " << std::left << std::setw(12) << ""
            << " Lap vs analytic:  Linf=" << linfLP
            << "   RelL2=" << relL2LP << "\n";
  std::cout << "    (RelL2 := ||error||_2 / ||reference||_2 over the whole grid)\n";

  // component-wise tables at one interior point
  print_point_curlcurl(N, L, CCx, CCy, CCz, CCxA, CCyA, CCzA, IDx, IDy, IDz, label);
  print_point_gd     (N, L, GDx, GDy, GDz, GDAx, GDAy, GDAz, label);
  print_point_lap    (N, L, LPx, LPy, LPz, LPAx, LPAy, LPAz, label);

  // return the main CC vs analytic stats for convergence summary
  ErrStats s; s.linf = linfA; s.l2abs = l2absA; s.l2rel = relL2A;
  return s;
}

} // anon

// ------------------------------ Registration ------------------------------
namespace CurlCurlE {

struct Variant {
  const char* name;
  std::function<void(cGradDivEStencil*, double, double, double)>   build_grad_div;
  std::function<void(cLaplacianStencil*, double, double, double)>  build_lap;
  std::function<void(cCurlCurlEStencil*, double, double, double)>  build_curl_curl;
};

int Run(const std::vector<std::string>&) {
  const double L = 1.0;
  const double a = 2.0*M_PI, b = 3.0*M_PI, c = 5.0*M_PI;

  const Variant variants[] = {
    { "2nd-compact",
      [](cGradDivEStencil G[3], double dx,double dy,double dz){ SecondOrder::InitGradDivEBStencils_compact(G,dx,dy,dz); },
      [](cLaplacianStencil* Ls, double dx,double dy,double dz){ SecondOrder::InitLaplacianStencil(Ls,dx,dy,dz); },
      [](cCurlCurlEStencil CC[3], double dx,double dy,double dz){ SecondOrder::InitCurlCurlEStencils_compact(CC,dx,dy,dz); }
    },
    { "2nd-wide",
      [](cGradDivEStencil G[3], double dx,double dy,double dz){ SecondOrder::InitGradDivEBStencils_wide(G,dx,dy,dz); },
      [](cLaplacianStencil* Ls, double dx,double dy,double dz){ SecondOrder::InitLaplacianStencil(Ls,dx,dy,dz); },
      [](cCurlCurlEStencil CC[3], double dx,double dy,double dz){ SecondOrder::InitCurlCurlEStencils_wide(CC,dx,dy,dz); }
    },
    { "4th",
      [](cGradDivEStencil G[3], double dx,double dy,double dz){ FourthOrder::InitGradDivEBStencils(G,dx,dy,dz); },
      [](cLaplacianStencil* Ls, double dx,double dy,double dz){ FourthOrder::InitLaplacianStencil(Ls,dx,dy,dz); },
      [](cCurlCurlEStencil CC[3], double dx,double dy,double dz){ FourthOrder::InitCurlCurlEStencils(CC,dx,dy,dz); }
    },
    { "6th",
      [](cGradDivEStencil G[3], double dx,double dy,double dz){ SixthOrder::InitGradDivEBStencils(G,dx,dy,dz); },
      [](cLaplacianStencil* Ls, double dx,double dy,double dz){ SixthOrder::InitLaplacianStencil(Ls,dx,dy,dz); },
      [](cCurlCurlEStencil CC[3], double dx,double dy,double dz){ SixthOrder::InitCurlCurlEStencils(CC,dx,dy,dz); }
    },
    { "8th",
      [](cGradDivEStencil G[3], double dx,double dy,double dz){ EighthOrder::InitGradDivEBStencils(G,dx,dy,dz); },
      [](cLaplacianStencil* Ls, double dx,double dy,double dz){ EighthOrder::InitLaplacianStencil(Ls,dx,dy,dz); },
      [](cCurlCurlEStencil CC[3], double dx,double dy,double dz){ EighthOrder::InitCurlCurlEStencils(CC,dx,dy,dz); }
    }
  };

  // For the convergence study we’ll use this fixed sweep.
  // (If you want CLI-driven N, we can parse args and override this vector.)
  std::vector<int> Ns = {16, 24, 32, 48, 64};

  std::cout << "\n=== Corner curl_curl(E): analytic vs numerical, identity, and GD/Lap verification ===\n"
            << "Domain L=" << L << ", wavenumbers a=2π, b=3π, c=5π\n";

  // --- Collect per-order errors (for CC vs analytic) to build the summary table at the end.
  struct Row { double linf=0.0, l2=0.0; };
  std::vector<Row> second_rows, fourth_rows, sixth_rows, eighth_rows;
  auto push_row = [](std::vector<Row>& vec, const ErrStats& s){
    Row r; r.linf = s.linf; r.l2 = s.l2abs; vec.push_back(r);
  };

  for (size_t t = 0; t < Ns.size(); ++t) {
    const int N = Ns[t];

    // Eye-catching separator per N (easier to scan logs)
    std::cout << "\n======== N = " << N << " ============================================================\n";

    // Run every variant; record 2nd-compact, 4th, 6th, 8th into the summary vectors.
    ErrStats s2c, s4, s6, s8;

    for (size_t v=0; v<sizeof(variants)/sizeof(variants[0]); ++v) {
      ErrStats res = run_one_order(variants[v].name,
                                   variants[v].build_grad_div,
                                   variants[v].build_lap,
                                   variants[v].build_curl_curl,
                                   N, L, a, b, c);

      // Map variants to the four summary columns.
      const std::string name = variants[v].name;
      if      (name == "2nd-compact") s2c = res;
      else if (name == "4th")         s4  = res;
      else if (name == "6th")         s6  = res;
      else if (name == "8th")         s8  = res;

      // Note: we still run "2nd-wide" for diagnostics, but we do not include it in the summary table.
    }

    push_row(second_rows, s2c);
    push_row(fourth_rows, s4);
    push_row(sixth_rows,  s6);
    push_row(eighth_rows, s8);

    std::cout << "=====================================================================================\n";
  }

  // ---------------- Convergence summary (L_inf & L2 RMS of CC vs analytic) ----------------
  auto safe_ord = [](double e_prev, double e_curr, int N_prev, int N_curr)->double{
    if (e_prev<=0.0 || e_curr<=0.0) return 0.0;
    const double rN = double(N_curr)/double(N_prev);
    return std::log(e_prev/e_curr) / std::log(rN);
  };

  auto print_row = [&](int idx){
    const int N = Ns[idx];
    const Row& s2 = second_rows[idx];
    const Row& s4 = fourth_rows[idx];
    const Row& s6 = sixth_rows[idx];
    const Row& s8 = eighth_rows[idx];

    // Observed orders relative to the previous N (h ~ 1/N)
    double o2_inf=0, o2_l2=0, o4_inf=0, o4_l2=0, o6_inf=0, o6_l2=0, o8_inf=0, o8_l2=0;
    if (idx>0) {
      o2_inf = safe_ord(second_rows[idx-1].linf, s2.linf, Ns[idx-1], N);
      o2_l2  = safe_ord(second_rows[idx-1].l2,   s2.l2,   Ns[idx-1], N);
      o4_inf = safe_ord(fourth_rows[idx-1].linf, s4.linf, Ns[idx-1], N);
      o4_l2  = safe_ord(fourth_rows[idx-1].l2,   s4.l2,   Ns[idx-1], N);
      o6_inf = safe_ord(sixth_rows[idx-1].linf,  s6.linf, Ns[idx-1], N);
      o6_l2  = safe_ord(sixth_rows[idx-1].l2,    s6.l2,   Ns[idx-1], N);
      o8_inf = safe_ord(eighth_rows[idx-1].linf, s8.linf, Ns[idx-1], N);
      o8_l2  = safe_ord(eighth_rows[idx-1].l2,   s8.l2,   Ns[idx-1], N);
    }

    std::cout << std::setw(5) << N << " | "
              << std::scientific << std::setprecision(3)
              << std::setw(12) << s2.linf  << " "
              << std::fixed << std::setprecision(2) << std::setw(4) << o2_inf << "   "
              << std::scientific << std::setprecision(3) << std::setw(12) << s2.l2 << "   " 
              << std::fixed << std::setprecision(2) << std::setw(5) << o2_l2 << " | "

              << std::scientific << std::setprecision(3) << std::setw(12) << s4.linf
              << std::fixed << std::setprecision(2) << std::setw(5) << o4_inf << "   "
              << std::scientific << std::setprecision(3) << std::setw(12) << s4.l2 << "   " 
              << std::fixed << std::setprecision(2) << std::setw(5) << o4_l2 << " | "

              << std::scientific << std::setprecision(3) << std::setw(13) << s6.linf
              << std::fixed << std::setprecision(2) << std::setw(5) << o6_inf << "   "
              << std::scientific << std::setprecision(3) << std::setw(12) << s6.l2 << "   " 
              << std::fixed << std::setprecision(2) << std::setw(5) << o6_l2 << " | "

              << std::scientific << std::setprecision(3) << std::setw(13) << s8.linf
              << std::fixed << std::setprecision(2) << std::setw(5) << o8_inf << "   "
              << std::scientific << std::setprecision(3) << std::setw(12) << s8.l2 << "   " 
              << std::fixed << std::setprecision(2) << std::setw(5) << o8_l2
              << " |\n";
  };

  // Header & table
  std::cout << "--------------------------------------------------------------------------------------------------------------\n"
            << "   N  |  Second (L_inf)  Ord   Second (L2)     Ord  |  Fourth (L_inf)  Ord   Fourth (L2)     Ord  |"
            << "   Sixth (L_inf)  Ord    Sixth (L2)     Ord  |  Eighth (L_inf)  Ord   Eighth (L2)     Ord\n"
            << "--------------------------------------------------------------------------------------------------------------\n";
  for (size_t i=0; i<Ns.size(); ++i) print_row((int)i);
  std::cout << "--------------------------------------------------------------------------------------------------------------\n";

  return 0;
}

} // namespace CurlCurlE

// Register in the harness (matches framework’s macro & force-link pattern)
REGISTER_STENCIL_TEST(CurlCurlE,
  "curl_curl_e",
  "Direct curl_curl(E) vs analytic and (grad_div−laplacian); also verify GD and Lap vs analytic.");

// Force-link shim (called from test_force_link_all.cpp)
namespace CurlCurlE {
  void ForceLinkAllTests() {}
}

