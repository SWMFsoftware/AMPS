/***************************************************************************************
 * curl_curl_e.cpp  (C++11/C++14 compatible)
 * -------------------------------------------------------------------------------------
 * PURPOSE
 *   Validation & convergence tests for the discrete operator
 *
 *       curl_curl(E) = ∇×(∇×E) = ∇(∇·E) − ∇²E
 *
 *   evaluated at CORNER nodes (collocated stencils) at 2nd/4th/6th/8th order.
 *   For second order, both grad_div variants are tested:
 *     • SecondOrder::InitCurlCurlEStencils_compact
 *     • SecondOrder::InitCurlCurlEStencils_wide
 *
 * TESTS
 *   • Accuracy vs analytic curl_curl(E) on a smooth periodic field.
 *   • Convergence vs N; also prints a component-wise comparison at one interior point.
 *
 * FIELD & DOMAIN
 *   • Periodic cube [0,L]^3 with N×N×N CORNER samples (wrap indexing).
 *   • E(x,y,z) = { sin(ax) cos(by) cos(cz),  cos(ax) sin(by) cos(cz),  cos(ax) cos(by) sin(cz) },
 *     with a=2π, b=3π, c=5π (unequal to avoid trivial cancellations).
 *
 * USAGE
 *   ./stencil_tests curl_curl_e
 *   Runs 2nd-compact, 2nd-wide, 4th, 6th, 8th; prints Linf/RelL2 vs analytic.
 ***************************************************************************************/

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <utility>

#include "test_harness.h"
#include "test_register.h"
#include "test_force_link_all.h"
#include "pic.h"  // cStencil, cGradDivEStencil, cLaplacianStencil, cCurlCurlEStencil, ExportStencil

using namespace PIC::FieldSolver::Electromagnetic::ECSIM::Stencil;

namespace {

// ---------------- small helpers ----------------
inline int wrap(int i, int N){ int r = i % N; return (r<0)? r+N : r; }
struct Vec3 { double x,y,z; };

// Analytic field at CORNERS
inline Vec3 analyticE(double x, double y, double z, double a, double b, double c){
  Vec3 v;
  v.x = std::sin(a*x)*std::cos(b*y)*std::cos(c*z);
  v.y = std::cos(a*x)*std::sin(b*y)*std::cos(c*z);
  v.z = std::cos(a*x)*std::cos(b*y)*std::sin(c*z);
  return v;
}

// Analytic curl_curl(E) via identity: ∇×∇×E = ∇(∇·E) − ∇²E
inline Vec3 analyticCurlCurlE(double x,double y,double z,double a,double b,double c){
  const double sx = std::sin(a*x), cx = std::cos(a*x);
  const double sy = std::sin(b*y), cy = std::cos(b*y);
  const double sz = std::sin(c*z), cz = std::cos(c*z);

  const double Ex = sx*cy*cz;
  const double Ey = cx*sy*cz;
  const double Ez = cx*cy*sz;

  // divE = (a+b+c) * cx*cy*cz
  const double K = (a + b + c);

  // grad(divE): divE = K * cx*cy*cz
  const double dDiv_dx = K * (-a*sx)*cy*cz;
  const double dDiv_dy = K * cx*(-b*sy)*cz;
  const double dDiv_dz = K * cx*cy*(-c*sz);

  // Laplacian: λ = a^2+b^2+c^2 (separable eigenfunction)
  const double lam = (a*a + b*b + c*c);
  const double lapx = -lam * Ex;
  const double lapy = -lam * Ey;
  const double lapz = -lam * Ez;

  Vec3 out;
  out.x = dDiv_dx - lapx;
  out.y = dDiv_dy - lapy;
  out.z = dDiv_dz - lapz;
  return out;
}

// Apply an exported (integerized) stencil to scalar field F at corner (i,j,k)
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

// Pretty print component-wise comparison
static void print_point_comparison(int N, double L,
                                   const std::vector<double>& CCx,
                                   const std::vector<double>& CCy,
                                   const std::vector<double>& CCz,
                                   const std::vector<double>& CCxA,
                                   const std::vector<double>& CCyA,
                                   const std::vector<double>& CCzA,
                                   const char* flavor_label)
{
  const double dx = L/N;
  // interior point with nontrivial values
  const int ii = N/4, jj = N/3, kk = (2*N)/5;
  const int Nx=N, Ny=N, Nz=N;
  const size_t idx = (size_t)kk*Ny*Nx + (size_t)jj*Nx + (size_t)ii;
  const double x = ii*dx, y = jj*dx, z = kk*dx;

  auto line = [&](const char* name, double a_, double n_){
    double err = std::abs(n_ - a_);
    std::cout << "    " << std::left << std::setw(14) << name
              << std::right << std::scientific << std::setprecision(8)
              << std::setw(16) << a_ << "   "
              << std::setw(16) << n_ << "   "
              << std::setprecision(3) << std::setw(8) << err << "\n";
  };

  std::cout << "\n[" << flavor_label << "] Component-wise comparison at interior point:\n"
            << "  Grid: N=" << N << ", (i,j,k)=(" << ii << "," << jj << "," << kk << ")"
            << ", (x,y,z)=(" << std::fixed << std::setprecision(6)
            << x << ", " << y << ", " << z << ")\n"
            << "  -----------------------------------------------------------------------------\n"
            << "    Component         Analytic               Numerical               AbsErr\n"
            << "  -----------------------------------------------------------------------------\n";

  line("(CurlCurlE)_x", CCxA[idx], CCx[idx]);
  line("(CurlCurlE)_y", CCyA[idx], CCy[idx]);
  line("(CurlCurlE)_z", CCzA[idx], CCz[idx]);
}

struct ErrStats { double linf=0.0, l2rel=0.0; };

// Generic runner for one order/variant
template<typename BuildGradDiv, typename BuildLap, typename BuildCurlCurl>
ErrStats run_one_order(const char* label,
                   BuildGradDiv build_grad_div,
                   BuildLap     build_lap,
                   BuildCurlCurl build_curl_curl,
                   int N, double L,
                   double a, double b, double c,
                   std::vector<double>* out_CCx = nullptr,
                   std::vector<double>* out_CCy = nullptr,
                   std::vector<double>* out_CCz = nullptr,
                   std::vector<double>* out_CCxA = nullptr,
                   std::vector<double>* out_CCyA = nullptr,
                   std::vector<double>* out_CCzA = nullptr)
{
  const double dx=L/N, dy=L/N, dz=L/N;
  const int Nx=N, Ny=N, Nz=N;
  const size_t NT = (size_t)Nx*Ny*Nz;

  // initialize E at CORNERS
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

  // export taps like curl_b.cpp does
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

  // numeric curl_curl(E) at corners
  std::vector<double> CCx(NT), CCy(NT), CCz(NT);

  // accumulate errors vs analytic
  double linfA=0.0, l2A=0.0, l2refA=0.0;

  for (int k=0; k<Nz; ++k){
    const double z=k*dz;
    for (int j=0; j<Ny; ++j){
      const double y=j*dy;
      for (int i=0; i<Nx; ++i){
        const double x=i*dx;
        const size_t id = (size_t)k*Ny*Nx + (size_t)j*Nx + (size_t)i;

        // direct CC rows
        const double vx =
            apply_exported(CCxEx, Ex,i,j,k,Nx,Ny,Nz)
          + apply_exported(CCxEy, Ey,i,j,k,Nx,Ny,Nz)
          + apply_exported(CCxEz, Ez,i,j,k,Nx,Ny,Nz);
        const double vy =
            apply_exported(CCyEx, Ex,i,j,k,Nx,Ny,Nz)
          + apply_exported(CCyEy, Ey,i,j,k,Nx,Ny,Nz)
          + apply_exported(CCyEz, Ez,i,j,k,Nx,Ny,Nz);
        const double vz =
            apply_exported(CCzEx, Ex,i,j,k,Nx,Ny,Nz)
          + apply_exported(CCzEy, Ey,i,j,k,Nx,Ny,Nz)
          + apply_exported(CCzEz, Ez,i,j,k,Nx,Ny,Nz);

        CCx[id]=vx; CCy[id]=vy; CCz[id]=vz;

        // analytic
        const Vec3 CCa = analyticCurlCurlE(x,y,z,a,b,c);
        const double eax=vx-CCa.x, eay=vy-CCa.y, eaz=vz-CCa.z;
        linfA = std::max(linfA, std::max(std::abs(eax), std::max(std::abs(eay), std::abs(eaz))));
        l2A  += eax*eax + eay*eay + eaz*eaz;
        l2refA += CCa.x*CCa.x + CCa.y*CCa.y + CCa.z*CCa.z;
      }
    }
  }

  // analytic arrays for the pretty printer
  std::vector<double> CCxA(NT), CCyA(NT), CCzA(NT);
  for (int k=0; k<Nz; ++k)
    for (int j=0; j<Ny; ++j)
      for (int i=0; i<Nx; ++i) {
        const size_t id = (size_t)k*Ny*Nx + (size_t)j*Nx + (size_t)i;
        const Vec3 CCa = analyticCurlCurlE(i*dx, j*dy, k*dz, a,b,c);
        CCxA[id]=CCa.x; CCyA[id]=CCa.y; CCzA[id]=CCa.z;
      }

  if (out_CCx)  *out_CCx  = CCx;
  if (out_CCy)  *out_CCy  = CCy;
  if (out_CCz)  *out_CCz  = CCz;
  if (out_CCxA) *out_CCxA = CCxA;
  if (out_CCyA) *out_CCyA = CCyA;
  if (out_CCzA) *out_CCzA = CCzA;

  // norms
  const double relL2 = (l2refA>0.0)? std::sqrt(l2A/l2refA) : std::sqrt(l2A);

  ErrStats s;
  s.linf  = linfA;
  s.l2rel = relL2;
  (void)label; // label printed by driver
  return s;
}

} // anon

// ------------------------------ Registration ------------------------------
namespace CurlCurlE {

struct Variant { 
  const char* name;
  std::function<void(cGradDivEStencil*, double, double, double)> build_grad_div;
  std::function<void(cLaplacianStencil*, double, double, double)> build_lap;
  std::function<void(cCurlCurlEStencil*, double, double, double)> build_curl_curl;
};

int Run(const std::vector<std::string>&) {
  const double L = 1.0;
  // Unequal wavenumbers to avoid trivial cancellations (mirrors curl_b.cpp)
  const double a = 2.0*M_PI;
  const double b = 3.0*M_PI; 
  const double c = 5.0*M_PI;

  // Define all variants
  const Variant variants[] = {
    {
      "2nd-compact",
      [](cGradDivEStencil G[3], double dx, double dy, double dz) {
        SecondOrder::InitGradDivEBStencils_compact(G, dx, dy, dz);
      },
      [](cLaplacianStencil* Ls, double dx, double dy, double dz) {
        SecondOrder::InitLaplacianStencil(Ls, dx, dy, dz);
      },
      [](cCurlCurlEStencil CC[3], double dx, double dy, double dz) {
        SecondOrder::InitCurlCurlEStencils_compact(CC, dx, dy, dz);
      }
    },
    {
      "2nd-wide",
      [](cGradDivEStencil G[3], double dx, double dy, double dz) {
        SecondOrder::InitGradDivEBStencils_wide(G, dx, dy, dz);
      },
      [](cLaplacianStencil* Ls, double dx, double dy, double dz) {
        SecondOrder::InitLaplacianStencil(Ls, dx, dy, dz);
      },
      [](cCurlCurlEStencil CC[3], double dx, double dy, double dz) {
        SecondOrder::InitCurlCurlEStencils_wide(CC, dx, dy, dz);
      }
    },
    {
      "4th",
      [](cGradDivEStencil G[3], double dx, double dy, double dz) {
        FourthOrder::InitGradDivEBStencils(G, dx, dy, dz);
      },
      [](cLaplacianStencil* Ls, double dx, double dy, double dz) {
        FourthOrder::InitLaplacianStencil(Ls, dx, dy, dz);
      },
      [](cCurlCurlEStencil CC[3], double dx, double dy, double dz) {
        FourthOrder::InitCurlCurlEStencils(CC, dx, dy, dz);
      }
    },
    {
      "6th",
      [](cGradDivEStencil G[3], double dx, double dy, double dz) {
        SixthOrder::InitGradDivEBStencils(G, dx, dy, dz);
      },
      [](cLaplacianStencil* Ls, double dx, double dy, double dz) {
        SixthOrder::InitLaplacianStencil(Ls, dx, dy, dz);
      },
      [](cCurlCurlEStencil CC[3], double dx, double dy, double dz) {
        SixthOrder::InitCurlCurlEStencils(CC, dx, dy, dz);
      }
    },
    {
      "8th",
      [](cGradDivEStencil G[3], double dx, double dy, double dz) {
        EighthOrder::InitGradDivEBStencils(G, dx, dy, dz);
      },
      [](cLaplacianStencil* Ls, double dx, double dy, double dz) {
        EighthOrder::InitLaplacianStencil(Ls, dx, dy, dz);
      },
      [](cCurlCurlEStencil CC[3], double dx, double dy, double dz) {
        EighthOrder::InitCurlCurlEStencils(CC, dx, dy, dz);
      }
    }
  };

  const int num_variants = sizeof(variants) / sizeof(variants[0]);
  std::vector<int> Ns = {16, 24, 32, 48, 64};
  
  struct Row { double linf, l2; };
  std::vector<Row> err[5];  // up to 5 variants
  for (int f = 0; f < num_variants; ++f) {
    err[f].resize(Ns.size());
  }

  std::cout << "\n=== Corner curl_curl(E) Convergence (2nd: compact/wide, 4th, 6th, 8th) ===\n"
            << "Domain: L=" << L << " (cube), wavenumbers: a=2π, b=3π, c=5π\n";

  // Run all variants and resolutions
  for (size_t t = 0; t < Ns.size(); ++t) {
    const int current_N = Ns[t];
    std::cout << "\nN = " << current_N << "\n";
    
    for (int f = 0; f < num_variants; ++f) {
      std::vector<double> CCx, CCy, CCz, CCxA, CCyA, CCzA;
      
      ErrStats s = run_one_order(
        variants[f].name,
        variants[f].build_grad_div,
        variants[f].build_lap,
        variants[f].build_curl_curl,
        current_N, L, a, b, c,
        &CCx, &CCy, &CCz, &CCxA, &CCyA, &CCzA
      );
      
      err[f][t] = Row{s.linf, s.l2rel};
      
      std::cout << "  " << std::left << std::setw(12) << variants[f].name
                << " Linf=" << std::scientific << std::setprecision(3) << s.linf
                << "   L2=" << s.l2rel << std::fixed << "\n";
      
      print_point_comparison(current_N, L, CCx, CCy, CCz, CCxA, CCyA, CCzA, variants[f].name);
    }
  }

  // Combined convergence table (format mirrors curl_b.cpp)
  std::cout << "\n" << std::string(150, '-') << "\n";
  std::cout << "   N  ";
  for (int f = 0; f < num_variants; ++f) {
    std::cout << " | " << std::setw(12) << variants[f].name << " (L_inf)  Ord   "
              << std::setw(12) << variants[f].name << " (L2)     Ord";
  }
  std::cout << "\n" << std::string(150, '-') << "\n";

  for (size_t t = 0; t < Ns.size(); ++t) {
    std::cout << std::setw(5) << Ns[t] << " ";
    for (int f = 0; f < num_variants; ++f) {
      double p_inf = 0.0, p_l2 = 0.0;
      if (t > 0) {
        const double rN = double(Ns[t]) / double(Ns[t-1]);
        const double den = std::log(rN);
        const double num_inf = err[f][t-1].linf / err[f][t].linf;
        const double num_l2  = err[f][t-1].l2  / err[f][t].l2;
        if (den > 0 && num_inf > 0) p_inf = std::log(num_inf) / den;
        if (den > 0 && num_l2  > 0) p_l2  = std::log(num_l2 ) / den;
      }
      std::cout << " | "
        << std::scientific << std::setprecision(3) << std::setw(12) << err[f][t].linf << " "
        << std::fixed << std::setprecision(2) << std::setw(5) << (t ? p_inf : 0.0) << "  "
        << std::scientific << std::setprecision(3) << std::setw(12) << err[f][t].l2 << " "
        << std::fixed << std::setprecision(2) << std::setw(5) << (t ? p_l2 : 0.0);
    }
    std::cout << "\n";
  }
  std::cout << std::string(150, '-') << "\n";

  return 0;
}

} // namespace CurlCurlE

// Register in the harness (matches framework’s macro & force-link pattern)
REGISTER_STENCIL_TEST(CurlCurlE,
  "curl_curl_e",
  "Build curl_curl(E) at 2nd (compact & wide) / 4th / 6th / 8th; compare with analytic and convergence.");

// Force-link shim (called from test_force_link_all.cpp)
namespace CurlCurlE {
  void ForceLinkAllTests() {}
}

