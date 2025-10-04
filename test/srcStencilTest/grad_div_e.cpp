/**
 * @file    grad_div_e.cpp
 * @brief   Convergence test for ∇(∇·E) stencils (SecondOrder/compact, SecondOrder/wide, FourthOrder),
 *          with optional Tecplot ASCII volume export (Ordered 3D, BLOCK).
 *
 * =============================================================================
 * 1) PURPOSE
 * =============================================================================
 * Verifies correctness and order of accuracy of discrete operators approximating
 *     G(E) = ∇(∇·E) ,  E=(Ex,Ey,Ez),
 * using three stencil families:
 *   (A) SecondOrder/Compact  : minimal centered stencils (3-pt D2, mixed via nested D1)
 *   (B) SecondOrder/Wide     : second-order face/edge-averaged construction
 *   (C) FourthOrder          : 5-point 4th-order D2; mixed via composition of 4th-order D1
 *
 * The program runs a mesh-refinement study (periodic domain) and reports L∞ and relative L²
 * errors with observed orders (log₂ error ratios). Optionally writes Tecplot volumes.
 *
 * =============================================================================
 * 2) ANALYTIC TEST (PERIODIC)
 * =============================================================================
 *   Ex = sin(ax) cos(by) cos(cz)
 *   Ey = cos(ax) sin(by) cos(cz)
 *   Ez = cos(ax) cos(by) sin(cz)
 *   div E = (a + b + c) cos(ax) cos(by) cos(cz)
 *   ∇(div E) =
 *     [ -a(a+b+c) sin(ax) cos(by) cos(cz),
 *       -b(a+b+c) cos(ax) sin(by) cos(cz),
 *       -c(a+b+c) cos(ax) cos(by) sin(cz) ].
 *
 * =============================================================================
 * 3) GRID & NUMERICS
 * =============================================================================
 * - Domain [0,Lx)×[0,Ly)×[0,Lz), uniform Nx×Ny×Nz (here cubic, Nx=Ny=Nz=N).
 * - Nodes at (i+½)dx etc. (periodic wrap; interior order preserved).
 * - Each operator exports 3 rows S[r], each with 3 column stencils acting on (Ex,Ey,Ez).
 * - We pre-export sparse taps **once per row** and apply them in the inner loops
 *   (faster and avoids const-qualification issues).
 *
 * =============================================================================
 * 4) OUTPUT
 * =============================================================================
 * - Table per operator:
 *     N    L_inf error  Order    L2_rel error  Order
 * - Optional single pointwise comparison at one interior node (first grid by default).
 * - Optional Tecplot file on the finest grid: variables
 *   x,y,z,Ex,Ey,Ez,Gx_num,Gy_num,Gz_num,Gx_ana,Gy_ana,Gz_ana,Err_mag
 *
 * =============================================================================
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <iomanip>
#include <iostream>
#include <fstream>

#include "pic.h"

using namespace PIC::FieldSolver::Electromagnetic::ECSIM::Stencil;

// ---------- Analytic field and ∇(∇·E) ----------
struct Vec3 { double x,y,z; };

static inline Vec3 analyticE(double x, double y, double z, double a, double b, double c) {
  return { std::sin(a*x)*std::cos(b*y)*std::cos(c*z),
           std::cos(a*x)*std::sin(b*y)*std::cos(c*z),
           std::cos(a*x)*std::cos(b*y)*std::sin(c*z) };
}

static inline Vec3 analyticGradDivE(double x, double y, double z, double a, double b, double c) {
  const double s = (a + b + c);
  return { -a*s*std::sin(a*x)*std::cos(b*y)*std::cos(c*z),
           -b*s*std::cos(a*x)*std::sin(b*y)*std::cos(c*z),
           -c*s*std::cos(a*x)*std::cos(b*y)*std::sin(c*z) };
}

// ---------- Utility: periodic wrap ----------
static inline int wrap(int idx, int N) {
  int r = idx % N; if (r < 0) r += N; return r;
}

// ---------- Apply a single exported stencil to a scalar field ----------
static inline double apply_exported(const cStencil::cStencilData& S,
                                    const std::vector<double>& F,
                                    int i, int j, int k, int Nx, int Ny, int Nz)
{
  double acc = 0.0;
  for (int n = 0; n < S.Length; ++n) {
    const int ii = wrap(i + S.Data[n].i, Nx);
    const int jj = wrap(j + S.Data[n].j, Ny);
    const int kk = wrap(k + S.Data[n].k, Nz);
    const size_t idx = (size_t)kk*Ny*Nx + (size_t)jj*Nx + (size_t)ii;
    acc += S.Data[n].a * F[idx];
  }
  return acc;
}

// ---------- TECplot export (Ordered 3D, BLOCK packing) ----------
static void write_tecplot_block(
    const std::string& filename,
    int Nx, int Ny, int Nz,
    double dx, double dy, double dz,
    const std::vector<double>& Ex,
    const std::vector<double>& Ey,
    const std::vector<double>& Ez,
    const std::vector<double>& Gx_num,
    const std::vector<double>& Gy_num,
    const std::vector<double>& Gz_num,
    const std::vector<double>& Gx_ana,
    const std::vector<double>& Gy_ana,
    const std::vector<double>& Gz_ana)
{
  std::ofstream out(filename.c_str());
  if (!out) {
    std::cerr << "ERROR: cannot open Tecplot file: " << filename << "\n";
    return;
  }

  out << "TITLE = \"GradDivE test\"\n";
  out << "VARIABLES = \"x\" \"y\" \"z\" "
         "\"Ex\" \"Ey\" \"Ez\" "
         "\"Gx_num\" \"Gy_num\" \"Gz_num\" "
         "\"Gx_ana\" \"Gy_ana\" \"Gz_ana\" "
         "\"Err_mag\"\n";
  out << "ZONE T=\"volume\", I=" << Nx << ", J=" << Ny << ", K=" << Nz
      << ", DATAPACKING=BLOCK\n";

  auto write_block = [&](auto get_value) {
    out.setf(std::ios::scientific);
    out << std::setprecision(8);
    for (int k=0; k<Nz; ++k)
      for (int j=0; j<Ny; ++j)
        for (int i=0; i<Nx; ++i)
          out << get_value(i,j,k) << "\n";
  };

  // x, y, z
  write_block([&](int i,int, int){ return (i + 0.5) * dx; });
  write_block([&](int, int j,int){ return (j + 0.5) * dy; });
  write_block([&](int, int, int k){ return (k + 0.5) * dz; });

  auto at = [&](const std::vector<double>& A, int i,int j,int k)->double {
    const size_t idx = (size_t)k*Ny*Nx + (size_t)j*Nx + (size_t)i;
    return A[idx];
  };

  // Ex, Ey, Ez
  write_block([&](int i,int j,int k){ return at(Ex,i,j,k); });
  write_block([&](int i,int j,int k){ return at(Ey,i,j,k); });
  write_block([&](int i,int j,int k){ return at(Ez,i,j,k); });

  // Gx_num, Gy_num, Gz_num
  write_block([&](int i,int j,int k){ return at(Gx_num,i,j,k); });
  write_block([&](int i,int j,int k){ return at(Gy_num,i,j,k); });
  write_block([&](int i,int j,int k){ return at(Gz_num,i,j,k); });

  // Gx_ana, Gy_ana, Gz_ana
  write_block([&](int i,int j,int k){ return at(Gx_ana,i,j,k); });
  write_block([&](int i,int j,int k){ return at(Gy_ana,i,j,k); });
  write_block([&](int i,int j,int k){ return at(Gz_ana,i,j,k); });

  // Err_mag = ||G_num - G_ana||
  write_block([&](int i,int j,int k){
    const double dx_ = at(Gx_num,i,j,k) - at(Gx_ana,i,j,k);
    const double dy_ = at(Gy_num,i,j,k) - at(Gy_ana,i,j,k);
    const double dz_ = at(Gz_num,i,j,k) - at(Gz_ana,i,j,k);
    return std::sqrt(dx_*dx_ + dy_*dy_ + dz_*dz_);
  });

  out.close();
  std::cout << "Tecplot: wrote " << filename << " ("
            << Nx << "x" << Ny << "x" << Nz << ")\n";
}

// ---------- Driver: evaluate one operator flavor and compute errors ----------
enum class OpFlavor { SecondCompact, SecondWide, FourthOrder };

static void build_stencil_rows(OpFlavor f, cGradDivEStencil S[3], double dx, double dy, double dz) {
  switch (f) {
    case OpFlavor::SecondCompact:
      SecondOrder::InitGradDivEBStencils_compact(S, dx, dy, dz);
      break;
    case OpFlavor::SecondWide:
      SecondOrder::InitGradDivEBStencils_wide(S, dx, dy, dz);
      break;
    case OpFlavor::FourthOrder:
      FourthOrder::InitGradDivEBStencils(S, dx, dy, dz);
      break;
  }
}

static std::string flavor_name(OpFlavor f) {
  switch (f) {
    case OpFlavor::SecondCompact: return "SecondOrder_Compact";
    case OpFlavor::SecondWide:    return "SecondOrder_Wide";
    case OpFlavor::FourthOrder:   return "FourthOrder";
  }
  return "Unknown";
}

struct ErrStats { double linf, l2; };

static ErrStats run_case(OpFlavor flavor, int N, double Lx, double Ly, double Lz,
                         double a, double b, double c,
                         // Optional array outputs for pointwise print or Tecplot
                         std::vector<double>* out_Gx = nullptr,
                         std::vector<double>* out_Gy = nullptr,
                         std::vector<double>* out_Gz = nullptr,
                         std::vector<double>* out_Gx_true = nullptr,
                         std::vector<double>* out_Gy_true = nullptr,
                         std::vector<double>* out_Gz_true = nullptr,
                         bool tecplot_dump = false)
{
  const int Nx = N, Ny = N, Nz = N;
  const double dx = Lx / Nx, dy = Ly / Ny, dz = Lz / Nz;

  // Build stencils
  cGradDivEStencil S[3];
  build_stencil_rows(flavor, S, dx, dy, dz);

  // Pre-export stencils once (per row, per component)
  cStencil::cStencilData SX[3], SY[3], SZ[3];
  for (int r = 0; r < 3; ++r) {
    S[r].Ex.ExportStencil(&SX[r]);
    S[r].Ey.ExportStencil(&SY[r]);
    S[r].Ez.ExportStencil(&SZ[r]);
  }

  // Sample E and analytic G on the grid
  std::vector<double> Ex(Nx*Ny*Nz), Ey(Nx*Ny*Nz), Ez(Nx*Ny*Nz);
  std::vector<double> Gx_true(Nx*Ny*Nz), Gy_true(Nx*Ny*Nz), Gz_true(Nx*Ny*Nz);

  for (int k=0; k<Nz; ++k) {
    const double z = (k + 0.5) * dz;
    for (int j=0; j<Ny; ++j) {
      const double y = (j + 0.5) * dy;
      for (int i=0; i<Nx; ++i) {
        const double x = (i + 0.5) * dx;
        const size_t idx = (size_t)k*Ny*Nx + (size_t)j*Nx + (size_t)i;

        const Vec3 E = analyticE(x,y,z, a,b,c);
        Ex[idx] = E.x; Ey[idx] = E.y; Ez[idx] = E.z;

        const Vec3 G = analyticGradDivE(x,y,z, a,b,c);
        Gx_true[idx] = G.x; Gy_true[idx] = G.y; Gz_true[idx] = G.z;
      }
    }
  }

  // Apply GradDiv rows using exported taps
  std::vector<double> Gx(Nx*Ny*Nz), Gy(Nx*Ny*Nz), Gz(Nx*Ny*Nz);
  for (int k=0; k<Nz; ++k)
    for (int j=0; j<Ny; ++j)
      for (int i=0; i<Nx; ++i) {
        const size_t idx = (size_t)k*Ny*Nx + (size_t)j*Nx + (size_t)i;

        Gx[idx] = apply_exported(SX[0], Ex, i,j,k, Nx,Ny,Nz)
                + apply_exported(SY[0], Ey, i,j,k, Nx,Ny,Nz)
                + apply_exported(SZ[0], Ez, i,j,k, Nx,Ny,Nz);

        Gy[idx] = apply_exported(SX[1], Ex, i,j,k, Nx,Ny,Nz)
                + apply_exported(SY[1], Ey, i,j,k, Nx,Ny,Nz)
                + apply_exported(SZ[1], Ez, i,j,k, Nx,Ny,Nz);

        Gz[idx] = apply_exported(SX[2], Ex, i,j,k, Nx,Ny,Nz)
                + apply_exported(SY[2], Ey, i,j,k, Nx,Ny,Nz)
                + apply_exported(SZ[2], Ez, i,j,k, Nx,Ny,Nz);
      }

  // Optional arrays to caller
  if (out_Gx)      *out_Gx      = Gx;
  if (out_Gy)      *out_Gy      = Gy;
  if (out_Gz)      *out_Gz      = Gz;
  if (out_Gx_true) *out_Gx_true = Gx_true;
  if (out_Gy_true) *out_Gy_true = Gy_true;
  if (out_Gz_true) *out_Gz_true = Gz_true;

  // Error norms
  double linf = 0.0, l2_num = 0.0, l2_den = 0.0;
  for (size_t idx=0; idx<Gx.size(); ++idx) {
    const double ex = std::abs(Gx[idx]-Gx_true[idx]);
    const double ey = std::abs(Gy[idx]-Gy_true[idx]);
    const double ez = std::abs(Gz[idx]-Gz_true[idx]);
    linf = std::max(linf, std::max(ex, std::max(ey, ez)));

    l2_num += (Gx[idx]-Gx_true[idx])*(Gx[idx]-Gx_true[idx])
            + (Gy[idx]-Gy_true[idx])*(Gy[idx]-Gy_true[idx])
            + (Gz[idx]-Gz_true[idx])*(Gz[idx]-Gz_true[idx]);
    l2_den += Gx_true[idx]*Gx_true[idx]
            + Gy_true[idx]*Gy_true[idx]
            + Gz_true[idx]*Gz_true[idx];
  }

  // Tecplot export (optional)
  if (tecplot_dump) {
    const std::string fname = "graddivE_" + flavor_name(flavor) + "_N" + std::to_string(N) + ".dat";
    write_tecplot_block(fname, Nx,Ny,Nz, dx,dy,dz,
                        Ex,Ey,Ez, Gx,Gy,Gz, Gx_true,Gy_true,Gz_true);
  }

  return ErrStats{ linf, (l2_den>0) ? std::sqrt(l2_num/l2_den) : std::sqrt(l2_num) };
}

// ---------- Single pointwise comparison (utility) ----------
static void print_point_comparison(int N, double Lx,double Ly,double Lz,
                                   double a,double b,double c,
                                   const std::vector<double>& Gx,
                                   const std::vector<double>& Gy,
                                   const std::vector<double>& Gz,
                                   const std::vector<double>& Gx_true,
                                   const std::vector<double>& Gy_true,
                                   const std::vector<double>& Gz_true,
                                   const char* flavor_label)
{
  const int Nx=N, Ny=N, Nz=N;
  const double dx=Lx/Nx, dy=Ly/Ny, dz=Lz/Nz;
  const int ii=Nx/2, jj=Ny/2, kk=Nz/2;
  const size_t idx=(size_t)kk*Ny*Nx + (size_t)jj*Nx + (size_t)ii;
  const double x=(ii+0.5)*dx, y=(jj+0.5)*dy, z=(kk+0.5)*dz;

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

  line("(GradDivE)_x", Gx_true[idx], Gx[idx]);
  line("(GradDivE)_y", Gy_true[idx], Gy[idx]);
  line("(GradDivE)_z", Gz_true[idx], Gz[idx]);
}

// ---------- Combined convergence table across all flavors ----------
// ---------- Combined convergence table across all flavors ----------
// Computes errors for all Ns and all stencil flavors, prints the per-N,
// per-stencil component-wise blocks inline (so you can see progress),
// and then prints the COMBINED TABLE ONCE at the very end.
static void print_convergence_combined(const std::vector<int>& Ns,
                                       double Lx, double Ly, double Lz,
                                       double a, double b, double c,
                                       bool tecplot_on_finest = true)
{
  const OpFlavor flavors[3] = {
    OpFlavor::SecondCompact,
    OpFlavor::SecondWide,
    OpFlavor::FourthOrder
  };
  const char* headers[3] = {
    "Second/Compact",
    "Second/Wide",
    "Fourth"
  };

  struct Row { double linf = 0.0, l2 = 0.0; };
  // Buffer all errors; err[f][t] corresponds to flavor f at Ns[t]
  std::vector<Row> err[3];
  for (int f=0; f<3; ++f) err[f].resize(Ns.size());

  std::cout << "\n=== GradDivE Convergence (Combined run) ===\n"
            << "Domain Lx=" << Lx << ", Ly=" << Ly << ", Lz=" << Lz
            << ", wave numbers a=" << a << ", b=" << b << ", c=" << c << "\n";

  // --- Compute everything first (and print component-wise blocks as requested) ---
  for (size_t t=0; t<Ns.size(); ++t) {
    const int N = Ns[t];
    const bool dump_tp = tecplot_on_finest && (t == Ns.size()-1);

    // For each flavor, compute errors + capture arrays for pointwise print
    struct Arrs {
      std::vector<double> Gx,Gy,Gz, GxT,GyT,GzT;
    } out[3];

    for (int f=0; f<3; ++f) {
      ErrStats s = run_case(flavors[f], N, Lx,Ly,Lz, a,b,c,
                            &out[f].Gx, &out[f].Gy, &out[f].Gz,
                            &out[f].GxT,&out[f].GyT,&out[f].GzT,
                            dump_tp);
      err[f][t] = Row{s.linf, s.l2};
    }

    // After finishing this N for all flavors, print each component-wise block
    for (int f=0; f<3; ++f) {
      print_point_comparison(N, Lx,Ly,Lz, a,b,c,
                             out[f].Gx, out[f].Gy, out[f].Gz,
                             out[f].GxT,out[f].GyT,out[f].GzT,
                             headers[f]);
    }
  }

  // --- Now print the combined table ONCE at the end ---
  std::cout << "\n---------------------------------------------------------------------------------------------------------------------------------------------------------\n";
  std::cout << "   N  |     " << std::setw(14) << headers[0] << " (L_inf)  Ord  "
                    << std::setw(14) << headers[0] << " (L2)     Ord  |  "
                    << std::setw(14) << headers[1] << " (L_inf)  Ord  "
                    << std::setw(14) << headers[1] << " (L2)     Ord  |  "
                    << std::setw(14) << headers[2] << " (L_inf)  Ord  "
                    << std::setw(14) << headers[2] << " (L2)     Ord\n";
  std::cout << "---------------------------------------------------------------------------------------------------------------------------------------------------------\n";

  for (size_t t=0; t<Ns.size(); ++t) {
    std::cout << std::setw(5) << Ns[t] << " ";
    for (int f=0; f<3; ++f) {
      double p_inf = 0.0, p_l2 = 0.0;
      if (t > 0) {
        // per-flavor observed orders (ratio vs previous N for the same flavor)
        p_inf = std::log2(err[f][t-1].linf / err[f][t].linf);
        p_l2  = std::log2(err[f][t-1].l2   / err[f][t].l2);
      }
      std::cout << " | "
                << std::scientific << std::setprecision(3) << std::setw(12) << err[f][t].linf << " "
                << std::fixed      << std::setprecision(2) << std::setw(5)  << (t==0 ? 0.0 : p_inf) << "  "
                << std::scientific << std::setprecision(3) << std::setw(12) << err[f][t].l2   << " "
                << std::fixed      << std::setprecision(2) << std::setw(5)  << (t==0 ? 0.0 : p_l2);
    }
    std::cout << "\n";
  }

  std::cout << "---------------------------------------------------------------------------------------------------------------------------------------------------------\n";
}

int test_grad_div_e() {
  // Domain and wavenumbers (periodic)
  const double Lx = 1.0, Ly = 1.0, Lz = 1.0;
  const double a = 2.0*M_PI, b = 2.0*M_PI, c = 2.0*M_PI;

  // Grid sizes (ensure compatibility with 4th-order 5-point stencils)
  std::vector<int> Ns = {16, 24, 32, 48, 64};

  // Combined table; now prints component-wise blocks for each stencil and each N.
  print_convergence_combined(Ns, Lx,Ly,Lz, a,b,c, /*tecplot_on_finest=*/true);

  std::cout << "\nTecplot files (finest grid per operator):\n"
               "  graddivE_SecondOrder_Compact_N64.dat\n"
               "  graddivE_SecondOrder_Wide_N64.dat\n"
               "  graddivE_FourthOrder_N64.dat\n";
  return 0;
}

