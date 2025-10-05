/**
 * @file    curl_b.cpp
 * @brief   Convergence test for corner curl(B) stencils:
 *          2nd-order (edge-based), 2nd-order (face-based), 4th-order, 6th-order.
 *
 * Analytic B (periodic, smooth):
 *   Bx =  sin(ax) cos(by) cos(cz)
 *   By =  cos(ax) sin(by) cos(cz)
 *   Bz =  cos(ax) cos(by) sin(cz)
 *
 * Curl(B) (continuous):
 *   (curl B)_x =  ∂Bz/∂y − ∂By/∂z
 *   (curl B)_y =  ∂Bx/∂z − ∂Bz/∂x
 *   (curl B)_z =  ∂By/∂x − ∂Bx/∂y
 *
 * Discrete target locations:
 *   • B components are sampled at cell centers: (i+1/2, j+1/2, k+1/2).
 *   • curl(B) is evaluated at cell corners:     (i,     j,     k    ).
 *     All stencils below are constructed for the (-1/2,-1/2,-1/2) corner relative
 *     to the base cell (consistent with AMPS ECSIM staggering). Periodic wrap is used.
 *
 * Output:
 *   • L∞ and relative L2 errors over a refinement set N in {16, 24, 32, 48, 64}.
 *   • Tecplot file on the finest grid for each variant.
 *   • Component-wise interior-point comparison (numeric vs analytic), like GradDivE.
 *
 * IMPORTANT: Avoid a=b=c for this analytic field.
 *   If a=b=c, the field is a gradient of ψ=-(1/a)cos(ax)cos(by)cos(cz) ⇒ curl(B)=0 identically.
 *   We therefore use a=2π, b=3π, c=5π so the analytic curl is non-zero and the comparison is meaningful.
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "test_harness.h"
#include "test_force_link_all.h"
#include "pic.h"                  // cStencil, cCurlBStencil, ExportStencil API

namespace EMSt = PIC::FieldSolver::Electromagnetic::ECSIM::Stencil;

// -------------------- utilities --------------------
static inline int wrap(int idx, int N) { int r = idx % N; if (r < 0) r += N; return r; }

struct V3 { double x,y,z; };

static inline V3 analyticB(double x, double y, double z, double a, double b, double c) {
  return {
    std::sin(a*x)*std::cos(b*y)*std::cos(c*z),
    std::cos(a*x)*std::sin(b*y)*std::cos(c*z),
    std::cos(a*x)*std::cos(b*y)*std::sin(c*z)
  };
}

static inline V3 analyticCurlB(double x, double y, double z, double a, double b, double c) {
  const double sxa = std::sin(a*x), cxa = std::cos(a*x);
  const double syb = std::sin(b*y), cyb = std::cos(b*y);
  const double szc = std::sin(c*z), czc = std::cos(c*z);

  const double dBz_dy = -b * cxa * syb * szc;
  const double dBy_dz = -c * cxa * syb * szc;
  const double dBx_dz = -c * sxa * cyb * szc;
  const double dBz_dx = -a * sxa * cyb * szc;
  const double dBy_dx = -a * sxa * syb * czc;
  const double dBx_dy = -b * sxa * syb * czc;

  return { dBz_dy - dBy_dz, dBx_dz - dBz_dx, dBy_dx - dBx_dy };
}

// apply an exported (integerized) stencil to a scalar field F at corner (i,j,k)
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

// Pretty print: component-wise numeric vs analytic at a central interior point
static void print_point_comparison(int N, double Lx, double Ly, double Lz,
                                   const std::vector<double>& Cx,
                                   const std::vector<double>& Cy,
                                   const std::vector<double>& Cz,
                                   const std::vector<double>& CxA,
                                   const std::vector<double>& CyA,
                                   const std::vector<double>& CzA,
                                   const char* flavor_label)
{
  const int Nx=N, Ny=N, Nz=N;
  const double dx=Lx/Nx, dy=Ly/Ny, dz=Lz/Nz;
  const int ii=Nx/2, jj=Ny/2, kk=Nz/2;  // a central interior point
  const size_t idx=(size_t)kk*Ny*Nx + (size_t)jj*Nx + (size_t)ii;
  const double x=(ii+0.0)*dx, y=(jj+0.0)*dy, z=(kk+0.0)*dz; // corners at integers

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

  line("(CurlB)_x", CxA[idx], Cx[idx]);
  line("(CurlB)_y", CyA[idx], Cy[idx]);
  line("(CurlB)_z", CzA[idx], Cz[idx]);
}

// -------------------- test driver --------------------
namespace CurlB {

using InitFn = void (*)(EMSt::cCurlBStencil*, double,double,double);
struct Variant { const char* name; InitFn init; };

static const Variant kVariants[] = {
  { "2nd-order (edge-based)",
    EMSt::SecondOrder::InitCurlBStencils_edge_based },
  { "2nd-order (face-based)",
    EMSt::SecondOrder::InitCurlBStencils_face_based },
  { "4th-order",
    EMSt::FourthOrder::InitCurlBStencils },
  { "6th-order",
    EMSt::SixthOrder::InitCurlBStencils }
};

struct ErrStats { double linf=0.0, l2rel=0.0; };

static ErrStats run_one(const Variant& V,
                        int N, double Lx, double Ly, double Lz,
                        double a, double b, double c,
                        bool tecplot_dump=false,
                        // optional outputs for interior-point comparison
                        std::vector<double>* out_Cx = nullptr,
                        std::vector<double>* out_Cy = nullptr,
                        std::vector<double>* out_Cz = nullptr,
                        std::vector<double>* out_CxA = nullptr,
                        std::vector<double>* out_CyA = nullptr,
                        std::vector<double>* out_CzA = nullptr)
{
  const int Nx=N, Ny=N, Nz=N;
  const size_t NTc = (size_t)Nx*Ny*Nz;     // centers
  const size_t NTk = (size_t)Nx*Ny*Nz;     // corners
  const double dx=Lx/Nx, dy=Ly/Ny, dz=Lz/Nz;

  // cell-centered B fields
  std::vector<double> Bx(NTc), By(NTc), Bz(NTc);

  for (int k=0; k<Nz; ++k)
    for (int j=0; j<Ny; ++j)
      for (int i=0; i<Nx; ++i) {
        const size_t idc = (size_t)k*Ny*Nx + (size_t)j*Nx + (size_t)i;
        const double x=(i+0.5)*dx, y=(j+0.5)*dy, z=(k+0.5)*dz;
        const V3 Bval = analyticB(x,y,z, a,b,c);
        Bx[idc]=Bval.x; By[idc]=Bval.y; Bz[idc]=Bval.z;
      }

  // build stencils for this variant
  EMSt::cCurlBStencil Rows[3];
  V.init(Rows, dx,dy,dz);

  // export integerized taps once per row/component
  cStencil::cStencilData SxBx, SxBy, SxBz;
  cStencil::cStencilData SyBx, SyBy, SyBz;
  cStencil::cStencilData SzBx, SzBy, SzBz;

  Rows[0].Bx.ExportStencil(&SxBx);
  Rows[0].By.ExportStencil(&SxBy);
  Rows[0].Bz.ExportStencil(&SxBz);

  Rows[1].Bx.ExportStencil(&SyBx);
  Rows[1].By.ExportStencil(&SyBy);
  Rows[1].Bz.ExportStencil(&SyBz);

  Rows[2].Bx.ExportStencil(&SzBx);
  Rows[2].By.ExportStencil(&SzBy);
  Rows[2].Bz.ExportStencil(&SzBz);

  // numeric curl(B) at corners (i,j,k) -> position (i*dx, j*dy, k*dz)
  std::vector<double> Cx(NTk), Cy(NTk), Cz(NTk);
  for (int k=0; k<Nz; ++k)
    for (int j=0; j<Ny; ++j)
      for (int i=0; i<Nx; ++i) {
        const size_t idk = (size_t)k*Ny*Nx + (size_t)j*Nx + (size_t)i;
        // row 0 (curl B)_x
        double vx = 0.0;
        vx += apply_exported(SxBx, Bx, i,j,k, Nx,Ny,Nz);
        vx += apply_exported(SxBy, By, i,j,k, Nx,Ny,Nz);
        vx += apply_exported(SxBz, Bz, i,j,k, Nx,Ny,Nz);
        Cx[idk] = vx;

        // row 1 (curl B)_y
        double vy = 0.0;
        vy += apply_exported(SyBx, Bx, i,j,k, Nx,Ny,Nz);
        vy += apply_exported(SyBy, By, i,j,k, Nx,Ny,Nz);
        vy += apply_exported(SyBz, Bz, i,j,k, Nx,Ny,Nz);
        Cy[idk] = vy;

        // row 2 (curl B)_z
        double vz = 0.0;
        vz += apply_exported(SzBx, Bx, i,j,k, Nx,Ny,Nz);
        vz += apply_exported(SzBy, By, i,j,k, Nx,Ny,Nz);
        vz += apply_exported(SzBz, Bz, i,j,k, Nx,Ny,Nz);
        Cz[idk] = vz;
      }

  // analytic curl(B) at the corners (x=i*dx, y=j*dy, z=k*dz)
  std::vector<double> CxA(NTk), CyA(NTk), CzA(NTk);
  for (int k=0; k<Nz; ++k)
    for (int j=0; j<Ny; ++j)
      for (int i=0; i<Nx; ++i) {
        const size_t idk = (size_t)k*Ny*Nx + (size_t)j*Nx + (size_t)i;
        const double x=i*dx, y=j*dy, z=k*dz;
        const V3 cu = analyticCurlB(x,y,z, a,b,c);
        CxA[idk]=cu.x; CyA[idk]=cu.y; CzA[idk]=cu.z;
      }

  // norms
  double linf=0.0, l2num=0.0, l2den=0.0;
  for (size_t t=0; t<NTk; ++t) {
    const double ex = Cx[t] - CxA[t];
    const double ey = Cy[t] - CyA[t];
    const double ez = Cz[t] - CzA[t];
    linf = std::max(linf, std::max(std::abs(ex), std::max(std::abs(ey), std::abs(ez))));
    l2num += ex*ex + ey*ey + ez*ez;
    l2den += CxA[t]*CxA[t] + CyA[t]*CyA[t] + CzA[t]*CzA[t];
  }

  // optional array copies for pointwise print
  if (out_Cx)  { *out_Cx  = Cx; }
  if (out_Cy)  { *out_Cy  = Cy; }
  if (out_Cz)  { *out_Cz  = Cz; }
  if (out_CxA) { *out_CxA = CxA; }
  if (out_CyA) { *out_CyA = CyA; }
  if (out_CzA) { *out_CzA = CzA; }

  // optional Tecplot on the finest grid
  if (tecplot_dump) {
    const std::string fn = std::string("curlB_")
                         + V.name + "_N" + std::to_string(N) + ".dat";
    std::ofstream out(fn.c_str());
    if (out) {
      out << "TITLE = \"curl(B) vs analytic\"\n";
      out << "VARIABLES = \"x\" \"y\" \"z\" "
             "\"Bx\" \"By\" \"Bz\" "
             "\"Cx_num\" \"Cy_num\" \"Cz_num\" "
             "\"Cx_ana\" \"Cy_ana\" \"Cz_ana\" "
             "\"Err_mag\"\n";
      out << "ZONE T=\"vol\", I="<<Nx<<", J="<<Ny<<", K="<<Nz<<", DATAPACKING=BLOCK\n";
      auto block = [&](auto f){
        out.setf(std::ios::scientific); out<<std::setprecision(8);
        for (int kk=0; kk<Nz; ++kk)
          for (int jj=0; jj<Ny; ++jj)
            for (int ii=0; ii<Nx; ++ii) out<<f(ii,jj,kk)<<"\n";
      };
      block([&](int i,int,  int){ return (i+0.0)*dx; });
      block([&](int,  int j,int){ return (j+0.0)*dy; });
      block([&](int,  int,  int k){ return (k+0.0)*dz; });

      auto dump = [&](const std::vector<double>& A){ for (size_t t=0;t<NTk;++t) out<<A[t]<<"\n"; };
      // keep B blank in this Tecplot for brevity
      std::vector<double> zeros(NTk, 0.0);
      dump(zeros); dump(zeros); dump(zeros);

      dump(Cx); dump(Cy); dump(Cz);
      dump(CxA); dump(CyA); dump(CzA);

      for (size_t t=0; t<NTk; ++t) {
        const double ex = Cx[t]-CxA[t], ey = Cy[t]-CyA[t], ez = Cz[t]-CzA[t];
        out << std::sqrt(ex*ex + ey*ey + ez*ez) << "\n";
      }
    }
  }

  ErrStats s;
  s.linf  = linf;
  s.l2rel = (l2den>0.0) ? std::sqrt(l2num/l2den) : std::sqrt(l2num);
  return s;
}

static int Run(const std::vector<std::string>& args) {
  // domain
  double Lx=1.0, Ly=1.0, Lz=1.0;

  // Use unequal wavenumbers so curl(B) ≠ 0
  double a=2.0*M_PI, b=3.0*M_PI, c=5.0*M_PI;

  std::vector<int> Ns = {16, 24, 32, 48, 64};
  bool tecplot_on_finest = true;

  for (size_t i=0; i<args.size(); ++i) {
    if (args[i]=="--no-tecplot") tecplot_on_finest=false;
  }

  struct Row { double linf, l2; };
  const int F = (int)(sizeof(kVariants)/sizeof(kVariants[0]));
  std::vector<Row> err[F];
  for (int f=0; f<F; ++f) err[f].resize(Ns.size());

  std::cout << "\n=== Corner curl(B) Convergence (edge/face 2nd, 4th, 6th) ===\n"
            << "Domain: Lx="<<Lx<<", Ly="<<Ly<<", Lz="<<Lz
            << "  wave numbers: a="<<a<<", b="<<b<<", c="<<c<<"\n";

  for (size_t t=0; t<Ns.size(); ++t) {
    const int N = Ns[t];
    const bool dump_tp = tecplot_on_finest && (t == Ns.size()-1);
    std::cout << "\nN = " << N << "\n";
    for (int f=0; f<F; ++f) {
      // capture arrays for the component-wise interior-point comparison
      std::vector<double> Cx, Cy, Cz, CxA, CyA, CzA;

      const ErrStats s = run_one(kVariants[f], N, Lx,Ly,Lz, a,b,c, dump_tp,
                                 &Cx, &Cy, &Cz, &CxA, &CyA, &CzA);
      err[f][t] = Row{s.linf, s.l2rel};

      std::cout << "  " << std::left << std::setw(22) << kVariants[f].name
                << " Linf=" << std::scientific << std::setprecision(3) << s.linf
                << "   L2=" << s.l2rel << std::fixed << "\n";

      // Detailed interior-point comparison (like GradDivE)
      print_point_comparison(N, Lx, Ly, Lz, Cx, Cy, Cz, CxA, CyA, CzA, kVariants[f].name);
    }
  }

  // combined convergence table
  std::cout << "\n------------------------------------------------------------------------------------------------------------------------------------------------------\n";
  std::cout << "   N  | ";
  for (int f=0; f<F; ++f)
    std::cout << std::setw(22) << kVariants[f].name << " (L_inf)  Ord   "
              << std::setw(22) << kVariants[f].name << " (L2)     Ord"
              << (f==F-1 ? "" : "  |  ");
  std::cout << "\n------------------------------------------------------------------------------------------------------------------------------------------------------\n";

  for (size_t t=0; t<Ns.size(); ++t) {
    std::cout << std::setw(5) << Ns[t] << " ";
    for (int f=0; f<F; ++f) {
      double p_inf=0.0, p_l2=0.0;
      if (t>0) {
        const double rN = double(Ns[t]) / double(Ns[t-1]);         // refinement ratio
        const double den = std::log(rN);                            // log(h_{t-1}/h_t) = log(N_t/N_{t-1})
        const double num_inf = err[f][t-1].linf / err[f][t].linf;
        const double num_l2  = err[f][t-1].l2   / err[f][t].l2;
        if (den > 0 && num_inf > 0) p_inf = std::log(num_inf) / den;
        if (den > 0 && num_l2  > 0) p_l2  = std::log(num_l2 ) / den;     
      }
      std::cout << " | "
        << std::scientific << std::setprecision(3) << std::setw(12) << err[f][t].linf << " "
        << std::fixed      << std::setprecision(2) << std::setw(5)  << (t? p_inf:0.0) << "  "
        << std::scientific << std::setprecision(3) << std::setw(12) << err[f][t].l2   << " "
        << std::fixed      << std::setprecision(2) << std::setw(5)  << (t? p_l2 :0.0);
    }
    std::cout << "\n";
  }
  std::cout << "------------------------------------------------------------------------------------------------------------------------------------------------------\n";

  if (tecplot_on_finest) {
    std::cout << "\nTecplot files written on finest grid:\n";
    for (int f=0; f<F; ++f)
      std::cout << "  curlB_" << kVariants[f].name << "_N" << Ns.back() << ".dat\n";
  }
  return 0;
}

} // namespace CurlB

// Auto-register in the harness
static TestHarness::AutoRegister _auto_reg_curlb(
  "curl_b",
  [](const std::vector<std::string>& args){ return CurlB::Run(args); },
  "Corner curl(B) stencils: build, apply, and convergence (2nd, 4th, 6th order) + component-wise comparison.");

// ---- Per-suite force-link shim so this TU is retained even with dead-stripping ----
namespace CurlB {
  void ForceLinkAllTests() {}
}

