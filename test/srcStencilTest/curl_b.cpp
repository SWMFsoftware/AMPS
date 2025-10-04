/***************************************************************************************
 * curl_b.cpp — Convergence test for corner curl(B) stencils (2nd & 4th order)
 * -------------------------------------------------------------------------------------
 * ENHANCED VERSION: Added numerical curl calculation using face-center averaging
 ***************************************************************************************/

#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <map>

#include "test_harness.h"

static constexpr bool kPrintDebug = false; // master switch for verbose/debug output


namespace StencilTests {
namespace CurlB {

// Shorthand to the production stencil namespace for readability in this TU.
namespace EMSt = PIC::FieldSolver::Electromagnetic::ECSIM::Stencil;

// ----------------------------
// Small uniform grid container
// ----------------------------
struct Grid {
  int Nx, Ny, Nz;       // number of cells in each direction
  double dx, dy, dz;    // cell sizes (uniform)
  std::vector<double> Bx, By, Bz; // cell-centered fields

  explicit Grid(int nx,int ny,int nz):Nx(nx),Ny(ny),Nz(nz) {
    dx = 1.0 / Nx;  dy = 1.0 / Ny;  dz = 1.0 / Nz;
    Bx.assign(Nx*Ny*Nz, 0.0);
    By.assign(Nx*Ny*Nz, 0.0);
    Bz.assign(Nx*Ny*Nz, 0.0);
  }
  inline int idx(int i,int j,int k) const { return (k*Ny + j)*Nx + i; }
  
  // Safe accessor with bounds checking
  inline double getBx(int i, int j, int k) const {
    if (i<0 || i>=Nx || j<0 || j>=Ny || k<0 || k>=Nz) return 0.0;
    return Bx[idx(i,j,k)];
  }
  inline double getBy(int i, int j, int k) const {
    if (i<0 || i>=Nx || j<0 || j>=Ny || k<0 || k>=Nz) return 0.0;
    return By[idx(i,j,k)];
  }
  inline double getBz(int i, int j, int k) const {
    if (i<0 || i>=Nx || j<0 || j>=Ny || k<0 || k>=Nz) return 0.0;
    return Bz[idx(i,j,k)];
  }
};

// ----------------------------------------
// Analytic field and exact curl at corners
// ----------------------------------------

// Fill Bx,By,Bz at cell centers with:  Bx=y^3, By=z^3, Bz=x^3
static void FillCubic(Grid& g) {
  for (int k=0;k<g.Nz;k++) {
    const double zc = (k + 0.5)*g.dz;
    for (int j=0;j<g.Ny;j++) {
      const double yc = (j + 0.5)*g.dy;
      for (int i=0;i<g.Nx;i++) {
        const double xc = (i + 0.5)*g.dx;
        g.Bx[g.idx(i,j,k)] = yc*yc*yc; // y^3
        g.By[g.idx(i,j,k)] = zc*zc*zc; // z^3
        g.Bz[g.idx(i,j,k)] = xc*xc*xc; // x^3
      }
    }
  }
}

// Exact curl(B) at a corner location (x=i*dx, y=j*dy, z=k*dz)
static inline std::array<double,3> CurlExactCorner(double x,double y,double z) {
  // curl = (-3 z^2, -3 x^2, -3 y^2)
  return { -3.0*z*z, -3.0*x*x, -3.0*y*y };
}

// -----------------------------------------------------------------------
// NEW: Numerical curl calculation at cell corners from cell-centered B
// -----------------------------------------------------------------------
// Corner at (i,j,k) is at physical location (i*dx, j*dy, k*dz).
// Cell centers are at ((i+0.5)*dx, (j+0.5)*dy, (k+0.5)*dz).
//
// To compute curl at corner (i,j,k), we need derivatives at that corner.
// We'll use a 2nd-order accurate centered difference approach:
//
// For (curl B)_x = dBz/dy - dBy/dz at corner (i,j,k):
//   - dBz/dy: Use cells straddling the corner in y-direction
//   - dBy/dz: Use cells straddling the corner in z-direction
//
// Strategy: Average the 4 surrounding cell-centered values to get face values,
// then take differences between faces.
//
// Corner (i,j,k) is surrounded by 8 cells with indices:
//   (i-1,j-1,k-1), (i,j-1,k-1), (i-1,j,k-1), (i,j,k-1),
//   (i-1,j-1,k),   (i,j-1,k),   (i-1,j,k),   (i,j,k)
//

// Helper to print the stencil used for face-averaging method
static void PrintFaceAvgStencil(const Grid& g, int i, int j, int k, double dx, double dy, double dz) {
  if (!kPrintDebug) { return; }

  std::printf("\n========== FACE-AVERAGING STENCIL at corner (i,j,k)=(%d,%d,%d) ==========\n", i, j, k);
  std::printf("Corner location: (x,y,z) = (%.6g, %.6g, %.6g)\n", i*dx, j*dy, k*dz);
  std::printf("\nThis stencil uses 8 surrounding cell centers to compute curl via face averaging.\n");
  std::printf("Cell centers are at (i+0.5, j+0.5, k+0.5) in index space.\n\n");
  
  // List the 8 surrounding cells
  std::printf("--- 8 Surrounding Cell Centers ---\n");
  int cells[8][3] = {
    {i-1, j-1, k-1}, {i, j-1, k-1}, {i-1, j, k-1}, {i, j, k-1},
    {i-1, j-1, k},   {i, j-1, k},   {i-1, j, k},   {i, j, k}
  };
  for (int n = 0; n < 8; ++n) {
    int ci = cells[n][0], cj = cells[n][1], ck = cells[n][2];
    double px = (ci+0.5)*dx, py = (cj+0.5)*dy, pz = (ck+0.5)*dz;
    std::printf("  Cell %d: (%d,%d,%d) at (%.6g,%.6g,%.6g) | Bx=%.6g, By=%.6g, Bz=%.6g\n", 
                n, ci, cj, ck, px, py, pz,
                g.getBx(ci,cj,ck), g.getBy(ci,cj,ck), g.getBz(ci,cj,ck));
  }
  
  // Calculate actual numerical values for derivatives
  // dBz/dy
  double Bz_yplus  = 0.25 * (g.getBz(i-1,j,k-1) + g.getBz(i,j,k-1) + 
                              g.getBz(i-1,j,k)   + g.getBz(i,j,k));
  double Bz_yminus = 0.25 * (g.getBz(i-1,j-1,k-1) + g.getBz(i,j-1,k-1) + 
                              g.getBz(i-1,j-1,k)   + g.getBz(i,j-1,k));
  double dBz_dy = (Bz_yplus - Bz_yminus) / dy;
  
  // dBy/dz
  double By_zplus  = 0.25 * (g.getBy(i-1,j-1,k) + g.getBy(i,j-1,k) + 
                              g.getBy(i-1,j,k)   + g.getBy(i,j,k));
  double By_zminus = 0.25 * (g.getBy(i-1,j-1,k-1) + g.getBy(i,j-1,k-1) + 
                              g.getBy(i-1,j,k-1)   + g.getBy(i,j,k-1));
  double dBy_dz = (By_zplus - By_zminus) / dz;
  
  // dBx/dz
  double Bx_zplus  = 0.25 * (g.getBx(i-1,j-1,k) + g.getBx(i,j-1,k) + 
                              g.getBx(i-1,j,k)   + g.getBx(i,j,k));
  double Bx_zminus = 0.25 * (g.getBx(i-1,j-1,k-1) + g.getBx(i,j-1,k-1) + 
                              g.getBx(i-1,j,k-1)   + g.getBx(i,j,k-1));
  double dBx_dz = (Bx_zplus - Bx_zminus) / dz;
  
  // dBz/dx
  double Bz_xplus  = 0.25 * (g.getBz(i,j-1,k-1) + g.getBz(i,j,k-1) + 
                              g.getBz(i,j-1,k)   + g.getBz(i,j,k));
  double Bz_xminus = 0.25 * (g.getBz(i-1,j-1,k-1) + g.getBz(i-1,j,k-1) + 
                              g.getBz(i-1,j-1,k)   + g.getBz(i-1,j,k));
  double dBz_dx = (Bz_xplus - Bz_xminus) / dx;
  
  // dBy/dx
  double By_xplus  = 0.25 * (g.getBy(i,j-1,k-1) + g.getBy(i,j,k-1) + 
                              g.getBy(i,j-1,k)   + g.getBy(i,j,k));
  double By_xminus = 0.25 * (g.getBy(i-1,j-1,k-1) + g.getBy(i-1,j,k-1) + 
                              g.getBy(i-1,j-1,k)   + g.getBy(i-1,j,k));
  double dBy_dx = (By_xplus - By_xminus) / dx;
  
  // dBx/dy
  double Bx_yplus  = 0.25 * (g.getBx(i-1,j,k-1) + g.getBx(i,j,k-1) + 
                              g.getBx(i-1,j,k)   + g.getBx(i,j,k));
  double Bx_yminus = 0.25 * (g.getBx(i-1,j-1,k-1) + g.getBx(i,j-1,k-1) + 
                              g.getBx(i-1,j-1,k)   + g.getBx(i,j-1,k));
  double dBx_dy = (Bx_yplus - Bx_yminus) / dy;
  
  // Show curl_x = dBz/dy - dBy/dz
  std::printf("\n--- Component: (curl B)_x = dBz/dy - dBy/dz ---\n");
  std::printf("  dBz/dy:\n");
  std::printf("    Bz at y+ face (y=%.6g): avg of cells (%d,%d,%d), (%d,%d,%d), (%d,%d,%d), (%d,%d,%d)\n",
              j*dy + dy/2, i-1,j,k-1, i,j,k-1, i-1,j,k, i,j,k);
  std::printf("    Bz at y- face (y=%.6g): avg of cells (%d,%d,%d), (%d,%d,%d), (%d,%d,%d), (%d,%d,%d)\n",
              j*dy - dy/2, i-1,j-1,k-1, i,j-1,k-1, i-1,j-1,k, i,j-1,k);
  std::printf("    => dBz/dy = (%.10g - %.10g) / %.6g = %.10g\n", Bz_yplus, Bz_yminus, dy, dBz_dy);
  std::printf("    Contribution to curl_x:\n");
  std::printf("      Bz[%d,%d,%d]: +%.6g  (from y+ avg)\n", i-1,j,k-1, 0.25/dy);
  std::printf("      Bz[%d,%d,%d]: +%.6g\n", i,j,k-1, 0.25/dy);
  std::printf("      Bz[%d,%d,%d]: +%.6g\n", i-1,j,k, 0.25/dy);
  std::printf("      Bz[%d,%d,%d]: +%.6g\n", i,j,k, 0.25/dy);
  std::printf("      Bz[%d,%d,%d]: %.6g  (from y- avg)\n", i-1,j-1,k-1, -0.25/dy);
  std::printf("      Bz[%d,%d,%d]: %.6g\n", i,j-1,k-1, -0.25/dy);
  std::printf("      Bz[%d,%d,%d]: %.6g\n", i-1,j-1,k, -0.25/dy);
  std::printf("      Bz[%d,%d,%d]: %.6g\n", i,j-1,k, -0.25/dy);
  
  std::printf("  dBy/dz:\n");
  std::printf("    By at z+ face (z=%.6g): avg of cells (%d,%d,%d), (%d,%d,%d), (%d,%d,%d), (%d,%d,%d)\n",
              k*dz + dz/2, i-1,j-1,k, i,j-1,k, i-1,j,k, i,j,k);
  std::printf("    By at z- face (z=%.6g): avg of cells (%d,%d,%d), (%d,%d,%d), (%d,%d,%d), (%d,%d,%d)\n",
              k*dz - dz/2, i-1,j-1,k-1, i,j-1,k-1, i-1,j,k-1, i,j,k-1);
  std::printf("    => dBy/dz = (%.10g - %.10g) / %.6g = %.10g\n", By_zplus, By_zminus, dz, dBy_dz);
  std::printf("    Contribution to curl_x:\n");
  std::printf("      By[%d,%d,%d]: %.6g  (from z+ avg, negative sign)\n", i-1,j-1,k, -0.25/dz);
  std::printf("      By[%d,%d,%d]: %.6g\n", i,j-1,k, -0.25/dz);
  std::printf("      By[%d,%d,%d]: %.6g\n", i-1,j,k, -0.25/dz);
  std::printf("      By[%d,%d,%d]: %.6g\n", i,j,k, -0.25/dz);
  std::printf("      By[%d,%d,%d]: +%.6g  (from z- avg, negative sign)\n", i-1,j-1,k-1, 0.25/dz);
  std::printf("      By[%d,%d,%d]: +%.6g\n", i,j-1,k-1, 0.25/dz);
  std::printf("      By[%d,%d,%d]: +%.6g\n", i-1,j,k-1, 0.25/dz);
  std::printf("      By[%d,%d,%d]: +%.6g\n", i,j,k-1, 0.25/dz);
  
  // Show curl_y = dBx/dz - dBz/dx
  std::printf("\n--- Component: (curl B)_y = dBx/dz - dBz/dx ---\n");
  std::printf("  dBx/dz:\n");
  std::printf("    Bx at z+ face (z=%.6g): avg of cells (%d,%d,%d), (%d,%d,%d), (%d,%d,%d), (%d,%d,%d)\n",
              k*dz + dz/2, i-1,j-1,k, i,j-1,k, i-1,j,k, i,j,k);
  std::printf("    Bx at z- face (z=%.6g): avg of cells (%d,%d,%d), (%d,%d,%d), (%d,%d,%d), (%d,%d,%d)\n",
              k*dz - dz/2, i-1,j-1,k-1, i,j-1,k-1, i-1,j,k-1, i,j,k-1);
  std::printf("    => dBx/dz = (%.10g - %.10g) / %.6g = %.10g\n", Bx_zplus, Bx_zminus, dz, dBx_dz);
  std::printf("    => Contribution: coeff = +1/dz = %.6g\n", 1.0/dz);
  
  std::printf("  dBz/dx:\n");
  std::printf("    Bz at x+ face (x=%.6g): avg of cells (%d,%d,%d), (%d,%d,%d), (%d,%d,%d), (%d,%d,%d)\n",
              i*dx + dx/2, i,j-1,k-1, i,j,k-1, i,j-1,k, i,j,k);
  std::printf("    Bz at x- face (x=%.6g): avg of cells (%d,%d,%d), (%d,%d,%d), (%d,%d,%d), (%d,%d,%d)\n",
              i*dx - dx/2, i-1,j-1,k-1, i-1,j,k-1, i-1,j-1,k, i-1,j,k);
  std::printf("    => dBz/dx = (%.10g - %.10g) / %.6g = %.10g\n", Bz_xplus, Bz_xminus, dx, dBz_dx);
  std::printf("    => Contribution: coeff = -1/dx = %.6g\n", -1.0/dx);
  
  // Show curl_z = dBy/dx - dBx/dy
  std::printf("\n--- Component: (curl B)_z = dBy/dx - dBx/dy ---\n");
  std::printf("  dBy/dx:\n");
  std::printf("    By at x+ face (x=%.6g): avg of cells (%d,%d,%d), (%d,%d,%d), (%d,%d,%d), (%d,%d,%d)\n",
              i*dx + dx/2, i,j-1,k-1, i,j,k-1, i,j-1,k, i,j,k);
  std::printf("    By at x- face (x=%.6g): avg of cells (%d,%d,%d), (%d,%d,%d), (%d,%d,%d), (%d,%d,%d)\n",
              i*dx - dx/2, i-1,j-1,k-1, i-1,j,k-1, i-1,j-1,k, i-1,j,k);
  std::printf("    => dBy/dx = (%.10g - %.10g) / %.6g = %.10g\n", By_xplus, By_xminus, dx, dBy_dx);
  std::printf("    => Contribution: coeff = +1/dx = %.6g\n", 1.0/dx);
  
  std::printf("  dBx/dy:\n");
  std::printf("    Bx at y+ face (y=%.6g): avg of cells (%d,%d,%d), (%d,%d,%d), (%d,%d,%d), (%d,%d,%d)\n",
              j*dy + dy/2, i-1,j,k-1, i,j,k-1, i-1,j,k, i,j,k);
  std::printf("    Bx at y- face (y=%.6g): avg of cells (%d,%d,%d), (%d,%d,%d), (%d,%d,%d), (%d,%d,%d)\n",
              j*dy - dy/2, i-1,j-1,k-1, i,j-1,k-1, i-1,j-1,k, i,j-1,k);
  std::printf("    => dBx/dy = (%.10g - %.10g) / %.6g = %.10g\n", Bx_yplus, Bx_yminus, dy, dBx_dy);
  std::printf("    => Contribution: coeff = -1/dy = %.6g\n", -1.0/dy);
  
  std::printf("====================================\n\n");
}

static std::array<double,3> CurlNumericalFaceAvg(const Grid& g, int i, int j, int k, bool print_stencil = false) {
  const double dx = g.dx, dy = g.dy, dz = g.dz;
  
  if (print_stencil) {
    PrintFaceAvgStencil(g, i, j, k, dx, dy, dz);
  }
  
  // Corner is at (i*dx, j*dy, k*dz)
  // The 8 surrounding cell centers are at (i±0.5, j±0.5, k±0.5) in index space
  
  // ===== (curl B)_x = dBz/dy - dBy/dz =====
  
  // dBz/dy: Need Bz at y-faces (at y = j*dy ± dy/2)
  // y+ face (j+0.5 in index space, j*dy+dy/2 in physical): average 4 cells in xz plane at j,j+1
  // y- face (j-0.5 in index space, j*dy-dy/2 in physical): average 4 cells in xz plane at j-1,j
  double Bz_yplus  = 0.25 * (g.getBz(i-1,j,k-1) + g.getBz(i,j,k-1) + 
                              g.getBz(i-1,j,k)   + g.getBz(i,j,k));
  double Bz_yminus = 0.25 * (g.getBz(i-1,j-1,k-1) + g.getBz(i,j-1,k-1) + 
                              g.getBz(i-1,j-1,k)   + g.getBz(i,j-1,k));
  double dBz_dy = (Bz_yplus - Bz_yminus) / dy;
  
  // dBy/dz: Need By at z-faces (at z = k*dz ± dz/2)
  // z+ face (k+0.5 in index space): average 4 cells in xy plane at k,k+1
  // z- face (k-0.5 in index space): average 4 cells in xy plane at k-1,k
  double By_zplus  = 0.25 * (g.getBy(i-1,j-1,k) + g.getBy(i,j-1,k) + 
                              g.getBy(i-1,j,k)   + g.getBy(i,j,k));
  double By_zminus = 0.25 * (g.getBy(i-1,j-1,k-1) + g.getBy(i,j-1,k-1) + 
                              g.getBy(i-1,j,k-1)   + g.getBy(i,j,k-1));
  double dBy_dz = (By_zplus - By_zminus) / dz;
  
  double curl_x = dBz_dy - dBy_dz;
  
  // ===== (curl B)_y = dBx/dz - dBz/dx =====
  
  // dBx/dz: Need Bx at z-faces
  double Bx_zplus  = 0.25 * (g.getBx(i-1,j-1,k) + g.getBx(i,j-1,k) + 
                              g.getBx(i-1,j,k)   + g.getBx(i,j,k));
  double Bx_zminus = 0.25 * (g.getBx(i-1,j-1,k-1) + g.getBx(i,j-1,k-1) + 
                              g.getBx(i-1,j,k-1)   + g.getBx(i,j,k-1));
  double dBx_dz = (Bx_zplus - Bx_zminus) / dz;
  
  // dBz/dx: Need Bz at x-faces (at x = i*dx ± dx/2)
  // x+ face (i+0.5 in index space): average 4 cells in yz plane at i,i+1
  // x- face (i-0.5 in index space): average 4 cells in yz plane at i-1,i
  double Bz_xplus  = 0.25 * (g.getBz(i,j-1,k-1) + g.getBz(i,j,k-1) + 
                              g.getBz(i,j-1,k)   + g.getBz(i,j,k));
  double Bz_xminus = 0.25 * (g.getBz(i-1,j-1,k-1) + g.getBz(i-1,j,k-1) + 
                              g.getBz(i-1,j-1,k)   + g.getBz(i-1,j,k));
  double dBz_dx = (Bz_xplus - Bz_xminus) / dx;
  
  double curl_y = dBx_dz - dBz_dx;
  
  // ===== (curl B)_z = dBy/dx - dBx/dy =====
  
  // dBy/dx: Need By at x-faces
  double By_xplus  = 0.25 * (g.getBy(i,j-1,k-1) + g.getBy(i,j,k-1) + 
                              g.getBy(i,j-1,k)   + g.getBy(i,j,k));
  double By_xminus = 0.25 * (g.getBy(i-1,j-1,k-1) + g.getBy(i-1,j,k-1) + 
                              g.getBy(i-1,j-1,k)   + g.getBy(i-1,j,k));
  double dBy_dx = (By_xplus - By_xminus) / dx;
  
  // dBx/dy: Need Bx at y-faces
  double Bx_yplus  = 0.25 * (g.getBx(i-1,j,k-1) + g.getBx(i,j,k-1) + 
                              g.getBx(i-1,j,k)   + g.getBx(i,j,k));
  double Bx_yminus = 0.25 * (g.getBx(i-1,j-1,k-1) + g.getBx(i,j-1,k-1) + 
                              g.getBx(i-1,j-1,k)   + g.getBx(i,j-1,k));
  double dBx_dy = (Bx_yplus - Bx_yminus) / dy;
  
  double curl_z = dBy_dx - dBx_dy;
  
  return {curl_x, curl_y, curl_z};
}

// --------------------------------------------------
// Utility: parse N-list override (e.g., "N=16,24,32")
// --------------------------------------------------
static std::vector<int> ParseNList(const std::vector<std::string>& args) {
  std::vector<int> Ns = {16, 24, 32, 48, 64}; // default refinement ladder
  for (const auto& a : args) {
    if (a.rfind("N=",0) == 0) {
      Ns.clear();
      std::stringstream ss(a.substr(2));
      while (ss.good()) {
        int v; char sep;
        if (ss >> v) Ns.push_back(v);
        if (!(ss >> sep)) break;
      }
    }
  }
  return Ns;
}

// -----------------------------------------------------------------
// Export helper: get the compact tables for each contributing field
// -----------------------------------------------------------------
using Compact = ::cStencil::cStencilData;

static inline void ExportCurlComponentTables(EMSt::cCurlBStencil& C, Compact out[3],double scale) {
  C.Bx.ExportStencil(&out[0],scale);
  C.By.ExportStencil(&out[1],scale);
  C.Bz.ExportStencil(&out[2],scale);
}

// --------------------------------------------------------
// NEW: Print stencil in readable format
// --------------------------------------------------------
static void PrintStencilTables(const Compact S[3][3], const char* variantName) {
  const char* comp_names[3] = {"(curl B)_x", "(curl B)_y", "(curl B)_z"};
  const char* field_names[3] = {"Bx", "By", "Bz"};
  
  std::printf("\n========== STENCIL TABLES: %s ==========\n", variantName);
  
  for (int c = 0; c < 3; ++c) {
    std::printf("\n--- Component: %s ---\n", comp_names[c]);
    for (int f = 0; f < 3; ++f) {
      std::printf("  Contribution from %s:\n", field_names[f]);
      const Compact& St = S[c][f];
      if (St.Length == 0) {
        std::printf("    (empty)\n");
        continue;
      }
      for (int n = 0; n < St.Length; ++n) {
        std::printf("    [%2d] coeff=%+.8e  offset=(%+2d,%+2d,%+2d)\n",
                    n, St.Data[n].a, St.Data[n].i, St.Data[n].j, St.Data[n].k);
      }
    }
  }
  std::printf("====================================\n\n");
}

// --------------------------------------------------------
// Compute maximum reach (padding) from a 3×3 set of tables
// --------------------------------------------------------
static int StencilMaxPad(const Compact S[3][3]) {
  int pad = 0;
  for (int comp=0; comp<3; ++comp) {
    for (int f=0; f<3; ++f) {
      for (int n=0; n<S[comp][f].Length; ++n) {
        const int di = std::abs(S[comp][f].Data[n].i);
        const int dj = std::abs(S[comp][f].Data[n].j);
        const int dk = std::abs(S[comp][f].Data[n].k);
        pad = std::max(pad, std::max(di, std::max(dj, dk)));
      }
    }
  }
  return pad;
}

// ------------------------------------------------------
// Apply one curl component's tables at base (i,j,k) cell
// ------------------------------------------------------
static inline double ApplyComponentAt(
    const Grid& g, int i,int j,int k,
    const Compact& SBx, const Compact& SBy, const Compact& SBz)
{
  auto accum = [&](const Compact& S, const std::vector<double>& F) {
    double s = 0.0;
    for (int n=0; n<S.Length; ++n) {
      const int ii = i + S.Data[n].i;
      const int jj = j + S.Data[n].j;
      const int kk = k + S.Data[n].k;
      s += S.Data[n].a * F[g.idx(ii,jj,kk)];
    }
    return s;
  };

  return accum(SBx, g.Bx) + accum(SBy, g.By) + accum(SBz, g.Bz);
}

// ------------------------------------------------------
// Apply component with detailed printing (verbose version)
// ------------------------------------------------------
static inline double ApplyComponentAtVerbose(
    const Grid& g, int i, int j, int k,
    const Compact& SBx, const Compact& SBy, const Compact& SBz,
    const char* comp_name, const char* formula)
{
  const bool _print = kPrintDebug;

  if(_print) std::printf("\n--- Component: %s = %s ---\n", comp_name, formula);
  
  auto accum_verbose = [&](const Compact& S, const std::vector<double>& F, const char* field_name) {
    double s = 0.0;
    if (S.Length == 0) {
      if(_print) std::printf("  No contribution from %s\n", field_name);
      return s;
    }
    
    if(_print) std::printf("  Contribution from %s (indices as i,j,k):\n", field_name);
    
    for (int n = 0; n < S.Length; ++n) {
      const int ii = i + S.Data[n].i;
      const int jj = j + S.Data[n].j;
      const int kk = k + S.Data[n].k;
      if(_print) std::printf("    [%d,%d,%d] coeff=%+.6g F=%e\n", ii, jj, kk, S.Data[n].a,F[g.idx(ii, jj, kk)]);
      s += S.Data[n].a * F[g.idx(ii, jj, kk)];
    }
    
    if(_print) std::printf("    Total from %s: %.10g\n", field_name, s);
    return s;
  };

  double result = 0.0;
  result += accum_verbose(SBx, g.Bx, "Bx");
  result += accum_verbose(SBy, g.By, "By");
  result += accum_verbose(SBz, g.Bz, "Bz");
  
  if(_print) std::printf("  ===> Final %s = %.10g\n", comp_name, result);
  
  return result;
}

// ------------------------
// Simple error accumulators
// ------------------------
struct Err { double l2=0, linf=0; long n=0; };
static inline void Accum(Err& e, const double num[3], const double ex[3]) {
  for (int c=0;c<3;c++) {
    const double d = num[c] - ex[c];
    e.l2 += d*d;
    e.linf = std::max(e.linf, std::abs(d));
    e.n++;
  }
}
static inline void Finish(Err& e) {
  e.l2 = std::sqrt(e.l2 / std::max<long>(1, e.n));
}

// ---------------------
// The test entry point
// ---------------------
static int Run(const std::vector<std::string>& args) {
  namespace EMSt = PIC::FieldSolver::Electromagnetic::ECSIM::Stencil;

  // Variant descriptor: name + builder function
  using InitFn = void (*)(EMSt::cCurlBStencil*, double,double,double);
  struct Variant { const char* name; InitFn init; };

  // Three requested variants
  const Variant variants[] = {
    { "2nd-order (edge-based)",
      PIC::FieldSolver::Electromagnetic::ECSIM::Stencil::SecondOrder::InitCurlBStencils_edge_based },
    { "2nd-order (face-based)",
      PIC::FieldSolver::Electromagnetic::ECSIM::Stencil::SecondOrder::InitCurlBStencils_face_based },
    { "4th-order",
      PIC::FieldSolver::Electromagnetic::ECSIM::Stencil::FourthOrder::InitCurlBStencils }
  };

  // Grid sizes to test
  const auto Ns = ParseNList(args);

  // Loop over stencil variants
  for (const auto& V : variants) {
    // (1) Build stencils with unit spacings
    EMSt::cCurlBStencil Curl[3];
    V.init(Curl, 1.0, 1.0, 1.0);

    // (2) Export compact tables (unit scale for display)
    Compact S_unit[3][3];
    for (int c=0;c<3;c++) {
      ExportCurlComponentTables(Curl[c], S_unit[c], 1.0);
    }

    // (3) Print stencil tables
    if (kPrintDebug) PrintStencilTables(S_unit, V.name);

    // (4) Determine interior pad
    const int pad = StencilMaxPad(S_unit);

    std::printf("\n=== curl_b variant: %s ===\n", V.name);
    std::printf("Interior padding required: %d\n\n", pad);
    std::printf("%6s %10s %16s %16s\n", "N", "h", "L2", "Linf");

    double prevL2 = -1.0;
    for (int N : Ns) {
      // Build grid and fill analytic cell-centered B
      Grid g(N, N, N);
      FillCubic(g);

      // Export and scale compact tables
      Compact S[3][3];
      for (int c=0;c<3;c++) {
        ExportCurlComponentTables(Curl[c], S[c], N);
      }

      Err e; // L2/Linf accumulator

      // Sweep interior corners
      for (int k = pad; k < g.Nz - pad; ++k) {
        const double z = k * g.dz;
        for (int j = pad; j < g.Ny - pad; ++j) {
          const double y = j * g.dy;
          for (int i = pad; i < g.Nx - pad; ++i) {
            const double x = i * g.dx;

            // Exact curl at this corner
            const auto ex = CurlExactCorner(x,y,z);

            // Numerical curl (stencil-based)
            double num[3];
            num[0] = ApplyComponentAt(g, i,j,k, S[0][0], S[0][1], S[0][2]);
            num[1] = ApplyComponentAt(g, i,j,k, S[1][0], S[1][1], S[1][2]);
            num[2] = ApplyComponentAt(g, i,j,k, S[2][0], S[2][1], S[2][2]);

            // Accumulate vector error
            Accum(e, num, ex.data());
          }
        }
      }

      // Finalize norms and print
      Finish(e);
      std::printf("%6d %10.4g %16.3e %16.3e\n", N, g.dx, e.l2, e.linf);

      // Observed order
      if (prevL2 > 0.0) {
        const double p = std::log(prevL2 / e.l2) / std::log(2.0);
        std::printf("      observed L2 order ≈ %.2f\n", p);
      }
      prevL2 = e.l2;

      // ----------------------------------------------------------------------
      // Component-wise comparison at representative interior corner
      // ----------------------------------------------------------------------
      const int ic = std::min(std::max(N/2, pad), N-1-pad);
      const int jc = std::min(std::max(N/2, pad), N-1-pad);
      const int kc = std::min(std::max(N/2, pad), N-1-pad);
      const double xc = ic * g.dx, yc = jc * g.dy, zc = kc * g.dz;

      const auto exC = CurlExactCorner(xc, yc, zc);
      
      // Compute stencil-based numerical curl (quiet)
      double numStencil[3];
      numStencil[0] = ApplyComponentAt(g, ic,jc,kc, S[0][0], S[0][1], S[0][2]);
      numStencil[1] = ApplyComponentAt(g, ic,jc,kc, S[1][0], S[1][1], S[1][2]);
      numStencil[2] = ApplyComponentAt(g, ic,jc,kc, S[2][0], S[2][1], S[2][2]);
      
      // For the first grid size, optionally print a verbose breakdown
      if (kPrintDebug && (N == Ns[0])) {
        std::printf("\n========== STENCIL APPLICATION: %s, N=%d ==========%s", V.name, N, "\n");
        std::printf("Corner: (i,j,k)=(%d,%d,%d), Physical: (x,y,z)=(%.6g,%.6g,%.6g)\n", ic, jc, kc, xc, yc, zc);
        (void)ApplyComponentAtVerbose(g, ic,jc,kc, S[0][0], S[0][1], S[0][2], "(curl B)_x", "dBz/dy - dBy/dz");
        (void)ApplyComponentAtVerbose(g, ic,jc,kc, S[1][0], S[1][1], S[1][2], "(curl B)_y", "dBx/dz - dBz/dx");
        (void)ApplyComponentAtVerbose(g, ic,jc,kc, S[2][0], S[2][1], S[2][2], "(curl B)_z", "dBy/dx - dBx/dy");
      }

      
      // Stencil-based numerical curl (verbose for first grid)
      if (kPrintDebug && (N == Ns[0])) {
        std::printf("\n========== STENCIL APPLICATION: %s, N=%d ==========\n", V.name, N);
        std::printf("Corner: (i,j,k)=(%d,%d,%d), Physical: (x,y,z)=(%.6g,%.6g,%.6g)\n", 
                    ic, jc, kc, xc, yc, zc);
        
        numStencil[0] = ApplyComponentAtVerbose(g, ic,jc,kc, S[0][0], S[0][1], S[0][2],
                                                 "(curl B)_x", "dBz/dy - dBy/dz");
        numStencil[1] = ApplyComponentAtVerbose(g, ic,jc,kc, S[1][0], S[1][1], S[1][2],
                                                 "(curl B)_y", "dBx/dz - dBz/dx");
        numStencil[2] = ApplyComponentAtVerbose(g, ic,jc,kc, S[2][0], S[2][1], S[2][2],
                                                 "(curl B)_z", "dBy/dx - dBx/dy");
        std::printf("====================================\n\n");
      } else {
        numStencil[0] = ApplyComponentAt(g, ic,jc,kc, S[0][0], S[0][1], S[0][2]);
        numStencil[1] = ApplyComponentAt(g, ic,jc,kc, S[1][0], S[1][1], S[1][2]);
        numStencil[2] = ApplyComponentAt(g, ic,jc,kc, S[2][0], S[2][1], S[2][2]);
      }
      
      // Face-averaged numerical curl (print stencil for first grid only)
      const auto numFaceAvg = CurlNumericalFaceAvg(g, ic, jc, kc, (N == Ns[0]));

      auto rel = [](double a, double b){ 
        const double d = std::max(std::abs(b), 1e-14); 
        return std::abs(a-b)/d; 
      };

      std::puts("\n    Component-wise comparison at representative interior corner:");
      std::puts("    Component-wise comparison at representative interior corner:");
      std::printf("      Grid: N=%d, h=%.6g\n", N, g.dx);
      std::printf("      Corner: (i,j,k)=(%d,%d,%d), (x,y,z)=(%.6g, %.6g, %.6g)\n", 
                  ic, jc, kc, xc, yc, zc);
      std::puts("      ----------------------------------------------------------------------------------------------");
      std::puts("      Component      Analytic          Stencil           RelErr      FaceAvg           RelErr");
      std::puts("      ----------------------------------------------------------------------------------------------");
      
      const char* nm[3] = {"(curl B)_x", "(curl B)_y", "(curl B)_z"};
      for (int c=0;c<3;c++) {
        std::printf("      %-12s % .8e  % .8e  %.3e  % .8e  %.3e\n",
                    nm[c], exC[c], numStencil[c], rel(numStencil[c], exC[c]),
                    numFaceAvg[c], rel(numFaceAvg[c], exC[c]));
      }
      std::puts("");
    } // N
  } // variants

  return 0;
}

} // namespace CurlB

// ---------------------
// Harness registration
// ---------------------
static TestHarness::AutoRegister _auto_reg_curlb(
  "curl_b",
  [](const std::vector<std::string>& args){ return CurlB::Run(args); },
  "Corner curl(B) stencils: build, apply, and convergence (2nd vs 4th order) + component-wise comparison.");

void ForceLinkAllTests() {}

} // namespace StencilTests
