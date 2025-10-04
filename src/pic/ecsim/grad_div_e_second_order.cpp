/**
 * @file    grad_div_e_second_order.cpp
 * @brief   Second-order ∇(∇·E) stencil builders with curl(B)-style signature.
 *
 * API
 * ---
 *   struct cGradDivEStencil { cStencil Ex, Ey, Ez; };
 *
 *   // Compact (nested centered first/second differences)
 *   void PIC::FieldSolver::Electromagnetic::ECSIM::Stencil::SecondOrder
 *        ::InitGradDivEBStencils_compact(cGradDivEStencil* S, double dx, double dy, double dz);
 *
 *   // Wide (symmetric face/edge construction + rotations)
 *   void PIC::FieldSolver::Electromagnetic::ECSIM::Stencil::SecondOrder
 *        ::InitGradDivEBStencils_wide   (cGradDivEStencil* S, double dx, double dy, double dz);
 *
 * Semantics
 * ---------
 *   S has length 3 and holds the rows of G = ∇(∇·E):
 *     S[0] (Gx): [ dxx, dxy, dxz ]
 *     S[1] (Gy): [ dxy, dyy, dyz ]
 *     S[2] (Gz): [ dxz, dyz, dzz ]
 *
 * Notes
 * -----
 * - Both are formally **second order** on uniform grids.
 * - Metric scaling is applied internally:
 *     dxx *= 1/dx^2, dxy *= 1/(dx*dy), …, dzz *= 1/dz^2
 */

#include "../pic.h"

namespace PIC {
namespace FieldSolver {
namespace Electromagnetic {
namespace ECSIM {
namespace Stencil {

// -----------------------------------------------------------------------------
// Public row container (curl-style)
// -----------------------------------------------------------------------------
//struct cGradDivEStencil {
//  cStencil Ex, Ey, Ez;  // columns acting on (Ex, Ey, Ez) for one output row
//};

namespace SecondOrder {

// -----------------------------------------------------------------------------
// Common helpers (used by both flavors)
// -----------------------------------------------------------------------------
static inline void ExportToRows(const cStencil& dxx, const cStencil& dyy, const cStencil& dzz,
                                const cStencil& dxy, const cStencil& dxz, const cStencil& dyz,
                                cGradDivEStencil* S) {
  S[0].Ex = dxx; S[0].Ey = dxy; S[0].Ez = dxz;
  S[1].Ex = dxy; S[1].Ey = dyy; S[1].Ez = dyz;
  S[2].Ex = dxz; S[2].Ey = dyz; S[2].Ez = dzz;

  // Tidy up possible duplicate entries
  S[0].Ex.Simplify(); S[0].Ey.Simplify(); S[0].Ez.Simplify();
  S[1].Ex.Simplify(); S[1].Ey.Simplify(); S[1].Ez.Simplify();
  S[2].Ex.Simplify(); S[2].Ey.Simplify(); S[2].Ez.Simplify();
}

static inline void ScaleMetric(cStencil& dxx, cStencil& dyy, cStencil& dzz,
                               cStencil& dxy, cStencil& dxz, cStencil& dyz,
                               double dx, double dy, double dz) {
  const double i_dx2  = 1.0/(dx*dx);
  const double i_dy2  = 1.0/(dy*dy);
  const double i_dz2  = 1.0/(dz*dz);
  const double i_dxdy = 1.0/(dx*dy);
  const double i_dxdz = 1.0/(dx*dz);
  const double i_dydz = 1.0/(dy*dz);

  dxx *= i_dx2;  dyy *= i_dy2;  dzz *= i_dz2;
  dxy *= i_dxdy; dxz *= i_dxdz; dyz *= i_dydz;
}

// -----------------------------------------------------------------------------
// Helpers dedicated to the WIDE builder
// -----------------------------------------------------------------------------
namespace Helper_Wide {

// Build 2×2×2 corner average around (0,0,0) shifted to octant (I,J,K)∈{0,1}³
inline void BuildCorner(cStencil& st, int I, int J, int K) {
  for (int di = -1; di <= 0; ++di)
    for (int dj = -1; dj <= 0; ++dj)
      for (int dk = -1; dk <= 0; ++dk)
        st.add(1.0/8.0, I + di, J + dj, K + dk);
}

// Build the 12 edge averages (each a 2×2 slab)
inline void BuildEdges(cStencil edges[12]) {
  struct R { int imin,imax, jmin,jmax, kmin,kmax; };
  static const R rng[12] = {
    { 0,0,-1,0,-1,0},{ 0,0, 0,1,-1,0},{ 0,0, 0,1, 0,1},{ 0,0,-1,0, 0,1},
    {-1,0, 0,0,-1,0},{ 0,1, 0,0,-1,0},{ 0,1, 0,0, 0,1},{-1,0, 0,0, 0,1},
    {-1,0,-1,0, 0,0},{ 0,1,-1,0, 0,0},{ 0,1, 0,1, 0,0},{-1,0, 0,1, 0,0}
  };
  for (int f = 0; f < 12; ++f) {
    edges[f] = cStencil();
    for (int i = rng[f].imin; i <= rng[f].imax; ++i)
      for (int j = rng[f].jmin; j <= rng[f].jmax; ++j)
        for (int k = rng[f].kmin; k <= rng[f].kmax; ++k)
          edges[f].add(1.0/4.0, i, j, k);
  }
}

} // namespace Helper_Wide


// ============================================================================
// Compact (nested centered first/second differences)
// ============================================================================
void InitGradDivEBStencils_compact(cGradDivEStencil* S, double dx, double dy, double dz) {
  // Centered first derivatives (unscaled)
  cStencil Dx, Dy, Dz;
  Dx.add(+0.5, +1, 0, 0); Dx.add(-0.5, -1, 0, 0);
  Dy.add(+0.5,  0,+1, 0); Dy.add(-0.5,  0,-1, 0);
  Dz.add(+0.5,  0, 0,+1); Dz.add(-0.5,  0, 0,-1);

  // Classic 3-pt second derivatives (unscaled)
  cStencil Dxx, Dyy, Dzz;
  Dxx.add(+1.0, +1, 0, 0); Dxx.add(+1.0, -1, 0, 0); Dxx.add(-2.0, 0, 0, 0);
  Dyy.add(+1.0,  0,+1, 0); Dyy.add(+1.0,  0,-1, 0); Dyy.add(-2.0, 0, 0, 0);
  Dzz.add(+1.0,  0, 0,+1); Dzz.add(+1.0,  0, 0,-1); Dzz.add(-2.0, 0, 0, 0);

  // Mixed derivatives via nested centered first differences (unscaled)
  // Dxy = ∂x(∂/∂y) : outer ∂/∂x across i±1 planes of the Dy stencil
  cStencil Dxy; Dxy.AddShifted(Dy, +1, 0, 0, +0.5);
               Dxy.AddShifted(Dy, -1, 0, 0, -0.5);

  // Dxz = ∂x(∂/∂z)
  cStencil Dxz; Dxz.AddShifted(Dz, +1, 0, 0, +0.5);
               Dxz.AddShifted(Dz, -1, 0, 0, -0.5);

  // Dyz = ∂y(∂/∂z)
  cStencil Dyz; Dyz.AddShifted(Dz, 0, +1, 0, +0.5);
               Dyz.AddShifted(Dz, 0, -1, 0, -0.5);

  // Apply metric scaling
  ScaleMetric(Dxx, Dyy, Dzz, Dxy, Dxz, Dyz, dx, dy, dz);

  // Export to rows
  ExportToRows(Dxx, Dyy, Dzz, Dxy, Dxz, Dyz, S);
}


// ============================================================================
// Wide (symmetric face/edge construction + rotations)
// ============================================================================
void InitGradDivEBStencils_wide(cGradDivEStencil* S, double dx, double dy, double dz) {
  using namespace Helper_Wide;

  // ---- Corner and edge building blocks
  cStencil Corner[2][2][2];
  for (int I=0; I<2; ++I)
    for (int J=0; J<2; ++J)
      for (int K=0; K<2; ++K)
        BuildCorner(Corner[I][J][K], I, J, K);

  cStencil Edge[12]; BuildEdges(Edge);

  // ---- Pure seconds: build along x then average and rotate
  cStencil Dxx_face0, Dxx, Dyy, Dzz;

  // central second difference in x using an x-edge average
  Dxx_face0.AddShifted(Edge[0], +1, 0, 0, +1.0);
  Dxx_face0.add(-2.0, 0, 0, 0);
  Dxx_face0.AddShifted(Edge[0], -1, 0, 0, +1.0);

  // average that face over the four subfaces around the point
  for (int j=0; j<2; ++j)
    for (int k=0; k<2; ++k)
      Dxx.AddShifted(Dxx_face0, 0, j, k, 0.25);

  Dyy = Dxx; Dyy.SwitchAxes(0,1);
  Dzz = Dxx; Dzz.SwitchAxes(0,2);

  // ---- Mixed seconds
  // dxy: difference of ∂/∂x taken on adjacent y-faces
  cStencil Dxy_face, Dxy;
  {
    // ∂/∂x @ y-face j=0
    cStencil dEx_y0 = Corner[1][0][0]; dEx_y0.AddShifted(Corner[0][0][0], 0,0,0, -1.0);
    // ∂/∂x @ y-face j=1
    cStencil dEx_y1 = Corner[1][1][0]; dEx_y1.AddShifted(Corner[0][1][0], 0,0,0, -1.0);

    // centered ∂/∂y of (∂/∂x): (j=+1) − (j=0)
    Dxy_face = dEx_y1; Dxy_face.AddShifted(dEx_y0, 0,0,0, -1.0);

    // average across the two z positions to center at the node
    for (int k=0; k<2; ++k) Dxy.AddShifted(Dxy_face, 0, 0, k, 0.5);
  }

  // dxz: difference of ∂/∂x taken on adjacent z-faces
  cStencil Dxz_face, Dxz;
  {
    // ∂/∂x @ z-face k=0
    cStencil dEx_z0 = Corner[1][0][0]; dEx_z0.AddShifted(Corner[0][0][0], 0,0,0, -1.0);
    // ∂/∂x @ z-face k=1
    cStencil dEx_z1 = Corner[1][0][1]; dEx_z1.AddShifted(Corner[0][0][1], 0,0,0, -1.0);

    // centered ∂/∂z of (∂/∂x): (k=+1) − (k=0)
    Dxz_face = dEx_z1; Dxz_face.AddShifted(dEx_z0, 0,0,0, -1.0);

    // average across the two y positions to center at the node
    for (int j=0; j<2; ++j) Dxz.AddShifted(Dxz_face, 0, j, 0, 0.5);
  }

  // dyz: obtain by rotation from dxy template (swap x↔z)
  cStencil Dyz = Dxy; Dyz.SwitchAxes(0,2);

  // ---- Metric scaling
  ScaleMetric(Dxx, Dyy, Dzz, Dxy, Dxz, Dyz, dx, dy, dz);

  // ---- Export to rows
  ExportToRows(Dxx, Dyy, Dzz, Dxy, Dxz, Dyz, S);
}

} // namespace SecondOrder

} // namespace Stencil
} // namespace ECSIM
} // namespace Electromagnetic
} // namespace FieldSolver
} // namespace PIC

