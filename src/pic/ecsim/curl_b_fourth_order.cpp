/**
 * @file    curlb_stencil_fourth_order.cpp
 * @brief   Fourth-order curl(B) at cell corners from cell-centered B.
 *
 * DETAILED DERIVATION (compact)
 * =============================
 * Target: compute (∇×B) at corner (i-1/2, j-1/2, k-1/2) with B at cell centers (i,j,k).
 * Components:
 *   (curl B)_x = ∂Bz/∂y - ∂By/∂z   ; (curl B)_y = ∂Bx/∂z - ∂Bz/∂x   ; (curl B)_z = ∂By/∂x - ∂Bx/∂y
 *
 * For a partial like ∂Bz/∂y|corner to be O(h^4), we:
 *  (A) Form **4th-order edge values** of Bz at y-indices {j+1, j, j-1, j-2} on the edge through the corner:
 *      Each edge point is at (i-1/2, y_index, k-1/2).
 *      To reach it from cell centers, apply **4th-order midpoint interpolation** independently in x and z:
 *        1D midpoint weights:  w^(1/2) = [-1, 9, 9, -1] / 16  at indices {m-2, m-1, m, m+1}.
 *      Tensor the two directions (x and z): 4×4 = 16 contributing cells per edge value.
 *  (B) Take a **4th-order centered difference** at y-midpoint (j-1/2):
 *        f'(j-1/2) ≈ [ f_{j+1} - 8 f_{j} + 8 f_{j-1} - f_{j-2} ] / (12 h_y)
 * Repeat cyclically for the other partials. This yields O(h^4) at the corner. Many offsets overlap and
 * will be merged by S.Simplify().
 *
 * Usage:
 *   PIC::FieldSolver::Electromagnetic::ECSIM::Stencil::cCurlBStencil Curl4[3];
 *   PIC::FieldSolver::ECSIM::Stencil::FourthOrder::InitCurlBStencils(Curl4, dx, dy, dz);
 *   // Curl4[comp].{Bx,By,Bz} are the contributions to (curl B)_comp.
 */
#include "../pic.h"

namespace PIC {
namespace FieldSolver {
namespace ECSIM {
namespace Stencil {
namespace FourthOrder {

//using cStencil = PIC::FieldSolver::Electromagnetic::ECSIM::Stencil::cStencil;
using CurlComp = PIC::FieldSolver::Electromagnetic::ECSIM::Stencil::cCurlBStencil;

static inline void midpoint_w(double w[4]) { w[0]=-1.0/16.0; w[1]=9.0/16.0; w[2]=9.0/16.0; w[3]=-1.0/16.0; }
static inline void dmid_w(double w[4])     { w[0]=+1.0/12.0;  w[1]=-8.0/12.0; w[2]=+8.0/12.0; w[3]=-1.0/12.0; }

static void addFourthOrderPartial(cStencil& S, int axisDer, double inv_h, const char* sym) {
  S = cStencil();
  const cFrac f(-1,2); S.SetSymbol(sym); S.SetBase(f,f,f);

  const int axA = (axisDer + 1) % 3; // first transverse
  const int axB = (axisDer + 2) % 3; // second transverse

  double wA[4], wB[4], wd[4]; midpoint_w(wA); midpoint_w(wB); dmid_w(wd);
  const int o[4]  = {-2,-1, 0, 1};   // transverse offsets for midpoint interp
  const int ed[4] = {+1, 0,-1,-2};   // derivative sample offsets along axisDer

  for (int E=0; E<4; ++E) {
    const double wD = wd[E] * inv_h;
    for (int a=0; a<4; ++a) {
      for (int b=0; b<4; ++b) {
        const double coeff = wD * wA[a] * wB[b];
        int d[3] = {0,0,0};        // (di,dj,dk) in x,y,z
        d[axisDer] = ed[E];
        d[axA]     = o[a];
        d[axB]     = o[b];
        S.add(coeff, d[0], d[1], d[2]);
      }
    }
  }
  S.Simplify();
}

void InitCurlBStencils(CurlComp CurlBStencilSecondOrder[3], double dx, double dy, double dz) {
  // (curl B)_x = +∂Bz/∂y − ∂By/∂z
  CurlBStencilSecondOrder[0].Bx = cStencil(); CurlBStencilSecondOrder[0].Bx.SetSymbol("curl_Bx_from_Bx_4o");
  addFourthOrderPartial(CurlBStencilSecondOrder[0].Bz, /*d/dy*/1, 1.0/dy, "curl_Bx_from_Bz_4o");
  addFourthOrderPartial(CurlBStencilSecondOrder[0].By, /*d/dz*/2, 1.0/dz, "curl_Bx_from_By_4o");
  CurlBStencilSecondOrder[0].By *= -1.0;

  // (curl B)_y = +∂Bx/∂z − ∂Bz/∂x
  CurlBStencilSecondOrder[1].By = cStencil(); CurlBStencilSecondOrder[1].By.SetSymbol("curl_By_from_By_4o");
  addFourthOrderPartial(CurlBStencilSecondOrder[1].Bx, /*d/dz*/2, 1.0/dz, "curl_By_from_Bx_4o");
  addFourthOrderPartial(CurlBStencilSecondOrder[1].Bz, /*d/dx*/0, 1.0/dx, "curl_By_from_Bz_4o");
  CurlBStencilSecondOrder[1].Bz *= -1.0;

  // (curl B)_z = +∂By/∂x − ∂Bx/∂y
  CurlBStencilSecondOrder[2].Bz = cStencil(); CurlBStencilSecondOrder[2].Bz.SetSymbol("curl_Bz_from_Bz_4o");
  addFourthOrderPartial(CurlBStencilSecondOrder[2].By, /*d/dx*/0, 1.0/dx, "curl_Bz_from_By_4o");
  addFourthOrderPartial(CurlBStencilSecondOrder[2].Bx, /*d/dy*/1, 1.0/dy, "curl_Bz_from_Bx_4o");
  CurlBStencilSecondOrder[2].Bx *= -1.0;
}

} // namespace FourthOrder
} // namespace Stencil
} // namespace ECSIM
} // namespace FieldSolver
} // namespace PIC

