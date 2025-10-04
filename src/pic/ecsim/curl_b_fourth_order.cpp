/**
 * @file    curl_b_fourth_order.cpp
 * @brief   Fourth-order curl(B) at cell corners from cell-centered B.
 *
 * Geometry & Scheme (summary)
 * ---------------------------
 * B is cell-centered (i,j,k). We construct (∇×B) at the corner (i-1/2,j-1/2,k-1/2).
 *
 * For a partial like ∂Bz/∂y|corner to be O(h^4):
 *   (A) 4th-order midpoint interpolation in the two *transverse* axes to move Bz
 *       from centers to the corner-directed half-indices using
 *         w_mid = [-1, 9, 9, -1] / 16 on offsets {-2,-1,0,1}.
 *   (B) 4th-order *half-index* derivative along the differentiation axis using
 *         w_dmid = (1/(24*h)) * [-1, 27, -27, 1] on offsets {+1, 0, -1, -2}.
 *
 * Tensoring (A)×(A)×(B) yields 4th-order accuracy at the corner. Repeat axis permutations
 * for other partials and combine with curl signs:
 *   (curl B)_x = +∂Bz/∂y − ∂By/∂z
 *   (curl B)_y = +∂Bx/∂z − ∂Bz/∂x
 *   (curl B)_z = +∂By/∂x − ∂Bx/∂y
 */

#include "../pic.h"  // brings in cStencil, cFrac, and cCurlBStencil

namespace PIC {
namespace FieldSolver {
namespace Electromagnetic {
namespace ECSIM {
namespace Stencil {
namespace FourthOrder {

// Use the project’s curl container type so the linker symbol matches exactly.
using CurlComp = PIC::FieldSolver::Electromagnetic::ECSIM::Stencil::cCurlBStencil;

// ---- 4-pt 4th-order midpoint interpolation weights: {m-2,m-1,m,m+1} → value at (m-1/2)
static inline void midpoint_w(double w[4]) {
  w[0] = -1.0/16.0;
  w[1] =  9.0/16.0;
  w[2] =  9.0/16.0;
  w[3] = -1.0/16.0;
}

/**
 * 4-pt, 4th-order derivative evaluated at a half-index (midpoint),
 * using integer samples at offsets {+1, 0, -1, -2} relative to the corner-directed midpoint.
 * When multiplied by inv_h, this yields (1/(24*h)) * [-1, 27, -27, 1].
 * This corrects the earlier centered-at-integer formula.
 */
static inline void dmid_w(double w[4]) {
  w[0] = -1.0/24.0; // +1
  w[1] = +27.0/24.0; //  0
  w[2] = -27.0/24.0; // -1
  w[3] = + 1.0/24.0; // -2
}

/**
 * Build one partial derivative at the corner by:
 *  - midpoint interpolation in the two transverse axes (4th order)
 *  - half-index derivative in the differentiation axis (4th order)
 *
 * axisDer: 0→d/dx, 1→d/dy, 2→d/dz
 * inv_h:   1/dx, 1/dy, or 1/dz accordingly
 * sym:     symbol string for debugging/printing
 */
static void addFourthOrderPartial(cStencil& S, int axisDer, double inv_h, const char* sym) {
  S = cStencil();
  const cFrac f(0,1);
  S.SetSymbol(sym);
  S.SetBase(f,f,f); // base at integer cell center; offsets are integers

  const int axA = (axisDer + 1) % 3; // first transverse
  const int axB = (axisDer + 2) % 3; // second transverse

  double wA[4], wB[4], wd[4];
  midpoint_w(wA);
  midpoint_w(wB);
  dmid_w(wd);

  // Transverse midpoint interpolation toward the (-1/2) corner
  const int o[4]  = {-2, -1,  0,  1};

  // Half-index derivative samples along derivative axis, aligned with dmid_w
  // order {+1, 0, -1, -2}
  const int ed[4] = {+1,  0, -1, -2};

  for (int E = 0; E < 4; ++E) {
    const double wD = wd[E] * inv_h;
    for (int a = 0; a < 4; ++a) {
      for (int b = 0; b < 4; ++b) {
        const double coeff = wD * wA[a] * wB[b];

        int d[3] = {0,0,0};  // (di,dj,dk) in x,y,z
        d[axisDer] = ed[E];
        d[axA]     = o[a];
        d[axB]     = o[b];

        S.add(coeff, d[0], d[1], d[2]);
      }
    }
  }
  S.Simplify();
}

/**
 * Initialize curl(B) stencils for all three components at the (-1/2,-1/2,-1/2) corner.
 * Other corners are obtained by index shifts at application time.
 *
 * (curl B)_x = +∂Bz/∂y − ∂By/∂z
 * (curl B)_y = +∂Bx/∂z − ∂Bz/∂x
 * (curl B)_z = +∂By/∂x − ∂Bx/∂y
 */
void InitCurlBStencils(CurlComp* CurlBStencilSecondOrder, double dx, double dy, double dz) {
  // (curl B)_x
  CurlBStencilSecondOrder[0].Bx = cStencil(); // empty holder for symmetry
  CurlBStencilSecondOrder[0].Bx.SetSymbol("curl_Bx_from_Bx_4o");
  addFourthOrderPartial(CurlBStencilSecondOrder[0].Bz, /*d/dy*/1, 1.0/dy, "curl_Bx_from_Bz_4o");
  addFourthOrderPartial(CurlBStencilSecondOrder[0].By, /*d/dz*/2, 1.0/dz, "curl_Bx_from_By_4o");
  CurlBStencilSecondOrder[0].By *= -1.0;

  // (curl B)_y
  CurlBStencilSecondOrder[1].By = cStencil(); // empty holder
  CurlBStencilSecondOrder[1].By.SetSymbol("curl_By_from_By_4o");
  addFourthOrderPartial(CurlBStencilSecondOrder[1].Bx, /*d/dz*/2, 1.0/dz, "curl_By_from_Bx_4o");
  addFourthOrderPartial(CurlBStencilSecondOrder[1].Bz, /*d/dx*/0, 1.0/dx, "curl_By_from_Bz_4o");
  CurlBStencilSecondOrder[1].Bz *= -1.0;

  // (curl B)_z
  CurlBStencilSecondOrder[2].Bz = cStencil(); // empty holder
  CurlBStencilSecondOrder[2].Bz.SetSymbol("curl_Bz_from_Bz_4o");
  addFourthOrderPartial(CurlBStencilSecondOrder[2].By, /*d/dx*/0, 1.0/dx, "curl_Bz_from_By_4o");
  addFourthOrderPartial(CurlBStencilSecondOrder[2].Bx, /*d/dy*/1, 1.0/dy, "curl_Bz_from_Bx_4o");
  CurlBStencilSecondOrder[2].Bx *= -1.0;
}

} // namespace FourthOrder
} // namespace Stencil
} // namespace ECSIM
} // namespace Electromagnetic
} // namespace FieldSolver
} // namespace PIC

