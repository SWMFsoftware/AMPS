/*
================================================================================
 example_units.cpp — End-to-end SI↔PIC normalization demonstration & checks
--------------------------------------------------------------------------------
WHAT THIS PROGRAM DOES
  • Builds a single normalization regime (Factors) from chosen SI scales (ℓ, u, m).
  • Converts two physically motivated input states (Case 1 and Case 2) from SI →
    PIC-normalized units using the identical formulas used in FLEKS/fluid.
  • Prints the normalized values and compares them to target outputs you provided
    (Ux_no, dt_no, Bz_no, Ey_no). Soft warnings are emitted if values deviate
    beyond tight tolerances (intended as quick regressions, not hard asserts).

NORMALIZATION IN A NUTSHELL
  Inputs:  ℓ_SI [m], u_SI [m/s], m_SI [kg]
  CGS:     L0=100ℓ_SI [cm], U0=100u_SI [cm/s], M0=m_SI [g]
  Derived: ρ0=M0/L0^3, B0=√(ρ0)U0, Q0=√(M0L0)U0, P0=ρ0U0^2, J0=Q0U0/L0^3, T0=L0/U0
  Factors: multiply SI by Si2NoX to get normalized; inverses are No2SiX.
  Time:    Si2NoT = Si2NoL/Si2NoV = 1/T0, so dt_no = dt_SI/T0 and vice versa.

WHY THESE SCALES?
  We back-solve the scales from Case 1 targets to make the example transparent:
   • From Ux_no = v_SI/u_SI ⇒ u_SI = 26.4 km/s / 8.8e-3 = 3.0e6 m/s.
   • From dt_no = dt_SI·(u_SI/ℓ_SI) = 600 with dt_SI=0.2 s ⇒ ℓ_SI = 1.0e3 m.
   • Mass scale m_SI is chosen so that Bz_no (for 325 nT) matches the target,
     using B0=√(ρ0)U0 with ρ0=M0/L0^3; this yields m_SI≈1.0898097712911516e-3 kg.
  With these scales, Case 2 automatically lands on the requested normalized values
  (Ux_no=1e-2, dt_no=600, Bz_no≈9.5791e-3, Ey_no≈9.5791e-5).

EXPECTED OUTPUT (abridged)
  Case 1:
    Ux_no ≈ -8.800000e-03, dt_no ≈ 6.000000e+02,
    Bz_no ≈ 1.0377357e-02,  Ey_no ≈ -9.1320740e-05
  Case 2:
    Ux_no ≈  1.000000e-02,  dt_no ≈ 6.000000e+02,
    Bz_no ≈  9.5790988e-03, Ey_no ≈  9.5790987e-05

HOW E IS CHECKED
  We do not solve Maxwell here; instead we use the convective form consistent
  with the chosen normalization: Ey_no = ±Ux_no * Bz_no (sign differs by case).

BUILD & RUN
  c++ -O2 -std=c++17 example_units.cpp -o example_units
  ./example_units

INTEGRATION NOTES
  • The conversions here are precisely those in pic_units_normalization.hpp.
  • For production tests, swap soft warnings with <cassert> or a test framework.
  • To compare to code output, ensure the same (ℓ, u, m) are used at runtime.
================================================================================
*/



#include <iostream>
#include <iomanip>
#include <cmath>
#include "pic_units_normalization.h"
using namespace picunits;


int main(){
  std::cout.setf(std::ios::scientific);
  std::cout << std::setprecision(9);


  // Unified normalization scales chosen to satisfy both cases
  // Deduced from Case 1 targets (see header comments):
  const double uSI = 3.0e6; // [m/s]
  const double lSI = 1.0e3; // [m]
  const double mSI = 1.0898097712911516e-3; // [kg]


  auto F = build({lSI, uSI, mSI});

  auto almost_eq = [](double a, double b, double tol){
    return std::abs(a-b) <= tol * std::max(1.0, std::abs(b));
  };


  // ===============================
  // Case 1 (from user description)
  // ===============================
  const double U1_SI[3] = {-26.4e3, 0.0, 0.0}; // m/s
  const double B1_SI[3] = {0.0, 0.0, 325e-9}; // Tesla
  const double dt1_SI = 0.2; // s

  double U1_no[3], B1_no[3];
  si2no_v3(U1_SI, U1_no, F);
  si2no_B3(B1_SI, B1_no, F);
  const double dt1_no = si2no_t(dt1_SI, F);
  const double Ey1_no = (U1_no[0] * B1_no[2]);


  std::cout << "\n\nCase 1 (targets in parentheses)\n";
  std::cout << " Ux_no = " << U1_no[0] << " ( -8.800000e-03 )\n";
  std::cout << " dt_no = " << dt1_no << " ( 6.000000e+02 )\n";
  std::cout << " Bz_no = " << B1_no[2] << " ( 1.0377357e-02 )\n";
  std::cout << " Ey_no = " << Ey1_no << " ( -9.1320740e-05 )\n";


  if(!almost_eq(U1_no[0], -8.8e-3, 1e-9)) std::cerr << "WARN: Case1 Ux_no mismatch\n";
  if(!almost_eq(dt1_no, 600.0, 1e-12)) std::cerr << "WARN: Case1 dt_no mismatch\n";
  if(!almost_eq(B1_no[2], 0.010377357, 1e-9)) std::cerr << "WARN: Case1 Bz_no mismatch\n";
  if(!almost_eq(Ey1_no, -9.1320740e-05, 1e-9))std::cerr << "WARN: Case1 Ey_no mismatch\n";


  // ===============================
  // Case 2 (from user description)
  // ===============================
  const double U2_SI[3] = {+30.0e3, 0.0, 0.0}; // m/s
  const double B2_SI[3] = {0.0, 0.0, 300e-9}; // Tesla
  const double dt2_SI = 0.2; // s


  double U2_no[3], B2_no[3];
  si2no_v3(U2_SI, U2_no, F);
  si2no_B3(B2_SI, B2_no, F);
  const double dt2_no = si2no_t(dt2_SI, F);
  const double Ey2_no = (U2_no[0] * B2_no[2]);


  std::cout << "\n\nCase 2 (targets in parentheses)\n";
  std::cout << std::setprecision(10);
  std::cout << " Ux_no = " << U2_no[0] << " ( 1.000000e-02 )\n";
  std::cout << " dt_no = " << dt2_no << " ( 6.000000e+02 )\n";
  std::cout << " Bz_no = " << B2_no[2] << " ( 9.5790988e-03 )\n";
  std::cout << " Ey_no = " << Ey2_no << " ( 9.5790987e-05 )\n";


  if(!almost_eq(U2_no[0], 1.0e-2, 1e-12)) std::cerr << "WARN: Case2 Ux_no mismatch\n";
  if(!almost_eq(dt2_no, 600.0, 1e-12)) std::cerr << "WARN: Case2 dt_no mismatch\n";
  if(!almost_eq(B2_no[2], 9.5790988e-03, 1e-10)) std::cerr << "WARN: Case2 Bz_no mismatch\n";
  if(!almost_eq(Ey2_no, 9.5790987e-05, 1e-10)) std::cerr << "WARN: Case2 Ey_no mismatch\n";

  return 0;
}
