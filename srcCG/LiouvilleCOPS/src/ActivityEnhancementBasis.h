// Versioned, fully normalized real spherical harmonics for the optional
// localized H2O activity enhancement.  The normative equations, coordinate
// convention, coefficient order, and relationship to the legacy srcCG basis
// are in HARMONIC_BASIS.md.
//
// The mathematical component remains independent of cSourceModel.  Stage A3
// composes it through ActivityEnhancement.* only when a validated optional
// coefficient file is enabled; the default source path does not call it.
#ifndef ACTIVITY_ENHANCEMENT_BASIS_H
#define ACTIVITY_ENHANCEMENT_BASIS_H

#include <cstddef>
#include <vector>

namespace ActivityEnhancementBasis {

constexpr const char* BasisId="liouville_cops.real_orthonormal_cs_v1";
constexpr int MinSupportedDegree=0;
constexpr int MaxSupportedDegree=8;

struct Mode {
  int Degree;
  int Order;
};

// Validate 0 <= minDegree <= maxDegree <= 8 and return the inclusive number
// of modes: (maxDegree+1)^2-minDegree^2.  Invalid inputs throw
// std::invalid_argument instead of terminating the application.
std::size_t CoefficientCount(int minDegree,int maxDegree);

// Normative increasing-degree, then increasing signed-order sequence.
std::vector<Mode> Modes(int minDegree,int maxDegree);

// Zero-based location of (degree,order) in Modes(minDegree,maxDegree).
std::size_t CoefficientIndex(int degree,int order,int minDegree,int maxDegree);

// Evaluate all modes at a nonzero Cartesian vector in the mesh/body frame.
// Only direction matters.  At x=y=0 the implementation sets longitude to
// zero; every |m|>0 function is zero at a pole.  Nonfinite/zero coordinates
// and unsupported degrees throw std::invalid_argument.
void EvaluateCartesian(double x,double y,double z,int minDegree,int maxDegree,
                       double* output,std::size_t outputSize);

// Owning convenience overload with the same ordering and validation.
std::vector<double> EvaluateCartesian(double x,double y,double z,
                                      int minDegree=5,int maxDegree=8);

} // namespace ActivityEnhancementBasis

#endif
