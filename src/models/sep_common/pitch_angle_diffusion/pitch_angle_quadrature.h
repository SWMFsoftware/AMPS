#ifndef SEP_COMMON_PITCH_ANGLE_DIFFUSION_PITCH_ANGLE_QUADRATURE_H
#define SEP_COMMON_PITCH_ANGLE_DIFFUSION_PITCH_ANGLE_QUADRATURE_H

// ============================================================================
// Internal header (not part of the public API): adaptive Gauss-Kronrod (7-15)
// quadrature shared by pitch_angle_diffusion.cpp (normalization integrals)
// and pitch_angle_kernels.cpp (finite-volume transmissibilities).
//
// Global adaptive bisection of the segment with the largest error estimate
// (the QUADPACK QAG strategy).  The error of a segment is the conservative
// |K15 - G7|; the total is accepted when it is below max(absTol, relTol |I|)
// of the *whole* integral, so an integrable algebraic endpoint behaviour
// (for example t^(2/3)) is resolved by repeated bisection toward it instead
// of failing a local relative test.  Accumulation in long double.  The
// integrand must be finite inside every segment; nonintegrable endpoints are
// removed by the callers' substitutions.  Failure to reach the tolerance
// within the segment budget is returned, never hidden.
// ============================================================================

#include <algorithm>
#include <cmath>
#include <vector>

namespace SEP {
namespace PitchAngleDiffusion {
namespace Quadrature {

constexpr long double kXgk[8] = {
    0.991455371120812639206854697526329L, 0.949107912342758524526189684047851L,
    0.864864423359769072789712788640926L, 0.741531185599394439863864773280788L,
    0.586087235467691130294144845693013L, 0.405845151377397166906606412076961L,
    0.207784955007898467600689403773245L, 0.000000000000000000000000000000000L};
constexpr long double kWgk[8] = {
    0.022935322010529224963732008058970L, 0.063092092629978553290700663189204L,
    0.104790010322250183839876322541518L, 0.140653259715525918745189590510238L,
    0.169004726639267902826583426598550L, 0.190350578064785409913256402421014L,
    0.204432940075298892414161999234649L, 0.209482141084727828012999174891714L};
constexpr long double kWg[4] = {
    0.129484966168869693270611432679082L, 0.279705391489276667901467771423780L,
    0.381830050505118944950369775488975L, 0.417959183673469387755102040816327L};

struct Segment {
  long double a, b, value, error;
};

template <class F>
inline Segment KronrodSegment(const F& f, long double a, long double b) {
  const long double centre = 0.5L * (a + b);
  const long double half = 0.5L * (b - a);
  const long double fc = f(centre);
  long double kronrod = kWgk[7] * fc;
  long double gauss = kWg[3] * fc;
  for (int j = 0; j < 7; ++j) {
    const long double x = half * kXgk[j];
    const long double sum = f(centre - x) + f(centre + x);
    kronrod += kWgk[j] * sum;
    if (j % 2 == 1) gauss += kWg[j / 2] * sum;
  }
  Segment s;
  s.a = a;
  s.b = b;
  s.value = kronrod * half;
  s.error = std::fabs((kronrod - gauss) * half);
  return s;
}

// Integrate f over the ordered breakpoints (at least two).  Returns false if
// the tolerance max(absTol, relTol |I|) is not reached within maxSegments.
template <class F>
inline bool AdaptiveIntegrate(const F& f, const std::vector<long double>& points,
                       long double relTol, long double absTol,
                       int maxSegments, long double* value,
                       long double* error) {
  std::vector<Segment> segments;
  for (std::size_t i = 0; i + 1 < points.size(); ++i)
    if (points[i + 1] > points[i])
      segments.push_back(KronrodSegment(f, points[i], points[i + 1]));
  for (;;) {
    long double total = 0.0L, totalError = 0.0L;
    std::size_t worst = 0;
    for (std::size_t i = 0; i < segments.size(); ++i) {
      total += segments[i].value;
      totalError += segments[i].error;
      if (segments[i].error > segments[worst].error) worst = i;
    }
    if (!std::isfinite(static_cast<double>(total))) return false;
    if (totalError <= std::max(absTol, relTol * std::fabs(total))) {
      *value = total;
      *error = totalError;
      return true;
    }
    if (static_cast<int>(segments.size()) >= maxSegments) return false;
    const Segment bad = segments[worst];
    const long double middle = 0.5L * (bad.a + bad.b);
    if (!(middle > bad.a && middle < bad.b)) return false;
    segments[worst] = KronrodSegment(f, bad.a, middle);
    segments.push_back(KronrodSegment(f, middle, bad.b));
  }
}

}  // namespace Quadrature
}  // namespace PitchAngleDiffusion
}  // namespace SEP

#endif  // SEP_COMMON_PITCH_ANGLE_DIFFUSION_PITCH_ANGLE_QUADRATURE_H
