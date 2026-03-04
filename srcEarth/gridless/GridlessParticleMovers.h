//======================================================================================
// GridlessParticleMovers.h
//======================================================================================
// PURPOSE
//   Shared particle movers for gridless energetic-particle tools (cutoff rigidity,
//   density/spectrum, and other backtracing-based diagnostics).
//
// MOTIVATION
//   Historically, each gridless tool carried a local copy of a relativistic Boris
//   pusher (and related vector helpers). This led to code drift and made it harder
//   to introduce new movers consistently.
//
//   This header+source pair factors the movers into a single implementation that
//   can be used from multiple gridless solvers.
//
// DESIGN PRINCIPLES
//   - DO NOT perform any coordinate transformations here.
//     The movers operate strictly in the coordinate system used by the caller.
//     (In the current gridless tools, that is typically GSM in SI units.)
//   - Keep the mover API independent of any particular field model.
//     The caller supplies a field-evaluator callback object that returns B(x).
//   - Keep the original classic Boris as an option for regression continuity.
//   - Add a higher-accuracy option that reduces spatial field-sampling error:
//       BORIS_MIDPOINT: evaluate B at x_{n+1/2}.
//
// UNITS
//   - Position x: meters [m]
//   - Momentum p: kg*m/s
//   - Charge q: Coulombs [C]
//   - Rest mass m0: kg
//   - Magnetic field B: Tesla [T]
//   - dt: seconds [s]
//
// NOTE ABOUT COMMENTS
//   The gridless solvers contain extensive in-file documentation and historical
//   notes. When integrating this header into those solvers, keep those comments
//   in place. This module adds additional comments here, but does not remove
//   solver-level documentation.
//======================================================================================

#ifndef _AMPS_GRIDLESS_PARTICLE_MOVERS_H_
#define _AMPS_GRIDLESS_PARTICLE_MOVERS_H_

#include <cmath>
#include <string>

// The gridless solvers already rely on AMPS constants and relativistic helpers.
// We keep this include lightweight: SpeedOfLight is used in gamma evaluation.
#include "constants.h"

//------------------------------------------------------------------------------
// Shared small vector type
//------------------------------------------------------------------------------
// We intentionally expose the same V3 struct name used historically in the
// gridless solvers. This allows those compilation units to include this header
// and remove their duplicated V3 definitions without rewriting the whole file.
//
// If other parts of the code base define a different V3, do not include this
// header there. This header is meant for the gridless solvers only.
struct V3 { double x,y,z; };

static inline V3 add(const V3&a,const V3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline V3 mul(double s,const V3&a){return {s*a.x,s*a.y,s*a.z};}
static inline double dot(const V3&a,const V3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline V3 cross(const V3&a,const V3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm(const V3&a){return std::sqrt(dot(a,a));}
static inline V3 unit(const V3&a){double n=norm(a); return (n>0)?mul(1.0/n,a):V3{0,0,0};}

//------------------------------------------------------------------------------
// Field evaluator interface
//------------------------------------------------------------------------------
// We do NOT want movers to know anything about Tsyganenko/IGRF/dipole.
// The solvers already have an internal field evaluator class; they can adapt
// it to this interface with a small wrapper.
class IGridlessFieldEvaluator {
public:
  virtual ~IGridlessFieldEvaluator() = default;

  // Evaluate magnetic field at position x (meters), returning B in Tesla.
  virtual void GetB_T(const V3& x_m, V3& B_T) const = 0;
};

//------------------------------------------------------------------------------
// Mover selection
//------------------------------------------------------------------------------
// Keep the classic Boris as the baseline and add a higher-accuracy variant.
enum class MoverType {
  BORIS,          // Classic relativistic Boris, B sampled at x_n
  BORIS_MIDPOINT  // Boris with midpoint B sampling: B(x_{n+1/2})
};

//------------------------------------------------------------------------------
// Process-wide default mover
//------------------------------------------------------------------------------
// See header for intended precedence and usage.
extern MoverType gDefaultMover;


//------------------------------------------------------------------------------
// Default mover (shared across gridless tools)
//------------------------------------------------------------------------------
// Many AMPS entry points currently do not have direct access to the parsed
// AmpsParam object at the exact place where stepping occurs (e.g., deep inside
// worker loops). To support a simple CLI override *today* (and a parser-based
// configuration *later*), we provide a process-wide default mover.
//
// Intended precedence (recommended):
//   CLI option (explicit)  >  input-file setting  >  compiled default
//
// Tools can call SetDefaultMoverType(...) during initialization (after parsing
// CLI and/or input file). Then stepping code can use GetDefaultMoverType().
//
// NOTE: This is deliberately minimal and does not introduce any dependency on
// MPI or thread primitives. If you run multiple different solvers in the same
// process with different mover choices, you should pass MoverType explicitly
// instead of using this global default.
void SetDefaultMoverType(MoverType m);
MoverType GetDefaultMoverType();

// Helper utilities for parsing/printing mover names in CLI/help text.
// Supported names (case-insensitive):
//   "BORIS" and "BORIS_MIDPOINT" (also accepts "MIDPOINT" as shorthand).
bool ParseMoverType(const std::string& s, MoverType& out);
const char* MoverTypeToString(MoverType m);

//------------------------------------------------------------------------------
// Public API
//------------------------------------------------------------------------------
// The dispatcher is what solvers should call. They can select the mover type
// at runtime (later wired via parser / input file), but for now the solvers can
// keep a local variable that chooses which mover to use.
void StepParticle(MoverType mover,
                  V3& x_m, V3& p_SI,
                  double q_C, double m0_kg,
                  double dt,
                  const IGridlessFieldEvaluator& field);

// Convenience overload: step using the process-wide default mover.
// This is what most gridless solvers should use until the parser wiring lands.
inline void StepParticle(V3& x_m, V3& p_SI,
                         double q_C, double m0_kg,
                         double dt,
                         const IGridlessFieldEvaluator& field) {
  StepParticle(GetDefaultMoverType(), x_m, p_SI, q_C, m0_kg, dt, field);
}

// Optional: expose the individual movers for debugging/experimentation.
void BorisStep(V3& x_m, V3& p_SI,
               double q_C, double m0_kg,
               double dt,
               const IGridlessFieldEvaluator& field);

void BorisMidpointStep(V3& x_m, V3& p_SI,
                       double q_C, double m0_kg,
                       double dt,
                       const IGridlessFieldEvaluator& field);

#endif
