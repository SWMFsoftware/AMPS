//======================================================================================
// GridlessParticleMovers.h
//======================================================================================
#ifndef _AMPS_GRIDLESS_PARTICLE_MOVERS_H_
#define _AMPS_GRIDLESS_PARTICLE_MOVERS_H_

#include <cmath>
#include <string>
#include "constants.h"

struct V3 { double x,y,z; };

static inline V3 add(const V3&a,const V3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline V3 sub(const V3&a,const V3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline V3 mul(double s,const V3&a){return {s*a.x,s*a.y,s*a.z};}
static inline double dot(const V3&a,const V3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline V3 cross(const V3&a,const V3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm(const V3&a){return std::sqrt(dot(a,a));}
static inline V3 unit(const V3&a){double n=norm(a); return (n>0)?mul(1.0/n,a):V3{0,0,0};}

class IGridlessFieldEvaluator {
public:
  virtual ~IGridlessFieldEvaluator() = default;
  virtual void GetB_T(const V3& x_m, V3& B_T) const = 0;
};

enum class MoverType {
  BORIS,
  RK2,
  RK4,
  RK6
};

extern MoverType gDefaultMover;

void SetDefaultMoverType(MoverType m);
MoverType GetDefaultMoverType();
bool ParseMoverType(const std::string& s, MoverType& out);
const char* MoverTypeToString(MoverType m);

void StepParticle(MoverType mover,
                  V3& x_m, V3& p_SI,
                  double q_C, double m0_kg,
                  double dt,
                  const IGridlessFieldEvaluator& field);

inline void StepParticle(V3& x_m, V3& p_SI,
                         double q_C, double m0_kg,
                         double dt,
                         const IGridlessFieldEvaluator& field) {
  StepParticle(GetDefaultMoverType(), x_m, p_SI, q_C, m0_kg, dt, field);
}

// Return false when the integrator detects that the trajectory entered the inner
// sphere during the step (including any intermediate RK stage or segment between
// consecutive stage positions). In that case x/p are advanced to the point at
// which the crossing was detected or to the last stage state available.
bool StepParticleChecked(MoverType mover,
                         V3& x_m, V3& p_SI,
                         double q_C, double m0_kg,
                         double dt,
                         const IGridlessFieldEvaluator& field,
                         double rInner_m);

inline bool StepParticleChecked(V3& x_m, V3& p_SI,
                                double q_C, double m0_kg,
                                double dt,
                                const IGridlessFieldEvaluator& field,
                                double rInner_m) {
  return StepParticleChecked(GetDefaultMoverType(), x_m, p_SI, q_C, m0_kg, dt, field, rInner_m);
}

void BorisStep(V3& x_m, V3& p_SI,
               double q_C, double m0_kg,
               double dt,
               const IGridlessFieldEvaluator& field);

void RK2Step(V3& x_m, V3& p_SI,
             double q_C, double m0_kg,
             double dt,
             const IGridlessFieldEvaluator& field);

void RK4Step(V3& x_m, V3& p_SI,
             double q_C, double m0_kg,
             double dt,
             const IGridlessFieldEvaluator& field);

void RK6Step(V3& x_m, V3& p_SI,
             double q_C, double m0_kg,
             double dt,
             const IGridlessFieldEvaluator& field);

#endif
