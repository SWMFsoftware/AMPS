#ifndef SEP_CORONAL_CME_INTERFACE_BALANCE_H
#define SEP_CORONAL_CME_INTERFACE_BALANCE_H

#include "sep_coronal_cme/vector_math.h"
#include "sep_status.h"

#include <string>
#include <vector>

namespace SEP { namespace CoronalCME {

enum class InterfaceKind { OpenClosedSeparatrix, GenericOpenOpen };
enum class BalancePolicy { DiagnosticKinematic, BoundedApproximation,
                           StationaryTangentialDiscontinuity };
enum class BalanceRepresentation { SharpOneSided, FiniteWidthVolume };
enum class StateOrigin { AnalyticComposite, SolverProduced, Imported };

struct PrimitiveState {
  double densityKgM3=0.0,pressurePa=0.0,fastSpeedMPerS=0.0;
  Vec3 velocityMPerS,magneticFieldT;
};
struct InterfaceTolerances {
  double magneticAbsoluteT=0.0,magneticRelative=0.0;
  double speedAbsoluteMPerS=0.0,speedMach=0.0;
  double massFluxAbsolute=0.0,massFluxRelative=0.0,massFluxFloor=1e-30;
  double tractionAbsolutePa=0.0,tractionRelative=0.0,tractionFloorPa=1e-30;
};
struct SharpBalance {
  Vec3 tractionJumpPa;
  double tractionNormPa=0.0,signedMassFluxJumpKgM2S=0.0;
  bool magneticGatePassed=false,kinematicGatePassed=false;
  bool massFluxGatePassed=false,tractionGatePassed=false,policyPassed=false;
};
struct VolumeBalance { Vec3 residualNPerM3; double normNPerM3=0.0; };

Core::Result<SharpBalance> EvaluateSharpInterface(InterfaceKind kind,
    BalancePolicy policy,StateOrigin originA,StateOrigin originB,
    const PrimitiveState& sideA,const PrimitiveState& sideB,Vec3 normal,
    double interfaceNormalSpeedMPerS,const InterfaceTolerances& tolerances,
    bool hasUncertaintyAndThreeLevelConvergence=false);
Core::Result<VolumeBalance> EvaluateVolumeMomentumResidual(Vec3 timeDerivative,
    Vec3 fluxDivergence,Vec3 gravityForce,Vec3 declaredForce);

} }  // namespace SEP::CoronalCME
#endif
