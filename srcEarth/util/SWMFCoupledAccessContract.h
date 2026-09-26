#ifndef _SRC_EARTH_UTIL_SWMF_COUPLED_ACCESS_CONTRACT_H_
#define _SRC_EARTH_UTIL_SWMF_COUPLED_ACCESS_CONTRACT_H_

//======================================================================================
// SWMFCoupledAccessContract.h
//======================================================================================
//
// Dependency-free contracts used by Roadmap Step 10, "Coupled cutoff/access
// production path".
//
// The live SWMF bridge depends on AMPS, MPI, and the SWMF coupler.  Scheduling,
// product naming, and outer-boundary classification do not.  Keeping those decisions
// in this small header gives the numerical tests an independent way to exercise the
// exact production rules without mocking the complete coupled executable.
//
// Three mistakes are deliberately made impossible here:
//
//   1. Cadence is based on callback count instead of authoritative simulation time.
//      CadenceGate accepts only finite, non-negative time, never advances on a skipped
//      callback, and rejects an in-process clock rollback.  A fresh gate represents a
//      restarted process and therefore evaluates the restart epoch reproducibly.
//
//   2. Two different SWMF states overwrite one another.  BuildProductSuffix() binds
//      the exact simulation time and complete snapshot ID into every coupled product
//      name.  Callback/rank/thread counts are intentionally absent because they are
//      execution details and would make restart/layout comparisons ambiguous.
//
//   3. Missing AMR data are interpreted as physical particle access.  BOX means that
//      the six computational faces are the physical escape surface.  SHUE means that
//      the Shue magnetopause (plus the configured nightside X cap) is physical; hitting
//      another computational face while still inside that surface is MESH_UNAVAILABLE,
//      never ALLOWED.
//
// All distances are SI metres except the conventional Shue r0 parameter, which is
// expressed in Earth radii.  Dynamic pressure is nPa and IMF Bz is nT in the published
// Shue et al. parameterization.
//======================================================================================

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace Earth {
namespace SWMFCoupledAccess {

//--------------------------------------------------------------------------------------
// Collective cadence state
//--------------------------------------------------------------------------------------

enum class CadenceAction {
  Run,
  Skip,
  RejectStaleTime
};

struct CadenceDecision {
  CadenceAction action{CadenceAction::Skip};
  double simulationTime_s{0.0};
  double lastCompletedTime_s{-1.0};
  double nextDueTime_s{0.0};
  std::string reason;
};

inline double CadenceToleranceSeconds(double cadence_s) {
  // This is only a floating-point comparison allowance.  It is ten orders of
  // magnitude smaller than a one-second cadence and is never added to the physical
  // time stored in product provenance.
  return 1.0e-10*std::max(1.0,std::fabs(cadence_s));
}

class CadenceGate {
 public:
  CadenceGate()=default;

  bool HasCompletedEpoch() const { return hasCompletedEpoch_; }
  double LastCompletedTimeSeconds() const {
    return hasCompletedEpoch_ ? lastCompletedTime_s_ : -1.0;
  }

  CadenceDecision Evaluate(double simulationTime_s,double cadence_s) const {
    if (!std::isfinite(simulationTime_s) || simulationTime_s<0.0)
      throw std::invalid_argument(
          "SWMF coupled cadence requires finite non-negative simulation time");
    if (!std::isfinite(cadence_s))
      throw std::invalid_argument("SWMF coupled cadence must be finite");

    CadenceDecision out;
    out.simulationTime_s=simulationTime_s;
    out.lastCompletedTime_s=LastCompletedTimeSeconds();

    if (!hasCompletedEpoch_) {
      out.action=CadenceAction::Run;
      out.nextDueTime_s=simulationTime_s;
      out.reason="first complete SWMF snapshot in this process";
      return out;
    }

    const double tolerance_s=CadenceToleranceSeconds(cadence_s);
    const double elapsed_s=simulationTime_s-lastCompletedTime_s_;
    if (elapsed_s < -tolerance_s) {
      out.action=CadenceAction::RejectStaleTime;
      out.nextDueTime_s=lastCompletedTime_s_;
      out.reason="SWMF simulation time moved backward after a completed product";
      return out;
    }

    // A non-positive cadence is the documented explicit request to evaluate every
    // callback.  An exact duplicate time is still skipped inside one process: it is
    // the same received state, not a new physical epoch.  A restarted process owns a
    // fresh gate, so that same epoch is evaluated again for the restart parity test.
    if (elapsed_s<=tolerance_s) {
      out.action=CadenceAction::Skip;
      out.nextDueTime_s=(cadence_s>0.0) ? lastCompletedTime_s_+cadence_s
                                       : lastCompletedTime_s_;
      out.reason="duplicate callback for the last completed SWMF epoch";
      return out;
    }

    if (cadence_s<=0.0 || elapsed_s+tolerance_s>=cadence_s) {
      out.action=CadenceAction::Run;
      out.nextDueTime_s=simulationTime_s;
      out.reason=(cadence_s<=0.0) ? "run-every-callback cadence" :
                                   "requested physical cadence reached";
      return out;
    }

    out.action=CadenceAction::Skip;
    out.nextDueTime_s=lastCompletedTime_s_+cadence_s;
    out.reason="requested physical cadence has not been reached";
    return out;
  }

  void CommitCompleted(double simulationTime_s) {
    if (!std::isfinite(simulationTime_s) || simulationTime_s<0.0)
      throw std::invalid_argument(
          "cannot commit a non-finite or negative SWMF product epoch");
    if (hasCompletedEpoch_ &&
        simulationTime_s+CadenceToleranceSeconds(1.0)<lastCompletedTime_s_)
      throw std::logic_error("cannot commit an SWMF epoch older than the last product");
    hasCompletedEpoch_=true;
    lastCompletedTime_s_=simulationTime_s;
  }

  // Used only when the owning executable deliberately begins a new restart
  // generation.  Normal callbacks must never reset the gate to bypass cadence.
  void ResetForRestart() {
    hasCompletedEpoch_=false;
    lastCompletedTime_s_=-1.0;
  }

 private:
  bool hasCompletedEpoch_{false};
  double lastCompletedTime_s_{-1.0};
};

//--------------------------------------------------------------------------------------
// Deterministic product identity
//--------------------------------------------------------------------------------------

inline std::string SanitizeSnapshotId(const std::string& snapshotId) {
  if (snapshotId.empty())
    throw std::invalid_argument("coupled cutoff product requires a snapshot ID");
  std::string out;
  out.reserve(snapshotId.size());
  for (std::string::const_iterator it=snapshotId.begin();it!=snapshotId.end();++it) {
    const unsigned char c=static_cast<unsigned char>(*it);
    out.push_back((std::isalnum(c) || c=='-' || c=='_') ? static_cast<char>(c) : '_');
  }
  return out;
}

inline std::string SimulationTimeToken(double simulationTime_s) {
  if (!std::isfinite(simulationTime_s) || simulationTime_s<0.0)
    throw std::invalid_argument(
        "coupled product suffix requires finite non-negative simulation time");
  std::ostringstream raw;
  raw.imbue(std::locale::classic());
  raw << std::fixed << std::setprecision(9) << simulationTime_s;
  std::string token=raw.str();
  const std::size_t minimumWidth=20; // 10 integer + decimal point + 9 fractional digits
  if (token.size()<minimumWidth) token.insert(0,minimumWidth-token.size(),'0');
  return token;
}

inline std::string BuildProductSuffix(double simulationTime_s,
                                      const std::string& snapshotId) {
  return ".swmf_t"+SimulationTimeToken(simulationTime_s)+"s_sid"+
         SanitizeSnapshotId(snapshotId);
}

//--------------------------------------------------------------------------------------
// BOX/Shue physical escape versus unavailable mesh
//--------------------------------------------------------------------------------------

enum class OuterBoundaryKind {
  Box,
  Shue
};

enum class PositionDisposition {
  Interior,
  PhysicalEscape,
  MeshUnavailable,
  Invalid
};

enum class SegmentDisposition {
  NoExit,
  PhysicalEscape,
  MeshUnavailable,
  Invalid
};

struct Box {
  std::array<double,3> minimum_m{{0.0,0.0,0.0}};
  std::array<double,3> maximum_m{{0.0,0.0,0.0}};
};

struct ShueParameters {
  double r0_Re{0.0};
  double alpha{0.0};
  double earthRadius_m{0.0};
  double tailCapX_m{0.0};
  bool r0WasAuto{false};
  bool alphaWasAuto{false};
};

struct BoundaryPolicy {
  OuterBoundaryKind kind{OuterBoundaryKind::Box};
  Box computationalBox;
  ShueParameters shue;
};

struct SegmentResult {
  SegmentDisposition disposition{SegmentDisposition::NoExit};
  double fraction{1.0};
  std::array<double,3> position_m{{0.0,0.0,0.0}};
  std::string surface{"NONE"};
};

inline bool IsFinitePoint(const std::array<double,3>& x) {
  return std::isfinite(x[0]) && std::isfinite(x[1]) && std::isfinite(x[2]);
}

inline void ValidateBox(const Box& box) {
  for (int d=0;d<3;++d) {
    if (!std::isfinite(box.minimum_m[d]) || !std::isfinite(box.maximum_m[d]) ||
        !(box.maximum_m[d]>box.minimum_m[d]))
      throw std::invalid_argument(
          "SWMF cutoff computational box must have finite ordered bounds");
  }
}

inline bool InsideBox(const Box& box,const std::array<double,3>& x,
                      double tolerance_m=0.0) {
  if (!IsFinitePoint(x)) return false;
  for (int d=0;d<3;++d) {
    if (x[d]<box.minimum_m[d]-tolerance_m ||
        x[d]>box.maximum_m[d]+tolerance_m) return false;
  }
  return true;
}

inline std::string UpperTrimmed(std::string value) {
  std::size_t first=0,last=value.size();
  while (first<last && std::isspace(static_cast<unsigned char>(value[first]))) ++first;
  while (last>first && std::isspace(static_cast<unsigned char>(value[last-1]))) --last;
  value=value.substr(first,last-first);
  std::transform(value.begin(),value.end(),value.begin(),
                 [](unsigned char c){ return static_cast<char>(std::toupper(c)); });
  return value;
}

inline double ParseStrictPositive(const std::string& token,const char* label) {
  std::size_t used=0;
  double value=0.0;
  try { value=std::stod(token,&used); }
  catch (const std::exception&) {
    throw std::invalid_argument(std::string(label)+" must be AUTO or a positive number");
  }
  while (used<token.size() &&
         std::isspace(static_cast<unsigned char>(token[used]))) ++used;
  if (used!=token.size() || !std::isfinite(value) || !(value>0.0))
    throw std::invalid_argument(std::string(label)+" must be AUTO or a positive number");
  return value;
}

// Shue et al. (1998) empirical magnetopause coefficients.  These formulas provide an
// independently checkable analytic surface for the Step-10 boundary tests; they do not
// infer a magnetopause from an arbitrary numerical contour in the MHD solution.
inline double ShueAutoR0Re(double dynamicPressure_nPa,double imfBz_nT) {
  if (!std::isfinite(dynamicPressure_nPa) || !(dynamicPressure_nPa>0.0) ||
      !std::isfinite(imfBz_nT))
    throw std::invalid_argument("SHUE AUTO requires PDYN>0 nPa and finite IMF_BZ");
  return (10.22+1.29*std::tanh(0.184*(imfBz_nT+8.14)))*
         std::pow(dynamicPressure_nPa,-1.0/6.6);
}

inline double ShueAutoAlpha(double dynamicPressure_nPa,double imfBz_nT) {
  if (!std::isfinite(dynamicPressure_nPa) || !(dynamicPressure_nPa>0.0) ||
      !std::isfinite(imfBz_nT))
    throw std::invalid_argument("SHUE AUTO requires PDYN>0 nPa and finite IMF_BZ");
  return (0.58-0.007*imfBz_nT)*(1.0+0.024*std::log(dynamicPressure_nPa));
}

inline ShueParameters ResolveShueParameters(
    const std::string& r0Token,const std::string& alphaToken,
    double dynamicPressure_nPa,double imfBz_nT,double earthRadius_m,
    double tailCapX_m) {
  if (!std::isfinite(earthRadius_m) || !(earthRadius_m>0.0) ||
      !std::isfinite(tailCapX_m))
    throw std::invalid_argument("invalid Earth radius or Shue tail cap");

  ShueParameters out;
  out.earthRadius_m=earthRadius_m;
  out.tailCapX_m=tailCapX_m;

  const std::string r0=UpperTrimmed(r0Token.empty() ? "AUTO" : r0Token);
  const std::string alpha=UpperTrimmed(alphaToken.empty() ? "AUTO" : alphaToken);
  out.r0WasAuto=(r0=="AUTO");
  out.alphaWasAuto=(alpha=="AUTO");
  out.r0_Re=out.r0WasAuto ? ShueAutoR0Re(dynamicPressure_nPa,imfBz_nT) :
                            ParseStrictPositive(r0Token,"SHUE_R0");
  out.alpha=out.alphaWasAuto ? ShueAutoAlpha(dynamicPressure_nPa,imfBz_nT) :
                              ParseStrictPositive(alphaToken,"SHUE_ALPHA");
  if (!std::isfinite(out.r0_Re) || !(out.r0_Re>0.0) ||
      !std::isfinite(out.alpha) || !(out.alpha>0.0))
    throw std::invalid_argument("resolved Shue parameters must be finite and positive");
  return out;
}

inline void ValidateBoundaryPolicy(const BoundaryPolicy& policy) {
  ValidateBox(policy.computationalBox);
  if (policy.kind==OuterBoundaryKind::Shue) {
    if (!std::isfinite(policy.shue.r0_Re) || !(policy.shue.r0_Re>0.0) ||
        !std::isfinite(policy.shue.alpha) || !(policy.shue.alpha>0.0) ||
        !std::isfinite(policy.shue.earthRadius_m) ||
        !(policy.shue.earthRadius_m>0.0) ||
        !std::isfinite(policy.shue.tailCapX_m))
      throw std::invalid_argument("invalid Shue boundary parameters");
  }
}

inline double ShueRadiusMeters(const ShueParameters& shue,double cosTheta) {
  cosTheta=std::max(-1.0,std::min(1.0,cosTheta));
  const double denominator=1.0+cosTheta;
  if (denominator<=64.0*std::numeric_limits<double>::epsilon())
    return std::numeric_limits<double>::infinity();
  return shue.r0_Re*shue.earthRadius_m*
         std::pow(2.0/denominator,shue.alpha);
}

inline double ShueSignedMarginMeters(const ShueParameters& shue,
                                     const std::array<double,3>& x) {
  if (!IsFinitePoint(x)) return -std::numeric_limits<double>::infinity();
  const double r=std::sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  if (r==0.0) return shue.r0_Re*shue.earthRadius_m;
  return ShueRadiusMeters(shue,x[0]/r)-r;
}

inline PositionDisposition ClassifyPosition(const BoundaryPolicy& policy,
                                            const std::array<double,3>& x,
                                            bool meshSampleAvailable=true,
                                            double tolerance_m=0.0) {
  ValidateBoundaryPolicy(policy);
  if (!IsFinitePoint(x)) return PositionDisposition::Invalid;

  const bool insideBox=InsideBox(policy.computationalBox,x,tolerance_m);
  if (policy.kind==OuterBoundaryKind::Box) {
    if (!insideBox) return PositionDisposition::PhysicalEscape;
    return meshSampleAvailable ? PositionDisposition::Interior :
                                 PositionDisposition::MeshUnavailable;
  }

  // The negative-X face is the explicit nightside cap for the open Shue surface.
  // Other computational faces are merely limits of available numerical data.
  if (x[0]<=policy.shue.tailCapX_m+tolerance_m)
    return PositionDisposition::PhysicalEscape;
  if (ShueSignedMarginMeters(policy.shue,x)<=tolerance_m)
    return PositionDisposition::PhysicalEscape;
  if (!insideBox || !meshSampleAvailable)
    return PositionDisposition::MeshUnavailable;
  return PositionDisposition::Interior;
}

inline std::array<double,3> InterpolatePoint(const std::array<double,3>& a,
                                             const std::array<double,3>& b,
                                             double fraction) {
  std::array<double,3> out{{0.0,0.0,0.0}};
  for (int d=0;d<3;++d) out[d]=a[d]+fraction*(b[d]-a[d]);
  return out;
}

inline bool FirstBoxExitFraction(const Box& box,const std::array<double,3>& a,
                                 const std::array<double,3>& b,double& fraction,
                                 std::string& face) {
  if (!InsideBox(box,a,0.0)) return false;
  double best=2.0;
  const char* bestFace="NONE";
  for (int d=0;d<3;++d) {
    const double delta=b[d]-a[d];
    if (delta>0.0) {
      const double t=(box.maximum_m[d]-a[d])/delta;
      if (t>=0.0 && t<=1.0 && t<best) {
        best=t; bestFace=(d==0 ? "XMAX" : (d==1 ? "YMAX" : "ZMAX"));
      }
    }
    else if (delta<0.0) {
      const double t=(box.minimum_m[d]-a[d])/delta;
      if (t>=0.0 && t<=1.0 && t<best) {
        best=t; bestFace=(d==0 ? "XMIN" : (d==1 ? "YMIN" : "ZMIN"));
      }
    }
  }
  if (best>1.0 || InsideBox(box,b,0.0)) return false;
  fraction=best;
  face=bestFace;
  return true;
}

inline bool FirstShueExitFraction(const ShueParameters& shue,
                                  const std::array<double,3>& a,
                                  const std::array<double,3>& b,
                                  double& fraction) {
  double previous=ShueSignedMarginMeters(shue,a);
  if (!(previous>0.0)) { fraction=0.0; return true; }

  // A particle chord can be long when the user selects a large fixed time step.
  // Search fixed subintervals for the first sign change instead of assuming the
  // nonlinear Shue margin is monotone over the entire chord.
  const int brackets=64;
  double lo=0.0,hi=0.0;
  bool found=false;
  for (int n=1;n<=brackets;++n) {
    const double t=static_cast<double>(n)/static_cast<double>(brackets);
    const double current=ShueSignedMarginMeters(shue,InterpolatePoint(a,b,t));
    if (current<=0.0 && previous>0.0) {
      lo=static_cast<double>(n-1)/static_cast<double>(brackets);
      hi=t;
      found=true;
      break;
    }
    previous=current;
  }
  if (!found) return false;

  for (int iteration=0;iteration<64;++iteration) {
    const double mid=0.5*(lo+hi);
    if (ShueSignedMarginMeters(shue,InterpolatePoint(a,b,mid))>0.0) lo=mid;
    else hi=mid;
  }
  fraction=hi;
  return true;
}

inline SegmentResult ClassifySegment(const BoundaryPolicy& policy,
                                     const std::array<double,3>& a,
                                     const std::array<double,3>& b) {
  ValidateBoundaryPolicy(policy);
  SegmentResult out;
  if (!IsFinitePoint(a) || !IsFinitePoint(b)) {
    out.disposition=SegmentDisposition::Invalid;
    out.surface="NONFINITE";
    return out;
  }

  double boxFraction=2.0;
  std::string boxFace="NONE";
  const bool crossesBox=FirstBoxExitFraction(
      policy.computationalBox,a,b,boxFraction,boxFace);

  if (policy.kind==OuterBoundaryKind::Box) {
    if (!crossesBox) return out;
    out.disposition=SegmentDisposition::PhysicalEscape;
    out.fraction=boxFraction;
    out.position_m=InterpolatePoint(a,b,boxFraction);
    out.surface=boxFace;
    return out;
  }

  double physicalFraction=2.0;
  std::string physicalSurface="NONE";
  double shueFraction=2.0;
  if (FirstShueExitFraction(policy.shue,a,b,shueFraction)) {
    physicalFraction=shueFraction;
    physicalSurface="SHUE";
  }
  // The negative-X computational face is also the declared physical tail cap.
  if (crossesBox && boxFace=="XMIN" && boxFraction<physicalFraction) {
    physicalFraction=boxFraction;
    physicalSurface="XMIN_TAIL_CAP";
  }

  if (physicalFraction<=1.0 &&
      (!crossesBox || physicalFraction<=boxFraction+1.0e-13)) {
    out.disposition=SegmentDisposition::PhysicalEscape;
    out.fraction=physicalFraction;
    out.position_m=InterpolatePoint(a,b,physicalFraction);
    out.surface=physicalSurface;
    return out;
  }
  if (crossesBox) {
    out.disposition=SegmentDisposition::MeshUnavailable;
    out.fraction=boxFraction;
    out.position_m=InterpolatePoint(a,b,boxFraction);
    out.surface=boxFace;
    return out;
  }
  return out;
}

inline const char* OuterBoundaryKindName(OuterBoundaryKind kind) {
  return kind==OuterBoundaryKind::Shue ? "SHUE" : "BOX";
}

//--------------------------------------------------------------------------------------
// Machine-readable coupled cutoff/access artifact manifest
//--------------------------------------------------------------------------------------

inline std::string JsonEscape(const std::string& text) {
  std::ostringstream out;
  for (std::string::const_iterator it=text.begin();it!=text.end();++it) {
    switch (*it) {
      case '\\': out << "\\\\"; break;
      case '"': out << "\\\""; break;
      case '\n': out << "\\n"; break;
      case '\r': out << "\\r"; break;
      case '\t': out << "\\t"; break;
      default: out << *it; break;
    }
  }
  return out.str();
}

inline std::string BuildAccessManifestJson(
    const std::string& result,const std::string& snapshotId,
    const std::string& contentFingerprint,const std::string& meshRevision,
    const std::string& epochUTC,double simulationTime_s,
    const std::string& outputMode,const std::string& boundaryPolicy,
    const std::string& outputSuffix,const std::vector<std::string>& artifacts,
    const std::string& message) {
  if (result!="PASS" && result!="FAILED")
    throw std::invalid_argument("coupled access manifest result must be PASS or FAILED");
  if (snapshotId.empty() || contentFingerprint.empty() || meshRevision.empty() ||
      epochUTC.empty())
    throw std::invalid_argument(
        "coupled access manifest requires complete snapshot provenance");
  if (!std::isfinite(simulationTime_s) || simulationTime_s<0.0)
    throw std::invalid_argument(
        "coupled access manifest requires finite non-negative simulation time");
  const std::string canonicalOutputMode=UpperTrimmed(outputMode);
  if (canonicalOutputMode!="POINTS" && canonicalOutputMode!="TRAJECTORY" &&
      canonicalOutputMode!="SHELLS")
    throw std::invalid_argument(
        "coupled access manifest output mode must be POINTS, TRAJECTORY, or SHELLS");
  const std::string canonicalBoundary=UpperTrimmed(boundaryPolicy);
  if (canonicalBoundary!="BOX" && canonicalBoundary!="SHUE")
    throw std::invalid_argument(
        "coupled access manifest boundary policy must be BOX or SHUE");
  if (outputSuffix!=BuildProductSuffix(simulationTime_s,snapshotId))
    throw std::invalid_argument(
        "coupled access manifest suffix does not match its time/snapshot identity");
  if (result=="PASS" && artifacts.empty())
    throw std::invalid_argument(
        "PASS coupled access manifest requires at least one closed artifact");
  for (std::size_t n=0;n<artifacts.size();++n) {
    if (artifacts[n].empty())
      throw std::invalid_argument("coupled access manifest contains an empty artifact");
    if (std::find(artifacts.begin(),artifacts.begin()+n,artifacts[n])!=
        artifacts.begin()+n)
      throw std::invalid_argument("coupled access manifest contains duplicate artifacts");
  }
  std::ostringstream out;
  out.imbue(std::locale::classic());
  out << "{\n"
      << "  \"schema\": \"sep-in-geospace/swmf-coupled-access/v1\",\n"
      << "  \"RESULT\": \"" << result << "\",\n"
      << "  \"phase_1_interpretation\": \"INSTANTANEOUS_QUASI_STATIC\",\n"
      << "  \"snapshot_id\": \"" << JsonEscape(snapshotId) << "\",\n"
      << "  \"content_fingerprint\": \"" << JsonEscape(contentFingerprint) << "\",\n"
      << "  \"mesh_revision\": \"" << JsonEscape(meshRevision) << "\",\n"
      << "  \"epoch_utc\": \"" << JsonEscape(epochUTC) << "\",\n"
      << "  \"simulation_time_s\": " << std::setprecision(17)
      << simulationTime_s << ",\n"
      << "  \"output_mode\": \"" << JsonEscape(canonicalOutputMode) << "\",\n"
      << "  \"outer_boundary_policy\": \"" << JsonEscape(canonicalBoundary) << "\",\n"
      << "  \"output_suffix\": \"" << JsonEscape(outputSuffix) << "\",\n"
      << "  \"artifacts\": [";
  for (std::size_t n=0;n<artifacts.size();++n) {
    if (n!=0) out << ", ";
    out << "\"" << JsonEscape(artifacts[n]) << "\"";
  }
  out << "],\n"
      << "  \"message\": \"" << JsonEscape(message) << "\"\n"
      << "}\n";
  return out.str();
}

} // namespace SWMFCoupledAccess
} // namespace Earth

#endif // _SRC_EARTH_UTIL_SWMF_COUPLED_ACCESS_CONTRACT_H_
