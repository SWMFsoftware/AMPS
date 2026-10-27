#ifndef _SRC_GEOSPACE_GPU_ENGINE_FIELD_H_
#define _SRC_GEOSPACE_GPU_ENGINE_FIELD_H_

//======================================================================================
// EngineField.h -- magnetic-field functors used by the batched trajectory engine.
//
// A field functor provides
//   void GetB(const Vec3<R>& x,Vec3<R>& B)   x and B in internal units
//   void ResetAvailability() / bool Available()
//
// Availability reproduces the Mode3D sticky "field sample unavailable" status of
// cMode3DMeshFieldEval (AMR coverage hole or incomplete interpolation row).  Gridless
// fields are always available.
//
// DeviceField implements kFieldNone, kFieldDipole and kFieldMesh and is device code.
// HostField (EngineHostField.h) wraps a host evaluator for Geopack/Tsyganenko models.
//======================================================================================

#include "EngineMath.h"
#include "EngineMesh.h"

#include <type_traits>

namespace Geospace {
namespace GPU {

// Constants of gridless/DipoleInterface.h (Geospace::GridlessMode::Dipole).
struct EgDipoleConstants {
  static constexpr double kPI=3.14159265358979323846;
  static constexpr double kReKm=6371.2;
  static constexpr double kRe_m=kReKm*1000.0;
  static constexpr double kMu0Over4Pi=1e-7;
  static constexpr double kBeqRe=3.12e-5;
  static constexpr double kM_E_Am2=kBeqRe*kRe_m*kRe_m*kRe_m/kMu0Over4Pi;
};

// Dipole moment components prepared on the host exactly as the legacy evaluator does
// (gridless: M*m_hat from Dipole::MakeParams; Mode3D: moment*sin/cos(tilt)).
struct EgDipole {
  double mx,my,mz;   // SI [A m^2]
};

// Accumulator of the Mode3D DIPOLE mesh-interpolation error diagnostic.
struct EgDipoleErrorStats {
  unsigned long long count;
  double sum;
  double maxValue;
  double maxLocation[3];
};

EG_HD void ResetDipoleErrorStats(EgDipoleErrorStats& s) {
  s.count=0; s.sum=0.0; s.maxValue=-1.0;
  s.maxLocation[0]=s.maxLocation[1]=s.maxLocation[2]=0.0;
}

// SI dipole formula of Dipole::GetB_Tesla / EvaluateConfiguredDipoleField_.
// Returns false for the cases the Mode3D strict evaluator rejects.
EG_HD bool DipoleFieldSI(const EgDipole& d,const double x_m[3],double B_T[3],bool strict) {
  const double x=x_m[0],y=x_m[1],z=x_m[2];
  const double r2=x*x+y*y+z*z;
  if (strict) {
    if (!(r2>0.0) || !std::isfinite(r2)) return false;
  }
  else if (r2<=0.0) { B_T[0]=B_T[1]=B_T[2]=0.0; return true; }

  const double r=std::sqrt(r2);
  const double r3=r2*r;
  const double r5=r3*r2;
  const double mdotr=d.mx*x+d.my*y+d.mz*z;
  const double c=EgDipoleConstants::kMu0Over4Pi;
  B_T[0]=c*(3.0*x*mdotr/r5-d.mx/r3);
  B_T[1]=c*(3.0*y*mdotr/r5-d.my/r3);
  B_T[2]=c*(3.0*z*mdotr/r5-d.mz/r3);
  if (strict) return std::isfinite(B_T[0]) && std::isfinite(B_T[1]) && std::isfinite(B_T[2]);
  return true;
}

template <class R,class U>
class DeviceField {
 public:
  EG_HD DeviceField(const EgContext& cx,const EgMeshView<R>* mesh,EgRowStencil<R>* row,
                    EgRowStencil<R>* fineRow,int lastNode)
      : kind_(cx.fieldKind), strict_(cx.dipoleStrict!=0),
        sampleError_(cx.sampleDipoleError!=0), mesh_(mesh), row_(row), fineRow_(fineRow),
        lastNode_(lastNode), available_(true), error_(kErrNone) {
    dip_.mx=cx.dipoleMx; dip_.my=cx.dipoleMy; dip_.mz=cx.dipoleMz;
    if constexpr (!U::kSI) {
      // B_int = (mu0/4pi) M/(B_unit L_unit^3) (3 x (m.x)/r^5 - m/r^3), x in L_unit.
      const double s=EgDipoleConstants::kMu0Over4Pi/
          (U::kField_T*U::kLength_m*U::kLength_m*U::kLength_m);
      cm_[0]=R(s*dip_.mx); cm_[1]=R(s*dip_.my); cm_[2]=R(s*dip_.mz);
    }
    for (int d=0;d<3;d++) {
      domMin_[d]=R(ToLength<U>(cx.fieldDomainMin_m[d]));
      domMax_[d]=R(ToLength<U>(cx.fieldDomainMax_m[d]));
      if constexpr (!std::is_same<R,double>::value) {
        // Float: a sample within float resolution of a domain face can round onto the
        // face, where the AMR search already reports "outside the mesh" while the
        // configured-box test would still say "inside" and flag an AMR coverage hole
        // (INVALID_FIELD).  The legacy tracer (double) has this ambiguity only for samples
        // exactly on the face; in float the band is a few ulps wide, so it is treated
        // as outside.  The double build keeps the exact legacy test.
        const R tolMin=R(8)*std::numeric_limits<R>::epsilon()*std::fabs(domMin_[d]);
        const R tolMax=R(8)*std::numeric_limits<R>::epsilon()*std::fabs(domMax_[d]);
        domMin_[d]+=tolMin;
        domMax_[d]-=tolMax;
      }
    }
    ResetDipoleErrorStats(stats_);
  }

  EG_HD void ResetAvailability() { available_=true; }
  EG_HD bool Available() const { return available_; }
  EG_HD int Error() const { return error_; }
  EG_HD int LastNode() const { return lastNode_; }
  EG_HD const EgDipoleErrorStats& Stats() const { return stats_; }
  // First field sample that was unavailable (diagnostics): kind 1 = no AMR leaf inside
  // the configured domain, 2 = empty interpolation row.
  EG_HD int UnavailableKind() const { return unavailableKind_; }
  EG_HD double UnavailableX(int d) const { return FromLength<U>(double(unavailableX_[d])); }

  EG_HD void GetB(const Vec3<R>& x,Vec3<R>& B) {
    if (kind_==kFieldNone) {
      B.x=R(0); B.y=R(0); B.z=R(0);
      return;
    }
    if (kind_==kFieldDipole) {
      if constexpr (U::kSI) {
        double xs[3]={x.x,x.y,x.z},bs[3];
        if (!DipoleFieldSI(dip_,xs,bs,strict_)) {
          error_=kErrDipoleField;
          B.x=B.y=B.z=R(0);
          return;
        }
        B.x=bs[0]; B.y=bs[1]; B.z=bs[2];
      }
      else {
        // Normalized units: same formula with the moment pre-scaled to internal units,
        // evaluated in R so that a float build does no double arithmetic here.
        const R r2=x.x*x.x+x.y*x.y+x.z*x.z;
        if (strict_ && (!(r2>R(0)) || !efinite(r2))) {
          error_=kErrDipoleField;
          B.x=B.y=B.z=R(0);
          return;
        }
        if (!strict_ && r2<=R(0)) { B.x=B.y=B.z=R(0); return; }
        const R r=std::sqrt(r2);
        const R r3=r2*r;
        const R r5=r3*r2;
        const R mdotr=cm_[0]*x.x+cm_[1]*x.y+cm_[2]*x.z;
        B.x=R(3)*x.x*mdotr/r5-cm_[0]/r3;
        B.y=R(3)*x.y*mdotr/r5-cm_[1]/r3;
        B.z=R(3)*x.z*mdotr/r5-cm_[2]/r3;
      }
      return;
    }

    // kFieldMesh: cMode3DMeshFieldEval::GetB_T production path.
    const R xa[3]={x.x,x.y,x.z};
    EgMesh<R> mesh(*mesh_);
    const int node=mesh.Find(xa,lastNode_);
    lastNode_=node;
    if (node<0) {
      const bool insideConfiguredBox=
          xa[0]>=domMin_[0] && xa[0]<=domMax_[0] &&
          xa[1]>=domMin_[1] && xa[1]<=domMax_[1] &&
          xa[2]>=domMin_[2] && xa[2]<=domMax_[2];
      if (insideConfiguredBox) { available_=false; NoteUnavailable(xa,1); }
      B.x=B.y=B.z=R(0);
      return;
    }
    R b[3]={R(0),R(0),R(0)};
    int err=kErrNone;
    const bool ok=mesh.Interpolate(xa,node,b,*row_,*fineRow_,err);
    if (err!=kErrNone) error_=err;
    if (!ok) {
      available_=false;
      NoteUnavailable(xa,2);
      B.x=B.y=B.z=R(0);
      return;
    }
    B.x=b[0]; B.y=b[1]; B.z=b[2];
    if (sampleError_) RecordDipoleError(xa,b);
  }

 private:
  // RecordDipoleFieldError_() + cDipoleFieldErrorLocalStats_::Add()
  EG_HD void RecordDipoleError(const R xa[3],const R b[3]) {
    double xs[3],bm[3],bref[3]={0.0,0.0,0.0};
    for (int d=0;d<3;d++) { xs[d]=FromLength<U>(xa[d]); bm[d]=FromField<U>(b[d]); }
    if (!DipoleFieldSI(dip_,xs,bref,true)) return;
    const double referenceNorm2=bref[0]*bref[0]+bref[1]*bref[1]+bref[2]*bref[2];
    if (!(referenceNorm2>std::numeric_limits<double>::min())) return;
    const double dBx=bm[0]-bref[0];
    const double dBy=bm[1]-bref[1];
    const double dBz=bm[2]-bref[2];
    const double differenceNorm2=dBx*dBx+dBy*dBy+dBz*dBz;
    const double relativeError=std::sqrt(differenceNorm2/referenceNorm2);
    if (!std::isfinite(relativeError) || relativeError<0.0) return;
    ++stats_.count;
    stats_.sum+=relativeError;
    if (relativeError>stats_.maxValue) {
      stats_.maxValue=relativeError;
      stats_.maxLocation[0]=xs[0];
      stats_.maxLocation[1]=xs[1];
      stats_.maxLocation[2]=xs[2];
    }
  }

  int kind_;
  bool strict_;
  bool sampleError_;
  EgDipole dip_;
  R cm_[3];
  R domMin_[3],domMax_[3];
  const EgMeshView<R>* mesh_;
  EgRowStencil<R>* row_;
  EgRowStencil<R>* fineRow_;
  int lastNode_;
  bool available_;
  int error_;
  int unavailableKind_=0;
  R unavailableX_[3]={R(0),R(0),R(0)};
  EG_HD void NoteUnavailable(const R xa[3],int kind) {
    if (unavailableKind_!=0) return;
    unavailableKind_=kind;
    unavailableX_[0]=xa[0]; unavailableX_[1]=xa[1]; unavailableX_[2]=xa[2];
  }
  EgDipoleErrorStats stats_;
};

} // namespace GPU
} // namespace Geospace

#endif
