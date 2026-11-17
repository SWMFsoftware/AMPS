#include "sep_polarized_wkb_transport.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <exception>
#include <limits>
#include <set>
#include <sstream>

namespace SEP { namespace Turbulence { namespace Polarized { namespace WKB { namespace {
using Transport::Status;using Transport::StatusCode;
const double LightSpeedMPerS=299792458.0;
Status Error(const std::string& message){return Status::Error(StatusCode::InvalidArgument,message);}
Status Unresolved(const std::string& message){return Status::Error(StatusCode::UnresolvedCoefficient,message);}
bool Positive(double x){return std::isfinite(x)&&x>0;}
long double Dot(const Vector& a,const Vector& b) {
  long double value=0;for(int d=0;d<3;++d)value+=static_cast<long double>(a[d])*b[d];return value;
}
long double Norm(const Vector& a){return std::sqrt(Dot(a,a));}
bool Finite(const Vector& a){for(double x:a)if(!std::isfinite(x))return false;return true;}
Vector Unit(const Vector& a) {
  const auto n=Norm(a);Vector b;for(int d=0;d<3;++d)b[d]=static_cast<double>(a[d]/n);return b;
}
// These bounds describe roundoff of the explicitly shown arithmetic, not an
// accepted physical-domain error. Each scale retains the residual's own units.
bool Roundoff(long double residual,long double scale) {
  return std::isfinite(residual)&&std::isfinite(scale)&&
      std::fabs(residual)<=256*std::numeric_limits<double>::epsilon()*scale;
}
Status Represent(long double value,double* output,bool nonnegative=false) {
  const double rounded=static_cast<double>(value);
  if(!std::isfinite(value)||!std::isfinite(rounded)||(value!=0&&rounded==0)||
      (nonnegative&&value<0))return Unresolved("WKB energy/action is outside resolved finite floating-point range");
  *output=rounded;return Status::Ok();
}
Status GridValid(const Grid& grid) {
  CenteredLogGrid centered;centered.centersPerM=grid.centersPerM;centered.logWidth=grid.logWidth;
  const auto resolved=ResolveCenteredLogGrid(centered);
  if(!resolved.status.ok())return resolved.status;
  return resolved.grid.edgesPerM==grid.edgesPerM?Status::Ok():Error("WKB requires the exact declared fixed k edges");
}
Status BinsValid(const FourBins& values,std::size_t bins) {
  for(const auto& component:values) {
    if(component.size()!=bins)return Error("WKB four-component shape differs from its fixed grid");
    for(double e:component)if(!std::isfinite(e)||e<0)return Error("WKB magnetic energy/density is negative or nonfinite");
  }
  return Status::Ok();
}
long double ActionFactor(const Grid& grid,std::size_t j,long double va) {
  const long double low=grid.edgesPerM[j],width=static_cast<long double>(grid.edgesPerM[j+1])-low;
  // For the SAME flat-per-linear-k reconstruction consumed by scattering:
  // A_bin=integral 2*g_B/(k*v_A) dk = 2*Ebin*ln(khi/klo)/(v_A*Delta k).
  return 2*std::log1p(width/low)/(va*width);
}
long double Sum(const std::vector<double>& a){long double sum=0;for(double x:a)sum+=x;return sum;}
long double Action(const Grid& grid,const std::vector<double>& a,long double va) {
  long double value=0;for(std::size_t j=0;j<a.size();++j)value+=a[j]*ActionFactor(grid,j,va);return value;
}
std::string ReservoirIdentity(const InitialReservoir& reservoir) {
  // A boundary reservoir is a fixed external value object, not a Snapshot.
  // This versioned canonical byte encoding detects accidental wire mutation;
  // it neither duplicates SnapshotIdentity nor authenticates historical origin.
  // Numeric fields are already finite/shape-validated before this helper runs.
  std::uint64_t hash=UINT64_C(1469598103934665603);
  auto byte=[&](unsigned char value){hash^=value;hash*=UINT64_C(1099511628211);};
  auto integer=[&](std::uint64_t value){for(int i=0;i<8;++i){byte(static_cast<unsigned char>(value&255));value>>=8;}};
  auto number=[&](double value){std::uint64_t bits=0;std::memcpy(&bits,&value,sizeof(bits));integer(bits);};
  auto text=[&](const std::string& value){integer(value.size());for(unsigned char c:value)byte(c);};
  static_assert(sizeof(double)==sizeof(std::uint64_t)&&std::numeric_limits<double>::is_iec559,
      "WKB reservoir identity requires IEEE binary64");
  integer(reservoir.fieldLine);integer(static_cast<std::uint64_t>(reservoir.side));number(reservoir.initialEpochS);
  text(reservoir.initialSnapshotIdentity);text(reservoir.configurationFingerprint);text(reservoir.provenance);
  number(reservoir.grid.logWidth);integer(reservoir.grid.centersPerM.size());
  for(double value:reservoir.grid.centersPerM)number(value);
  integer(reservoir.grid.edgesPerM.size());for(double value:reservoir.grid.edgesPerM)number(value);
  for(const auto& component:reservoir.magneticDensityJPerM3) {
    integer(component.size());for(double value:component)number(value);
  }
  std::ostringstream result;result<<"initial-wave-reservoir-v1:"<<std::hex<<hash;return result.str();
}
Status ReservoirValid(const InitialReservoir& reservoir) {
  auto status=GridValid(reservoir.grid);if(!status.ok())return status;
  if(reservoir.initialSnapshotIdentity.empty()||reservoir.configurationFingerprint.empty()||reservoir.provenance.empty()||!std::isfinite(reservoir.initialEpochS)||
      (reservoir.side!=Side::Lower&&reservoir.side!=Side::Upper))return Error("initial boundary reservoir lacks explicit provenance/side/epoch");
  status=BinsValid(reservoir.magneticDensityJPerM3,reservoir.grid.centersPerM.size());if(!status.ok())return status;
  return reservoir.identity==ReservoirIdentity(reservoir)?Status::Ok():Error("fixed initial wave reservoir content identity differs from its captured values");
}
Status FaceValid(const FaceInterval& face,double mu) {
  if(!face.integralsKnown||!face.noCharacteristicReversalKnown)
    return Status::Error(StatusCode::UnsupportedConfiguration,"WKB face needs declared interval integrals and no unresolved characteristic reversal");
  if(!Positive(mu)||!Positive(face.initialAreaM2)||!Positive(face.finalAreaM2)||
      !Positive(face.initialMassDensityKgM3)||!Positive(face.finalMassDensityKgM3))return Error("invalid WKB face area/density");
  const Vector* fields[2]={&face.initialMagneticFieldT,&face.finalMagneticFieldT};
  const Vector* velocity[2]={&face.initialPlasmaVelocityMPerS,&face.finalPlasmaVelocityMPerS};
  const Vector* grid[2]={&face.initialGridVelocityMPerS,&face.finalGridVelocityMPerS};
  const Vector* normal[2]={&face.initialNormal,&face.finalNormal};
  const double density[2]={face.initialMassDensityKgM3,face.finalMassDensityKgM3};
  for(int end=0;end<2;++end) {
    if(!Finite(*fields[end])||!Finite(*velocity[end])||!Finite(*grid[end])||!(Norm(*fields[end])>0)||Norm(*velocity[end])>=LightSpeedMPerS)
      return Error("WKB face has invalid B/U/grid primitives");
    if(!Finite(*normal[end])||!Roundoff(Norm(*normal[end])-1,1))return Error("invalid WKB end-cap unit normal");
    const auto b=Unit(*fields[end]);const auto polarity=Dot(b,*normal[end]);
    if(polarity==0)return Error("WKB tube end cap carries no magnetic flux");
    const auto va=Norm(*fields[end])/std::sqrt(static_cast<long double>(mu)*density[end]);
    if(!std::isfinite(va)||va>=LightSpeedMPerS)return Error("WKB face Alfven speed is outside the physical frame");
    for(int branch=0;branch<2;++branch) {
      const int sigma=branch==0?1:-1;
      const long double q=face.plasmaSweptVolumeM3-face.gridSweptVolumeM3+sigma*face.alfvenSweptVolumeM3;
      const long double a=face.relativeFlowOverAlfvenAreaTimeM2S+sigma*face.signedAreaTimeM2S;
      Vector relative;for(int d=0;d<3;++d)relative[d]=(*velocity[end])[d]-(*grid[end])[d];
      const auto speed=Dot(relative,*normal[end])+sigma*va*polarity;
      if(!std::isfinite(q)||!std::isfinite(a)||q*speed<0||q*a<0||(q==0&&a!=0))
        return Error("WKB face integrated energy/action characteristic disagrees with its declared sign");
    }
  }
  for(double value:{face.plasmaSweptVolumeM3,face.gridSweptVolumeM3,face.alfvenSweptVolumeM3,
      face.relativeMassTransportKg,face.relativeFlowOverAlfvenAreaTimeM2S,face.signedAreaTimeM2S})
    if(!std::isfinite(value))return Error("nonfinite WKB face transport integral");
  return Status::Ok();
}
// Stable analytic integrals of the crossing energy/action, including the exact
// r=-1 energy and r=0 action limits. No epsilon k floor or almost-equal edge
// rule is inserted. std::expm1 keeps the narrow-bin integral well conditioned.
long double Exprel(long double x){return x==0?1:std::expm1(x)/x;}
long double CrossingIntegral(long double low,long double high,long double edge,long double power,bool inverseK) {
  if(high<=low)return 0;
  const auto span=std::log1p((high-low)/low);
  return (inverseK?1:low)*std::exp(power*std::log1p((low-edge)/edge))*span*
      Exprel((power+(inverseK?0:1))*span);
}
Status AccountIncrement(double increment,Transport::ExactNonnegativeSum* exact,
    double* nearest,double* remainder,double* bound,bool validateOnly=false) {
  Transport::ExactSumApproximation display;
  auto status=Transport::ApproximateExactNonnegative(*exact,&display);if(!status.ok())return status;
  if(validateOnly)return *nearest==display.nearest&&*remainder==display.signedRemainder&&*bound==display.absoluteRoundingBound?
      Status::Ok():Error("external spectral diagnostics do not match their exact integer authority");
  Transport::ExactNonnegativeSum candidate;
  status=Transport::AddExactNonnegative(*exact,increment,&candidate);if(!status.ok())return status;
  status=Transport::ApproximateExactNonnegative(candidate,&display);if(!status.ok())return status;
  *exact=candidate;*nearest=display.nearest;*remainder=display.signedRemainder;*bound=display.absoluteRoundingBound;
  return Status::Ok();
}
} // namespace

ReservoirResult CaptureInitialReservoir(const Snapshot& initial,std::uint64_t line,Side side) {
  ReservoirResult result;result.status=ValidateSnapshot(initial);if(!result.status.ok())return result;
  if(initial.revision!=1) {result.status=Error("incoming reservoir must be captured from declared initialization, not an evolved state");return result;}
  const Cell* selected=nullptr;
  for(const auto& cell:initial.cells)if(cell.fieldLine==line&&(!selected||side==Side::Upper))selected=&cell;
  if(!selected){result.status=Error("initial boundary field-line namespace is absent");return result;}
  auto& reservoir=result.reservoir;reservoir.grid=initial.grid;reservoir.fieldLine=line;reservoir.side=side;
  reservoir.initialEpochS=initial.epochS;reservoir.initialSnapshotIdentity=initial.identity;
  reservoir.configurationFingerprint=initial.configurationFingerprint;
  reservoir.provenance="fixed declared INITIAL magnetic J/m^3 per bin; incoming characteristics only";
  for(int c=0;c<4;++c)for(double energy:selected->magneticEnergyJ[c]) {
    double density;result.status=Represent(static_cast<long double>(energy)/selected->volumeM3,&density,true);
    if(!result.status.ok())return result;
    reservoir.magneticDensityJPerM3[c].push_back(density);
  }
  reservoir.identity=ReservoirIdentity(reservoir);result.status=ReservoirValid(reservoir);return result;
}

BoundaryFluxResult EvaluateBoundaryFlux(const InitialReservoir& reservoir,const FourBins& interior,
    const FaceInterval& face,double mu) {
  BoundaryFluxResult result;result.status=ReservoirValid(reservoir);if(!result.status.ok())return result;
  result.status=FaceValid(face,mu);if(!result.status.ok())return result;
  const auto bins=reservoir.grid.centersPerM.size();result.status=BinsValid(interior,bins);if(!result.status.ok())return result;
  const int normal=reservoir.side==Side::Lower?-1:1;
  for(int c=0;c<4;++c) {
    const int branch=c/2,sigma=branch==0?1:-1;
    const long double q=normal*(face.plasmaSweptVolumeM3-face.gridSweptVolumeM3+sigma*face.alfvenSweptVolumeM3);
    const long double actionQ=normal*(face.relativeFlowOverAlfvenAreaTimeM2S+sigma*face.signedAreaTimeM2S);
    result.flux.incoming[branch]=q<0;
    for(std::size_t j=0;j<bins;++j) {
      const double density=q<0?reservoir.magneticDensityJPerM3[c][j]:interior[c][j];
      const long double low=reservoir.grid.edgesPerM[j],width=static_cast<long double>(reservoir.grid.edgesPerM[j+1])-low;
      double energy,action,pressure;
      result.status=Represent(q*density,&energy);if(!result.status.ok())return result;
      result.status=Represent(actionQ*2*std::log1p(width/low)/width*density,&action);if(!result.status.ok())return result;
      result.status=Represent(normal*static_cast<long double>(face.plasmaSweptVolumeM3)*density,&pressure);if(!result.status.ok())return result;
      result.flux.outwardMagneticJ[c].push_back(energy);result.flux.outwardActionJs[c].push_back(action);
      result.flux.outwardPressureWorkJ[c].push_back(pressure);
    }
  }
  result.status=Status::Ok();return result;
}

ExternalAccountResult BeginExternalAccount(const Snapshot& initial) {
  ExternalAccountResult result;result.status=ValidateSnapshot(initial);if(!result.status.ok())return result;
  if(initial.revision!=1){result.status=Error("external spectral account must start at declared initialization");return result;}
  result.account.originSnapshotIdentity=initial.identity;result.account.lastSnapshotIdentity=initial.identity;return result;
}

SpectralCellResult PrepareSpectralWork(const Grid& grid,const FourBins& energy,double va0,double va1,
    double theta,double strain,double gradVA) {
  SpectralCellResult result;result.status=GridValid(grid);if(!result.status.ok())return result;
  result.status=BinsValid(energy,grid.centersPerM.size());if(!result.status.ok())return result;
  if(!Positive(va0)||!Positive(va1)||va0>=LightSpeedMPerS||va1>=LightSpeedMPerS||
      !std::isfinite(theta)||!std::isfinite(strain)||!std::isfinite(gradVA)) {
    result.status=Error("spectral WKB source lacks finite physical speeds/integrated derivatives");return result;
  }
  const std::size_t bins=grid.centersPerM.size();
  const long double workLog=-0.5L*theta,work=std::exp(workLog),q=std::log(static_cast<long double>(va1)/va0);
  if(!std::isfinite(work)||!(work>0)){result.status=Unresolved("WKB compression factor is unresolved");return result;}
  for(int c=0;c<4;++c) {
    const int sigma=c/2==0?1:-1;
    const long double shift=-static_cast<long double>(strain)-sigma*static_cast<long double>(gradVA);
    const auto inverseShift=std::exp(-shift);
    if(!std::isfinite(inverseShift)||!(inverseShift>0)){result.status=Unresolved("WKB ray shift exceeds numerical range");return result;}
    if(shift!=0&&inverseShift==1&&Sum(energy[c])>0){result.status=Unresolved("nonzero WKB ray shift is below edge resolution");return result;}
    std::vector<long double> finalEnergy(bins,0);
    for(std::size_t j=0;j<bins;++j) {
      const long double low=grid.edgesPerM[j]*inverseShift,high=grid.edgesPerM[j+1]*inverseShift;
      if(!std::isfinite(low)||!std::isfinite(high)){result.status=Unresolved("WKB ray preimage exceeds numerical range");return result;}
      auto found=std::upper_bound(grid.edgesPerM.begin(),grid.edgesPerM.end(),low);
      std::size_t old=found==grid.edgesPerM.begin()?0:static_cast<std::size_t>(found-grid.edgesPerM.begin()-1);
      for(;old<bins&&grid.edgesPerM[old]<high;++old) {
        const auto left=std::max(low,static_cast<long double>(grid.edgesPerM[old]));
        const auto right=std::min(high,static_cast<long double>(grid.edgesPerM[old+1]));
        if(right>left)finalEnergy[j]+=work*energy[c][old]*(right-left)/
            (static_cast<long double>(grid.edgesPerM[old+1])-grid.edgesPerM[old]);
      }
    }
    long double escapedB[2]={0,0},escapedA[2]={0,0};
    if(shift!=0) {
      const int side=shift<0?0:1;const long double edge=side==0?grid.edgesPerM.front():grid.edgesPerM.back();
      const auto threshold=edge*inverseShift;
      for(std::size_t old=0;old<bins;++old) {
        const long double low=side==0?grid.edgesPerM[old]:std::max(static_cast<long double>(grid.edgesPerM[old]),threshold);
        const long double high=side==0?std::min(static_cast<long double>(grid.edgesPerM[old+1]),threshold):grid.edgesPerM[old+1];
        if(high<=low)continue;
        const long double density=energy[c][old]/(static_cast<long double>(grid.edgesPerM[old+1])-grid.edgesPerM[old]);
        // Crossing fraction f=ln(kedge/k0)/shift. Work and VA are evaluated
        // at f, so escaped waves receive no fictitious later on-grid work.
        escapedB[side]+=density*CrossingIntegral(low,high,edge,-workLog/shift,false);
        escapedA[side]+=2*density/va0*CrossingIntegral(low,high,edge,(shift+q-workLog)/shift,true);
      }
    }
    long double finalSum=0;
    for(long double value:finalEnergy) {
      double rounded;result.status=Represent(value,&rounded,true);if(!result.status.ok())return result;
      result.magneticEnergyAfterJ[c].push_back(rounded);finalSum+=value;
    }
    auto& ledger=result.ledger[c];const auto initial=Sum(energy[c]);
    const long double physicalWork=theta==0?0:finalSum+escapedB[0]+escapedB[1]-initial;
    if(initial>0&&theta!=0&&(physicalWork==0||physicalWork*theta>=0)) {
      result.status=Unresolved("WKB compression work is below resolved signed energy precision");return result;
    }
    double* fields[]={&ledger.magneticBeforeJ,&ledger.magneticAfterJ,&ledger.lowerSpectralOutgoingMagneticJ,
        &ledger.upperSpectralOutgoingMagneticJ,&ledger.magneticCompressionWorkJ,&ledger.actionBeforeJs,
        &ledger.actionAfterJs,&ledger.lowerSpectralOutgoingActionJs,&ledger.upperSpectralOutgoingActionJs};
    const long double values[]={initial,Sum(result.magneticEnergyAfterJ[c]),escapedB[0],escapedB[1],physicalWork,
        Action(grid,energy[c],va0),Action(grid,result.magneticEnergyAfterJ[c],va1),escapedA[0],escapedA[1]};
    for(std::size_t i=0;i<9;++i){result.status=Represent(values[i],fields[i]);if(!result.status.ok())return result;}
    ledger.totalBackgroundCompressionReceiptJ=-2*ledger.magneticCompressionWorkJ;
    ledger.totalBackgroundMechanicalReceiptJ=ledger.totalBackgroundCompressionReceiptJ;
    ledger.energyRoundoffResidualJ=static_cast<double>(values[1]+values[2]+values[3]-values[0]-values[4]);
    ledger.actionDiscretizationResidualJs=static_cast<double>(values[6]+values[7]+values[8]-values[5]);
  }
  result.status=Status::Ok();return result;
}

// PrepareInterval is defined below the separately testable boundary/source
// stages. Keeping publication outside this module prevents partial native
// writes if a later cell, budget, geometry or characteristic check fails.

StageResult PrepareInterval(const Snapshot& before,const Interval& interval) {
  StageResult result;result.status=ValidateSnapshot(before);if(!result.status.ok())return result;
  if(!Positive(interval.durationS)||interval.reflectionEnabled||interval.cascadeEnabled||!interval.zeroExteriorSpectralIncoming||
      (interval.geometryMode!=GeometryMode::MaterialTube&&interval.geometryMode!=GeometryMode::ControlledExplicitAle)) {
    result.status=Status::Error(StatusCode::UnsupportedConfiguration,
        "controlled WKB requires explicit grid mode/positive interval, zero k-exterior inflow, reflection OFF and cascade OFF");return result;
  }
  if(!Roundoff(static_cast<long double>(interval.target.epochS)-before.epochS-interval.durationS,
      std::fabs(static_cast<long double>(interval.target.epochS))+std::fabs(before.epochS)+interval.durationS)||
      interval.externalBefore.originSnapshotIdentity.empty()||interval.externalBefore.lastSnapshotIdentity!=before.identity) {
    result.status=Error("WKB target time/external account does not bind the immutable input realization");return result;
  }
  for(int c=0;c<4;++c) {
    for(double value:{interval.externalBefore.lowerMagneticJ[c],interval.externalBefore.upperMagneticJ[c],
        interval.externalBefore.lowerActionJs[c],interval.externalBefore.upperActionJs[c]})
      if(!std::isfinite(value)||value<0){result.status=Error("invalid external spectral account");return result;}
    for(double value:{interval.externalBefore.lowerMagneticRoundoffJ[c],interval.externalBefore.upperMagneticRoundoffJ[c],
        interval.externalBefore.lowerActionRoundoffJs[c],interval.externalBefore.upperActionRoundoffJs[c]})
      if(!std::isfinite(value)){result.status=Error("invalid external spectral rounding ledger");return result;}
  }
  try {
    result.externalAfter=interval.externalBefore;
    for(int c=0;c<4;++c) {
      auto& a=result.externalAfter;
      result.status=AccountIncrement(0,&a.lowerMagneticExact[c],&a.lowerMagneticJ[c],&a.lowerMagneticRoundoffJ[c],&a.lowerMagneticDiagnosticBoundJ[c],true);
      if(!result.status.ok())return result;
      result.status=AccountIncrement(0,&a.upperMagneticExact[c],&a.upperMagneticJ[c],&a.upperMagneticRoundoffJ[c],&a.upperMagneticDiagnosticBoundJ[c],true);
      if(!result.status.ok())return result;
      result.status=AccountIncrement(0,&a.lowerActionExact[c],&a.lowerActionJs[c],&a.lowerActionRoundoffJs[c],&a.lowerActionDiagnosticBoundJs[c],true);
      if(!result.status.ok())return result;
      result.status=AccountIncrement(0,&a.upperActionExact[c],&a.upperActionJs[c],&a.upperActionRoundoffJs[c],&a.upperActionDiagnosticBoundJs[c],true);
      if(!result.status.ok())return result;
    }
    auto updated=before.cells;std::set<std::uint64_t> coveredLines;std::size_t coveredCells=0;
    const auto bins=before.grid.centersPerM.size();const double mu=before.vacuumPermeabilityHPerM;
    for(const auto& tube:interval.tubes) {
      if(!coveredLines.insert(tube.fieldLine).second){result.status=Error("duplicate WKB field-line interval");return result;}
      std::vector<std::size_t> indices;
      for(std::size_t i=0;i<before.cells.size();++i)if(before.cells[i].fieldLine==tube.fieldLine)indices.push_back(i);
      const auto count=indices.size();
      if(count==0||tube.cells.size()!=count||tube.faces.size()!=count+1) {
        result.status=Error("WKB interval does not cover the complete physical tube");return result;
      }
      coveredCells+=count;
      if(tube.lower.fieldLine!=tube.fieldLine||tube.upper.fieldLine!=tube.fieldLine||
          tube.lower.side!=Side::Lower||tube.upper.side!=Side::Upper||
          tube.lower.grid.edgesPerM!=before.grid.edgesPerM||tube.upper.grid.edgesPerM!=before.grid.edgesPerM||
          tube.lower.configurationFingerprint!=before.configurationFingerprint||tube.upper.configurationFingerprint!=before.configurationFingerprint||
          tube.lower.initialSnapshotIdentity!=interval.externalBefore.originSnapshotIdentity||
          tube.upper.initialSnapshotIdentity!=interval.externalBefore.originSnapshotIdentity) {
        result.status=Error("WKB boundary reservoirs differ from the fixed declared initialization/grid");return result;
      }
      for(const auto& face:tube.faces) {
        result.status=FaceValid(face,mu);if(!result.status.ok())return result;
        if(interval.geometryMode==GeometryMode::MaterialTube) {
          for(int d=0;d<3;++d)if(!Roundoff(face.initialPlasmaVelocityMPerS[d]-face.initialGridVelocityMPerS[d],
              std::fabs(face.initialPlasmaVelocityMPerS[d])+std::fabs(face.initialGridVelocityMPerS[d]))||
              !Roundoff(face.finalPlasmaVelocityMPerS[d]-face.finalGridVelocityMPerS[d],
              std::fabs(face.finalPlasmaVelocityMPerS[d])+std::fabs(face.finalGridVelocityMPerS[d]))) {
            result.status=Error("material WKB vertices must move with the full plasma velocity");return result;
          }
          if(!Roundoff(face.plasmaSweptVolumeM3-face.gridSweptVolumeM3,
              std::fabs(face.plasmaSweptVolumeM3)+std::fabs(face.gridSweptVolumeM3))||
              face.relativeMassTransportKg!=0||face.relativeFlowOverAlfvenAreaTimeM2S!=0) {
            result.status=Error("material WKB cannot retain Eulerian U convection in relative face fluxes");return result;
          }
        }
      }
      // A fixed labeled magnetic flux is checked at both time endpoints and
      // every face. Areas are physical inputs; no conical/constant-area guess
      // or seed-field fallback occurs in this controlled operator.
      const long double flux=tube.faces.front().initialAreaM2*
          Dot(tube.faces.front().initialMagneticFieldT,tube.faces.front().initialNormal);
      for(const auto& face:tube.faces)for(int end=0;end<2;++end) {
        const auto measured=(end==0?face.initialAreaM2:face.finalAreaM2)*
            Dot(end==0?face.initialMagneticFieldT:face.finalMagneticFieldT,end==0?face.initialNormal:face.finalNormal);
        if(!Roundoff(measured-flux,std::fabs(measured)+std::fabs(flux))) {
          result.status=Error("WKB tube magnetic flux is not conserved by the supplied face geometry");return result;
        }
      }
      std::vector<FourBins> density(count);
      for(std::size_t i=0;i<count;++i) {
        const auto& cell=before.cells[indices[i]];const auto& input=tube.cells[i];
        if(cell.segment!=i||!input.inputsKnown||!input.fixedFieldDirectionKnown||!input.lateralTubeMaterialKnown||
            !Finite(input.initialMagneticFieldT)||!Finite(input.finalMagneticFieldT)||
            !(Norm(input.initialMagneticFieldT)>0)||!(Norm(input.finalMagneticFieldT)>0)||
            !Positive(input.finalVolumeM3)||!Positive(input.finalMassDensityKgM3)||
            !Positive(input.initialCellMassKg)||!Positive(input.finalCellMassKg)) {
          result.status=Status::Error(StatusCode::UnsupportedConfiguration,
              "WKB needs all contiguous cells, explicit ideal interval primitives and material side-wall/fixed-direction capabilities");return result;
        }
        const auto b=Unit(input.initialMagneticFieldT),endB=Unit(input.finalMagneticFieldT);
        Vector cross{{b[1]*endB[2]-b[2]*endB[1],b[2]*endB[0]-b[0]*endB[2],b[0]*endB[1]-b[1]*endB[0]}};
        if(!Roundoff(Norm(input.initialMagneticFieldT)-cell.meanMagneticFieldT,cell.meanMagneticFieldT)||
            !Roundoff(Norm(cross),1)||Dot(b,endB)<=0) {
          result.status=Error("WKB interval B differs from input or rotates beyond its declared fixed-direction stage");return result;
        }
        long double theta=0,strain=0;
        for(int d=0;d<3;++d)for(int e=0;e<3;++e) {
          const double g=input.integratedVelocityGradient[3*d+e];
          if(!std::isfinite(g)){result.status=Error("WKB full velocity-gradient integral is unavailable");return result;}
          if(d==e)theta+=g;
          strain+=static_cast<long double>(b[d])*g*b[e];
        }
        if(!Finite(input.integratedAlfvenSpeedGradient)||!std::isfinite(input.relativeFlowDotGradLogB)||
            !std::isfinite(input.relativeFlowDotGradLogRho)||!std::isfinite(input.lateralGridSweptVolumeM3)) {
          result.status=Error("WKB background advection/gradient or geometric integral is nonfinite");return result;
        }
        const auto logB=std::log(Norm(input.finalMagneticFieldT)/Norm(input.initialMagneticFieldT));
        const auto logRho=std::log(static_cast<long double>(input.finalMassDensityKgM3)/cell.massDensityKgM3);
        if(!Roundoff(logB-strain+theta+input.relativeFlowDotGradLogB,1+std::fabs(logB)+std::fabs(strain)+std::fabs(theta)+std::fabs(input.relativeFlowDotGradLogB))||
            !Roundoff(logRho+theta+input.relativeFlowDotGradLogRho,1+std::fabs(logRho)+std::fabs(theta)+std::fabs(input.relativeFlowDotGradLogRho))) {
          result.status=Error("supplied WKB background violates ideal induction/continuity in its declared grid frame");return result;
        }
        const auto swept=tube.faces[i+1].gridSweptVolumeM3-tube.faces[i].gridSweptVolumeM3+input.lateralGridSweptVolumeM3;
        const auto deltaV=static_cast<long double>(input.finalVolumeM3)-cell.volumeM3;
        if(!Roundoff(deltaV-swept,cell.volumeM3+input.finalVolumeM3+std::fabs(swept)+
            std::fabs(tube.faces[i+1].gridSweptVolumeM3)+std::fabs(tube.faces[i].gridSweptVolumeM3))) {
          result.status=Error("WKB supplied geometry violates its discrete geometric conservation law");return result;
        }
        const long double initialMass=input.initialCellMassKg,finalMass=input.finalCellMassKg;
        const auto transportedMass=tube.faces[i+1].relativeMassTransportKg-tube.faces[i].relativeMassTransportKg;
        if(!Roundoff(finalMass-initialMass+transportedMass,std::fabs(finalMass)+std::fabs(initialMass)+
            std::fabs(tube.faces[i+1].relativeMassTransportKg)+std::fabs(tube.faces[i].relativeMassTransportKg))) {
          result.status=Error("WKB interval mass flux/end states violate source-free continuity");return result;
        }
        for(int c=0;c<4;++c)for(double energy:cell.magneticEnergyJ[c]) {
          double value;result.status=Represent(static_cast<long double>(energy)/cell.volumeM3,&value,true);
          if(!result.status.ok())return result;
          density[i][c].push_back(value);
        }
      }
      const auto lower=EvaluateBoundaryFlux(tube.lower,density.front(),tube.faces.front(),mu);
      if(!lower.status.ok()){result.status=lower.status;return result;}
      const auto upper=EvaluateBoundaryFlux(tube.upper,density.back(),tube.faces.back(),mu);
      if(!upper.status.ok()){result.status=upper.status;return result;}
      std::vector<FourBins> faceEnergy(count+1),faceAction(count+1);
      for(std::size_t f=0;f<=count;++f)for(int c=0;c<4;++c) {
        const int sigma=c/2==0?1:-1;const auto& face=tube.faces[f];
        const long double q=face.plasmaSweptVolumeM3-face.gridSweptVolumeM3+sigma*face.alfvenSweptVolumeM3;
        const long double actionQ=face.relativeFlowOverAlfvenAreaTimeM2S+sigma*face.signedAreaTimeM2S;
        for(std::size_t j=0;j<bins;++j) {
          double magnetic,action;
          if(f==0){magnetic=-lower.flux.outwardMagneticJ[c][j];action=-lower.flux.outwardActionJs[c][j];}
          else if(f==count){magnetic=upper.flux.outwardMagneticJ[c][j];action=upper.flux.outwardActionJs[c][j];}
          else {
            const double upstream=q>=0?density[f-1][c][j]:density[f][c][j];
            result.status=Represent(q*upstream,&magnetic);if(!result.status.ok())return result;
            const long double low=before.grid.edgesPerM[j],width=static_cast<long double>(before.grid.edgesPerM[j+1])-low;
            result.status=Represent(actionQ*2*std::log1p(width/low)/width*upstream,&action);if(!result.status.ok())return result;
          }
          faceEnergy[f][c].push_back(magnetic);faceAction[f][c].push_back(action);
        }
      }
      for(std::size_t i=0;i<count;++i) {
        const auto& old=before.cells[indices[i]];const auto& input=tube.cells[i];FourBins advected;
        for(int c=0;c<4;++c) {
          const int sigma=c/2==0?1:-1;
          const auto qleft=tube.faces[i].plasmaSweptVolumeM3-tube.faces[i].gridSweptVolumeM3+sigma*tube.faces[i].alfvenSweptVolumeM3;
          const auto qright=tube.faces[i+1].plasmaSweptVolumeM3-tube.faces[i+1].gridSweptVolumeM3+sigma*tube.faces[i+1].alfvenSweptVolumeM3;
          const long double donor=(std::max(-static_cast<long double>(qleft),0.0L)+std::max(static_cast<long double>(qright),0.0L))/old.volumeM3;
          result.maximumSpatialDonorFraction=std::max(result.maximumSpatialDonorFraction,static_cast<double>(donor));
          if(donor>1){result.status=Unresolved("WKB interval exceeds the positive donor-cell CFL; supply a resolved subinterval");return result;}
          for(std::size_t j=0;j<bins;++j) {
            const long double incoming=std::max(static_cast<long double>(faceEnergy[i][c][j]),0.0L)+
                std::max(-static_cast<long double>(faceEnergy[i+1][c][j]),0.0L);
            // Both neighbors consume the identical stored face flux. A
            // marginal-CFL roundoff deficit is rejected, never clipped.
            const long double outgoing=std::max(-static_cast<long double>(faceEnergy[i][c][j]),0.0L)+
                std::max(static_cast<long double>(faceEnergy[i+1][c][j]),0.0L);
            const long double retained=static_cast<long double>(old.magneticEnergyJ[c][j])-outgoing;
            if(retained<0){result.status=Unresolved("WKB donor flux exceeds resolved cell energy at the CFL boundary");return result;}
            double value;result.status=Represent(retained+incoming,&value,true);
            if(!result.status.ok())return result;
            advected[c].push_back(value);
          }
        }
        const auto b=Unit(input.initialMagneticFieldT);long double theta=0,strain=0;
        for(int d=0;d<3;++d)for(int e=0;e<3;++e) {
          const auto g=input.integratedVelocityGradient[3*d+e];if(d==e)theta+=g;
          strain+=static_cast<long double>(b[d])*g*b[e];
        }
        const auto nextVA=Norm(input.finalMagneticFieldT)/std::sqrt(static_cast<long double>(mu)*input.finalMassDensityKgM3);
        const auto ray=PrepareSpectralWork(before.grid,advected,old.alfvenSpeedMPerS,static_cast<double>(nextVA),
            static_cast<double>(theta),static_cast<double>(strain),static_cast<double>(Dot(b,input.integratedAlfvenSpeedGradient)));
        if(!ray.status.ok()){result.status=ray.status;return result;}
        auto& cell=updated[indices[i]];cell.magneticEnergyJ=ray.magneticEnergyAfterJ;
        cell.volumeM3=input.finalVolumeM3;cell.meanMagneticFieldT=static_cast<double>(Norm(input.finalMagneticFieldT));
        cell.massDensityKgM3=input.finalMassDensityKgM3;cell.alfvenSpeedMPerS=static_cast<double>(nextVA);cell.domain=input.finalDomain;
        for(int c=0;c<4;++c) {
          auto& total=result.ledger[c];const auto& part=ray.ledger[c];
          // Every accepted cell crossing enters the exact external authority
          // BEFORE floating diagnostic aggregation. Otherwise a small escape
          // could disappear while summing a line's widely separated energies,
          // even if the subsequent inter-interval accumulator were exact.
          auto& account=result.externalAfter;
          result.status=AccountIncrement(part.lowerSpectralOutgoingMagneticJ,&account.lowerMagneticExact[c],
              &account.lowerMagneticJ[c],&account.lowerMagneticRoundoffJ[c],&account.lowerMagneticDiagnosticBoundJ[c]);
          if(!result.status.ok())return result;
          result.status=AccountIncrement(part.upperSpectralOutgoingMagneticJ,&account.upperMagneticExact[c],
              &account.upperMagneticJ[c],&account.upperMagneticRoundoffJ[c],&account.upperMagneticDiagnosticBoundJ[c]);
          if(!result.status.ok())return result;
          result.status=AccountIncrement(part.lowerSpectralOutgoingActionJs,&account.lowerActionExact[c],
              &account.lowerActionJs[c],&account.lowerActionRoundoffJs[c],&account.lowerActionDiagnosticBoundJs[c]);
          if(!result.status.ok())return result;
          result.status=AccountIncrement(part.upperSpectralOutgoingActionJs,&account.upperActionExact[c],
              &account.upperActionJs[c],&account.upperActionRoundoffJs[c],&account.upperActionDiagnosticBoundJs[c]);
          if(!result.status.ok())return result;
          total.magneticBeforeJ+=static_cast<double>(Sum(old.magneticEnergyJ[c]));total.magneticAfterJ+=part.magneticAfterJ;
          total.actionBeforeJs+=static_cast<double>(Action(before.grid,old.magneticEnergyJ[c],old.alfvenSpeedMPerS));
          total.actionAfterJs+=part.actionAfterJs;
          total.lowerSpectralOutgoingMagneticJ+=part.lowerSpectralOutgoingMagneticJ;
          total.upperSpectralOutgoingMagneticJ+=part.upperSpectralOutgoingMagneticJ;
          total.lowerSpectralOutgoingActionJs+=part.lowerSpectralOutgoingActionJs;total.upperSpectralOutgoingActionJs+=part.upperSpectralOutgoingActionJs;
          total.magneticCompressionWorkJ+=part.magneticCompressionWorkJ;
          // Material side walls have no relative advective wave flux, yet
          // their physical U*P pressure-work flux is generally nonzero. This
          // is the same declared frozen donor-pressure quadrature as the cap
          // traces; its time/splitting error remains a convergence quantity.
          double pressure;result.status=Represent(Sum(old.magneticEnergyJ[c])/old.volumeM3*
              input.lateralGridSweptVolumeM3,&pressure);if(!result.status.ok())return result;
          total.lateralPressureWorkOutJ+=pressure;total.pressureWorkBoundaryOutJ+=pressure;
        }
      }
      for(int c=0;c<4;++c)for(const auto* boundary:{&lower.flux,&upper.flux}) {
        auto& total=result.ledger[c];
        for(std::size_t j=0;j<bins;++j) {
          const double energy=boundary->outwardMagneticJ[c][j],action=boundary->outwardActionJs[c][j];
          if(energy<0){total.spatialIncomingMagneticJ-=energy;total.spatialIncomingActionJs-=action;}
          else {total.spatialOutgoingMagneticJ+=energy;total.spatialOutgoingActionJs+=action;}
          total.pressureWorkBoundaryOutJ+=boundary->outwardPressureWorkJ[c][j];
        }
      }
    }
    if(coveredCells!=before.cells.size()){result.status=Error("WKB interval omitted an authoritative field line");return result;}
    for(int c=0;c<4;++c) {
      auto& a=result.externalAfter;auto& l=result.ledger[c];
      l.totalBackgroundCompressionReceiptJ=-2*l.magneticCompressionWorkJ;
      l.totalBackgroundMechanicalReceiptJ=l.totalBackgroundCompressionReceiptJ-l.pressureWorkBoundaryOutJ;
      l.energyRoundoffResidualJ=l.magneticAfterJ-l.magneticBeforeJ-l.spatialIncomingMagneticJ+l.spatialOutgoingMagneticJ+
          l.lowerSpectralOutgoingMagneticJ+l.upperSpectralOutgoingMagneticJ-l.magneticCompressionWorkJ;
      l.actionDiscretizationResidualJs=l.actionAfterJs-l.actionBeforeJs-l.spatialIncomingActionJs+l.spatialOutgoingActionJs+
          l.lowerSpectralOutgoingActionJs+l.upperSpectralOutgoingActionJs;
      for(double value:{a.lowerMagneticJ[c],a.upperMagneticJ[c],a.lowerActionJs[c],a.upperActionJs[c],
          a.lowerMagneticRoundoffJ[c],a.upperMagneticRoundoffJ[c],a.lowerActionRoundoffJs[c],a.upperActionRoundoffJs[c],
          l.magneticBeforeJ,l.magneticAfterJ,l.magneticCompressionWorkJ,l.totalBackgroundCompressionReceiptJ,
          l.pressureWorkBoundaryOutJ,l.totalBackgroundMechanicalReceiptJ,l.actionBeforeJs,l.actionAfterJs,
          l.energyRoundoffResidualJ,l.actionDiscretizationResidualJs})
        if(!std::isfinite(value)){result.status=Unresolved("WKB global energy/action account exceeds finite representation");return result;}
    }
    result.candidate=PrepareEvolutionCandidate(before,updated,interval.target);
    if(!result.candidate.status.ok()){result.status=result.candidate.status;return result;}
    result.externalAfter.lastSnapshotIdentity=result.candidate.snapshot->identity;
    result.status=Status::Ok();return result;
  } catch(const std::exception& error) {
    result.status=Error(std::string("WKB staging failed without publication: ")+error.what());return result;
  }
}

} } } }
