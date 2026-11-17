#include "sep_polarized_wkb_validation.h"
#include "sep_polarized_wkb_transport.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace SEP { namespace Testing { namespace {
namespace PW=Turbulence::Polarized;namespace W=PW::WKB;
const double Mu0=1.256637062e-6,B0=5e-9,Rho0=5e6*1.67262192369e-27;
double VA(){return B0/std::sqrt(Mu0*Rho0);}
bool Near(double a,double b,double relative=1e-11) {
  return std::isfinite(a)&&std::isfinite(b)&&std::fabs(static_cast<long double>(a)-b)<=
      relative*std::max(std::fabs(static_cast<long double>(a)),std::fabs(static_cast<long double>(b)));
}
Result Finish(bool pass,const std::string& message) {
  Result result;result.status=pass?Status::Pass:Status::Fail;result.message=message;
  result.configuration.push_back("capability=controlled-ideal-WKB-prepare-only;native_publication=false;production_qualification=false");
  result.configuration.push_back("state=4-independent-magnetic-J-per-k-bin;reconstruction=Ebin/Delta-k");
  result.configuration.push_back("incoming=initial-J/m3/bin;zero-k-incoming;escaped-energy/action=external-account;reflection/cascade=OFF");
  return result;
}
PW::Initialization Initial(std::size_t spatial,std::size_t spectral) {
  PW::Initialization input;const auto grid=PW::ResolveCenteredLogGrid(spectral,1e-10,1e-8);
  input.grid.centersPerM=grid.grid.centersPerM;input.grid.logWidth=grid.grid.logWidth;input.powerLawIndex=0;
  // Deliberately unequal manufactured components verify independent transport;
  // these are operator fixtures, not a change to native initial L/R=0.5.
  input.plus.known=input.minus.known=true;input.plus.leftFraction=.2;input.minus.leftFraction=.7;
  input.plus.provenance=input.minus.provenance="explicit manufactured component-independence fixture";
  input.vacuumPermeabilityHPerM=Mu0;input.epochS=0;input.validUntilS=1e12;input.geometryGeneration=1;
  input.backgroundIdentity="controlled ideal-MHD input";input.configurationFingerprint="controlled-WKB;reflection=off;cascade=off";
  input.initializerProvenance="explicit SI manufactured data; not observational or production qualification";
  for(std::size_t i=0;i<spatial;++i) {
    PW::InitialCell cell;cell.fieldLine=0;cell.segment=i;cell.volumeM3=1e8;
    cell.meanMagneticFieldT=B0;cell.massDensityKgM3=Rho0;cell.ePlusJ=1e-20;cell.eMinusJ=2e-20;
    cell.domain.authority=Transport::ResonantCollision::DomainAuthority::ValidationStudy;
    cell.domain.provenance="controlled WKB equation verification; no selected physical error limits";
    input.cells.push_back(cell);
  }
  return input;
}
W::FaceInterval ConstantFace(double area,double u,double grid,double va,double b,double rho,double dt) {
  W::FaceInterval face;face.initialNormal=face.finalNormal={{1,0,0}};
  face.initialAreaM2=face.finalAreaM2=area;
  face.initialMagneticFieldT=face.finalMagneticFieldT={{b,0,0}};
  face.initialMassDensityKgM3=face.finalMassDensityKgM3=rho;
  face.initialPlasmaVelocityMPerS=face.finalPlasmaVelocityMPerS={{u,0,0}};
  face.initialGridVelocityMPerS=face.finalGridVelocityMPerS={{grid,0,0}};
  const double polarity=b>0?1:-1;
  face.plasmaSweptVolumeM3=area*u*dt;face.gridSweptVolumeM3=area*grid*dt;
  face.alfvenSweptVolumeM3=area*polarity*va*dt;face.relativeMassTransportKg=rho*area*(u-grid)*dt;
  face.relativeFlowOverAlfvenAreaTimeM2S=area*(u-grid)/va*dt;face.signedAreaTimeM2S=area*polarity*dt;
  face.integralsKnown=face.noCharacteristicReversalKnown=true;return face;
}
W::FourBins Density(const PW::Cell& cell) {
  auto result=cell.magneticEnergyJ;for(auto& component:result)for(double& e:component)e/=cell.volumeM3;return result;
}
Result IncomingBoundary() {
  const auto prepared=PW::PrepareInitialization(Initial(2,16));
  if(!prepared.status.ok())return Finish(false,prepared.status.message);
  const auto& state=*prepared.snapshot;const auto lower=W::CaptureInitialReservoir(state,0,W::Side::Lower);
  const auto upper=W::CaptureInitialReservoir(state,0,W::Side::Upper);
  if(!lower.status.ok()||!upper.status.ok())return Finish(false,"initial reservoir capture failed");
  auto interior=Density(state.cells.back());for(auto& c:interior)for(double& v:c)v*=3;
  const double dt=10,u=2e5;
  const auto face=ConstantFace(1,u,0,VA(),B0,Rho0,dt);
  const auto in=W::EvaluateBoundaryFlux(lower.reservoir,interior,face,Mu0);
  const auto out=W::EvaluateBoundaryFlux(upper.reservoir,interior,face,Mu0);
  bool pass=in.status.ok()&&out.status.ok()&&in.flux.incoming[0]&&in.flux.incoming[1]&&
      !out.flux.incoming[0]&&!out.flux.incoming[1];
  if(pass)for(int c=0;c<4;++c)for(std::size_t j=0;j<16;++j) {
    const double speed=u+(c/2==0?VA():-VA());
    pass=pass&&Near(in.flux.outwardMagneticJ[c][j],-speed*dt*lower.reservoir.magneticDensityJPerM3[c][j])&&
        Near(out.flux.outwardMagneticJ[c][j],speed*dt*interior[c][j]);
  }
  const auto moving=ConstantFace(1,1e5,1.3e5,VA(),B0,Rho0,dt);
  const auto classified=W::EvaluateBoundaryFlux(upper.reservoir,interior,moving,Mu0);
  const auto inwardB=ConstantFace(1,1e5,1.3e5,VA(),-B0,Rho0,dt);
  const auto reversed=W::EvaluateBoundaryFlux(upper.reservoir,interior,inwardB,Mu0);
  pass=pass&&classified.status.ok()&&reversed.status.ok()&&!classified.flux.incoming[0]&&classified.flux.incoming[1]&&
      reversed.flux.incoming[0]&&!reversed.flux.incoming[1];
  // Public wire/scratch structs are mutable before an outer owner pins them.
  // Changed densities or provenance must not silently redefine a fixed INITIAL
  // reservoir while retaining its original snapshot/configuration metadata.
  auto changed=lower.reservoir;changed.magneticDensityJPerM3[0][0]*=2;
  pass=pass&&!W::EvaluateBoundaryFlux(changed,interior,face,Mu0).status.ok();
  changed=lower.reservoir;changed.initialEpochS=1;
  pass=pass&&!W::EvaluateBoundaryFlux(changed,interior,face,Mu0).status.ok();
  changed=lower.reservoir;changed.side=W::Side::Upper;
  pass=pass&&!W::EvaluateBoundaryFlux(changed,interior,face,Mu0).status.ok();
  PW::EvolutionMetadata later;later.epochS=1;later.validUntilS=state.validUntilS;later.geometryGeneration=state.geometryGeneration;
  later.backgroundIdentity=state.backgroundIdentity;later.operationProvenance="manufactured no-op evolution must not refresh initial reservoirs";
  const auto evolved=PW::PrepareEvolutionCandidate(state,state.cells,later);
  pass=pass&&evolved.status.ok()&&!W::CaptureInitialReservoir(*evolved.snapshot,0,W::Side::Lower).status.ok();
  return Finish(pass,"all four captured initial intensities supply only signed incoming characteristics; outflow has no hidden Eminus reset, polarity/grid motion classify independently, and mutated/evolved initial reservoirs reject");
}
Result ExactRayAndExternalAccount() {
  const auto prepared=PW::PrepareInitialization(Initial(1,16));if(!prepared.status.ok())return Finish(false,prepared.status.message);
  const auto& grid=prepared.snapshot->grid;const auto bins=grid.centersPerM.size();const double h=grid.logWidth;
  W::FourBins energy;for(int c=0;c<4;++c){energy[c].assign(bins,0);energy[c][7]=(c+1)*1e-20;}
  const auto extension=W::PrepareSpectralWork(grid,energy,VA(),VA()*std::exp(h),0,h,0);
  bool pass=extension.status.ok();double largest=0;
  if(pass)for(int c=0;c<4;++c) {
    const double e=(c+1)*1e-20;largest=std::max(largest,std::fabs(extension.magneticEnergyAfterJ[c][6]/e-1));
    pass=pass&&Near(extension.ledger[c].magneticAfterJ,e)&&extension.ledger[c].magneticCompressionWorkJ==0&&
        Near(extension.ledger[c].actionBeforeJs,extension.ledger[c].actionAfterJs);
  }
  // A whole initial low-k bin leaves the lower spectral edge during an ideal
  // isotropic expansion. Independent power integration evaluates the outgoing
  // energy at each crossing, not at final time: integral x^(-3/2) dx.
  for(int c=0;c<4;++c){energy[c].assign(bins,0);energy[c][0]=(c+1)*1e-20;}
  const auto expansion=W::PrepareSpectralWork(grid,energy,VA(),VA()*std::exp(-h/2),3*h,h,0);
  pass=pass&&expansion.status.ok();
  if(expansion.status.ok())for(int c=0;c<4;++c) {
    const double e=(c+1)*1e-20,lo=grid.edgesPerM[0],hi=grid.edgesPerM[1];
    const double escaped=e/(hi-lo)*2*lo*(1-std::sqrt(lo/hi));
    pass=pass&&Near(expansion.ledger[c].lowerSpectralOutgoingMagneticJ,escaped)&&
        Near(expansion.ledger[c].lowerSpectralOutgoingActionJs,expansion.ledger[c].actionBeforeJs)&&
        Near(expansion.ledger[c].magneticAfterJ+escaped+expansion.ledger[c].totalBackgroundCompressionReceiptJ/2,e)&&
        expansion.ledger[c].upperSpectralOutgoingMagneticJ==0;
  }
  return Finish(pass&&largest<1e-12,"frozen ideal strain shifts exact k characteristics without mixing helicities; escaped magnetic/total-wave action is evaluated at crossing and kept separate from signed compression work");
}
// Three-point Gauss integration is exact for this factored quartic. This
// reference integrates the continuous manufactured profile independently of
// the production overlap remapper and its piecewise-constant representation.
long double ProfileIntegral(long double low,long double high,long double a,long double b) {
  low=std::max(low,a);high=std::min(high,b);if(high<=low)return 0;
  const long double midpoint=(low+high)/2,half=(high-low)/2,node=std::sqrt(3.0L/5);
  auto value=[&](long double x){return (x-a)*(x-a)*(b-x)*(b-x);};
  return half*(8.0L/9*value(midpoint)+5.0L/9*(value(midpoint-half*node)+value(midpoint+half*node)));
}
double SpectralError(std::size_t bins,double* workError,double* actionError) {
  const auto grid=PW::ResolveCenteredLogGrid(bins,1e-10,1e-8).grid;
  W::FourBins energy;const double scale=1e-20,k0=1e-9,horizon=.25;
  const long double normalization=std::pow(2.0L,5)/30;
  for(int c=0;c<4;++c)for(std::size_t j=0;j<bins;++j)
    energy[c].push_back(static_cast<double>((c+1)*scale*ProfileIntegral(grid.edgesPerM[j]/k0,grid.edgesPerM[j+1]/k0,2,4)/normalization));
  const auto initial=energy;double va=VA();const auto steps=static_cast<std::size_t>(std::ceil(horizon/(.35*grid.logWidth)));
  const double stretch=horizon/steps;long double accumulatedWork[4]={0,0,0,0};
  for(std::size_t step=0;step<steps;++step) {
    const double next=va*std::exp(-stretch/2);
    const auto result=W::PrepareSpectralWork(grid,energy,va,next,3*stretch,stretch,0);
    if(!result.status.ok())return std::numeric_limits<double>::infinity();
    for(int c=0;c<4;++c)accumulatedWork[c]+=result.ledger[c].magneticCompressionWorkJ;
    energy=result.magneticEnergyAfterJ;va=next;
  }
  long double error=0,total=0,initialTotal=0,finalTotal=0,work=0,action=0;
  for(int c=0;c<4;++c) {
    work+=accumulatedWork[c];
    for(std::size_t j=0;j<bins;++j) {
      const long double expected=(c+1)*scale*std::exp(-1.5L*horizon)*ProfileIntegral(
          grid.edgesPerM[j]*std::exp(horizon)/k0,grid.edgesPerM[j+1]*std::exp(horizon)/k0,2,4)/normalization;
      error+=std::fabs(energy[c][j]-expected);total+=expected;initialTotal+=initial[c][j];finalTotal+=energy[c][j];
      action+=2*static_cast<long double>(energy[c][j])*std::log(static_cast<long double>(grid.edgesPerM[j+1])/grid.edgesPerM[j])/
          (va*(static_cast<long double>(grid.edgesPerM[j+1])-grid.edgesPerM[j]));
    }
  }
  *workError=static_cast<double>(std::fabs(finalTotal-initialTotal-work)/initialTotal);
  // Independent CONTINUOUS action oracle. On a geometric grid this particular
  // homogeneous remap preserves its discrete1/k moment up to roundoff, so
  // demanding monotone decay of roundoff itself is not a physical refinement
  // test. Instead compare with the exact integral of quartic(x)/x.
  auto antiderivative=[](long double x) {
    const long double a=2,b=4;
    return std::pow(x,4)/4-2*(a+b)*std::pow(x,3)/3+(a*a+4*a*b+b*b)*x*x/2-
        2*a*b*(a+b)*x+a*a*b*b*std::log(x);
  };
  const long double exactAction=2*10*scale/(VA()*k0)*(antiderivative(4)-antiderivative(2))/normalization;
  *actionError=static_cast<double>(std::fabs(action-exactAction)/exactAction);
  return static_cast<double>(error/total);
}
Result RayConvergence() {
  double w[3]={0,0,0},a[3]={0,0,0},error[3];const std::size_t bins[3]={64,128,256};
  for(int i=0;i<3;++i)error[i]=SpectralError(bins[i],&w[i],&a[i]);
  const double order=std::log(error[1]/error[2])/std::log(2.0);
  const bool pass=std::isfinite(error[0])&&error[1]<error[0]&&error[2]<error[1]&&error[2]<.12&&order>.65&&
      w[2]<1e-12&&a[2]<a[1];
  auto result=Finish(pass,"actual Ebin/Delta-k WKB remapping converges to an independently integrated compact quartic spectrum; physical work and numerical action-projection error remain distinct");
  result.metrics.push_back({"fine_continuous_spectrum_L1",error[2],.12,"<","relative"});
  result.metrics.push_back({"observed_k_refinement_order",order,.65,">","dimensionless"});
  result.metrics.push_back({"magnetic_work_balance_error",w[2],1e-12,"<","relative"});
  result.metrics.push_back({"fine_continuous_action_error",a[2],a[1],"<","relative"});return result;
}
// Four-point Gauss integration exactly integrates x^power times the quartic
// for power0,1,2. The production spatial operator uses no such quadrature.
long double SpatialIntegral(long double low,long double high,int power) {
  const long double a=.8L,b=1.2L;low=std::max(low,a);high=std::min(high,b);if(high<=low)return 0;
  const long double middle=(low+high)/2,half=(high-low)/2;
  const long double root=std::sqrt(6.0L/5),outer=std::sqrt((3+2*root)/7),inner=std::sqrt((3-2*root)/7);
  const long double wo=(18-std::sqrt(30.0L))/36,wi=(18+std::sqrt(30.0L))/36;
  auto f=[&](long double x){return std::pow(x,power)*(x-a)*(x-a)*(b-x)*(b-x);};
  return half*(wo*(f(middle-half*outer)+f(middle+half*outer))+wi*(f(middle-half*inner)+f(middle+half*inner)));
}
long double LogProfileActionIntegral(long double a,long double b) {
  auto f=[&](long double x){return std::pow(x,4)/4-2*(a+b)*std::pow(x,3)/3+
      (a*a+4*a*b+b*b)*x*x/2-2*a*b*(a+b)*x+a*a*b*b*std::log(x);};
  return f(b)-f(a);
}
double StaticGradientError(std::size_t count,double* actionError,double* energyError) {
  // Static ideal equilibrium: U=0, radial curl-free B~r^-2, rho~r^-6,
  // constant gas pressure and no gravity. Thus VA=VA0*r/r0; this is not the
  // nonideal canonical Leblanc wind. All true masses/volumes are analytic.
  const double r0=1e11,k0=1e-9,amplitude=1e-20,aRate=VA()/r0,horizon=.1/aRate;
  auto initialization=Initial(count,count);
  for(std::size_t i=0;i<count;++i) {
    const double lo=.5+static_cast<double>(i)/count,hi=.5+static_cast<double>(i+1)/count,mid=(lo+hi)/2;
    auto& cell=initialization.cells[i];cell.volumeM3=r0*(hi*hi*hi-lo*lo*lo)/3;
    cell.meanMagneticFieldT=B0/(mid*mid);cell.massDensityKgM3=Rho0/std::pow(mid,6);
    cell.ePlusJ=cell.eMinusJ=0;
  }
  const auto declared=PW::PrepareInitialization(initialization);if(!declared.status.ok())return std::numeric_limits<double>::infinity();
  auto cells=declared.snapshot->cells;const auto& grid=declared.snapshot->grid;
  const auto spatialNorm=SpatialIntegral(.8L,1.2L,2),spectralNorm=std::pow(2.0L,5)/30;
  for(std::size_t i=0;i<count;++i)for(int c=0;c<4;++c)for(std::size_t j=0;j<count;++j) {
    const long double lo=.5L+static_cast<long double>(i)/count,hi=.5L+static_cast<long double>(i+1)/count;
    cells[i].magneticEnergyJ[c][j]=static_cast<double>((c+1)*amplitude*SpatialIntegral(lo,hi,2)/spatialNorm*
        ProfileIntegral(grid.edgesPerM[j]/k0,grid.edgesPerM[j+1]/k0,2,4)/spectralNorm);
  }
  PW::EvolutionMetadata initialData;initialData.epochS=0;initialData.validUntilS=1e12;initialData.geometryGeneration=1;
  initialData.backgroundIdentity="static radial exact ideal equilibrium";
  initialData.operationProvenance="manufactured compact Cauchy data for equation verification, not a physical injection event";
  auto current=PW::PrepareEvolutionCandidate(*declared.snapshot,cells,initialData).snapshot;
  if(!current)return std::numeric_limits<double>::infinity();
  auto external=W::BeginExternalAccount(*declared.snapshot).account;external.lastSnapshotIdentity=current->identity;
  const auto lower=W::CaptureInitialReservoir(*declared.snapshot,0,W::Side::Lower).reservoir;
  const auto upper=W::CaptureInitialReservoir(*declared.snapshot,0,W::Side::Upper).reservoir;
  double limit=std::numeric_limits<double>::infinity();
  for(std::size_t i=0;i<count;++i) {
    const double hi=.5+static_cast<double>(i+1)/count;
    limit=std::min(limit,.3*current->cells[i].volumeM3/(VA()*hi*hi*hi));
  }
  const auto steps=static_cast<std::size_t>(std::ceil(horizon/limit));const double dt=horizon/steps;
  long double roundoff=0;
  for(std::size_t step=0;step<steps;++step) {
    W::Interval interval;interval.durationS=dt;interval.geometryMode=W::GeometryMode::MaterialTube;
    interval.zeroExteriorSpectralIncoming=true;interval.externalBefore=external;
    interval.target.epochS=(step+1)*dt;interval.target.validUntilS=1e12;interval.target.geometryGeneration=1;
    interval.target.backgroundIdentity="static radial exact ideal equilibrium";
    interval.target.operationProvenance="explicit material U=grid=0 WKB donor-cell/ray interval";
    W::TubeInterval tube;tube.fieldLine=0;tube.lower=lower;tube.upper=upper;
    for(std::size_t f=0;f<=count;++f) {
      const double r=.5+static_cast<double>(f)/count;
      tube.faces.push_back(ConstantFace(r*r,0,0,VA()*r,B0/(r*r),Rho0/std::pow(r,6),dt));
    }
    for(std::size_t i=0;i<count;++i) {
      const double lo=.5+static_cast<double>(i)/count,hi=.5+static_cast<double>(i+1)/count;
      W::CellInterval cell;cell.initialMagneticFieldT=cell.finalMagneticFieldT={{current->cells[i].meanMagneticFieldT,0,0}};
      cell.finalVolumeM3=current->cells[i].volumeM3;cell.finalMassDensityKgM3=current->cells[i].massDensityKgM3;
      cell.initialCellMassKg=cell.finalCellMassKg=Rho0*r0*(std::pow(lo,-3)-std::pow(hi,-3))/3;
      cell.integratedAlfvenSpeedGradient={{aRate*dt,0,0}};cell.finalDomain=current->cells[i].domain;
      cell.inputsKnown=cell.fixedFieldDirectionKnown=cell.lateralTubeMaterialKnown=true;tube.cells.push_back(cell);
    }
    interval.tubes.push_back(tube);const auto stage=W::PrepareInterval(*current,interval);
    if(!stage.status.ok())return std::numeric_limits<double>::infinity();
    for(const auto& component:stage.ledger)roundoff+=std::fabs(component.energyRoundoffResidualJ);
    current=stage.candidate.snapshot;external=stage.externalAfter;
  }
  long double error=0,total=0,action=0;
  for(std::size_t i=0;i<count;++i)for(int c=0;c<4;++c)for(std::size_t j=0;j<count;++j) {
    const double sigma=c/2==0?1:-1,stretch=std::exp(sigma*aRate*horizon);
    const long double lo=.5L+static_cast<long double>(i)/count,hi=.5L+static_cast<long double>(i+1)/count;
    const long double exact=(c+1)*amplitude*SpatialIntegral(lo/stretch,hi/stretch,2)/spatialNorm*
        ProfileIntegral(grid.edgesPerM[j]*stretch/k0,grid.edgesPerM[j+1]*stretch/k0,2,4)/spectralNorm;
    const auto energy=current->cells[i].magneticEnergyJ[c][j];error+=std::fabs(energy-exact);total+=exact;
    action+=2*static_cast<long double>(energy)*std::log(static_cast<long double>(grid.edgesPerM[j+1])/grid.edgesPerM[j])/
        (current->cells[i].alfvenSpeedMPerS*(static_cast<long double>(grid.edgesPerM[j+1])-grid.edgesPerM[j]));
  }
  const long double exactAction=2*10*amplitude/(VA()*k0)*SpatialIntegral(.8L,1.2L,1)/spatialNorm*
      LogProfileActionIntegral(2,4)/spectralNorm;
  *actionError=static_cast<double>(std::fabs(action-exactAction)/exactAction);*energyError=static_cast<double>(roundoff/(10*amplitude));
  return static_cast<double>(error/total);
}
Result StaticGradientConvergence() {
  double action[3]={0,0,0},energy[3]={0,0,0},error[3];const std::size_t resolution[3]={32,64,128};
  for(int i=0;i<3;++i)error[i]=StaticGradientError(resolution[i],&action[i],&energy[i]);
  const double order=std::log(error[1]/error[2])/std::log(2.0);
  const bool pass=std::isfinite(error[0])&&error[1]<error[0]&&error[2]<error[1]&&error[2]<.15&&order>.65&&
      action[2]<action[1]&&action[2]<.05&&energy[2]<1e-10;
  auto result=Finish(pass,"full spatial-plus-logk evolution converges to independent radial ideal-MHD characteristics r~exp(sigma*a*t), k~exp(-sigma*a*t), with constant packet magnetic energy/action");
  result.metrics.push_back({"fine_spatial_spectral_L1",error[2],.15,"<","relative"});
  result.metrics.push_back({"observed_joint_refinement_order",order,.65,">","dimensionless"});
  result.metrics.push_back({"fine_continuous_action_error",action[2],.05,"<","relative"});
  result.metrics.push_back({"accumulated_energy_roundoff",energy[2],1e-10,"<","relative"});return result;
}
struct ControlledFixture {
  std::shared_ptr<const PW::Snapshot> state;
  W::InitialReservoir lower,upper;W::ExternalSpectralAccount external;
};
ControlledFixture UniformFixture(std::size_t cells=8,std::size_t bins=16) {
  ControlledFixture f;const auto initialized=PW::PrepareInitialization(Initial(cells,bins));
  if(!initialized.status.ok())return f;
  f.state=initialized.snapshot;f.lower=W::CaptureInitialReservoir(*f.state,0,W::Side::Lower).reservoir;
  f.upper=W::CaptureInitialReservoir(*f.state,0,W::Side::Upper).reservoir;
  f.external=W::BeginExternalAccount(*f.state).account;return f;
}
W::Interval StaticAleInterval(const ControlledFixture& fixture) {
  W::Interval interval;const auto& state=*fixture.state;const auto n=state.cells.size();const double dt=300;
  interval.durationS=dt;interval.zeroExteriorSpectralIncoming=true;interval.externalBefore=fixture.external;
  interval.geometryMode=W::GeometryMode::ControlledExplicitAle;
  interval.target.epochS=state.epochS+dt;interval.target.validUntilS=state.validUntilS;
  interval.target.geometryGeneration=state.geometryGeneration+1;interval.target.backgroundIdentity="uniform physical medium; explicitly moving ALE grid";
  interval.target.operationProvenance="controlled moving-grid FV+WKB; not native production";
  W::TubeInterval tube;tube.fieldLine=0;tube.lower=fixture.lower;tube.upper=fixture.upper;
  for(std::size_t i=0;i<=n;++i) {
    const double grid=.1*VA()*std::sin(2*3.14159265358979323846*i/n);
    tube.faces.push_back(ConstantFace(1,0,grid,VA(),B0,Rho0,dt));
  }
  for(std::size_t i=0;i<n;++i) {
    W::CellInterval cell;cell.initialMagneticFieldT=cell.finalMagneticFieldT={{B0,0,0}};
    cell.finalVolumeM3=state.cells[i].volumeM3+tube.faces[i+1].gridSweptVolumeM3-tube.faces[i].gridSweptVolumeM3;
    cell.finalMassDensityKgM3=Rho0;cell.finalDomain=state.cells[i].domain;
    cell.initialCellMassKg=Rho0*state.cells[i].volumeM3;cell.finalCellMassKg=Rho0*cell.finalVolumeM3;
    cell.inputsKnown=cell.fixedFieldDirectionKnown=cell.lateralTubeMaterialKnown=true;tube.cells.push_back(cell);
  }
  interval.tubes.push_back(tube);return interval;
}
Result AleGeometricConservation() {
  const auto fixture=UniformFixture();if(!fixture.state)return Finish(false,"ALE fixture initialization failed");
  const auto interval=StaticAleInterval(fixture);const auto staged=W::PrepareInterval(*fixture.state,interval);
  if(!staged.status.ok())return Finish(false,staged.status.message);
  bool pass=staged.candidate.snapshot&&!staged.qualifiedForProduction;
  double maximum=0;
  for(std::size_t i=0;i<fixture.state->cells.size();++i)for(int c=0;c<4;++c)for(std::size_t j=0;j<16;++j) {
    const double old=fixture.state->cells[i].magneticEnergyJ[c][j]/fixture.state->cells[i].volumeM3;
    const double now=staged.candidate.snapshot->cells[i].magneticEnergyJ[c][j]/staged.candidate.snapshot->cells[i].volumeM3;
    if(old>0)maximum=std::max(maximum,std::fabs(now/old-1));
  }
  for(const auto& l:staged.ledger)pass=pass&&l.magneticCompressionWorkJ==0&&l.totalBackgroundCompressionReceiptJ==0;
  auto result=Finish(pass&&maximum<1e-12,"arbitrary longitudinal grid motion in a static physical medium preserves each uniform magnetic bin density via the supplied discrete GCL and introduces no compression work");
  result.metrics.push_back({"maximum_ALE_uniform_density_error",maximum,1e-12,"<","relative"});return result;
}
Result RejectedIntervals() {
  const auto fixture=UniformFixture();if(!fixture.state)return Finish(false,"invalid-interval fixture initialization failed");
  const auto original=fixture.state->identity;auto interval=StaticAleInterval(fixture);
  interval.tubes[0].cells.back().finalVolumeM3*=1.01;
  auto result=W::PrepareInterval(*fixture.state,interval);bool pass=!result.status.ok()&&!result.candidate.snapshot;
  interval=StaticAleInterval(fixture);interval.tubes[0].cells.back().integratedVelocityGradient[0]=.1;
  result=W::PrepareInterval(*fixture.state,interval);pass=pass&&!result.status.ok()&&!result.candidate.snapshot;
  interval=StaticAleInterval(fixture);for(auto& face:interval.tubes[0].faces) {
    face.alfvenSweptVolumeM3*=100;face.signedAreaTimeM2S*=100;
  }
  result=W::PrepareInterval(*fixture.state,interval);pass=pass&&!result.status.ok()&&!result.candidate.snapshot;
  interval=StaticAleInterval(fixture);interval.reflectionEnabled=true;
  result=W::PrepareInterval(*fixture.state,interval);pass=pass&&!result.status.ok();
  interval=StaticAleInterval(fixture);interval.zeroExteriorSpectralIncoming=false;
  result=W::PrepareInterval(*fixture.state,interval);pass=pass&&!result.status.ok();
  pass=pass&&fixture.state->identity==original&&PW::ValidateSnapshot(*fixture.state).ok()&&
      fixture.external.lastSnapshotIdentity==original;
  return Finish(pass,"malformed ideal induction/geometry, unresolved donor CFL and unsupported operators reject before candidate publication; input wave arrays and external account remain unchanged");
}
Result MaterialSidePressure() {
  const auto fixture=UniformFixture();if(!fixture.state)return Finish(false,"material pressure fixture failed");
  const auto& state=*fixture.state;const double h=.01,H=1e-5,dt=h/H,stretch=std::exp(h);
  W::Interval interval;interval.durationS=dt;interval.geometryMode=W::GeometryMode::MaterialTube;
  interval.zeroExteriorSpectralIncoming=true;interval.externalBefore=fixture.external;
  interval.target.epochS=dt;interval.target.validUntilS=state.validUntilS;interval.target.geometryGeneration=2;
  interval.target.backgroundIdentity="explicit incompressible material affine deformation";
  interval.target.operationProvenance="uniform-pressure cap/side cancellation under declared donor pressure quadrature";
  W::TubeInterval tube;tube.fieldLine=0;tube.lower=fixture.lower;tube.upper=fixture.upper;
  const double length=state.cells.front().volumeM3;
  for(std::size_t f=0;f<=state.cells.size();++f) {
    const double x=f*length;auto face=ConstantFace(1,H*x,H*x,VA(),B0,Rho0,dt);
    face.finalAreaM2=1/stretch;face.finalMagneticFieldT={{B0*stretch,0,0}};
    face.finalPlasmaVelocityMPerS=face.finalGridVelocityMPerS={{H*x*stretch,0,0}};
    face.signedAreaTimeM2S=-std::expm1(-h)/H;tube.faces.push_back(face);
  }
  for(const auto& old:state.cells) {
    W::CellInterval cell;cell.initialMagneticFieldT={{B0,0,0}};cell.finalMagneticFieldT={{B0*stretch,0,0}};
    cell.finalVolumeM3=old.volumeM3;cell.finalMassDensityKgM3=Rho0;
    cell.initialCellMassKg=cell.finalCellMassKg=Rho0*old.volumeM3;cell.integratedVelocityGradient[0]=h;cell.integratedVelocityGradient[4]=-h;
    cell.lateralGridSweptVolumeM3=-h*old.volumeM3;cell.finalDomain=old.domain;
    cell.inputsKnown=cell.fixedFieldDirectionKnown=cell.lateralTubeMaterialKnown=true;tube.cells.push_back(cell);
  }
  interval.tubes.push_back(tube);const auto staged=W::PrepareInterval(state,interval);
  if(!staged.status.ok())return Finish(false,staged.status.message);
  bool pass=true;double largest=0;
  for(const auto& l:staged.ledger) {
    pass=pass&&Near(l.lateralPressureWorkOutJ,-h*l.magneticBeforeJ)&&l.magneticCompressionWorkJ==0;
    largest=std::max(largest,std::fabs(l.totalBackgroundMechanicalReceiptJ)/l.magneticBeforeJ);
  }
  auto result=Finish(pass&&largest<1e-12,"material incompressible extension includes side-wall U*P work; uniform-pressure cap and side contributions cancel without fictitious background mechanical work");
  result.metrics.push_back({"uniform_pressure_mechanical_receipt",largest,1e-12,"<","relative"});return result;
}
Result ExactExternalArithmetic() {
  using Transport::ExactNonnegativeSum;using Transport::ExactSumApproximation;
  ExactNonnegativeSum total,next;
  bool pass=true;for(double value:{1.0,std::ldexp(1.0,-54),std::ldexp(1.0,-108)}) {
    pass=pass&&Transport::AddExactNonnegative(total,value,&next).ok();total=next;
  }
  // Independent integer locations: double1 is bit1074 in2^-1074 units;
  // the two smaller exact increments occupy bits1020 and966 respectively.
  auto bit=[](const ExactNonnegativeSum& s,int n){return (s.words[n/64]>>(n%64)&1)!=0;};
  pass=pass&&bit(total,1074)&&bit(total,1020)&&bit(total,966);
  ExactSumApproximation display;pass=pass&&Transport::ApproximateExactNonnegative(total,&display).ok()&&
      display.nearest==1&&display.signedRemainder==std::ldexp(1.0,-54)&&display.absoluteRoundingBound>=std::ldexp(1.0,-108);
  ExactNonnegativeSum reversed;
  for(double value:{std::ldexp(1.0,-108),std::ldexp(1.0,-54),1.0}) {
    pass=pass&&Transport::AddExactNonnegative(reversed,value,&next).ok();reversed=next;
  }
  pass=pass&&reversed.words==total.words;
  ExactNonnegativeSum maximum;pass=pass&&Transport::AddExactNonnegative(maximum,std::numeric_limits<double>::max(),&maximum).ok();
  const auto saved=total;
  pass=pass&&!Transport::AddExactNonnegative(maximum,std::numeric_limits<double>::denorm_min(),&total).ok()&&total.words==saved.words;
  ExactNonnegativeSum tiny;
  pass=pass&&Transport::AddExactNonnegative(tiny,std::numeric_limits<double>::denorm_min(),&tiny).ok()&&
      Transport::AddExactNonnegative(tiny,std::numeric_limits<double>::denorm_min(),&tiny).ok()&&
      Transport::ApproximateExactNonnegative(tiny,&display).ok()&&display.nearest==2*std::numeric_limits<double>::denorm_min();
  // Independent geometric sum across every subnormal-to-unit exponent:
  // sum_(n=-1074)^0 2^n = 2-2^-1074, hence integer bits0..1074 are ALL one.
  // This stresses every carry boundary without calling the production encoder
  // to construct its own expected integer or floating rounding oracle.
  ExactNonnegativeSum broad;
  for(int exponent=-1074;exponent<=0;++exponent)
    pass=pass&&Transport::AddExactNonnegative(broad,std::ldexp(1.0,exponent),&broad).ok();
  for(std::size_t word=0;word<broad.words.size();++word) {
    const std::uint64_t expected=word<16?UINT64_MAX:word==16?(UINT64_C(1)<<51)-1:0;
    pass=pass&&broad.words[word]==expected;
  }
  pass=pass&&Transport::ApproximateExactNonnegative(broad,&display).ok()&&display.nearest==2&&
      display.signedRemainder==-std::numeric_limits<double>::denorm_min()&&display.absoluteRoundingBound==0;
  // Binary64 ties-to-even must work on both parity sides of a unit-scale ULP.
  ExactNonnegativeSum tie;
  pass=pass&&Transport::AddExactNonnegative(tie,1,&tie).ok()&&
      Transport::AddExactNonnegative(tie,std::ldexp(1.0,-53),&tie).ok()&&
      Transport::ApproximateExactNonnegative(tie,&display).ok()&&display.nearest==1&&
      display.signedRemainder==std::ldexp(1.0,-53)&&display.absoluteRoundingBound==0;
  pass=pass&&Transport::AddExactNonnegative(tie,std::ldexp(1.0,-52),&tie).ok()&&
      Transport::ApproximateExactNonnegative(tie,&display).ok()&&display.nearest==1+std::ldexp(1.0,-51)&&
      display.signedRemainder==-std::ldexp(1.0,-53)&&display.absoluteRoundingBound==0;
  const auto unchanged=tie;
  pass=pass&&!Transport::AddExactNonnegative(tie,-1,&tie).ok()&&tie.words==unchanged.words&&
      !Transport::AddExactNonnegative(tie,std::numeric_limits<double>::infinity(),&tie).ok()&&tie.words==unchanged.words;
  return Finish(pass,"exact binary64 external account retains multiscale increments/subnormals across1075 exponents independently of order; ties-to-even diagnostics and bounded residuals are explicit, invalid/overflow additions reject transactionally");
}
Descriptor Make(const char* id,const char* name,TestCallback callback) {
  Descriptor d;d.id=id;d.name=name;d.description=name;d.group="polarized-wkb";d.initialization=InitializationLevel::None;
  d.runtime=RuntimeClass::Routine;d.seedCapability=SeedCapability::Deterministic;d.repeatSupported=true;
  d.seedPolicy="deterministic explicit SI ideal-MHD manufactured intervals";
  d.supportedBuildModes="identical native registry and dependency-light sanitizer callback";
  d.stateIsolation="private immutable snapshots, boundary inputs, staged candidates and external accounts";d.callback=callback;return d;
}
} // namespace
std::vector<Descriptor> PolarizedWkbDescriptors() {
  return {Make("WKB01","Signed incoming four-spectrum reservoir",IncomingBoundary),
          Make("WKB02","Exact ideal wave work/ray escape accounting",ExactRayAndExternalAccount),
          Make("WKB03","Independent compact-spectrum WKB convergence",RayConvergence),
          Make("WKB04","Joint spatial/ray ideal-MHD characteristic convergence",StaticGradientConvergence),
          Make("WKB05","ALE geometric conservation of all four spectra",AleGeometricConservation),
          Make("WKB06","Atomic rejection of invalid WKB intervals",RejectedIntervals),
          Make("WKB07","Material side-wall pressure-work conservation",MaterialSidePressure),
          Make("WKB08","Exact multi-scale external spectral account",ExactExternalArithmetic)};
}
} }
