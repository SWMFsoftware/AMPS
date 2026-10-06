#include "sep_corona_swcme/piston_solver.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace SEP { namespace CoronaSwcme { namespace {

struct Derived {
  std::vector<double> density,pressure,magnetic1,magnetic2,radialMagnetic,
      totalPressure,
      viscosity,signal;
};

struct Rate {
  std::vector<double> x,u,e,b1,b2;
  double pistonPower=0;
  double outerPower=0;
  double bodyForcePower=0;
  double heatingPower=0;
  double magneticSourcePower=0;
  double radialMagneticBoundaryPower=0;
};

template<class T> Core::Result<T> Numerical(const std::string& message) {
  return Core::Result<T>::Failure(Core::StatusCode::NumericalFailure,message);
}

double FrameSpeed(const PlanarPistonInput& in) {
  // A uniform Cartesian translation is a true equilibrium.  A constant radial
  // spherical speed is divergent and must not be removed as a frame motion.
  return in.geometry==PistonTubeGeometry::Planar?
      in.initialVelocityMPerS:0;
}

double FaceArea(const PlanarPistonInput& in,double position) {
  return in.geometry==PistonTubeGeometry::Planar?in.areaM2:
      in.solidAngleSr*position*position;
}

double CellVolume(const PlanarPistonInput& in,double left,double right) {
  return in.geometry==PistonTubeGeometry::Planar?in.areaM2*(right-left):
      in.solidAngleSr*(right*right*right-left*left*left)/3;
}

double CellRadius(const PlanarPistonInput& in,double left,double right) {
  if(in.geometry==PistonTubeGeometry::Planar)return 0.5*(left+right);
  // Volume centroid of r over a spherical sector.  The quotient remains well
  // conditioned for qualified positive shell widths at current resolutions.
  return 0.75*(std::pow(right,4)-std::pow(left,4))/
      (std::pow(right,3)-std::pow(left,3));
}

double PhysicalPosition(const PlanarPistonInput& in,double stored,double time) {
  return stored+FrameSpeed(in)*(time-in.startS);
}

PistonHistory DefaultOuterBoundary(const PlanarPistonInput& in) {
  return [in](double time) {
    if(!std::isfinite(time))return Core::Result<CoronalCME::KinematicValue>::
        Failure(Core::StatusCode::OutOfDomain,
            "outer material boundary received nonfinite time");
    return Core::Result<CoronalCME::KinematicValue>::Success({
        in.leftPositionM+in.columnLengthM+
            in.initialVelocityMPerS*(time-in.startS),
        in.initialVelocityMPerS,0});
  };
}

Core::Result<Derived> Derive(const PlanarPistonInput& in,
    const std::vector<double>& mass,const std::vector<double>& x,
    const std::vector<double>& u,const std::vector<double>& e,
    const std::array<std::vector<double>,2>& transverseInvariant,
    const std::vector<double>& radialFluxInvariant) {
  const int n=static_cast<int>(mass.size());
  Derived d;
  d.density.resize(n);d.pressure.resize(n);d.magnetic1.resize(n);
  d.magnetic2.resize(n);
  d.radialMagnetic.resize(n);
  d.totalPressure.resize(n);d.viscosity.resize(n);d.signal.resize(n);
  for(int i=0;i<n;++i) {
    const double width=x[i+1]-x[i];
    if(!(std::isfinite(width)&&width>0&&std::isfinite(e[i])&&e[i]>0))
      return Numerical<Derived>("piston tube cell "+std::to_string(i)+
          " has width="+std::to_string(width)+" m and e="+
          std::to_string(e[i])+" J/kg");
    const double volume=CellVolume(in,x[i],x[i+1]);
    if(!(std::isfinite(volume)&&volume>0))return Numerical<Derived>(
        "piston tube produced a nonpositive physical cell volume");
    d.density[i]=mass[i]/volume;
    d.pressure[i]=(in.gammaAdiabatic-1)*d.density[i]*e[i];
    if(!(std::isfinite(d.density[i])&&d.density[i]>0&&
        std::isfinite(d.pressure[i])&&d.pressure[i]>0))
      return Numerical<Derived>("planar piston produced nonphysical thermodynamics");
    const double radialFactor=in.geometry==PistonTubeGeometry::Planar?1:
        CellRadius(in,x[i],x[i+1]);
    d.magnetic1[i]=transverseInvariant[0][i]*d.density[i]*radialFactor;
    d.magnetic2[i]=transverseInvariant[1][i]*d.density[i]*radialFactor;
    d.radialMagnetic[i]=in.geometry==PistonTubeGeometry::Planar?
        radialFluxInvariant[i]:radialFluxInvariant[i]/(radialFactor*radialFactor);
    const double magneticSquared=d.magnetic1[i]*d.magnetic1[i]+
        d.magnetic2[i]*d.magnetic2[i];
    const double magneticPressure=magneticSquared/
        (2*in.magneticPermeabilityNPerA2);
    d.totalPressure[i]=d.pressure[i]+magneticPressure;
    d.signal[i]=std::sqrt((in.gammaAdiabatic*d.pressure[i]+
        magneticSquared/in.magneticPermeabilityNPerA2)/
        d.density[i]);
    if(!std::isfinite(magneticSquared)||
        !std::isfinite(d.radialMagnetic[i])||
        !std::isfinite(d.totalPressure[i])||
        !(std::isfinite(d.signal[i])&&d.signal[i]>0))
      return Numerical<Derived>(
          "piston tube produced nonphysical magnetic pressure or signal speed");
    const double jump=u[i+1]-u[i];
    const bool linearActive=jump<0&&-jump/d.signal[i]>=
        in.linearViscosityActivation;
    d.viscosity[i]=jump<0?d.density[i]*(in.quadraticViscosity*jump*jump+
        (linearActive?in.linearViscosity*d.signal[i]*std::abs(jump):0)):0;
  }
  return Core::Result<Derived>::Success(std::move(d));
}

Core::Result<Rate> EvaluateRate(const PlanarPistonInput& in,
    const PistonHistory& piston,const PistonHistory& outerBoundary,
    const PistonSourceHistory& source,double time,
    const std::vector<double>& mass,const std::vector<double>& x,
    const std::vector<double>& u,const std::vector<double>& e,
    const std::array<std::vector<double>,2>& transverseInvariant,
    const std::vector<double>& radialFluxInvariant) {
  using Return=Core::Result<Rate>;
  const auto boundary=piston(time);
  if(!boundary.ok())return Return::Failure(
      boundary.status.code,boundary.status.message);
  const auto outer=outerBoundary(time);
  if(!outer.ok())return Return::Failure(outer.status.code,outer.status.message);
  const auto d=Derive(in,mass,x,u,e,transverseInvariant,radialFluxInvariant);
  if(!d.ok())return Return::Failure(d.status.code,d.status.message);
  Rate r;
  r.x=u;
  for(double& velocity:r.x)velocity-=FrameSpeed(in);
  const int cells=static_cast<int>(mass.size());
  r.u.assign(cells+1,0);
  r.e.resize(cells);
  r.b1.assign(cells,0);
  r.b2.assign(cells,0);
  r.x.front()=boundary.value.firstDerivative-FrameSpeed(in);
  r.u.front()=boundary.value.secondDerivative;
  r.x.back()=outer.value.firstDerivative-FrameSpeed(in);
  r.u.back()=outer.value.secondDerivative;

  // The nodal mass is the half-sum of adjacent immutable cell masses.  With
  // the same pressure-plus-Q on the cell work and nodal force, all interior
  // mechanical work cancels pairwise.  Q therefore transfers kinetic energy
  // to heat instead of acting as an unledgered source.
  for(int node=1;node<cells;++node) {
    const double nodalMass=0.5*(mass[node-1]+mass[node]);
    const double left=d.value.totalPressure[node-1]+d.value.viscosity[node-1];
    const double right=d.value.totalPressure[node]+d.value.viscosity[node];
    r.u[node]=FaceArea(in,x[node])*(left-right)/nodalMass;
    if(in.geometry==PistonTubeGeometry::RadialSpherical) {
      double acceleration=0;
      for(int cell:{node-1,node}) {
        const double magneticSquared=d.value.magnetic1[cell]*
            d.value.magnetic1[cell]+d.value.magnetic2[cell]*
            d.value.magnetic2[cell];
        acceleration+=magneticSquared/(in.magneticPermeabilityNPerA2*
          d.value.density[cell]*CellRadius(in,x[cell],x[cell+1]));
      }
      r.u[node]-=0.5*acceleration;
    }
  }

  // Source accelerations act on the same staggered nodal masses as the
  // pressure force.  Sampling at physical nodes gives a resolution-independent
  // Eulerian source law; exact cancellation is expected only for profiles
  // whose discrete pressure gradient represents the continuum derivative.
  // More general maintained atmospheres must demonstrate convergence rather
  // than manufacturing a mesh-dependent cancelling force.
  std::vector<PistonVolumeSource> nodeSource(cells+1),cellSource(cells);
  if(source) {
    for(int node=0;node<=cells;++node) {
      const auto value=source(PhysicalPosition(in,x[node],time),time);
      if(!value.ok())return Return::Failure(
          value.status.code,value.status.message);
      nodeSource[node]=value.value;
      if(!std::isfinite(nodeSource[node].radialAccelerationMPerS2))
        return Numerical<Rate>("piston source returned nonfinite acceleration");
    }
    for(int cell=0;cell<cells;++cell) {
      const double radius=PhysicalPosition(in,
          CellRadius(in,x[cell],x[cell+1]),time);
      const auto value=source(radius,time);
      if(!value.ok())return Return::Failure(
          value.status.code,value.status.message);
      cellSource[cell]=value.value;
      if(!std::isfinite(cellSource[cell].heatingWPerM3)||
          !std::isfinite(cellSource[cell].transverseInvariantRate[0])||
          !std::isfinite(cellSource[cell].transverseInvariantRate[1]))
        return Numerical<Rate>("piston source returned a nonfinite cell term");
    }
    for(int node=1;node<cells;++node) {
      const double nodalMass=0.5*(mass[node-1]+mass[node]);
      r.u[node]+=nodeSource[node].radialAccelerationMPerS2;
      r.bodyForcePower+=nodalMass*
          nodeSource[node].radialAccelerationMPerS2*u[node];
    }
  }
  for(int cell=0;cell<cells;++cell) {
    const double totalPressure=d.value.pressure[cell]+d.value.viscosity[cell];
    const double volumeRate=FaceArea(in,x[cell+1])*u[cell+1]-
        FaceArea(in,x[cell])*u[cell];
    r.e[cell]=-totalPressure*volumeRate/mass[cell];
    if(source) {
      const double volume=CellVolume(in,x[cell],x[cell+1]);
      r.e[cell]+=cellSource[cell].heatingWPerM3*volume/mass[cell];
      r.b1[cell]=cellSource[cell].transverseInvariantRate[0];
      r.b2[cell]=cellSource[cell].transverseInvariantRate[1];
      r.heatingPower+=cellSource[cell].heatingWPerM3*volume;
      const double radialFactor=in.geometry==PistonTubeGeometry::Planar?1:
          CellRadius(in,x[cell],x[cell+1]);
      const double factor=d.value.density[cell]*radialFactor;
      r.magneticSourcePower+=volume/in.magneticPermeabilityNPerA2*
          (d.value.magnetic1[cell]*factor*r.b1[cell]+
           d.value.magnetic2[cell]*factor*r.b2[cell]);
    }
  }

  const double firstPressure=d.value.totalPressure.front()+d.value.viscosity.front();
  const double lastPressure=d.value.totalPressure.back()+d.value.viscosity.back();
  const double leftNodeMass=0.5*mass.front();
  const double rightNodeMass=0.5*mass.back();
  // Endpoint half masses are part of the staggered kinetic quadrature.  Their
  // actuator inertia therefore belongs in the external-work ledger.  It
  // vanishes with refinement and is exactly zero for constant-speed tests.
  const double leftSource=source?nodeSource.front().radialAccelerationMPerS2:0;
  const double rightSource=source?nodeSource.back().radialAccelerationMPerS2:0;
  // The actuator supplies only the force not already provided by the volume
  // source.  Source work on endpoint half masses is ledgered separately, so
  // pressure + actuator + source exactly recover d(kinetic energy)/dt.
  r.pistonPower=(FaceArea(in,x.front())*firstPressure+
      leftNodeMass*(boundary.value.secondDerivative-leftSource))*
      boundary.value.firstDerivative;
  r.outerPower=(-FaceArea(in,x.back())*lastPressure+
      rightNodeMass*(outer.value.secondDerivative-rightSource))*
      outer.value.firstDerivative;
  r.bodyForcePower+=leftNodeMass*leftSource*u.front()+
      rightNodeMass*rightSource*u.back();
  // A material surface threaded by B_r is not a closed magnetic-energy
  // surface.  The radial Maxwell stress contributes -B_r^2/(2 mu0) to the
  // material-boundary work even though pressure and tension cancel from the
  // bulk radial Lorentz force.  Keeping it separate avoids treating passive
  // B_r as a spurious compression pressure.
  r.radialMagneticBoundaryPower=
      -FaceArea(in,x.front())*d.value.radialMagnetic.front()*
          d.value.radialMagnetic.front()/(2*in.magneticPermeabilityNPerA2)*
          boundary.value.firstDerivative+
      FaceArea(in,x.back())*d.value.radialMagnetic.back()*
          d.value.radialMagnetic.back()/(2*in.magneticPermeabilityNPerA2)*
          outer.value.firstDerivative;
  return Return::Success(std::move(r));
}

double Energy(const PlanarPistonInput& in,const std::vector<double>& mass,
    const std::vector<double>& x,const std::vector<double>& u,
    const std::vector<double>& e,
    const std::array<std::vector<double>,2>& transverseInvariant,
    const std::vector<double>& radialFluxInvariant) {
  double energy=0;
  const int cells=static_cast<int>(mass.size());
  for(int i=0;i<cells;++i) {
    energy+=mass[i]*e[i];
    const double volume=CellVolume(in,x[i],x[i+1]);
    const double density=mass[i]/volume;
    const double radialFactor=in.geometry==PistonTubeGeometry::Planar?1:
        CellRadius(in,x[i],x[i+1]);
    const double magnetic1=transverseInvariant[0][i]*density*radialFactor;
    const double magnetic2=transverseInvariant[1][i]*density*radialFactor;
    const double radialMagnetic=in.geometry==PistonTubeGeometry::Planar?
        radialFluxInvariant[i]:radialFluxInvariant[i]/(radialFactor*radialFactor);
    energy+=(magnetic1*magnetic1+magnetic2*magnetic2+
        radialMagnetic*radialMagnetic)/
        (2*in.magneticPermeabilityNPerA2)*
        volume;
  }
  energy+=0.25*mass.front()*u.front()*u.front();
  for(int i=1;i<cells;++i)
    energy+=0.25*(mass[i-1]+mass[i])*u[i]*u[i];
  energy+=0.25*mass.back()*u.back()*u.back();
  return energy;
}

} // namespace

Core::Result<std::unique_ptr<PlanarPistonSolver>> PlanarPistonSolver::Create(
    PlanarPistonInput input,PistonHistory piston,PistonSourceHistory source,
    PistonHistory outerBoundary) {
  PistonInitialState state;
  if(input.cells>0) {
    state.nodeVelocityMPerS.assign(input.cells+1,input.initialVelocityMPerS);
    state.cellDensityKgM3.assign(input.cells,input.initialDensityKgM3);
    state.cellPressurePa.assign(input.cells,input.initialPressurePa);
    state.cellTransverseMagneticFieldT.assign(input.cells,
        input.initialTransverseMagneticFieldT);
    state.cellTransverseMagneticField2T.assign(input.cells,
        input.initialTransverseMagneticField2T);
    state.cellRadialMagneticFieldT.assign(input.cells,
        input.initialRadialMagneticFieldT);
  }
  if(piston) {
    const auto boundary=piston(input.startS);
    if(boundary.ok()&&!state.nodeVelocityMPerS.empty())
      state.nodeVelocityMPerS.front()=boundary.value.firstDerivative;
  }
  return CreateInitialized(input,std::move(piston),std::move(state),
      std::move(source),std::move(outerBoundary));
}

Core::Result<std::unique_ptr<PlanarPistonSolver>>
PlanarPistonSolver::CreateInitialized(PlanarPistonInput input,
    PistonHistory piston,PistonInitialState state,PistonSourceHistory source,
    PistonHistory outerBoundary) {
  using Return=Core::Result<std::unique_ptr<PlanarPistonSolver>>;
  // Existing one-component manufactured fixtures lie in a declared tangent
  // plane.  Supplying no second component means physical zero, not missing
  // state; production projected-IMF initialization supplies both explicitly.
  if(state.cellTransverseMagneticField2T.empty()&&input.cells>0)
    state.cellTransverseMagneticField2T.assign(input.cells,0);
  if(state.cellRadialMagneticFieldT.empty()&&input.cells>0)
    state.cellRadialMagneticFieldT.assign(input.cells,0);
  const bool validMeasure=input.geometry==PistonTubeGeometry::Planar?
      input.areaM2>0:input.solidAngleSr>0&&input.leftPositionM>0;
  if(!piston||!std::isfinite(input.startS)||!std::isfinite(input.endS)||
      !(input.endS>input.startS)||!(input.columnLengthM>0)||
      !validMeasure||!(input.initialDensityKgM3>0)||
      !(input.initialPressurePa>0)||!std::isfinite(input.initialVelocityMPerS)||
      !std::isfinite(input.initialTransverseMagneticFieldT)||
      !std::isfinite(input.initialTransverseMagneticField2T)||
      !std::isfinite(input.initialRadialMagneticFieldT)||
      !(input.magneticPermeabilityNPerA2>0)||
      !(input.gammaAdiabatic>1)||input.cells<32||!(input.cfl>0&&input.cfl<=0.8)||
      !(input.quadraticViscosity>=0)||!(input.linearViscosity>=0)||
      !(input.linearViscosityActivation>=0&&
        input.linearViscosityActivation<1)||
      !(input.shockThreshold>0))return Return::Failure(
          Core::StatusCode::InvalidConfiguration,
      "planar piston input is incomplete or outside stable bounds");
  if(state.nodeVelocityMPerS.size()!=static_cast<std::size_t>(input.cells+1)||
      state.cellDensityKgM3.size()!=static_cast<std::size_t>(input.cells)||
      state.cellPressurePa.size()!=static_cast<std::size_t>(input.cells)||
      state.cellTransverseMagneticFieldT.size()!=
          static_cast<std::size_t>(input.cells)||
      state.cellTransverseMagneticField2T.size()!=
          static_cast<std::size_t>(input.cells)||
      state.cellRadialMagneticFieldT.size()!=
          static_cast<std::size_t>(input.cells))return Return::Failure(
              Core::StatusCode::InvalidConfiguration,
              "initialized piston state dimensions do not match the tube");
  const auto initial=piston(input.startS);
  if(!initial.ok())return Return::Failure(initial.status.code,initial.status.message);
  if(!std::isfinite(initial.value.value)||
      std::abs(initial.value.value-input.leftPositionM)>
          64*std::numeric_limits<double>::epsilon()*
          std::max(1.0,std::abs(input.leftPositionM)))return Return::Failure(
              Core::StatusCode::InvalidConfiguration,
              "piston initial position does not match the material column");
  if(!outerBoundary)outerBoundary=DefaultOuterBoundary(input);
  const auto initialOuter=outerBoundary(input.startS);
  if(!initialOuter.ok())return Return::Failure(
      initialOuter.status.code,initialOuter.status.message);
  const double expectedOuter=input.leftPositionM+input.columnLengthM;
  if(!std::isfinite(initialOuter.value.value)||
      std::abs(initialOuter.value.value-expectedOuter)>
          64*std::numeric_limits<double>::epsilon()*
          std::max(1.0,std::abs(expectedOuter)))return Return::Failure(
              Core::StatusCode::InvalidConfiguration,
              "outer boundary initial position does not match the material column");
  std::unique_ptr<PlanarPistonSolver> out(new PlanarPistonSolver);
  out->input_=input;out->piston_=std::move(piston);
  out->outerBoundary_=std::move(outerBoundary);
  out->source_=std::move(source);out->timeS_=input.startS;
  out->cellMass_.resize(input.cells);
  out->x_.resize(input.cells+1);
  out->u_=std::move(state.nodeVelocityMPerS);
  out->e_.resize(input.cells);
  for(int i=0;i<=input.cells;++i)
    out->x_[i]=input.leftPositionM+input.columnLengthM*i/input.cells;
  out->transverseInvariant_[0].resize(input.cells);
  out->transverseInvariant_[1].resize(input.cells);
  out->radialFluxInvariant_.resize(input.cells);
  for(int i=0;i<input.cells;++i) {
    if(!(std::isfinite(state.cellDensityKgM3[i])&&
        state.cellDensityKgM3[i]>0&&std::isfinite(state.cellPressurePa[i])&&
        state.cellPressurePa[i]>0&&
        std::isfinite(state.cellTransverseMagneticFieldT[i])&&
        std::isfinite(state.cellTransverseMagneticField2T[i])&&
        std::isfinite(state.cellRadialMagneticFieldT[i])))
      return Return::Failure(Core::StatusCode::InvalidConfiguration,
          "initialized piston cell state is nonphysical");
    out->cellMass_[i]=state.cellDensityKgM3[i]*
        CellVolume(input,out->x_[i],out->x_[i+1]);
    out->e_[i]=state.cellPressurePa[i]/
        ((input.gammaAdiabatic-1)*state.cellDensityKgM3[i]);
    const double radialFactor=input.geometry==PistonTubeGeometry::Planar?1:
        CellRadius(input,out->x_[i],out->x_[i+1]);
    out->transverseInvariant_[0][i]=state.cellTransverseMagneticFieldT[i]/
        (state.cellDensityKgM3[i]*radialFactor);
    out->transverseInvariant_[1][i]=state.cellTransverseMagneticField2T[i]/
        (state.cellDensityKgM3[i]*radialFactor);
    out->radialFluxInvariant_[i]=state.cellRadialMagneticFieldT[i]*
        (input.geometry==PistonTubeGeometry::Planar?1:radialFactor*radialFactor);
  }
  for(double velocity:out->u_)if(!std::isfinite(velocity))return Return::Failure(
      Core::StatusCode::InvalidConfiguration,
      "initialized piston node velocity is nonfinite");
  if(std::abs(out->u_.front()-initial.value.firstDerivative)>
      1e-12*std::max(1.0,std::abs(initial.value.firstDerivative))||
      std::abs(out->u_.back()-initialOuter.value.firstDerivative)>
      1e-12*std::max(1.0,std::abs(initialOuter.value.firstDerivative)))
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "initialized node velocities violate piston or outer boundary data");
  out->initialEnergyJ_=Energy(input,out->cellMass_,out->x_,out->u_,out->e_,
      out->transverseInvariant_,out->radialFluxInvariant_);
  return Return::Success(std::move(out));
}

Core::Status PlanarPistonSolver::AdvanceTo(double target) {
  const auto fail=[](Core::StatusCode code,const std::string& message) {
    return Core::Status::Failure(code,message);
  };
  if(!std::isfinite(target)||target<timeS_||target>input_.endS)
    return fail(Core::StatusCode::OutOfDomain,
        "planar piston target time is outside forward support");
  auto x=x_,u=u_,e=e_;
  auto b=transverseInvariant_;
  double time=timeS_,pistonWork=pistonWorkJ_,outerWork=outerBoundaryWorkJ_;
  double bodyWork=bodyForceWorkJ_,heating=volumeHeatingJ_,
      magneticWork=magneticSourceWorkJ_;
  double radialBoundaryWork=radialMagneticBoundaryWorkJ_;
  while(time<target) {
    const auto d=Derive(input_,cellMass_,x,u,e,b,radialFluxInvariant_);
    if(!d.ok())return d.status;
    double dt=target-time;
    for(std::size_t i=0;i<cellMass_.size();++i) {
      const double width=x[i+1]-x[i];
      dt=std::min(dt,input_.cfl*width/
          (d.value.signal[i]+std::abs(u[i+1]-u[i])+1e-300));
      const double volumeRate=FaceArea(input_,x[i+1])*u[i+1]-
          FaceArea(input_,x[i])*u[i];
      const double thermalRate=-(d.value.pressure[i]+d.value.viscosity[i])*
          volumeRate/cellMass_[i];
      if(thermalRate<0)dt=std::min(dt,0.5*e[i]/(-thermalRate));
    }
    if(!(std::isfinite(dt)&&dt>0))return fail(
        Core::StatusCode::NumericalFailure,"planar piston CFL step is invalid");
    if(dt<1e-12*std::max(1.0,std::abs(time)))return fail(
        Core::StatusCode::NumericalFailure,
        "piston tube time step collapsed at t="+std::to_string(time)+
        " s with dt="+std::to_string(dt)+" s");
    const auto k1=EvaluateRate(input_,piston_,outerBoundary_,source_,time,
        cellMass_,x,u,e,b,radialFluxInvariant_);
    if(!k1.ok())return k1.status;
    const auto stage=[&](const Rate& rate,double fraction,
        std::vector<double>* xs,std::vector<double>* us,
        std::vector<double>* es,std::array<std::vector<double>,2>* bs)->
        Core::Status {
      *xs=x;*us=u;*es=e;
      *bs=b;
      for(std::size_t i=0;i<xs->size();++i) {
        (*xs)[i]+=fraction*dt*rate.x[i];
        (*us)[i]+=fraction*dt*rate.u[i];
      }
      for(std::size_t i=0;i<es->size();++i)
        (*es)[i]+=fraction*dt*rate.e[i];
      for(std::size_t i=0;i<es->size();++i) {
        (*bs)[0][i]+=fraction*dt*rate.b1[i];
        (*bs)[1][i]+=fraction*dt*rate.b2[i];
      }
      const double stageTime=time+fraction*dt;
      const auto boundary=piston_(stageTime);
      if(!boundary.ok())return boundary.status;
      const auto outer=outerBoundary_(stageTime);
      if(!outer.ok())return outer.status;
      const double frame=FrameSpeed(input_)*
          (stageTime-input_.startS);
      xs->front()=boundary.value.value-frame;
      us->front()=boundary.value.firstDerivative;
      xs->back()=outer.value.value-frame;
      us->back()=outer.value.firstDerivative;
      return Core::Status::Success();
    };

    // Classical RK4 is used rather than two-stage explicit RK.  In the
    // linear, viscosity-free limit the centered Lagrangian acoustic operator
    // has imaginary eigenvalues; RK2 amplifies them slightly on every step,
    // producing a refinement-worsening mode.  RK4 has a finite imaginary-axis
    // stability interval and preserves the deliberately undamped linear-wave
    // acceptance case at the selected CFL.
    std::vector<double> xs,us,es;
    std::array<std::vector<double>,2> bs;
    auto staged=stage(k1.value,0.5,&xs,&us,&es,&bs);
    if(!staged.ok())return staged;
    const auto k2=EvaluateRate(input_,piston_,outerBoundary_,source_,
        time+0.5*dt,
        cellMass_,xs,us,es,bs,radialFluxInvariant_);
    if(!k2.ok())return k2.status;
    staged=stage(k2.value,0.5,&xs,&us,&es,&bs);
    if(!staged.ok())return staged;
    const auto k3=EvaluateRate(input_,piston_,outerBoundary_,source_,
        time+0.5*dt,
        cellMass_,xs,us,es,bs,radialFluxInvariant_);
    if(!k3.ok())return k3.status;
    staged=stage(k3.value,1,&xs,&us,&es,&bs);
    if(!staged.ok())return staged;
    const auto k4=EvaluateRate(input_,piston_,outerBoundary_,source_,time+dt,
        cellMass_,xs,us,es,bs,radialFluxInvariant_);
    if(!k4.ok())return k4.status;
    for(std::size_t i=0;i<x.size();++i) {
      x[i]+=dt*(k1.value.x[i]+2*k2.value.x[i]+2*k3.value.x[i]+
          k4.value.x[i])/6;
      u[i]+=dt*(k1.value.u[i]+2*k2.value.u[i]+2*k3.value.u[i]+
          k4.value.u[i])/6;
    }
    for(std::size_t i=0;i<e.size();++i)e[i]+=dt*(k1.value.e[i]+
        2*k2.value.e[i]+2*k3.value.e[i]+k4.value.e[i])/6;
    for(std::size_t i=0;i<e.size();++i) {
      b[0][i]+=dt*(k1.value.b1[i]+2*k2.value.b1[i]+2*k3.value.b1[i]+
          k4.value.b1[i])/6;
      b[1][i]+=dt*(k1.value.b2[i]+2*k2.value.b2[i]+2*k3.value.b2[i]+
          k4.value.b2[i])/6;
    }
    const auto boundary=piston_(time+dt);
    if(!boundary.ok())return boundary.status;
    const auto outer=outerBoundary_(time+dt);
    if(!outer.ok())return outer.status;
    const double frame=FrameSpeed(input_)*
        (time+dt-input_.startS);
    x.front()=boundary.value.value-frame;u.front()=boundary.value.firstDerivative;
    x.back()=outer.value.value-frame;
    u.back()=outer.value.firstDerivative;
    pistonWork+=dt*(k1.value.pistonPower+2*k2.value.pistonPower+
        2*k3.value.pistonPower+k4.value.pistonPower)/6;
    outerWork+=dt*(k1.value.outerPower+2*k2.value.outerPower+
        2*k3.value.outerPower+k4.value.outerPower)/6;
    bodyWork+=dt*(k1.value.bodyForcePower+2*k2.value.bodyForcePower+
        2*k3.value.bodyForcePower+k4.value.bodyForcePower)/6;
    heating+=dt*(k1.value.heatingPower+2*k2.value.heatingPower+
        2*k3.value.heatingPower+k4.value.heatingPower)/6;
    magneticWork+=dt*(k1.value.magneticSourcePower+
        2*k2.value.magneticSourcePower+2*k3.value.magneticSourcePower+
        k4.value.magneticSourcePower)/6;
    radialBoundaryWork+=dt*(k1.value.radialMagneticBoundaryPower+
        2*k2.value.radialMagneticBoundaryPower+
        2*k3.value.radialMagneticBoundaryPower+
        k4.value.radialMagneticBoundaryPower)/6;
    time+=dt;
    const auto physical=Derive(input_,cellMass_,x,u,e,b,radialFluxInvariant_);
    if(!physical.ok())return physical.status;
  }
  x_=std::move(x);u_=std::move(u);e_=std::move(e);
  transverseInvariant_=std::move(b);timeS_=time;
  pistonWorkJ_=pistonWork;outerBoundaryWorkJ_=outerWork;
  bodyForceWorkJ_=bodyWork;volumeHeatingJ_=heating;
  magneticSourceWorkJ_=magneticWork;
  radialMagneticBoundaryWorkJ_=radialBoundaryWork;
  return Core::Status::Success();
}

Core::Status PlanarPistonSolver::AppendAmbient(PistonAppendState appended) {
  const auto fail=[](Core::StatusCode code,const std::string& message) {
    return Core::Status::Failure(code,message);
  };
  const std::size_t added=appended.nodePositionM.empty()?0:
      appended.nodePositionM.size()-1;
  if(added==0||!appended.outerBoundary||
      appended.state.nodeVelocityMPerS.size()!=added+1||
      appended.state.cellDensityKgM3.size()!=added||
      appended.state.cellPressurePa.size()!=added||
      appended.state.cellTransverseMagneticFieldT.size()!=added)
    return fail(Core::StatusCode::InvalidConfiguration,
        "ambient append dimensions or outer history are incomplete");
  if(appended.state.cellTransverseMagneticField2T.empty())
    appended.state.cellTransverseMagneticField2T.assign(added,0);
  if(appended.state.cellRadialMagneticFieldT.empty())
    appended.state.cellRadialMagneticFieldT.assign(added,0);
  if(appended.state.cellTransverseMagneticField2T.size()!=added||
      appended.state.cellRadialMagneticFieldT.size()!=added)
    return fail(Core::StatusCode::InvalidConfiguration,
        "ambient append magnetic dimensions are incomplete");

  const double frame=FrameSpeed(input_)*(timeS_-input_.startS);
  const auto close=[](double a,double b) {
    return std::abs(a-b)<=128*std::numeric_limits<double>::epsilon()*
        std::max({1.0,std::abs(a),std::abs(b)});
  };
  if(!close(appended.nodePositionM.front(),x_.back()+frame)||
      !close(appended.state.nodeVelocityMPerS.front(),u_.back()))
    return fail(Core::StatusCode::InvalidConfiguration,
        "ambient append does not share the committed outer material node");
  const auto newOuter=appended.outerBoundary(timeS_);
  if(!newOuter.ok())return newOuter.status;
  if(!close(newOuter.value.value,appended.nodePositionM.back())||
      !close(newOuter.value.firstDerivative,
          appended.state.nodeVelocityMPerS.back())||
      !std::isfinite(newOuter.value.secondDerivative))
    return fail(Core::StatusCode::InvalidConfiguration,
        "ambient append outer history does not match its last node");

  auto mass=cellMass_,x=x_,u=u_,e=e_;
  auto b=transverseInvariant_;
  auto radial=radialFluxInvariant_;
  const double energyBefore=Energy(input_,mass,x,u,e,b,radial);
  for(std::size_t node=1;node<appended.nodePositionM.size();++node) {
    const double stored=appended.nodePositionM[node]-frame;
    if(!std::isfinite(stored)||!(stored>x.back())||
        !std::isfinite(appended.state.nodeVelocityMPerS[node]))
      return fail(Core::StatusCode::InvalidConfiguration,
          "ambient append nodes are nonfinite or not strictly outward");
    x.push_back(stored);u.push_back(appended.state.nodeVelocityMPerS[node]);
  }
  for(std::size_t cell=0;cell<added;++cell) {
    const double rho=appended.state.cellDensityKgM3[cell];
    const double pressure=appended.state.cellPressurePa[cell];
    const double bt1=appended.state.cellTransverseMagneticFieldT[cell];
    const double bt2=appended.state.cellTransverseMagneticField2T[cell];
    const double br=appended.state.cellRadialMagneticFieldT[cell];
    if(!(std::isfinite(rho)&&rho>0&&std::isfinite(pressure)&&pressure>0&&
        std::isfinite(bt1)&&std::isfinite(bt2)&&std::isfinite(br)))
      return fail(Core::StatusCode::InvalidConfiguration,
          "ambient append contains a nonphysical cell primitive");
    const std::size_t index=cellMass_.size()+cell;
    const double volume=CellVolume(input_,x[index],x[index+1]);
    const double radius=CellRadius(input_,x[index],x[index+1]);
    const double radialFactor=input_.geometry==PistonTubeGeometry::Planar?1:radius;
    mass.push_back(rho*volume);
    e.push_back(pressure/((input_.gammaAdiabatic-1)*rho));
    b[0].push_back(bt1/(rho*radialFactor));
    b[1].push_back(bt2/(rho*radialFactor));
    radial.push_back(br*(input_.geometry==PistonTubeGeometry::Planar?
        1:radius*radius));
  }
  const auto physical=Derive(input_,mass,x,u,e,b,radial);
  if(!physical.ok())return physical.status;
  const double energyAfter=Energy(input_,mass,x,u,e,b,radial);
  if(!(std::isfinite(energyAfter)&&energyAfter>energyBefore))return fail(
      Core::StatusCode::NumericalFailure,
      "ambient append did not add finite positive material energy");

  cellMass_=std::move(mass);x_=std::move(x);u_=std::move(u);e_=std::move(e);
  transverseInvariant_=std::move(b);radialFluxInvariant_=std::move(radial);
  outerBoundary_=std::move(appended.outerBoundary);
  appendedEnergyJ_+=energyAfter-energyBefore;
  return Core::Status::Success();
}

std::unique_ptr<PlanarPistonSolver> PlanarPistonSolver::Clone() const {
  return std::unique_ptr<PlanarPistonSolver>(new PlanarPistonSolver(*this));
}

Core::Result<std::vector<PlanarCellState>> PlanarPistonSolver::Cells() const {
  using Return=Core::Result<std::vector<PlanarCellState>>;
  const auto d=Derive(input_,cellMass_,x_,u_,e_,transverseInvariant_,
      radialFluxInvariant_);
  if(!d.ok())return Return::Failure(d.status.code,d.status.message);
  std::vector<PlanarCellState> out(cellMass_.size());
  const double frame=FrameSpeed(input_)*(timeS_-input_.startS);
  for(std::size_t i=0;i<cellMass_.size();++i)out[i]={
      CellRadius(input_,x_[i],x_[i+1])+frame,
      x_[i+1]-x_[i],cellMass_[i],d.value.density[i],d.value.pressure[i],
      e_[i],d.value.viscosity[i],d.value.magnetic1[i],d.value.magnetic2[i],
      d.value.radialMagnetic[i],d.value.totalPressure[i],
      0.5*(u_[i]+u_[i+1])};
  return Return::Success(std::move(out));
}

std::vector<double> PlanarPistonSolver::NodePositionsM() const {
  std::vector<double> physical=x_;
  const double frame=FrameSpeed(input_)*(timeS_-input_.startS);
  for(double& position:physical)position+=frame;
  return physical;
}

Core::Result<PlanarShockState> PlanarPistonSolver::DetectShock() const {
  using Return=Core::Result<PlanarShockState>;
  const auto cells=Cells();
  if(!cells.ok())return Return::Failure(cells.status.code,cells.status.message);
  int last=-1;
  PlanarShockState out;
  bool insideZone=false;
  for(int i=0;i<static_cast<int>(cellMass_.size());++i) {
    const double ratio=cells.value[i].artificialPressurePa/
        cells.value[i].totalPressurePa;
    out.maximumArtificialPressureRatio=std::max(
        out.maximumArtificialPressureRatio,ratio);
    const bool active=ratio>input_.shockThreshold;
    if(active&&!insideZone)++out.shockZoneCount;
    insideZone=active;
    if(active)last=i;
  }
  if(last<0)return Return::Success(out);
  int first=last;
  while(first>0&&cells.value[first-1].artificialPressurePa/
      cells.value[first-1].totalPressurePa>input_.shockThreshold)--first;
  out.present=true;out.firstShockCell=first;out.lastShockCell=last;
  double weighted=0,weight=0;
  for(int i=first;i<=last;++i) {
    weighted+=cells.value[i].artificialPressurePa*cells.value[i].centerM;
    weight+=cells.value[i].artificialPressurePa;
    out.selectedZoneMaximumArtificialPressureRatio=std::max(
        out.selectedZoneMaximumArtificialPressureRatio,
        cells.value[i].artificialPressurePa/cells.value[i].totalPressurePa);
  }
  out.radiusM=weighted/weight;
  // Shock-zone admission and two-sided state extraction are deliberately
  // separate.  A newly formed shock can be located by Q/P before six-cell
  // plateaus exist; consumers that require RH states must also check
  // statesAvailable and may not fill them from the dissipative zone.
  constexpr int gap=2,window=6;
  if(first-gap-window<0||
      last+gap+window>=static_cast<int>(cellMass_.size()))
    return Return::Success(out);
  auto mean=[&](int begin,int end,auto field) {
    double sum=0;for(int i=begin;i<end;++i)sum+=field(cells.value[i]);
    return sum/(end-begin);
  };
  const int d1=first-gap,d0=d1-window,u0=last+gap+1,u1=u0+window;
  out.downstreamDensityKgM3=mean(d0,d1,[](const auto& s){return s.densityKgM3;});
  out.downstreamPressurePa=mean(d0,d1,[](const auto& s){return s.pressurePa;});
  out.downstreamVelocityMPerS=mean(d0,d1,[](const auto& s){return s.velocityMPerS;});
  out.downstreamTransverseMagneticFieldT=mean(d0,d1,
      [](const auto& s){return s.transverseMagneticFieldT;});
  out.downstreamTransverseMagneticField2T=mean(d0,d1,
      [](const auto& s){return s.transverseMagneticField2T;});
  out.downstreamRadialMagneticFieldT=mean(d0,d1,
      [](const auto& s){return s.radialMagneticFieldT;});
  out.upstreamDensityKgM3=mean(u0,u1,[](const auto& s){return s.densityKgM3;});
  out.upstreamPressurePa=mean(u0,u1,[](const auto& s){return s.pressurePa;});
  out.upstreamVelocityMPerS=mean(u0,u1,[](const auto& s){return s.velocityMPerS;});
  out.upstreamTransverseMagneticFieldT=mean(u0,u1,
      [](const auto& s){return s.transverseMagneticFieldT;});
  out.upstreamTransverseMagneticField2T=mean(u0,u1,
      [](const auto& s){return s.transverseMagneticField2T;});
  out.upstreamRadialMagneticFieldT=mean(u0,u1,
      [](const auto& s){return s.radialMagneticFieldT;});
  const double denominator=out.downstreamDensityKgM3-out.upstreamDensityKgM3;
  if(!(denominator>0))return Return::Success(out);
  out.speedMPerS=(out.downstreamDensityKgM3*out.downstreamVelocityMPerS-
      out.upstreamDensityKgM3*out.upstreamVelocityMPerS)/denominator;
  out.compressionRatio=out.downstreamDensityKgM3/out.upstreamDensityKgM3;
  out.statesAvailable=std::isfinite(out.speedMPerS)&&out.speedMPerS>0;
  return Return::Success(out);
}

Core::Result<PlanarEnergyLedger> PlanarPistonSolver::EnergyLedger() const {
  using Return=Core::Result<PlanarEnergyLedger>;
  const auto physical=Derive(input_,cellMass_,x_,u_,e_,transverseInvariant_,
      radialFluxInvariant_);
  if(!physical.ok())return Return::Failure(
      physical.status.code,physical.status.message);
  PlanarEnergyLedger out;
  out.initialEnergyJ=initialEnergyJ_;
  out.currentEnergyJ=Energy(input_,cellMass_,x_,u_,e_,transverseInvariant_,
      radialFluxInvariant_);
  out.pistonWorkJ=pistonWorkJ_;out.outerBoundaryWorkJ=outerBoundaryWorkJ_;
  out.bodyForceWorkJ=bodyForceWorkJ_;out.volumeHeatingJ=volumeHeatingJ_;
  out.magneticSourceWorkJ=magneticSourceWorkJ_;
  out.radialMagneticBoundaryWorkJ=radialMagneticBoundaryWorkJ_;
  out.appendedEnergyJ=appendedEnergyJ_;
  out.residualJ=out.currentEnergyJ-out.initialEnergyJ-
      out.pistonWorkJ-out.outerBoundaryWorkJ-out.bodyForceWorkJ-
      out.volumeHeatingJ-out.magneticSourceWorkJ-out.appendedEnergyJ;
  out.residualJ-=out.radialMagneticBoundaryWorkJ;
  return Return::Success(out);
}

} } // namespace SEP::CoronaSwcme
