/*
 * mover.cpp
 *
 *  Created on: May 16, 2020
 *      Author: vtenishe
 */
#include <algorithm>
#include <math.h>

#include "sep.h"
#include "amps2swmf.h"

bool SEP::AccountTransportCoefficient=true;
SEP::fParticleMover SEP::ParticleMoverPtr=ParticleMover_FTE;
double SEP::MaxTurbulenceLevel=0.1;
bool SEP::MaxTurbulenceEnforceLimit=false;

//set the lower limit of the mean free path being the local Larmor radius of the particle
bool SEP::LimitMeanFreePath=false;

bool SEP::LimitScatteringUpcomingWave=false;

//set the numerical limit on the number of simulated scattering events
bool SEP::NumericalScatteringEventMode=false;
double SEP::NumericalScatteringEventLimiter=-1.0;

// Apply adiabatic cooling only if the flag is set
bool SEP::AccountAdiabaticCoolingFlag=true;

// Historical field-line implementations remain private implementation material;
// public selection is owned exclusively by the three-entry production registry.
//===================================================================================================================================
int SEP::ParticleMover_Droge_2009_AJ(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  namespace PB = PIC::ParticleBuffer;
  namespace FL = PIC::FieldLine;

  PIC::ParticleBuffer::byte *ParticleData;
  double W[2],mu,AbsB,absB2,vParallel,vNormal,v,DivAbsB,vParallelInit,vNormalInit;
  double FieldLineCoord,Lmax,vAlfven;
  int iFieldLine,spec;
  FL::cFieldLineSegment *Segment;

  static int nCallCnt=0;
  nCallCnt++;

  ParticleData=PB::GetParticleDataPointer(ptr);

  FieldLineCoord=PB::GetFieldLineCoord(ParticleData);
  iFieldLine=PB::GetFieldLineId(ParticleData);
  spec=PB::GetI(ParticleData);

  //velocity is in the frame moving with solar wind
  vParallel=PB::GetVParallel(ParticleData);
  vNormal=PB::GetVNormal(ParticleData);

  vParallelInit=vParallel,vNormalInit=vNormal;

  /*double  ee=Relativistic::Speed2E(sqrt(vNormal*vNormal+vParallel*vParallel),PIC::MolecularData::GetMass(spec));
  ee*=J2MeV;

  if (ee>200) {
    double ss=0.0;

    ss+=23;
  }*/

  //determine the segment of the particle location
  Segment=FL::FieldLinesAll[iFieldLine].GetSegment(FieldLineCoord);

  //double AbsBDeriv;
  double vSolarWind[3],vSolarWindParallel;
  double FieldLineCoord_init=FieldLineCoord;

  //get the new value of 'mu'
  double D,dD_dmu;

  double mu_init=mu;
  double time_counter=0.0;
  double dt=dtTotal;
  double dmu=0.0;
  double delta;

  bool first_pass_flag=true;
  static long int loop_cnt=0;

  SEP::Diffusion::cD_SA D_SA;
  SEP::Diffusion::cD_mu_mu D_mu_mu;
  static SEP::Diffusion::cD_x_x<SEP::Diffusion::cD_mu_mu> D_x_x;
  static SEP::Diffusion::cD_mu_mu_Jokopii1966AJ<100,100> D_mu_mu_Jokopii1966AJ;

  SEP::Diffusion::cD_x_x<SEP::Diffusion::cD_mu_mu> *D_x_x_ptr=&D_x_x;
  SEP::Diffusion::cD_SA *D_SA_ptr=&D_SA;
  SEP::Diffusion::cD_mu_mu *D_mu_mu_TwoWaves_ptr=&D_mu_mu;
  SEP::Diffusion::cDiffusionCoeffcient *D_mu_mu_ptr=&D_mu_mu_Jokopii1966AJ;

  double *B0,*B1,B[3],r2;
  double *W0,*W1;
  double *x0,*x1;
  double w0,w1;
  double PlasmaDensity0,PlasmaDensity1,PlasmaDensity,PlasmaDensityPrev0,PlasmaDensityPrev1,PlasmaDensityPrev;
  double NuPlus,NuMinus;

  auto Interpolate = [&] () {
    double x[3];
    Segment->GetCartesian(x, FieldLineCoord);

    D_SA_ptr->SetLocation(x);
    D_SA_ptr->Init();

    D_mu_mu_TwoWaves_ptr->SetLocation(x);
    D_mu_mu_TwoWaves_ptr->Init();

    D_mu_mu_ptr->SetLocation(x);
    D_mu_mu_ptr->Init(spec);


    D_x_x_ptr->SetLocation(x);
    D_x_x_ptr->Init(spec);

    Segment=FL::FieldLinesAll[iFieldLine].GetSegment(FieldLineCoord);
    if (Segment==NULL) return false;

    FL::cFieldLineVertex* VertexBegin=Segment->GetBegin();
    FL::cFieldLineVertex* VertexEnd=Segment->GetEnd();

    absB2=0.0;

    //get the magnetic field and the plasma waves at the corners of the segment
    B0=VertexBegin->GetDatum_ptr(FL::DatumAtVertexMagneticField);
    B1=VertexEnd->GetDatum_ptr(FL::DatumAtVertexMagneticField);

    W0=VertexBegin->GetDatum_ptr(FL::DatumAtVertexPlasmaWaves);
    W1=VertexEnd->GetDatum_ptr(FL::DatumAtVertexPlasmaWaves);

    VertexBegin->GetDatum(FL::DatumAtVertexPlasmaDensity,&PlasmaDensity0);
    VertexEnd->GetDatum(FL::DatumAtVertexPlasmaDensity,&PlasmaDensity1);

    VertexBegin->GetDatum(FL::DatumAtVertexPrevious::DatumAtVertexPlasmaDensity,&PlasmaDensityPrev0);
    VertexEnd->GetDatum(FL::DatumAtVertexPrevious::DatumAtVertexPlasmaDensity,&PlasmaDensityPrev1);

    x0=VertexBegin->GetX();
    x1=VertexEnd->GetX();

    //determine the interpolation coefficients
    w1=fmod(FieldLineCoord,1);
    w0=1.0-w1;

    for (int idim=0;idim<3;idim++) {
      double t;

      B[idim]=w0*B0[idim]+w1*B1[idim];
      absB2+=B[idim]*B[idim];

      t=w0*x0[idim]+w1*x1[idim];
      r2+=t*t;
    }

    W[0]=w0*W0[0]+w1*W1[0];
    W[1]=w0*W0[1]+w1*W1[1];
    PlasmaDensity=(w0*PlasmaDensity0+w1*PlasmaDensity1)*PIC::CPLR::SWMF::MeanPlasmaAtomicMass;
    PlasmaDensityPrev=(w0*PlasmaDensityPrev0+w1*PlasmaDensityPrev1)*PIC::CPLR::SWMF::MeanPlasmaAtomicMass;

    AbsB=sqrt(absB2);
    vAlfven=AbsB/sqrt(VacuumPermeability*PlasmaDensity);

    D_SA_ptr->SetW(W);
    D_SA_ptr->SetVelAlfven(vAlfven);
    D_SA_ptr->SetAbsB(AbsB);

    D_mu_mu_TwoWaves_ptr->SetW(W);
    D_mu_mu_TwoWaves_ptr->SetVelAlfven(vAlfven);
    D_mu_mu_TwoWaves_ptr->SetAbsB(AbsB);

    D_mu_mu_ptr->SetW(W);
    D_mu_mu_ptr->SetVelAlfven(vAlfven);
    D_mu_mu_ptr->SetAbsB(AbsB);

    return true;
  };


  v=sqrt(vParallel*vParallel+vNormal*vNormal);
  mu=vParallel/v;

  if (v>0.99*SpeedOfLight) {
    double t=0.99*SpeedOfLight/v;

    v=0.99*SpeedOfLight;
    vParallel*=t;
    vNormal*=t;
  }

  double ParticleStatWeight=node->block->GetLocalParticleWeight(spec);
  ParticleStatWeight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr);

  if (Interpolate()==false) exit(__LINE__,__FILE__"Error: the local coorsinate is outside of the field line");

  double dtSubStep=dtTotal;
  double dD_mu_mu_dmu_Plus,dD_mu_mu_dmu_Minus;
  bool FastParticleFlag=false;

  if (Interpolate()==false) exit(__LINE__,__FILE__"Error: the local coorsinate is outside of the field line");

  double speed=sqrt(vParallel*vParallel+vNormal*vNormal);
  mu=vParallel/speed;

  D_SA_ptr->SetVelocity(speed,mu);
  D_mu_mu_TwoWaves_ptr->SetVelocity(speed,mu);

  D_mu_mu_ptr->SetVelocity(speed,mu);

  double t0=SEP::Diffusion::AccelerationModelVelocitySwitchFactor*vAlfven;

  if (vNormal*vNormal+vParallel*vParallel>t0*t0) {
    //fast particle
    FastParticleFlag=true;
  }
  else {
    FastParticleFlag=false;

    double dD_mu_mu_dMu=D_mu_mu_TwoWaves_ptr->GetdDdMuSolarFrame();

    if (SEP::Diffusion::muTimeStepVariationLimitFlag==false) {
      if (fabs(dD_mu_mu_dMu)*dtSubStep>0.1) dtSubStep=0.1/fabs(dD_mu_mu_dMu);
    }

    if (std::isfinite(dD_mu_mu_dMu)==false) {
      dD_mu_mu_dMu=D_mu_mu_TwoWaves_ptr->GetdDdMuSolarFrame();
      exit(__LINE__,__FILE__,"Error: NAN is found");
    }
  }

  //integrate particle trajectory
  double DivVsw=0.0,ds;

  while (time_counter<dtTotal) {
    loop_cnt++;

    if (Interpolate()==false) break;

    //determine the which method should be used
    double MeanFreePath;

    D_x_x_ptr->SetVelocity(speed);
    MeanFreePath=D_x_x_ptr->GetMeanFreePath(FieldLineCoord,Segment,iFieldLine);


#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
    if (AMPS2SWMF::MagneticFieldLineUpdate::SecondCouplingFlag==true) {
      DivVsw=-log(PlasmaDensity/PlasmaDensityPrev)/(AMPS2SWMF::MagneticFieldLineUpdate::LastCouplingTime-AMPS2SWMF::MagneticFieldLineUpdate::LastLastCouplingTime);
    }
#else
    DivVsw=-log(PlasmaDensity/PlasmaDensityPrev)/dtTotal;
#endif

    double t0=SEP::Diffusion::AccelerationModelVelocitySwitchFactor*vAlfven;
    if (vNormal*vNormal+vParallel*vParallel>t0*t0) {
      //fast particle
      FastParticleFlag=true;
      dtSubStep=dtTotal-time_counter;
    }

    if ((FastParticleFlag==false)&&(SEP::Diffusion::AccelerationType==SEP::Diffusion::AccelerationTypeScattering)) {
      D_mu_mu_TwoWaves_ptr->SetVelocity(speed,mu);
      NuPlus=fabs(speed*mu)/D_mu_mu_TwoWaves_ptr->D_mu_mu_Plus.GetLambda();
      NuMinus=fabs(speed*mu)/D_mu_mu_TwoWaves_ptr->D_mu_mu_Minus.GetLambda();
    }

    double MovingTime,ScatteringTime;
    bool ScatteringFlag;


    //set the numerical limit on the number of simulated scattering events
    extern bool NumericalScatteringEventMode;
    extern double NumericalScatteringEventLimiter;

    //decide is scattering occured
    if ((FastParticleFlag==false)&&(SEP::Diffusion::AccelerationType==SEP::Diffusion::AccelerationTypeScattering)) {
      ScatteringTime=-log(rnd())/(NuPlus+NuMinus);

      if (time_counter+ScatteringTime<dtTotal) {
        //scattering occured
        ScatteringFlag=true;

        MovingTime=ScatteringTime;
        time_counter+=ScatteringTime;
      }
      else {
        //no scattering
        ScatteringFlag=false;

        MovingTime=dtTotal-time_counter;
        time_counter=dtTotal;
      }
    }
    else {
      ScatteringFlag=false;

      if (time_counter+dtSubStep<dtTotal) {
        MovingTime=dtSubStep;
        time_counter+=dtSubStep;
      }
      else {
        MovingTime=dtTotal-time_counter;
        time_counter=dtTotal;
      }
    }

    //determine the new particle pitch angle and location
    double L,AbsBDeriv;

    AbsBDeriv = (pow(B1[0]*B1[0] + B1[1]*B1[1] + B1[2]*B1[2], 0.5) -
        pow(B0[0]*B0[0] + B0[1]*B0[1] + B0[2]*B0[2], 0.5)) /  FL::FieldLinesAll[iFieldLine].GetSegmentLength(FieldLineCoord);

    L=-Vector3D::Length(B)/AbsBDeriv;
    mu+=(1.0-mu*mu)/(2.0*L)*MovingTime;

    if (mu<-1.0+muLimit) mu=-1.0+muLimit;
    if (mu>1.0-muLimit) mu=1.0-muLimit;

    //determine the shift of a particle position
    ds=MovingTime*speed*mu;

    //increment particle momentum
    double p=Relativistic::Speed2Momentum(speed,PIC::MolecularData::GetMass(spec));
    p*=exp(DivVsw*MovingTime/3.0);
    speed=Relativistic::Momentum2Speed(p,PIC::MolecularData::GetMass(spec));

    //limit scattering only with the incoming wave (if vParallel>0, then scatter only of the wave movinf with -vAlfven, or if vParallel<0, them scatter on the wave moveing with +vAlfven)
    if (LimitScatteringUpcomingWave==true) {
      if (mu>=0.0) NuPlus=0.0;
      else NuMinus=0.0;
    }


    if ((isfinite(speed)==false)||(isfinite(mu)==false)) {
      exit(__LINE__,__FILE__,"Error: NaN found");
    }

    //model scattering
    switch (SEP::Diffusion::AccelerationType) {
    case SEP::Diffusion::AccelerationTypeScattering:
      if (FastParticleFlag==false) {
        if (ScatteringFlag==true) {
          SEP::Diffusion::WaveScatteringModel(vAlfven,NuPlus,NuMinus,speed,mu);

          if ((isfinite(speed)==false)||(isfinite(mu)==false)) {
            exit(__LINE__,__FILE__,"Error: NaN found");
          }
        }
      }
      else {
        double muNew,pNew;
        double dMu;
        double x[3];
        Segment->GetCartesian(x, FieldLineCoord);

        D_mu_mu_ptr->SetVelocity(speed,mu);
        dMu=D_mu_mu_ptr->Get_dMu(MovingTime);


        if(fabs(dMu)<1.0) { // (MeanFreePath>ds) {
          //Mean free path is "large" -> integrate the pich angle evalution
          D_mu_mu_ptr->SetVelocity(speed,mu);
          D_SA_ptr->SetVelocity(speed,mu);

          muNew=D_mu_mu_ptr->DistributeMu(MovingTime);
          pNew=D_SA_ptr->DistributeP(MovingTime);

          if ((isfinite(muNew)==false)||(isfinite(pNew)==false)) {
            exit(__LINE__,__FILE__,"Error: NaN found");
          }

          mu=muNew;

          D_SA_ptr->Convert2Velocity();
          speed=D_SA_ptr->speed;
        }
        else {
          // Mean Free path is 'small" -> assume multiple scattering during the particle moving step
          D_x_x_ptr->SetVelocity(speed);
          ds=D_x_x_ptr->Get_ds(MovingTime,FieldLineCoord,Segment,iFieldLine);
          mu=-1.0+muLimit+rnd()*2.0*(1.0-muLimit);
        }
      }
      break;
    case SEP::Diffusion::AccelerationTypeDiffusion:
    {
      double muNew,pNew;

      if (MeanFreePath>ds) {
        //Mean free path is "large" -> integrate the pich angle evalution
        D_mu_mu_TwoWaves_ptr->SetVelocity(speed,mu);
        D_SA_ptr->SetVelocity(speed,mu);

        muNew=D_mu_mu_TwoWaves_ptr->DistributeMu(MovingTime);
        pNew=D_SA_ptr->DistributeP(MovingTime);

        if ((isfinite(muNew)==false)||(isfinite(pNew)==false)) {
          exit(__LINE__,__FILE__,"Error: NaN found");
        }

        mu=D_mu_mu_TwoWaves_ptr->mu;

        D_SA_ptr->Convert2Velocity();
        speed=D_SA_ptr->speed;
      }
      else {
        // Mean Free path is 'small" -> assume multiple scattering during the particle moving step
        D_x_x_ptr->SetVelocity(speed);
        ds=D_x_x_ptr->Get_ds(MovingTime,FieldLineCoord,Segment,iFieldLine);
        mu=-1.0+muLimit+rnd()*2.0*(1.0-muLimit);
      }
    }
    break;
    default:
      exit(__LINE__,__FILE__,"Error: the oprion is not recognized");
    }
  }

  //update the particle location
  FieldLineCoord=FL::FieldLinesAll[iFieldLine].move(FieldLineCoord,ds);

  //get the segment of the new particle location
  if ((Segment=FL::FieldLinesAll[iFieldLine].GetSegment(FieldLineCoord))==NULL) {
    //the particle left the computational domain
    int code=_PARTICLE_DELETED_ON_THE_FACE_;

    //call the function that process particles that leaved the coputational domain
    switch (code) {
    case _PARTICLE_DELETED_ON_THE_FACE_:
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;

    default:
      exit(__LINE__,__FILE__,"Error: not implemented");
    }
  }


  vParallel=speed*mu;
  vNormal=speed*sqrt(1.0-mu*mu);

  //set the new values of the normal and parallel particle velocities
  PB::SetVParallel(vParallel,ParticleData);
  PB::SetVNormal(vNormal,ParticleData);

  if (std::isfinite(vParallel)==false) exit(__LINE__,__FILE__);
  if (std::isfinite(vNormal)==false) exit(__LINE__,__FILE__);

  //set the new particle coordinate
  PB::SetFieldLineCoord(FieldLineCoord,ParticleData);

  //attach the particle to the temporaty list
  // Step 5 guarantees segment attachment at compile time, so the mover can
  // commit directly to the destination segment without a mesh-cell branch.
  {
    long int temp=Segment->tempFirstParticleIndex.exchange(ptr);
      PIC::ParticleBuffer::SetNext(temp,ParticleData);
      PIC::ParticleBuffer::SetPrev(-1,ParticleData);

      if (temp!=-1) PIC::ParticleBuffer::SetPrev(ptr,temp);
  }

  return _PARTICLE_MOTION_FINISHED_;
}


//===================================================================================================================================



int SEP::ParticleMover_FTE(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  namespace PB = PIC::ParticleBuffer;
  namespace FL = PIC::FieldLine;

  PIC::ParticleBuffer::byte *ParticleData;
  double Speed,AbsB,L,vParallel,vNormal,DivAbsB,vParallelInit,vNormalInit;
  double FieldLineCoord,xCartesian[3];
  int iFieldLine,spec;
  FL::cFieldLineSegment *Segment;

  static int ncall=0;
  ncall++;

  if (ncall==290946) {
	  double rr=0.0;
  }

  ParticleData=PB::GetParticleDataPointer(ptr);

  FieldLineCoord=PB::GetFieldLineCoord(ParticleData);
  iFieldLine=PB::GetFieldLineId(ParticleData);
  spec=PB::GetI(ParticleData);

  //velocity is in the frame moving with solar wind
  vParallel=PB::GetVParallel(ParticleData);
  vNormal=PB::GetVNormal(ParticleData);

  if ((Segment=FL::FieldLinesAll[iFieldLine].GetSegment(FieldLineCoord))==NULL) {
    exit(__LINE__,__FILE__,"Error: cannot find the segment");
  }

  double TimeCounter=0.0,dt;
  double D_mumu,ds;
  double energy,vnew[3],l[3],x[3],rHelio,v,mu,dmu,dmu_mean;

  Segment->GetCartesian(x,FieldLineCoord);
  rHelio=Vector3D::Length(x);
  Speed=sqrt(vNormal*vNormal+vParallel*vParallel);

  v=sqrt(vParallel*vParallel+vNormal*vNormal);

  // ---------------------------------------------------------------------------
  // Numerical/physical sanity check.  The focused-transport mover is formulated
  // in terms of the pitch-angle cosine mu=v_parallel/v.  A zero-speed particle
  // would therefore immediately generate a division by zero and later could be
  // written back to the particle buffer as v_parallel=v_normal=0.  Such a particle
  // is not a valid SEP macro-particle; it usually indicates an injection or an
  // earlier mover/coupling error.  Stop with a detailed message instead of
  // silently propagating NaNs or zero velocities into the turbulence coupling.
  // ---------------------------------------------------------------------------
  const double MinFocusedTransportSpeed=1.0; // [m/s], only a numerical floor
  if ((isfinite(v)==false)||(v<=MinFocusedTransportSpeed)) {
    char msg[512];
    sprintf(msg,
        "Error: ParticleMover_FTE received a particle with invalid/zero speed. "
        "v=%e m/s, vParallel=%e m/s, vNormal=%e m/s, fieldLine=%i, coord=%e. "
        "A focused-transport particle must have non-zero speed because mu=vParallel/v.",
        v,vParallel,vNormal,iFieldLine,FieldLineCoord);
    exit(__LINE__,__FILE__,msg);
  }

  mu=vParallel/v;

  //get D_mu_mu and evaluate the time substep
  D_mumu=QLT::calculateDmuMu(v,mu,rHelio);

  // Mean absolute value of a standard normal deviate is sqrt(2/pi)=0.79788456.
  // The previous value was smaller by a factor of ten, which made the adaptive
  // substep estimate inconsistent with the actual stochastic pitch-angle kick.
  dmu_mean=sqrt(2.0*D_mumu*dtTotal)*0.7978845608028654;

  if (dmu_mean>0.2) {
    double t=0.2/dmu_mean;

    if (t<0.2) {
      //the subtime step should be too small -> switch Parker eq instead
      return SEP::ParticleMover_Parker_MeanFreePath(ptr,dtTotal,node);
    }

    dt=dtTotal*t*t;
  }
  else {
    dt=dtTotal;
  }

  FL::cFieldLineSegment *LastSegment=NULL;
  double B;

  while (TimeCounter<dtTotal) {
    // Clip the local substep to the remaining part of the AMPS particle time step.
    // Without this clipping, the last substep can advance the particle and the
    // turbulence-coupling path integral beyond dtTotal.
    double dt_step=std::min(dt,dtTotal-TimeCounter);
    if (dt_step<=0.0) break;

    D_mumu=QLT::calculateDmuMu(v,mu,rHelio);

    // Store the velocity components used during this explicit transport substep.
    // The pitch angle is updated below for the next substep, but the spatial
    // displacement over the current substep is computed with the old velocity.
    // These same old components must be passed to the wave-particle coupling
    // accumulator so that the resonant streaming source is consistent with the
    // path actually traveled by the particle.
    double vParallel_substep=vParallel;
    double vNormal_substep=vNormal;

    ds=vParallel_substep*dt_step;
    double FieldLineCoordStart=FieldLineCoord;

    //calculate L
    if (Segment!=LastSegment) {
      auto b=Segment->GetBegin();
      auto e=Segment->GetEnd();
      double *B0,*B1,dB=0.0;

      B0=b->GetMagneticField();
      B1=e->GetMagneticField();
      B=0.0;

      for (int i=0;i<3;i++) {
        double t;

	t=0.5*(B0[i]+B1[i]);
	B+=t*t;

	t=B1[i]-B0[i];
	dB+=t*t;
      }

      if ((dB>0.0)&&(B>0.0)) {
        L=-Segment->GetLength()*sqrt(B/dB);
      }
      else {
        // Uniform |B| over the segment gives an infinite focusing length and
        // therefore no deterministic focusing contribution.
        L=1.0e100;
      }
      LastSegment=Segment;
    }

    dmu=0.0;
    if ((isfinite(L)==true)&&(fabs(L)>1.0e-100)) {
      dmu=-(1.0-mu*mu)/(2.0*L)*v*dt_step;
    }
    if ((D_mumu>0.0)&&(isfinite(D_mumu)==true)) {
      dmu+=sqrt(2.0*D_mumu*dt_step)*Vector3D::Distribution::Normal();
    }
    mu+=dmu;

    // Reflect the pitch-angle cosine at the physical boundaries mu=+/-1.  The
    // earlier implementation wrapped mu by subtracting one, e.g. 1.1 -> 0.1,
    // which is not a reflecting boundary and produces an artificial large-angle
    // scattering event.  The reflection below preserves the distance beyond the
    // boundary: 1.1 -> 0.9 and -1.1 -> -0.9.
    while ((-1.0>mu)||(mu>1.0)) {
      if (mu>1.0) mu=2.0-mu;
      if (mu<-1.0) mu=-2.0-mu;
    }

    if (mu>1.0-1.0e-12) mu=1.0-1.0e-12;
    if (mu<-1.0+1.0e-12) mu=-1.0+1.0e-12;


    vParallel=mu*v;
    FieldLineCoord=FL::FieldLinesAll[iFieldLine].move(FieldLineCoord,ds,Segment);

    if (Segment==NULL) {
      //the particle has left the simulation, and it is need to be deleted
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;
    }

    // Feed the manager-based wave-particle coupling with the actual path traveled
    // during this focused-transport substep.  This is the missing link that makes
    // the default FTE mover compatible with both the integrated and the
    // wave-number-resolved turbulence models: the coupling manager can update
    // E_+(k) and E_-(k) only after the mover fills G_+(k) and G_-(k).
    //
    // The normalization time passed to AccumulateParticleFluxForWaveCoupling() is
    // dtTotal, not dt_step.  The coupling arrays represent the time-averaged
    // streaming source over the full AMPS particle step.  Since this mover can split
    // dtTotal into several substeps, each substep contributes only the fraction of
    // the full-step path/time it actually covers.
    if (SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag &&
        SEP::AlfvenTurbulence_Kolmogorov::ParticleCouplingMode &&
        fabs(ds)>0.0) {
      SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::AccumulateParticleFluxForWaveCoupling(
          iFieldLine,ptr,dtTotal,
          vParallel_substep,vNormal_substep,
          FieldLineCoordStart,FieldLineCoord,ds);
    }



    TimeCounter+=dt_step;

    Segment->GetCartesian(x,FieldLineCoord);
    rHelio=Vector3D::Length(x);
  }

  //set the new values of the normal and parallel particle velocities
  // Use a guarded square-root argument because roundoff can make 1-mu^2 slightly
  // negative when |mu| is extremely close to one.
  vNormal=sqrt(std::max(0.0,1.0-mu*mu))*v;
  if ((isfinite(vNormal)==false)||(isfinite(vParallel)==false)||
      (sqrt(vParallel*vParallel+vNormal*vNormal)<=MinFocusedTransportSpeed)) {
    char msg[512];
    sprintf(msg,
        "Error: ParticleMover_FTE produced an invalid/zero final velocity. "
        "vParallel=%e m/s, vNormal=%e m/s, mu=%e, speed=%e m/s, fieldLine=%i, coord=%e",
        vParallel,vNormal,mu,sqrt(vParallel*vParallel+vNormal*vNormal),iFieldLine,FieldLineCoord);
    exit(__LINE__,__FILE__,msg);
  }

  PB::SetVParallel(vParallel,ParticleData);
  PB::SetVNormal(vNormal,ParticleData);

  //set the new particle coordinate
  PB::SetFieldLineCoord(FieldLineCoord,ParticleData);

  //attach the particle to the temporaty list
  // Step 5 guarantees segment attachment at compile time, so the mover can
  // commit directly to the destination segment without a mesh-cell branch.
  {
    long int temp=Segment->tempFirstParticleIndex.exchange(ptr);
    PIC::ParticleBuffer::SetNext(temp,ParticleData);
    PIC::ParticleBuffer::SetPrev(-1,ParticleData);

    if (temp!=-1) PIC::ParticleBuffer::SetPrev(ptr,temp);
  }

  return _PARTICLE_MOTION_FINISHED_;
}

int SEP::ParticleMover_Parker_MeanFreePath(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  namespace PB = PIC::ParticleBuffer;
  namespace FL = PIC::FieldLine;

  PIC::ParticleBuffer::byte *ParticleData;
  double Speed,mu,AbsB,L,vParallel,vNormal,DivAbsB,vParallelInit,vNormalInit;
  double FieldLineCoord,xCartesian[3];
  int iFieldLine,spec;
  FL::cFieldLineSegment *Segment;

  ParticleData=PB::GetParticleDataPointer(ptr);

  FieldLineCoord=PB::GetFieldLineCoord(ParticleData);
  iFieldLine=PB::GetFieldLineId(ParticleData);
  spec=PB::GetI(ParticleData);

  //velocity is in the frame moving with solar wind
  vParallel=PB::GetVParallel(ParticleData);
  vNormal=PB::GetVNormal(ParticleData);

  if ((Segment=FL::FieldLinesAll[iFieldLine].GetSegment(FieldLineCoord))==NULL) {
    exit(__LINE__,__FILE__,"Error: cannot find the segment");
  }

  double TimeCounter=0.0,dt;
  double MeanFreePath,ds;
  double energy,vnew[3],l[3],x[3],rHelio,dxx;

  Segment->GetCartesian(x,FieldLineCoord);
  rHelio=Vector3D::Length(x);
  Speed=sqrt(vNormal*vNormal+vParallel*vParallel);

  while (TimeCounter<dtTotal) {
    //get the value of the backgound magnetic field
    double AbsB;

    switch (_PIC_COUPLER_MODE_) {
    case _PIC_COUPLER_MODE__SWMF_:
      AbsB=SEP::FieldLineData::GetAbsB(FieldLineCoord,Segment,iFieldLine);
      break;
    default:
      AbsB=SEP::ParkerSpiral::GetAbsB(rHelio);
    }

    switch (SEP::Scattering::MeanFreePathMode) {
    case SEP::Scattering::MeanFreePathMode_QLT:
      MeanFreePath=QLT::calculateMeanFreePath(rHelio,Speed);
      break;
    case SEP::Scattering::MeanFreePathMode_QLT1:
      MeanFreePath=QLT1::calculateMeanFreePath(rHelio,Speed,AbsB);
      break;
    case SEP::Scattering::MeanFreePathMode_Tenishev2005AIAA:
      energy=Relativistic::Speed2E(Speed,PIC::MolecularData::GetMass(spec));

      MeanFreePath=SEP::Scattering::Tenishev2005AIAA::lambda0*
        pow(energy/GeV2J,SEP::Scattering::Tenishev2005AIAA::alpha)*
        pow(rHelio/_AU_,SEP::Scattering::Tenishev2005AIAA::beta);
      break;
    case SEP::Scattering::MeanFreePathMode_Chen2024AA:
       energy=Relativistic::Speed2E(Speed,PIC::MolecularData::GetMass(spec));
       dxx=SEP::Diffusion::Chen2024AA::GetDxx(rHelio,energy);
       MeanFreePath=3.0*dxx/Speed; //Eq, 15, Liu-2024-arXiv;
       break;

    default:
      exit(__LINE__,__FILE__,"Error: the oprion is unknown");
    }


    if ((SEP::Offset::MeanFreePath!=-1)&&(SEP::Sampling::MeanFreePath::active_flag==true)) {
      *((double*)(ParticleData+SEP::Offset::MeanFreePath))=MeanFreePath;
    }



    ds=-MeanFreePath*log(rnd());
    dt=ds/fabs(vParallel);


    if (TimeCounter+dt<dtTotal) {
      //scattering occured begore the end of the time interval
      FieldLineCoord=FL::FieldLinesAll[iFieldLine].move(FieldLineCoord,ds,Segment);

      if (Segment==NULL) {
        //the particle has left the simulation, and it is need to be deleted
        PIC::ParticleBuffer::DeleteParticle(ptr);
        return _PARTICLE_LEFT_THE_DOMAIN_;
      }


      //simulate scattering of the particle
      Vector3D::Distribution::Uniform(vnew,Speed);

      Segment->GetDir(l);
      Vector3D::GetComponents(vParallel,vNormal,vnew,l);
    }
    else {
      //scattering does not occur before the enf of the simulated time interval
      dt=dtTotal-TimeCounter;
      ds=vParallel*dt;

      FieldLineCoord=FL::FieldLinesAll[iFieldLine].move(FieldLineCoord,ds,Segment);

      if (Segment==NULL) {
        //the particle has left the simulation, and it is need to be deleted
        PIC::ParticleBuffer::DeleteParticle(ptr);
        return _PARTICLE_LEFT_THE_DOMAIN_;
      }

    }


    TimeCounter+=dt;

    Segment->GetCartesian(x,FieldLineCoord);
    rHelio=Vector3D::Length(x);
  }

  //set the new values of the normal and parallel particle velocities
  PB::SetVParallel(vParallel,ParticleData);
  PB::SetVNormal(vNormal,ParticleData);

  //set the new particle coordinate
  PB::SetFieldLineCoord(FieldLineCoord,ParticleData);

  //attach the particle to the temporaty list
  // Step 5 guarantees segment attachment at compile time, so the mover can
  // commit directly to the destination segment without a mesh-cell branch.
  {
    long int temp=Segment->tempFirstParticleIndex.exchange(ptr);
    PIC::ParticleBuffer::SetNext(temp,ParticleData);
    PIC::ParticleBuffer::SetPrev(-1,ParticleData);

    if (temp!=-1) PIC::ParticleBuffer::SetPrev(ptr,temp);
  }

  return _PARTICLE_MOTION_FINISHED_;
}


int SEP::ParticleMover_Parker_Dxx(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  namespace PB = PIC::ParticleBuffer;
  namespace FL = PIC::FieldLine;

  PIC::ParticleBuffer::byte *ParticleData;
  double Speed,mu,AbsB,L,vParallel,vNormal,DivAbsB,vParallelInit,vNormalInit;
  double FieldLineCoord,xCartesian[3];
  int iFieldLine,spec;
  FL::cFieldLineSegment *Segment;

  double totalTraversedPath=0.0;

  ParticleData=PB::GetParticleDataPointer(ptr);

  FieldLineCoord=PB::GetFieldLineCoord(ParticleData);
  iFieldLine=PB::GetFieldLineId(ParticleData);
  spec=PB::GetI(ParticleData);

  //velocity is in the frame moving with solar wind
  vParallel=PB::GetVParallel(ParticleData);
  vNormal=PB::GetVNormal(ParticleData);

  if ((Segment=FL::FieldLinesAll[iFieldLine].GetSegment(FieldLineCoord))==NULL) {
    exit(__LINE__,__FILE__,"Error: cannot find the segment");
  }

  // Save initial values for Parker flux sampling
  double s_init = FieldLineCoord;
  double dtTotal_saved = dtTotal;
  FL::cFieldLineSegment *segment_start = Segment;

  double TimeCounter=0.0,dt;
  double MeanFreePath,ds,Dxx,dDxx_ds;
  double energy,vnew[3],l[3],x[3],r;

  Segment->GetCartesian(x,FieldLineCoord);
  r=Vector3D::Length(x);
  Speed=sqrt(vNormal*vNormal+vParallel*vParallel);

  // Helper function to calculate mean free path (adapted from ParticleMover_Parker_MeanFreePath)
  auto CalculateMeanFreePath = [&] (int spec, double rHelio, double Speed, double AbsB) {
    double MeanFreePath, dxx, energy;

    switch (SEP::Scattering::MeanFreePathMode) {
    case SEP::Scattering::MeanFreePathMode_QLT:
      MeanFreePath = QLT::calculateMeanFreePath(rHelio, Speed);
      break;
    case SEP::Scattering::MeanFreePathMode_QLT1:
      MeanFreePath = QLT1::calculateMeanFreePath(rHelio, Speed, AbsB);
      break;
    case SEP::Scattering::MeanFreePathMode_Tenishev2005AIAA:
      energy = Relativistic::Speed2E(Speed, PIC::MolecularData::GetMass(spec));
      MeanFreePath = SEP::Scattering::Tenishev2005AIAA::lambda0 *
        pow(energy/GeV2J, SEP::Scattering::Tenishev2005AIAA::alpha) *
        pow(rHelio/_AU_, SEP::Scattering::Tenishev2005AIAA::beta);
      break;
    case SEP::Scattering::MeanFreePathMode_Chen2024AA:
      energy = Relativistic::Speed2E(Speed, PIC::MolecularData::GetMass(spec));
      dxx = SEP::Diffusion::Chen2024AA::GetDxx(rHelio, energy);
      MeanFreePath = 3.0 * dxx / Speed; // Eq. 15, Liu-2024-arXiv
      break;
    default:
      exit(__LINE__, __FILE__, "Error: the option is unknown");
    }

    return MeanFreePath;
  };

  // Helper function to calculate Dxx from mean free path
  auto CalculateDxx = [&] (double MeanFreePath, double Speed) {
    return MeanFreePath * Speed / 3.0;  // Standard relation: Dxx = λ * v / 3
  };

  // Helper function to calculate d(Dxx)/ds using finite differences
  auto CalculateDxxGradient = [&] (double currentFieldLineCoord, FL::cFieldLineSegment* currentSegment,
                                   int fieldLineId, double currentSpeed, int particleSpec) {
    double dDxx_ds = 0.0;
    double delta_s = 0.01 * currentSegment->GetLength(); // Small displacement for finite difference

    // Get current position and properties
    double x_current[3], r_current;
    currentSegment->GetCartesian(x_current, currentFieldLineCoord);
    r_current = Vector3D::Length(x_current);

    // Get magnetic field at current position
    double AbsB_current;
    switch (_PIC_COUPLER_MODE_) {
    case _PIC_COUPLER_MODE__SWMF_:
      AbsB_current = SEP::FieldLineData::GetAbsB(currentFieldLineCoord, currentSegment, fieldLineId);
      break;
    default:
      AbsB_current = SEP::ParkerSpiral::GetAbsB(r_current);
    }

    double MeanFreePath_current = CalculateMeanFreePath(particleSpec, r_current, currentSpeed, AbsB_current);
    double Dxx_current = CalculateDxx(MeanFreePath_current, currentSpeed);

    // Try forward difference
    double FieldLineCoord_forward = FL::FieldLinesAll[fieldLineId].move(currentFieldLineCoord, delta_s, currentSegment);
    FL::cFieldLineSegment* Segment_forward = FL::FieldLinesAll[fieldLineId].GetSegment(FieldLineCoord_forward);

    if (Segment_forward != NULL) {
      double x_forward[3], r_forward;
      Segment_forward->GetCartesian(x_forward, FieldLineCoord_forward);
      r_forward = Vector3D::Length(x_forward);

      double AbsB_forward;
      switch (_PIC_COUPLER_MODE_) {
      case _PIC_COUPLER_MODE__SWMF_:
        AbsB_forward = SEP::FieldLineData::GetAbsB(FieldLineCoord_forward, Segment_forward, fieldLineId);
        break;
      default:
        AbsB_forward = SEP::ParkerSpiral::GetAbsB(r_forward);
      }

      double MeanFreePath_forward = CalculateMeanFreePath(particleSpec, r_forward, currentSpeed, AbsB_forward);
      double Dxx_forward = CalculateDxx(MeanFreePath_forward, currentSpeed);

      // Try backward difference
      double FieldLineCoord_backward = FL::FieldLinesAll[fieldLineId].move(currentFieldLineCoord, -delta_s, currentSegment);
      FL::cFieldLineSegment* Segment_backward = FL::FieldLinesAll[fieldLineId].GetSegment(FieldLineCoord_backward);

      if (Segment_backward != NULL) {
        double x_backward[3], r_backward;
        Segment_backward->GetCartesian(x_backward, FieldLineCoord_backward);
        r_backward = Vector3D::Length(x_backward);

        double AbsB_backward;
        switch (_PIC_COUPLER_MODE_) {
        case _PIC_COUPLER_MODE__SWMF_:
          AbsB_backward = SEP::FieldLineData::GetAbsB(FieldLineCoord_backward, Segment_backward, fieldLineId);
          break;
        default:
          AbsB_backward = SEP::ParkerSpiral::GetAbsB(r_backward);
        }

        double MeanFreePath_backward = CalculateMeanFreePath(particleSpec, r_backward, currentSpeed, AbsB_backward);
        double Dxx_backward = CalculateDxx(MeanFreePath_backward, currentSpeed);

        // Central difference
        dDxx_ds = (Dxx_forward - Dxx_backward) / (2.0 * delta_s);
      } else {
        // Forward difference only
        dDxx_ds = (Dxx_forward - Dxx_current) / delta_s;
      }
    } else {
      // Try backward difference only
      double FieldLineCoord_backward = FL::FieldLinesAll[fieldLineId].move(currentFieldLineCoord, -delta_s, currentSegment);
      FL::cFieldLineSegment* Segment_backward = FL::FieldLinesAll[fieldLineId].GetSegment(FieldLineCoord_backward);

      if (Segment_backward != NULL) {
        double x_backward[3], r_backward;
        Segment_backward->GetCartesian(x_backward, FieldLineCoord_backward);
        r_backward = Vector3D::Length(x_backward);

        double AbsB_backward;
        switch (_PIC_COUPLER_MODE_) {
        case _PIC_COUPLER_MODE__SWMF_:
          AbsB_backward = SEP::FieldLineData::GetAbsB(FieldLineCoord_backward, Segment_backward, fieldLineId);
          break;
        default:
          AbsB_backward = SEP::ParkerSpiral::GetAbsB(r_backward);
        }

        double MeanFreePath_backward = CalculateMeanFreePath(particleSpec, r_backward, currentSpeed, AbsB_backward);
        double Dxx_backward = CalculateDxx(MeanFreePath_backward, currentSpeed);

        // Backward difference
        dDxx_ds = (Dxx_current - Dxx_backward) / delta_s;
      } else {
        // Cannot calculate gradient - set to zero
        dDxx_ds = 0.0;
      }
    }

    return dDxx_ds;
  };

  while (TimeCounter < dtTotal) {
    // Get current magnetic field magnitude
    double AbsB;
    switch (_PIC_COUPLER_MODE_) {
    case _PIC_COUPLER_MODE__SWMF_:
      AbsB = SEP::FieldLineData::GetAbsB(FieldLineCoord, Segment, iFieldLine);
      break;
    default:
      AbsB = SEP::ParkerSpiral::GetAbsB(r);
    }

    // Calculate mean free path using the same method as ParticleMover_Parker_MeanFreePath
    MeanFreePath = CalculateMeanFreePath(spec, r, Speed, AbsB);

    // Calculate Dxx from mean free path
    Dxx = CalculateDxx(MeanFreePath, Speed);

    // Calculate gradient of Dxx
    dDxx_ds = CalculateDxxGradient(FieldLineCoord, Segment, iFieldLine, Speed, spec);

    // Store Dxx in particle data if sampling is active
    if ((SEP::Offset::MeanFreePath != -1) && (SEP::Sampling::MeanFreePath::active_flag == true)) {
      *((double*)(ParticleData + SEP::Offset::MeanFreePath)) = MeanFreePath;
    }

    // Determine time step based on scattering time scale
    double scattering_time = MeanFreePath / Speed;  // Characteristic scattering time scale
    dt = -scattering_time * log(rnd());  // Exponential distribution for scattering events

    // Limit time step to not exceed remaining simulation time
    if (TimeCounter + dt > dtTotal) {
      dt = dtTotal - TimeCounter;
    }

    dt=dtTotal;

    // Calculate particle displacement accounting for diffusion and its gradient
    // Based on the solution to: ds/dt = dDxx/ds + sqrt(2*Dxx) * η(t)
    // where η(t) is white noise

    double stochastic_displacement = sqrt(2.0 * Dxx * dt) * Vector3D::Distribution::Normal();
    double deterministic_displacement = dDxx_ds * dt;

    ds = deterministic_displacement + stochastic_displacement;

    double debug_effective_speed=ds/dt;
    if (fabs(debug_effective_speed)>Speed) ds*=Speed*dt/fabs(ds);

    totalTraversedPath+=ds;


    // Move particle along field line
    FieldLineCoord = FL::FieldLinesAll[iFieldLine].move(FieldLineCoord, ds, Segment);

    if (Segment == NULL) {
      // The particle has left the simulation domain - sample flux before deletion
      double s_final = FieldLineCoord;

if (SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag) SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::AccumulateParticleFluxForWaveCoupling(
    iFieldLine, //int field_line_idx,
    ptr, //long int particle_index,
    dtTotal_saved, //double dt,
    Speed, //double speed,
    s_init, //double s_start,
    s_final, //double s_finish,
    totalTraversedPath
);



      // Now delete the particle
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;
    }

    // Update particle position and radial distance
    Segment->GetCartesian(x, FieldLineCoord);
    r = Vector3D::Length(x);

    // Check if scattering occurred (based on whether we used the full scattering time)
    if (TimeCounter + scattering_time * (-log(rnd())) < dtTotal) {
      // Scattering event occurred - randomize velocity direction
      Vector3D::Distribution::Uniform(vnew, Speed);
      Segment->GetDir(l);
      Vector3D::GetComponents(vParallel, vNormal, vnew, l);
    }

    TimeCounter += dt;
  }

  // Sample Parker flux using the final position
  double s_final = FieldLineCoord;


if (SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag) SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::AccumulateParticleFluxForWaveCoupling(
    iFieldLine, //int field_line_idx,
    ptr, //long int particle_index,
    dtTotal_saved, //double dt,
    Speed, //double speed,
    s_init, //double s_start,
    s_final, //double s_finish,
    totalTraversedPath
);





  // Set the new values of the normal and parallel particle velocities
  PB::SetVParallel(vParallel, ParticleData);
  PB::SetVNormal(vNormal, ParticleData);

  // Set the new particle coordinate
  PB::SetFieldLineCoord(FieldLineCoord, ParticleData);

  // Attach the particle to the temporary list
  // Step 5 guarantees segment attachment at compile time, so the mover can
  // commit directly to the destination segment without a mesh-cell branch.
  {
    long int temp = Segment->tempFirstParticleIndex.exchange(ptr);
      PIC::ParticleBuffer::SetNext(temp, ParticleData);
      PIC::ParticleBuffer::SetPrev(-1, ParticleData);

      if (temp != -1) PIC::ParticleBuffer::SetPrev(ptr, temp);
  }

  return _PARTICLE_MOTION_FINISHED_;
}

int SEP::ParticleMover_Tenishev_2005_FL(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  namespace PB = PIC::ParticleBuffer;
  namespace FL = PIC::FieldLine;

  PIC::ParticleBuffer::byte *ParticleData;
  double mu,AbsB,L,vParallel,vNormal,v,DivAbsB,vParallelInit,vNormalInit;
  double FieldLineCoord,xCartesian[3];
  int iFieldLine,spec;
  FL::cFieldLineSegment *Segment;

  ParticleData=PB::GetParticleDataPointer(ptr);

  FieldLineCoord=PB::GetFieldLineCoord(ParticleData);
  iFieldLine=PB::GetFieldLineId(ParticleData);
  spec=PB::GetI(ParticleData);

  //velocity is in the frame moving with solar wind
  vParallel=PB::GetVParallel(ParticleData);
  vNormal=PB::GetVNormal(ParticleData);

  //shift location of the particle
  FieldLineCoord=FL::FieldLinesAll[iFieldLine].move(FieldLineCoord,dtTotal*vParallel);

  //get the segment of the new particle location
  if ((Segment=FL::FieldLinesAll[iFieldLine].GetSegment(FieldLineCoord))==NULL) {
    //the particle left the computational domain
    int code=_PARTICLE_DELETED_ON_THE_FACE_;

    //call the function that process particles that leaved the coputational domain
    switch (code) {
    case _PARTICLE_DELETED_ON_THE_FACE_:
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;

    default:
      exit(__LINE__,__FILE__,"Error: not implemented");
    }
  }


  //simulate scattering of the particle
  if (SEP::Scattering::Tenishev2005AIAA::status==SEP::Scattering::Tenishev2005AIAA::_enabled) {
    double Speed,p,energy,lambda;

    FL::FieldLinesAll[iFieldLine].GetCartesian(xCartesian,FieldLineCoord);

    Speed=sqrt(vParallel*vParallel+vNormal*vNormal);
    energy=Relativistic::Speed2E(Speed,PIC::MolecularData::GetMass(spec));

    lambda=SEP::Scattering::Tenishev2005AIAA::lambda0*
      pow(energy/GeV2J,SEP::Scattering::Tenishev2005AIAA::alpha)*
      pow(Vector3D::Length(xCartesian)/_AU_,SEP::Scattering::Tenishev2005AIAA::beta);

    //the prabability of scattering event during the current time step
    p=1.0-exp(-dtTotal*Speed/lambda);

    if (p>rnd()) {
      //scattering occured
      double vnew[3],l[3];

      Vector3D::Distribution::Uniform(vnew,Speed);
      Segment->GetDir(l);

      Vector3D::GetComponents(vParallel,vNormal,vnew,l);
    }
  }

  //set the new values of the normal and parallel particle velocities
  PB::SetVParallel(vParallel,ParticleData);
  PB::SetVNormal(vNormal,ParticleData);

  //set the new particle coordinate
  PB::SetFieldLineCoord(FieldLineCoord,ParticleData);

  //attach the particle to the temporaty list
  // Step 5 guarantees segment attachment at compile time, so the mover can
  // commit directly to the destination segment without a mesh-cell branch.
  {
    long int temp=Segment->tempFirstParticleIndex.exchange(ptr);
      PIC::ParticleBuffer::SetNext(temp,ParticleData);
      PIC::ParticleBuffer::SetPrev(-1,ParticleData);

      if (temp!=-1) PIC::ParticleBuffer::SetPrev(ptr,temp);
  }

  return _PARTICLE_MOTION_FINISHED_;
}


//=============================================================================================================
void SEP::GetTransportCoefficients (double& dP,double& dLogP,double& dmu,double v,double mu,PIC::FieldLine::cFieldLineSegment *Segment,double FieldLineCoord,double dt,int iFieldLine,double& vSolarWindParallel) {
  namespace PB = PIC::ParticleBuffer;
  namespace FL = PIC::FieldLine;

  //calculate B and L
  double B[3],B0[3],B1[3], AbsBDeriv;
  double L,AbsB;

  FL::FieldLinesAll[iFieldLine].GetMagneticField(B0, (int)FieldLineCoord);
  FL::FieldLinesAll[iFieldLine].GetMagneticField(B,       FieldLineCoord);
  FL::FieldLinesAll[iFieldLine].GetMagneticField(B1, (int)FieldLineCoord+1-1E-7);
  AbsB   = pow(B[0]*B[0] + B[1]*B[1] + B[2]*B[2], 0.5);

  AbsBDeriv = (pow(B1[0]*B1[0] + B1[1]*B1[1] + B1[2]*B1[2], 0.5) -
    pow(B0[0]*B0[0] + B0[1]*B0[1] + B0[2]*B0[2], 0.5)) /  FL::FieldLinesAll[iFieldLine].GetSegmentLength(FieldLineCoord);

  L=-Vector3D::Length(B)/AbsBDeriv;

  #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
  if (::AMPS2SWMF::MagneticFieldLineUpdate::SecondCouplingFlag==false) {
    dLogP=0.0,dP=0.0;
    dmu=(1.0-mu*mu)/2.0*v/L*dt;
    return;
  }
  #else
  dLogP=0.0,dP=0.0;
  dmu=(1.0-mu*mu)/2.0*v/L*dt;
  return;
  #endif

  //calculate dVsw_dz
  double vSolarWind[3],vSW1,vSW0,dVz_dz;

  FL::FieldLinesAll[iFieldLine].GetPlasmaVelocity(vSolarWind,(int)FieldLineCoord);
  vSW0=Vector3D::DotProduct(vSolarWind,B0)/Vector3D::Length(B0);

  FL::FieldLinesAll[iFieldLine].GetPlasmaVelocity(vSolarWind,(int)FieldLineCoord+1-1E-7);
  vSW1=Vector3D::DotProduct(vSolarWind,B1)/Vector3D::Length(B1);

  dVz_dz=(vSW1-vSW0)/FL::FieldLinesAll[iFieldLine].GetSegmentLength(FieldLineCoord);

  //calculate div(vSW) : Dln(Rho)=-div(vSW)*dt
  double PlasmaDensityCurrent,PlasmaDensityOld,DivVsw,PlasmaDensityCurrentParticle=0.0,PlasmaDensityOldParticle=0.0;
  auto Vertex0=Segment->GetBegin();
  auto Vertex1=Segment->GetEnd();

  double weight0=1.0-(FieldLineCoord-floor(FieldLineCoord));
  double weight1=1.0-weight0;

  Vertex0->GetDatum(FL::DatumAtVertexPlasmaDensity,&PlasmaDensityCurrent);
  Vertex0->GetDatum(FL::DatumAtVertexPrevious::DatumAtVertexPlasmaDensity,&PlasmaDensityOld);

  PlasmaDensityCurrentParticle=weight0*PlasmaDensityCurrent;
  PlasmaDensityOldParticle=weight0*PlasmaDensityOld;

  Vertex1->GetDatum(FL::DatumAtVertexPlasmaDensity,&PlasmaDensityCurrent);
  Vertex1->GetDatum(FL::DatumAtVertexPrevious::DatumAtVertexPlasmaDensity,&PlasmaDensityOld);

  PlasmaDensityCurrentParticle+=weight1*PlasmaDensityCurrent;
  PlasmaDensityOldParticle+=weight1*PlasmaDensityOld;

  if ((PlasmaDensityCurrentParticle==0.0)||(PlasmaDensityOldParticle==0.0)) {
    dLogP=0.0,dP=0.0,dmu=0.0;
    return;
  }

  #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
  DivVsw=-log(PlasmaDensityCurrentParticle/PlasmaDensityOldParticle)/(AMPS2SWMF::MagneticFieldLineUpdate::LastCouplingTime-AMPS2SWMF::MagneticFieldLineUpdate::LastLastCouplingTime);
  #else
  DivVsw=-log(PlasmaDensityCurrentParticle/PlasmaDensityOld)/dt;
  #endif


  if (isfinite(DivVsw)==false) {
    dLogP=0.0,dP=0.0,dmu=0.0;
    return;
   }

  double mu2=mu*mu;

  if (v>=SpeedOfLight) v=0.99*SpeedOfLight;

  dLogP=-((1.0-mu2)/2.0*(DivVsw-dVz_dz)+mu2*dVz_dz)*dt;
  dP=Relativistic::Speed2Momentum(v,_H__MASS_)*dLogP;

  dmu=((1.0-mu2)/2.0*(v/L+mu*(DivVsw-3.0*dVz_dz)))*dt;
}


int SEP::ParticleMover_He_2011_AJ(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  namespace PB = PIC::ParticleBuffer;
  namespace FL = PIC::FieldLine;

  PIC::ParticleBuffer::byte *ParticleData;
  double mu,AbsB,L,vParallel,vNormal,v,DivAbsB,vParallelInit,vNormalInit;
  double FieldLineCoord;
  int iFieldLine,spec;
  FL::cFieldLineSegment *Segment;

  ParticleData=PB::GetParticleDataPointer(ptr);

  FieldLineCoord=PB::GetFieldLineCoord(ParticleData);
  iFieldLine=PB::GetFieldLineId(ParticleData);
  spec=PB::GetI(ParticleData);

  //velocity is in the frame moving with solar wind
  vParallel=PB::GetVParallel(ParticleData);
  vNormal=PB::GetVNormal(ParticleData);

  vParallelInit=vParallel,vNormalInit=vNormal;

  //determine the segment of the particle location
  Segment=FL::FieldLinesAll[iFieldLine].GetSegment(FieldLineCoord);

  //get the new value of 'mu'
  double D,dD_dmu;

  double mu_init=mu;
  double time_counter=0.0;
  double dt=dtTotal;
  double dmu=0.0,dv;
  double delta;

  bool first_pass_flag=true;
  bool first_transport_coeffcient=true;

  //calculate B and L
  double B[3],B0[3],B1[3], AbsBDeriv;

  FL::FieldLinesAll[iFieldLine].GetMagneticField(B0, (int)FieldLineCoord);
  FL::FieldLinesAll[iFieldLine].GetMagneticField(B,       FieldLineCoord);
  FL::FieldLinesAll[iFieldLine].GetMagneticField(B1, (int)FieldLineCoord+1-1E-7);
  AbsB   = pow(B[0]*B[0] + B[1]*B[1] + B[2]*B[2], 0.5);

  AbsBDeriv = (pow(B1[0]*B1[0] + B1[1]*B1[1] + B1[2]*B1[2], 0.5) -
    pow(B0[0]*B0[0] + B0[1]*B0[1] + B0[2]*B0[2], 0.5)) /  FL::FieldLinesAll[iFieldLine].GetSegmentLength(FieldLineCoord);

  //calculate solarwind velocity,particle velocity and mu in the frame moving with solar wind
  double vSolarWind[3],vSolarWindParallel;

  FL::FieldLinesAll[iFieldLine].GetPlasmaVelocity(vSolarWind,FieldLineCoord);
  vSolarWindParallel=Vector3D::DotProduct(vSolarWind,B)/AbsB;


  v=sqrt(vParallel*vParallel+vNormal*vNormal);

  if (v>=0.99*SpeedOfLight) {
    double t=0.99*SpeedOfLight/v;

    v*=t;
    vParallel*=t;
    vNormal*=t;
  }


/*
if (Relativistic::Speed2E(v,_H__MASS_)>100.0*MeV2J) {
cout << "found" << endl;
}*/


  mu=vParallel/v;

  static long int ncall=0,loop_cnt=0;
  ncall++;

  while (time_counter<dtTotal) {
    if (time_counter+dt>dtTotal) dt=dtTotal-time_counter;

    loop_cnt++;
    dmu=0.0;

    if (_PIC_DEBUGGER_MODE_== _PIC_DEBUGGER_MODE_ON_) {
      if (isfinite(mu)==false) exit(__LINE__,__FILE__);
    }

    int iR,iE,iMu;

    if (SEP::Diffusion::GetPitchAngleDiffusionCoefficient!=NULL) {
      SEP::Diffusion::GetPitchAngleDiffusionCoefficient(D,dD_dmu,mu,vParallel,vNormal,spec,FieldLineCoord,Segment);

      if (SEP::Diffusion::PitchAngleDifferentialMode==SEP::Diffusion::PitchAngleDifferentialModeNumerical) {
        double t,mu_plus,mu_minus,D_plus,D_minus,D_mu_mu_numerical;

        mu_plus=mu+SEP::Diffusion::muNumericalDifferentiationStep;
        if (mu_plus>1.0) mu_plus=1.0;

        mu_minus=mu-SEP::Diffusion::muNumericalDifferentiationStep;
        if (mu_minus<-1.0) mu_minus=-1.0;

        if (mu_plus*mu_minus<0.0) {
           if (fabs(mu_minus)<fabs(mu_plus)) {
             mu_minus=0.0;
           }
           else {
             mu_plus=0.0;
           }
        }

        SEP::Diffusion::GetPitchAngleDiffusionCoefficient(D_plus,t,mu_plus,vParallel,vNormal,spec,FieldLineCoord,Segment);
        SEP::Diffusion::GetPitchAngleDiffusionCoefficient(D_minus,t,mu_minus,vParallel,vNormal,spec,FieldLineCoord,Segment);

        D_mu_mu_numerical=(D_plus-D_minus)/(mu_plus-mu_minus);

        dD_dmu=D_mu_mu_numerical;
      }
      else if (SEP::Diffusion::PitchAngleDifferentialMode!=SEP::Diffusion::PitchAngleDifferentialModeAnalytical) {
        exit(__LINE__,__FILE__,"Error: the option is unknown");
      }

      if (first_pass_flag==true) {
        delta=sqrt(2.0*D*dt)*Vector3D::Distribution::Normal();
        if (isfinite(delta)==false) exit(__LINE__,__FILE__);

        //sample Dmumu
        double x[3],e,speed,ParticleWeight;

        speed=sqrt(vNormal*vNormal+vParallel*vParallel);
        if (speed>0.99*SpeedOfLight) speed=0.99*SpeedOfLight;

        e=Relativistic::Speed2E(speed,PIC::MolecularData::GetMass(spec));
        iE=log(e/SEP::Sampling::PitchAngle::emin)/SEP::Sampling::PitchAngle::dLogE;

        if (iE>=SEP::Sampling::PitchAngle::nEnergySamplingIntervals) iE=SEP::Sampling::PitchAngle::nEnergySamplingIntervals-1;
        if (iE<0)iE=0;

        iMu=(int)((mu+1.0)/SEP::Sampling::PitchAngle::dMu);
        if (iMu>=SEP::Sampling::PitchAngle::nMuIntervals) iMu=SEP::Sampling::PitchAngle::nMuIntervals-1;

        FL::FieldLinesAll[iFieldLine].GetCartesian(x,FieldLineCoord);
        iR=(int)(Vector3D::Length(x)/SEP::Sampling::PitchAngle::dR);
        if (iR>=SEP::Sampling::PitchAngle::nRadiusIntervals) iR=SEP::Sampling::PitchAngle::nRadiusIntervals-1;

        ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
        ParticleWeight*=PB::GetIndividualStatWeightCorrection(ParticleData);

        SEP::Sampling::PitchAngle::DmumuSamplingTable(0,iMu,iE,iR,iFieldLine)+=D*ParticleWeight;
        SEP::Sampling::PitchAngle::DmumuSamplingTable(1,iMu,iE,iR,iFieldLine)+=ParticleWeight;

        if (ModelEquation==ModelEquationParker) {
          return ParticleMover_ParkerEquation(ptr,dtTotal,node);
        }

        if (fabs(dD_dmu*dt)>0.05) {
          dt=0.05/fabs(dD_dmu);

          if (dt/dtTotal<TimeStepRatioSwitch_FTE2PE) {
            return ParticleMover_ParkerEquation(ptr,dtTotal,node);
          }
        }

        if (sqrt(2.0*D*dt)>0.05) { // fabs(delta)>0.1) {
          double t=dt*pow(0.05/fabs(delta),2);

          if (t<dt) dt=t;

          if (dt/dtTotal<TimeStepRatioSwitch_FTE2PE) {
            return ParticleMover_ParkerEquation(ptr,dtTotal,node);
          }
        }

        first_pass_flag=false;
      }

      delta=sqrt(2.0*D*dt)*Vector3D::Distribution::Normal();
      dmu+=delta;
      dmu+=dD_dmu*dt; //IMPORTANT. It is actually should be '+' [Dresing, 2012 Arxive; Droge-2009-AJ]

      mu+=dmu;
      dmu=0.0;

  if (mu>1.0) {
    double d=mu-1.0;
//    mu=1.0-d;

    mu=(mu>=1.01) ? -1.0+2*rnd() : 1.0-d;
  }
  else if (mu<-1.0) {
    double d=1.0+mu;
    //mu=-1.0-d;

    mu=(mu<-1.01) ? -1.0+2*rnd() : -1.0-d;
  }

    }


    double dlogp,dp;
    if (v>=SpeedOfLight) v=0.99*SpeedOfLight;

    if (AccountTransportCoefficient==true) {
      GetTransportCoefficients(dp,dlogp,dmu,v,mu,Segment,FieldLineCoord,dt,iFieldLine,vSolarWindParallel);
    }
    else {
      dp=0.0,dlogp=0.0;
    }


    FieldLineCoord=FL::FieldLinesAll[iFieldLine].move(FieldLineCoord,dt*vParallel);


    if ((first_transport_coeffcient==true)&&(SEP::Diffusion::GetPitchAngleDiffusionCoefficient!=NULL)) {
      first_transport_coeffcient=false;

      double ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
      ParticleWeight*=PB::GetIndividualStatWeightCorrection(ParticleData);

      SEP::Sampling::PitchAngle::DmumuSamplingTable(2,iMu,iE,iR,iFieldLine)+=dmu/dt*ParticleWeight;
      SEP::Sampling::PitchAngle::DmumuSamplingTable(3,iMu,iE,iR,iFieldLine)+=dp/dt*ParticleWeight;
    }


    if (_PIC_DEBUGGER_MODE_== _PIC_DEBUGGER_MODE_ON_) {
      if ((isfinite(dp)==false)||(isfinite(dlogp)==false)) {
        exit(__LINE__,__FILE__);
      }
    }

    double p=Relativistic::Speed2Momentum(v,_H__MASS_);

    p*=exp(dlogp);
    v=Relativistic::Momentum2Speed(p,_H__MASS_);

    if (v>=0.99*SpeedOfLight) {
      v=0.99*SpeedOfLight;
    }

    vParallel=mu*v;
    vNormal=sqrt(1.0-mu*mu)*v;

    if (_PIC_DEBUGGER_MODE_== _PIC_DEBUGGER_MODE_ON_) {
      if ((isfinite(mu)==false)||(isfinite(v)==false)) {
        exit(__LINE__,__FILE__);
      }
    }

    //get the segment of the new particle location
    if ((FieldLineCoord<0.0) || ((Segment=FL::FieldLinesAll[iFieldLine].GetSegment(FieldLineCoord))==NULL)) {
      //the particle left the computational domain
      int code=_PARTICLE_DELETED_ON_THE_FACE_;

      //call the function that process particles that leaved the coputational domain
      switch (code) {
      case _PARTICLE_DELETED_ON_THE_FACE_:
        PIC::ParticleBuffer::DeleteParticle(ptr);
        return _PARTICLE_LEFT_THE_DOMAIN_;

      default:
        exit(__LINE__,__FILE__,"Error: not implemented");
      }
    }

    time_counter+=dt;
  }


  //set the new values of the normal and parallel particle velocities
  vParallel=mu*v;
  vNormal=sqrt(1.0-mu*mu)*v;

  PB::SetVParallel(vParallel,ParticleData);
  PB::SetVNormal(vNormal,ParticleData);

  //set the new particle coordinate
  PB::SetFieldLineCoord(FieldLineCoord,ParticleData);

  //attach the particle to the temporaty list
  // Step 5 guarantees segment attachment at compile time, so the mover can
  // commit directly to the destination segment without a mesh-cell branch.
  {
    long int temp=Segment->tempFirstParticleIndex.exchange(ptr);
      PIC::ParticleBuffer::SetNext(temp,ParticleData);
      PIC::ParticleBuffer::SetPrev(-1,ParticleData);

      if (temp!=-1) PIC::ParticleBuffer::SetPrev(ptr,temp);
  }

  return _PARTICLE_MOTION_FINISHED_;
}


//=============================================================================================================
//Sokolov-2004-AJ
int SEP::ParticleMover_ParkerEquation(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  namespace PB = PIC::ParticleBuffer;
  namespace FL = PIC::FieldLine;

  PIC::ParticleBuffer::byte *ParticleData;
  double mu,AbsB,L,vParallel,vNormal,v,DivAbsB,vParallelInit,vNormalInit;
  double FieldLineCoord;
  int iFieldLine,spec;
  FL::cFieldLineSegment *Segment;

  ParticleData=PB::GetParticleDataPointer(ptr);

  FieldLineCoord=PB::GetFieldLineCoord(ParticleData);
  iFieldLine=PB::GetFieldLineId(ParticleData);
  spec=PB::GetI(ParticleData);

  //velocity is in the frame moving with solar wind
  vParallel=PB::GetVParallel(ParticleData);
  vNormal=PB::GetVNormal(ParticleData);

  v=sqrt(vParallel*vParallel+vNormal*vNormal);

  vParallelInit=vParallel,vNormalInit=vNormal;

  //determine the segment of the particle location
  Segment=FL::FieldLinesAll[iFieldLine].GetSegment(FieldLineCoord);


  auto GetTransportCoefficients = [&] (double& dLogP,double v,FL::cFieldLineSegment *Segment,double FieldLineCoord,double dt,int iFieldLine,double& vSolarWindParallel) {
    //calculate div(vSW) : Dln(Rho)=-div(vSW)*dt
    double PlasmaDensityCurrent,PlasmaDensityOld,DivVsw,PlasmaDensityCurrentParticle=0.0,PlasmaDensityOldParticle=0.0;
    auto Vertex0=Segment->GetBegin();
    auto Vertex1=Segment->GetEnd();

    double weight0=1.0-(FieldLineCoord-floor(FieldLineCoord));
    double weight1=1.0-weight0;

    Vertex0->GetDatum(FL::DatumAtVertexPlasmaDensity,&PlasmaDensityCurrent);
    Vertex0->GetDatum(FL::DatumAtVertexPrevious::DatumAtVertexPlasmaDensity,&PlasmaDensityOld);

    PlasmaDensityCurrentParticle=weight0*PlasmaDensityCurrent;
    PlasmaDensityOldParticle=weight0*PlasmaDensityOld;

    Vertex1->GetDatum(FL::DatumAtVertexPlasmaDensity,&PlasmaDensityCurrent);
    Vertex1->GetDatum(FL::DatumAtVertexPrevious::DatumAtVertexPlasmaDensity,&PlasmaDensityOld);

    PlasmaDensityCurrentParticle+=weight1*PlasmaDensityCurrent;
    PlasmaDensityOldParticle+=weight1*PlasmaDensityOld;

    if ((PlasmaDensityCurrentParticle==0.0)||(PlasmaDensityOldParticle==0.0)) {
      dLogP=0.0;
      return;
    }

    #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
    DivVsw=-log(PlasmaDensityCurrentParticle/PlasmaDensityOldParticle)/(AMPS2SWMF::MagneticFieldLineUpdate::LastCouplingTime-AMPS2SWMF::MagneticFieldLineUpdate::LastLastCouplingTime);
    #else
    DivVsw=-log(PlasmaDensityCurrentParticle/PlasmaDensityOld)/dtTotal;
    #endif



    if (isfinite(DivVsw)==false) {
      dLogP=0.0;
      return;
    }

    dLogP=-(DivVsw)*dt/3.0;
  };


  //get the new value of 'mu'
  double D;
  double mu_init=mu;
  double time_counter=0.0;
  double dt=dtTotal;
  double dmu=0.0,dv;
  double delta,vSolarWindParallel;

  bool first_pass_flag=true;
  bool first_transport_coeffcient=true;

  v=sqrt(vParallel*vParallel+vNormal*vNormal);

  if (v>=SpeedOfLight) {
    double t=0.99*SpeedOfLight/v;

    v*=t;
    vParallel*=t;
    vNormal*=t;
  }

  mu=vParallel/v;

  static long int ncall=0;
  ncall++;

  if (_PIC_DEBUGGER_MODE_== _PIC_DEBUGGER_MODE_ON_) {
    if (isfinite(mu)==false) exit(__LINE__,__FILE__);
  }

  double dx,dDxx_dx;

  dx=0.0;

  while (time_counter<dtTotal) {
    if (time_counter+dt>dtTotal) dt=dtTotal-time_counter;

    if (SEP::Diffusion::GetPitchAngleDiffusionCoefficient!=NULL) {
      SEP::Diffusion::GetDxx(D,dDxx_dx,v,spec,FieldLineCoord,Segment,iFieldLine);
      delta=sqrt(2.0*D*dt)*Vector3D::Distribution::Normal();

      if (_PIC_DEBUGGER_MODE_== _PIC_DEBUGGER_MODE_ON_) {
        if (isfinite(delta)==false) exit(__LINE__,__FILE__);
      }

      if (first_pass_flag==true) {
        if (fabs(dDxx_dx*dt)>0.5*Segment->GetLength()) {
          double dt_new=0.5*Segment->GetLength()/fabs(dDxx_dx);

          delta*=sqrt(dt_new/dt);
          dt=dt_new;
        }

        if (fabs(delta)>0.5*Segment->GetLength()) {
          double dt_new=dt*pow(0.5*Segment->GetLength()/delta,2);

          if (dt<dt_new) dt=dt_new;
        }

        first_pass_flag=false;
      }

      delta=sqrt(2.0*D*dt)*Vector3D::Distribution::Normal();
      dx+=delta;
      dx-=dDxx_dx*dt;
    }


    double dp,dlogp;

    if (v>=SpeedOfLight) v=0.99*SpeedOfLight;

    GetTransportCoefficients(dlogp,v,Segment,FieldLineCoord,dt,iFieldLine,vSolarWindParallel);
    FieldLineCoord=FL::FieldLinesAll[iFieldLine].move(FieldLineCoord,dx);

    double p=Relativistic::Speed2Momentum(v,_H__MASS_);
    p*=exp(dlogp);

    v=Relativistic::Momentum2Speed(p,_H__MASS_);

    if (v>=0.99*SpeedOfLight) {
      v=0.99*SpeedOfLight;
    }

    vParallel=v;
    vNormal=0.0;

    if ((isfinite(mu)==false)||(isfinite(v)==false)) {
      exit(__LINE__,__FILE__);
    }

    //get the segment of the new particle location
    if ((FieldLineCoord<0.0) || ((Segment=FL::FieldLinesAll[iFieldLine].GetSegment(FieldLineCoord))==NULL)) {
      //the particle left the computational domain
      int code=_PARTICLE_DELETED_ON_THE_FACE_;

      //call the function that process particles that leaved the coputational domain
      switch (code) {
      case _PARTICLE_DELETED_ON_THE_FACE_:
        PIC::ParticleBuffer::DeleteParticle(ptr);
        return _PARTICLE_LEFT_THE_DOMAIN_;

      default:
        exit(__LINE__,__FILE__,"Error: not implemented");
      }
    }

    time_counter+=dt;
  }


  //set the new values of the normal and parallel particle velocities
  mu=-1.0+2.0*rnd();
  vParallel=mu*v;
  vNormal=sqrt(1.0-mu*mu)*v;

  PB::SetVParallel(vParallel,ParticleData);
  PB::SetVNormal(vNormal,ParticleData);

  //set the new particle coordinate
  PB::SetFieldLineCoord(FieldLineCoord,ParticleData);

  /*
  //determine the final location of the particle in 3D
  double xFinal[3];
  FL::FieldLinesAll[iFieldLine].GetCartesian(xFinal,FieldLineCoord);

  if (Vector3D::Length(xFinal)>=_AU_) {
    PIC::ParticleBuffer::DeleteParticle(ptr);
    return _PARTICLE_LEFT_THE_DOMAIN_;
  }
  */


  // Field-line attachment is the only supported production representation.
  {
    PIC::ParticleBuffer::SetPrev(-1,ParticleData);

    long int tempFirstParticleIndex;

    tempFirstParticleIndex=atomic_exchange(&Segment->tempFirstParticleIndex,ptr);
    if (tempFirstParticleIndex!=-1) PIC::ParticleBuffer::SetPrev(ptr,tempFirstParticleIndex);
    PIC::ParticleBuffer::SetNext(tempFirstParticleIndex,ParticleData);
  }

  return _PARTICLE_MOTION_FINISHED_;
}

//=========================================================================================================
int SEP::ParticleMover_MeanFreePathScattering(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  namespace PB = PIC::ParticleBuffer;
  namespace FL = PIC::FieldLine;

  PIC::ParticleBuffer::byte *ParticleData;
  double mu,AbsB,L,vParallel,vNormal,v,DivAbsB,vParallelInit,vNormalInit;
  double FieldLineCoord;
  int iFieldLine,spec;
  FL::cFieldLineSegment *Segment;

  ParticleData=PB::GetParticleDataPointer(ptr);

  FieldLineCoord=PB::GetFieldLineCoord(ParticleData);
  iFieldLine=PB::GetFieldLineId(ParticleData);
  spec=PB::GetI(ParticleData);

  //velocity is in the frame moving with solar wind
  vParallel=PB::GetVParallel(ParticleData);
  vNormal=PB::GetVNormal(ParticleData);

  vParallelInit=vParallel,vNormalInit=vNormal;

  //determine the segment of the particle location
  Segment=FL::FieldLinesAll[iFieldLine].GetSegment(FieldLineCoord);





  //get the new value of 'mu'
  double D,dD_dmu;

  double mu_init=mu;
  double time_counter=0.0;
  double dt=dtTotal;
  double dmu=0.0,dv;
  double delta,vSolarWindParallel;

  bool first_pass_flag=true;
  bool first_transport_coeffcient=true;

  v=sqrt(vParallel*vParallel+vNormal*vNormal);

if (v>=SpeedOfLight) {
double t=0.99*SpeedOfLight/v;

v*=t;
vParallel*=t;
vNormal*=t;
}

  mu=vParallel/v;

static long int ncall=0;

ncall++;

if (ncall==2369588) {
double d33=0.0;

d33+=34;
cout << d33 << endl;
}



  while (time_counter<dtTotal) {
    if (time_counter+dt>dtTotal) dt=dtTotal-time_counter;

    dmu=0.0;


if (isfinite(mu)==false) exit(__LINE__,__FILE__);


int iR,iE,iMu;





double dp,dlogp;



if (v>=SpeedOfLight) v=0.99*SpeedOfLight;


double MeanFreePath=SEP::Diffusion::GetMeanFreePath(v,spec,FieldLineCoord,Segment,iFieldLine);

if (rnd()<1.0-exp(-dt*v/MeanFreePath)) {
  //scattering occured
  time_counter+=dt;

  //determine the new location of the particle
  //1. determine the time before the scattering and push the particle forward
  double MaxPathLength=v*dt*fabs(mu);
  double PathLength;
  double dt_before_scattering;

  if (MaxPathLength>0.0) {
    PathLength=-MeanFreePath*log(1.0-rnd()*(1.0-exp(-MaxPathLength/MeanFreePath)));
    dt_before_scattering=PathLength/(v*fabs(mu));

    FieldLineCoord=FL::FieldLinesAll[iFieldLine].move(
        FieldLineCoord,dt_before_scattering*(vParallel+vSolarWindParallel));

    if ((FieldLineCoord<0.0) || ((Segment=FL::FieldLinesAll[iFieldLine].GetSegment(FieldLineCoord))==NULL)) {
      //the particle left the computational domain
      int code=_PARTICLE_DELETED_ON_THE_FACE_;

      //call the function that process particles that leaved the coputational domain
      switch (code) {
      case _PARTICLE_DELETED_ON_THE_FACE_:
        PIC::ParticleBuffer::DeleteParticle(ptr);
        return _PARTICLE_LEFT_THE_DOMAIN_;

      default:
        exit(__LINE__,__FILE__,"Error: not implemented");
      }
    }

    //push the particle after scattering
    mu=-1.0+rnd()*2.0;
    vParallel=mu*v;
    vNormal=v*sqrt(1.0-mu*mu);

    FieldLineCoord=FL::FieldLinesAll[iFieldLine].move(
        FieldLineCoord,(vParallel+vSolarWindParallel)*(dt-dt_before_scattering));

    if ((FieldLineCoord<0.0) || ((Segment=FL::FieldLinesAll[iFieldLine].GetSegment(FieldLineCoord))==NULL)) {
      //the particle left the computational domain
      int code=_PARTICLE_DELETED_ON_THE_FACE_;

      //call the function that process particles that leaved the coputational domain
      switch (code) {
      case _PARTICLE_DELETED_ON_THE_FACE_:
        PIC::ParticleBuffer::DeleteParticle(ptr);
        return _PARTICLE_LEFT_THE_DOMAIN_;

      default:
        exit(__LINE__,__FILE__,"Error: not implemented");
      }
    }


    continue;
  }
}



    GetTransportCoefficients(dp,dlogp,dmu,v,mu,Segment,FieldLineCoord,dt,iFieldLine,vSolarWindParallel);
    FieldLineCoord=FL::FieldLinesAll[iFieldLine].move(FieldLineCoord,dt*(vParallel+vSolarWindParallel));


/*    if (first_transport_coeffcient==true) {
      first_transport_coeffcient=false;

      double ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
      ParticleWeight*=PB::GetIndividualStatWeightCorrection(ParticleData);

      SEP::Sampling::PitchAngle::DmumuSamplingTable(2,iMu,iE,iR,iFieldLine)+=dmu/dt*ParticleWeight;
      SEP::Sampling::PitchAngle::DmumuSamplingTable(3,iMu,iE,iR,iFieldLine)+=dp/dt*ParticleWeight;
    }*/


double p=Relativistic::Speed2Momentum(v,_H__MASS_);
//p+=dp;

p*=exp(dlogp);

v=Relativistic::Momentum2Speed(p,_H__MASS_);

if (v>=0.99*SpeedOfLight) {
v=0.99*SpeedOfLight;
}

//    v+=dv;
    mu+=dmu;
    dmu=0.0;


      if (mu>0.999) mu=0.999;
      if (mu<-0.999) mu=-0.999;

  vParallel=mu*v;
  vNormal=sqrt(1.0-mu*mu)*v;

if ((isfinite(mu)==false)||(isfinite(v)==false)) {
  exit(__LINE__,__FILE__);
}

    //get the segment of the new particle location
    if ((FieldLineCoord<0.0) || ((Segment=FL::FieldLinesAll[iFieldLine].GetSegment(FieldLineCoord))==NULL)) {
      //the particle left the computational domain
      int code=_PARTICLE_DELETED_ON_THE_FACE_;

      //call the function that process particles that leaved the coputational domain
      switch (code) {
      case _PARTICLE_DELETED_ON_THE_FACE_:
        PIC::ParticleBuffer::DeleteParticle(ptr);
        return _PARTICLE_LEFT_THE_DOMAIN_;

      default:
        exit(__LINE__,__FILE__,"Error: not implemented");
      }
    }

    time_counter+=dt;
  }


  //set the new values of the normal and parallel particle velocities
  vParallel=mu*v;
  vNormal=sqrt(1.0-mu*mu)*v;

  PB::SetVParallel(vParallel,ParticleData);
  PB::SetVNormal(vNormal,ParticleData);

  //set the new particle coordinate
  PB::SetFieldLineCoord(FieldLineCoord,ParticleData);

  // Field-line attachment is the only supported production representation.
    {
    PIC::ParticleBuffer::SetPrev(-1,ParticleData);

    long int tempFirstParticleIndex;

    tempFirstParticleIndex=atomic_exchange(&Segment->tempFirstParticleIndex,ptr);
    if (tempFirstParticleIndex!=-1) PIC::ParticleBuffer::SetPrev(ptr,tempFirstParticleIndex);
    PIC::ParticleBuffer::SetNext(tempFirstParticleIndex,ParticleData);
  }

  return _PARTICLE_MOTION_FINISHED_;
}
