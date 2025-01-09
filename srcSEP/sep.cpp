
#include "sep.h"

int SEP::Offset::Momentum=-1;
int SEP::Offset::CosPitchAngle=-1;
int SEP::Offset::p_par=-1;
int SEP::Offset::p_norm=-1;
int SEP::Offset::RadialLocation=-1;
int SEP::Offset::MeanFreePath=-1;

//in the case the model is run as a part of the SWMF, FreezeTimeSimulationMHD  is the sumulation time starting which the control of the
//model run is not returned to the SWMF and the sumulation continues with AMPS only and "freezed" MHD solar wind
double SEP::FreezeSolarWindModelTime=-1.0;

//composition table of the GCR composition
cCompositionGroupTable *SEP::CompositionGroupTable=NULL;
int *SEP::CompositionGroupTableIndex=NULL;
int SEP::nCompositionGroups=0;

cInternalSphericalData* SEP::InnerBoundary=NULL;

//the type of the equations that is soleved 
int SEP::ModelEquation=SEP::ModelEquationFTE;

//types of the differentiation of the pitch angle diffusion coeffcient
int SEP::Diffusion::PitchAngleDifferentialMode=SEP::Diffusion::PitchAngleDifferentialModeAnalytical; 

//parameters of the scattering model
int SEP::Scattering::Tenishev2005AIAA::status=SEP::Scattering::Tenishev2005AIAA::_disabled;
double SEP::Scattering::Tenishev2005AIAA::alpha=0.0;
double SEP::Scattering::Tenishev2005AIAA::beta=0.0;
double SEP::Scattering::Tenishev2005AIAA::lambda0=0.4*_AU_;

//the limit to switch from solving FTE to the Parker Equation when the D_{\mu\mu} is to high
double SEP::TimeStepRatioSwitch_FTE2PE=-1.0;

//min/max particle number limit during a run 
int SEP::MinParticleLimit=10,SEP::MaxParticleLimit=20;

//title that will be printed inn Tecplot output file (simuation time)
void SEP::TecplotFileTitle(char* title) {
  sprintf(title,"time=%e",PIC::SimulationTime::Get());
}  


void SEP::Init() {
  //title that will be printed inn Tecplot output file (simuation time)
  PIC::FieldLine::UserDefinedTecplotFileTitle=TecplotFileTitle;

  //composition of the GCRs
  nCompositionGroups=1;
  CompositionGroupTable=new cCompositionGroupTable[nCompositionGroups];
  CompositionGroupTableIndex=new int[PIC::nTotalSpecies];

  for (int spec=0;spec<PIC::nTotalSpecies;spec++) CompositionGroupTableIndex[spec]=0; //all simulated model species are hydrogen
  CompositionGroupTable[0].FistGroupSpeciesNumber=0;
  CompositionGroupTable[0].nModelSpeciesGroup=PIC::nTotalSpecies;

  CompositionGroupTable[0].minVelocity=Relativistic::E2Speed(SEP::BoundingBoxInjection::minEnergy,PIC::MolecularData::GetMass(0));
  CompositionGroupTable[0].maxVelocity=Relativistic::E2Speed(SEP::BoundingBoxInjection::maxEnergy,PIC::MolecularData::GetMass(0));

  CompositionGroupTable[0].GroupVelocityStep=(CompositionGroupTable[0].maxVelocity-CompositionGroupTable[0].minVelocity)/CompositionGroupTable[0].nModelSpeciesGroup;

  if (PIC::ThisThread==0) {
    cout << "$PREFIX: Composition Group Velocity and Energy Characteristics:\nspec\tmin Velocity [m/s]\tmax Velovity[m/s]\t min Energy[eV]\tmax Energy[eV]" << endl;

    for (int s=0;s<PIC::nTotalSpecies;s++) {
      double minV,maxV,minE,maxE,mass;

      mass=PIC::MolecularData::GetMass(s);

      minV=::SEP::CompositionGroupTable[0].GetMinVelocity(s);
      maxV=::SEP::CompositionGroupTable[0].GetMaxVelocity(s);

      //convert velocity into energy and distribute energy of a new particles
      minE=Relativistic::Speed2E(minV,mass);
      maxE=Relativistic::Speed2E(maxV,mass);
      
      cout << s << "\t" << minV << "\t" << maxV << "\t" << minE*J2eV << "\t" <<  maxE*J2eV << endl;
    }
  }
  
  //init source models of SEP and GCR
  if (_PIC_EARTH_GCR__MODE_==_PIC_MODE_ON_) BoundingBoxInjection::GCR::Init();
} 


void SEP::RequestParticleData() {
  long int offset;

  switch (_SEP_MOVER_) {
  case _SEP_MOVER_HE_2019_AJL_:
    PIC::ParticleBuffer::RequestDataStorage(offset,sizeof(double));
    Offset::Momentum=offset;

    PIC::ParticleBuffer::RequestDataStorage(offset,sizeof(double));
    Offset::CosPitchAngle=offset;
    break;
  case _SEP_MOVER_BOROVIKOV_2019_ARXIV_:
    PIC::ParticleBuffer::RequestDataStorage(offset,sizeof(double));
    Offset::p_par=offset;

    PIC::ParticleBuffer::RequestDataStorage(offset,sizeof(double));
    Offset::p_norm=offset;
    break;
  }


  //request the memory to store particle's distance from the magnetic field line 
  PIC::ParticleBuffer::RequestDataStorage(offset,sizeof(double));
  Offset::RadialLocation=offset;   

  //request the memory to store particle's mean free path for sample  
  PIC::ParticleBuffer::RequestDataStorage(offset,sizeof(double));
  Offset::MeanFreePath=offset; 
}
