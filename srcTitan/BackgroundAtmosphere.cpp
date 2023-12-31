//====================================================
//$Id$
//====================================================
//define specific functions for modeling collisions with the backgound atmosphere

#include <iostream>
#include <fstream>
#include <sstream>
#include "pic.h"


//collison cross section of the model particles with the background atmosphere
double PIC::MolecularCollisions::BackgroundAtmosphere::GetCollisionCrossSectionBackgoundAtmosphereParticle(int spec,int BackgroundSpecieNumber,PIC::ParticleBuffer::byte *modelParticleData,PIC::ParticleBuffer::byte *BackgroundAtmosphereParticleData,double TranslationalEnergy,double cr2) {
  /* const static double BackgroundSpecieDiameter[2]={1.95E-10,4.0E-10};

  const static double RefDiameter_O=2.52E-10;

  return Pi*pow(0.5*(RefDiameter_O+BackgroundSpecieDiameter[BackgroundSpecieNumber]),2);*/

  return 3.0E-19;
}


double PIC::MolecularCollisions::BackgroundAtmosphere::GetSigmaCrMax(int spec,int BackgroundSpecieNumber,PIC::ParticleBuffer::byte *modelParticleData) {
  return 10000.0*3.0E-19;
}

//distribute the direction of the relative velocity of particles after collision
double PIC::MolecularCollisions::BackgroundAtmosphere::GetCollisionScatteringAngle(double* Vrel,double TranslationalEnergy,int spec,int BackgroundSpecieNumber) {

  exit(__LINE__,__FILE__,"Error: not implemented");

  return 0.0;
}


double GetBackgroundMolecularMass(int BackgroundSpecieNumber) {
  double res;

  switch (BackgroundSpecieNumber) {
  case _N2_BACKGROUND_SPEC_ :
    res=_MASS_(_N2_);
    break;
  default:
    exit(__LINE__,__FILE__,"BackgroundSpecieNumber is out of range");
  }

  return res;
}

void PIC::MolecularCollisions::BackgroundAtmosphere::GenerateBackgoundAtmosphereParticle(PIC::ParticleBuffer::byte *BackgroundAtmosphereParticleData,int BackgroundSpecieNumber,PIC::Mesh::cDataCenterNode *cell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  int idim;
  double *xmin,*xmax,*xMiddle,x[3],v[3],beta;

  static const double GlobalNeutalTemperature=147.0;

  //generate positions of the background particle in the cell
  xmin=node->xmin;
  xmax=node->xmax;
  xMiddle=cell->GetX();

  x[0]=xMiddle[0]+(xmax[0]-xmin[0])/_BLOCK_CELLS_X_*(rnd()-0.5);
  x[1]=xMiddle[1]+(xmax[1]-xmin[1])/_BLOCK_CELLS_Y_*(rnd()-0.5);
  x[2]=xMiddle[2]+(xmax[2]-xmin[2])/_BLOCK_CELLS_Z_*(rnd()-0.5);

  PIC::ParticleBuffer::SetX(x,BackgroundAtmosphereParticleData);

  //generate velocity vector for a particle representing the bacground atmosphere
  beta=GetBackgroundMolecularMass(BackgroundSpecieNumber)/(2*Kbol*GlobalNeutalTemperature);

  for (idim=0;idim<3;idim++) v[idim]=cos(2*Pi*rnd())*sqrt(-log(rnd())/beta);

  PIC::ParticleBuffer::SetV(v,BackgroundAtmosphereParticleData);
  PIC::ParticleBuffer::SetI(BackgroundSpecieNumber,BackgroundAtmosphereParticleData);
}

double PIC::MolecularCollisions::BackgroundAtmosphere::GetBackgroundLocalNumberDensity(int BackgroundSpecieNumber,double *x) {
  double r,res=0.0;

  static const double r0=_RADIUS_(_TITAN_);   3788.921E3;
  static const double n0=1.51e15;
  static const double T=147.0;

  r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  res=(r>r0) ? n0*exp(-(r-r0) / ((Kbol*T*r*r0)/(GravityConstant*_MASS_(_TITAN_)*GetBackgroundMolecularMass(BackgroundSpecieNumber))) ) : n0; //   n0 exp [-(r – r0)/{ k T r r0/(G MT mN2)}]


  return res;
}


double PIC::MolecularCollisions::BackgroundAtmosphere::GetCellMeanBackgroundNumberDensity(int BackgroundSpecieNumber,PIC::Mesh::cDataCenterNode *cell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  double x[3];
  cell->GetX(x);

  return GetBackgroundLocalNumberDensity(BackgroundSpecieNumber,x);
}

double PIC::MolecularCollisions::BackgroundAtmosphere::GetCellMaximumBackgroundNumberDensity(int BackgroundSpecieNumber,PIC::Mesh::cDataCenterNode *cell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  return GetCellMeanBackgroundNumberDensity(BackgroundSpecieNumber,cell,node);
}

double PIC::MolecularCollisions::BackgroundAtmosphere::GetCellLocalBackgroundNumberDensity(double x[3],int BackgroundSpecieNumber,PIC::Mesh::cDataCenterNode *cell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  return GetBackgroundLocalNumberDensity(BackgroundSpecieNumber,x);
}

int GetTotalNumberBackgroundSpecies() {
  return 1;
}

//acceptance croterion for the simulated ions
bool PIC::MolecularCollisions::BackgroundAtmosphere::KeepConditionModelParticle(PIC::ParticleBuffer::byte *ModelParticleData) {
  double *v;

  static double const VelocityLimit=100.0;

  v=PIC::ParticleBuffer::GetV(ModelParticleData);

  return (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]>VelocityLimit*VelocityLimit) ? true : false;
}

//calculate the stopping power
double PIC::MolecularCollisions::StoppingPowerModel::GetStoppingPower(double *x,double *v,int spec,PIC::Mesh::cDataCenterNode *cell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  double mass,charge;

  struct cStoppingPowerData {
    double B,T0,q,s;
  };

  static cStoppingPowerData SeParam={1.3798e-16,0.05,0.92,0.6};
  static cStoppingPowerData SnParamLow={2.665e-11,4700.0,0.99,1.68};
  static cStoppingPowerData SnParamHigh={1.57e-11,4000,0.97,1.65};
  static double EaEnergyThrehold=1000.0;

  mass=PIC::MolecularData::GetMass(spec);
  charge=PIC::MolecularData::GetElectricCharge(spec);

  if (charge==0.0) return 0.0;

  double Ea,nu,Se,Sn;

  Ea=mass*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])/2.0;
  Ea = Ea/eV2J;

  nu=(4.0*mass*_N__MASS_/(pow(mass+_N__MASS_,2.0)));

  Se=(SeParam.B/pow(Ea,SeParam.s))*
      (pow(SeParam.T0,2.0-SeParam.q)-pow(SeParam.T0+nu*Ea,1.0-SeParam.q)*
      (SeParam.T0+(-1.0+SeParam.q)*nu*Ea))/
      ((1.0-SeParam.q)*(2.0-SeParam.q));

  if (Ea<EaEnergyThrehold) {
    Sn=(SnParamLow.B/pow(Ea,SnParamLow.s))*
        (pow(SnParamLow.T0,2.0-SnParamLow.q)-pow(SnParamLow.T0+nu*Ea,1.0-SnParamLow.q)*
        (SnParamLow.T0+(-1.0+SnParamLow.q)*nu*Ea))/
        ((1.0-SnParamLow.q)*(2.0-SnParamLow.q));
  }
  else {
    Sn=(SnParamHigh.B/pow(Ea,SnParamHigh.s))*
        (pow(SnParamHigh.T0,2.0-SnParamHigh.q)-pow(SnParamHigh.T0+nu*Ea,1.0-SnParamHigh.q)*
        (SnParamHigh.T0+(-1.0+SnParamHigh.q)*nu*Ea))/
        ((1.0-SnParamHigh.q)*(2.0-SnParamHigh.q));
  }

  //get the totan background species number density
  double nTotalNumberDensity=0.0;
  for (int nBackgroundSpec=0;nBackgroundSpec<GetTotalNumberBackgroundSpecies();nBackgroundSpec++) {
    nTotalNumberDensity+=PIC::MolecularCollisions::BackgroundAtmosphere::GetBackgroundLocalNumberDensity(nBackgroundSpec,x);
  }
 
  return 2.0*(Se+Sn)*eV2J*1.0e-4*nTotalNumberDensity;
}
