#include "sep.h"
//analytic model of a shock wave (Tenishev-2005-AIAA-4928

double SEP::ParticleSource::ShockWave::Tenishev2005::rShock = 0.0;
bool SEP::ParticleSource::ShockWave::Tenishev2005::InitFlag=false;
double SEP::ParticleSource::ShockWave::Tenishev2005::MinFieldLineHeliocentricDistance=-1.0;  


void SEP::ParticleSource::ShockWave::Tenishev2005::Init() {
  InitFlag=true;

  //determine the  initial location of the shock that is the minimum helpocentric distance of the beginning of the simulated field lines  
  //loop through all simulated field lines
  for (int i = 0; i < PIC::FieldLine::nFieldLine; i++) {
    double r, *x;

    x = PIC::FieldLine::FieldLinesAll[i].GetFirstSegment()->GetBegin()->GetX();
    r = Vector3D::Length(x);

    if ((MinFieldLineHeliocentricDistance < 0.0) || (MinFieldLineHeliocentricDistance > r))
      MinFieldLineHeliocentricDistance = r;
  }

  rShock=MinFieldLineHeliocentricDistance;
}
  

double SEP::ParticleSource::ShockWave::Tenishev2005::GetShockSpeed() {
  double r = rShock / _AU_;
  double res;

  if (r < 0.1)
    res = 1800.0;
  else if (r < 0.15)
    res = 1800.0 + (1500.0 - 1800.0) / (0.15 - 0.1) * (r - 0.1);
  else if (r < 0.3)
    res = 1500.0;
  else if (r < 0.5)
    res = 1500.0 + (1100.0 - 1500.0) / (0.5 - 0.1) * (r - 0.3);
  else if (r < 1.3)
    res = 1100 + (900.0 - 1100.0) / (1.3 - 0.5) * (r - 0.5);
  else
    res = 900.0;

  return res*1.0E3;
}

void SEP::ParticleSource::ShockWave::Tenishev2005::UpdateShockLocation() {
  static double simulation_time_last = 0.0;

  if (InitFlag==false) Init();

  if (simulation_time_last != PIC::SimulationTime::TimeCounter) {
    double speed = GetShockSpeed();
    rShock += speed * (PIC::SimulationTime::TimeCounter - simulation_time_last);
    simulation_time_last=PIC::SimulationTime::TimeCounter;
  }
}

double SEP::ParticleSource::ShockWave::Tenishev2005::GetInjectionRate() {
  double res;

  if (InitFlag==false) Init(); 

  return pow(MinFieldLineHeliocentricDistance/rShock,2);
}

int SEP::ParticleSource::ShockWave::Tenishev2005::GetInjectionLocation(int iFieldLine) {
  double r, *x, r2 = rShock * rShock;
  int iSegment;
  PIC::FieldLine::cFieldLineSegment *Segment =
      PIC::FieldLine::FieldLinesAll[iFieldLine].GetFirstSegment();

  UpdateShockLocation();
  x=Segment->GetBegin()->GetX(); 

  if (r2 <= Vector3D::DotProduct(x, x)) return -1; 

  for (iSegment=0; Segment != NULL; iSegment++,Segment = Segment->GetNext()) {
    x = Segment->GetBegin()->GetX();
    if (Vector3D::DotProduct(x, x) < r2) return iSegment; 
  }

  return -1;
}

