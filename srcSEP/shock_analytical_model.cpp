#include "sep.h"
//analytic model of a shock wave (Tenishev-2005-AIAA-4928

double SEP::ParticleSource::ShockWave::Tenishev2005::rShock = 0.0;

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

  return res;
}

void SEP::ParticleSource::ShockWave::Tenishev2005::UpdateShockLocation() {
  static double simulation_time_last = 0.0;

  if (simulation_time_last != PIC::SimulationTime::TimeCounter) {
    double speed = GetShockSpeed();
    rShock += speed * (PIC::SimulationTime::TimeCounter - simulation_time_last);
  }
}

double SEP::ParticleSource::ShockWave::Tenishev2005::GetInjectionRate() {
  static bool init_flag;
  static double rmin;
  double res;

  if (init_flag == false) {
    //find the valiue of rmin
    init_flag = true;

    rmin = -1.0;

    //loop through all simulated field lines
    for (int i = 0; i < PIC::FieldLine::nFieldLine; i++) {
      double r, *x;

      x = PIC::FieldLine::FieldLinesAll[i].GetFirstSegment()->GetBegin()->GetX();
      r = Vector3D::Length(x);

      if ((rmin < 0.0) || (rmin > r))
        rmin = r;
    }
  }

  return pow(rmin / rShock, 2);
}

PIC::FieldLine::cFieldLineSegment* SEP::ParticleSource::ShockWave::Tenishev2005::GetInjectionLocation(int iFieldLine) {
  PIC::FieldLine::cFieldLineSegment *Segment =
      PIC::FieldLine::FieldLinesAll[iFieldLine].GetFirstSegment();
  double r, *x, r2 = rShock * rShock;

  UpdateShockLocation();
  x=Segment->GetBegin()->GetX(); 

  if (r2 <= Vector3D::DotProduct(x, x))
    return NULL;

  for (; Segment != NULL; Segment = Segment->GetNext()) {
    x = Segment->GetBegin()->GetX();

    if (Vector3D::DotProduct(x, x) < r2)
      return Segment;
  }

  return NULL;

}

