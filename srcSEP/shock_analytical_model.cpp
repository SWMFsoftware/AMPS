#include "sep.h"
//analytic model of a shock wave (Tenishev-2005-AIAA-4928

#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
#include "../srcInterface/amps2swmf.h"
#endif

double SEP::ParticleSource::ShockWave::Tenishev2005::rShock =1.0E-5*_SUN__RADIUS_;
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

  #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
  if (_PIC_FIELD_LINE_MODE_==_PIC_MODE_OFF_) { 
    rShock=AMPS2SWMF::Heliosphere::rMin;
  }
  #endif 
}
  
double SEP::ParticleSource::ShockWave::Tenishev2005::GetCompressionRatio() {
  double r = rShock / _AU_;
  double res;

  if (r<0.04) {
    res=1.7;
  }
  else {
    res=2.0+(1.4-2.0)/(0.14-0.04)*(r-0.04);
    if (res<1.0) res=1.0;
  } 

  return res;
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

double SEP::ParticleSource::ShockWave::Tenishev2005::GetSolarWindDensity() {
  double t=_AU_/rShock;
  return 5.0E6*t*t;
}

double SEP::ParticleSource::ShockWave::Tenishev2005::GetInjectionRate() {
  double s,density,efficientcy,res;

  if (InitFlag==false) Init(); 

  s=SEP::ParticleSource::ShockWave::Tenishev2005::GetCompressionRatio();

  density=GetSolarWindDensity();
  efficientcy=(s-1.0)/s;
  res=density*efficientcy;

  return res;
}

int SEP::ParticleSource::ShockWave::Tenishev2005::GetInjectionLocation(int iFieldLine,double &S,double *xInjection) {
  double r, *x, r2 = rShock * rShock;
  int iSegment;
  PIC::FieldLine::cFieldLineSegment *Segment =
      PIC::FieldLine::FieldLinesAll[iFieldLine].GetFirstSegment();


  UpdateShockLocation();
  x = Segment->GetBegin()->GetX();

  S = -1.0;
  if (r2 <= Vector3D::DotProduct(x, x)) return -1;

  for (iSegment = 0; Segment != NULL; iSegment++, Segment = Segment->GetNext()) {
    double *x1 = Segment->GetBegin()->GetX();
    double r1_2 = Vector3D::DotProduct(x1, x1);

    if (r1_2 < r2) {
        double *x2 = Segment->GetEnd()->GetX();
        double r2_2 = Vector3D::DotProduct(x2, x2);

        if (r2_2 >= r2) {
            // Solve quadratic equation |x1 + t(x2-x1)|^2 = r2
            double d[3];
            double a = 0.0, b = 0.0, c = -r2;

            for (int i = 0; i < 3; i++) {
                d[i] = x2[i] - x1[i];
                a += d[i] * d[i];           // |d|^2
                b += 2.0 * x1[i] * d[i];    // 2(x1·d)
                c += x1[i] * x1[i];         // |x1|^2
            }

            // Solve at^2 + bt + c = 0
            double discriminant = b*b - 4*a*c;

            if (discriminant < 0) {
                // Use segment start point since it's inside shock
                for (int i = 0; i < 3; i++) {
                   xInjection[i] = x1[i];
                }

		S=iSegment; 
                return iSegment;
            }

            // Two real solutions
            double t1 = (-b - sqrt(discriminant))/(2*a);
            double t2 = (-b + sqrt(discriminant))/(2*a);

            // Find valid t
            if (t1 >= 0 && t1 <= 1) {  // t1 is in [0,1]
                for (int i = 0; i < 3; i++) {
                    xInjection[i] = x1[i] + t1 * d[i];
                }

		S=t1+iSegment; 
                return iSegment;
            }
            else if (t2 >= 0 && t2 <= 1) {  // t2 is in [0,1]
                for (int i = 0; i < 3; i++) {
                    xInjection[i] = x1[i] + t2 * d[i];
                }

		S=t2+iSegment; 
                return iSegment;
            }
            else if (t1 < 0 && t2 > 1) {  // Solutions bracket the segment
                for (int i = 0; i < 3; i++) {
                    xInjection[i] = x1[i];
                }

		S=iSegment;
                return iSegment;
            }
        }
    }
}


return -1;
}


