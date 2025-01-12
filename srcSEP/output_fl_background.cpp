#include "sep.h"
#include "amps2swmf.h"

void SEP::FieldLine::OutputBackgroundData(char* fname, int iFieldLine) {
  namespace FL = PIC::FieldLine;
  FILE *fout=fopen(fname,"w");
  double s=0.0,*x0,*x1;
  FL::cFieldLineSegment *Segment=FL::FieldLinesAll[iFieldLine].GetFirstSegment();

  fprintf(fout,"VARIABLES=\"S[AU]\", \"Heliocentric Distance [AU]\", \"AbsB (Parker) [T]\"");

  if (_PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_) {
    fprintf(fout,", \"AbsB [T]\", \"Number Density\", \"Speed\", \"Temp\", \"Pressure\", \"w+\", \"w-\"");
  }

  fprintf(fout,"\n");

  auto OutputVertex = [&] (double s,double r,FL::cFieldLineVertex* v) {
    double AbsB_Parker;

    AbsB_Parker=SEP::ParkerSpiral::GetAbsB(r); 
    fprintf(fout,"%e %e %e ",s/_AU_,r/_AU_,AbsB_Parker);

    if (_PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_) {
      double pv[3],n,t,p,*W,AbsB;

      AbsB=Vector3D::Length(v->GetMagneticField());
      v->GetPlasmaVelocity(pv);
      v->GetPlasmaDensity(n);
      v->GetPlasmaTemperature(t);
      v->GetPlasmaPressure(p);
      W=v->GetDatum_ptr(FL::DatumAtVertexPlasmaWaves);

      fprintf(fout,"%e %e %e %e %e %e %e ", AbsB,n,Vector3D::Length(pv),t,p,W[0],W[1]);
    }

    fprintf(fout,"\n");
  };  

  x0=Segment->GetBegin()->GetX(); 

  for (;Segment!=NULL;Segment=Segment->GetNext()) {
    x1=Segment->GetBegin()->GetX();
    s+=Vector3D::Distance(x0,x1);
    x0=x1;

    OutputVertex(s,Vector3D::Length(x1),Segment->GetBegin());

    if (Segment->GetNext()==NULL) {
       x1=Segment->GetEnd()->GetX();
       s+=Vector3D::Distance(x0,x1);
       OutputVertex(s,Vector3D::Length(x1),Segment->GetEnd());
    }
  }

  fclose(fout);
}
