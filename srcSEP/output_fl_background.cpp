#include "sep.h"
#include "amps2swmf.h"

void SEP::FieldLine::OutputBackgroundData(char* fname, int iFieldLine) {
  namespace FL = PIC::FieldLine;
  FILE *fout=fopen(fname,"w");
  double s=0.0,*x0,*x1;
  FL::cFieldLineSegment *Segment=FL::FieldLinesAll[iFieldLine].GetFirstSegment();

  fprintf(fout,"VARIABLES=\"S[AU]\", \"Heliocentric Distance [AU]\", \"AbsB [T]\", \"AbsB (Parker) [T]\"");

  auto OutputVertex = [&] (double s,double r,FL::cFieldLineVertex* v) {
    fprintf(fout,"%e %e ",s/_AU_,r/_AU_);

    double AbsB,AbsB_Parker; 

    AbsB=Vector3D::Length(v->GetDatum_ptr(FL::DatumAtVertexMagneticField)); 
    AbsB_Parker=SEP::ParkerSpiral::GetAbsB(r);

    fprintf(fout,"%e %e \n",AbsB,AbsB_Parker); 

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
