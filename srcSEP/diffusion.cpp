//SEP diffusion models



#include "sep.h"


//========= Roux2004AJ (LeRoux-2004-AJ) =============================
void SEP::Diffusion::Roux2004AJ::GetPitchAngleDiffusionCoefficient(double& D,double &dD_dmu,double mu,double vParallel,double vNorm,int spec,double FieldLineCoord,PIC::FieldLine::cFieldLineSegment *Segment) {
  static double c=7.0*Pi/8.0*2.0E3/(0.01*_AU_)*(0.2*0.2);

  D=c*(1.0-mu*mu);
  dD_dmu=-2.0*c*mu;
} 

//========= Borovokov_2019_ARXIV (Borovikov-2019-ARXIV, Eq. 6.11)==== 
void SEP::Diffusion::Borovokov_2019_ARXIV::GetPitchAngleDiffusionCoefficient(double& D,double &dD_dmu,double mu,double vParallel,double vNorm,int spec,double FieldLineCoord,PIC::FieldLine::cFieldLineSegment *Segment) {
  namespace FL = PIC::FieldLine;
  namespace MD = PIC::MolecularData;

  FL::cFieldLineVertex* VertexBegin=Segment->GetBegin();
  FL::cFieldLineVertex* VertexEnd=Segment->GetEnd();

  double *B0,*B1;
  double *W0,*W1;
  double *x0,*x1;
  double w0,w1;

  //get the magnetic field and the plasma waves at the corners of the segment
  B0=VertexBegin->GetDatum_ptr(FL::DatumAtVertexMagneticField);
  B1=VertexEnd->GetDatum_ptr(FL::DatumAtVertexMagneticField);

  W0=VertexBegin->GetDatum_ptr(FL::DatumAtVertexPlasmaWaves); 
  W1=VertexEnd->GetDatum_ptr(FL::DatumAtVertexPlasmaWaves);

  x0=VertexBegin->GetX();
  x1=VertexEnd->GetX();
  
  //determine the interpolation coefficients
  w1=fmod(FieldLineCoord,1); 
  w0=1.0-w1;

  //calculate lambda
  double lambda,r2=0.0,absB2=0.0,SummedW=0.0;
  int idim;

  for (idim=0;idim<3;idim++) {
    double t;

    t=w0*B0[idim]+w1*B1[idim];
    absB2+=t*t;

    t=w0*x0[idim]+w1*x1[idim];
    r2+=t*t;
  }


   SummedW=w0*(W0[0]+W0[1])+w1*(W1[0]+W1[1]); 

   double Lmax=0.03*sqrt(r2);
   double rLarmor=MD::GetMass(spec)*Relativistic::E2Speed(1.0*GeV2J,MD::GetMass(spec))/fabs(MD::GetElectricCharge(spec)*sqrt(absB2));
   double speed=sqrt(vParallel*vParallel+vNorm*vNorm);

   lambda=0.5*absB2/(VacuumPermeability*SummedW)*pow(Lmax*Lmax*rLarmor*Relativistic::Speed2E(speed,MD::GetMass(spec))*J2GeV,1.0/3.0); 
 
   if (lambda!=0) {
     double mu2=mu*mu;
     
     mu=fabs(mu);

     D=speed/lambda*(1.0-mu2)*pow(mu,2.0/3.0);
     dD_dmu=speed/lambda*(2.0/3.0/pow(mu,1.0/3.0)-8.0/3.0*pow(mu,5.0/3.0));  
   }
   else {
     D=0.0,dD_dmu=0.0;
   }
} 

//========= Jokopii1966AJ (Jokopii-1966-AJ) =============================
void SEP::Diffusion::Jokopii1966AJ::GetPitchAngleDiffusionCoefficient(double& D,double &dD_dmu,double mu,double vParallel,double vNorm,int spec,double FieldLineCoord,PIC::FieldLine::cFieldLineSegment *Segment) {
  namespace FL = PIC::FieldLine;
  namespace MD = PIC::MolecularData;

  FL::cFieldLineVertex* VertexBegin=Segment->GetBegin();
  FL::cFieldLineVertex* VertexEnd=Segment->GetEnd();

  double *B0,*B1;
  double *W0,*W1;
  double *x0,*x1;
  double w0,w1;

  //reference values of k_max and k_min at 1 AU
  //k_max and k_min are scaled with B, which is turne is scaled with 1/R^2
  const double k_ref_min=1.0E-10;
  const double k_ref_max=1.0E-7;

   

  //get the magnetic field and the plasma waves at the corners of the segment
  B0=VertexBegin->GetDatum_ptr(FL::DatumAtVertexMagneticField);
  B1=VertexEnd->GetDatum_ptr(FL::DatumAtVertexMagneticField);

  W0=VertexBegin->GetDatum_ptr(FL::DatumAtVertexPlasmaWaves);
  W1=VertexEnd->GetDatum_ptr(FL::DatumAtVertexPlasmaWaves);

  x0=VertexBegin->GetX();
  x1=VertexEnd->GetX();

  //determine the interpolation coefficients
  w1=fmod(FieldLineCoord,1);
  w0=1.0-w1;

  double absB2=0.0,r2=0.0;
  int idim;

  for (idim=0;idim<3;idim++) {
    double t;

    t=w0*B0[idim]+w1*B1[idim];
    absB2+=t*t;

    t=w0*x0[idim]+w1*x1[idim]; 
    r2+=t*t;
  }
  
  double k_min,k_max;

  double t=_AU_*_AU_/r2;

  k_min=t*k_ref_min;
  k_max=t*k_ref_max;
  
  double SummW=w0*(W0[0]+W0[1])+w1*(W1[0]+W1[1]);
  double omega,k,P;

  omega=fabs(MD::GetElectricCharge(spec))*sqrt(absB2)/MD::GetMass(spec);

  double kmin=omega/Relativistic::E2Speed(400.0*MeV2J,MD::GetMass(spec));
  double C;

//  C=SummW*VacuumPermeability/(3.0*(pow(k_min,-2.0/3.0)-pow(k_max,-2.0/3.0))/2.0);

  C=0.05*absB2/(3.0*(pow(k_min,-2.0/3.0)-pow(k_max,-2.0/3.0))/2.0);

  k=omega/fabs(vParallel);

  //C/=1.0E7;

  P=C/pow(k,5.0/3.0); 
 
  double c=Pi/4.0*omega*k*P/absB2;

  D=c*(1.0-mu*mu);
  dD_dmu=-c*2*mu;
}





