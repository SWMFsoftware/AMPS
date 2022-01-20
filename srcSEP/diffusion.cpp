//SEP diffusion models


#include <cmath>
#include <ctgmath>

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
   double rLarmor=PIC::MolecularData::GetMass(spec)*Relativistic::E2Speed(1.0*GeV2J,PIC::MolecularData::GetMass(spec))/fabs(PIC::MolecularData::GetElectricCharge(spec)*sqrt(absB2));
   double speed=sqrt(vParallel*vParallel+vNorm*vNorm);

   lambda=0.5*absB2/(VacuumPermeability*SummedW)*pow(Lmax*Lmax*rLarmor*Relativistic::Speed2E(speed,PIC::MolecularData::GetMass(spec))*J2GeV,1.0/3.0); 
 
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


//========= Florinskiy =============================
double SEP::Diffusion::Florinskiy::P_s_plus(double k_parallel,double delta_B_s_2) {
  double res,t;

  static double C=2.0*tgamma(gamma_div_two)/(sqrtPi*tgamma(gamma_div_two-0.5)*k_0_s); 

  if (k_parallel<0.0) return 0.0;
  
  t=k_parallel/k_0_s;
  res=C*delta_B_s_2*pow(1.0+t*t,-gamma_div_two); 

  return res;
}

double SEP::Diffusion::Florinskiy::P_s_minus(double k_parallel,double delta_B_s_2) {
  double res,t;

  static double C=2.0*tgamma(gamma_div_two)/(sqrtPi*tgamma(gamma_div_two-0.5)*k_0_s);

  if (k_parallel>0.0) return 0.0;

  t=k_parallel/k_0_s;
  res=C*delta_B_s_2*pow(1.0+t*t,-gamma_div_two);

  return res;
}

  


void SEP::Diffusion::Florinskiy::GetB(double *B,PIC::InterpolationRoutines::CellCentered::cStencil& Stencil) {
  int idim;

  for (idim=0;idim<3;idim++) B[idim]=0.0;

  for (int iStencil=0;iStencil<Stencil.Length;iStencil++) {
    double *ptr=(double*)(Stencil.cell[iStencil]->GetAssociatedDataBufferPointer()+PIC::CPLR::SWMF::MagneticFieldOffset);

    for (idim=0;idim<3;idim++) B[idim]+=Stencil.Weight[iStencil]*ptr[idim];
  }
} 

double SEP::Diffusion::Florinskiy::GetD_mu_mu(double *x,double *v,int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *Node) {
  double r_s,w_plus,w_minus,r_2D_A,sigma_2D_c;
  double D_2D_mu_mu;

  
  //calculate plasma density and magnetic field 
  double B[3]={0.0,0.0,0.0},PlasmaDensity=0.0,omega_minus=0.0,omega_plus=0.0;
  PIC::InterpolationRoutines::CellCentered::cStencil Stencil;
  int idim;

  PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(x,Node,Stencil);
  
  for (int iStencil=0;iStencil<Stencil.Length;iStencil++) {
    double *ptr=(double*)(Stencil.cell[iStencil]->GetAssociatedDataBufferPointer()+PIC::CPLR::SWMF::MagneticFieldOffset);

    for (idim=0;idim<3;idim++) B[idim]+=Stencil.Weight[iStencil]*ptr[idim];

    PlasmaDensity+=(*((double*)(Stencil.cell[iStencil]->GetAssociatedDataBufferPointer()+PIC::CPLR::SWMF::PlasmaNumberDensityOffset)))*Stencil.Weight[iStencil]; 
    omega_minus+=(*((double*)(Stencil.cell[iStencil]->GetAssociatedDataBufferPointer()+PIC::CPLR::SWMF::AlfvenWaveI01Offset)))*Stencil.Weight[iStencil]; 
    omega_plus+=(*(1+(double*)(Stencil.cell[iStencil]->GetAssociatedDataBufferPointer()+PIC::CPLR::SWMF::AlfvenWaveI01Offset)))*Stencil.Weight[iStencil];
  }

  //calculate the Alfven speed
  double B2=Vector3D::DotProduct(B,B);
  double mu=Vector3D::DotProduct(v,B)/sqrt(B2);
  double v_abs=Vector3D::Length(v);

  double vAlfven=sqrt(B2/(VacuumPermeability*PlasmaDensity*_MASS_(_H_))); 
  double Omega=fabs(PIC::MolecularData::GetElectricCharge(spec))*sqrt(B2)/PIC::MolecularData::GetMass(spec); 

  //calculate D_mu_mu
  double delta_B_s_minus_2,delta_B_s_plus_2; 

  delta_B_s_plus_2=( (sigma_c_2D*(1.0-r_s)+r_s)*(omega_minus+omega_plus)+(omega_plus-omega_minus) ) /2.0;
  delta_B_s_minus_2=( (sigma_c_2D*(1.0-r_s)-r_s)*(omega_minus+omega_plus)-(omega_plus-omega_minus) ) /2.0;


  
  auto GetD_mu_mu = [&] (double mu) {
    double res,t0,t1;

    t0=fabs(v_abs*mu-vAlfven);
    t1=fabs(v_abs*mu+vAlfven);

    res=Pi*Omega*Omega*(1.0-mu*mu)/(4.0*B2)*
    (P_s_plus(Omega/t0,delta_B_s_plus_2)/t0+P_s_minus(Omega/t1,delta_B_s_minus_2)/t0);  

    return res;
  };

  double D_mu_mu=GetD_mu_mu(mu);
  
  //calculate d(D_mu_mu)/d(mu)
  double mu_step=2.0/100.0; 
  
  double mu_plus=mu+mu_step; 
  double mu_minus=mu-mu_step;
 
  double dD_mu_mu__d_mu=(GetD_mu_mu(mu_plus)-GetD_mu_mu(mu_minus))/(mu_plus-mu_minus);
  
  return D_mu_mu;
} 





