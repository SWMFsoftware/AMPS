/*
 * Test of the gas-surface interaction models of src/models/surface (Surface.h).
 *
 * Usage: amps                          the checks (report gas-surface-interaction.report)
 *        amps unset-wall-temperature   a diffuse reflection without a wall temperature: must stop with a message
 *        amps unset-parameter          a CLL reflection without coefficients: must stop with a message
 *
 * Every check calls the model through Surface::ParticleInteractionProcessor() (the run-time selection of the model of
 * the species, Surface::Model::Set() or the runtime input block 'physics: gas_surface_interaction') with a particle of
 * the AMPS particle buffer and a face whose frame (ExternalNormal, e0Orthogonal, e1Orthogonal) is skewed with respect to
 * the coordinate axes.  The incident velocity has components along both tangents and the normal.  Sampled moments are
 * compared with their theoretical values within 5 standard errors (the errors are estimated from the samples); the
 * deterministic models (specular, CLL with zero coefficients, the Schamberg speed and cone) are checked to 1e-12.
 * With a = sqrt(2 k Tw/m), v the incident velocity, v_n = v.n < 0, v_t the tangential components, u the re-emitted one:
 *   diffuse      <u.n> = sqrt(pi)/2 a, <(u.n)^2> = a^2, <u_t> = 0, <u_t^2> = a^2/2 each, <u^2> = 2a^2 (all species)
 *   specular     u = v - 2 v_n n exactly
 *   maxwell      <u.n> = alpha sqrt(pi)/2 a + (1-alpha)|v_n|, <u_t> = (1-alpha) v_t, <u_t^2> = alpha a^2/2 + (1-alpha) v_t^2
 *   cll          alpha_n = sigma_t = 1 is diffuse (the moments above), alpha_n = sigma_t = 0 is specular (exactly);
 *                otherwise (Rice normal part with nu = sqrt(1-alpha_n)|v_n|, sigma^2 = alpha_n a^2/2):
 *                <u.n> = sigma sqrt(pi/2) L_{1/2}(-nu^2/(2 sigma^2)), <(u.n)^2> = nu^2 + 2 sigma^2,
 *                <(u.n)^4> = nu^4 + 8 sigma^2 nu^2 + 8 sigma^4, <u_t> = (1-sigma_t) v_t,
 *                Var(u_t) = sigma_t (2-sigma_t) a^2/2 each (also sigma_t > 1: back-scattering), per species;
 *                reciprocity: with the incident molecules of a gas in equilibrium with the wall (the flux Maxwellian at
 *                Tw), the re-emitted ones have the diffuse moments, and the joint histograms of (|v_n|, u.n) and of
 *                (v.e0, u.e0) are symmetric under (x, y) -> (y, x) and (x, y) -> (-y, -x) (chi-square test)
 *   quasi-specular CLL   the CLL moments with small coefficients
 *   nocilla      s = U_r cos(theta_r)/c_r: <u.n> = c_r <x>, <(u.n)^2> = c_r^2 <x^2> with
 *                <x> = (s e^{-s^2} + (1/2+s^2) sqrt(pi)(1+erf s))/(e^{-s^2} + s sqrt(pi)(1+erf s)),
 *                <x^2> = (e^{-s^2}(1+s^2) + sqrt(pi)(1+erf s)(3s/2+s^3))/(e^{-s^2} + s sqrt(pi)(1+erf s)),
 *                <u.t> = U_r sin(theta_r) along the incident tangential direction t, <u.(n x t)> = 0, variances c_r^2/2;
 *                the same for another incident speed; U_r = 0, T_r = Tw is diffuse
 *   nocilla with accommodation coefficients (Hurlbut & Sherman; Cook & Hoffbauer 1997, NASA TP-3701, Eqs. 1, 2, 10,
 *                15-17): the solved parameters reproduce A = <u.n>, Q = <u.t> and V2 = <u^2> exactly; sampled
 *                <u.n> = A, <u.t> = Q, <u^2> = V2, the tangential and binormal variances c_r^2/2; a drift toward the wall
 *                (Lambda < 0, also the sampler with s < 0); sigma' = sigma = epsilon = 1 is diffuse; inconsistent
 *                coefficients re-emit the mean velocity and are counted
 *   schamberg    |u| = V_r = sqrt((1-alpha)|v|^2 + alpha 4kTw/m) exactly; every direction within phi0 of the axis;
 *                the axis at sin(theta_r) = sin(theta_i)^nu in the plane of incidence; <u/|u|> = Phi(phi0) axis, Phi by an
 *                independent quadrature of the beam ~ cos(pi phi/(2 phi0)) (and Schamberg's closed form, Eq. 16 of
 *                Mostaza Prieto et al. 2014); nu = 1, phi0 = 0, alpha = 0 is specular (exactly); nu = inf, phi0 = 90,
 *                alpha = 1 is the cosine law (<cos> = 2/3, <cos^2> = 1/2) at sqrt(4kTw/m); a beam that reaches below
 *                the surface re-emits only above it
 *   adsorption   the particle is removed
 *   input        the block 'physics: gas_surface_interaction' installs a valid model and rejects invalid ones with a
 *                message (and installs nothing then)
 */

#include "pic.h"
#include "Surface.h"
#include "src/general/input/amps_input.h"

#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <functional>
#include <string>
#include <vector>

namespace {

FILE *Report=NULL;
int nFailed=0,nPassed=0;

void Check(bool pass,const char *name,const char *format,...) {
  char message[2000];
  va_list args;

  va_start(args,format);
  vsnprintf(message,sizeof(message),format,args);
  va_end(args);

  //one line per check (a diagnostic of the input reader quotes the offending line on a second line)
  for (char *c=message;*c!='\0';c++) if (*c=='\n') *c='|';

  fprintf(Report,"%s %s %s\n",pass ? "PASS" : "FAIL",name,message);
  printf("%s %s %s\n",pass ? "PASS" : "FAIL",name,message);
  fflush(Report);

  if (pass) nPassed++;
  else nFailed++;
}

void Info(const char *format,...) {
  char message[2000];
  va_list args;

  va_start(args,format);
  vsnprintf(message,sizeof(message),format,args);
  va_end(args);

  fprintf(Report,"INFO %s\n",message);
  printf("INFO %s\n",message);
}

const double Tw=300.0;
const long int N=1000000;    //samples per statistical check
const double Z=5.0;          //the limit of the statistical checks in standard errors

CutCell::cTriangleFace Face;
double *n=Face.ExternalNormal,*e0=Face.e0Orthogonal,*e1=Face.e1Orthogonal;
long int ParticlePtr[PIC::nTotalSpecies];

double Dot(const double *a,const double *b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
void Cross(double *c,const double *a,const double *b) {c[0]=a[1]*b[2]-a[2]*b[1],c[1]=a[2]*b[0]-a[0]*b[2],c[2]=a[0]*b[1]-a[1]*b[0];}

double MostProbableSpeed(int spec,double T) {return sqrt(2.0*Kbol*T/PIC::MolecularData::GetMass(spec));}

//the velocity with the components (c0,c1,cn) in the face frame
void FromFrame(double *v,double c0,double c1,double cn) {for (int i=0;i<3;i++) v[i]=c0*e0[i]+c1*e1[i]+cn*n[i];}

//a statistic: mean and standard error of the mean
struct cStat {
  double s=0.0,s2=0.0;
  long int n=0;

  void Add(double x) {s+=x,s2+=x*x,n++;}
  double Mean() const {return s/n;}
  double Error() const {double m=Mean(); return sqrt(std::max(s2/n-m*m,0.0)/n);}
};

//compare a sampled mean with theory: |mean-theory| <= Z standard errors (and an absolute floor for exact zeros)
void CheckMean(const char *name,const char *what,const cStat &st,double theory,double scale) {
  double err=st.Error(),dev=st.Mean()-theory;
  double z=(err>0.0) ? dev/err : ((fabs(dev)<=1.0E-12*scale) ? 0.0 : 1.0E300);

  Check(fabs(z)<=Z,name,"%s: sampled %.6e, theory %.6e, deviation %+.2f standard errors",what,st.Mean()/scale,theory/scale,z);
}

//one interaction through the run-time dispatcher
int Interact(int spec,double *v) {
  double x[3]={0.0,0.0,0.0};

  return Surface::ParticleInteractionProcessor(ParticlePtr[spec],x,v,&Face,NULL);
}

//the moments of the re-emitted velocity in the face frame
struct cMoments {
  cStat un,un2,un4,ut0,ut02,ut1,ut12,u2,dt,db,dt2,db2;
  long int nonpositive=0;
};

cMoments Sample(int spec,const double *vin,long int nSamples,const double *t=NULL) {
  cMoments m;
  double v[3],b[3];

  if (t!=NULL) Cross(b,n,t);

  for (long int i=0;i<nSamples;i++) {
    memcpy(v,vin,3*sizeof(double));
    Interact(spec,v);

    double un=Dot(v,n),u0=Dot(v,e0),u1=Dot(v,e1);

    m.un.Add(un),m.un2.Add(un*un),m.un4.Add(un*un*un*un);
    m.ut0.Add(u0),m.ut02.Add(u0*u0),m.ut1.Add(u1),m.ut12.Add(u1*u1);
    m.u2.Add(Dot(v,v));
    if (t!=NULL) {double a=Dot(v,t),c=Dot(v,b); m.dt.Add(a),m.db.Add(c),m.dt2.Add(a*a),m.db2.Add(c*c);}
    if (un<=0.0) m.nonpositive++;
  }

  return m;
}

//the diffuse moments at the temperature T
void CheckDiffuseMoments(const char *name,const cMoments &m,double a) {
  CheckMean(name,"<u.n>/a",m.un,sqrt(Pi)/2.0*a,a);
  CheckMean(name,"<(u.n)^2>/a^2",m.un2,a*a,a*a);
  CheckMean(name,"<u.e0>/a",m.ut0,0.0,a);
  CheckMean(name,"<u.e1>/a",m.ut1,0.0,a);
  CheckMean(name,"<(u.e0)^2>/a^2",m.ut02,0.5*a*a,a*a);
  CheckMean(name,"<(u.e1)^2>/a^2",m.ut12,0.5*a*a,a*a);
  CheckMean(name,"<u^2>/a^2",m.u2,2.0*a*a,a*a);
  Check(m.nonpositive==0,name,"every re-emitted particle leaves the wall (u.n > 0): %ld do not",m.nonpositive);
}

//the CLL moments (v: incident, in the face frame)
void CheckCLLMoments(const char *name,const cMoments &m,double a,double alpha_n,double sigma_t,double v0,double v1,double vn) {
  const double nu=sqrt(1.0-alpha_n)*fabs(vn),sig2=0.5*alpha_n*a*a,alpha_t=sigma_t*(2.0-sigma_t);

  //the mean of the Rice distribution: sigma sqrt(pi/2) L_{1/2}(x), x = -nu^2/(2 sigma^2),
  //L_{1/2}(x) = e^{x/2} ((1-x) I0(-x/2) - x I1(-x/2))
  if (sig2>0.0) {
    const double x=-nu*nu/(2.0*sig2),y=-0.5*x;
    const double L=exp(-y)*((1.0-x)*std::cyl_bessel_i(0.0,y)-x*std::cyl_bessel_i(1.0,y));
    CheckMean(name,"<u.n>/a (Rice mean)",m.un,sqrt(sig2)*sqrt(0.5*Pi)*L,a);
  }

  CheckMean(name,"<(u.n)^2>/a^2 = ((1-alpha_n) v_n^2 + alpha_n a^2)/a^2",m.un2,nu*nu+2.0*sig2,a*a);
  CheckMean(name,"<(u.n)^4>/a^4",m.un4,pow(nu,4)+8.0*sig2*nu*nu+8.0*sig2*sig2,pow(a,4));
  CheckMean(name,"<u.e0>/a = (1-sigma_t) v.e0/a",m.ut0,(1.0-sigma_t)*v0,a);
  CheckMean(name,"<u.e1>/a = (1-sigma_t) v.e1/a",m.ut1,(1.0-sigma_t)*v1,a);
  CheckMean(name,"<(u.e0)^2>/a^2",m.ut02,pow((1.0-sigma_t)*v0,2)+0.5*alpha_t*a*a,a*a);
  CheckMean(name,"<(u.e1)^2>/a^2",m.ut12,pow((1.0-sigma_t)*v1,2)+0.5*alpha_t*a*a,a*a);
  Check(m.nonpositive==0,name,"every re-emitted particle leaves the wall (u.n > 0): %ld do not",m.nonpositive);
}

//the moments of the normal speed ratio of the Nocilla re-emission
double NocillaMeanX(double s) {
  const double E=erfc(-s);   //1+erf(s), accurate also for s < 0
  return (s*exp(-s*s)+(0.5+s*s)*sqrt(Pi)*E)/(exp(-s*s)+s*sqrt(Pi)*E);
}

double NocillaMeanX2(double s) {
  const double E=erfc(-s);   //1+erf(s), accurate also for s < 0
  return (exp(-s*s)*(1.0+s*s)+sqrt(Pi)*E*(1.5*s+s*s*s))/(exp(-s*s)+s*sqrt(Pi)*E);
}

//the mean cosine (and mean cos^2) of the Schamberg beam by an independent composite Simpson quadrature
void BeamQuadrature(double phi0,double *mcos,double *mcos2) {
  const int m=20000;
  double h=phi0/m,s0=0.0,s1=0.0,s2=0.0;

  for (int i=0;i<=m;i++) {
    double phi=i*h,w=((i==0)||(i==m)) ? 1.0 : ((i%2==1) ? 4.0 : 2.0);
    double f=cos(0.5*Pi*phi/phi0)*sin(phi)*w;

    s0+=f,s1+=f*cos(phi),s2+=f*cos(phi)*cos(phi);
  }

  *mcos=s1/s0,*mcos2=s2/s0;
}

//select a model with its parameters through the runtime input block (the path of the applications)
bool ApplyBlock(const std::string &text,std::string *error) {
  AMPS::Input::Document doc;
  AMPS::Input::Status st=AMPS::Input::ParseText(text,"test-block.in",".",&doc);

  if (!st.ok()) {*error=st.message; return false;}
  return Surface::RuntimeInput::Apply(doc.root.blocks.front(),error);
}

std::string SpeciesBlock(const char *species,const std::string &keys) {
  return std::string("    #begin species: ")+species+"\n"+keys+"    #end species\n";
}

}  //namespace

//----------------------------------------------------------------------------------------------------------------
int main(int argc,char **argv) {
  const std::string mode=(argc>1) ? argv[1] : "checks";

  PIC::Init_BeforeParser();
  PIC::ParticleBuffer::Init(100);

  //the face frame (as CutCell::cTriangleFace::SetFaceNodes builds it from two skewed edges)
  {
    double a[3]={1.0,0.2,-0.3},b[3]={0.1,1.0,0.4},l;

    Cross(n,a,b);
    l=sqrt(Dot(n,n));
    for (int i=0;i<3;i++) n[i]/=l;
    l=sqrt(Dot(a,a));
    for (int i=0;i<3;i++) e0[i]=a[i]/l;
    Cross(e1,n,e0);
  }

  for (int s=0;s<PIC::nTotalSpecies;s++) {
    ParticlePtr[s]=PIC::ParticleBuffer::GetNewParticle();
    PIC::ParticleBuffer::SetI(s,ParticlePtr[s]);
  }

  //the modes that must stop the run
  if (mode=="unset-wall-temperature") {
    double v[3];

    FromFrame(v,700.0,-400.0,-900.0);
    Surface::Model::Set(-1,Surface::Model::Diffuse);
    Interact(0,v);
    printf("ERROR: the diffuse reflection ran without a wall temperature\n");
    MPI_Finalize();
    return 0;
  }

  if (mode=="unset-parameter") {
    double v[3];

    FromFrame(v,700.0,-400.0,-900.0);
    Surface::Temperature::Isothremal::Temp=Tw;
    Surface::Model::Set(-1,Surface::Model::CLL);
    Interact(0,v);
    printf("ERROR: the CLL reflection ran without coefficients\n");
    MPI_Finalize();
    return 0;
  }

  if (mode!="checks") {
    printf("usage: %s [unset-wall-temperature|unset-parameter]\n",argv[0]);
    MPI_Finalize();
    return 2;
  }

  Report=fopen("gas-surface-interaction.report","w");
  Info("species: %d; samples per check %ld; limit %.1f standard errors; frame n=(%.6f,%.6f,%.6f) e0=(%.6f,%.6f,%.6f) e1=(%.6f,%.6f,%.6f)",
      PIC::nTotalSpecies,N,Z,n[0],n[1],n[2],e0[0],e0[1],e0[2],e1[0],e1[1],e1[2]);
  Check((fabs(Dot(n,n)-1.0)<1e-14)&&(fabs(Dot(e0,e0)-1.0)<1e-14)&&(fabs(Dot(e1,e1)-1.0)<1e-14)&&(fabs(Dot(n,e0))<1e-14)&&(fabs(Dot(n,e1))<1e-14)&&(fabs(Dot(e0,e1))<1e-14),
      "frame","the face frame is orthonormal");

  //the compile-time choice is the initial model of every species
  {
    bool ok=true;
    for (int s=0;s<PIC::nTotalSpecies;s++) ok&=(Surface::Model::Selected[s]==_SURFACE__GAS_SURFACE_INTERACTION_MODEL_);
    Check(ok,"compile-time","the initial model of every species is _SURFACE__GAS_SURFACE_INTERACTION_MODEL_ (%s)",Surface::Model::Name(_SURFACE__GAS_SURFACE_INTERACTION_MODEL_));
  }

  Surface::Temperature::Isothremal::Temp=Tw;

  //the incident velocity of each species: components along both tangents and the normal, in units of its a
  double vin[PIC::nTotalSpecies][3],vc[PIC::nTotalSpecies][3],aw[PIC::nTotalSpecies];

  for (int s=0;s<PIC::nTotalSpecies;s++) {
    aw[s]=MostProbableSpeed(s,Tw);
    vc[s][0]=1.33*aw[s],vc[s][1]=-0.76*aw[s],vc[s][2]=-1.71*aw[s];
    FromFrame(vin[s],vc[s][0],vc[s][1],vc[s][2]);
  }

  //------------------------------------------------------------------------- specular
  {
    Surface::Model::Set(-1,Surface::Model::Specular);
    double maxerr=0.0;

    for (int s=0;s<PIC::nTotalSpecies;s++) {
      double v[3],expect[3];

      memcpy(v,vin[s],sizeof(v));
      int code=Interact(s,v);
      FromFrame(expect,vc[s][0],vc[s][1],-vc[s][2]);
      for (int i=0;i<3;i++) maxerr=std::max(maxerr,fabs(v[i]-expect[i])/aw[s]);
      Check(code==_PARTICLE_REJECTED_ON_THE_FACE_,"specular","species %s: the particle is re-emitted",PIC::MolecularData::GetChemSymbol(s));
    }

    Check(maxerr<1.0E-12,"specular","u = v - 2 (v.n) n: largest error %.2e a",maxerr);
  }

  //------------------------------------------------------------------------- adsorption
  {
    Surface::Model::Set(-1,Surface::Model::Adsorption);
    double v[3];

    memcpy(v,vin[0],sizeof(v));
    Check(Interact(0,v)==_PARTICLE_DELETED_ON_THE_FACE_,"adsorption","the particle is removed (_PARTICLE_DELETED_ON_THE_FACE_)");
  }

  //------------------------------------------------------------------------- diffuse (every species: its own a)
  Surface::Model::Set(-1,Surface::Model::Diffuse);

  for (int s=0;s<PIC::nTotalSpecies;s++) {
    char name[100];

    snprintf(name,sizeof(name),"diffuse[%s]",PIC::MolecularData::GetChemSymbol(s));
    CheckDiffuseMoments(name,Sample(s,vin[s],N),aw[s]);
  }

  //------------------------------------------------------------------------- maxwell (per-species coefficients)
  {
    const double alpha[3]={0.75,0.3,1.0};
    std::string error;

    Surface::Model::Set(-1,Surface::Model::Maxwell);
    for (int s=0;s<PIC::nTotalSpecies;s++) Surface::MaxwellReflection::SetAccommodationCoefficient(s,alpha[s%3],&error);

    for (int s=0;s<PIC::nTotalSpecies;s++) {
      char name[100];
      const double al=alpha[s%3],a=aw[s];
      cMoments m=Sample(s,vin[s],N);

      snprintf(name,sizeof(name),"maxwell[%s,alpha=%g]",PIC::MolecularData::GetChemSymbol(s),al);
      CheckMean(name,"<u.n>/a",m.un,al*sqrt(Pi)/2.0*a+(1.0-al)*fabs(vc[s][2]),a);
      CheckMean(name,"<(u.n)^2>/a^2",m.un2,al*a*a+(1.0-al)*vc[s][2]*vc[s][2],a*a);
      CheckMean(name,"<u.e0>/a",m.ut0,(1.0-al)*vc[s][0],a);
      CheckMean(name,"<u.e1>/a",m.ut1,(1.0-al)*vc[s][1],a);
      CheckMean(name,"<(u.e0)^2>/a^2",m.ut02,0.5*al*a*a+(1.0-al)*vc[s][0]*vc[s][0],a*a);
    }
  }

  //------------------------------------------------------------------------- CLL
  {
    std::string error;

    Surface::Model::Set(-1,Surface::Model::CLL);

    //alpha_n = sigma_t = 1: diffuse
    Surface::CLL::SetCoefficients(-1,1.0,1.0,&error);
    for (int s=0;s<PIC::nTotalSpecies;s++) {
      char name[100];

      snprintf(name,sizeof(name),"cll-limit-diffuse[%s]",PIC::MolecularData::GetChemSymbol(s));
      CheckDiffuseMoments(name,Sample(s,vin[s],N),aw[s]);
    }

    //alpha_n = sigma_t = 0: specular, exactly
    Surface::CLL::SetCoefficients(-1,0.0,0.0,&error);
    double maxerr=0.0;

    for (int s=0;s<PIC::nTotalSpecies;s++) for (int k=0;k<1000;k++) {
      double v[3],expect[3];

      memcpy(v,vin[s],sizeof(v));
      Interact(s,v);
      FromFrame(expect,vc[s][0],vc[s][1],-vc[s][2]);
      for (int i=0;i<3;i++) maxerr=std::max(maxerr,fabs(v[i]-expect[i])/aw[s]);
    }

    Check(maxerr<1.0E-12,"cll-limit-specular","alpha_n = sigma_t = 0 gives u = v - 2 (v.n) n: largest error %.2e a (1000 samples per species)",maxerr);

    //intermediate coefficients, different for each species (the second one back-scatters: sigma_t > 1)
    const double an[3]={0.3,0.6,0.1},st[3]={0.8,1.5,0.35};

    for (int s=0;s<PIC::nTotalSpecies;s++) Surface::CLL::SetCoefficients(s,an[s%3],st[s%3],&error);

    for (int s=0;s<PIC::nTotalSpecies;s++) {
      char name[100];

      snprintf(name,sizeof(name),"cll[%s,alpha_n=%g,sigma_t=%g]",PIC::MolecularData::GetChemSymbol(s),an[s%3],st[s%3]);
      CheckCLLMoments(name,Sample(s,vin[s],N),aw[s],an[s%3],st[s%3],vc[s][0],vc[s][1],vc[s][2]);
    }

    //reciprocity: the incident molecules of a gas in equilibrium with the wall (the flux Maxwellian at Tw toward the wall)
    {
      const int s=0,nb=12;
      const double a=aw[s],wmax=3.0*a,tmax=3.0*a;
      std::vector<double> Hn(nb*nb,0.0),Ht(nb*nb,0.0);
      cMoments m;
      double v[3];
      CutCell::cTriangleFace *f=&Face;

      Surface::CLL::SetCoefficients(s,0.3,0.8,&error);

      for (long int i=0;i<4*N;i++) {
        //the incident velocity: the negative of a diffusely emitted one (Surface::DiffuseReflection::Reemit at Tw)
        Surface::DiffuseReflection::Reemit(v,Tw,s,f);
        for (int k=0;k<3;k++) v[k]=-v[k];

        double win=-Dot(v,n),tin=Dot(v,e0);
        Interact(s,v);
        double wout=Dot(v,n),tout=Dot(v,e0);

        m.un.Add(wout),m.un2.Add(wout*wout),m.ut0.Add(tout),m.ut02.Add(tout*tout),m.u2.Add(Dot(v,v));
        m.ut1.Add(Dot(v,e1)),m.ut12.Add(Dot(v,e1)*Dot(v,e1));
        if (wout<=0.0) m.nonpositive++;

        int i0=(int)(win/wmax*nb),i1=(int)(wout/wmax*nb);
        if ((i0<nb)&&(i1<nb)) Hn[i0+nb*i1]+=1.0;

        int j0=(int)floor((tin+tmax)/(2.0*tmax)*nb),j1=(int)floor((tout+tmax)/(2.0*tmax)*nb);
        if ((j0>=0)&&(j0<nb)&&(j1>=0)&&(j1<nb)) Ht[j0+nb*j1]+=1.0;
      }

      CheckDiffuseMoments("cll-reciprocity-equilibrium",m,a);

      //detailed balance of the normal part: H(win in I, wout in J) = H(win in J, wout in I)
      double chi2=0.0;
      int dof=0;

      for (int i=0;i<nb;i++) for (int j=i+1;j<nb;j++) {
        double x=Hn[i+nb*j],y=Hn[j+nb*i];
        if (x+y>20.0) chi2+=(x-y)*(x-y)/(x+y),dof++;
      }

      Check(chi2<dof+Z*sqrt(2.0*dof),"cll-reciprocity-normal","the joint histogram of (|v.n|, u.n) is symmetric: chi2 = %.1f for %d degrees of freedom",chi2,dof);

      //detailed balance of the tangential part: H(tin in I, tout in J) = H(tin in -J, tout in -I)
      chi2=0.0,dof=0;

      for (int i=0;i<nb;i++) for (int j=0;j<nb;j++) {
        int ii=nb-1-j,jj=nb-1-i;   //the bins -J and -I

        if (i+nb*j<ii+nb*jj) {
          double x=Ht[i+nb*j],y=Ht[ii+nb*jj];
          if (x+y>20.0) chi2+=(x-y)*(x-y)/(x+y),dof++;
        }
      }

      Check(chi2<dof+Z*sqrt(2.0*dof),"cll-reciprocity-tangential","the joint histogram of (v.e0, u.e0) is symmetric under (x,y)->(-y,-x): chi2 = %.1f for %d degrees of freedom",chi2,dof);
    }

    //the coefficients are checked
    {error.clear(); bool ok=Surface::CLL::SetCoefficients(0,1.2,0.5,&error); Check((ok==false)&&(error.empty()==false),"cll-parameters","alpha_n = 1.2 is rejected: %s",error.c_str());}
    {error.clear(); bool ok=Surface::CLL::SetCoefficients(0,0.5,2.5,&error); Check((ok==false)&&(error.empty()==false),"cll-parameters","sigma_t = 2.5 is rejected: %s",error.c_str());}
    Check((Surface::CLL::NormalPartKineticEnergyAccommodationCoefficient[0]==0.3)&&(Surface::CLL::TangentialMomentumAccommodationCoefficient[0]==0.8),
        "cll-parameters","a rejected value does not change the coefficients");
  }

  //------------------------------------------------------------------------- quasi-specular CLL
  {
    std::string error;
    const double an=0.05,st=0.15;

    Surface::Model::Set(-1,Surface::Model::QuasiSpecularCLL);
    Surface::QuasiSpecularReflection::CLL::SetCoefficients(-1,an,st,&error);

    for (int s=0;s<PIC::nTotalSpecies;s++) {
      char name[100];
      cMoments m=Sample(s,vin[s],N);

      snprintf(name,sizeof(name),"quasi-specular-cll[%s,alpha_n=%g,sigma_t=%g]",PIC::MolecularData::GetChemSymbol(s),an,st);
      CheckCLLMoments(name,m,aw[s],an,st,vc[s][0],vc[s][1],vc[s][2]);

      //a narrow lobe about the specular direction: the rms deviation from the specular velocity is small
      double dev2=m.un2.Mean()-2.0*fabs(vc[s][2])*m.un.Mean()+vc[s][2]*vc[s][2]+
          m.ut02.Mean()-2.0*vc[s][0]*m.ut0.Mean()+vc[s][0]*vc[s][0]+m.ut12.Mean()-2.0*vc[s][1]*m.ut1.Mean()+vc[s][1]*vc[s][1];
      double vmag=sqrt(vc[s][0]*vc[s][0]+vc[s][1]*vc[s][1]+vc[s][2]*vc[s][2]);
      Check(sqrt(dev2)<0.35*vmag,name,"quasi-specular: the rms deviation from the specular velocity is %.3f |v| (< 0.35)",sqrt(dev2)/vmag);
    }

    {error.clear(); bool ok=Surface::QuasiSpecularReflection::CLL::SetCoefficients(0,0.1,1.5,&error); Check((ok==false)&&(error.empty()==false),"quasi-specular-cll-parameters","sigma_t = 1.5 (back-scattering) is rejected: %s",error.c_str());}
  }

  //------------------------------------------------------------------------- Nocilla
  {
    std::string error;

    //the sampler of the normal speed ratio x ~ x exp(-(x-s)^2)
    const double svals[4]={0.0,0.5,2.47,6.0};

    for (double s : svals) {
      cStat x1,x2;
      char name[100];

      for (long int i=0;i<N;i++) {
        double x=Surface::QuasiSpecularReflection::Nocilla::SampleFluxNormalSpeedRatio(s);
        x1.Add(x),x2.Add(x*x);
      }

      snprintf(name,sizeof(name),"nocilla-normal-sampler[s=%g]",s);
      CheckMean(name,"<x>",x1,NocillaMeanX(s),1.0);
      CheckMean(name,"<x^2>",x2,NocillaMeanX2(s),1.0);
    }

    //the re-emission: U_r, theta_r, T_r per species
    const double Ur[3]={1500.0,400.0,0.0},thr[3]={30.0,70.0,0.0},Tr[3]={300.0,150.0,Tw};

    Surface::Model::Set(-1,Surface::Model::QuasiSpecularNocilla);
    for (int s=0;s<PIC::nTotalSpecies;s++) Surface::QuasiSpecularReflection::Nocilla::SetParameters(s,Ur[s%3],thr[s%3]*Pi/180.0,Tr[s%3],&error);

    for (int s=0;s<PIC::nTotalSpecies;s++) {
      const double cr=MostProbableSpeed(s,Tr[s%3]),th=thr[s%3]*Pi/180.0,sr=Ur[s%3]*cos(th)/cr;
      double t[3],l=sqrt(vc[s][0]*vc[s][0]+vc[s][1]*vc[s][1]);

      for (int i=0;i<3;i++) t[i]=(vc[s][0]*e0[i]+vc[s][1]*e1[i])/l;

      for (int speed=0;speed<2;speed++) {
        char name[150];
        double v[3];

        //two incident speeds in the same direction: the re-emission does not depend on the incident speed
        for (int i=0;i<3;i++) v[i]=vin[s][i]*((speed==0) ? 1.0 : 3.0);
        cMoments m=Sample(s,v,N,t);

        snprintf(name,sizeof(name),"nocilla[%s,U_r=%g,theta_r=%g,T_r=%g,|v|x%d]",PIC::MolecularData::GetChemSymbol(s),Ur[s%3],thr[s%3],Tr[s%3],(speed==0) ? 1 : 3);
        CheckMean(name,"<u.n>/c_r",m.un,cr*NocillaMeanX(sr),cr);
        CheckMean(name,"<(u.n)^2>/c_r^2",m.un2,cr*cr*NocillaMeanX2(sr),cr*cr);
        CheckMean(name,"<u.t>/c_r (forward drift along the incident tangential direction t)",m.dt,Ur[s%3]*sin(th),cr);
        CheckMean(name,"<u.(n x t)>/c_r",m.db,0.0,cr);
        CheckMean(name,"<(u.t)^2>/c_r^2",m.dt2,pow(Ur[s%3]*sin(th),2)+0.5*cr*cr,cr*cr);
        CheckMean(name,"<(u.(n x t))^2>/c_r^2",m.db2,0.5*cr*cr,cr*cr);
        Check(m.nonpositive==0,name,"every re-emitted particle leaves the wall (u.n > 0): %ld do not",m.nonpositive);
      }

      if (Ur[s%3]==0.0) {
        char name[100];

        snprintf(name,sizeof(name),"nocilla-limit-diffuse[%s]",PIC::MolecularData::GetChemSymbol(s));
        CheckDiffuseMoments(name,Sample(s,vin[s],N),aw[s]);
      }
    }

    {error.clear(); bool ok=Surface::QuasiSpecularReflection::Nocilla::SetParameters(0,1000.0,95.0*Pi/180.0,300.0,&error); Check((ok==false)&&(error.empty()==false),"nocilla-parameters","theta_r = 95 deg is rejected: %s",error.c_str());}
    {error.clear(); bool ok=Surface::QuasiSpecularReflection::Nocilla::SetParameters(0,1000.0,0.5,0.0,&error); Check((ok==false)&&(error.empty()==false),"nocilla-parameters","T_r = 0 is rejected: %s",error.c_str());}
  }

  //------------------------------------------------------------------------- Nocilla with accommodation coefficients
  //(Hurlbut & Sherman; the equations of Cook & Hoffbauer 1997, NASA TP-3701: <u>.t = (1-sigma) v_t (Eq. 1),
  //<u>.n = (1-sigma')|v_n| + sigma' sqrt(pi k Tw/2m) (Eq. 2), <u^2> = (1-epsilon) v^2 + epsilon 4kTw/m (Eq. 10))
  {
    namespace NC=Surface::QuasiSpecularReflection::Nocilla;
    std::string error;

    //the sampler of x ~ x exp(-(x-s)^2) for a drift toward the wall (s < 0)
    for (double s : {-0.5,-3.0}) {
      cStat x1,x2;
      char name[100];

      for (long int i=0;i<N;i++) {
        double x=NC::SampleFluxNormalSpeedRatio(s);
        x1.Add(x),x2.Add(x*x);
      }

      snprintf(name,sizeof(name),"nocilla-normal-sampler[s=%g]",s);
      CheckMean(name,"<x>",x1,NocillaMeanX(s),1.0);
      CheckMean(name,"<x^2>",x2,NocillaMeanX2(s),1.0);
    }

    //the moments of the library (scaled erfc for Lambda < 0) against the formulas of this test
    {
      double maxerr=0.0;

      for (double L : {-5.0,-2.0,-0.3,0.0,0.7,2.47,8.0}) {
        double m1,m2;

        NC::FluxNormalMoments(L,&m1,&m2);
        maxerr=std::max(maxerr,std::max(fabs(m1/NocillaMeanX(L)-1.0),fabs(m2/NocillaMeanX2(L)-1.0)));
      }

      Check(maxerr<1.0E-9,"nocilla-hs-moments","<x>, <x^2> of x exp(-(x-Lambda)^2) for Lambda = -5..8: largest relative difference %.2e",maxerr);
    }

    //per species: (sigma', sigma, epsilon); the second one has Lambda < 0 (a drift toward the wall), the third is near-specular
    const double sn[3]={0.6,1.0,0.2},st[3]={0.7,1.0,0.3},ep[3]={0.5,0.5,0.1};

    Surface::Model::Set(-1,Surface::Model::QuasiSpecularNocilla);
    for (int s=0;s<PIC::nTotalSpecies;s++) NC::SetAccommodationCoefficients(s,sn[s%3],st[s%3],ep[s%3],&error);

    for (int s=0;s<PIC::nTotalSpecies;s++) {
      const double mass=PIC::MolecularData::GetMass(s),a=aw[s];
      const double vn=fabs(vc[s][2]),vt=sqrt(vc[s][0]*vc[s][0]+vc[s][1]*vc[s][1]),v2=vn*vn+vt*vt;
      const double A=(1.0-sn[s%3])*vn+sn[s%3]*sqrt(Pi*Kbol*Tw/(2.0*mass)),Q=(1.0-st[s%3])*vt,V2=(1.0-ep[s%3])*v2+ep[s%3]*4.0*Kbol*Tw/mass;
      double Lambda,c_r,A1,Q1,t[3];
      char name[200];

      int code=NC::SolveAccommodation(vn,vt,sn[s%3],st[s%3],ep[s%3],Tw,mass,&Lambda,&c_r,&A1,&Q1);
      snprintf(name,sizeof(name),"nocilla-hs[%s,sigma'=%g,sigma=%g,epsilon=%g]",PIC::MolecularData::GetChemSymbol(s),sn[s%3],st[s%3],ep[s%3]);

      //the solution reproduces the three moments (Cook & Hoffbauer Eqs. 15-17), with the moment formulas of this test
      double m1=NocillaMeanX(Lambda),m2=NocillaMeanX2(Lambda);
      Check((code==0)&&(fabs(A1-A)<1e-12*a)&&(fabs(Q1-Q)<1e-12*a)&&(fabs(c_r*m1-A)<1.0E-9*a)&&(fabs(c_r*c_r*(m2+1.0)+Q*Q-V2)<1.0E-9*a*a),name,
          "the solution Lambda=%.6f, T_r=%.3f K gives <u.n> = c_r <x> = A and <u^2> = c_r^2 (<x^2>+1) + Q^2 = V2 (Eqs. 15-17)",
          Lambda,mass*c_r*c_r/(2.0*Kbol));
      if (s%3==1) Check(Lambda<0.0,name,"Lambda = %.4f < 0: the drift points toward the wall (theta_u > 90 deg)",Lambda);

      for (int i=0;i<3;i++) t[i]=(vc[s][0]*e0[i]+vc[s][1]*e1[i])/vt;
      cMoments m=Sample(s,vin[s],N,t);

      CheckMean(name,"<u.n>/a = A/a (Eq. 2)",m.un,A,a);
      CheckMean(name,"<u.t>/a = Q/a (Eq. 1)",m.dt,Q,a);
      CheckMean(name,"<u.(n x t)>/a",m.db,0.0,a);
      CheckMean(name,"<u^2>/a^2 = V2/a^2 (Eq. 10)",m.u2,V2,a*a);
      CheckMean(name,"<(u.(n x t))^2>/a^2 = c_r^2/2 (the temperature T_r)",m.db2,0.5*c_r*c_r,a*a);
      CheckMean(name,"<(u.t)^2>/a^2 = Q^2 + c_r^2/2",m.dt2,Q*Q+0.5*c_r*c_r,a*a);
      Check(m.nonpositive==0,name,"every re-emitted particle leaves the wall (u.n > 0): %ld do not",m.nonpositive);
    }

    //the diffuse limit: sigma' = sigma = epsilon = 1
    NC::SetAccommodationCoefficients(-1,1.0,1.0,1.0,&error);
    for (int s=0;s<PIC::nTotalSpecies;s++) {
      char name[100];

      snprintf(name,sizeof(name),"nocilla-hs-limit-diffuse[%s]",PIC::MolecularData::GetChemSymbol(s));
      CheckDiffuseMoments(name,Sample(s,vin[s],N),aw[s]);
    }

    //inconsistent coefficients (sigma' = sigma = 0 keeps the mean velocity, epsilon = 0.5 asks for less energy than
    //that): the molecule is re-emitted with the mean velocity A n + Q t, and counted
    {
      const long int before=NC::nInconsistentReemissions;
      double v[3],expect[3],maxerr=0.0;

      NC::SetAccommodationCoefficients(0,0.0,0.0,0.5,&error);
      memcpy(v,vin[0],sizeof(v));
      Interact(0,v);
      FromFrame(expect,vc[0][0],vc[0][1],-vc[0][2]);
      for (int i=0;i<3;i++) maxerr=std::max(maxerr,fabs(v[i]-expect[i])/aw[0]);
      Check((maxerr<1.0E-12)&&(NC::nInconsistentReemissions==before+1),"nocilla-hs-inconsistent",
          "sigma'=sigma=0, epsilon=0.5: re-emitted with the mean velocity (error %.1e a) and counted (%ld)",maxerr,NC::nInconsistentReemissions-before);
    }

    {error.clear(); bool ok=NC::SetAccommodationCoefficients(0,1.5,0.5,0.5,&error); Check((ok==false)&&(error.empty()==false),"nocilla-hs-parameters","sigma' = 1.5 is rejected: %s",error.c_str());}
    {error.clear(); bool ok=NC::SetAccommodationCoefficients(0,0.5,0.5,-0.1,&error); Check((ok==false)&&(error.empty()==false),"nocilla-hs-parameters","epsilon = -0.1 is rejected: %s",error.c_str());}
  }

  //------------------------------------------------------------------------- Schamberg
  {
    std::string error;
    Surface::Model::Set(-1,Surface::Model::QuasiSpecularSchamberg);

    //the beam-shape function against an independent quadrature
    const double phis[5]={1.0E-4,0.3,Pi/4.0,1.2,0.5*Pi};

    for (double p : phis) {
      double mc,mc2;

      BeamQuadrature(p,&mc,&mc2);
      Check(fabs(Surface::QuasiSpecularReflection::Schamberg::BeamShapeFunction(p)-mc)<1.0E-8,"schamberg-beam-shape",
          "Phi(%.4f) = %.10f, quadrature of the beam %.10f",p,Surface::QuasiSpecularReflection::Schamberg::BeamShapeFunction(p),mc);
    }

    Check((fabs(Surface::QuasiSpecularReflection::Schamberg::BeamShapeFunction(0.0)-1.0)<1e-12)&&(fabs(Surface::QuasiSpecularReflection::Schamberg::BeamShapeFunction(0.5*Pi)-2.0/3.0)<1e-8),
        "schamberg-beam-shape","Phi(0) = 1, Phi(pi/2) = 2/3 (the cosine law)");

    //a beam that stays above the surface: theta_i = 51.8 deg, nu = 2 -> theta_r = 38.0 deg, phi0 = 30 deg
    const double al[3]={0.5,0.0,0.9},nu[3]={2.0,1.0,0.7},phi0[3]={30.0,10.0,20.0};

    for (int s=0;s<PIC::nTotalSpecies;s++) {
      const double mass=PIC::MolecularData::GetMass(s),p0=phi0[s%3]*Pi/180.0;
      char name[150];
      double v[3],axis[3],b1[3],b2[3],t[3];
      double vmag=sqrt(Dot(vin[s],vin[s])),lt=sqrt(vc[s][0]*vc[s][0]+vc[s][1]*vc[s][1]);
      double sin_r=pow(lt/vmag,nu[s%3]),cos_r=sqrt(1.0-sin_r*sin_r);
      double Vr=sqrt((1.0-al[s%3])*vmag*vmag+al[s%3]*4.0*Kbol*Tw/mass);

      Surface::QuasiSpecularReflection::Schamberg::SetParameters(s,al[s%3],nu[s%3],p0,&error);
      for (int i=0;i<3;i++) t[i]=(vc[s][0]*e0[i]+vc[s][1]*e1[i])/lt;
      for (int i=0;i<3;i++) axis[i]=cos_r*n[i]+sin_r*t[i],b1[i]=-sin_r*n[i]+cos_r*t[i];
      Cross(b2,n,t);

      double speederr=0.0,maxangle=0.0;
      cStat ca,cb1,cb2;

      for (long int i=0;i<N;i++) {
        memcpy(v,vin[s],sizeof(v));
        Interact(s,v);

        double u=sqrt(Dot(v,v));
        speederr=std::max(speederr,fabs(u-Vr)/Vr);
        double c=Dot(v,axis)/u;
        maxangle=std::max(maxangle,acos(std::min(c,1.0)));
        ca.Add(c),cb1.Add(Dot(v,b1)/u),cb2.Add(Dot(v,b2)/u);
      }

      double mc,mc2;
      BeamQuadrature(p0,&mc,&mc2);

      snprintf(name,sizeof(name),"schamberg[%s,alpha=%g,nu=%g,phi0=%g]",PIC::MolecularData::GetChemSymbol(s),al[s%3],nu[s%3],phi0[s%3]);
      Check(speederr<1.0E-12,name,"|u| = sqrt((1-alpha)|v|^2 + alpha 4kTw/m) = %.6e m/s: largest relative error %.2e",Vr,speederr);
      Check(maxangle<=p0*(1.0+1.0E-12),name,"every direction is within phi0 of the axis (sin(theta_r) = sin(theta_i)^nu = %.6f): largest angle %.6f deg",sin_r,maxangle*180.0/Pi);
      CheckMean(name,"<cos(angle to the axis)> = Phi(phi0)",ca,mc,1.0);
      CheckMean(name,"<direction . b1> (in the plane of incidence)",cb1,0.0,1.0);
      CheckMean(name,"<direction . (n x t)>",cb2,0.0,1.0);
    }

    //specular limit: nu = 1, phi0 = 0, alpha = 0
    Surface::QuasiSpecularReflection::Schamberg::SetParameters(-1,0.0,1.0,0.0,&error);
    double maxerr=0.0;

    for (int s=0;s<PIC::nTotalSpecies;s++) {
      double v[3],expect[3];

      memcpy(v,vin[s],sizeof(v));
      Interact(s,v);
      FromFrame(expect,vc[s][0],vc[s][1],-vc[s][2]);
      for (int i=0;i<3;i++) maxerr=std::max(maxerr,fabs(v[i]-expect[i])/aw[s]);
    }

    Check(maxerr<1.0E-12,"schamberg-limit-specular","nu = 1, phi0 = 0, alpha = 0 gives u = v - 2 (v.n) n: largest error %.2e a",maxerr);

    //quasi-diffuse limit: nu = inf, phi0 = 90 deg, alpha = 1: the cosine law at the speed sqrt(4kTw/m)
    Surface::QuasiSpecularReflection::Schamberg::SetParameters(-1,1.0,HUGE_VAL,0.5*Pi,&error);

    for (int s=0;s<PIC::nTotalSpecies;s++) {
      char name[100];
      cStat c1,c2,ct0,ct1,sp;
      double v[3];

      for (long int i=0;i<N;i++) {
        memcpy(v,vin[s],sizeof(v));
        Interact(s,v);
        double u=sqrt(Dot(v,v));
        c1.Add(Dot(v,n)/u),c2.Add(pow(Dot(v,n)/u,2)),ct0.Add(Dot(v,e0)/u),ct1.Add(Dot(v,e1)/u),sp.Add(u);
      }

      snprintf(name,sizeof(name),"schamberg-limit-cosine-law[%s]",PIC::MolecularData::GetChemSymbol(s));
      CheckMean(name,"<cos(angle to n)> = 2/3",c1,2.0/3.0,1.0);
      CheckMean(name,"<cos^2(angle to n)> = 1/2",c2,0.5,1.0);
      CheckMean(name,"<direction . e0>",ct0,0.0,1.0);
      CheckMean(name,"<direction . e1>",ct1,0.0,1.0);
      Check(fabs(sp.Mean()-sqrt(4.0*Kbol*Tw/PIC::MolecularData::GetMass(s)))<1.0E-9*sp.Mean(),name,"the speed is sqrt(4kTw/m) = sqrt(2) a");
    }

    //a beam that reaches below the surface (grazing incidence, wide beam): only directions above the surface
    Surface::QuasiSpecularReflection::Schamberg::SetParameters(-1,0.2,1.0,60.0*Pi/180.0,&error);
    {
      double v[3],vg[3];
      long int below=0;

      FromFrame(vg,3.0*aw[0],0.5*aw[0],-0.2*aw[0]);
      for (long int i=0;i<N;i++) {
        memcpy(v,vg,sizeof(v));
        Interact(0,v);
        if (Dot(v,n)<=0.0) below++;
      }

      Check(below==0,"schamberg-truncated-beam","with theta_r + phi0 > 90 deg every re-emitted particle leaves the wall: %ld do not",below);
    }

    {error.clear(); bool ok=Surface::QuasiSpecularReflection::Schamberg::SetParameters(0,0.5,0.0,0.3,&error); Check((ok==false)&&(error.empty()==false),"schamberg-parameters","nu = 0 is rejected: %s",error.c_str());}
    {error.clear(); bool ok=Surface::QuasiSpecularReflection::Schamberg::SetParameters(0,0.5,2.0,2.0,&error); Check((ok==false)&&(error.empty()==false),"schamberg-parameters","phi0 = 2 rad is rejected: %s",error.c_str());}
  }

  //------------------------------------------------------------------------- the runtime input block
  {
    std::string error;
    const char *sp=PIC::MolecularData::GetChemSymbol(0);
    std::string all;

    //a valid block: CLL with per-species coefficients
    std::string text="#begin physics: gas_surface_interaction\n  model = cll\n  wall_temperature_k = 250\n";
    for (int s=0;s<PIC::nTotalSpecies;s++) text+=SpeciesBlock(PIC::MolecularData::GetChemSymbol(s),
        "      normal_energy_accommodation = "+std::to_string(0.1+0.1*s)+"\n      tangential_momentum_accommodation = "+std::to_string(0.5+0.2*s)+"\n");
    text+="#end physics\n";

    bool ok=ApplyBlock(text,&error);
    bool installed=ok&&(Surface::Temperature::Isothremal::Temp==250.0);
    for (int s=0;s<PIC::nTotalSpecies;s++) installed&=(Surface::Model::Selected[s]==Surface::Model::CLL)&&
        (fabs(Surface::CLL::NormalPartKineticEnergyAccommodationCoefficient[s]-(0.1+0.1*s))<1e-12)&&(fabs(Surface::CLL::TangentialMomentumAccommodationCoefficient[s]-(0.5+0.2*s))<1e-12);
    Check(installed,"input","a valid block installs the model, the wall temperature and the per-species coefficients%s%s",ok ? "" : ": ",error.c_str());

    //a valid block in the dotted form: Nocilla, angle in degrees
    text="#begin physics: gas_surface_interaction\n  model = nocilla\n";
    for (int s=0;s<PIC::nTotalSpecies;s++) {
      std::string p=std::string("  species.")+PIC::MolecularData::GetChemSymbol(s)+".";
      text+=p+"drift_speed_m_s = 1200\n"+p+"drift_angle_deg = 45\n"+p+"temperature_k = 200\n";
    }
    text+="#end physics\n";
    ok=ApplyBlock(text,&error);
    Check(ok&&(Surface::Model::Selected[0]==Surface::Model::QuasiSpecularNocilla)&&(fabs(Surface::QuasiSpecularReflection::Nocilla::DriftAngle[0]-Pi/4.0)<1e-15),
        "input","the dotted keys 'species.SPECIES.KEY' are accepted; the angle is converted from degrees%s%s",ok ? "" : ": ",error.c_str());

    //the Nocilla model with accommodation coefficients
    text="#begin physics: gas_surface_interaction\n  model = nocilla\n  wall_temperature_k = 280\n";
    for (int s=0;s<PIC::nTotalSpecies;s++) text+=SpeciesBlock(PIC::MolecularData::GetChemSymbol(s),
        "      normal_momentum_accommodation = 0.6\n      tangential_momentum_accommodation = 0.7\n      energy_accommodation = 0.5\n");
    text+="#end physics\n";
    ok=ApplyBlock(text,&error);
    Check(ok&&(Surface::QuasiSpecularReflection::Nocilla::Mode[0]==Surface::QuasiSpecularReflection::Nocilla::AccommodationCoefficients)&&
        (Surface::QuasiSpecularReflection::Nocilla::NormalMomentumAccommodation[0]==0.6)&&(Surface::Temperature::Isothremal::Temp==280.0),
        "input","the Nocilla accommodation coefficients are accepted and installed%s%s",ok ? "" : ": ",error.c_str());

    //invalid blocks: each must be rejected with a message, and leave the installed state unchanged
    struct cCase {const char *what,*block,*expect;};
    std::string s1=std::string("    #begin species: ")+sp+"\n      accommodation = 0.5\n    #end species\n";
    //Nocilla parameters of every species: both sets (direct and accommodation coefficients), or the accommodation ones
    std::string nocillaMixed,nocillaHS;
    for (int s=0;s<PIC::nTotalSpecies;s++) {
      nocillaHS+=SpeciesBlock(PIC::MolecularData::GetChemSymbol(s),"      normal_momentum_accommodation = 0.6\n      tangential_momentum_accommodation = 0.7\n      energy_accommodation = 0.5\n");
      nocillaMixed+=SpeciesBlock(PIC::MolecularData::GetChemSymbol(s),"      normal_momentum_accommodation = 0.6\n      tangential_momentum_accommodation = 0.7\n      energy_accommodation = 0.5\n"
          "      drift_speed_m_s = 1000\n      drift_angle_deg = 45\n      temperature_k = 600\n");
    }

    //every species has the accommodation coefficient, the first one out of range
    for (int s=0;s<PIC::nTotalSpecies;s++) all+=SpeciesBlock(PIC::MolecularData::GetChemSymbol(s),(s==0) ? "      accommodation = 1.5\n" : "      accommodation = 0.5\n");
    std::vector<std::pair<std::string,std::string>> bad={
      {"an unknown model","#begin physics: gas_surface_interaction\n  model = lambertian\n#end physics\n"},
      {"no model","#begin physics: gas_surface_interaction\n  wall_temperature_k = 300\n#end physics\n"},
      {"diffuse without the wall temperature","#begin physics: gas_surface_interaction\n  model = diffuse\n#end physics\n"},
      {"a missing parameter","#begin physics: gas_surface_interaction\n  model = maxwell\n  wall_temperature_k = 300\n#end physics\n"},
      {"a parameter of another model","#begin physics: gas_surface_interaction\n  model = specular\n"+s1+"#end physics\n"},
      {"an out-of-range value","#begin physics: gas_surface_interaction\n  model = maxwell\n  wall_temperature_k = 300\n"+all+"#end physics\n"},
      {"a value that is not a number","#begin physics: gas_surface_interaction\n  model = maxwell\n  wall_temperature_k = 300\n"+SpeciesBlock(sp,"      accommodation = half\n")+"#end physics\n"},
      {"an unknown species","#begin physics: gas_surface_interaction\n  model = specular\n  #begin species: Xe\n    accommodation = 0.5\n  #end species\n#end physics\n"},
      {"a key given twice (block and dotted form)","#begin physics: gas_surface_interaction\n  model = maxwell\n  wall_temperature_k = 300\n  species."+std::string(sp)+".accommodation = 0.5\n"+s1+"#end physics\n"},
      {"a negative wall temperature","#begin physics: gas_surface_interaction\n  model = specular\n  wall_temperature_k = -3\n#end physics\n"},
      {"an unknown key","#begin physics: gas_surface_interaction\n  model = specular\n  accomodation = 0.5\n#end physics\n"},
      {"an unknown child block","#begin physics: gas_surface_interaction\n  model = specular\n  #begin wall\n    t = 1\n  #end wall\n#end physics\n"},
      {"Nocilla with both parameter sets","#begin physics: gas_surface_interaction\n  model = nocilla\n  wall_temperature_k = 300\n"+nocillaMixed+"#end physics\n"},
      {"Nocilla accommodation coefficients without the wall temperature","#begin physics: gas_surface_interaction\n  model = nocilla\n"+nocillaHS+"#end physics\n"},
    };

    for (auto &b : bad) {
      int before[PIC::nTotalSpecies];
      double beforeT=Surface::Temperature::Isothremal::Temp,beforeU=Surface::QuasiSpecularReflection::Nocilla::DriftSpeed[0];
      int beforeMode=Surface::QuasiSpecularReflection::Nocilla::Mode[0];

      for (int s=0;s<PIC::nTotalSpecies;s++) before[s]=Surface::Model::Selected[s];
      error.clear();
      ok=ApplyBlock(b.second,&error);

      bool unchanged=(Surface::Temperature::Isothremal::Temp==beforeT)&&(Surface::QuasiSpecularReflection::Nocilla::DriftSpeed[0]==beforeU)&&
          (Surface::QuasiSpecularReflection::Nocilla::Mode[0]==beforeMode);
      for (int s=0;s<PIC::nTotalSpecies;s++) unchanged&=(Surface::Model::Selected[s]==before[s]);
      Check((ok==false)&&(error.empty()==false)&&unchanged,"input-rejected","%s is rejected with a message, nothing installed: %s",b.first.c_str(),error.c_str());
    }
  }

  fprintf(Report,"INFO %d checks passed, %d failed\n",nPassed,nFailed);
  printf("INFO %d checks passed, %d failed\n",nPassed,nFailed);
  fclose(Report);

  MPI_Finalize();
  return (nFailed==0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
