// Float-vs-double experiment for legacy-tracer-style backtracing in a centered dipole.
// Normalized units: length Re, speed c, time Re/c, momentum m_p c, charge e, and
// magnetic field B_u = m_p c/(e Re).  In these units dp/dt = q (v x B), R[GV]=p*0.938.
//
// Template parameters:
//   S = precision of the state (x,p) and of all field/mover arithmetic
//   A = precision of the accumulators (trace time, path length)
//   X = precision of the position only (lets x be double while p/B are float)
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <vector>
#include <string>

static const double kRe=6371.0e3, kC=299792458.0, kE=1.602176634e-19;
static const double kMp=1.007276466621*1.66053906660e-27;
static const double kMpc2_GeV=kMp*kC*kC/kE*1e-9;          // 0.938272 GeV
static const double kTimeUnit=kRe/kC;                      // s
static const double kBunit=kMp*kC/(kE*kRe);                // T
static const double kBeq=3.12e-5;                          // T (gridless/DipoleInterface.h)

template<class T> struct V { T x,y,z; };
template<class T> inline V<T> operator+(V<T>a,V<T>b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
template<class T> inline V<T> operator-(V<T>a,V<T>b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
template<class T> inline V<T> operator*(T s,V<T>a){return {s*a.x,s*a.y,s*a.z};}
template<class T> inline T dot(V<T>a,V<T>b){return a.x*b.x+a.y*b.y+a.z*b.z;}
template<class T> inline V<T> cross(V<T>a,V<T>b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
template<class T,class U> inline V<T> cvt(V<U>a){return {T(a.x),T(a.y),T(a.z)};}

enum Mover { BORIS, HC4, RK4 };
enum Term { ALLOWED=0, INNER=1, TIMELIMIT=2, STEPLIMIT=3, DISTLIMIT=4 };

struct Config {
  double box=35.0, rInner=1.0;            // Re
  double tMax_s=600.0, dtMax_s=1.0;       // s
  double distMax_Re=400.0;                // Re
  long maxSteps=300000L;           // default MAX_STEPS
  double gyro=0.15;
  Mover mover=BORIS;
  double dtScale=1.0;                // calibration perturbation of the step
};

// B in normalized units, moment along +z as in gridless/DipoleInterface.h.
template<class S> inline V<S> Dipole(V<S> x) {
  const S B0=S(kBeq/kBunit);              // equatorial surface field
  const S r2=dot(x,x), r=std::sqrt(r2), r5=r2*r2*r;
  const S mz=B0;                           // |m| mu0/4pi / Re^3 in B units
  const S mdotr=mz*x.z;
  const S inv3=S(1)/(r2*r);
  return { S(3)*x.x*mdotr/r5, S(3)*x.y*mdotr/r5, S(3)*x.z*mdotr/r5-mz*inv3 };
}

template<class S> inline V<S> Vel(V<S> p){ const S g=std::sqrt(S(1)+dot(p,p)); return (S(1)/g)*p; }

template<class S> inline void Rotate(V<S>& p,V<S> B,S q,S dt){
  const S g=std::sqrt(S(1)+dot(p,p));
  const V<S> t=(q*dt/(S(2)*g))*B;
  const S t2=dot(t,t);
  const V<S> s=(S(2)/(S(1)+t2))*t;
  const V<S> pp=p+cross(p,t);
  p=p+cross(pp,s);
}

// Legacy BorisStep: B at x_n, half drift, rotate, half drift.
template<class S,class X> inline void BorisStep(V<X>& x,V<S>& p,S q,S dt){
  const V<S> B=Dipole<S>(cvt<S>(x));
  x=x+cvt<X>((S(0.5)*dt)*Vel(p));
  Rotate(p,B,q,dt);
  x=x+cvt<X>((S(0.5)*dt)*Vel(p));
}
template<class S,class X> inline void H2(V<X>& x,V<S>& p,S q,S dt){
  x=x+cvt<X>((S(0.5)*dt)*Vel(p));
  const V<S> B=Dipole<S>(cvt<S>(x));
  Rotate(p,B,q,dt);
  x=x+cvt<X>((S(0.5)*dt)*Vel(p));
}
template<class S,class X> inline void HC4Step(V<X>& x,V<S>& p,S q,S dt){
  const S c=S(1.25992104989487316477), w1=S(1)/(S(2)-c), w0=-c/(S(2)-c);
  H2(x,p,q,w1*dt); H2(x,p,q,w0*dt); H2(x,p,q,w1*dt);
}
template<class S,class X> inline void RK4Step(V<X>& x,V<S>& p,S q,S dt){
  auto f=[&](V<S> xx,V<S> pp,V<S>& dx,V<S>& dp){ const V<S> v=Vel(pp); dx=v; dp=q*cross(v,Dipole<S>(xx)); };
  const V<S> x0=cvt<S>(x); V<S> k1x,k1p,k2x,k2p,k3x,k3p,k4x,k4p;
  f(x0,p,k1x,k1p);
  f(x0+(S(0.5)*dt)*k1x,p+(S(0.5)*dt)*k1p,k2x,k2p);
  f(x0+(S(0.5)*dt)*k2x,p+(S(0.5)*dt)*k2p,k3x,k3p);
  f(x0+dt*k3x,p+dt*k3p,k4x,k4p);
  const S h=dt/S(6);
  x=x+cvt<X>(h*(k1x+S(2)*k2x+S(2)*k3x+k4x));
  p=p+h*(k1p+S(2)*k2p+S(2)*k3p+k4p);
}

struct Result { int term; long steps; double t_s; double pSpread; double dtMin_s=1e300; };
static double gDtMinGlobal=1e300; static double gSteps=0,gTime=0;

// Structured (upper-bound-only) timestep; segment checks against inner sphere; box escape.
template<class S,class A,class X>
Result Trace(const Config& cfg,V<double> x0,V<double> dir,double R_GV,double q){
  V<X> x=cvt<X>(x0);
  const double pmag=R_GV/kMpc2_GeV;                  // units of m_p c for |q|=1
  V<S> p=cvt<S>((pmag)*dir);
  const S qs=S(q);
  const A tMax=A(cfg.tMax_s/kTimeUnit), dtMax=A(cfg.dtMax_s/kTimeUnit);
  A t=0, s=0; long n=0;
  const double p0=std::sqrt(dot(cvt<double>(p),cvt<double>(p)));
  double pmin=p0,pmax=p0;
  const double r2in=cfg.rInner*cfg.rInner;
  while (true) {
    if (n>=cfg.maxSteps) return {STEPLIMIT,n,double(t)*kTimeUnit,(pmax-pmin)/p0};
    if (!(t<tMax)) return {TIMELIMIT,n,double(t)*kTimeUnit,(pmax-pmin)/p0};
    if (cfg.distMax_Re>0 && !(s<A(cfg.distMax_Re))) return {DISTLIMIT,n,double(t)*kTimeUnit,(pmax-pmin)/p0};
    // dt: min(DT_TRACE, gyro/Omega, remaining)
    const V<S> Bh=Dipole<S>(cvt<S>(x));
    const S g=std::sqrt(S(1)+dot(p,p));
    const S om=std::fabs(qs)*std::sqrt(dot(Bh,Bh))/g;
    A dt=A(double(dtMax)*cfg.dtScale);
    if (om>S(0)) dt=std::min(dt,A(S(cfg.gyro*cfg.dtScale)/om));
    dt=std::min(dt,tMax-t);
    if (double(dt)*kTimeUnit<gDtMinGlobal && double(dt)<double(tMax-t)) gDtMinGlobal=double(dt)*kTimeUnit;
    const V<X> xb=x;
    switch (cfg.mover) {
      case BORIS: BorisStep<S,X>(x,p,qs,S(dt)); break;
      case HC4: { const int nc=4; for(int i=0;i<nc;++i){ HC4Step<S,X>(x,p,qs,S(dt)/S(nc)); } } break;
      case RK4: RK4Step<S,X>(x,p,qs,S(dt)); break;
    }
    // segment / inner sphere (double arithmetic on the converted endpoints)
    const V<double> a=cvt<double>(xb), b=cvt<double>(x), d=b-a;
    const double dd=dot(d,d);
    double tc=dd>0? -dot(a,d)/dd : 0.0; tc=std::max(0.0,std::min(1.0,tc));
    const V<double> qc=a+tc*d;
    const A t_old=t;
    t=t+dt; s=s+A(std::sqrt(dd)); ++n;
    (void)t_old;
    const double pm=std::sqrt(dot(cvt<double>(p),cvt<double>(p))); pmin=std::min(pmin,pm); pmax=std::max(pmax,pm);
    if (dot(qc,qc)<=r2in) return {INNER,n,double(t)*kTimeUnit,(pmax-pmin)/p0};
    if (std::fabs(b.x)>cfg.box||std::fabs(b.y)>cfg.box||std::fabs(b.z)>cfg.box)
      return {ALLOWED,n,double(t)*kTimeUnit,(pmax-pmin)/p0};
  }
}

static double RigidityFromMeV(double E){ const double m=kMpc2_GeV*1e3; return std::sqrt(E*(E+2*m))*1e-3; }
static double Stormer(double latDeg,double rRe){ const double R0=0.299792458*0.25*kBeq*kRe; double c=std::cos(latDeg*M_PI/180); return R0*c*c*c*c/(rRe*rRe); }

template<class S,class A,class X>
double UpperScan(const Config& cfg,V<double> x0,V<double> arrival,int nScan,long& traces){
  const V<double> v0=(-1.0)*arrival;
  const double Rmin=RigidityFromMeV(1.0), Rmax=RigidityFromMeV(20000.0);
  std::vector<double> g(nScan);
  for(int i=0;i<nScan;++i) g[i]=std::exp(std::log(Rmin)+(std::log(Rmax)-std::log(Rmin))*i/(nScan-1.0));
  auto allowed=[&](double R){ ++traces; return Trace<S,A,X>(cfg,x0,v0,R,1.0).term==ALLOWED; };
  if (!allowed(g.back())) return -1;
  for (int i=nScan-2;i>=0;--i) if (!allowed(g[i])) { double lo=g[i],hi=g[i+1];
      while (hi-lo>1e-3) { double m=0.5*(lo+hi); if (allowed(m)) hi=m; else lo=m; } return hi; }
  return Rmin;
}

struct VariantName { const char* name; };

template<class S,class A,class X>
void RunCutoffs(const char* label,const Config& cfg,std::vector<double>& out,double& sec){
  auto t0=std::chrono::steady_clock::now();
  long traces=0;
  const double lats[]={0,30,45,60}, alts[]={500,9000};
  out.clear();
  for(double alt:alts) for(double lat:lats){
    const double r=1.0+alt*1e3/kRe, la=lat*M_PI/180;
    V<double> x{r*std::cos(la),0,r*std::sin(la)}; V<double> up=(1.0/r)*x;
    out.push_back(UpperScan<S,A,X>(cfg,x,(-1.0)*up,100,traces));
  }
  sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
  (void)label;
}

// Directional transmission at one point: equal solid angle grid, energies around Rc_vert.
template<class S,class A,class X>
void RunTransmission(const Config& cfg,double latDeg,double altKm,const std::vector<double>& R,
                     std::vector<int>& states,std::vector<double>& T,double& sec,double& maxSpread){
  auto t0=std::chrono::steady_clock::now();
  const int nMu=24,nPhi=48; const double r=1.0+altKm*1e3/kRe, la=latDeg*M_PI/180;
  V<double> x{r*std::cos(la),0,r*std::sin(la)};
  states.assign(R.size()*nMu*nPhi,0); T.assign(R.size(),0); maxSpread=0;
  for(size_t ie=0;ie<R.size();++ie){ int nAllowed=0;
    for(int i=0;i<nMu;++i){ const double mu=-1+2*(i+0.5)/nMu, st=std::sqrt(1-mu*mu);
      for(int j=0;j<nPhi;++j){ const double ph=2*M_PI*(j+0.5)/nPhi; V<double> d{st*std::cos(ph),st*std::sin(ph),mu};
        Result res=Trace<S,A,X>(cfg,x,(-1.0)*d,R[ie],1.0); gSteps+=res.steps; gTime+=res.t_s;
        states[(ie*nMu+i)*nPhi+j]=res.term; if(res.term==ALLOWED) ++nAllowed;
        if (res.term!=TIMELIMIT && res.term!=STEPLIMIT) maxSpread=std::max(maxSpread,res.pSpread);
      }}
    T[ie]=double(nAllowed)/(nMu*nPhi);
  }
  sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
}

int main(int argc,char** argv){
  Config cfg;
  std::string mv=argc>1?argv[1]:"BORIS";
  cfg.mover= mv=="HC4"?HC4: (mv=="RK4"?RK4:BORIS);
  printf("mover=%s  Bunit=%.4e T  time unit=%.4e s  B_eq'=%.4f  units check R0=%.4f GV\n",
         mv.c_str(),kBunit,kTimeUnit,kBeq/kBunit,0.299792458*0.25*kBeq*kRe);

  // ---- A: vertical UPPER_SCAN cutoffs ------------------------------------------------
  std::vector<double> dd,fd,ff,xd; double sdd,sfd,sff,sxd;
  RunCutoffs<double,double,double>("DD",cfg,dd,sdd);
  RunCutoffs<float,double,float>("FD",cfg,fd,sfd);
  RunCutoffs<float,double,double>("XD",cfg,xd,sxd);
  RunCutoffs<float,float,float>("FF",cfg,ff,sff);
  std::vector<double> dp; double sdp; Config cp=cfg; cp.dtScale=1.0+1.0e-7;
  RunCutoffs<double,double,double>("DDp",cp,dp,sdp);
  printf("\nA) vertical UPPER_SCAN cutoff [GV] (N=100 log grid + 1 MV bisection)\n");
  printf("%6s %5s %9s %9s %9s %9s %9s %9s | rel(DDp-DD) rel(FD-DD) rel(XD-DD) rel(FF-DD)\n","alt","lat","Stormer","DD","DDp","FD","XD","FF");
  const double lats[]={0,30,45,60}, alts[]={500,9000}; int k=0;
  for(double alt:alts) for(double lat:lats){ const double r=1+alt*1e3/kRe;
    auto rel=[&](double a){ return dd[k]>0? (a-dd[k])/dd[k] : 0.0; };
    printf("%6.0f %5.0f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f | %+10.2e %+10.2e %+10.2e %+10.2e\n",alt,lat,Stormer(lat,r),dd[k],dp[k],fd[k],xd[k],ff[k],rel(dp[k]),rel(fd[k]),rel(xd[k]),rel(ff[k])); ++k; }
  printf("CPU s: DD=%.2f FD=%.2f XD=%.2f FF=%.2f\n",sdd,sfd,sxd,sff);

  // ---- B: directional transmission in the penumbra ------------------------------------
  const double lat=45, alt=9000, rc=Stormer(lat,1+alt*1e3/kRe);
  std::vector<double> R; for(int i=0;i<12;++i) R.push_back(rc*(0.4+0.15*i));
  std::vector<int> sD,sF,sX; std::vector<double> TD,TF,TX; double tD,tF,tX,spD,spF,spX;
  gDtMinGlobal=1e300; gSteps=0; gTime=0;
  RunTransmission<double,double,double>(cfg,lat,alt,R,sD,TD,tD,spD);
  printf("\nE) adaptive vs global fixed dt on the DD transmission workload: steps=%.3e, traced time=%.3e s, min gyro dt=%.3e s\n"
         "   a fixed dt equal to the smallest gyro-limited step would need %.3e steps (%.0fx)\n",
         gSteps,gTime,gDtMinGlobal,gTime/gDtMinGlobal,gTime/gDtMinGlobal/gSteps);
  RunTransmission<float,double,float>(cfg,lat,alt,R,sF,TF,tF,spF);
  RunTransmission<float,double,double>(cfg,lat,alt,R,sX,TX,tX,spX);
  std::vector<int> sP; std::vector<double> TP; double tP,spP; Config cp2=cfg; cp2.dtScale=1.0+1.0e-7;
  RunTransmission<double,double,double>(cp2,lat,alt,R,sP,TP,tP,spP);
  printf("\nB) 24x48 directional transmission at lat=%g alt=%g km (Rc_vert=%.4f GV)\n",lat,alt,rc);
  printf("%9s %8s %8s %8s %8s %9s %9s %9s\n","R/Rc","T_DD","T_DDp","T_FD","T_XD","mism_DDp","mism_FD","mism_XD");
  long mF=0,mX=0,mP=0,nTL=0;
  for(size_t ie=0;ie<R.size();++ie){ int a=0,b=0,c=0; for(int j=0;j<24*48;++j){ size_t q=ie*1152+j;
      if((sD[q]==ALLOWED)!=(sF[q]==ALLOWED)) ++a; if((sD[q]==ALLOWED)!=(sX[q]==ALLOWED)) ++b; if((sD[q]==ALLOWED)!=(sP[q]==ALLOWED)) ++c; if(sD[q]==TIMELIMIT||sD[q]==STEPLIMIT) ++nTL; }
    mF+=a; mX+=b; mP+=c; printf("%9.3f %8.4f %8.4f %8.4f %8.4f %9d %9d %9d\n",R[ie]/rc,TD[ie],TP[ie],TF[ie],TX[ie],c,a,b); }
  printf("total direction mismatches vs DD: DDp=%ld FD=%ld XD=%ld of %zu ; DD time/step-limited=%ld\n",mP,mF,mX,R.size()*1152,nTL);
  printf("max |p| relative spread on resolved traces: DD=%.2e FD=%.2e XD=%.2e (TRAP_ENERGY_REL_TOL default 1e-4)\n",spD,spF,spX);
  printf("CPU s: DD=%.2f FD=%.2f XD=%.2f\n",tD,tF,tX);

  // ---- C: accumulator stall and long-orbit invariants ---------------------------------
  Config c2=cfg; c2.tMax_s=1200; c2.distMax_Re=0; c2.maxSteps=100000000L;
  const double R10=RigidityFromMeV(10.0); V<double> x4{4,0,0}; V<double> d4{0,std::cos(1.0),std::sin(1.0)};
  Result a=Trace<double,double,double>(c2,x4,d4,R10,1), b=Trace<float,double,float>(c2,x4,d4,R10,1), c=Trace<float,float,float>(c2,x4,d4,R10,1);
  printf("\nC) trapped 10 MeV proton at 4 Re, tMax=1200 s, no distance cap\n");
  printf("   DD: term=%d steps=%ld t=%.3f s dp/p=%.2e\n",a.term,a.steps,a.t_s,a.pSpread);
  printf("   FD: term=%d steps=%ld t=%.3f s dp/p=%.2e\n",b.term,b.steps,b.t_s,b.pSpread);
  printf("   FF: term=%d steps=%ld t=%.3f s dp/p=%.2e  (term 2=TIME,3=STEP)\n",c.term,c.steps,c.t_s,c.pSpread);
  // ---- D: time accumulator alone (sum a near-Earth gyro step until 1200 s) -------------
  { const double dts=5.0e-5/kTimeUnit, tEnd=1200.0/kTimeUnit; float tf=0; double td=0; long n=0;
    while (td<tEnd && n<100000000L) { tf+=float(dts); td+=dts; ++n; }
    printf("\nD) %ld steps of dt=5e-5 s: double t=%.3f s, float t=%.3f s (float stalls when dt<ulp(t)/2)\n",n,td*kTimeUnit,double(tf)*kTimeUnit); }
  return 0;
}
