#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>

#include "Liouville.h"

void cLiouville::Init(const cInput& input) {
  SurfaceBC=input.GetOption("SurfaceBC",{"flat","rough"},"flat");
  SpacecraftVelocity=input.GetBool("SpacecraftVelocity",true);
  SubdivisionAccuracy=input.GetDouble("SubdivisionAccuracy",0.05);
  MaxSubdivisionLevel=input.GetInt("MaxSubdivisionLevel",8);

  FlightTimeRotation=input.GetOption("FlightTimeRotation",{"none","exact"},"none");
  nSpeedNodes=input.GetInt("SpeedNodes",16);
  PathTolerance=input.GetDouble("PathTolerance",0.5);

  //Gauss-Legendre nodes and weights on [0,xmax] (Newton iteration on the Legendre polynomial)
  double xmax=input.GetDouble("SpeedRangeMax",5.5);
  int N=nSpeedNodes;

  SpeedNode.resize(N),SpeedWeight.resize(N);

  for (int i=0;i<N;i++) {
    double z=cos(Const::Pi*(i+0.75)/(N+0.5)),pp=0.0,dz;

    do {
      double p1=1.0,p2=0.0;
      for (int j=1;j<=N;j++) {
        double p3=p2;
        p2=p1;
        p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      pp=N*(z*p1-p2)/(z*z-1.0);
      dz=p1/pp;
      z-=dz;
    }
    while (fabs(dz)>1.0E-15);

    SpeedNode[i]=0.5*xmax*(1.0-z);
    SpeedWeight[i]=xmax/((1.0-z*z)*pp*pp);
  }

  printf("Liouville: surface BC=%s, spacecraft velocity=%s, subdivision: element size <= %g x distance, max level %i\n",
      SurfaceBC.c_str(),SpacecraftVelocity ? "on" : "off",SubdivisionAccuracy,MaxSubdivisionLevel);
  printf("Liouville: flight-time rotation=%s (speed nodes %i on x=beta v in [0,%g], path tolerance %g m)\n",FlightTimeRotation.c_str(),nSpeedNodes,xmax,PathTolerance);
}

//rotation of vector v by angle phi about the unit axis a (Rodrigues)
static inline Vec3 Rotate(const Vec3& v,const Vec3& a,double phi) {
  double c=cos(phi),s=sin(phi);
  return v*c+Cross(a,v)*s+a*(Dot(a,v)*(1.0-c));
}

bool PointInsideSurface(const cMesh& mesh,const cBVH& bvh,const Vec3& x) {
  //count the intersections along a fixed oblique direction (avoids alignment with mesh edges)
  Vec3 d=Normalized(Vec3(0.5773502691896,0.5773502691897,0.5773502691898));
  int cnt=0,ignore=-1;
  double t,tmin=0.0;

  while (true) {
    int iface=bvh.FirstHit(x,d,tmin,1.0E30,t,ignore);
    if (iface<0) break;

    cnt++,tmin=t,ignore=iface;
  }

  return (cnt%2)==1;
}

void cLiouville::EvaluateKernels(const cMesh& mesh,const cBVH& bvh,const cGaugeModel& gauges,const cObservationGeometry& g,cLiouvilleResult& res) const {
  const int nFaces=(int)mesh.Faces.size();
  const Vec3 xg=g.xSC;
  const bool flat=(SurfaceBC=="flat");

  std::vector<double> Kn(nFaces,0.0),KPhi(nFaces,0.0),KV(nFaces,0.0),KnBelow(nFaces,0.0),dist(nFaces,0.0);
  const Vec3 Zsc=g.RG.Boresight;   //RG boresight = +Z of the spacecraft (ROS_ROSINA_COPS is co-aligned with ROS_SPACECRAFT)

  res.GaugeInsideNucleus=PointInsideSurface(mesh,bvh,xg);
  res.VdotB_RG=SpacecraftVelocity ? Dot(g.vSC,g.RG.Boresight) : 0.0;

  #pragma omp parallel for schedule(dynamic,64)
  for (int iface=0;iface<nFaces;iface++) {
    const cTriangle& f=mesh.Faces[iface];
    Vec3 a=mesh.Nodes[f.node[0]],b=mesh.Nodes[f.node[1]],c=mesh.Nodes[f.node[2]];

    dist[iface]=PointTriangleDistance(xg,a,b,c);

    //the gauge must be in front of the plane of the element
    if (Dot(f.Normal,xg-a)<=0.0) continue;

    //subdivision level: sub-element size <= SubdivisionAccuracy x distance to the element
    int level=0;
    double h=f.MaxEdge;

    while ((level<MaxSubdivisionLevel)&&(h>SubdivisionAccuracy*dist[iface])) h*=0.5,level++;

    int n=1<<level;
    double subArea=f.Area/((double)n*n);
    Vec3 ea=(b-a)/n,eb=(c-a)/n;
    double kn=0.0,kphi=0.0,kv=0.0,knb=0.0;

    auto Sample=[&](const Vec3& p) {
      Vec3 d=p-xg;          //from the gauge to the sub-element
      double t;

      //line of sight: the first surface crossed by the ray from the gauge must be this element
      if (bvh.FirstHit(xg,d,0.0,1.0+1.0E-9,t)!=iface) return;

      double r2=Dot(d,d),r=sqrt(r2);
      Vec3 l=-d/r;           //molecule velocity direction (surface -> gauge)
      double cosTheta=Dot(l,f.Normal);

      if (cosTheta<=0.0) return;

      double w=subArea/r2*(flat ? cosTheta/Const::Pi : 1.0/(2.0*Const::Pi));

      double wng=w*gauges.NGWeight(l,g);

      kn+=wng;
      if (Dot(l,Zsc)>0.0) knb+=wng;

      double acc=gauges.RGAcceptance(l,g);
      if (acc>0.0) {
        kphi+=w*acc*(-Dot(l,g.RG.Boresight));
        kv+=w*acc;
      }
    };

    //regular subdivision into n^2 congruent sub-triangles (equivalent to 'level' midpoint refinements)
    for (int i=0;i<n;i++) for (int j=0;i+j<n;j++) {
      Sample(a+ea*(i+1.0/3.0)+eb*(j+1.0/3.0));
      if (i+j<n-1) Sample(a+ea*(i+2.0/3.0)+eb*(j+2.0/3.0));
    }

    Kn[iface]=kn,KPhi[iface]=kphi,KV[iface]=kv,KnBelow[iface]=knb;
  }

  res.Altitude=*std::min_element(dist.begin(),dist.end());
  res.Kernels.clear();

  for (int i=0;i<nFaces;i++) if ((Kn[i]!=0.0)||(KPhi[i]!=0.0)||(KV[i]!=0.0)) res.Kernels.push_back({i,Kn[i],KPhi[i],KV[i],KnBelow[i]});

  res.nVisibleElements=(int)res.Kernels.size();
}

void cLiouville::ApplySource(const cMesh& mesh,const cGaugeModel& gauges,const std::vector<cSurfaceState>& state,double TemperatureRG,cLiouvilleResult& res) const {
  res.TemperatureRG=TemperatureRG;
  double nTot=0.0,nBelow=0.0;

  for (int s=0;s<Const::nSpecies;s++) {
    double n=0.0,phi=0.0,Q=0.0;

    for (auto& k : res.Kernels) {
      const cSurfaceState& st=state[k.Face];
      double beta=sqrt(Const::SpeciesMass[s]/(2.0*Const::kBoltzmann*st.Temperature));

      n+=st.Flux[s]*beta*Const::SqrtPi/2.0*k.Kn;
      nBelow+=st.Flux[s]*beta*Const::SqrtPi/2.0*k.KnBelowDeck;
      phi+=st.Flux[s]*(k.KPhi+res.VdotB_RG*beta*Const::SqrtPi/2.0*k.KV);
    }

    for (size_t i=0;i<mesh.Faces.size();i++) Q+=state[i].Flux[s]*mesh.Faces[i].Area;

    res.n[s]=n;
    res.PhiIn[s]=phi;
    res.pNG[s]=gauges.NGPressure(s,n);
    res.pRG[s]=gauges.RGPressure(s,phi,TemperatureRG);
    res.TotalProduction[s]=Q;
    nTot+=n;
  }

  res.NGBelowDeckFraction=(nTot>0.0) ? nBelow/nTot : 0.0;
}

//=====================================================================================================
//Rotation of the nucleus during the flight of the molecules (FlightTimeRotation=exact)
//
//Molecules move on straight lines in the inertial frame while the nucleus rotates with the angular velocity
//omega (constant during a flight of <= tens of seconds). Work in the body-fixed frame at the reception time t
//(an inertial frame frozen at t). A molecule received at t with speed v left the body-fixed point p at t-tau:
//   p'=R(-omega tau) p,  tau=|x_g-p'|/v (iterated),  u=(x_g-p')/L,  L=|x_g-p'|,  n'=R(-omega tau) n.
//The emitting surface moves with u_s=omega x p'; the molecule velocity relative to the surface is w=v u-u_s and the
//distribution function is Galilean invariant: f_gauge(v u)=f_surface(w). Changing the integration variable from
//the arrival direction to the emitting area at fixed v gives (Jacobian of the retarded-time mapping)
//   dOmega = dS (w.n') / [v L^2 (1-u.u_s/v)],
//so that, with the emitted distributions f_flat=(2 F beta^4/pi) exp(-beta^2 w^2) (w.n'>0) and
//f_rough=F h(w)/(2 pi w^3 cos(theta_w)),
//   flat : dn = dS/L^2 * (w.n') v/(1-u.u_s/v) * (2 beta^4/pi) F exp(-beta^2 w^2) dv
//   rough: dn = dS/L^2 * |w|    v/(1-u.u_s/v) * (  beta^4/pi) F exp(-beta^2 w^2) dv
//which reduce to the static expressions for omega=0. The flux into the RG aperture is dn (-(v u-V).b), and the
//gauge acceptances use the direction of the velocity relative to the spacecraft, v u-V.
//Obstruction is tested along the image of the straight inertial path in the rotating body frame,
//   q(s)=R(+omega tau (1-s)) [p'+s (x_g-p')],  0<=s<=1  (q(0)=p, q(1)=x_g),
//with straight segments of sagitta <= PathTolerance. The source (flux, temperature) is evaluated at the emission
//time t-tau: the Sun direction in the body frame is R(+omega tau) s(t); the cast-shadow state is taken at t.
//The speed integral is evaluated with the Gauss-Legendre rule in x=beta(t) v.
void cLiouville::EvaluateExact(const cMesh& mesh,const cBVH& bvh,const cGaugeModel& gauges,const cSourceModel& source,
    const std::vector<cSurfaceState>& state,const cObservationGeometry& g,double TemperatureRG,cLiouvilleResult& res) const {
  const int nFaces=(int)mesh.Faces.size();
  const Vec3 xg=g.xSC;
  const bool flat=(SurfaceBC=="flat");
  const double w=Length(g.Omega);
  const Vec3 axis=(w>0.0) ? g.Omega/w : Vec3(0.0,0.0,1.0);
  const Vec3 Zsc=g.RG.Boresight;
  const Vec3 Vsc=SpacecraftVelocity ? g.vSC : Vec3(0.0,0.0,0.0);
  const double VdotB=Dot(Vsc,g.RG.Boresight);
  const Vec3 bRG=g.RG.Boresight;

  std::vector<double> dist(nFaces,0.0);
  std::vector<cElementKernelExact> K(nFaces);
  std::vector<double> nF(nFaces*Const::nSpecies,0.0),phiF(nFaces*Const::nSpecies,0.0),nBelow(nFaces,0.0);

  res.GaugeInsideNucleus=PointInsideSurface(mesh,bvh,xg);
  res.VdotB_RG=VdotB;
  res.TemperatureRG=TemperatureRG;

  #pragma omp parallel for schedule(dynamic,64)
  for (int iface=0;iface<nFaces;iface++) {
    const cTriangle& f=mesh.Faces[iface];
    Vec3 a=mesh.Nodes[f.node[0]],b=mesh.Nodes[f.node[1]],c=mesh.Nodes[f.node[2]];
    cElementKernelExact& k=K[iface];

    k.Face=iface;
    for (int s=0;s<Const::nSpecies;s++) k.Gn[s]=0.0,k.GPhi[s]=0.0;

    dist[iface]=PointTriangleDistance(xg,a,b,c);

    //the rotation during the flight changes the orientation by <<0.1 rad: elements clearly behind their plane cannot contribute
    Vec3 dc=xg-f.Centroid;
    if (Dot(f.Normal,dc)<-0.1*Length(dc)) continue;

    int level=0;
    double h=f.MaxEdge;

    while ((level<MaxSubdivisionLevel)&&(h>SubdivisionAccuracy*dist[iface])) h*=0.5,level++;

    int n=1<<level;
    double subArea=f.Area/((double)n*n);
    Vec3 ea=(b-a)/n,eb=(c-a)/n;

    //state of the element at the reception time (temperature defines the speed nodes; shadow state)
    const cSurfaceState& st=state[iface];
    double betaT[Const::nSpecies];

    for (int s=0;s<Const::nSpecies;s++) betaT[s]=sqrt(Const::SpeciesMass[s]/(2.0*Const::kBoltzmann*st.Temperature));

    auto Sample=[&](const Vec3& p) {
      for (int s=0;s<Const::nSpecies;s++) for (int iv=0;iv<nSpeedNodes;iv++) {
        double v=SpeedNode[iv]/betaT[s],dv=SpeedWeight[iv]/betaT[s];
        double tau=Length(xg-p)/v;
        Vec3 pp=p;

        for (int it=0;it<6;it++) {
          pp=Rotate(p,axis,-w*tau);
          tau=Length(xg-pp)/v;
        }

        Vec3 d=xg-pp;
        double L=Length(d);
        Vec3 u=d/L;
        Vec3 nn=Rotate(f.Normal,axis,-w*tau);
        Vec3 us=Cross(g.Omega,pp);           //velocity of the emitting surface
        Vec3 wv=u*v-us;                      //molecule velocity relative to the surface
        double wn=Dot(wv,nn);

        if (wn<=0.0) continue;

        //obstruction along the image of the straight inertial path in the rotating body frame
        double phi=w*tau;
        int nseg=(phi*L>4.0*PathTolerance) ? (int)ceil(sqrt(phi*L/(4.0*PathTolerance))) : 1;
        bool blocked=false;

        if (nseg>256) nseg=256;

        Vec3 q0=p;
        for (int is=1;(is<=nseg)&&(!blocked);is++) {
          double sg=(double)is/nseg;
          Vec3 q1=(is==nseg) ? xg : Rotate(pp+d*sg,axis,phi*(1.0-sg));

          if (bvh.AnyHit(q0,q1-q0,1.0E-9,1.0,iface)) blocked=true;
          q0=q1;
        }

        if (blocked) continue;

        //source at the emission time
        Vec3 sunE=Rotate(g.SunDirection,axis,phi);
        double cosSZA=Dot(f.Normal,sunE),F[Const::nSpecies],T;

        source.EvaluateElement(iface,cosSZA,(cosSZA<=0.0)||st.Shadow,g.HeliocentricDistance,F,T);

        double beta=sqrt(Const::SpeciesMass[s]/(2.0*Const::kBoltzmann*T)),b4=pow(beta,4);
        double w2=Dot(wv,wv),jac=v/(1.0-Dot(u,us)/v);
        double dnUnit=subArea/(L*L)*jac*exp(-beta*beta*w2)*dv*(flat ? wn*2.0*b4/Const::Pi : sqrt(w2)*b4/Const::Pi);

        //velocity relative to the spacecraft
        Vec3 vrel=u*v-Vsc;
        Vec3 urel=Normalized(vrel);
        double wng=gauges.NGWeight(urel,g),acc=gauges.RGAcceptance(urel,g);
        double gn=dnUnit*wng,gphi=(acc>0.0) ? dnUnit*acc*(-Dot(vrel,bRG)) : 0.0;

        //per unit flux at the reception time (kernels), and with the flux at the emission time (model values)
        double Ft=st.Flux[s];
        double ratio=(Ft>0.0) ? F[s]/Ft : 1.0;

        k.Gn[s]+=gn*ratio,k.GPhi[s]+=gphi*ratio;
        nF[iface*Const::nSpecies+s]+=F[s]*gn;
        phiF[iface*Const::nSpecies+s]+=F[s]*gphi;
        if (Dot(urel,Zsc)>0.0) nBelow[iface]+=F[s]*gn;
      }
    };

    for (int i=0;i<n;i++) for (int j=0;i+j<n;j++) {
      Sample(a+ea*(i+1.0/3.0)+eb*(j+1.0/3.0));
      if (i+j<n-1) Sample(a+ea*(i+2.0/3.0)+eb*(j+2.0/3.0));
    }
  }

  res.Altitude=*std::min_element(dist.begin(),dist.end());
  res.Kernels.clear(),res.KernelsExact.clear();

  double nTot=0.0,nB=0.0;

  for (int s=0;s<Const::nSpecies;s++) {
    double nn=0.0,phi=0.0,Q=0.0;

    for (int i=0;i<nFaces;i++) nn+=nF[i*Const::nSpecies+s],phi+=phiF[i*Const::nSpecies+s];
    for (int i=0;i<nFaces;i++) Q+=state[i].Flux[s]*mesh.Faces[i].Area;

    res.n[s]=nn;
    res.PhiIn[s]=phi;
    res.pNG[s]=gauges.NGPressure(s,nn);
    res.pRG[s]=gauges.RGPressure(s,phi,TemperatureRG);
    res.TotalProduction[s]=Q;
    nTot+=nn;
  }

  for (int i=0;i<nFaces;i++) {
    nB+=nBelow[i];

    bool nz=false;
    for (int s=0;s<Const::nSpecies;s++) if ((K[i].Gn[s]!=0.0)||(K[i].GPhi[s]!=0.0)) nz=true;
    if (nz) res.KernelsExact.push_back(K[i]);
  }

  res.NGBelowDeckFraction=(nTot>0.0) ? nB/nTot : 0.0;
  res.nVisibleElements=(int)res.KernelsExact.size();
}
