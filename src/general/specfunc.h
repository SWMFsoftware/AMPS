//$Id$
//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
#ifndef SPECFUNC
#define SPECFUNC

#include "mpi.h"

#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "rnd.h"
#include "constants.h"

#define _STDOUT_ERRORLOG_MODE__ON_   0
#define _STDOUT_ERRORLOG_MODE__OFF_  1
#define _STDOUT_ERRORLOG_MODE_ _STDOUT_ERRORLOG_MODE__ON_

extern int ThisThread;
extern int TotalThreadsNumber;

/* gamma and error functions used from g++ math library
double erf(double);
double gam(double);
*/

long int nint(double);

void PrintErrorLog(const char*);
void PrintErrorLog(long int,const char*,const char*);

void StampSignature(char*);
extern int ExitErrorCode;
void exit(long int,const char*,const char* =NULL);
void PrintLineMark(long int,char*,char* =NULL);

template<class T>
void PrintLineMark(long int nline ,char* fname ,T code,const char *msg=NULL) {
  long thread;
  char *buffer=new char[sizeof(T)*TotalThreadsNumber];

#ifdef MPI_ON
  MPI_Gather((char*)&code,sizeof(T),MPI_CHAR,buffer,sizeof(T),MPI_CHAR,0,MPI_GLOBAL_COMMUNICATOR);
#else
  char *ptr;
  int i;  

  for (ptr=(char*)&code,i=0;i<sizeof(T);i++,ptr++) buffer[i]=*ptr;  
#endif

  if (ThisThread==0) {
    std::cout << "linemark: line=" << nline << ", file=" << fname;
    if (msg!=NULL) std::cout << ", msg=" << msg;
    std::cout << ": code=";

    for (thread=0;thread<TotalThreadsNumber;thread++) std::cout << "  " << ((T*)buffer)[thread]; 

    std::cout << std::endl;
  }

  delete [] buffer;
}


#ifdef DIM
template<class TMesh>
bool GetGradient(double* gradQ,double cellQ,double* Q,long int ncell,TMesh &grid) {
  int counter,pface;
  long int neib;
  double dx2,x[4][3],x0[3],df[4];
  double A[3][3],aa[3][3],af[3],detaa;

  switch (DIM) {
  case 1:
    grid.cell[ncell].GetCellCenter(x0);
    counter=0,gradQ[0]=0.0,dx2=0.0;

    for (pface=0;pface<DIM+1;pface++) if ((neib=grid.cell[ncell].neighbour_cellno[pface])>=0) {
      counter++;

      grid.cell[neib].GetCellCenter(x[pface]);
      gradQ[0]+=(Q[pface]-cellQ)*(x[pface][0]-x0[0]);
      dx2=pow(x[pface][0]-x0[0],2);
    }

    gradQ[0]/=(counter!=0) ? dx2 : 1.0; 

    if (counter!=DIM+1) return false;
    break;
  case 2:
    grid.cell[ncell].GetCellCenter(x0);
    counter=0,gradQ[0]=0.0,gradQ[1]=0.0;

    for (pface=0;pface<DIM+1;pface++) if ((neib=grid.cell[ncell].neighbour_cellno[pface])>=0) {
      counter++;

      grid.cell[neib].GetCellCenter(x[pface]);
      A[pface][0]=x[pface][0]-x0[0],A[pface][1]=x[pface][1]-x0[1];
      df[pface]=Q[pface]-cellQ;
    }

    if (counter!=DIM+1) return false;

    aa[0][0]=pow(A[0][0],2)+pow(A[1][0],2)+pow(A[2][0],2);
    aa[0][1]=A[0][0]*A[0][1]+A[1][0]*A[1][1]+A[2][0]*A[2][1];
    aa[1][0]=aa[0][1];
    aa[1][1]=pow(A[0][1],2)+pow(A[1][1],2)+pow(A[2][1],2);

    af[0]=A[0][0]*df[0]+A[1][0]*df[1]+A[2][0]*df[2];
    af[1]=A[0][1]*df[0]+A[1][1]*df[1]+A[2][1]*df[2];

    detaa=aa[0][0]*aa[1][1]-aa[1][0]*aa[0][1];

    gradQ[0]=(af[0]*aa[1][1]-af[1]*aa[0][1])/detaa;
    gradQ[1]=(aa[0][0]*af[1]-aa[1][0]*af[0])/detaa;

    break;
  default:
    printf("$PREFIX:Error: GetGradient. DIM=%i is not implemented\n",DIM);
    exit(__LINE__,__FILE__);
  }

  return true;
}
#endif

//=========================================================
//calculation of CRC-32
class CRC32 {
private:
  unsigned long crc_accum,crc_table[256];

  //generate the table of CRC remainders for all possible bytes 
  void generare_crc_table() { 
    int i, j;
    unsigned long crc_accum;

    for (i=0;i<256;i++) { 
      crc_accum=((unsigned long)i<<24);
      for (j=0;j<8;j++) crc_accum=(crc_accum&0x80000000L) ? (crc_accum<<1)^0x04c11db7L : crc_accum<<1;
      crc_table[i]=crc_accum; 
    }
  };

public: 

  CRC32 () {
    crc_accum=0;
    generare_crc_table();
  };

  //update the CRC on the data block one byte at a time
  template <class T> 
  void add(T* buffer, long int size) {
    char *data_blk_ptr=(char*)buffer;
    long int data_blk_size=size*sizeof(T); 
    long int i,j;

    for (j=0;j<data_blk_size;j++) { 
      i=((int)(crc_accum>>24)^ *data_blk_ptr++)&0xff;
      crc_accum=(crc_accum<<8)^crc_table[i]; 
    }
  } 

  template <class T>
  void add(T t) {
    add(&t,1);
  } 

  void clear() {
    crc_accum=0;
  }

  unsigned long checksum() { 
    return crc_accum;
  }

  void PrintChecksum(long int nline,const char* fname) {
    char message[1000];
    
    sprintf(message," line=%ld, file=%s",nline,fname);
    PrintChecksum(message);
  }

  void PrintChecksumThread(long int nline,const char* fname,int ThisThread=-1) {
    char message[1000];

    sprintf(message," line=%ld, file=%s",nline,fname);
    if (ThisThread!=-1)  sprintf(message,"%s, thread=%i",message,ThisThread);

    printf("$PREFIX:CRC32 checksum=0x%lx, message=%s\n",checksum(),message);
  }

  void PrintChecksum(const char* message=NULL) {
    unsigned long int *buffer=new unsigned long int[TotalThreadsNumber];
    long int thread;

    buffer[0]=checksum();

#ifdef MPI_ON
    unsigned long int bufferRecv[TotalThreadsNumber];

    MPI_Gather(buffer,1,MPI_UNSIGNED_LONG,bufferRecv,1,MPI_UNSIGNED_LONG,0,MPI_GLOBAL_COMMUNICATOR);
    memcpy(buffer,bufferRecv,TotalThreadsNumber*sizeof(unsigned long int));
#endif

    if (ThisThread==0) {
      CRC32 cumulativeSignature;
      cumulativeSignature.add(buffer,TotalThreadsNumber);

      if (message!=NULL) printf("$PREFIX:CRC32 checksum, cumulativeSignature=0x%lx, message=%s:\n",cumulativeSignature.checksum(),message);
      else printf("$PREFIX:CRC32 checksum, cumulativeSignature=0x%lx:\n",cumulativeSignature.checksum());

      for (thread=0;thread<TotalThreadsNumber;thread++) printf("$PREFIX:thread=%ld, sum=0x%lx\n",thread,buffer[thread]);
    }

    delete [] buffer;
  }

  void PrintChecksumSingleThread(const char* message=NULL) {
    int ThisThread=0;

    #ifdef MPI_ON
    MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&ThisThread);
    #endif

    if (message!=NULL) printf("$PREFIX:CRC32 checksum=0x%lx, message=%s (thread=%i):\n",checksum(),message,ThisThread);
    else printf("$PREFIX:CRC32 checksum=0x%lx  (thread=%i):\n",checksum(),ThisThread);
  }
};


//=========================================================
//Vector Rotations
namespace VectorRotation {
  inline void Along_Z_direction(double *l,double angle) {
    double temp[2],cosAngle,sinAngle;

    cosAngle=cos(angle);
    sinAngle=sin(angle);

    temp[0]=cosAngle*l[0]-sinAngle*l[1];
    temp[1]=sinAngle*l[0]+cosAngle*l[1];

    memcpy(l,temp,2*sizeof(double));
  }
}


//=========================================================
//Vector Operations
namespace Vector3D {
  inline void CrossProduct(double *res,double *a,double *b) {
    res[0]=a[1]*b[2]-a[2]*b[1];
    res[1]=a[2]*b[0]-a[0]*b[2];
    res[2]=a[0]*b[1]-a[1]*b[0];
  }

  inline double Length(double *x) {
    return sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  }

  inline void Copy(double *target,double *source,int length=3) {
    memcpy(target,source,length*sizeof(double));
  }

  inline double DotProduct(double *a,double *b) {
    int i;
    double res=0.0;

    for (i=0;i<3;i++) res+=a[i]*b[i];
    return res;
  }

  inline void Orthogonalize(double *PrimaryVector,double *OrthogonalVector) {
     double l2=0.0,c=0.0;
     int i;

     //get the dot products of the vectors 
     for (i=0;i<3;i++) {
        l2+=PrimaryVector[i]*PrimaryVector[i];
        c+=PrimaryVector[i]*OrthogonalVector[i];
     }

     //get the orthogonal component of 'OrthogonalVector'
     c/=l2;

     for (i=0;i<3;i++) {
       OrthogonalVector[i]-=c*PrimaryVector[i];
      }
   }

  inline double ParallelComponentLength(double *Vector,double *Axis) {
    int i;
    double l=0.0,c=0.0;

    for (i=0;i<3;i++) {
      l+=pow(Axis[i],2);
      c+=Vector[i]*Axis[i];
    }

    return c/sqrt(l);
  }
 

  //determine an orthogonal frame of rederence: z is input; e1 and e2 and orthogonal to z and form a right handed frame of reference
  inline void GetNormFrame(double *e0,double *e1,double *z) {
    double l,e[3];
    int idim;

    l=Length(z);
    for (idim=0;idim<3;idim++) e[idim]=z[idim]/l;

    //get e0
    if (fabs(e[0])>1.0E-5) e0[0]=-e[1],e0[1]=e[0],e0[2]=0.0;
    else e0[0]=0.0,e0[1]=e[2],e0[2]=-e[1];

    l=Length(e0);
    for (idim=0;idim<3;idim++) e0[idim]/=l;

    //get e1: e1=z \times e0
    CrossProduct(e1,z,e0);
  }

  inline void GetRandomNormFrame(double *e0,double *e1,double *z) {
    double e0temp[3],e1temp[3],theta,sin_theta,cos_theta;

    GetNormFrame(e0temp,e1temp,z);
    theta=rnd()*PiTimes2;
    sin_theta=sin(theta);
    cos_theta=cos(theta);

    for (int idim=0;idim<3;idim++) {
      e0[idim]=cos_theta*e0temp[idim]+sin_theta*e1temp[idim];
      e1[idim]=-sin_theta*e0temp[idim]+cos_theta*e1temp[idim];
    }
  }

  inline void Normalize(double *x,double NewLength=1.0) {
    double l=NewLength/Length(x);

    for (int idim=0;idim<3;idim++) x[idim]*=l;
  }

  //distribute the vector direction
  namespace Distribution {
    //uniform distribution of the
    inline void Uniform(double *a,double length=1.0) {
      for (int i=0;i<3;i++) a[i]=sqrt(-log(rnd()))*cos(PiTimes2*rnd());
      Vector3D::Normalize(a,length);
    }

    namespace Circle {
       inline void Uniform(double *a, double Radius) {
         double phi,r;

         phi=PiTimes2*rnd();
         r=sqrt(rnd())*Radius;

         a[0]=r*sin(phi);
         a[1]=r*cos(phi);
       }

       inline void Uniform(double *x, double *e0,double *e1,double Radius) {
         double xPlane_e0_e1[2];

         //location of the generated point in the e0-e1 plane
         Distribution::Circle::Uniform(xPlane_e0_e1,Radius);

         //a global coordinate of the point
         for (int idim=0;idim<3;idim++) x[idim]=xPlane_e0_e1[0]*e0[idim]+xPlane_e0_e1[1]*e1[idim];
       }

      inline void Uniform(double *x, double *e0,double *e1,double *xCicleCenterLocation,double Radius) {
        double xPlane_e0_e1[2];

        //location of the generated point in the e0-e1 plane
        Distribution::Circle::Uniform(xPlane_e0_e1,Radius);

        //a global coordinate of the point
        for (int idim=0;idim<3;idim++) x[idim]=xCicleCenterLocation[idim]+xPlane_e0_e1[0]*e0[idim]+xPlane_e0_e1[1]*e1[idim];
      }

      inline void Gaussian(double *x, double *e0,double *e1,double *xCicleCenterLocation,double Radius,double ProbabilityAtRadius=1.0E-10) {
        double A,phi,r,xPlane_e0_e1[2];

        //distribute the radius of the point f=exp(-A*t^2)
        A=-log(ProbabilityAtRadius)/sqrt(Radius);

        do {
          r=Radius*rnd();
        }
        while (exp(-A*r*r)<rnd());

        //distribute the angle
        phi=PiTimes2*rnd();

        //calculate the vector
        xPlane_e0_e1[0]=r*sin(phi);
        xPlane_e0_e1[1]=r*cos(phi);

        for (int idim=0;idim<3;idim++) x[idim]=xCicleCenterLocation[idim]+xPlane_e0_e1[0]*e0[idim]+xPlane_e0_e1[1]*e1[idim];
      }
    }

  }
}
  

//=========================================================
//Relativistic functions
namespace Relativistic {
  inline double Speed2E(double Speed,double mass) {
    double beta2,gamma2;

    //relativistic kinetic energy: KE=(gamma-1) \times m_0 \times c^2
    beta2=pow(Speed/SpeedOfLight,2);
    gamma2=1.0/(1.0-beta2);
    return mass*SpeedOfLight*SpeedOfLight*beta2*gamma2/(1.0+sqrt(gamma2)); //(gamma-1)=(gamma^2-1)/(gamma+1)
  }

  inline double E2Speed(double E,double mass) {
    double mc2;

    mc2=mass*SpeedOfLight*SpeedOfLight;
    return SpeedOfLight*sqrt(E*(E+2.0*mc2))/(E+mc2);
  }

  inline double GetGamma(double *v) {
    return 1.0/sqrt(1.0-(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])/(SpeedOfLight*SpeedOfLight));
  }

  inline double Speed2Momentum(double Speed,double mass) {
    return mass*Speed/sqrt(1.0-(Speed*Speed)/(SpeedOfLight*SpeedOfLight));
  }

  inline double Momentum2Speed(double Momentum,double mass) {
    return Momentum*SpeedOfLight/sqrt(Momentum*Momentum+pow(mass*SpeedOfLight,2));
  }

  inline double Momentum2Speed(double* MomentumVector,double mass) {
    return Momentum2Speed(Vector3D::Length(MomentumVector),mass);
  }

  inline double Momentum2Energy(double Momentum,double mass) {
    double t=mass*SpeedOfLight;

    return t*SpeedOfLight*sqrt(1.0+pow(Momentum/t,2));
  }

  //get gyro frequency
  inline double GetGyroFrequency(double *v,double ParticleRestMass, double ElectricCharge, double* B) {
    return fabs(ElectricCharge)*Vector3D::Length(B)/(PiTimes2*ParticleRestMass*GetGamma(v));
  }

  inline double GetGyroRadius(double *v,double ParticleRestMass,double ElectricCharge,double *B) {
    double c,VelocityPerp=0.0,absB,absB2=0.0;
    int idim;

    for (idim=0;idim<3;idim++) absB2+=B[idim]*B[idim];

    c=Vector3D::DotProduct(v,B)/absB2;
    absB=sqrt(absB2);

    for (idim=0;idim<3;idim++) VelocityPerp+=pow(v[idim]-c*B[idim],2);

    VelocityPerp=sqrt(VelocityPerp);
    return ParticleRestMass*GetGamma(v)*VelocityPerp/(fabs(ElectricCharge)*absB);
  }

  namespace GuidingCeneter {
    inline double GetGamma(double v_parallel,double p_perp,double mass) {
      return sqrt((1.0+pow(p_perp/(mass*SpeedOfLight),2))/(1.0-pow(v_parallel/SpeedOfLight,2)));
    }

    inline double GetEnergy(double v_parallel,double p_perp,double mass) {
      double gamma,t2;

      gamma=GetGamma(v_parallel,p_perp,mass);
      t2=pow(v_parallel/SpeedOfLight,2);

      return (t2+pow(p_perp/(mass*SpeedOfLight),2))/((1.0-t2)*(gamma+1));
    }
  }

}


//===========================================================================================
//functions that can be used for the code debugging
 namespace Debugger {
   //save data into debugger stream
   void SaveDataIntoStream(void* data,int length,const char* msg);

   template <class T>
   void SaveDataIntoStream(T data,const char* msg);
 }


#endif
   
