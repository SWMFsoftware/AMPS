//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//=======================================================================
//$Id$
//=======================================================================
//the definition and functions for the internal surfaces installed into the mesh


#ifndef _AMR_INTERNAL_SURFACE_SPHERE_
#define _AMR_INTERNAL_SURFACE_SPHERE_

#include "math.h"
#include <string>


#include "meshAMRdef.h"
#include "mpichannel.h"
#include "rnd.h"
#include "specfunc.h"


//=======================================================================
//the class describes the data that defines the spherical internal boundary; the class may contains the user defined data
class cInternalSphericalData : public cAMRexit
#if _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ == _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ON_
, public cInternalSphericalData_UserDefined
#endif
{
//protected:
public:
  double OriginPosition[3],Radius;

  static long int nZenithSurfaceElements,nAzimuthalSurfaceElements;
  static double dAzimuthalAngle;

  static double dCosZenithAngle;
  static double dZenithAngle;

  typedef void (*fPrintVariableList)(FILE*);
  fPrintVariableList PrintVariableList;

  typedef void (*fPrintTitle)(FILE*);
  fPrintTitle PrintTitle;

  char TitleMessage[400];

  typedef void (*fPrintDataStateVector)(FILE* fout,long int nZenithPoint,long int nAzimuthalPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalSphericalData *Sphere,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads);
  fPrintDataStateVector PrintDataStateVector;

  typedef double (*fLocalResolution)(double *);
  fLocalResolution localResolution;

  #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
  long int Temp_ID;
  #endif


  void cleanDataBuffer() {
    OriginPosition[0]=0.0,OriginPosition[1]=0.0,OriginPosition[2]=0.0;
    Radius=1.0;

    PrintVariableList=NULL,PrintDataStateVector=NULL,PrintTitle=NULL,localResolution=NULL;

    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    Temp_ID=0;
    #endif
  }

  cInternalSphericalData ()
#if _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ == _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ON_
  : cInternalSphericalData_UserDefined()
#endif
  {
    cleanDataBuffer();
    SetGeneralSurfaceMeshParameters(nZenithSurfaceElements,nAzimuthalSurfaceElements);

    //set up the default TitleMessage
    TitleMessage[0]=0;
  }

  static void SetGeneralSurfaceMeshParameters(long int nZenithElements,long int nAzimuthalElements) {
    nZenithSurfaceElements=nZenithElements,nAzimuthalSurfaceElements=nAzimuthalElements;

    dAzimuthalAngle=2.0*Pi/nAzimuthalSurfaceElements;

    switch (_INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_) {
    case _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_COSINE_DISTRIBUTION_: 
      dCosZenithAngle=2.0/nZenithSurfaceElements;
      break;
    case  _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_UNIFORM_DISTRIBUTION_: 
      dZenithAngle=Pi/nZenithSurfaceElements;
      break;
    default:
      ::exit(__LINE__,__FILE__,"Error: the option is not defiend");
    }
  }

  void SetSphereGeometricalParameters(double *x0,double r) {
     for (int idim=0;idim<3;idim++) OriginPosition[idim]=x0[idim];
     Radius=r;
  }

  void GetSphereGeometricalParameters(double* &x0,double &r) {
     x0=OriginPosition;
     r=Radius;
  }

  inline long int GetLocalSurfaceElementNumber(long int nZenithElement,long int nAzimuthalElement) {
    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    if ((nZenithElement<0)||(nZenithElement>=nZenithSurfaceElements)||(nAzimuthalElement<0)||(nAzimuthalElement>=nAzimuthalSurfaceElements)) {
      char msg[1000];

      sprintf(msg,"Error: 'nZenithElement' or 'nAzimuthalElement' are outside of the range\nnZenithElement=%lo\nnAzimuthalElement=%li\n",nZenithElement,nAzimuthalElement);
      exit(__LINE__,__FILE__,msg);
    }
    #endif

    return nZenithElement+nZenithSurfaceElements*nAzimuthalElement;
  }

  inline void GetSurfaceElementIndex(int &nZenithElement,int &nAzimuthalElement,int nSurfaceElement) {
    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    if ((nSurfaceElement<0)||(nSurfaceElement>=nZenithSurfaceElements*nAzimuthalSurfaceElements)) ::exit(__LINE__,__FILE__,"Error: 'nSurfaceElement' is out of range");
    #endif

    nAzimuthalElement=(int)(nSurfaceElement/nZenithSurfaceElements); 
    nZenithElement=nSurfaceElement-nZenithSurfaceElements*nAzimuthalElement; 
  }

  long int GetTotalSurfaceElementsNumber() {return nZenithSurfaceElements*nAzimuthalSurfaceElements;}


  double GetSurfaceElementArea(int nZenithElement,int nAzimuthalElement) {
    double cosZenithAngleMin,cosZenithAngleMax;

    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    if ((nZenithElement<0)||(nZenithElement>=nZenithSurfaceElements)||(nAzimuthalElement<0)||(nAzimuthalElement>=nAzimuthalSurfaceElements)) exit(__LINE__,__FILE__,"Error: 'nZenithElement' or 'nAzimuthalElement' are outside of the range ");
    #endif

    switch (_INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_) {
    case  _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_COSINE_DISTRIBUTION_: 
      cosZenithAngleMin=1.0-dCosZenithAngle*nZenithElement;
      cosZenithAngleMax=1.0-dCosZenithAngle*(1+nZenithElement);
      break;
    case  _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_UNIFORM_DISTRIBUTION_: 
      cosZenithAngleMin=cos(nZenithElement*dZenithAngle);
      cosZenithAngleMax=cos((1+nZenithElement)*dZenithAngle);
      break;
    default:
      exit(__LINE__,__FILE__,"Error: the option is not defined");
    }

    return (cosZenithAngleMin-cosZenithAngleMax)*dAzimuthalAngle*pow(Radius,2);
  }

  double GetSurfaceElementArea(int SurfaceElementNumber) {
    int nZenithElement,nAzimuthalElement;

    GetSurfaceElementIndex(nZenithElement,nAzimuthalElement,SurfaceElementNumber);
    return GetSurfaceElementArea(nZenithElement,nAzimuthalElement);
  } 


  inline void GetSurfaceElementNormal(double *norm,int nZenithElement,int nAzimuthalElement) {
    double ZenithAngle,AzimuthalAngle;

    AzimuthalAngle=(nAzimuthalElement+0.5)*dAzimuthalAngle;

    switch (_INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_) {
    case  _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_COSINE_DISTRIBUTION_: 
      ZenithAngle=acos(1.0-dCosZenithAngle*(nZenithElement+0.5));
      break;
    case _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_UNIFORM_DISTRIBUTION_: 
      ZenithAngle=(0.5+nZenithElement)*dZenithAngle; 
      break;
    default:
      exit(__LINE__,__FILE__,"Error: the option is not defined");
    }

    norm[0]=sin(ZenithAngle)*cos(AzimuthalAngle);
    norm[1]=sin(ZenithAngle)*sin(AzimuthalAngle);
    norm[2]=cos(ZenithAngle);
  }

  inline void GetSurfaceElementNormal(double *norm,int SurfaceElementNumber) {
    int nZenithElement,nAzimuthalElement;

    GetSurfaceElementIndex(nZenithElement,nAzimuthalElement,SurfaceElementNumber);
    GetSurfaceElementNormal(norm,nZenithElement,nAzimuthalElement);
  }



  inline void GetSurfaceElementMiddlePoint(double *x,int nZenithElement,int nAzimuthalElement) {
    int idim;

    GetSurfaceElementNormal(x,nZenithElement,nAzimuthalElement);
    for (idim=0;idim<3;idim++) x[idim]=x[idim]*Radius+OriginPosition[idim];
  }


  inline void GetSurfaceElementMiddlePoint(double *x,int SurfaceElementNumber) {
    int nZenithElement,nAzimuthalElement;

    GetSurfaceElementIndex(nZenithElement,nAzimuthalElement,SurfaceElementNumber);
    GetSurfaceElementMiddlePoint(x,nZenithElement,nAzimuthalElement);
  }

  inline void GetSurfaceElementRandomDirection(double *x,int nZenithElement,int nAzimuthalElement) {
    double ZenithAngle,AzimuthalAngle;
    double c0,c1;


    AzimuthalAngle=(nAzimuthalElement+rnd())*dAzimuthalAngle;

    switch (_INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_) {
    case  _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_COSINE_DISTRIBUTION_: 
      c0=1.0-dCosZenithAngle*nZenithElement;
      c1=1.0-dCosZenithAngle*(nZenithElement+1);
      break;
    case  _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_UNIFORM_DISTRIBUTION_: 
      c0=cos(nZenithElement*dZenithAngle);
      c1=cos((1+nZenithElement)*dZenithAngle);
      break;
    default:
      exit(__LINE__,__FILE__,"Error: the option is not defined");
    }

    ZenithAngle=acos(c0+rnd()*(c1-c0));

    x[0]=sin(ZenithAngle)*cos(AzimuthalAngle);
    x[1]=sin(ZenithAngle)*sin(AzimuthalAngle);
    x[2]=cos(ZenithAngle);
  }

  inline void GetSurfaceElementRandomPoint(double *x,int nZenithElement,int nAzimuthalElement) {
    GetSurfaceElementRandomDirection(x,nZenithElement,nAzimuthalElement);

    x[0]=Radius*x[0]+OriginPosition[0];
    x[1]=Radius*x[1]+OriginPosition[1];
    x[2]=Radius*x[2]+OriginPosition[2];
  }

  inline void GetSurfaceElementRandomPoint(double *x,int SurfaceElementNumber) {
    int nZenithElement,nAzimuthalElement;

    GetSurfaceElementIndex(nZenithElement,nAzimuthalElement,SurfaceElementNumber);
    GetSurfaceElementRandomPoint(x,nZenithElement,nAzimuthalElement);
  } 

  inline void GetSurfaceElementProjectionIndex(double *x,double &ZenithAngle,long int &nZenithElement,double &AzimuthalAngle,long int &nAzimuthalElement) {
    double r,r2,xNormalized[3];
    int idim;

    for (r2=0.0,idim=0;idim<3;idim++) r2+=pow(x[idim]-OriginPosition[idim],2);
    for (r=sqrt(r2),idim=0;idim<3;idim++) xNormalized[idim]=(x[idim]-OriginPosition[idim])/r;

    if ((r=pow(xNormalized[0],2)+pow(xNormalized[1],2))>0.0) {
      AzimuthalAngle=acos(xNormalized[0]/sqrt(r));
      if (xNormalized[1]<0.0) AzimuthalAngle=2.0*Pi-AzimuthalAngle;
    }
    else AzimuthalAngle=0.0;

    ZenithAngle=acos(xNormalized[2]);

    switch (_INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_) {
    case  _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_COSINE_DISTRIBUTION_: 
      nZenithElement=(long int)((1.0-xNormalized[2])/dCosZenithAngle);
      break;
    case _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_UNIFORM_DISTRIBUTION_: 
      nZenithElement=(long int)(ZenithAngle/dZenithAngle);
      break;
    default:
      exit(__LINE__,__FILE__,"Error: the option is not defined");
    }

    nAzimuthalElement=(long int)(AzimuthalAngle/dAzimuthalAngle);

    if (nZenithElement==nZenithSurfaceElements) --nZenithElement;
    if (nAzimuthalElement==nAzimuthalSurfaceElements) --nAzimuthalElement;

    #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
    static unsigned long int FunctionCallCounter=0;

    if ((nZenithElement<0)||(nZenithElement>=nZenithSurfaceElements)||(nAzimuthalElement<0)||(nAzimuthalElement>=nAzimuthalSurfaceElements)) {
      printf("$PREFIX:Error: out of range\n nZenithElement=%ld, nAzimuthalElement=%ld (%s@%i)\n",nZenithElement,nAzimuthalElement,__FILE__,__LINE__);
      printf("$PREFIX:x=%e, %e, %e\nCallCounter=%ld\n",x[0],x[1],x[2],FunctionCallCounter);
      exit(__LINE__,__FILE__,"Error: 'nZenithElement' or 'nAzimuthalElement' are outside of the range ");
    }

    FunctionCallCounter++;
    #endif
  }

  inline void GetSurfaceElementProjectionIndex(double *x,long int &nZenithElement,long int &nAzimuthalElement) {
    double ZenithAngle,AzimuthalAngle;
    GetSurfaceElementProjectionIndex(x,ZenithAngle,nZenithElement,AzimuthalAngle,nAzimuthalElement);
  }

  inline void etSurfaceElementProjectionAngle(double *x,double &ZenithAngle,double &AzimuthalAngle) {
    long int nZenithElement,nAzimuthalElement;
    GetSurfaceElementProjectionIndex(x,ZenithAngle,nZenithElement,AzimuthalAngle,nAzimuthalElement);
  }

  inline void GetSurfaceNormal(double *x,double iZenithPoint,double  iAzimutalPoint) {
    double ZenithAngle,cosZenithAngle,sinZenithAngle,AzimuthalAngle;

    switch (_INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_MODE_) {
    case _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_COSINE_DISTRIBUTION_: 
      cosZenithAngle=1.0-dCosZenithAngle*iZenithPoint;

      if (cosZenithAngle<-1.0) cosZenithAngle=-1.0;
      if (cosZenithAngle>1.0) cosZenithAngle=1.0; 

      sinZenithAngle=sqrt(1.0-cosZenithAngle*cosZenithAngle);
      break;
    case _INTERNAL_BOUNDARY_SPHERE_ZENITH_ANGLE_UNIFORM_DISTRIBUTION_: 
      ZenithAngle=dZenithAngle*iZenithPoint;

      cosZenithAngle=cos(ZenithAngle);
      sinZenithAngle=sin(ZenithAngle);
      break;
    default:
      exit(__LINE__,__FILE__,"Error: the option is not defined");
    }

    AzimuthalAngle=dAzimuthalAngle*iAzimutalPoint;

    x[0]=sinZenithAngle*cos(AzimuthalAngle);
    x[1]=sinZenithAngle*sin(AzimuthalAngle);
    x[2]=cosZenithAngle;
  }

  inline void GetSurfaceCoordinate(double *x,double iZenithPoint,double  iAzimutalPoint) {
    GetSurfaceNormal(x,iZenithPoint,iAzimutalPoint);

    x[0]=Radius*x[0]+OriginPosition[0];
    x[1]=Radius*x[1]+OriginPosition[1];
    x[2]=Radius*x[2]+OriginPosition[2];
  }

  inline void GetSurfaceLonLatNormal(double &lon,double &lat,double iZenithPoint,double  iAzimutalPoint) {
    double x[3],r,t;
    int idim;

    GetSurfaceCoordinate(x,iZenithPoint,iAzimutalPoint);

    for (idim=0;idim<3;idim++) x[idim]-=OriginPosition[idim];
    r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);

    t=x[2]/r;

    if (t>=1.0) lat=90.0;
    else if (t<=-1.0) lat=-90.0;
    else lat=180.0/Pi*asin(t);

    if (isfinite(lat)==false) {
      char msg[1000];

      sprintf(msg,"Error: something wrong with the surface coordinate(x=%e, %e, %e, r=%e, t=%e, lat=%e, iZenithPoint=%e,iAzimutalPoint=%e\n",
             x[0],x[1],x[2],r,t,lat,iZenithPoint,iAzimutalPoint);
      exit(__LINE__,__FILE__,msg);
    }

    lon=180.0/Pi*dAzimuthalAngle*iAzimutalPoint;
  }

  void PrintSurfaceData(const char *fname,int nDataSet, bool PrintStateVectorFlag=true) {
    long int iZenith,iAzimuthal;
    FILE *fout=NULL,*fout2d=NULL;
    double x[3];

    CMPI_channel pipe(1000000);
    int ThisThread=0,nTotalThreads=1;
    MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&ThisThread);
    MPI_Comm_size(MPI_GLOBAL_COMMUNICATOR,&nTotalThreads);

    if (ThisThread==0) {
      // The caller controls the base path length, so construct the companion
      // filename dynamically rather than appending into a 300-byte buffer.
      const std::string fname2d=std::string(fname)+".2d.dat";

      fout=fopen(fname,"w");
      fout2d=fopen(fname2d.c_str(),"w");
      pipe.openRecvAll();

      //print the output file title
      if (PrintTitle!=NULL) {
        PrintTitle(fout);
        fprintf(fout,"\n");
      }
      else {
        if (TitleMessage[0]!=0) fprintf(fout,"%s\n",TitleMessage);
      }

      //print the variable list
      fprintf(fout,"VARIABLES=\"X\", \"Y\", \"Z\"");
      fprintf(fout2d,"VARIABLES=\"Lon\", \"Lat\"");


      if (PrintStateVectorFlag==true) {
        if (PrintVariableList==NULL) exit(__LINE__,__FILE__,"Error: PrintVariableList is not defined");
        PrintVariableList(fout);
        PrintVariableList(fout2d);
      }

      //print the number of variables and blocks
      fprintf(fout,"\nZONE N=%ld, E=%ld, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n",(nZenithSurfaceElements+1)*nAzimuthalSurfaceElements,nZenithSurfaceElements*nAzimuthalSurfaceElements);
      fprintf(fout2d,"\nZONE I=%ld, J=%ld, DATAPACKING=POINT\n",nAzimuthalSurfaceElements,nZenithSurfaceElements+1);
    }
    else pipe.openSend(0);

    //interpolate and print the state vector
    long int InterpolationList[nAzimuthalSurfaceElements],InterpolationListLength=0;

    for (iZenith=0;iZenith<nZenithSurfaceElements+1;iZenith++) for (iAzimuthal=0;iAzimuthal<nAzimuthalSurfaceElements;iAzimuthal++) {
      GetSurfaceCoordinate(x,iZenith,iAzimuthal);

      if (ThisThread==0) {
        fprintf(fout,"%e %e %e ",x[0],x[1],x[2]);

        double lon,lat;
        GetSurfaceLonLatNormal(lon,lat,iZenith,iAzimuthal);
        fprintf(fout2d,"%e %e ",lon,lat);
      }

      if (PrintStateVectorFlag==true) {
        if (PrintDataStateVector==NULL) exit(__LINE__,__FILE__,"Error: PrintDataStateVector is not defined");

        //prepare the interpolation stencil
        InterpolationListLength=0;

        if (iZenith==0) {
          InterpolationList[InterpolationListLength++]=GetLocalSurfaceElementNumber(0,iAzimuthal);
          InterpolationList[InterpolationListLength++]=GetLocalSurfaceElementNumber(0,((iAzimuthal>0) ? iAzimuthal-1 : nAzimuthalSurfaceElements-1));
        }
        else if (iZenith==nZenithSurfaceElements) {
          InterpolationList[InterpolationListLength++]=GetLocalSurfaceElementNumber(nZenithSurfaceElements-1,iAzimuthal);
          InterpolationList[InterpolationListLength++]=GetLocalSurfaceElementNumber(nZenithSurfaceElements-1,((iAzimuthal>0) ? iAzimuthal-1 : nAzimuthalSurfaceElements-1));
        }
        else {
          int iA,iZ,A[2],Z[2];

          Z[0]=iZenith-1,Z[1]=iZenith;

          A[0]=(iAzimuthal!=0) ? iAzimuthal-1 : nAzimuthalSurfaceElements-1;
          A[1]=iAzimuthal;

          for (iA=0;iA<2;iA++) for (iZ=0;iZ<2;iZ++) InterpolationList[InterpolationListLength++]=GetLocalSurfaceElementNumber(Z[iZ],A[iA]);
        }


        PrintDataStateVector(fout,iZenith,iAzimuthal,InterpolationList,InterpolationListLength,this,nDataSet,&pipe,ThisThread,nTotalThreads);
        PrintDataStateVector(fout2d,iZenith,iAzimuthal,InterpolationList,InterpolationListLength,this,nDataSet,&pipe,ThisThread,nTotalThreads);
      }

      if (ThisThread==0) {
        fprintf(fout,"\n");
        fprintf(fout2d,"\n");
      }
    }

    //close the pipe
    if (ThisThread==0) pipe.closeRecvAll();
    else pipe.closeSend();

    //print the connectivity list
    long int iAzimuthalMax,iAzimuthalMin;
    long int nd0,nd1,nd2,nd3;

    if (ThisThread==0) {
      for (iZenith=0;iZenith<nZenithSurfaceElements;iZenith++) for (iAzimuthal=0;iAzimuthal<nAzimuthalSurfaceElements;iAzimuthal++) {
        iAzimuthalMax=(iAzimuthal+1!=nAzimuthalSurfaceElements) ? iAzimuthal+1 : 0;
        iAzimuthalMin=iAzimuthal;

        nd0=1+iAzimuthalMin+iZenith*nAzimuthalSurfaceElements;
        nd1=1+iAzimuthalMax+iZenith*nAzimuthalSurfaceElements;
        nd2=1+iAzimuthalMax+(iZenith+1)*nAzimuthalSurfaceElements;
        nd3=1+iAzimuthalMin+(iZenith+1)*nAzimuthalSurfaceElements;

        fprintf(fout,"%ld %ld %ld %ld\n",nd0,nd1,nd2,nd3);
      }

      fclose(fout);
      fclose(fout2d);
    }
  }


   void PrintSurfaceMesh(const char *fname) {PrintSurfaceData(fname,0,false);}

   //=======================================================================
   //intersection of a block with the sphere
   int BlockIntersection(double *xBlockMin,double *xBlockMax,double EPS) {
     int idim,i,j,k,nCounter;
     double x[3];

     //if all corners of the block are within the sphere -> the block is entirely within the sphere
     for (nCounter=0,i=0;i<2;i++) {
       x[0]=((i==0) ? xBlockMin[0] : xBlockMax[0])-OriginPosition[0];

       for (j=0;j<2;j++) {
         x[1]=((j==0) ? xBlockMin[1] : xBlockMax[1])-OriginPosition[1];

         for (k=0;k<2;k++) {
           x[2]=((k==0) ? xBlockMin[2] : xBlockMax[2])-OriginPosition[2];


           if (x[0]*x[0]+x[1]*x[1]+x[2]*x[2]<=pow(Radius+EPS,2)) nCounter++;
           else goto Not_Inside_Sphere;

         }
       }
     }


Not_Inside_Sphere:
     if (nCounter==8) return _AMR_BLOCK_OUTSIDE_DOMAIN_;

     /* The squared distance from a point to an axis-aligned box is the sum
      * of squared coordinate gaps to the nearest interval endpoint.  The
      * sphere intersects (or is wholly contained in) the box exactly when
      * this distance is <= R^2.  Unlike face/edge-only searches, this also
      * detects a small sphere completely enclosed by a large block. */
     double MinimumDistance2=0.0;
     for (idim=0;idim<3;idim++) {
       double gap=0.0;
       if (OriginPosition[idim]<xBlockMin[idim])
         gap=xBlockMin[idim]-OriginPosition[idim];
       else if (OriginPosition[idim]>xBlockMax[idim])
         gap=OriginPosition[idim]-xBlockMax[idim];
       MinimumDistance2+=gap*gap;
     }

     return (MinimumDistance2<=pow(Radius+EPS,2)) ?
         _AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_ : _AMR_BLOCK_INSIDE_DOMAIN_;
   }


 private:

   inline double _angle_aux_2(double x, double y) {
     // auxilary angle first arising in volume cut by 2 planes
     // valid for the unit sphere
     double x2 = x * x, y2 = y * y;
     double misc = 1.0 - x2 * y2 - x2 - y2;

     if(misc == 0.0) return 0.5 * Pi;
     double corr = (misc < 0.0) ? Pi : 0.0;
     return corr + atan(2.0 * x * y * sqrt(1.0 - x2 - y2) / misc );
   }

   inline double _angle_aux_3(double x, double y, double z) {
     // auxilary angle first arising in volume cut by 3 planes
     // valid for the unit sphere
     double x2 = x * x, y2 = y * y, z2 = z * z;
     double sqrtxz = sqrt(1.0 - x2 - z2), sqrtyz = sqrt(1.0 - y2 - z2);
     double misc = x * y -  sqrtxz * sqrtyz;

     if(misc == 0.0) return 0.5 * Pi;
     return atan( (y * sqrtxz + x * sqrtyz) / misc);
   }

   inline double _volume_3(double x, double y, double z) {
     // volume cut by 3 planes, x >= 0, y >= 0, z >= 0
     // valid for the unit sphere
     double x2 = x * x, y2 = y * y, z2 = z * z;
     double sqrtxy = sqrt(1.0 - x2 - y2);
     double sqrtxz = sqrt(1.0 - x2 - z2);
     double sqrtyz = sqrt(1.0 - z2 - y2);

     return (x * y * sqrtxy + y * z * sqrtyz + x * z * sqrtxz) / 3. - x * y * z + (Pi -
         x*(3. - x2)*(0.5*Pi + _angle_aux_3(y,z,x)) -
         y*(3. - y2)*(0.5*Pi + _angle_aux_3(x,z,y)) -
         z*(3. - z2)*(0.5*Pi + _angle_aux_3(x,y,z)) -
         _angle_aux_2(x,y) - _angle_aux_2(x,z) - _angle_aux_2(y,z)) / 6.;
   }
   

 public:

   double GetRemainedBlockVolume(double *xBlockMinInit,double *xBlockMaxInit,double EPS,int *IntersectionStatus) {
     int BlockIntersectionCode;
     double res=0.0;
     int idim;

     /*
      * Return the part of an axis-aligned block that remains in the AMPS
      * computational domain, i.e. the block volume outside the excluded
      * sphere.  All input coordinates, Radius, EPS, and the returned volume
      * use one consistent caller-selected length unit (the result is unit^3).
      * For an intersected block we translate by OriginPosition and normalize
      * every coordinate and width by R.  Inclusion/exclusion of _volume_3()
      * then supplies the overlap with the unit sphere and multiplication by
      * R^3 restores physical volume:
      *
      *   V_remaining = R^3 (V_normalized_box - V_box_intersect_unit_sphere).
      *
      * The closed-form octant primitive requires a positive radius.  A
      * non-positive radius is invalid geometry and is rejected through the
      * normal AMPS fatal path rather than fabricating a volume.
      */
     if (IntersectionStatus==NULL) {
       exit(__LINE__,__FILE__,"Error: sphere/block intersection status pointer is null");
       return 0.0;
     }

     if (!(Radius>0.0) || !isfinite(Radius)) {
       exit(__LINE__,__FILE__,"Error: sphere radius must be finite and positive");
       *IntersectionStatus=_AMR_BLOCK_INSIDE_DOMAIN_;
       return 0.0;
     }

     BlockIntersectionCode=BlockIntersection(xBlockMinInit,xBlockMaxInit,EPS);
     *IntersectionStatus=BlockIntersectionCode;

     if (BlockIntersectionCode==_AMR_BLOCK_OUTSIDE_DOMAIN_) {
       // Every block corner lies within the excluded sphere.
       return 0.0;
     }

     if (BlockIntersectionCode==_AMR_BLOCK_INSIDE_DOMAIN_) {
       // The sphere does not intersect the block: its entire volume remains.
       for (res=1.0,idim=0;idim<3;idim++) res*=xBlockMaxInit[idim]-xBlockMinInit[idim];
       return res;
     }


     // Split at a sphere-centered coordinate plane.  This maps every child to
     // one octant, where the analytic primitive is valid, and is translation
     // invariant for a sphere whose center is not the global origin.
     for(idim = 0; idim < 3; idim++) {
       if((xBlockMinInit[idim]-OriginPosition[idim]) *
          (xBlockMaxInit[idim]-OriginPosition[idim]) < 0.0) {
         // Preserve the two unsplit extents.  A single "middle" vector used
         // as both a maximum and a minimum collapses those dimensions to zero
         // and loses seven octants when a box encloses the sphere.  The left
         // child's maximum begins as the original maximum; the right child's
         // minimum begins as the original minimum.  Only the selected axis is
         // replaced by the sphere-centered plane.
         double xLeftMax[3],xRightMin[3];

         memcpy(xLeftMax,xBlockMaxInit,3*sizeof(double));
         memcpy(xRightMin,xBlockMinInit,3*sizeof(double));
         xLeftMax[idim]=OriginPosition[idim];
         xRightMin[idim]=OriginPosition[idim];
         res += GetRemainedBlockVolume(xBlockMinInit,xLeftMax,EPS,IntersectionStatus);
         res += GetRemainedBlockVolume(xRightMin,xBlockMaxInit,EPS,IntersectionStatus);

         // The status describes the caller's original block, not whichever
         // recursive child happened to be evaluated last.
         *IntersectionStatus=_AMR_BLOCK_INTERSECT_DOMAIN_BOUNDARY_;
         return res;
       }
     }

     // reflect coordinates if needed (for convenience)
     // also normalize both positions and widths to the sphere radius.
     double x[3], d[3], R3 = Radius*Radius*Radius;

     for(idim = 0; idim < 3; idim++) {
       d[idim]=(xBlockMaxInit[idim]-xBlockMinInit[idim])/Radius;

       const double xmin=xBlockMinInit[idim]-OriginPosition[idim];
       const double xmax=xBlockMaxInit[idim]-OriginPosition[idim];
       if ((xmin < 0.) && (xmax <= 0.)) x[idim]=-xmax/Radius;
       else x[idim]=xmin/Radius;
     }

     // now the cell lays fully in the first octant
     // and x are the coordinates of the vertex closest to the origin

     // check the trivial cases: cube whithin or outside the sphere
     if(x[0]*x[0] + x[1]*x[1] + x[2]*x[2] >= 1.0) return R3 * d[0] * d[1] *d[2];        // outside

     if ((x[0]+d[0])*(x[0]+d[0]) + (x[1]+d[1])*(x[1]+d[1]) + (x[2]+d[2])*(x[2]+d[2]) < 1.0) return 0.0;// within

     // add volumes cut by three plains with weight -1, 0 or +1
     for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++) for (int k = 0; k < 2; k++) {
       if((x[0]+i*d[0])*(x[0]+i*d[0]) + (x[1]+j*d[1])*(x[1]+j*d[1]) + (x[2]+k*d[2])*(x[2]+k*d[2]) < 1.0) {
         res +=(((i+j+k)%2==0)?1.:-1.)*_volume_3(x[0]+i*d[0],x[1]+j*d[1],x[2]+k*d[2]);
       }
     }

     return R3 * (d[0]*d[1]*d[2]-res);
   }

   double GetRemainedBlockVolumeNumerical(double *xBlockMinInit,double *xBlockMaxInit,double EPS,int *IntersectionStatus) {
     int BlockIntersectionCode;
     double res=0.0,TotalResult=0.0;
     int idim;

     static const int nLevelMax=6;

     struct cLevelData {
       double xSubBlockMin[3],xSubBlockMax[3],dx,dy;
       int i,j;
     };

     cLevelData LevelData[nLevelMax+1];
     cLevelData *levelDataPtr;
     int nLevel=0;

     //init the level 0 data
     for (levelDataPtr=LevelData,idim=0;idim<3;idim++) levelDataPtr->xSubBlockMax[idim]=xBlockMaxInit[idim],levelDataPtr->xSubBlockMin[idim]=xBlockMinInit[idim];

FunctionBegins:
     res=0.0;
     levelDataPtr=LevelData+nLevel;

     BlockIntersectionCode=BlockIntersection(levelDataPtr->xSubBlockMin,levelDataPtr->xSubBlockMax,EPS);
     *IntersectionStatus=BlockIntersectionCode;

     if (BlockIntersectionCode==_AMR_BLOCK_OUTSIDE_DOMAIN_) {
       if (nLevel==0) return 0.0;
       else {
         res=0.0;
         goto LevelProcessingDone;
       }
     }

     if (BlockIntersectionCode==_AMR_BLOCK_INSIDE_DOMAIN_) {
       for (res=1.0,idim=0;idim<3;idim++) res*=levelDataPtr->xSubBlockMax[idim]-levelDataPtr->xSubBlockMin[idim];

       #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
       if (res<0.0) exit(__LINE__,__FILE__,"Error: out of range");
       #endif

       if (nLevel==0) return res;
       else {
         goto LevelProcessingDone;
       }
     }

     //the block is intersected by the sphere
     if (nLevel<nLevelMax) {
       levelDataPtr->dx=(levelDataPtr->xSubBlockMax[0]-levelDataPtr->xSubBlockMin[0])/2.0;
       levelDataPtr->dy=(levelDataPtr->xSubBlockMax[1]-levelDataPtr->xSubBlockMin[1])/2.0;

       (levelDataPtr+1)->xSubBlockMin[2]=levelDataPtr->xSubBlockMin[2];
       (levelDataPtr+1)->xSubBlockMax[2]=levelDataPtr->xSubBlockMax[2];

       for (levelDataPtr->i=0;levelDataPtr->i<2;levelDataPtr->i++) {
         (levelDataPtr+1)->xSubBlockMin[0]=levelDataPtr->xSubBlockMin[0]+levelDataPtr->i*levelDataPtr->dx;
         (levelDataPtr+1)->xSubBlockMax[0]=levelDataPtr->xSubBlockMin[0]+(levelDataPtr->i+1)*levelDataPtr->dx;

         for (levelDataPtr->j=0;levelDataPtr->j<2;levelDataPtr->j++) {
           (levelDataPtr+1)->xSubBlockMin[1]=levelDataPtr->xSubBlockMin[1]+levelDataPtr->j*levelDataPtr->dy;
           (levelDataPtr+1)->xSubBlockMax[1]=levelDataPtr->xSubBlockMin[1]+(levelDataPtr->j+1)*levelDataPtr->dy;

           nLevel++;
           res=0.0;
           goto FunctionBegins;

//           res+=GetRemainedBlockVolume(xSubBlockMin,xSubBlockMax,EPS,&BlockIntersectionCode,nLevel+1);

LevelProcessingDone:
           TotalResult+=res;
           res=0.0;
           nLevel--;
           levelDataPtr=LevelData+nLevel;
         }
       }
     }
     else {
       double SegmentSplittingTime[4],t1,t2,x0[3],l,xy2,R2;
       double a,b,c,d,sqrt_d;
       int i,nSegments=1;

       double dii=0.5*(levelDataPtr->xSubBlockMax[0]-levelDataPtr->xSubBlockMin[0]);
       double djj=0.5*(levelDataPtr->xSubBlockMax[1]-levelDataPtr->xSubBlockMin[1]);
       int ii,jj;
       bool IntersectionFound=false;

       l=levelDataPtr->xSubBlockMax[2]-levelDataPtr->xSubBlockMin[2];

       R2=Radius*Radius;

       for (ii=0;ii<3;ii++) for (jj=0;jj<3;jj++) {
         x0[0]=levelDataPtr->xSubBlockMin[0]+ii*dii;
         x0[1]=levelDataPtr->xSubBlockMin[1]+jj*djj;
         x0[2]=levelDataPtr->xSubBlockMin[2];


         a=l*l;
         b=2.0*l*(x0[2]-OriginPosition[2]);
         for (c=-R2,idim=0;idim<3;idim++) c+=pow(x0[idim]-OriginPosition[idim],2);

         d=b*b-4.0*a*c;

         if (d<=0.0) {
//         if (pow(x0[0]-OriginPosition[0],2)+pow(x0[1]-OriginPosition[1],2)+pow(x0[2]-OriginPosition[2]+l*0.5,2)>R2) res+=1;

//         res=(levelDataPtr->xSubBlockMax[0]-levelDataPtr->xSubBlockMin[0])*(levelDataPtr->xSubBlockMax[1]-levelDataPtr->xSubBlockMin[1])*(levelDataPtr->xSubBlockMax[2]-levelDataPtr->xSubBlockMin[2]);
         }
         else {
           if (b<0.0) a*=-1.0,b*=-1.0,c*=-1.0;

           sqrt_d=sqrt(d);
           t1=-(b+sqrt_d)/(2.0*a);
           t2=-2.0*c/(b+sqrt_d);
           nSegments=0;

           SegmentSplittingTime[0]=0.0;
           if ((t1>0.0)&&(t1<1.0)) SegmentSplittingTime[++nSegments]=t1;

           if ((t2>0.0)&&(t2<1.0)) {
             if (t2>SegmentSplittingTime[nSegments]) SegmentSplittingTime[++nSegments]=t2;
             else {
               SegmentSplittingTime[2]=SegmentSplittingTime[1];
               SegmentSplittingTime[1]=t2;
               nSegments++;
             }
           }

           SegmentSplittingTime[++nSegments]=1.0;
           xy2=pow(x0[0]-OriginPosition[0],2)+pow(x0[1]-OriginPosition[1],2);

           for (i=0;i<nSegments;i++) {
             t1=SegmentSplittingTime[i],t2=SegmentSplittingTime[i+1];

             if (xy2+pow(x0[2]-OriginPosition[2]+l*0.5*(t1+t2),2)>R2) {
               res+=t2-t1;
               IntersectionFound=true;
             }
           }
         }
       }

       res/=9.0;
       if (IntersectionFound==false) res=pow(EPS,3.0);
       else res*=(levelDataPtr->xSubBlockMax[0]-levelDataPtr->xSubBlockMin[0])*(levelDataPtr->xSubBlockMax[1]-levelDataPtr->xSubBlockMin[1])*(levelDataPtr->xSubBlockMax[2]-levelDataPtr->xSubBlockMin[2]);

       #if _AMR_DEBUGGER_MODE_ == _AMR_DEBUGGER_MODE_ON_
       if (res<0.0) {
         exit(__LINE__,__FILE__,"Error: out of range");
       }
       #endif
     }

     if (nLevel!=0) goto LevelProcessingDone;

     if (TotalResult<0.0) exit(__LINE__,__FILE__,"Error: out of range");

     return TotalResult;
   }
};

//set up the surface parameters of the sphere
void inline SetGeneralSphereSurfaceMeshParameters(long int nZenithElements,long int nAzimuthalElements) {
  cInternalSphericalData::nAzimuthalSurfaceElements=nAzimuthalElements,cInternalSphericalData::nZenithSurfaceElements=nZenithElements;
  cInternalSphericalData::dAzimuthalAngle=2.0*Pi/cInternalSphericalData::nAzimuthalSurfaceElements;

  cInternalSphericalData::dCosZenithAngle=2.0/cInternalSphericalData::nZenithSurfaceElements;
  cInternalSphericalData::dZenithAngle=Pi/cInternalSphericalData::nZenithSurfaceElements;
}


#endif
