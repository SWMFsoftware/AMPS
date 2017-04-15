//$Id$
//the set of function for processing of the cut-cells and cut-faces of the internal irregular surfaces

#include "pic.h"

int PIC::Mesh::IrregularSurface::nCutFaceInformationCopyAttempts=0; 

void PIC::Mesh::IrregularSurface::InitExternalNormalVector() {
  //calculate external normals to the faces
  double xStart[3],xFinish[3],l,l0,*norm;
  long int nface;
  cTriangleFace *fcptr;
  int idim,iIntersections;

  const double angleCosMin=cos(85.0/180.0*3.141592654);

  int nStartFace,nFinishFace,nTotalThreads,ThisThread,nFaceThread;

  //init the ray tracing module if needed
  PIC::RayTracing::Init();

  //check whether external normal vectors are already have been determined for the surface trianguletion
  unsigned long int TriangulationSignature;
  char fname[_MAX_STRING_LENGTH_PIC_];
  FILE *fExternalVectorFile=NULL;

  TriangulationSignature=CutCell::GetTriangulationSignature();
  sprintf(fname,"amr.sig=0x%lx.TriangulationExternalNormalVector.bin",TriangulationSignature);

  fExternalVectorFile=fopen(fname,"r");

  if (fExternalVectorFile!=NULL) {
    //the binary file containing the external normal information exists -> load it
    if (PIC::ThisThread==0) {
      printf("$PREFIX: Binary file with the external normals (%s) is found: loading.... ",fname);
      fflush(stdout);
    }

    for (nface=0;nface<nBoundaryTriangleFaces;nface++) fread(BoundaryTriangleFaces[nface].ExternalNormal,sizeof(double),3,fExternalVectorFile);

    fclose(fExternalVectorFile);
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    if (PIC::ThisThread==0) {
      printf("done\n");
      fflush(stdout);
    }

    return;
  }

  if (PIC::ThisThread==0) {
    printf("$PREFIX: Binary file with the external normals is not found: generating.... ");
    fflush(stdout);
  }

  MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&ThisThread);
  MPI_Comm_size(MPI_GLOBAL_COMMUNICATOR,&nTotalThreads);

  nFaceThread=nBoundaryTriangleFaces/nTotalThreads;
  nStartFace=nFaceThread*ThisThread;
  nFinishFace=nStartFace+nFaceThread;
  if (ThisThread==nTotalThreads-1) nFinishFace=nBoundaryTriangleFaces;

  //evaluate the distance size of the domain
  double lDomain=2.0*sqrt(pow(PIC::Mesh::mesh.xGlobalMax[0]-PIC::Mesh::mesh.xGlobalMin[0],2)+
      pow(PIC::Mesh::mesh.xGlobalMax[1]-PIC::Mesh::mesh.xGlobalMin[1],2)+
      pow(PIC::Mesh::mesh.xGlobalMax[2]-PIC::Mesh::mesh.xGlobalMin[2],2));

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#pragma omp parallel for schedule(dynamic,1) default (none) shared(nStartFace,nFinishFace,BoundaryTriangleFaces,lDomain) \
  private (nface,fcptr,norm,idim,xFinish,xStart,l,l0,iIntersections)
#endif
  for (nface=nStartFace;nface<nFinishFace;nface++) {
    fcptr=BoundaryTriangleFaces+nface;

    //get the center point and the normal of the face
    fcptr->GetCenterPosition(xStart);
    norm=fcptr->ExternalNormal;

    if (Vector3D::Length(norm)<0.1) {
      if (fcptr->SurfaceArea>1.0E-15) exit(__LINE__,__FILE__,"Error: something is wrong with the face normal. The length of the normal is zero, while the surface area of the face is not zero");
      continue;
    }

    do {
      for (idim=0,l=0.0,l0=0.0;idim<3;idim++) {
        xFinish[idim]=sqrt(-2.0*log(rnd()))*cos(2.0*3.1415926*rnd());
        l+=pow(xFinish[idim],2);
        l0+=xFinish[idim]*norm[idim];
      }

      l=sqrt(l);
    }
    while (fabs(l0)/l<angleCosMin);

    if (l0<0.0) l*=-1.0;
    for (idim=0;idim<3;idim++) xFinish[idim]=lDomain*xFinish[idim]/l+xStart[idim];

    //count face intersections
    iIntersections=PIC::RayTracing::CountFaceIntersectionNumber(xStart,xFinish,fcptr->MeshFileID,fcptr);

    if (iIntersections%2!=0) {
      //the norm has to be reversed
      for (idim=0;idim<3;idim++) fcptr->ExternalNormal[idim]*=-1.0;
    }
  }

  //collect the surface normals
  double *sendBuffer=new double[3*2*nFaceThread];
  int thread,cnt;

  for (thread=0;thread<nTotalThreads;thread++) {
    nStartFace=nFaceThread*thread;
    nFinishFace=nStartFace+nFaceThread;
    if (thread==nTotalThreads-1) nFinishFace=nBoundaryTriangleFaces;

    if (thread==ThisThread) {
      for (nface=nStartFace,cnt=0;nface<nFinishFace;nface++,cnt++) memcpy(sendBuffer+3*cnt,BoundaryTriangleFaces[nface].ExternalNormal,3*sizeof(double));
    }

    MPI_Bcast(sendBuffer,3*(nFinishFace-nStartFace),MPI_DOUBLE,thread,MPI_GLOBAL_COMMUNICATOR);

    if (thread!=ThisThread) {
      for (nface=nStartFace,cnt=0;nface<nFinishFace;nface++,cnt++) memcpy(BoundaryTriangleFaces[nface].ExternalNormal,sendBuffer+3*cnt,3*sizeof(double));
    }
  }

  delete [] sendBuffer;

  //save a file with the external normals for the future use
  if (PIC::ThisThread==0) {
    fExternalVectorFile=fopen(fname,"w");

    for (nface=0;nface<nBoundaryTriangleFaces;nface++) fwrite(BoundaryTriangleFaces[nface].ExternalNormal,sizeof(double),3,fExternalVectorFile);

    fclose(fExternalVectorFile);
  }


  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread==0) {
    printf("done\n");
    fflush(stdout);
  }
}

//check weather a point (x0) in insed the domain:
//if the number if interasections of the ray (x=x0+l*t) is even than the point is within the domain; otherwise the point is outsede the domain
//l -> is a random ray (intersection search) direction
bool PIC::Mesh::IrregularSurface::CheckPointInsideDomain_default(double *x,PIC::Mesh::IrregularSurface::cTriangleFace* SurfaceTriangulation,int nSurfaceTriangulation,bool ParallelCheck,double EPS) {
  int nface,nfaceStart,nfaceFinish,iIntersections;
  double l,xFinish[3];
  int idim;
  bool flag=true;

  if (SurfaceTriangulation==NULL) return true;

  //evaluate the distance size of the domain
  double lDomain=2.0*sqrt(pow(PIC::Mesh::mesh.xGlobalMax[0]-PIC::Mesh::mesh.xGlobalMin[0],2)+
      pow(PIC::Mesh::mesh.xGlobalMax[1]-PIC::Mesh::mesh.xGlobalMin[1],2)+
      pow(PIC::Mesh::mesh.xGlobalMax[2]-PIC::Mesh::mesh.xGlobalMin[2],2));

  //distribute ditrction of the search
  for (l=0.0,idim=0;idim<3;idim++) {
    xFinish[idim]=sqrt(-2.0*log(rnd()))*cos(2.0*3.1415926*rnd());
    l+=pow(xFinish[idim],2);
  }

  for (l=sqrt(l),idim=0;idim<3;idim++) xFinish[idim]=lDomain*xFinish[idim]/l+x[idim];

  //xFinish is outside of the domain -> the point outside of the surface
  //calculate the number of the face intersections between 'x' and 'xFinish'

  iIntersections=PIC::RayTracing::CountFaceIntersectionNumber(x,xFinish,-1);

  return (iIntersections%2==0) ? true : false;

 /*   SearchDirection[idim]/=l;


  static bool initflag=false;
  static int ThisThread,nTotalThreads;

  if (initflag==false) {
    MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&ThisThread);
    MPI_Comm_size(MPI_GLOBAL_COMMUNICATOR,&nTotalThreads);
    initflag=true;
  }

  if (ParallelCheck==true) {
    nfaceStart=ThisThread*(nSurfaceTriangulation/nTotalThreads);
    nfaceFinish=(ThisThread+1)*(nSurfaceTriangulation/nTotalThreads);
    if (ThisThread==nTotalThreads-1) nfaceFinish=nSurfaceTriangulation;
  }
  else nfaceStart=0,nfaceFinish=nSurfaceTriangulation;

  bool flagbuffer[nTotalThreads];

  do {
    //distribute ditrction of the search
    for (l=0.0,idim=0;idim<3;idim++) {
      SearchDirection[idim]=sqrt(-2.0*log(rnd()))*cos(2.0*3.1415926*rnd());
      l+=pow(SearchDirection[idim],2);
    }

    if (ParallelCheck==true) MPI_Bcast(SearchDirection,3,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);

    for (l=sqrt(l),idim=0;idim<3;idim++) SearchDirection[idim]/=l;
    iIntersections=0;
    flag=true;

    //find intersections with the faces on the mesh
    for (nface=nfaceStart;nface<nfaceFinish;nface++) {
      if (SurfaceTriangulation[nface].RayIntersection(x,SearchDirection,EPS)==true) iIntersections++;

      for (l=0.0,idim=0;idim<3;idim++) l+=pow(SurfaceTriangulation[nface].ExternalNormal[idim]*SearchDirection[idim],2);
      if (l<1.0E-10) {
        flag=false;
        break;
      }

    }

    if (ParallelCheck==true) {
      MPI_Gather(&flag,sizeof(bool),MPI_CHAR,flagbuffer,sizeof(bool),MPI_CHAR,0,MPI_GLOBAL_COMMUNICATOR);
      if (ThisThread==0) for (int thread=1;thread<nTotalThreads;thread++) if (flagbuffer[thread]==false) flag=false;
      MPI_Bcast(&flag,sizeof(bool),MPI_CHAR,0,MPI_GLOBAL_COMMUNICATOR);
    }
  }
  while (flag==false);

  if (ParallelCheck==true) {
    int t;

    MPI_Allreduce(&iIntersections,&t,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
    iIntersections=t;
  }


  return (2*(iIntersections/2)==iIntersections) ? true : false;
*/}


//copy information about the triangular cut faces into the nighbouring blocks 
void PIC::Mesh::IrregularSurface::CopyCutFaceInformation(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  static unsigned int *BoundaryTriangleMap=NULL;
  static const int nResetBoundaryTriangleMapTestCounter=100000000;
  static unsigned int BoundaryTriangleMapTestCounter=0;

  static int CutFaceDescriptorTablePointer=-1;
  static CutCell::cTriangleFaceDescriptor *CutFaceDescriptorTable=NULL;
  static int CutFaceDescriptorTableLength=(int)(0.5*CutCell::nBoundaryTriangleFaces);     

  struct CutFaceData {
    void ResetTempPointerTable(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *nd) {
      if (nd->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) {
        int i;
        cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;

        for (i=0;i<(1<<DIM);i++) if ((downNode=nd->downNode[i])!=NULL) {
          ResetTempPointerTable(downNode);
        }
      }
      else nd->neibFirstTriangleCutFace_temp=NULL;
    } 

    void SetTempCutFacePointers(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *nd) {
      if (nd->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) {
        int i;
        cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;

        for (i=0;i<(1<<DIM);i++) if ((downNode=nd->downNode[i])!=NULL) {
          SetTempCutFacePointers(downNode);
        }
      }
      else if (nd->neibFirstTriangleCutFace_temp!=NULL) {
        //add the pointers to the neibFirstTriangleCutFace list
        CutCell::cTriangleFaceDescriptor *t;

        if (nd->neibFirstTriangleCutFace==NULL) {
          nd->neibFirstTriangleCutFace=nd->neibFirstTriangleCutFace_temp;
          nd->neibFirstTriangleCutFace_temp=NULL;
        }
        else {
          t=nd->neibFirstTriangleCutFace;
          while (t->next!=NULL) t=t->next; 

          //link the lists
          t->next=nd->neibFirstTriangleCutFace_temp;
          nd->neibFirstTriangleCutFace_temp->prev=t;

          nd->neibFirstTriangleCutFace_temp=NULL;
        }
      }
    }

    //add cut faces from block NodeFrom to NodeTo 
    void AddCutFaceData(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *To,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *From) { 
      int iTestNeib,iTestNode;
      CutCell::cTriangleFaceDescriptor *t;

      //save all cut-faces from the neib node
      if (++BoundaryTriangleMapTestCounter==nResetBoundaryTriangleMapTestCounter) {
        //reset the counters
        BoundaryTriangleMapTestCounter=1;
        for (int i=0;i<CutCell::nBoundaryTriangleFaces;i++) BoundaryTriangleMap[i]=0;
      }

      for (iTestNode=0;iTestNode<3;iTestNode++) {
        switch (iTestNode) {
        case 0: 
          t=To->FirstTriangleCutFace;
          break;
        case 1: 
          t=To->neibFirstTriangleCutFace;
          break;
        case 2:
          t=To->neibFirstTriangleCutFace_temp;
          break;
        default:
          exit(__LINE__,__FILE__,"Error: something is wrong");
        }

        for (;t!=NULL;t=t->next) BoundaryTriangleMap[t->TriangleFace->Temp_ID]=BoundaryTriangleMapTestCounter;
      }

      //check whether all cut-faces from startNode are regiteres in neibNode
      for (iTestNode=0;iTestNode<2;iTestNode++) {
        t=(iTestNode==0) ? From->FirstTriangleCutFace : From->neibFirstTriangleCutFace;

        for (;t!=NULL;t=t->next) if (BoundaryTriangleMap[t->TriangleFace->Temp_ID]!=BoundaryTriangleMapTestCounter) {
          CutCell::cTriangleFaceDescriptor *NewDescriptor;

          if (CutFaceDescriptorTablePointer==-1) {
            CutFaceDescriptorTable=new CutCell::cTriangleFaceDescriptor[CutFaceDescriptorTableLength];
            CutFaceDescriptorTablePointer=CutFaceDescriptorTableLength-1;
          }

          NewDescriptor=CutFaceDescriptorTable+CutFaceDescriptorTablePointer;
          CutFaceDescriptorTablePointer--;

          NewDescriptor->TriangleFace=t->TriangleFace;
          NewDescriptor->prev=NULL;
          NewDescriptor->next=To->neibFirstTriangleCutFace_temp;

          if (To->neibFirstTriangleCutFace_temp!=NULL) To->neibFirstTriangleCutFace_temp->prev=NewDescriptor;
          To->neibFirstTriangleCutFace_temp=NewDescriptor;

          BoundaryTriangleMap[t->TriangleFace->Temp_ID]=BoundaryTriangleMapTestCounter;
        }
      }
    }

 
  } ProcessCutFaceData;
     

  if (startNode==PIC::Mesh::mesh.rootTree) {
    BoundaryTriangleMapTestCounter=1;
    BoundaryTriangleMap=new unsigned int [CutCell::nBoundaryTriangleFaces];  
         
    for (int i=0;i<CutCell::nBoundaryTriangleFaces;i++) {
      BoundaryTriangleMap[i]=0,CutCell::BoundaryTriangleFaces[i].Temp_ID=i;
    }

    //allocate the buffer with the new descriptors if needed 
    if (CutFaceDescriptorTable==NULL) {
      CutFaceDescriptorTable=new CutCell::cTriangleFaceDescriptor[CutFaceDescriptorTableLength];
      CutFaceDescriptorTablePointer=CutFaceDescriptorTableLength-1;
    }


    ProcessCutFaceData.ResetTempPointerTable(PIC::Mesh::mesh.rootTree); 
  }

  //move to the botton of the tree
  if (startNode->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) {
    int i;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;
    double c;

    for (i=0;i<(1<<DIM);i++) if ((downNode=startNode->downNode[i])!=NULL) {
      CopyCutFaceInformation(downNode); 
    }
  }
  else if ((startNode->FirstTriangleCutFace!=NULL)||(startNode->neibFirstTriangleCutFace!=NULL)) {
    const int ProcessedNeibBlockTableLength=60;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* ProcessedNeibBlockTable[ProcessedNeibBlockTableLength];
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* neibNode;
    int i,j,k,iFace,iEdge,iCorner,iNeib,ProcessedNodeCounter=0;
    bool found;

    //connection through corners
    for (iCorner=0;iCorner<(1<<DIM);iCorner++) if ((neibNode=startNode->GetNeibCorner(iCorner))!=NULL) {
      for (found=false,iNeib=0;iNeib<ProcessedNodeCounter;iNeib++) if (neibNode==ProcessedNeibBlockTable[iNeib]) {
        found=true;
        break;
      }

      if (found==false) {
        //the block is not processed yet -> copy cut face information 
        ProcessedNeibBlockTable[ProcessedNodeCounter++]=neibNode;
        ProcessCutFaceData.AddCutFaceData(neibNode,startNode);  
      }
    }
  
    //connection through edges 
    for (iEdge=0;iEdge<12;iEdge++) for (i=0;i<2;i++) if ((neibNode=startNode->GetNeibEdge(iEdge,i))!=NULL) {
      for (found=false,iNeib=0;iNeib<ProcessedNodeCounter;iNeib++) if (neibNode==ProcessedNeibBlockTable[iNeib]) {
        found=true;
        break;
      }

      if (found==false) {
        //the block is not processed yet -> copy cut face information
        ProcessedNeibBlockTable[ProcessedNodeCounter++]=neibNode;
        ProcessCutFaceData.AddCutFaceData(neibNode,startNode);
      }
    }

    //connection through faces   
    for (iFace=0;iFace<6;iFace++) for (i=0;i<2;i++) for (j=0;j<2;j++) if ((neibNode=startNode->GetNeibFace(iFace,i,j))!=NULL) {
      for (found=false,iNeib=0;iNeib<ProcessedNodeCounter;iNeib++) if (neibNode==ProcessedNeibBlockTable[iNeib]) {
        found=true;
        break;
      }

      if (found==false) {
        //the block is not processed yet -> copy cut face information
        ProcessedNeibBlockTable[ProcessedNodeCounter++]=neibNode;
        ProcessCutFaceData.AddCutFaceData(neibNode,startNode);
      }
    }

    if (ProcessedNodeCounter>ProcessedNeibBlockTableLength) exit(__LINE__,__FILE__,"Error: the counting is wrong"); 
  } 

 





/*
    int iNeib,jNeib,kNeib;
    CutCell::cTriangleFaceDescriptor *t,*tNeib;

    //scan through all neibours
    int iNeibNode,jNeibNode,kNeibNode;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *neibNode;

    static const int iNeibNodeMin=-1,iNeibNodeMax=1,jNeibNodeMin=-1,jNeibNodeMax=1,kNeibNodeMin=-1,kNeibNodeMax=1;

    for (iNeibNode=iNeibNodeMin;iNeibNode<=iNeibNodeMax;iNeibNode++) { 
      for (jNeibNode=jNeibNodeMin;jNeibNode<=jNeibNodeMax;jNeibNode++) {
        for (kNeibNode=kNeibNodeMin;kNeibNode<=kNeibNodeMax;kNeibNode++)  {
          if (abs(iNeibNode)+abs(jNeibNode)+abs(kNeibNode)!=0) {
            if ((neibNode=PIC::Mesh::mesh.getNeibNode(iNeibNode,jNeibNode,kNeibNode,startNode))!=NULL) {
              int iTestNeib,iTestNode;

              //save all cut-faces from the neib node
              if (++BoundaryTriangleMapTestCounter==nResetBoundaryTriangleMapTestCounter) {
                //reset the counters 
                BoundaryTriangleMapTestCounter=1;
                for (int i=0;i<CutCell::nBoundaryTriangleFaces;i++) BoundaryTriangleMap[i]=0; 
              }

              for (iTestNode=0;iTestNode<2;iTestNode++) {
                t=(iTestNode==0) ? neibNode->FirstTriangleCutFace : neibNode->neibFirstTriangleCutFace;

                for (;t!=NULL;t=t->next) BoundaryTriangleMap[t->TriangleFace->Temp_ID]=BoundaryTriangleMapTestCounter;
              }     

              //check whether all cut-faces from startNode are regiteres in neibNode
              for (iTestNode=0;iTestNode<2;iTestNode++) {
                t=(iTestNode==0) ? startNode->FirstTriangleCutFace : startNode->neibFirstTriangleCutFace; 

                for (;t!=NULL;t=t->next) if (BoundaryTriangleMap[t->TriangleFace->Temp_ID]!=BoundaryTriangleMapTestCounter) {
                  CutCell::cTriangleFaceDescriptor *NewDescriptor;

                  if (CutFaceDescriptorTablePointer==-1) {
                    CutFaceDescriptorTable=new CutCell::cTriangleFaceDescriptor[CutFaceDescriptorTableLength];
                    CutFaceDescriptorTablePointer=CutFaceDescriptorTableLength-1;  
                  }

                  NewDescriptor=CutFaceDescriptorTable+CutFaceDescriptorTablePointer;
                  CutFaceDescriptorTablePointer--;

                  NewDescriptor->TriangleFace=t->TriangleFace;
                  NewDescriptor->prev=NULL;
                  NewDescriptor->next=neibNode->neibFirstTriangleCutFace;

                  if (neibNode->neibFirstTriangleCutFace!=NULL) neibNode->neibFirstTriangleCutFace->prev=NewDescriptor;
                  neibNode->neibFirstTriangleCutFace=NewDescriptor;   

                  BoundaryTriangleMap[t->TriangleFace->Temp_ID]=BoundaryTriangleMapTestCounter; 
                }
              }

            }
          }
        }
      }
    } 
  }
*/

  if (startNode==PIC::Mesh::mesh.rootTree) {
    delete [] BoundaryTriangleMap;

    BoundaryTriangleMapTestCounter=0;
    BoundaryTriangleMap=NULL;

    ProcessCutFaceData.SetTempCutFacePointers(PIC::Mesh::mesh.rootTree);
  }
}
 

//========================================================================================================================
double PIC::Mesh::IrregularSurface::GetClosestDistance(double *x) {
  double xClosestPoint[3];
  int iClosestTriangularFace;

  return GetClosestDistance(x,xClosestPoint,iClosestTriangularFace);
}

double PIC::Mesh::IrregularSurface::GetClosestDistance(double *x,double *xClosestPoint,int& iClosestTriangularFace) {
  double xFace[3],c,*ExternNormal,Altitude=-1.0,l[3],xIntersection[3],xIntersectionLocal[3],IntersectionTime,t;
  int iFace,idim,iPoint;
  int nThreadsOpenMP=1;

  //determine whether the point is insde the triangulated surface
  if (CutCell::CheckPointInsideDomain(x,CutCell::BoundaryTriangleFaces,CutCell::nBoundaryTriangleFaces,false,0.0)==false) return -1.0;

  //loop through the cut-faces to detemine the closest distance to the surface
#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
   #pragma omp parallel shared(nThreadsOpenMP)
   {
     #pragma omp single
     {
       nThreadsOpenMP=omp_get_num_threads();
     }
   }

   double AltitudeTable[nThreadsOpenMP];
   double **xClosestPointTable=new double* [nThreadsOpenMP];
   int iClosestTriangularFaceTable[nThreadsOpenMP];


   xClosestPointTable[0]=new double [3*nThreadsOpenMP];
   for (int i=1;i<nThreadsOpenMP;i++) xClosestPointTable[i]=xClosestPointTable[0]+3*i;

   #pragma omp parallel default(none) private (iPoint,xIntersectionLocal,xIntersection,IntersectionTime,xFace,iFace,ExternNormal,c,idim,t,l) shared (x,PIC::Mesh::mesh,iClosestTriangularFaceTable,xClosestPointTable,AltitudeTable,CutCell::nBoundaryTriangleFaces,CutCell::BoundaryTriangleFaces,nThreadsOpenMP)
   {
   int iThreadOpenMP=omp_get_thread_num();

   AltitudeTable[iThreadOpenMP]=-1.0;
   for (iFace=0;iFace<CutCell::nBoundaryTriangleFaces;iFace++) if (iFace/nThreadsOpenMP==iThreadOpenMP) {
#else
   int iThreadOpenMP=0;
   double AltitudeTable[1]={-1.0};
   double **xClosestPointTable=new double* [1];
   int iClosestTriangularFaceTable[1];

   xClosestPointTable[0]=new double[3];

   for (iFace=0;iFace<CutCell::nBoundaryTriangleFaces;iFace++) {
#endif

    //the external point has to be pointed in the direction of the point of test
    ExternNormal=CutCell::BoundaryTriangleFaces[iFace].ExternalNormal;

    CutCell::BoundaryTriangleFaces[iFace].GetCenterPosition(xFace);
    for (c=0.0,idim=0;idim<DIM;idim++) c+=(l[idim]=(x[idim]-xFace[idim]))*ExternNormal[idim];

    if (c<0.0) continue;

    //the extermal normal of the face is derected toward the tested point ==>
    //evaluate the distance to the face at the center, corners of the face, and in the direction normal to the surface
    t=Vector3D::Length(l);
    if ((t<AltitudeTable[iThreadOpenMP])||(AltitudeTable[iThreadOpenMP]<0.0)) {
      AltitudeTable[iThreadOpenMP]=t;
      memcpy(xClosestPointTable[iThreadOpenMP],xFace,3*sizeof(double));
      iClosestTriangularFaceTable[iThreadOpenMP]=iFace;
    }

    //the closest point to the face is that along the normal check intersecion of the line along the normal with the surface element
    for (idim=0;idim<DIM;idim++) l[idim]=-ExternNormal[idim];

    if (CutCell::BoundaryTriangleFaces[iFace].RayIntersection(x,l,IntersectionTime,xIntersectionLocal,xIntersection,PIC::Mesh::mesh.EPS)==true) {
      for (c=0.0,idim=0;idim<DIM;idim++) c+=pow(x[idim]-xIntersection[idim],2);

      t=sqrt(c);

      if (t<AltitudeTable[iThreadOpenMP]) {
        AltitudeTable[iThreadOpenMP]=t;
        memcpy(xClosestPointTable[iThreadOpenMP],xIntersection,3*sizeof(double));
        iClosestTriangularFaceTable[iThreadOpenMP]=iFace;
      }
    }
    else {
      //in case there is no intersection of the line that is along the  normal -> get the distances to the corners to the face
      for (iPoint=0;iPoint<DIM;iPoint++) {
        //get the corner point on the face
        switch (iPoint) {
        case 0:
          memcpy(xFace,CutCell::BoundaryTriangleFaces[iFace].x0Face,3*sizeof(double));
          break;
        case 1:
          memcpy(xFace,CutCell::BoundaryTriangleFaces[iFace].x1Face,3*sizeof(double));
          break;
        case 2:
          memcpy(xFace,CutCell::BoundaryTriangleFaces[iFace].x2Face,3*sizeof(double));
          break;
        default:
          exit(__LINE__,__FILE__,"Error: something went wrong");
        }

        for (c=0.0,idim=0;idim<DIM;idim++) c+=pow(x[idim]-xFace[idim],2);

        t=sqrt(c);

        if (t<AltitudeTable[iThreadOpenMP]) {
          AltitudeTable[iThreadOpenMP]=t;
          memcpy(xClosestPointTable[iThreadOpenMP],xFace,3*sizeof(double));
          iClosestTriangularFaceTable[iThreadOpenMP]=iFace;
        }
      }
    }
  }
#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  }
#endif

  //collect altitude information from all OpenMP threads
  Altitude=AltitudeTable[0];
  memcpy(xClosestPoint,xClosestPointTable[0],3*sizeof(double));
  iClosestTriangularFace=iClosestTriangularFaceTable[0];

  for (int thread=1;thread<nThreadsOpenMP;thread++) if (Altitude=AltitudeTable[thread]) {
    Altitude=AltitudeTable[thread];
    memcpy(xClosestPoint,xClosestPointTable[thread],3*sizeof(double));
    iClosestTriangularFace=iClosestTriangularFaceTable[thread];
  }

  delete [] xClosestPointTable[0];
  delete [] xClosestPointTable;

  return Altitude;
}



