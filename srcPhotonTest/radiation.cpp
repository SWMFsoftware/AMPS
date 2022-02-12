#include "radiation.h"


//particle weight == energy carried by a model particle [eV]


long int Radiation::PhotonFreqOffset=-1;
int Radiation::MaterialTemperatureOffset=-1;

void Radiation::Init() {
  //request a place in a particle's tate vector
  PIC::ParticleBuffer::RequestDataStorage(PhotonFreqOffset,sizeof(double));
  PIC::IndividualModelSampling::RequestStaticCellData.push_back(RequestStaticCellData);

  //register output functions
  PIC::Mesh::AddVaraibleListFunction(PrintVariableList);
  PIC::Mesh::PrintDataCenterNode.push_back(PrintData);
  PIC::Mesh::InterpolateCenterNode.push_back(Interpolate);

  //particle interaction with the boundary of the domain
  PIC::Mover::ProcessOutsideDomainParticles=ProcessParticlesBoundaryIntersection;

  //particle produced in thermal radiation
  PIC::BC::UserDefinedParticleInjectionFunction=ThermalRadiation::InjectParticles;
}

void Radiation::PrintVariableList(FILE* fout,int DataSetNumber) {
  fprintf(fout,", \"Material Temparature\"");
}

void Radiation::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {
  bool gather_output_data=false;

  if (pipe==NULL) gather_output_data=true;
  else if (pipe->ThisThread==CenterNodeThread) gather_output_data=true;

  struct cDataExchengeBuffer {
    double MaterialTemparature;
  } buffer;

  if (gather_output_data==true) {
    buffer.MaterialTemparature=*((double*)(MaterialTemperatureOffset+CenterNode->GetAssociatedDataBufferPointer()));
  } 

  if ((PIC::ThisThread==0)||(pipe==NULL)) { 
    if ((CenterNodeThread!=0)&&(pipe!=NULL)) { 
      pipe->recv((char*)&buffer,sizeof(cDataExchengeBuffer),CenterNodeThread);
    }

    fprintf(fout,"%e ",buffer.MaterialTemparature);
  }
  else {
    pipe->send((char*)&buffer,sizeof(cDataExchengeBuffer));
  }
}

void Radiation::Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode) { 
  int i;
  double c,InterpolatedMaterialTemparature=0.0;

  //interpolate the sampled data
  for (i=0;i<nInterpolationCoeficients;i++) {
    c=InterpolationCoeficients[i];
    InterpolatedMaterialTemparature+=c*(*(double*)(MaterialTemperatureOffset+InterpolationList[i]->GetAssociatedDataBufferPointer()));
  }

  //store interpolated data
  *(double*)(MaterialTemperatureOffset+CenterNode->GetAssociatedDataBufferPointer())=InterpolatedMaterialTemparature; 
}

int Radiation::RequestStaticCellData(int offset) {
  MaterialTemperatureOffset=offset;
  return sizeof(double);
}

//scatter model particles from the boundaries of the computational domain 
int Radiation::ProcessParticlesBoundaryIntersection(long int ptr,double* xInit,double* vInit,int nIntersectionFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) { //(long int ptr,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  int res;

  //particles are reflected from all boundaries exept that with face:
  if (nIntersectionFace==0) {
    res=_PARTICLE_DELETED_ON_THE_FACE_;
  }
  else {
    //the particle is reflected back into the domain  
    double v[3];
   
    PIC::ParticleBuffer::GetV(v,ptr);

    switch (nIntersectionFace) {
    case 1:
      v[0]*=-1.0;
      break;
    case 2:case 3:
      v[1]*=-1.0;
      break;
    case 4:case 5:
      v[2]*=-1.0;
    } 

    PIC::ParticleBuffer::SetV(v,ptr);
    res=_PARTICLE_REJECTED_ON_THE_FACE_;
  }

  return res;
}

bool Radiation::Injection::BoundingBoxParticleInjectionIndicator(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  bool ExternalFaces[6];

  if (PIC::Mesh::mesh->ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) { 
    if (ExternalFaces[0]==true) return true;
  } 

  return false;
}

long int Radiation::Injection::BoundingBoxInjection(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  bool ExternalFaces[6];
  int nface,idim;
  long int newParticle;
  PIC::ParticleBuffer::byte *newParticleData;
  long int nInjectedParticles=0;
  double c0,c1,TimeCounter,ModelParticlesInjectionRate,ParticleWeight,LocalTimeStep,x[3],v[3],vr,theta;
  double x0[3],e0[3],e1[3];

  if (PIC::Mesh::mesh->ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    if (ExternalFaces[0]==true) {
      //inject new particles 
      TimeCounter=0.0;
      ModelParticlesInjectionRate;

      ParticleWeight=startNode->block->GetLocalParticleWeight(spec);
      LocalTimeStep=startNode->block->GetLocalTimeStep(spec);

      PIC::Mesh::mesh->GetBlockFaceCoordinateFrame_3D(x0,e0,e1,0,startNode);
   
      ModelParticlesInjectionRate=1.0;//kev* cm^2 /ns
      ModelParticlesInjectionRate*=Vector3D::Length(e0)*Vector3D::Length(e1); //kev*/ns 
//      ModelParticlesInjectionRate/=_NANO_; // kev/s 



     double T0=1.0; //equilibrium temeprature
     double I0=Opasity::GetSigma(T0)*Material::RadiationConstant*SpeedOfLight_cm*pow(T0,4)/(4*Pi);
     double ModelParticleDensity=I0/ParticleWeight;  
     double ModelParticleFlux=ModelParticleDensity*SpeedOfLight_cm/2.0;


      ModelParticlesInjectionRate=ModelParticleFlux*Vector3D::Length(e0)*Vector3D::Length(e1);;


      while ((TimeCounter+=-log(rnd())/ModelParticlesInjectionRate)<LocalTimeStep) {
        //generate the new particle position on the face
        for (idim=0,c0=rnd(),c1=rnd();idim<DIM;idim++) x[idim]=x0[idim]+c0*e0[idim]+c1*e1[idim];

        //generate a particle
        newParticle=PIC::ParticleBuffer::GetNewParticle();
        newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
        nInjectedParticles++;

        //generate particles' velocity
        double theta=asin(rnd()); 
        v[0]=SpeedOfLight_cm*cos(theta); 
        
     //   v[0]=sqrt(rnd())*SpeedOfLight; 
        vr=sqrt(SpeedOfLight_cm*SpeedOfLight_cm-v[0]*v[0]);

        theta=PiTimes2*rnd(); 
        v[1]=vr*sin(theta); 
        v[2]=vr*cos(theta);

        PIC::ParticleBuffer::SetX(x,newParticleData);
        PIC::ParticleBuffer::SetV(v,newParticleData);
        PIC::ParticleBuffer::SetI(spec,newParticleData);
        PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,newParticleData);

        //inject the particle into the system
        _PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_(newParticle,LocalTimeStep-TimeCounter,startNode);
      }
    }
  } 

  return nInjectedParticles;
}

void Radiation::IC::Set() {
  int iNode,i,j,k;
  PIC::Mesh::cDataBlockAMR *block;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  PIC::Mesh::cDataCenterNode *CenterNode;

  for (int iNode=0;iNode<PIC::DomainBlockDecomposition::nLocalBlocks;iNode++) { 
    node=PIC::DomainBlockDecomposition::BlockTable[iNode];
    block=node->block;

    if (block==NULL) continue;

    for (i=0;i<_BLOCK_CELLS_X_;i++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (k=0;k<_BLOCK_CELLS_Z_;k++) {
      CenterNode=block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k));
 
      *(double*)(MaterialTemperatureOffset+CenterNode->GetAssociatedDataBufferPointer())=T;
    }
  } 
}




long int Radiation::Injection::BoundingBoxInjection(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  return BoundingBoxInjection(0,startNode);
}


int Radiation::Mover(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  double v[3],x[3];
  int idim; 
  PIC::ParticleBuffer::byte *ParticleData;

  ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  PIC::ParticleBuffer::GetX(x,ParticleData);
  PIC::ParticleBuffer::GetV(v,ParticleData);

  for (idim=0;idim<3;idim++) x[idim]+=dtTotal*v[idim];

  if (x[0]<PIC::Mesh::mesh->xGlobalMin[0]) {
    //the particle existed from the domain 
    PIC::ParticleBuffer::DeleteParticle(ptr);
    return _PARTICLE_DELETED_ON_THE_FACE_; 
  }

  //check if the particle ger scattered from the walls 
  for (idim=0;idim<3;idim++) {
     double delta;

    if (x[idim]<PIC::Mesh::mesh->xGlobalMin[idim]) {
      delta=PIC::Mesh::mesh->xGlobalMin[idim]-x[idim];

      x[idim]=PIC::Mesh::mesh->xGlobalMin[idim]+delta;
      v[idim]*=-1.0;
    }
    
    if (x[idim]>PIC::Mesh::mesh->xGlobalMax[idim]) {
      delta=x[idim]-PIC::Mesh::mesh->xGlobalMax[idim];

      x[idim]=PIC::Mesh::mesh->xGlobalMax[idim]-delta;
      v[idim]*=-1.0;
    }
  }

  //find a new location of the particle and add it to the processed particle list 
  node=PIC::Mesh::mesh->findTreeNode(x,node);

  if (node==NULL) {
    PIC::ParticleBuffer::DeleteParticle(ptr);
    return _PARTICLE_DELETED_ON_THE_FACE_;
  }

  PIC::ParticleBuffer::SetX(x,ParticleData);
  PIC::ParticleBuffer::SetV(v,ParticleData);

  //attach the particle to the temporaty list
  int i,j,k;
  long int tempFirstCellParticle,*tempFirstCellParticlePtr;

  if (PIC::Mesh::mesh->FindCellIndex(x,i,j,k,node,false)==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located");

  tempFirstCellParticlePtr=node->block->tempParticleMovingListTable+i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k);
  tempFirstCellParticle=(*tempFirstCellParticlePtr);

  PIC::ParticleBuffer::SetNext(tempFirstCellParticle,ParticleData);
  PIC::ParticleBuffer::SetPrev(-1,ParticleData);

  if (tempFirstCellParticle!=-1) PIC::ParticleBuffer::SetPrev(ptr,tempFirstCellParticle);
  *tempFirstCellParticlePtr=ptr;

  return _PARTICLE_MOTION_FINISHED_;
}  

 

//absorption of radiation 
void Radiation::Absorption(long int ptr,long int& FirstParticleCell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  int i,j,k,spec;
  PIC::Mesh::cDataCenterNode *cell;
  double w,x[3],ParticleWeight,MaterialTemparature;
  double LocalTimeStep,GlobalWeight,ParticleWeightCorrection;

  PIC::ParticleBuffer::byte *ParticleData;

  ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  PIC::ParticleBuffer::GetX(x,ParticleData);
  spec=PIC::ParticleBuffer::GetI(ParticleData);
  ParticleWeightCorrection=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);

  LocalTimeStep=node->block->GetLocalTimeStep(spec);
  GlobalWeight=node->block->GetLocalParticleWeight(spec);
  
  if (PIC::Mesh::mesh->FindCellIndex(x,i,j,k,node,false)==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located");
  cell=node->block->GetCenterNode(i,j,k);

if (true) {
  MaterialTemparature= *(double*)(MaterialTemperatureOffset+cell->GetAssociatedDataBufferPointer());

  double dU,sigma=Opasity::GetSigma(MaterialTemparature); 
  w=exp(-sigma*SpeedOfLight_cm*LocalTimeStep); 

  if (w*ParticleWeightCorrection<1.0E-2) w=0.0; 

  dU=(1.0-w)*ParticleWeightCorrection*GlobalWeight; //GJ 
  MaterialTemparature+=dU/(Material::Density*cell->Measure)/Material::SpecificHeat; //keV 
  *(double*)(MaterialTemperatureOffset+cell->GetAssociatedDataBufferPointer())=MaterialTemparature;
}
else {
  double T0=*(double*)(MaterialTemperatureOffset+cell->GetAssociatedDataBufferPointer());
  
    double dU,sigma=Opasity::GetSigma(T0);
  w=exp(-sigma*SpeedOfLight_cm*LocalTimeStep);

  if (w*ParticleWeightCorrection<1.0E-2) w=0.0;

  dU=(1.0-w)*ParticleWeightCorrection*GlobalWeight; //GJ
  double T1=T0+dU/(Material::Density*cell->Measure)/Material::SpecificHeat; //keV
  
  double J,delta;

int cnt=4;

  do {
    J=-1+ParticleWeightCorrection*GlobalWeight/
      ((Material::Density*cell->Measure)/(Material::SpecificHeat)*
      400/pow(T0,4)*SpeedOfLight_cm*LocalTimeStep*exp(-100/pow(T0,3)*SpeedOfLight_cm*LocalTimeStep));  

double c0=ParticleWeightCorrection*GlobalWeight;
double c1=(Material::Density*cell->Measure)*(Material::SpecificHeat);
double c2=400/pow(T0,4)*SpeedOfLight_cm*LocalTimeStep;
double c3=exp(-100/pow(T0,3)*SpeedOfLight_cm*LocalTimeStep);

J=-1+c0/c1*c2*c3; 



    delta=(T0-T1)+ParticleWeightCorrection*GlobalWeight/
       ((Material::Density*cell->Measure)/(Material::SpecificHeat)*
        exp(-100/pow(T0,3)*SpeedOfLight_cm*LocalTimeStep)); 

double d0=ParticleWeightCorrection*GlobalWeight;
double d1=(Material::Density*cell->Measure)*(Material::SpecificHeat);
double d2=exp(-100/pow(T0,3)*SpeedOfLight_cm*LocalTimeStep);

delta=(T0-T1)+d0/d1*d2; 

    delta/=J; 

    T1=T0-0.1*delta;
   
   }
   while (cnt-->0);
}



  if (w==0.0) {
    //remove the particle from the simulation
    PIC::ParticleBuffer::DeleteParticle(ptr); 
  }
  else {
    //keep the particle in the simulation 
    PIC::ParticleBuffer::SetNext(FirstParticleCell,ptr);
    PIC::ParticleBuffer::SetPrev(-1,ptr);

    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(w*ParticleWeightCorrection,ParticleData);
    
    if (FirstParticleCell!=-1) PIC::ParticleBuffer::SetPrev(ptr,FirstParticleCell);
    FirstParticleCell=ptr;
  }
}        


long int Radiation::ThermalRadiation::InjectParticles() {
  int LocalCellNumber,i,j,k;
  long int res=0;

  for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::Mesh::mesh->ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    PIC::Mesh::cDataBlockAMR *block=node->block;

    if (!block) continue;

    for (k=0;k<_BLOCK_CELLS_Z_;k++) {
      for (j=0;j<_BLOCK_CELLS_Y_;j++) {
        for (i=0;i<_BLOCK_CELLS_X_;i++) {
          LocalCellNumber=_getCenterNodeLocalNumber(i,j,k);
          res+=InjectParticles(0,i,j,k,block->GetCenterNode(LocalCellNumber),block,node);
        }
      }
    }
  }

  return res;
} 

long int Radiation::ThermalRadiation::InjectParticles(int spec,int i,int j,int k,PIC::Mesh::cDataCenterNode* cell, PIC::Mesh::cDataBlockAMR* block, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) { 
  double MaterialTemparature;
  double anpart,dU,dU_kev,dU_GJ,ParticleWeight,WeightCorrectionFactor=1.0,LocalTimeStep,x[3],v[3];
  int npart;
  long int ptr;

  MaterialTemparature=*(double*)(MaterialTemperatureOffset+cell->GetAssociatedDataBufferPointer()); 

  ParticleWeight=node->block->GetLocalParticleWeight(spec);
  LocalTimeStep=node->block->GetLocalTimeStep(spec);

//  dU=Opasity::GetSigma(MaterialTemparature)*Material::RadiationConstant*SpeedOfLight*100.0*pow(MaterialTemparature,4)/(4.0*Pi); //GJ/g/ns  


  double sigma=100/pow(MaterialTemparature,3); //cm^2/g
  double a=0.01372; //GJ/(cm^3*keV^4)
  double rho=3.0; //g/cm^3
  double SpecificHeat=0.1; //GJ/g/keV

  dU=sigma*a*pow(Radiation::SpeedOfLight_cm,2)*pow(MaterialTemparature,4)/(4.0*Pi); //GJ/g/ns 
  dU_GJ=dU*LocalTimeStep; //GJ/g 
  dU=dU_GJ*cell->Measure*Material::Density; //GJ

//  dU_kev=dU*1.0E9*J2eV/1.0E3;
  anpart=dU/ParticleWeight;

/*
  dU=100.0/pow(MaterialTemparature,3)* 0.01372 * 29.98 *  pow(MaterialTemparature,4)/(4.0*Pi); //GJ/g/ns 


  dU_GJ=dU*LocalTimeStep/_NANO_; //GJ/g 

  dU=dU_GJ*cell->Measure*1.0E6*Material::Density; //GJ

  dU_kev=dU*1.0E9*J2eV/1.0E3; 
  anpart=dU_kev/ParticleWeight;
*/

  if (anpart>100.0) {
    WeightCorrectionFactor=anpart/100.0;
    anpart=100.0;
  }


  npart=(int)anpart;

  if (rnd()<npart-anpart) npart++;

  if (npart!=0) MaterialTemparature-=dU_GJ/Material::SpecificHeat; //keV 
  *(double*)(MaterialTemperatureOffset+cell->GetAssociatedDataBufferPointer())=MaterialTemparature;

  for (int ii=0;ii<npart;ii++) {
    x[0]=node->xmin[0]+(node->xmax[0]-node->xmin[0])*(i+rnd())/_BLOCK_CELLS_X_; 
    x[1]=node->xmin[1]+(node->xmax[1]-node->xmin[1])*(j+rnd())/_BLOCK_CELLS_Y_;
    x[2]=node->xmin[2]+(node->xmax[2]-node->xmin[2])*(k+rnd())/_BLOCK_CELLS_Z_;

    Vector3D::Distribution::Uniform(v,SpeedOfLight); 

    ptr=PIC::ParticleBuffer::GetNewParticle(block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)],true);

    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,ptr);
    PIC::ParticleBuffer::SetX(x,ptr);
    PIC::ParticleBuffer::SetV(v,ptr);
    PIC::ParticleBuffer::SetI(spec,ptr);
  }  

  return npart;
}






