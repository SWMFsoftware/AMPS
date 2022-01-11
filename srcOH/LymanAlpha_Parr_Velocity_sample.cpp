/*
This is a function meant to sample to Parallel velocity of neutral hydrogen along a LOS

Inputs: XYZ of a star in question
Outputs: A data file with the Parallel components of Velocity of hydrogen along the line of sight 

structure:
    Who knows yet


Plan for improvements:
1. Get it to compile using the Moon_Sample_Veloc... as a template
   a. I will be using a generic Coordinate in XYZ to as the target at first
  
2. Make it  use multiple XYZ coordinates in a text file perhaphs
3. Be able to give it an RA, Dec and have it translate to XYZ


*/




#include "pic.h"



OH::Sampling::LymanAlpha::cVelocitySampleBuffer *OH::Sampling::LymanAlpha::SampleBuffer=NULL;
int OH::Sampling::LymanAlpha::nTotalSampleDirections=0.0;
int OH::Sampling::LymanAlpha::nZenithPoints=0;


void OH::Sampling::LymanAlpha::Init() {


  //This is written for a single direction  (X,Y,Z) = (1,1,1)
  int nSampleDirections=1;  
  SampleBuffer=new cVelocitySampleBuffer[1];
  SampleBuffer[nTotalSampleDirections].lGSE[0]=1;
  SampleBuffer[nTotalSampleDirections].lGSE[1]=1;
  SampleBuffer[nTotalSampleDirections].lGSE[2]=1;

  nTotalSampleDirections++;
 
  
}

void OH::Sampling::LymanAlpha::Sampling() {


  //Initializing the varibles needed
  int Direction_i;
  int idim;
  double dl;
  //the ratio between the step of the integration procedure and the local cell size
  static const double IntegrationStep2CellSizeRatio=0.3;

  //initalizing the coordinates needed in the loop
  double startCoord[3]={0,0,0};
  //Adjust this for the coordinate that you need it to be
  double endCoord[3]={1,1,1};
  

  //Leaving this as a for loop because I would like to implement multiple directions 
  for (Direction_i=0;Direction_i<nTotalSampleDirections;Direction_i++) {

    //temporary sampling buffer and the buffer for the particle data
    cVelocitySampleBuffer tempSamplingBuffer;
    char tempParticleData[PIC::ParticleBuffer::ParticleDataLength];

    memcpy(&tempSamplingBuffer,SampleBuffer+Direction_i,sizeof(cVelocitySampleBuffer));

        
    //determine the limits of integration and initalize the nodes
    double IntegrationPathLength,xStart[3],xFinish[3];
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* xNode;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *xFinishNode;
    
    PIC::ColumnIntegration::FindIntegrationLimits(startCoord,endCoord,IntegrationPathLength,xStart,xNode,xFinish,xFinishNode);
    memcpy(startCoord,xStart,3*sizeof(double));

    while (xNode!=NULL) {
      dl=IntegrationStep2CellSizeRatio*xNode->GetCharacteristicCellSize();

      if (xNode->Thread==PIC::ThisThread) {
        //find the cell for sampling
        int i,j,k;
        long int ncell,ptr;
        PIC::ParticleBuffer::byte *ParticleData;
        int spec;
        double *v,*x,LocalParticleWeight,VelocityLineOfSight,HeliocentricRadialVelocity,c,cellMeasure,rHeliocentric,Speed;

        ncell=PIC::Mesh::mesh->FindCellIndex(startCoord,i,j,k,xNode);
        cellMeasure=xNode->block->GetCenterNode(ncell)->Measure;
//        dl=IntegrationStep2CellSizeRatio*xNode->GetCharacteristicCellSize();

        //sample particles
        ptr=xNode->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

	//cycle through the particles
        while (ptr!=-1) {
	  //Get the paricle data
          ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
          memcpy((void*)tempParticleData,(void*)ParticleData,PIC::ParticleBuffer::ParticleDataLength);

          spec=PIC::ParticleBuffer::GetI((PIC::ParticleBuffer::byte*)tempParticleData);
          v=PIC::ParticleBuffer::GetV((PIC::ParticleBuffer::byte*)tempParticleData);
          x=PIC::ParticleBuffer::GetX((PIC::ParticleBuffer::byte*)tempParticleData);

          LocalParticleWeight=xNode->block->GetLocalParticleWeight(spec)/cellMeasure;
          LocalParticleWeight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection((PIC::ParticleBuffer::byte*)tempParticleData);

          //calculate velocity component in the direction of the line of sight
	  //rx*vx + ry*vy + rz*vz = v||r
          for (idim=0,VelocityLineOfSight=0.0;idim<DIM;idim++) VelocityLineOfSight+=(v[idim])*endCoord[idim];
          i=(int)((VelocityLineOfSight+maxVelocityLimit)/VelocityBinWidth);
          if ((i>=0.0)&&(i<nVelocitySamplePoints)) tempSamplingBuffer.VelocityLineOfSight[spec][i]+=LocalParticleWeight*dl;

          
         

	  //Next particle
          ptr=PIC::ParticleBuffer::GetNext((PIC::ParticleBuffer::byte*)tempParticleData);
        }
      }   

      //next position along the line of sight
      for (idim=0;idim<3;idim++) startCoord[idim]+=dl*endCoord[idim];

      //break the loop if the new point is within the Moon or outside of the computational domain
      if (startCoord[0]*startCoord[0]+startCoord[1]*startCoord[1]+startCoord[2]*startCoord[2]
	  <endCoord[0]*endCoord[0]+endCoord[1]*endCoord[1]+endCoord[2]*endCoord[2]) break;

      xNode=PIC::Mesh::mesh->findTreeNode(startCoord,xNode);
    } 


    memcpy(SampleBuffer+Direction_i,&tempSamplingBuffer,sizeof(cVelocitySampleBuffer));
  }

  
}

void OH::Sampling::LymanAlpha::OutputSampledData(int DataOutputFileNumber) {

#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  char fname[_MAX_STRING_LENGTH_PIC_];
  int nZone,s,i,cnt;
  FILE *fout=NULL,*fMAP=NULL;

  const int mpiZoneExchangeBufferLength=PIC::nTotalSpecies*(nVelocitySamplePoints);
  double mpiZoneExchangeBuffer[mpiZoneExchangeBufferLength];

  //Need to redo this whole section for my problem --later problem
 
  if (PIC::ThisThread==0) {
    sprintf(fname,"%s/pic.LymanAlpha_Parr_Velocity_Dist.out=%i.dat",PIC::OutputDataFileDirectory,DataOutputFileNumber);
    fout=fopen(fname,"w");
    fprintf(fout,"VARIABLES=\"v\"");

    for (s=0;s<PIC::nTotalSpecies;s++) fprintf(fout,", \"f Velocity along the Line of Sight (%s) (in respect to the SUN)\" ",
					       PIC::MolecularData::GetChemSymbol(s);
    fprintf(fout,"\n");


    //create a map of the zones
    sprintf(fname,"%s/pic.Moon.VelocityDistribution.LymanAlpha.Map.out=%i.dat",PIC::OutputDataFileDirectory,DataOutputFileNumber);
    fMAP=fopen(fname,"w");

    fprintf(fMAP,"VARIABLES=\"X\", \"Y\", \"Z\", \"nZone\"");
    for (s=0;s<PIC::nTotalSpecies;s++) fprintf(fMAP,", \"Mean Velocity along the Line of Sight (%s) (in respect to the Earth)\"",
        PIC::MolecularData::GetChemSymbol(s));
    fprintf(fMAP,"\n");

    fprintf(fMAP,"ZONE T=\"Velocity Distribution Zone Map\"\n");
    fprintf(fMAP,"X=1, Y=1, Z=1, ZONETYPE=Ordered\n");
    fprintf(fMAP,"DATAPACKING=POINT\n");
  }



  for (nZone=0;nZone<nTotalSampleDirections;nZone++) {
    double XYZ[6]
   
    
    XYZ= SampleBuffer[nZone].lGSE;
    

 
    if (PIC::ThisThread==0) {
      fprintf(fout,"ZONE T=\"X=%e , Y=%e , Z=%e \",I=%i\n",XYZ[0],XYZ[1],XYZ[2],nVelocitySamplePoints);
      
    }

    //combine all sampled data on the root processor and print them into a file
    for (cnt=0,s=0;s<PIC::nTotalSpecies;s++) {
      for (i=0;i<nVelocitySamplePoints;i++) {
        mpiZoneExchangeBuffer[cnt++]=SampleBuffer[nZone].VelocityLineOfSight[s][i];
       }
    }


    if (PIC::ThisThread==0) {
      int thread;

      for (thread=1;thread<PIC::nTotalThreads;thread++) {
        MPI_Status status;

        MPI_Recv(mpiZoneExchangeBuffer,mpiZoneExchangeBufferLength,MPI_DOUBLE,thread,0,MPI_GLOBAL_COMMUNICATOR,&status);

        for (cnt=0,s=0;s<PIC::nTotalSpecies;s++) {
          for (i=0;i<nVelocitySamplePoints;i++) {
            SampleBuffer[nZone].VelocityLineOfSight[s][i]+=mpiZoneExchangeBuffer[cnt++];
        }
      }

      //normalize the distribution function
      for (s=0;s<PIC::nTotalSpecies;s++) {
        double w;

        for (w=0.0,i=0;i<nVelocitySamplePoints;i++) w+=SampleBuffer[nZone].VelocityLineOfSight[s][i];
        if (w>0.0) for (i=0;i<nVelocitySamplePoints;i++) SampleBuffer[nZone].VelocityLineOfSight[s][i]/=w*VelocityBinWidth;


      }


      //output sample data into a file
      for (i=0;i<nVelocitySamplePoints;i++) {
        fprintf(fout,"%e  ",-maxVelocityLimit+VelocityBinWidth*(i+0.5));

        for (s=0;s<PIC::nTotalSpecies;s++) fprintf(fout,"  %e  ", SampleBuffer[nZone].VelocityLineOfSight[s][i]);
        fprintf(fout,"\n");
      }


      fprintf(fMAP,"\n");
    }
    else MPI_Send(mpiZoneExchangeBuffer,mpiZoneExchangeBufferLength,MPI_DOUBLE,0,0,MPI_GLOBAL_COMMUNICATOR);

    //clear the sampling buffer
    for (s=0;s<PIC::nTotalSpecies;s++) {
      for (i=0;i<nVelocitySamplePoints;i++) SampleBuffer[nZone].VelocityLineOfSight[s][i]=0.0;
    }

    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  }


  //close the file
  if (PIC::Mesh::mesh->ThisThread==0) {
    fclose(fout);
    fclose(fMAP);
  }
#endif

  
}


