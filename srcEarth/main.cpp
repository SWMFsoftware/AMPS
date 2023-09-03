
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <time.h>
#include <cstdlib>


#include <sys/time.h>
#include <sys/resource.h>

//$Id$



//the particle class
#include "pic.h"
#include "constants.h"
#include "Earth.h"

#include "GeopackInterface.h"
#include "T96Interface.h"
#include "T05Interface.h"

int nZenithElements=200;
int nAzimuthalElements=200;

void amps_init();
void amps_init_mesh();
void amps_time_step();

#define _NIGHTLY_TEST__CUTOFF_  0
#define _NIGHTLY_TEST__LEGACY_ 1

#ifndef _NIGHTLY_TEST_
#define _NIGHTLY_TEST_ _NIGHTLY_TEST__CUTOFF_
#endif


#undef _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ 
#define _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_  _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ON_ 


//output directional cutoff rigidity and min energy 
void PrintDirectionalRegidityCutoffTitle(FILE* fout) {
  fprintf(fout,"TITLE=\"%s\"","DirectionalRegidityCutoff");
}

void PrintDirectionalRegidityCutoffVariableList(FILE* fout) {
  if (Earth::ModelMode!=Earth::CutoffRigidityMode) return;

  for (int i=0;i<PIC::nTotalSpecies;i++) fprintf(fout,", \"Cutoff Rigidity[%i]\", \"Min Energy(spec=%i)[MeV]\"",i,i);
}

void PrintDirectionalRigidityCutoffDataStateVector(FILE* fout,long int nZenithPoint,long int nAzimuthalPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalSphericalData *Sphere,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads) {
  int nInterpolationElement,nSurfaceElement,iZenith,iAzimuth;
  double InterpolationNormalization=0.0,InterpolationCoefficient;

  double CutoffRigidity=0.0;
  double InterpolatedInjectedParticleNumber=0.0,normInterpolatedInjectedParticleNumber=0.0;
  int InjectedParticleNumber;

  if (Earth::ModelMode!=Earth::CutoffRigidityMode) return;

  for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
    double t;

    CutoffRigidity=-1.0;

    for (nInterpolationElement=0;nInterpolationElement<SurfaceElementsInterpolationListLength;nInterpolationElement++) {
      nSurfaceElement=SurfaceElementsInterpolationList[nInterpolationElement];
      Sphere->GetSurfaceElementIndex(iZenith,iAzimuth,nSurfaceElement);

      t=Sphere->minRigidity[spec][nSurfaceElement];
 
      if (PIC::ThisThread!=0) {
        pipe->send(t);
      }
      else {
        if ((t>=0.0) && ((CutoffRigidity<0.0)||(t<CutoffRigidity)) ) CutoffRigidity=t;

        for (int thread=1;thread<PIC::nTotalThreads;thread++) {
          t=pipe->recv<double>(thread);

          if ((t>=0.0) && ((CutoffRigidity<0.0)||(t<CutoffRigidity)) ) CutoffRigidity=t;
        }
      }
    }

    if (PIC::ThisThread==0) {
      double momentum=-1.0,energy=-1.0; 

      if (CutoffRigidity>0.0) {
        momentum=CutoffRigidity*fabs(PIC::MolecularData::GetElectricCharge(spec))*1.0E9/SpeedOfLight;
        energy=Relativistic::Momentum2Energy(momentum,PIC::MolecularData::GetMass(spec))*J2MeV; 
      }

      fprintf(fout," %e  %e ",CutoffRigidity,energy);
    }
  }
}




void SampleIndividualLocations(int nMaxIterations) {
  int IterationCounter=0,localParticleGenerationFlag=0,globalParticleGenerationFlag;

  bool ShortTrajectoryFound=false;
  int ShortTrajectoryIndex,ShortTrajectoryIndexAll;


    PIC::SamplingMode=_TEMP_DISABLED_SAMPLING_MODE_;

  //estimate the total flux and rigidity in a set of the defined locations
  if (Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength!=0) {
    int LacalParticleNumber,GlobalParticleNumber;
    int nIngectedParticlePerIteration=Earth::CutoffRigidity::IndividualLocations::nTotalTestParticlesPerLocations/Earth::CutoffRigidity::IndividualLocations::nParticleInjectionIterations;
    int nTotalInjectedParticles=0;

    if (nIngectedParticlePerIteration==0) nIngectedParticlePerIteration=Earth::CutoffRigidity::IndividualLocations::nTotalTestParticlesPerLocations;

    if (PIC::ThisThread==0) {
      cout << "nTotalTestParticlesPerLocations=" << Earth::CutoffRigidity::IndividualLocations::nTotalTestParticlesPerLocations << endl;
      cout << "nParticleInjectionIterations=" << Earth::CutoffRigidity::IndividualLocations::nParticleInjectionIterations << endl;
      cout << "nIngectedParticlePerIteration=" << nIngectedParticlePerIteration << endl;
    }

    PIC::Mover::BackwardTimeIntegrationMode=_PIC_MODE_ON_;
    Earth::CutoffRigidity::DomainBoundaryParticleProperty::EnableSampleParticleProperty=true;

    if (Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleDomainBoundaryParticleProperty==true) Earth::CutoffRigidity::DomainBoundaryParticleProperty::Allocate(std::max(1,Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength));

    if (Earth::CutoffRigidity::IndividualLocations::CutoffRigidityTable.IsAllocated()==false) {
      Earth::CutoffRigidity::IndividualLocations::CutoffRigidityTable.init(PIC::nTotalSpecies,std::max(1,Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength));
    }

    Earth::CutoffRigidity::IndividualLocations::CutoffRigidityTable=-1.0;

    //allocate spherical objects for sampling directional rigidity cutoff
    cInternalSphericalData::SetGeneralSurfaceMeshParameters(150,150);
    Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable=new cInternalSphericalData[Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength];

    for (int i=0;i<Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength;i++) {
      double x[3]={0.0,0.0,0.0};

      Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[i].SetSphereGeometricalParameters(x,1.0);
      Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[i].TotalSurfaceElementNumber=Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[i].GetTotalSurfaceElementsNumber();  

      Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[i].PrintDataStateVector=PrintDirectionalRigidityCutoffDataStateVector;
      Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[i].PrintVariableList=PrintDirectionalRegidityCutoffVariableList; 
      Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[i].PrintTitle=PrintDirectionalRegidityCutoffTitle;

      double *xp=Earth::CutoffRigidity::IndividualLocations::xTestLocationTable[i];

      sprintf(Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[i].TitleMessage,
        "Directional rigidity cutoff: x=%e,%e,%e[m] (%e,%e,%e)R_Earth",
         xp[0],xp[1],xp[2],xp[0]/_RADIUS_(_EARTH_),xp[1]/_RADIUS_(_EARTH_),xp[2]/_RADIUS_(_EARTH_)); 


      Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[i].minRigidity=new double* [PIC::nTotalSpecies];
      Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[i].minRigidity[0]=new double[PIC::nTotalSpecies*Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[i].GetTotalSurfaceElementsNumber()];

      for (int spec=1;spec<PIC::nTotalSpecies;spec++) {
        Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[i].minRigidity[spec]=Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[i].minRigidity[spec-1]+Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[i].GetTotalSurfaceElementsNumber();
      }

      for (int ii=0;ii<PIC::nTotalSpecies*Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[i].TotalSurfaceElementNumber;ii++) {
        Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[i].minRigidity[0][ii]=-1.0;
      }
    }

    int nInjectionRounds=0;
    int iInjectedParticle;
    
    double dRlog=log(Earth::CutoffRigidity::IndividualLocations::MaxInjectionRigidityLimit/Earth::CutoffRigidity::IndividualLocations::MinInjectionRigidityLimit)/
        Earth::CutoffRigidity::IndividualLocations::nRigiditySearchIntervals;
        
    do {
      //reset the partilce generation flag
      localParticleGenerationFlag=0;
      
      //redefine nTotalTestParticlesPerLocations in case of Earth::CutoffRigidity::IndividualLocations::InjectionMode==Earth::CutoffRigidity::IndividualLocations::_rigidity_mesh
      if (Earth::CutoffRigidity::IndividualLocations::InjectionMode==Earth::CutoffRigidity::IndividualLocations::_rigidity_grid_injection) {
        Earth::CutoffRigidity::IndividualLocations::nTotalTestParticlesPerLocations=Earth::CutoffRigidity::IndividualLocations::nRigiditySearchIntervals*
            Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[0].GetTotalSurfaceElementsNumber();
      }

      //Inject new particles
      if (nTotalInjectedParticles<Earth::CutoffRigidity::IndividualLocations::nTotalTestParticlesPerLocations) {
        nTotalInjectedParticles+=nIngectedParticlePerIteration;
        nInjectionRounds++;
        
        int nParticles2Inject=nIngectedParticlePerIteration;
       
        if (nTotalInjectedParticles>Earth::CutoffRigidity::IndividualLocations::nTotalTestParticlesPerLocations) {
          nParticles2Inject=Earth::CutoffRigidity::IndividualLocations::nTotalTestParticlesPerLocations-(nTotalInjectedParticles-nIngectedParticlePerIteration);
        }

        //inject the new portion of the particles
        for (int spec=0;spec<PIC::nTotalSpecies;spec++) for (int iLocation=0;iLocation<Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength;iLocation++) {
          double x[3],v[3];
          int idim,iCell,jCell,kCell;
          cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode;

          for (idim=0;idim<3;idim++) x[idim]=Earth::CutoffRigidity::IndividualLocations::xTestLocationTable[iLocation][idim];
          startNode=PIC::Mesh::mesh->findTreeNode(x);

          if (startNode->Thread==PIC::ThisThread) {
            //generate a new particle velocity
            double mass,speed,energy,rigidity,momentum;

            static double logMinEnergyLimit=log(Earth::CutoffRigidity::IndividualLocations::MinEnergyLimit);
            static double logMaxEnergyLimit=log(Earth::CutoffRigidity::IndividualLocations::MaxEnergyLimit);

            mass=PIC::MolecularData::GetMass(spec);

            if (PIC::Mesh::mesh->FindCellIndex(x,iCell,jCell,kCell,startNode,false)==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located");

            for (int iNewParticle=0;iNewParticle<nParticles2Inject;iNewParticle++) {
              double momentum,rigidity;
              int iR,t;
              long int nZenithElement,nAzimuthalElement;
              
              iInjectedParticle=(nInjectionRounds-1)*nIngectedParticlePerIteration+iNewParticle;

              switch(Earth::CutoffRigidity::IndividualLocations::InjectionMode) {
              case Earth::CutoffRigidity::IndividualLocations::_rigidity_grid_injection:                
                //determine the velocity direction
                iR=iInjectedParticle/Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[iLocation].GetTotalSurfaceElementsNumber();
                t=iInjectedParticle%Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[iLocation].GetTotalSurfaceElementsNumber();
                
                nZenithElement=t/Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[iLocation].nAzimuthalSurfaceElements;
                nAzimuthalElement=t%Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[iLocation].nAzimuthalSurfaceElements;
                
                Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[iLocation].GetSurfaceElementRandomDirection(v,nZenithElement,nAzimuthalElement);
                
                //determine the velocity vector
                rigidity=Earth::CutoffRigidity::IndividualLocations::MinInjectionRigidityLimit*exp(dRlog*(iR+rnd()));

                if (rigidity>Earth::CutoffRigidity::IndividualLocations::MaxInjectionRigidityLimit) {
                  exit(__LINE__,__FILE__,"Error: the rigidity value exeeds the limit");
                }
       
                momentum=rigidity*fabs(PIC::MolecularData::GetElectricCharge(spec))*1.0E9/SpeedOfLight;
                speed=Relativistic::Momentum2Speed(momentum,mass);
                
                Vector3D::MultiplyScalar(speed,v);
                break;
                
              case Earth::CutoffRigidity::IndividualLocations::_energy_injection:
                energy=exp(logMinEnergyLimit+rnd()*(logMaxEnergyLimit-logMinEnergyLimit));

                speed=Relativistic::E2Speed(energy,mass);
                Vector3D::Distribution::Uniform(v,speed);
                Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[iLocation].GetSurfaceElementProjectionIndex(v,nZenithElement,nAzimuthalElement);
                break;
              case Earth::CutoffRigidity::IndividualLocations::_rigidity_injection:
                //def:          rigidity=(charge>0.0) ? momentum*SpeedOfLight/charge/1.0E9 : 0.0; //cutoff rigidity in SI -> GV (Moraal-2013-SSR)
                rigidity=Earth::CutoffRigidity::IndividualLocations::MinInjectionRigidityLimit+
                         rnd()*(Earth::CutoffRigidity::IndividualLocations::MaxInjectionRigidityLimit-Earth::CutoffRigidity::IndividualLocations::MinInjectionRigidityLimit);

                if (rigidity>Earth::CutoffRigidity::IndividualLocations::MaxInjectionRigidityLimit) {
                  exit(__LINE__,__FILE__,"Error: the rigidity value exeeds the limit");
                }

                
                momentum=rigidity*fabs(PIC::MolecularData::GetElectricCharge(spec))*1.0E9/SpeedOfLight;
                speed=Relativistic::Momentum2Speed(momentum,mass);
                Vector3D::Distribution::Uniform(v,speed);
                Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[iLocation].GetSurfaceElementProjectionIndex(v,nZenithElement,nAzimuthalElement);
                break;
              }


              //at this point, 'v' point into the direction from where the particles came. Hence, the velocity of the particles sgould be -v
              Vector3D::MultiplyScalar(-1.0,v);
   

              //generate a new particle
              long int newParticle=PIC::ParticleBuffer::GetNewParticle(startNode->block->FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)]);
              PIC::ParticleBuffer::byte *newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);

              PIC::ParticleBuffer::SetV(v,newParticleData);
              PIC::ParticleBuffer::SetX(x,newParticleData);
              PIC::ParticleBuffer::SetI(spec,newParticleData);

              //set the particle generation flag
              localParticleGenerationFlag=1;

              //apply condition of tracking the particle
              if (_PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_) {
                PIC::ParticleTracker::InitParticleID(newParticleData);
                PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x,v,spec,newParticleData,(void*)startNode);
              }

              *((int*)(newParticleData+Earth::CutoffRigidity::ParticleDataOffset::OriginLocationIndex))=iLocation;
              *((double*)(newParticleData+Earth::CutoffRigidity::ParticleDataOffset::OriginalSpeed))=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

              //set the initial value of the integrate path length
              if (Earth::CutoffRigidity::IntegratedPathLengthOffset!=-1) {
                *((double*)(newParticleData+Earth::CutoffRigidity::IntegratedPathLengthOffset))=0.0;
              }

              //set initial value for the integration time 
              if (Earth::CutoffRigidity::IntegratedTimeOffset!=-1) {
                *((double*)(newParticleData+Earth::CutoffRigidity::IntegratedTimeOffset))=0.0;
              }

              //set the velocity direction ID
              if (Earth::CutoffRigidity::ParticleDataOffset::OriginalVelocityDirectionID!=-1) {
                int id;
               // long int nZenithElement,nAzimuthalElement;

                //Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[iLocation].GetSurfaceElementProjectionIndex(v,nZenithElement,nAzimuthalElement);
                id=Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[iLocation].GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement); 

                *((int*)(newParticleData+Earth::CutoffRigidity::ParticleDataOffset::OriginalVelocityDirectionID))=id;
              }
   

              //set up the particle rigidity
              if (Earth::CutoffRigidity::InitialRigidityOffset!=-1) {
                double momentum,charge,rigidity;

                charge=fabs(PIC::MolecularData::GetElectricCharge(spec));

                momentum=Relativistic::Speed2Momentum(speed,mass);
                rigidity=(charge>0.0) ? momentum*SpeedOfLight/charge/1.0E9 : 0.0; //cutoff rigidity in SI -> GV (Moraal-2013-SSR) 

                //rigidity=(charge>0.0) ? momentum/charge : 0.0;


                if (Earth::CutoffRigidity::IndividualLocations::InjectionMode==Earth::CutoffRigidity::IndividualLocations::_rigidity_grid_injection) {
                  if (rigidity>Earth::CutoffRigidity::IndividualLocations::MaxInjectionRigidityLimit) {
                    exit(__LINE__,__FILE__,"Error: the rigidity value exeeds the limit");
                  }
                }

                *((double*)(newParticleData+Earth::CutoffRigidity::InitialRigidityOffset))=rigidity;
              }

              //save the original location of the particle
              if (Earth::CutoffRigidity::InitialLocationOffset!=-1) {
                memcpy(newParticleData+Earth::CutoffRigidity::InitialLocationOffset,x,3*sizeof(double));
              }


            }
          }
        }
      }


      //preform the next iteration
      amps_time_step();

    //search for particles that left geospace
    if (Earth::GeospaceFlag::offset!=-1) {
      //loop through all cells
      for (int iLocalNode=0;iLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;iLocalNode++) {
        cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::DomainBlockDecomposition::BlockTable[iLocalNode];
        PIC::Mesh::cDataBlockAMR *block;
        PIC::Mesh::cDataCenterNode *CenterNode;
        long int ParticleList,ptr;
        int i,j,k,nd;
        char *offset;

        block=node->block;
        if ((block=node->block)==NULL) continue;

        for (k=0;k<_BLOCK_CELLS_Z_;k++) {
          for (j=0;j<_BLOCK_CELLS_Y_;j++) {
            for (i=0;i<_BLOCK_CELLS_X_;i++) {
              nd=PIC::Mesh::mesh->getCenterNodeLocalNumber(i,j,k);

              if ((CenterNode=node->block->GetCenterNode(nd))==NULL) continue;
              offset=CenterNode->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::CenterNodeAssociatedDataOffsetBegin+PIC::CPLR::DATAFILE::MULTIFILE::CurrDataFileOffset;

              if (*((double*)(offset+Earth::GeospaceFlag::offset))==0.0) {
                //the cell is outside of the geospace: search for particles and removed those found

                ParticleList=block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

                while (ParticleList!=-1) {
                  ptr=ParticleList;
                  ParticleList=PIC::ParticleBuffer::GetNext(ParticleList);

                  //account for the particle in calculation of the rigidity
                  double x[3],v[3];

                  PIC::ParticleBuffer::GetV(v,ptr);
                  PIC::ParticleBuffer::GetX(x,ptr);

                  Earth::CutoffRigidity::ProcessOutsideDomainParticles(ptr,x,v,-1,node);
                  PIC::ParticleBuffer::DeleteParticle(ptr);
                }

                block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]=-1;
              }
            }
          }
        }
      }
    }

        //search for particles that spend too much time in the simulation
        ShortTrajectoryFound=false;

      for (int iLocalNode=0;iLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;iLocalNode++) {
        cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::DomainBlockDecomposition::BlockTable[iLocalNode];
        PIC::Mesh::cDataBlockAMR *block;
        PIC::Mesh::cDataCenterNode *CenterNode;
        long int ParticleList,ptr;
        int i,j,k,nd;
        char *offset;
        double x[3];

        block=node->block;
        if ((block=node->block)==NULL) continue;


        for (int k=0;k<_BLOCK_CELLS_Z_;k++) {
          for (int j=0;j<_BLOCK_CELLS_Y_;j++) {
            for (int i=0;i<_BLOCK_CELLS_X_;i++) {
              long int ParticleList=block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

              while (ParticleList!=-1) {
                ptr=ParticleList;
                ParticleList=PIC::ParticleBuffer::GetNext(ParticleList);

                double t,l,v[3];
                PIC::ParticleBuffer::byte *ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);

                t=*((double*)(ParticleData+Earth::CutoffRigidity::IntegratedTimeOffset));
                PIC::ParticleBuffer::GetX(x,ptr);

                if (t*t*Vector3D::DotProduct(v,v)>Earth::CutoffRigidity::MaxIntegrationLength*Earth::CutoffRigidity::MaxIntegrationLength) {
                  PIC::ParticleBuffer::DeleteParticle(ptr,block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]);
                }
                else {
                 *((double*)(ParticleData+Earth::CutoffRigidity::IntegratedTimeOffset))+=PIC::ParticleWeightTimeStep::GlobalTimeStep[PIC::ParticleBuffer::GetI(ptr)];
                 if (Earth::CutoffRigidity::SearchShortTrajectory==true) ShortTrajectoryFound=true;
              }
            }
         }
       }
     }
   }

      static int LastDataOutputFileNumber=-1;

      if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
        PIC::RequiredSampleLength*=2;
        if (PIC::RequiredSampleLength>50000) PIC::RequiredSampleLength=50000;


        LastDataOutputFileNumber=PIC::DataOutputFileNumber;
        if (PIC::Mesh::mesh->ThisThread==0) cout << "The new sample length is " << PIC::RequiredSampleLength << endl;
      }


      //get the total number of particles in the system
      LacalParticleNumber=PIC::ParticleBuffer::GetAllPartNum();
      MPI_Allreduce(&LacalParticleNumber,&GlobalParticleNumber,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

      //determine whether any particle has been generated during the current iteration
      MPI_Allreduce(&localParticleGenerationFlag,&globalParticleGenerationFlag,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

      if (globalParticleGenerationFlag!=0) {
        //at least one particle has been generated -> reset the iteration counter
        IterationCounter=0;
      }
      else {
        //increment the iteration counter
        IterationCounter++;
      }


    ShortTrajectoryIndex=(ShortTrajectoryFound==true) ? 1 : 0;
    MPI_Allreduce(&ShortTrajectoryIndex,&ShortTrajectoryIndexAll,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

    if (PIC::ThisThread==0) cout << "IterationCounter=" << IterationCounter << ", nMaxIterations=" << nMaxIterations << ", GlobalParticleNumber=" << GlobalParticleNumber << ", ShortTrajectoryIndexAll=" << ShortTrajectoryIndexAll << endl;


    }
    while (((GlobalParticleNumber!=0)&&(IterationCounter<nMaxIterations)) || (ShortTrajectoryIndexAll>0));

    //delete all particles that still present in the system
    PIC::ParticleBuffer::DeleteAllParticles();

    //print out the results 
    //output directional regidity cutoff
    for (int iLocation=0;iLocation<Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength;iLocation++) {
      char fname[200];

      sprintf(fname,"%s/directional-rigidity-cutoff--point-%i.dat",PIC::OutputDataFileDirectory,iLocation);
      Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[iLocation].PrintSurfaceData(fname,0,true);
    }

    int nel=Earth::CutoffRigidity::IndividualLocations::CutoffRigidityTable.GetElementNumber();
    double *buff=new double [PIC::nTotalThreads*nel];
    
    MPI_Gather(Earth::CutoffRigidity::IndividualLocations::CutoffRigidityTable.GetBufferPointer(),nel,MPI_DOUBLE,buff,nel,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);

if (PIC::nTotalSpecies!=1) exit(__LINE__,__FILE__,"Error: heen to be generalized from more then 1 species");

    for (int el=0;el<nel;el++) {
      for (int thread=1;thread<PIC::nTotalThreads;thread++) {
        if ((buff[el]<=0.0)||(buff[el]>buff[el+thread*nel])) buff[el]=buff[el+thread*nel];
      }
    }

    if (PIC::ThisThread==0) {
      memcpy(Earth::CutoffRigidity::IndividualLocations::CutoffRigidityTable.GetBufferPointer(),buff,nel*sizeof(double));
    }

    delete [] buff;

    MPI_Bcast(Earth::CutoffRigidity::IndividualLocations::CutoffRigidityTable.GetBufferPointer(),nel,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR); 
  
    if (PIC::ThisThread==0) {
      ofstream fout("res.dat"); 

      cout << "iloc\tspec\tcutoff rigidity\n";
      fout << "iloc\tspec\tcutoff rigidity\n";

      for (int iloc=0;iloc<Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength;iloc++) {
        for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
          cout << iloc << "\t" << spec << "\t" << Earth::CutoffRigidity::IndividualLocations::CutoffRigidityTable(spec,iloc) << endl; 
          fout << iloc << "\t" << spec << "\t" << Earth::CutoffRigidity::IndividualLocations::CutoffRigidityTable(spec,iloc) << endl;
        }
      }

      if (Earth::CutoffRigidity::IndividualLocations::InjectionMode==Earth::CutoffRigidity::IndividualLocations::_energy_injection) {
        cout << "iloc\tspec\tflux\n";
        fout << "iloc\tspec\tflux\n";

        for (int iloc=0;iloc<Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength;iloc++) {
          for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
           cout << iloc << "\t" << spec << "\t" <<
             Earth::CutoffRigidity::IndividualLocations::SampledFluxTable(spec,iloc)*4.0*Pi/Earth::CutoffRigidity::IndividualLocations::nTotalTestParticlesPerLocations << endl; 

           fout << iloc << "\t" << spec << "\t" <<
             Earth::CutoffRigidity::IndividualLocations::SampledFluxTable(spec,iloc)*4.0*Pi/Earth::CutoffRigidity::IndividualLocations::nTotalTestParticlesPerLocations << endl;
          }
        }
      }

      fout.close();
    }


    //determine the flux and eneregy spectra of the energetic particles in the poins of the observation
    if (Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleDomainBoundaryParticleProperty==true) {
      double v[3],KineticEnergy,Speed,DiffFlux,dSurface,norm;
      int offset,iTestsLocation,spec,i,j,iface,iTable,jTable,Index,iBit,iByte,iE;

      const int nTotalEnergySpectrumIntervals=25;
      const double logMinEnergyLimit=log10(Earth::CutoffRigidity::IndividualLocations::MinEnergyLimit);
      const double logMaxEnergyLimit=log10(Earth::CutoffRigidity::IndividualLocations::MaxEnergyLimit);
      const double dE=(logMaxEnergyLimit-logMinEnergyLimit)/nTotalEnergySpectrumIntervals;

      //allocate the data buffers
      //TotalFlux[iLocation][spec]
      //EnergySpectrum[iLocation][spec][iEnergyInterval]

      array_3d<double> EnergySpectrum(Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength,PIC::nTotalSpecies,nTotalEnergySpectrumIntervals);
      array_2d<double> TotalFlux(Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength,PIC::nTotalSpecies);

      EnergySpectrum=0.0;
      TotalFlux=0.0;

      //calculate the flux and energey spectrum
      for (iTestsLocation=0;iTestsLocation<Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength;iTestsLocation++) {
        for (spec=0;spec<PIC::nTotalSpecies;spec++) {
          for (iface=0;iface<6;iface++) {
            //surface area of the element of the surface mesh that covers the boundary of the computational domain
            double lx,ly,lz;

            switch (iface) {
            case 0:case 1:
              ly=(PIC::Mesh::mesh->xGlobalMax[1]-PIC::Mesh::mesh->xGlobalMin[1])/Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleMaskNumberPerSpatialDirection;
              lz=(PIC::Mesh::mesh->xGlobalMax[2]-PIC::Mesh::mesh->xGlobalMin[2])/Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleMaskNumberPerSpatialDirection;

              dSurface=ly*lz;
              break;
            case 2:case 3:
              lx=(PIC::Mesh::mesh->xGlobalMax[0]-PIC::Mesh::mesh->xGlobalMin[0])/Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleMaskNumberPerSpatialDirection;
              lz=(PIC::Mesh::mesh->xGlobalMax[2]-PIC::Mesh::mesh->xGlobalMin[2])/Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleMaskNumberPerSpatialDirection;

              dSurface=lx*lz;
              break;
            case 4:case 5:
              lx=(PIC::Mesh::mesh->xGlobalMax[0]-PIC::Mesh::mesh->xGlobalMin[0])/Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleMaskNumberPerSpatialDirection;
              ly=(PIC::Mesh::mesh->xGlobalMax[1]-PIC::Mesh::mesh->xGlobalMin[1])/Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleMaskNumberPerSpatialDirection;

              dSurface=ly*lz;
            }

            //loop through the mesh that covers face 'iface' on the computational domain
            for (iTable=0;iTable<Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleMaskNumberPerSpatialDirection;iTable++) {
              for (jTable=0;jTable<Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleMaskNumberPerSpatialDirection;jTable++) {
                Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleTable(spec,iTestsLocation,iface,iTable,jTable).Gather();

                for (iByte=0;iByte<Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleTable(spec,iTestsLocation,iface,iTable,jTable).FlagTableLength[0];iByte++) for (iBit=0;iBit<8;iBit++) {
                  Index=iBit+8*iByte;

                  if (Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleTable(spec,iTestsLocation,iface,iTable,jTable).Test(Index)==true) {
                    //at least one particle that corrsponds to 'Index' has been detected. Add a contribution of such particles to the total energy spectrum and flux as observed at the point of the observation 'iTestsLocation'

                    Earth::CutoffRigidity::DomainBoundaryParticleProperty::ConvertVelocityVectorIndex2Velocity(spec,v,iface,Index);
                    Speed=Vector3D::Length(v);
                    KineticEnergy=Relativistic::Speed2E(Speed,PIC::MolecularData::GetMass(spec));

                    //probability density that particles has velocity 'v'
                    DiffFlux=GCR_BADAVI2011ASR::Hydrogen::GetDiffFlux(KineticEnergy);

                    //determine the contributio of the particles into the 'observed' flux and energy spectrum
                    TotalFlux(iTestsLocation,spec)+=DiffFlux*dSurface;

                    //determine contribution of the particles to the energy flux
                    iE=(log10(KineticEnergy)-logMinEnergyLimit)/dE;
                    if (iE<0) iE=0;
                    if (iE>=nTotalEnergySpectrumIntervals) iE=nTotalEnergySpectrumIntervals-1;

                    EnergySpectrum(iTestsLocation,spec,iE)+=DiffFlux*dSurface;
                  }
                }
              }
            }
          }

          //normalize the energy spectrum
          for (iE=0,norm=0.0;iE<nTotalEnergySpectrumIntervals;iE++) norm+=EnergySpectrum(iTestsLocation,spec,iE)*dE;
          if (norm>0) for (iE=0;iE<nTotalEnergySpectrumIntervals;iE++) EnergySpectrum(iTestsLocation,spec,iE)/=norm;
        }
      }


      //output sampled particles flux and energy spectrum
      if (PIC::ThisThread==0) {
        //sampled energy spectrum
        FILE *fout;
        int spec,iTestsLocation,iE;
        char fname[400];

        for (spec=0;spec<PIC::nTotalSpecies;spec++) {
          sprintf(fname,"%s/EnergySpectrum[s=%i].dat",PIC::OutputDataFileDirectory,spec);
          fout=fopen(fname,"w");

          fprintf(fout,"VARIABLES=\"log10(Kinetic Energy[MeV]\"");
          for (iTestsLocation=0;iTestsLocation<Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength;iTestsLocation++) fprintf(fout,", \"Spectrum (iTestsLocation=%i)\"",iTestsLocation);
          fprintf(fout,"\n");

          for (iE=0;iE<nTotalEnergySpectrumIntervals;iE++) {
            double log10e=logMinEnergyLimit+iE*dE;
            double e=pow(10,log10e);

            e*=J2MeV;
            fprintf(fout,"%e  ",e);

            for (iTestsLocation=0;iTestsLocation<Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength;iTestsLocation++) fprintf(fout,"%e  ",EnergySpectrum(iTestsLocation,spec,iE));
            fprintf(fout,"\n");
          }

          fclose(fout);
        }

        //The total energetic particle flux
        for (spec=0;spec<PIC::nTotalSpecies;spec++) for (iTestsLocation=0;iTestsLocation<Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength;iTestsLocation++) {
          printf("spec=%i, iTestsLocation=%i: Flux=%e\n",spec,iTestsLocation,TotalFlux(iTestsLocation,spec));
        }
      }
    }
    else {
      //non-uniform distribution of the injected particles
      if (PIC::ThisThread==0) printf("%i, %s: Error: not implemented\n",__LINE__,__FILE__); //exit(__LINE__,__FILE__,"Error: not implemented");
      return; 
    }

  }


/*
  //output directional regidity cutoff
  for (int iLocation=0;iLocation<Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength;iLocation++) {
    char fname[200];

    sprintf(fname,"directional-rigidity-cutoff--point-%i.dat",iLocation); 
    Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[iLocation].PrintSurfaceData(fname,0,true);
  }
*/

  //release sampling buffers
 // Earth::CutoffRigidity::DomainBoundaryParticleProperty::Deallocate();
}



void SampleSphericalMaplLocations(double Radius,int nMaxIterations) {
  int IterationCounter=0,localParticleGenerationFlag=0,globalParticleGenerationFlag;
  int iLocation;
  double x[3]={0.0,0.0,0.0},v[3];

  PIC::Mover::BackwardTimeIntegrationMode=_PIC_MODE_ON_;


  int nTotalInjectedParticlePerPoint=25;

//  Earth::CutoffRigidity::IndividualLocations::nTotalTestParticlesPerLocations=2000;

  cInternalSphericalData Sphere;
  Sphere.SetGeneralSurfaceMeshParameters(nZenithElements,nAzimuthalElements);
  Sphere.SetSphereGeometricalParameters(x,Radius);

  //  Earth::CutoffRigidity::DomainBoundaryParticleProperty::Allocate(nZenithElements*nAzimuthalElements);

  double Speed,DiffFlux,dSurface,norm,KineticEnergy;
  int offset,iTestsLocation,spec,i,j,iface,iTable,jTable,Index,iBit,iByte,iE;


  const int nTotalEnergySpectrumIntervals=50;
  const double logMinEnergyLimit=log10(Earth::CutoffRigidity::IndividualLocations::MinEnergyLimit);
  const double logMaxEnergyLimit=log10(Earth::CutoffRigidity::IndividualLocations::MaxEnergyLimit);
  const double dE=(logMaxEnergyLimit-logMinEnergyLimit)/nTotalEnergySpectrumIntervals;


  array_2d<double> TotalFlux(nZenithElements*nAzimuthalElements,PIC::nTotalSpecies);
  array_3d<double> EnergySpectrum(nZenithElements*nAzimuthalElements,PIC::nTotalSpecies,nTotalEnergySpectrumIntervals);

  //set the values of the buffers to zero
  TotalFlux=0.0;
  EnergySpectrum=0.0;


  //estimate the total flux and rigidity in a set of the defined locations
  int LacalParticleNumber,GlobalParticleNumber=0;
  int nIngectedParticlePerIteration=Earth::CutoffRigidity::IndividualLocations::nTotalTestParticlesPerLocations/Earth::CutoffRigidity::IndividualLocations::nParticleInjectionIterations;
  int nTotalInjectedParticles=0;
  bool ShortTrajectoryFound=false;
  int ShortTrajectoryIndex,ShortTrajectoryIndexAll;

  if (nIngectedParticlePerIteration==0) nIngectedParticlePerIteration=Earth::CutoffRigidity::IndividualLocations::nTotalTestParticlesPerLocations;

  PIC::Mover::BackwardTimeIntegrationMode=_PIC_MODE_ON_;
  Earth::CutoffRigidity::DomainBoundaryParticleProperty::EnableSampleParticleProperty=true;

  //allocate the sampling buffer
  if (Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleDomainBoundaryParticleProperty==true) {
    Earth::CutoffRigidity::DomainBoundaryParticleProperty::Allocate(nAzimuthalElements*nAzimuthalElements);
  }

  if (Earth::CutoffRigidity::CutoffRigidityTable.IsAllocated()==true) {
    Earth::CutoffRigidity::CutoffRigidityTable.Deallocate(); 
    Earth::CutoffRigidity::InjectedParticleMap.Deallocate();
    Earth::CutoffRigidity::MaxEnergyInjectedParticles.Deallocate();
  }

  Earth::CutoffRigidity::MaxEnergyInjectedParticles.init(PIC::nTotalSpecies,nZenithElements*nAzimuthalElements);
  Earth::CutoffRigidity::MaxEnergyInjectedParticles=0.0;

  Earth::CutoffRigidity::CutoffRigidityTable.init(PIC::nTotalSpecies,nZenithElements*nAzimuthalElements);
  Earth::CutoffRigidity::CutoffRigidityTable=-1.0;

  Earth::CutoffRigidity::InjectedParticleMap.init(PIC::nTotalSpecies,nZenithElements*nAzimuthalElements);
  Earth::CutoffRigidity::InjectedParticleMap=0;

  if (PIC::ThisThread==0) {
    cout << "nIngectedParticlePerIteration=" << nIngectedParticlePerIteration << endl;
    cout << "nTotalTestParticlesPerLocations=" << Earth::CutoffRigidity::IndividualLocations::nTotalTestParticlesPerLocations << endl;
    cout << "nParticleInjectionIterations=" << Earth::CutoffRigidity::IndividualLocations::nParticleInjectionIterations << endl;
  }

  //deallocate the individual location sampling buffer
  Earth::CutoffRigidity::IndividualLocations::CutoffRigidityTable.Deallocate();


  do { //while there are particles in the system
    localParticleGenerationFlag=0;


    int sumTotalInjectedParticles;
    MPI_Allreduce(&nTotalInjectedParticles,&sumTotalInjectedParticles,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

    if (sumTotalInjectedParticles<Earth::CutoffRigidity::IndividualLocations::nTotalTestParticlesPerLocations*nZenithElements*nAzimuthalElements*PIC::nTotalSpecies) {
      //inject new particles in case it is needed
      localParticleGenerationFlag=1;

      //inject the new portion of the particles
      for (int spec=0;spec<PIC::nTotalSpecies;spec++) for (int iZenithElement=0;iZenithElement<nZenithElements;iZenithElement++) for (int iAzimutalElement=0;iAzimutalElement<nAzimuthalElements;iAzimutalElement++) {
        int idim,iCell,jCell,kCell;
        cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode;

        iLocation=Sphere.GetLocalSurfaceElementNumber(iZenithElement,iAzimutalElement);

        //generate a new particle velocity
        double mass,speed,energy,rigidity,momentum;
        int nd;

        static double logMinEnergyLimit=log(Earth::CutoffRigidity::IndividualLocations::MinEnergyLimit);
        static double logMaxEnergyLimit=log(Earth::CutoffRigidity::IndividualLocations::MaxEnergyLimit);

        mass=PIC::MolecularData::GetMass(spec);

        for (int iNewParticle=0;iNewParticle<nIngectedParticlePerIteration;iNewParticle++) {
          Sphere.GetSurfaceElementRandomPoint(x,iZenithElement,iAzimutalElement);
          startNode=PIC::Mesh::mesh->findTreeNode(x);

          if (startNode->Thread!=PIC::ThisThread) continue;

          if ((nd=PIC::Mesh::mesh->FindCellIndex(x,iCell,jCell,kCell,startNode,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located");

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
          PIC::Mesh::cDataCenterNode *cell;

          cell=startNode->block->GetCenterNode(nd);

          if (cell==NULL) exit(__LINE__,__FILE__,"Error: the cell measure is not initialized");
          if (cell->Measure<=0.0) exit(__LINE__,__FILE__,"Error: the cell measure is not initialized");
#endif

          //generate energy of the new particle
//          energy=exp(logMinEnergyLimit+rnd()*(logMaxEnergyLimit-logMinEnergyLimit));

 //         speed=Relativistic::E2Speed(energy,mass);
 //         Vector3D::Distribution::Uniform(v,speed);


              double momentum,rigidity;

              switch(Earth::CutoffRigidity::IndividualLocations::InjectionMode) {
              case Earth::CutoffRigidity::IndividualLocations::_energy_injection:
                energy=exp(logMinEnergyLimit+rnd()*(logMaxEnergyLimit-logMinEnergyLimit));

                speed=Relativistic::E2Speed(energy,mass);
                Vector3D::Distribution::Uniform(v,speed);
                break;
              case Earth::CutoffRigidity::IndividualLocations::_rigidity_injection:
                rigidity=Earth::CutoffRigidity::IndividualLocations::MinInjectionRigidityLimit+
                         rnd()*(Earth::CutoffRigidity::IndividualLocations::MaxInjectionRigidityLimit-Earth::CutoffRigidity::IndividualLocations::MinInjectionRigidityLimit);


  //def:          rigidity=(charge>0.0) ? momentum*SpeedOfLight/charge/1.0E9 : 0.0; //cutoff rigidity in SI -> GV (Moraal-2013-SSR)  


                momentum=rigidity*fabs(PIC::MolecularData::GetElectricCharge(spec))*1.0E9/SpeedOfLight;
                speed=Relativistic::Momentum2Speed(momentum,mass);
                Vector3D::Distribution::Uniform(v,speed);
                break;
              }


          //generate a new particle
          long int newParticle=PIC::ParticleBuffer::GetNewParticle(startNode->block->FirstCellParticleTable[iCell+_BLOCK_CELLS_X_*(jCell+_BLOCK_CELLS_Y_*kCell)]);
          PIC::ParticleBuffer::byte *newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);

          nTotalInjectedParticles++;
          Earth::CutoffRigidity::InjectedParticleMap(spec,iLocation)=1+Earth::CutoffRigidity::InjectedParticleMap(spec,iLocation); 

          if (x[0]*v[0]+x[1]*v[1]+x[2]*v[2]<0.0) v[0]=-v[0],v[1]=-v[1],v[2]=-v[2]; 

          PIC::ParticleBuffer::SetV(v,newParticleData);
          PIC::ParticleBuffer::SetX(x,newParticleData);
          PIC::ParticleBuffer::SetI(spec,newParticleData);

          if (energy>Earth::CutoffRigidity::MaxEnergyInjectedParticles(spec,iLocation)) Earth::CutoffRigidity::MaxEnergyInjectedParticles(spec,iLocation)=energy*J2MeV;

          //set the particle generation flag
          localParticleGenerationFlag=1;

          //apply condition of tracking the particle
          if (_PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_) {
            PIC::ParticleTracker::InitParticleID(newParticleData);
            PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x,v,spec,newParticleData,(void*)startNode);
          }

          *((int*)(newParticleData+Earth::CutoffRigidity::ParticleDataOffset::OriginLocationIndex))=iLocation;
          *((double*)(newParticleData+Earth::CutoffRigidity::ParticleDataOffset::OriginalSpeed))=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

          //set the initial value of the integrate path length
          if (Earth::CutoffRigidity::IntegratedPathLengthOffset!=-1) {
            *((double*)(newParticleData+Earth::CutoffRigidity::IntegratedPathLengthOffset))=0.0;
          }

          //set initial value for the integration time
          if (Earth::CutoffRigidity::IntegratedTimeOffset!=-1) {
            *((double*)(newParticleData+Earth::CutoffRigidity::IntegratedTimeOffset))=0.0;
          }


          //set up the particle rigidity
          if (Earth::CutoffRigidity::InitialRigidityOffset!=-1) {
            double momentum,charge,rigidity;

            charge=fabs(PIC::MolecularData::GetElectricCharge(spec));

            momentum=Relativistic::Speed2Momentum(speed,mass);
            rigidity=(charge>0.0) ? momentum*SpeedOfLight/charge/1.0E9 : 0.0; //cutoff rigidity in SI -> GV (Moraal-2013-SSR)  

            *((double*)(newParticleData+Earth::CutoffRigidity::InitialRigidityOffset))=rigidity;
          }

          //save the original location of the particle
          if (Earth::CutoffRigidity::InitialLocationOffset!=-1) {
            memcpy(newParticleData+Earth::CutoffRigidity::InitialLocationOffset,x,3*sizeof(double));
          }

        }  //all particle are allocated for the currect iteration

      } //loop throught he sphere
    } // generate new particles if needed


    //preform the next iteration
    amps_time_step();

    //search for particles that left geospace 
    if (Earth::GeospaceFlag::offset!=-1) {
      //loop through all cells  
      for (int iLocalNode=0;iLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;iLocalNode++) {
        cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::DomainBlockDecomposition::BlockTable[iLocalNode]; 
        PIC::Mesh::cDataBlockAMR *block;
        PIC::Mesh::cDataCenterNode *CenterNode;
        long int ParticleList,ptr;
        int i,j,k,nd;
        char *offset;

        block=node->block;
        if ((block=node->block)==NULL) continue;

        for (k=0;k<_BLOCK_CELLS_Z_;k++) {
          for (j=0;j<_BLOCK_CELLS_Y_;j++) {
            for (i=0;i<_BLOCK_CELLS_X_;i++) {
              nd=PIC::Mesh::mesh->getCenterNodeLocalNumber(i,j,k);

              if ((CenterNode=node->block->GetCenterNode(nd))==NULL) continue;
              offset=CenterNode->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::CenterNodeAssociatedDataOffsetBegin+PIC::CPLR::DATAFILE::MULTIFILE::CurrDataFileOffset;

              if (*((double*)(offset+Earth::GeospaceFlag::offset))==0.0) {
                //the cell is outside of the geospace: search for particles and removed those found 
                 
                ParticleList=block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

                while (ParticleList!=-1) { 
                  ptr=ParticleList;
                  ParticleList=PIC::ParticleBuffer::GetNext(ParticleList);

                  //account for the particle in calculation of the rigidity
                  double x[3],v[3]; 

                  PIC::ParticleBuffer::GetV(v,ptr);
                  PIC::ParticleBuffer::GetX(x,ptr);

                  Earth::CutoffRigidity::ProcessOutsideDomainParticles(ptr,x,v,-1,node);
                  PIC::ParticleBuffer::DeleteParticle(ptr);
                }

                block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]=-1; 
              }
            }
          }
        }

        //search for particles that spend too much time in the simulation 
        ShortTrajectoryFound=false;

        for (k=0;k<_BLOCK_CELLS_Z_;k++) {
          for (j=0;j<_BLOCK_CELLS_Y_;j++) {
            for (i=0;i<_BLOCK_CELLS_X_;i++) {
              ParticleList=block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

              while (ParticleList!=-1) {
                ptr=ParticleList;
                ParticleList=PIC::ParticleBuffer::GetNext(ParticleList);

                double t,l,v[3];
                PIC::ParticleBuffer::byte *ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);

                t=*((double*)(ParticleData+Earth::CutoffRigidity::IntegratedTimeOffset));
                PIC::ParticleBuffer::GetX(x,ptr);

                if (t*t*Vector3D::DotProduct(v,v)>Earth::CutoffRigidity::MaxIntegrationLength*Earth::CutoffRigidity::MaxIntegrationLength) {
                  PIC::ParticleBuffer::DeleteParticle(ptr,block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]); 
                } 
                else {
                 *((double*)(ParticleData+Earth::CutoffRigidity::IntegratedTimeOffset))+=PIC::ParticleWeightTimeStep::GlobalTimeStep[PIC::ParticleBuffer::GetI(ptr)];
                 if (Earth::CutoffRigidity::SearchShortTrajectory==true) ShortTrajectoryFound=true;
              }
            }
         }
       }
     } 
   }
}

    static int LastDataOutputFileNumber=-1;

    if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
      PIC::RequiredSampleLength*=2;
      if (PIC::RequiredSampleLength>50000) PIC::RequiredSampleLength=50000;


      LastDataOutputFileNumber=PIC::DataOutputFileNumber;
      if (PIC::Mesh::mesh->ThisThread==0) cout << "The new sample length is " << PIC::RequiredSampleLength << endl;
    }


    //get the total number of particles in the system
    LacalParticleNumber=PIC::ParticleBuffer::GetAllPartNum();
    MPI_Allreduce(&LacalParticleNumber,&GlobalParticleNumber,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

    //determine whether any particle has been generated during the current iteration
    MPI_Allreduce(&localParticleGenerationFlag,&globalParticleGenerationFlag,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

    if (globalParticleGenerationFlag!=0) {
      //at least one particle has been generated -> reset the iteration counter
      IterationCounter=0;
    }
    else {
      //increment the iteration counter
      IterationCounter++;
    }

    ShortTrajectoryIndex=(ShortTrajectoryFound==true) ? 1 : 0;
    MPI_Allreduce(&ShortTrajectoryIndex,&ShortTrajectoryIndexAll,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

    if (PIC::ThisThread==0) cout << "IterationCounter=" << IterationCounter << ", nMaxIterations=" << nMaxIterations << ", GlobalParticleNumber=" << GlobalParticleNumber << ", ShortTrajectoryIndexAll=" << ShortTrajectoryIndexAll << endl;

  }
  while (((GlobalParticleNumber!=0)&&(IterationCounter<nMaxIterations))||(ShortTrajectoryIndexAll>0));


  //calculate the flux and energey spectrum
  /*
#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#pragma omp parallel for schedule(dynamic,1) default (none)  \
    private (iSphereIndex,spec,iface,dSurface,iTable,jTable,iByte,iBit,Index,Speed,KineticEnergy,DiffFlux,iE,norm,v) \
    shared (iSphereIndexMin,iSphereIndexMax,Earth::CutoffRigidity::DomainBoundaryParticleProperty::dX,Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleMaskNumberPerSpatialDirection, \
        Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleTable,TotalFlux,EnergySpectrum)
#endif
   */



  TotalFlux=0.0;
  EnergySpectrum=0.0;


  if (Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleDomainBoundaryParticleProperty==true) for (int iZenithElement=0;iZenithElement<nZenithElements;iZenithElement++) for (int iAzimutalElement=0;iAzimutalElement<nAzimuthalElements;iAzimutalElement++) {
    iLocation=Sphere.GetLocalSurfaceElementNumber(iZenithElement,iAzimutalElement);

    for (spec=0;spec<PIC::nTotalSpecies;spec++) {
      for (iface=0;iface<6;iface++) {
        //surface area of the element of the surface mesh that covers the boundary of the computational domain
        double lx,ly,lz;

        switch (iface) {
        case 0:case 1:
          ly=(PIC::Mesh::mesh->xGlobalMax[1]-PIC::Mesh::mesh->xGlobalMin[1])/Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleMaskNumberPerSpatialDirection;
          lz=(PIC::Mesh::mesh->xGlobalMax[2]-PIC::Mesh::mesh->xGlobalMin[2])/Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleMaskNumberPerSpatialDirection;

          dSurface=ly*lz;
          break;
        case 2:case 3:
          lx=(PIC::Mesh::mesh->xGlobalMax[0]-PIC::Mesh::mesh->xGlobalMin[0])/Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleMaskNumberPerSpatialDirection;
          lz=(PIC::Mesh::mesh->xGlobalMax[2]-PIC::Mesh::mesh->xGlobalMin[2])/Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleMaskNumberPerSpatialDirection;

          dSurface=lx*lz;
          break;
        case 4:case 5:
          lx=(PIC::Mesh::mesh->xGlobalMax[0]-PIC::Mesh::mesh->xGlobalMin[0])/Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleMaskNumberPerSpatialDirection;
          ly=(PIC::Mesh::mesh->xGlobalMax[1]-PIC::Mesh::mesh->xGlobalMin[1])/Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleMaskNumberPerSpatialDirection;

          dSurface=ly*lz;
        }


        //loop through the mesh that covers face 'iface' on the computational domain
        for (iTable=0;iTable<Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleMaskNumberPerSpatialDirection;iTable++) {
          for (jTable=0;jTable<Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleMaskNumberPerSpatialDirection;jTable++) {
            Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleTable(spec,iLocation,iface,iTable,jTable).Gather();

            for (iByte=0;iByte<Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleTable(spec,iLocation,iface,iTable,jTable).FlagTableLength[0];iByte++) for (iBit=0;iBit<8;iBit++) {
              Index=iBit+8*iByte;

              if (Earth::CutoffRigidity::DomainBoundaryParticleProperty::SampleTable(spec,iLocation,iface,iTable,jTable).Test(Index)==true) {
                //at least one particle that corrsponds to 'Index' has been detected. Add a contribution of such particles to the total energy spectrum and flux as observed at the point of the observation 'iTestsLocation'

                Earth::CutoffRigidity::DomainBoundaryParticleProperty::ConvertVelocityVectorIndex2Velocity(spec,v,iface,Index);
                Speed=Vector3D::Length(v);
                KineticEnergy=Relativistic::Speed2E(Speed,PIC::MolecularData::GetMass(spec));

                //probability density that particles has velocity 'v'
                DiffFlux=GCR_BADAVI2011ASR::Hydrogen::GetDiffFlux(KineticEnergy);

                //determine the contributio of the particles into the 'observed' flux and energy spectrum
                TotalFlux(iLocation,spec)+=DiffFlux*dSurface;

                //determine contribution of the particles to the energy flux
                iE=(log10(KineticEnergy)-logMinEnergyLimit)/dE;
                if (iE<0) iE=0;
                if (iE>=nTotalEnergySpectrumIntervals) iE=nTotalEnergySpectrumIntervals-1;

                EnergySpectrum(iLocation,spec,iE)+=DiffFlux*dSurface;
              }
            }
          }
        }
      }


    }
  }


  //output the calculated map
  CMPI_channel pipe(1000000);
  FILE *fout2d_total_flux,*fout2d_rigidity,**fout2d_spectrum;

  if (PIC::ThisThread==0) {
    //sampled energy spectrum
    int spec,iTestsLocation,iE;
    char fname[400];

    pipe.openRecvAll();

    fout2d_spectrum=new FILE* [PIC::nTotalSpecies];

    //open the output file
    sprintf(fname,"%s/CutoffRigidityMap[R=%e].dat",PIC::OutputDataFileDirectory,Radius);
    fout2d_rigidity=fopen(fname,"w");
    fprintf(fout2d_rigidity,"VARIABLES=\"Lon\", \"Lat\"");

    sprintf(fname,"%s/TotalFluxMap[R=%e].dat",PIC::OutputDataFileDirectory,Radius);
    fout2d_total_flux=fopen(fname,"w");
    fprintf(fout2d_total_flux,"VARIABLES=\"Lon\", \"Lat\"");

    for (spec=0;spec<PIC::nTotalSpecies;spec++) {
      fprintf(fout2d_rigidity,",  \"Cutoff Rigidity [GV] (s=%i)\", \"Injected Particle Number (s=%i)\", \"Max energy injected particles (s=%i)\"",spec,spec,spec);
      fprintf(fout2d_total_flux,",  \"Total Flux [1/(s*m^2)] (s=%i)\"",spec);
    }

    fprintf(fout2d_rigidity,"\nZONE I=%i, J=%i, DATAPACKING=POINT\n",nAzimuthalElements+1,nZenithElements+1);
    fprintf(fout2d_total_flux,"\nZONE I=%i, J=%i, DATAPACKING=POINT\n",nAzimuthalElements+1,nZenithElements+1);

    //files to output the energy spectrum
    for (spec=0;spec<PIC::nTotalSpecies;spec++) {
      sprintf(fname,"%s/EnergySpectrumMap[R=%e].s=%i.dat",PIC::OutputDataFileDirectory,Radius,spec);
      fout2d_spectrum[spec]=fopen(fname,"w");
      fprintf(fout2d_spectrum[spec],"VARIABLES=\"Lon\", \"Lat\"");

      for (iE=0;iE<nTotalEnergySpectrumIntervals;iE++) {
        fprintf(fout2d_spectrum[spec],",  \"(%e[MeV]<E<%e[MeV])\"",pow(10.0,iE*dE+logMinEnergyLimit)*J2MeV,pow(10.0,(iE+1)*dE+logMinEnergyLimit)*J2MeV);
      }

      fprintf(fout2d_spectrum[spec],"\nZONE I=%i, J=%i, DATAPACKING=POINT\n",nAzimuthalElements+1,nZenithElements+1);
    }
  }
  else {
    pipe.openSend(0);
  }

  //interpolate and print the state vector
  long int InterpolationList[nZenithElements*nAzimuthalElements],InterpolationListLength=0;
  int AzimuthalShift=nAzimuthalElements/2;

  for (int iZenith=0;iZenith<nZenithElements+1;iZenith++) for (int iAzimuthalIn=0;iAzimuthalIn<nAzimuthalElements+1;iAzimuthalIn++) { 
//for (int iAzimuthal=0;iAzimuthal<nAzimuthalElements;iAzimuthal++) {

    int iAzimuthal=iAzimuthalIn+AzimuthalShift;
    bool shft_flag=true;

    if (iAzimuthal>=nAzimuthalElements) {
      iAzimuthal-=nAzimuthalElements;
      shft_flag=false;
    }

    if (iAzimuthalIn==nAzimuthalElements) {
      shft_flag=false;
      iAzimuthal=+AzimuthalShift;
    } 
   
    Sphere.GetSurfaceCoordinate(x,iZenith,iAzimuthal);

    if (PIC::ThisThread==0) {
      double lon,lat;
      Sphere.GetSurfaceLonLatNormal(lon,lat,iZenith,iAzimuthal);

      if (shft_flag==true) lon-=360.0;

      fprintf(fout2d_rigidity,"%e %e ",lon,lat);
      fprintf(fout2d_total_flux,"%e %e ",lon,lat);

      for (int spec=0;spec<PIC::nTotalSpecies;spec++) fprintf(fout2d_spectrum[spec],"%e %e ",lon,lat);
    }


    //prepare the interpolation stencil
    InterpolationListLength=0;

    if (iZenith==0) {
      InterpolationList[InterpolationListLength++]=Sphere.GetLocalSurfaceElementNumber(0,iAzimuthal);
      InterpolationList[InterpolationListLength++]=Sphere.GetLocalSurfaceElementNumber(0,((iAzimuthal>0) ? iAzimuthal-1 : nAzimuthalElements-1));
    }
    else if (iZenith==nZenithElements) {
      InterpolationList[InterpolationListLength++]=Sphere.GetLocalSurfaceElementNumber(nZenithElements-1,iAzimuthal);
      InterpolationList[InterpolationListLength++]=Sphere.GetLocalSurfaceElementNumber(nZenithElements-1,((iAzimuthal>0) ? iAzimuthal-1 : nAzimuthalElements-1));
    }
    else {
      int iA,iZ,A[2],Z[2];

      Z[0]=iZenith-1,Z[1]=iZenith;

      A[0]=(iAzimuthal!=0) ? iAzimuthal-1 : nAzimuthalElements-1;
      A[1]=iAzimuthal;

      for (iA=0;iA<2;iA++) for (iZ=0;iZ<2;iZ++) InterpolationList[InterpolationListLength++]=Sphere.GetLocalSurfaceElementNumber(Z[iZ],A[iA]);
    }

    //prepare and print the interpolated value of the cutoff rigidity
    //loop throught all species
    for (spec=0;spec<PIC::nTotalSpecies;spec++) {
      //loop therough elements of the interpolation stencil
      double minElementRigidity=Earth::CutoffRigidity::CutoffRigidityTable(spec,InterpolationList[0]);
      double norm=0.0,flux_total=0.0,sum_area=0.0;
      array_1d<double> LocalEnergySpectrum(nTotalEnergySpectrumIntervals);

      double interpolated_number_injected_partilces=0.0;
      int number_injected_partilces=0;
      double max_injected_particle_energy=0.0;

      LocalEnergySpectrum=0.0;

      for (int el=0;el<InterpolationListLength;el++) {
        //loop through all MPI processes

        if (PIC::ThisThread==0) {
          //this process will output data
          double t;

          //the totalsurface area
          sum_area+=Sphere.GetSurfaceElementArea(InterpolationList[el]);
          flux_total+=TotalFlux(InterpolationList[el],spec);
          number_injected_partilces=Earth::CutoffRigidity::InjectedParticleMap(spec,InterpolationList[el]);

          max_injected_particle_energy=Earth::CutoffRigidity::MaxEnergyInjectedParticles(spec,InterpolationList[el]);

          //collect cutoff regidiry from other MPI processes
          for (int thread=1;thread<PIC::nTotalThreads;thread++) {
            number_injected_partilces+=pipe.recv<int>(thread);

            pipe.recv(t,thread);
            if (t>max_injected_particle_energy) max_injected_particle_energy=t;

            pipe.recv(t,thread);
            if ( (t>0.0) && ((minElementRigidity<0.0)||(t<minElementRigidity)) )  minElementRigidity=t;
          }

          interpolated_number_injected_partilces+=number_injected_partilces*Sphere.GetSurfaceElementArea(InterpolationList[el]);

          //recieve the energy specgtrum
          for (iE=0,norm=0.0;iE<nTotalEnergySpectrumIntervals;iE++)  LocalEnergySpectrum(iE)+=EnergySpectrum(InterpolationList[el],spec,iE);

          for (int thread=1;thread<PIC::nTotalThreads;thread++) for (iE=0,norm=0.0;iE<nTotalEnergySpectrumIntervals;iE++) {
            pipe.recv(t,thread);

            LocalEnergySpectrum(iE)+=t;
          }
        }
        else {
          //send the number of injected model particles 
          pipe.send(Earth::CutoffRigidity::InjectedParticleMap(spec,InterpolationList[el]));

          //send the max energy of the injected particles 
          pipe.send(Earth::CutoffRigidity::MaxEnergyInjectedParticles(spec,InterpolationList[el]));

          //send cutoff rigidity to the "root" MPI process
          pipe.send(Earth::CutoffRigidity::CutoffRigidityTable(spec,InterpolationList[el]));

          //send ebergy spectrum to the root
          for (iE=0,norm=0.0;iE<nTotalEnergySpectrumIntervals;iE++) pipe.send(EnergySpectrum(InterpolationList[el],spec,iE));

        }
      }

      //the interpolated value is calcualted by the root MPI process
      if (PIC::ThisThread==0) {

        //normalize the energy spectrum
        for (iE=0,norm=0.0;iE<nTotalEnergySpectrumIntervals;iE++) norm+=LocalEnergySpectrum(iE)*(pow((iE+1)*dE+logMinEnergyLimit,10)-pow(iE*dE+logMinEnergyLimit,10));
        if (norm>0.0) for (iE=0;iE<nTotalEnergySpectrumIntervals;iE++) LocalEnergySpectrum(iE)/=norm;


        fprintf(fout2d_rigidity,"  %e  %e  %e",minElementRigidity,interpolated_number_injected_partilces/sum_area,max_injected_particle_energy);
        fprintf(fout2d_total_flux,"  %e",flux_total/sum_area);

        for (iE=0;iE<nTotalEnergySpectrumIntervals;iE++) fprintf(fout2d_spectrum[spec],"  %e",LocalEnergySpectrum(iE));
      }
    }

    //the cutoff rigidity is computed for all speces at the given point in the map
    if (PIC::ThisThread==0) {
      fprintf(fout2d_rigidity,"\n");
      fprintf(fout2d_total_flux,"\n");

      for (spec=0;spec<PIC::nTotalSpecies;spec++) fprintf(fout2d_spectrum[spec],"\n");
    }
  }

  //close the pipe nad the file
  if (ThisThread==0) {
    fclose(fout2d_rigidity);
    fclose(fout2d_total_flux);

    for (spec=0;spec<PIC::nTotalSpecies;spec++) fclose(fout2d_spectrum[spec]);

    pipe.closeRecvAll();
  }
  else pipe.closeSend();



  //delete all particles that still present in the system
  PIC::ParticleBuffer::DeleteAllParticles();
//  Earth::CutoffRigidity::CutoffRigidityTable.Deallocate();
  Earth::CutoffRigidity::CutoffRigidityTable=0.0;
}


void CutoffRigidityCalculation(int nMaxIterations) {
  //disable sampling
  PIC::Sampling::RuntimeSamplingSwitch=false;

  //estimate the total flux and rigidity at a sphere

  if (Earth::RigidityCalculationMode==Earth::_sphere) {
    if (Earth::RigidityCalculationSphereRadius==0.0) {
      exit(__LINE__,__FILE__,"The radius was not set");
    }
    else {
      SampleSphericalMaplLocations(Earth::RigidityCalculationSphereRadius,nMaxIterations);
    }
  }
  else {
    SampleIndividualLocations(nMaxIterations);
  }


  /*
  //start forward integration
  //enable sampling
  PIC::Sampling::RuntimeSamplingSwitch=true;
  PIC::ParticleBuffer::DeleteAllParticles();
  Earth::ForwardParticleModeling(nMaxIterations);

  //estimate the cutoff rigidity and energy spectrum in individual locations
  PIC::ParticleBuffer::DeleteAllParticles();
  SampleIndividualLocations(nMaxIterations);
  */
}

void CutoffRigidityCalculation_Legacy(int nTotalIterations) {
  int LastDataOutputFileNumber=PIC::DataOutputFileNumber;

  if (Earth::CutoffRigidity::DomainBoundaryParticleProperty::SamplingParameters::ActiveFlag==true) {
    PIC::Mover::BackwardTimeIntegrationMode=_PIC_MODE_ON_;
    Earth::CutoffRigidity::DomainBoundaryParticleProperty::EnableSampleParticleProperty=true;

    //particles will be injected only in the near Earth's region
    Earth::BoundingBoxInjection::BoundaryInjectionMode=false;
    Earth::CutoffRigidity::ParticleInjector::ParticleInjectionMode=true;

    for (long int niter=0;(niter<nTotalIterations)&&(LastDataOutputFileNumber<Earth::CutoffRigidity::DomainBoundaryParticleProperty::SamplingParameters::LastActiveOutputCycleNumber);niter++) {
      amps_time_step();

      if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
        PIC::RequiredSampleLength*=2;
        if (PIC::RequiredSampleLength>50000) PIC::RequiredSampleLength=50000;


        LastDataOutputFileNumber=PIC::DataOutputFileNumber;
        if (PIC::Mesh::mesh->ThisThread==0) cout << "The new sample length is " << PIC::RequiredSampleLength << endl;
      }
    }

    Earth::CutoffRigidity::DomainBoundaryParticleProperty::Gather();
    Earth::CutoffRigidity::DomainBoundaryParticleProperty::SmoothSampleTable();
    PIC::Mover::BackwardTimeIntegrationMode=_PIC_MODE_OFF_;

    //partices will be injected from the boundary of the domain
    Earth::BoundingBoxInjection::BoundaryInjectionMode=true;
    Earth::CutoffRigidity::ParticleInjector::ParticleInjectionMode=false;
    Earth::CutoffRigidity::DomainBoundaryParticleProperty::ApplyInjectionPhaseSpaceLimiting=true;
    Earth::CutoffRigidity::DomainBoundaryParticleProperty::EnableSampleParticleProperty=false;
  }

  //time step with the forward integration
  for (long int niter=0;niter<nTotalIterations;niter++) {
    amps_time_step();
    
    if (PIC::Mesh::mesh->ThisThread==0) {
      time_t TimeValue=time(NULL);
      tm *ct=localtime(&TimeValue);
      printf(": (%i/%i %i:%i:%i), Iteration: %ld  (current sample length:%ld, %ld interations to the next output)\n",
       ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec,niter,
       PIC::RequiredSampleLength,
       PIC::RequiredSampleLength-PIC::CollectingSampleCounter);
    }

     if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
       PIC::RequiredSampleLength*=2;
       if (PIC::RequiredSampleLength>50000) PIC::RequiredSampleLength=50000;


       LastDataOutputFileNumber=PIC::DataOutputFileNumber;
       if (PIC::Mesh::mesh->ThisThread==0) cout << "The new sample length is " << PIC::RequiredSampleLength << endl;
     }
  }
}


int main(int argc,char **argv) {
  static int LastDataOutputFileNumber=0;

  Earth::CutoffRigidity::SampleRigidityMode=true;

  //read post-compile input file  
  if (PIC::PostCompileInputFileName!="") {
    Earth::Parser::ReadFile(PIC::PostCompileInputFileName,Earth::Parser::_reading_mode_pre_init);
  } 

  amps_init_mesh();

  if (PIC::PostCompileInputFileName!="") {
    Earth::Parser::ReadFile(PIC::PostCompileInputFileName,Earth::Parser::_reading_mode_pre_init);
  }

  Earth::CutoffRigidity::Init_BeforeParser();
  Earth::CutoffRigidity::AllocateCutoffRigidityTable();

  amps_init();

/*
  //read post-compile input file
  if (PIC::PostCompileInputFileName!="") {
    Earth::Parser::ReadFile(PIC::PostCompileInputFileName,Earth::Parser::_reading_mode_post_init);
  }
*/



  if (Earth::ModelMode==Earth::BoundaryInjectionMode) {
    int nIterations,nTotalIterations=100000001;
    double et;

    if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_)  nTotalIterations=3200;

    //load parameters of T05 model
    if (Earth::BackgroundMagneticFieldT05Data!="") {
      T05::LoadDataFile(Earth::BackgroundMagneticFieldT05Data.c_str()); 
    }    

    //init the simulation start time 
    str2et_c(Exosphere::SimulationStartTimeString,&et); 

    int iEtInterval=-1;

    for (int i=0;i<T05::Data.size()-1;i++) {  
      if ((T05::Data[i].et<=et)&&(et<=T05::Data[i+1].et)) {
        iEtInterval=i;
        break;
      }
    }

    if (iEtInterval==-1) exit(__LINE__,__FILE__,"Error: the time intervais was not found");

    T05::SetSolarWindPressure_nano(T05::Data[iEtInterval].Pdyn);
    T05::SetDST_nano(T05::Data[iEtInterval].SYMH);
    T05::SetBYIMF_nano(T05::Data[iEtInterval].BYGSM);
    T05::SetBZIMF_nano(T05::Data[iEtInterval].BZGSM);

    T05::SetW1(T05::Data[iEtInterval].W1);
    T05::SetW2(T05::Data[iEtInterval].W2);
    T05::SetW3(T05::Data[iEtInterval].W3);
    T05::SetW4(T05::Data[iEtInterval].W4);
    T05::SetW5(T05::Data[iEtInterval].W5);
    T05::SetW6(T05::Data[iEtInterval].W6);

    Earth::InitMagneticField();

  //read post-compile input file
  if (PIC::PostCompileInputFileName!="") {
    Earth::Parser::ReadFile(PIC::PostCompileInputFileName,Earth::Parser::_reading_mode_post_init);
  }


    //time step
    for (long int niter=0;niter<nTotalIterations;niter++) {
      static int LastDataOutputFileNumber=-1;

      PIC::TimeStep();

      et+=PIC::ParticleWeightTimeStep::GlobalTimeStep[0];
      int iEtInterval_prev=iEtInterval;
      
      for (;iEtInterval<T05::Data.size()-1;iEtInterval++) {  
        if ((T05::Data[iEtInterval].et<=et)&&(et<=T05::Data[iEtInterval+1].et)) {
          break;
        }
      }

      if (iEtInterval_prev!=iEtInterval) {
        //update magnetic field model
        T05::SetSolarWindPressure_nano(T05::Data[iEtInterval].Pdyn);
        T05::SetDST_nano(T05::Data[iEtInterval].SYMH);
        T05::SetBYIMF_nano(T05::Data[iEtInterval].BYGSM);
        T05::SetBZIMF_nano(T05::Data[iEtInterval].BZGSM);
        T05::SetIMF_nano(T05::Data[iEtInterval].BXGSM,T05::Data[iEtInterval].BYGSM,T05::Data[iEtInterval].BZGSM);
    
        T05::SetW1(T05::Data[iEtInterval].W1);
        T05::SetW2(T05::Data[iEtInterval].W2);
        T05::SetW3(T05::Data[iEtInterval].W3);
        T05::SetW4(T05::Data[iEtInterval].W4);
        T05::SetW5(T05::Data[iEtInterval].W5);
        T05::SetW6(T05::Data[iEtInterval].W6);

        Earth::InitMagneticField();
      }


      if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
        PIC::RequiredSampleLength*=2;
        if (PIC::RequiredSampleLength>1000) PIC::RequiredSampleLength=1000;


        LastDataOutputFileNumber=PIC::DataOutputFileNumber;
        if (PIC::Mesh::mesh->ThisThread==0) cout << "The new sample length is " << PIC::RequiredSampleLength << endl;
      }

      if (PIC::Mesh::mesh->ThisThread==0) {
        time_t TimeValue=time(NULL);
        tm *ct=localtime(&TimeValue);

        printf(": (%i/%i %i:%i:%i), Iteration: %ld  (current sample length:%ld, %ld interations to the next output)\n",ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec,niter,PIC::RequiredSampleLength,PIC::RequiredSampleLength-PIC::CollectingSampleCounter);
      }
    }
  }
  else if (Earth::ModelMode==Earth::CutoffRigidityMode) {
    if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_) {
      //execute the nightly test routine

      int niter=(Earth::CutoffRigidity::nTotalIterations!=-1) ? Earth::CutoffRigidity::nTotalIterations : 10000;

      switch (_NIGHTLY_TEST_) {
      case _NIGHTLY_TEST__CUTOFF_:
        if (PIC::ThisThread==0) cout << "_NIGHTLY_TEST_=_NIGHTLY_TEST__CUTOFF_" << endl;

        CutoffRigidityCalculation(niter);
        break;
      case _NIGHTLY_TEST__LEGACY_:
        if (PIC::ThisThread==0) cout << "_NIGHTLY_TEST_=_NIGHTLY_TEST__LEGACY_" << endl;

        CutoffRigidityCalculation_Legacy(0);
        break;
       }

      //output the particle statistics of the test run
      char fname[300];
      sprintf(fname,"%s/test_Earth.dat",PIC::OutputDataFileDirectory);
      PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);
    }
    else {
      int niter=(Earth::CutoffRigidity::nTotalIterations!=-1) ? Earth::CutoffRigidity::nTotalIterations : 10000;

      CutoffRigidityCalculation(niter);
    }
  }
  else exit(__LINE__,__FILE__,"Error: the option is not recognized");


  //finish the run
  MPI_Finalize();
  cout << "End of the run:" << PIC::nTotalSpecies << endl;
  
  return EXIT_SUCCESS;
}
