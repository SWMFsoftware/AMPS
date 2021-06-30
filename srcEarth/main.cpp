
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


#include <sys/time.h>
#include <sys/resource.h>

#include "T05Interface.h"
#include "T96Interface.h"

//$Id$



//the particle class
#include "pic.h"
#include "constants.h"
#include "Earth.h"

extern int nZenithElements;
extern int nAzimuthalElements;

void amps_init();
void amps_init_mesh();
void amps_time_step();

#define _NIGHTLY_TEST__CUTOFF_  0
#define _NIGHTLY_TEST__LEGACY_ 1

#ifndef _NIGHTLY_TEST_
#define _NIGHTLY_TEST_ _NIGHTLY_TEST__CUTOFF_
#endif


void SampleIndividualLocations(int nMaxIterations) {
  int IterationCounter=0,localParticleGenerationFlag=0,globalParticleGenerationFlag;

  //estimate the total flux and rigidity in a set of the defined locations
  if (Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength!=0) {
    int LacalParticleNumber,GlobalParticleNumber;
    int nIngectedParticlePerIteration=Earth::CutoffRigidity::IndividualLocations::nTotalTestParticlesPerLocations/Earth::CutoffRigidity::IndividualLocations::nParticleInjectionIterations;
    int nTotalInjectedParticles=0;

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

    do {
      //reset the partilce generation flag
      localParticleGenerationFlag=0;

      //Inject new particles
      if (nTotalInjectedParticles<Earth::CutoffRigidity::IndividualLocations::nTotalTestParticlesPerLocations) {
        nTotalInjectedParticles+=nIngectedParticlePerIteration;

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

            if (PIC::Mesh::mesh->fingCellIndex(x,iCell,jCell,kCell,startNode,false)==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located");

            for (int iNewParticle=0;iNewParticle<nIngectedParticlePerIteration;iNewParticle++) {
              energy=exp(logMinEnergyLimit+rnd()*(logMaxEnergyLimit-logMinEnergyLimit));

              speed=Relativistic::E2Speed(energy,mass);
              Vector3D::Distribution::Uniform(v,speed);

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

              //set up the particle rigidity
              if (Earth::CutoffRigidity::InitialRigidityOffset!=-1) {
                double momentum,charge,rigidity;

                charge=fabs(PIC::MolecularData::GetElectricCharge(spec));

                momentum=Relativistic::Speed2Momentum(speed,mass);
                rigidity=(charge>0.0) ? momentum/charge : 0.0;

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


/*      //limit the integration time (remove particle moving in trapped orbits)
      if ((true)&&(Earth::CutoffRigidity::InitialLocationOffset!=-1)) {
        //determine min altitude;
        static double rMax=-1.0;
        double r;
        int idim;

        if (rMax<0.0) {
          for (int iLocation=0;iLocation<Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength;iLocation++) {
            r=Vector3D::Length(Earth::CutoffRigidity::IndividualLocations::xTestLocationTable[iLocation]);

            if ((rMax<0.0)||(r>rMax)) rMax=r;
          }
        }

        //estimation of the domain length
        double IntegratioinLengthMax=sqrt(pow(PIC::Mesh::mesh->xGlobalMax[0]-PIC::Mesh::mesh->xGlobalMin[0],2)+
            pow(PIC::Mesh::mesh->xGlobalMax[1]-PIC::Mesh::mesh->xGlobalMin[1],2)+
            pow(PIC::Mesh::mesh->xGlobalMax[2]-PIC::Mesh::mesh->xGlobalMin[2],2));

        //increase the limit of the total integrated oath length
        IntegratioinLengthMax/=4.0;

        IntegratioinLengthMax=2.0*Pi*rMax;

        //loop through all blocks
        for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh->BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
          int ptr,next;
          PIC::Mesh::cDataBlockAMR *block=node->block;
          double dt,l,v;
          PIC::ParticleBuffer::byte *ParticleData;

          if (sqrt(pow(node->xmax[0]+node->xmin[0],2)+pow(node->xmax[1]+node->xmin[1],2)+pow(node->xmax[2]+node->xmin[2],2))/2.0>rMax) continue;

          if (block!=NULL) for (int i=0;i<_BLOCK_CELLS_X_;i++) for (int j=0;j<_BLOCK_CELLS_Y_;j++) for (int k=0;k<_BLOCK_CELLS_Z_;k++) {
             ptr=block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

             while (ptr!=-1) {
               //loop through all particles and determin which are trapped and does not contribute to the incoming radiation
               ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);

               v=Vector3D::Length(PIC::ParticleBuffer::GetV(ParticleData));
               next=PIC::ParticleBuffer::GetNext(ParticleData);

               dt=block->GetLocalTimeStep(PIC::ParticleBuffer::GetI(ParticleData));

               //update estimation of the total integration path length
               l=*((double*)(ParticleData+Earth::CutoffRigidity::IntegratedPathLengthOffset));
               l+=v*dt;
               *((double*)(ParticleData+Earth::CutoffRigidity::IntegratedPathLengthOffset))=l;


               if (l>IntegratioinLengthMax) {
                 //the particle is below the minimum altitude
                 PIC::ParticleBuffer::DeleteParticle(ptr,block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]);
               }

               ptr=next;
             }

          }
        }
      }*/

      //preform the next iteration
      amps_time_step();

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

      if (PIC::ThisThread==0) {
        cout << "IterationCounter=" << IterationCounter << ", GlobalParticleNumber=" << GlobalParticleNumber << endl;
      }

    }
    while ((GlobalParticleNumber!=0)&&(IterationCounter<nMaxIterations));

    //delete all particles that still present in the system
    PIC::ParticleBuffer::DeleteAllParticles();

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
      exit(__LINE__,__FILE__,"Error: not implemented");
    }

  }

  //release sampling buffers
 // Earth::CutoffRigidity::DomainBoundaryParticleProperty::Deallocate();
}



void SampleSphericalMaplLocations(double Radius,int nMaxIterations) {
  int IterationCounter=0,localParticleGenerationFlag=0,globalParticleGenerationFlag;
  int iLocation;
  double x[3]={0.0,0.0,0.0},v[3];

  PIC::Mover::BackwardTimeIntegrationMode=_PIC_MODE_ON_;


  int nTotalInjectedParticlePerPoint=25;

  Earth::CutoffRigidity::IndividualLocations::nTotalTestParticlesPerLocations=2000;

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

          if ((nd=PIC::Mesh::mesh->fingCellIndex(x,iCell,jCell,kCell,startNode,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located");

#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
          PIC::Mesh::cDataCenterNode *cell;

          cell=startNode->block->GetCenterNode(nd);

          if (cell==NULL) exit(__LINE__,__FILE__,"Error: the cell measure is not initialized");
          if (cell->Measure<=0.0) exit(__LINE__,__FILE__,"Error: the cell measure is not initialized");
#endif

          //generate energy of the new particle
          energy=exp(logMinEnergyLimit+rnd()*(logMaxEnergyLimit-logMinEnergyLimit));

          speed=Relativistic::E2Speed(energy,mass);
          Vector3D::Distribution::Uniform(v,speed);

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

    if (PIC::ThisThread==0) cout << "IterationCounter=" << IterationCounter << ", nMaxIterations=" << nMaxIterations << ", GlobalParticleNumber=" << GlobalParticleNumber << endl;

  }
  while ((GlobalParticleNumber!=0)&&(IterationCounter<nMaxIterations));


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
  SampleSphericalMaplLocations(_EARTH__RADIUS_+500.0E3,nMaxIterations);

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


  //TEST T96/T05
  if ((_PIC_COUPLER_MODE_==_PIC_COUPLER_MODE__T05_)||(_PIC_COUPLER_MODE_==_PIC_COUPLER_MODE__T96_)) {
    //SI unit are used
    double x[3]={4.0*_RADIUS_(_EARTH_),1.0E-5*_RADIUS_(_EARTH_),-1.0E-5*_RADIUS_(_EARTH_)};
    double B[3];

    switch (_PIC_COUPLER_MODE_) {
    case _PIC_COUPLER_MODE__T05_: 
      T05::Init(Exosphere::SimulationStartTimeString,NULL);

      T05::SetSolarWindPressure(2.97*_NANO_); 
      T05::SetDST(-10.0*_NANO_);
      T05::SetBYIMF(-16.600000*_NANO_);
      T05::SetBZIMF(-8.200000*_NANO_);
      T05::SetW(1.022,0.874,0.807,1.657000,2.737,4.542000);

      T05::GetMagneticField(B,x);
      break;
    case _PIC_COUPLER_MODE__T96_:
      T96::Init(Exosphere::SimulationStartTimeString,NULL);

      T96::SetSolarWindPressure(Earth::T96::solar_wind_pressure);
      T96::SetDST(Earth::T96::dst);
      T96::SetBYIMF(Earth::T96::by);
      T96::SetBZIMF(Earth::T96::bz);
    } 
  
    if (PIC::ThisThread==0) {
      cout << "T05 test: x=" << x[0] << ", " << x[1] << ", " << x[2] << ",  B=" << B[0] << ", " << B[1] << ",  " << B[2] << endl;
    }
  }


  Earth::CutoffRigidity::SampleRigidityMode=true;

  amps_init_mesh();

  Earth::CutoffRigidity::Init_BeforeParser();
  Earth::CutoffRigidity::AllocateCutoffRigidityTable();

  amps_init();


  if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_) {
    //execute the nightly test routine

    switch (_NIGHTLY_TEST_) {
    case _NIGHTLY_TEST__CUTOFF_:
      if (PIC::ThisThread==0) cout << "_NIGHTLY_TEST_=_NIGHTLY_TEST__CUTOFF_" << endl;

      CutoffRigidityCalculation(10000);
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
    CutoffRigidityCalculation(10000);
  }




  //finish the run
  MPI_Finalize();
  cout << "End of the run:" << PIC::nTotalSpecies << endl;
  
  return EXIT_SUCCESS;
}
