
//$Id$

/*
 * ColumnIntegral.cpp
 *
 *  Created on: Oct 9, 2015
 *      Author: vtenishe
 */

//calcualtion of the coma column density and brightness of the dust


#include "pic.h"
#include "global.h"
#include "Exosphere.h"
#include "Dust.h"

//==============================================================================================================================
//sampling data offset
int Comet::Sampling::SamplingDataOffset=-1;


int Exosphere::ColumnIntegral::GetVariableList(char *vlist) {
  int s,nvars=0;

  if (vlist!=NULL) sprintf(vlist,"");


  for (int iCase=0;iCase<ElectricallyChargedDust::ScatteringEfficientcy::Mie::Fink2012Icarus::nCases;iCase++) {
    nvars++;
    if (vlist!=NULL) sprintf(vlist,"%s, \"Dust Density*a^2*Qmie(%i)/r2\"",vlist,iCase);
  }


  //column density of the loaded neutral profile
  for (s=0;s<Comet::CometData::nNeutrals;s++) {
    nvars++;
    if (vlist!=NULL) sprintf(vlist,"%s, \"Column Density[s=%i]\"",vlist,s);
  }

  //column density of the modeled species
  for (s=0;s<PIC::nTotalSpecies;s++) {
    nvars++;
    if (vlist!=NULL) sprintf(vlist,"%s, \"Modeled Column Density[s=%i]\"",vlist,s);
  }

  return nvars;
}

void Exosphere::ColumnIntegral::CoulumnDensityIntegrant(double *res,int resLength,double* x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  int i,j,k,nd,cnt=0,spec;
  double NumberDensity,wtot,DustBrightness=0.0;
  PIC::Mesh::cDataCenterNode *cell;
  char *CellData;

  PIC::InterpolationRoutines::CellCentered::cStencil* Stencil;

  //deterine the interpolation stencil
  if (_PIC_COUPLER__INTERPOLATION_MODE_ == _PIC_COUPLER__INTERPOLATION_MODE__CELL_CENTERED_LINEAR_) {
    Stencil=PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(x,node);
  }

  //init the return values
  for (i=0;i<resLength;i++) res[i]=0.0;

  nd=PIC::Mesh::mesh->FindCellIndex(x,i,j,k,node);

  if (node->block==NULL) return;
  if ((cell=node->block->GetCenterNode(nd))==NULL) return;
  CellData=cell->GetAssociatedDataBufferPointer();

  for (int iCase=0;iCase<ElectricallyChargedDust::ScatteringEfficientcy::Mie::Fink2012Icarus::nCases;iCase++) {
    DustBrightness=0.0;

    for (spec=0;spec<PIC::nTotalSpecies;spec++) if (_DUST_SPEC_<=spec && spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups) {
      wtot=cell->GetCompleteSampleCellParticleWeight(spec);

      if (wtot>0.0) {
        DustBrightness+=
            ((double*)(Comet::Sampling::SamplingDataOffset+PIC::Mesh::completedCellSampleDataPointerOffset+CellData))
            [iCase+(spec-_DUST_SPEC_)*ElectricallyChargedDust::ScatteringEfficientcy::Mie::Fink2012Icarus::nCases]/
            wtot*cell->GetNumberDensity(spec);
      }

    }

    res[cnt++]=DustBrightness;
  }


  //calculate column density dust to the background species
  for (int s=0;s<Comet::CometData::nNeutrals;s++) {
    if (_PIC_COUPLER__INTERPOLATION_MODE_ == _PIC_COUPLER__INTERPOLATION_MODE__CELL_CENTERED_LINEAR_) {
      for (int iStencil=0;iStencil<Stencil->Length;iStencil++) {
        res[cnt]+=Stencil->Weight[iStencil]*Comet::CometData::GetNeutralsMassDensity(s,Stencil->cell[iStencil]);
      }
    }
    else res[cnt]=Comet::CometData::GetNeutralsMassDensity(s,nd,node);

    ++cnt;
  }

  //calculate the column density of the modeled species
  for (int s=0;s<PIC::nTotalSpecies;s++) {
    if (_PIC_COUPLER__INTERPOLATION_MODE_ == _PIC_COUPLER__INTERPOLATION_MODE__CELL_CENTERED_LINEAR_) {
      for (int iStencil=0;iStencil<Stencil->Length;iStencil++) {
        res[cnt]+=Stencil->Weight[iStencil]*Stencil->cell[iStencil]->GetNumberDensity(s);
      }
    }
    else res[cnt]=node->block->GetCenterNode(nd)->GetNumberDensity(s);

    ++cnt;
  }


  if (cnt!=resLength) exit(__LINE__,__FILE__,"Error: something is wrong with the variable counting (cnt!=resLength)");

}





//==============================================================================================================================
//sampling of the particle data
void Comet::Sampling::SampleParticleData(char *ParticleData,double LocalParticleWeight,char  *SamplingBuffer,int spec) {

  if ((_DUST_SPEC_<=spec) && (spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups)) {
    //the particle is a dust grain
    double GrainRadius;
    int iRadiusGroup;

    iRadiusGroup=spec-_DUST_SPEC_;
    GrainRadius=ElectricallyChargedDust::GetGrainRadius((PIC::ParticleBuffer::byte*)ParticleData);

    for (int iCase=0;iCase<ElectricallyChargedDust::ScatteringEfficientcy::Mie::Fink2012Icarus::nCases;iCase++) {
      ((double*)(SamplingDataOffset+SamplingBuffer))[iCase+iRadiusGroup*ElectricallyChargedDust::ScatteringEfficientcy::Mie::Fink2012Icarus::nCases]+=
          Pi*pow(GrainRadius,2)*ElectricallyChargedDust::ScatteringEfficientcy::Mie::Fink2012Icarus::Get(GrainRadius,1.095E-6,iCase)*LocalParticleWeight;
    }
  }
}

//request data for the particle sampling
int Comet::Sampling::RequestSamplingData(int offset) {
  SamplingDataOffset=offset;

  return ElectricallyChargedDust::ScatteringEfficientcy::Mie::Fink2012Icarus::nCases*ElectricallyChargedDust::GrainVelocityGroup::nGroups*sizeof(double);
}

//==============================================================================================================================
//functions that are passed into the core for sampling of the model data (empty) and output of the column integrals
void Comet::Sampling::SampleModelData() {}

void Comet::Sampling::PrintBrightnessMap(int DataOutputFileNumber) {
  double halfAngleRange,r,theta;
  int iTestPoint,idim,iMap;

  if (_COMET_PRINT_BRIGHTNESS_MAP_MODE_==_PIC_MODE_OFF_) return;

  const int nCalcualtedMaps=4;

  for (iTestPoint=0;iTestPoint<ElectricallyChargedDust::Sampling::FluxMap::SampleLocations.size();iTestPoint++) {
    for (r=0.0,idim=0;idim<3;idim++) r+=pow(ElectricallyChargedDust::Sampling::FluxMap::SampleLocations[iTestPoint].x[idim],2);

    r=sqrt(r);
    theta=atan(2.0*_RADIUS_(_TARGET_)/r);

    for (iMap=0;iMap<nCalcualtedMaps;iMap++) {
      halfAngleRange=theta*pow(2.0,iMap);
      if (halfAngleRange>Pi/2.0) halfAngleRange=Pi/2.0;

      PrintBrightnessMap(halfAngleRange,iTestPoint,DataOutputFileNumber);
    }
  }
}

void Comet::Sampling::PrintBrightnessMap(double halfAngleRange,int iTestPoint,int DataOutputFileNumber) {
  int iZenithPoint,iAzimuthPoint,idim,i;
  double x0[3],e0[3],e1[3],e2[3],l[3],theta,phi;
  double ZenithAngle,AzimuthAngle;
  double cosPhi,sinPhi,sinTheta;
  FILE *fout=NULL;
  char fname[_MAX_STRING_LENGTH_PIC_];


  const int nZenithPoints=100;
  const int nAzimuthPoints=100;

  const double dZenitAngle=2.0*halfAngleRange/(nZenithPoints-1);
  const double dAzimuthAngle=2.0*halfAngleRange/(nAzimuthPoints-1);

  const int StateVectorLength=ElectricallyChargedDust::ScatteringEfficientcy::Mie::Fink2012Icarus::nCases+Comet::CometData::nNeutrals+PIC::nTotalSpecies;
  double StateVector[StateVectorLength];

  //open output file
  if (PIC::ThisThread==0) {
    sprintf(fname,"%s/Comet.ColumnIntegrals.MapAngularRange=%e.SamplePoint=%i.out=%i.dat",PIC::OutputDataFileDirectory,halfAngleRange/Pi*180.0,iTestPoint,DataOutputFileNumber);
    fout=fopen(fname,"w");

    fprintf(fout,"VARIABLES=\"Lon\", \"Lat\", \"Nucleus Projection\" ");

    for (int iCase=0;iCase<ElectricallyChargedDust::ScatteringEfficientcy::Mie::Fink2012Icarus::nCases;iCase++)
      fprintf(fout,", \"Dust Density * a^2 *Qmie(a)[case=%i] Integral\"",iCase);

    for (int s=0;s<Comet::CometData::nNeutrals;s++) fprintf(fout,", \"Column Density[s=%i]\"",s);

    for (int s=0;s<PIC::nTotalSpecies;s++) fprintf(fout,", \"Modeled Column Density[s=%i]\"",s);


    fprintf(fout,"\nZONE I=%ld, J=%ld, DATAPACKING=POINT\n",nAzimuthPoints,nZenithPoints);
  }

  //determine the frame of reference for the map
  for (idim=0;idim<3;idim++) {
    e0[idim]=-ElectricallyChargedDust::Sampling::FluxMap::SampleLocations[iTestPoint].e0[idim];
    e1[idim]=-ElectricallyChargedDust::Sampling::FluxMap::SampleLocations[iTestPoint].e1[idim];
    e2[idim]=-ElectricallyChargedDust::Sampling::FluxMap::SampleLocations[iTestPoint].e2[idim];

    x0[idim]=ElectricallyChargedDust::Sampling::FluxMap::SampleLocations[iTestPoint].x[idim];
  }

  //create the map
  for (iZenithPoint=0;iZenithPoint<nZenithPoints;iZenithPoint++) {
    ZenithAngle=halfAngleRange-dZenitAngle*iZenithPoint;

    for (iAzimuthPoint=0;iAzimuthPoint<nAzimuthPoints;iAzimuthPoint++) {
      AzimuthAngle=-halfAngleRange+dAzimuthAngle*iAzimuthPoint;

      cosPhi=cos(AzimuthAngle);
      sinPhi=sin(AzimuthAngle);
      sinTheta=sin(ZenithAngle);

      for (idim=0;idim<3;idim++) l[idim]=cosPhi*e0[idim]+sinPhi*e1[idim]+sinTheta*e2[idim];

      //calculate and output of the column integral
      PIC::ColumnIntegration::GetCoulumnIntegral(StateVector,StateVectorLength,x0,l,Exosphere::ColumnIntegral::CoulumnDensityIntegrant);

      //determine the nucleus projection
      double t,NucleusProjectionCode=-1.0;
      int iStartFace,iFinishFace;

      iStartFace=(PIC::Mesh::IrregularSurface::nBoundaryTriangleFaces/PIC::nTotalThreads)*PIC::ThisThread;
      iFinishFace=(PIC::Mesh::IrregularSurface::nBoundaryTriangleFaces/PIC::nTotalThreads)*(PIC::ThisThread+1);

      if (iFinishFace>PIC::Mesh::IrregularSurface::nBoundaryTriangleFaces) iFinishFace=PIC::Mesh::IrregularSurface::nBoundaryTriangleFaces;

      for (int i=iStartFace;i<iFinishFace;i++) {
        if (CutCell::BoundaryTriangleFaces[i].RayIntersection(x0,l,t,0.0)==true) {
          NucleusProjectionCode=1.0;
        }
      }

      //determine the minimum intersection time calculated by all processors
      double Buffer[PIC::nTotalThreads];
      MPI_Gather(&NucleusProjectionCode,1,MPI_DOUBLE,Buffer,1,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);

      //output the column integrals value
      if (PIC::ThisThread==0) {
        int cnt=0;

        for (int thread=0;thread<PIC::nTotalThreads;thread++) if (Buffer[thread]>0.0) NucleusProjectionCode=1.0;

        fprintf(fout,"%e %e %e ",ZenithAngle*180.0/Pi,AzimuthAngle*180.0/Pi,NucleusProjectionCode);

        //output the dust brightness
        for (int iCase=0;iCase<ElectricallyChargedDust::ScatteringEfficientcy::Mie::Fink2012Icarus::nCases;iCase++) fprintf(fout," %e",  StateVector[cnt++]);

        //output column densities of the background neutral species
        for (int s=0;s<Comet::CometData::nNeutrals;s++) fprintf(fout," %e",  StateVector[cnt++]);

        //column densities ofthe modeled species
        for (int s=0;s<PIC::nTotalSpecies;s++) fprintf(fout," %e",  StateVector[cnt++]);

        //finish the line
        fprintf(fout," \n");
      }

    }
  }

  if (PIC::ThisThread==0) fclose(fout);
}

















