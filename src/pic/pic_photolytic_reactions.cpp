//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//====================================================
//$Id$
//====================================================
//the model of photolytic reactions

#include "pic.h"

double *PIC::ChemicalReactions::PhotolyticReactions::ConstantTotalLifeTime=NULL;

void PIC::ChemicalReactions::PhotolyticReactions::Init() {
  //only one particle transformation model can be used
  if (_PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ == _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ON_) exit(__LINE__,__FILE__,"Error: only one particle transformation model can be used");

  // Initialization is intentionally idempotent.  This table is process-global
  // configuration state, and native applications may initialize a module once
  // before parsing and again after selecting a run-time chemistry deck.  Reuse
  // the allocation so a second call neither leaks the first table nor leaves
  // stale lifetimes active.  A negative value is the documented "reaction
  // disabled" sentinel; physical lifetimes accepted below are finite seconds
  // strictly greater than zero.
  if (ConstantTotalLifeTime==NULL) {
    ConstantTotalLifeTime=new double [PIC::nTotalSpecies];
  }

  for (int i=0;i<PIC::nTotalSpecies;i++) {
    ConstantTotalLifeTime[i]=-1.0;
  }
}


double PIC::ChemicalReactions::PhotolyticReactions::TotalLifeTime_default(double *x,int spec,long int ptr,bool &ReactionAllowedFlag,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  // Invalid species IDs and an uninitialized table are typed as unavailable,
  // not dereferenced.  This preserves the scheduler's no-transformation
  // contract and makes malformed run-time input fail closed without inventing
  // a chemical rate.  NaN and infinity are not lifetimes in SI seconds.
  if ((ConstantTotalLifeTime==NULL)||(spec<0)||(spec>=PIC::nTotalSpecies)) {
    ReactionAllowedFlag=false;
    return -1.0;
  }

  const double lifetime=ConstantTotalLifeTime[spec];
  ReactionAllowedFlag=(std::isfinite(lifetime)&&(lifetime>0.0));
  return lifetime;
}

int PIC::ChemicalReactions::PhotolyticReactions::PhotolyticReaction(double *x,long int ptr,int &spec,double &TimeInterval,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  int code=_PHOTOLYTIC_REACTIONS_NO_TRANSPHORMATION_;
  double p,lifetime;
  bool flag;

  // dt is an elapsed time in seconds.  Zero is an exact no-op; negative,
  // infinite, and NaN intervals are invalid domains.  The legacy API has no
  // error return distinct from "no transformation", so invalid domains fail
  // closed and leave TimeInterval unchanged.  Check the species before the
  // configured lifetime callback because the default callback indexes a
  // species-sized table.
  if ((spec<0)||(spec>=PIC::nTotalSpecies)||
      (!std::isfinite(TimeInterval))||(TimeInterval<=0.0)) {
    return _PHOTOLYTIC_REACTIONS_NO_TRANSPHORMATION_;
  }

  lifetime=_PIC_PHOTOLYTIC_REACTIONS__TOTAL_LIFETIME_(x,spec,ptr,flag,node);
  if ((flag==false)||(!std::isfinite(lifetime))||(lifetime<=0.0)) {
    return _PHOTOLYTIC_REACTIONS_NO_TRANSPHORMATION_;
  }

  // For a constant first-order rate k=1/tau, survival is exp(-dt/tau).
  // expm1 avoids cancellation in 1-exp(-dt/tau) when dt << tau.  Conditional
  // event time inversion uses t=-tau*log(1-U*p), where U is uniform on [0,1)
  // and p is the event probability over the requested interval.
  p=-expm1(-TimeInterval/lifetime);

  if (rnd()<p) {
    //determine the time offset when the reaction had occured
    TimeInterval=-lifetime*log1p(-rnd()*p);

    return _PHOTOLYTIC_REACTION_OCCURES_;
  }

  return code;
}

void PIC::ChemicalReactions::PhotolyticReactions::InitProductStatWeight() {
  bool NewWeightFound=false;
  int s,sSource,sProduct,nProducts,i,iReaction;
  bool UnknownWeightTableBefore[PIC::nTotalSpecies],UnknownWeightTableAfter[PIC::nTotalSpecies];
  int nUnknownWeightsBefore;

  double StatWeight[PIC::nTotalSpecies];

  for (s=0,nUnknownWeightsBefore=0;s<PIC::nTotalSpecies;s++) {
    StatWeight[s]=PIC::ParticleWeightTimeStep::GlobalParticleWeight[s];
    UnknownWeightTableBefore[s]=false;

    if (StatWeight[s]==0) UnknownWeightTableBefore[s]=true,nUnknownWeightsBefore++;
  }

  memcpy(UnknownWeightTableAfter,UnknownWeightTableBefore,PIC::nTotalSpecies*sizeof(bool));

  if (nUnknownWeightsBefore!=0) {
    do {
      NewWeightFound=false;

      for (sSource=0;sSource<PIC::nTotalSpecies;sSource++) if (UnknownWeightTableBefore[sSource]==false) {
        //check if weight of products of this reaction is unknown
        for (iReaction=0;iReaction<nTotalUnimolecularReactions;iReaction++) if (UnimoleculecularReactionDescriptor[iReaction].SourceSpecie==sSource) {
          nProducts=UnimoleculecularReactionDescriptor[iReaction].nProducts;

          for (i=0;i<nProducts;i++) {
            sProduct=UnimoleculecularReactionDescriptor[iReaction].ProductSpecies[i];

            if (UnknownWeightTableBefore[sProduct]==true) {
              //the weight of the ptoduct species is unknown yet
              StatWeight[sProduct]+=UnimoleculecularReactionDescriptor[iReaction].ReactionRate*StatWeight[sSource];
              UnknownWeightTableAfter[sProduct]=false;
              NewWeightFound=true;
            }
          }

        }
      }

      memcpy(UnknownWeightTableBefore,UnknownWeightTableAfter,PIC::nTotalSpecies*sizeof(bool));
    }
    while (NewWeightFound==true);
  }

  //update the stat weight
  for (s=0;s<PIC::nTotalSpecies;s++) PIC::ParticleWeightTimeStep::GlobalParticleWeight[s]=StatWeight[s];
}

//=======================================================================================================
//default photolytic reaction processor
//by default particle will stay in the system
void PIC::ChemicalReactions::PhotolyticReactions::PhotolyticReactionProcessor_default(long int ptr,long int& FirstParticleCell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  PIC::ParticleBuffer::SetNext(FirstParticleCell,ptr);
  PIC::ParticleBuffer::SetPrev(-1,ptr);

  if (FirstParticleCell!=-1) PIC::ParticleBuffer::SetPrev(ptr,FirstParticleCell);
  FirstParticleCell=ptr;
}

//=======================================================================================================
//execute the photolytic reaction model
void PIC::ChemicalReactions::PhotolyticReactions::ExecutePhotochemicalModel() {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  int i,j,k;
  long int oldFirstCellParticle,newFirstCellParticle,p,pnext;


#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#pragma omp parallel for schedule(dynamic,1) default (none)  \
    private (k,j,i,node,oldFirstCellParticle,newFirstCellParticle,p,pnext) \
    shared (DomainBlockDecomposition::nLocalBlocks,PIC::DomainBlockDecomposition::BlockTable)
#endif

for (int nLocalNode=0;nLocalNode<DomainBlockDecomposition::nLocalBlocks;nLocalNode++) for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++)  for (i=0;i<_BLOCK_CELLS_X_;i++) {
    double StartTime=MPI_Wtime();

    node=DomainBlockDecomposition::BlockTable[nLocalNode];
    if (node->block!=NULL) {
       oldFirstCellParticle=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
       newFirstCellParticle=-1;

       if (oldFirstCellParticle!=-1) {
         p=oldFirstCellParticle;

         while (p!=-1) {
           pnext=PIC::ParticleBuffer::GetNext(p);
           _PIC_PHOTOLYTIC_REACTIONS__REACTION_PROCESSOR_(p,newFirstCellParticle,node);
           p=pnext;
         }

         node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]=newFirstCellParticle;
      }

    if (_PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_) {
      node->ParallelLoadMeasure+=MPI_Wtime()-StartTime;
    }

    }
  }
}
