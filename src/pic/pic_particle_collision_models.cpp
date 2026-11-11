//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//$Id$
//modeling of collision between the model particles

#include "pic.h"

#include <limits>
#include <vector>

//sampling offset of the collision frequentcy
int PIC::MolecularCollisions::ParticleCollisionModel::CollsionFrequentcySampling::SamplingBufferOffset=-1;

//the limit of the collisions per a single particles
double PIC::MolecularCollisions::ParticleCollisionModel::CollisionLimitingThrehold=20.0;

PIC::MolecularCollisions::ParticleCollisionModel::cNTCValidationDiagnostic
    PIC::MolecularCollisions::ParticleCollisionModel::NTCValidationDiagnostic={};

//a user-defined function to accept a cell for modeling collisions
PIC::MolecularCollisions::ParticleCollisionModel::fDoSimulateCellCollisions PIC::MolecularCollisions::ParticleCollisionModel::DoSimulateCellCollisions=NULL;

namespace {

// Route every collision algorithm through one cross-section adapter.  The
// historical MolecularModels callback carries only species and velocities,
// which is sufficient for HS/VHS/VSS laws.  Application databases can also
// require cell state (for example, a temperature-dependent neutral law), so
// the user-defined route preserves the block and cell already available to
// NTC and MF instead of hiding them in mutable global state.
inline double EvaluateCollisionCrossSection(int species0,double *velocity0,
    int species1,double *velocity1,PIC::Mesh::cDataBlockAMR *block,
    PIC::Mesh::cDataCenterNode *cell) {
#if _PIC__PARTICLE_COLLISION_MODEL__CROSS_SECTION_ == _PIC__PARTICLE_COLLISION_MODEL__CROSS_SECTION__USER_DEFINED_
  return PIC::MolecularCollisions::ParticleCollisionModel::UserDefined::
      GetTotalCrossSection(velocity0,velocity1,species0,species1,block,cell);
#elif _PIC__PARTICLE_COLLISION_MODEL__CROSS_SECTION_ == _PIC__PARTICLE_COLLISION_MODEL__CROSS_SECTION__HS_
  return PIC::MolecularData::MolecularModels::GetTotalCrossSection(
      species0,velocity0,species1,velocity1);
#else
#error Unknown model-particle collision cross-section mode
#endif
}

} // namespace

void PIC::MolecularCollisions::ParticleCollisionModel::
ResetNTCValidationDiagnostic(bool enable) {
#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  // The focused Stage-4 validator uses one execution thread per MPI rank.
  // Reject an opt-in hybrid request until the counters and extrema below are
  // implemented as thread-local records with an explicit thread reduction;
  // silently racing a validation diagnostic would produce false evidence.
  if (enable)
    exit(__LINE__,__FILE__,
        "general NTC validation diagnostic requires an MPI-only build");
#endif
  NTCValidationDiagnostic=cNTCValidationDiagnostic{};
  NTCValidationDiagnostic.Enabled=enable;
  NTCValidationDiagnostic.MinimumUpdateProbability=
      std::numeric_limits<double>::infinity();
  NTCValidationDiagnostic.MaximumUpdateProbability=
      -std::numeric_limits<double>::infinity();
}

PIC::MolecularCollisions::ParticleCollisionModel::cNTCValidationDiagnostic
PIC::MolecularCollisions::ParticleCollisionModel::
GetNTCValidationDiagnostic() {
  return NTCValidationDiagnostic;
}



//request the sampling data
int PIC::MolecularCollisions::ParticleCollisionModel::RequestSamplingData(int offset) {
  int SamplingLength=0;

  if (_PIC__PARTICLE_COLLISION_MODEL__SAMPLE_COLLISION_FREQUENTCY_MODE__ == _PIC_MODE_ON_) {
    CollsionFrequentcySampling::SamplingBufferOffset = 
        offset + SamplingLength;
    SamplingLength += 
        sizeof(double) * (PIC::nTotalSpecies * PIC::nTotalSpecies);
  }

  return SamplingLength;
}

void PIC::MolecularCollisions::ParticleCollisionModel::PrintVariableList(FILE* fout,int DataSetNumber) {
  int s;

  for (s=0;s<PIC::nTotalSpecies;s++) fprintf(fout,", \"Coll Freq(%s) [m^{-3} s^{-1}]\"",PIC::MolecularData::GetChemSymbol(s));
}

void PIC::MolecularCollisions::ParticleCollisionModel::Interpolate(PIC::Mesh::cDataCenterNode** InterpolationList,double *InterpolationCoeficients,int nInterpolationCoeficients,PIC::Mesh::cDataCenterNode *CenterNode)  {
  double CollisionFrequentcy[PIC::nTotalSpecies*PIC::nTotalSpecies];
  int i,s;
  double *SamplingBuffer;

  for (i=0;i<PIC::nTotalSpecies*PIC::nTotalSpecies;i++) CollisionFrequentcy[i]=0.0;

  for (i=0;i<nInterpolationCoeficients;i++) {
    SamplingBuffer=(double*)(InterpolationList[i]->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset+CollsionFrequentcySampling::SamplingBufferOffset);

    for (s=0;s<PIC::nTotalSpecies*PIC::nTotalSpecies;s++) CollisionFrequentcy[s]+=SamplingBuffer[s]*InterpolationCoeficients[i];
  }

  SamplingBuffer=(double*)(CenterNode->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset+CollsionFrequentcySampling::SamplingBufferOffset);
  for (s=0;s<PIC::nTotalSpecies*PIC::nTotalSpecies;s++) SamplingBuffer[s]=CollisionFrequentcy[s]/((PIC::LastSampleLength!=0) ? PIC::LastSampleLength : 1);
}

void PIC::MolecularCollisions::ParticleCollisionModel::PrintData(FILE* fout,int DataSetNumber,CMPI_channel *pipe,int CenterNodeThread,PIC::Mesh::cDataCenterNode *CenterNode) {
  double t;
  double *SamplingBuffer=(double*)(CenterNode->GetAssociatedDataBufferPointer()+PIC::Mesh::completedCellSampleDataPointerOffset+CollsionFrequentcySampling::SamplingBufferOffset);
  int s;

  bool gather_output_data=false;

  if (pipe==NULL) gather_output_data=true;
  else if (pipe->ThisThread==CenterNodeThread) gather_output_data=true;


  for (s=0;s<PIC::nTotalSpecies;s++) {
    if (gather_output_data==true) { //(pipe->ThisThread==CenterNodeThread) {
      t= *(SamplingBuffer+CollsionFrequentcySampling::Offset(DataSetNumber,s));
    }

    if ((PIC::ThisThread==0)||(pipe==NULL)) { //if (pipe->ThisThread==0) {
      if ((CenterNodeThread!=0)&&(pipe!=NULL)) pipe->recv(t,CenterNodeThread);

      fprintf(fout,"%e ",t);
    }
    else pipe->send(t);
  }
}

void PIC::MolecularCollisions::ParticleCollisionModel::Init() {
  if (_PIC__PARTICLE_COLLISION_MODEL__SAMPLE_COLLISION_FREQUENTCY_MODE__ == _PIC_MODE_ON_) {
    //request sampling buffer and particle fields
    PIC::IndividualModelSampling::RequestSamplingData.push_back(RequestSamplingData);

    //print out of the otuput file
    PIC::Mesh::AddVaraibleListFunction(PrintVariableList);
    PIC::Mesh::PrintDataCenterNode.push_back(PrintData);
    PIC::Mesh::InterpolateCenterNode.push_back(Interpolate);
  }
}


//Non-time Counter collision procedure
//Modeling collision in a cell
/*!
 * Simplified version of the Non-time Counter (NTC) collision model implementation
 * This version assumes:
 * 1. All species have the same time step
 * 2. All species have the same particle weight
 * 3. No individual particle weight corrections (_INDIVIDUAL_PARTICLE_WEIGHT_MODE_ is off)
 */
void PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_ntc_simplified(
    int i, int j, int k, 
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  // This optimized kernel assumes equal statistical carriers and selects
  // model particles uniformly.  It therefore cannot be repaired merely by
  // applying a correction to the post-collision update probability.  Fail at
  // entry, before inspecting a node or particle, so every caller receives an
  // explicit setup error instead of a biased collision rate.  The validated
  // correction-aware implementation is ModelCellCollisions_ntc().
  SetExitErrorCode(__LINE__,1);
  exit(__LINE__,__FILE__,
      "simplified NTC does not support individual particle-weight corrections; use general NTC");
#endif
    
  const int SigmaCrMax_nTest = 100;
  const double SigmaCrMax_SafetyMargin = 1.3;

  // Structure to store particle data
  struct cParticleDataList {
    double *v;  // Direct pointer to velocity state vector
    PIC::ParticleBuffer::byte *ParticleData;
  };

  // Get block and check validity
  PIC::Mesh::cDataBlockAMR *block = node->block;
  if (block == NULL) return;

  // Get cell and check validity
  PIC::Mesh::cDataCenterNode *cell = block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k));
  if (cell == NULL) return;

  // Get first particle in cell and cell measure
  long int FirstCellParticle = block->FirstCellParticleTable[i + _BLOCK_CELLS_X_*(j + _BLOCK_CELLS_Y_*k)];
  double cellMeasure = cell->Measure;

  // Get sampling data buffer if sampling mode is on
  char *SamplingData = NULL;
  if (_PIC__PARTICLE_COLLISION_MODEL__SAMPLE_COLLISION_FREQUENTCY_MODE__ == _PIC_MODE_ON_) {
    SamplingData = cell->GetAssociatedDataBufferPointer() + PIC::Mesh::collectingCellSampleDataPointerOffset;
  }

  if (FirstCellParticle == -1) return;

  // Count particles per species
  int nParticleNumber[PIC::nTotalSpecies] = {0};
  for (long int ptr = FirstCellParticle; ptr != -1; ptr = PIC::ParticleBuffer::GetNext(ptr)) {
    nParticleNumber[PIC::ParticleBuffer::GetI(ptr)]++;
  }

  // Find maximum number of particles among all species
  int ParticleDataListLength = 0;
  for (int s = 0; s < PIC::nTotalSpecies; s++) {
    if (ParticleDataListLength < nParticleNumber[s]) {
      ParticleDataListLength = nParticleNumber[s];
    }
  }

  // Allocate particle lists
  cParticleDataList *s0ParticleDataList = new cParticleDataList[ParticleDataListLength];
  cParticleDataList *s1ParticleDataList = new cParticleDataList[ParticleDataListLength];

  // Loop through first species s0
  for (int s0 = 0; s0 < PIC::nTotalSpecies; s0++) {
    if (nParticleNumber[s0] == 0 || PIC::MolecularData::GetSpecieType(s0) != _PIC_SPECIE_TYPE__GAS_) {
      continue;
    }

    double m0 = PIC::MolecularData::GetMass(s0);
    double LocalTimeStep = block->GetLocalTimeStep(s0);  // Same for all species
    double LocalParticleWeight = block->GetLocalParticleWeight(s0);  // Same for all species

    // Populate s0 particle list
    cParticleDataList *s0List = s0ParticleDataList;
    int cnt = 0;
    for (long int ptr = FirstCellParticle; ptr != -1; ptr = PIC::ParticleBuffer::GetNext(ptr)) {
      if (PIC::ParticleBuffer::GetI(ptr) == (unsigned)s0) {
        PIC::ParticleBuffer::byte *ParticleData = PIC::ParticleBuffer::GetParticleDataPointer(ptr);
        s0List[cnt].ParticleData = ParticleData;
        s0List[cnt].v = PIC::ParticleBuffer::GetV(ParticleData);
        cnt++;
      }
    }

    // Loop through second species s1
    for (int s1 = s0; s1 < PIC::nTotalSpecies; s1++) {
      if (nParticleNumber[s1] == 0 || PIC::MolecularData::GetSpecieType(s1) != _PIC_SPECIE_TYPE__GAS_) {
        continue;
      }

      // Check if collisions should be simulated for this pair
      if (DoSimulateCellCollisions != NULL) {
        if (!DoSimulateCellCollisions(cell, s0, s1)) continue;
      }

      double m1 = PIC::MolecularData::GetMass(s1);

      // Populate s1 particle list (if different from s0)
      cParticleDataList *s1List = (s0 == s1) ? s0List : s1ParticleDataList;
      if (s0 != s1) {
        cnt = 0;
        for (long int ptr = FirstCellParticle; ptr != -1; ptr = PIC::ParticleBuffer::GetNext(ptr)) {
          if (PIC::ParticleBuffer::GetI(ptr) == (unsigned)s1) {
            PIC::ParticleBuffer::byte *ParticleData = PIC::ParticleBuffer::GetParticleDataPointer(ptr);
            s1List[cnt].ParticleData = ParticleData;
            s1List[cnt].v = PIC::ParticleBuffer::GetV(ParticleData);
            cnt++;
          }
        }
      }

      // Calculate maximum collision cross section
      double SigmaCrMax = 0.0;
      for (int ntest = 0; ntest < SigmaCrMax_nTest; ntest++) {
        int s0ptr = (int)(rnd()*nParticleNumber[s0]);
        int s1ptr = (int)(rnd()*nParticleNumber[s1]);
        
        double cr = Vector3D::Distance(s1List[s1ptr].v,s0List[s0ptr].v); 
        double SigmaCr = cr * EvaluateCollisionCrossSection(
            s0,s0List[s0ptr].v,s1,s1List[s1ptr].v,block,cell);
        
        if (SigmaCr > SigmaCrMax) SigmaCrMax = SigmaCr;
      }
      
      SigmaCrMax *= SigmaCrMax_SafetyMargin;

      // Calculate number of collision pairs
      double ancoll;
      if (s0 == s1) {
        ancoll = 0.5 * nParticleNumber[s0] * (nParticleNumber[s0] /*- 1*/) * 
                 LocalParticleWeight * SigmaCrMax * LocalTimeStep / cellMeasure;
      } else {
        ancoll = nParticleNumber[s0] * nParticleNumber[s1] * 
                 LocalParticleWeight * SigmaCrMax * LocalTimeStep / cellMeasure;
      }

      // Apply unbiased collision-trial thinning. CollisionLimitingThrehold is
      // dimensionless and means retained trials per model particle per step;
      // the cap is therefore L*max(N0,N1), not N0*N1. If A prospective
      // trials exceed cap C, run C trials and multiply every accepted sample
      // by A/C. This preserves the expected physical rate while preventing a
      // long sequence of collisions from evolving velocities beyond the
      // single pre-loop sigma*g majorant. A non-positive/non-finite limit is
      // invalid configuration rather than a silent no-collision model.
      double CollisionLimitingFactor = 1.0;
      if (!(CollisionLimitingThrehold > 0.0) ||
          !isfinite(CollisionLimitingThrehold))
        exit(__LINE__,__FILE__,
            "CollisionLimitingThrehold must be positive and finite");
      const double maxCollisionTrials = CollisionLimitingThrehold *
          max(nParticleNumber[s0],nParticleNumber[s1]);
      if (maxCollisionTrials < ancoll) {
        CollisionLimitingFactor = ancoll / maxCollisionTrials;
        ancoll /= CollisionLimitingFactor;
      }

      // Determine actual number of collisions
      long int ncoll = (long int)ancoll;
      if (rnd() < (ancoll - ncoll)) ncoll++;

      // Simulate collisions
      while (ncoll-- > 0) {
        // Select random particle pairs
        int s0ptr = (int)(rnd() * nParticleNumber[s0]);
        int s1ptr;
        do {
          s1ptr = (int)(rnd() * nParticleNumber[s1]);
        } while (s0 == s1 && s0ptr == s1ptr);

        double *v0 = s0List[s0ptr].v;
        double *v1 = s1List[s1ptr].v;

        // Calculate relative velocity
        double cr = 0.0;
        double vrel[3], vcm[3];

        for (int idim = 0; idim < 3; idim++) {
          vrel[idim] = v1[idim] - v0[idim];
          vcm[idim] = (m1*v1[idim] + m0*v0[idim])/(m1 + m0);
          cr += vrel[idim] * vrel[idim];
        }
        cr = sqrt(cr);

        // Check if collision occurs
        double SigmaCr = cr * EvaluateCollisionCrossSection(
            s0,v0,s1,v1,block,cell);
        if (rnd() * SigmaCrMax >= SigmaCr) continue;

        //model the internal energy exchange
        if (_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_) {
          double crInit = cr;
	  bool UpdateFlag[2]={true,true};

          PIC::IDF::RedistributeEnergy(s0List[s0ptr].ParticleData,s1List[s1ptr].ParticleData, cr, UpdateFlag, cell);

          for (int idim = 0; idim < 3; idim++) {
            vrel[idim] = (cr > 1.0E-10) ? vrel[idim] * cr / crInit : 0.0;
          }
        }


        // Perform collision - directly modifies v0 and v1 which point to particle state
        // The generated molecular model owns both sigma(g) and its angular
        // law. HS/VHS remain isotropic; VSS consumes the configured alpha.
        if (!PIC::MolecularCollisions::VelocityScattering::Scatter(
                s0,v0,m0,s1,v1,m1,vcm,vrel))
          exit(__LINE__,__FILE__,
              "Invalid or unsupported molecular angular-scattering model");

        // Sample collision frequency if sampling is enabled
        if (_PIC__PARTICLE_COLLISION_MODEL__SAMPLE_COLLISION_FREQUENTCY_MODE__ == _PIC_MODE_ON_) {
          // Sample for first species
          int CollFreqOffset_s0 = CollsionFrequentcySampling::SamplingBufferOffset + 
                                 sizeof(double) * CollsionFrequentcySampling::Offset(s0, s1);
          *((double*)(SamplingData + CollFreqOffset_s0)) += LocalParticleWeight / 
                                                           LocalTimeStep / cellMeasure * 
                                                           CollisionLimitingFactor;

          // Sample for second species (if different from first)
            int CollFreqOffset_s1 = CollsionFrequentcySampling::SamplingBufferOffset + 
                                   sizeof(double) * CollsionFrequentcySampling::Offset(s1, s0);
            *((double*)(SamplingData + CollFreqOffset_s1)) += LocalParticleWeight / 
                                                             LocalTimeStep / cellMeasure * 
                                                             CollisionLimitingFactor;
        }
      }
    }
  }

  // Clean up
  delete [] s1ParticleDataList;
  delete [] s0ParticleDataList;
}

void PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_ntc(int i,int j,int k, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  int SigmaCrMax_nTest=100;
  double SigmaCrMax_SafetyMargin=1.05;

  struct cParticleDataList {
    double vel[3];
    // The base block weight is common to one species bucket.  Consequently
    // selecting a particle in proportion to its physical carrier W*c is
    // equivalent to selecting it in proportion to c.  Keep both c and its
    // prefix sum with the compact velocity state so every NTC trial can draw
    // the statistically correct pair without another particle-buffer read.
    double WeightCorrection;
    double CumulativeWeightCorrection;
    bool ValueChangedFlag;
    PIC::ParticleBuffer::byte *ParticleData;
  };

  // Retain one independently reusable set of species buckets per execution
  // thread.  The old implementation allocated two arrays for every occupied
  // cell and repeatedly traversed the cell linked list for each species pair.
  // CollisionModelManager may dispatch cells through OpenMP, so thread_local
  // storage is required even though the current comet validation uses one
  // OpenMP thread per MPI rank.  clear() preserves capacity: after a thread
  // has encountered its largest cell, later cells do not call the allocator.
  static thread_local std::vector<cParticleDataList>
      ParticleDataBySpecies[PIC::nTotalSpecies];

  int s,s0,s1;

  int nParticleNumber[PIC::nTotalSpecies],cnt;
  long int FirstCellParticle,ptr;
  double cellMeasure;
  int thread;
  PIC::Mesh::cDataBlockAMR *block;

  // s0List and s1List point into the thread-local species buckets above.
  cParticleDataList *s0List=NULL,*s1List=NULL;

  //sampling data
  PIC::Mesh::cDataCenterNode *cell;
  char *SamplingData=NULL;

  block=node->block;
  if (block==NULL) return;

  FirstCellParticle=block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
  cell=block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k));
  if (cell==NULL) return;

  //simulate particle's collisions
  cellMeasure=cell->Measure;

  if (_PIC__PARTICLE_COLLISION_MODEL__SAMPLE_COLLISION_FREQUENTCY_MODE__ == _PIC_MODE_ON_) {
    SamplingData=cell->GetAssociatedDataBufferPointer()+PIC::Mesh::collectingCellSampleDataPointerOffset;
  }

  if (FirstCellParticle!=-1) {
    //simulate collision in the cell
    // 1. Bucket every particle by species in one linked-list traversal.  The
    // velocity is copied because trial collisions operate on a compact array;
    // ParticleData retains the authoritative AMPS record for final updates.
    // Weight-correction extrema and sums are also species invariants within
    // this cell collision call, so calculate them once rather than once for
    // every (s0,s1) pair.
    double minWeightCorrection[PIC::nTotalSpecies];
    double sumWeightCorrection[PIC::nTotalSpecies];
    for (s=0;s<PIC::nTotalSpecies;s++) {
      ParticleDataBySpecies[s].clear();
      nParticleNumber[s]=0;
      minWeightCorrection[s]=1.0;
      sumWeightCorrection[s]=0.0;
    }

    for (ptr=FirstCellParticle;ptr!=-1;
         ptr=PIC::ParticleBuffer::GetNext(ptr)) {
      const int species=PIC::ParticleBuffer::GetI(ptr);
      // This kernel has always collided gas species only. Dust corrections
      // may span a steep size distribution and are not a gas probability
      // measure: even adding a tiny grain weight to a large one can lose all
      // precision. Skip non-gas records before building sums/majorants, just
      // as the species-pair loop below does, preserving the separate grain
      // dynamics and avoiding a new rejection of otherwise valid dust data.
      if (PIC::MolecularData::GetSpecieType(species)!=_PIC_SPECIE_TYPE__GAS_)
        continue;
      PIC::ParticleBuffer::byte *ParticleData=
          PIC::ParticleBuffer::GetParticleDataPointer(ptr);
      std::vector<cParticleDataList> &list=ParticleDataBySpecies[species];
      cParticleDataList entry;
      PIC::ParticleBuffer::GetV(entry.vel,ParticleData);
      entry.ParticleData=ParticleData;
      entry.ValueChangedFlag=false;

      // GetIndividualStatWeightCorrection() returns one in correction-off
      // builds.  Using the same path in both configurations keeps the NTC
      // estimator identical algebraically and limits the special handling to
      // the particle-record accessor.  A non-positive/non-finite correction
      // cannot define a probability measure and must fail before it reaches a
      // collision majorant or an update-probability denominator.
      entry.WeightCorrection=
          PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);
      if (!(entry.WeightCorrection>0.0) ||
          !std::isfinite(entry.WeightCorrection))
        exit(__LINE__,__FILE__,
            "general NTC requires a positive finite particle-weight correction");

      const double nextCorrectionSum=
          sumWeightCorrection[species]+entry.WeightCorrection;
      if (!(nextCorrectionSum>sumWeightCorrection[species]) ||
          !std::isfinite(nextCorrectionSum))
        exit(__LINE__,__FILE__,
            "general NTC particle-weight-correction sum is not representable");
      sumWeightCorrection[species]=nextCorrectionSum;
      entry.CumulativeWeightCorrection=nextCorrectionSum;
      list.push_back(entry);

      if (nParticleNumber[species]==0 ||
          entry.WeightCorrection<minWeightCorrection[species])
        minWeightCorrection[species]=entry.WeightCorrection;
      nParticleNumber[species]++;
    }

    // Draw from a species bucket with probability c_i/sum(c).  The physical
    // carrier of entry i is q_i=W_species*c_i, and W_species is common to the
    // bucket, so this is exactly a q-weighted draw.  The binary search operates
    // on strictly increasing prefix sums checked during bucket construction.
    auto SelectCorrectionWeightedParticle=[](cParticleDataList *list,
        int count,double correctionSum) -> int {
      if (list==NULL || count<=0 || !(correctionSum>0.0) ||
          !std::isfinite(correctionSum))
        exit(__LINE__,__FILE__,
            "invalid species bucket in general NTC weighted selection");

      double target=rnd()*correctionSum;
      // rnd() is expected to be in [0,1), but map a rounded upper endpoint to
      // the final positive bin rather than returning an out-of-range index.
      if (!(target<correctionSum)) return count-1;

      int left=0,right=count-1;
      while (left<right) {
        const int middle=left+(right-left)/2;
        if (target<list[middle].CumulativeWeightCorrection) right=middle;
        else left=middle+1;
      }
      return left;
    };

    //loop through the first species
    double m0,m1,am;
    double LocalParticleWeight_s0,LocalParticleWeight_s1;
    double LocalTimeStep_s0,LocalTimeStep_s1;
    double minParticleWeightCorrection_s0,minParticleWeightCorrection_s1,sumWeightCorrection_s0,sumWeightCorrection_s1;

    for (s0=0;s0<PIC::nTotalSpecies;s0++) if ((PIC::MolecularData::GetSpecieType(s0)==_PIC_SPECIE_TYPE__GAS_)&&(nParticleNumber[s0]!=0)) {
      m0=PIC::MolecularData::GetMass(s0);
      LocalParticleWeight_s0=block->GetLocalParticleWeight(s0);
      LocalTimeStep_s0=block->GetLocalTimeStep(s0);
      minParticleWeightCorrection_s0=minWeightCorrection[s0];
      sumWeightCorrection_s0=sumWeightCorrection[s0];

      // The species may already have appeared as s1 for a lower-index pair.
      // Its compact velocities contain those committed updates; only the
      // per-pair writeback flags need resetting before it becomes s0.
      s0List=ParticleDataBySpecies[s0].data();
      for (cnt=0;cnt<nParticleNumber[s0];cnt++)
        s0List[cnt].ValueChangedFlag=false;

      //loop through the second speices
      for (s1=s0;s1<PIC::nTotalSpecies;s1++) if ((PIC::MolecularData::GetSpecieType(s1)==_PIC_SPECIE_TYPE__GAS_)&&(nParticleNumber[s1]!=0)) {
        // One model particle cannot collide with itself.  The former repeated
        // uniform draw never terminated for this case when a non-zero trial
        // count was produced.
        if (s0==s1 && nParticleNumber[s0]<2) continue;

        //call a user-defined funstion to verify whether collision should be modeled in the cell
        if (DoSimulateCellCollisions!=NULL) {
          if (DoSimulateCellCollisions(cell,s0,s1)==false) continue;
        }

        m1=PIC::MolecularData::GetMass(s1);
        am=m0+m1;
        LocalParticleWeight_s1=block->GetLocalParticleWeight(s1);
        LocalTimeStep_s1=block->GetLocalTimeStep(s1);
        sumWeightCorrection_s1=sumWeightCorrection[s1];

        // Select the already-populated species bucket.  For a distinct s1,
        // reset writeback flags just as the historical repopulation loop did;
        // velocities remain current because every accepted collision updates
        // this same bucket before it is committed to ParticleData.
        if (s0==s1) s1List=s0List,minParticleWeightCorrection_s1=minParticleWeightCorrection_s0,sumWeightCorrection_s1=sumWeightCorrection_s0;
        else {
          s1List=ParticleDataBySpecies[s1].data();
          minParticleWeightCorrection_s1=minWeightCorrection[s1];
          for (cnt=0;cnt<nParticleNumber[s1];cnt++)
            s1List[cnt].ValueChangedFlag=false;
        }

        //simulate collsions between the pair of species
        double SigmaCrMax=0.0,SigmaCr,ancoll;
        long int ncoll;
        double v0[3],v1[3],cr;

        //1.Evaluate the maximum value of SigmaCr
        for (int ntest=0;ntest<SigmaCrMax_nTest;ntest++) {
          // This loop constructs an envelope; it is not an estimator of the
          // collision integral.  Probe model particles uniformly so even a
          // low-correction, high-relative-speed tail can establish SigmaCrMax.
          // Conditioning this envelope-only probe on distinct same-species
          // particles does not alter the independently weighted NTC trials
          // below.  The one-particle case was rejected before entering here.
          const int s0Test=static_cast<int>(rnd()*nParticleNumber[s0]);
          int s1Test;
          do s1Test=static_cast<int>(rnd()*nParticleNumber[s1]);
          while (s0==s1 && s0Test==s1Test);
          memcpy(v0,s0List[s0Test].vel,3*sizeof(double));
          memcpy(v1,s1List[s1Test].vel,3*sizeof(double));
          cr=sqrt(pow(v1[0]-v0[0],2)+pow(v1[1]-v0[1],2)+pow(v1[2]-v0[2],2));

          SigmaCr=cr*EvaluateCollisionCrossSection(s0,v0,s1,v1,block,cell);

          if (SigmaCr>SigmaCrMax) SigmaCrMax=SigmaCr;
        }

        SigmaCrMax*=SigmaCrMax_SafetyMargin;

        //2.Evaluate the prospective number of collisions
        double maxLocalTimeStep,minLocalParticleWeight;

        maxLocalTimeStep=max(LocalTimeStep_s0,LocalTimeStep_s1);
        minLocalParticleWeight=min(LocalParticleWeight_s0*minParticleWeightCorrection_s0,LocalParticleWeight_s1*minParticleWeightCorrection_s1);

        const double totalPhysicalWeight_s0=
            sumWeightCorrection_s0*LocalParticleWeight_s0;
        const double totalPhysicalWeight_s1=
            sumWeightCorrection_s1*LocalParticleWeight_s1;

        // Let q_i=W_s*c_i be the physical weight of particle i and q_min the
        // smallest carrier in this species pair.  Cross-species NTC uses
        //
        //   N_trial = sum_i(q_i) sum_j(q_j) sigma_g,max dt_max
        //             / (q_min V).
        //
        // Selecting i and j with probabilities q_i/sum(q) makes one ordered
        // pair's trial rate q_i*q_j/q_min.  Updating each selected particle
        // with probability q_min/q_i (and its local-dt factor) then recovers
        // the correct q_i*q_j collision-integral contribution.  Uniform pair
        // selection is biased whenever c is correlated with velocity.
        //
        // For one species, the factor 1/2 converts ordered to unordered pairs.
        // Independent weighted draws deliberately include i==j in N_trial;
        // those zero-measure self trials are discarded below.  This preserves
        // the exact sum over i!=j without conditional-redraw normalization.
        if (s0==s1)
          ancoll=0.5*totalPhysicalWeight_s0*totalPhysicalWeight_s0*
              SigmaCrMax*LocalTimeStep_s0/
              (minLocalParticleWeight*cellMeasure);
        else
          ancoll=totalPhysicalWeight_s0*totalPhysicalWeight_s1*
              SigmaCrMax*maxLocalTimeStep/
              (minLocalParticleWeight*cellMeasure);

        //3. Collision limiting: use the same per-particle cap and unbiased
        //sample multiplier as the simplified NTC path above. The prior
        //expression ancoll/N0*N1 also had left-to-right precedence, producing
        //ancoll*N1/N0 rather than ancoll/(N0*N1); the explicit cap removes
        //that ambiguity and gives this configuration variable its stated
        //meaning across selector implementations.
        double CollisionLimitingFactor=1.0;

        if (!(CollisionLimitingThrehold>0.0) ||
            !isfinite(CollisionLimitingThrehold))
          exit(__LINE__,__FILE__,
              "CollisionLimitingThrehold must be positive and finite");
        const double maxCollisionTrials=CollisionLimitingThrehold*
            max(nParticleNumber[s0],nParticleNumber[s1]);
        if (maxCollisionTrials<ancoll) {
          CollisionLimitingFactor=ancoll/maxCollisionTrials;
          ancoll/=CollisionLimitingFactor;
        }

        if (!(ancoll>=0.0) || !std::isfinite(ancoll) ||
            ancoll>static_cast<double>(std::numeric_limits<long int>::max()))
          exit(__LINE__,__FILE__,
              "invalid prospective collision count in general NTC");

        ncoll=(long int)ancoll;
        ancoll-=ncoll;
        if (rnd()<ancoll) ncoll++;

        //4. simulate collisions
        int s0ptr,s1ptr,idim;
        double vrel[3],vcm[3];

        while (ncoll-->0) {
          if (NTCValidationDiagnostic.Enabled)
            NTCValidationDiagnostic.TrialCount++;

          s0ptr=SelectCorrectionWeightedParticle(s0List,
              nParticleNumber[s0],sumWeightCorrection_s0);
          s1ptr=SelectCorrectionWeightedParticle(s1List,
              nParticleNumber[s1],sumWeightCorrection_s1);

          // Do not redraw a same-species self pair.  Redrawing conditions the
          // second-particle distribution and destroys the independent
          // q_i*q_j normalization used in ancoll.  Discarding it instead
          // removes exactly the nonphysical i==j term and also eliminates the
          // historical one-particle infinite loop.
          if (s0==s1 && s0ptr==s1ptr) {
            if (NTCValidationDiagnostic.Enabled)
              NTCValidationDiagnostic.SelfPairRejectionCount++;
            continue;
          }

          memcpy(v0,s0List[s0ptr].vel,3*sizeof(double));
          memcpy(v1,s1List[s1ptr].vel,3*sizeof(double));

          for (cr=0.0,idim=0;idim<3;idim++) {
            vrel[idim]=v1[idim]-v0[idim];
            cr+=vrel[idim]*vrel[idim];
	    vcm[idim] = (m1*v1[idim] + m0*v0[idim])/(m1 + m0);
          }

          cr=sqrt(cr);
          SigmaCr=cr*EvaluateCollisionCrossSection(s0,v0,s1,v1,block,cell);

          if (rnd()*SigmaCrMax>=SigmaCr) continue;
          if (NTCValidationDiagnostic.Enabled)
            NTCValidationDiagnostic.AcceptedCollisionCount++;

          //determine weather the properties of particles were updated
          double pUpdate_s0,pUpdate_s1;
          bool UpdateFlag[2];

          pUpdate_s0=minLocalParticleWeight/
              (LocalParticleWeight_s0*s0List[s0ptr].WeightCorrection)*
              LocalTimeStep_s0/maxLocalTimeStep;
          pUpdate_s1=minLocalParticleWeight/
              (LocalParticleWeight_s1*s1List[s1ptr].WeightCorrection)*
              LocalTimeStep_s1/maxLocalTimeStep;

          // q_min<=q_selected and dt_species<=dt_max imply p in [0,1]
          // analytically.  Reject non-finite or materially out-of-range values
          // so an invalid correction cannot silently become an always/never
          // update.  Clamp only the few ulps permitted for roundoff in the
          // products and ratios above.
          const double probabilityTolerance=
              64.0*std::numeric_limits<double>::epsilon();
          auto ValidateUpdateProbability=[&](double probability) -> double {
            const bool invalid=!std::isfinite(probability) ||
                probability < -probabilityTolerance ||
                probability > 1.0+probabilityTolerance;
            if (NTCValidationDiagnostic.Enabled) {
              NTCValidationDiagnostic.UpdateProbabilityCount++;
              if (std::isfinite(probability)) {
                NTCValidationDiagnostic.MinimumUpdateProbability=std::min(
                    NTCValidationDiagnostic.MinimumUpdateProbability,
                    probability);
                NTCValidationDiagnostic.MaximumUpdateProbability=std::max(
                    NTCValidationDiagnostic.MaximumUpdateProbability,
                    probability);
              }
              if (invalid)
                NTCValidationDiagnostic.InvalidUpdateProbabilityCount++;
            }
            if (invalid)
              exit(__LINE__,__FILE__,
                  "general NTC collision update probability is outside [0,1]");
            if (probability<0.0) {
              if (NTCValidationDiagnostic.Enabled)
                NTCValidationDiagnostic.RoundoffClampCount++;
              return 0.0;
            }
            if (probability>1.0) {
              if (NTCValidationDiagnostic.Enabled)
                NTCValidationDiagnostic.RoundoffClampCount++;
              return 1.0;
            }
            return probability;
          };
          pUpdate_s0=ValidateUpdateProbability(pUpdate_s0);
          pUpdate_s1=ValidateUpdateProbability(pUpdate_s1);

          UpdateFlag[0]=(rnd()<pUpdate_s0) ? true :false;
          UpdateFlag[1]=(rnd()<pUpdate_s1) ? true :false;
          if (NTCValidationDiagnostic.Enabled) {
            if (UpdateFlag[0]) NTCValidationDiagnostic.ParticleUpdateCount++;
            if (UpdateFlag[1]) NTCValidationDiagnostic.ParticleUpdateCount++;
          }

          //model the internal energy exchange
          if (_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_) {
            double crInit=cr;

            PIC::IDF::RedistributeEnergy(s0List[s0ptr].ParticleData,s1List[s1ptr].ParticleData,cr,UpdateFlag,cell);
            for (int idim=0;idim<3;idim++) {
              vrel[idim] = (cr > 1.0E-10) ? vrel[idim]*cr/crInit : 0.0;
            }
          }

          //the collision is considered to be true
          if (!PIC::MolecularCollisions::VelocityScattering::Scatter(
                  s0,v0,m0,s1,v1,m1,vcm,vrel))
            exit(__LINE__,__FILE__,
                "Invalid or unsupported molecular angular-scattering model");

          //update the velocities in the lists
          if (UpdateFlag[0]==true) {
            s0List[s0ptr].ValueChangedFlag=true;
            memcpy(s0List[s0ptr].vel,v0,3*sizeof(double));

            if (_PIC__PARTICLE_COLLISION_MODEL__SAMPLE_COLLISION_FREQUENTCY_MODE__ == _PIC_MODE_ON_) {
              int CollFreqOffset=CollsionFrequentcySampling::SamplingBufferOffset+sizeof(double)*CollsionFrequentcySampling::Offset(s0,s1);

              *((double*)(SamplingData+CollFreqOffset))+=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(s0List[s0ptr].ParticleData)*LocalParticleWeight_s0/LocalTimeStep_s0/cellMeasure*CollisionLimitingFactor; // calculate collision frequency taking into account the collision limiting factor
            }
          }

          if (UpdateFlag[1]==true) {
            s1List[s1ptr].ValueChangedFlag=true;
            memcpy(s1List[s1ptr].vel,v1,3*sizeof(double));

            if (_PIC__PARTICLE_COLLISION_MODEL__SAMPLE_COLLISION_FREQUENTCY_MODE__ == _PIC_MODE_ON_) {
              int CollFreqOffset=CollsionFrequentcySampling::SamplingBufferOffset+sizeof(double)*CollsionFrequentcySampling::Offset(s1,s0);

              *((double*)(SamplingData+CollFreqOffset))+=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(s1List[s1ptr].ParticleData)*LocalParticleWeight_s1/LocalTimeStep_s1/cellMeasure*CollisionLimitingFactor; // calculate collision frequency taking into account the collision limiting factor
            }
          }
        }

        //update the velocities of the species 's1'
        if (s0!=s1) {
          for (cnt=0;cnt<nParticleNumber[s1];cnt++) if (s1List[cnt].ValueChangedFlag==true) {
            PIC::ParticleBuffer::SetV(s1List[cnt].vel,s1List[cnt].ParticleData);
          }
        }
      }

      //update velocities of species 's0'
      for (cnt=0;cnt<nParticleNumber[s0];cnt++) if (s0List[cnt].ValueChangedFlag==true) {
        PIC::ParticleBuffer::SetV(s0List[cnt].vel,s0List[cnt].ParticleData);
      }
    }
  }
}

//NTC: Manager
void PIC::MolecularCollisions::ParticleCollisionModel::ntc() {
  CollisionModelManager(ModelCellCollisions_ntc);
}

void PIC::MolecularCollisions::ParticleCollisionModel::CollisionModelManager(void (*fCellCollision)(int, int, int, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)) {
  int thread;

  //global offsets
  int LoadBalancingMeasureOffset=PIC::Mesh::cDataBlockAMR_static_data::LoadBalancingMeasureOffset;

  //simulate particle's collisions
  #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  if (_PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_) {
    //reset the balancing counters
    for (int nLocalNode=0;nLocalNode<DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
      for (int thread=0;thread<PIC::nTotalThreadsOpenMP;thread++) {
        cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=DomainBlockDecomposition::BlockTable[nLocalNode];

        if (node->block!=NULL) {
          *(thread+(double*)(node->block->GetAssociatedDataBufferPointer()+
              LoadBalancingMeasureOffset))=0.0;
        }
      }
    }
  } //_PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_

  #pragma omp parallel for schedule(dynamic,1) 
  #endif  //_COMPILATION_MODE_
  for (int CellCounter=0;CellCounter<DomainBlockDecomposition::nLocalBlocks*_BLOCK_CELLS_Z_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_;CellCounter++) {
    int nLocalNode,ii=CellCounter,i,j,k;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;

    nLocalNode=ii/(_BLOCK_CELLS_Z_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_);
    ii-=nLocalNode*_BLOCK_CELLS_Z_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_;

    k=ii/(_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_);
    ii-=k*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_;

    j=ii/_BLOCK_CELLS_X_;
    ii-=j*_BLOCK_CELLS_X_;

    i=ii;

    double StartTime=MPI_Wtime();

    node=DomainBlockDecomposition::BlockTable[nLocalNode];
    auto block=node->block;

    if (block==NULL) continue;  

    fCellCollision(i,j,k,node); 

    #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
    thread=omp_get_thread_num(); 
    #else
    thread=0;
    #endif

    
    if (_PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_) {
      switch (_COMPILATION_MODE_) {
      case _COMPILATION_MODE__MPI_:
        node->ParallelLoadMeasure+=MPI_Wtime()-StartTime;
        break;
      case _COMPILATION_MODE__HYBRID_:
      {
        int iThread_OpenMP=0;

        #if _COMPILATION_MODE_==_COMPILATION_MODE__HYBRID_
        iThread_OpenMP=omp_get_thread_num();
        #endif

        *(iThread_OpenMP+(double*)(block->GetAssociatedDataBufferPointer()+LoadBalancingMeasureOffset))+=MPI_Wtime()-StartTime;
      }
      break;
      default:
        exit(__LINE__,__FILE__,"The option is unknown");
      }
    }
  }

  if ((_COMPILATION_MODE_==_COMPILATION_MODE__HYBRID_)&&(_PIC_DYNAMIC_LOAD_BALANCING_MODE_==_PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_)) {
    //sum-up the balancing measure
    for (int nLocalNode=0;nLocalNode<DomainBlockDecomposition::nLocalBlocks;nLocalNode++) for (int thread=0;thread<PIC::nTotalThreadsOpenMP;thread++) {
      auto node=DomainBlockDecomposition::BlockTable[nLocalNode];
      node->ParallelLoadMeasure+=*(thread+(double*)(node->block->GetAssociatedDataBufferPointer()+LoadBalancingMeasureOffset));
    }
  }
}


//Majorant Frequency collision procedure
//Modeling collision in a cell
void PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_mf_simplified(int i, int j, int k,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) { 
  const int SigmaCrMax_nTest = 100;
  const double SigmaCrMax_SafetyMargin = 1.05;

  // Structure to store particle data
  struct cParticleDataList {
    double *v;  // Direct pointer to velocity state vector
    PIC::ParticleBuffer::byte *ParticleData;
  };

  // Get block and check validity
  PIC::Mesh::cDataBlockAMR *block = node->block;
  if (block == NULL) return;

  // Get cell and check validity
  PIC::Mesh::cDataCenterNode *cell = block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k));
  if (cell == NULL) return;

  double cellMeasure=cell->Measure;

  // Get first particle in cell and cell measure
  long int FirstCellParticle = block->FirstCellParticleTable[i + _BLOCK_CELLS_X_*(j + _BLOCK_CELLS_Y_*k)];
  // Get sampling data buffer if sampling mode is on
  char *SamplingData = NULL;
  if (_PIC__PARTICLE_COLLISION_MODEL__SAMPLE_COLLISION_FREQUENTCY_MODE__ == _PIC_MODE_ON_) {
    SamplingData = cell->GetAssociatedDataBufferPointer() + PIC::Mesh::collectingCellSampleDataPointerOffset;
  }

  // If no particles in cell, return
  if (FirstCellParticle == -1) return;

  // Count particles per species
  int nParticleNumber[PIC::nTotalSpecies] = {0};
  for (long int ptr = FirstCellParticle; ptr != -1; ptr = PIC::ParticleBuffer::GetNext(ptr)) {
    nParticleNumber[PIC::ParticleBuffer::GetI(ptr)]++;
  }

  // Find maximum number of particles among all species
  int ParticleDataListLength = 0;
  for (int s = 0; s < PIC::nTotalSpecies; s++) {
    if (ParticleDataListLength < nParticleNumber[s]) {
      ParticleDataListLength = nParticleNumber[s];
    }
  }

  // Allocate particle lists
  cParticleDataList *s0ParticleDataList = new cParticleDataList[ParticleDataListLength];
  cParticleDataList *s1ParticleDataList = new cParticleDataList[ParticleDataListLength];

  // Loop through first species s0
  for (int s0 = 0; s0 < PIC::nTotalSpecies; s0++) {
    if (nParticleNumber[s0] == 0 || PIC::MolecularData::GetSpecieType(s0) != _PIC_SPECIE_TYPE__GAS_) {
      continue;
    }

    double m0 = PIC::MolecularData::GetMass(s0);
    double LocalTimeStep = block->GetLocalTimeStep(s0);  // Same for all species
    double LocalParticleWeight = block->GetLocalParticleWeight(s0);  // Same for all species
    
    // Populate s0 particle list
    cParticleDataList *s0List = s0ParticleDataList;
    int cnt = 0;
    for (long int ptr = FirstCellParticle; ptr != -1; ptr = PIC::ParticleBuffer::GetNext(ptr)) {
      if (PIC::ParticleBuffer::GetI(ptr) == (unsigned)s0) {
        PIC::ParticleBuffer::byte *ParticleData = PIC::ParticleBuffer::GetParticleDataPointer(ptr);
        s0List[cnt].ParticleData = ParticleData;
        s0List[cnt].v = PIC::ParticleBuffer::GetV(ParticleData);
        cnt++;
      }
    }

    // Loop through second species s1
    for (int s1 = s0; s1 < PIC::nTotalSpecies; s1++) {
      if (nParticleNumber[s1] == 0 || PIC::MolecularData::GetSpecieType(s1) != _PIC_SPECIE_TYPE__GAS_) {
        continue;
      }

      // Check if collisions should be simulated for this pair
      if (DoSimulateCellCollisions != NULL) {
        if (!DoSimulateCellCollisions(cell, s0, s1)) continue;
      }

      double m1 = PIC::MolecularData::GetMass(s1);
      
      // Populate s1 particle list (if different from s0)
      cParticleDataList *s1List = (s0 == s1) ? s0List : s1ParticleDataList;
      if (s0 != s1) {
        cnt = 0;
        for (long int ptr = FirstCellParticle; ptr != -1; ptr = PIC::ParticleBuffer::GetNext(ptr)) {
          if (PIC::ParticleBuffer::GetI(ptr) == (unsigned)s1) {
            PIC::ParticleBuffer::byte *ParticleData = PIC::ParticleBuffer::GetParticleDataPointer(ptr);
            s1List[cnt].ParticleData = ParticleData;
            s1List[cnt].v = PIC::ParticleBuffer::GetV(ParticleData);
            cnt++;
          }
        }
      }

      // Calculate maximum collision cross section
      double SigmaCrMax = 0.0;
      double SigmaCrMean = 0.0;
      
      for (int ntest = 0; ntest < SigmaCrMax_nTest; ntest++) {
        int s0ptr = (int)(rnd()*nParticleNumber[s0]);
        int s1ptr = (int)(rnd()*nParticleNumber[s1]);
        
        double cr = Vector3D::Distance(s1List[s1ptr].v,s0List[s0ptr].v);
        double SigmaCr = cr * EvaluateCollisionCrossSection(
            s0,s0List[s0ptr].v,s1,s1List[s1ptr].v,block,cell);
        
        SigmaCrMean += SigmaCr;
        if (SigmaCr > SigmaCrMax) SigmaCrMax = SigmaCr;
      }
      
      SigmaCrMax *= SigmaCrMax_SafetyMargin;
      SigmaCrMean *= SigmaCrMax_SafetyMargin / SigmaCrMax_nTest;

      // Calculate majorant collision frequency
      double MajorantFrequency;
      if (s0 == s1) {
        MajorantFrequency = 0.5 * nParticleNumber[s0] * (nParticleNumber[s0] /*-1*/) * 
                           LocalParticleWeight * SigmaCrMax / cellMeasure;
      } else {
        MajorantFrequency = nParticleNumber[s0] * nParticleNumber[s1] * 
                           LocalParticleWeight * SigmaCrMax / cellMeasure;
      }

      // Apply collision limiting if needed
      double CollisionLimitingFactor = 1.0;
      if (LocalTimeStep * MajorantFrequency > CollisionLimitingThrehold * 
          max(nParticleNumber[s0], nParticleNumber[s1])) {
        CollisionLimitingFactor = LocalTimeStep * MajorantFrequency / 
                                (CollisionLimitingThrehold * max(nParticleNumber[s0], nParticleNumber[s1]));
        MajorantFrequency /= CollisionLimitingFactor;
      }

      // Simulate collisions using Majorant Frequency scheme
      double TimeCounter = 0.0;

      while ((TimeCounter += -log(rnd())/MajorantFrequency) < LocalTimeStep) {
        // Get particle pointers
        int s0ptr = (int)(rnd() * nParticleNumber[s0]);
        int s1ptr;
        do {
          s1ptr = (int)(rnd() * nParticleNumber[s1]);
        } while (s0 == s1 && s0ptr == s1ptr);

        double *v0 = s0List[s0ptr].v;
        double *v1 = s1List[s1ptr].v;

        // Calculate relative velocity
        double cr = 0.0;
        double vrel[3],vcm[3];

        for (int idim = 0; idim < 3; idim++) {
          vrel[idim] = v1[idim] - v0[idim];
          cr += vrel[idim] * vrel[idim];
	  vcm[idim] = (m1*v1[idim] + m0*v0[idim])/(m1 + m0);
        }
        cr = sqrt(cr);

        // Check if collision occurs
        double SigmaCr = cr * EvaluateCollisionCrossSection(
            s0,v0,s1,v1,block,cell);
        if (rnd() * SigmaCrMax >= SigmaCr) continue;

	//model the internal energy exchange
        if (_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_) {
          double crInit = cr;
          bool UpdateFlag[2]={true,true};

          PIC::IDF::RedistributeEnergy(s0List[s0ptr].ParticleData,s1List[s1ptr].ParticleData, cr, UpdateFlag, cell);

          for (int idim = 0; idim < 3; idim++) {
            vrel[idim] = (cr > 1.0E-10) ? vrel[idim] * cr / crInit : 0.0;
          }
        }

        // Perform collision - directly modifies v0 and v1 which point to particle state
        if (!PIC::MolecularCollisions::VelocityScattering::Scatter(
                s0,v0,m0,s1,v1,m1,vcm,vrel))
          exit(__LINE__,__FILE__,
              "Invalid or unsupported molecular angular-scattering model");

        // Sample collision frequency if sampling is enabled
        if (_PIC__PARTICLE_COLLISION_MODEL__SAMPLE_COLLISION_FREQUENTCY_MODE__ == _PIC_MODE_ON_) {
          // Sample for first species
          int CollFreqOffset_s0 = CollsionFrequentcySampling::SamplingBufferOffset + 
                                 sizeof(double) * CollsionFrequentcySampling::Offset(s0, s1);
          *((double*)(SamplingData + CollFreqOffset_s0)) += LocalParticleWeight / 
                                                           LocalTimeStep / cellMeasure * 
                                                           CollisionLimitingFactor;

          // Sample for second species (if different from first)
          int CollFreqOffset_s1 = CollsionFrequentcySampling::SamplingBufferOffset + 
                                   sizeof(double) * CollsionFrequentcySampling::Offset(s1, s0);
          *((double*)(SamplingData + CollFreqOffset_s1)) += LocalParticleWeight / 
                                                             LocalTimeStep / cellMeasure * 
                                                             CollisionLimitingFactor;
        }
      }

      // No need to update particle buffer - velocities are modified directly
    }
  }

  // Clean up
  delete [] s1ParticleDataList;
  delete [] s0ParticleDataList;
}


void PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_mf(int i,int j,int k, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  int SigmaCrMax_nTest=100;
  double SigmaCrMax_SafetyMargin=1.05;

  struct cParticleDataList {
    double vel[3];
    bool ValueChangedFlag;
    PIC::ParticleBuffer::byte *ParticleData;
  };


  int s,s0,s1;

  int nParticleNumber[PIC::nTotalSpecies],nMaxSpecParticleNumber,cnt;
  //  long int FirstCellParticleTable[_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_],
  long int FirstCellParticle,ptr;
  double cellMeasure;
  PIC::Mesh::cDataBlockAMR *block;

  //particle lists
  int ParticleDataListLength=0;
  cParticleDataList *s0ParticleDataList; //=new cParticleDataList[ParticleDataListLength];
  cParticleDataList *s1ParticleDataList; //=new cParticleDataList[ParticleDataListLength];  //s0ParticleDataList,s1ParticleDataList are used to store the actual data
  cParticleDataList *s0List=NULL,*s1List=NULL; //s0List and s1List are used to access s0ParticleDataList and s1ParticleDataList

  //sampling data
  PIC::Mesh::cDataCenterNode *cell;
  char *SamplingData;

  //global offsets
  int LoadBalancingMeasureOffset=PIC::Mesh::cDataBlockAMR_static_data::LoadBalancingMeasureOffset;

  //simulate particle's collisions
  block=node->block;
  if (block==NULL) return;

  FirstCellParticle=block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
  cell=block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k));

  if (cell==NULL) return;
  cellMeasure=cell->Measure;

  if (_PIC__PARTICLE_COLLISION_MODEL__SAMPLE_COLLISION_FREQUENTCY_MODE__ == _PIC_MODE_ON_) {
    SamplingData=cell->GetAssociatedDataBufferPointer()+PIC::Mesh::collectingCellSampleDataPointerOffset;
  }

  if (FirstCellParticle!=-1) {
    //simulate collision in the cell
    //1. coalcualte the number of the model particles in the cell and re-allocate the particle lists if needed
    for (s=0;s<PIC::nTotalSpecies;s++) nParticleNumber[s]=0;
    for (ptr=FirstCellParticle;ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) nParticleNumber[PIC::ParticleBuffer::GetI(ptr)]++;
    for (s=0;s<PIC::nTotalSpecies;s++) if (ParticleDataListLength<nParticleNumber[s]) ParticleDataListLength=nParticleNumber[s];

    s0ParticleDataList=new cParticleDataList[ParticleDataListLength];
    s1ParticleDataList=new cParticleDataList[ParticleDataListLength];

    //loop through the first species
    double m0,m1,am;
    double LocalParticleWeight_s0,LocalParticleWeight_s1;
    double LocalTimeStep_s0,LocalTimeStep_s1;
    double minParticleWeightCorrection_s0,minParticleWeightCorrection_s1,sumWeightCorrection_s0,sumWeightCorrection_s1;
    PIC::ParticleBuffer::byte *ParticleData;

    for (s0=0;s0<PIC::nTotalSpecies;s0++) if ((PIC::MolecularData::GetSpecieType(s0)==_PIC_SPECIE_TYPE__GAS_)&&(nParticleNumber[s0]!=0)) {

    m0=PIC::MolecularData::GetMass(s0);
    LocalParticleWeight_s0=block->GetLocalParticleWeight(s0);
    LocalTimeStep_s0=block->GetLocalTimeStep(s0);
    minParticleWeightCorrection_s0=1,sumWeightCorrection_s0=0.0;

    //populate the particle list
    s0List=s0ParticleDataList;

    for (cnt=0,ptr=FirstCellParticle;ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) if (PIC::ParticleBuffer::GetI(ptr)==(unsigned)s0) {
      ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);

      PIC::ParticleBuffer::GetV(s0ParticleDataList[cnt].vel,ParticleData);
      s0ParticleDataList[cnt].ParticleData=ParticleData;
      s0ParticleDataList[cnt].ValueChangedFlag=false;

      if (_INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_) {
        double wc = PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);

        if (cnt == 0) {
          minParticleWeightCorrection_s0 = wc;
        } else if (wc < minParticleWeightCorrection_s0) {
          minParticleWeightCorrection_s0 = wc;
        }

        sumWeightCorrection_s0 += wc;
      } else {
        sumWeightCorrection_s0 += 1.0;
      }

      cnt++;
    }

    //loop through the second speices
    for (s1=s0;s1<PIC::nTotalSpecies;s1++) if ((PIC::MolecularData::GetSpecieType(s1)==_PIC_SPECIE_TYPE__GAS_)&&(nParticleNumber[s1]!=0)) {
      //call a user-defined funstion to verify whether collision should be modeled in the cell

      if (DoSimulateCellCollisions!=NULL) {
        if (DoSimulateCellCollisions(cell,s0,s1)==false) continue;
      }

      m1=PIC::MolecularData::GetMass(s1);
      am=m0+m1;
      LocalParticleWeight_s1=block->GetLocalParticleWeight(s1);
      LocalTimeStep_s1=block->GetLocalTimeStep(s1);
      sumWeightCorrection_s1=0.0;

      //populate the list
      if (s0==s1) {
        s1List=s0List,minParticleWeightCorrection_s1=minParticleWeightCorrection_s0,sumWeightCorrection_s1=sumWeightCorrection_s0;
      }
      else {
        s1List=s1ParticleDataList;

        for (cnt=0,ptr=FirstCellParticle;ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) if (PIC::ParticleBuffer::GetI(ptr)==(unsigned)s1) {
          ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);

          PIC::ParticleBuffer::GetV(s1ParticleDataList[cnt].vel,ParticleData);
          s1ParticleDataList[cnt].ParticleData=ParticleData;
          s1ParticleDataList[cnt].ValueChangedFlag=false;

          if (_INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_) {
            double wc = PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);

            if (cnt == 0) {
              minParticleWeightCorrection_s1 = wc;
            }
            else if (wc < minParticleWeightCorrection_s1) {
              minParticleWeightCorrection_s1 = wc;
            }

            sumWeightCorrection_s1 += wc;
          }
          else {
            sumWeightCorrection_s1 += 1.0;
          }

          cnt++;
        }
      }

      //simulate collsions between the pair of species
      double SigmaCrMax=0.0,SigmaCr,TimeCounter,SigmaCrMean=0.0;
      double v0[3],v1[3],cr;

      //1.Evaluate the maximum value of SigmaCr
      for (int ntest=0;ntest<SigmaCrMax_nTest;ntest++) {
        memcpy(v0,s0List[(int)(rnd()*nParticleNumber[s0])].vel,3*sizeof(double));
        memcpy(v1,s1List[(int)(rnd()*nParticleNumber[s1])].vel,3*sizeof(double));
        cr=sqrt(pow(v1[0]-v0[0],2)+pow(v1[1]-v0[1],2)+pow(v1[2]-v0[2],2));

        SigmaCr=cr*EvaluateCollisionCrossSection(s0,v0,s1,v1,block,cell);

        if (SigmaCr>SigmaCrMax) SigmaCrMax=SigmaCr;
          SigmaCrMean+=SigmaCr;
        }

        SigmaCrMax*=SigmaCrMax_SafetyMargin;
        SigmaCrMean*=SigmaCrMax_SafetyMargin/SigmaCrMax_nTest;

        //2.Evaluate the prospective number of collisions
        double maxLocalTimeStep,minLocalParticleWeight;
        double MajorantFrequency;

        maxLocalTimeStep=max(LocalTimeStep_s0,LocalTimeStep_s1);
        minLocalParticleWeight=min(LocalParticleWeight_s0*minParticleWeightCorrection_s0,LocalParticleWeight_s1*minParticleWeightCorrection_s1);

        if (s0==s1) {
          MajorantFrequency=0.5*sumWeightCorrection_s0/minParticleWeightCorrection_s0*
           (sumWeightCorrection_s0/minParticleWeightCorrection_s0 /*-1.0*/)*
            LocalParticleWeight_s0*minParticleWeightCorrection_s0*SigmaCrMax/cellMeasure;
        }
        else {
          MajorantFrequency=(sumWeightCorrection_s0*LocalParticleWeight_s0)*(sumWeightCorrection_s1*LocalParticleWeight_s1)*
            SigmaCrMax/minLocalParticleWeight/cellMeasure;
        }

        if (MajorantFrequency*maxLocalTimeStep<1.0E-5) continue;

        //3. Collision Limiting
        double CollisionLimitingFactor=1.0;

        // Computation of collision limiting factor
        if (CollisionLimitingThrehold*max(nParticleNumber[s0],nParticleNumber[s1])<MajorantFrequency/SigmaCrMax*SigmaCrMean*maxLocalTimeStep) {
          CollisionLimitingFactor=(MajorantFrequency/SigmaCrMax*SigmaCrMean*maxLocalTimeStep)/(CollisionLimitingThrehold*max(nParticleNumber[s0],nParticleNumber[s1]));
          MajorantFrequency/=CollisionLimitingFactor;
        }


        //4. simulate collisions
        int s0ptr,s1ptr,idim;
        double vrel[3],vcm[3];

        //set the initial value of the counter
        TimeCounter=0.0;

        while ((TimeCounter+=-log(rnd())/MajorantFrequency)<maxLocalTimeStep) {
          s0ptr=(int)((int)(rnd()*nParticleNumber[s0]));

          if (s0!=s1) s1ptr=(int)((int)(rnd()*nParticleNumber[s1]));
          else {
            do {
              s1ptr=(int)((int)(rnd()*nParticleNumber[s1]));
            }
            while (s0ptr==s1ptr);
          }

          memcpy(v0,s0List[s0ptr].vel,3*sizeof(double));
          memcpy(v1,s1List[s1ptr].vel,3*sizeof(double));

          for (cr=0.0,idim=0;idim<3;idim++) {
            vrel[idim]=v1[idim]-v0[idim];
            cr+=vrel[idim]*vrel[idim];
	    vcm[idim] = (m1*v1[idim] + m0*v0[idim])/(m1 + m0);
          }

          cr=sqrt(cr);
          SigmaCr=cr*EvaluateCollisionCrossSection(s0,v0,s1,v1,block,cell);

          if (rnd()*SigmaCrMax>=SigmaCr) continue;

          //determine weather the properties of particles were updated
          double pUpdate_s0,pUpdate_s1;
          bool UpdateFlag[2];

          pUpdate_s0=minLocalParticleWeight/LocalParticleWeight_s0 * LocalTimeStep_s0/maxLocalTimeStep;
          pUpdate_s1=minLocalParticleWeight/LocalParticleWeight_s1 * LocalTimeStep_s1/maxLocalTimeStep;

          if (_INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_) {
            pUpdate_s0 /= PIC::ParticleBuffer::GetIndividualStatWeightCorrection(s0List[s0ptr].ParticleData);
            pUpdate_s1 /= PIC::ParticleBuffer::GetIndividualStatWeightCorrection(s1List[s1ptr].ParticleData);
          }

          UpdateFlag[0]=(rnd()<pUpdate_s0) ? true :false;
          UpdateFlag[1]=(rnd()<pUpdate_s1) ? true :false;

          //model the internal energy exchange
          if (_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_) {
            double crInit = cr;

            PIC::IDF::RedistributeEnergy(s0List[s0ptr].ParticleData,s1List[s1ptr].ParticleData,cr,UpdateFlag,cell);

            for (int idim = 0; idim < 3; idim++) {
              vrel[idim] = (cr > 1.0E-10) ? vrel[idim] * cr / crInit : 0.0;
            }
          }

          //the collision is considered to be true
          if (!PIC::MolecularCollisions::VelocityScattering::Scatter(
                  s0,v0,m0,s1,v1,m1,vcm,vrel))
            exit(__LINE__,__FILE__,
                "Invalid or unsupported molecular angular-scattering model");

          //update the velocities in the lists
          if (UpdateFlag[0]==true) {
            s0List[s0ptr].ValueChangedFlag=true;
            memcpy(s0List[s0ptr].vel,v0,3*sizeof(double));

            if (_PIC__PARTICLE_COLLISION_MODEL__SAMPLE_COLLISION_FREQUENTCY_MODE__ == _PIC_MODE_ON_) {
              int CollFreqOffset=CollsionFrequentcySampling::SamplingBufferOffset+sizeof(double)*CollsionFrequentcySampling::Offset(s0,s1);

              *((double*)(SamplingData+CollFreqOffset))+=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(s0List[s0ptr].ParticleData)*LocalParticleWeight_s0/LocalTimeStep_s0/cellMeasure*CollisionLimitingFactor; // calculate collision frequency taking into account the collision limiting factor
            }
         }

         if (UpdateFlag[1]==true) {
            s1List[s1ptr].ValueChangedFlag=true;
            memcpy(s1List[s1ptr].vel,v1,3*sizeof(double));

            if (_PIC__PARTICLE_COLLISION_MODEL__SAMPLE_COLLISION_FREQUENTCY_MODE__ == _PIC_MODE_ON_) {
              int CollFreqOffset = CollsionFrequentcySampling::SamplingBufferOffset 
                 + sizeof(double)*CollsionFrequentcySampling::Offset(s1, s0);
 
 
              *((double*)(SamplingData+CollFreqOffset)) 
                += PIC::ParticleBuffer::GetIndividualStatWeightCorrection(s1List[s1ptr].ParticleData)
                * LocalParticleWeight_s1 / LocalTimeStep_s1 / cellMeasure 
                * CollisionLimitingFactor; // calculate collision frequency taking into account the collision limiting factor
            }
          }
        }

        //update the velocities of the species 's1'
        if (s0!=s1) {
          for (cnt=0;cnt<nParticleNumber[s1];cnt++) if (s1List[cnt].ValueChangedFlag==true) {
            PIC::ParticleBuffer::SetV(s1List[cnt].vel,s1List[cnt].ParticleData);
          }
        }
      }

      //update velocities of species 's0'
      for (cnt=0;cnt<nParticleNumber[s0];cnt++) if (s0List[cnt].ValueChangedFlag==true) {
        PIC::ParticleBuffer::SetV(s0List[cnt].vel,s0List[cnt].ParticleData);
      }
    }

    delete [] s0ParticleDataList;
    delete [] s1ParticleDataList;
  } //end simulation of collisons in a cell
}

void PIC::MolecularCollisions::ParticleCollisionModel::mf() {
  CollisionModelManager(ModelCellCollisions_mf);
}




//Majorant Frequency collision procedure
void PIC::MolecularCollisions::ParticleCollisionModel::mf_Yinsi() {
  CollisionModelManager(ModelCellCollisions_mf_Yinsi);
}


void PIC::MolecularCollisions::ParticleCollisionModel::ModelCellCollisions_mf_Yinsi(int i,int j,int k, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  int SigmaCrMax_nTest=100;
  double SigmaCrMax_SafetyMargin=1.05;

  // CollisionLimitingThrehold is shared configuration for all selectors.
  // Earlier code overwrote it with 40 on every Yinsi cell call, making a
  // later NTC/MF result depend on selector execution order.  The Yinsi
  // equations below already consume the common dimensionless limit; do not
  // mutate process-global configuration from a cell-local physics kernel.

  struct cParticleDataList {
    double vel[3];
    bool ValueChangedFlag;
    PIC::ParticleBuffer::byte *ParticleData;
  };


  int s,s0,s1;
  int nParticleNumber[PIC::nTotalSpecies],nMaxSpecParticleNumber,cnt;
  long int FirstCellParticle,ptr;
  double cellMeasure;
  int thread;

  //particle lists
  int ParticleDataListLength=0;
  cParticleDataList *s0ParticleDataList=NULL; //=new cParticleDataList[ParticleDataListLength];
  cParticleDataList *s1ParticleDataList=NULL; //=new cParticleDataList[ParticleDataListLength];  //s0ParticleDataList,s1ParticleDataList are used to store the actual data
  cParticleDataList *s0List=NULL,*s1List=NULL; //s0List and s1List are used to access s0ParticleDataList and s1ParticleDataList

  //sampling data
  PIC::Mesh::cDataCenterNode *cell;
  char *SamplingData;
  PIC::Mesh::cDataBlockAMR *block;

  block=node->block;
  if (block==NULL) return;

  cell=block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k));
  if (cell==NULL) return;

  FirstCellParticle=block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

  cell=block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k));
  cellMeasure=cell->Measure;

  if (_PIC__PARTICLE_COLLISION_MODEL__SAMPLE_COLLISION_FREQUENTCY_MODE__ == _PIC_MODE_ON_) {
    SamplingData=cell->GetAssociatedDataBufferPointer()+PIC::Mesh::collectingCellSampleDataPointerOffset;
  }

  if (FirstCellParticle!=-1) {
    //simulate collision in the cell
    //1. coalcualte the number of the model particles in the cell and re-allocate the particle lists if needed
    for (s=0;s<PIC::nTotalSpecies;s++) nParticleNumber[s]=0;
    for (ptr=FirstCellParticle;ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) nParticleNumber[PIC::ParticleBuffer::GetI(ptr)]++;
    for (ParticleDataListLength=0,s=0;s<PIC::nTotalSpecies;s++) if (ParticleDataListLength<nParticleNumber[s]) ParticleDataListLength=nParticleNumber[s];

    s0ParticleDataList=new cParticleDataList[ParticleDataListLength];
    s1ParticleDataList=new cParticleDataList[ParticleDataListLength];

    //loop through the first species
    double m0,m1,am;
    double LocalParticleWeight_s0,LocalParticleWeight_s1;
    double LocalTimeStep_s0,LocalTimeStep_s1;
    double minParticleWeightCorrection_s0,minParticleWeightCorrection_s1,sumWeightCorrection_s0,sumWeightCorrection_s1;
    PIC::ParticleBuffer::byte *ParticleData;

    for (s0=0;s0<PIC::nTotalSpecies;s0++) if ((PIC::MolecularData::GetSpecieType(s0)==_PIC_SPECIE_TYPE__GAS_)&&(nParticleNumber[s0]!=0)) {

    m0=PIC::MolecularData::GetMass(s0);
    LocalParticleWeight_s0=block->GetLocalParticleWeight(s0);
    LocalTimeStep_s0=block->GetLocalTimeStep(s0);
    minParticleWeightCorrection_s0=1,sumWeightCorrection_s0=0.0;

    //populate the particle list
    s0List=s0ParticleDataList;

    for (cnt=0,ptr=FirstCellParticle;ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) if (PIC::ParticleBuffer::GetI(ptr)==(unsigned)s0) {
      ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);

      PIC::ParticleBuffer::GetV(s0ParticleDataList[cnt].vel,ParticleData);
      s0ParticleDataList[cnt].ParticleData=ParticleData;
      s0ParticleDataList[cnt].ValueChangedFlag=false;

      if (_INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_) {
        double wc = PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);

        if (cnt == 0) {
          minParticleWeightCorrection_s0 = wc;
        } else if (wc < minParticleWeightCorrection_s0) {
          minParticleWeightCorrection_s0 = wc;
        }

        sumWeightCorrection_s0 += wc;
      } else {
        sumWeightCorrection_s0 += 1.0;
      }

      cnt++;
    }

    //loop through the second speices
    for (s1=s0;s1<PIC::nTotalSpecies;s1++) if ((PIC::MolecularData::GetSpecieType(s1)==_PIC_SPECIE_TYPE__GAS_)&&(nParticleNumber[s1]!=0)) {
      //call a user-defined function to verify whether collision should be modeled in the cell

      if (DoSimulateCellCollisions!=NULL) {
        if (DoSimulateCellCollisions(cell,s0,s1)==false) continue;
      }

      m1=PIC::MolecularData::GetMass(s1);
      am=m0+m1;
      LocalParticleWeight_s1=block->GetLocalParticleWeight(s1);
      LocalTimeStep_s1=block->GetLocalTimeStep(s1);
      sumWeightCorrection_s1=0.0;

      //populate the list
      if (s0==s1) s1List=s0List,minParticleWeightCorrection_s1=minParticleWeightCorrection_s0,sumWeightCorrection_s1=sumWeightCorrection_s0;
      else {
        s1List=s1ParticleDataList;

        for (cnt=0,ptr=FirstCellParticle;ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) if (PIC::ParticleBuffer::GetI(ptr)==(unsigned)s1) {
          ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);

          PIC::ParticleBuffer::GetV(s1ParticleDataList[cnt].vel,ParticleData);
          s1ParticleDataList[cnt].ParticleData=ParticleData;
          s1ParticleDataList[cnt].ValueChangedFlag=false;

          if (_INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_) {
            double wc=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);

            if (cnt==0) {
              minParticleWeightCorrection_s1=wc;
            } else if (wc<minParticleWeightCorrection_s1) {
              minParticleWeightCorrection_s1=wc;
            }

            sumWeightCorrection_s1+=wc;
          } else {
            sumWeightCorrection_s1+=1.0;
          }

          cnt++;
        }
      }

      //simulate collsions between the pair of species
      double SigmaCrMax=0.0,SigmaCr,TimeCounter,SigmaCrMean=0.0;
      double v0[3],v1[3],cr;

      //1.Evaluate the maximum value of SigmaCr
      for (int ntest=0;ntest<SigmaCrMax_nTest;ntest++) {
        memcpy(v0,s0List[(int)(rnd()*nParticleNumber[s0])].vel,3*sizeof(double));
        memcpy(v1,s1List[(int)(rnd()*nParticleNumber[s1])].vel,3*sizeof(double));
        cr=sqrt(pow(v1[0]-v0[0],2)+pow(v1[1]-v0[1],2)+pow(v1[2]-v0[2],2));

        SigmaCr=cr*EvaluateCollisionCrossSection(s0,v0,s1,v1,block,cell);

        if (SigmaCr>SigmaCrMax) SigmaCrMax=SigmaCr;
          SigmaCrMean+=SigmaCr;
        }

        SigmaCrMax*=SigmaCrMax_SafetyMargin;
        SigmaCrMean*=SigmaCrMax_SafetyMargin/SigmaCrMax_nTest;

        //2.Evaluate the prospective number of collisions
        double maxLocalTimeStep,minLocalParticleWeight[2];
        double MajorantFrequency[2];

        maxLocalTimeStep=max(LocalTimeStep_s0,LocalTimeStep_s1);
        minLocalParticleWeight[0]=LocalParticleWeight_s0*minParticleWeightCorrection_s0;
        minLocalParticleWeight[1]=LocalParticleWeight_s1*minParticleWeightCorrection_s1;

        if (s0==s1) {
          MajorantFrequency[0]=0.5*(nParticleNumber[s0] /*-1*/)*
             (sumWeightCorrection_s0*LocalParticleWeight_s0)*SigmaCrMax/cellMeasure;
          //N_a* (sum weight)/cellVol*vel
          // 0.5* n0* (n0-1)* w0 *sigma*relV/ cellV
          MajorantFrequency[1]=MajorantFrequency[0];
        }
        else {
          MajorantFrequency[0]=nParticleNumber[s0]*(sumWeightCorrection_s1*LocalParticleWeight_s1)*SigmaCrMax/cellMeasure;
          // n0*n1*w1*sig*rel/cellV
          // no*(sum weight)*rel/CellV
          // n0=sum(weightCorrection0)/minWeightCorr0  
          //one particle undergoes one collision,change v once equivalent to
          //two particles with half weight undergoes two collisions, each changes v once
          //two particles with 0.1 and 0.9 weight correction, each should undergo v change once

          MajorantFrequency[1]=nParticleNumber[s1]*(sumWeightCorrection_s0*LocalParticleWeight_s0)*SigmaCrMax/cellMeasure;
          // n1*(n0*w0*sig*rel/cellV)
          //(n0*w0*sig*rel/cellV) =1/s=same as analytical
          // n1=sum(weightCorrection1)/minWeightCorr1
        }

        //3. Collision Limiting
        double CollisionLimitingFactor[2]={1.0,1.0};

        // Computation of collision limiting factor
        if (maxLocalTimeStep*MajorantFrequency[0]>CollisionLimitingThrehold*nParticleNumber[s0]){
          //printf("maxLocalTimeStep:%e,MajorantFrequency:%e,prod:%e, s0:%d, s1:%d\n",maxLocalTimeStep,MajorantFrequency,maxLocalTimeStep*MajorantFrequency,s0,s1);
          CollisionLimitingFactor[0] = maxLocalTimeStep*MajorantFrequency[0]/CollisionLimitingThrehold/nParticleNumber[s0];
          MajorantFrequency[0] /=CollisionLimitingFactor[0]; 
        }

        if (maxLocalTimeStep*MajorantFrequency[1]>CollisionLimitingThrehold*nParticleNumber[s1]){
          //printf("maxLocalTimeStep:%e,MajorantFrequency:%e,prod:%e, s0:%d, s1:%d\n",maxLocalTimeStep,MajorantFrequency,maxLocalTimeStep*MajorantFrequency,s0,s1);
          CollisionLimitingFactor[1] = maxLocalTimeStep*MajorantFrequency[1]/CollisionLimitingThrehold/nParticleNumber[s1];
          MajorantFrequency[1] /=CollisionLimitingFactor[1]; 
        }

        //4. simulate collisions
        int s0ptr,s1ptr,idim;
        double vrel[3],vcm[3];

        //set the initial value of the counter
        TimeCounter=0.0;
        double MajorFreq=max(MajorantFrequency[0],MajorantFrequency[1]);
 
        while ((TimeCounter+=-log(rnd())/MajorFreq)<maxLocalTimeStep) {
          s0ptr=(int)((int)(rnd()*nParticleNumber[s0]));

          if (s0!=s1) s1ptr=(int)((int)(rnd()*nParticleNumber[s1]));
          else {
            do {
              s1ptr=(int)((int)(rnd()*nParticleNumber[s1]));
            }
            while (s0ptr==s1ptr);
          }

          memcpy(v0,s0List[s0ptr].vel,3*sizeof(double));
          memcpy(v1,s1List[s1ptr].vel,3*sizeof(double));

          for (cr=0.0,idim=0;idim<3;idim++) {
            vrel[idim]=v1[idim]-v0[idim];
            cr+=vrel[idim]*vrel[idim];
	    vcm[idim] = (m1*v1[idim] + m0*v0[idim])/(m1 + m0);
          }

          cr=sqrt(cr);
          SigmaCr=cr*EvaluateCollisionCrossSection(s0,v0,s1,v1,block,cell);

          if (rnd()*SigmaCrMax>=SigmaCr) continue;

          //determine weather the properties of particles were updated
          double pUpdate_s0,pUpdate_s1;
          bool UpdateFlag[2];

          pUpdate_s0=MajorantFrequency[0]/MajorFreq * LocalTimeStep_s0/maxLocalTimeStep;
          pUpdate_s1=MajorantFrequency[1]/MajorFreq * LocalTimeStep_s1/maxLocalTimeStep;

          if (_INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_) {
            exit(__LINE__,__FILE__,"Error: not implemented");
            //pUpdate_s0/=meanWeightCorr[0]/minParticleWeightCorrection_s0;//min weight correction increases the majorant frequency
            //pUpdate_s1/=meanWeightCorr[1]/minParticleWeightCorrection_s1;
          }
              
          UpdateFlag[0]=(rnd()<pUpdate_s0) ? true :false;
          UpdateFlag[1]=(rnd()<pUpdate_s1) ? true :false;

          //model the internal energy exchange
          if (_PIC_INTERNAL_DEGREES_OF_FREEDOM_MODE_ == _PIC_MODE_ON_) {
            double crInit = cr;

            PIC::IDF::RedistributeEnergy(s0List[s0ptr].ParticleData,s1List[s1ptr].ParticleData, cr, UpdateFlag, cell);
          
            for (int idim = 0; idim < 3; idim++) {
              vrel[idim] = (cr > 1.0E-10) ? vrel[idim] * cr / crInit : 0.0;
            }
          }

          //the collision is considered to be true
          if (!PIC::MolecularCollisions::VelocityScattering::Scatter(
                  s0,v0,m0,s1,v1,m1,vcm,vrel))
            exit(__LINE__,__FILE__,
                "Invalid or unsupported molecular angular-scattering model");

          //update the velocities in the lists
          if (UpdateFlag[0]==true) {
            s0List[s0ptr].ValueChangedFlag=true;
            memcpy(s0List[s0ptr].vel,v0,3*sizeof(double));

            if (_PIC__PARTICLE_COLLISION_MODEL__SAMPLE_COLLISION_FREQUENTCY_MODE__ == _PIC_MODE_ON_) {
              int CollFreqOffset=CollsionFrequentcySampling::SamplingBufferOffset+sizeof(double)*CollsionFrequentcySampling::Offset(s0,s1);

              *((double*)(SamplingData+CollFreqOffset))+=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(s0List[s0ptr].ParticleData)*LocalParticleWeight_s0/LocalTimeStep_s0/cellMeasure*CollisionLimitingFactor[0]; // calculate collision frequency taking into account the collision limiting factor
            }
          }

          if (UpdateFlag[1]==true) {
            s1List[s1ptr].ValueChangedFlag=true;
            memcpy(s1List[s1ptr].vel,v1,3*sizeof(double));

            if (_PIC__PARTICLE_COLLISION_MODEL__SAMPLE_COLLISION_FREQUENTCY_MODE__ == _PIC_MODE_ON_) {
              int CollFreqOffset = 
                CollsionFrequentcySampling::SamplingBufferOffset + 
                sizeof(double) * CollsionFrequentcySampling::Offset(s1, s0);

              *((double*)(SamplingData + CollFreqOffset)) += 
                 PIC::ParticleBuffer::GetIndividualStatWeightCorrection(s1List[s1ptr].ParticleData) * 
                 LocalParticleWeight_s1 / LocalTimeStep_s1 / cellMeasure * 
                 CollisionLimitingFactor[1]; // calculate collision frequency taking into account the collision limiting factor
            } //_PIC__PARTICLE_COLLISION_MODEL__SAMPLE_COLLISION_FREQUENTCY_MODE__
          }
        }

        //update the velocities of the species 's1'
        if (s0!=s1) {
          for (cnt=0;cnt<nParticleNumber[s1];cnt++) if (s1List[cnt].ValueChangedFlag==true) {
            PIC::ParticleBuffer::SetV(s1List[cnt].vel,s1List[cnt].ParticleData);
          }
        }
      }//for (s1=s0;s1<PIC::nTotalSpecies;s1++) 

      //update velocities of species 's0'
      for (cnt=0;cnt<nParticleNumber[s0];cnt++) if (s0List[cnt].ValueChangedFlag==true) {
        PIC::ParticleBuffer::SetV(s0List[cnt].vel,s0List[cnt].ParticleData);
      }
  
    }
  
    delete [] s1ParticleDataList;
    delete [] s0ParticleDataList; 
  } //end simulation of collisons in a cell
}
 
