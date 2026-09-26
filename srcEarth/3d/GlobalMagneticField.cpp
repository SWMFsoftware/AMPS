//======================================================================================
// GlobalMagneticField.cpp
//======================================================================================
//
// Compact global cell-centered B/E storage for Mode3D backward trajectory tracing.
// Step 9 additionally retains u for an SWMF snapshot so the frozen live B/u state can
// be exported and replayed without reconstructing either quantity.
//
// Only owner-rank interior cells are packed.  MPI_Allreduce then creates identical
// compact arrays on every process.  The AMR tree itself is already globally replicated
// by AMPS, so node->Temp_ID plus local interior cell indices provide a complete global
// address without allocating remote cDataBlockAMR objects or reconstructing ghost-cell
// buffers.  Field evaluation is performed with cRowStencil, whose entries identify the
// physical owning leaf and interior (i,j,k) even for unallocated remote blocks.
//======================================================================================

#include "GlobalMagneticField.h"
#include "../util/SWMFSnapshotContract.h"

#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <mpi.h>

namespace Earth {
namespace Mode3D {
namespace GlobalMagneticField {
namespace {

// Compact arrays replicated on every MPI rank.  They are assembled before entering
// the parallel trajectory calculation and remain read-only while worker threads are
// active, so interpolation requires no mutex.
std::vector<double> GlobalMagneticField_;
std::vector<double> GlobalElectricField_;
// Step 9 retains the exact bulk velocity used to derive E.  Earlier code kept only
// B and E, which was sufficient for tracing but could not prove that a standalone
// replay used the same coupled B/u state.  This array is populated only when a plasma-
// velocity offset is supplied and remains immutable with the other compact arrays.
std::vector<double> GlobalPlasmaVelocity_;
std::vector<int> GlobalCellPresence_;
long int GlobalUsedLeafBlocks_=0;
long int GlobalInteriorCellCount_=0;
bool GlobalFieldsReady_=false;
Earth::Field::SnapshotMetadata GlobalSnapshotMetadata_;
unsigned long long GlobalSnapshotGeneration_=0;
bool GlobalPlasmaVelocityAvailable_=false;
std::string GlobalContentFingerprint_;
std::string GlobalMeshRevision_;
double GlobalSourceSimulationTime_s_=0.0;

// A product batch holds a logical lease on the published generation.  The SWMF
// coupler invokes PT callbacks serially, so later receives remain queued by the
// scheduler until the callback returns; this explicit lease also prevents any AMPS
// helper from clearing/reassembling the arrays from inside the callback.
bool GlobalFrozenBatchActive_=false;
std::string GlobalFrozenBatchSnapshotId_;
unsigned long long GlobalFrozenBatchGeneration_=0;

void RequireNoFrozenBatch_(const char* operation) {
  if (!GlobalFrozenBatchActive_) return;
  std::ostringstream message;
  message << "Cannot " << (operation!=NULL ? operation : "modify compact fields")
          << " while SWMF snapshot '" << GlobalFrozenBatchSnapshotId_
          << "' is frozen for a product batch; the next receive must remain queued.";
  throw std::runtime_error(message.str());
}

std::string SafeTag_(const char* diagnosticTag) {
  return (diagnosticTag!=NULL && diagnosticTag[0]!='\0') ?
         std::string(diagnosticTag) : std::string("Mode3D::GlobalMagneticField");
}

long int InteriorCellsPerBlock_() {
  return static_cast<long int>(_BLOCK_CELLS_X_) *
         static_cast<long int>(_BLOCK_CELLS_Y_) *
         static_cast<long int>(_BLOCK_CELLS_Z_);
}

long int InteriorCellIndex_(int i,int j,int k) {
  return static_cast<long int>(i) +
         static_cast<long int>(_BLOCK_CELLS_X_) *
         (static_cast<long int>(j) +
          static_cast<long int>(_BLOCK_CELLS_Y_) * static_cast<long int>(k));
}

long int GlobalCellIndex_(cAMRNode* node,int i,int j,int k) {
  return node->Temp_ID*InteriorCellsPerBlock_()+InteriorCellIndex_(i,j,k);
}

// Temp_ID is scratch storage used by several AMPS algorithms.  Reset every tree node
// before assigning field-array IDs so stale values from a previous mesh operation or
// snapshot can never alias a valid global-field row.  Resetting non-leaf and unused
// nodes is important because a row-stencil consistency error must be detected rather
// than accidentally indexing an old field location.
void ResetTreeTempIds_(cAMRNode* node) {
  if (node==NULL) return;

  node->Temp_ID=-1;
  if (node->block!=NULL) node->block->Temp_ID=-1;

  if (node->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) {
    for (int i=0;i<(1<<_MESH_DIMENSION_);i++) {
      ResetTreeTempIds_(node->downNode[i]);
    }
  }
}

// Assign deterministic dense IDs by traversing the globally replicated tree in fixed
// child order.  Because every rank has the same tree topology, a given physical leaf
// obtains the same Temp_ID on every rank without communication.
void AssignGlobalLeafTempIds_(cAMRNode* node,long int& nUsedLeafBlocks) {
  if (node==NULL) return;

  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if (node->IsUsedInCalculationFlag==true) {
      node->Temp_ID=nUsedLeafBlocks++;
      if (node->block!=NULL) node->block->Temp_ID=node->Temp_ID;
    }

    return;
  }

  for (int i=0;i<(1<<_MESH_DIMENSION_);i++) {
    AssignGlobalLeafTempIds_(node->downNode[i],nUsedLeafBlocks);
  }
}

void CollectUsedLeafNodes_(cAMRNode* node,std::vector<cAMRNode*>& nodes) {
  if (node==NULL) return;

  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if ((node->IsUsedInCalculationFlag==true) && (node->Temp_ID>=0)) {
      nodes.push_back(node);
    }
    return;
  }

  for (int i=0;i<(1<<_MESH_DIMENSION_);i++) {
    CollectUsedLeafNodes_(node->downNode[i],nodes);
  }
}

// MPI count parameters are int in the interface used by AMPS.  Reduce large arrays in
// bounded chunks and use MPI_IN_PLACE so no second full-sized receive array is needed.
void AllreduceDoubleVectorInPlace_(std::vector<double>& data) {
  const long int n=static_cast<long int>(data.size());
  const long int chunkMax=100000000;

  for (long int offset=0;offset<n;offset+=chunkMax) {
    const int chunk=static_cast<int>(std::min(chunkMax,n-offset));
    MPI_Allreduce(MPI_IN_PLACE,&data[static_cast<size_t>(offset)],chunk,
                  MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  }
}

void AllreduceIntVectorInPlace_(std::vector<int>& data) {
  const long int n=static_cast<long int>(data.size());
  const long int chunkMax=100000000;

  for (long int offset=0;offset<n;offset+=chunkMax) {
    const int chunk=static_cast<int>(std::min(chunkMax,n-offset));
    MPI_Allreduce(MPI_IN_PLACE,&data[static_cast<size_t>(offset)],chunk,
                  MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  }
}

// A rank-zero filesystem failure must become the same exception on every rank. If
// only rank zero threw, its peers could enter the next collective or trajectory batch
// and hang indefinitely. The text is broadcast as well as the status so the common
// failure retains the useful underlying I/O diagnostic.
void RequireCollectiveRootWriteSuccess_(int rootSucceeded,
                                        std::string rootError,
                                        const char* context) {
  MPI_Bcast(&rootSucceeded,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);
  int messageLength=(PIC::ThisThread==0) ?
      static_cast<int>(rootError.size()) : 0;
  MPI_Bcast(&messageLength,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);
  if (messageLength<0)
    throw std::runtime_error("Invalid collective SWMF write error length");
  if (PIC::ThisThread!=0)
    rootError.assign(static_cast<std::size_t>(messageLength),'\0');
  if (messageLength>0)
    MPI_Bcast(&rootError[0],messageLength,MPI_CHAR,0,MPI_GLOBAL_COMMUNICATOR);
  if (!rootSucceeded) {
    std::ostringstream message;
    message << (context!=NULL ? context : "rank-zero SWMF write") << " failed";
    if (!rootError.empty()) message << ": " << rootError;
    throw std::runtime_error(message.str());
  }
}

// Pack only authoritative owner-rank interior cells.  Nonowner blocks can be present
// because AMPS keeps a local neighbor layer, but those copies may contain stale ghost
// data and must never contribute to the global snapshot.
long int PackOwnedInteriorFields_(
    const std::vector<cAMRNode*>& nodes,
    long int magneticFieldDataOffset,
    long int electricFieldDataOffset,
    long int plasmaVelocityDataOffset,
    bool deriveElectricFromVelocity,
    std::vector<double>& magneticField,
    std::vector<double>& electricField,
    std::vector<double>& plasmaVelocity,
    std::vector<int>& presence) {

  long int nPacked=0;

  for (std::vector<cAMRNode*>::const_iterator it=nodes.begin();it!=nodes.end();++it) {
    cAMRNode* node=*it;

    if (node->Thread!=PIC::ThisThread) continue;
    if (node->block==NULL) continue;

    for (int i=0;i<_BLOCK_CELLS_X_;i++) {
      for (int j=0;j<_BLOCK_CELLS_Y_;j++) {
        for (int k=0;k<_BLOCK_CELLS_Z_;k++) {
          PIC::Mesh::cDataCenterNode* cell=
              node->block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k));

          if (cell==NULL) continue;

          const long int cellIndex=GlobalCellIndex_(node,i,j,k);
          const long int vectorIndex=3*cellIndex;
          char* data=cell->GetAssociatedDataBufferPointer();

          std::memcpy(&magneticField[static_cast<size_t>(vectorIndex)],
                      data+magneticFieldDataOffset,3*sizeof(double));

          // PIC::CPLR stores the imported SWMF quantities in AMPS SI units (tesla and
          // m/s). Retain u even for the released magnetic-only mode: it is required for
          // provenance, replay, and the independent E=-u×B diagnostic, but it does not
          // make E available to a trajectory unless the experimental mode is explicit.
          double velocity[3]={0.0,0.0,0.0};
          if (plasmaVelocityDataOffset>=0) {
            std::memcpy(velocity,data+plasmaVelocityDataOffset,3*sizeof(double));
            std::memcpy(&plasmaVelocity[static_cast<size_t>(vectorIndex)],
                        velocity,3*sizeof(double));
          }

          if (electricFieldDataOffset>=0) {
            // Standalone DATAFILE path: E was initialized explicitly in the same
            // cell-associated buffer as B.
            std::memcpy(&electricField[static_cast<size_t>(vectorIndex)],
                        data+electricFieldDataOffset,3*sizeof(double));
          }
          else if (deriveElectricFromVelocity) {
            // Experimental SWMF path: derive E only after the caller explicitly opts
            // in. The released Phase-1 default leaves the valid compact E array zero
            // and advertises electricFieldAvailable=false in its snapshot metadata.
            double b[3],electric[3];
            std::memcpy(b,data+magneticFieldDataOffset,3*sizeof(double));
            Earth::SWMFSnapshot::DeriveElectricField(velocity,b,electric);
            std::memcpy(&electricField[static_cast<size_t>(vectorIndex)],
                        electric,3*sizeof(double));
          }
          else {
            // No electric-field source was requested.  The vector was initialized to
            // zero, so no write is required.  Keeping an explicitly valid zero array
            // simplifies future movers that request both fields.
          }

          // Reject bad source data before MPI reduction.  Allowing NaN/Inf into the
          // compact arrays would contaminate interpolation on every rank and could be
          // mistaken downstream for a physical forbidden trajectory.
          const double* packedB=&magneticField[static_cast<size_t>(vectorIndex)];
          const double* packedE=&electricField[static_cast<size_t>(vectorIndex)];
          const double* packedU=(plasmaVelocityDataOffset>=0) ?
              &plasmaVelocity[static_cast<size_t>(vectorIndex)] : NULL;
          if (!Earth::Field::FiniteVector3(packedB) ||
              !Earth::Field::FiniteVector3(packedE) ||
              (plasmaVelocityDataOffset>=0 &&
               !Earth::Field::FiniteVector3(packedU))) {
            std::ostringstream msg;
            msg << "[Mode3D::GlobalMagneticField] non-finite owner-cell field at "
                << "Temp_ID=" << node->Temp_ID
                << ", cell=(" << i << ',' << j << ',' << k << ").";
            const std::string text=msg.str();
            exit(__LINE__,__FILE__,text.c_str());
          }

          presence[static_cast<size_t>(cellIndex)]=1;
          nPacked++;
        }
      }
    }
  }

  return nPacked;
}

// Re-read every authoritative owner cell after the collective reduction and compare it
// with the compact value that trajectories will see. This is a deliberately exact gate
// for B and u (ordinary equality treats +/-0 as the same physical value). It catches a
// wrong offset, an owner/ghost mix-up, or a coupler buffer that changed while the
// snapshot was being assembled. Derived E is recomputed with the shared contract so a
// sign/convention drift cannot pass merely because both arrays are finite.
void ValidateOwnerCellParity_(
    const std::vector<cAMRNode*>& nodes,
    long int magneticFieldDataOffset,
    long int electricFieldDataOffset,
    long int plasmaVelocityDataOffset,
    bool deriveElectricFromVelocity,
    const std::string& tag) {
  long int localMismatches=0;
  for (std::vector<cAMRNode*>::const_iterator it=nodes.begin();it!=nodes.end();++it) {
    cAMRNode* node=*it;
    if (node->Thread!=PIC::ThisThread || node->block==NULL) continue;
    for (int i=0;i<_BLOCK_CELLS_X_;++i) {
      for (int j=0;j<_BLOCK_CELLS_Y_;++j) {
        for (int k=0;k<_BLOCK_CELLS_Z_;++k) {
          PIC::Mesh::cDataCenterNode* cell=
              node->block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k));
          if (cell==NULL) continue;
          const long int index=GlobalCellIndex_(node,i,j,k);
          const std::size_t vectorIndex=static_cast<std::size_t>(3*index);
          char* source=cell->GetAssociatedDataBufferPointer();
          double b[3];
          std::memcpy(b,source+magneticFieldDataOffset,3*sizeof(double));
          for (int d=0;d<3;++d)
            if (b[d]!=GlobalMagneticField_[vectorIndex+d]) ++localMismatches;

          double velocity[3]={0.0,0.0,0.0};
          if (plasmaVelocityDataOffset>=0) {
            std::memcpy(velocity,source+plasmaVelocityDataOffset,3*sizeof(double));
            for (int d=0;d<3;++d)
              if (velocity[d]!=GlobalPlasmaVelocity_[vectorIndex+d])
                ++localMismatches;
          }

          double expectedE[3]={0.0,0.0,0.0};
          if (electricFieldDataOffset>=0)
            std::memcpy(expectedE,source+electricFieldDataOffset,3*sizeof(double));
          else if (deriveElectricFromVelocity)
            Earth::SWMFSnapshot::DeriveElectricField(velocity,b,expectedE);
          for (int d=0;d<3;++d)
            if (expectedE[d]!=GlobalElectricField_[vectorIndex+d])
              ++localMismatches;
        }
      }
    }
  }

  long int globalMismatches=0;
  MPI_Allreduce(&localMismatches,&globalMismatches,1,MPI_LONG,MPI_SUM,
                MPI_GLOBAL_COMMUNICATOR);
  if (globalMismatches!=0) {
    std::ostringstream message;
    message << "[" << tag << "] compact SWMF owner-cell parity failed for "
            << globalMismatches
            << " component(s); field offsets, ownership, or receive stability are invalid.";
    exit(__LINE__,__FILE__,message.str().c_str());
  }
}

void ValidateMeshAndOffset_(const std::string& tag,long int magneticFieldDataOffset) {
  if (PIC::Mesh::mesh==NULL) {
    const std::string msg="["+tag+"] compact global field assembly called before PIC::Mesh::mesh is initialized.";
    exit(__LINE__,__FILE__,msg.c_str());
  }

  if (PIC::Mesh::mesh->rootTree==NULL) {
    const std::string msg="["+tag+"] compact global field assembly called before the AMR root tree exists.";
    exit(__LINE__,__FILE__,msg.c_str());
  }

  if (magneticFieldDataOffset<0) {
    const std::string msg="["+tag+"] compact global field assembly received a negative magnetic-field data offset.";
    exit(__LINE__,__FILE__,msg.c_str());
  }
}

// Validate one row-stencil element and return its compact global cell index.  A bad
// Temp_ID after assembly means some other code reused Temp_ID while the snapshot was
// active; continuing would read an unrelated field cell and silently corrupt particle
// trajectories, so this is intentionally fatal.
long int CheckedGlobalCellIndex_(cAMRNode* node,int i,int j,int k,const char* fieldName) {
  if (node==NULL) {
    std::string msg="[Mode3D::GlobalMagneticField] NULL AMR node in ";
    msg+=fieldName;
    msg+=" row stencil.";
    exit(__LINE__,__FILE__,msg.c_str());
  }

  if ((i<0)||(i>=_BLOCK_CELLS_X_) ||
      (j<0)||(j>=_BLOCK_CELLS_Y_) ||
      (k<0)||(k>=_BLOCK_CELLS_Z_)) {
    std::ostringstream msg;
    msg << "[Mode3D::GlobalMagneticField] non-interior row-stencil index for "
        << fieldName << ": (" << i << "," << j << "," << k << ").";
    const std::string text=msg.str();
    exit(__LINE__,__FILE__,text.c_str());
  }

  if ((node->Temp_ID<0)||(node->Temp_ID>=GlobalUsedLeafBlocks_)) {
    std::ostringstream msg;
    msg << "[Mode3D::GlobalMagneticField] invalid node->Temp_ID=" << node->Temp_ID
        << " while evaluating " << fieldName
        << ". Temp_ID must not be reused after compact global fields are assembled.";
    const std::string text=msg.str();
    exit(__LINE__,__FILE__,text.c_str());
  }

  const long int cellIndex=GlobalCellIndex_(node,i,j,k);

  if ((cellIndex<0)||(cellIndex>=GlobalInteriorCellCount_) ||
      (GlobalCellPresence_[static_cast<size_t>(cellIndex)]!=1)) {
    std::ostringstream msg;
    msg << "[Mode3D::GlobalMagneticField] missing compact global " << fieldName
        << " value for Temp_ID=" << node->Temp_ID
        << ", cell=(" << i << "," << j << "," << k << ").";
    const std::string text=msg.str();
    exit(__LINE__,__FILE__,text.c_str());
  }

  return cellIndex;
}

bool GetCellCenteredField_(cAMRNode* node,int i,int j,int k,double* field,
                           const std::vector<double>& storage,const char* fieldName) {
  if ((field==NULL)||(GlobalFieldsReady_==false)) return false;

  const long int cellIndex=CheckedGlobalCellIndex_(node,i,j,k,fieldName);
  const long int vectorIndex=3*cellIndex;

  field[0]=storage[static_cast<size_t>(vectorIndex+0)];
  field[1]=storage[static_cast<size_t>(vectorIndex+1)];
  field[2]=storage[static_cast<size_t>(vectorIndex+2)];
  return true;
}

bool InterpolateField_(const double* x,cAMRNode* node,double* field,
                       const std::vector<double>& storage,const char* fieldName) {
  if ((x==NULL)||(field==NULL)||(GlobalFieldsReady_==false)) return false;

  double xLocal[3]={x[0],x[1],x[2]};
  cAMRNode* interpolationNode=node;

  if (interpolationNode==NULL) {
    interpolationNode=PIC::Mesh::mesh->findTreeNode(xLocal);
  }

  if (interpolationNode==NULL) return false;

  // The row-only overload intentionally does not require interpolationNode->block.
  // AMPS constructs the same geometric linear/multiblock/blended stencil as the
  // legacy pointer path, but records owning nodes and interior cell indices.
  PIC::InterpolationRoutines::cRowStencil row;
  PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(
      xLocal,interpolationNode,row);

  if (row.Length<=0) return false;

  field[0]=field[1]=field[2]=0.0;

  for (int s=0;s<row.Length;s++) {
    const PIC::InterpolationRoutines::cRowStencil::cElement& e=row.Element[s];
    const long int cellIndex=CheckedGlobalCellIndex_(e.node,e.i,e.j,e.k,fieldName);
    const long int vectorIndex=3*cellIndex;

    field[0]+=e.Weight*storage[static_cast<size_t>(vectorIndex+0)];
    field[1]+=e.Weight*storage[static_cast<size_t>(vectorIndex+1)];
    field[2]+=e.Weight*storage[static_cast<size_t>(vectorIndex+2)];
  }

  return true;
}

Earth::SWMFSnapshot::Snapshot BuildPortableSWMFSnapshot_(
    const Earth::Field::SnapshotMetadata& metadata,
    const std::vector<cAMRNode*>& nodes) {
  // Build records only from the already-reduced compact arrays.  Never revisit the
  // mutable SWMF cell buffers here: the exported file must be exactly the generation
  // used by trajectories, even if the coupler receives another state later.
  Earth::SWMFSnapshot::Snapshot snapshot;
  snapshot.epochUTC=metadata.epochUTC;
  snapshot.simulationTime_s=GlobalSourceSimulationTime_s_;
  snapshot.electricFieldMode=metadata.electricFieldAvailable ?
      Earth::SWMFSnapshot::kExperimentalIdealMhdMode :
      Earth::SWMFSnapshot::kMagneticOnlyMode;
  snapshot.usedLeafBlocks=GlobalUsedLeafBlocks_;
  snapshot.blockCellsX=_BLOCK_CELLS_X_;
  snapshot.blockCellsY=_BLOCK_CELLS_Y_;
  snapshot.blockCellsZ=_BLOCK_CELLS_Z_;
  snapshot.domain=metadata.domain;
  snapshot.cells.reserve(static_cast<std::size_t>(GlobalInteriorCellCount_));

  const std::size_t expectedVectorValues=
      static_cast<std::size_t>(3*GlobalInteriorCellCount_);
  if (GlobalMagneticField_.size()!=expectedVectorValues ||
      GlobalPlasmaVelocity_.size()!=expectedVectorValues ||
      GlobalCellPresence_.size()!=
          static_cast<std::size_t>(GlobalInteriorCellCount_))
    throw std::runtime_error(
        "SWMF portable snapshot requested from incomplete compact arrays");

  for (std::vector<cAMRNode*>::const_iterator it=nodes.begin();it!=nodes.end();++it) {
    cAMRNode* node=*it;
    // Temp_ID is shared AMPS scratch state. Validate it before using it as a compact
    // vector index so an intervening mesh algorithm becomes a clean stale-state
    // failure rather than an out-of-bounds read during export.
    if (node==NULL || node->Temp_ID<0 || node->Temp_ID>=GlobalUsedLeafBlocks_)
      throw std::runtime_error(
          "SWMF portable snapshot found an invalid AMR block/Temp_ID mapping");
    const double dx[3]={
      (node->xmax[0]-node->xmin[0])/_BLOCK_CELLS_X_,
      (node->xmax[1]-node->xmin[1])/_BLOCK_CELLS_Y_,
      (node->xmax[2]-node->xmin[2])/_BLOCK_CELLS_Z_};
    for (int i=0;i<_BLOCK_CELLS_X_;++i) {
      for (int j=0;j<_BLOCK_CELLS_Y_;++j) {
        for (int k=0;k<_BLOCK_CELLS_Z_;++k) {
          const long int cellIndex=GlobalCellIndex_(node,i,j,k);
          if (cellIndex<0 || cellIndex>=GlobalInteriorCellCount_ ||
              GlobalCellPresence_[static_cast<std::size_t>(cellIndex)]!=1)
            throw std::runtime_error(
                "Cannot export SWMF snapshot: compact cell coverage is incomplete");

          Earth::SWMFSnapshot::Cell cell;
          cell.blockId=node->Temp_ID;
          cell.i=i; cell.j=j; cell.k=k;
          cell.position_m[0]=node->xmin[0]+(i+0.5)*dx[0];
          cell.position_m[1]=node->xmin[1]+(j+0.5)*dx[1];
          cell.position_m[2]=node->xmin[2]+(k+0.5)*dx[2];
          for (int d=0;d<3;++d) {
            cell.magneticField_T[d]=GlobalMagneticField_[
                static_cast<std::size_t>(3*cellIndex+d)];
            cell.plasmaVelocity_m_s[d]=GlobalPlasmaVelocity_[
                static_cast<std::size_t>(3*cellIndex+d)];
          }
          snapshot.cells.push_back(cell);
        }
      }
    }
  }
  Earth::SWMFSnapshot::FinalizeIdentity(snapshot);
  Earth::SWMFSnapshot::Validate(snapshot,true);
  return snapshot;
}

Earth::Field::SnapshotMetadata LegacySnapshotMetadata_(
    const std::string& tag,
    bool electricFieldAvailable,
    bool derivedElectricField) {
  Earth::Field::SnapshotMetadata metadata;
  metadata.sourceId=tag+":LEGACY_CELL_BUFFER";
  metadata.modelName="LEGACY_CELL_BUFFER";
  metadata.epochUTC="UNSPECIFIED";
  metadata.frame=Earth::Field::CoordinateFrame::GSM;
  metadata.interpolation=derivedElectricField ?
      Earth::Field::InterpolationMode::CellCenteredLinearDerivedElectric :
      Earth::Field::InterpolationMode::CellCenteredLinear;
  metadata.magneticFieldAvailable=true;
  metadata.electricFieldAvailable=electricFieldAvailable;
  metadata.immutableDuringBatch=true;
  metadata.valid=true;
  metadata.validityMessage=
      "legacy caller did not supply physical snapshot provenance";
  std::ostringstream state;
  state << "legacy|electric=" << (electricFieldAvailable ? 1 : 0)
        << "|derived=" << (derivedElectricField ? 1 : 0);
  metadata.snapshotId=Earth::Field::MakeSnapshotId(
      metadata.sourceId,metadata.epochUTC,state.str());
  return metadata;
}

// Generic view over one published compact-array generation.  It does not own/copy the
// potentially large arrays; generation checks prevent it from ever sampling a later
// replacement under an earlier snapshot identity.
class CompactFieldSnapshotView_ : public Earth::Field::IFieldSnapshot {
public:
  CompactFieldSnapshotView_(const Earth::Field::SnapshotMetadata& metadata,
                            unsigned long long generation)
      : metadata_(metadata),generation_(generation) {}

  const Earth::Field::SnapshotMetadata& Metadata() const override {
    return metadata_;
  }

  Earth::Field::FieldSample Sample(
      const Earth::Field::FieldQuery& query) const override {
    Earth::Field::FieldSample sample;
    sample.snapshotId=metadata_.snapshotId;
    sample.interpolation=metadata_.interpolation;

    sample.status=Earth::Field::ValidateMetadata(metadata_,&sample.message);
    if (!sample.ok()) return sample;

    if (!GlobalFieldsReady_ || generation_!=GlobalSnapshotGeneration_ ||
        metadata_.snapshotId!=GlobalSnapshotMetadata_.snapshotId) {
      sample.status=Earth::Field::FieldSampleStatus::StaleEpoch;
      sample.message="compact field snapshot was superseded by another assembly";
      return sample;
    }

    sample.status=Earth::Field::ValidateQuery(metadata_,query,&sample.message);
    if (!sample.ok()) return sample;

    if (PIC::Mesh::mesh==NULL) {
      sample.status=Earth::Field::FieldSampleStatus::SourceUnavailable;
      sample.message="AMPS mesh is unavailable for compact field interpolation";
      return sample;
    }

    double position[3]={query.position_m[0],query.position_m[1],query.position_m[2]};
    cAMRNode* node=PIC::Mesh::mesh->findTreeNode(position);
    if (node==NULL) {
      sample.status=Earth::Field::FieldSampleStatus::OutsideDomain;
      sample.message="field query position is outside the used AMR tree";
      return sample;
    }

    if (!InterpolateField_(position,node,sample.magneticField_T,
                           GlobalMagneticField_,"magnetic field")) {
      sample.status=Earth::Field::FieldSampleStatus::InterpolationFailure;
      sample.message="failed to construct the magnetic-field interpolation row";
      return sample;
    }
    if (metadata_.electricFieldAvailable &&
        !InterpolateField_(position,node,sample.electricField_V_m,
                           GlobalElectricField_,"electric field")) {
      sample.status=Earth::Field::FieldSampleStatus::InterpolationFailure;
      sample.message="failed to construct the electric-field interpolation row";
      return sample;
    }
    if (!Earth::Field::FiniteVector3(sample.magneticField_T) ||
        (metadata_.electricFieldAvailable &&
         !Earth::Field::FiniteVector3(sample.electricField_V_m))) {
      sample.status=Earth::Field::FieldSampleStatus::NonFiniteValue;
      sample.message="compact field interpolation returned a non-finite value";
      return sample;
    }

    sample.status=Earth::Field::FieldSampleStatus::Valid;
    sample.message.clear();
    return sample;
  }

private:
  const Earth::Field::SnapshotMetadata metadata_;
  const unsigned long long generation_;
};

class CompactFieldProviderView_ : public Earth::Field::IFieldProvider {
public:
  std::string SourceId() const override {
    return GlobalFieldsReady_ ? GlobalSnapshotMetadata_.sourceId :
                                std::string("MODE3D:UNAVAILABLE");
  }

  std::shared_ptr<const Earth::Field::IFieldSnapshot> CreateSnapshot(
      const Earth::Field::SnapshotRequest& request) override {
    if (!GlobalFieldsReady_)
      throw std::runtime_error(
          "Mode3D compact field provider requested before snapshot assembly");
    if (!request.epochUTC.empty() &&
        request.epochUTC!=GlobalSnapshotMetadata_.epochUTC)
      throw std::runtime_error(
          "Mode3D field-provider request epoch differs from published snapshot");

    // requestId is diagnostic only and cannot affect physical snapshot identity.
    return std::shared_ptr<const Earth::Field::IFieldSnapshot>(
        new CompactFieldSnapshotView_(
            GlobalSnapshotMetadata_,GlobalSnapshotGeneration_));
  }
};

} // anonymous namespace

long int DataFileMagneticFieldDataOffset() {
  return PIC::CPLR::DATAFILE::CenterNodeAssociatedDataOffsetBegin +
         PIC::CPLR::DATAFILE::MULTIFILE::CurrDataFileOffset +
         PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
}

long int DataFileElectricFieldDataOffset() {
  return PIC::CPLR::DATAFILE::CenterNodeAssociatedDataOffsetBegin +
         PIC::CPLR::DATAFILE::MULTIFILE::CurrDataFileOffset +
         PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
}

MaterializationStats AssembleCellCenteredFieldsForCutoff(
    const char* diagnosticTag,
    long int magneticFieldDataOffset,
    long int electricFieldDataOffset,
    long int plasmaVelocityDataOffset,
    bool verbose) {

  const std::string tag=SafeTag_(diagnosticTag);
  const bool electricAvailable=
      (electricFieldDataOffset>=0 || plasmaVelocityDataOffset>=0);
  const Earth::Field::SnapshotMetadata metadata=LegacySnapshotMetadata_(
      tag,electricAvailable,
      electricFieldDataOffset<0 && plasmaVelocityDataOffset>=0);
  return AssembleCellCenteredFieldsForCutoff(
      diagnosticTag,magneticFieldDataOffset,electricFieldDataOffset,
      plasmaVelocityDataOffset,metadata,verbose);
}

MaterializationStats AssembleCellCenteredFieldsForCutoff(
    const char* diagnosticTag,
    long int magneticFieldDataOffset,
    long int electricFieldDataOffset,
    long int plasmaVelocityDataOffset,
    const Earth::Field::SnapshotMetadata& metadata,
    bool verbose,
    double sourceSimulationTime_s) {

  const std::string tag=SafeTag_(diagnosticTag);
  RequireNoFrozenBatch_("assemble a replacement field snapshot");

  // Invalidate first, before validating the incoming offsets/metadata. A malformed
  // replacement is not permission to keep serving the previous epoch. Arrays are left
  // allocated for efficient refill, but every public sampler now reports unavailable
  // until all ownership/content gates below publish a new generation.
  GlobalFieldsReady_=false;
  GlobalPlasmaVelocityAvailable_=false;
  GlobalContentFingerprint_.clear();
  GlobalMeshRevision_.clear();
  GlobalSnapshotMetadata_=Earth::Field::SnapshotMetadata();
  ValidateMeshAndOffset_(tag,magneticFieldDataOffset);
  if (!std::isfinite(sourceSimulationTime_s) || sourceSimulationTime_s<0.0) {
    const std::string msg="["+tag+
        "] source simulation time must be finite and nonnegative.";
    exit(__LINE__,__FILE__,msg.c_str());
  }

  std::string metadataError;
  const Earth::Field::FieldSampleStatus metadataStatus=
      Earth::Field::ValidateMetadata(metadata,&metadataError);
  if (metadataStatus!=Earth::Field::FieldSampleStatus::Valid) {
    std::ostringstream msg;
    msg << "[" << tag << "] invalid field-snapshot metadata: status="
        << Earth::Field::FieldSampleStatusName(metadataStatus)
        << ", detail=" << metadataError << ".";
    const std::string text=msg.str();
    exit(__LINE__,__FILE__,text.c_str());
  }

  const bool deriveElectricFromVelocity=
      metadata.electricFieldAvailable && electricFieldDataOffset<0 &&
      plasmaVelocityDataOffset>=0;
  const bool electricSourceAvailable=
      (electricFieldDataOffset>=0 || deriveElectricFromVelocity);
  if (metadata.electricFieldAvailable && !electricSourceAvailable) {
    const std::string msg="["+tag+
        "] metadata advertises E, but no E or plasma-velocity source was supplied.";
    exit(__LINE__,__FILE__,msg.c_str());
  }
  const Earth::Field::InterpolationMode expectedInterpolation=
      deriveElectricFromVelocity ?
      Earth::Field::InterpolationMode::CellCenteredLinearDerivedElectric :
      Earth::Field::InterpolationMode::CellCenteredLinear;
  if (metadata.interpolation!=expectedInterpolation) {
    const std::string msg="["+tag+
        "] metadata interpolation mode does not match compact-array assembly.";
    exit(__LINE__,__FILE__,msg.c_str());
  }

  MaterializationStats stats;

  // The old generation was invalidated before incoming validation. Assembly is
  // expected to run before worker threads; the new state remains unavailable until
  // the final atomic publication below.
  GlobalSourceSimulationTime_s_=sourceSimulationTime_s;

  // The explicit reset is required because Temp_ID is shared scratch storage in AMPS.
  // IDs are then reassigned deterministically over the complete global tree.
  ResetTreeTempIds_(PIC::Mesh::mesh->rootTree);
  long int nUsedLeafBlocks=0;
  AssignGlobalLeafTempIds_(PIC::Mesh::mesh->rootTree,nUsedLeafBlocks);

  std::vector<cAMRNode*> nodes;
  nodes.reserve(static_cast<size_t>(std::max(static_cast<long int>(0),nUsedLeafBlocks)));
  CollectUsedLeafNodes_(PIC::Mesh::mesh->rootTree,nodes);

  if (static_cast<long int>(nodes.size())!=nUsedLeafBlocks) {
    std::ostringstream msg;
    msg << "[" << tag << "] inconsistent used-leaf traversal: assigned "
        << nUsedLeafBlocks << " Temp_ID values but collected " << nodes.size()
        << " leaves.";
    const std::string text=msg.str();
    exit(__LINE__,__FILE__,text.c_str());
  }

  const long int nInteriorCells=nUsedLeafBlocks*InteriorCellsPerBlock_();

  GlobalUsedLeafBlocks_=nUsedLeafBlocks;
  GlobalInteriorCellCount_=nInteriorCells;

  // The AMR topology is invariant during a standalone multi-snapshot run, so these
  // compact arrays have the same dimensions at every epoch.  Resize only when the
  // required size actually changes and otherwise clear the existing storage in
  // place.  This makes snapshot batching reuse both the distributed AMR blocks and
  // the compact replicated field buffers instead of asking the allocator to rebuild
  // three large vectors for every field update.
  //
  // Temp_ID is still reset/reassigned above because it is shared AMPS scratch state;
  // caching it across products would be unsafe.  Buffer reuse is independent of that
  // correctness guard and preserves the deterministic leaf-to-cell mapping.
  const std::size_t nVectorValues=static_cast<std::size_t>(3*nInteriorCells);
  const std::size_t nPresenceValues=static_cast<std::size_t>(nInteriorCells);
  if (GlobalMagneticField_.size()!=nVectorValues)
    GlobalMagneticField_.resize(nVectorValues);
  if (GlobalElectricField_.size()!=nVectorValues)
    GlobalElectricField_.resize(nVectorValues);
  if (plasmaVelocityDataOffset>=0) {
    if (GlobalPlasmaVelocity_.size()!=nVectorValues)
      GlobalPlasmaVelocity_.resize(nVectorValues);
  }
  else GlobalPlasmaVelocity_.clear();
  if (GlobalCellPresence_.size()!=nPresenceValues)
    GlobalCellPresence_.resize(nPresenceValues);
  std::fill(GlobalMagneticField_.begin(),GlobalMagneticField_.end(),0.0);
  std::fill(GlobalElectricField_.begin(),GlobalElectricField_.end(),0.0);
  if (plasmaVelocityDataOffset>=0)
    std::fill(GlobalPlasmaVelocity_.begin(),GlobalPlasmaVelocity_.end(),0.0);
  std::fill(GlobalCellPresence_.begin(),GlobalCellPresence_.end(),0);

  const long int nPackedLocal=PackOwnedInteriorFields_(
      nodes,magneticFieldDataOffset,electricFieldDataOffset,
      plasmaVelocityDataOffset,deriveElectricFromVelocity,
      GlobalMagneticField_,GlobalElectricField_,
      GlobalPlasmaVelocity_,GlobalCellPresence_);

  // Replicate only the compact physical fields and one integer validity entry per
  // interior cell.  No cDataBlockAMR, center-node state vector, or ghost cell is
  // allocated by this operation.
  AllreduceDoubleVectorInPlace_(GlobalMagneticField_);
  AllreduceDoubleVectorInPlace_(GlobalElectricField_);
  if (plasmaVelocityDataOffset>=0)
    AllreduceDoubleVectorInPlace_(GlobalPlasmaVelocity_);
  AllreduceIntVectorInPlace_(GlobalCellPresence_);

  long int nPackedGlobal=0;
  MPI_Allreduce((void*)&nPackedLocal,&nPackedGlobal,1,MPI_LONG,MPI_SUM,
                MPI_GLOBAL_COMMUNICATOR);

  long int nMissing=0,nDuplicate=0;
  for (long int c=0;c<nInteriorCells;c++) {
    const int count=GlobalCellPresence_[static_cast<size_t>(c)];

    if (count==0) nMissing++;
    else if (count>1) nDuplicate++;

    // Average defensively before reporting duplicate ownership.  A duplicate is still
    // fatal below, but this keeps the arrays finite for debugger inspection.
    if (count>1) {
      const double inv=1.0/static_cast<double>(count);
      for (int d=0;d<3;d++) {
        GlobalMagneticField_[static_cast<size_t>(3*c+d)]*=inv;
        GlobalElectricField_[static_cast<size_t>(3*c+d)]*=inv;
        if (plasmaVelocityDataOffset>=0)
          GlobalPlasmaVelocity_[static_cast<size_t>(3*c+d)]*=inv;
      }
    }

    const double* reducedB=&GlobalMagneticField_[static_cast<size_t>(3*c)];
    const double* reducedE=&GlobalElectricField_[static_cast<size_t>(3*c)];
    const double* reducedU=(plasmaVelocityDataOffset>=0) ?
        &GlobalPlasmaVelocity_[static_cast<size_t>(3*c)] : NULL;
    if (!Earth::Field::FiniteVector3(reducedB) ||
        !Earth::Field::FiniteVector3(reducedE) ||
        (plasmaVelocityDataOffset>=0 && !Earth::Field::FiniteVector3(reducedU))) {
      std::ostringstream msg;
      msg << "[" << tag << "] non-finite reduced field at compact cell " << c << ".";
      const std::string text=msg.str();
      exit(__LINE__,__FILE__,text.c_str());
    }
  }

  stats.usedLeafBlocks=nUsedLeafBlocks;
  stats.ownerInteriorCells=nPackedGlobal;
  stats.expectedInteriorCells=nInteriorCells;
  stats.missingInteriorCells=nMissing;
  stats.duplicateInteriorCells=nDuplicate;
  stats.magneticFieldBytes=static_cast<long int>(GlobalMagneticField_.size()*sizeof(double));
  stats.electricFieldBytes=static_cast<long int>(GlobalElectricField_.size()*sizeof(double));
  stats.plasmaVelocityBytes=(plasmaVelocityDataOffset>=0) ?
      static_cast<long int>(GlobalPlasmaVelocity_.size()*sizeof(double)) : 0;
  stats.electricFieldReadFromBuffer=(electricFieldDataOffset>=0);
  stats.electricFieldDerivedFromVelocity=deriveElectricFromVelocity;
  stats.plasmaVelocityAvailable=(plasmaVelocityDataOffset>=0);

  if ((nPackedGlobal!=nInteriorCells)||(nMissing!=0)||(nDuplicate!=0)) {
    std::ostringstream msg;
    msg << "[" << tag << "] compact global field gather is inconsistent: ownerCells="
        << nPackedGlobal << ", expected=" << nInteriorCells
        << ", missing=" << nMissing << ", duplicate=" << nDuplicate
        << ". Every used interior cell must be supplied by exactly one AMPS owner rank.";
    const std::string text=msg.str();
    exit(__LINE__,__FILE__,text.c_str());
  }

  ValidateOwnerCellParity_(nodes,magneticFieldDataOffset,electricFieldDataOffset,
                           plasmaVelocityDataOffset,deriveElectricFromVelocity,tag);
  stats.ownerCellParityValidated=true;

  // A live SWMF generation receives a content-derived identity from the exact reduced
  // B/u arrays.  This makes the identity independent of MPI/OpenMP decomposition and
  // of the callback counter, while changing it for any physical cell value, epoch,
  // topology, or domain change.  Standalone empirical snapshots have no u array and
  // retain their established driver-derived Step-3 identity.
  Earth::Field::SnapshotMetadata publishedMetadata=metadata;
  GlobalPlasmaVelocityAvailable_=(plasmaVelocityDataOffset>=0);
  GlobalContentFingerprint_.clear();
  GlobalMeshRevision_.clear();
  if (GlobalPlasmaVelocityAvailable_) {
    const Earth::SWMFSnapshot::Snapshot portable=
        BuildPortableSWMFSnapshot_(metadata,nodes);
    publishedMetadata.snapshotId=portable.snapshotId;
    GlobalContentFingerprint_=portable.contentFingerprint;
    GlobalMeshRevision_=portable.meshRevision;

    const std::string::size_type dash=portable.contentFingerprint.rfind('-');
    const std::string hex=(dash==std::string::npos) ? portable.contentFingerprint :
                                                    portable.contentFingerprint.substr(dash+1);
    const unsigned long long localFingerprint=
        static_cast<unsigned long long>(std::strtoull(hex.c_str(),NULL,16));
    unsigned long long minFingerprint=0,maxFingerprint=0;
    MPI_Allreduce((void*)&localFingerprint,&minFingerprint,1,
                  MPI_UNSIGNED_LONG_LONG,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce((void*)&localFingerprint,&maxFingerprint,1,
                  MPI_UNSIGNED_LONG_LONG,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
    if (minFingerprint!=maxFingerprint) {
      const std::string msg="["+tag+
          "] SWMF snapshot epoch/domain/content identity differs across MPI ranks.";
      exit(__LINE__,__FILE__,msg.c_str());
    }
  }

  stats.snapshotId=publishedMetadata.snapshotId;
  stats.contentFingerprint=GlobalContentFingerprint_;
  stats.meshRevision=GlobalMeshRevision_;

  // Publish metadata and generation only after ownership, finite-value, topology, and
  // cross-rank identity gates pass.  A consumer can never observe new identity with
  // old arrays or an identity assembled from different states on different ranks.
  GlobalSnapshotMetadata_=publishedMetadata;
  ++GlobalSnapshotGeneration_;
  GlobalFieldsReady_=true;

  if ((verbose==true)&&(PIC::ThisThread==0)) {
    const double mib=1024.0*1024.0;
    std::cout << "[" << tag << "] Prepared compact global cell-centered fields:"
              << " usedLeafBlocks=" << stats.usedLeafBlocks
              << ", interiorCells=" << stats.expectedInteriorCells
              << ", B=" << stats.magneticFieldBytes/mib << " MiB"
              << ", E=" << stats.electricFieldBytes/mib << " MiB"
              << ", u=" << stats.plasmaVelocityBytes/mib << " MiB"
              << ", snapshot=" << stats.snapshotId
              << ", content=" << (stats.contentFingerprint.empty() ?
                                    std::string("not-applicable") :
                                    stats.contentFingerprint)
              << ", meshRevision=" << (stats.meshRevision.empty() ?
                                         std::string("not-applicable") :
                                         stats.meshRevision)
              << ", epoch=" << metadata.epochUTC
              << ", tSim_s=" << sourceSimulationTime_s
              << ", frame=" << Earth::Field::CoordinateFrameName(metadata.frame)
              << ", ESource=";

    if (stats.electricFieldReadFromBuffer) std::cout << "cell buffer";
    else if (stats.electricFieldDerivedFromVelocity) std::cout << "-VxB";
    else std::cout << "zero";

    std::cout << ", ownerParity=PASS. No nonlocal AMR blocks were allocated.\n";
    std::cout.flush();
  }

  return stats;
}

MaterializationStats MaterializeCellCenteredMagneticFieldForCutoff(
    const char* diagnosticTag,long int magneticFieldDataOffset,bool verbose) {
  return AssembleCellCenteredFieldsForCutoff(
      diagnosticTag,magneticFieldDataOffset,-1,-1,verbose);
}

bool GetCellCenteredMagneticField(cAMRNode* node,int i,int j,int k,double* B) {
  return GetCellCenteredField_(node,i,j,k,B,GlobalMagneticField_,"magnetic field");
}

bool GetCellCenteredElectricField(cAMRNode* node,int i,int j,int k,double* E) {
  return GetCellCenteredField_(node,i,j,k,E,GlobalElectricField_,"electric field");
}

bool GetCellCenteredPlasmaVelocity(cAMRNode* node,int i,int j,int k,double* velocity) {
  if (!GlobalPlasmaVelocityAvailable_) return false;
  return GetCellCenteredField_(node,i,j,k,velocity,GlobalPlasmaVelocity_,
                               "plasma velocity");
}

bool InterpolateMagneticField(const double* x,cAMRNode* node,double* B) {
  return InterpolateField_(x,node,B,GlobalMagneticField_,"magnetic field");
}

bool InterpolateElectricField(const double* x,cAMRNode* node,double* E) {
  return InterpolateField_(x,node,E,GlobalElectricField_,"electric field");
}

bool GlobalFieldsReady() {
  return GlobalFieldsReady_;
}

long int GlobalCellCount() {
  return GlobalInteriorCellCount_;
}

void ClearGlobalFields() {
  RequireNoFrozenBatch_("clear compact field arrays");
  GlobalFieldsReady_=false;
  GlobalUsedLeafBlocks_=0;
  GlobalInteriorCellCount_=0;
  GlobalMagneticField_.clear();
  GlobalElectricField_.clear();
  GlobalPlasmaVelocity_.clear();
  GlobalCellPresence_.clear();
  GlobalPlasmaVelocityAvailable_=false;
  GlobalContentFingerprint_.clear();
  GlobalMeshRevision_.clear();
  GlobalSourceSimulationTime_s_=0.0;
  GlobalSnapshotMetadata_=Earth::Field::SnapshotMetadata();
  ++GlobalSnapshotGeneration_;
}

const Earth::Field::SnapshotMetadata& CurrentSnapshotMetadata() {
  return GlobalSnapshotMetadata_;
}

std::shared_ptr<const Earth::Field::IFieldSnapshot> CurrentSnapshot() {
  if (!GlobalFieldsReady_)
    return std::shared_ptr<const Earth::Field::IFieldSnapshot>();
  return std::shared_ptr<const Earth::Field::IFieldSnapshot>(
      new CompactFieldSnapshotView_(
          GlobalSnapshotMetadata_,GlobalSnapshotGeneration_));
}

std::shared_ptr<Earth::Field::IFieldProvider> CurrentFieldProvider() {
  return std::shared_ptr<Earth::Field::IFieldProvider>(
      new CompactFieldProviderView_());
}

void BeginFrozenFieldBatch(const std::string& expectedSnapshotId) {
  if (!GlobalFieldsReady_)
    throw std::runtime_error("Cannot freeze an unavailable compact field snapshot");
  if (GlobalFrozenBatchActive_)
    throw std::runtime_error("A compact field snapshot batch is already frozen");
  if (expectedSnapshotId.empty() ||
      expectedSnapshotId!=GlobalSnapshotMetadata_.snapshotId)
    throw std::runtime_error(
        "Cannot freeze compact fields under a stale or empty snapshot identity");
  GlobalFrozenBatchActive_=true;
  GlobalFrozenBatchSnapshotId_=expectedSnapshotId;
  GlobalFrozenBatchGeneration_=GlobalSnapshotGeneration_;
}

void EndFrozenFieldBatch(const std::string& expectedSnapshotId) {
  if (!GlobalFrozenBatchActive_)
    throw std::runtime_error("No compact field snapshot batch is frozen");
  if (expectedSnapshotId!=GlobalFrozenBatchSnapshotId_ ||
      expectedSnapshotId!=GlobalSnapshotMetadata_.snapshotId ||
      GlobalFrozenBatchGeneration_!=GlobalSnapshotGeneration_)
    throw std::runtime_error(
        "Compact field generation changed while a product batch was frozen");
  GlobalFrozenBatchActive_=false;
  GlobalFrozenBatchSnapshotId_.clear();
  GlobalFrozenBatchGeneration_=0;
}

bool FrozenFieldBatchActive() {
  return GlobalFrozenBatchActive_;
}

bool PlasmaVelocityAvailable() {
  return GlobalFieldsReady_ && GlobalPlasmaVelocityAvailable_;
}

std::string CurrentContentFingerprint() {
  return GlobalFieldsReady_ ? GlobalContentFingerprint_ : std::string();
}

std::string CurrentMeshRevision() {
  return GlobalFieldsReady_ ? GlobalMeshRevision_ : std::string();
}

std::string ExportCurrentSWMFSnapshot(const std::string& fileName) {
  // Validate locally inside a catch boundary, then reduce the outcome. This is needed
  // even though normal states are replicated: a changed Temp_ID or rank-local memory
  // error must make every rank stop before rank zero enters filesystem I/O.
  int localValidationSucceeded=1;
  std::string localValidationError;
  std::vector<cAMRNode*> nodes;
  Earth::SWMFSnapshot::Snapshot snapshot;
  try {
    if (fileName.empty())
      throw std::invalid_argument(
          "SWMF snapshot export requires a non-empty file name");
    if (!GlobalFieldsReady_)
      throw std::runtime_error(
          "SWMF snapshot export requested before compact field assembly");
    if (!GlobalFrozenBatchActive_ ||
        GlobalFrozenBatchSnapshotId_!=GlobalSnapshotMetadata_.snapshotId)
      throw std::runtime_error(
          "SWMF snapshot export requires the active generation to be frozen first");
    if (!GlobalPlasmaVelocityAvailable_)
      throw std::runtime_error(
          "SWMF snapshot export requires the frozen plasma-velocity array");
    if (PIC::Mesh::mesh==NULL || PIC::Mesh::mesh->rootTree==NULL)
      throw std::runtime_error(
          "SWMF snapshot export requires an initialized AMR tree");

    nodes.reserve(static_cast<std::size_t>(GlobalUsedLeafBlocks_));
    CollectUsedLeafNodes_(PIC::Mesh::mesh->rootTree,nodes);
    if (static_cast<long int>(nodes.size())!=GlobalUsedLeafBlocks_)
      throw std::runtime_error(
          "SWMF snapshot export detected changed AMR Temp_ID/topology state");

    snapshot=BuildPortableSWMFSnapshot_(GlobalSnapshotMetadata_,nodes);
    if (snapshot.snapshotId!=GlobalSnapshotMetadata_.snapshotId ||
        snapshot.contentFingerprint!=GlobalContentFingerprint_ ||
        snapshot.meshRevision!=GlobalMeshRevision_)
      throw std::runtime_error(
          "SWMF snapshot export identity differs from the published trajectory state");
  }
  catch (const std::exception& error) {
    localValidationSucceeded=0;
    localValidationError=error.what();
  }
  catch (...) {
    localValidationSucceeded=0;
    localValidationError="unknown exception";
  }
  int allRanksValidated=0;
  MPI_Allreduce(&localValidationSucceeded,&allRanksValidated,1,MPI_INT,MPI_MIN,
                MPI_GLOBAL_COMMUNICATOR);
  if (!allRanksValidated) {
    std::ostringstream message;
    message << "SWMF snapshot export validation failed on at least one MPI rank";
    if (!localValidationError.empty()) message << ": " << localValidationError;
    throw std::runtime_error(message.str());
  }

  // The compact arrays are replicated and byte-identical on all ranks. Only rank zero
  // writes the portable file, then broadcasts success/failure before any rank can
  // proceed to a product. This avoids both parallel filesystem races and MPI hangs
  // caused by a one-rank I/O failure.
  int writeSucceeded=1;
  std::string writeError;
  if (PIC::ThisThread==0) {
    try { Earth::SWMFSnapshot::Write(snapshot,fileName); }
    catch (const std::exception& error) {
      writeSucceeded=0;
      writeError=error.what();
    }
    catch (...) {
      writeSucceeded=0;
      writeError="unknown exception";
    }
  }
  RequireCollectiveRootWriteSuccess_(
      writeSucceeded,writeError,"SWMF snapshot export");
  return snapshot.contentFingerprint;
}

MaterializationStats ImportSWMFSnapshot(const std::string& fileName,
                                        const std::string& expectedEpochUTC,
                                        bool enableExperimentalDerivedElectric,
                                        bool verbose) {
  RequireNoFrozenBatch_("import a replacement SWMF snapshot");

  // Fail closed before opening the new file. If parsing or validation fails, callers
  // cannot catch the exception and continue tracing against the old generation.
  GlobalFieldsReady_=false;
  GlobalPlasmaVelocityAvailable_=false;
  GlobalContentFingerprint_.clear();
  GlobalMeshRevision_.clear();
  GlobalSnapshotMetadata_=Earth::Field::SnapshotMetadata();
  const int localMeshReady=
      (PIC::Mesh::mesh!=NULL && PIC::Mesh::mesh->rootTree!=NULL) ? 1 : 0;
  int allRanksMeshReady=0;
  MPI_Allreduce((void*)&localMeshReady,&allRanksMeshReady,1,MPI_INT,MPI_MIN,
                MPI_GLOBAL_COMMUNICATOR);
  if (!allRanksMeshReady)
    throw std::runtime_error(
        "SWMF snapshot replay requested before the Mode3D AMR tree is initialized "
        "on every MPI rank");

  // Every rank reads the shared artifact because the compact arrays are replicated.
  // Convert a rank-local filesystem/parser error into a collective failure before any
  // rank touches topology or enters a later collective.
  Earth::SWMFSnapshot::Snapshot snapshot;
  int localReadSucceeded=1;
  std::string localReadError;
  try { snapshot=Earth::SWMFSnapshot::Read(fileName); }
  catch (const std::exception& error) {
    localReadSucceeded=0;
    localReadError=error.what();
  }
  catch (...) {
    localReadSucceeded=0;
    localReadError="unknown exception";
  }
  int allRanksReadSucceeded=0;
  MPI_Allreduce(&localReadSucceeded,&allRanksReadSucceeded,1,MPI_INT,MPI_MIN,
                MPI_GLOBAL_COMMUNICATOR);
  if (!allRanksReadSucceeded) {
    std::ostringstream message;
    message << "SWMF snapshot replay failed on at least one MPI rank";
    if (!localReadError.empty()) message << ": " << localReadError;
    throw std::runtime_error(message.str());
  }

  // Read() has recomputed the fixed-prefix fingerprint. Compare its physical-state
  // word across ranks before allocating replay arrays; a per-rank filesystem view can
  // otherwise produce different but individually valid snapshots.
  const std::string::size_type fingerprintDash=
      snapshot.contentFingerprint.rfind('-');
  const std::string fingerprintHex=(fingerprintDash==std::string::npos) ?
      snapshot.contentFingerprint :
      snapshot.contentFingerprint.substr(fingerprintDash+1);
  const unsigned long long localFingerprint=static_cast<unsigned long long>(
      std::strtoull(fingerprintHex.c_str(),NULL,16));
  unsigned long long minimumFingerprint=0,maximumFingerprint=0;
  MPI_Allreduce((void*)&localFingerprint,&minimumFingerprint,1,
                MPI_UNSIGNED_LONG_LONG,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce((void*)&localFingerprint,&maximumFingerprint,1,
                MPI_UNSIGNED_LONG_LONG,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  if (minimumFingerprint!=maximumFingerprint)
    throw std::runtime_error(
        "SWMF replay content fingerprint differs across MPI ranks");
  long int nUsedLeafBlocks=0;
  int localReplayValidationSucceeded=1;
  std::string localReplayValidationError;
  try {
    if (expectedEpochUTC.empty() || snapshot.epochUTC!=expectedEpochUTC)
      throw std::runtime_error(
          "SWMF snapshot replay epoch differs from #BACKGROUND_FIELD/EPOCH");
    if (snapshot.blockCellsX!=_BLOCK_CELLS_X_ ||
        snapshot.blockCellsY!=_BLOCK_CELLS_Y_ ||
        snapshot.blockCellsZ!=_BLOCK_CELLS_Z_)
      throw std::runtime_error(
          "SWMF snapshot replay block dimensions differ from the compiled AMPS mesh");
    if (Earth::SWMFSnapshot::UsesExperimentalIdealMhdElectricField(snapshot)!=
        enableExperimentalDerivedElectric)
      throw std::runtime_error(
          "SWMF snapshot electric_field_mode differs from "
          "SWMF_DERIVED_ELECTRIC_FIELD in the replay input");

    ResetTreeTempIds_(PIC::Mesh::mesh->rootTree);
    AssignGlobalLeafTempIds_(PIC::Mesh::mesh->rootTree,nUsedLeafBlocks);
    std::vector<cAMRNode*> nodes;
    nodes.reserve(static_cast<std::size_t>(nUsedLeafBlocks));
    CollectUsedLeafNodes_(PIC::Mesh::mesh->rootTree,nodes);
    if (nUsedLeafBlocks!=snapshot.usedLeafBlocks ||
        static_cast<long int>(nodes.size())!=snapshot.usedLeafBlocks)
      throw std::runtime_error(
          "SWMF snapshot replay AMR leaf count differs from the current Mode3D mesh");

    GlobalUsedLeafBlocks_=nUsedLeafBlocks;
    GlobalInteriorCellCount_=nUsedLeafBlocks*InteriorCellsPerBlock_();
    const std::size_t nVectorValues=
        static_cast<std::size_t>(3*GlobalInteriorCellCount_);
    GlobalMagneticField_.assign(nVectorValues,0.0);
    GlobalElectricField_.assign(nVectorValues,0.0);
    GlobalPlasmaVelocity_.assign(nVectorValues,0.0);
    GlobalCellPresence_.assign(static_cast<std::size_t>(GlobalInteriorCellCount_),0);

    std::vector<Earth::SWMFSnapshot::Cell> cells=
        Earth::SWMFSnapshot::CanonicalCells(snapshot);
    for (std::vector<Earth::SWMFSnapshot::Cell>::iterator it=cells.begin();
         it!=cells.end();++it) {
      cAMRNode* node=nodes[static_cast<std::size_t>(it->blockId)];
      if (node->Temp_ID!=it->blockId)
        throw std::runtime_error(
            "SWMF snapshot replay encountered a noncanonical block ID");
      const double dx[3]={
        (node->xmax[0]-node->xmin[0])/_BLOCK_CELLS_X_,
        (node->xmax[1]-node->xmin[1])/_BLOCK_CELLS_Y_,
        (node->xmax[2]-node->xmin[2])/_BLOCK_CELLS_Z_};
      const double expectedPosition[3]={
        node->xmin[0]+(it->i+0.5)*dx[0],
        node->xmin[1]+(it->j+0.5)*dx[1],
        node->xmin[2]+(it->k+0.5)*dx[2]};
      for (int d=0;d<3;++d) {
        const double tolerance=64.0*std::numeric_limits<double>::epsilon()*
            std::max(1.0,std::max(std::fabs(expectedPosition[d]),std::fabs(dx[d])));
        if (std::fabs(it->position_m[d]-expectedPosition[d])>tolerance)
          throw std::runtime_error(
              "SWMF snapshot replay cell centre differs from the current AMR mesh");
        // Replace with the actual current centre after the tolerance diagnostic. The
        // mesh-revision check below then requires exact topology/geometry identity, not
        // merely a matching leaf count or coordinates that happened to be close.
        it->position_m[d]=expectedPosition[d];
      }

      const long int cellIndex=GlobalCellIndex_(node,it->i,it->j,it->k);
      if (cellIndex<0 || cellIndex>=GlobalInteriorCellCount_ ||
          GlobalCellPresence_[static_cast<std::size_t>(cellIndex)]!=0)
        throw std::runtime_error(
            "SWMF snapshot replay has invalid/duplicate cell mapping");
      double electric[3]={0.0,0.0,0.0};
      if (enableExperimentalDerivedElectric)
        Earth::SWMFSnapshot::DeriveElectricField(
            it->plasmaVelocity_m_s,it->magneticField_T,electric);
      for (int d=0;d<3;++d) {
        GlobalMagneticField_[static_cast<std::size_t>(3*cellIndex+d)]=
            it->magneticField_T[d];
        GlobalPlasmaVelocity_[static_cast<std::size_t>(3*cellIndex+d)]=
            it->plasmaVelocity_m_s[d];
        GlobalElectricField_[static_cast<std::size_t>(3*cellIndex+d)]=electric[d];
      }
      GlobalCellPresence_[static_cast<std::size_t>(cellIndex)]=1;
    }

    for (long int cell=0;cell<GlobalInteriorCellCount_;++cell)
      if (GlobalCellPresence_[static_cast<std::size_t>(cell)]!=1)
        throw std::runtime_error(
            "SWMF snapshot replay did not populate every AMR cell");

    Earth::SWMFSnapshot::Snapshot currentMesh;
    currentMesh.usedLeafBlocks=nUsedLeafBlocks;
    currentMesh.blockCellsX=_BLOCK_CELLS_X_;
    currentMesh.blockCellsY=_BLOCK_CELLS_Y_;
    currentMesh.blockCellsZ=_BLOCK_CELLS_Z_;
    currentMesh.domain=snapshot.domain;
    currentMesh.cells.swap(cells);
    if (Earth::SWMFSnapshot::ComputeMeshRevision(currentMesh)!=snapshot.meshRevision)
      throw std::runtime_error(
          "SWMF snapshot mesh_revision differs from the current Mode3D AMR geometry");
  }
  catch (const std::exception& error) {
    localReplayValidationSucceeded=0;
    localReplayValidationError=error.what();
  }
  catch (...) {
    localReplayValidationSucceeded=0;
    localReplayValidationError="unknown exception";
  }
  int allRanksReplayValidated=0;
  MPI_Allreduce(&localReplayValidationSucceeded,&allRanksReplayValidated,1,
                MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  if (!allRanksReplayValidated) {
    std::ostringstream message;
    message << "SWMF snapshot topology/content replay validation failed on at least "
               "one MPI rank";
    if (!localReplayValidationError.empty())
      message << ": " << localReplayValidationError;
    throw std::runtime_error(message.str());
  }

  Earth::Field::SnapshotMetadata metadata;
  metadata.snapshotId=snapshot.snapshotId;
  metadata.sourceId="STANDALONE_MODE3D:SWMF_SNAPSHOT";
  metadata.modelName="SWMF";
  metadata.epochUTC=snapshot.epochUTC;
  metadata.frame=Earth::Field::CoordinateFrame::GSM;
  metadata.interpolation=enableExperimentalDerivedElectric ?
      Earth::Field::InterpolationMode::CellCenteredLinearDerivedElectric :
      Earth::Field::InterpolationMode::CellCenteredLinear;
  metadata.domain=snapshot.domain;
  metadata.magneticFieldAvailable=true;
  metadata.electricFieldAvailable=enableExperimentalDerivedElectric;
  metadata.immutableDuringBatch=true;
  metadata.valid=true;

  GlobalSnapshotMetadata_=metadata;
  GlobalPlasmaVelocityAvailable_=true;
  GlobalContentFingerprint_=snapshot.contentFingerprint;
  GlobalMeshRevision_=snapshot.meshRevision;
  GlobalSourceSimulationTime_s_=snapshot.simulationTime_s;
  ++GlobalSnapshotGeneration_;
  GlobalFieldsReady_=true;

  MaterializationStats stats;
  stats.usedLeafBlocks=nUsedLeafBlocks;
  stats.ownerInteriorCells=GlobalInteriorCellCount_;
  stats.expectedInteriorCells=GlobalInteriorCellCount_;
  stats.magneticFieldBytes=static_cast<long int>(GlobalMagneticField_.size()*sizeof(double));
  stats.electricFieldBytes=static_cast<long int>(GlobalElectricField_.size()*sizeof(double));
  stats.plasmaVelocityBytes=static_cast<long int>(GlobalPlasmaVelocity_.size()*sizeof(double));
  stats.electricFieldDerivedFromVelocity=enableExperimentalDerivedElectric;
  stats.plasmaVelocityAvailable=true;
  stats.snapshotId=snapshot.snapshotId;
  stats.contentFingerprint=snapshot.contentFingerprint;
  stats.meshRevision=snapshot.meshRevision;
  // Replay has no live owner buffer to re-read. Its equivalent gates are the exact
  // canonical cell mapping and mesh-revision checks above, so do not mislabel them as
  // an owner/coupler parity measurement.
  stats.ownerCellParityValidated=false;

  if (verbose && PIC::ThisThread==0) {
    std::cout << "[Mode3D::GlobalMagneticField] Imported frozen SWMF snapshot: file="
              << fileName << ", cells=" << GlobalInteriorCellCount_
              << ", snapshot=" << snapshot.snapshotId
              << ", content=" << snapshot.contentFingerprint << ".\n";
    std::cout.flush();
  }
  return stats;
}

long int RedefineGlobalMagneticField(
    const char* diagnosticTag,
    void (*fieldCallback)(double*,double*),
    bool verbose) {

  const std::string tag=SafeTag_(diagnosticTag);
  if (fieldCallback==NULL) return 0;
  RequireNoFrozenBatch_("redefine the active magnetic field");

  if ((PIC::Mesh::mesh==NULL)||(PIC::Mesh::mesh->rootTree==NULL)) {
    const std::string msg="["+tag+"] global magnetic-field redefinition called before the AMR tree is initialized.";
    exit(__LINE__,__FILE__,msg.c_str());
  }

  GlobalFieldsReady_=false;
  ResetTreeTempIds_(PIC::Mesh::mesh->rootTree);

  long int nUsedLeafBlocks=0;
  AssignGlobalLeafTempIds_(PIC::Mesh::mesh->rootTree,nUsedLeafBlocks);

  std::vector<cAMRNode*> nodes;
  nodes.reserve(static_cast<size_t>(nUsedLeafBlocks));
  CollectUsedLeafNodes_(PIC::Mesh::mesh->rootTree,nodes);

  const long int nInteriorCells=nUsedLeafBlocks*InteriorCellsPerBlock_();
  GlobalUsedLeafBlocks_=nUsedLeafBlocks;
  GlobalInteriorCellCount_=nInteriorCells;
  GlobalMagneticField_.assign(static_cast<size_t>(3*nInteriorCells),0.0);
  GlobalElectricField_.assign(static_cast<size_t>(3*nInteriorCells),0.0);
  GlobalPlasmaVelocity_.clear();
  GlobalCellPresence_.assign(static_cast<size_t>(nInteriorCells),1);

  double x[3],b[3];

  for (std::vector<cAMRNode*>::const_iterator it=nodes.begin();it!=nodes.end();++it) {
    cAMRNode* node=*it;
    const double dx[3]={
      (node->xmax[0]-node->xmin[0])/_BLOCK_CELLS_X_,
      (node->xmax[1]-node->xmin[1])/_BLOCK_CELLS_Y_,
      (node->xmax[2]-node->xmin[2])/_BLOCK_CELLS_Z_};

    for (int i=0;i<_BLOCK_CELLS_X_;i++) {
      for (int j=0;j<_BLOCK_CELLS_Y_;j++) {
        for (int k=0;k<_BLOCK_CELLS_Z_;k++) {
          x[0]=node->xmin[0]+(i+0.5)*dx[0];
          x[1]=node->xmin[1]+(j+0.5)*dx[1];
          x[2]=node->xmin[2]+(k+0.5)*dx[2];
          fieldCallback(x,b);

          if (!Earth::Field::FiniteVector3(b)) {
            std::ostringstream msg;
            msg << "[" << tag << "] callback returned non-finite B at x=("
                << x[0] << ',' << x[1] << ',' << x[2] << ").";
            const std::string text=msg.str();
            exit(__LINE__,__FILE__,text.c_str());
          }

          const long int c=GlobalCellIndex_(node,i,j,k);
          GlobalMagneticField_[static_cast<size_t>(3*c+0)]=b[0];
          GlobalMagneticField_[static_cast<size_t>(3*c+1)]=b[1];
          GlobalMagneticField_[static_cast<size_t>(3*c+2)]=b[2];
        }
      }
    }
  }

  // A callback replacement is a new physical generation.  Preserve useful domain and
  // epoch provenance when it exists, but never leave an earlier snapshot ID attached
  // to new array values.  E is unavailable after replacement because its old values
  // were explicitly reset to zero and are no longer physically synchronized with B.
  const std::string previousId=GlobalSnapshotMetadata_.snapshotId;
  Earth::Field::SnapshotMetadata replacement=GlobalSnapshotMetadata_;
  replacement.sourceId=tag+":CALLBACK_REDEFINE";
  replacement.modelName="CALLBACK";
  if (replacement.epochUTC.empty()) replacement.epochUTC="UNSPECIFIED";
  if (replacement.frame==Earth::Field::CoordinateFrame::Unknown)
    replacement.frame=Earth::Field::CoordinateFrame::GSM;
  replacement.interpolation=Earth::Field::InterpolationMode::CellCenteredLinear;
  replacement.magneticFieldAvailable=true;
  replacement.electricFieldAvailable=false;
  replacement.immutableDuringBatch=true;
  replacement.valid=true;
  replacement.validityMessage.clear();
  replacement.snapshotId=Earth::Field::MakeSnapshotId(
      replacement.sourceId,replacement.epochUTC,
      previousId+"|callback-redefine");
  GlobalSnapshotMetadata_=replacement;
  GlobalPlasmaVelocityAvailable_=false;
  GlobalContentFingerprint_.clear();
  GlobalMeshRevision_.clear();
  GlobalSourceSimulationTime_s_=0.0;
  ++GlobalSnapshotGeneration_;
  GlobalFieldsReady_=true;

  if ((verbose==true)&&(PIC::ThisThread==0)) {
    std::cout << "[" << tag << "] Replaced compact global B field from callback:"
              << " usedLeafBlocks=" << nUsedLeafBlocks
              << ", interiorCells=" << nInteriorCells
              << ". E was reset to zero; no AMR blocks were allocated.\n";
    std::cout.flush();
  }

  return nInteriorCells;
}

long int RedefineAllAllocatedMagneticField(
    const char* diagnosticTag,
    long int magneticFieldDataOffset,
    void (*fieldCallback)(double*,double*),
    bool allocateMissingBlocks,
    bool verbose) {

  // Preserve the old function signature for source compatibility.  The two arguments
  // below described writes into replicated AMPS cell buffers; compact-array operation
  // deliberately ignores them and never allocates a nonlocal block.
  (void)magneticFieldDataOffset;
  (void)allocateMissingBlocks;
  return RedefineGlobalMagneticField(diagnosticTag,fieldCallback,verbose);
}

} // namespace GlobalMagneticField
} // namespace Mode3D
} // namespace Earth
