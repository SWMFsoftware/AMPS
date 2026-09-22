//======================================================================================
// GlobalMagneticField.cpp
//======================================================================================
//
// Compact global cell-centered B/E storage for Mode3D backward trajectory tracing.
//
// Only owner-rank interior cells are packed.  MPI_Allreduce then creates identical
// compact arrays on every process.  The AMR tree itself is already globally replicated
// by AMPS, so node->Temp_ID plus local interior cell indices provide a complete global
// address without allocating remote cDataBlockAMR objects or reconstructing ghost-cell
// buffers.  Field evaluation is performed with cRowStencil, whose entries identify the
// physical owning leaf and interior (i,j,k) even for unallocated remote blocks.
//======================================================================================

#include "GlobalMagneticField.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
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
std::vector<int> GlobalCellPresence_;
long int GlobalUsedLeafBlocks_=0;
long int GlobalInteriorCellCount_=0;
bool GlobalFieldsReady_=false;
Earth::Field::SnapshotMetadata GlobalSnapshotMetadata_;
unsigned long long GlobalSnapshotGeneration_=0;

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

void Cross_(const double* a,const double* b,double* c) {
  c[0]=a[1]*b[2]-a[2]*b[1];
  c[1]=a[2]*b[0]-a[0]*b[2];
  c[2]=a[0]*b[1]-a[1]*b[0];
}

// Pack only authoritative owner-rank interior cells.  Nonowner blocks can be present
// because AMPS keeps a local neighbor layer, but those copies may contain stale ghost
// data and must never contribute to the global snapshot.
long int PackOwnedInteriorFields_(
    const std::vector<cAMRNode*>& nodes,
    long int magneticFieldDataOffset,
    long int electricFieldDataOffset,
    long int plasmaVelocityDataOffset,
    std::vector<double>& magneticField,
    std::vector<double>& electricField,
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

          if (electricFieldDataOffset>=0) {
            // Standalone DATAFILE path: E was initialized explicitly in the same
            // cell-associated buffer as B.
            std::memcpy(&electricField[static_cast<size_t>(vectorIndex)],
                        data+electricFieldDataOffset,3*sizeof(double));
          }
          else if (plasmaVelocityDataOffset>=0) {
            // SWMF path: E is not stored as an independent field.  Reconstruct the
            // cell-centered ideal-MHD field from the imported plasma velocity and B.
            double velocity[3],b[3],vxB[3];
            std::memcpy(velocity,data+plasmaVelocityDataOffset,3*sizeof(double));
            std::memcpy(b,data+magneticFieldDataOffset,3*sizeof(double));
            Cross_(velocity,b,vxB);

            electricField[static_cast<size_t>(vectorIndex+0)]=-vxB[0];
            electricField[static_cast<size_t>(vectorIndex+1)]=-vxB[1];
            electricField[static_cast<size_t>(vectorIndex+2)]=-vxB[2];
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
          if (!Earth::Field::FiniteVector3(packedB) ||
              !Earth::Field::FiniteVector3(packedE)) {
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
    bool verbose) {

  const std::string tag=SafeTag_(diagnosticTag);
  ValidateMeshAndOffset_(tag,magneticFieldDataOffset);

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

  const bool electricSourceAvailable=
      (electricFieldDataOffset>=0 || plasmaVelocityDataOffset>=0);
  if (metadata.electricFieldAvailable && !electricSourceAvailable) {
    const std::string msg="["+tag+
        "] metadata advertises E, but no E or plasma-velocity source was supplied.";
    exit(__LINE__,__FILE__,msg.c_str());
  }
  const Earth::Field::InterpolationMode expectedInterpolation=
      (electricFieldDataOffset<0 && plasmaVelocityDataOffset>=0) ?
      Earth::Field::InterpolationMode::CellCenteredLinearDerivedElectric :
      Earth::Field::InterpolationMode::CellCenteredLinear;
  if (metadata.interpolation!=expectedInterpolation) {
    const std::string msg="["+tag+
        "] metadata interpolation mode does not match compact-array assembly.";
    exit(__LINE__,__FILE__,msg.c_str());
  }

  MaterializationStats stats;

  // Invalidate the previous snapshot before touching Temp_ID or resizing arrays.  A
  // concurrent field evaluation during assembly would otherwise observe a mixture of
  // old indices and new storage.  Assembly is expected to run before worker threads.
  GlobalFieldsReady_=false;

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
  if (GlobalCellPresence_.size()!=nPresenceValues)
    GlobalCellPresence_.resize(nPresenceValues);
  std::fill(GlobalMagneticField_.begin(),GlobalMagneticField_.end(),0.0);
  std::fill(GlobalElectricField_.begin(),GlobalElectricField_.end(),0.0);
  std::fill(GlobalCellPresence_.begin(),GlobalCellPresence_.end(),0);

  const long int nPackedLocal=PackOwnedInteriorFields_(
      nodes,magneticFieldDataOffset,electricFieldDataOffset,
      plasmaVelocityDataOffset,GlobalMagneticField_,GlobalElectricField_,
      GlobalCellPresence_);

  // Replicate only the compact physical fields and one integer validity entry per
  // interior cell.  No cDataBlockAMR, center-node state vector, or ghost cell is
  // allocated by this operation.
  AllreduceDoubleVectorInPlace_(GlobalMagneticField_);
  AllreduceDoubleVectorInPlace_(GlobalElectricField_);
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
      }
    }

    const double* reducedB=&GlobalMagneticField_[static_cast<size_t>(3*c)];
    const double* reducedE=&GlobalElectricField_[static_cast<size_t>(3*c)];
    if (!Earth::Field::FiniteVector3(reducedB) ||
        !Earth::Field::FiniteVector3(reducedE)) {
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
  stats.electricFieldReadFromBuffer=(electricFieldDataOffset>=0);
  stats.electricFieldDerivedFromVelocity=
      ((electricFieldDataOffset<0)&&(plasmaVelocityDataOffset>=0));
  stats.snapshotId=metadata.snapshotId;

  if ((nPackedGlobal!=nInteriorCells)||(nMissing!=0)||(nDuplicate!=0)) {
    std::ostringstream msg;
    msg << "[" << tag << "] compact global field gather is inconsistent: ownerCells="
        << nPackedGlobal << ", expected=" << nInteriorCells
        << ", missing=" << nMissing << ", duplicate=" << nDuplicate
        << ". Every used interior cell must be supplied by exactly one AMPS owner rank.";
    const std::string text=msg.str();
    exit(__LINE__,__FILE__,text.c_str());
  }

  // Publish metadata and generation only after ownership, completeness, and finite-
  // value gates pass.  A consumer can never observe new identity with old arrays.
  GlobalSnapshotMetadata_=metadata;
  ++GlobalSnapshotGeneration_;
  GlobalFieldsReady_=true;

  if ((verbose==true)&&(PIC::ThisThread==0)) {
    const double mib=1024.0*1024.0;
    std::cout << "[" << tag << "] Prepared compact global cell-centered fields:"
              << " usedLeafBlocks=" << stats.usedLeafBlocks
              << ", interiorCells=" << stats.expectedInteriorCells
              << ", B=" << stats.magneticFieldBytes/mib << " MiB"
              << ", E=" << stats.electricFieldBytes/mib << " MiB"
              << ", snapshot=" << stats.snapshotId
              << ", epoch=" << metadata.epochUTC
              << ", frame=" << Earth::Field::CoordinateFrameName(metadata.frame)
              << ", ESource=";

    if (stats.electricFieldReadFromBuffer) std::cout << "cell buffer";
    else if (stats.electricFieldDerivedFromVelocity) std::cout << "-VxB";
    else std::cout << "zero";

    std::cout << ". No nonlocal AMR blocks were allocated.\n";
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
  GlobalFieldsReady_=false;
  GlobalUsedLeafBlocks_=0;
  GlobalInteriorCellCount_=0;
  GlobalMagneticField_.clear();
  GlobalElectricField_.clear();
  GlobalCellPresence_.clear();
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

long int RedefineGlobalMagneticField(
    const char* diagnosticTag,
    void (*fieldCallback)(double*,double*),
    bool verbose) {

  const std::string tag=SafeTag_(diagnosticTag);
  if (fieldCallback==NULL) return 0;

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
