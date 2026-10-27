//======================================================================================
// Mode3DMeshExport.cpp -- see Mode3DMeshExport.h.
//
// This file is compiled with the AMPS C++ compiler as part of the application library
// (it needs pic.h for the AMR tree); it never runs in a device kernel.
//======================================================================================

#include "Mode3DMeshExport.h"

#include "BatchTracer.h"
#include "EngineGlue.h"
#include "Fiber.h"
#include "../3d/GlobalMagneticField.h"

#include "pic.h"

#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <mpi.h>

namespace Geospace {
namespace GPU {
namespace Mode3DMesh {

namespace {

typedef cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> Node;

// Upload bookkeeping.  gResidentGeneration is read lock-free by IsCurrent() on every
// batched mesh trajectory; the remaining fields are touched only under gUploadMutex.
std::mutex gUploadMutex;
std::atomic<unsigned long long> gResidentGeneration(0);
unsigned long long gResidentChecksum=0;
const void* gResidentMesh=nullptr;

bool EnvFlag(const char* name) {
  const char* e=std::getenv(name);
  return e!=nullptr && std::string(e)=="1";
}

long EnvLong(const char* name,long fallback) {
  if (const char* e=std::getenv(name)) {
    const long v=std::atol(e);
    if (v>0) return v;
  }
  return fallback;
}

int WorldRank() {
  int initialized=0,finalized=0;
  MPI_Initialized(&initialized);
  if (initialized) MPI_Finalized(&finalized);
  if (!initialized || finalized) return 0;
  int rank=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  return rank;
}

// FNV-1a over raw bytes; used only to detect a changed tree/Temp_ID assignment for
// an unchanged field generation (defensive: the application never does this while a
// generation is published, but a silent mismatch would corrupt every trajectory).
struct Fnv {
  unsigned long long h=1469598103934665603ULL;
  void Add(const void* p,std::size_t n) {
    const unsigned char* c=static_cast<const unsigned char*>(p);
    for (std::size_t i=0;i<n;++i) { h^=c[i]; h*=1099511628211ULL; }
  }
  template <class T> void AddVec(const std::vector<T>& v) {
    if (!v.empty()) Add(v.data(),v.size()*sizeof(T));
  }
};

// Flattened copy of the AMR tree in the layout of EgMeshHost (EnginePool.h).
struct FlatMesh {
  std::vector<Node*> nodes;
  std::unordered_map<const Node*,int> index;
  std::vector<double> xmin,xmax;
  std::vector<int> ixMin,sizeIndex,up,down,level,minL,maxL;
  std::vector<unsigned char> used,presence;
  std::vector<long long> tempId;
  unsigned long long checksum=0;
  EgMeshHost host{};
};

// Pre-order depth-first traversal through downNode[0..7], exactly the set of nodes
// cMeshAMRgeneric::findTreeNode(int*) can visit (it follows downNode[] while non-NULL
// and upNode[] while non-NULL).  The traversal order only defines the engine node
// numbering; the search result (a node) is independent of it.
void CollectNodes(Node* root,std::vector<Node*>& out) {
  std::vector<Node*> stack;
  if (root!=nullptr) stack.push_back(root);
  while (!stack.empty()) {
    Node* n=stack.back();
    stack.pop_back();
    out.push_back(n);
    for (int c=(1<<_MESH_DIMENSION_)-1;c>=0;--c)
      if (n->downNode[c]!=nullptr) stack.push_back(n->downNode[c]);
  }
}

void Flatten(FlatMesh& f,const Geospace::Mode3D::GlobalMagneticField::CompactMagneticFieldView& v) {
#if _MESH_DIMENSION_ != 3
  throw std::logic_error("Mode3DMeshExport: the engine mesh port requires a 3-D AMPS mesh");
#endif
  auto* mesh=PIC::Mesh::mesh;
  CollectNodes(mesh->rootTree,f.nodes);
  const std::size_t nn=f.nodes.size();
  if (nn>static_cast<std::size_t>(0x7fffffff))
    throw std::runtime_error("Mode3DMeshExport: AMR tree has more than INT_MAX nodes");
  f.index.reserve(nn*2);
  for (std::size_t i=0;i<nn;++i) f.index.emplace(f.nodes[i],static_cast<int>(i));

  f.xmin.resize(3*nn); f.xmax.resize(3*nn); f.ixMin.resize(3*nn);
  f.sizeIndex.resize(nn); f.up.resize(nn); f.down.resize(8*nn);
  f.level.resize(nn); f.minL.resize(nn); f.maxL.resize(nn);
  f.used.resize(nn); f.tempId.resize(nn);

  auto IndexOf=[&](const Node* p) -> int {
    if (p==nullptr) return -1;
    auto it=f.index.find(p);
    if (it==f.index.end())
      throw std::runtime_error("Mode3DMeshExport: tree link to a node that is not "
                               "reachable from rootTree");
    return it->second;
  };

  for (std::size_t i=0;i<nn;++i) {
    const Node* n=f.nodes[i];
    for (int d=0;d<3;++d) {
      f.xmin[3*i+d]=n->xmin[d];
      f.xmax[3*i+d]=n->xmax[d];
      f.ixMin[3*i+d]=n->xMinGlobalIndex[d];
    }
    f.sizeIndex[i]=n->NodeGeometricSizeIndex;
    f.up[i]=IndexOf(n->upNode);
    for (int c=0;c<8;++c) f.down[8*i+c]=IndexOf(n->downNode[c]);
    // char -> int preserves the value; InitStencil compares these levels only with
    // each other (==, <), so the comparison results are unchanged.
    f.level[i]=static_cast<int>(n->RefinmentLevel);
    f.minL[i]=static_cast<int>(n->minNeibRefinmentLevel);
    f.maxL[i]=static_cast<int>(n->maxNeibRefinmentLevel);
    f.used[i]=n->IsUsedInCalculationFlag ? 1 : 0;
    f.tempId[i]=static_cast<long long>(n->Temp_ID);
  }

  // CheckedGlobalCellIndex_ accepts a cell only when GlobalCellPresence_[c]==1 (the
  // owner count after the Allreduce).  Map that predicate exactly to 0/1 so a count
  // such as 257 cannot alias to 1 in an unsigned char.
  const std::size_t nc=static_cast<std::size_t>(v.interiorCells);
  f.presence.resize(nc);
  for (std::size_t c=0;c<nc;++c) f.presence[c]=(v.cellPresence[c]==1) ? 1 : 0;

  EgMeshHost& h=f.host;
  h.nNodes=static_cast<int>(nn);
  h.root=nn>0 ? 0 : -1;
  h.nCells[0]=_BLOCK_CELLS_X_; h.nCells[1]=_BLOCK_CELLS_Y_; h.nCells[2]=_BLOCK_CELLS_Z_;
  for (int d=0;d<3;++d) {
    h.xGlobalMin[d]=mesh->xGlobalMin[d];
    h.xGlobalMax[d]=mesh->xGlobalMax[d];
    h.dxMaxRefinement[d]=mesh->dx_max_refinment[d];
  }
  h.xmin=f.xmin.data(); h.xmax=f.xmax.data(); h.ixMin=f.ixMin.data();
  h.sizeIndex=f.sizeIndex.data(); h.up=f.up.data(); h.down=f.down.data();
  h.level=f.level.data(); h.minNeibLevel=f.minL.data(); h.maxNeibLevel=f.maxL.data();
  h.used=f.used.data(); h.tempId=f.tempId.data();
  h.usedLeafBlocks=static_cast<long long>(v.usedLeafBlocks);
  h.interiorCells=static_cast<long long>(v.interiorCells);
  // Not copied here: BatchTracer::SetMesh copies B into the host mirror and the device
  // pool.  The pointer stays valid for the published generation.
  h.B=v.magneticField;
  h.presence=f.presence.data();

  Fnv fnv;
  fnv.Add(h.nCells,sizeof(h.nCells));
  fnv.Add(h.xGlobalMin,sizeof(h.xGlobalMin));
  fnv.Add(h.xGlobalMax,sizeof(h.xGlobalMax));
  fnv.Add(h.dxMaxRefinement,sizeof(h.dxMaxRefinement));
  fnv.AddVec(f.xmin); fnv.AddVec(f.xmax); fnv.AddVec(f.ixMin); fnv.AddVec(f.sizeIndex);
  fnv.AddVec(f.up); fnv.AddVec(f.down); fnv.AddVec(f.level); fnv.AddVec(f.minL);
  fnv.AddVec(f.maxL); fnv.AddVec(f.used); fnv.AddVec(f.tempId);
  f.checksum=fnv.h;
}

// Deterministic generator for the self-check sample points (splitmix64).
struct SplitMix {
  unsigned long long s;
  explicit SplitMix(unsigned long long seed) : s(seed) {}
  unsigned long long Next() {
    unsigned long long z=(s+=0x9E3779B97F4A7C15ULL);
    z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL;
    z=(z^(z>>27))*0x94D049BB133111EBULL;
    return z^(z>>31);
  }
  double Uniform() { return (Next()>>11)*(1.0/9007199254740992.0); }
};

bool SameBits(double a,double b) { return std::memcmp(&a,&b,sizeof(double))==0; }

// Bitwise comparison of the engine mesh port with the AMPS/legacy routines (see the
// header).  Runs on the calling (main) thread after the upload; read-only with respect
// to every AMPS/legacy structure, except that InterpolateMagneticField() may exit()
// on a compact-array consistency error exactly as it would during tracing.
void SelfCheck(const FlatMesh& f,
               const Geospace::Mode3D::GlobalMagneticField::CompactMagneticFieldView& v,
               const char* label) {
#if defined(GEOSPACE_PRECISION_FLOAT)
  // Float/normalized build: bitwise identity is not defined.  Report instead how the
  // float engine field compares with AMPS (double): node-search disagreements, samples
  // where AMPS interpolates but the engine reports the field unavailable (or vice
  // versa), and the maximum relative |B| difference.
  const auto t0=std::chrono::steady_clock::now();
  HostMesh hm;
  hm.Set(f.host);
  const EgMesh<Real> em(hm.View());
  EgRowStencil<Real> row,fineRow;
  auto* mesh=PIC::Mesh::mesh;
  const long nSamples=EnvLong("AMPS_EARTH_MESH_VERIFY_SAMPLES",60000);
  SplitMix rng(0x5EED0000ULL+v.generation);
  long long nodeDiff=0,unavailable=0,extra=0,compared=0,errors=0;
  double maxRel=0.0;
  for (long s=0;s<nSamples;++s) {
    double x[3];
    for (int d=0;d<3;++d) {
      const double lo=f.host.xGlobalMin[d],hi=f.host.xGlobalMax[d],w=hi-lo;
      x[d]=lo+w*rng.Uniform();
    }
    Node* ref=mesh->findTreeNode(x);
    const Real xi[3]={Real(ToLength<Units>(x[0])),Real(ToLength<Units>(x[1])),
                      Real(ToLength<Units>(x[2]))};
    const int engIdx=em.Find(xi,-1);
    const int refIdx=(ref==nullptr) ? -1 : f.index.at(ref);
    if (engIdx!=refIdx) ++nodeDiff;
    if (ref==nullptr || engIdx<0) continue;
    int err=kErrNone;
    Real Be[3]={0,0,0};
    const bool okEngine=em.Interpolate(xi,engIdx,Be,row,fineRow,err);
    if (err!=kErrNone) { ++errors; continue; }
    double Bref[3]={0.0,0.0,0.0};
    const bool okRef=Geospace::Mode3D::GlobalMagneticField::InterpolateMagneticField(x,ref,Bref);
    if (okRef && !okEngine) { ++unavailable; continue; }
    if (!okRef && okEngine) { ++extra; continue; }
    if (!okRef) continue;
    ++compared;
    double dn=0.0,rn=0.0;
    for (int d=0;d<3;++d) {
      const double be=FromField<Units>(double(Be[d]));
      dn+=(be-Bref[d])*(be-Bref[d]);
      rn+=Bref[d]*Bref[d];
    }
    if (rn>0.0) maxRel=std::max(maxRel,std::sqrt(dn/rn));
  }
  const double dt=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
  if (WorldRank()==0 || unavailable>0 || errors>0)
    std::printf("[earth-engine][mesh-verify][float] %s: rank=%d samples=%ld compared=%lld "
                "node-differs=%lld engine-unavailable=%lld engine-extra=%lld errors=%lld "
                "max-rel-|dB|=%.3e (%.2fs)\n",label!=nullptr ? label : "",WorldRank(),
                nSamples,compared,nodeDiff,unavailable,extra,errors,maxRel,dt);
  std::fflush(stdout);
#else
  const auto t0=std::chrono::steady_clock::now();
  HostMesh hm;
  hm.Set(f.host);
  const EgMesh<Real> em(hm.View());
  EgRowStencil<Real> row,fineRow;
  auto* mesh=PIC::Mesh::mesh;

  std::vector<int> leaves;
  for (std::size_t i=0;i<f.nodes.size();++i)
    if (f.down[8*i+0]<0 && f.down[8*i+1]<0 && f.down[8*i+2]<0 && f.down[8*i+3]<0 &&
        f.down[8*i+4]<0 && f.down[8*i+5]<0 && f.down[8*i+6]<0 && f.down[8*i+7]<0)
      leaves.push_back(static_cast<int>(i));

  const long nSamples=EnvLong("AMPS_EARTH_MESH_VERIFY_SAMPLES",60000);
  SplitMix rng(0x5EED0000ULL+v.generation);
  long long findMismatch=0,hintMismatch=0,interpMismatch=0,interpolated=0,nullNodes=0;
  int printed=0;
  const int rank=WorldRank();

  auto Report=[&](const char* what,const double x[3],const std::string& detail) {
    if (printed>=10) return;
    ++printed;
    std::printf("[earth-engine][mesh-verify][rank %d] %s mismatch at x=(%.17g,%.17g,%.17g): %s\n",
                rank,what,x[0],x[1],x[2],detail.c_str());
  };

  for (long s=0;s<nSamples;++s) {
    double x[3];
    if (s%3==0 || leaves.empty()) {
      // Uniform over the global AMR box enlarged by 1% (also exercises NULL results).
      for (int d=0;d<3;++d) {
        const double lo=f.host.xGlobalMin[d],hi=f.host.xGlobalMax[d],w=hi-lo;
        x[d]=lo-0.01*w+1.02*w*rng.Uniform();
      }
    }
    else {
      // Leaf-relative special positions: faces, near-face cell-centre lines, the
      // one-cell transition bands of InitStencil, and points just outside the leaf.
      const int leaf=leaves[rng.Next()%leaves.size()];
      for (int d=0;d<3;++d) {
        const double n=static_cast<double>(f.host.nCells[d]);
        const double table[16]={-1.0e-7,0.0,1.0e-12,0.25/n,0.5/n,0.75/n,1.0/n,1.25/n,
                                rng.Uniform(),1.0-1.25/n,1.0-1.0/n,1.0-0.5/n,1.0-0.25/n,
                                1.0-1.0e-12,1.0,1.0+1.0e-7};
        const double frac=table[rng.Next()%16];
        const double lo=f.xmin[3*leaf+d],hi=f.xmax[3*leaf+d];
        x[d]=lo+frac*(hi-lo);
      }
    }

    // 1. findTreeNode(x) from the root.
    Node* ref=mesh->findTreeNode(x);
    const int refIdx=(ref==nullptr) ? -1 : f.index.at(ref);
    const int engIdx=em.Find(x,-1);
    if (engIdx!=refIdx) {
      ++findMismatch;
      Report("find",x,"engine="+std::to_string(engIdx)+" amps="+std::to_string(refIdx));
    }

    // 2. findTreeNode(x,hint) with a random start node (the legacy evaluator passes lastNode_).
    const int hint=static_cast<int>(rng.Next()%f.nodes.size());
    Node* refHint=mesh->findTreeNode(x,f.nodes[hint]);
    const int refHintIdx=(refHint==nullptr) ? -1 : f.index.at(refHint);
    const int engHintIdx=em.Find(x,hint);
    if (engHintIdx!=refHintIdx) {
      ++hintMismatch;
      Report("find(hint)",x,"engine="+std::to_string(engHintIdx)+" amps="+
             std::to_string(refHintIdx));
    }

    if (ref==nullptr) { ++nullNodes; continue; }

    // 3. InterpolateMagneticField(x,node) -- the exact call of cMode3DMeshFieldEval.
    int err=kErrNone;
    Real Be[3]={0,0,0};
    const bool okEngine=em.Interpolate(x,refIdx,Be,row,fineRow,err);
    if (err!=kErrNone) {
      ++interpMismatch;
      Report("interpolate",x,"engine reported mesh error "+std::to_string(err)+
             " (AMPS would terminate in CheckedGlobalCellIndex_/AddCell)");
      continue;
    }
    double Bref[3]={0.0,0.0,0.0};
    const bool okRef=Geospace::Mode3D::GlobalMagneticField::InterpolateMagneticField(x,ref,Bref);
    ++interpolated;
    bool same=(okEngine==okRef);
    if (same && okRef)
      same=SameBits(Be[0],Bref[0]) && SameBits(Be[1],Bref[1]) && SameBits(Be[2],Bref[2]);
    if (!same) {
      ++interpMismatch;
      char buf[320];
      std::snprintf(buf,sizeof(buf),"engine ok=%d B=(%.17g,%.17g,%.17g) amps ok=%d "
                    "B=(%.17g,%.17g,%.17g) rowLength=%d",okEngine?1:0,Be[0],Be[1],Be[2],
                    okRef?1:0,Bref[0],Bref[1],Bref[2],row.Length);
      Report("interpolate",x,buf);
    }
  }

  const double dt=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
  const bool bad=(findMismatch+hintMismatch+interpMismatch)>0;
  if (rank==0 || bad) {
    std::printf("[earth-engine][mesh-verify] %s: rank=%d generation=%llu nodes=%d leaves=%zu "
                "samples=%ld null=%lld interpolated=%lld find-mismatch=%lld "
                "hint-mismatch=%lld interp-mismatch=%lld (%.2fs)\n",
                label!=nullptr ? label : "",rank,v.generation,f.host.nNodes,leaves.size(),
                nSamples,nullNodes,interpolated,findMismatch,hintMismatch,interpMismatch,dt);
  }
  std::fflush(stdout);
#endif
}

} // namespace

bool IsCurrent() {
  const unsigned long long g=gResidentGeneration.load(std::memory_order_acquire);
  return g!=0 && g==Geospace::Mode3D::GlobalMagneticField::CurrentSnapshotGeneration() &&
         Geospace::Mode3D::GlobalMagneticField::GlobalFieldsReady();
}

unsigned long long ResidentGeneration() {
  return gResidentGeneration.load(std::memory_order_acquire);
}

bool EnsureUploaded(const char* label) {
  if (!BatchTracer::Enabled()) return false;
  // A fiber can never upload (trajectories of other fibers may be in flight); it may
  // only observe whether the resident copy is current.
  if (Fiber::InFiber()) return IsCurrent();

  std::lock_guard<std::mutex> lock(gUploadMutex);
  const Geospace::Mode3D::GlobalMagneticField::CompactMagneticFieldView view=
      Geospace::Mode3D::GlobalMagneticField::GetCompactMagneticFieldView();
  if (!view.ready || PIC::Mesh::mesh==nullptr || PIC::Mesh::mesh->rootTree==nullptr ||
      view.magneticField==nullptr || view.cellPresence==nullptr)
    return false;
  if (view.magneticFieldValues!=3*view.interiorCells ||
      view.cellPresenceValues!=view.interiorCells)
    throw std::runtime_error("Mode3DMeshExport: compact field arrays have inconsistent sizes");

  const auto t0=std::chrono::steady_clock::now();
  FlatMesh f;
  Flatten(f,view);
  if (gResidentGeneration.load(std::memory_order_acquire)==view.generation &&
      gResidentChecksum==f.checksum && gResidentMesh==PIC::Mesh::mesh)
    return true;

  // The generation is published only after the new values are visible; clear it first
  // so a concurrent IsCurrent() never accepts a half-replaced copy.
  gResidentGeneration.store(0,std::memory_order_release);
  // The flattened arrays stay resident (moving the vectors keeps their buffers, so the
  // pointers in f.host remain valid): CPU double/SI builds alias them, and the compact
  // B array of GlobalMagneticField, instead of copying (BatchTracer::SetMesh).
  static FlatMesh resident;
  resident=std::move(f);
  BatchTracer::Get().SetMesh(resident.host,view.generation);
  gResidentChecksum=resident.checksum;
  gResidentMesh=PIC::Mesh::mesh;
  gResidentGeneration.store(view.generation,std::memory_order_release);

  if (EnvFlag("AMPS_EARTH_ENGINE_VERBOSE") && WorldRank()==0) {
    const double dt=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    std::printf("[earth-engine] mesh upload (%s): generation=%llu nodes=%d usedLeafBlocks=%lld "
                "interiorCells=%lld (%.3fs)\n",label!=nullptr ? label : "",view.generation,
                resident.host.nNodes,resident.host.usedLeafBlocks,resident.host.interiorCells,dt);
    std::fflush(stdout);
  }

  if (VerifyEnabled() || EnvFlag("AMPS_EARTH_MESH_VERIFY")) SelfCheck(resident,view,label);
  return true;
}

} // namespace Mode3DMesh
} // namespace GPU
} // namespace Geospace
