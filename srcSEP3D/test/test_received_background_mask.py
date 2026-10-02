"""Check native halo evidence against the actual AMPS packing mask.

Run from any directory: python3 srcSEP3D/test/test_received_background_mask.py
The fixture compiles the production layer-mask generator and the entire native
CaptureRuntimeMeshBackground function. Only mesh allocation and scalar wire
bytes are fixtures. No MPI behavior or successful native exchange is claimed.
"""
from pathlib import Path
import subprocess
import sys
import tempfile

APP_ROOT=Path(__file__).resolve().parents[1]
AMPS_ROOT=APP_ROOT.parent
sys.path.insert(0,str(APP_ROOT/'test'))
from test_runtime_background_coupling import cpp_block
h=(AMPS_ROOT/'src/pic/pic.h').read_text()
mask=(AMPS_ROOT/'src/pic/pic_block_send_mask.cpp').read_text()
app=(APP_ROOT/'main_lib.cpp').read_text()
mask_namespace = cpp_block(h, r'namespace BlockElementSendMask')
definitions = [cpp_block(mask, r'void PIC::Mesh::BlockElementSendMask::Set\([^\n]+\)'),
               cpp_block(mask, r'int PIC::Mesh::BlockElementSendMask::CornerNode::GetSize\(\)'),
               cpp_block(mask, r'int PIC::Mesh::BlockElementSendMask::CenterNode::GetSize\(\)'),
               cpp_block(mask, r'void PIC::Mesh::BlockElementSendMask::InitLayerBlockBasic\([^\n]+\)'),
               cpp_block(mask, r'void PIC::Mesh::BlockElementSendMask::InitLayerBlock\([^\n]+\)')]
capture = cpp_block(app, r'void CaptureRuntimeMeshBackground\([^\n]+\n[^\n]+\)')
fixture = r'''
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>
#define _TARGET_HOST_
#define _TARGET_DEVICE_
#define _CUDA_MANAGED_
#define DIM 3
#define _BLOCK_CELLS_X_ 4
#define _BLOCK_CELLS_Y_ 4
#define _BLOCK_CELLS_Z_ 4
#define _TOTAL_BLOCK_CELLS_X_ 4
#define _TOTAL_BLOCK_CELLS_Y_ 4
#define _TOTAL_BLOCK_CELLS_Z_ 4
void exit(int,const char*,const char*) { throw std::runtime_error("AMPS error"); }
template<class A,class B,class C> struct cAmpsMesh;
template<class B> struct cTreeNodeAMR {
  B* block=nullptr;
  int Thread=1,RefinmentLevel=0,face=1;
  bool IsUsedInCalculationFlag=true;
  double xmin[3]={40,0,0},xmax[3]={44,4,4};
  cTreeNodeAMR *nextBranchBottomNode=nullptr,*neighbor=nullptr;
  template<class M> cTreeNodeAMR* GetNeibFace(int f,int,int,M*) {return f==face?neighbor:nullptr;}
  template<class M> cTreeNodeAMR* GetNeibEdge(int,int,M*) {return nullptr;}
  template<class M> cTreeNodeAMR* GetNeibCorner(int,M*) {return nullptr;}
};
namespace PIC {
int ThisThread=0,nTotalThreads=2;
namespace Mesh {
struct cDataCenterNode { double value=0; };
struct cDataCornerNode {};
struct cDataBlockAMR {
  cDataCenterNode cells[64];
  cDataCenterNode* GetCenterNode(int n) {return &cells[n];}
};
cAmpsMesh<cDataCornerNode,cDataCenterNode,cDataBlockAMR>* mesh=nullptr;
int _getCenterNodeLocalNumber(int i,int j,int k) { return i+4*(j+4*k); }
int _getCornerNodeLocalNumber(int i,int j,int k) { return i+5*(j+5*k); }
@MASK_NAMESPACE@
}
}
template<class A,class B,class C> struct cAmpsMesh {
  cTreeNodeAMR<C>* BranchBottomNodeList=nullptr;
  struct Exchange {
    int* GlobalSendTable=nullptr;
    cTreeNodeAMR<C>*** RecvNodeTable=nullptr;
    unsigned char** RecvCenterNodePackingTable=nullptr;
    int BlockCenterNodeSendMaskLength=0;
  } ParallelBlockDataExchangeData;
  int getCenterNodeLocalNumber(int i,int j,int k) {return i+4*(j+4*k);}
};
int PIC::Mesh::BlockElementSendMask::CommunicationDepthSmall=1;
int PIC::Mesh::BlockElementSendMask::CommunicationDepthLarge=2;
@MASK_DEFINITIONS@
namespace SEP3D {
namespace Core {
struct Vec3 {
  double x=0,y=0,z=0;
  Vec3()=default;
  Vec3(double a,double b,double c):x(a),y(b),z(c){}
  Vec3 operator-(const Vec3& b)const {return {x-b.x,y-b.y,z-b.z};}
  double Norm()const {return std::sqrt(x*x+y*y+z*z);}
};
}
}
struct MockStatus {std::string message;};
struct Sample {bool valid=true;double value=120;MockStatus status;};
struct Metadata {double epochS=120;unsigned long long generation=3;};
struct Snapshot {
  std::vector<Sample> points{{true,120}};
  Metadata tag;
  const std::vector<Sample>& samples()const {return points;}
  const Metadata& metadata()const {return tag;}
} snapshot;
struct Provider {
  Metadata tag;
  const Metadata* PreparedMetadata()const{return &tag;}
  Sample Evaluate(const SEP3D::Core::Vec3&)const{return {true,120};}
} prepared;
Snapshot* gInstalledBackground=&snapshot;
Provider* gRuntimeBackgroundProvider=&prepared;
bool gNativeAmpsBackgroundReady=true;
struct Options {double innerRadiusM=1,outerRadiusM=100;SEP3D::Core::Vec3 coordinateOriginM;};
struct Config {Options values;const Options& options()const{return values;}} config;
const Config& Configuration(){return config;}
PIC::Mesh::cDataCenterNode owner;
struct AmpsCellReference {PIC::Mesh::cDataCenterNode* cell;SEP3D::Core::Vec3 positionM;};
std::vector<AmpsCellReference> CollectOwnedPhysicalCells(){return {{&owner,{10,0,0}}};}
// Byte comparison and wire transport are scalar fixtures. The production mask
// generator and the entire production capture/representative loop are intact.
bool MeshBackgroundBytesMatch(PIC::Mesh::cDataCenterNode* cell,const Sample& expected) {
  return std::memcmp(&cell->value,&expected.value,sizeof(double))==0;
}
@CAPTURE@
int main(){
  using namespace PIC::Mesh;
  cAmpsMesh<cDataCornerNode,cDataCenterNode,cDataBlockAMR> grid;
  mesh=&grid;
  cDataBlockAMR block;
  cTreeNodeAMR<cDataBlockAMR> remote,local;
  remote.block=&block;remote.neighbor=&local;local.Thread=0;
  grid.BranchBottomNodeList=&remote;
  owner.value=120;
  std::vector<unsigned char> center(BlockElementSendMask::CenterNode::GetSize());
  std::vector<unsigned char> corner(BlockElementSendMask::CornerNode::GetSize());
  int sends[4]={0,0,1,0};
  cTreeNodeAMR<cDataBlockAMR>* receivedNodes[1]={&remote};
  cTreeNodeAMR<cDataBlockAMR>** recvTables[2]={nullptr,receivedNodes};
  unsigned char* masks[2]={nullptr,center.data()};
  grid.ParallelBlockDataExchangeData={sends,recvTables,masks,static_cast<int>(center.size())};
  int verified=0;
  for(int face=0;face<6;++face){
    remote.face=face;
    for(auto& cell:block.cells)cell.value=0;
    BlockElementSendMask::InitLayerBlock(&remote,0,center.data(),corner.data());
    int received=0;
    for(int k=0;k<4;++k)for(int j=0;j<4;++j)for(int i=0;i<4;++i){
      if(!BlockElementSendMask::CenterNode::Test(i,j,k,center.data()))continue;
      block.GetCenterNode(grid.getCenterNodeLocalNumber(i,j,k))->value=120;
      ++received;
    }
    bool owned=false,ghosts=false,provider=false;unsigned long long count=0;
    CaptureRuntimeMeshBackground(&owned,&ghosts,&provider,&count);
    if(received!=16||!owned||!provider||count!=1)throw std::runtime_error("fixture contract failed");
    const bool selectedWasReceived=BlockElementSendMask::CenterNode::Test(0,0,0,center.data());
    if(!ghosts)throw std::runtime_error("a current received face falsely failed");
    ++verified;
    // A genuine stale transfer must remain a failure; fixing representative
    // selection must not weaken the primitive/derivative byte comparison.
    for(auto& cell:block.cells)cell.value=0;
    CaptureRuntimeMeshBackground(&owned,&ghosts,&provider,&count);
    if(ghosts||count!=1)throw std::runtime_error("stale received data was accepted");
    std::cout<<"face="<<face<<" received_cells="<<received
             <<" selected_(0,0,0)_received="<<selectedWasReceived
             <<" current_received_match=1 stale_received_rejected=1\n";
  }
  if(verified!=6)throw std::runtime_error("missing face coverage");
  // A null mask represents a full-block receive, not a zero-length transfer.
  masks[1]=nullptr;for(auto& cell:block.cells)cell.value=120;
  bool owned,ghosts,provider;unsigned long long count;
  CaptureRuntimeMeshBackground(&owned,&ghosts,&provider,&count);
  if(!ghosts||count!=1)throw std::runtime_error("unmasked receive failed");
  std::cout<<"PASS: six current face receives, six stale-data rejections and one full-block receive\n";
}
'''
fixture = fixture.replace('@MASK_NAMESPACE@', mask_namespace)
fixture = fixture.replace('@MASK_DEFINITIONS@', '\n'.join(definitions))
fixture = fixture.replace('@CAPTURE@', capture)
with tempfile.TemporaryDirectory(prefix='amps-ghost-mask-') as temp:
    source=Path(temp)/'fixture.cpp'
    source.write_text(fixture)
    executable=Path(temp)/'fixture'
    subprocess.run(['g++','-std=c++17','-O1','-Wall','-Wextra',str(source),'-o',str(executable)],check=True)
    run=subprocess.run([str(executable)],capture_output=True,text=True,check=True)
    print(run.stdout,end='')

