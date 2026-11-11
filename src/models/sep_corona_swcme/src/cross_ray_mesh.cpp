#include "sep_corona_swcme/cross_ray_mesh.h"

#include "sep_coronal_cme/configuration_parser.h"
#include "sep_coronal_cme/constants.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <locale>
#include <map>
#include <new>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace SEP { namespace CoronaSwcme { namespace CrossRay { namespace {

using CoronalCME::Norm;
using CoronalCME::Vec3;
using GridIndex=std::array<std::size_t,3>;
constexpr double kRoundoff=128*std::numeric_limits<double>::epsilon();

struct BuildFailure : std::runtime_error {
  Core::StatusCode code;
  BuildFailure(Core::StatusCode c,const std::string& message):
      std::runtime_error(message),code(c){}
};

void Require(bool condition,const std::string& message,
    Core::StatusCode code=Core::StatusCode::InvalidConfiguration) {
  if(!condition)throw BuildFailure(code,message);
}

Vec3 Position(CurvilinearPoint point,const SectorBasis& basis) {
  const double c=std::sqrt((1-point.mu)*(1+point.mu));
  return point.radiusM*(c*std::cos(point.phiRad)*basis.x+
      c*std::sin(point.phiRad)*basis.y+point.mu*basis.z);
}

bool ValidDigest(const std::string& value) {
  return value.size()==64&&std::all_of(value.begin(),value.end(),[](char c) {
    return (c>='0'&&c<='9')||(c>='a'&&c<='f');
  });
}

void CheckEdges(const std::vector<double>& values) {
  Require(values.size()>=2,"mesh axis requires at least two ordered edges");
  for(std::size_t i=0;i<values.size();++i) {
    Require(std::isfinite(values[i]),"mesh coordinate is nonfinite");
    if(i)Require(values[i]>values[i-1],"mesh coordinates must increase strictly");
  }
}

std::string Identity(const StaticMeshInput& input,const std::string& event) {
  std::ostringstream bytes;bytes.imbue(std::locale::classic());
  bytes<<"static-spherical-shared-mesh-v1\n"<<event<<'\n'<<std::hexfloat;
  for(const auto* axis:{&input.radiusEdgesM,&input.muEdges,&input.phiEdgesRad}) {
    bytes<<axis->size()<<':';for(double value:*axis)bytes<<value<<',';bytes<<'\n';
  }
  for(const auto vector:{input.basis.x,input.basis.y,input.basis.z})
    bytes<<vector.x<<','<<vector.y<<','<<vector.z<<'\n';
  bytes<<input.periodicPhi;
  return CoronalCME::ComputeContentChecksum(bytes.str());
}

} // namespace

Core::Result<std::shared_ptr<const StaticSphericalMesh>> StaticSphericalMesh::Create(
    StaticMeshInput input,std::string event,double epoch,std::uint64_t generation) {
  using Return=Core::Result<std::shared_ptr<const StaticSphericalMesh>>;
  try {
    Require(ValidDigest(event)&&std::isfinite(epoch)&&generation>0,
        "mesh requires an owning event SHA-256, finite epoch and positive generation");
    CheckEdges(input.radiusEdgesM);CheckEdges(input.muEdges);CheckEdges(input.phiEdgesRad);
    const std::size_t nr=input.radiusEdgesM.size()-1,nm=input.muEdges.size()-1,
        np=input.phiEdgesRad.size()-1;
    Require(input.maximumCells>0&&nr<=input.maximumCells/nm&&
        nr*nm<=input.maximumCells/np,"mesh exceeds its explicit cell budget");
    Require(!(nm==1&&input.muEdges.front()==-1&&input.muEdges.back()==1),
        "one cell cannot span both coordinate poles; angular subdivision required");
    if(input.periodicPhi)Require(np>=3&&
        std::abs((input.phiEdgesRad.back()-input.phiEdgesRad.front())-
            2*CoronalCME::Constants::kPi)<=kRoundoff*2*CoronalCME::Constants::kPi,
        "periodic longitude requires a full circle and at least three cells");
    std::shared_ptr<StaticSphericalMesh> out(new StaticSphericalMesh);
    out->input_=std::move(input);out->eventFingerprint_=std::move(event);
    out->epochS_=epoch;out->generation_=generation;
    const auto& in=out->input_;
    const auto cellIndex=[&](std::size_t r,std::size_t m,std::size_t p) {
      return (r*nm+m)*np+p;
    };
    out->cells_.reserve(nr*nm*np);
    for(std::size_t r=0;r<nr;++r)for(std::size_t m=0;m<nm;++m)
      for(std::size_t p=0;p<np;++p) {
        MeshCell cell;
        cell.bounds={in.radiusEdgesM[r],in.radiusEdgesM[r+1],
            in.muEdges[m],in.muEdges[m+1],in.phiEdgesRad[p],in.phiEdgesRad[p+1]};
        const auto measured=MeasureSphericalSector(cell.bounds,in.basis);
        if(!measured.ok())throw BuildFailure(measured.status.code,measured.status.message);
        cell.measure=measured.value;
        for(std::size_t f=0;f<6;++f)
          cell.faces[f].collapsedGeometricFace=cell.measure.faces[f].physicalAreaM2==0;
        out->cells_.push_back(std::move(cell));
      }
    std::map<GridIndex,std::size_t> nodeIds;
    std::map<std::pair<std::size_t,std::size_t>,std::size_t> edgeIds;
    const auto coordinate=[&](const GridIndex& grid) {
      return CurvilinearPoint{in.radiusEdgesM[grid[0]],in.muEdges[grid[1]],
          in.phiEdgesRad[grid[2]]};
    };
    const auto node=[&](GridIndex grid) {
      const auto physical=coordinate(grid);
      // Only exact coordinate identifications are collapsed. No distance
      // tolerance can merge distinct physical vertices near a pole or seam.
      if(in.radiusEdgesM[grid[0]]==0){grid[1]=0;grid[2]=0;}
      else if(std::abs(in.muEdges[grid[1]])==1)grid[2]=0;
      else if(in.periodicPhi&&grid[2]==np)grid[2]=0;
      const auto found=nodeIds.find(grid);
      if(found!=nodeIds.end())return found->second;
      const std::size_t id=out->nodes_.size();
      out->nodes_.push_back(Position(physical,in.basis));nodeIds.emplace(grid,id);
      return id;
    };
    const auto edge=[&](std::size_t first,std::size_t last,
        CurvilinearPoint a,CurvilinearPoint b) {
      Require(first!=last,"noncollapsed face contains a zero topological edge",
          Core::StatusCode::DataIntegrityFailure);
      const auto key=std::minmax(first,last);
      const std::pair<std::size_t,std::size_t> ownedKey{key.first,key.second};
      const int orientation=first<last?1:-1;
      const auto found=edgeIds.find(ownedKey);
      if(found!=edgeIds.end())return OrientedMeshEdge{found->second,orientation};
      MeshEdge value;value.firstNode=key.first;value.lastNode=key.second;
      if(a.radiusM!=b.radiusM) {
        value.axis=MeshAxis::Radial;
        if(a.radiusM==0){a.mu=b.mu;a.phiRad=b.phiRad;}
        else {b.mu=a.mu;b.phiRad=a.phiRad;}
        value.arcLengthM=std::abs(b.radiusM-a.radiusM);
      }else if(a.mu!=b.mu) {
        value.axis=MeshAxis::Mu;
        // After a polar corner collapses, its stored azimuth is arbitrary.
        // The remaining nonpolar endpoint identifies the actual meridian;
        // interpolating two different pole azimuths would invent a spiral.
        if(std::abs(a.mu)==1)a.phiRad=b.phiRad;
        else b.phiRad=a.phiRad;
        value.arcLengthM=a.radiusM*std::abs(std::asin(b.mu)-std::asin(a.mu));
      }else {
        value.axis=MeshAxis::Phi;
        value.arcLengthM=a.radiusM*std::sqrt((1-a.mu)*(1+a.mu))*
            std::abs(b.phiRad-a.phiRad);
      }
      Require(std::isfinite(value.arcLengthM)&&value.arcLengthM>0,
          "curved edge length is unresolved",Core::StatusCode::NumericalFailure);
      if(orientation>0){value.start=a;value.end=b;}
      else {value.start=b;value.end=a;}
      const double scale=std::max(a.radiusM,b.radiusM);
      Require(Norm(Position(value.start,in.basis)-out->nodes_[value.firstNode])<=
          kRoundoff*scale&&
          Norm(Position(value.end,in.basis)-out->nodes_[value.lastNode])<=kRoundoff*scale,
          "edge curve does not meet its canonical nodes",Core::StatusCode::DataIntegrityFailure);
      const std::size_t id=out->edges_.size();out->edges_.push_back(value);
      edgeIds.emplace(ownedKey,id);return OrientedMeshEdge{id,orientation};
    };
    const auto face=[&](MeshAxis axis,std::size_t left,std::size_t right,
        std::size_t leftSide,std::size_t rightSide,
        const std::array<GridIndex,4>& corners) {
      Require(left!=NoMeshIndex||right!=NoMeshIndex,"face has no owning cell");
      const auto& raw=left!=NoMeshIndex?out->cells_[left].measure.faces[leftSide]:
          out->cells_[right].measure.faces[rightSide];
      if(raw.physicalAreaM2==0)return;
      if(left!=NoMeshIndex&&right!=NoMeshIndex) {
        const auto& other=out->cells_[right].measure.faces[rightSide];
        const double scale=raw.physicalAreaM2+other.physicalAreaM2;
        Require(std::abs(raw.physicalAreaM2-other.physicalAreaM2)<=kRoundoff*scale&&
            Norm(raw.outwardAreaVectorM2+other.outwardAreaVectorM2)<=kRoundoff*scale,
            "neighbors disagree on a shared face measure",Core::StatusCode::DataIntegrityFailure);
      }
      MeshFace value;value.axis=axis;value.leftCell=left;value.rightCell=right;
      value.measure=raw;
      if(left==NoMeshIndex)value.measure.outwardAreaVectorM2=
          -1.0*value.measure.outwardAreaVectorM2;
      std::vector<CurvilinearPoint> coordinates;
      for(const auto& corner:corners) {
        const auto id=node(corner);
        if(value.positiveNodeLoop.empty()||value.positiveNodeLoop.back()!=id) {
          value.positiveNodeLoop.push_back(id);coordinates.push_back(coordinate(corner));
        }
      }
      if(value.positiveNodeLoop.front()==value.positiveNodeLoop.back()) {
        value.positiveNodeLoop.pop_back();coordinates.pop_back();
      }
      Require(value.positiveNodeLoop.size()>=3,
          "finite curved face needs further angular subdivision",
          Core::StatusCode::UnsupportedCapability);
      for(std::size_t k=0;k<value.positiveNodeLoop.size();++k) {
        const std::size_t next=(k+1)%value.positiveNodeLoop.size();
        value.positiveEdgeLoop.push_back(edge(value.positiveNodeLoop[k],
            value.positiveNodeLoop[next],coordinates[k],coordinates[next]));
      }
      const std::size_t id=out->faces_.size();out->faces_.push_back(std::move(value));
      if(left!=NoMeshIndex)out->cells_[left].faces[leftSide]={id,1,false};
      if(right!=NoMeshIndex)out->cells_[right].faces[rightSide]={id,-1,false};
    };
    // The parameter tuple (r,mu,phi) is left-handed. The explicit loops below
    // therefore use (phi,mu), (r,phi), (mu,r) to point in +r,+mu,+phi.
    for(std::size_t r=0;r<=nr;++r)for(std::size_t m=0;m<nm;++m)
      for(std::size_t p=0;p<np;++p)
        face(MeshAxis::Radial,r?cellIndex(r-1,m,p):NoMeshIndex,
            r<nr?cellIndex(r,m,p):NoMeshIndex,1,0,
            {{{r,m,p},{r,m,p+1},{r,m+1,p+1},{r,m+1,p}}});
    for(std::size_t r=0;r<nr;++r)for(std::size_t m=0;m<=nm;++m)
      for(std::size_t p=0;p<np;++p)
        face(MeshAxis::Mu,m?cellIndex(r,m-1,p):NoMeshIndex,
            m<nm?cellIndex(r,m,p):NoMeshIndex,3,2,
            {{{r,m,p},{r+1,m,p},{r+1,m,p+1},{r,m,p+1}}});
    const std::size_t phiPlanes=in.periodicPhi?np:np+1;
    for(std::size_t r=0;r<nr;++r)for(std::size_t m=0;m<nm;++m)
      for(std::size_t p=0;p<phiPlanes;++p) {
        const std::size_t left=p?cellIndex(r,m,p-1):
            (in.periodicPhi?cellIndex(r,m,np-1):NoMeshIndex);
        face(MeshAxis::Phi,left,p<np?cellIndex(r,m,p):NoMeshIndex,5,4,
            {{{r,m,p},{r,m+1,p},{r+1,m+1,p},{r+1,m,p}}});
      }
    for(auto& cell:out->cells_) {
      std::map<std::size_t,int> boundaryOfBoundary;
      double area=0;
      for(const auto& use:cell.faces) {
        if(use.collapsedGeometricFace) {
          Require(use.index==NoMeshIndex,"collapsed face gained finite flux ownership");
          continue;
        }
        Require(use.index<out->faces_.size()&&std::abs(use.outwardOrientation)==1,
            "mesh cell has an unowned physical face",Core::StatusCode::DataIntegrityFailure);
        const auto& value=out->faces_[use.index];
        area+=value.measure.physicalAreaM2;
        cell.sharedAreaClosureResidualM2=cell.sharedAreaClosureResidualM2+
            use.outwardOrientation*value.measure.outwardAreaVectorM2;
        for(const auto& e:value.positiveEdgeLoop)
          boundaryOfBoundary[e.index]+=use.outwardOrientation*e.orientation;
      }
      for(const auto& coefficient:boundaryOfBoundary)Require(coefficient.second==0,
          "cell boundary-of-boundary is nonzero; CT topology is inconsistent",
          Core::StatusCode::DataIntegrityFailure);
      Require(Norm(cell.sharedAreaClosureResidualM2)<=kRoundoff*area,
          "shared-face measures fail static geometric closure",Core::StatusCode::NumericalFailure);
    }
    out->fingerprint_=Identity(in,out->eventFingerprint_);
    return Return::Success(std::move(out));
  }catch(const BuildFailure& error) {
    return Return::Failure(error.code,error.what());
  }catch(const std::bad_alloc&) {
    return Return::Failure(Core::StatusCode::NumericalFailure,"static mesh allocation failed");
  }catch(const std::length_error&) {
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "static mesh capacity is not representable on this host");
  }
}

Core::Result<EdgeSample> StaticSphericalMesh::SampleEdge(std::size_t index,double s) const {
  using Return=Core::Result<EdgeSample>;
  if(index>=edges_.size()||!std::isfinite(s)||s<0||s>1)return Return::Failure(
      Core::StatusCode::OutOfDomain,"edge index/parameter is outside static mesh support");
  const auto& edge=edges_[index];const auto& b=input_.basis;
  const auto& a=edge.start;const auto& z=edge.end;
  EdgeSample out;
  if(edge.axis==MeshAxis::Radial) {
    CurvilinearPoint p=a;p.radiusM=a.radiusM+s*(z.radiusM-a.radiusM);
    out.positionHciM=Position(p,b);
    p.radiusM=1;out.derivativeM=(z.radiusM-a.radiusM)*Position(p,b);
  }else if(edge.axis==MeshAxis::Mu) {
    const double start=std::asin(a.mu),delta=std::asin(z.mu)-start;
    const double latitude=start+s*delta;
    const Vec3 equator=std::cos(a.phiRad)*b.x+std::sin(a.phiRad)*b.y;
    out.positionHciM=a.radiusM*(std::cos(latitude)*equator+std::sin(latitude)*b.z);
    out.derivativeM=a.radiusM*delta*(-std::sin(latitude)*equator+std::cos(latitude)*b.z);
  }else {
    const double delta=z.phiRad-a.phiRad,phi=a.phiRad+s*delta;
    const double c=std::sqrt((1-a.mu)*(1+a.mu));
    out.positionHciM=a.radiusM*(c*std::cos(phi)*b.x+c*std::sin(phi)*b.y+a.mu*b.z);
    out.derivativeM=a.radiusM*c*delta*(-std::sin(phi)*b.x+std::cos(phi)*b.y);
  }
  return Return::Success(out);
}

} } } // namespace SEP::CoronaSwcme::CrossRay
