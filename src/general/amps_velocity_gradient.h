#ifndef AMPS_VELOCITY_GRADIENT_H
#define AMPS_VELOCITY_GRADIENT_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace AMPS { namespace Coupling {

// Spatial reconstruction only: x is Cartesian m, U is Cartesian m/s, and all
// samples must describe the same epoch/frame. No density continuity, radial
// symmetry, field-line tangent derivative or temporal extrapolation supplies
// a missing transverse derivative. The caller owns epoch/frame validation.
class VelocityGradient {
 public:
  struct Sample {
    std::array<double,3> x,u;
    std::uint64_t line;
  };
  struct Result {
    bool valid=false;
    std::array<double,9> gradient{{0,0,0,0,0,0,0,0,0}}; // dU_i/dx_j [1/s]
    double condition=0.0;
    std::size_t neighbors=0;
    std::string message;
  };
  // These are numerical stencil/roundoff controls, not physical-domain limits.
  // With cond(R)<=1e6, roundoff amplification is bounded independently of the
  // spatial truncation error. Four distinct producer lines require genuinely
  // transverse data; curvature of one densely sampled line is not a 3-D grid.
  enum {MinimumLines=4,InitialNeighbors=32,MaximumNeighbors=512};
  explicit VelocityGradient(const std::vector<Sample>& samples):samples_(samples) {
    if(samples_.size()<4)throw std::invalid_argument("velocity gradient requires a three-dimensional spatial stencil");
    std::vector<std::size_t> order(samples_.size());
    for(std::size_t i=0;i<samples_.size();++i) {
      order[i]=i;
      for(int d=0;d<3;++d)if(!std::isfinite(samples_[i].x[d])||!std::isfinite(samples_[i].u[d]))
        throw std::invalid_argument("velocity gradient received nonfinite spatial samples");
    }
    // Balanced median kd-tree bounds nearest-neighbor lookup without an O(N^2)
    // all-vertices scan. Build once between imports; const queries are read-only.
    nodes_.reserve(samples_.size());root_=Build(order,0,order.size(),0);
  }
  Result Evaluate(std::size_t index) const {
    Result result;
    if(index>=samples_.size()){result.message="velocity gradient query index is invalid";return result;}
    const std::size_t maximum=std::min<std::size_t>(MaximumNeighbors,samples_.size());
    for(std::size_t count=std::min<std::size_t>(InitialNeighbors,maximum);;) {
      Heap heap;Nearest(root_,samples_[index].x,count,&heap);
      std::vector<std::size_t> stencil;
      while(!heap.empty()){stencil.push_back(heap.top().second);heap.pop();}
      std::reverse(stencil.begin(),stencil.end());
      result=Fit(index,stencil);
      if(result.valid||count==maximum)return result;
      count=std::min(maximum,count*2);
    }
  }
 private:
  using Row=std::array<long double,3>;
  using Heap=std::priority_queue<std::pair<double,std::size_t> >;
  struct Node {std::size_t sample;int left,right,axis;};
  std::vector<Sample> samples_;
  std::vector<Node> nodes_;
  int root_=-1;
  int Build(std::vector<std::size_t>& order,std::size_t begin,std::size_t end,int depth) {
    if(begin==end)return -1;
    const int axis=depth%3;const std::size_t middle=begin+(end-begin)/2;
    std::nth_element(order.begin()+begin,order.begin()+middle,order.begin()+end,
      [&](std::size_t a,std::size_t b){return samples_[a].x[axis]<samples_[b].x[axis]||
        (samples_[a].x[axis]==samples_[b].x[axis]&&a<b);});
    const int node=static_cast<int>(nodes_.size());nodes_.push_back({order[middle],-1,-1,axis});
    const int left=Build(order,begin,middle,depth+1),right=Build(order,middle+1,end,depth+1);
    nodes_[node].left=left;nodes_[node].right=right;return node;
  }
  void Nearest(int node,const std::array<double,3>& x,std::size_t count,Heap* heap) const {
    if(node<0)return;
    const auto& n=nodes_[node];const auto& point=samples_[n.sample].x;
    double distance=0.0;for(int d=0;d<3;++d)distance+=(point[d]-x[d])*(point[d]-x[d]);
    const std::pair<double,std::size_t> candidate(distance,n.sample);
    if(heap->size()<count)heap->push(candidate);
    else if(candidate<heap->top()){heap->pop();heap->push(candidate);}
    const double delta=x[n.axis]-point[n.axis];
    Nearest(delta<0?n.left:n.right,x,count,heap);
    if(heap->size()<count||delta*delta<=heap->top().first)
      Nearest(delta<0?n.right:n.left,x,count,heap);
  }
  Result Fit(std::size_t index,const std::vector<std::size_t>& stencil) const {
    Result result;result.neighbors=stencil.size();
    std::set<std::uint64_t> lines;std::vector<Row> a,rhs;
    const auto& center=samples_[index];
    for(auto neighbor:stencil) {
      const auto& s=samples_[neighbor];Row dx,du;
      long double r2=0.0L;
      for(int d=0;d<3;++d){dx[d]=static_cast<long double>(s.x[d])-center.x[d];
        du[d]=static_cast<long double>(s.u[d])-center.u[d];r2+=dx[d]*dx[d];}
      if(r2==0.0L) {
        if(du[0]!=0||du[1]!=0||du[2]!=0){result.message="coincident velocity samples disagree";return result;}
        continue;
      }
      // Solve min sum |(delta U-G*delta x)/|delta x||^2. Centering removes
      // bulk-flow offsets and inverse-distance weighting is explicit. Rows
      // have no coordinate units; RHS already has derivative units [1/s].
      const long double radius=std::sqrt(r2);
      for(int d=0;d<3;++d){dx[d]/=radius;du[d]/=radius;}
      a.push_back(dx);rhs.push_back(du);lines.insert(s.line);
    }
    if(lines.size()<MinimumLines||a.size()<3){result.message="velocity stencil lacks transverse producer lines";return result;}
    int permutation[3]={0,1,2};
    // Column-pivoted Householder QR avoids squaring the stencil condition
    // number through normal equations. Apply the same orthogonal transforms
    // to all three velocity components, preserving full Cartesian information.
    for(int k=0;k<3;++k) {
      int pivot=k;long double largest=-1.0L;
      for(int j=k;j<3;++j){long double norm=0;for(std::size_t i=k;i<a.size();++i)norm+=a[i][j]*a[i][j];
        if(norm>largest){largest=norm;pivot=j;}}
      if(!(largest>0)){result.message="velocity stencil is rank deficient";return result;}
      if(pivot!=k){for(auto& row:a)std::swap(row[k],row[pivot]);std::swap(permutation[k],permutation[pivot]);}
      const long double alpha=-std::copysign(std::sqrt(largest),a[k][k]);
      std::vector<long double> v(a.size()-k);long double squared=0;
      for(std::size_t i=k;i<a.size();++i){v[i-k]=a[i][k]-(i==static_cast<std::size_t>(k)?alpha:0);squared+=v[i-k]*v[i-k];}
      if(!(squared>0)){result.message="velocity QR reflector is unresolved";return result;}
      for(int j=k;j<3;++j){long double dot=0;for(std::size_t i=k;i<a.size();++i)dot+=v[i-k]*a[i][j];
        for(std::size_t i=k;i<a.size();++i)a[i][j]-=2*dot*v[i-k]/squared;}
      for(int component=0;component<3;++component){long double dot=0;for(std::size_t i=k;i<a.size();++i)dot+=v[i-k]*rhs[i][component];
        for(std::size_t i=k;i<a.size();++i)rhs[i][component]-=2*dot*v[i-k]/squared;}
      a[k][k]=alpha;
    }
    long double inverse[3][3]={{0}},norm=0,inverseNorm=0;
    for(int column=0;column<3;++column)for(int i=2;i>=0;--i) {
      long double value=i==column?1.0L:0.0L;
      for(int j=i+1;j<3;++j)value-=a[i][j]*inverse[j][column];
      if(a[i][i]==0){result.message="velocity stencil has a singular triangular factor";return result;}
      inverse[i][column]=value/a[i][i];
    }
    for(int j=0;j<3;++j){long double sum=0,isum=0;for(int i=0;i<3;++i){if(i<=j)sum+=std::abs(a[i][j]);isum+=std::abs(inverse[i][j]);}
      norm=std::max(norm,sum);inverseNorm=std::max(inverseNorm,isum);}
    result.condition=static_cast<double>(norm*inverseNorm);
    if(!std::isfinite(result.condition)||result.condition>1.0e6){result.message="velocity stencil is numerically ill-conditioned";return result;}
    for(int component=0;component<3;++component) {
      long double value[3]={0};
      for(int i=2;i>=0;--i){value[i]=rhs[i][component];for(int j=i+1;j<3;++j)value[i]-=a[i][j]*value[j];value[i]/=a[i][i];}
      for(int d=0;d<3;++d){const double entry=static_cast<double>(value[d]);
        if(!std::isfinite(entry)){result.message="velocity derivative is nonfinite";return result;}
        result.gradient[3*component+permutation[d]]=entry;}
    }
    result.valid=true;return result;
  }
};
} }
#endif
