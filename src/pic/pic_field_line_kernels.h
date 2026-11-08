#ifndef _AMPS_PIC_FIELD_LINE_KERNELS_H_
#define _AMPS_PIC_FIELD_LINE_KERNELS_H_

/* Dependency-light geometry and conservation kernels for field-line data. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace PIC { namespace FieldLineKernels {
enum class Status { Success,InvalidArgument,OutsideDomain,DuplicateIndex };

inline Status LinearInterpolate(const double* begin,const double* end,
                                std::size_t length,double s,double* output) {
  if (!begin || !end || !output || length==0 || !std::isfinite(s))
    return Status::InvalidArgument;
  if (s<0.0 || s>1.0) return Status::OutsideDomain;
  for (std::size_t i=0;i<length;++i) {
    if (!std::isfinite(begin[i]) || !std::isfinite(end[i]))
      return Status::InvalidArgument;
    output[i]=(1.0-s)*begin[i]+s*end[i];
  }
  return Status::Success;
}

inline Status InterpolatePolylineByArc(const std::vector<double>& xyz,
                                       double arc,double output[3]) {
  if (!output || xyz.size()<6 || xyz.size()%3!=0 || !std::isfinite(arc))
    return Status::InvalidArgument;
  if (arc<0.0) return Status::OutsideDomain;
  double passed=0.0;
  for (std::size_t v=0;v+3<xyz.size();v+=3) {
    const double dx=xyz[v+3]-xyz[v],dy=xyz[v+4]-xyz[v+1],dz=xyz[v+5]-xyz[v+2];
    const double length=std::sqrt(dx*dx+dy*dy+dz*dz);
    if (!std::isfinite(length) || length<=0.0) return Status::InvalidArgument;
    if (arc<=passed+length || v+6==xyz.size()) {
      if (arc>passed+length) return Status::OutsideDomain;
      const double s=(arc-passed)/length;
      return LinearInterpolate(&xyz[v],&xyz[v+3],3,s,output);
    }
    passed+=length;
  }
  return Status::OutsideDomain;
}

inline Status ConservativeRepartition1D(const std::vector<double>& old_edges,
                                        const std::vector<double>& old_total,
                                        const std::vector<double>& new_edges,
                                        std::vector<double>& new_total) {
  if (old_edges.size()!=old_total.size()+1 || old_total.empty() ||
      new_edges.size()<2 || old_edges.front()!=new_edges.front() ||
      old_edges.back()!=new_edges.back()) return Status::InvalidArgument;
  for (std::size_t i=1;i<old_edges.size();++i)
    if (!std::isfinite(old_edges[i]) || old_edges[i]<=old_edges[i-1] ||
        !std::isfinite(old_total[i-1])) return Status::InvalidArgument;
  for (std::size_t i=1;i<new_edges.size();++i)
    if (!std::isfinite(new_edges[i]) || new_edges[i]<=new_edges[i-1])
      return Status::InvalidArgument;
  new_total.assign(new_edges.size()-1,0.0);
  for (std::size_t n=0;n<new_total.size();++n)
    for (std::size_t o=0;o<old_total.size();++o) {
      const double left=std::max(new_edges[n],old_edges[o]);
      const double right=std::min(new_edges[n+1],old_edges[o+1]);
      if (right>left)
        new_total[n]+=old_total[o]*(right-left)/(old_edges[o+1]-old_edges[o]);
    }
  return Status::Success;
}

struct WeightedParticle { double s,v_parallel,v_normal,weight; };
inline Status SplitParticle(const WeightedParticle& parent,std::size_t count,
                            std::vector<WeightedParticle>& children) {
  if (count<2 || !std::isfinite(parent.s) || !std::isfinite(parent.v_parallel) ||
      !std::isfinite(parent.v_normal) || !std::isfinite(parent.weight) ||
      parent.weight<=0.0) return Status::InvalidArgument;
  children.assign(count,parent);
  for (auto& child:children) child.weight=parent.weight/count;
  return Status::Success;
}

inline Status ExchangeIndexedValues(const std::vector<std::size_t>& index,
                                    const std::vector<double>& packed,
                                    std::size_t item_count,std::size_t width,
                                    std::vector<double>& unpacked) {
  if (width==0 || index.size()!=item_count || packed.size()!=item_count*width)
    return Status::InvalidArgument;
  unpacked.assign(item_count*width,std::numeric_limits<double>::quiet_NaN());
  std::vector<bool> seen(item_count,false);
  for (std::size_t item=0;item<item_count;++item) {
    if (index[item]>=item_count) return Status::OutsideDomain;
    if (seen[index[item]]) return Status::DuplicateIndex;
    seen[index[item]]=true;
    for (std::size_t component=0;component<width;++component)
      unpacked[index[item]*width+component]=packed[item*width+component];
  }
  return Status::Success;
}
}} // namespace PIC::FieldLineKernels
#endif
