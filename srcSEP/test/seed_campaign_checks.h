#ifndef SRCSEP_TEST_SEED_CAMPAIGN_CHECKS_H
#define SRCSEP_TEST_SEED_CAMPAIGN_CHECKS_H

#include "sep_test_registry.h"
#include <algorithm>
#include <map>
#include <sstream>
#include <stdexcept>

namespace SEP { namespace Testing {

// An external check of the callback's actual stochastic observable. Merely
// changing result.seed would not satisfy these assertions. Execution order is
// deliberately reversed for the second campaign rather than being hidden by
// Select's canonical sorting. Statistical physics gates still run normally.
inline bool CheckAppliedCampaignSeeds(const Registry& registry,
    std::vector<const Descriptor*> selected,std::ostream& out) {
  CampaignOptions options;options.seedProvided=true;options.campaignSeed=123;
  options.repeatCount=2;std::ostringstream quiet;
  const auto first=registry.RunCampaign(selected,options,quiet);
  std::reverse(selected.begin(),selected.end());
  const auto reordered=registry.RunCampaign(selected,options,quiet);
  options.campaignSeed=456;
  const auto different=registry.RunCampaign(selected,options,quiet);
  bool ok=first.ExitCode()==0&&reordered.ExitCode()==0&&different.ExitCode()==0&&
      first.results.size()==selected.size()&&reordered.results.size()==selected.size()&&
      different.results.size()==selected.size();
  using Pair=std::pair<std::uint64_t,double>;
  auto values=[](const Summary& summary) {
    std::map<std::string,std::vector<Pair>> result;
    for(const auto& aggregate:summary.results) {
      for(const auto& invocation:aggregate.invocations) {
        if(!invocation.result||!invocation.result->hasSeed||
           invocation.result->status!=Status::Pass)
          throw std::runtime_error("seeded invocation lacks successful applied seed");
        const auto& child=*invocation.result;bool found=false;
        for(const auto& metric:child.metrics)if(metric.name=="actual_seeded_observable") {
          result[aggregate.id].push_back({child.seed,metric.value});found=true;
        }
        if(!found)throw std::runtime_error("seed test requires a reported actual stochastic observable");
      }
    }
    return result;
  };
  try {
    const auto a=values(first),b=values(reordered),c=values(different);
    ok=ok&&a==b;
    for(const auto& entry:a) {
      const auto& x=entry.second;const auto& y=c.at(entry.first);
      ok=ok&&x.size()==2&&y.size()==2;
      if(x.size()==2&&y.size()==2) {
        ok=ok&&x[0].first!=x[1].first&&x[0].second!=x[1].second;
        for(unsigned replicate=0;replicate<2;++replicate)
          ok=ok&&x[replicate].first!=y[replicate].first&&
              x[replicate].second!=y[replicate].second;
      }
    }
    std::string error;
    ok=ok&&WriteJsonSummary(first,"seeded-campaign-A.json",&error)&&
        WriteJsonSummary(reordered,"seeded-campaign-reordered.json",&error)&&
        WriteJsonSummary(different,"seeded-campaign-B.json",&error)&&
        WriteJUnitSummary(first,"seeded-campaign-A.xml",&error);
    if(!error.empty())out<<error<<'\n';
  } catch(const std::exception& exception) {
    out<<"seeded campaign check: "<<exception.what()<<'\n';ok=false;
  }
  out<<(ok?"PASS":"FAIL")<<" APPLIED-SEEDS: actual observables change with campaign/replicate keys and remain identical under reversed execution order\n";
  return ok;
}

} }
#endif
