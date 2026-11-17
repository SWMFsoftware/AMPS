#include "sep_cli_campaign_validation.h"
#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)

#include <mpi.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

// The application may use a communicator supplied by its enclosing AMPS
// driver. Do not substitute MPI_COMM_WORLD, or rank-local fault provenance
// would refer to a different process group in a coupled embedding.
extern MPI_Comm MPI_GLOBAL_COMMUNICATOR;

namespace SEP { namespace Testing { namespace {

struct FaultInput {
  bool enabled=false;
  std::uint32_t rank=0,replicate=0;
  std::string exactText;
};

std::uint32_t Unsigned(const std::string& text) {
  if(text.empty())throw std::runtime_error("CLIEX01 fault input has an empty unsigned value");
  std::uint64_t value=0;
  for(char c:text) {
    if(c<'0'||c>'9')throw std::runtime_error("CLIEX01 fault input requires unsigned decimal values");
    const unsigned digit=static_cast<unsigned>(c-'0');
    if(value>(std::numeric_limits<std::uint32_t>::max()-digit)/10)
      throw std::runtime_error("CLIEX01 fault input integer overflows uint32");
    value=value*10+digit;
  }
  return static_cast<std::uint32_t>(value);
}

FaultInput ReadFaultInput(const ExecutionContext& context,int ranks) {
  FaultInput fault;
  if(context.inputPath.empty())return fault;
  std::ifstream input(context.inputPath.c_str(),std::ios::binary);
  if(!input)throw std::runtime_error("CLIEX01 cannot read explicit fault-test input");
  std::ostringstream exact;exact<<input.rdbuf();
  if(input.bad())throw std::runtime_error("CLIEX01 fault-test input read failed");
  fault.exactText=exact.str();std::istringstream lines(fault.exactText);
  std::string line;
  if(!std::getline(lines,line)||line!="srcsep-cli-fault-v1")
    throw std::runtime_error("CLIEX01 requires srcsep-cli-fault-v1 input schema");
  std::set<std::string> keys;
  while(std::getline(lines,line)) {
    const auto equal=line.find('=');
    if(equal==std::string::npos||line.find('=',equal+1)!=std::string::npos)
      throw std::runtime_error("CLIEX01 fault input requires exactly key=value lines");
    const auto key=line.substr(0,equal);
    if(!keys.insert(key).second)throw std::runtime_error("CLIEX01 duplicate fault-test input key");
    const auto value=Unsigned(line.substr(equal+1));
    if(key=="fault_rank")fault.rank=value;
    else if(key=="fault_replicate")fault.replicate=value;
    else throw std::runtime_error("CLIEX01 unknown fault-test input key: "+key);
  }
  if(keys.size()!=2||!keys.count("fault_rank")||!keys.count("fault_replicate"))
    throw std::runtime_error("CLIEX01 fault-test input must declare rank and replicate");
  if(fault.rank>=static_cast<std::uint32_t>(ranks)||fault.replicate>=context.replicateCount)
    throw std::runtime_error("CLIEX01 fault-test input refers to a rank or replicate outside this campaign");
  fault.enabled=true;return fault;
}

void Commit(const std::string& path,const std::string& data) {
  // Each invocation owns its directory. The temporary name is never shared by
  // ranks or replicates; flush/check/rename makes a report reference complete.
  const auto temporary=path+".tmp";
  std::ofstream output(temporary.c_str(),std::ios::binary|std::ios::trunc);
  output<<data;output.flush();
  if(!output.good())throw std::runtime_error("CLIEX01 cannot flush artifact: "+path);
  output.close();
  if(output.fail()||std::rename(temporary.c_str(),path.c_str())!=0)
    throw std::runtime_error("CLIEX01 cannot commit artifact: "+path);
}

Result Probe() {
  const auto& context=GetExecutionContext();
  int rank=0,ranks=1;
  MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&rank);
  MPI_Comm_size(MPI_GLOBAL_COMMUNICATOR,&ranks);
  const auto fault=ReadFaultInput(context,ranks);
  if(context.artifactDirectory.empty())
    throw std::runtime_error("CLIEX01 requires the native campaign artifact directory");
  const auto seed=TestSeed(177,UINT64_C(0x434c4945583031));
  Transport::KeyedRandomStream random(seed,37,UINT64_C(0x434c4950524f4245),0);
  const double uniform=random.UniformOpen01(),normal=random.Normal01();
  const auto state=random.SaveState();
  std::ostringstream draw;
  draw<<std::setprecision(17)<<"{\"schema\":\"srcsep-cli-probe-v1\",\"rank\":"<<rank
      <<",\"replicate\":"<<context.replicateIndex<<",\"seed\":"<<seed
      <<",\"uniform\":"<<uniform<<",\"normal\":"<<normal
      <<",\"integer_state\":"<<state.integerState<<"}\n";
  const auto artifact=context.artifactDirectory+"/CLIEX01_draw.json";
  Commit(artifact,draw.str());
  // This deliberately undeclared scratch file tests --keep-test-artifacts.
  // It is neither a physics output nor a permanent declared artifact.
  Commit(".cli-probe-unlisted.tmp","CLIEX01 owned scratch marker\n");
  Result result;
  const bool requested=fault.enabled&&fault.rank==static_cast<std::uint32_t>(rank)&&
      fault.replicate==context.replicateIndex;
  const bool valid=uniform>0.0&&uniform<1.0&&std::isfinite(normal);
  result.status=requested||!valid?Status::Fail:Status::Pass;
  result.message=requested?"explicit infrastructure fault requested on this rank and replicate":
      valid?"actual keyed draws and transactional rank-local artifact completed":
      "keyed generator returned an invalid actual uniform/normal draw";
  result.hasSeed=true;result.seed=seed;result.artifacts.push_back(artifact);
  // A relative reference outside owned scratch must be retained in place,
  // never moved by the harness as though it were a temporary output.
  Commit(context.artifactDirectory+"/CLIEX01_immutable_reference.txt","CLIEX01 immutable reference sentinel\n");
  result.artifacts.push_back("../CLIEX01_immutable_reference.txt");
  result.configuration.push_back("evidence_level=infrastructure-probe;not-physics-validation");
  result.configuration.push_back("fault_source="+(fault.enabled?context.inputPath:std::string("none")));
  result.configuration.push_back("rank="+std::to_string(rank));
  result.configuration.push_back("rng_integer_state="+std::to_string(state.integerState));
  if(fault.enabled) {
    const auto provenance=context.artifactDirectory+"/CLIEX01_fault_input.txt";
    Commit(provenance,fault.exactText);result.artifacts.push_back(provenance);
    result.configuration.push_back("fault_rank="+std::to_string(fault.rank));
    result.configuration.push_back("fault_replicate="+std::to_string(fault.replicate));
  }
  result.metrics.push_back({"actual_uniform_draw",uniform,0.0,">","dimensionless"});
  result.metrics.push_back({"actual_normal_draw",normal,0.0,"observed","dimensionless"});
  result.metrics.push_back({"explicit_fault_requested",requested?1.0:0.0,0.0,"==","count"});
  result.metrics.push_back({"assertion_failures",valid?0.0:1.0,0.0,"==","count"});
  std::cout<<"CLIEX01 rank="<<rank<<" replicate="<<context.replicateIndex<<" committed actual random draws\n";
  return result;
}

} // namespace

std::vector<Descriptor> CliCampaignDescriptors() {
  Descriptor descriptor;
  descriptor.id="CLIEX01";descriptor.name="Native CLI campaign infrastructure probe";
  descriptor.group="cli-infrastructure";
  descriptor.description="Exercise applied keyed seeds, per-rank artifacts, repeats and explicitly requested rank-local fault reduction; no physics qualification.";
  descriptor.initialization=InitializationLevel::None;
  descriptor.supportedBuildModes="native AMPS MPI campaign; no mesh or particle initialization";
  descriptor.runtime=RuntimeClass::Routine;
  descriptor.seedPolicy="TestSeed schema1; default177 with explicit probe purpose; actual uniform and normal draws";
  descriptor.seedCapability=SeedCapability::Context;descriptor.repeatSupported=true;
  descriptor.stateIsolation="private keyed generator; harness-owned rank/replicate scratch; optional exact fault input archived";
  descriptor.tags.push_back("infrastructure");descriptor.callback=Probe;
  return {descriptor};
}

} }
