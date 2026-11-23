#include "sep_test_registry.h"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <limits>
#include <cstdio>
#include <unistd.h>

namespace SEP {
namespace Testing {

namespace {

// The standalone executable installs this record once, after command-line
// validation and before registry dispatch. Component callbacks receive a const
// view, preventing one case from redirecting a later case's input or evidence.
ExecutionContext executionContext;
struct SeedUse {std::uint64_t original,purpose,actual;};
std::vector<SeedUse> seedUses;
void HashWord(std::uint64_t value,std::uint64_t& hash) {
  for(unsigned byte=0;byte<8;++byte) {
    hash^=static_cast<unsigned char>(value>>(8*byte));
    hash*=UINT64_C(1099511628211);
  }
}

}  // namespace

void SetExecutionContext(const ExecutionContext& context) {
  executionContext = context;
}

const ExecutionContext& GetExecutionContext() {
  return executionContext;
}
std::uint64_t DeriveTestCaseSeed(std::uint64_t campaign,const std::string& id,
                                std::uint32_t replicate) {
  std::uint64_t hash=UINT64_C(14695981039346656037);
  HashWord(UINT64_C(0x5445535443415345),hash);HashWord(TestSeedSchema,hash);
  HashWord(campaign,hash);
  for(unsigned char value:id) {hash^=static_cast<unsigned char>(std::toupper(value));hash*=UINT64_C(1099511628211);}
  HashWord(replicate,hash);return hash;
}
std::uint64_t TestSeed(std::uint64_t original,std::uint64_t purpose) {
  if(!executionContext.seedActive)return original;
  for(const auto& use:seedUses)
    if(use.original==original&&use.purpose==purpose)return use.actual;
  std::uint64_t hash=UINT64_C(14695981039346656037);
  HashWord(UINT64_C(0x5445535453454544),hash);HashWord(TestSeedSchema,hash);
  HashWord(executionContext.caseSeed,hash);HashWord(original,hash);HashWord(purpose,hash);
  SeedUse use;use.original=original;use.purpose=purpose;use.actual=hash;seedUses.push_back(use);
  return hash;
}
namespace {

std::string FoldCase(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) {
                   return static_cast<char>(std::tolower(c));
                 });
  return value;
}

bool DescriptorIdLess(const Descriptor& left, const Descriptor& right) {
  return left.id < right.id;
}

bool DescriptorPointerIdLess(const Descriptor* left, const Descriptor* right) {
  return left->id < right->id;
}

std::string JsonEscape(const std::string& value) {
  std::ostringstream out;
  for (std::size_t i = 0; i < value.size(); ++i) {
    const unsigned char c = static_cast<unsigned char>(value[i]);
    switch (c) {
      case '\\': out << "\\\\"; break;
      case '"': out << "\\\""; break;
      case '\b': out << "\\b"; break;
      case '\f': out << "\\f"; break;
      case '\n': out << "\\n"; break;
      case '\r': out << "\\r"; break;
      case '\t': out << "\\t"; break;
      default:
        if (c < 0x20U) {
          out << "\\u" << std::hex << std::setw(4) << std::setfill('0')
              << static_cast<unsigned>(c) << std::dec << std::setfill(' ');
        }
        else out << static_cast<char>(c);
    }
  }
  return out.str();
}

std::string XmlEscape(const std::string& value) {
  std::string result;
  for (std::size_t i = 0; i < value.size(); ++i) {
    switch (value[i]) {
      case '&': result += "&amp;"; break;
      case '<': result += "&lt;"; break;
      case '>': result += "&gt;"; break;
      case '"': result += "&quot;"; break;
      case '\'': result += "&apos;"; break;
      default: result += value[i];
    }
  }
  return result;
}

class AtomicReport {
 public:
  std::ofstream stream;
  std::string temporary;
  ~AtomicReport(){if(!temporary.empty())std::remove(temporary.c_str());}
  bool Open(const std::string& path,std::string* error) {
    if(path.empty()) {if(error)*error="structured report path is empty";return false;}
    const auto pattern=path+".tmp.XXXXXX";
    std::vector<char> name(pattern.begin(),pattern.end());name.push_back('\0');
    const int fd=mkstemp(name.data());
    if(fd<0) {if(error)*error="cannot create temporary structured report: "+path;return false;}
    close(fd);temporary=name.data();stream.open(temporary.c_str());
    if(!stream.good()) {if(error)*error="cannot open temporary structured report: "+path;return false;}
    return true;
  }
  bool Commit(const std::string& path,std::string* error) {
    stream.flush();stream.close();
    if(stream.fail()) {if(error)*error="failed while writing structured report: "+path;return false;}
    if(std::rename(temporary.c_str(),path.c_str())!=0) {
      if(error)*error="cannot atomically commit structured report: "+path;
      return false;
    }
    temporary.clear();return true;
  }
};

void WriteJsonFloatingPoint(std::ostream& out, double value) {
  // JSON has no NaN or infinity tokens. Preserve non-finite failure evidence as
  // a quoted spelling instead of emitting invalid JSON or silently coercing it
  // to zero. Consumers should treat number-or-string as an explicit diagnostic
  // value, while ordinary successful metrics remain JSON numbers.
  if (std::isfinite(value)) out << std::setprecision(17) << value;
  else if (std::isnan(value)) out << "\"nan\"";
  else out << (value > 0.0 ? "\"+infinity\"" : "\"-infinity\"");
}
void WriteJsonResult(std::ostream& out,const Result& result) {
  out<<"{\"id\": \""<<JsonEscape(result.id)<<"\", \"status\": \""<<StatusName(result.status)
      <<"\", \"message\": \""<<JsonEscape(result.message)<<"\", \"elapsed_seconds\": ";
  WriteJsonFloatingPoint(out,result.elapsedSeconds);out<<", \"seed\": ";
  if(result.hasSeed)out<<result.seed;else out<<"null";
  out<<", \"configuration\": [";
  for(std::size_t j=0;j<result.configuration.size();++j) {
    if(j)out<<", ";
    out<<'"'<<JsonEscape(result.configuration[j])<<'"';
  }
  out<<"], \"metrics\": [";
  for(std::size_t j=0;j<result.metrics.size();++j) {
    if(j)out<<", ";
    const auto& metric=result.metrics[j];
    out<<"{\"name\": \""<<JsonEscape(metric.name)<<"\", \"value\": ";
    WriteJsonFloatingPoint(out,metric.value);out<<", \"tolerance\": ";
    WriteJsonFloatingPoint(out,metric.tolerance);
    out<<", \"comparison\": \""<<JsonEscape(metric.comparison)<<"\", \"units\": \""
        <<JsonEscape(metric.units)<<"\"}";
  }
  out<<"], \"artifacts\": [";
  for(std::size_t j=0;j<result.artifacts.size();++j) {
    if(j)out<<", ";
    out<<'"'<<JsonEscape(result.artifacts[j])<<'"';
  }
  out<<"], \"invocations\": [";
  for(std::size_t j=0;j<result.invocations.size();++j) {
    if(j)out<<", ";
    out<<"{\"replicate_index\": "<<result.invocations[j].replicateIndex<<", \"result\": ";
    if(result.invocations[j].result)WriteJsonResult(out,*result.invocations[j].result);
    else out<<"null";
    out<<'}';
  }
  out<<"]}";
}

}  // namespace

int Summary::ExitCode() const {
  if (errors != 0) return 2;
  if (failed != 0) return 1;
  return 0;
}

void Summary::Add(const Result& result) {
  results.push_back(result);
  if (result.status == Status::Pass) ++passed;
  else if (result.status == Status::Fail) ++failed;
  else if (result.status == Status::Skip) ++skipped;
  else ++errors;
}

const char* StatusName(Status status) {
  switch (status) {
    case Status::Pass: return "PASS";
    case Status::Fail: return "FAIL";
    case Status::Skip: return "SKIP";
    case Status::Error: return "ERROR";
  }

  return "ERROR";
}

const char* InitializationName(InitializationLevel level) {
  switch (level) {
    case InitializationLevel::None: return "none";
    case InitializationLevel::FieldLineModel: return "field-line-model";
    case InitializationLevel::NativeStorage: return "native-storage";
  }

  return "unknown";
}

const char* RuntimeClassName(RuntimeClass runtime) {
  return runtime == RuntimeClass::Routine ? "routine" : "extended";
}

RefinementOrderEstimate EstimateRefinementOrder(
    double coarseError, double fineError,
    double coarseResolution, double fineResolution) {
  RefinementOrderEstimate result;
  // Zero error cannot define a logarithmic order, even though it may indicate
  // an exact special case.  Refinement tests must choose a nontrivial solution
  // and positive norms so an accidentally exact/empty calculation is visible.
  if (!std::isfinite(coarseError) || !std::isfinite(fineError) ||
      !std::isfinite(coarseResolution) ||
      !std::isfinite(fineResolution) || coarseError <= 0.0 ||
      fineError <= 0.0 || coarseResolution <= fineResolution ||
      fineResolution <= 0.0) {
    result.message =
        "refinement order requires finite positive errors and h_coarse>h_fine>0";
    return result;
  }

  result.coarseToFineErrorRatio = coarseError / fineError;
  result.coarseToFineResolutionRatio =
      coarseResolution / fineResolution;
  result.observedOrder = std::log(result.coarseToFineErrorRatio) /
      std::log(result.coarseToFineResolutionRatio);
  if (!std::isfinite(result.observedOrder)) {
    result.message = "refinement order is not finite";
    return result;
  }
  result.valid = true;
  result.message = "observed refinement order calculated";
  return result;
}

Registry::Registry(std::vector<Descriptor> descriptors)
    : descriptors_(std::move(descriptors)) {
  // Validate once, before a selection can trigger expensive model setup.  A
  // duplicate ID is a programming error: accepting it would make --test ID
  // depend on construction/link order and therefore break reproducibility.
  std::set<std::string> foldedIds;
  for (const Descriptor& descriptor : descriptors_) {
    if (descriptor.id.empty() || descriptor.name.empty() ||
        descriptor.group.empty() || descriptor.description.empty() ||
        descriptor.supportedBuildModes.empty() ||
        descriptor.seedPolicy.empty() || descriptor.stateIsolation.empty() ||
        !descriptor.callback) {
      throw std::invalid_argument(
          "component-test descriptors require complete metadata and a callback");
    }

    const std::string foldedId = FoldCase(descriptor.id);
    if (!foldedIds.insert(foldedId).second) {
      throw std::invalid_argument("duplicate component-test ID: " + descriptor.id);
    }
    std::set<std::string> tags;
    for(const auto& tag:descriptor.tags)
      if(tag.empty()||!tags.insert(FoldCase(tag)).second)
        throw std::invalid_argument("test tags must be nonempty and unique: "+descriptor.id);
  }

  std::sort(descriptors_.begin(), descriptors_.end(), DescriptorIdLess);
}

std::vector<const Descriptor*> Registry::Select(
    const std::vector<std::string>& ids,
    const std::vector<std::string>& groups,
    bool allRoutine,const std::vector<std::string>& expressions) const {
  std::set<std::string> requestedIds;
  std::set<std::string> requestedGroups;
  for (const std::string& id : ids) requestedIds.insert(FoldCase(id));
  for (const std::string& group : groups) requestedGroups.insert(FoldCase(group));
  bool allComplete=false;
  std::set<std::string> requestedTags;
  for(const auto& expression:expressions) {
    std::string normalized=expression;
    std::replace(normalized.begin(),normalized.end(),',',' ');
    std::istringstream tokens(normalized);std::string token;bool any=false;
    while(tokens>>token) {
      any=true;token=FoldCase(token);
      if(token=="all") {allComplete=true;continue;}
      if(token.compare(0,3,"id:")==0) {requestedIds.insert(token.substr(3));continue;}
      if(token.compare(0,6,"group:")==0) {requestedGroups.insert(token.substr(6));continue;}
      if(token.compare(0,4,"tag:")==0) {requestedTags.insert(token.substr(4));continue;}
      bool id=false,group=false,tag=false;
      for(const auto& descriptor:descriptors_) {
        id=id||FoldCase(descriptor.id)==token;group=group||FoldCase(descriptor.group)==token;
        for(const auto& item:descriptor.tags)tag=tag||FoldCase(item)==token;
      }
      if(id)requestedIds.insert(token);
      else if(group)requestedGroups.insert(token);
      else if(tag)requestedTags.insert(token);
      else throw std::invalid_argument("unknown component-test expression: "+token);
    }
    if(!any)throw std::invalid_argument("empty component-test expression");
  }

  // Diagnose every unknown selector before initialization.  Silent empty group
  // matches are especially dangerous because they could otherwise turn a test
  // command into a successful no-op.
  for (const std::string& requestedId : requestedIds) {
    bool found = false;
    for (const Descriptor& descriptor : descriptors_) {
      if (FoldCase(descriptor.id) == requestedId) found = true;
    }
    if (!found) throw std::invalid_argument("unknown component-test ID: " + requestedId);
  }

  for (const std::string& requestedGroup : requestedGroups) {
    bool found = false;
    for (const Descriptor& descriptor : descriptors_) {
      if (FoldCase(descriptor.group) == requestedGroup) found = true;
    }
    if (!found) {
      throw std::invalid_argument("unknown component-test group: " + requestedGroup);
    }
  }
  for(const auto& requested:requestedTags) {
    bool found=false;
    for(const auto& descriptor:descriptors_)for(const auto& tag:descriptor.tags)
      found=found||FoldCase(tag)==requested;
    if(!found)throw std::invalid_argument("unknown component-test tag: "+requested);
  }

  std::vector<const Descriptor*> selected;
  for (const Descriptor& descriptor : descriptors_) {
    const bool selectedById =
        requestedIds.count(FoldCase(descriptor.id)) != 0;
    const bool selectedByGroup =
        requestedGroups.count(FoldCase(descriptor.group)) != 0;
    const bool selectedByAll =
        allComplete||(allRoutine && descriptor.runtime == RuntimeClass::Routine);
    bool selectedByTag=false;
    for(const auto& tag:descriptor.tags)selectedByTag=selectedByTag||requestedTags.count(FoldCase(tag))!=0;

    if (selectedById || selectedByGroup || selectedByAll||selectedByTag) {
      selected.push_back(&descriptor);
    }
  }

  std::sort(selected.begin(), selected.end(), DescriptorPointerIdLess);
  selected.erase(std::unique(selected.begin(), selected.end()), selected.end());
  return selected;
}

Result Registry::RunOne(const Descriptor& descriptor) const {
  seedUses.clear();
  const std::chrono::steady_clock::time_point start =
      std::chrono::steady_clock::now();
  Result result;

  try {
    result = descriptor.callback();
  }
  catch (const std::exception& exception) {
    result.status = Status::Error;
    result.message = std::string("unexpected exception: ") + exception.what();
  }
  catch (...) {
    result.status = Status::Error;
    result.message = "unexpected non-standard exception";
  }

  result.id = descriptor.id;
  result.elapsedSeconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  if (result.message.empty()) result.message = "test returned no diagnostic message";

  // A legacy callback may return PASS while retaining an old hidden failure
  // counter.  Treat any positive/non-finite assertion_failures metric as FAIL,
  // making it impossible for a converted diagnostic to hide an assertion from
  // the process exit code and structured reports.
  if (result.status == Status::Pass) {
    for (std::size_t i = 0; i < result.metrics.size(); ++i) {
      if (result.metrics[i].name == "assertion_failures" &&
          (!std::isfinite(result.metrics[i].value) ||
           result.metrics[i].value > 0.0)) {
        result.status = Status::Fail;
        result.message = "callback reported PASS with assertion failures: " +
                         result.message;
        break;
      }
    }
  }
  if(executionContext.seedActive) {
    result.configuration.push_back("test_seed_schema="+std::to_string(TestSeedSchema));
    result.configuration.push_back("test_campaign_seed="+std::to_string(executionContext.campaignSeed));
    result.configuration.push_back("test_case_seed="+std::to_string(executionContext.caseSeed));
    for(const auto& use:seedUses)result.configuration.push_back(
        "test_stream_seed default="+std::to_string(use.original)+" purpose="+
        std::to_string(use.purpose)+" actual="+std::to_string(use.actual));
    if(descriptor.seedCapability==SeedCapability::Context&&seedUses.empty()&&
       result.status==Status::Pass) {
      result.status=Status::Error;
      result.message="context-seeded callback did not apply TestSeed to an actual stream";
    }
  }
  return result;
}

Summary Registry::RunCampaign(const std::vector<const Descriptor*>& selected,
    const CampaignOptions& options,std::ostream& out,
    const std::function<Result(const Descriptor&)>& runner) const {
  if(options.repeatCount==0)throw std::invalid_argument("test-repeat must be positive");
  if(!selected.empty()&&options.repeatCount>std::numeric_limits<std::size_t>::max()/selected.size())
    throw std::invalid_argument("test selection times repetition count overflows campaign storage");
  const ExecutionContext saved=executionContext;
  struct Restore {ExecutionContext context;~Restore(){executionContext=context;seedUses.clear();}} restore{saved};
  Summary summary;bool stopped=false;
  const auto severity=[](Status status){return status==Status::Error?3:status==Status::Fail?2:status==Status::Skip?1:0;};
  for(const Descriptor* descriptor:selected) {
    if(!descriptor)throw std::invalid_argument("campaign selection contains null descriptor");
    if(stopped) {
      summary.notRun.push_back(descriptor->id+"#0.."+std::to_string(options.repeatCount-1));
      summary.notRunInvocationCount+=options.repeatCount;continue;
    }
    Result aggregate;aggregate.id=descriptor->id;aggregate.status=Status::Pass;
    for(std::uint32_t replicate=0;replicate<options.repeatCount;++replicate) {
      executionContext=saved;seedUses.clear();
      executionContext.testId=descriptor->id;executionContext.campaignSeed=options.campaignSeed;
      executionContext.caseSeed=DeriveTestCaseSeed(options.campaignSeed,descriptor->id,replicate);
      executionContext.replicateIndex=replicate;executionContext.replicateCount=options.repeatCount;
      executionContext.seedActive=(options.seedProvided||options.repeatCount>1)&&
          descriptor->seedCapability==SeedCapability::Context;
      if(!options.artifactRoot.empty())executionContext.artifactDirectory=options.artifactRoot+
          "/"+descriptor->id+"/replicate-"+std::to_string(replicate);
      Result result;
      const bool unsupported=(options.seedProvided&&descriptor->seedCapability==SeedCapability::FixedContract)||
         (options.repeatCount>1&&!descriptor->repeatSupported);
      if(unsupported) {
        result.status=Status::Error;result.id=descriptor->id;
        result.message="requested seed/repeat is unsupported by this fixed or unconverted test contract";
      } else result=runner?runner(*descriptor):RunOne(*descriptor);
      Invocation invocation;invocation.replicateIndex=replicate;
      invocation.result=std::make_shared<Result>(result);aggregate.invocations.push_back(invocation);
      if(aggregate.invocations.size()==1) {
        aggregate.message=result.message;aggregate.hasSeed=result.hasSeed;aggregate.seed=result.seed;
        aggregate.configuration=result.configuration;aggregate.metrics=result.metrics;aggregate.artifacts=result.artifacts;
      } else {
        aggregate.hasSeed=aggregate.hasSeed&&result.hasSeed&&aggregate.seed==result.seed;
        aggregate.metrics.clear();aggregate.artifacts.insert(aggregate.artifacts.end(),result.artifacts.begin(),result.artifacts.end());
      }
      if(severity(result.status)>severity(aggregate.status))aggregate.status=result.status;
      aggregate.elapsedSeconds+=result.elapsedSeconds;
      if(unsupported||(options.failFast&&(result.status==Status::Fail||result.status==Status::Error))) {
        // An unsupported contract cannot execute any requested replicate.
        // Report one authoritative ERROR, then retain the other requested
        // invocations as a compact range even without global fail-fast. This
        // avoids allocating billions of identical contract errors while still
        // continuing to the next selected descriptor when fail-fast is off.
        stopped=options.failFast;
        const std::uint64_t remaining=options.repeatCount-replicate-1;
        if(remaining) {
          summary.notRun.push_back(descriptor->id+"#"+std::to_string(replicate+1)+".."+
              std::to_string(options.repeatCount-1));
          summary.notRunInvocationCount+=remaining;
        }
        break;
      }
    }
    if(!aggregate.invocations.empty()) {
      if(options.repeatCount>1)aggregate.configuration.push_back("test_replicates_executed="+
          std::to_string(aggregate.invocations.size())+" requested="+std::to_string(options.repeatCount));
      summary.Add(aggregate);PrintResult(aggregate,out);
    }
  }
  out<<"Component-test summary: PASS="<<summary.passed<<" FAIL="<<summary.failed
      <<" SKIP="<<summary.skipped<<" ERROR="<<summary.errors<<'\n';
  return summary;
}

Summary Registry::Run(const std::vector<const Descriptor*>& selected,
                      std::ostream& out) const {
  Summary summary;
  for (const Descriptor* descriptor : selected) {
    if (!descriptor) {
      Result result;
      result.status = Status::Error;
      result.id = "<null>";
      result.message = "selection contains a null test descriptor";
      summary.Add(result);
      PrintResult(result, out);
      continue;
    }

    Result result = RunOne(*descriptor);
    summary.Add(result);
    PrintResult(result, out);
  }

  out << "Component-test summary: PASS=" << summary.passed
      << " FAIL=" << summary.failed << " SKIP=" << summary.skipped
      << " ERROR=" << summary.errors << '\n';
  return summary;
}

bool WriteJsonSummary(const Summary& summary,const std::string& path,std::string* error) {
  AtomicReport transaction;
  if(!transaction.Open(path,error))return false;
  auto& out=transaction.stream;
  out<<"{\n  \"schema\": \"srcsep-component-tests-v1\",\n"
      <<"  \"test_seed_schema\": "<<TestSeedSchema<<",\n"
      <<"  \"exit_code\": "<<summary.ExitCode()<<",\n"
      <<"  \"not_run_invocation_count\": "<<summary.notRunInvocationCount<<",\n"
      <<"  \"totals\": {\"passed\": "<<summary.passed<<", \"failed\": "<<summary.failed
      <<", \"skipped\": "<<summary.skipped<<", \"errors\": "<<summary.errors<<"},\n"
      <<"  \"not_run\": [";
  for(std::size_t i=0;i<summary.notRun.size();++i) {
    if(i)out<<", ";
    out<<'"'<<JsonEscape(summary.notRun[i])<<'"';
  }
  out<<"],\n  \"results\": [\n";
  for(std::size_t i=0;i<summary.results.size();++i) {
    out<<"    ";WriteJsonResult(out,summary.results[i]);
    out<<(i+1==summary.results.size()?"\n":",\n");
  }
  out<<"  ]\n}\n";
  return transaction.Commit(path,error);
}

bool WriteJUnitSummary(const Summary& summary,const std::string& path,std::string* error) {
  std::vector<std::pair<std::string,const Result*>> executions;
  unsigned failures=0,errors=0,skipped=0;
  for(const auto& aggregate:summary.results) {
    if(aggregate.invocations.empty())executions.emplace_back(aggregate.id,&aggregate);
    else for(const auto& invocation:aggregate.invocations) {
      if(!invocation.result) {if(error)*error="null invocation in structured report";return false;}
      executions.emplace_back(aggregate.id+"[replicate="+std::to_string(invocation.replicateIndex)+"]",invocation.result.get());
    }
  }
  for(const auto& execution:executions) {
    failures+=execution.second->status==Status::Fail;
    errors+=execution.second->status==Status::Error;
    skipped+=execution.second->status==Status::Skip;
  }
  AtomicReport transaction;
  if(!transaction.Open(path,error))return false;
  auto& out=transaction.stream;
  out<<"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
      <<"<testsuite name=\"srcSEP component tests\" tests=\""<<executions.size()
      <<"\" failures=\""<<failures<<"\" errors=\""<<errors<<"\" skipped=\""<<skipped<<"\">\n";
  for(const auto& execution:executions) {
    const Result& result=*execution.second;
    out<<"  <testcase classname=\"srcSEP."<<XmlEscape(result.id)<<"\" name=\""
        <<XmlEscape(execution.first)<<"\" time=\""<<std::setprecision(17)<<result.elapsedSeconds<<"\">\n";
    if(result.status==Status::Fail)out<<"    <failure message=\""<<XmlEscape(result.message)<<"\"/>\n";
    else if(result.status==Status::Error)out<<"    <error message=\""<<XmlEscape(result.message)<<"\"/>\n";
    else if(result.status==Status::Skip)out<<"    <skipped message=\""<<XmlEscape(result.message)<<"\"/>\n";
    out<<"    <system-out>";
    if(result.hasSeed)out<<"seed="<<result.seed<<"\n";
    for(const auto& setting:result.configuration)out<<XmlEscape(setting)<<"\n";
    for(const auto& metric:result.metrics)out<<XmlEscape(metric.name)<<'='<<metric.value<<"\n";
    for(const auto& artifact:result.artifacts)out<<"artifact="<<XmlEscape(artifact)<<"\n";
    out<<"</system-out>\n  </testcase>\n";
  }
  if(!summary.notRun.empty()) {
    out<<"  <system-out>";
    for(const auto& id:summary.notRun)out<<"not_run="<<XmlEscape(id)<<"\n";
    out<<"</system-out>\n";
  }
  out<<"</testsuite>\n";
  return transaction.Commit(path,error);
}

void Registry::PrintList(std::ostream& out) const {
  out << "Registered srcSEP standalone component tests\n";
  out << "ID | group | class | initialization | build modes | description\n";
  for (const Descriptor& descriptor : descriptors_) {
    out << descriptor.id << " | " << descriptor.group << " | "
        << RuntimeClassName(descriptor.runtime) << " | "
        << InitializationName(descriptor.initialization) << " | "
        << descriptor.supportedBuildModes << " | "
        << descriptor.description << '\n';
    out << "  seed: " << descriptor.seedPolicy << '\n'
        << "  state/isolation: " << descriptor.stateIsolation << '\n';
    out<<"  seed-capability: "<<(descriptor.seedCapability==SeedCapability::Context?"context":
        descriptor.seedCapability==SeedCapability::Deterministic?"deterministic":"fixed-or-unconverted")
       <<"; repeat-supported: "<<(descriptor.repeatSupported?"yes":"no")<<'\n';
    out<<"  tags:";
    for(const auto& tag:descriptor.tags)out<<' '<<tag;
    out<<'\n';
  }
}

void PrintResult(const Result& result, std::ostream& out) {
  out << '[' << result.id << "] " << StatusName(result.status)
      << " duration_s=" << std::fixed << std::setprecision(6)
      << result.elapsedSeconds << " message=\"" << result.message << "\"\n";

  if (result.hasSeed) out << "  seed=" << result.seed << '\n';
  for (const std::string& setting : result.configuration) {
    out << "  configuration " << setting << '\n';
  }
  for (const Metric& metric : result.metrics) {
    out << "  metric " << metric.name << '=' << std::setprecision(12)
        << metric.value << ' ' << metric.comparison << ' ' << metric.tolerance;
    if (!metric.units.empty()) out << ' ' << metric.units;
    out << '\n';
  }
  for (const std::string& artifact : result.artifacts) {
    out << "  artifact " << artifact << '\n';
  }
}

}  // namespace Testing
}  // namespace SEP
