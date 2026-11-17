#ifndef SEP_COMMON_SEP_TEST_REGISTRY_H
#define SEP_COMMON_SEP_TEST_REGISTRY_H

#include <cstdint>
#include <functional>
#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

namespace SEP {
namespace Testing {

// A test result is deliberately richer than a Boolean.  The standalone runner
// must distinguish a scientific assertion failure from an unavailable
// prerequisite or an unexpected execution error, because those outcomes have
// different meanings in automated campaigns and in developer diagnostics.
enum class Status { Pass, Fail, Skip, Error };

// The registry tells main.cpp how far initialization must progress before a
// callback is safe.  Step 1 needs only tests that are independent of model state
// and tests that consume the initialized AMPS field-line model.  Adding an enum
// value later is preferable to hiding extra initialization inside a callback.
enum class InitializationLevel { None, FieldLineModel };

// Routine tests are included by --all-tests and by the bounded make test target.
// Extended tests remain discoverable and individually selectable, but are kept
// out of routine CI because their Monte-Carlo sample counts are intentionally
// large.
enum class RuntimeClass { Routine, Extended };
enum class SeedCapability { FixedContract, Deterministic, Context };
const std::uint64_t TestSeedSchema=1;

struct Metric {
  Metric() = default;
  Metric(const std::string& metricName, double metricValue,
         double metricTolerance, const std::string& metricComparison,
         const std::string& metricUnits)
      : name(metricName), value(metricValue), tolerance(metricTolerance),
        comparison(metricComparison), units(metricUnits) {}

  std::string name;
  double value = 0.0;
  double tolerance = 0.0;
  std::string comparison;
  std::string units;
};

// A refinement study should report the convergence order implied by its data,
// rather than compare a hand-derived fine/coarse ratio whose meaning changes
// with the refinement factor.  For errors e_c and e_f evaluated at resolutions
// h_c > h_f, the observed order is
//
//   p = log(e_c/e_f) / log(h_c/h_f).
//
// The estimate retains both ratios so structured test evidence is sufficient
// to reproduce p.  Invalid and degenerate inputs return valid=false with a
// diagnostic; callers must fail the test instead of silently reporting p=0.
struct RefinementOrderEstimate {
  bool valid = false;
  double observedOrder = 0.0;
  double coarseToFineErrorRatio = 0.0;
  double coarseToFineResolutionRatio = 0.0;
  std::string message;
};

RefinementOrderEstimate EstimateRefinementOrder(double coarseError,
                                                 double fineError,
                                                 double coarseResolution,
                                                 double fineResolution);

struct Result;
struct Invocation {
  std::uint32_t replicateIndex=0;
  std::shared_ptr<Result> result;
};
struct Result {
  Status status = Status::Error;
  std::string id;
  std::string message;
  double elapsedSeconds = 0.0;

  // Zero is a valid deterministic seed, so hasSeed records whether the seed is
  // meaningful rather than overloading a numeric sentinel.
  bool hasSeed = false;
  std::uint64_t seed = 0;
  std::vector<std::string> configuration;
  std::vector<Metric> metrics;
  std::vector<std::string> artifacts;
  // Stable aggregate ID remains unique in results[]. Every executed replicate
  // retains its complete result here, so older readers cannot deduplicate away
  // distinct seeds, failures or artifacts. Child results have no nested runs.
  std::vector<Invocation> invocations;
};

// Optional file context for end-to-end cases executed by the linked
// application. Ordinary component tests leave both paths empty. The command
// line owns path selection, main.cpp installs the immutable context before any
// callback runs, and validation callbacks consume it without reading process
// environment variables or inventing a second CLI.
struct ExecutionContext {
  std::string inputPath;
  std::string artifactDirectory;
  // Optional metadata was appended after the original two path fields. Keep
  // its empty default explicit so existing two-field aggregate initializers
  // remain valid under -Werror=missing-field-initializers. This preserves the
  // aggregate API and leaves ordinary non-campaign clients without a test ID.
  std::string testId{};
  bool seedActive=false;
  std::uint64_t campaignSeed=0,caseSeed=0;
  std::uint32_t replicateIndex=0,replicateCount=1;
};

void SetExecutionContext(const ExecutionContext& context);
const ExecutionContext& GetExecutionContext();
// A versioned, order/rank-independent case/replicate key. TestSeed preserves
// legacy defaults when no seed/repeat control is requested, and otherwise
// derives a purpose-tagged actual seed. Distinct old seed constants remain
// distinct domains. It records the actual mapping, not only report metadata.
std::uint64_t DeriveTestCaseSeed(std::uint64_t campaignSeed,
    const std::string& canonicalId,std::uint32_t replicateIndex);
std::uint64_t TestSeed(std::uint64_t defaultSeed,std::uint64_t purpose=0);

using TestCallback = std::function<Result()>;

struct Descriptor {
  std::string id;
  std::string name;
  std::string group;
  std::string description;
  InitializationLevel initialization = InitializationLevel::None;
  std::string supportedBuildModes;
  RuntimeClass runtime = RuntimeClass::Routine;
  std::string seedPolicy;
  std::string stateIsolation;
  SeedCapability seedCapability=SeedCapability::FixedContract;
  bool repeatSupported=false;
  std::vector<std::string> tags;
  TestCallback callback;
};
struct CampaignOptions {
  bool seedProvided=false;
  std::uint64_t campaignSeed=0;
  std::uint32_t repeatCount=1;
  bool failFast=false;
  std::string artifactRoot;
};

struct Summary {
  std::vector<Result> results;
  int passed = 0;
  int failed = 0;
  int skipped = 0;
  int errors = 0;
  std::vector<std::string> notRun;
  std::uint64_t notRunInvocationCount=0;

  // PASS and SKIP produce a successful command status.  A requested SKIP is
  // still displayed explicitly and is never relabelled PASS; FAIL and ERROR
  // return nonzero so Make and batch schedulers cannot hide a broken test.
  int ExitCode() const;

  // Add is the only supported aggregation path.  Keeping status counters next
  // to the stored result prevents a text reporter and a structured reporter
  // from silently describing different outcomes.
  void Add(const Result& result);
};

class Registry {
 public:
  explicit Registry(std::vector<Descriptor> descriptors);

  const std::vector<Descriptor>& Descriptors() const { return descriptors_; }

  // IDs and groups are matched case-insensitively for consistency with the
  // existing SEP CLI value policy.  The returned descriptors are de-duplicated
  // and sorted by canonical ID, independent of selector order or link order.
  std::vector<const Descriptor*> Select(const std::vector<std::string>& ids,
                                        const std::vector<std::string>& groups,
                                        bool allRoutine,
                                        const std::vector<std::string>& expressions={}) const;

  Result RunOne(const Descriptor& descriptor) const;
  Summary Run(const std::vector<const Descriptor*>& selected,
              std::ostream& out) const;
  // Runner optionally supplies native fixture/seed scope and MPI outcome
  // reduction. It must return the globally authoritative result before this
  // loop applies fail-fast; canonical registry code has no PIC/MPI dependency.
  Summary RunCampaign(const std::vector<const Descriptor*>& selected,
      const CampaignOptions& options,std::ostream& out,
      const std::function<Result(const Descriptor&)>& runner={}) const;
  void PrintList(std::ostream& out) const;

 private:
  std::vector<Descriptor> descriptors_;
};

const char* StatusName(Status status);
const char* InitializationName(InitializationLevel level);
const char* RuntimeClassName(RuntimeClass runtime);
void PrintResult(const Result& result, std::ostream& out);

// Structured reports contain the same complete Result objects printed to the
// terminal, including seeds, effective configuration, metrics, and artifacts.
// A write error is returned to the caller and must change the process status;
// test evidence is never best-effort when the operator requested a report.
bool WriteJsonSummary(const Summary& summary, const std::string& path,
                      std::string* error);
bool WriteJUnitSummary(const Summary& summary, const std::string& path,
                       std::string* error);

}  // namespace Testing
}  // namespace SEP

#endif  // SEP_COMMON_SEP_TEST_REGISTRY_H
