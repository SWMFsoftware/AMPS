//======================================================================================
// EngineGlue.cpp -- verification helpers (see EngineGlue.h).
//======================================================================================

#include "EngineGlue.h"

#include <atomic>
#include <cstdio>
#include <mutex>
#include <set>

namespace Geospace {
namespace GPU {

namespace {

std::atomic<long long> gVerified(0);
std::atomic<long long> gMismatches(0);
std::mutex gWarnMutex;
std::set<std::string> gWarned;

bool SameBits(double a,double b) { return std::memcmp(&a,&b,sizeof(double))==0; }

struct VerifyReporter {
  ~VerifyReporter() {
    if (gVerified.load()>0) {
      std::printf("[earth-engine][verify] trajectories compared=%lld bitwise mismatches=%lld\n",
                  gVerified.load(),gMismatches.load());
      std::fflush(stdout);
    }
  }
} gReporter;

} // namespace

bool VerifyEnabled() {
  static int enabled=-1;
  if (enabled<0) {
    const char* e=std::getenv("AMPS_EARTH_ENGINE_VERIFY");
    enabled=(e!=nullptr && std::string(e)=="1") ? 1 : 0;
  }
  return enabled==1;
}

void VerifyCompare(const char* where,double R_GV,const Geospace::Trajectory::Result& e,
                   const Geospace::Trajectory::Result& l) {
  ++gVerified;
  std::string diff;
  auto D=[&](const char* name,double a,double b) {
    if (!SameBits(a,b)) {
      char buf[256];
      std::snprintf(buf,sizeof(buf)," %s engine=%.17g legacy=%.17g;",name,a,b);
      diff+=buf;
    }
  };
  auto I=[&](const char* name,long long a,long long b) {
    if (a!=b) {
      char buf[160];
      std::snprintf(buf,sizeof(buf)," %s engine=%lld legacy=%lld;",name,a,b);
      diff+=buf;
    }
  };
  I("termination",static_cast<int>(e.termination),static_cast<int>(l.termination));
  I("steps",e.steps,l.steps);
  D("traceTime_s",e.traceTime_s,l.traceTime_s);
  D("traceDistance_m",e.traceDistance_m,l.traceDistance_m);
  I("exit.valid",e.exitState.valid,l.exitState.valid);
  for (int d=0;d<3;++d) {
    D("exit.x",e.exitState.x_exit_m[d],l.exitState.x_exit_m[d]);
    D("exit.p",e.exitState.p_exit_SI[d],l.exitState.p_exit_SI[d]);
  }
  D("exit.cosAlpha",e.exitState.cosAlpha,l.exitState.cosAlpha);
  D("exit.t",e.exitState.traceTimeAtExit_s,l.exitState.traceTimeAtExit_s);
  I("mirrorPoints",e.mirrorPoints,l.mirrorPoints);
  I("bounceCycles",e.bounceCycles,l.bounceCycles);
  I("driftRevolutions",e.driftRevolutions,l.driftRevolutions);
  D("driftAngle",e.driftAngle_rad,l.driftAngle_rad);
  D("minimumOuterMargin",e.minimumOuterMargin_m,l.minimumOuterMargin_m);
  D("momentumRelativeSpread",e.momentumRelativeSpread,l.momentumRelativeSpread);
  I("trapMechanism",e.trapMechanism,l.trapMechanism);
  D("driftMatchedBinFraction",e.driftMatchedBinFraction,l.driftMatchedBinFraction);
  I("driftTurnStates",static_cast<long long>(e.driftTurnStates.size()),
    static_cast<long long>(l.driftTurnStates.size()));
  if (!diff.empty()) {
    const long long n=++gMismatches;
    if (n<=20) {
      std::printf("[earth-engine][verify][%s] R=%.17g mismatch:%s\n",where,R_GV,diff.c_str());
      std::fflush(stdout);
    }
  }
}

void WarnOnce(const char* key,const std::string& message) {
  std::lock_guard<std::mutex> lock(gWarnMutex);
  if (gWarned.insert(key).second) {
    std::fprintf(stderr,"[earth-engine] WARNING: %s\n",message.c_str());
    std::fflush(stderr);
  }
}

} // namespace GPU
} // namespace Geospace
