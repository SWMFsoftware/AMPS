/*
 * Deterministic Stage-6 LET and spectral-hardness tests.
 *
 * These tests are intentionally independent of Geant4. They protect the
 * numerical boundary that transport hands to the result writer: explicit
 * non-clamping logarithmic LET bins, exact partial-bin threshold integration,
 * and covariance-aware ratio/Fieller statistics. Stochastic step and surface
 * scoring are exercised by test_stage6_transport.py.
 */

#include "Stage6Physics.hh"
#include "StrictJson.hh"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace {
int failures=0;
void Check(bool condition,const std::string& description){
  if(!condition){ ++failures; std::cerr<<"FAIL: "<<description<<'\n'; }
}
bool Near(double actual,double expected,double relative=1.0e-12){
  return std::abs(actual-expected)<=relative*std::max({1.0e-300,std::abs(actual),
                                                       std::abs(expected)});
}
}

int main(int argc,char** argv){
  if(argc!=2){ std::cerr<<"usage: stage6_core_tests EVIDENCE_DIR\n"; return 2; }
  namespace fs=std::filesystem;
  const fs::path evidence=argv[1];
  fs::create_directories(evidence);
  shield_sim::Json::Array comparisons;

  shield_sim::LogarithmicLETBinning bins(1.0e-3,1.0e4,20);
  Check(bins.BinCount()==140,"seven decades at 20 bins/decade gives 140 bins");
  Check(bins.Locate(9.0e-4)==-1,"LET underflow is not clamped");
  Check(bins.Locate(1.0e4)==139,"final LET edge belongs to the final bin");
  Check(bins.Locate(1.0e4+1.0)==140,"LET overflow is not clamped");
  Check(Near(bins.Lower(20),1.0e-2) && Near(bins.Upper(20),
        std::pow(10.0,-2.0+1.0/20.0)),"logarithmic LET edges are exact by decade");

  // Constant differential fluence over [5,20] has exactly 10*phi above the
  // strict 10-MeV threshold; no histogram-center approximation is involved.
  shield_sim::PowerLawFluenceBin constant{5.0,20.0,10.0,3.0,0.0};
  const double constantIntegral=shield_sim::IntegratePowerLawAbove(constant,10.0);
  Check(Near(constantIntegral,30.0),"constant partial-bin threshold integral");
  comparisons.emplace_back(shield_sim::Json::Object{{"actual",constantIntegral},
    {"expected",30.0},{"id","constant_partial_bin_integral"},{"unit","cm-2"}});

  // For phi(E)=C*(E/Eref)^-2 the primitive is -C*Eref^2/E.
  shield_sim::PowerLawFluenceBin inverseSquare{20.0,200.0,20.0,5.0,-2.0};
  const double inverseExpected=5.0*20.0*20.0*(1.0/50.0-1.0/200.0);
  const double inverseActual=shield_sim::IntegratePowerLawAbove(inverseSquare,50.0);
  Check(Near(inverseActual,inverseExpected),"power-law partial-bin threshold integral");
  comparisons.emplace_back(shield_sim::Json::Object{{"actual",inverseActual},
    {"expected",inverseExpected},{"id","power_law_partial_bin_integral"},{"unit","cm-2"}});

  // Piecewise-constant spectrum: 10--100 contains 90 units and 100--200
  // contains 200 units, so H100/10=200/290 exactly.
  const auto hardness=shield_sim::FoldH100Over10({
    {10.0,100.0,10.0,1.0,0.0},{100.0,200.0,100.0,2.0,0.0}});
  const double expectedHardness=200.0/290.0;
  Check(hardness.defined && Near(hardness.value,expectedHardness),
        "exact analytic H100/10 fold across threshold edges");
  comparisons.emplace_back(shield_sim::Json::Object{{"actual",hardness.value},
    {"expected",expectedHardness},{"id","analytic_H100_10_fold"},{"unit","1"}});

  // Three correlated event pairs (N,D)=(2,3),(4,5),(6,7). The production
  // estimator must use ratio-of-means 12/15, retain covariance, and produce a
  // finite Fieller interval because the denominator is resolved from zero.
  shield_sim::RatioMoments ratio;
  for(const auto pair:{std::pair<double,double>{2.0,3.0},{4.0,5.0},{6.0,7.0}}){
    ++ratio.histories;
    ratio.numeratorSum+=pair.first; ratio.denominatorSum+=pair.second;
    ratio.numeratorSumSquares+=pair.first*pair.first;
    ratio.denominatorSumSquares+=pair.second*pair.second;
    ratio.numeratorDenominatorSumProducts+=pair.first*pair.second;
  }
  const auto estimate=shield_sim::EstimateRatioOfMeans(ratio);
  Check(estimate.status==shield_sim::RatioStatus::Available &&
        Near(estimate.value,0.8) && estimate.covarianceOfMeans>0.0 &&
        estimate.interval95Low<=estimate.value && estimate.interval95High>=estimate.value,
        "correlated ratio-of-means and finite Fieller interval");
  comparisons.emplace_back(shield_sim::Json::Object{{"actual",estimate.value},
    {"expected",0.8},{"id","correlated_ratio_of_means"},{"unit","1"}});

  shield_sim::RatioMoments zero;
  zero.histories=10;
  Check(shield_sim::EstimateRatioOfMeans(zero).status==
        shield_sim::RatioStatus::ZeroDenominator,"zero denominator is undefined");
  shield_sim::RatioMoments unresolved;
  unresolved.histories=2;
  unresolved.numeratorSum=1.0; unresolved.denominatorSum=2.0;
  unresolved.numeratorSumSquares=1.0; unresolved.denominatorSumSquares=4.0;
  unresolved.numeratorDenominatorSumProducts=2.0;
  Check(shield_sim::EstimateRatioOfMeans(unresolved).status==
        shield_sim::RatioStatus::StatisticallyUnresolvedDenominator,
        "Fieller rejects a denominator confidence set containing zero");

  shield_sim::Json document(shield_sim::Json::Object{
    {"comparisons",shield_sim::Json(std::move(comparisons))},
    {"schema_version","shieldSim.validation_evidence.v1"},
    {"test_case","stage6_core"}});
  std::ofstream out(evidence/"stage6_core.json");
  out<<document.DumpCanonical()<<'\n';
  if(failures){ std::cerr<<failures<<" Stage-6 core test(s) failed\n"; return 1; }
  std::cout<<"All Stage-6 core tests passed\n";
  return 0;
}
