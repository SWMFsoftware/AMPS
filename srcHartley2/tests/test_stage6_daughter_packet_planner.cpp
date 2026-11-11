#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <vector>

#include "srcHartley2/DaughterPacketPlanner.h"

namespace Packet=Comet::DaughterPacketPlanner;

namespace {

void Require(bool condition,const char *message) {
  if (!condition) throw std::runtime_error(message);
}

template<class Exception,class Function>
void RequireThrows(Function function,const char *message) {
  bool threw=false;
  try { function(); }
  catch (const Exception&) { threw=true; }
  Require(threw,message);
}

double ProductSum(const Packet::cPlan &plan,std::size_t product) {
  double result=0.0;
  for (std::size_t packet=0;packet<plan.PacketCount;packet++)
    result+=plan.Correction(packet,product);
  return result;
}

void TestIdentityAndExactSplits() {
  const std::vector<double> corrections={0.7,10.25,1.0/3.0};
  for (std::size_t packets:{1U,2U,4U,8U}) {
    const Packet::cPlan plan=Packet::BuildPlan(
        corrections,packets,8,24,23);
    Require(plan.ConcreteProductCount==3 && plan.PacketCount==packets,
        "packet plan dimensions are incorrect");
    Require(plan.TotalDaughterCount==3*packets &&
        plan.RequiredAdditionalParticles==3*packets-1,
        "packet transaction allocation is incorrect");
    for (std::size_t product=0;product<corrections.size();product++)
      Require(ProductSum(plan,product)==corrections[product],
          "packet corrections do not reconstruct the concrete product");
  }
}

void TestRepeatedProductsCapacityAndSinks() {
  // Three entries model CO + H + H.  Repeated hydrogen products must remain
  // separate columns in every complete packet.
  const Packet::cPlan repeated=Packet::BuildPlan(
      {2.0,0.5,0.5},8,8,24,23);
  Require(repeated.TotalDaughterCount==24 &&
      repeated.RequiredAdditionalParticles==23,
      "repeated-product packet count is incorrect");
  Require(ProductSum(repeated,1)==0.5 && ProductSum(repeated,2)==0.5,
      "repeated concrete products were aggregated");

  const Packet::cPlan sink=Packet::BuildPlan({},1,8,24,0);
  Require(sink.TotalDaughterCount==0 &&
      sink.RequiredAdditionalParticles==0 && sink.Corrections.empty(),
      "one-packet sink did not remain empty");

  RequireThrows<std::invalid_argument>([]() {
    Packet::BuildPlan({},2,8,24,0);
  },"split sink was accepted");
  RequireThrows<std::length_error>([]() {
    Packet::BuildPlan({1.0,1.0,1.0},8,8,23,23);
  },"transaction hard limit was ignored");
  RequireThrows<std::length_error>([]() {
    Packet::BuildPlan({1.0,1.0,1.0},8,8,24,22);
  },"capacity shortage was ignored");
  RequireThrows<std::length_error>([]() {
    Packet::BuildPlan({1.0},9,8,24,23);
  },"packet maximum was ignored");
  RequireThrows<std::invalid_argument>([]() {
    Packet::BuildPlan({std::numeric_limits<double>::infinity()},2,8,24,23);
  },"non-finite correction was accepted");
}

} // namespace

int main() {
  try {
    TestIdentityAndExactSplits();
    TestRepeatedProductsCapacityAndSinks();
  }
  catch (const std::exception &exception) {
    std::fprintf(stderr,"Stage-6 daughter-packet planner failed: %s\n",
        exception.what());
    return 1;
  }
  std::printf("PASS Stage-6 daughter-packet planner identity=1 exact-sum=1 "
      "packet-counts=1,2,4,8 repeated-products=1 capacity=1 sink=1\n");
  return 0;
}
