#include "test_harness.h"
#include "pic_callback_output_kernels.h"
#include <cstdio>
#include <functional>
#include <limits>
#include <string>
#include <vector>
namespace C=PIC::CallbackOutputKernels;

AMPS_TEST_CASE("PIC-CALLBACK-RUNTIME-ORDER-001","runtime callbacks execute once in registration order") {
  std::vector<int> order;
  std::vector<std::function<void()>> callbacks;
  callbacks.emplace_back([&](){order.push_back(3);}); callbacks.emplace_back([&](){order.push_back(1);}); callbacks.emplace_back([&](){order.push_back(4);});
  AMPS_EXPECT_TRUE(C::InvokeInRegistrationOrder(callbacks));
  AMPS_EXPECT_TRUE(order==std::vector<int>({3,1,4}));
}
AMPS_TEST_CASE("PIC-CALLBACK-PARTICLE-ONCE-001","user-particle callback boundary invokes the selected function exactly once") {
  int count=0,value=7;
  AMPS_EXPECT_TRUE(C::InvokeExactlyOnce([&](int increment){++count; value+=increment;},5));
  AMPS_EXPECT_EQ(count,1); AMPS_EXPECT_EQ(value,12);
}
AMPS_TEST_CASE("PIC-TRAJECTORY-FORMAT-001","trajectory base record has deterministic variable order and formatting") {
  FILE* stream=std::tmpfile(); AMPS_EXPECT_TRUE(stream!=nullptr);
  const double x[3]={1.0,-2.5,0.125},v[3]={4.0,5.0,-6.0};
  AMPS_EXPECT_TRUE(C::WriteTrajectoryBaseRecord(stream,x,2,9.0,v,11.5));
  std::fflush(stream); std::rewind(stream); char buffer[512]={0};
  AMPS_EXPECT_TRUE(std::fgets(buffer,sizeof(buffer),stream)!=nullptr); std::fclose(stream);
  AMPS_EXPECT_TRUE(std::string(buffer)=="1.000000e+00  -2.500000e+00  1.250000e-01  2  9.000000e+00   4.000000e+00  5.000000e+00  -6.000000e+00 1.150000e+01");
}
AMPS_TEST_CASE("PIC-CALLBACK-OUTPUT-INVALID-001","null callbacks and nonfinite trajectory fields are rejected") {
  void (*callback)()=nullptr; AMPS_EXPECT_TRUE(!C::InvokeExactlyOnce(callback));
  FILE* stream=std::tmpfile(); const double x[3]={0,0,0},v[3]={0,0,std::numeric_limits<double>::quiet_NaN()};
  AMPS_EXPECT_TRUE(!C::WriteTrajectoryBaseRecord(stream,x,0,0.0,v,0.0)); std::fclose(stream);
}
