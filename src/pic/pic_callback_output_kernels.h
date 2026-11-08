#ifndef _AMPS_PIC_CALLBACK_OUTPUT_KERNELS_H_
#define _AMPS_PIC_CALLBACK_OUTPUT_KERNELS_H_

/* Checked callback dispatch and deterministic trajectory base formatting. */
#include <cmath>
#include <cstdio>
#include <functional>
#include <type_traits>
#include <utility>

namespace PIC { namespace CallbackOutputKernels {

inline const char* TrajectoryBaseVariables() {
  return "VARIABLES=\"x\", \"y\", \"z\", \"spec\", \"Speed\", \"vx\", \"vy\", \"vz\", \"Kinetic Energy [eV]\"";
}

inline bool WriteTrajectoryBaseVariables(FILE* output) {
  return output!=nullptr && std::fprintf(output,"%s",TrajectoryBaseVariables())>=0;
}

template<class Container>
inline bool InvokeInRegistrationOrder(Container& callbacks) {
  for (auto& callback:callbacks) {
    if constexpr (std::is_pointer_v<std::decay_t<decltype(callback)>>) {
      if (callback==nullptr) return false;
    }
    std::invoke(callback);
  }
  return true;
}

template<class Function,class... Args>
inline bool InvokeExactlyOnce(Function&& callback,Args&&... args) {
  if constexpr (std::is_pointer_v<std::decay_t<Function>>) {
    if (callback==nullptr) return false;
  }
  std::invoke(std::forward<Function>(callback),std::forward<Args>(args)...);
  return true;
}

inline bool WriteTrajectoryBaseRecord(FILE* output,const double x[3],int species,
                                      double speed,const double velocity[3],
                                      double kinetic_energy_eV) {
  if (!output || !x || !velocity || species<0 || !std::isfinite(speed) ||
      !std::isfinite(kinetic_energy_eV)) return false;
  for (int d=0;d<3;++d)
    if (!std::isfinite(x[d]) || !std::isfinite(velocity[d])) return false;
  return std::fprintf(output,"%e  %e  %e  %i  %e   %e  %e  %e %e",
      x[0],x[1],x[2],species,speed,velocity[0],velocity[1],velocity[2],
      kinetic_energy_eV)>=0;
}
}} // namespace PIC::CallbackOutputKernels
#endif
