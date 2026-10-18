//Generated COPS products live outside the application: <AMPS>/cops_output.
//Absolute output paths outside the application remain overrides. Input search never moves old data.
#ifndef LIOUVILLE_COPS_PATHS_H
#define LIOUVILLE_COPS_PATHS_H

#include <cstdlib>
#include <filesystem>
#include <cstdio>
#include <string>

namespace CopsPaths {
namespace fs = std::filesystem;

inline fs::path ApplicationRoot() {
  //On Linux this also handles invoking the binary through PATH or a symlink.
  std::error_code error;
  fs::path executable=fs::read_symlink("/proc/self/exe",error);
  return error ? fs::current_path() : executable.parent_path();
}

inline fs::path OutputRoot() {
  const char* configured=std::getenv("COPS_OUTPUT_ROOT");
  if (configured && *configured) return fs::absolute(configured);
  return ApplicationRoot().parent_path().parent_path()/"cops_output";
}

inline bool IsWithin(const fs::path& path,const fs::path& root) {
  auto iPath=path.begin(),iRoot=root.begin();
  for (;iRoot!=root.end();++iRoot,++iPath)
    if (iPath==path.end() || *iPath!=*iRoot) return false;
  return true;
}

inline std::string OutputPath(const std::string& name) {
  fs::path path(name);
  fs::path candidate=(path.is_absolute() ? path : OutputRoot()/path).lexically_normal();
  std::error_code error;
  fs::path application=fs::weakly_canonical(ApplicationRoot(),error);
  if (error) application=fs::absolute(ApplicationRoot()).lexically_normal();
  error.clear();
  fs::path resolved=fs::weakly_canonical(candidate,error);
  if (error) resolved=fs::absolute(candidate).lexically_normal();
  if (IsWithin(resolved,application)) {
    std::fprintf(stderr,"Error: generated output must be outside LiouvilleCOPS: %s\n",resolved.c_str());
    std::exit(EXIT_FAILURE);
  }
  return candidate.string();
}

inline std::string ProductInput(const std::string& name) {
  fs::path path(name);
  if (path.is_absolute() || name.rfind("./",0)==0 || name.rfind("../",0)==0)
    return fs::absolute(path).lexically_normal().string();
  fs::path primary=OutputRoot()/path;
  if (fs::exists(primary)) return primary.string();
  if (fs::exists(path)) return fs::absolute(path).string();
  fs::path legacy=ApplicationRoot()/path;
  return (fs::exists(legacy) ? legacy : primary).string();
}
}
#endif
