#ifndef SRCMOON_I17_AUTHORIZATION_H_
#define SRCMOON_I17_AUTHORIZATION_H_

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <string>

namespace Moon {
namespace I17Authorization {

// Native accident guard for the reviewed December-only adapter. The Python
// runner is the authority for complete JSON approval, actual SHA-256 values,
// unchanged committed review, model scope, sources/build and one-time release.
// It writes this small fixed key/value document only AFTER those checks and
// selected-only staging. AMPS has no duplicate JSON or cryptographic library:
// this helper checks the complete sidecar schema/status/path/hash-token shape
// before RunI17 opens a science CSV. It does not authenticate a signature or
// claim observational qualification. Units do not enter this control record;
// all paths identify exact staged files in the existing campaign directory.
inline bool Validate(const std::string& authorizationFile,
    const std::string& historyFile,const std::string& freezeFile,
    std::string* error) {
  if (error==nullptr) return false;
  error->clear();
  const auto reject=[&](const std::string& reason) {
    *error="I17 hold-out authorization: "+reason;
    return false;
  };
  if (authorizationFile.empty() || historyFile.empty() || freezeFile.empty())
    return reject("explicit authorization/history/review paths required");
  std::ifstream input(authorizationFile);
  if (!input) return reject("native release sidecar unavailable; use the reviewed I17 runner");
  const std::set<std::string> required={"schema","selection","review_status",
      "history_file","history_sha256","freeze_file","freeze_sha256",
      "validation_file","validation_sha256","source_identity_sha256",
      "physics_identity_sha256","review_token"};
  std::map<std::string,std::string> fields;
  std::string line;
  while (std::getline(input,line)) {
    const auto at=line.find('=');
    if (at==std::string::npos || at==0 || at+1==line.size())
      return reject("empty/malformed key=value record");
    const std::string key=line.substr(0,at),value=line.substr(at+1);
    if (!required.count(key) || fields.count(key))
      return reject("unknown or duplicate field: "+key);
    fields.emplace(key,value);
  }
  if (!input.eof() || fields.size()!=required.size())
    return reject("incomplete or unreadable native release record");
  if (fields["schema"]!="srcMoon-I17-native-authorization-v1" ||
      fields["selection"]!="holdout" || fields["review_status"]!="APPROVED_FROZEN")
    return reject("schema, reviewed state or selected stream disagrees");
  const auto samePath=[](const std::string& first,const std::string& second) {
    const std::filesystem::path a(first),b(second);
    return a.is_absolute() && b.is_absolute() && a.lexically_normal()==b.lexically_normal();
  };
  if (!samePath(fields["history_file"],historyFile) ||
      !samePath(fields["freeze_file"],freezeFile) ||
      !samePath(fields["validation_file"],
          (std::filesystem::path(authorizationFile).parent_path()/
              "holdout_release_validation.json").string()))
    return reject("release names different staged history/review/validation files");
  const auto hexadecimal=[](const std::string& value) {
    return !value.empty() && std::all_of(value.begin(),value.end(),[](char c) {
      return (c>='0' && c<='9') || (c>='a' && c<='f');
    });
  };
  for (const char* key:{"history_sha256","freeze_sha256","validation_sha256",
                       "source_identity_sha256","physics_identity_sha256"})
    if (fields[key].size()!=64 || !hexadecimal(fields[key]))
      return reject(std::string("invalid SHA-256 identity field: ")+key);
  const std::string prefix=fields["freeze_sha256"]+"/";
  if (fields["review_token"].compare(0,prefix.size(),prefix)!=0)
    return reject("review token differs from frozen review identity");
  const std::string blob=fields["review_token"].substr(prefix.size());
  if ((blob.size()!=40 && blob.size()!=64) || !hexadecimal(blob))
    return reject("review token lacks a committed Git blob identity");
  for (const char* key:{"history_file","freeze_file","validation_file"}) {
    std::error_code code;
    if (!std::filesystem::is_regular_file(fields[key],code) || code)
      return reject(std::string("staged file is unavailable: ")+key);
  }
  return true;
}

} // namespace I17Authorization
} // namespace Moon
#endif
