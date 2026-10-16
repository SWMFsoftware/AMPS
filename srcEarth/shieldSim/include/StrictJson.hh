#ifndef SHIELDSIM_STRICT_JSON_HH
#define SHIELDSIM_STRICT_JSON_HH

/*
 * StrictJson.hh
 * -------------
 * A deliberately small JSON value and parser used for shieldSim scientific
 * input files.  Keeping this code in-tree avoids making the Geant4 executable
 * depend on whichever JSON package happens to be installed on a workstation.
 *
 * The parser is intentionally stricter than many convenience parsers:
 * duplicate object keys, trailing tokens, invalid UTF-16 surrogate pairs,
 * non-finite numbers, and malformed escape sequences are errors.  Objects use
 * std::map, so DumpCanonical() is deterministic and can be content-hashed.
 */

#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

namespace shield_sim {

class JsonError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

class Json {
public:
  using Array  = std::vector<Json>;
  using Object = std::map<std::string, Json>;
  using Value  = std::variant<std::nullptr_t, bool, double, std::string, Array, Object>;

  Json() : value_(nullptr) {}
  Json(std::nullptr_t) : value_(nullptr) {}
  Json(bool value) : value_(value) {}
  Json(double value) : value_(value) {}
  Json(int value) : value_(static_cast<double>(value)) {}
  Json(std::string value) : value_(std::move(value)) {}
  Json(const char* value) : value_(std::string(value)) {}
  Json(Array value) : value_(std::move(value)) {}
  Json(Object value) : value_(std::move(value)) {}

  static Json Parse(const std::string& text, const std::string& source = "<string>");
  static Json ParseFile(const std::string& path);

  bool IsNull()   const;
  bool IsBool()   const;
  bool IsNumber() const;
  bool IsString() const;
  bool IsArray()  const;
  bool IsObject() const;

  bool               AsBool(const std::string& context = "value") const;
  double             AsNumber(const std::string& context = "value") const;
  const std::string& AsString(const std::string& context = "value") const;
  const Array&       AsArray(const std::string& context = "value") const;
  const Object&      AsObject(const std::string& context = "value") const;
  Array&             AsArray(const std::string& context = "value");
  Object&            AsObject(const std::string& context = "value");

  bool Contains(const std::string& key) const;
  const Json& At(const std::string& key, const std::string& context = "object") const;
  Json& At(const std::string& key, const std::string& context = "object");

  // Compact serialization with lexicographically sorted object keys.  The
  // output is stable across runs and is the exact byte sequence hashed by the
  // material, geometry, and case provenance code.
  std::string DumpCanonical() const;

private:
  Value value_;
};

// Reject fields that are not part of a versioned schema.  This is central to
// the "no silent default at a scientific boundary" Stage-1 requirement: a
// misspelled field must fail rather than be ignored.
void RejectUnknownFields(const Json::Object& object,
                         const std::vector<std::string>& allowed,
                         const std::string& context);

} // namespace shield_sim

#endif
