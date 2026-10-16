#include "StrictJson.hh"

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>

namespace shield_sim {
namespace {

class Parser {
public:
  Parser(const std::string& text, std::string source)
    : text_(text), source_(std::move(source)) {}

  Json ParseDocument() {
    SkipWhitespace();
    Json result = ParseValue();
    SkipWhitespace();
    if (position_ != text_.size()) Fail("unexpected trailing token");
    return result;
  }

private:
  [[noreturn]] void Fail(const std::string& message) const {
    std::size_t line = 1;
    std::size_t column = 1;
    for (std::size_t i = 0; i < position_ && i < text_.size(); ++i) {
      if (text_[i] == '\n') { ++line; column = 1; }
      else { ++column; }
    }
    std::ostringstream out;
    out << source_ << ':' << line << ':' << column << ": " << message;
    throw JsonError(out.str());
  }

  void SkipWhitespace() {
    while (position_ < text_.size()) {
      const char c = text_[position_];
      if (c == ' ' || c == '\t' || c == '\n' || c == '\r') ++position_;
      else break;
    }
  }

  char Peek() const {
    if (position_ >= text_.size()) Fail("unexpected end of input");
    return text_[position_];
  }

  bool Consume(char wanted) {
    if (position_ < text_.size() && text_[position_] == wanted) {
      ++position_;
      return true;
    }
    return false;
  }

  void Require(char wanted, const std::string& description) {
    if (!Consume(wanted)) Fail("expected " + description);
  }

  Json ParseValue() {
    SkipWhitespace();
    const char c = Peek();
    if (c == '{') return ParseObject();
    if (c == '[') return ParseArray();
    if (c == '"') return Json(ParseString());
    if (c == 't') return ParseLiteral("true", Json(true));
    if (c == 'f') return ParseLiteral("false", Json(false));
    if (c == 'n') return ParseLiteral("null", Json(nullptr));
    if (c == '-' || (c >= '0' && c <= '9')) return ParseNumber();
    Fail("expected a JSON value");
  }

  Json ParseLiteral(const char* literal, Json result) {
    for (const char* p = literal; *p; ++p) {
      if (position_ >= text_.size() || text_[position_] != *p)
        Fail(std::string("expected '") + literal + "'");
      ++position_;
    }
    return result;
  }

  static void AppendUtf8(std::string& out, std::uint32_t codepoint) {
    if (codepoint <= 0x7f) {
      out.push_back(static_cast<char>(codepoint));
    } else if (codepoint <= 0x7ff) {
      out.push_back(static_cast<char>(0xc0 | (codepoint >> 6)));
      out.push_back(static_cast<char>(0x80 | (codepoint & 0x3f)));
    } else if (codepoint <= 0xffff) {
      out.push_back(static_cast<char>(0xe0 | (codepoint >> 12)));
      out.push_back(static_cast<char>(0x80 | ((codepoint >> 6) & 0x3f)));
      out.push_back(static_cast<char>(0x80 | (codepoint & 0x3f)));
    } else {
      out.push_back(static_cast<char>(0xf0 | (codepoint >> 18)));
      out.push_back(static_cast<char>(0x80 | ((codepoint >> 12) & 0x3f)));
      out.push_back(static_cast<char>(0x80 | ((codepoint >> 6) & 0x3f)));
      out.push_back(static_cast<char>(0x80 | (codepoint & 0x3f)));
    }
  }

  std::uint32_t ParseHex4() {
    std::uint32_t result = 0;
    for (int i = 0; i < 4; ++i) {
      if (position_ >= text_.size()) Fail("incomplete Unicode escape");
      const char c = text_[position_++];
      result <<= 4;
      if (c >= '0' && c <= '9') result += c - '0';
      else if (c >= 'a' && c <= 'f') result += 10 + c - 'a';
      else if (c >= 'A' && c <= 'F') result += 10 + c - 'A';
      else Fail("invalid hexadecimal digit in Unicode escape");
    }
    return result;
  }

  std::string ParseString() {
    Require('"', "a string");
    std::string out;
    while (position_ < text_.size()) {
      const unsigned char c = static_cast<unsigned char>(text_[position_++]);
      if (c == '"') return out;
      if (c < 0x20) Fail("unescaped control character in string");
      if (c != '\\') {
        out.push_back(static_cast<char>(c));
        continue;
      }
      if (position_ >= text_.size()) Fail("incomplete string escape");
      const char e = text_[position_++];
      switch (e) {
        case '"': out.push_back('"'); break;
        case '\\': out.push_back('\\'); break;
        case '/':  out.push_back('/');  break;
        case 'b':  out.push_back('\b'); break;
        case 'f':  out.push_back('\f'); break;
        case 'n':  out.push_back('\n'); break;
        case 'r':  out.push_back('\r'); break;
        case 't':  out.push_back('\t'); break;
        case 'u': {
          std::uint32_t cp = ParseHex4();
          if (cp >= 0xd800 && cp <= 0xdbff) {
            if (position_ + 2 > text_.size() || text_[position_] != '\\' ||
                text_[position_ + 1] != 'u')
              Fail("high UTF-16 surrogate is not followed by a low surrogate");
            position_ += 2;
            const std::uint32_t low = ParseHex4();
            if (low < 0xdc00 || low > 0xdfff)
              Fail("invalid low UTF-16 surrogate");
            cp = 0x10000 + ((cp - 0xd800) << 10) + (low - 0xdc00);
          } else if (cp >= 0xdc00 && cp <= 0xdfff) {
            Fail("unpaired low UTF-16 surrogate");
          }
          AppendUtf8(out, cp);
          break;
        }
        default: Fail("invalid string escape");
      }
    }
    Fail("unterminated string");
  }

  Json ParseNumber() {
    const std::size_t begin = position_;
    Consume('-');
    if (Consume('0')) {
      if (position_ < text_.size() && text_[position_] >= '0' && text_[position_] <= '9')
        Fail("leading zero in number");
    } else {
      if (position_ >= text_.size() || text_[position_] < '1' || text_[position_] > '9')
        Fail("invalid number");
      while (position_ < text_.size() && text_[position_] >= '0' && text_[position_] <= '9')
        ++position_;
    }
    if (Consume('.')) {
      if (position_ >= text_.size() || text_[position_] < '0' || text_[position_] > '9')
        Fail("fraction requires at least one digit");
      while (position_ < text_.size() && text_[position_] >= '0' && text_[position_] <= '9')
        ++position_;
    }
    if (position_ < text_.size() && (text_[position_] == 'e' || text_[position_] == 'E')) {
      ++position_;
      if (position_ < text_.size() && (text_[position_] == '+' || text_[position_] == '-'))
        ++position_;
      if (position_ >= text_.size() || text_[position_] < '0' || text_[position_] > '9')
        Fail("exponent requires at least one digit");
      while (position_ < text_.size() && text_[position_] >= '0' && text_[position_] <= '9')
        ++position_;
    }
    const std::string token = text_.substr(begin, position_ - begin);
    std::size_t used = 0;
    double value = 0.0;
    try { value = std::stod(token, &used); }
    catch (...) { Fail("number is outside the supported double range"); }
    if (used != token.size() || !std::isfinite(value))
      Fail("number must be finite");
    return Json(value);
  }

  Json ParseArray() {
    Require('[', "'['");
    Json::Array array;
    SkipWhitespace();
    if (Consume(']')) return Json(std::move(array));
    for (;;) {
      array.push_back(ParseValue());
      SkipWhitespace();
      if (Consume(']')) break;
      Require(',', "',' or ']'");
      SkipWhitespace();
    }
    return Json(std::move(array));
  }

  Json ParseObject() {
    Require('{', "'{'");
    Json::Object object;
    SkipWhitespace();
    if (Consume('}')) return Json(std::move(object));
    for (;;) {
      if (Peek() != '"') Fail("object key must be a string");
      const std::string key = ParseString();
      SkipWhitespace();
      Require(':', "':' after object key");
      SkipWhitespace();
      Json value = ParseValue();
      if (!object.emplace(key, std::move(value)).second)
        Fail("duplicate object key '" + key + "'");
      SkipWhitespace();
      if (Consume('}')) break;
      Require(',', "',' or '}'");
      SkipWhitespace();
    }
    return Json(std::move(object));
  }

  const std::string& text_;
  std::string source_;
  std::size_t position_ = 0;
};

std::string Escape(const std::string& value) {
  std::ostringstream out;
  out << '"';
  for (const unsigned char c : value) {
    switch (c) {
      case '"': out << "\\\""; break;
      case '\\': out << "\\\\"; break;
      case '\b': out << "\\b"; break;
      case '\f': out << "\\f"; break;
      case '\n': out << "\\n"; break;
      case '\r': out << "\\r"; break;
      case '\t': out << "\\t"; break;
      default:
        if (c < 0x20) {
          out << "\\u" << std::hex << std::setw(4) << std::setfill('0')
              << static_cast<unsigned>(c) << std::dec;
        } else {
          out << static_cast<char>(c);
        }
    }
  }
  out << '"';
  return out.str();
}

std::string Dump(const Json& value) {
  if (value.IsNull()) return "null";
  if (value.IsBool()) return value.AsBool() ? "true" : "false";
  if (value.IsNumber()) {
    std::ostringstream out;
    out << std::setprecision(std::numeric_limits<double>::max_digits10)
        << value.AsNumber();
    return out.str();
  }
  if (value.IsString()) return Escape(value.AsString());
  if (value.IsArray()) {
    std::string out = "[";
    const auto& array = value.AsArray();
    for (std::size_t i = 0; i < array.size(); ++i) {
      if (i) out += ',';
      out += Dump(array[i]);
    }
    return out + ']';
  }
  std::string out = "{";
  bool first = true;
  for (const auto& pair : value.AsObject()) {
    if (!first) out += ',';
    first = false;
    out += Escape(pair.first) + ':' + Dump(pair.second);
  }
  return out + '}';
}

[[noreturn]] void TypeError(const std::string& context, const std::string& wanted) {
  throw JsonError(context + " must be " + wanted);
}

} // namespace

Json Json::Parse(const std::string& text, const std::string& source) {
  return Parser(text, source).ParseDocument();
}

Json Json::ParseFile(const std::string& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) throw JsonError("cannot open JSON file '" + path + "'");
  std::ostringstream text;
  text << input.rdbuf();
  if (!input.good() && !input.eof()) throw JsonError("cannot read JSON file '" + path + "'");
  return Parse(text.str(), path);
}

bool Json::IsNull() const { return std::holds_alternative<std::nullptr_t>(value_); }
bool Json::IsBool() const { return std::holds_alternative<bool>(value_); }
bool Json::IsNumber() const { return std::holds_alternative<double>(value_); }
bool Json::IsString() const { return std::holds_alternative<std::string>(value_); }
bool Json::IsArray() const { return std::holds_alternative<Array>(value_); }
bool Json::IsObject() const { return std::holds_alternative<Object>(value_); }

bool Json::AsBool(const std::string& context) const {
  if (!IsBool()) TypeError(context, "a boolean");
  return std::get<bool>(value_);
}
double Json::AsNumber(const std::string& context) const {
  if (!IsNumber()) TypeError(context, "a number");
  return std::get<double>(value_);
}
const std::string& Json::AsString(const std::string& context) const {
  if (!IsString()) TypeError(context, "a string");
  return std::get<std::string>(value_);
}
const Json::Array& Json::AsArray(const std::string& context) const {
  if (!IsArray()) TypeError(context, "an array");
  return std::get<Array>(value_);
}
const Json::Object& Json::AsObject(const std::string& context) const {
  if (!IsObject()) TypeError(context, "an object");
  return std::get<Object>(value_);
}
Json::Array& Json::AsArray(const std::string& context) {
  if (!IsArray()) TypeError(context, "an array");
  return std::get<Array>(value_);
}
Json::Object& Json::AsObject(const std::string& context) {
  if (!IsObject()) TypeError(context, "an object");
  return std::get<Object>(value_);
}

bool Json::Contains(const std::string& key) const {
  if (!IsObject()) return false;
  return std::get<Object>(value_).count(key) != 0;
}
const Json& Json::At(const std::string& key, const std::string& context) const {
  const auto& object = AsObject(context);
  const auto found = object.find(key);
  if (found == object.end()) throw JsonError(context + " is missing required field '" + key + "'");
  return found->second;
}
Json& Json::At(const std::string& key, const std::string& context) {
  auto& object = AsObject(context);
  const auto found = object.find(key);
  if (found == object.end()) throw JsonError(context + " is missing required field '" + key + "'");
  return found->second;
}

std::string Json::DumpCanonical() const { return Dump(*this); }

void RejectUnknownFields(const Json::Object& object,
                         const std::vector<std::string>& allowed,
                         const std::string& context) {
  const std::set<std::string> known(allowed.begin(), allowed.end());
  for (const auto& pair : object) {
    if (!known.count(pair.first))
      throw JsonError(context + " contains unknown field '" + pair.first + "'");
  }
}

} // namespace shield_sim
