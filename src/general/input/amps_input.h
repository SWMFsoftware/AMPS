// ============================================================================
// src/general/input/amps_input.h
//
// Generic AMPS runtime input reader (the grammar of ``amps.in``).
//
// One runtime input format configures every AMPS application: planets, moons,
// comets, spacecraft re-entry, SEP and GCR transport.  The file is a tree of
// named blocks.  This header implements only the *grammar*; it never assigns
// physical meaning to a key.  Each block is handed to exactly one owner (the
// AMPS core, the mesh, a model library under src/models, or an application),
// and that owner validates its own keys, values and units.  The complete
// specification, including the block kinds and their owners, is
// src/INPUT_FILE.md; this header must stay consistent with it.
//
// Grammar (summary; see src/INPUT_FILE.md for the normative text)
// ---------------------------------------------------------------
//   !  starts a comment that runs to the end of the physical line.
//   \  as the last non-blank character before a comment joins the next
//      physical line to this one (one space is inserted between fragments).
//   #begin KIND            opens a block; KIND = [a-z][a-z0-9_]* (case-
//   #begin KIND: NAME      insensitive, stored lower case); NAME =
//                          [A-Za-z0-9_-]+ (case-sensitive) identifies one
//                          instance, e.g. ``#begin observer: earth``.
//   #end KIND              closes the innermost open block, which must have
//                          the same KIND.
//   #include PATH          inserts another file here; PATH may be bare,
//   #include "PATH"        "quoted" or <angled>; a relative path is resolved
//   #include <PATH>        from the directory of the *including* file.
//   KEY = VALUE            one assignment inside a block.  KEY =
//                          [A-Za-z][A-Za-z0-9_.-]* (case-sensitive; dots are
//                          allowed so model libraries can keep dotted
//                          vocabularies); VALUE is the trimmed remainder and
//                          must not be empty.  Values may contain blanks
//                          ("400 km/s", "0, 0, 1").
//
// Structural rules enforced here (every violation is fatal and reported as
// ``file:line: reason`` followed by the offending logical line):
//   * an assignment must be inside a block;
//   * a key may appear only once in one block;
//   * two sibling blocks may not share the same KIND and NAME;
//   * every block opened in a file must be closed in the same file (an
//     included file is therefore always block-balanced; it may still be
//     included inside an open block to supply that block's assignments);
//   * unknown ``#`` directives, include cycles, includes nested deeper than
//     kMaximumIncludeDepth, missing files and dangling continuations fail.
//
// Ownership checks (unknown block kinds, unknown keys, required keys, value
// domains) belong to the block owners.  ``RequireOnlyBlocks`` is the common
// helper with which an owner rejects child blocks it does not understand; by
// the AMPS policy an unknown block is always an error, never ignored.
//
// Design constraints
//   * Header-only, C++17, standard library only: model libraries under
//     src/models must stay usable without PIC, MPI or application headers,
//     and they add only ``-I<AMPS>/src/general/input`` to use this reader.
//   * Deterministic: the canonical text produced by ``ToCanonicalText`` is a
//     pure function of the parsed tree (used for MPI broadcast, provenance
//     and configuration fingerprints).
//   * Every assignment keeps its source file and physical line, so an owner
//     can report a domain error at the exact place the user wrote it.
// ============================================================================

#ifndef AMPS_GENERAL_INPUT_AMPS_INPUT_H
#define AMPS_GENERAL_INPUT_AMPS_INPUT_H

#include <cctype>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace AMPS {
namespace Input {

// Maximum nesting of #include directives.  The bound turns a pathological
// (but acyclic) include chain into a diagnosed error instead of unbounded
// recursion; real inputs use one or two levels.
constexpr std::size_t kMaximumIncludeDepth = 64;

// ----------------------------------------------------------------------------
// Result type
// ----------------------------------------------------------------------------

// Outcome of a reader or owner operation.  ``message`` already contains the
// location prefix ``file:line:`` when the error refers to user input.
struct Status {
  bool ok_ = true;
  std::string message;

  bool ok() const { return ok_; }
  static Status OK() { return Status(); }
  static Status Error(std::string text) {
    Status status;
    status.ok_ = false;
    status.message = std::move(text);
    return status;
  }
};

// ----------------------------------------------------------------------------
// Parsed tree
// ----------------------------------------------------------------------------

// Where a block or an assignment was written.  ``line`` is the 1-based
// physical line of the first fragment of the logical line.
struct Location {
  std::string file;
  std::size_t line = 0;

  std::string Prefix() const { return file + ":" + std::to_string(line) + ": "; }
};

// One ``KEY = VALUE`` statement.  ``text`` is the logical line after comment
// removal and continuation joining; owners quote it in diagnostics.
struct Assignment {
  std::string key;
  std::string value;
  Location where;
  std::string text;
};

// One ``#begin KIND[: NAME] ... #end KIND`` block.  The document root is a
// Block with an empty kind; it holds only child blocks.
struct Block {
  std::string kind;
  std::string name;
  Location begin;
  std::vector<Assignment> assignments;
  std::vector<Block> blocks;

  // "kind" or "kind: name" exactly as a user would write it after #begin.
  std::string Label() const {
    if (kind.empty()) return "<document>";
    return name.empty() ? kind : kind + ": " + name;
  }

  // Assignment with this key, or nullptr.  Keys are unique within a block.
  const Assignment* Find(const std::string& key) const {
    for (const Assignment& item : assignments)
      if (item.key == key) return &item;
    return nullptr;
  }

  // First child block with this kind (and name, when ``name`` is given).
  const Block* FindBlock(const std::string& childKind,
                         const std::string& childName = std::string()) const {
    for (const Block& child : blocks)
      if (child.kind == childKind &&
          (childName.empty() || child.name == childName))
        return &child;
    return nullptr;
  }

  // Every child block of one kind, in file order (e.g. all observers).
  std::vector<const Block*> BlocksOfKind(const std::string& childKind) const {
    std::vector<const Block*> result;
    for (const Block& child : blocks)
      if (child.kind == childKind) result.push_back(&child);
    return result;
  }
};

// A parsed input: the root block plus every file that was read (root first,
// then includes in the order they were expanded) for provenance receipts.
struct Document {
  Block root;
  std::vector<std::string> files;
};

// Diagnostic for an error at a user-visible location, in the uniform
// ``file:line: reason`` form, optionally followed by the offending line.
inline Status ErrorAt(const Location& where, const std::string& reason,
                      const std::string& text = std::string()) {
  std::string message = where.Prefix() + reason;
  if (!text.empty()) message += "\n  line: " + text;
  return Status::Error(message);
}

// ----------------------------------------------------------------------------
// Implementation details (lexing and recursive file expansion)
// ----------------------------------------------------------------------------
namespace detail {

inline std::string Trim(const std::string& text) {
  std::size_t first = 0;
  while (first < text.size() &&
         std::isspace(static_cast<unsigned char>(text[first])))
    ++first;
  std::size_t last = text.size();
  while (last > first &&
         std::isspace(static_cast<unsigned char>(text[last - 1])))
    --last;
  return text.substr(first, last - first);
}

inline std::string Lower(std::string text) {
  for (char& c : text)
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return text;
}

// KIND: lower-case letter followed by lower-case letters, digits, '_'.
inline bool ValidKind(const std::string& text) {
  if (text.empty() || !std::islower(static_cast<unsigned char>(text[0])))
    return false;
  for (char c : text)
    if (!(std::islower(static_cast<unsigned char>(c)) ||
          std::isdigit(static_cast<unsigned char>(c)) || c == '_'))
      return false;
  return true;
}

// NAME: letters, digits, '_' and '-'.  Dots are excluded so that the dotted
// qualified name of a key (see QualifiedName) is unambiguous.
inline bool ValidName(const std::string& text) {
  if (text.empty()) return false;
  for (char c : text)
    if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-'))
      return false;
  return true;
}

// KEY: a letter followed by letters, digits, '_', '.', '-'.
inline bool ValidKey(const std::string& text) {
  if (text.empty() || !std::isalpha(static_cast<unsigned char>(text[0])))
    return false;
  for (char c : text)
    if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_' ||
          c == '.' || c == '-'))
      return false;
  return true;
}

// One logical line (comments removed, continuations joined) with the
// location of its first physical fragment.
struct LogicalLine {
  Location where;
  std::string text;
};

// Split raw file text into logical lines.  A '!' starts a comment on every
// physical line; a trailing '\' (after comment removal) continues the logical
// line.  A continuation still open at end of text is an error.
inline Status SplitLogicalLines(const std::string& text,
                                const std::string& file,
                                std::vector<LogicalLine>* lines) {
  std::istringstream input(text);
  std::string physical;
  std::size_t number = 0;
  bool continuing = false;
  LogicalLine pending;
  while (std::getline(input, physical)) {
    ++number;
    if (!physical.empty() && physical.back() == '\r') physical.pop_back();
    const std::size_t comment = physical.find('!');
    if (comment != std::string::npos) physical.erase(comment);
    std::string fragment = Trim(physical);
    bool continues = false;
    if (!fragment.empty() && fragment.back() == '\\') {
      continues = true;
      fragment = Trim(fragment.substr(0, fragment.size() - 1));
    }
    if (!continuing) {
      pending = LogicalLine();
      pending.where = Location{file, number};
    }
    if (!fragment.empty()) {
      if (!pending.text.empty()) pending.text += ' ';
      pending.text += fragment;
    }
    continuing = continues;
    if (!continuing && !pending.text.empty()) lines->push_back(pending);
  }
  if (continuing)
    return ErrorAt(pending.where,
                   "line continuation '\\' at end of file has no next line",
                   pending.text);
  return Status::OK();
}

inline Status ReadWholeFile(const std::filesystem::path& path,
                            std::string* text) {
  std::ifstream stream(path, std::ios::binary);
  if (!stream) return Status::Error("cannot open input file: " + path.string());
  std::ostringstream buffer;
  buffer << stream.rdbuf();
  if (stream.bad())
    return Status::Error("cannot read input file: " + path.string());
  *text = buffer.str();
  return Status::OK();
}

// Extract the include target from the text after "#include".
inline Status IncludeTarget(const LogicalLine& line, const std::string& rest,
                            std::string* target) {
  std::string body = Trim(rest);
  if (body.empty())
    return ErrorAt(line.where, "#include requires a file name", line.text);
  if (body.front() == '"' || body.front() == '<') {
    const char close = body.front() == '"' ? '"' : '>';
    const std::size_t end = body.find(close, 1);
    if (end == std::string::npos)
      return ErrorAt(line.where, "unterminated #include file name", line.text);
    if (!Trim(body.substr(end + 1)).empty())
      return ErrorAt(line.where, "unexpected text after the #include file name",
                     line.text);
    *target = body.substr(1, end - 1);
  } else {
    *target = body;
  }
  if (target->empty())
    return ErrorAt(line.where, "#include requires a file name", line.text);
  return Status::OK();
}

// Recursive parser state.  ``stack`` holds the currently open blocks (root
// first); ``includeChain`` the canonical paths of the files being expanded,
// for cycle detection.
struct Parser {
  Document* document = nullptr;
  std::vector<Block*> stack;
  std::vector<std::string> includeChain;

  Status ParseText(const std::string& text, const std::string& displayName,
                   const std::filesystem::path& baseDirectory) {
    std::vector<LogicalLine> lines;
    Status status = SplitLogicalLines(text, displayName, &lines);
    if (!status.ok()) return status;
    // Blocks opened in this file must be closed in this file.
    const std::size_t entryDepth = stack.size();
    for (const LogicalLine& line : lines) {
      status = ParseLine(line, baseDirectory, entryDepth);
      if (!status.ok()) return status;
    }
    if (stack.size() != entryDepth) {
      const Block* open = stack.back();
      return ErrorAt(open->begin, "block '" + open->Label() +
                     "' is not closed by '#end " + open->kind +
                     "' before the end of this file");
    }
    return Status::OK();
  }

  Status ParseFile(const std::filesystem::path& path, const LogicalLine* from) {
    std::error_code error;
    std::filesystem::path canonical = std::filesystem::weakly_canonical(path, error);
    if (error) canonical = path;
    const std::string key = canonical.string();
    for (const std::string& open : includeChain)
      if (open == key)
        return from ? ErrorAt(from->where, "#include cycle through " + key,
                              from->text)
                    : Status::Error("#include cycle through " + key);
    if (includeChain.size() >= kMaximumIncludeDepth)
      return from ? ErrorAt(from->where, "#include nesting exceeds " +
                            std::to_string(kMaximumIncludeDepth) + " files",
                            from->text)
                  : Status::Error("#include nesting is too deep");
    std::string text;
    Status status = ReadWholeFile(path, &text);
    if (!status.ok())
      return from ? ErrorAt(from->where, status.message, from->text) : status;
    document->files.push_back(path.string());
    includeChain.push_back(key);
    status = ParseText(text, path.string(), path.parent_path());
    includeChain.pop_back();
    return status;
  }

  Status ParseLine(const LogicalLine& line,
                   const std::filesystem::path& baseDirectory,
                   std::size_t entryDepth) {
    const std::string& text = line.text;
    if (text.front() == '#') return ParseDirective(line, baseDirectory, entryDepth);

    const std::size_t equals = text.find('=');
    if (equals == std::string::npos)
      return ErrorAt(line.where, "expected 'key = value', '#begin', '#end' "
                     "or '#include'", text);
    const std::string key = Trim(text.substr(0, equals));
    const std::string value = Trim(text.substr(equals + 1));
    if (!ValidKey(key))
      return ErrorAt(line.where, "invalid key '" + key + "' (a key starts with "
                     "a letter and contains letters, digits, '_', '.', '-')",
                     text);
    if (value.empty())
      return ErrorAt(line.where, "key '" + key + "' has no value", text);
    Block* owner = stack.back();
    if (owner->kind.empty())
      return ErrorAt(line.where, "assignment '" + key + "' is outside any "
                     "block; place it between '#begin KIND' and '#end KIND'",
                     text);
    if (const Assignment* previous = owner->Find(key))
      return ErrorAt(line.where, "duplicate key '" + key + "' in block '" +
                     owner->Label() + "' (first set at " +
                     previous->where.file + ":" +
                     std::to_string(previous->where.line) + ")", text);
    owner->assignments.push_back(Assignment{key, value, line.where, text});
    return Status::OK();
  }

  Status ParseDirective(const LogicalLine& line,
                        const std::filesystem::path& baseDirectory,
                        std::size_t entryDepth) {
    const std::string& text = line.text;
    std::size_t end = 1;
    while (end < text.size() &&
           std::isalpha(static_cast<unsigned char>(text[end])))
      ++end;
    const std::string directive = Lower(text.substr(1, end - 1));
    const std::string rest = text.substr(end);

    if (directive == "include") {
      std::string target;
      Status status = IncludeTarget(line, rest, &target);
      if (!status.ok()) return status;
      std::filesystem::path path(target);
      if (path.is_relative()) path = baseDirectory / path;
      return ParseFile(path, &line);
    }

    if (directive == "begin") {
      std::string body = Trim(rest);
      std::string kind = body, name;
      const std::size_t colon = body.find(':');
      if (colon != std::string::npos) {
        kind = Trim(body.substr(0, colon));
        name = Trim(body.substr(colon + 1));
        if (!ValidName(name))
          return ErrorAt(line.where, "invalid block name '" + name +
                         "' (letters, digits, '_' and '-' only)", text);
      }
      kind = Lower(kind);
      if (!ValidKind(kind))
        return ErrorAt(line.where, "invalid block kind '" + kind +
                       "' (a lower-case letter followed by letters, digits "
                       "and '_')", text);
      Block* parent = stack.back();
      for (const Block& sibling : parent->blocks)
        if (sibling.kind == kind && sibling.name == name)
          return ErrorAt(line.where, "duplicate block '" + sibling.Label() +
                         "' in '" + parent->Label() + "' (first opened at " +
                         sibling.begin.file + ":" +
                         std::to_string(sibling.begin.line) + ")", text);
      Block child;
      child.kind = kind;
      child.name = name;
      child.begin = line.where;
      parent->blocks.push_back(std::move(child));
      stack.push_back(&parent->blocks.back());
      return Status::OK();
    }

    if (directive == "end") {
      const std::string kind = Lower(Trim(rest));
      if (kind.empty())
        return ErrorAt(line.where, "'#end' must name the block kind it "
                       "closes, e.g. '#end mesh'", text);
      if (stack.size() <= entryDepth || stack.back()->kind.empty())
        return ErrorAt(line.where, "'#end " + kind + "' has no matching "
                       "'#begin " + kind + "' in this file", text);
      if (stack.back()->kind != kind)
        return ErrorAt(line.where, "'#end " + kind + "' does not match the "
                       "open block '" + stack.back()->Label() + "' (opened at " +
                       stack.back()->begin.file + ":" +
                       std::to_string(stack.back()->begin.line) + ")", text);
      stack.pop_back();
      return Status::OK();
    }

    return ErrorAt(line.where, "unknown directive '#" + directive +
                   "' (only #begin, #end and #include exist; comments start "
                   "with '!')", text);
  }
};

inline void AppendCanonical(const Block& block, std::size_t depth,
                            std::string* out) {
  const std::string indent(2 * depth, ' ');
  for (const Assignment& item : block.assignments)
    *out += indent + item.key + " = " + item.value + '\n';
  for (const Block& child : block.blocks) {
    *out += indent + "#begin " + child.Label() + '\n';
    AppendCanonical(child, depth + 1, out);
    *out += indent + "#end " + child.kind + '\n';
  }
}

}  // namespace detail

// ----------------------------------------------------------------------------
// Public entry points
// ----------------------------------------------------------------------------

// Parse a file and everything it includes.
inline Status ParseFile(const std::string& path, Document* document) {
  if (document == nullptr) return Status::Error("input document output is null");
  Document candidate;
  detail::Parser parser;
  parser.document = &candidate;
  parser.stack.push_back(&candidate.root);
  Status status = parser.ParseFile(std::filesystem::path(path), nullptr);
  if (!status.ok()) return status;
  *document = std::move(candidate);
  return Status::OK();
}

// Parse text held in memory.  ``displayName`` is used in diagnostics and
// ``baseDirectory`` resolves relative #include paths (tests use it to parse
// fixtures without writing files).
inline Status ParseText(const std::string& text, const std::string& displayName,
                        const std::string& baseDirectory, Document* document) {
  if (document == nullptr) return Status::Error("input document output is null");
  Document candidate;
  candidate.files.push_back(displayName);
  detail::Parser parser;
  parser.document = &candidate;
  parser.stack.push_back(&candidate.root);
  Status status = parser.ParseText(text, displayName,
                                   std::filesystem::path(baseDirectory));
  if (!status.ok()) return status;
  *document = std::move(candidate);
  return Status::OK();
}

// Canonical text of a block's contents: assignments first (file order), then
// child blocks, two-space indentation, no comments, includes expanded.  It
// is valid input text again, so ParseText(ToCanonicalText(b)) reproduces the
// same tree (locations aside).  Used for broadcast and fingerprints.
inline std::string ToCanonicalText(const Block& block) {
  std::string out;
  detail::AppendCanonical(block, 0, &out);
  return out;
}

// Flatten a block subtree into dotted keys relative to ``block``: a key ``k``
// in child ``#begin ambient`` becomes ``ambient.k``; a named child
// ``#begin observer: earth`` contributes ``observer.earth.``.  Model
// libraries whose vocabulary is dotted (e.g. the reduced shock front, SWCME)
// receive their block in this form, with locations preserved.
inline void Flatten(const Block& block, const std::string& prefix,
                    std::vector<Assignment>* out) {
  for (const Assignment& item : block.assignments) {
    Assignment copy = item;
    copy.key = prefix + item.key;
    out->push_back(std::move(copy));
  }
  for (const Block& child : block.blocks) {
    std::string childPrefix = prefix + child.kind + ".";
    if (!child.name.empty()) childPrefix += child.name + ".";
    Flatten(child, childPrefix, out);
  }
}

inline std::vector<Assignment> Flatten(const Block& block) {
  std::vector<Assignment> out;
  Flatten(block, std::string(), &out);
  return out;
}

// Owner helper: every child block of ``parent`` must have a kind in
// ``allowedKinds``.  AMPS policy: a block whose owner is not part of this
// executable is an error, never silently ignored.
inline Status RequireOnlyBlocks(const Block& parent,
                                const std::set<std::string>& allowedKinds,
                                const std::string& ownerDescription) {
  for (const Block& child : parent.blocks) {
    if (allowedKinds.count(child.kind) == 0) {
      std::string allowed;
      for (const std::string& kind : allowedKinds)
        allowed += (allowed.empty() ? "" : ", ") + kind;
      return ErrorAt(child.begin, "block '" + child.Label() + "' is not "
                     "handled by " + ownerDescription + " (allowed here: " +
                     (allowed.empty() ? std::string("none") : allowed) + ")");
    }
  }
  return Status::OK();
}

// Owner helper: the block must not contain assignments (container blocks
// such as a document root or a pure grouping block).
inline Status RequireNoAssignments(const Block& block) {
  if (!block.assignments.empty()) {
    const Assignment& first = block.assignments.front();
    return ErrorAt(first.where, "block '" + block.Label() + "' holds only "
                   "blocks; key '" + first.key + "' is not allowed here",
                   first.text);
  }
  return Status::OK();
}

}  // namespace Input
}  // namespace AMPS

#endif  // AMPS_GENERAL_INPUT_AMPS_INPUT_H
