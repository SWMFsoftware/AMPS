// ============================================================================
// src/general/input/test_amps_input.cpp
//
// Dependency-free regression test of the generic AMPS runtime input grammar
// (amps_input.h, specified in src/INPUT_FILE.md).  Every fixture is written
// by hand from the specification; no expected value is produced by the reader
// itself.  Build and run (from the AMPS root):
//
//   make -C src/general/input test
//
// Output: one "[INP-xx] PASS|FAIL message" line per test (the format read by
// tools/sep_test_orchestrator.py), then a summary; exit status 1 on failure.
// ============================================================================

#include "amps_input.h"

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <random>
#include <string>
#include <vector>

namespace {

namespace In = AMPS::Input;
namespace fs = std::filesystem;

struct Outcome {
  bool pass;
  std::string message;
};

Outcome Pass(const std::string& message = "ok") { return {true, message}; }
Outcome Fail(const std::string& message) { return {false, message}; }

// Parse text that must fail and check that the diagnostic names the expected
// physical line and contains the expected phrase.
Outcome ExpectError(const std::string& text, const std::string& lineTag,
                    const std::string& phrase) {
  In::Document document;
  const In::Status status = In::ParseText(text, "fixture.in", ".", &document);
  if (status.ok()) return Fail("accepted invalid input: " + phrase);
  if (status.message.find(lineTag) == std::string::npos)
    return Fail("diagnostic lacks location '" + lineTag + "': " + status.message);
  if (status.message.find(phrase) == std::string::npos)
    return Fail("diagnostic lacks '" + phrase + "': " + status.message);
  return Pass();
}

// INP-01: nested blocks, named blocks, comments, continuation, values with
// blanks, and locations of every element.
Outcome Grammar() {
  const std::string text =
      "! leading comment\n"
      "#begin mesh            ! comment after a directive\n"
      "  maximum_level = 7\n"
      "  #begin refinement: solar\n"
      "    model = radial-smoothstep\n"
      "    surface_cell_size_m = \\   ! continued\n"
      "      1.5e9\n"
      "  #end refinement\n"
      "#end mesh\n"
      "#BEGIN Physics: swcme\n"
      "  ambient.wind_speed = 400 km/s\n"
      "#End PHYSICS\n";
  In::Document document;
  const In::Status status = In::ParseText(text, "fixture.in", ".", &document);
  if (!status.ok()) return Fail(status.message);
  const In::Block& root = document.root;
  if (root.blocks.size() != 2) return Fail("expected two top-level blocks");
  const In::Block* mesh = root.FindBlock("mesh");
  if (!mesh || mesh->begin.line != 2) return Fail("mesh block or its line");
  const In::Assignment* level = mesh->Find("maximum_level");
  if (!level || level->value != "7" || level->where.line != 3)
    return Fail("maximum_level value/line");
  const In::Block* solar = mesh->FindBlock("refinement", "solar");
  if (!solar || solar->Label() != "refinement: solar") return Fail("named block");
  const In::Assignment* cell = solar->Find("surface_cell_size_m");
  if (!cell || cell->value != "1.5e9" || cell->where.line != 6)
    return Fail("continuation must join to the first physical line");
  const In::Block* swcme = root.FindBlock("physics", "swcme");
  if (!swcme) return Fail("kind is case-insensitive, name case-sensitive");
  const In::Assignment* wind = swcme->Find("ambient.wind_speed");
  if (!wind || wind->value != "400 km/s") return Fail("value with blanks");
  return Pass("nested/named blocks, comments, continuation, locations");
}

// INP-02: structural errors with file:line diagnostics.
Outcome StructuralErrors() {
  struct Case { std::string text, lineTag, phrase; };
  const std::vector<Case> cases = {
      {"key = 1\n", "fixture.in:1:", "outside any block"},
      {"#begin run\nx = 1\nx = 2\n#end run\n", "fixture.in:3:", "duplicate key 'x'"},
      {"#begin run\n#end mesh\n", "fixture.in:2:", "does not match"},
      {"#begin run\nx = 1\n", "fixture.in:1:", "is not closed"},
      {"#end run\n", "fixture.in:1:", "no matching"},
      {"#begin run\n#end\n", "fixture.in:2:", "must name the block kind"},
      {"#begin o: a\n#end o\n#begin o: a\n#end o\n", "fixture.in:3:", "duplicate block 'o: a'"},
      {"#begin run\nx =\n#end run\n", "fixture.in:2:", "has no value"},
      {"#begin run\n1x = 2\n#end run\n", "fixture.in:2:", "invalid key"},
      {"#begin Run-Control\n#end run\n", "fixture.in:1:", "invalid block kind"},
      {"#begin o: a.b\n#end o\n", "fixture.in:1:", "invalid block name"},
      {"#section begin: sep3d\n", "fixture.in:1:", "unknown directive"},
      {"#begin run\nno equals sign\n#end run\n", "fixture.in:2:", "expected 'key = value'"},
      {"#begin run\nx = \\\n", "fixture.in:2:", "continuation"},
  };
  for (const Case& item : cases) {
    const Outcome outcome = ExpectError(item.text, item.lineTag, item.phrase);
    if (!outcome.pass) return outcome;
  }
  return Pass(std::to_string(cases.size()) + " structural errors diagnosed");
}

// INP-03: #include resolution relative to the including file, include inside
// an open block, balanced-file rule, cycles and missing files.
Outcome Includes() {
  const fs::path root = fs::temp_directory_path() / "amps_input_test_includes";
  std::error_code ignored;
  fs::remove_all(root, ignored);
  fs::create_directories(root / "parts");
  auto write = [](const fs::path& path, const std::string& text) {
    std::ofstream(path) << text;
  };
  write(root / "amps.in",
        "#include \"parts/mesh.in\"\n"
        "#begin physics: background\n"
        "  #include <parts/model_keys.in>\n"
        "#end physics\n");
  write(root / "parts/mesh.in", "#begin mesh\n  maximum_level = 4\n#end mesh\n");
  write(root / "parts/model_keys.in", "model = reduced-shock-front\n");
  In::Document document;
  In::Status status = In::ParseFile((root / "amps.in").string(), &document);
  if (!status.ok()) return Fail(status.message);
  if (document.files.size() != 3) return Fail("expected 3 files in provenance");
  const In::Block* background = document.root.FindBlock("physics", "background");
  if (!background || !background->Find("model"))
    return Fail("include inside an open block must supply its assignments");
  if (background->Find("model")->where.file.find("model_keys.in") == std::string::npos)
    return Fail("assignment must keep the included file as its location");

  write(root / "unbalanced.in", "#begin mesh\n#include \"parts/close.in\"\n");
  write(root / "parts/close.in", "#end mesh\n");
  status = In::ParseFile((root / "unbalanced.in").string(), &document);
  if (status.ok() || status.message.find("no matching") == std::string::npos)
    return Fail("an included file must not close an outer block");

  write(root / "a.in", "#include \"b.in\"\n");
  write(root / "b.in", "#include \"a.in\"\n");
  status = In::ParseFile((root / "a.in").string(), &document);
  if (status.ok() || status.message.find("cycle") == std::string::npos)
    return Fail("include cycle not diagnosed");

  write(root / "missing.in", "#include \"nowhere.in\"\n");
  status = In::ParseFile((root / "missing.in").string(), &document);
  if (status.ok() || status.message.find("missing.in:1:") == std::string::npos)
    return Fail("missing include must be reported at the #include line");
  fs::remove_all(root, ignored);
  return Pass("relative/angled/quoted includes, balance, cycle, missing file");
}

// INP-04: canonical text round trip, flattening and ownership helpers.
Outcome CanonicalAndOwnership() {
  const std::string text =
      "#begin physics: shock_front\n"
      "  schema = v1\n"
      "  #begin ambient\n    density = 7\n  #end ambient\n"
      "  #begin observer: earth\n    radius = 1\n  #end observer\n"
      "#end physics\n"
      "#begin photolysis\n#end photolysis\n";
  In::Document document;
  In::Status status = In::ParseText(text, "fixture.in", ".", &document);
  if (!status.ok()) return Fail(status.message);
  const std::string canonical = In::ToCanonicalText(document.root);
  const std::string expected =
      "#begin physics: shock_front\n"
      "  schema = v1\n"
      "  #begin ambient\n    density = 7\n  #end ambient\n"
      "  #begin observer: earth\n    radius = 1\n  #end observer\n"
      "#end physics\n"
      "#begin photolysis\n#end photolysis\n";
  if (canonical != expected) return Fail("canonical text:\n" + canonical);
  In::Document again;
  status = In::ParseText(canonical, "canonical", ".", &again);
  if (!status.ok() || In::ToCanonicalText(again.root) != canonical)
    return Fail("canonical text does not round-trip");

  const std::vector<In::Assignment> flat =
      In::Flatten(*document.root.FindBlock("physics", "shock_front"));
  if (flat.size() != 3 || flat[0].key != "schema" ||
      flat[1].key != "ambient.density" || flat[2].key != "observer.earth.radius")
    return Fail("flattened keys");
  if (flat[1].where.line != 4) return Fail("flattening must keep locations");

  status = In::RequireOnlyBlocks(document.root, {"physics", "mesh"}, "this test");
  if (status.ok() || status.message.find("fixture.in:10:") == std::string::npos ||
      status.message.find("'photolysis'") == std::string::npos)
    return Fail("unknown block must be fatal at its #begin line");
  status = In::RequireNoAssignments(*document.root.FindBlock("physics"));
  if (status.ok()) return Fail("RequireNoAssignments accepted a key");
  return Pass("canonical round trip, flattening, unknown-block policy");
}

}  // namespace

int main(int argc, char** argv) {
  struct Test {
    std::string id;
    const char* description;
    std::function<Outcome()> function;
  };
  std::vector<Test> tests = {
      {"INP-01", "nested/named blocks, comments, continuations, and locations",
       Grammar},
      {"INP-02", "structural-error diagnostics and physical line numbers",
       StructuralErrors},
      {"INP-03", "relative includes, balance, cycles, and missing files",
       Includes},
      {"INP-04", "canonical round trip, flattening, and ownership policy",
       CanonicalAndOwnership},
  };

  // Keep this dependency-light program's selection contract aligned with the
  // shared C++ harness used by configured AMPS tests.  Exact repeated --case
  // arguments let the global runner prove that every case is independently
  // executable.  Deterministic reverse/shuffle modes expose hidden mutable
  // state and order dependencies; the default remains the historical order.
  std::vector<std::string> selected_ids;
  std::string filter;
  std::string order = "normal";
  unsigned int shuffle_seed = 0;
  for (int i = 1; i < argc; ++i) {
    const std::string argument = argv[i];
    if (argument == "--case" && i + 1 < argc) {
      selected_ids.push_back(argv[++i]);
    } else if (argument == "--filter" && i + 1 < argc) {
      filter = argv[++i];
    } else if (argument == "--order" && i + 1 < argc) {
      order = argv[++i];
    } else if (argument == "--seed" && i + 1 < argc) {
      try {
        shuffle_seed = static_cast<unsigned int>(std::stoul(argv[++i]));
      } catch (const std::exception&) {
        std::fprintf(stderr, "invalid --seed value\n");
        return 2;
      }
    } else {
      std::fprintf(stderr,
                   "usage: %s [--case ID]... [--filter SUBSTRING] "
                   "[--order normal|reverse|shuffle] [--seed N]\n",
                   argv[0]);
      return 2;
    }
  }
  if (order != "normal" && order != "reverse" && order != "shuffle") {
    std::fprintf(stderr, "unknown test order '%s'\n", order.c_str());
    return 2;
  }
  if (!selected_ids.empty() && !filter.empty()) {
    std::fprintf(stderr, "--case and --filter are mutually exclusive\n");
    return 2;
  }

  // Reject a misspelled exact ID even when another selected case would run;
  // silently accepting a partial selection would create false coverage.
  for (const std::string& requested : selected_ids) {
    const bool known = std::any_of(tests.begin(), tests.end(),
        [&requested](const Test& test) { return test.id == requested; });
    if (!known) {
      std::fprintf(stderr, "unknown input-test case '%s'\n", requested.c_str());
      return 2;
    }
  }

  // Exact selection begins in repeated --case order, matching the shared
  // harness contract.  This is important because the global runner has
  // already materialized reverse/shuffle order before invoking the adapter;
  // filtering the fixed registry order would silently undo that evidence.
  std::vector<Test> execution;
  if (!selected_ids.empty()) {
    for (const std::string& requested : selected_ids) {
      const auto match = std::find_if(tests.begin(), tests.end(),
          [&requested](const Test& test) { return test.id == requested; });
      execution.push_back(*match);
    }
  } else {
    for (const Test& test : tests)
      if (filter.empty() || test.id.find(filter) != std::string::npos)
        execution.push_back(test);
  }
  if (order == "reverse") {
    std::reverse(execution.begin(), execution.end());
  } else if (order == "shuffle") {
    std::mt19937 generator(shuffle_seed);
    std::shuffle(execution.begin(), execution.end(), generator);
  }

  int failed = 0;
  for (const auto& test : execution) {
    std::printf("[CASE] RUN  %s %s\n", test.id.c_str(), test.description);
    std::fflush(stdout);
    const Outcome outcome = test.function();
    std::printf("[%s] %s %s\n", test.id.c_str(), outcome.pass ? "PASS" : "FAIL",
                outcome.message.c_str());
    if (!outcome.pass) ++failed;
  }
  if (execution.empty()) {
    std::fprintf(stderr, "no input tests matched the requested selection\n");
    return 2;
  }
  std::printf("amps_input tests: total=%zu pass=%zu fail=%d\n", execution.size(),
              execution.size() - static_cast<std::size_t>(failed), failed);
  return failed == 0 ? 0 : 1;
}
