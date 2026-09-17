#ifndef SWCME_TEST_TECPLOT_PARSER_HPP
#define SWCME_TEST_TECPLOT_PARSER_HPP

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace swcme_test::tecplot {

// ParsedZone preserves Tecplot's two data organizations instead of flattening
// them: BLOCK values remain grouped by variable, while POINT values remain
// grouped by row.  OUT01 and OUT07 can therefore verify the layout presented
// to an external reader as well as its aggregate value count.
struct ParsedZone {
  std::vector<std::vector<double>> blocks;
  std::vector<std::vector<double>> rows;
  std::vector<std::array<std::size_t,3>> triangles;
  std::size_t nodes=0;
  std::size_t elements=0;
  std::size_t ni=0;
  std::size_t nj=0;
  std::size_t nk=0;
};

// ParsedDocument retains the exact external title and ordered variable tokens
// so test callers can make product-specific assertions after grammar parsing.
struct ParsedDocument {
  std::string title;
  std::vector<std::string> variables;
  std::vector<ParsedZone> zones;
};

// These entry points expose only the validation-side parser.  They do not
// include or call CheckedTextFile, production formatting constants, or writer
// helpers, preserving OUT01's independent-consumer boundary when OUT07 reuses
// the parser for demonstration artifacts.
bool parse_surface_text(const std::string& text,ParsedDocument& document,
                        std::string& error);
bool parse_bundle_text(const std::string& text,ParsedDocument& document,
                       std::string& error);
bool parse_face_text(const std::string& text,ParsedDocument& document,
                     std::string& error);

// Demo-only POINT products have smaller schemas than the canonical 25-field
// 3-D outputs.  The caller supplies an independently declared expected title,
// ordered variables, zone title, and count key ("I" or "N"); the same strict
// numeric-row and exact-EOF checks are then applied.
bool parse_point_text(const std::string& text,const std::string& expected_title,
                      const std::vector<std::string>& expected_variables,
                      const std::string& expected_zone,
                      const std::string& count_key,ParsedDocument& document,
                      std::string& error);

}  // namespace swcme_test::tecplot

#endif
