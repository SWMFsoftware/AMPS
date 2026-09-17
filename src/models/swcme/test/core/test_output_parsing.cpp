#include "test_framework.hpp"
#include "tecplot_parser.hpp"

#include <swcme3d.hpp>

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <locale>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace swcme_test::tecplot {

// OUT01 owns this schema instead of importing a string or table from the
// production writer.  That deliberate duplication is the independence
// boundary: a writer-side omission, reordering, or unit change cannot silently
// update the parser's expectation and make the validation pass.
const std::vector<std::string> kExpectedVariables={
    "X[m]","Y[m]","Z[m]","n[m^-3]","Vx[m/s]","Vy[m/s]","Vz[m/s]",
    "Bx[T]","By[T]","Bz[T]","divVsw[s^-1]","rc[-]","Vsh_n[m/s]",
    "nx[-]","ny[-]","nz[-]","area[m^2]","rc_mean[-]",
    "Vsh_n_mean[m/s]","tnx[-]","tny[-]","tnz[-]","cx[m]","cy[m]",
    "cz[m]"};

// This line-oriented reader is intentionally built from standard C++ text
// primitives, not CheckedTextFile or any writer helper.  It accepts the small
// Tecplot subset SWCME promises, rejects unknown assignments and trailing
// records, and retains the first line-specific diagnostic for test output.
class IndependentTecplotParser {
public:
  explicit IndependentTecplotParser(const std::string& text) {
    std::istringstream input(text);
    input.imbue(std::locale::classic());
    std::string line;
    while (std::getline(input,line)) {
      if (!line.empty() && line.back()=='\r') line.pop_back();
      lines_.push_back(line);
    }
  }

  const std::string& error() const { return error_; }

  bool parse_surface(ParsedDocument& document) {
    if (!parse_header(
            "Shock surface (cell metrics + nodal rc)",kExpectedVariables,
            document))
      return false;
    ParsedZone surface;
    if (!parse_surface_zone(surface)) return false;
    document.zones.push_back(surface);
    return require_end_of_file();
  }

  bool parse_bundle(ParsedDocument& document) {
    if (!parse_header("SW+CME dataset",kExpectedVariables,document))
      return false;

    ParsedZone surface;
    if (!parse_surface_zone(surface)) return false;
    document.zones.push_back(surface);

    ParsedZone nodal;
    if (!parse_zone_assignments(
            {{"T","surface_nodal"},{"N",std::to_string(surface.nodes)},
             {"E",std::to_string(surface.elements)},{"F","FEPOINT"},
             {"ET","TRIANGLE"}}))
      return false;
    nodal.nodes=surface.nodes;
    nodal.elements=surface.elements;
    if (!parse_point_rows(nodal.nodes,kExpectedVariables.size(),nodal.rows) ||
        !parse_connectivity(nodal.elements,nodal.nodes,nodal.triangles))
      return false;
    document.zones.push_back(nodal);

    ParsedZone volume;
    std::map<std::string,std::string> volume_fields;
    if (!take_zone_assignments(volume_fields) ||
        !require_exact_keys(volume_fields,{"T","I","J","K","DATAPACKING"}) ||
        !require_value(volume_fields,"T","volume_box") ||
        !require_value(volume_fields,"DATAPACKING","POINT") ||
        !parse_positive_size(volume_fields["I"],volume.ni,"volume I") ||
        !parse_positive_size(volume_fields["J"],volume.nj,"volume J") ||
        !parse_positive_size(volume_fields["K"],volume.nk,"volume K"))
      return false;
    std::size_t volume_plane=0;
    std::size_t volume_rows=0;
    if (!checked_product(volume.ni,volume.nj,volume_plane,"volume I*J") ||
        !checked_product(volume_plane,volume.nk,volume_rows,"volume I*J*K"))
      return false;
    if (!parse_point_rows(volume_rows,kExpectedVariables.size(),volume.rows))
      return false;
    document.zones.push_back(volume);

    ParsedZone face;
    std::map<std::string,std::string> face_fields;
    std::size_t face_rows=0;
    if (!take_zone_assignments(face_fields) ||
        !parse_face_zone_fields(face_fields,face) ||
        !checked_product(face.ni,face.nj,face_rows,"bundle face I*J") ||
        !parse_point_rows(face_rows,kExpectedVariables.size(),face.rows))
      return false;
    document.zones.push_back(face);
    return require_end_of_file();
  }

  bool parse_face(ParsedDocument& document) {
    if (!parse_header("Box face (minX)",kExpectedVariables,document))
      return false;
    ParsedZone face;
    std::map<std::string,std::string> fields;
    std::size_t face_rows=0;
    if (!take_zone_assignments(fields) || !parse_face_zone_fields(fields,face) ||
        !checked_product(face.ni,face.nj,face_rows,"standalone face I*J") ||
        !parse_point_rows(face_rows,kExpectedVariables.size(),face.rows))
      return false;
    document.zones.push_back(face);
    return require_end_of_file();
  }

  bool parse_point(const std::string& expected_title,
                   const std::vector<std::string>& expected_variables,
                   const std::string& expected_zone,
                   const std::string& count_key,ParsedDocument& document) {
    if (count_key!="I" && count_key!="N")
      return fail("POINT count key must be I or N");
    if (!parse_header(expected_title,expected_variables,document)) return false;
    std::map<std::string,std::string> fields;
    if (!take_zone_assignments(fields) ||
        !require_exact_keys(fields,{"T",count_key,"F"}) ||
        !require_value(fields,"T",expected_zone) ||
        !require_value(fields,"F","POINT"))
      return false;
    ParsedZone zone;
    std::size_t row_count=0;
    if (!parse_positive_size(fields[count_key],row_count,"POINT row count") ||
        !parse_point_rows(row_count,expected_variables.size(),zone.rows))
      return false;
    if (count_key=="I") zone.ni=row_count;
    else zone.nodes=row_count;
    document.zones.push_back(zone);
    return require_end_of_file();
  }

private:
  static std::string trim(const std::string& value) {
    const std::size_t first=value.find_first_not_of(" \t");
    if (first==std::string::npos) return "";
    const std::size_t last=value.find_last_not_of(" \t");
    return value.substr(first,last-first+1);
  }

  bool fail(const std::string& message) {
    if (error_.empty())
      error_="line "+std::to_string(cursor_+1)+": "+message;
    return false;
  }

  bool take_line(std::string& line) {
    if (cursor_>=lines_.size()) return fail("unexpected end of file");
    line=lines_[cursor_++];
    if (line.empty()) return fail("unexpected empty record");
    return true;
  }

  // TITLE is parsed as grammar rather than compared as raw bytes so harmless
  // whitespace around '=' is accepted while missing quotes, trailing content,
  // or a changed semantic title remains a failure.
  bool parse_quoted_assignment(const std::string& line,const std::string& key,
                               std::string& value) {
    std::size_t position=0;
    while (position<line.size() && (line[position]==' ' || line[position]=='\t'))
      ++position;
    if (line.compare(position,key.size(),key)!=0) return fail("expected "+key);
    position+=key.size();
    while (position<line.size() && (line[position]==' ' || line[position]=='\t'))
      ++position;
    if (position>=line.size() || line[position++]!='=')
      return fail("missing '=' after "+key);
    while (position<line.size() && (line[position]==' ' || line[position]=='\t'))
      ++position;
    if (position>=line.size() || line[position++]!='\"')
      return fail("missing opening quote after "+key);
    const std::size_t closing=line.find('\"',position);
    if (closing==std::string::npos) return fail("unterminated quoted "+key);
    value=line.substr(position,closing-position);
    if (!trim(line.substr(closing+1)).empty())
      return fail("trailing content after "+key);
    return true;
  }

  bool parse_header(const std::string& expected_title,
                    const std::vector<std::string>& expected_variables,
                    ParsedDocument& document) {
    std::string line;
    if (!take_line(line) ||
        !parse_quoted_assignment(line,"TITLE",document.title))
      return false;
    if (document.title!=expected_title) return fail("unexpected TITLE value");
    if (!take_line(line) || !parse_variable_list(line,document.variables))
      return false;
    if (document.variables!=expected_variables)
      return fail("VARIABLES names, units, or order differ from validation schema");
    return true;
  }

  // VARIABLES uses a strict quoted comma-list.  No splitting helper from the
  // emitter is shared, and the expected unit-qualified tokens above are owned
  // by the test, making format drift observable.
  bool parse_variable_list(const std::string& line,
                           std::vector<std::string>& variables) {
    std::size_t position=0;
    while (position<line.size() && (line[position]==' ' || line[position]=='\t'))
      ++position;
    const std::string keyword="VARIABLES";
    if (line.compare(position,keyword.size(),keyword)!=0)
      return fail("expected VARIABLES record");
    position+=keyword.size();
    while (position<line.size() && (line[position]==' ' || line[position]=='\t'))
      ++position;
    if (position>=line.size() || line[position++]!='=')
      return fail("missing '=' after VARIABLES");
    for (;;) {
      while (position<line.size() && (line[position]==' ' || line[position]=='\t'))
        ++position;
      if (position>=line.size() || line[position++]!='\"')
        return fail("VARIABLES item is not quoted");
      const std::size_t closing=line.find('\"',position);
      if (closing==std::string::npos)
        return fail("unterminated VARIABLES item");
      variables.push_back(line.substr(position,closing-position));
      position=closing+1;
      while (position<line.size() && (line[position]==' ' || line[position]=='\t'))
        ++position;
      if (position==line.size()) break;
      if (line[position++]!=',') return fail("expected VARIABLES comma");
    }
    return !variables.empty() || fail("empty VARIABLES list");
  }

  // Zone declarations are tokenized without depending on assignment order.
  // Quoted zone titles are unwrapped, duplicate/unknown keys are rejected, and
  // a trailing comma is accepted because Tecplot permits the following
  // VARLOCATION record to continue the declaration.
  bool take_zone_assignments(std::map<std::string,std::string>& fields) {
    std::string line;
    if (!take_line(line)) return false;
    const std::string prefix="ZONE";
    std::size_t position=line.find_first_not_of(" \t");
    if (position==std::string::npos ||
        line.compare(position,prefix.size(),prefix)!=0 ||
        (position+prefix.size()<line.size() &&
         line[position+prefix.size()]!=' ' &&
         line[position+prefix.size()]!='\t'))
      return fail("expected ZONE declaration");
    position+=prefix.size();
    std::string current;
    bool quoted=false;
    for (; position<=line.size(); ++position) {
      const char character=position<line.size() ? line[position] : ',';
      if (character=='\"') quoted=!quoted;
      if (character==',' && !quoted) {
        const std::string assignment=trim(current);
        current.clear();
        if (assignment.empty()) continue;
        const std::size_t equals=assignment.find('=');
        if (equals==std::string::npos)
          return fail("ZONE item lacks '='");
        const std::string key=trim(assignment.substr(0,equals));
        std::string value=trim(assignment.substr(equals+1));
        if (value.size()>=2 && value.front()=='\"' && value.back()=='\"')
          value=value.substr(1,value.size()-2);
        if (key.empty() || value.empty() || fields.count(key)!=0)
          return fail("invalid or duplicate ZONE assignment");
        fields[key]=value;
      } else {
        current.push_back(character);
      }
    }
    if (quoted) return fail("unterminated quote in ZONE declaration");
    return true;
  }

  bool require_exact_keys(const std::map<std::string,std::string>& fields,
                          const std::vector<std::string>& keys) {
    if (fields.size()!=keys.size()) return fail("wrong ZONE assignment count");
    for (const std::string& key : keys)
      if (fields.count(key)==0) return fail("missing ZONE field "+key);
    return true;
  }

  bool require_value(const std::map<std::string,std::string>& fields,
                     const std::string& key,const std::string& expected) {
    const auto item=fields.find(key);
    return (item!=fields.end() && item->second==expected) ||
        fail("unexpected ZONE value for "+key);
  }

  bool parse_positive_size(const std::string& text,std::size_t& value,
                           const std::string& label) {
    const char* begin=text.data();
    const char* end=begin+text.size();
    const std::from_chars_result result=std::from_chars(begin,end,value);
    if (result.ec!=std::errc() || result.ptr!=end || value==0)
      return fail("invalid positive count for "+label);
    return true;
  }

  // Malformed validation input must not wrap its declared cardinality and make
  // a truncated file appear complete.  This independent checked multiply is
  // intentionally local to the parser rather than reusing OUT04 production
  // validation, which operates before serialization rather than after it.
  bool checked_product(std::size_t left,std::size_t right,std::size_t& product,
                       const std::string& label) {
    if (right!=0 && left>std::numeric_limits<std::size_t>::max()/right)
      return fail("overflowed declared count for "+label);
    product=left*right;
    return true;
  }

  bool parse_zone_assignments(
      const std::map<std::string,std::string>& expected) {
    std::map<std::string,std::string> actual;
    if (!take_zone_assignments(actual)) return false;
    return actual==expected || fail("ZONE declaration differs from schema");
  }

  bool parse_surface_zone(ParsedZone& zone) {
    std::map<std::string,std::string> fields;
    if (!take_zone_assignments(fields) ||
        !require_exact_keys(fields,{"T","N","E","ZONETYPE","DATAPACKING"}) ||
        !require_value(fields,"T","surface_cells") ||
        !require_value(fields,"ZONETYPE","FETRIANGLE") ||
        !require_value(fields,"DATAPACKING","BLOCK") ||
        !parse_positive_size(fields["N"],zone.nodes,"surface N") ||
        !parse_positive_size(fields["E"],zone.elements,"surface E"))
      return false;

    std::string locations;
    if (!take_line(locations)) return false;
    std::string compact;
    for (char character : locations)
      if (character!=' ' && character!='\t') compact.push_back(character);
    if (compact!="VARLOCATION=([1-3,12-13]=NODAL,[4-11,14-25]=CELLCENTERED)")
      return fail("unexpected VARLOCATION contract");

    zone.blocks.resize(kExpectedVariables.size());
    for (std::size_t variable=0; variable<zone.blocks.size(); ++variable) {
      const bool nodal=variable<3 || variable==11 || variable==12;
      if (!parse_block(nodal ? zone.nodes : zone.elements,
                       zone.blocks[variable]))
        return false;
    }
    return parse_connectivity(zone.elements,zone.nodes,zone.triangles);
  }

  bool parse_face_zone_fields(const std::map<std::string,std::string>& fields,
                              ParsedZone& zone) {
    return require_exact_keys(fields,{"T","I","J","DATAPACKING"}) &&
        require_value(fields,"T","box_face_minX") &&
        require_value(fields,"DATAPACKING","POINT") &&
        parse_positive_size(fields.at("I"),zone.ni,"face I") &&
        parse_positive_size(fields.at("J"),zone.nj,"face J");
  }

  // SWCME's BLOCK writer intentionally wraps at eight values.  Checking each
  // physical line as well as the aggregate count catches missing separators,
  // merged variables, and an early next-zone header instead of merely proving
  // that some total number of numeric tokens exists.
  bool parse_block(std::size_t expected_values,std::vector<double>& values) {
    std::size_t remaining=expected_values;
    while (remaining>0) {
      std::vector<double> line_values;
      if (!parse_numeric_line(line_values)) return false;
      const std::size_t expected_on_line=remaining<8 ? remaining : 8;
      if (line_values.size()!=expected_on_line)
        return fail("wrong BLOCK line width");
      values.insert(values.end(),line_values.begin(),line_values.end());
      remaining-=expected_on_line;
    }
    return true;
  }

  bool parse_point_rows(std::size_t expected_rows,std::size_t expected_width,
                        std::vector<std::vector<double>>& rows) {
    rows.reserve(expected_rows);
    for (std::size_t row=0; row<expected_rows; ++row) {
      std::vector<double> values;
      if (!parse_numeric_line(values)) return false;
      if (values.size()!=expected_width) return fail("wrong POINT row width");
      rows.push_back(values);
    }
    return true;
  }

  bool parse_numeric_line(std::vector<double>& values) {
    std::string line;
    if (!take_line(line)) return false;
    std::istringstream input(line);
    input.imbue(std::locale::classic());
    double value=0.0;
    while (input>>value) {
      if (!std::isfinite(value)) return fail("nonfinite numeric value");
      values.push_back(value);
    }
    input.clear();
    std::string residue;
    input>>residue;
    if (!residue.empty()) return fail("nonnumeric data token");
    return !values.empty() || fail("empty numeric record");
  }

  bool parse_connectivity(
      std::size_t expected_rows,std::size_t node_count,
      std::vector<std::array<std::size_t,3>>& triangles) {
    triangles.reserve(expected_rows);
    for (std::size_t row=0; row<expected_rows; ++row) {
      std::string line;
      if (!take_line(line)) return false;
      std::istringstream input(line);
      input.imbue(std::locale::classic());
      long long a=0,b=0,c=0;
      std::string extra;
      if (!(input>>a>>b>>c) || (input>>extra))
        return fail("connectivity row must contain exactly three integers");
      if (a<1 || b<1 || c<1 ||
          static_cast<unsigned long long>(a)>node_count ||
          static_cast<unsigned long long>(b)>node_count ||
          static_cast<unsigned long long>(c)>node_count ||
          a==b || a==c || b==c)
        return fail("connectivity index is invalid");
      triangles.push_back({static_cast<std::size_t>(a),
                           static_cast<std::size_t>(b),
                           static_cast<std::size_t>(c)});
    }
    return true;
  }

  bool require_end_of_file() {
    return cursor_==lines_.size() || fail("extra record after final zone");
  }

  std::vector<std::string> lines_;
  std::size_t cursor_=0;
  std::string error_;
};

std::string read_text(const std::filesystem::path& path) {
  std::ifstream input(path,std::ios::binary);
  return std::string(std::istreambuf_iterator<char>(input),
                     std::istreambuf_iterator<char>());
}

bool parse_surface_text(const std::string& text,ParsedDocument& document,
                        std::string& error) {
  IndependentTecplotParser parser(text);
  const bool valid=parser.parse_surface(document);
  error=parser.error();
  return valid;
}

bool parse_bundle_text(const std::string& text,ParsedDocument& document,
                       std::string& error) {
  IndependentTecplotParser parser(text);
  const bool valid=parser.parse_bundle(document);
  error=parser.error();
  return valid;
}

bool parse_face_text(const std::string& text,ParsedDocument& document,
                     std::string& error) {
  IndependentTecplotParser parser(text);
  const bool valid=parser.parse_face(document);
  error=parser.error();
  return valid;
}

bool parse_point_text(const std::string& text,const std::string& expected_title,
                      const std::vector<std::string>& expected_variables,
                      const std::string& expected_zone,
                      const std::string& count_key,ParsedDocument& document,
                      std::string& error) {
  IndependentTecplotParser parser(text);
  const bool valid=parser.parse_point(
      expected_title,expected_variables,expected_zone,count_key,document);
  error=parser.error();
  return valid;
}

bool serialized_near(double actual,double expected) {
  // Production uses %.9e, i.e. ten significant decimal digits.  One part in
  // 1e9 safely covers final-digit rounding while remaining far tighter than a
  // physical-model tolerance that could conceal a wrong field or unit.
  const double scale=std::max(std::abs(expected),1.0e-300);
  return std::isfinite(actual) && std::abs(actual-expected)<=1.0e-9*scale;
}

std::vector<std::array<std::size_t,3>> expected_triangles(
    const swcme3d::ShockMesh& mesh) {
  std::vector<std::array<std::size_t,3>> triangles;
  triangles.reserve(mesh.tri_i.size());
  for (std::size_t element=0; element<mesh.tri_i.size(); ++element)
    triangles.push_back({static_cast<std::size_t>(mesh.tri_i[element]),
                         static_cast<std::size_t>(mesh.tri_j[element]),
                         static_cast<std::size_t>(mesh.tri_k[element])});
  return triangles;
}

swcme3d::Params output_params() {
  // The deterministic sphere is comfortably outside the solar-wind floor, has
  // no finite-cap edge cases, and keeps OUT01 focused on serialization rather
  // than nonlinear shock or geometry failure behavior.
  swcme3d::Params params;
  params.shape=swcme3d::ShockShape::Sphere;
  params.kinematics_mode=swcme::kinematics::Mode::Ballistic;
  params.r0_Rs=40.0;
  params.V0_sh_kms=1200.0;
  params.V_sw_kms=400.0;
  params.cme_dir[0]=0.73;
  params.cme_dir[1]=-0.41;
  params.cme_dir[2]=0.547;
  params.region_mode=swcme::regions::Mode::ShockOnly;
  params.shock_acceleration_mode=swcme::acceleration::Mode::Source;
  return params;
}

void remove_if_present(const std::filesystem::path& path) {
  std::error_code error;
  (void)std::filesystem::remove(path,error);
}

}  // namespace swcme_test::tecplot

// OUT01: validate the bytes consumed by external Tecplot readers rather than
// trusting successful writer statuses.  Three public products are written
// through the production filesystem backend and then parsed by the independent
// grammar above: surface-only, the four-zone surface/volume/face bundle, and
// the standalone min-X face.
void test_out01(swcme_test::Context& context) {
  using namespace swcme_test::tecplot;
  std::cout << "OUT01 independent output parsing\n";

  const std::filesystem::path surface_path="output/OUT01_surface.dat";
  const std::filesystem::path bundle_path="output/OUT01_bundle.dat";
  const std::filesystem::path face_path="output/OUT01_face.dat";
  remove_if_present(surface_path);
  remove_if_present(bundle_path);
  remove_if_present(face_path);

  swcme3d::Model model(output_params());
  const swcme3d::StepState step=model.prepare_step(0.0);
  const swcme3d::ShockMesh mesh=model.build_shock_mesh(step,3,6);
  swcme3d::TriMetrics metrics;
  model.compute_triangle_metrics(mesh,metrics);
  const swcme3d::BoxSpec box=model.default_apex_box(step,0.02,2);

  const swcme::ModelStatus surface_status=
      model.write_shock_surface_center_metrics_tecplot_checked(
          mesh,metrics,surface_path.c_str());
  const swcme::ModelStatus bundle_status=
      model.write_tecplot_dataset_bundle_checked(
          mesh,metrics,step,box,bundle_path.c_str());
  const swcme::ModelStatus face_status=
      model.write_box_face_minX_tecplot_structured_checked(
          step,box,face_path.c_str());
  context.expect_true(surface_status.ok() && bundle_status.ok() &&
                          face_status.ok(),
                      "all deterministic OUT01 products commit successfully");

  const std::string surface_text=read_text(surface_path);
  const std::string bundle_text=read_text(bundle_path);
  const std::string face_text=read_text(face_path);

  ParsedDocument surface_document;
  std::string surface_error;
  const bool surface_valid=parse_surface_text(
      surface_text,surface_document,surface_error);
  context.expect_true(surface_valid,
                      "surface output satisfies independent schema: "+
                          surface_error);

  ParsedDocument bundle_document;
  std::string bundle_error;
  const bool bundle_valid=parse_bundle_text(
      bundle_text,bundle_document,bundle_error);
  context.expect_true(bundle_valid,
                      "bundle output satisfies independent schema: "+
                          bundle_error);

  ParsedDocument face_document;
  std::string face_error;
  const bool face_valid=parse_face_text(face_text,face_document,face_error);
  context.expect_true(face_valid,
                      "standalone face satisfies independent schema: "+
                          face_error);

  const std::vector<std::array<std::size_t,3>> topology=
      expected_triangles(mesh);
  if (surface_valid) {
    const ParsedZone& zone=surface_document.zones[0];
    context.expect_true(zone.nodes==mesh.x.size() &&
                            zone.elements==mesh.tri_i.size() &&
                            zone.triangles==topology,
                        "surface declarations and connectivity match the mesh");
    context.expect_true(serialized_near(zone.blocks[0][0],mesh.x[0]) &&
                            serialized_near(zone.blocks[11][0],mesh.rc[0]) &&
                            serialized_near(zone.blocks[16][0],metrics.area[0]),
                        "surface nodal and cell values retain formatting precision");
  }

  if (bundle_valid) {
    const ParsedZone& cell_zone=bundle_document.zones[0];
    const ParsedZone& nodal_zone=bundle_document.zones[1];
    const ParsedZone& volume_zone=bundle_document.zones[2];
    const ParsedZone& bundled_face=bundle_document.zones[3];
    context.expect_true(bundle_document.zones.size()==4 &&
                            cell_zone.triangles==topology &&
                            nodal_zone.triangles==topology,
                        "bundle contains exactly four zones with preserved topology");
    context.expect_true(volume_zone.ni==static_cast<std::size_t>(box.Ni) &&
                            volume_zone.nj==static_cast<std::size_t>(box.Nj) &&
                            volume_zone.nk==static_cast<std::size_t>(box.Nk) &&
                            volume_zone.rows.size()==8 &&
                            bundled_face.rows.size()==4,
                        "volume and bundled-face dimensions match the 2x2x2 request");
    context.expect_true(serialized_near(nodal_zone.rows[0][0],mesh.x[0]) &&
                            serialized_near(nodal_zone.rows[0][11],mesh.rc[0]) &&
                            serialized_near(nodal_zone.rows[0][13],mesh.n_hat_x[0]),
                        "surface POINT values match direct mesh fields");

    const double x=box.cx-box.hx;
    const double y=box.cy-box.hy;
    const double z=box.cz-box.hz;
    double n=0.0,Vx=0.0,Vy=0.0,Vz=0.0,Bx=0.0,By=0.0,Bz=0.0,div=0.0;
    model.evaluate_cartesian_with_B(
        step,&x,&y,&z,&n,&Vx,&Vy,&Vz,&Bx,&By,&Bz,1);
    model.compute_divV_radial(step,&x,&y,&z,&div,1,1e-3);
    const std::array<double,11> expected={x,y,z,n,Vx,Vy,Vz,Bx,By,Bz,div};
    bool volume_values_match=true;
    bool face_values_match=true;
    for (std::size_t field=0; field<expected.size(); ++field) {
      volume_values_match=volume_values_match &&
          serialized_near(volume_zone.rows[0][field],expected[field]);
      face_values_match=face_values_match &&
          serialized_near(bundled_face.rows[0][field],expected[field]);
    }
    context.expect_true(volume_values_match && face_values_match,
                        "volume and bundled face agree with direct API evaluation");
  }

  if (bundle_valid && face_valid) {
    context.expect_true(face_document.zones[0].ni==
                            static_cast<std::size_t>(box.Nj) &&
                            face_document.zones[0].nj==
                            static_cast<std::size_t>(box.Nk) &&
                            face_document.zones[0].rows==
                                bundle_document.zones[3].rows,
                        "standalone and bundled min-X faces are byte-value equivalent");
  }

  // Parser sensitivity probes are deliberately derived after successful
  // production parsing.  Each changes one independent contract dimension and
  // proves the validator would not bless extra records, unit drift, malformed
  // connectivity, or a widened POINT row.
  ParsedDocument corrupt_document;
  std::string corrupt_error;
  context.expect_true(!parse_surface_text(
                          surface_text+"0\n",corrupt_document,corrupt_error),
                      "independent parser rejects trailing records");

  std::string wrong_unit=surface_text;
  const std::size_t unit_position=wrong_unit.find("\"X[m]\"");
  if (unit_position!=std::string::npos)
    wrong_unit.replace(unit_position,6,"\"X[km]\"");
  corrupt_document=ParsedDocument{};
  corrupt_error.clear();
  context.expect_true(unit_position!=std::string::npos &&
                          !parse_surface_text(
                              wrong_unit,corrupt_document,corrupt_error),
                      "independent parser rejects variable-unit drift");

  std::string bad_connectivity=surface_text;
  if (!bad_connectivity.empty() && bad_connectivity.back()=='\n')
    bad_connectivity.pop_back();
  const std::size_t last_line=bad_connectivity.rfind('\n');
  if (last_line!=std::string::npos) {
    const std::size_t first_space=bad_connectivity.find(' ',last_line+1);
    if (first_space!=std::string::npos)
      bad_connectivity.replace(last_line+1,first_space-last_line-1,"0");
  }
  corrupt_document=ParsedDocument{};
  corrupt_error.clear();
  context.expect_true(last_line!=std::string::npos &&
                          !parse_surface_text(
                              bad_connectivity,corrupt_document,corrupt_error),
                      "independent parser rejects out-of-range connectivity");

  std::string wide_face=face_text;
  const std::string face_zone="DATAPACKING=POINT\n";
  const std::size_t face_zone_end=wide_face.find(face_zone);
  const std::size_t first_row_end=face_zone_end==std::string::npos
      ? std::string::npos : wide_face.find('\n',face_zone_end+face_zone.size());
  if (first_row_end!=std::string::npos) wide_face.insert(first_row_end," 0");
  corrupt_document=ParsedDocument{};
  corrupt_error.clear();
  context.expect_true(first_row_end!=std::string::npos &&
                          !parse_face_text(
                              wide_face,corrupt_document,corrupt_error),
                      "independent parser rejects incorrect POINT row width");

  remove_if_present(surface_path);
  remove_if_present(bundle_path);
  remove_if_present(face_path);
}
