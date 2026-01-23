#include "cli.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <algorithm>

// ----------- small utilities -----------
static inline std::string ltrim(const std::string& s) {
  size_t i=0;
  while (i<s.size() && std::isspace((unsigned char)s[i])) ++i;
  return s.substr(i);
}
static inline std::string rtrim(const std::string& s) {
  if (s.empty()) return s;
  size_t i=s.size();
  while (i>0 && std::isspace((unsigned char)s[i-1])) --i;
  return s.substr(0,i);
}
static inline std::string trim(const std::string& s) { return rtrim(ltrim(s)); }

static inline std::string tolower_str(std::string s) {
  for (auto& c: s) c = (char)std::tolower((unsigned char)c);
  return s;
}

// Strip comments: supports '#', ';', '//' (first occurrence wins).
static std::string StripComments(const std::string& line) {
  size_t cut = std::string::npos;

  auto upd = [&](size_t pos) {
    if (pos!=std::string::npos) cut = (cut==std::string::npos) ? pos : std::min(cut,pos);
  };

  upd(line.find('#'));
  upd(line.find(';'));
  upd(line.find("//"));

  if (cut==std::string::npos) return line;
  return line.substr(0,cut);
}

static bool ParseDoublesCSV(const std::string& rhs, std::vector<double>& out, std::vector<std::string>& outStr) {
  out.clear(); outStr.clear();
  std::string s = rhs;
  // Allow whitespace around commas
  std::stringstream ss(s);
  std::string tok;
  while (std::getline(ss, tok, ',')) {
    tok = trim(tok);
    if (tok.empty()) continue;

    char* end=nullptr;
    double v = std::strtod(tok.c_str(), &end);
    if (end && end!=tok.c_str()) {
      // accept if remainder is whitespace
      while (*end && std::isspace((unsigned char)*end)) ++end;
      if (*end=='\0') { out.push_back(v); continue; }
    }
    // non-numeric token
    outStr.push_back(tok);
  }
  return (!out.empty() || !outStr.empty());
}

static bool ToBool(const std::vector<double>& nums, const std::vector<std::string>& strs, bool defaultVal=false) {
  if (!strs.empty()) {
    std::string t = tolower_str(strs[0]);
    if (t=="1" || t=="true" || t=="yes" || t=="on") return true;
    if (t=="0" || t=="false"|| t=="no"  || t=="off") return false;
  }
  if (!nums.empty()) return (nums[0]!=0.0);
  return defaultVal;
}

static void ApplyKeyValue(TestConfig& cfg, const std::string& keyRaw,
                          const std::vector<double>& nums, const std::vector<std::string>& strs) {
  // Normalize key: lowercase; remove leading '-' ; '_' -> '-'
  std::string key = tolower_str(trim(keyRaw));
  while (!key.empty() && key[0]=='-') key.erase(key.begin());
  std::replace(key.begin(), key.end(), '_', '-');

  auto need3 = [&](double v[3]) {
    if (nums.size()>=3) { v[0]=nums[0]; v[1]=nums[1]; v[2]=nums[2]; return true; }
    return false;
  };

  // Mode switches
  if (key=="mode") {
    if (!strs.empty()) {
      std::string m = tolower_str(strs[0]);
      if (m=="particles" || m=="withparticles" || m=="with-particles") cfg.mode = TestConfig::Mode::WithParticles;
      else if (m=="no-particles" || m=="noparticles" || m=="field-only") cfg.mode = TestConfig::Mode::NoParticles;
      else if (m=="fieldonlyb" || m=="b") cfg.mode = TestConfig::Mode::FieldOnlyB;
      else if (m=="fieldonlye" || m=="e") cfg.mode = TestConfig::Mode::FieldOnlyE;
      return;
    }
    // numeric fallback: 1 -> particles, 0 -> no-particles
    cfg.mode = (ToBool(nums,strs,true) ? TestConfig::Mode::WithParticles : TestConfig::Mode::NoParticles);
    return;
  }
// Domain boundary condition type: bc=dirichlet|neumann  (applied on all 6 faces)
if (key=="bc" || key=="bc-type" || key=="domain-bc") {
  std::string v;
  if (!strs.empty()) v = tolower_str(strs[0]);
  else if (!nums.empty()) v = (nums[0]!=0.0 ? "dirichlet" : "neumann"); // allow bc=1/0 as shorthand
  v = tolower_str(trim(v));
  if (!v.empty()) {
    cfg.user_domain_bc = true;
    if (v=="d" || v=="dirichlet") cfg.domain_bc = TestConfig::DomainBCType::Dirichlet;
    else if (v=="n" || v=="neumann") cfg.domain_bc = TestConfig::DomainBCType::Neumann;
    else {
      std::fprintf(stderr,"[ConstFieldBC] WARNING: unknown bc='%s' in input; keeping default\n", v.c_str());
    }
  }
  return;
}

// Per-face domain boundary condition overrides:
//   bc-xmin=dirichlet|neumann, ..., bc-zmax=...
// Accepted synonyms (after normalization): bc-xmin, domain-bc-xmin, bc-type-xmin.
{
  auto parse_bc = [&](const std::vector<double>& numsL, const std::vector<std::string>& strsL, bool& ok)->TestConfig::DomainBCType {
    std::string v;
    if (!strsL.empty()) v = tolower_str(strsL[0]);
    else if (!numsL.empty()) v = (numsL[0]!=0.0 ? "dirichlet" : "neumann"); // allow 1/0 shorthand
    v = tolower_str(trim(v));
    ok = true;
    if (v=="d" || v=="dirichlet") return TestConfig::DomainBCType::Dirichlet;
    if (v=="n" || v=="neumann")   return TestConfig::DomainBCType::Neumann;
    ok = false;
    return cfg.domain_bc; // fallback; caller decides whether to apply
  };

  auto try_face = [&](const std::string& faceKey, int faceIdx)->bool {
    if (key==faceKey || key==("domain-bc-"+faceKey.substr(3)) || key==("bc-type-"+faceKey.substr(3))) {
      bool ok=false;
      TestConfig::DomainBCType t = parse_bc(nums, strs, ok);
      if (ok) {
        cfg.user_domain_bc_face[faceIdx] = true;
        cfg.domain_bc_face[faceIdx] = t;
      }
      else {
        std::fprintf(stderr,"[ConstFieldBC] WARNING: unknown %s in input; ignoring\n", key.c_str());
      }
      return true;
    }
    return false;
  };

  if (try_face("bc-xmin",0)) return;
  if (try_face("bc-xmax",1)) return;
  if (try_face("bc-ymin",2)) return;
  if (try_face("bc-ymax",3)) return;
  if (try_face("bc-zmin",4)) return;
  if (try_face("bc-zmax",5)) return;
}

  if (key=="particles" || key=="withparticles") {
    cfg.mode = (ToBool(nums,strs,true) ? TestConfig::Mode::WithParticles : cfg.mode);
    return;
  }
  if (key=="no-particles" || key=="noparticles") {
    if (ToBool(nums,strs,true)) cfg.mode = TestConfig::Mode::NoParticles;
    return;
  }

  // Domain
  if (key=="l" || key=="domain-l" || key=="domain" || key=="box") {
    if (nums.size()==1) { cfg.use_domain_L=true; cfg.domain_L[0]=cfg.domain_L[1]=cfg.domain_L[2]=nums[0]; return; }
    if (nums.size()>=3) { cfg.use_domain_L=true; cfg.domain_L[0]=nums[0]; cfg.domain_L[1]=nums[1]; cfg.domain_L[2]=nums[2]; return; }
    return;
  }

  // Background fields
  if (key=="b" || key=="b0") {
    double tmp[3];
    if (need3(tmp)) {
      cfg.userB=true;
      cfg.B0[0]=tmp[0]; cfg.B0[1]=tmp[1]; cfg.B0[2]=tmp[2];
      // If not in particle mode, interpret as field-only B init.
      if (cfg.mode != TestConfig::Mode::WithParticles) cfg.mode = TestConfig::Mode::FieldOnlyB;
    }
    return;
  }
  if (key=="e" || key=="e0") {
    double tmp[3];
    if (need3(tmp)) {
      cfg.userE=true;
      cfg.E0[0]=tmp[0]; cfg.E0[1]=tmp[1]; cfg.E0[2]=tmp[2];
      if (cfg.mode != TestConfig::Mode::WithParticles) cfg.mode = TestConfig::Mode::FieldOnlyE;
    }
    return;
  }

  // Background fields in *physical* units (SI). Values are converted to solver units in
  // FinalizeConfigUnits(). These keys are available both on the CLI (as -BnT/-EmVm/-EVm)
  // and in the input file (as BnT=..., EmVm=..., EVm=...).
  if (key=="bnt" || key=="b-nt") {
    double tmp[3];
    if (need3(tmp)) {
      cfg.userB_SI = true;
      cfg.B0_SI_T[0] = tmp[0]*1.0e-9;
      cfg.B0_SI_T[1] = tmp[1]*1.0e-9;
      cfg.B0_SI_T[2] = tmp[2]*1.0e-9;
      if (cfg.mode != TestConfig::Mode::WithParticles) cfg.mode = TestConfig::Mode::FieldOnlyB;
    }
    return;
  }
  if (key=="evm" || key=="e-vm") {
    double tmp[3];
    if (need3(tmp)) {
      cfg.userE_SI = true;
      cfg.E0_SI_Vm[0] = tmp[0];
      cfg.E0_SI_Vm[1] = tmp[1];
      cfg.E0_SI_Vm[2] = tmp[2];
      if (cfg.mode != TestConfig::Mode::WithParticles) cfg.mode = TestConfig::Mode::FieldOnlyE;
    }
    return;
  }
  if (key=="emvm" || key=="e-mvm") {
    double tmp[3];
    if (need3(tmp)) {
      cfg.userE_SI = true;
      cfg.E0_SI_Vm[0] = tmp[0]*1.0e-3;
      cfg.E0_SI_Vm[1] = tmp[1]*1.0e-3;
      cfg.E0_SI_Vm[2] = tmp[2]*1.0e-3;
      if (cfg.mode != TestConfig::Mode::WithParticles) cfg.mode = TestConfig::Mode::FieldOnlyE;
    }
    return;
  }

  // Stencil
  if (key=="stencil-order" || key=="stencilorder") {
    if (!nums.empty()) cfg.stencilOrder = (int)nums[0];
    return;
  }

  // Solar wind IC in test units
  if (key=="sw" || key=="solar-wind" || key=="solarwind") {
    cfg.mode = TestConfig::Mode::WithParticles;
    return;
  }
  if (key=="sw-rho" || key=="sw-rho0") {
    if (!nums.empty()) { cfg.mode = TestConfig::Mode::WithParticles; cfg.sw_rho0 = nums[0]; }
    return;
  }
  if (key=="sw-p" || key=="sw-p0") {
    if (!nums.empty()) { cfg.mode = TestConfig::Mode::WithParticles; cfg.sw_p0 = nums[0]; }
    return;
  }
  if (key=="sw-u" || key=="sw-u0") {
    if (nums.size()>=3) { cfg.mode = TestConfig::Mode::WithParticles; cfg.sw_u0[0]=nums[0]; cfg.sw_u0[1]=nums[1]; cfg.sw_u0[2]=nums[2]; }
    return;
  }
  if (key=="sw-no-round" || key=="sw-noround") {
    // sw-no-round=1 disables rounding
    if (ToBool(nums,strs,true)) cfg.sw_use_rounding = false;
    return;
  }

  // Units normalization scales used when converting SI inputs to normalized/code units.
  // These are ONLY used when the ECSIM build expects normalized units, but can be set in any build.
  if (key=="units-lm" || key=="units-l" || key=="units-length-m") {
    if (!nums.empty() && nums[0]>0.0) { cfg.units_user_set=true; cfg.units_lSI_m=nums[0]; }
    return;
  }
  if (key=="units-umps" || key=="units-u" || key=="units-speed-mps") {
    if (!nums.empty() && nums[0]>0.0) { cfg.units_user_set=true; cfg.units_uSI_mps=nums[0]; }
    return;
  }
  if (key=="units-mkg" || key=="units-m" || key=="units-mass-kg") {
    if (!nums.empty() && nums[0]>0.0) { cfg.units_user_set=true; cfg.units_mSI_kg=nums[0]; }
    return;
  }

  // Physical solar-wind inputs (SI). These are converted to solver units in FinalizeConfigUnits().
  if (key=="sw-ncm3" || key=="sw-n-cm3") {
    if (!nums.empty()) { cfg.mode = TestConfig::Mode::WithParticles; cfg.sw_has_ncm3=true; cfg.sw_n_cm3=nums[0]; }
    return;
  }
  if (key=="sw-tk" || key=="sw-t") {
    if (!nums.empty()) { cfg.mode = TestConfig::Mode::WithParticles; cfg.sw_has_TK=true; cfg.sw_TK=nums[0]; }
    return;
  }
  if (key=="sw-bnt") {
    if (nums.size()>=3) {
      cfg.mode = TestConfig::Mode::WithParticles;
      cfg.sw_has_BnT=true;
      cfg.sw_BnT[0]=nums[0]; cfg.sw_BnT[1]=nums[1]; cfg.sw_BnT[2]=nums[2];
      // Also treat as background-field specification in SI.
      cfg.userB_SI=true;
      cfg.B0_SI_T[0]=nums[0]*1.0e-9; cfg.B0_SI_T[1]=nums[1]*1.0e-9; cfg.B0_SI_T[2]=nums[2]*1.0e-9;
    }
    return;
  }
  if (key=="sw-u-kms" || key=="sw-ukms") {
    if (nums.size()>=3) {
      cfg.mode = TestConfig::Mode::WithParticles;
      cfg.sw_has_u_kms=true;
      cfg.sw_u_kms[0]=nums[0]; cfg.sw_u_kms[1]=nums[1]; cfg.sw_u_kms[2]=nums[2];
    }
    return;
  }
  if (key=="sw-u-mps" || key=="sw-umps") {
    if (nums.size()>=3) {
      cfg.mode = TestConfig::Mode::WithParticles;
      cfg.sw_has_u_mps=true;
      cfg.sw_u_mps[0]=nums[0]; cfg.sw_u_mps[1]=nums[1]; cfg.sw_u_mps[2]=nums[2];
    }
    return;
  }
  if (key=="sw-evxb") {
    cfg.mode = TestConfig::Mode::WithParticles;
    cfg.sw_evxb = ToBool(nums,strs,true) ? 1 : 0;
    return;
  }
  if (key=="sw-emvm") {
    if (nums.size()>=3) {
      cfg.mode = TestConfig::Mode::WithParticles;
      cfg.sw_has_EmVm=true;
      cfg.sw_E_mVm[0]=nums[0]; cfg.sw_E_mVm[1]=nums[1]; cfg.sw_E_mVm[2]=nums[2];
      // If an explicit E is provided, it should take precedence over E=-u×B.
      cfg.sw_evxb = 0;
    }
    return;
  }
  if (key=="sw-evm") {
    if (nums.size()>=3) {
      cfg.mode = TestConfig::Mode::WithParticles;
      cfg.sw_has_EVm=true;
      cfg.sw_E_Vm[0]=nums[0]; cfg.sw_E_Vm[1]=nums[1]; cfg.sw_E_Vm[2]=nums[2];
      cfg.sw_evxb = 0;
    }
    return;
  }

  // Particle count target
  if (key=="ppc" || key=="target-ppc") {
    if (!nums.empty() && nums[0]>0.0) {
      cfg.target_ppc = nums[0];
      cfg.user_target_ppc = true;
    }
    return;
  }

  // Unknown: ignore but warn (only rank 0 later; here we don't have MPI)
  std::fprintf(stderr,"[ConstFieldBC input] Warning: unknown key '%s'\n", keyRaw.c_str());
  exit(__LINE__,__FILE__,"Unknown option");
}

static void ApplyInputFile(TestConfig& cfg, const std::string& filename) {
  std::ifstream in(filename.c_str());
  if (!in.is_open()) {
    std::fprintf(stderr,"[ConstFieldBC] ERROR: cannot open input file '%s'\n", filename.c_str());
    return;
  }

  std::string line;
  int lineno=0;
  while (std::getline(in, line)) {
    ++lineno;
    line = StripComments(line);
    line = trim(line);
    if (line.empty()) continue;

    // Expect key=value
    size_t eq = line.find('=');
    if (eq==std::string::npos) continue;

    std::string key = trim(line.substr(0,eq));
    std::string rhs = trim(line.substr(eq+1));

    // Allow empty rhs for boolean flags: "sw-no-round="
    std::vector<double> nums;
    std::vector<std::string> strs;
    if (!rhs.empty()) ParseDoublesCSV(rhs, nums, strs);

    ApplyKeyValue(cfg, key, nums, strs);
  }
}

// -----------------------------------------------------------------------------
// Help
// -----------------------------------------------------------------------------
void PrintHelpAndExit(const char* prog) {
  std::printf(
    "\n"
    "ECSIM Constant-Field / Boundary-Condition (BC) Test\n"
    "---------------------------------------------------\n"
    "Regression/validation driver for the ECSIM field solver and BC handling.\n"
    "Use it to catch regressions in matrix/RHS assembly, BC short-circuit rows\n"
    "in GetStencil(), and MPI corner/center exchange.\n"
    "\n"
    "Configuration precedence:\n"
    "  1) Defaults\n"
    "  2) Input file (-i <file>)\n"
    "  3) Command line (overrides input file)\n"
    "\n"
    "Input file format:\n"
    "  key=value1,value2,value3,...\n"
    "  • comments allowed with '#', ';', or '//' anywhere on a line\n"
    "  • examples:  sw-ncm3=5   ; sw-TK=1e5   ; sw-u=0.05,0,0\n"
    "\n"
    "Usage:\n"
    "  %s [options]\n"
    "\n"
    "General:\n"
    "  -h, --help\n"
    "      Print this help.\n"
    "  -i <file>\n"
    "      Read parameters from an input file (key=value,...). CLI options override.\n"
    "\n"
    "Modes:\n"
    "  -no-particles\n"
    "      Field-only test (J=0, rho=0). Uniform E/B should remain uniform.\n"
    "  -particles\n"
    "      Enable particles. With -sw* options, initializes a uniform solar-wind-like\n"
    "      plasma (ions+electrons) with bulk drift and thermal pressure.\n"
    "\n"
    "Domain:\n"
    "  -L  L\n"
    "      Set cubic domain size L, centered at (0,0,0).\n"
    "  -L  Lx Ly Lz\n"
    "      Set domain size (Lx,Ly,Lz), centered at (0,0,0).\n"
    "\n"
    "Background fields (any mode):\n"
    "  -B  Bx By Bz\n"
    "      Set uniform background magnetic field (code units of this test).\n"
    "  -E  Ex Ey Ez\n"
    "      Set uniform background electric field  (code units of this test).\n"
    "  -sw-BnT Bx By Bz\n"
    "      Set uniform background B in nT (converted internally; applied as -B).\n"
    "\n"
    "Solar-wind plasma IC (implies -particles):\n"
    "  -sw | -solar-wind\n"
    "      Enable particles with default solar-wind-like parameters.\n"
    "  -sw-rho RHO\n"
    "      Background mass density parameter (driver code units).\n"
    "  -sw-p   P\n"
    "      Background scalar pressure parameter (driver code units).\n"
    "  -sw-u  ux uy uz\n"
    "      Bulk flow velocity (driver code units).\n"
    "  -ppc N\n"
    "      Target macro-particles per cell per species for the uniform IC (default 100).\n"
    "      The code sets a global particle weight so initial injection produces ~N ppc.\n"
    "  -sw-no-round\n"
    "      Disable stochastic rounding of particles-per-cell.\n"
    "\n"
    "Solar-wind plasma IC (physical inputs):\n"
    "  -sw-ncm3 N\n"
    "      Number density in cm^-3 (converted internally).\n"
    "  -sw-TK   T\n"
    "      Temperature in K (assumes Ti=Te=T; sets pressure from n,T).\n"
    "\n"
    "Options:\n"
    "  -stencil-order=N\n"
    "      FD stencil order for this test (default 2).\n"
    "\n"
    "Domain electromagnetic BC (global + per-face overrides):\n"
    "  -bc dirichlet|neumann\n"
    "      Set DomainBC type on all 6 faces (default: dirichlet).\n"
    "      Input-file keys (global): bc=..., domain-bc=..., bc-type=...\n"
    "\n"
    "  --bc-xmin dirichlet|neumann   --bc-xmax dirichlet|neumann\n"
    "  --bc-ymin dirichlet|neumann   --bc-ymax dirichlet|neumann\n"
    "  --bc-zmin dirichlet|neumann   --bc-zmax dirichlet|neumann\n"
    "      Override the BC on an individual face. These take precedence over -bc.\n"
    "      Input-file keys (per-face):\n"
    "        bc-xmin=..., bc-xmax=..., bc-ymin=..., bc-ymax=..., bc-zmin=..., bc-zmax=...\n"
    "\n"
    "      Precedence (most specific wins):\n"
    "        per-face > global, and CLI > input file > defaults\n"
    "\n"
    "Examples:\n"
    "  Field-only constant B:\n"
    "    %s -no-particles -B 0 5e-9 0\n"
    "\n"
    "  Solar wind (physical units) + constant B, ~100 ppc/spec:\n"
    "    %s -particles -sw-ncm3 5 -sw-TK 1e5 -sw-BnT 0 5 0 -sw-u 0.05 0 0 -ppc 100\n"
    "\n"
    "  Same, but driven by an input file (CLI overrides):\n"
    "    %s -i sw.in -ppc 200\n"
    "\n"
    "  Long box centered at origin:\n"
    "    %s -particles -L 128 16 16 -sw-ncm3 5 -sw-TK 1e5 -sw-BnT 0 5 0 -ppc 100\n"
    "\n",
    prog, prog, prog, prog, prog);
  std::exit(0);
}

static bool TryRead3(int& i, int argc, char** argv, double v[3]) {
  if (i + 3 >= argc) return false;
  char* end=nullptr;
  for (int k=0; k<3; ++k) {
    end = nullptr;
    v[k] = std::strtod(argv[i+1+k], &end);
    if (end==argv[i+1+k] || !end) return false;
  }
  i += 3;
  return true;
}

static bool TryRead1or3(int& i, int argc, char** argv, double v[3]) {
  if (i + 1 >= argc) return false;

  auto parse_double = [](const char* s, double& out)->bool {
    char* end=nullptr;
    out = std::strtod(s, &end);
    if (end==s) return false;               // no conversion
    while (*end==' ' || *end=='\t') ++end;  // tolerate trailing whitespace
    return (*end=='\0');                    // must consume full token
  };

  double a0=0.0;
  if (!parse_double(argv[i+1], a0)) return false;

  // Try to read two more doubles. If both parse cleanly, interpret as 3-vector.
  double a1=0.0, a2=0.0;
  bool have3 = false;
  if (i + 3 < argc) {
    if (parse_double(argv[i+2], a1) && parse_double(argv[i+3], a2)) {
      have3 = true;
    }
  }

  if (have3) {
    v[0]=a0; v[1]=a1; v[2]=a2;
    i += 3;
  }
  else {
    v[0]=a0; v[1]=a0; v[2]=a0;
    i += 1;
  }
  return true;
}


static bool ParseIntAfterEqOrNext(int& i, int argc, char** argv, const char* /*opt*/, int& outVal) {
  std::string a(argv[i]);
  auto pos = a.find('=');
  if (pos != std::string::npos) {
    const char* s = a.c_str() + pos + 1;
    char* end = nullptr;
    long v = std::strtol(s, &end, 10);
    if (end && end != s) { outVal = static_cast<int>(v); return true; }
    return false;
  } else {
    if (i + 1 >= argc) return false;
    char* end = nullptr;
    long v = std::strtol(argv[i+1], &end, 10);
    if (end && end != argv[i+1]) { outVal = static_cast<int>(v); ++i; return true; }
    return false;
  }
}

void ConfigureTestFromArgs(TestConfig& cfg,int argc, char** argv) {
  for (int i=1; i<argc; ++i) {
    std::string a(argv[i]);

    // Input file option is handled in ConfigureTestFromArgsWithInput(). Skip it here.
    if (a=="-i" || a.rfind("-i=",0)==0) { if (a=="-i" && i+1<argc) ++i; continue; }

    // ---- Mode selection ----
    if (a=="-particles") {
      cfg.mode = TestConfig::Mode::WithParticles;
      continue;
    }
    if (a=="-no-particles") {
      cfg.mode = TestConfig::Mode::NoParticles;
      continue;
    }

    
// ---- Domain boundary conditions ----
// Applies to all faces (Dirichlet or Neumann). Default is Dirichlet.
if (a=="-bc" || a=="--bc" || a.rfind("-bc=",0)==0 || a.rfind("--bc=",0)==0) {
  std::string v;
  if (a.rfind("=",0)!=std::string::npos) {
    v = a.substr(a.find('=')+1);
  }
  else {
    if (i+1>=argc) { std::fprintf(stderr,"[ConstFieldBC] ERROR: -bc requires a value (dirichlet|neumann)\n"); std::exit(1); }
    v = argv[++i];
  }
  v = tolower_str(trim(v));
  cfg.user_domain_bc = true;
  if (v=="d" || v=="dirichlet") cfg.domain_bc = TestConfig::DomainBCType::Dirichlet;
  else if (v=="n" || v=="neumann") cfg.domain_bc = TestConfig::DomainBCType::Neumann;
  else {
    std::fprintf(stderr,"[ConstFieldBC] ERROR: unknown BC type '%s' (expected dirichlet|neumann)\n", v.c_str());
    std::exit(1);
  }
  continue;
}

// ---- Per-face domain boundary conditions ----
// CLI overrides input file. Per-face overrides global -bc.
//   --bc-xmin dirichlet|neumann   (also accepts --bc-xmin=... and -bc-xmin/ -bc-xmin=...)
// Face index convention: 0:xmin,1:xmax,2:ymin,3:ymax,4:zmin,5:zmax.
{
  auto parse_bc_str = [&](const std::string& raw)->TestConfig::DomainBCType {
    std::string v = tolower_str(trim(raw));
    if (v=="d" || v=="dirichlet") return TestConfig::DomainBCType::Dirichlet;
    if (v=="n" || v=="neumann")   return TestConfig::DomainBCType::Neumann;
    std::fprintf(stderr,"[ConstFieldBC] ERROR: unknown BC type '%s' (expected dirichlet|neumann)\n", v.c_str());
    std::exit(1);
  };

  auto parse_face_opt = [&](const std::string& opt, int faceIdx)->bool {
    if (a==opt || a==("-"+opt) || a==("--"+opt) || a.rfind(opt+"=",0)==0 || a.rfind("-"+opt+"=",0)==0 || a.rfind("--"+opt+"=",0)==0) {
      std::string v;
      auto pos = a.find('=');
      if (pos!=std::string::npos) v = a.substr(pos+1);
      else {
        if (i+1>=argc) { std::fprintf(stderr,"[ConstFieldBC] ERROR: %s requires a value (dirichlet|neumann)\n", opt.c_str()); std::exit(1); }
        v = argv[++i];
      }
      cfg.user_domain_bc_face[faceIdx] = true;
      cfg.domain_bc_face[faceIdx] = parse_bc_str(v);
      return true;
    }
    return false;
  };

  if (parse_face_opt("bc-xmin",0)) continue;
  if (parse_face_opt("bc-xmax",1)) continue;
  if (parse_face_opt("bc-ymin",2)) continue;
  if (parse_face_opt("bc-ymax",3)) continue;
  if (parse_face_opt("bc-zmin",4)) continue;
  if (parse_face_opt("bc-zmax",5)) continue;
}
// ---- Particle (solar wind) controls ----
    // These options imply particles are enabled.
    if (a=="-sw" || a=="-solar-wind") {
      cfg.mode = TestConfig::Mode::WithParticles;
      // keep defaults sw_rho0/sw_p0/sw_u0
      continue;
    }
    if (a=="-sw-rho") {
      if (i+1>=argc) { std::printf("-sw-rho requires a value\n"); continue; }
      cfg.mode = TestConfig::Mode::WithParticles;
      cfg.sw_rho0 = std::atof(argv[++i]);
      continue;
    }
    if (a=="-sw-p") {
      if (i+1>=argc) { std::printf("-sw-p requires a value\n"); continue; }
      cfg.mode = TestConfig::Mode::WithParticles;
      cfg.sw_p0 = std::atof(argv[++i]);
      continue;
    }
    if (a=="-sw-u") {
      cfg.mode = TestConfig::Mode::WithParticles;
      if (!TryRead3(i, argc, argv, cfg.sw_u0)) {
        std::printf("-sw-u requires three values: ux uy uz\n");
      }
      continue;
    }
    if (a=="-sw-no-round") {
      cfg.sw_use_rounding = false;
      continue;
    }

    // ---- Unit normalization scales (SI) for converting physical inputs -> solver units ----
    // Accept either "-units-lm V" or "-units-lm=V" (similarly for -units-umps / -units-mkg).
    auto ParseDoubleAfterEqOrNext = [&](const std::string& opt, double& outv)->bool {
      if (a.rfind(opt+"=",0)==0) {
        const char* s = a.c_str() + (opt.size()+1);
        char* end=nullptr;
        outv = std::strtod(s,&end);
        return end!=s;
      }
      if (a==opt) {
        if (i+1>=argc) return false;
        const char* s = argv[++i];
        char* end=nullptr;
        outv = std::strtod(s,&end);
        return end!=s;
      }
      return false;
    };

    if (a=="-units-lm" || a.rfind("-units-lm=",0)==0) {
      double v;
      if (!ParseDoubleAfterEqOrNext("-units-lm", v) || v<=0.0) {
        std::printf("-units-lm requires a positive value in meters\n");
      } else { cfg.units_user_set=true; cfg.units_lSI_m=v; }
      continue;
    }
    if (a=="-units-umps" || a.rfind("-units-umps=",0)==0) {
      double v;
      if (!ParseDoubleAfterEqOrNext("-units-umps", v) || v<=0.0) {
        std::printf("-units-umps requires a positive value in m/s\n");
      } else { cfg.units_user_set=true; cfg.units_uSI_mps=v; }
      continue;
    }
    if (a=="-units-mkg" || a.rfind("-units-mkg=",0)==0) {
      double v;
      if (!ParseDoubleAfterEqOrNext("-units-mkg", v) || v<=0.0) {
        std::printf("-units-mkg requires a positive value in kg\n");
      } else { cfg.units_user_set=true; cfg.units_mSI_kg=v; }
      continue;
    }

    if (a=="-ppc" || a.rfind("-ppc=",0)==0) {
      // Target macro-particles per cell (per species) for the uniform IC.
      // Accept "-ppc N" or "-ppc=N".
      const char* s = nullptr;
      if (a=="-ppc") {
        if (i+1>=argc) { std::printf("-ppc requires a value\n"); continue; }
        s = argv[++i];
      }
      else {
        s = a.c_str() + 5; // after "-ppc="
      }

      char* end=nullptr;
      double v = std::strtod(s,&end);
      if (end==s) { std::printf("Invalid -ppc value: %s\n", s); continue; }
      if (v<=0.0) { std::printf("-ppc must be > 0\n"); continue; }

      cfg.target_ppc = v;
      cfg.user_target_ppc = true;
      continue;
    }


    if (a=="-sw-ncm3") {
      if (i+1>=argc) { std::printf("-sw-ncm3 requires a value\n"); continue; }
      cfg.mode = TestConfig::Mode::WithParticles;
      cfg.sw_has_ncm3 = true;
      cfg.sw_n_cm3 = std::atof(argv[++i]);
      continue;
    }
    if (a=="-sw-TK") {
      if (i+1>=argc) { std::printf("-sw-TK requires a value\n"); continue; }
      cfg.mode = TestConfig::Mode::WithParticles;
      cfg.sw_has_TK = true;
      cfg.sw_TK = std::atof(argv[++i]);
      continue;
    }
    if (a=="-sw-BnT") {
      // Background magnetic field in nT (physical). Stored as SI and converted later.
      cfg.sw_has_BnT = true;

      double bnT[3] = {0.0,0.0,0.0};
      bool ok = TryRead3(i, argc, argv, bnT);
      if (!ok) { bnT[0]=0.0; bnT[1]=5.0; bnT[2]=0.0; } // default 5 nT along +y
      cfg.sw_BnT[0]=bnT[0]; cfg.sw_BnT[1]=bnT[1]; cfg.sw_BnT[2]=bnT[2];

      // Also treat as a background-field specification in SI (Tesla).
      cfg.userB_SI = true;
      cfg.B0_SI_T[0] = bnT[0]*1.0e-9;
      cfg.B0_SI_T[1] = bnT[1]*1.0e-9;
      cfg.B0_SI_T[2] = bnT[2]*1.0e-9;

      continue;
    }

    if (a=="-sw-u-kms" || a=="-sw-ukms") {
      cfg.mode = TestConfig::Mode::WithParticles;
      cfg.sw_has_u_kms = true;
      if (!TryRead3(i, argc, argv, cfg.sw_u_kms)) {
        std::printf("-sw-u-kms requires three values: ux uy uz [km/s]\n");
      }
      continue;
    }

    if (a=="-sw-u-mps" || a=="-sw-umps") {
      cfg.mode = TestConfig::Mode::WithParticles;
      cfg.sw_has_u_mps = true;
      if (!TryRead3(i, argc, argv, cfg.sw_u_mps)) {
        std::printf("-sw-u-mps requires three values: ux uy uz [m/s]\n");
      }
      continue;
    }

    if (a.rfind("-sw-evxb",0)==0) {
      // Enable/disable convective electric field: E = -u x B (physical SI, then converted).
      int val = 1;
      if (!ParseIntAfterEqOrNext(i, argc, argv, "-sw-evxb", val)) {
        std::printf("Invalid -sw-evxb value; using default 1\n");
        val = 1;
      }
      cfg.mode = TestConfig::Mode::WithParticles;
      cfg.sw_evxb = (val!=0) ? 1 : 0;
      continue;
    }

    if (a=="-sw-EmVm") {
      cfg.mode = TestConfig::Mode::WithParticles;
      cfg.sw_has_EmVm = true;
      if (!TryRead3(i, argc, argv, cfg.sw_E_mVm)) {
        std::printf("-sw-EmVm requires three values: Ex Ey Ez [mV/m]\n");
      }
      cfg.sw_evxb = 0; // explicit E overrides E=-u×B
      continue;
    }

    if (a=="-sw-EVm") {
      cfg.mode = TestConfig::Mode::WithParticles;
      cfg.sw_has_EVm = true;
      if (!TryRead3(i, argc, argv, cfg.sw_E_Vm)) {
        std::printf("-sw-EVm requires three values: Ex Ey Ez [V/m]\n");
      }
      cfg.sw_evxb = 0;
      continue;
    }



// ---- Domain size ----
// -L L     -> cubic domain of size L, centered at (0,0,0)
// -L Lx Ly Lz -> anisotropic domain, centered at (0,0,0)
if (a=="-L") {
  cfg.use_domain_L = true;
  if (!TryRead1or3(i, argc, argv, cfg.domain_L)) {
    std::printf("-L requires 1 value (L) or 3 values (Lx Ly Lz)\n");
    cfg.use_domain_L = false;
  }
  continue;
}


    // ---- Background fields in *physical* units (SI) ----
    // These do NOT force particle mode; they are convenience setters.
    if (a=="-BnT") {
      double bnT[3] = {0.0,0.0,0.0};
      bool ok = TryRead3(i, argc, argv, bnT);
      if (!ok) { bnT[0]=0.0; bnT[1]=5.0; bnT[2]=0.0; }
      cfg.userB_SI = true;
      cfg.B0_SI_T[0]=bnT[0]*1.0e-9; cfg.B0_SI_T[1]=bnT[1]*1.0e-9; cfg.B0_SI_T[2]=bnT[2]*1.0e-9;
      if (cfg.mode != TestConfig::Mode::WithParticles) cfg.mode = TestConfig::Mode::FieldOnlyB;
      continue;
    }
    if (a=="-EVm") {
      double Evm[3] = {0.0,0.0,0.0};
      bool ok = TryRead3(i, argc, argv, Evm);
      if (!ok) { Evm[0]=1.0; Evm[1]=0.0; Evm[2]=0.0; }
      cfg.userE_SI = true;
      cfg.E0_SI_Vm[0]=Evm[0]; cfg.E0_SI_Vm[1]=Evm[1]; cfg.E0_SI_Vm[2]=Evm[2];
      if (cfg.mode != TestConfig::Mode::WithParticles) cfg.mode = TestConfig::Mode::FieldOnlyE;
      continue;
    }
    if (a=="-EmVm") {
      double EmVm[3] = {0.0,0.0,0.0};
      bool ok = TryRead3(i, argc, argv, EmVm);
      if (!ok) { EmVm[0]=0.0; EmVm[1]=0.0; EmVm[2]=0.0; }
      cfg.userE_SI = true;
      cfg.E0_SI_Vm[0]=EmVm[0]*1.0e-3; cfg.E0_SI_Vm[1]=EmVm[1]*1.0e-3; cfg.E0_SI_Vm[2]=EmVm[2]*1.0e-3;
      if (cfg.mode != TestConfig::Mode::WithParticles) cfg.mode = TestConfig::Mode::FieldOnlyE;
      continue;
    }
    // ---- Background fields ----
    // In particle mode, -B/-E set background fields but DO NOT switch to field-only.
    // In non-particle modes, -B/-E imply field-only initialization.
    if (a=="-B") {
      cfg.userB = TryRead3(i, argc, argv, cfg.B0);
      if (!cfg.userB) { cfg.B0[0]=0.0; cfg.B0[1]=1.0; cfg.B0[2]=0.0; }
      if (cfg.mode != TestConfig::Mode::WithParticles) cfg.mode = TestConfig::Mode::FieldOnlyB;
      continue;
    }
    if (a=="-E") {
      cfg.userE = TryRead3(i, argc, argv, cfg.E0);
      if (!cfg.userE) { cfg.E0[0]=1.0; cfg.E0[1]=0.0; cfg.E0[2]=0.0; }
      if (cfg.mode != TestConfig::Mode::WithParticles) cfg.mode = TestConfig::Mode::FieldOnlyE;
      continue;
    }

    // ---- Operator stencil order (if your solver reads g_TestStencilOrder) ----
    if (a.rfind("-stencil-order",0)==0) {
      int val = cfg.stencilOrder;
      if (!ParseIntAfterEqOrNext(i, argc, argv, "-stencil-order", val)) {
        std::printf("Invalid -stencil-order value; keeping default %d\n", cfg.stencilOrder);
      } else {
        cfg.stencilOrder = val;
      }
      continue;
    }


if (a=="-h" || a=="--help") { PrintHelpAndExit(argv[0]); }


    std::printf("Unknown option: %s (use -h for help)\n", argv[i]);
    exit(__LINE__,__FILE__,"Unknown option");
  }
}



// Wrapper: parse input file first (if -i is provided), then parse the CLI (CLI overrides).
void ConfigureTestFromArgsWithInput(TestConfig& cfg, int argc, char** argv) {
  cfg = TestConfig(); // defaults

  // 1) scan for -i / -i=FILE
  std::string infile;
  for (int i=1; i<argc; ++i) {
    std::string a(argv[i]);
    if (a=="-i") {
      if (i+1<argc) infile = argv[++i];
    }
    else if (a.rfind("-i=",0)==0) {
      infile = a.substr(3);
    }
  }

  if (!infile.empty()) {
    cfg.inputFile = infile;
    ApplyInputFile(cfg, infile);
  }

  // 2) apply CLI options on top
  ConfigureTestFromArgs(cfg, argc, argv);
}
