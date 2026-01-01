#ifndef _CONSTFIELDBC_CLI_H_
#define _CONSTFIELDBC_CLI_H_

#include <string>

// ---------------- CLI configuration ----------------
struct TestConfig {
  enum class Mode { WithParticles, NoParticles, FieldOnlyB, FieldOnlyE } mode = Mode::WithParticles;

  bool userB = false, userE = false;
  double B0[3] = {0.0, 0.0, 0.0};
  double E0[3] = {0.0, 0.0, 0.0};
  int stencilOrder = 2;

  // Domain size control (optional). If -L is provided, the domain is centered at (0,0,0)
  // with extents [-Lx/2,Lx/2] etc.
  bool   use_domain_L = false;
  double domain_L[3] = {32.0,16.0,8.0}; // defaults match legacy xmin/xmax

  // Particle (solar wind) initialization controls (used when mode==WithParticles).
  // NOTE: sw_rho0/sw_p0 are in the same normalized convention used by this legacy test
  //       (see rho_conv/p_conv in PrepopulateDomain()).
  double sw_rho0 = 1.0;              // background mass density parameter (before rho_conv)
  double sw_p0   = 4.5e-4;           // background scalar pressure parameter (before p_conv)
  double sw_u0[3] = {0.05,0.0,0.0};  // bulk flow velocity (test units)
  bool   sw_use_rounding = true;     // stochastic rounding for particle counts

  // Optional physical-unit convenience inputs (interpreted in SI):
  //   -sw-ncm3 : number density in cm^-3 (converted to m^-3)
  //   -sw-TK   : temperature in K (assumes Ti=Te=T)
  //   -sw-BnT  : background magnetic field in nT (converted to Tesla and stored in B0)
  bool   sw_has_ncm3 = false;
  double sw_n_cm3 = 0.0;

  bool   sw_has_TK = false;
  double sw_TK = 0.0;

  bool   sw_has_BnT = false;
  double sw_BnT[3] = {0.0,0.0,0.0};

  // Target macro-particles per cell (per species) for the uniform solar-wind-like IC.
  // Default: 100 ppc/spec. Override with -ppc <N>.
  double target_ppc = 100.0;
  bool   user_target_ppc = false;

  // Optional input file path (for bookkeeping; not required)
  std::string inputFile;
};

// Parse configuration. If '-i <file>' is present, the file is parsed first and then
// command-line options are applied on top (CLI overrides input file).
void ConfigureTestFromArgsWithInput(TestConfig& cfg, int argc, char** argv);

// For standalone help printing.
void PrintHelpAndExit(const char* prog);

#endif
