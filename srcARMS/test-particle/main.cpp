// mc_fieldline.cpp
// Monte Carlo test-particle model for electrons moving along a magnetic field line,
// with a simple guiding-center mover (motion along arc-length s, mirror force,
// optional constant E_parallel) and optional pitch-angle diffusion.
// Also writes Tecplot ASCII outputs (fieldline XY, pitch-angle distributions).
//
// Build:
//   g++ -O2 -std=c++17 mc_fieldline.cpp -o mc_fieldline
//
// -----------------------------------------------------------------------------
// INPUT FIELD LINE FORMAT
// -----------------------------------------------------------------------------
// Comment/header lines beginning with '#' are ignored.
// A bare integer line may appear giving the point count (and is ignored).
//
// Data lines support either:
//   (A) 4 columns: x y z Bmag
//   (B) 6 columns: x y z Bx By Bz        (|B| computed)
//   (C) 7+ cols :  x y z Bx By Bz Bmag   (first 7 used)
//
// Units are assumed SI (meters, Tesla), but the code does not enforce units.
//
// -----------------------------------------------------------------------------
// PHYSICS MODEL (guiding-center along s)
// -----------------------------------------------------------------------------
// Each particle state (non-relativistic):
//   s        : arc-length coordinate along the line [m]
//   v_par    : parallel velocity along +s [m/s]
//   mu_mag   : magnetic moment = m v_perp^2 / (2 B)   [J/T]  (constant in deterministic step)
//
// Deterministic guiding-center update:
//   ds/dt      = v_par
//   dv_par/dt  = (q/m) E_par  -  (mu_mag/m) dB/ds
//
// where E_par is an optional *constant* parallel electric field [V/m] (CLI -Epar_Vm).
// The mirror term uses piecewise dB/ds computed from the loaded B(s).
//
// Pitch-angle diffusion (optional, CLI -D) is applied as diffusion in mu_cos = cos(alpha):
//   d(mu_cos) = sqrt(2 D dt) * N(0,1)
// and we re-split kinetic energy between v_par and v_perp while keeping speed |v| fixed
// during the scattering substep.
//
// INITIALIZATION / INJECTION
// -----------------------------------------------------------------------------
// The code injects N electrons at s = inj_sfrac * L, where L is the total field line length.
// Velocity is sampled from a 3D Maxwellian at temperature T_eV (CLI -T_eV).
// Optionally, -mu0 can fix the initial mu_cos for all particles.
// If -forward_only 1 is used, initial v_par is forced >= 0.
//
// OUTPUTS (Tecplot ASCII)
// -----------------------------------------------------------------------------
// 1) <prefix>_fieldline.dat      : s vs abs(x,y,z), B components, |B|
// 2) <prefix>_mc_mu_hist.dat     : final mu_cos histogram PDF
// 3) <prefix>_mc_mu_vs_s.dat     : binned mean/var(mu_cos) vs s
// 4) <prefix>_mc_mu_s_2d.dat     : 2D distribution f(s, mu_cos) (counts normalized per s-bin)
//
// -----------------------------------------------------------------------------
// CLI (run with -help)
// -----------------------------------------------------------------------------

#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

static constexpr double kElectronMass   = 9.1093837015e-31; // kg
static constexpr double kElectronCharge = -1.602176634e-19; // C (electron charge is negative)
static constexpr double kBoltzmann      = 1.380649e-23;     // J/K
static constexpr double kEvToJ          = 1.602176634e-19;  // J/eV

struct Vec3 { double x=0, y=0, z=0; };

static inline Vec3 add(const Vec3& a, const Vec3& b){ return {a.x+b.x, a.y+b.y, a.z+b.z}; }
static inline Vec3 sub(const Vec3& a, const Vec3& b){ return {a.x-b.x, a.y-b.y, a.z-b.z}; }
static inline Vec3 mul(const Vec3& a, double s){ return {a.x*s, a.y*s, a.z*s}; }
static inline double dot(const Vec3& a, const Vec3& b){ return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline double norm(const Vec3& a){ return std::sqrt(dot(a,a)); }
static inline Vec3 unit_or_zero(const Vec3& a){
  double n = norm(a);
  if(n<=0.0) return {0,0,0};
  return mul(a, 1.0/n);
}
static inline double clamp(double x, double lo, double hi){ return std::min(hi, std::max(lo, x)); }

struct FieldLinePoint {
  Vec3 r;          // position [m]
  Vec3 Bvec;       // [T] optional; may be (0,0,0) if not present
  double Bmag=0.0; // |B| [T]
  double s=0.0;    // arc-length from start [m]
};

struct InterpSample {
  Vec3 r;      // position at s
  Vec3 t_hat;  // tangent unit vector dr/ds at s (used as local b-hat if B direction absent)
  Vec3 Bvec;   // interpolated Bvec (may be 0)
  double Bmag; // interpolated |B|
  double dBds; // piecewise derivative d|B|/ds
};

struct FieldLine {
  std::vector<FieldLinePoint> p;

  double length() const { return p.empty() ? 0.0 : p.back().s; }

  bool load(const std::string& path, std::string* err=nullptr){
    std::ifstream in(path);
    if(!in){
      if(err) *err = "Cannot open file: " + path;
      return false;
    }

    std::string line;
    std::vector<std::array<double,7>> rows; // max columns we care about

    while(std::getline(in, line)){
      // trim leading spaces
      size_t i=0;
      while(i<line.size() && std::isspace(static_cast<unsigned char>(line[i]))) i++;
      if(i==line.size()) continue;
      if(line[i]=='#') continue;

      // skip a bare integer line (N points)
      {
        std::istringstream iss(line);
        long long maybeN;
        if((iss >> maybeN) && iss.eof()) continue;
      }

      std::istringstream iss(line);
      std::vector<double> cols;
      double v;
      while(iss >> v) cols.push_back(v);

      if(cols.size()==4){
        // x y z B
        std::array<double,7> r{};
        r[0]=cols[0]; r[1]=cols[1]; r[2]=cols[2];
        r[3]=0; r[4]=0; r[5]=0;
        r[6]=cols[3];
        rows.push_back(r);
      } else if(cols.size()>=6){
        // x y z Bx By Bz [B]
        std::array<double,7> r{};
        r[0]=cols[0]; r[1]=cols[1]; r[2]=cols[2];
        r[3]=cols[3]; r[4]=cols[4]; r[5]=cols[5];
        if(cols.size()>=7) r[6]=cols[6];
        else r[6]=std::sqrt(r[3]*r[3]+r[4]*r[4]+r[5]*r[5]);
        rows.push_back(r);
      } else {
        // ignore malformed
        continue;
      }
    }

    if(rows.size()<2){
      if(err) *err = "Parsed <2 data rows from file: " + path;
      return false;
    }

    p.clear();
    p.reserve(rows.size());
    for(const auto& r: rows){
      FieldLinePoint q;
      q.r    = {r[0], r[1], r[2]};
      q.Bvec = {r[3], r[4], r[5]};
      q.Bmag = r[6];
      p.push_back(q);
    }

    // arc-length
    p[0].s = 0.0;
    for(size_t k=1;k<p.size();k++){
      double ds = norm(sub(p[k].r, p[k-1].r));
      p[k].s = p[k-1].s + ds;
    }

    // If Bmag looks unusable, compute it from Bvec
    bool allZero=true;
    for(const auto& q: p){ if(q.Bmag!=0.0){ allZero=false; break; } }
    if(allZero){
      for(auto& q: p) q.Bmag = norm(q.Bvec);
    }

    // If the file provides only |B| (Bvec == 0 everywhere) but |B| is present,
    // reconstruct a consistent vector field direction using the geometric tangent:
    //   Bvec = |B| * t_hat, where t_hat ~ dr/ds.
    bool bvecAllZero = true;
    bool bmagAnyPos  = false;
    for(const auto& q: p){
      if(norm(q.Bvec) > 0.0) bvecAllZero = false;
      if(q.Bmag > 0.0) bmagAnyPos = true;
    }

    if(bvecAllZero && bmagAnyPos){
      const size_t n = p.size();
      for(size_t k=0;k<n;k++){
        Vec3 t{0,0,0};
        if(n>=2){
          if(k==0)          t = sub(p[1].r, p[0].r);
          else if(k==n-1)   t = sub(p[n-1].r, p[n-2].r);
          else {
            // centered difference for smoother direction
            t = sub(p[k+1].r, p[k-1].r);
          }
        }
        Vec3 th = unit_or_zero(t);
        p[k].Bvec = mul(th, p[k].Bmag);
      }
    }

    // Guard against zeros
    for(auto& q: p){
      if(!(q.Bmag>0.0)) q.Bmag = 0.0;
    }

    return true;
  }

  // Find bracket segment for s_query (returns i such that p[i].s <= s < p[i+1].s)
  size_t bracket(double s_query) const {
    if(s_query<=p.front().s) return 0;
    if(s_query>=p.back().s) return p.size()-2;
    auto it = std::lower_bound(p.begin(), p.end(), s_query,
      [](const FieldLinePoint& a, double s){ return a.s < s; });
    size_t j = static_cast<size_t>(it - p.begin());
    if(j==0) return 0;
    if(j>=p.size()) return p.size()-2;
    return j-1;
  }

  InterpSample sample(double s_query) const {
    InterpSample out{};
    if(p.size()<2){
      return out;
    }

    // Clamp s into range
    double s = clamp(s_query, p.front().s, p.back().s);
    size_t i = bracket(s);
    size_t j = i+1;

    double s0=p[i].s, s1=p[j].s;
    double t = (s1>s0) ? (s - s0)/(s1 - s0) : 0.0;

    out.r    = add(mul(p[i].r, 1.0-t), mul(p[j].r, t));
    out.Bvec = add(mul(p[i].Bvec,1.0-t), mul(p[j].Bvec,t));
    out.Bmag = p[i].Bmag*(1.0-t) + p[j].Bmag*t;

    // Tangent along the segment
    Vec3 dr = sub(p[j].r, p[i].r);
    double ds = (s1 - s0);
    out.t_hat = (ds>0.0) ? unit_or_zero(dr) : Vec3{0,0,0};

    // Piecewise dB/ds on the same segment
    if(ds>0.0) out.dBds = (p[j].Bmag - p[i].Bmag)/ds;
    else out.dBds = 0.0;

    // Robustness: ensure Bmag positive-ish
    if(!(out.Bmag>0.0)) out.Bmag = 0.0;

    return out;
  }

  // Tecplot XY exporter: X axis = s (relative distance), Y variables = abs(x,y,z,Bx,By,Bz,B)
  bool write_tecplot_xy(const std::string& out_path, std::string* err=nullptr) const {
    std::ofstream out(out_path);
    if(!out){
      if(err) *err = "Cannot write: " + out_path;
      return false;
    }

    out << "TITLE = \"Field line XY export\"\n";
    out << "VARIABLES = \"s_m\" \"absX_m\" \"absY_m\" \"absZ_m\" \"Bx_T\" \"By_T\" \"Bz_T\" \"B_T\"\n";
    out << "ZONE T=\"fieldline\", I=" << p.size() << ", F=POINT\n";
    out << std::setprecision(12) << std::scientific;

    for(const auto& q: p){
      out << q.s << " "
          << std::fabs(q.r.x) << " " << std::fabs(q.r.y) << " " << std::fabs(q.r.z) << " "
          << q.Bvec.x << " " << q.Bvec.y << " " << q.Bvec.z << " "
          << q.Bmag << "\n";
    }
    return true;
  }
};

struct Cli {
  std::string fieldline_path;
  std::string out_prefix = "out";

  // Injection / plasma parameters
  double inj_sfrac = 0.0;     // fraction of field line length [0..1]
  double T_eV = 100.0;        // injected electron temperature [eV]
  double Epar_Vm = 0.0;       // constant parallel electric field [V/m]

  // Monte Carlo / numerics
  int npart = 100000;
  double dt = 1e-3;           // [s]
  double tmax = 1.0;          // [s]
  double D = 0.0;             // [1/s] pitch-angle diffusion coefficient

  bool mu0_provided = false;
  double mu0 = 0.0;           // fixed initial mu_cos = cos(alpha)
  bool forward_only = false;

  uint64_t seed = 1;

  int nbins_mu = 80;
  int nbins_s  = 80;

  static void print_help(const char* exe){
    std::cout
      << "Usage:\n"
      << "  " << exe << " -fieldline <path> [options]\n\n"
      << "Required:\n"
      << "  -fieldline <path>         Input field line file.\n"
      << "                            Supports columns:\n"
      << "                              x y z B\n"
      << "                              x y z Bx By Bz\n"
      << "                              x y z Bx By Bz B\n\n"
      << "Key options (requested enhancements):\n"
      << "  -inj_sfrac <0..1>          Injection location as fraction of total field-line length.\n"
      << "                            0.0=start, 0.5=middle, 1.0=end (default: 0.0)\n"
      << "  -T_eV <eV>                 Injected electron temperature (Maxwellian) (default: 100 eV)\n\n"
      << "Guiding-center / acceleration:\n"
      << "  -Epar_Vm <V/m>             Constant parallel electric field along +s (default: 0)\n\n"
      << "Monte Carlo / numerics:\n"
      << "  -npart <int>               Number of particles (default: 100000)\n"
      << "  -dt <s>                    Time step (default: 1e-3)\n"
      << "  -tmax <s>                  Total time (default: 1.0)\n"
      << "  -D <1/s>                   Pitch-angle diffusion coeff for mu=cos(alpha) (default: 0)\n"
      << "  -mu0 <value>               Fix initial mu for all particles (overrides isotropic Maxwellian angles)\n"
      << "  -forward_only <0|1>         Force initial v_par >= 0 (default: 0)\n"
      << "  -seed <int>                RNG seed (default: 1)\n"
      << "  -nbins_mu <int>            mu histogram bins (default: 80)\n"
      << "  -nbins_s <int>             s bins for statistics and 2D distribution (default: 80)\n"
      << "  -out <prefix>              Output prefix (default: out)\n"
      << "  -help                      Show this help and exit\n\n"
      << "Outputs (Tecplot ASCII):\n"
      << "  <prefix>_fieldline.dat     Field line XY: s vs abs(x,y,z), B components and |B|\n"
      << "  <prefix>_mc_mu_hist.dat    Final mu histogram (PDF)\n"
      << "  <prefix>_mc_mu_vs_s.dat    mean/var(mu) binned vs s\n"
      << "  <prefix>_mc_mu_s_2d.dat    2D distribution f(s,mu) normalized per s-bin\n";
  }

  bool parse(int argc, char** argv, std::string* err=nullptr){
    if(argc<=1){
      if(err) *err = "No arguments. Use -help.";
      return false;
    }

    for(int i=1;i<argc;i++){
      std::string a = argv[i];
      auto need = [&](const std::string& opt)->std::string{
        if(i+1>=argc) throw std::runtime_error("Missing value after " + opt);
        return std::string(argv[++i]);
      };

      try{
        if(a=="-help" || a=="--help" || a=="-h"){
          print_help(argv[0]);
          std::exit(0);
        } else if(a=="-fieldline"){
          fieldline_path = need(a);
        } else if(a=="-out"){
          out_prefix = need(a);
        } else if(a=="-inj_sfrac"){
          inj_sfrac = std::stod(need(a));
        } else if(a=="-T_eV"){
          T_eV = std::stod(need(a));
        } else if(a=="-Epar_Vm"){
          Epar_Vm = std::stod(need(a));
        } else if(a=="-npart"){
          npart = std::stoi(need(a));
        } else if(a=="-dt"){
          dt = std::stod(need(a));
        } else if(a=="-tmax"){
          tmax = std::stod(need(a));
        } else if(a=="-D"){
          D = std::stod(need(a));
        } else if(a=="-mu0"){
          mu0 = std::stod(need(a));
          mu0_provided = true;
        } else if(a=="-forward_only"){
          forward_only = (std::stoi(need(a))!=0);
        } else if(a=="-seed"){
          seed = static_cast<uint64_t>(std::stoull(need(a)));
        } else if(a=="-nbins_mu"){
          nbins_mu = std::max(4, std::stoi(need(a)));
        } else if(a=="-nbins_s"){
          nbins_s = std::max(4, std::stoi(need(a)));
        } else {
          throw std::runtime_error("Unknown option: " + a);
        }
      } catch(const std::exception& e){
        if(err) *err = e.what();
        return false;
      }
    }

    if(fieldline_path.empty()){
      if(err) *err = "Missing -fieldline <path>.";
      return false;
    }
    if(!(inj_sfrac>=0.0 && inj_sfrac<=1.0)){
      if(err) *err = "-inj_sfrac must be in [0,1].";
      return false;
    }
    if(!(T_eV>0.0)){
      if(err) *err = "-T_eV must be > 0.";
      return false;
    }
    if(npart<=0){
      if(err) *err = "-npart must be > 0.";
      return false;
    }
    if(!(dt>0.0 && tmax>0.0)){
      if(err) *err = "-dt and -tmax must be > 0.";
      return false;
    }
    if(D<0.0){
      if(err) *err = "-D must be >= 0.";
      return false;
    }
    if(mu0_provided && (mu0<-1.0 || mu0>1.0)){
      if(err) *err = "-mu0 must be in [-1,1].";
      return false;
    }
    return true;
  }
};

struct Particle {
  double s=0.0;       // [m]
  double vpar=0.0;    // [m/s]
  double mu_mag=0.0;  // [J/T] = m v_perp^2 /(2B)
  bool alive=true;    // particle is on the field line
};

struct MCResult {
  std::vector<double> s_final;
  std::vector<double> mu_final; // mu_cos
};

static inline double vth_sigma_from_Te_eV(double T_eV){
  // For Maxwellian, each velocity component is N(0, sigma^2) with sigma^2 = kT/m.
  // Convert eV -> Joules using kT = T_eV * e (since 1 eV = e Joules).
  const double kT = T_eV * kEvToJ; // J
  return std::sqrt(kT / kElectronMass);
}

static inline Vec3 sample_maxwellian_velocity(std::mt19937_64& rng, double sigma){
  std::normal_distribution<double> gauss(0.0, sigma);
  return {gauss(rng), gauss(rng), gauss(rng)};
}

// For mu diffusion, keep speed constant (during scattering), update (vpar, mu_mag) at local B.
static inline void apply_pitch_angle_diffusion(Particle& part,
                                               double B_local,
                                               double dt,
                                               double D_mumu,
                                               std::mt19937_64& rng)
{
  if(D_mumu<=0.0) return;

  // Compute current v_perp from mu_mag and B
  const double vperp2 = std::max(0.0, 2.0*part.mu_mag*B_local / kElectronMass);
  const double v2 = part.vpar*part.vpar + vperp2;
  const double v  = std::sqrt(std::max(0.0, v2));
  if(v<=0.0) return;

  double mu_cos = part.vpar / v;
  mu_cos = clamp(mu_cos, -1.0, 1.0);

  std::normal_distribution<double> gauss01(0.0, 1.0);
  const double sigma = std::sqrt(2.0*D_mumu*dt);
  mu_cos = clamp(mu_cos + sigma*gauss01(rng), -1.0, 1.0);

  // Re-split speed between parallel/perp at *fixed* |v|
  part.vpar = mu_cos * v;
  const double vperp2_new = (1.0 - mu_cos*mu_cos) * v2;

  // Update mu_mag consistent with new v_perp at the current B.
  const double B = std::max(B_local, 1e-30);
  part.mu_mag = 0.5*kElectronMass*vperp2_new / B;
}

static MCResult run_guiding_center_mc(const FieldLine& fl, const Cli& cli){
  std::mt19937_64 rng(cli.seed);

  const double L = fl.length();
  const double s0 = cli.inj_sfrac * L;

  const int nsteps = static_cast<int>(std::ceil(cli.tmax / cli.dt));

  // Initialize particles injected at s0
  std::vector<Particle> part(static_cast<size_t>(cli.npart));

  // Local tangent direction (used as b-hat for defining pitch angle at injection)
  InterpSample samp0 = fl.sample(s0);
  Vec3 bhat = samp0.t_hat;
  if(norm(bhat)<=0.0){
    // fallback if tangent is degenerate
    bhat = {1,0,0};
  }

  const double sigma = vth_sigma_from_Te_eV(cli.T_eV); // std of each component

  // Uniform distribution for mu if mu0 not specified (isotropic in mu)
  std::uniform_real_distribution<double> uni01(0.0, 1.0);

  for(int i=0;i<cli.npart;i++){
    Particle p;
    p.s = s0;

    // Sample thermal velocity vector
    Vec3 v = sample_maxwellian_velocity(rng, sigma);
    double vmag = norm(v);

    // If user fixes mu0, reorient the parallel/perp split while keeping |v|
    double mu_cos;
    if(cli.mu0_provided){
      mu_cos = clamp(cli.mu0, -1.0, 1.0);
    } else {
      if(cli.forward_only){
        mu_cos = uni01(rng);           // [0,1]
      } else {
        mu_cos = 2.0*uni01(rng) - 1.0; // [-1,1]
      }
    }

    // Define v_par from mu_cos and |v|, independent of sampled direction
    // (this avoids needing full gyro-phase and perpendicular basis).
    double vpar = mu_cos * vmag;
    if(cli.forward_only && vpar<0.0) vpar = -vpar;

    // Perpendicular speed
    double vperp2 = std::max(0.0, vmag*vmag - vpar*vpar);

    // Magnetic moment
    double B = std::max(samp0.Bmag, 1e-30);
    p.mu_mag = 0.5*kElectronMass*vperp2 / B;
    p.vpar = vpar;

    part[static_cast<size_t>(i)] = p;
  }

  // Time advance
  const double q_over_m = kElectronCharge / kElectronMass;

  for(int n=0;n<nsteps;n++){
    for(auto& p : part){
      if(!p.alive) continue;

      // Sample local field
      InterpSample sOld = fl.sample(p.s);
      const double B_old = std::max(sOld.Bmag, 1e-30);

      // Optional scattering (acts at old position)
      apply_pitch_angle_diffusion(p, B_old, cli.dt, cli.D, rng);

      // Deterministic GC update
      const double dBds = sOld.dBds;
      const double dvpar = cli.dt * ( q_over_m*cli.Epar_Vm - (p.mu_mag/kElectronMass)*dBds );
      p.vpar += dvpar;

      // Update position
      p.s += cli.dt * p.vpar;

      // Delete particle if it leaves the field line (no reflection)
      if(p.s < 0.0 || p.s > L){
        p.alive = false;
        continue;
      }

      // Update mu_mag is constant in deterministic step (already stored),
      // but B changes => v_perp changes implicitly via mu_mag = m v_perp^2 /(2B).
      // We do not need to explicitly store v_perp, but it enters mu_cos at output.
    }
  }

  // Collect outputs: final mu_cos and final s (ONLY for particles that stayed on the field line)
  MCResult out;
  out.s_final.reserve(part.size());
  out.mu_final.reserve(part.size());

  for(const Particle& p : part){
    if(!p.alive) continue;
    if(p.s < 0.0 || p.s > L) continue; // should not happen, but keep it robust

    InterpSample sf = fl.sample(p.s);
    const double B = std::max(sf.Bmag, 1e-30);
    const double vperp2 = std::max(0.0, 2.0*p.mu_mag*B / kElectronMass);
    const double v2 = p.vpar*p.vpar + vperp2;
    const double v  = std::sqrt(std::max(0.0, v2));
    double mu_cos = (v>0.0) ? (p.vpar / v) : 0.0;
    mu_cos = clamp(mu_cos, -1.0, 1.0);

    out.s_final.push_back(p.s);
    out.mu_final.push_back(mu_cos);
  }

  return out;
}

static bool write_mu_hist_tecplot(const std::string& out_path,
                                  const std::vector<double>& mu,
                                  int nbins,
                                  std::string* err=nullptr)
{
  std::vector<double> counts(static_cast<size_t>(nbins), 0.0);

  for(double m : mu){
    double x = (m + 1.0) * 0.5; // [-1,1] -> [0,1]
    int b = static_cast<int>(std::floor(x * nbins));
    if(b<0) b=0;
    if(b>=nbins) b=nbins-1;
    counts[static_cast<size_t>(b)] += 1.0;
  }

  const double N  = static_cast<double>(mu.size());
  const double dm = 2.0 / nbins;

  std::ofstream out(out_path);
  if(!out){
    if(err) *err = "Cannot write: " + out_path;
    return false;
  }

  out << "TITLE = \"Pitch-angle distribution (mu histogram)\"\n";
  out << "VARIABLES = \"mu\" \"pdf\"\n";
  out << "ZONE T=\"mu_hist\", I=" << nbins << ", F=POINT\n";
  out << std::setprecision(12) << std::scientific;

  for(int b=0;b<nbins;b++){
    const double mu_center = -1.0 + (b + 0.5)*dm;
    const double pdf = (N>0.0) ? ((counts[static_cast<size_t>(b)]/N)/dm) : 0.0;
    out << mu_center << " " << pdf << "\n";
  }
  return true;
}

static bool write_mu_vs_s_tecplot(const std::string& out_path,
                                  const std::vector<double>& s,
                                  const std::vector<double>& mu,
                                  double smax,
                                  int nbins_s,
                                  std::string* err=nullptr)
{
  std::vector<double> sum_mu(static_cast<size_t>(nbins_s),0.0);
  std::vector<double> sum_mu2(static_cast<size_t>(nbins_s),0.0);
  std::vector<double> count(static_cast<size_t>(nbins_s),0.0);

  for(size_t i=0;i<s.size();i++){
    double x = clamp(s[i]/smax, 0.0, 1.0);
    int b = static_cast<int>(std::floor(x*nbins_s));
    if(b>=nbins_s) b=nbins_s-1;
    size_t bi = static_cast<size_t>(b);
    sum_mu[bi] += mu[i];
    sum_mu2[bi] += mu[i]*mu[i];
    count[bi] += 1.0;
  }

  std::ofstream out(out_path);
  if(!out){
    if(err) *err = "Cannot write: " + out_path;
    return false;
  }

  out << "TITLE = \"Pitch-angle moments vs distance\"\n";
  out << "VARIABLES = \"s_m\" \"mu_mean\" \"mu_var\" \"count\"\n";
  out << "ZONE T=\"mu_vs_s\", I=" << nbins_s << ", F=POINT\n";
  out << std::setprecision(12) << std::scientific;

  const double ds = smax/nbins_s;
  for(int b=0;b<nbins_s;b++){
    const double s_center = (b + 0.5)*ds;
    size_t bi = static_cast<size_t>(b);
    if(count[bi]>0.0){
      const double m = sum_mu[bi]/count[bi];
      const double v = sum_mu2[bi]/count[bi] - m*m;
      out << s_center << " " << m << " " << v << " " << count[bi] << "\n";
    } else {
      out << s_center << " " << 0.0 << " " << 0.0 << " " << 0.0 << "\n";
    }
  }
  return true;
}

static bool write_mu_s_2d_tecplot(const std::string& out_path,
                                  const std::vector<double>& s,
                                  const std::vector<double>& mu,
                                  double smax,
                                  int nbins_s,
                                  int nbins_mu,
                                  std::string* err=nullptr)
{
  // counts[sbin][mubin]
  std::vector<double> counts(static_cast<size_t>(nbins_s*nbins_mu), 0.0);
  std::vector<double> count_s(static_cast<size_t>(nbins_s), 0.0);

  auto idx = [&](int bs, int bm)->size_t{
    return static_cast<size_t>(bs*nbins_mu + bm);
  };

  for(size_t i=0;i<s.size();i++){
    double xs = clamp(s[i]/smax, 0.0, 1.0);
    int bs = static_cast<int>(std::floor(xs*nbins_s));
    if(bs>=nbins_s) bs=nbins_s-1;

    double xm = (clamp(mu[i], -1.0, 1.0) + 1.0)*0.5;
    int bm = static_cast<int>(std::floor(xm*nbins_mu));
    if(bm<0) bm=0;
    if(bm>=nbins_mu) bm=nbins_mu-1;

    counts[idx(bs,bm)] += 1.0;
    count_s[static_cast<size_t>(bs)] += 1.0;
  }

  const double ds = smax/nbins_s;
  const double dm = 2.0/nbins_mu;

  std::ofstream out(out_path);
  if(!out){
    if(err) *err = "Cannot write: " + out_path;
    return false;
  }

  out << "TITLE = \"2D distribution f(s,mu)\"\n";
  out << "VARIABLES = \"s_m\" \"mu\" \"pdf\"\n";
  // Tecplot ordered zone: I varies fastest, then J.
  // We'll use I=nbins_s (s) and J=nbins_mu (mu) so that rows correspond to mu bins.
  out << "ZONE T=\"mu_s_2d\", I=" << nbins_s << ", J=" << nbins_mu << ", F=POINT\n";
  out << std::setprecision(12) << std::scientific;

  for(int jm=0;jm<nbins_mu;jm++){
    const double mu_center = -1.0 + (jm + 0.5)*dm;
    for(int is=0;is<nbins_s;is++){
      const double s_center = (is + 0.5)*ds;

      const double cs = count_s[static_cast<size_t>(is)];
      // Normalize per s-bin so that for each s-bin: sum_j pdf(s,mu_j)*dm = 1
      double pdf = 0.0;
      if(cs>0.0){
        pdf = (counts[idx(is,jm)]/cs)/dm;
      }
      out << s_center << " " << mu_center << " " << pdf << "\n";
    }
  }

  return true;
}

int main(int argc, char** argv){
  Cli cli;
  std::string err;

  if(!cli.parse(argc, argv, &err)){
    std::cerr << "Error: " << err << "\n\n";
    Cli::print_help(argv[0]);
    return 1;
  }

  FieldLine fl;
  if(!fl.load(cli.fieldline_path, &err)){
    std::cerr << "Error: " << err << "\n";
    return 2;
  }

  if(fl.p.size()<2 || !(fl.length()>0.0)){
    std::cerr << "Error: field line is too short or invalid.\n";
    return 3;
  }

  // 1) Field line Tecplot output
  {
    const std::string out_path = cli.out_prefix + "_fieldline.dat";
    if(!fl.write_tecplot_xy(out_path, &err)){
      std::cerr << "Error: " << err << "\n";
      return 4;
    }
    std::cout << "Wrote: " << out_path << "\n";
  }

  // 2) Guiding-center Monte Carlo
  MCResult res = run_guiding_center_mc(fl, cli);

  // 3) Final mu histogram
  {
    const std::string out_path = cli.out_prefix + "_mc_mu_hist.dat";
    if(!write_mu_hist_tecplot(out_path, res.mu_final, cli.nbins_mu, &err)){
      std::cerr << "Error: " << err << "\n";
      return 5;
    }
    std::cout << "Wrote: " << out_path << "\n";
  }

  // 4) mu moments vs s
  {
    const std::string out_path = cli.out_prefix + "_mc_mu_vs_s.dat";
    if(!write_mu_vs_s_tecplot(out_path, res.s_final, res.mu_final, fl.length(), cli.nbins_s, &err)){
      std::cerr << "Error: " << err << "\n";
      return 6;
    }
    std::cout << "Wrote: " << out_path << "\n";
  }

  // 5) 2D distribution f(s,mu) (normalized per s-bin)
  {
    const std::string out_path = cli.out_prefix + "_mc_mu_s_2d.dat";
    if(!write_mu_s_2d_tecplot(out_path, res.s_final, res.mu_final, fl.length(), cli.nbins_s, cli.nbins_mu, &err)){
      std::cerr << "Error: " << err << "\n";
      return 7;
    }
    std::cout << "Wrote: " << out_path << "\n";
  }

  std::cout << "Done.\n";
  return 0;
}

