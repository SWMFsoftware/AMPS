// Stand-alone tests of srcHartley2/CometChemistryVariants.{h,cpp}: SHA-256,
// the OH lifetime table reader, linear interpolation with clamping, and the
// SA1982_VH rate formula (srcHartley2/CHEMISTRY_VARIANTS.md, sections 3-4).
//
// Usage:
//   test_chemistry_variants --self-test TABLE
//       built-in checks; prints PASS/FAIL lines and exits non-zero on failure
//   test_chemistry_variants --rates TABLE SCALE_S ION_FRACTION V1 [V2 ...]
//       prints "v ohtaus clamped k_d k_ion" with 17 significant digits for each
//       velocity; test_chemistry_variants.py compares these with its own,
//       independent evaluation of the legacy formula.
// No AMPS, PIC, or MPI code is involved.

#include "CometChemistryVariants.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

int Failures=0;

void Check(bool condition,const std::string &name) {
  std::printf("%s %s\n",condition ? "PASS" : "FAIL",name.c_str());
  if (!condition) Failures++;
}

bool Near(double a,double b,double relative) {
  return std::fabs(a-b)<=relative*std::fmax(std::fabs(a),std::fabs(b));
}

void ExpectParseError(const std::string &text,const std::string &fragment,
    const std::string &name) {
  Comet::ChemistryVariants::cLifetimeTable table;
  std::string error;
  const bool ok=Comet::ChemistryVariants::ParseLifetimeTable(text,"t.csv",&table,&error);
  Check(!ok && error.find(fragment)!=std::string::npos,
      name+" -> '"+error+"'");
}

int SelfTest(const std::string &path) {
  using namespace Comet::ChemistryVariants;
  // FIPS 180-4 example vectors.
  Check(Sha256Hex("")==
      "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
      "sha256 of the empty string");
  Check(Sha256Hex("abc")==
      "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
      "sha256 of 'abc'");
  Check(Sha256Hex("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq")==
      "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1",
      "sha256 of the 448-bit message");

  cLifetimeTable table;
  std::string error;
  Check(ReadLifetimeTable(path,&table,&error),"read "+path+" "+error);
  Check(table.VelocityKmS.size()==61,"61 rows");
  Check(table.TauFuvS==4.24e5,"tau_fuv_s = 4.24e5 s");
  Check(table.VelocityKmS.front()==-60.0 && table.VelocityKmS.back()==60.0,
      "grid -60..+60 km/s");
  Check(table.Sha256.size()==64,"sha256 recorded");

  // Grid points return the table value exactly; mid-points the mean.
  Check(InterpolateOhtaus(table,-40.0).Value==3.80 &&
      !InterpolateOhtaus(table,-40.0).Clamped,"grid point -40 km/s = 3.80");
  Check(Near(InterpolateOhtaus(table,-39.0).Value,0.5*(3.80+3.50),1e-15),
      "mid-point -39 km/s");
  Check(Near(InterpolateOhtaus(table,-39.7).Value,3.755,1e-14),
      "case 1, -39.7 km/s = 3.755");
  // Clamping at both ends, and exact ends are not clamps.
  Check(InterpolateOhtaus(table,-75.0).Value==2.85 &&
      InterpolateOhtaus(table,-75.0).Clamped,"clamp below -60 km/s");
  Check(InterpolateOhtaus(table,61.0).Value==2.85 &&
      InterpolateOhtaus(table,61.0).Clamped,"clamp above +60 km/s");
  Check(!InterpolateOhtaus(table,-60.0).Clamped &&
      !InterpolateOhtaus(table,60.0).Clamped,"table ends are not clamps");

  cSA1982Rates rates;
  Check(EvaluateSA1982(table,0.4,0.833e5,0.05,&rates,&error),"evaluate +0.4 km/s");
  Check(Near(rates.Dissociation1AU,1.0/4.24e5+1.0/(0.833e5*2.76),1e-14),
      "k_d at +0.4 km/s");
  Check(Near(rates.Ionization1AU,0.05*rates.Dissociation1AU,1e-15),
      "k_ion = 0.05 k_d");
  Check(!EvaluateSA1982(table,NAN,0.833e5,0.05,&rates,&error),"reject NaN velocity");
  Check(!EvaluateSA1982(table,0.0,0.0,0.05,&rates,&error),"reject zero scale");
  Check(!EvaluateSA1982(table,0.0,0.833e5,-0.1,&rates,&error),"reject negative fraction");

  // Malformed tables: each must fail with its own message.
  const std::string head="tau_fuv_s,4.24e5\nv_h_km_s,ohtaus_1e5_s\n";
  ExpectParseError(head+"0,1\n0,2\n","increase strictly","repeated velocity");
  ExpectParseError(head+"1,1\n0,2\n","increase strictly","decreasing velocity");
  ExpectParseError(head+"0,1\n2,nan\n","finite positive","NaN value");
  ExpectParseError(head+"0,1\n2,-1\n","finite positive","negative value");
  ExpectParseError(head+"0,1\n","at least two","one row");
  ExpectParseError("v_h_km_s,ohtaus_1e5_s\n0,1\n2,1\n","tau_fuv_s","missing tau_fuv_s");
  ExpectParseError("tau_fuv_s,0\nv_h_km_s,ohtaus_1e5_s\n0,1\n2,1\n",
      "finite positive number of seconds","zero tau_fuv_s");
  ExpectParseError("tau_fuv_s,1\nv,o\n0,1\n2,1\n","column header","wrong header");
  ExpectParseError(head+"0,1,3\n2,1\n","two comma-separated","three fields");
  ExpectParseError(head+"0,1\n2,1x\n","finite positive","trailing characters");
  ExpectParseError(head+"0,1\nfoo,1\n","t.csv:4","error names file and line");
  {
    cLifetimeTable missing;
    Check(!ReadLifetimeTable("/nonexistent/oh.csv",&missing,&error) &&
        error.find("cannot open")!=std::string::npos,"missing file");
  }
  // A non-uniform grid is accepted and interpolated on its own spacing.
  {
    cLifetimeTable custom;
    Check(ParseLifetimeTable("# c\n\ntau_fuv_s,1e5\nv_h_km_s,ohtaus_1e5_s\n-1,1\n0,2\n10,4\n",
        "c.csv",&custom,&error),"non-uniform grid accepted");
    Check(Near(InterpolateOhtaus(custom,5.0).Value,3.0,1e-15),"non-uniform interpolation");
  }

  std::printf("%s test_chemistry_variants failures=%d\n",
      Failures==0 ? "PASS" : "FAIL",Failures);
  return Failures==0 ? EXIT_SUCCESS : EXIT_FAILURE;
}

int Rates(int argc,char **argv) {
  using namespace Comet::ChemistryVariants;
  cLifetimeTable table;
  std::string error;
  if (!ReadLifetimeTable(argv[2],&table,&error)) {
    std::fprintf(stderr,"%s\n",error.c_str());
    return EXIT_FAILURE;
  }
  const double scale=std::atof(argv[3]);
  const double fraction=std::atof(argv[4]);
  for (int argument=5;argument<argc;argument++) {
    cSA1982Rates rates;
    if (!EvaluateSA1982(table,std::atof(argv[argument]),scale,fraction,&rates,&error)) {
      std::fprintf(stderr,"%s\n",error.c_str());
      return EXIT_FAILURE;
    }
    std::printf("%.17g %.17g %d %.17g %.17g\n",rates.VelocityKmS,rates.Ohtaus,
        rates.Clamped ? 1 : 0,rates.Dissociation1AU,rates.Ionization1AU);
  }
  std::printf("sha256 %s\n",table.Sha256.c_str());
  return EXIT_SUCCESS;
}

}  // namespace

int main(int argc,char **argv) {
  if (argc==3 && std::string(argv[1])=="--self-test") return SelfTest(argv[2]);
  if (argc>=6 && std::string(argv[1])=="--rates") return Rates(argc,argv);
  std::fprintf(stderr,"usage: %s --self-test TABLE | --rates TABLE SCALE_S ION_FRACTION V...\n",
      argv[0]);
  return EXIT_FAILURE;
}
