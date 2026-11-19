#include "ElectronImpact.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <map>
#include <sstream>
#include <utility>

namespace AMPS { namespace ChemicalDatabase { namespace ElectronImpact {
namespace {
constexpr double kBoltzmannJPerK=1.380649e-23;
constexpr double kElectronVoltJ=1.602176634e-19;
constexpr long double kCm3ToM3=1.0e-6L;
const char* const kElements[]={"", "H","He","Li","Be","B","C","N","O",
    "F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
    "Sc","Ti","V","Cr","Mn","Fe","Co","Ni"};

std::string Trim(const std::string& value) {
  const auto first=value.find_first_not_of(" \t\r\n");
  return first==std::string::npos ? std::string() :
      value.substr(first,value.find_last_not_of(" \t\r\n")-first+1);
}

template<class T> bool Number(const std::string& field,T* output) {
  // Classic locale makes native decimal points independent of the machine's
  // locale. Require complete consumption; NaN/Inf are rejected separately.
  std::istringstream input(Trim(field));input.imbue(std::locale::classic());
  T candidate;std::string extra;
  if (!(input>>candidate) || (input>>extra)) return false;
  *output=candidate;return true;
}

bool Index(NeutralSpecies species,std::size_t* index) {
  switch (species) {
    case NeutralSpecies::He:*index=0;return true;
    case NeutralSpecies::Na:*index=1;return true;
    case NeutralSpecies::Ne:*index=2;return true;
    case NeutralSpecies::Ar:*index=3;return true;
  }
  return false;
}

bool Positive(double value) { return std::isfinite(value) && value>0.0; }
} // namespace

const char* NeutralSpeciesName(NeutralSpecies species) {
  switch (species) {
    case NeutralSpecies::He:return "He";
    case NeutralSpecies::Na:return "Na";
    case NeutralSpecies::Ne:return "Ne";
    case NeutralSpecies::Ar:return "Ar";
  }
  return "unsupported";
}

bool NeutralSpeciesFromName(const std::string& name,NeutralSpecies* species,
                            std::string* error) {
  if (species==nullptr || error==nullptr) return false;
  error->clear();
  for (const auto value:{NeutralSpecies::He,NeutralSpecies::Na,
                         NeutralSpecies::Ne,NeutralSpecies::Ar}) {
    if (name==NeutralSpeciesName(value)) { *species=value;return true; }
  }
  *error="unsupported neutral electron-impact species '"+name+"' (He/Na/Ne/Ar only)";
  return false;
}

bool ElectronTemperatureEvFromK(double temperatureK,double* temperatureEv,
                                std::string* error) {
  if (temperatureEv==nullptr || error==nullptr) return false;
  error->clear();
  const double candidate=temperatureK*(kBoltzmannJPerK/kElectronVoltJ);
  if (!Positive(temperatureK) || !Positive(candidate)) {
    *error="electron temperature conversion requires finite positive representable K/eV";
    return false;
  }
  *temperatureEv=candidate;return true;
}

bool MaxwellianVoronovTable::LoadNativeTable(const std::string& path,
                                            std::string* error) {
  if (error==nullptr) return false;
  error->clear();std::ifstream input(path);
  if (!input) { *error="cannot open native Voronov cfit.dat '"+path+"'";return false; }
  std::map<std::pair<int,int>,Coefficient> rows;
  std::string line;std::size_t lineNumber=0;
  while (std::getline(input,line)) {
    ++lineNumber;
    if (!line.empty() && line.back()=='\r') line.pop_back();
    const auto fail=[&](const std::string& reason) {
      *error="native Voronov record "+std::to_string(lineNumber)+": "+reason;
      return false;
    };
    if (line.size()!=53) return fail("expected exact53-byte cfit.txt layout");
    for (std::size_t separator:{5,8,11,19,21,31,38,43,48})
      if (line[separator]!=' ') return fail("invalid fixed-field separator");
    Coefficient row;row.sourceLine=lineNumber;int minimumEv=0,maxKev=0;
    if (!Number(line.substr(6,2),&row.atomicNumber) ||
        !Number(line.substr(9,2),&row.electronCount) ||
        !Number(line.substr(12,7),&row.ionizationEnergyEv) ||
        !Number(line.substr(20,1),&row.p) || !Number(line.substr(22,9),&row.aCm3S) ||
        !Number(line.substr(32,6),&row.x) || !Number(line.substr(39,4),&row.k) ||
        !Number(line.substr(44,4),&minimumEv) ||
        !Number(line.substr(49,4),&maxKev)) return fail("malformed numeric field");
    row.temperatureMinEv=minimumEv;
    if (row.atomicNumber<1 || row.atomicNumber>28 || row.electronCount<1 ||
        row.electronCount>row.atomicNumber || (row.p!=0 && row.p!=1) ||
        !Positive(row.ionizationEnergyEv) || !Positive(row.aCm3S) ||
        !std::isfinite(row.x) || row.x<0 || !std::isfinite(row.k) || row.k<0 ||
        !Positive(row.temperatureMinEv) || !Positive(maxKev))
      return fail("invalid atomic stage, coefficient or native temperature bound");
    row.temperatureMaxEv=1000.0*maxKev;
    if (!Positive(row.temperatureMaxEv) || row.temperatureMaxEv<row.temperatureMinEv)
      return fail("native Tmax[keV] is below Tmin[eV]");
    std::string label;
    for (char character:line.substr(0,5)) if (character!=' ') label+=character;
    const int charge=row.atomicNumber-row.electronCount;
    const std::string expected=std::string(kElements[row.atomicNumber])+
        (charge==0 ? "" : std::to_string(charge)+"+");
    if (label!=expected) return fail("ion label disagrees with atomic number/electron count");
    if (!rows.emplace(std::make_pair(row.atomicNumber,row.electronCount),row).second)
      return fail("duplicate atomic-number/electron-count row");
  }
  if (!input.eof()) { *error="I/O failure reading native Voronov table";return false; }
  if (lineNumber!=406 || rows.size()!=406) {
    *error="native Voronov table requires all406 unique atomic-stage records";return false;
  }
  for (int z=1;z<=28;++z) for (int n=1;n<=z;++n)
    if (!rows.count({z,n})) { *error="native Voronov table is missing an atomic stage";return false; }
  std::array<Coefficient,4> selected;
  const int z[]={2,11,10,18};
  for (std::size_t index=0;index<selected.size();++index) selected[index]=rows.at({z[index],z[index]});
  // Transactional commit: failed reload leaves an earlier qualified table
  // intact; no partial species selection is ever published to callbacks.
  neutral_=selected;loaded_=true;return true;
}

bool MaxwellianVoronovTable::GetCoefficient(NeutralSpecies species,
    Coefficient* coefficient,std::string* error) const {
  if (coefficient==nullptr || error==nullptr) return false;
  error->clear();std::size_t index=0;
  if (!loaded_ || !Index(species,&index)) {
    *error="Voronov neutral table is uninitialized or species is unsupported";return false;
  }
  *coefficient=neutral_[index];return true;
}

bool MaxwellianVoronovTable::RateCoefficientM3S(NeutralSpecies species,
    double temperatureEv,double* coefficientM3S,std::string* error) const {
  if (coefficientM3S==nullptr || error==nullptr) return false;
  Coefficient row;
  if (!GetCoefficient(species,&row,error)) return false;
  if (!Positive(temperatureEv) || temperatureEv<row.temperatureMinEv ||
      temperatureEv>row.temperatureMaxEv) {
    std::ostringstream message;message << std::setprecision(17)
        << "electron temperature [eV] must be finite and in native inclusive ["
        << row.temperatureMinEv << ',' << row.temperatureMaxEv << "] for "
        << NeutralSpeciesName(species) << "; no Voronov extrapolation is allowed";
    *error=message.str();return false;
  }
  const long double u=static_cast<long double>(row.ionizationEnergyEv)/temperatureEv;
  // The published expression is evaluated once, in log space. In particular
  // the legacy single-precision cfit.f U>80 cutoff is not a physical boundary:
  // only actual representability underflow may return zero inside valid data.
  const long double logarithm=std::log(static_cast<long double>(row.aCm3S))+
      std::log1p(row.p*std::sqrt(u))-std::log(row.x+u)+row.k*std::log(u)-u+
      std::log(kCm3ToM3);
  const long double value=std::exp(logarithm);
  if (!std::isfinite(value) || value<0 || value>std::numeric_limits<double>::max()) {
    *error="Voronov coefficient is not representable as finite nonnegative m^3/s";return false;
  }
  *coefficientM3S=static_cast<double>(value);return true;
}

bool MaxwellianVoronovTable::IonizationRateS(NeutralSpecies species,double temperatureEv,
    double densityM3,double* rateS,std::string* error) const {
  if (rateS==nullptr || error==nullptr) return false;
  error->clear();
  if (!std::isfinite(densityM3) || densityM3<0) {
    *error="electron density [m^-3] must be finite and nonnegative";return false;
  }
  double coefficient=0;
  // Validate the selected species and temperature even at zero density; zero
  // plasma does not turn an invalid/unsupported physical configuration valid.
  if (!RateCoefficientM3S(species,temperatureEv,&coefficient,error)) return false;
  const long double rate=static_cast<long double>(densityM3)*coefficient;
  if (!std::isfinite(rate) || rate>std::numeric_limits<double>::max()) {
    *error="electron-impact rate [s^-1] overflows; review density/temperature inputs";return false;
  }
  *rateS=static_cast<double>(rate);return true;
}

} } } // namespace AMPS::ChemicalDatabase::ElectronImpact
