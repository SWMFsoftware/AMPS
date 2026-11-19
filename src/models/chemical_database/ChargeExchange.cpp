#include "ChargeExchange.h"
#include "DataIO.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <utility>

namespace AMPS { namespace ChemicalDatabase { namespace ChargeExchange {
namespace {
bool Fail(std::string* error,const std::string& message) {
  if (error) *error=message;
  return false;
}

// Transcribed from the NUMERICAL table/fit on original ORNL-6086/V1 A-32
// (PDF page48), visually checked against the archived original PDF. These
// are not coefficients recovered from plotted pixels or a rounded web fit.
// Appendix I defines ln(sigma[cm2])=A0/2+sum(Aj*Tj(x),j=1..8), with x the
// logarithmic energy coordinate in [-1,1]. In particular A2 is NEGATIVE.
constexpr std::array<long double,9> kHeliumCoefficients{{
  -92.7819L,-4.69073L,-8.99666L,-0.490484L,0.935276L,
  -0.0504853L,-0.0410750L,-0.0730136L,0.0189320L}};
constexpr double kHeliumMinimumEvPerU=99.0,kHeliumMaximumEvPerU=11000000.0;

// Morgan et al. (1985), original printed p984/PDF page15, Table2, Na
// column: capture into ALL hydrogen states. Blanks below1keV/u are absent
// measurements, not zeros; Table3's 2s-only values are a different quantity.
constexpr std::array<double,11> kSodiumEnergyEvPerU{{
  1000,2000,3000,5000,7000,10000,20000,50000,70000,100000,200000}};
constexpr std::array<long double,11> kSodiumSourceUnits{{
  43,68,75,66,53,24,3.5L,0.71L,0.47L,0.32L,0.19L}};
constexpr long double kSodiumUnitToM2=1.0e-20L; // 10^-16 cm2 * 10^-4 m2/cm2.

bool Known(ChannelId id) {
  return id==ChannelId::ProtonHeOrnlA32 || id==ChannelId::ProtonNaMorgan1985LogLog;
}
} // namespace

ReferenceFiles FrozenV1References(const std::string& directory) {
  const std::string prefix=directory+(directory.empty() || directory.back()!='/' ? "/" : "");
  ReferenceFiles result;
  result.heliumPublicationPath=prefix+"ornl_redbook_v1_osti.pdf";
  result.sodiumPublicationPath=prefix+"nist_morgan_1985_retry2.pdf";
  result.heliumPublicationSha256="3bb3fdb4be8f0f806d378083be8798cf3ce347714c3d06d58f90ee483b10c030";
  result.sodiumPublicationSha256="7c412cb93f827b3412b6da000e674c63237cacb0e201733a6f5a4f431359952e";
  return result;
}

bool Registry::LoadPublishedReferences(const ReferenceFiles& source,std::string* error) {
  if (!error) return false;
  error->clear();const auto reviewed=FrozenV1References("");
  if (source.schema!=reviewed.schema || source.heliumPublicationSha256!=reviewed.heliumPublicationSha256 ||
      source.sodiumPublicationSha256!=reviewed.sodiumPublicationSha256)
    return Fail(error,"total-capture v1 schema/source identity differs from the reviewed original publications");
  if (!ReactionData::VerifyFileSha256(source.heliumPublicationPath,source.heliumPublicationSha256,error) ||
      !ReactionData::VerifyFileSha256(source.sodiumPublicationPath,source.sodiumPublicationSha256,error)) return false;
  // Only constants whose publication bytes/law have been checked are exposed.
  // The reviewed sources define total rates, not product-state distributions;
  // this load never installs a particle callback or a plasma distribution.
  ReferenceFiles candidate=source;references_=std::move(candidate);loaded_=true;return true;
}

bool Registry::GetMetadata(ChannelId id,ChannelMetadata* output,std::string* error) const {
  if (!output || !error) return false;
  error->clear();
  if (!loaded_ || !Known(id)) return Fail(error,"unloaded/unsupported total-capture channel; Ne/Ar/alpha/full kinetics are not supplied");
  ChannelMetadata result;
  if (id==ChannelId::ProtonHeOrnlA32) {
    result.id="H+_He_total_capture_ORNL_A32";
    result.reactants="H+ + He";result.products="H* + He+; unresolved final states";
    result.source="Barnett1990 ORNL-6086/V1 A-32, Appendix I";
    result.sourceSha256=references_.heliumPublicationSha256;
    result.law="original degree8 Chebyshev ln(sigma[cm2]); A0/2 convention";
    result.minimumEnergyEvPerU=kHeliumMinimumEvPerU;
    result.maximumEnergyEvPerU=kHeliumMaximumEvPerU;
  }
  else {
    result.id="H+_Na3s_total_capture_Morgan1985_Table2_log_log";
    result.reactants="H+ + Na(3s)";result.products="H* + Na+; all H states unresolved";
    result.source="MorganEtAl1985 DOI10.1063/1.555752, Table2 Na column";
    result.sourceSha256=references_.sodiumPublicationSha256;
    result.law="named source-knot log-log interpolation; exact populated knots; no extrapolation";
    result.minimumEnergyEvPerU=kSodiumEnergyEvPerU.front();
    result.maximumEnergyEvPerU=kSodiumEnergyEvPerU.back();
  }
  *output=std::move(result);return true;
}

bool Registry::CrossSectionM2(ChannelId id,double energy,double* output,std::string* error) const {
  if (!output || !error) return false;
  ChannelMetadata metadata;
  if (!GetMetadata(id,&metadata,error)) return false;
  if (!std::isfinite(energy) || energy<metadata.minimumEnergyEvPerU || energy>metadata.maximumEnergyEvPerU) {
    std::ostringstream message;message<<std::setprecision(17)
      <<"channel '"<<metadata.id<<"' target-rest energy "<<energy
      <<" eV/u is nonfinite/outside its source domain ["<<metadata.minimumEnergyEvPerU
      <<','<<metadata.maximumEnergyEvPerU<<"] eV/u; supply supported inputs or obtain a separately reviewed source; no clamping, zero substitution or extrapolation is allowed";
    return Fail(error,message.str());
  }
  long double result=0;
  if (id==ChannelId::ProtonHeOrnlA32) {
    const long double x=(2*std::log(static_cast<long double>(energy))-
        std::log(static_cast<long double>(kHeliumMinimumEvPerU))-
        std::log(static_cast<long double>(kHeliumMaximumEvPerU)))/
        (std::log(static_cast<long double>(kHeliumMaximumEvPerU))-
         std::log(static_cast<long double>(kHeliumMinimumEvPerU)));
    // Clenshaw recurrence for the original A0/2 series. The independent test
    // oracle constructs T_j directly at high precision, including endpoints;
    // it does not copy the webpage's incorrectly indexed sample program.
    long double next=0,next2=0;
    for (int j=8;j>=1;--j) {
      const long double current=kHeliumCoefficients[j]+2*x*next-next2;
      next2=next;next=current;
    }
    result=std::exp(kHeliumCoefficients[0]/2+x*next-next2)*1.0e-4L;
  }
  else {
    const auto upper=std::lower_bound(kSodiumEnergyEvPerU.begin(),kSodiumEnergyEvPerU.end(),energy);
    const std::size_t i=static_cast<std::size_t>(upper-kSodiumEnergyEvPerU.begin());
    if (upper!=kSodiumEnergyEvPerU.end() && energy==*upper)
      result=kSodiumSourceUnits[i]*kSodiumUnitToM2;
    else {
      // The interval is strictly inside source support and its two published
      // values are positive. This named interpolation creates a model value,
      // not a new measurement or a normalized/adjusted reference table.
      const long double fraction=std::log(static_cast<long double>(energy)/kSodiumEnergyEvPerU[i-1])/
          std::log(static_cast<long double>(kSodiumEnergyEvPerU[i])/kSodiumEnergyEvPerU[i-1]);
      result=std::exp(std::log(kSodiumSourceUnits[i-1])+
          fraction*std::log(kSodiumSourceUnits[i]/kSodiumSourceUnits[i-1]))*kSodiumUnitToM2;
    }
  }
  if (!std::isfinite(result) || result<=0 || result>std::numeric_limits<double>::max())
    return Fail(error,"source total-capture cross section is not a representable positive m2 value");
  *output=static_cast<double>(result);return true;
}

bool Registry::ColdBeamRateS(ChannelId id,double energy,double speed,double density,
    double* output,std::string* error) const {
  if (!output || !error) return false;
  if (!std::isfinite(speed) || speed<=0 || !std::isfinite(density) || density<0)
    return Fail(error,"cold-beam source-energy query requires explicit finite positive relative speed [m/s] and nonnegative density [m^-3]");
  double sigma;
  // Validate source identity/energy even at zero density. The application must
  // supply a consistent energy-speed pair in its declared frame; this algebra
  // does not approve a mass/relativistic mapping or an unbounded thermal PDF.
  if (!CrossSectionM2(id,energy,&sigma,error)) return false;
  const long double result=static_cast<long double>(density)*sigma*speed;
  if (!std::isfinite(result) || result>std::numeric_limits<double>::max())
    return Fail(error,"cold-beam total-capture rate [s^-1] is not representable");
  *output=static_cast<double>(result);return true;
}

bool Registry::SodiumSourceKnots(std::vector<double>* energies,std::vector<double>* values,std::string* error) const {
  if (!energies || !values || !error) return false;
  error->clear();if (!loaded_) return Fail(error,"total-capture registry is not loaded");
  std::vector<double> energy(kSodiumEnergyEvPerU.begin(),kSodiumEnergyEvPerU.end()),sigma;
  for (auto value:kSodiumSourceUnits) sigma.push_back(static_cast<double>(value*kSodiumUnitToM2));
  *energies=std::move(energy);*values=std::move(sigma);return true;
}

} } } // namespace AMPS::ChemicalDatabase::ChargeExchange
