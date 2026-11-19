#include "Photochemistry.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <set>
#include <sstream>
#include <utility>

namespace AMPS { namespace ChemicalDatabase {
namespace {
bool Fail(std::string* error,const std::string& message) {
  if (error) *error=message;
  return false;
}

using Row=ReactionData::Row;
using Csv=ReactionData::Csv;
using ReactionData::ReadPinnedCsv;

bool Number(const std::string& field,double* value) {
  std::istringstream input(field);input.imbue(std::locale::classic());
  double candidate;std::string extra;
  if (!(input>>candidate) || (input>>extra) || !std::isfinite(candidate)) return false;
  *value=candidate;return true;
}

struct Definition {const char *id,*key,*parent,*type,*products,*label;};
const Definition kDefinitions[]={
 {"PHIDRATES:He:00:he_e","He","He","photoionization","He+ + e","He+"},
 {"PHIDRATES:Na:00:na_e","Na","Na","photoionization","Na+ + e","Na+"},
 {"PHIDRATES:Ne:00:ne_e","Ne","Ne","photoionization","Ne+ + e","Ne+"},
 {"PHIDRATES:Ar:00:ar_e","Ar","Ar","photoionization","Ar+ + e","Ar+"},
 {"PHIDRATES:H2O:00:total","H2O","H2O","total_photoabsorption","TOTAL","Total"},
 {"PHIDRATES:H2O:01:h_oh","H2O","H2O","photodissociation","H + OH","H/OH"},
 {"PHIDRATES:H2O:02:h2_o_1d","H2O","H2O","photodissociation","H2 + O(1D)","H2/O1D"},
 {"PHIDRATES:H2O:03:o_h_h","H2O","H2O","photodissociation","O + H + H","O/H/H"},
 {"PHIDRATES:H2O:04:oh_h_e","H2O","H2O","photodissociative_ionization","OH+ + H + e","OH+H"},
 {"PHIDRATES:H2O:05:o_h2_e","H2O","H2O","photodissociative_ionization","O+ + H2 + e","O+H2"},
 {"PHIDRATES:H2O:06:h_oh_e","H2O","H2O","photodissociative_ionization","H+ + OH + e","H+OH"},
 {"PHIDRATES:H2O:07:h2o_e","H2O","H2O","photoionization","H2O+ + e","H2O+"},
 {"PHIDRATES:OH_exp:00:total","OH_exp","OH","total_photoabsorption","TOTAL","Total"},
 {"PHIDRATES:OH_exp:01:o_3p_h","OH_exp","OH","photodissociation","O(3P) + H","O/H"},
 {"PHIDRATES:OH_exp:02:o_1s_h","OH_exp","OH","photodissociation","O(1S) + H","O1S/H"},
 {"PHIDRATES:OH_exp:03:oh_e","OH_exp","OH","photoionization","OH+ + e","OH+"},
 {"PHIDRATES:OH_exp:04:o_1d_h","OH_exp","OH","photodissociation","O(1D) + H","O1D/H"}};

std::map<std::string,unsigned> Nuclei(const std::string& name) {
  if (name=="H") return {{"H",1}};
  if (name=="H2") return {{"H",2}};
  if (name=="O") return {{"O",1}};
  if (name=="OH") return {{"H",1},{"O",1}};
  if (name=="H2O") return {{"H",2},{"O",1}};
  if (name=="He" || name=="Na" || name=="Ne" || name=="Ar") return {{name,1}};
  return {};
}

bool Products(Photochemistry::Channel* channel,std::string* error) {
  if (!channel->reactive) return true;
  std::size_t begin=0;int charge=0;std::map<std::string,unsigned> total;
  while (begin<channel->sourceProducts.size()) {
    const auto end=channel->sourceProducts.find(" + ",begin);
    std::string token=channel->sourceProducts.substr(begin,end==std::string::npos ? end : end-begin);
    Photochemistry::Product product;
    if (token=="e") {product.species="e";product.charge=-1;}
    else {
      if (!token.empty() && token.back()=='+') {product.charge=1;token.pop_back();}
      const auto state=token.find('(');
      if (state!=std::string::npos) {
        if (token.back()!=')' || token.substr(0,state)!="O") return Fail(error,"unknown D04 electronic-state notation");
        product.electronicState=token.substr(state+1,token.size()-state-2);token.erase(state);
      }
      product.species=token+(product.charge==1 ? "+" : "");product.nuclei=Nuclei(token);
      if (product.nuclei.empty()) return Fail(error,"unknown D04 product nuclei");
    }
    charge+=product.charge;
    for (const auto& item:product.nuclei) total[item.first]+=item.second;
    auto same=std::find_if(channel->products.begin(),channel->products.end(),[&](const Photochemistry::Product& p) {
      return p.species==product.species && p.electronicState==product.electronicState && p.charge==product.charge;});
    if (same==channel->products.end()) channel->products.push_back(product);
    else ++same->multiplicity;
    if (end==std::string::npos) break;
    begin=end+3;
  }
  if (charge!=0 || total!=Nuclei(channel->physicalParent)) return Fail(error,"D04 products violate neutral-parent charge/nucleus identity");
  return true;
}

bool SolarIndex(Photochemistry::SolarCase solar,std::size_t* index,std::string* error) {
  if (solar==Photochemistry::SolarCase::Quiet) {*index=0;return true;}
  if (solar==Photochemistry::SolarCase::Active) {*index=1;return true;}
  return Fail(error,"D04 solar case must be the explicit quiet or active endpoint");
}
} // namespace

namespace Photochemistry {
SourceFiles FrozenV1Sources(const std::string& directory) {
  const std::string prefix=directory+(directory.empty() || directory.back()!='/' ? "/" : "");
  SourceFiles source;
  source.channelsPath=prefix+"reaction_channels.csv";
  source.ratesPath=prefix+"photo_rates_1au.csv";
  source.thresholdsPath=prefix+"thresholds.csv";
  source.crossSectionsPath=prefix+"cross_sections.csv";
  source.channelsSha256="2d9fa1c8eb72d3bbba0b91ab4c9a912e6ebce2abc80b251096df2f41a4354061";
  source.ratesSha256="3aa54915a08258371dbe975949e6afb2ab74618164fadca8bceebe876c926847";
  source.thresholdsSha256="a3335b3b37c3d275a4f72e92ae8e5b1f6ecfcbc4e2d84735cb20694dc93b48f1";
  source.crossSectionsSha256="e900eb5514b616dcbbde7e3ca124b5c06d4d22a28a3cd429ed58a14c1933c860";
  return source;
}

bool Registry::LoadFrozenD04(const SourceFiles& source,std::string* error) {
  if (!error) return false;
  error->clear();const SourceFiles reviewed=FrozenV1Sources("");
  if (source.schema!=reviewed.schema || source.channelsSha256!=reviewed.channelsSha256 ||
      source.ratesSha256!=reviewed.ratesSha256 || source.thresholdsSha256!=reviewed.thresholdsSha256 ||
      source.crossSectionsSha256!=reviewed.crossSectionsSha256)
    return Fail(error,"D04 v1 schema/hash identity differs from the reviewed source version");
  Csv catalogue,rates,thresholds,cross;
  if (!ReadPinnedCsv(source.channelsPath,source.channelsSha256,6012,
      {"species","reaction_id","reaction_type","parent","products","phidrates_channel_label","source_reference","source_file"},&catalogue,error) ||
      !ReadPinnedCsv(source.ratesPath,source.ratesSha256,11460,
      {"species","reaction_id","reaction_type","products","rate_1au_s-1","excess_energy_eV","solar_activity","optical_depth","distance_au","source_reference","source_file","source_record_or_index"},&rates,error) ||
      !ReadPinnedCsv(source.thresholdsPath,source.thresholdsSha256,2948,
      {"species","reaction_id","products","threshold_wavelength_angstrom","threshold_wavelength_nm","threshold_eV","source_reference","source_file","source_record_or_index"},&thresholds,error) ||
      !ReadPinnedCsv(source.crossSectionsPath,source.crossSectionsSha256,1034075,
      {"species","reaction_id","reaction_type","parent","products","wavelength_nm","photon_energy_eV","cross_section_m2","threshold_eV","branch_fraction","source_reference","source_file","source_record_or_index"},&cross,error)) return false;
  if (catalogue.size()!=17 || rates.size()!=30 || thresholds.size()!=4 || cross.size()!=2305)
    return Fail(error,"D04 v1 complete native record counts are not17/30/4/2305");
  Registry candidate;std::map<std::string,const Definition*> definitions;
  for (const auto& d:kDefinitions) definitions[d.id]=&d;
  for (const auto& row:catalogue) {
    const auto found=definitions.find(row[1]);
    if (found==definitions.end()) return Fail(error,"unknown D04 reaction identity");
    const auto& d=*found->second;
    if (row[0]!=d.key || row[2]!=d.type || row[3]!=d.parent || row[4]!=d.products || row[5]!=d.label)
      return Fail(error,"D04 reaction metadata/state/multiplicity differs from reviewed definition");
    Channel channel;channel.speciesKey=row[0];channel.id=row[1];channel.reactionType=row[2];
    channel.physicalParent=row[3];channel.sourceProducts=row[4];channel.sourceChannelLabel=row[5];
    channel.sourceReference=row[6];channel.sourceFile=row[7];
    channel.reactive=channel.reactionType!="total_photoabsorption";
    channel.experimentalOh=channel.speciesKey=="OH_exp";
    if (channel.sourceFile!="raw/phidrates/"+channel.speciesKey+"/cross_sections/data/brnout_001.raw" ||
        !Products(&channel,error)) return Fail(error,"D04 source path or product identity mismatch");
    if (!candidate.byId_.emplace(channel.id,candidate.channels_.size()).second)
      return Fail(error,"duplicate D04 reaction identity");
    candidate.channels_.push_back(std::move(channel));
  }
  std::map<std::string,std::array<bool,2>> seenRates;
  for (const auto& row:rates) {
    const auto found=candidate.byId_.find(row[1]);
    if (found==candidate.byId_.end()) return Fail(error,"rate references unknown D04 channel");
    auto& c=candidate.channels_[found->second];double rate,energy,activity,tau,distance;
    if (!c.reactive || row[0]!=c.speciesKey || row[2]!=c.reactionType || row[3]!=c.sourceProducts || row[9]!=c.sourceReference ||
        !Number(row[4],&rate) || rate<=0 || !Number(row[5],&energy) || energy<0 ||
        !Number(row[6],&activity) || (activity!=0 && activity!=1) ||
        !Number(row[7],&tau) || tau!=0 || !Number(row[8],&distance) || distance!=1 || row[11].empty())
      return Fail(error,"invalid D04 rate, endpoint, units/domain or source metadata");
    const std::size_t index=activity==0 ? 0 : 1;
    if (seenRates[c.id][index]) return Fail(error,"duplicate D04 rate endpoint");
    const std::string solar=index==0 ? "0p000" : "1p000";
    if (row[10]!="raw/phidrates/"+c.speciesKey+"/solar_rates/activity_"+solar+"_tau_0p000/result.html")
      return Fail(error,"D04 rate source does not match its declared solar endpoint");
    seenRates[c.id][index]=true;c.rate1AuS[index]=rate;c.meanExcessEnergyEv[index]=energy;
    c.rateSourceFile[index]=row[10];c.rateSourceRecord[index]=row[11];
  }
  std::set<std::string> seenThresholds;
  for (const auto& row:thresholds) {
    const auto found=candidate.byId_.find(row[1]);double angstrom,nm,energy;
    if (found==candidate.byId_.end()) return Fail(error,"unknown D04 threshold identity");
    auto& c=candidate.channels_[found->second];
    if (!c.reactive || c.speciesKey=="H2O" || c.speciesKey=="OH_exp" ||
        row[0]!=c.speciesKey || row[2]!=c.sourceProducts || row[6]!=c.sourceReference || row[7]!=c.sourceFile || row[8].empty() ||
        !Number(row[3],&angstrom) || angstrom<=0 || !Number(row[4],&nm) || nm<=0 ||
        !Number(row[5],&energy) || energy<=0 || !seenThresholds.insert(c.id).second)
      return Fail(error,"invalid/duplicate atomic D04 threshold metadata");
    c.hasThreshold=true;c.thresholdEnergyEv=energy;c.thresholdWavelengthNm=nm;
  }
  for (const auto& row:cross) {
    const auto found=candidate.byId_.find(row[1]);double nm,energy,sigma;
    if (found==candidate.byId_.end()) return Fail(error,"unknown D04 cross-section identity");
    auto& c=candidate.channels_[found->second];
    if (row[0]!=c.speciesKey || row[2]!=c.reactionType || row[3]!=c.physicalParent || row[4]!=c.sourceProducts ||
        row[10]!=c.sourceReference || row[11]!=c.sourceFile || row[12].empty() ||
        !Number(row[5],&nm) || nm<=0 || !Number(row[6],&energy) || energy<=0 ||
        !Number(row[7],&sigma) || sigma<0)
      return Fail(error,"invalid D04 cross-section row/units/source");
    if (!row[8].empty()) {double threshold;if (!Number(row[8],&threshold) || !c.hasThreshold || threshold!=c.thresholdEnergyEv)
      return Fail(error,"D04 cross-section threshold metadata mismatch");}
    else if (c.hasThreshold) return Fail(error,"atomic D04 threshold unexpectedly missing");
    if (!row[9].empty()) {double fraction;if (!Number(row[9],&fraction) || fraction<0 || fraction>1)
      return Fail(error,"D04 source branch fraction is outside [0,1]");}
    auto& domain=c.tabulated;
    if (domain.sampleCount==0) {
      domain.minimumWavelengthNm=domain.maximumWavelengthNm=nm;
      domain.minimumPhotonEnergyEv=domain.maximumPhotonEnergyEv=energy;
    }
    else {
      domain.minimumWavelengthNm=std::min(domain.minimumWavelengthNm,nm);
      domain.maximumWavelengthNm=std::max(domain.maximumWavelengthNm,nm);
      domain.minimumPhotonEnergyEv=std::min(domain.minimumPhotonEnergyEv,energy);
      domain.maximumPhotonEnergyEv=std::max(domain.maximumPhotonEnergyEv,energy);
    }
    ++domain.sampleCount;
  }
  for (const auto& c:candidate.channels_)
    if (!c.tabulated.sampleCount || (c.reactive && (!seenRates[c.id][0] || !seenRates[c.id][1])))
      return Fail(error,"incomplete D04 channel/endpoints/sample inventory");
  candidate.loaded_=true;candidate.sources_=source;*this=std::move(candidate);return true;
}

bool Registry::GetChannel(const std::string& id,Channel* output,std::string* error) const {
  if (!output || !error) return false;
  error->clear();const auto found=byId_.find(id);
  if (!loaded_ || found==byId_.end()) return Fail(error,"unloaded or unknown D04 channel '"+id+"'");
  *output=channels_[found->second];return true;
}
bool Registry::ChannelsForParent(const std::string& key,std::vector<Channel>* output,std::string* error) const {
  if (!output || !error) return false;
  error->clear();std::vector<Channel> candidate;
  if (loaded_) for (const auto& c:channels_) if (c.reactive && c.speciesKey==key) candidate.push_back(c);
  if (candidate.empty()) return Fail(error,"unloaded/unsupported D04 species key '"+key+"'; experimental OH must be explicit OH_exp");
  *output=std::move(candidate);return true;
}
bool Registry::BranchFractions(const std::string& key,SolarCase solar,
    std::vector<std::pair<std::string,double>>* output,std::string* error) const {
  if (!output || !error) return false;
  std::size_t index;std::vector<Channel> channels;
  if (!SolarIndex(solar,&index,error) || !ChannelsForParent(key,&channels,error)) return false;
  long double total=0;for (const auto& c:channels) total+=c.rate1AuS[index];
  if (!(total>0)) return Fail(error,"D04 complete family has no reactive hazard");
  std::vector<std::pair<std::string,double>> result;
  for (const auto& c:channels) result.emplace_back(c.id,static_cast<double>(c.rate1AuS[index]/total));
  *output=std::move(result);return true;
}
bool Registry::SelectChannel(const std::string& key,SolarCase solar,double u,
    std::string* output,std::string* error) const {
  if (!output || !error) return false;
  error->clear();
  if (!std::isfinite(u) || u<0 || u>=1) return Fail(error,"conditional channel draw must be finite in [0,1)");
  std::size_t index;
  if (!SolarIndex(solar,&index,error)) return false;
  // Particle-event queries must not copy every source reference/product list.
  // Walk the immutable complete registry directly, preserving its frozen
  // source order and every positive hazard in the conditional distribution.
  long double total=0;
  if (loaded_) for (const auto& c:channels_)
    if (c.reactive && c.speciesKey==key) total+=c.rate1AuS[index];
  if (!(total>0)) return Fail(error,"unloaded/unsupported D04 reactive family '"+key+"'; experimental OH requires the explicit OH_exp key");
  const long double target=u*total;long double cumulative=0;
  for (const auto& c:channels_) if (c.reactive && c.speciesKey==key) {
    cumulative+=c.rate1AuS[index];
    if (target<cumulative) {*output=c.id;return true;}
  }
  return Fail(error,"conditional D04 selector did not span its complete hazard family");
}
bool Registry::ChannelRate1AuS(const std::string& id,SolarCase solar,double* output,std::string* error) const {
  if (!output || !error) return false;
  error->clear();
  std::size_t index;
  if (!SolarIndex(solar,&index,error)) return false;
  const auto found=byId_.find(id);
  if (!loaded_ || found==byId_.end()) return Fail(error,"unloaded/unknown D04 channel");
  const auto& c=channels_[found->second];
  if (!c.reactive) return Fail(error,"total_photoabsorption is not a reactive destruction channel");
  *output=c.rate1AuS[index];return true;
}
bool Registry::TotalDestructionRateS(const std::string& key,SolarCase solar,double distance,
    double visible,double* output,std::string* error) const {
  if (!output || !error) return false;
  error->clear();
  if (!std::isfinite(distance) || distance<=0 || !std::isfinite(visible) || visible<0 || visible>1)
    return Fail(error,"photon rate requires finite positive distance [AU] and normalized gas visibility in [0,1]");
  std::size_t index;
  if (!SolarIndex(solar,&index,error)) return false;
  long double rate=0;
  if (loaded_) for (const auto& c:channels_)
    if (c.reactive && c.speciesKey==key) rate+=c.rate1AuS[index];
  if (!(rate>0)) return Fail(error,"unloaded/unsupported D04 reactive family '"+key+"'; experimental OH requires the explicit OH_exp key");
  rate=rate*visible/static_cast<long double>(distance)/distance;
  if (!std::isfinite(rate) || rate<0 || rate>std::numeric_limits<double>::max())
    return Fail(error,"scaled photon rate is not a representable finite nonnegative s^-1 value");
  *output=static_cast<double>(rate);return true;
}

} } } // namespace AMPS::ChemicalDatabase::Photochemistry
