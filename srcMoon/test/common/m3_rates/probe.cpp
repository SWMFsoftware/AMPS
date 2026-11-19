// Thin PIC-independent adapter. All reaction data, rates, identities and
// selectors below call the production libraries; expected physics is evaluated
// independently by test.py/primary_reference.py, never duplicated in this file.
#include "Photochemistry.h"
#include "ChargeExchange.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace P=Moon::Species::Photochemistry;
namespace C=Moon::Species::ChargeExchange;
namespace D=Moon::Species::ReactionData;

std::string Json(const std::string& value) {
  std::ostringstream out;out<<'"';
  for (unsigned char c:value) {
    if (c=='"' || c=='\\') out<<'\\'<<char(c);
    else if (c<32) out<<"\\u00"<<std::hex<<std::setw(2)<<std::setfill('0')<<unsigned(c)<<std::dec;
    else out<<char(c);
  }
  out<<'"';return out.str();
}
double Number(const std::string& text) {
  std::size_t n=0;const double value=std::stod(text,&n);
  if (n!=text.size()) throw std::runtime_error("invalid complete numeric token");
  return value;
}
void Scalar(bool ok,double value,const std::string& error) {
  std::cout<<"{\"ok\":"<<(ok ? "true" : "false")<<",\"value\":"<<value
           <<",\"error\":"<<Json(error)<<"}\n";
}

int main(int argc,char** argv) {
  std::cout<<std::setprecision(17);
  if (argc!=3) return 2;
  std::string error;
  if (std::string(argv[1])=="sha") {
    std::string result="unchanged";const bool ok=D::FileSha256(argv[2],&result,&error);
    std::cout<<"{\"ok\":"<<(ok ? "true" : "false")<<",\"digest\":"<<Json(result)<<",\"error\":"<<Json(error)<<"}\n";
    return ok ? 0 : 2;
  }
  P::Registry photo;C::Registry capture;
  const bool isPhoto=std::string(argv[1])=="photo";
  if (!isPhoto && std::string(argv[1])!="cex") return 2;
  const auto photoSource=P::FrozenV1Sources(argv[2]);
  const auto captureSource=C::FrozenV1References(argv[2]);
  if (!(isPhoto ? photo.LoadFrozenD04(photoSource,&error) : capture.LoadPublishedReferences(captureSource,&error))) {
    std::cout<<"{\"initial_load_error\":"<<Json(error)<<"}\n";return 2;
  }
  std::string line;
  while (std::getline(std::cin,line)) {
    try {
      std::istringstream input(line);std::string command;input>>command;
      double value=17.125;bool ok=false;
      if (isPhoto) {
        if (command=="inventory") {
          std::cout<<"{\"catalogue\":"<<photo.catalogueSize()<<",\"reactive\":"<<photo.reactiveSize()<<",\"channels\":[";
          bool first=true;
          for (const auto& key:{"He","Na","Ne","Ar","H2O","OH_exp"}) {
            std::vector<P::Channel> channels;
            if (!photo.ChannelsForParent(key,&channels,&error)) return 2;
            // Nonreactive total rows are separately queried: they remain in
            // the catalogue but can never appear in an event/branch list.
            if (std::string(key)=="H2O" || std::string(key)=="OH_exp") {
              P::Channel total;
              if (!photo.GetChannel(std::string("PHIDRATES:")+key+":00:total",&total,&error)) return 2;
              channels.push_back(total);
            }
            for (const auto& c:channels) {
              if (!first) std::cout<<',';
              first=false;
              std::cout<<"{\"id\":"<<Json(c.id)<<",\"key\":"<<Json(c.speciesKey)
                <<",\"parent\":"<<Json(c.physicalParent)<<",\"parent_state\":"<<Json(c.parentState)
                <<",\"reactive\":"<<(c.reactive ? "true" : "false")
                <<",\"experimental_oh\":"<<(c.experimentalOh ? "true" : "false")
                <<",\"source_products\":"<<Json(c.sourceProducts)<<",\"source_label\":"<<Json(c.sourceChannelLabel)
                <<",\"source_reference\":"<<Json(c.sourceReference)<<",\"source_file\":"<<Json(c.sourceFile)
                <<",\"rates\":["<<c.rate1AuS[0]<<','<<c.rate1AuS[1]<<"],\"mean_energy_eV\":["
                <<c.meanExcessEnergyEv[0]<<','<<c.meanExcessEnergyEv[1]<<"],\"has_threshold\":"
                <<(c.hasThreshold ? "true" : "false")<<",\"threshold_eV\":"<<c.thresholdEnergyEv
                <<",\"sample_count\":"<<c.tabulated.sampleCount<<",\"wavelength_nm\":["
                <<c.tabulated.minimumWavelengthNm<<','<<c.tabulated.maximumWavelengthNm
                <<"],\"photon_energy_eV\":["<<c.tabulated.minimumPhotonEnergyEv<<','<<c.tabulated.maximumPhotonEnergyEv
                <<"],\"daughter_kinetics\":"<<(c.daughterKineticsAvailable ? "true" : "false")<<",\"products\":[";
              bool firstProduct=true;
              for (const auto& p:c.products) {
                if (!firstProduct) std::cout<<',';
                firstProduct=false;
                std::cout<<"{\"species\":"<<Json(p.species)<<",\"state\":"<<Json(p.electronicState)
                  <<",\"multiplicity\":"<<p.multiplicity<<",\"charge\":"<<p.charge<<",\"nuclei\":{";
                bool firstNucleus=true;
                for (const auto& item:p.nuclei) {
                  if (!firstNucleus) std::cout<<',';
                  firstNucleus=false;
                  std::cout<<Json(item.first)<<':'<<item.second;
                }
                std::cout<<"}}";
              }
              std::cout<<"]}";
            }
          }
          std::cout<<"]}\n";continue;
        }
        if (command=="rate") {
          std::string key,distance,visibility;int solar;
          input>>key>>solar>>distance>>visibility;
          ok=photo.TotalDestructionRateS(key,static_cast<P::SolarCase>(solar),Number(distance),Number(visibility),&value,&error);
        }
        else if (command=="channel-rate") {
          std::string id;int solar;input>>id>>solar;
          ok=photo.ChannelRate1AuS(id,static_cast<P::SolarCase>(solar),&value,&error);
        }
        else if (command=="branches") {
          std::string key;int solar;input>>key>>solar;
          std::vector<std::pair<std::string,double>> branches;
          ok=photo.BranchFractions(key,static_cast<P::SolarCase>(solar),&branches,&error);
          std::cout<<"{\"ok\":"<<(ok ? "true" : "false")<<",\"branches\":[";
          bool first=true;for (const auto& item:branches) {
            if (!first) std::cout<<',';
            first=false;
            std::cout<<'['<<Json(item.first)<<','<<item.second<<']';
          }
          std::cout<<"],\"error\":"<<Json(error)<<"}\n";continue;
        }
        else if (command=="select") {
          std::string key,u,id="unchanged";int solar;input>>key>>solar>>u;
          ok=photo.SelectChannel(key,static_cast<P::SolarCase>(solar),Number(u),&id,&error);
          std::cout<<"{\"ok\":"<<(ok ? "true" : "false")<<",\"id\":"<<Json(id)<<",\"error\":"<<Json(error)<<"}\n";continue;
        }
        else if (command=="grid") {
          std::string key;int solar,n;input>>key>>solar>>n;
          if (n<=0) return 2;
          std::map<std::string,int> counts;
          for (int i=0;i<n;++i) {
            std::string id;
            if (!photo.SelectChannel(key,static_cast<P::SolarCase>(solar),(i+0.5)/n,&id,&error)) return 2;
            ++counts[id];
          }
          std::cout<<"{\"counts\":{";bool first=true;
          for (const auto& c:counts) {if (!first) std::cout<<',';first=false;std::cout<<Json(c.first)<<':'<<c.second;}
          std::cout<<"}}\n";continue;
        }
        else if (command=="reload") {
          std::string field,path,hash,schema;input>>field>>std::quoted(path)>>hash>>schema;
          auto source=photoSource;
          if (field=="channels") {source.channelsPath=path;if (hash!="-") source.channelsSha256=hash;}
          else if (field=="rates") {source.ratesPath=path;if (hash!="-") source.ratesSha256=hash;}
          else if (field=="thresholds") {source.thresholdsPath=path;if (hash!="-") source.thresholdsSha256=hash;}
          else if (field=="cross") {source.crossSectionsPath=path;if (hash!="-") source.crossSectionsSha256=hash;}
          else return 2;
          if (schema!="-") source.schema=schema;
          ok=photo.LoadFrozenD04(source,&error);const std::string loadError=error;
          if (!photo.TotalDestructionRateS("Na",P::SolarCase::Quiet,1,1,&value,&error)) return 2;
          error=loadError;
        }
        else return 2;
      }
      else {
        if (command=="query" || command=="beam") {
          int id;std::string energy,speed,density;input>>id>>energy;
          if (command=="query") ok=capture.CrossSectionM2(static_cast<C::ChannelId>(id),Number(energy),&value,&error);
          else {input>>speed>>density;ok=capture.ColdBeamRateS(static_cast<C::ChannelId>(id),Number(energy),Number(speed),Number(density),&value,&error);}
        }
        else if (command=="metadata") {
          int id;input>>id;C::ChannelMetadata m;
          ok=capture.GetMetadata(static_cast<C::ChannelId>(id),&m,&error);
          std::cout<<"{\"ok\":"<<(ok ? "true" : "false")<<",\"id\":"<<Json(m.id)
            <<",\"law\":"<<Json(m.law)<<",\"energy_convention\":"<<Json(m.energyConvention)
            <<",\"domain\":["<<m.minimumEnergyEvPerU<<','<<m.maximumEnergyEvPerU
            <<"],\"source_sha256\":"<<Json(m.sourceSha256)<<",\"daughter_kinetics\":"
            <<(m.daughterKineticsAvailable ? "true" : "false")<<",\"particle_execution\":"
            <<(m.productionParticleExecutionEnabled ? "true" : "false")<<",\"error\":"<<Json(error)<<"}\n";continue;
        }
        else if (command=="knots") {
          std::vector<double> energies,values;
          if (!capture.SodiumSourceKnots(&energies,&values,&error)) return 2;
          std::cout<<"{\"knots\":[";
          for (std::size_t i=0;i<energies.size();++i) {if (i) std::cout<<',';std::cout<<'['<<energies[i]<<','<<values[i]<<']';}
          std::cout<<"]}\n";continue;
        }
        else if (command=="reload") {
          std::string field,path,hash,schema;input>>field>>std::quoted(path)>>hash>>schema;
          auto source=captureSource;
          if (field=="he") {source.heliumPublicationPath=path;if (hash!="-") source.heliumPublicationSha256=hash;}
          else if (field=="na") {source.sodiumPublicationPath=path;if (hash!="-") source.sodiumPublicationSha256=hash;}
          else return 2;
          if (schema!="-") source.schema=schema;
          ok=capture.LoadPublishedReferences(source,&error);const auto loadError=error;
          if (!capture.CrossSectionM2(C::ChannelId::ProtonNaMorgan1985LogLog,1000,&value,&error)) return 2;
          error=loadError;
        }
        else return 2;
      }
      Scalar(ok,value,error);
    }
    catch (const std::exception& e) {std::cerr<<e.what()<<'\n';return 2;}
  }
  return std::cin.bad() ? 2 : 0;
}
