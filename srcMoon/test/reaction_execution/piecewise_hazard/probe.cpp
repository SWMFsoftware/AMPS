// Thin numerical I/O adapter: every hazard integral/root is production code.
// The textual up:/down: notation creates adjacent native floating controls;
// it is fixture syntax only and is not a new application/data-file parser.
#include "srcMoon/species/PiecewiseHazard.h"
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
namespace H=Moon::Species::PiecewiseHazard;
using Real=long double;
Real Read() {
  std::string token;if(!(std::cin>>token))throw std::runtime_error("truncated fixture number");
  const bool up=token.rfind("up:",0)==0,down=token.rfind("down:",0)==0;
  const std::string number=token.substr(up?3:down?5:0);char* end=nullptr;
  Real value=std::strtold(number.c_str(),&end);
  if(!end||*end||end==number.c_str())throw std::runtime_error("malformed fixture number");
  if(up||down)value=std::nextafter(value,up?std::numeric_limits<Real>::infinity():-std::numeric_limits<Real>::infinity());
  return value;
}
bool Same(Real a,Real b){return a==b||(std::isnan(a)&&std::isnan(b));}
bool Same(const std::vector<H::Segment>& a,const std::vector<H::Segment>& b) {
  if(a.size()!=b.size())return false;
  for(std::size_t i=0;i<a.size();++i) {
    if(!Same(a[i].beginS,b[i].beginS)||!Same(a[i].endS,b[i].endS)||a[i].channels.size()!=b[i].channels.size())return false;
    for(std::size_t k=0;k<a[i].channels.size();++k) {
      const auto& x=a[i].channels[k];const auto& y=b[i].channels[k];
      if(x.reaction!=y.reaction||!Same(x.beginRateS,y.beginRateS)||!Same(x.endRateS,y.endRateS))return false;
    }
  }
  return true;
}
const char* Name(H::Status s) {
  switch(s) {
    case H::Status::Complete:return "Complete";
    case H::Status::Event:return "Event";
    case H::Status::IncompleteSupport:return "IncompleteSupport";
    case H::Status::InvalidInput:return "InvalidInput";
    case H::Status::NumericalError:return "NumericalError";
  }
  return "UNKNOWN";
}
int main() {
  try {
    std::cout<<std::setprecision(std::numeric_limits<Real>::max_digits10);
    std::cout<<"PLATFORM\t"<<std::numeric_limits<Real>::digits<<'\t'<<std::numeric_limits<Real>::epsilon()
        <<'\t'<<std::numeric_limits<Real>::max_exponent<<'\t'<<std::numeric_limits<Real>::min_exponent<<'\n';
    std::size_t count=0;if(!(std::cin>>count))throw std::runtime_error("missing fixture count");
    for(std::size_t i=0;i<count;++i) {
      std::string id,operation;if(!(std::cin>>id>>operation))throw std::runtime_error("truncated fixture identity");
      const Real begin=Read(),end=Read(),depth=Read();std::size_t maximum=0,n=0;
      if(!(std::cin>>maximum>>n))throw std::runtime_error("truncated fixture controls");
      std::vector<H::Segment> history;
      for(std::size_t s=0;s<n;++s) {
        H::Segment segment;segment.beginS=Read();segment.endS=Read();std::size_t channels=0;
        if(!(std::cin>>channels))throw std::runtime_error("truncated channel count");
        for(std::size_t k=0;k<channels;++k) {
          H::Channel channel;if(!(std::cin>>channel.reaction))throw std::runtime_error("truncated channel identity");
          channel.beginRateS=Read();channel.endRateS=Read();segment.channels.push_back(channel);
        }
        history.push_back(segment);
      }
      const auto before=history;std::string error;H::Status status;bool unchanged=false;
      H::Exposure exposure;exposure.integratedHazard=-13;exposure.channels={{999,-17}};
      H::EventLocation location;location.eventTimeS=-7;location.segmentIndex=999;
      location.integratedHazard=-19;location.totalInstantaneousRateS=-11;location.instantaneousRates={{999,-23}};
      if(operation=="integrate") {
        status=H::Integrate(history,begin,end,&exposure,&error);
        unchanged=exposure.integratedHazard==-13&&exposure.channels.size()==1&&
            exposure.channels[0].reaction==999&&exposure.channels[0].integratedHazard==-17;
      }
      else if(operation=="locate") {
        status=H::LocateEvent(history,begin,end,depth,H::NumericControls{maximum},&location,&error);
        unchanged=location.eventTimeS==-7&&location.segmentIndex==999&&location.integratedHazard==-19&&
            location.totalInstantaneousRateS==-11&&location.instantaneousRates.size()==1&&
            location.instantaneousRates[0].reaction==999&&location.instantaneousRates[0].rateS==-23;
      }
      else throw std::runtime_error("unknown fixture operation");
      std::cout<<"INPUT\t"<<id<<'\t'<<begin<<'\t'<<end<<'\t'<<depth<<'\n';
      std::cout<<"INPUTHEX\t"<<id<<'\t'<<std::hexfloat<<begin<<'\t'<<end<<'\t'<<depth<<std::defaultfloat<<'\n';
      std::cout<<"RESULT\t"<<id<<'\t'<<Name(status)<<'\t'<<unchanged<<'\t'<<Same(history,before)
          <<'\t'<<(operation=="integrate"?exposure.integratedHazard:location.integratedHazard)
          <<'\t'<<location.eventTimeS<<'\t'<<location.segmentIndex<<'\t'<<location.totalInstantaneousRateS<<'\n';
      std::cout<<"EVENTHEX\t"<<id<<'\t'<<std::hexfloat<<location.eventTimeS<<'\t'
          <<std::nextafter(location.eventTimeS,-std::numeric_limits<Real>::infinity())<<std::defaultfloat<<'\n';
      if(operation=="integrate")for(const auto& channel:exposure.channels)
        std::cout<<"EXPOSURE\t"<<id<<'\t'<<channel.reaction<<'\t'<<channel.integratedHazard<<'\n';
      else for(const auto& channel:location.instantaneousRates)
        std::cout<<"RATE\t"<<id<<'\t'<<channel.reaction<<'\t'<<channel.rateS<<'\n';
      std::cout<<"DIAGNOSTIC\t"<<id<<'\t'<<error<<'\n';
    }
    std::string extra;if(std::cin>>extra)throw std::runtime_error("trailing fixture tokens");
    return 0;
  }
  catch(const std::exception& error) {std::cerr<<"PROTOCOL ERROR: "<<error.what()<<'\n';return 2;}
}
