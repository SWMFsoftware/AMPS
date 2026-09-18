
#include "sep.h"
#include "adapters/swcme1d_adapter.h"

namespace {

// The legacy parser passes a block as strings only. Keep the matching source
// line numbers beside that block so canonical D02 diagnostics can identify the
// exact bad assignment without changing PIC's public parser API.
std::vector<std::size_t> gBlockLineNumbers;
std::string gInputFileName;

std::string TrimCopy(const std::string& input) {
  const std::string whitespace=" \t\r\n";
  const std::size_t begin=input.find_first_not_of(whitespace);
  if (begin==std::string::npos) return std::string();
  const std::size_t end=input.find_last_not_of(whitespace);
  return input.substr(begin,end-begin+1);
}

}  // namespace

void SEP::Parser::Scattering(vector<string>& StringVector) {
  string sub,s;
  int i;

  for (i=0;i<StringVector.size();i++) {
    s=StringVector[i];

    PIC::Parser::replace(s,"="," ");
    PIC::Parser::replace(s,"au","149598000.0E3");
  
    istringstream iss(s);
    iss >> sub;

    if (sub=="lambda0") {
      iss >> sub;
      SEP::Scattering::Tenishev2005AIAA::lambda0=PIC::Parser::Evaluate(sub);
    }
    else if (sub=="alpha") {
      iss >> sub;
      SEP::Scattering::Tenishev2005AIAA::alpha=PIC::Parser::Evaluate(sub); 
    }
    else if (sub=="beta") {
      iss >> sub;
      SEP::Scattering::Tenishev2005AIAA::beta=PIC::Parser::Evaluate(sub);
    }
    else {
      exit(__LINE__,__FILE__,"Error: unknown keyword");
    }
  }

  SEP::Scattering::Tenishev2005AIAA::status=SEP::Scattering::Tenishev2005AIAA::_enabled;
}

void SEP::Parser::SWCME1D(vector<string>& StringVector) {
  // Do not interpret values here. This application parser only preserves the
  // user's key, value, file, and line. Unit conversion, unknown/duplicate-key
  // rejection, preset layering, and physical validation are all owned by the
  // canonical src/models/swcme resolver.
  for (std::size_t i=0;i<StringVector.size();++i) {
    const std::string& row=StringVector[i];
    std::size_t separator=row.find('=');
    if (separator==std::string::npos)
      separator=row.find_first_of(" \t");
    if (separator==std::string::npos) {
      const std::string message="SWCME1D assignment requires key = value: "+row;
      exit(__LINE__,__FILE__,message.c_str());
    }
    SEP::SW1DAdapter::ParameterAssignment assignment;
    assignment.key=TrimCopy(row.substr(0,separator));
    assignment.value=TrimCopy(row.substr(separator+1));
    assignment.origin=gInputFileName;
    assignment.line=i+1<gBlockLineNumbers.size() ? gBlockLineNumbers[i+1] : 0;
    const SEP::SW1DAdapter::Status status=
        SEP::SW1DAdapter::StageInputAssignment(assignment);
    if (!status.ok()) exit(__LINE__,__FILE__,status.detail.c_str());
  }
}

void SEP::Parser::SelectCommand(vector<string>& StringVector) {
  string sub,s=StringVector[0];
  StringVector.erase(StringVector.begin());

  PIC::Parser::replace(s,"="," ");
  istringstream iss(s);

  iss >> sub;

  if (sub=="Scattering") {
    iss >> sub;

    if (sub=="on") {
      Scattering(StringVector);
    }
  }
  else if (sub=="SWCME1D" || sub=="SWCME") {
    // The optional word "on" is accepted for symmetry with Scattering. An
    // omitted word means enabled because the presence of the block is already
    // an explicit configuration action.
    std::string toggle;
    iss >> toggle;
    if (toggle.empty() || toggle=="on") SWCME1D(StringVector);
    else exit(__LINE__,__FILE__,"Error: SWCME1D block accepts only optional 'on'");
  }
  else {
    exit(__LINE__,__FILE__,"Error: unknown keyword");
  }

  StringVector.clear();
}


void SEP::Parser::ReadFile(string fname) {
  string str;
  vector<string> StringVector;
  std::size_t lineNumber=0;
  ifstream file (fname); //file just has some sentence


  if (!file) {
    exit(__LINE__,__FILE__,"Error: cannot open input file");
  }

  gInputFileName=fname;
  gBlockLineNumbers.clear();
  SEP::SW1DAdapter::ClearStagedInputAssignments();


  while (getline (file,str)) {
    ++lineNumber;
    //remove comments
    str=str.substr(0,str.find("!",0));
    PIC::Parser::replace(str,"\\"," ");
    PIC::Parser::trim(str);
    
    if (str=="") {
      //the input of the command is completed
      if (StringVector.size()!=0) {
        SelectCommand(StringVector);
        gBlockLineNumbers.clear();
      }
    }
    else {
      //the input of the command is not completed yet
      StringVector.push_back(str);
      gBlockLineNumbers.push_back(lineNumber);
      str.clear();
    }
  }
 
  if (!StringVector.empty()) {
    SelectCommand(StringVector);
    gBlockLineNumbers.clear();
  }
}
