#include "CaseSelection.h"
#include <cctype>
#include <sstream>
namespace Enceladus { namespace ForwardModel {
// CLI parsing has no MPI/PIC/filesystem side effects. Keeping option and ID
// validation here lets local inspection work on login nodes and allows all
// actual MPI ranks to reach a shared error decision before mesh allocation.
Status ParseCaseIdentity(const std::string& id,std::string* measurement,std::string* basis) {
  const auto delimiter=id.find("--");
  if(delimiter==std::string::npos||delimiter==0||delimiter+2==id.size()||id.find("--",delimiter+2)!=std::string::npos)
    return Status::Error("case ID must be MEASUREMENT--SOURCE_BASIS with exactly one '--'");
  const auto left=id.substr(0,delimiter),right=id.substr(delimiter+2);
  for(const auto& token:{left,right}) {
    if(token.front()=='-'||token.back()=='-')return Status::Error("case components cannot start or end with '-'");
    for(unsigned char c:token)if(!std::isalnum(c)&&c!='_'&&c!='-')
      return Status::Error("case components permit only letters, digits, '_' and single '-'");
  }
  if(measurement)*measurement=left;if(basis)*basis=right;
  return Status::Ok();
}
std::string RuntimeUsage(const std::string& p) {
  return "Usage: "+p+" [-i FILE|-input FILE|--input FILE|--input=FILE] "
    "[-case MEASUREMENT--SOURCE_BASIS|--case=ID] "
    "[--list-cases|--describe-case|--validate-input|--validate-runtime-input]\n"
    "Default run deck: input/enceladus.in\n"
    "Case data root: Config.pl -model-data-path=PATH; no fitting-strategy option.\n";
}
Status ParseRuntimeArguments(int argc,char** argv,RuntimeSelection* output) {
  RuntimeSelection c;bool selected=false;
  for(int i=1;i<argc;++i) {
    std::string a=argv[i],v;
    if(a=="--help"||a=="-h") {if(c.help)return Status::Error("duplicate help");c.help=true;continue;}
    if(a=="-i"||a=="-input"||a=="--input"||a.rfind("--input=",0)==0) {
      if(c.explicitInput)return Status::Error("input path was specified more than once");
      if(a.rfind("--input=",0)==0)v=a.substr(8);
      else if(i+1<argc)v=argv[++i];
      if(v.empty()||v[0]=='-')return Status::Error("input option requires a file");
      c.inputFile=v;c.explicitInput=true;continue;
    }
    if(a=="-case"||a=="--case"||a.rfind("--case=",0)==0) {
      if(selected)return Status::Error("case was specified more than once");
      if(a.rfind("--case=",0)==0)v=a.substr(7);
      else if(i+1<argc)v=argv[++i];
      auto s=ParseCaseIdentity(v,nullptr,nullptr);if(!s.success)return s;
      c.caseId=v;selected=true;continue;
    }
    bool* flag=a=="--list-cases"?&c.list:a=="--describe-case"?&c.describe:
      a=="--validate-input"?&c.validateLocal:a=="--validate-runtime-input"?&c.validateMPI:nullptr;
    if(!flag)return Status::Error("unknown CLI option: "+a);
    if(*flag)return Status::Error("duplicate CLI option: "+a);*flag=true;
  }
  if(int(c.list)+int(c.describe)+int(c.validateLocal)+int(c.validateMPI)>1)
    return Status::Error("list/describe/validation actions are mutually exclusive");
  if(c.list&&!c.caseId.empty())return Status::Error("--list-cases cannot select one case");
  if(c.describe&&c.caseId.empty())return Status::Error("--describe-case requires -case ID");
  *output=c;return Status::Ok();
}
} }
