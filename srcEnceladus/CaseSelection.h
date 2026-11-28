#ifndef ENCELADUS_CASE_SELECTION_H
#define ENCELADUS_CASE_SELECTION_H
#include "ForwardModel.h"

namespace Enceladus { namespace ForwardModel {
// Frontend syntax is dependency-free. Physics/source interpretation remains
// with the existing owned parsers, not with command-line string matching.
struct RuntimeSelection {
  std::string inputFile="input/enceladus.in",caseId;
  bool explicitInput=false,help=false,list=false,describe=false;
  bool validateLocal=false,validateMPI=false;
};
Status ParseRuntimeArguments(int argc,char** argv,RuntimeSelection*);
Status ParseCaseIdentity(const std::string&,std::string* measurement,std::string* basis);
std::string RuntimeUsage(const std::string& program);
} }
#endif
