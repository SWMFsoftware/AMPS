/* Command-line assembly and exact-weight folding for Stage-3 response cells. */

#include "ResponseMatrix.hh"
#include "StrictJson.hh"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void Usage(){
  std::cerr<<"usage:\n"
           <<"  shieldSimResponseTool assemble OUTPUT.json CELL.json [CELL.json ...]\n"
           <<"  shieldSimResponseTool fold MATRIX.json WEIGHTS.json OUTPUT.json\n";
}
}

int main(int argc,char** argv){
  try {
    if(argc>=4 && std::string(argv[1])=="assemble"){
      std::vector<std::string> cells;
      for(int i=3;i<argc;++i) cells.emplace_back(argv[i]);
      auto matrix=shield_sim::ResponseMatrix::Assemble(cells);
      matrix.Write(argv[2]);
      std::cout<<"response matrix written: "<<argv[2]<<"\n";
      return 0;
    }
    if(argc==5 && std::string(argv[1])=="fold"){
      auto matrix=shield_sim::ResponseMatrix::Load(argv[2]);
      auto folded=matrix.Fold(shield_sim::Json::ParseFile(argv[3]));
      std::ofstream output(argv[4],std::ios::binary|std::ios::trunc);
      if(!output) throw std::runtime_error("cannot write fold output");
      output<<folded.DumpCanonical()<<'\n';
      std::cout<<"response fold written: "<<argv[4]<<"\n";
      return 0;
    }
    Usage();
    return 2;
  } catch(const std::exception& error){
    std::cerr<<"shieldSimResponseTool: "<<error.what()<<'\n';
    return 1;
  }
}
