#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>

#include "Input.h"

static std::string Trim(const std::string& s) {
  size_t b=s.find_first_not_of(" \t\r\n");
  if (b==std::string::npos) return "";

  size_t e=s.find_last_not_of(" \t\r\n");
  return s.substr(b,e-b+1);
}

void cInput::Read(const char* fname) {
  std::ifstream fin(fname);
  std::string line;
  int nline=0;

  if (!fin) {
    fprintf(stderr,"Error: cannot open the input file '%s'\n",fname);
    exit(EXIT_FAILURE);
  }

  while (std::getline(fin,line)) {
    nline++;

    size_t c=line.find('#');
    if (c!=std::string::npos) line=line.substr(0,c);
    line=Trim(line);
    if (line.empty()) continue;

    size_t eq=line.find('=');
    if (eq==std::string::npos) {
      fprintf(stderr,"Error: line %i of '%s' is not in the 'Key = Value' format: %s\n",nline,fname,line.c_str());
      exit(EXIT_FAILURE);
    }

    std::string key=Trim(line.substr(0,eq));
    std::string value=Trim(line.substr(eq+1));

    Table[key].push_back(value);
  }
}

std::string cInput::GetString(const std::string& key) const {
  auto it=Table.find(key);

  if (it==Table.end()) {
    fprintf(stderr,"Error: the required key '%s' is not found in the input file\n",key.c_str());
    exit(EXIT_FAILURE);
  }

  Used[key]=true;
  Effective[key]={it->second.back()};
  return it->second.back();
}

std::string cInput::GetString(const std::string& key,const std::string& def) const {
  if (Has(key)) return GetString(key);
  Effective[key]={def};
  return def;
}

double cInput::GetDouble(const std::string& key) const {
  std::string s=GetString(key);
  char *end;
  double res=strtod(s.c_str(),&end);

  if (*end!=0) {
    fprintf(stderr,"Error: the value of '%s' (%s) is not a number\n",key.c_str(),s.c_str());
    exit(EXIT_FAILURE);
  }

  return res;
}

double cInput::GetDouble(const std::string& key,double def) const {
  if (Has(key)) return GetDouble(key);
  std::ostringstream value;
  value << std::setprecision(17) << def;
  Effective[key]={value.str()};
  return def;
}

int cInput::GetInt(const std::string& key,int def) const {
  if (Has(key)) return (int)GetDouble(key);
  Effective[key]={std::to_string(def)};
  return def;
}

bool cInput::GetBool(const std::string& key,bool def) const {
  if (!Has(key)) {
    Effective[key]={def ? "on" : "off"};
    return def;
  }

  std::string s=GetString(key);
  std::transform(s.begin(),s.end(),s.begin(),::tolower);

  if ((s=="on")||(s=="true")||(s=="yes")||(s=="1")) return true;
  if ((s=="off")||(s=="false")||(s=="no")||(s=="0")) return false;

  fprintf(stderr,"Error: the value of '%s' (%s) must be on/off\n",key.c_str(),s.c_str());
  exit(EXIT_FAILURE);
}

Vec3 cInput::GetVec3(const std::string& key,const Vec3& def) const {
  if (!Has(key)) {
    std::ostringstream value;
    value << std::setprecision(17) << def.x << " " << def.y << " " << def.z;
    Effective[key]={value.str()};
    return def;
  }

  std::istringstream ss(GetString(key));
  Vec3 v;

  if (!(ss >> v.x >> v.y >> v.z)) {
    fprintf(stderr,"Error: the value of '%s' must be three numbers\n",key.c_str());
    exit(EXIT_FAILURE);
  }

  return v;
}

std::string cInput::GetOption(const std::string& key,const std::vector<std::string>& options,const std::string& def) const {
  std::string s=GetString(key,def);

  for (auto& o : options) if (o==s) return s;

  fprintf(stderr,"Error: the value of '%s' (%s) is not one of:",key.c_str(),s.c_str());
  for (auto& o : options) fprintf(stderr," %s",o.c_str());
  fprintf(stderr,"\n");
  exit(EXIT_FAILURE);
}

const std::vector<std::string>& cInput::GetAll(const std::string& key) const {
  static const std::vector<std::string> empty;
  auto it=Table.find(key);

  Used[key]=true;
  Effective[key]=(it==Table.end()) ? empty : it->second;
  return (it==Table.end()) ? empty : it->second;
}

std::vector<std::string> cInput::UnusedKeys() const {
  std::vector<std::string> res;

  for (auto& t : Table) if (Used.count(t.first)==0) res.push_back(t.first);
  return res;
}
