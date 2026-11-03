#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <fnmatch.h>
#include <unistd.h>

#include "Input.h"

static std::string Trim(const std::string& s) {
  size_t b=s.find_first_not_of(" \t\r\n");
  if (b==std::string::npos) return "";

  size_t e=s.find_last_not_of(" \t\r\n");
  return s.substr(b,e-b+1);
}

static std::string Lower(std::string s) {
  std::transform(s.begin(),s.end(),s.begin(),[](unsigned char c) {return (char)std::tolower(c);});
  return s;
}

std::vector<std::string> cInput::MachineNames() {
  std::vector<std::string> res;
  auto add=[&res](std::string n) {
    n=Lower(Trim(n));
    if (!n.empty() && std::find(res.begin(),res.end(),n)==res.end()) res.push_back(n);
  };
  char host[256]={0};

  //$COPS_MACHINE replaces the detected names (to choose a machine's lines explicitly)
  if ((getenv("COPS_MACHINE")!=NULL)&&(Trim(getenv("COPS_MACHINE"))!="")) {
    add(getenv("COPS_MACHINE"));
    return res;
  }

  if (gethostname(host,sizeof(host)-1)==0) {
    std::string h(host);
    add(h.substr(0,h.find('.')));
    add(h);
  }

  if (getenv("SLURM_CLUSTER_NAME")!=NULL) add(getenv("SLURM_CLUSTER_NAME"));
  return res;
}

static std::string JoinNames(const std::vector<std::string>& v) {
  std::string res;
  for (auto& s : v) res+=(res.empty() ? "" : ", ")+s;
  return res;
}

void cInput::Read(const char* fname) {
  std::ifstream fin(fname);
  std::string line;
  int nline=0;
  const std::vector<std::string> names=MachineNames();
  std::map<std::string,std::vector<std::string> > generic,specific;
  std::map<std::string,std::string> patterns;   //machine patterns of the qualified lines, per key

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

    //machine qualifier: the last '@' of the line, after a blank
    size_t at=value.rfind('@');

    if ((at!=std::string::npos)&&(at>0)&&isspace((unsigned char)value[at-1])) {
      std::string list=value.substr(at+1);
      value=Trim(value.substr(0,at));
      for (auto& ch : list) if (ch==',') ch=' ';

      std::istringstream ss(list);
      std::string p;
      bool any=false,match=false;

      while (ss >> p) {
        any=true;
        patterns[key]+=(patterns[key].empty() ? "" : ", ")+p;
        for (auto& n : names) if (fnmatch(Lower(p).c_str(),n.c_str(),0)==0) match=true;
      }

      if (!any) {
        fprintf(stderr,"Error: line %i of '%s' has '@' without a machine name: %s\n",nline,fname,line.c_str());
        exit(EXIT_FAILURE);
      }

      if (match) specific[key].push_back(value);
      continue;
    }

    generic[key].push_back(value);
  }

  Table=generic;
  for (auto& s : specific) Table[s.first]=s.second;
  for (auto& p : patterns) if ((specific.count(p.first)==0)&&(generic.count(p.first)==0)) OtherMachines[p.first]=p.second;

  if (!patterns.empty()) {
    std::string applied;
    for (auto& s : specific) applied+=(applied.empty() ? "" : ", ")+s.first;
    printf("Input: machine %s; machine-specific values used for: %s\n",JoinNames(names).c_str(),
        applied.empty() ? "(none)" : applied.c_str());
  }
}

std::string cInput::GetString(const std::string& key) const {
  auto it=Table.find(key);

  if (it==Table.end()) {
    auto om=OtherMachines.find(key);

    if (om!=OtherMachines.end()) {
      fprintf(stderr,"Error: the required key '%s' is given only for the machines %s; this machine is %s "
          "(set COPS_MACHINE to choose one)\n",key.c_str(),om->second.c_str(),JoinNames(MachineNames()).c_str());
    }
    else fprintf(stderr,"Error: the required key '%s' is not found in the input file\n",key.c_str());

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
