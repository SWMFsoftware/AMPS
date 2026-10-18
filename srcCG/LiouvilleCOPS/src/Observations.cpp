#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>

#include "Observations.h"

//remove C and C++ comments
static std::string StripComments(const std::string& s) {
  std::string res;
  size_t i=0;

  while (i<s.size()) {
    if ((s[i]=='/')&&(i+1<s.size())&&(s[i+1]=='*')) {
      size_t e=s.find("*/",i+2);
      i=(e==std::string::npos) ? s.size() : e+2;
    }
    else if ((s[i]=='/')&&(i+1<s.size())&&(s[i+1]=='/')) {
      size_t e=s.find('\n',i);
      i=(e==std::string::npos) ? s.size() : e;
    }
    else if (s[i]=='"') {
      size_t e=s.find('"',i+1);
      if (e==std::string::npos) e=s.size()-1;
      res+=s.substr(i,e-i+1);
      i=e+1;
    }
    else res+=s[i++];
  }

  return res;
}

//return the text between the braces of the initializer of the array 'name'
static std::string ArrayBody(const std::string& text,const std::string& name) {
  size_t p=0;

  while ((p=text.find(name,p))!=std::string::npos) {
    size_t q=p+name.size();

    //the name must be followed by '[' (array declaration)
    while ((q<text.size())&&isspace(text[q])) q++;

    if ((q<text.size())&&(text[q]=='[')) {
      size_t eq=text.find('=',q),b=text.find('{',q);

      if ((eq!=std::string::npos)&&(b!=std::string::npos)&&(eq<b)) {
        size_t e=text.find('}',b);
        return text.substr(b+1,e-b-1);
      }
    }

    p=q;
  }

  fprintf(stderr,"Error: array '%s' is not found in the observation file\n",name.c_str());
  exit(EXIT_FAILURE);
}

static std::vector<double> ParseNumbers(const std::string& body) {
  std::vector<double> res;
  std::stringstream ss(body);
  std::string item;

  while (std::getline(ss,item,',')) {
    size_t b=item.find_first_not_of(" \t\r\n");
    if (b==std::string::npos) continue;

    char *end;
    double v=strtod(item.c_str()+b,&end);

    while ((*end!=0)&&isspace(*end)) end++;

    if (*end!=0) {
      fprintf(stderr,"Error: cannot parse the number '%s' in the observation file\n",item.c_str());
      exit(EXIT_FAILURE);
    }

    res.push_back(v);
  }

  return res;
}

void ReadObservations(const std::string& fname,std::vector<cObservation>& obs) {
  std::ifstream fin(fname);

  if (!fin) {
    fprintf(stderr,"Error: cannot open the observation file '%s'\n",fname.c_str());
    exit(EXIT_FAILURE);
  }

  std::stringstream buffer;
  buffer << fin.rdbuf();

  std::string text=StripComments(buffer.str());
  std::string tbody=ArrayBody(text,"ObservationTime");
  std::vector<std::string> times;

  for (size_t p=0;(p=tbody.find('"',p))!=std::string::npos;) {
    size_t e=tbody.find('"',p+1);
    times.push_back(tbody.substr(p+1,e-p-1));
    p=e+1;
  }

  std::vector<double> ng=ParseNumbers(ArrayBody(text,"NudeGaugeReferenceData"));
  std::vector<double> rg=ParseNumbers(ArrayBody(text,"RamGaugeReferenceData"));

  if ((times.size()!=ng.size())||(times.size()!=rg.size())) {
    fprintf(stderr,"Error: inconsistent lengths in the observation file: time=%zu, NG=%zu, RG=%zu\n",times.size(),ng.size(),rg.size());
    exit(EXIT_FAILURE);
  }

  obs.resize(times.size());

  for (size_t i=0;i<times.size();i++) {
    obs[i].UTC=times[i];
    obs[i].NG=ng[i];
    obs[i].RG=rg[i];
  }

  printf("Observations '%s': %zu points, %s .. %s\n",fname.c_str(),obs.size(),obs.front().UTC.c_str(),obs.back().UTC.c_str());
}
