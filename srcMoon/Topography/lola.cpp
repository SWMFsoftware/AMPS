//read tophological map derived from LRO/LOLA 
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string.h>
#include <list>
#include <utility>
#include <map>
#include <cmath>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <signal.h>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <iostream>
#include <string>
#include <algorithm>

using namespace std;

string fNameLON="ldem_4_float_LON.txt";
string fNameLAT="ldem_4_float_LAT.txt";
string fNameALT="ldem_4_float.txt";

vector<double> vLON,vLAT;
vector<vector<double> >vALT;

void ReadFileLON_LAT(string fName,vector<double>& vRes) {
  ifstream f(fName);
  string line;
  double  number;
  
  getline(f,line);
  std::istringstream iss(line);

  // Extract each float from the line
  while (iss >> number) { 
    vRes.push_back(number);
  }
 
  f.close();
}

void ReadFileALT(string fName,vector<vector<double> >& vRes) {
  ifstream f(fName);
  string line;
  double  number;

  for (int i=0;i<vRes.size();i++) {
    getline(f,line);
    std::istringstream iss(line);

    // Extract each float from the line
    while (iss >> number) {
      vRes[i].push_back(number);
    }
  }

  f.close();
}

void OutputTecplotMap() {
  FILE* fout=fopen("topologycam-map.dat","w");
  int iLon,iLonMax,iLat,iLatMax;

  iLonMax=vLON.size();
  iLatMax=vLAT.size();

  fprintf(fout,"\nZONE I=%ld, J=%ld, DATAPACKING=POINT\n",iLonMax,iLatMax); 

  for (iLon=0;iLon<iLonMax;iLon++) {
    for (iLat=0;iLat<iLatMax;iLat++) {
      fprintf(fout,"%e  %e  %e\n", vLON[iLon],vLAT[iLat],vALT[iLat][iLon]);
    }
  }

  fclose(fout);
}


int main () {
  

  //read the map	
  ReadFileLON_LAT(fNameLON,vLON);
  ReadFileLON_LAT(fNameLAT,vLAT);


  vALT.resize(vLAT.size());
  ReadFileALT(fNameALT,vALT);

  OutputTecplotMap();

  return 0;
}
