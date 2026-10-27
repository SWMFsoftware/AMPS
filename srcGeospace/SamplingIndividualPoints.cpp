
#include <iostream>
#include <string>

#include "pic.h"
#include "string_func.h"
#include "Earth.h"



#include "T96Interface.h"
#include "T05Interface.h"


string Geospace::IndividualPointSample::fname="";
int Geospace::IndividualPointSample::nSampleIntervals=0;
vector <Geospace::IndividualPointSample::cIndividualPointSampleElement> Geospace::IndividualPointSample::SamplePointTable;

double Geospace::IndividualPointSample::Emin;
double Geospace::IndividualPointSample::Emax;
double Geospace::IndividualPointSample::dE;
double Geospace::IndividualPointSample::LogEmin,Geospace::IndividualPointSample::LogEmax,Geospace::IndividualPointSample::dLogE;
int Geospace::IndividualPointSample::SamplingMode=Geospace::IndividualPointSample::ModeLinear;



//reading the datapoint file 
void Geospace::IndividualPointSample::Parser() {
  ifstream fin(fname);
  double x[3];
  string str,s;
  int nLocadedPoints=0;
 
  cIndividualPointSampleElement t;

  while (getline(fin,str)){
    FindAndReplaceAll(str,"("," ");
    FindAndReplaceAll(str,")"," ");
    FindAndReplaceAll(str,","," ");
    FindAndReplaceAll(str,";"," ");

    trim(str);
    
    CutFirstWord(s,str);
    t.x[0]=stod(s); 

    CutFirstWord(s,str);
    t.x[1]=stod(s);

    CutFirstWord(s,str);
    t.x[2]=stod(s);


    SamplePointTable.push_back(t); 
    nLocadedPoints++;    
  } 

  //init the points 
  for (auto& p : SamplePointTable) {
    p.EnergyDistribution.init(nSampleIntervals);
  }
} 

