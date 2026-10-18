//ROSINA/COPS reference observations
//The data are read from the C++ header used by srcCG (RosinaMeasurements.h): arrays
//ObservationTime[], NudeGaugeReferenceData[], RamGaugeReferenceData[]
#ifndef OBSERVATIONS_H
#define OBSERVATIONS_H

#include <string>
#include <vector>

struct cObservation {
  std::string UTC;   //UTC at the spacecraft
  double NG,RG;      //reference gauge values as stored in the data file (a value <=0 marks a data gap)
};

void ReadObservations(const std::string& fname,std::vector<cObservation>& obs);

#endif
