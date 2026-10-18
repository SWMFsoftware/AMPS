//parser of the "Key = Value" input file
//  - '#' starts a comment
//  - a key can be repeated (e.g. Kernel); all values are kept in the order of appearance
//  - keys are case sensitive
#ifndef INPUT_H
#define INPUT_H

#include <map>
#include <string>
#include <vector>

#include "Vec3.h"

class cInput {
public:
  void Read(const char* fname);

  bool Has(const std::string& key) const {return Table.count(key)!=0;}

  //return the last value of the key; exit with an error message when the key is missing and no default is given
  std::string GetString(const std::string& key) const;
  std::string GetString(const std::string& key,const std::string& def) const;
  double GetDouble(const std::string& key) const;
  double GetDouble(const std::string& key,double def) const;
  int GetInt(const std::string& key,int def) const;
  bool GetBool(const std::string& key,bool def) const;
  Vec3 GetVec3(const std::string& key,const Vec3& def) const;

  //return the value that must be one of 'options'
  std::string GetOption(const std::string& key,const std::vector<std::string>& options,const std::string& def) const;

  const std::vector<std::string>& GetAll(const std::string& key) const;

  //list of keys that were never requested (used to catch typos in the input file)
  std::vector<std::string> UnusedKeys() const;

  //Every value actually requested by the running model, including defaults.
  //Repeated keys retain their ordered values.  This is written to the run
  //manifest so a copied input file and the effective runtime configuration
  //remain separately auditable.
  const std::map<std::string,std::vector<std::string> >& EffectiveValues() const {return Effective;}

private:
  std::map<std::string,std::vector<std::string> > Table;
  mutable std::map<std::string,bool> Used;
  mutable std::map<std::string,std::vector<std::string> > Effective;
};

#endif
