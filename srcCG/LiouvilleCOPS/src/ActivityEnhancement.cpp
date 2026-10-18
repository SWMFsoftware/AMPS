#include "ActivityEnhancement.h"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

#include "ActivityEnhancementBasis.h"

namespace ActivityEnhancement {
namespace {

std::string Trim(const std::string& value) {
  const std::size_t first=value.find_first_not_of(" \t\r\n");
  if (first==std::string::npos) return "";
  const std::size_t last=value.find_last_not_of(" \t\r\n");
  return value.substr(first,last-first+1);
}

[[noreturn]] void FileError(const std::string& filename,int line,
                            const std::string& message) {
  std::ostringstream text;
  text << "invalid activity-enhancement file '" << filename << "'";
  if (line>0) text << " at line " << line;
  text << ": " << message;
  throw std::runtime_error(text.str());
}

double ParseFiniteDouble(const std::string& text,const std::string& filename,
                         int line,const std::string& field) {
  const std::string value=Trim(text);
  char* end=nullptr;
  errno=0;
  const double result=std::strtod(value.c_str(),&end);
  if (value.empty()||end==value.c_str()||*end!='\0'||errno==ERANGE||!std::isfinite(result))
    FileError(filename,line,field+" must be one finite floating-point value");
  return result;
}

int ParseInteger(const std::string& text,const std::string& filename,int line,
                 const std::string& field) {
  const double value=ParseFiniteDouble(text,filename,line,field);
  if (std::floor(value)!=value||value<-2147483648.0||value>2147483647.0)
    FileError(filename,line,field+" must be an integer");
  return static_cast<int>(value);
}

std::string Required(const std::map<std::string,std::pair<std::string,int> >& fields,
                     const std::string& key,const std::string& filename) {
  const auto item=fields.find(key);
  if (item==fields.end()) FileError(filename,0,"missing required field "+key);
  if (item->second.first.empty()) FileError(filename,item->second.second,key+" must not be empty");
  return item->second.first;
}

void RequireEqual(const std::string& actual,const std::string& expected,
                  const std::string& key,const std::string& filename,int line) {
  if (actual!=expected)
    FileError(filename,line,key+" is '"+actual+"'; expected '"+expected+"'");
}

} // anonymous namespace

Configuration ReadFile(const std::string& filename) {
  std::ifstream input(filename);
  if (!input) throw std::runtime_error("cannot open activity-enhancement file '"+filename+"'");

  // Singular metadata keys are stored with their source line for precise
  // diagnostics.  Coefficient is the sole repeatable key.
  const std::set<std::string> allowed={
      "Format","FormatVersion","BasisID","Species","BodyFrame",
      "AngleConvention","Ordering","Clipping","MinimumDegree",
      "MaximumDegree","Scale","ScaleInterpretation",
      "LegacyActivityIdentity","Regularization","Lambda","Coefficient"};
  std::map<std::string,std::pair<std::string,int> > fields;
  std::map<std::pair<int,int>,std::pair<double,int> > coefficients;
  std::string line;
  int lineNumber=0;

  while (std::getline(input,line)) {
    ++lineNumber;
    const std::size_t comment=line.find('#');
    if (comment!=std::string::npos) line=line.substr(0,comment);
    line=Trim(line);
    if (line.empty()) continue;
    const std::size_t equals=line.find('=');
    if (equals==std::string::npos)
      FileError(filename,lineNumber,"expected 'Key = Value'");
    const std::string key=Trim(line.substr(0,equals));
    const std::string value=Trim(line.substr(equals+1));
    if (!allowed.count(key)) FileError(filename,lineNumber,"unknown field "+key);

    if (key=="Coefficient") {
      std::istringstream values(value);
      int degree,order;
      std::string coefficientText,extra;
      if (!(values>>degree>>order>>coefficientText)||(values>>extra))
        FileError(filename,lineNumber,"Coefficient must be '<degree> <order> <value>'");
      const double coefficient=ParseFiniteDouble(coefficientText,filename,lineNumber,"coefficient");
      const std::pair<int,int> mode={degree,order};
      if (coefficients.count(mode))
        FileError(filename,lineNumber,"duplicate coefficient for (l,m)=("+
            std::to_string(degree)+","+std::to_string(order)+")");
      coefficients[mode]={coefficient,lineNumber};
    }
    else {
      if (fields.count(key)) FileError(filename,lineNumber,"duplicate field "+key);
      fields[key]={value,lineNumber};
    }
  }
  if (!input.eof()) throw std::runtime_error("read error in activity-enhancement file '"+filename+"'");

  Configuration result;
  result.FilePath=filename;
  result.Model=ModelId;
  result.Format=Required(fields,"Format",filename);
  result.FormatVersion=ParseInteger(Required(fields,"FormatVersion",filename),filename,
      fields.at("FormatVersion").second,"FormatVersion");
  result.BasisId=Required(fields,"BasisID",filename);
  result.Species=Required(fields,"Species",filename);
  result.BodyFrame=Required(fields,"BodyFrame",filename);
  result.AngleDefinition=Required(fields,"AngleConvention",filename);
  result.CoefficientOrdering=Required(fields,"Ordering",filename);
  result.Clipping=Required(fields,"Clipping",filename);
  result.MinimumDegree=ParseInteger(Required(fields,"MinimumDegree",filename),filename,
      fields.at("MinimumDegree").second,"MinimumDegree");
  result.MaximumDegree=ParseInteger(Required(fields,"MaximumDegree",filename),filename,
      fields.at("MaximumDegree").second,"MaximumDegree");
  result.Scale=ParseFiniteDouble(Required(fields,"Scale",filename),filename,
      fields.at("Scale").second,"Scale");
  result.ScaleMeaning=Required(fields,"ScaleInterpretation",filename);
  result.LegacyActivityIdentity=Required(fields,"LegacyActivityIdentity",filename);
  result.Regularization=Required(fields,"Regularization",filename);
  result.Lambda=ParseFiniteDouble(Required(fields,"Lambda",filename),filename,
      fields.at("Lambda").second,"Lambda");

  RequireEqual(result.Format,FileFormatId,"Format",filename,fields.at("Format").second);
  if (result.FormatVersion!=FileFormatVersion)
    FileError(filename,fields.at("FormatVersion").second,"unsupported FormatVersion "+
        std::to_string(result.FormatVersion));
  RequireEqual(result.BasisId,ActivityEnhancementBasis::BasisId,"BasisID",filename,
      fields.at("BasisID").second);
  RequireEqual(result.Species,"H2O","Species",filename,fields.at("Species").second);
  RequireEqual(result.AngleDefinition,AngleConvention,"AngleConvention",filename,
      fields.at("AngleConvention").second);
  RequireEqual(result.CoefficientOrdering,Ordering,"Ordering",filename,
      fields.at("Ordering").second);
  RequireEqual(result.Clipping,ClippingRule,"Clipping",filename,
      fields.at("Clipping").second);
  RequireEqual(result.ScaleMeaning,ScaleInterpretation,"ScaleInterpretation",filename,
      fields.at("ScaleInterpretation").second);
  if (result.MinimumDegree!=5)
    FileError(filename,fields.at("MinimumDegree").second,"MinimumDegree must be 5");
  if (result.MaximumDegree<4||result.MaximumDegree>8)
    FileError(filename,fields.at("MaximumDegree").second,"MaximumDegree must be in [4,8]");
  if (!(result.Scale>0.0))
    FileError(filename,fields.at("Scale").second,"Scale must be finite and positive");
  if (result.Lambda<0.0)
    FileError(filename,fields.at("Lambda").second,"Lambda must be finite and nonnegative");

  if (result.MaximumDegree>=5) {
    result.Coefficients.reserve(ActivityEnhancementBasis::CoefficientCount(5,result.MaximumDegree));
    for (int degree=5;degree<=result.MaximumDegree;degree++) {
      for (int order=-degree;order<=degree;order++) {
        const auto item=coefficients.find({degree,order});
        if (item==coefficients.end())
          FileError(filename,0,"missing coefficient for (l,m)=("+
              std::to_string(degree)+","+std::to_string(order)+")");
        result.Coefficients.push_back(item->second.first);
      }
    }
  }
  const std::size_t expected=result.Coefficients.size();
  if (coefficients.size()!=expected) {
    for (const auto& item : coefficients) {
      const int degree=item.first.first,order=item.first.second;
      if (degree<5||degree>result.MaximumDegree||order<-degree||order>degree)
        FileError(filename,item.second.second,"unexpected coefficient for (l,m)=("+
            std::to_string(degree)+","+std::to_string(order)+")");
    }
    FileError(filename,0,"coefficient count does not match the declared degree range");
  }
  return result;
}

void ValidateCompatibility(const Configuration& configuration,
                           const std::string& bodyFrame,
                           const std::string& legacyActivityIdentity) {
  if (!configuration.Enabled()) return;
  if (bodyFrame.empty()||configuration.BodyFrame!=bodyFrame)
    throw std::runtime_error("enhancement BodyFrame '"+configuration.BodyFrame+
        "' does not match active BodyFrame '"+bodyFrame+"'");
  if (configuration.LegacyActivityIdentity!=legacyActivityIdentity)
    throw std::runtime_error("enhancement LegacyActivityIdentity '"+
        configuration.LegacyActivityIdentity+"' does not match active baseline '"+
        legacyActivityIdentity+"'");
}

double EvaluateResidual(const Configuration& configuration,const Vec3& position) {
  if (!configuration.Enabled()) return 0.0;
  if (!configuration.HasResidualModes()) return 0.0;
  const std::size_t count=ActivityEnhancementBasis::CoefficientCount(5,configuration.MaximumDegree);
  if (configuration.Coefficients.size()!=count)
    throw std::runtime_error("activity-enhancement coefficient vector has an inconsistent size");

  // The maximum supported range contains 56 values.  Stack storage avoids a
  // heap allocation for every one of the O(100,000) mesh facets.
  std::array<double,56> basis;
  ActivityEnhancementBasis::EvaluateCartesian(position.x,position.y,position.z,
      5,configuration.MaximumDegree,basis.data(),count);
  double sum=0.0;
  for (std::size_t i=0;i<count;i++) sum+=configuration.Coefficients[i]*basis[i];
  return configuration.Scale*sum;
}

double Apply(const Configuration& configuration,const Vec3& position,
             double legacyLinear,double* residual,double* totalLinear) {
  const double delta=EvaluateResidual(configuration,position);
  const double total=legacyLinear+delta;
  if (residual!=nullptr) *residual=delta;
  if (totalLinear!=nullptr) *totalLinear=total;
  return std::max(0.0,total);
}

} // namespace ActivityEnhancement
