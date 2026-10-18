#include "ResponseMatrix.hh"

#include "ContentHash.hh"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace shield_sim {
namespace {

double Number(const Json::Object& object,const std::string& key,const std::string& context){
  const auto found=object.find(key);
  if(found==object.end()) throw JsonError(context+" missing '"+key+"'");
  const double value=found->second.AsNumber(context+"."+key);
  if(!std::isfinite(value)) throw JsonError(context+"."+key+" must be finite");
  return value;
}

std::string String(const Json::Object& object,const std::string& key,const std::string& context){
  const auto found=object.find(key);
  if(found==object.end()) throw JsonError(context+" missing '"+key+"'");
  const std::string value=found->second.AsString(context+"."+key);
  if(value.empty()) throw JsonError(context+"."+key+" must not be empty");
  return value;
}

int Integer(const Json::Object& object,const std::string& key,const std::string& context){
  const double value=Number(object,key,context);
  if(std::floor(value)!=value) throw JsonError(context+"."+key+" must be an integer");
  return static_cast<int>(value);
}

ResponseCell ParseCell(const Json& document,const std::string& context,
                       std::string& responseKey,Json& responseKeyDocument){
  const auto& root=document.AsObject(context);
  RejectUnknownFields(root,{"schema_version","absorber_id","case_hash","histories",
      "incident","normalization","random_seed","response_key","response_key_document",
      "score_covariance","scores","secondary_production_by_layer","validation_status"},context);
  if(String(root,"schema_version",context)!="shieldSim.response_cell.v1")
    throw JsonError(context+" is not shieldSim.response_cell.v1");
  const std::string key=String(root,"response_key",context);
  const Json& keyDocument=root.at("response_key_document");
  if(Sha256(keyDocument.DumpCanonical())!=key)
    throw JsonError(context+" response_key does not hash response_key_document");
  if(responseKey.empty()){
    responseKey=key;
    responseKeyDocument=keyDocument;
  } else if(key!=responseKey ||
            keyDocument.DumpCanonical()!=responseKeyDocument.DumpCanonical()){
    throw JsonError(context+" response key differs; incompatible cells cannot share a matrix");
  }
  const auto& incident=root.at("incident").AsObject(context+".incident");
  RejectUnknownFields(incident,{"A","Z","energy_coordinate","energy_coordinate_unit",
      "energy_coordinate_value","incident_id","mu","particle","phi_rad",
      "total_kinetic_energy_MeV"},context+".incident");
  ResponseCell cell;
  cell.incidentId=String(incident,"incident_id",context+".incident");
  cell.particle=String(incident,"particle",context+".incident");
  cell.atomicNumber=Integer(incident,"Z",context+".incident");
  cell.massNumber=Integer(incident,"A",context+".incident");
  const std::string energyCoordinate=String(incident,"energy_coordinate",context+".incident");
  const std::string energyUnit=String(incident,"energy_coordinate_unit",context+".incident");
  const double energyValue=Number(incident,"energy_coordinate_value",context+".incident");
  cell.totalKineticEnergyMeV=Number(incident,"total_kinetic_energy_MeV",context+".incident");
  cell.mu=Number(incident,"mu",context+".incident");
  cell.phiRad=Number(incident,"phi_rad",context+".incident");
  const double cellHistories=Number(root,"histories",context);
  if(cellHistories<=0.0 || std::floor(cellHistories)!=cellHistories)
    throw JsonError(context+" histories must be a positive integer");
  if(cell.atomicNumber<=0 || cell.massNumber<cell.atomicNumber ||
     cell.totalKineticEnergyMeV<=0.0 || cell.mu<=0.0 || cell.mu>1.0 ||
     cell.phiRad<0.0 || cell.phiRad>=2.0*3.14159265358979323846)
    throw JsonError(context+" has an invalid species/energy/direction coordinate");
  if((cell.particle=="proton" && (cell.atomicNumber!=1 || cell.massNumber!=1)) ||
     (cell.particle=="alpha" && (cell.atomicNumber!=2 || cell.massNumber!=4)) ||
     (cell.particle!="proton" && cell.particle!="alpha" && cell.particle!="ion"))
    throw JsonError(context+" has an inconsistent particle name/Z/A");
  double resolvedEnergy=0.0;
  if(energyCoordinate=="total_kinetic_energy"){
    if(energyUnit=="MeV") resolvedEnergy=energyValue;
    else if(energyUnit=="GeV") resolvedEnergy=1000.0*energyValue;
    else throw JsonError(context+" has an invalid total-energy unit");
  } else if(energyCoordinate=="kinetic_energy_per_nucleon"){
    if(energyUnit=="MeV/n") resolvedEnergy=cell.massNumber*energyValue;
    else if(energyUnit=="GeV/n") resolvedEnergy=1000.0*cell.massNumber*energyValue;
    else throw JsonError(context+" has an invalid per-nucleon energy unit");
  } else throw JsonError(context+" has an invalid energy coordinate");
  if(energyValue<=0.0 || std::abs(resolvedEnergy-cell.totalKineticEnergyMeV)>
     1.0e-12*std::max({1.0,std::abs(resolvedEnergy),std::abs(cell.totalKineticEnergyMeV)}))
    throw JsonError(context+" energy coordinate does not resolve to total_kinetic_energy_MeV");
  std::map<std::string,double> scoreSums;
  for(const auto& entry:root.at("scores").AsArray(context+".scores")){
    const auto& score=entry.AsObject(context+".scores[]");
    RejectUnknownFields(score,{"histories","mean_per_primary","sample_variance","score_id",
        "standard_error","sum","sum_squares","unit"},context+".scores[]");
    ResponseScore parsed;
    parsed.id=String(score,"score_id",context+".scores[]");
    parsed.unit=String(score,"unit",context+".scores[]");
    parsed.meanPerPrimary=Number(score,"mean_per_primary",context+".scores[]");
    parsed.standardError=Number(score,"standard_error",context+".scores[]");
    parsed.histories=Number(score,"histories",context+".scores[]");
    const double variance=Number(score,"sample_variance",context+".scores[]");
    const double sum=Number(score,"sum",context+".scores[]");
    const double sumSquares=Number(score,"sum_squares",context+".scores[]");
    if(parsed.standardError<0.0 || parsed.histories<=0.0 || parsed.histories!=cellHistories)
      throw JsonError(context+" score uncertainty/history is invalid");
    const double reconstructedMean=sum/parsed.histories;
    const double reconstructedVariance=parsed.histories>1.0
        ? std::max(0.0,(sumSquares-sum*sum/parsed.histories)/(parsed.histories-1.0)) : 0.0;
    const double reconstructedError=std::sqrt(reconstructedVariance/parsed.histories);
    auto inconsistent=[](double a,double b){
      return std::abs(a-b)>1.0e-12*std::max({1.0,std::abs(a),std::abs(b)});
    };
    if(variance<0.0 || sumSquares<0.0 || inconsistent(parsed.meanPerPrimary,reconstructedMean) ||
       inconsistent(variance,reconstructedVariance) ||
       inconsistent(parsed.standardError,reconstructedError))
      throw JsonError(context+" score moments are not internally consistent");
    scoreSums[parsed.id]=sum;
    if(!cell.scores.emplace(parsed.id,std::move(parsed)).second)
      throw JsonError(context+" has duplicate score_id");
  }
  for(const auto& entry:root.at("score_covariance").AsArray(context+".score_covariance")){
    const auto& covariance=entry.AsObject(context+".score_covariance[]");
    RejectUnknownFields(covariance,{"covariance_of_mean","histories","sample_covariance",
        "score_id_i","score_id_j","sum_products","unit"},context+".score_covariance[]");
    const std::string first=String(covariance,"score_id_i",context+".score_covariance[]");
    const std::string second=String(covariance,"score_id_j",context+".score_covariance[]");
    if(first>second || !cell.scores.count(first) || !cell.scores.count(second))
      throw JsonError(context+" covariance pair is not a canonical pair of known scores");
    const double histories=Number(covariance,"histories",context+".score_covariance[]");
    const double sample=Number(covariance,"sample_covariance",context+".score_covariance[]");
    const double covarianceOfMean=Number(covariance,"covariance_of_mean",
                                          context+".score_covariance[]");
    const double sumProducts=Number(covariance,"sum_products",context+".score_covariance[]");
    (void)String(covariance,"unit",context+".score_covariance[]");
    const double reconstructedSample=histories>1.0
        ? (sumProducts-scoreSums.at(first)*scoreSums.at(second)/histories)/(histories-1.0)
        : 0.0;
    if(histories!=cellHistories ||
       std::abs(sample-reconstructedSample)>
         1.0e-12*std::max({1.0,std::abs(sample),std::abs(reconstructedSample)}) ||
       std::abs(covarianceOfMean-sample/histories)>
         1.0e-12*std::max({1.0,std::abs(covarianceOfMean),std::abs(sample/histories)}))
      throw JsonError(context+" covariance history/mean scaling is inconsistent");
    if(!cell.meanCovariance.emplace(std::make_pair(first,second),covarianceOfMean).second)
      throw JsonError(context+" has a duplicate covariance pair");
  }
  const std::size_t expectedPairs=cell.scores.size()*(cell.scores.size()+1)/2;
  if(cell.meanCovariance.size()!=expectedPairs)
    throw JsonError(context+" score_covariance must contain the complete upper triangle");
  for(const auto& score:cell.scores){
    const double diagonal=cell.meanCovariance.at({score.first,score.first});
    if(std::abs(diagonal-score.second.standardError*score.second.standardError)>
       1.0e-12*std::max({1.0,std::abs(diagonal),
                         score.second.standardError*score.second.standardError}))
      throw JsonError(context+" covariance diagonal differs from score standard_error squared");
  }
  return cell;
}

Json CellToJson(const ResponseCell& cell){
  Json::Array scores;
  for(const auto& item:cell.scores)
    scores.emplace_back(Json::Object{{"histories",item.second.histories},
      {"mean_per_primary",item.second.meanPerPrimary},{"score_id",item.second.id},
      {"standard_error",item.second.standardError},{"unit",item.second.unit}});
  Json::Array covariance;
  for(const auto& item:cell.meanCovariance)
    covariance.emplace_back(Json::Object{{"covariance_of_mean",item.second},
      {"score_id_i",item.first.first},{"score_id_j",item.first.second}});
  return Json(Json::Object{
    {"incident",Json(Json::Object{{"A",cell.massNumber},{"Z",cell.atomicNumber},
      {"incident_id",cell.incidentId},{"mu",cell.mu},{"particle",cell.particle},
      {"phi_rad",cell.phiRad},{"total_kinetic_energy_MeV",cell.totalKineticEnergyMeV}})},
    {"score_covariance",Json(std::move(covariance))},
    {"scores",Json(std::move(scores))}});
}

std::string EscapeCsv(const std::string& value){
  if(value.find_first_of(",\"\n") == std::string::npos) return value;
  std::string result="\"";
  for(char c:value){ if(c=='\"') result+='\"'; result+=c; }
  return result+'\"';
}

std::string HashFile(const std::filesystem::path& path){
  std::ifstream input(path,std::ios::binary);
  if(!input) throw std::runtime_error("cannot hash file '"+path.string()+"'");
  std::ostringstream bytes;
  bytes<<input.rdbuf();
  if(!input.good() && !input.eof())
    throw std::runtime_error("failed reading file for hash '"+path.string()+"'");
  return Sha256(bytes.str());
}

} // namespace

ResponseMatrix ResponseMatrix::Assemble(const std::vector<std::string>& paths){
  if(paths.empty()) throw JsonError("response-matrix assembly requires at least one cell");
  ResponseMatrix matrix;
  std::set<std::string> ids;
  std::set<std::tuple<std::string,int,int,double,double,double>> coordinates;
  for(const std::string& path:paths){
    ResponseCell cell=ParseCell(Json::ParseFile(path),path,matrix.responseKey_,
                                matrix.responseKeyDocument_);
    if(!ids.insert(cell.incidentId).second)
      throw JsonError("duplicate incident_id '"+cell.incidentId+"' during matrix assembly");
    if(!coordinates.emplace(cell.particle,cell.atomicNumber,cell.massNumber,
                            cell.totalKineticEnergyMeV,cell.mu,cell.phiRad).second)
      throw JsonError("duplicate species/energy/direction coordinate during matrix assembly");
    matrix.cells_.push_back(std::move(cell));
  }
  std::sort(matrix.cells_.begin(),matrix.cells_.end(),[](const auto& a,const auto& b){
    return a.incidentId<b.incidentId;
  });
  Json canonical=matrix.CanonicalDocument();
  matrix.matrixHash_=Sha256(canonical.DumpCanonical());
  return matrix;
}

ResponseMatrix ResponseMatrix::AssembleSelected(
    const std::vector<std::string>& paths,
    const std::vector<std::string>& scoreIds){
  if(scoreIds.empty())
    throw JsonError("selected response-matrix assembly requires at least one score_id");
  std::set<std::string> selected;
  for(const std::string& id:scoreIds){
    if(id.empty()) throw JsonError("selected response-matrix score_id must not be empty");
    if(!selected.insert(id).second)
      throw JsonError("duplicate selected response-matrix score_id '"+id+"'");
  }

  if(paths.empty()) throw JsonError("response-matrix assembly requires at least one cell");
  ResponseMatrix matrix;
  std::set<std::string> ids;
  std::set<std::tuple<std::string,int,int,double,double,double>> coordinates;
  for(const std::string& path:paths){
    // ParseCell validates the unabridged score moments, complete covariance
    // triangle, response key, species, and coordinates. Filter immediately
    // afterward so only one full cell is resident at a time; large transport
    // campaigns must not construct or hash a temporary all-score matrix.
    ResponseCell cell=ParseCell(Json::ParseFile(path),path,matrix.responseKey_,
                                matrix.responseKeyDocument_);
    if(!ids.insert(cell.incidentId).second)
      throw JsonError("duplicate incident_id '"+cell.incidentId+
                      "' during matrix assembly");
    if(!coordinates.emplace(cell.particle,cell.atomicNumber,cell.massNumber,
                            cell.totalKineticEnergyMeV,cell.mu,cell.phiRad).second)
      throw JsonError("duplicate species/energy/direction coordinate during matrix assembly");
    for(const std::string& id:selected){
      if(!cell.scores.count(id))
        throw JsonError("selected score_id '"+id+"' is absent from response cell '"+
                        cell.incidentId+"'");
    }
    for(auto score=cell.scores.begin();score!=cell.scores.end();){
      if(!selected.count(score->first)) score=cell.scores.erase(score);
      else ++score;
    }
    for(auto covariance=cell.meanCovariance.begin();covariance!=cell.meanCovariance.end();){
      if(!selected.count(covariance->first.first) ||
         !selected.count(covariance->first.second))
        covariance=cell.meanCovariance.erase(covariance);
      else ++covariance;
    }
    matrix.cells_.push_back(std::move(cell));
  }
  std::sort(matrix.cells_.begin(),matrix.cells_.end(),[](const auto& a,const auto& b){
    return a.incidentId<b.incidentId;
  });
  // Score selection changes the canonical matrix content, so weights built
  // for the unabridged matrix cannot be reused accidentally.
  matrix.matrixHash_=Sha256(matrix.CanonicalDocument().DumpCanonical());
  return matrix;
}

ResponseMatrix ResponseMatrix::Load(const std::string& path){
  const Json document=Json::ParseFile(path);
  const auto& root=document.AsObject("response matrix");
  RejectUnknownFields(root,{"schema_version","cells","fold_contract","response_key",
      "response_key_document","validation_status","matrix_hash","csv_path","csv_sha256"},
      "response matrix");
  if(String(root,"schema_version","response matrix")!="shieldSim.response_matrix.v1")
    throw JsonError("response matrix has unsupported schema_version");
  ResponseMatrix result;
  result.responseKey_=String(root,"response_key","response matrix");
  result.responseKeyDocument_=root.at("response_key_document");
  if(Sha256(result.responseKeyDocument_.DumpCanonical())!=result.responseKey_)
    throw JsonError("response matrix response_key does not hash response_key_document");
  result.matrixHash_=String(root,"matrix_hash","response matrix");
  std::string ignoredKey;
  Json ignoredDocument;
  // Matrix cells have a compact internal representation rather than complete
  // response-cell documents, so parse them directly and strictly here.
  std::set<std::string> ids;
  std::set<std::tuple<std::string,int,int,double,double,double>> coordinates;
  for(const auto& item:root.at("cells").AsArray("response matrix.cells")){
    const auto& object=item.AsObject("response matrix.cells[]");
    RejectUnknownFields(object,{"incident","score_covariance","scores"},
                        "response matrix.cells[]");
    const auto& incident=object.at("incident").AsObject("response matrix.cells[].incident");
    RejectUnknownFields(incident,{"A","Z","incident_id","mu","particle","phi_rad",
        "total_kinetic_energy_MeV"},"response matrix.cells[].incident");
    ResponseCell cell;
    cell.incidentId=String(incident,"incident_id","matrix incident");
    cell.particle=String(incident,"particle","matrix incident");
    cell.atomicNumber=Integer(incident,"Z","matrix incident");
    cell.massNumber=Integer(incident,"A","matrix incident");
    cell.totalKineticEnergyMeV=Number(incident,"total_kinetic_energy_MeV","matrix incident");
    cell.mu=Number(incident,"mu","matrix incident");
    cell.phiRad=Number(incident,"phi_rad","matrix incident");
    if(!ids.insert(cell.incidentId).second) throw JsonError("duplicate matrix incident_id");
    if(cell.atomicNumber<=0 || cell.massNumber<cell.atomicNumber ||
       cell.totalKineticEnergyMeV<=0.0 || cell.mu<=0.0 || cell.mu>1.0 ||
       cell.phiRad<0.0 || cell.phiRad>=2.0*3.14159265358979323846 ||
       !coordinates.emplace(cell.particle,cell.atomicNumber,cell.massNumber,
                            cell.totalKineticEnergyMeV,cell.mu,cell.phiRad).second)
      throw JsonError("matrix has an invalid or duplicate incident coordinate");
    if((cell.particle=="proton" && (cell.atomicNumber!=1 || cell.massNumber!=1)) ||
       (cell.particle=="alpha" && (cell.atomicNumber!=2 || cell.massNumber!=4)) ||
       (cell.particle!="proton" && cell.particle!="alpha" && cell.particle!="ion"))
      throw JsonError("matrix has an inconsistent particle name/Z/A");
    for(const auto& scoreItem:object.at("scores").AsArray("matrix scores")){
      const auto& score=scoreItem.AsObject("matrix score");
      RejectUnknownFields(score,{"histories","mean_per_primary","score_id","standard_error","unit"},
                          "matrix score");
      ResponseScore parsed;
      parsed.id=String(score,"score_id","matrix score");
      parsed.unit=String(score,"unit","matrix score");
      parsed.meanPerPrimary=Number(score,"mean_per_primary","matrix score");
      parsed.standardError=Number(score,"standard_error","matrix score");
      parsed.histories=Number(score,"histories","matrix score");
      if(parsed.standardError<0.0 || parsed.histories<=0.0 ||
         std::floor(parsed.histories)!=parsed.histories)
        throw JsonError("matrix score uncertainty/history is invalid");
      if(!cell.scores.emplace(parsed.id,std::move(parsed)).second)
        throw JsonError("matrix cell has duplicate score_id");
    }
    for(const auto& covarianceItem:object.at("score_covariance").AsArray("matrix covariance")){
      const auto& covariance=covarianceItem.AsObject("matrix covariance entry");
      RejectUnknownFields(covariance,{"covariance_of_mean","score_id_i","score_id_j"},
                          "matrix covariance entry");
      const std::string first=String(covariance,"score_id_i","matrix covariance entry");
      const std::string second=String(covariance,"score_id_j","matrix covariance entry");
      const double value=Number(covariance,"covariance_of_mean","matrix covariance entry");
      if(first>second || !cell.scores.count(first) || !cell.scores.count(second) ||
         !cell.meanCovariance.emplace(std::make_pair(first,second),value).second)
        throw JsonError("matrix has an invalid or duplicate covariance pair");
    }
    if(cell.meanCovariance.size()!=cell.scores.size()*(cell.scores.size()+1)/2)
      throw JsonError("matrix covariance must contain the complete upper triangle");
    for(const auto& score:cell.scores){
      const double diagonal=cell.meanCovariance.at({score.first,score.first});
      if(std::abs(diagonal-score.second.standardError*score.second.standardError)>
         1.0e-12*std::max({1.0,std::abs(diagonal),
                           score.second.standardError*score.second.standardError}))
        throw JsonError("matrix covariance diagonal differs from standard_error squared");
    }
    result.cells_.push_back(std::move(cell));
  }
  const std::string declared=result.matrixHash_;
  result.matrixHash_.clear();
  const std::string calculated=Sha256(result.CanonicalDocument().DumpCanonical());
  result.matrixHash_=declared;
  if(calculated!=declared) throw JsonError("response matrix content hash mismatch");
  return result;
}

Json ResponseMatrix::CanonicalDocument() const {
  Json::Array cells;
  for(const auto& cell:cells_) cells.emplace_back(CellToJson(cell));
  Json::Object root{
    {"cells",Json(std::move(cells))},
    {"fold_contract",Json(Json::Object{
      {"angular_interpolation","declared_in_response_weights; validation default bilinear_mu_phi"},
      {"energy_interpolation","declared_in_response_weights; validation default linear_response_in_log_energy"},
      {"extrapolation","prohibited"},
      {"nonlinear_outputs","fold primitive numerator/denominator separately"},
      {"normalization","per generated primary"}})},
    {"response_key",responseKey_},{"response_key_document",responseKeyDocument_},
    {"schema_version","shieldSim.response_matrix.v1"},
    {"validation_status","implemented_stage3_transport_raw_not_physically_validated"}
  };
  return Json(std::move(root));
}

void ResponseMatrix::Write(const std::string& path) const {
  namespace fs=std::filesystem;
  fs::path manifest(path);
  fs::path csv=manifest;
  csv.replace_extension(".csv");
  std::ofstream table(csv,std::ios::binary|std::ios::trunc);
  if(!table) throw std::runtime_error("cannot write response matrix CSV '"+csv.string()+"'");
  table<<"incident_id,particle,Z,A,total_kinetic_energy_MeV,mu,phi_rad,score_id,unit,mean_per_primary,standard_error,histories\n";
  table<<std::setprecision(17);
  for(const auto& cell:cells_) for(const auto& score:cell.scores)
    table<<EscapeCsv(cell.incidentId)<<','<<EscapeCsv(cell.particle)<<','
         <<cell.atomicNumber<<','<<cell.massNumber<<','<<cell.totalKineticEnergyMeV<<','
         <<cell.mu<<','<<cell.phiRad<<','<<EscapeCsv(score.first)<<','
         <<EscapeCsv(score.second.unit)<<','<<score.second.meanPerPrimary<<','
         <<score.second.standardError<<','<<score.second.histories<<'\n';
  table.close();
  if(!table) throw std::runtime_error("failed writing response matrix CSV '"+csv.string()+"'");
  Json::Object root=CanonicalDocument().AsObject();
  root["csv_path"]=csv.filename().string();
  root["csv_sha256"]=HashFile(csv);
  root["matrix_hash"]=matrixHash_.empty()?Sha256(CanonicalDocument().DumpCanonical()):matrixHash_;
  std::ofstream output(manifest,std::ios::binary|std::ios::trunc);
  if(!output) throw std::runtime_error("cannot write response matrix manifest '"+path+"'");
  output<<Json(std::move(root)).DumpCanonical()<<'\n';
}

Json ResponseMatrix::Fold(const Json& document) const {
  const auto& root=document.AsObject("response weights");
  if(String(root,"schema_version","response weights")!="shieldSim.response_weights.v1")
    throw JsonError("response weights have unsupported schema_version");
  RejectUnknownFields(root,{"schema_version","matrix_hash","normalization","quadrature_policy",
                             "interpolation_policy","weights","source_weight_covariance"},
                      "response weights");
  if(String(root,"matrix_hash","response weights")!=matrixHash_)
    throw JsonError("response weights name a different matrix_hash");
  const std::string interpolation=String(root,"interpolation_policy","response weights");
  if(interpolation!="exact_nodes" && interpolation!="linear_response_in_log_energy_bilinear_mu_phi")
    throw JsonError("unsupported response interpolation policy");
  const auto& weightsJson=root.at("weights").AsArray("response weights.weights");
  if(weightsJson.empty()) throw JsonError("response weights must not be empty");
  std::map<std::string,std::size_t> cellIndex;
  for(std::size_t i=0;i<cells_.size();++i) cellIndex[cells_[i].incidentId]=i;
  std::vector<double> weights(cells_.size(),0.0);
  std::set<std::string> seen;
  double total=0.0;
  for(const auto& item:weightsJson){
    const auto& entry=item.AsObject("response weight");
    RejectUnknownFields(entry,{"incident_id","weight"},"response weight");
    const std::string id=String(entry,"incident_id","response weight");
    if(!seen.insert(id).second) throw JsonError("duplicate response weight incident_id");
    const auto found=cellIndex.find(id);
    if(found==cellIndex.end())
      throw JsonError("response weight requests an absent node; extrapolation is prohibited");
    const double weight=Number(entry,"weight","response weight");
    if(weight<0.0) throw JsonError("response weights must be nonnegative");
    weights[found->second]=weight;
    total+=weight;
  }
  if(std::abs(total-1.0)>1.0e-12)
    throw JsonError("response weights must sum to one within 1e-12");

  std::vector<std::vector<double>> covariance(cells_.size(),
                                               std::vector<double>(cells_.size(),0.0));
  if(root.count("source_weight_covariance")){
    const auto& rows=root.at("source_weight_covariance").AsArray("source covariance");
    if(rows.size()!=cells_.size()) throw JsonError("source covariance dimension differs from matrix");
    for(std::size_t i=0;i<rows.size();++i){
      const auto& row=rows[i].AsArray("source covariance row");
      if(row.size()!=cells_.size()) throw JsonError("source covariance must be square");
      for(std::size_t j=0;j<row.size();++j){
        covariance[i][j]=row[j].AsNumber("source covariance element");
        if(!std::isfinite(covariance[i][j]) ||
           std::abs(covariance[i][j]-rows[j].AsArray()[i].AsNumber())>1.0e-12)
          throw JsonError("source covariance must be finite and symmetric");
      }
      if(covariance[i][i]<0.0) throw JsonError("source covariance diagonal must be nonnegative");
    }
    // Symmetry and nonnegative diagonal are necessary but not sufficient for a
    // covariance matrix.  This tolerance-aware Cholesky factorization rejects
    // matrices with a negative eigen-direction while accepting legitimate
    // positive-semidefinite matrices (including singular normalization
    // covariances whose row sums are zero).
    std::vector<std::vector<double>> lower(cells_.size(),
                                            std::vector<double>(cells_.size(),0.0));
    constexpr double covarianceTolerance=1.0e-12;
    for(std::size_t i=0;i<cells_.size();++i){
      for(std::size_t j=0;j<=i;++j){
        double residual=covariance[i][j];
        for(std::size_t k=0;k<j;++k) residual-=lower[i][k]*lower[j][k];
        if(i==j){
          if(residual < -covarianceTolerance)
            throw JsonError("source covariance must be positive semidefinite");
          lower[i][j]=std::sqrt(std::max(0.0,residual));
        } else if(lower[j][j]>covarianceTolerance){
          lower[i][j]=residual/lower[j][j];
        } else if(std::abs(residual)>covarianceTolerance){
          throw JsonError("source covariance must be positive semidefinite");
        }
      }
    }
  }

  std::set<std::string> scoreIds;
  std::map<std::string,std::string> scoreUnits;
  for(const auto& cell:cells_) for(const auto& score:cell.scores){
    scoreIds.insert(score.first);
    const auto existing=scoreUnits.find(score.first);
    if(existing!=scoreUnits.end() && existing->second!=score.second.unit)
      throw JsonError("response score unit differs between cells");
    scoreUnits[score.first]=score.second.unit;
  }
  Json::Array results;
  for(const std::string& id:scoreIds){
    std::vector<double> response(cells_.size(),0.0),se(cells_.size(),0.0);
    for(std::size_t i=0;i<cells_.size();++i){
      const auto found=cells_[i].scores.find(id);
      if(found!=cells_[i].scores.end()){
        response[i]=found->second.meanPerPrimary;
        se[i]=found->second.standardError;
      }
    }
    double mean=0.0,matrixVariance=0.0,sourceVariance=0.0;
    for(std::size_t i=0;i<cells_.size();++i){
      mean+=weights[i]*response[i];
      matrixVariance+=weights[i]*weights[i]*se[i]*se[i];
      for(std::size_t j=0;j<cells_.size();++j)
        sourceVariance+=response[i]*covariance[i][j]*response[j];
    }
    if(sourceVariance<-1.0e-12)
      throw JsonError("source covariance produces a negative propagated variance");
    sourceVariance=std::max(0.0,sourceVariance);
    results.emplace_back(Json::Object{{"combined_standard_uncertainty",
      std::sqrt(matrixVariance+sourceVariance)},{"matrix_statistical_standard_uncertainty",
      std::sqrt(matrixVariance)},{"score_id",id},{"source_standard_uncertainty",
      std::sqrt(sourceVariance)},{"unit",scoreUnits[id]},{"value",mean}});
  }
  Json::Array resultCovariance;
  for(auto first=scoreIds.begin();first!=scoreIds.end();++first){
    for(auto second=first;second!=scoreIds.end();++second){
      std::vector<double> responseFirst(cells_.size(),0.0);
      std::vector<double> responseSecond(cells_.size(),0.0);
      double matrixCovariance=0.0;
      for(std::size_t i=0;i<cells_.size();++i){
        const auto firstScore=cells_[i].scores.find(*first);
        const auto secondScore=cells_[i].scores.find(*second);
        if(firstScore!=cells_[i].scores.end()) responseFirst[i]=firstScore->second.meanPerPrimary;
        if(secondScore!=cells_[i].scores.end()) responseSecond[i]=secondScore->second.meanPerPrimary;
        const auto covarianceEntry=cells_[i].meanCovariance.find({*first,*second});
        if(covarianceEntry!=cells_[i].meanCovariance.end())
          matrixCovariance+=weights[i]*weights[i]*covarianceEntry->second;
      }
      double sourceCovariance=0.0;
      for(std::size_t i=0;i<cells_.size();++i) for(std::size_t j=0;j<cells_.size();++j)
        sourceCovariance+=responseFirst[i]*covariance[i][j]*responseSecond[j];
      resultCovariance.emplace_back(Json::Object{
        {"combined_covariance",matrixCovariance+sourceCovariance},
        {"matrix_statistical_covariance",matrixCovariance},
        {"score_id_i",*first},{"score_id_j",*second},
        {"source_covariance",sourceCovariance},
        {"unit",scoreUnits[*first]+" * "+scoreUnits[*second]}});
    }
  }
  return Json(Json::Object{{"interpolation_policy",interpolation},
    {"matrix_hash",matrixHash_},{"normalization",String(root,"normalization","response weights")},
    {"quadrature_policy",String(root,"quadrature_policy","response weights")},
    {"result_covariance",Json(std::move(resultCovariance))},
    {"results",Json(std::move(results))},{"schema_version","shieldSim.response_fold.v1"},
    {"validation_status","implemented_stage3_transport_raw_not_physically_validated"}});
}

} // namespace shield_sim
